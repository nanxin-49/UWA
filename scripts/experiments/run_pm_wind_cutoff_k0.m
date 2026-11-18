function receipt = run_pm_wind_cutoff_k0(package_root, mode)
%RUN_PM_WIND_CUTOFF_K0 Default preflight; explicit 'execute' runs surface K0.
% Standalone allowlisted package. No PS, BIE, FullFK or production calls.
arguments
    package_root (1,:) char
    mode (1,:) char {mustBeMember(mode,{'preflight','execute','test-only'})} = 'preflight'
end
pm_sensitivity_io('absolute',package_root);
manifest=jsondecode(fileread(fullfile(package_root,'manifest.json')));
pm_sensitivity_io('verify',package_root,manifest.files);
binding=pm_sensitivity_io('sha',fullfile(package_root,'manifest.json'));
cfg=jsondecode(fileread(fullfile(package_root,'config.json')));
official=strcmp(cfg.kind,'K0_SURFACE_ONLY');
assert(official || strcmp(cfg.kind,'TEST_ONLY'));
assert(isequal(cfg.U_PM_mps(:),[5;6;8]) && isequal(cfg.N(:),[7;12;19;32;40]));
assert(isequal(cfg.seeds(:),(260001:260003)') && cfg.L_m==160);
assert(cfg.sample_count==8192 && strcmp(cfg.family,'legacy12_gaussian_extension_v1'));
if official
    assert(numel(cfg.cases)==45 && strcmp(cfg.wind_label,'U_PM'));
    expected=[];
    for seed=cfg.seeds(:).'
        for u=cfg.U_PM_mps(:).'
            for nn=cfg.N(:).', expected=[expected;seed u nn]; end %#ok<AGROW>
        end
    end
    actual=[[cfg.cases.seed]' [cfg.cases.U_PM_mps]' [cfg.cases.N]'];
    assert(isequal(actual,expected),'Frozen 45-case order changed');
end
for seed=cfg.seeds(:).'
    anchor=readtable(fullfile(package_root,'inputs',sprintf('seed%d.csv',seed)));
    d=build_pm_wind_cutoff_surface(anchor,seed,6,12,256);
    assert(isequal(d.coefficients.cos_depth_coeff_m,anchor.cos_coeff_m));
    assert(isequal(d.coefficients.sin_depth_coeff_m,anchor.sin_coeff_m));
end
receipt=struct('status','PREFLIGHT_PASS','case_count',numel(cfg.cases),'failed_checks',0, ...
    'manifest_sha256',binding,'acoustic_calls',0,'execution_mode','MATLAB_surface_only');
if strcmp(mode,'preflight'), disp(jsonencode(receipt)); return; end
assert((official && strcmp(mode,'execute')) || (~official && strcmp(mode,'test-only')),'Test/official mode mismatch');
if official
    assert(isunix && strcmp(version('-release'),'2025b'),'Official K0 requires Ubuntu MATLAB R2025b');
    assert(maxNumCompThreads==1,'K0 must use one computational thread');
end
out=fullfile(package_root,'output'); if ~isfolder(out), mkdir(out); end
seal=fullfile(out,'checkpoint.json');
if isfile(seal)
    previous=jsondecode(fileread(seal));
    assert(strcmp(previous.manifest_sha256,binding));
    pm_sensitivity_io('verify',out,previous.artifacts_sha256);
    for j=1:numel(cfg.cases)
        c=cfg.cases(j); id=sprintf('seed%d_U%d_N%02d',c.seed,c.U_PM_mps,c.N);
        cp=jsondecode(fileread(fullfile(out,id,'checkpoint.json')));
        assert(strcmp(cp.manifest_sha256,binding) && strcmp(cp.status,'AUDITED_PASS'));
        pm_sensitivity_io('verify',fullfile(out,id),cp.artifacts_sha256);
    end
    receipt=previous; receipt.action='VERIFIED_NO_REPEAT'; disp(jsonencode(receipt)); return
end
t0=tic; cpu0=cputime; completed=0;
for j=1:numel(cfg.cases)
    c=cfg.cases(j); id=sprintf('seed%d_U%d_N%02d',c.seed,c.U_PM_mps,c.N);
    case_dir=fullfile(out,id); committed=fullfile(case_dir,'checkpoint.json');
    if isfile(committed)
        cp=jsondecode(fileread(committed)); assert(strcmp(cp.manifest_sha256,binding));
        assert(strcmp(cp.status,'AUDITED_PASS')); pm_sensitivity_io('verify',case_dir,cp.artifacts_sha256);
        completed=completed+1; continue
    end
    if ~isfolder(case_dir), mkdir(case_dir); end
    generated=fullfile(case_dir,'generated.json');
    if ~isfile(generated)
        marker=fullfile(case_dir,'generation_started.json');
        assert(~isfile(marker),'Interrupted uncommitted generation: inspect this case; never auto-repeat');
        pm_sensitivity_io('write',marker,struct('manifest_sha256',binding,'case',c));
        generation_log=fullfile(case_dir,'generation.log');
        tc=tic;
        try
            pm_sensitivity_io('log',generation_log,sprintf('K0 %s; %s; %s',id,version,computer));
            anchor=readtable(fullfile(package_root,'inputs',sprintf('seed%d.csv',c.seed)));
            data=build_pm_wind_cutoff_surface(anchor,c.seed,c.U_PM_mps,c.N,cfg.sample_count);
            save(fullfile(case_dir,'surface.mat'),'data','-v7');
            writetable(data.coefficients,fullfile(case_dir,'coefficients.csv'));
            writetable(data.samples,fullfile(case_dir,'surface.csv'));
            writetable(data.spectra,fullfile(case_dir,'spectra.csv'));
            pm_sensitivity_io('write',fullfile(case_dir,'config.json'), ...
                struct('case',c,'run_config',cfg,'extension_seed',data.extension_seed, ...
                'wind_definition',data.wind_definition,'sign_convention',data.sign_convention));
            pm_sensitivity_io('write',fullfile(case_dir,'metrics.json'),data.metrics);
            pm_sensitivity_io('log',generation_log,sprintf('generation_elapsed_s=%.6f; acoustic_calls=0',toc(tc)));
        catch err
            pm_sensitivity_io('log',generation_log,getReport(err,'extended')); rethrow(err)
        end
        paths={'generation_started.json','generation.log','surface.mat','coefficients.csv','surface.csv','spectra.csv','config.json','metrics.json'};
        gen=struct('status','GENERATED_NOT_AUDITED','manifest_sha256',binding, ...
            'elapsed_s',toc(tc),'files',pm_sensitivity_io('manifest',case_dir,paths));
        pm_sensitivity_io('write',generated,gen);
    end
    gen=jsondecode(fileread(generated)); assert(strcmp(gen.manifest_sha256,binding));
    pm_sensitivity_io('verify',case_dir,gen.files);
    attempts=dir(fullfile(case_dir,'audit_attempt_*.log')); attempt=numel(attempts)+1;
    logname=sprintf('audit_attempt_%03d.log',attempt);
    audit_log=fullfile(case_dir,logname); ta=tic;
    pm_sensitivity_io('log',audit_log,sprintf('Independent saved-data audit %s attempt=%d',id,attempt));
    try
        audit=audit_pm_wind_cutoff_k0_case(case_dir);
        pm_sensitivity_io('log',audit_log,sprintf('audit_elapsed_s=%.6f failed_checks=%d',toc(ta),audit.failed_checks));
    catch err
        pm_sensitivity_io('log',audit_log,getReport(err,'extended'));
        pm_sensitivity_io('write',fullfile(case_dir,sprintf('failure_%03d.json',attempt)), ...
            struct('status','FAILED_STOPPED','error',err.message,'manifest_sha256',binding,'next_allowed_step','Inspect saved case; no automatic generation retry'));
        rethrow(err)
    end
    paths=[{gen.files.path},{'generated.json','audit.json',logname}];
    cp=struct('status','AUDITED_PASS','case',c,'manifest_sha256',binding, ...
        'completed_steps',{{'generation','saved_data_audit','sampling_audit'}}, ...
        'artifacts_sha256',pm_sensitivity_io('manifest',case_dir,paths), ...
        'failed_or_unrun_checks',{{}},'audit_status_and_path','audit.json', ...
        'resource_usage',struct('generation_elapsed_s',gen.elapsed_s,'audit_elapsed_s',toc(ta),'peak_RAM','see session resource.tsv'), ...
        'matlab_release',version('-release'),'execution_mode','surface_only_serial', ...
        'next_allowed_step','Next K0 case only; K1 disabled');
    pm_sensitivity_io('write',committed,cp); completed=completed+1;
end
% Rebuild summary into a unique attempt. Previously committed cases are read only.
listing=dir(fullfile(out,'summary_attempt_*')); summary=fullfile(out,sprintf('summary_attempt_%03d',numel(listing)+1));
mkdir(summary); summarylog=fullfile(summary,'summary.log');
summary_receipt=render_pm_wind_cutoff_k0(out,summary,cfg);
pm_sensitivity_io('log',summarylog,sprintf('summary failed_checks=%d',summary_receipt.failed_checks));
files=dir(fullfile(summary,'*')); names={files(~[files.isdir]).name};
% Avoid absolute paths in portable checkpoint manifests.
[~,summary_name]=fileparts(summary);
relative=cellfun(@(p)[summary_name '/' p],names,UniformOutput=false);
for j=1:numel(cfg.cases)
    c=cfg.cases(j); relative{end+1}=sprintf('seed%d_U%d_N%02d/checkpoint.json',c.seed,c.U_PM_mps,c.N); %#ok<AGROW>
end
receipt=struct('status','K0_COMPLETE_SURFACE_ONLY','case_count',completed,'manifest_sha256',binding, ...
    'completed_steps',{{'all_case_audits','summary_csv','summary_figures'}}, ...
    'artifacts_sha256',pm_sensitivity_io('manifest',out,relative), ...
    'elapsed_this_invocation_s',toc(t0),'cpu_this_invocation_s',cputime-cpu0, ...
    'failed_checks',0,'acoustic_calls',0,'matlab_release',version('-release'), ...
    'next_allowed_step','Stop; user reviews K0 before K1 authorization');
pm_sensitivity_io('write',seal,receipt); disp(jsonencode(receipt));
end
