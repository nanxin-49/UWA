function receipt=closeout_r2_saved_seed(case_root,seed,mode)
%CLOSEOUT_R2_SAVED_SEED Independently audit/seal existing arrays, never solve.
arguments
    case_root (1,:) char
    seed (1,1) double {mustBeMember(seed,[260002 260003])}
    mode (1,:) char
end
assert(ismember(mode,{'numeric','figures','seal'}));
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
task=fullfile(case_root,sprintf('seed%d',seed)); assert(isfolder(task));
pre=jsondecode(fileread(fullfile(case_root,'source_preflight.json'))); assert(pre.failed_checks==0);
record=pre.sources([pre.sources.seed]==seed); source=record.source; source_work=record.work;
oldpath=path; pathguard=onCleanup(@()path(oldpath));
addpath(fullfile(source,'code'),'-begin'); addpath(fullfile(root,'scripts','reporting'),'-begin');
addpath(fullfile(root,'scripts','validation'),'-begin');
threads=maxNumCompThreads(1); tg=onCleanup(@()maxNumCompThreads(threads));
guards=fullfile(task,'solve_sentinels'); if ~isfolder(guards), mkdir(guards); end
for name={'solve_helmholtz_bie_halfplane_vertical','solve_helmholtz_bie_halfplane_cached_vertical', ...
        'closeout_s3a_pm_fullfk','run_pe_1d_surface_reflection_validation','run_r2_piston_pm_job','run_r2_pm_seed'}
    p=fullfile(guards,[name{1} '.m']);
    if ~isfile(p), fid=fopen(p,'w'); assert(fid>=0); fprintf(fid,'function varargout=%s(varargin)\nerror(''r2closeout:NoSolve'',''Saved arrays only'');\nend\n',name{1}); fclose(fid); end
end
addpath(guards,'-begin'); rehash;
for name={'solve_helmholtz_bie_halfplane_vertical','closeout_s3a_pm_fullfk','run_r2_pm_seed'}
    assert(startsWith(which(name{1}),guards));
end
assert(strcmp(r1_bie_manual_io('sha',fullfile(source,'manifest.json')),record.manifest_sha256));
r2_pm_io('verify_manifest',source);
for name={'audit_r1_bie_manual_stage','audit_piston_pm_r1','audit_r2_pm_final','summarize_r2_pm_fields'}
    assert(strcmp(r1_bie_manual_io('sha',which(name{1})),r1_bie_manual_io('sha',fullfile(source,'code',[name{1} '.m']))), ...
        'r2closeout:CodeDrift','Frozen numerical audit dependency changed');
end
runtime=fullfile(task,'runtime.json');
if ~isfile(runtime)
    [~,sys]=memory; r1_bie_manual_io('write',runtime,struct('MATLAB_version',version,'mode','existing local MCP', ...
        'RAM_baseline_bytes',sys.PhysicalMemory.Available,'reserve_fraction',.3,'single_thread',true,'new_solves_allowed',false));
end
run=r1_bie_manual_io('read',runtime); local_guard(run.RAM_baseline_bytes);
timer=tic; diary(fullfile(task,[mode '_diary.log'])); dg=onCleanup(@()diary('off'));
work=fullfile(task,'work'); ids={'bie_coarse','bie_main','bie_self','bie_space','bie_window'};
try
switch mode
    case 'numeric'
        assert(~isfolder(work)); mkdir(work);
        for name={'common.mat','comparison.mat','model_metrics.csv','receiver_fields_and_errors.csv','policy.json'}
            copyfile(fullfile(source_work,name{1}),fullfile(work,name{1}));
        end
        mkdir(fullfile(work,'input')); copyfile(fullfile(source,'input','*'),fullfile(work,'input'));
        q=load(fullfile(work,'common.mat'),'cfg','x','z','profile');
        assert(q.cfg.seed==seed && q.cfg.source.total_width_m==.1 && q.cfg.source.frequency_hz==4000);
        assert(isequal(q.x,(-205:205)'*.1953125) && isequal(q.z,3*ones(411,1)));
        assert(isequal(q.profile,readtable(fullfile(source,'input','coefficients.csv'))));
        recipes=[60 10 96;60 12 96;60 12 64;60 14 96;70 12 96]; stages=cell(5,1);
        for j=1:5
            local_guard(run.RAM_baseline_bytes); dest=fullfile(work,'jobs',ids{j}); mkdir(dest);
            copyfile(fullfile(source,'jobs',ids{j},'data.mat'),fullfile(dest,'data.mat'));
            copyfile(fullfile(source,'jobs',ids{j},'points.csv'),fullfile(dest,'points.csv'));
            assert(strcmp(r1_bie_manual_io('sha',fullfile(dest,'data.mat')),r1_bie_manual_io('sha',fullfile(source,'jobs',ids{j},'data.mat'))));
            tick=tic; stages{j}=audit_r1_bie_manual_stage(dest,fullfile(work,'common.mat'),recipes(j,:));
            stages{j}.elapsed_s=toc(tick); assert(stages{j}.failed_checks==0);
            r1_bie_manual_io('write',fullfile(task,'numeric_checkpoint.json'),struct('completed_BIE_stage',j-1,'seed',seed, ...
                'new_acoustic_solves',0,'elapsed_s',toc(timer),'audit_path',['work/jobs/' ids{j} '/audit.json']));
        end
        comprehensive=audit_piston_pm_r1(work); diary(fullfile(task,'numeric_diary.log'));
        final=audit_r2_pm_final(work); summary=summarize_r2_pm_fields({work},fullfile(work,'spatial_summary'));
        assert(comprehensive.failed_checks==0 && final.failed_checks==0);
        receipt=struct('status','FULL_INDEPENDENT_SAVED_NUMERIC_PASS','seed',seed,'failed_checks',0,'stages',{stages}, ...
            'comprehensive',comprehensive,'reference',final,'spatial_summary',summary.status);
    case 'figures'
        n=r1_bie_manual_io('read',fullfile(task,'numeric_receipt.json')); assert(n.failed_checks==0);
        original_path=fullfile(task,'original_figure_data_audit');
        if isfile(fullfile(original_path,'receipt.json'))
            original=r1_bie_manual_io('read',fullfile(original_path,'receipt.json'));
            assert(original.failed_checks==0 && original.checked_error_curves==12);
        else
            original=audit_r2_final_figures({work},fullfile(source_work,'figures_final'),original_path,'original-data');
        end
        render_r2_stable_saved({work},fullfile(work,'figures_final'));
        final=audit_r2_final_figures({work},fullfile(work,'figures_final'),fullfile(task,'figure_audit'));
        assert(original.failed_checks==0 && final.failed_checks==0);
        receipt=struct('status','FIXED_LAYOUT_AND_STRICT_FIG_PNG_PASS','seed',seed,'failed_checks',0, ...
            'original_data_audit',original,'new_figure_audit',final,'original_status_preserved',record.original_status);
    case 'seal'
        assert(~isfile(fullfile(task,'committed.json')));
        n=r1_bie_manual_io('read',fullfile(task,'numeric_receipt.json'));
        f=r1_bie_manual_io('read',fullfile(task,'figures_receipt.json'));
        assert(n.failed_checks==0 && f.failed_checks==0 && strcmp(n.status,'FULL_INDEPENDENT_SAVED_NUMERIC_PASS') && ...
            strcmp(f.status,'FIXED_LAYOUT_AND_STRICT_FIG_PNG_PASS'));
        for j=1:5
            assert(strcmp(r1_bie_manual_io('sha',fullfile(work,'jobs',ids{j},'data.mat')), ...
                r1_bie_manual_io('sha',fullfile(source,'jobs',ids{j},'data.mat'))));
        end
        for name={'common.mat','comparison.mat','model_metrics.csv','receiver_fields_and_errors.csv','policy.json'}
            assert(strcmp(r1_bie_manual_io('sha',fullfile(work,name{1})),r1_bie_manual_io('sha',fullfile(source_work,name{1}))));
        end
        code={'closeout_r2_saved_seed','render_r2_stable_saved','audit_r2_final_figures', ...
            'audit_r1_bie_manual_stage','audit_piston_pm_r1','audit_r2_pm_final','summarize_r2_pm_fields'};
        versions=struct('path',{},'sha256',{}); mkdir(fullfile(task,'code'));
        for j=1:numel(code)
            p=which(code{j}); versions(j)=struct('path',p,'sha256',r1_bie_manual_io('sha',p)); copyfile(p,fullfile(task,'code',[code{j} '.m']));
        end
        r1_bie_manual_io('write',fullfile(task,'code_versions.json'),versions);
        receipt=struct('status','FORMAL_LOCAL_SAVED_CLOSEOUT_PASS','seed',seed,'failed_checks',0, ...
            'original_remote_status',record.original_status,'original_remote_exit_code',record.original_real_exit_code, ...
            'original_final_seal_present',record.original_final_seal,'qualification','independent saved data and strict FIG/PNG closeout', ...
            'source_manifest_sha256',record.manifest_sha256,'policy_sha256',record.policy_sha256, ...
            'numeric_receipt_sha256',r1_bie_manual_io('sha',fullfile(task,'numeric_receipt.json')), ...
            'figures_receipt_sha256',r1_bie_manual_io('sha',fullfile(task,'figures_receipt.json')), ...
            'new_acoustic_solves',0,'SSH_connections',0,'next_allowed_step','Only paired saved-data R1+R2 summary; no R3/broadband authorization');
        r1_bie_manual_io('write',fullfile(task,'closeout_receipt.json'),receipt);
        r1_bie_manual_io('write',fullfile(task,'seal_plan.json'),struct('signature',[record.manifest_sha256 ':LOCAL-SAVED-CLOSEOUT'], ...
            'note','After this MATLAB function returns, external integrity tool includes completed diaries and every local artifact in committed.json'));
end
local_guard(run.RAM_baseline_bytes); receipt.elapsed_s=toc(timer); receipt.new_acoustic_solves=0; receipt.SSH_connections=0;
if ~strcmp(mode,'seal'), r1_bie_manual_io('write',fullfile(task,[mode '_receipt.json']),receipt); end
fprintf('R2_SAVED_CLOSEOUT seed=%d mode=%s status=%s elapsed_s=%.6f\n',seed,mode,receipt.status,receipt.elapsed_s);
catch ME
    r1_bie_manual_io('write',fullfile(task,[mode '_failure.json']),struct('status','FAIL','error',ME.message,'identifier',ME.identifier,'stack',ME.stack));
    fprintf(2,'%s\n',getReport(ME,'extended','hyperlinks','off')); rethrow(ME);
end
end
function local_guard(baseline)
[~,sys]=memory; assert(sys.PhysicalMemory.Available>=.3*baseline,'r2closeout:RAM','Local saved-data RAM reserve breached');
end
