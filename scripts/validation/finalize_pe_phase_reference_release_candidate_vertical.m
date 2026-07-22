function release = finalize_pe_phase_reference_release_candidate_vertical()
%FINALIZE_PE_PHASE_REFERENCE_RELEASE_CANDIDATE_VERTICAL Produce PASS/FAIL/INCOMPLETE.
root=fileparts(fileparts(fileparts(mfilename('fullpath')))); addpath(root);
run_meta=pe_phase_release_run_meta_vertical(root);
run_dir=fullfile(root,'results','validation','pe_phase_release_candidate',run_meta.run_id);
if ~exist(run_dir,'dir'), mkdir(run_dir); end
checks=struct('name',{},'passed',{},'value',{},'limit',{},'note',{});
missing={};
paths=local_paths(root);
names=fieldnames(paths);
for ii=1:numel(names), if ~isfile(paths.(names{ii})), missing{end+1}=paths.(names{ii}); end, end %#ok<AGROW>
if isempty(missing)
    try
        phase=load(paths.phase,'result'); assert_pe_phase_release_artifact_vertical(phase.result,run_meta.run_id,'phase validation');
        checks(end+1)=local_check('phase_reference_validation',phase.result.pass.all,double(phase.result.pass.all),1,'all phase, migration, and compact research gates');
        checks(end+1)=local_check('nominal_delay_4ms',phase.result.pass.analytic_delay,1e3*phase.result.analytic.relative_delay_error_s,1e3/4000,'F=65 error <= one delay sample');
        checks(end+1)=local_check('migration_full_closure',phase.result.pass.migration,max(struct2array(rmfield(phase.result.migration,'meta'))),1e-12,'mean, C/P, complex and augmented EVD');
        pa=load(paths.phase_audit,'audit'); assert_pe_phase_release_artifact_vertical(pa.audit,run_meta.run_id,'uniform phase audit');
        checks(end+1)=local_check('uniform_phase_audit',pa.audit.passed,double(pa.audit.passed),1,'independent angular-spectrum sign audit');
        ad=load(paths.adjoint,'validation'); assert_pe_phase_release_artifact_vertical(ad.validation,run_meta.run_id,'adjoint validation');
        checks(end+1)=local_check('adjoint_full',ad.validation.pass.all,double(ad.validation.pass.all),1,'exact adjoint, projection, dense/FFT, F9 and F64');
        checks(end+1)=local_check('f9_split_floor',ad.validation.f9.comparison.pass_split_floor,max([ad.validation.f9.comparison.C_to_split_floor_ratio,ad.validation.f9.comparison.P_Cnorm_to_split_floor_ratio]),1.25,'analytic/sample versus split floor');
        checks(end+1)=local_check('f64_split_floor',ad.validation.f64.comparison.pass_split_floor,max([ad.validation.f64.comparison.C_to_split_floor_ratio,ad.validation.f64.comparison.P_Cnorm_to_split_floor_ratio]),1.25,'analytic/sample versus split floor');
        u5=load(paths.u5,'result'); u8=load(paths.u8,'result');
        assert_pe_phase_release_artifact_vertical(u5.result,run_meta.run_id,'U5 validation');
        assert_pe_phase_release_artifact_vertical(u8.result,run_meta.run_id,'U8 validation');
        checks(end+1)=local_check('u5_full',u5.result.pass.all,double(u5.result.pass.all),1,'full conditional validation');
        checks(end+1)=local_check('u8_full',u8.result.pass.all,double(u8.result.pass.all),1,'full conditional validation');
        m5=load(paths.model5,'model'); m8=load(paths.model8,'model');
        assert_pe_phase_release_artifact_vertical(m5.model,run_meta.run_id,'U5 model');
        assert_pe_phase_release_artifact_vertical(m8.model,run_meta.run_id,'U8 model');
        models_ok=startsWith(m5.model.schema_version,'2.') && startsWith(m8.model.schema_version,'2.') && ...
            strcmp(m5.model.phase_reference_meta.target_reference,'direct_dsp') && ...
            strcmp(m8.model.phase_reference_meta.target_reference,'direct_dsp') && ...
            m5.model.reference_delay_s==0 && m8.model.reference_delay_s==0;
        checks(end+1)=local_check('conditional_model_schema',models_ok,double(models_ok),1,'schema 2.x, direct_dsp, zero common shift');
        lib=load(paths.library,'library'); assert_pe_phase_release_artifact_vertical(lib.library,run_meta.run_id,'conditional library');
        checks(end+1)=local_check('conditional_library',strcmp(lib.library.schema_version,'2.0.0-discrete'),double(strcmp(lib.library.schema_version,'2.0.0-discrete')),1,'exact U5/U8 nodes');
        pub=load(paths.public,'result'); cp=load(paths.cached_public,'result');
        assert_pe_phase_release_artifact_vertical(pub.result,run_meta.run_id,'public modes');
        assert_pe_phase_release_artifact_vertical(cp.result,run_meta.run_id,'cached/public');
        checks(end+1)=local_check('public_regression',pub.result.pass,double(pub.result.pass),1,'scalar direct-only/direct-plus-reflect and unchanged default');
        checks(end+1)=local_check('cached_public_f64',cp.result.pass,cp.result.max_abs_total_error,1e-10,'F64 double consistency');
        comm=load(paths.comm,'validation'); assert_pe_phase_release_artifact_vertical(comm.validation,run_meta.run_id,'two-node communication');
        [comm_ok,comm_max]=local_comm_check(comm.validation);
        checks(end+1)=local_check('two_node_communication',comm_ok,comm_max,1e-10,'four sources per node, direct_dsp, finite, paired latent comparison');
        atlas=load(paths.atlas,'atlas'); assert_pe_phase_release_artifact_vertical(atlas.atlas,run_meta.run_id,'PE atlas');
        atlas_ok=atlas.atlas.checks.all_pass && local_atlas_files(fileparts(paths.atlas));
        checks(end+1)=local_check('propagation_atlas',atlas_ok,double(atlas_ok),1,'13 PNG, MP4, MAT, manifest, summary; same run_id');
        [parse_ok,parse_count]=local_parse_check(root);
        checks(end+1)=local_check('matlab_static_parse',parse_ok,parse_count,0,'no checkcode parse errors in RC sources');
    catch exception
        missing{end+1}=sprintf('Evaluation error: %s',exception.message);
    end
end
if ~isempty(missing), decision='INCOMPLETE';
elseif all([checks.passed]), decision='PASS'; else, decision='FAIL'; end
phase_reference_meta=struct('target_reference','direct_dsp', ...
    'relative_delay_s',(100+3-(100-3))/1500);
schema_version='1.0.0'; validation_run_meta=run_meta;
release=struct('schema_version',schema_version,'phase_reference_meta',phase_reference_meta, ...
    'validation_run_meta',validation_run_meta,'decision',decision,'checks',checks, ...
    'missing_or_errors',{missing},'recommend_merge',strcmp(decision,'PASS'), ...
    'roles',struct('public_pe','general propagation entry', ...
    'cached_forward','fixed-path regression oracle', ...
    'adjoint_projection','exact fast receiver realization for fixed environment', ...
    'analytic_fft_cp','receiver second-order statistics without receiver Monte Carlo', ...
    'conditional_generator','large communication Monte Carlo after physical/statistical validation'));
save(fullfile(run_dir,'pe_phase_release_candidate_decision.mat'),'release', ...
    'schema_version','phase_reference_meta','validation_run_meta','-v7.3');
local_write_report(release,fullfile(run_dir,'pe_phase_release_candidate_report.md'));
local_write_report(release,fullfile(root,'reports','pe_phase_reference_release_candidate_report.md'));
fprintf('PE phase release candidate %s: %s\n',run_meta.run_id,decision);
end

function p=local_paths(r)
p=struct( ...
 'phase',fullfile(r,'results','validation','pe_channel_phase_reference','pe_channel_phase_reference_validation.mat'), ...
 'phase_audit',fullfile(r,'results','validation','pe_phase_convention_uniform','pe_phase_convention_audit.mat'), ...
 'adjoint',fullfile(r,'results','validation','adjoint_pe_receiver_projection','adjoint_pe_receiver_projection_full.mat'), ...
 'u5',fullfile(r,'results','validation','u5_conditional_channel_f64','u5_conditional_channel_validation_f64_full.mat'), ...
 'u8',fullfile(r,'results','validation','u8_conditional_channel_f64','u8_conditional_channel_validation_f64_full.mat'), ...
 'model5',fullfile(r,'results','validation','u5_conditional_channel_f64','u5_conditional_channel_model_f64_full.mat'), ...
 'model8',fullfile(r,'results','validation','u8_conditional_channel_f64','u8_conditional_channel_model_f64_full.mat'), ...
 'library',fullfile(r,'results','validation','conditional_channel_library','conditional_channel_library_u5_u8_f64.mat'), ...
 'comm',fullfile(r,'results','validation','two_node_communication','two_node_communication_validation_full.mat'), ...
 'public',fullfile(r,'results','validation','validate_public_channel_modes_vertical_result.mat'), ...
 'cached_public',fullfile(r,'results','validation','cached_joint_kstat_pe_receiver','cached_public_consistency_double.mat'), ...
 'atlas',fullfile(r,'results','visualization','pe_propagation_atlas','pe_propagation_atlas_data.mat'));
end

function c=local_check(name,pass,value,limit,note), c=struct('name',name,'passed',logical(pass),'value',double(value),'limit',double(limit),'note',note); end

function [ok,mx]=local_comm_check(v)
ok=strcmp(v.mode,'full') && v.channel_count_per_source==32 && strcmp(v.channel_phase_reference,'direct_dsp'); mx=0;
for key={'u5','u8'}
    for name={'kdomain_pe','joint_cached_pe','stats_full','stats_99_9'}
        r=v.results.(key{1}).(name{1}); q=max(r.tap_meta.ifft_closure_error);
        mx=max(mx,q); ok=ok && all(isfinite([r.BER,r.SER,r.effective_snr_db_mean])) && q<=1e-10;
    end
    ok=ok && isfield(v.rank_pair.(key{1}),'ber_difference_low_minus_full');
end
end

function ok=local_atlas_files(folder)
ok=isfile(fullfile(folder,'pe_propagation_atlas_manifest.md')) && ...
    isfile(fullfile(folder,'pe_propagation_atlas_summary.txt'));
for ii=1:13, d=dir(fullfile(folder,sprintf('%02d_*.png',ii))); ok=ok&&~isempty(d)&&max([d.bytes])>1024; end
d=dir(fullfile(folder,'14_*.mp4')); ok=ok&&~isempty(d)&&max([d.bytes])>1024;
end

function [ok,count]=local_parse_check(root)
files={'pe_phase_release_run_meta_vertical.m','assert_pe_phase_release_artifact_vertical.m', ...
    fullfile('scripts','validation','prepare_pe_phase_reference_release_candidate_vertical.m'), ...
    fullfile('scripts','validation','audit_phase_reference_artifacts_vertical.m'), ...
    fullfile('scripts','validation','finalize_pe_phase_reference_release_candidate_vertical.m'), ...
    fullfile('scripts','reporting','generate_pe_propagation_atlas_vertical.m'), ...
    fullfile('scripts','reporting','plot_pe_propagation_atlas_vertical.m')};
count=0;
for ii=1:numel(files)
    messages=checkcode(fullfile(root,files{ii}),'-id');
    if ~isempty(messages), count=count+sum(contains({messages.message},'Parse error','IgnoreCase',true)); end
end
ok=count==0;
end

function local_write_report(r,path)
parent=fileparts(path); if ~exist(parent,'dir'), mkdir(parent); end
fid=fopen(path,'w','n','UTF-8'); cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'# PE carrier-phase release candidate\n\n');
fprintf(fid,'- Decision: **%s**\n- run_id: `%s`\n- revision: `%s`\n- fingerprint: `%s`\n- created: `%s`\n\n', ...
    r.decision,r.validation_run_meta.run_id,r.validation_run_meta.code_revision, ...
    r.validation_run_meta.code_fingerprint,r.validation_run_meta.created_at);
if ~isempty(r.missing_or_errors)
    fprintf(fid,'## Missing or evaluation errors\n\n'); for ii=1:numel(r.missing_or_errors), fprintf(fid,'- %s\n',r.missing_or_errors{ii}); end
end
fprintf(fid,'\n## Acceptance checks\n\n| Check | Pass | Value | Limit | Note |\n|---|:---:|---:|---:|---|\n');
for ii=1:numel(r.checks), c=r.checks(ii); fprintf(fid,'| `%s` | %d | %.6g | %.6g | %s |\n',c.name,c.passed,c.value,c.limit,c.note); end
fprintf(fid,'\n## Architecture decision\n\n');
if strcmp(r.decision,'PASS'), fprintf(fid,'Recommend integration into the main branch.\n\n');
elseif strcmp(r.decision,'FAIL'), fprintf(fid,'Do not integrate until the failed numerical or semantic gates are corrected.\n\n');
else, fprintf(fid,'No integration recommendation: required formal artifacts are incomplete.\n\n'); end
fprintf(fid,'- Public PE: general propagation entry.\n- Cached forward: fixed-path high-confidence regression oracle.\n- Adjoint projection: exact receiver realization generation in a fixed environment.\n- Analytic FFT C/P: receiver second-order statistics without receiver Monte Carlo.\n- Conditional generator: large communication Monte Carlo after physical/statistical validation.\n');
end
