function receipt=validate_k1_pm_preparation(package_root,output_dir)
%VALIDATE_K1_PM_PREPARATION Geometry/anchor/state/gate tests with solve sentinels.
arguments
    package_root (1,:) char
    output_dir (1,:) char
end
pm_sensitivity_io('absolute',package_root); pm_sensitivity_io('absolute',output_dir);
assert(~isfolder(output_dir)); mkdir(output_dir); timer=tic;
oldpath=path; pg=onCleanup(@()path(oldpath)); addpath(fullfile(package_root,'code'),'-begin');
threads=maxNumCompThreads(1); tg=onCleanup(@()maxNumCompThreads(threads));
manifest=k1_pm_io('read',fullfile(package_root,'manifest.json'));
assert(strcmp(manifest.scope,'K1_PM_MASTER'));
pm_sensitivity_io('verify',package_root,manifest.files);
guard_dir=fullfile(output_dir,'solve_sentinels'); mkdir(guard_dir);
names={'solve_helmholtz_bie_halfplane_vertical','solve_helmholtz_bie_halfplane_cached_vertical', ...
    'closeout_s3a_pm_fullfk','run_k1_piston_pm_job','baffled_strip_piston_field_vertical', ...
    'piston_independent_angular_spectrum_vertical','run_pe_1d_surface_reflection_validation'};
for j=1:numel(names)
    fid=fopen(fullfile(guard_dir,[names{j} '.m']),'w'); assert(fid>=0);
    fprintf(fid,'function varargout=%s(varargin)\nerror(''k1:NoAcoustics'',''Preparation forbids acoustic calls'');\nend\n',names{j}); fclose(fid);
end
addpath(guard_dir,'-begin'); rehash;
for j=1:numel(names), assert(startsWith(which(names{j}),guard_dir)); end
cfg=k1_pm_io('read',fullfile(package_root,'config.json')); checks=0;
assert(cfg.max_processes==1 && cfg.threads==1 && ~cfg.parpool && cfg.RAM_reserve_fraction==.3); checks=checks+1;
anchor_test=fullfile(output_dir,'anchor_recheck'); mkdir(anchor_test);
copyfile(fullfile(package_root,'anchor'),fullfile(anchor_test,'anchor'));
movefile(fullfile(anchor_test,'anchor','anchor_proof.json'),fullfile(anchor_test,'input_anchor_proof.json'));
proof=prepare_k1_pm_anchor(anchor_test);
assert(strcmp(proof.status,'STRICT_R1_K0_ANCHOR_MATCH_PASS') && proof.failed_checks==0); checks=checks+1;
for j=1:numel(cfg.cases)
    c=cfg.cases(j); root=fullfile(package_root,'cases',c.case_id);
    result=run_k1_pm_case(root);
    assert(strcmp(result.status,'DRY_RUN_INPUT_CHECKS_PASS') && result.new_acoustic_solves==0 && strcmp(result.case_id,c.case_id)); checks=checks+1;
end
% The physical formula bodies are unchanged from the validated R2 adapter.
reference=fileread(fullfile(package_root,'provenance','r2_reference_model_job.m'));
current=fileread(fullfile(package_root,'code','run_k1_piston_pm_job.m'));
assert(strcmp(extractAfter(reference,'function r=local_ps'),extractAfter(current,'function r=local_ps'))); checks=checks+1;
start='recipe=[128 .05 .1875;'; stop='sx=linspace(-50,50,2001)';
assert(strcmp(extractBetween(reference,start,stop),extractBetween(current,start,stop))); checks=checks+1;
% Synthetic checkpoint payloads exercise restart logic; they are not physics.
folder=fullfile(output_dir,'state_fixture'); mkdir(folder); sig='TEST_ONLY_INPUTS_PREDECESSOR';
assert(strcmp(k1_pm_io('resume_state',folder,sig),'SOLVE')); checks=checks+1;
fid=fopen(fullfile(folder,'points.csv'),'w'); fprintf(fid,'TEST_ONLY\n'); fclose(fid);
synthetic_payload=1; save(fullfile(folder,'data.mat'),'synthetic_payload','-v7');
k1_pm_io('write',fullfile(folder,'solve_invocation.json'),struct('invocations',1,'stage','C0'));
local_reject(@()k1_pm_io('resume_state',folder,sig)); checks=checks+1;
k1_pm_io('seal',folder,'solved.json',sig,{'data.mat','points.csv','solve_invocation.json'},'SOLVED_UNAUDITED');
assert(strcmp(k1_pm_io('resume_state',folder,sig),'AUDIT_ONLY')); checks=checks+1;
k1_pm_io('write',fullfile(folder,'audit.json'),struct('status','PASS','failed_checks',0));
k1_pm_io('write',fullfile(folder,'gates.json'),struct('status','PASS'));
k1_pm_io('seal',folder,'committed.json',sig,{'data.mat','points.csv','solved.json','solve_invocation.json','audit.json','gates.json'},'AUDITED_PASS');
assert(strcmp(k1_pm_io('resume_state',folder,sig),'SKIP_AUDITED')); checks=checks+1;
local_reject(@()k1_pm_io('resume_state',folder,'CHANGED_PREDECESSOR')); checks=checks+1;
fid=fopen(fullfile(folder,'points.csv'),'a'); fprintf(fid,'TAMPER'); fclose(fid);
local_reject(@()k1_pm_io('resume_state',folder,sig)); checks=checks+1;
interrupted=fullfile(output_dir,'unsealed_invocation'); mkdir(interrupted);
k1_pm_io('write',fullfile(interrupted,'solve_invocation.json'),struct('invocations',1));
local_reject(@()k1_pm_io('resume_state',interrupted,sig)); checks=checks+1;
% Preserve approved v2 nonmonotonic diagnostic and all strict numeric gates.
policy=k1_pm_io('read',fullfile(package_root,'anchor','policy.json'));
offset=[4e-6 2e-6 2e-6 4.4e-6 2.1e-6]; fields=arrayfun(@(v)ones(411,1)+v,offset,'UniformOutput',false);
g=r1_bie_continuation_gates(fields,(-205:205)'*.1953125,3*ones(411,1),policy);
assert(strcmp(g.status,'PASS') && g.grid_diagnostic.grid_nonmonotonic && strcmp(g.legacy_policy_status,'FAIL')); checks=checks+1;
fields{4}=fields{4}+1e-3; g=r1_bie_continuation_gates(fields,(-205:205)'*.1953125,3*ones(411,1),policy);
assert(strcmp(g.status,'FAIL')); checks=checks+1;
bad=policy; bad.gates.sensitivity=1e-3; local_reject(@()r1_bie_continuation_io('policy',bad)); checks=checks+1;
% Check unit source/profile/receiver drift rejection before any physics.
clone=fullfile(output_dir,'tampered_case'); copyfile(fullfile(package_root,'cases',cfg.cases(1).case_id),clone);
fid=fopen(fullfile(clone,'input','coefficients.csv'),'a'); fprintf(fid,' '); fclose(fid);
local_reject(@()run_k1_pm_case(clone)); checks=checks+1;
for rel={'../outside','C:/outside','/outside'}
    local_reject(@()k1_pm_io('path',folder,rel{1})); checks=checks+1;
end
pm_sensitivity_io('verify',package_root,manifest.files);
receipt=struct('status','K1_PREPARATION_TESTS_PASS','failed_checks',0,'checks',checks, ...
    'new_acoustic_solves',0,'acoustic_calls',0,'solve_sentinels',numel(names),'case_dry_runs',6, ...
    'anchor_recheck','STRICT_MATCH','physical_formula_regions_identical',true,'matlab_release',version('-release'), ...
    'elapsed_s',toc(timer),'input_manifest_sha256',k1_pm_io('sha',fullfile(package_root,'manifest.json')));
pm_sensitivity_io('write',fullfile(output_dir,'receipt.json'),receipt);
pm_sensitivity_io('log',fullfile(output_dir,'preparation.log'),jsonencode(receipt));
end

function local_reject(f)
rejected=false;
try
    f();
catch
    rejected=true;
end
assert(rejected,'Expected preparation rejection did not occur');
end
