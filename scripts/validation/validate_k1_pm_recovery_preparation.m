function receipt=validate_k1_pm_recovery_preparation(master,out)
%VALIDATE_K1_PM_RECOVERY_PREPARATION Real snapshot dry-run, no acoustic calls.
assert(~isfolder(out)); mkdir(out); t=tic;
oldpath=path; pg=onCleanup(@()path(oldpath)); addpath(fullfile(master,'code'),'-begin');
threads=maxNumCompThreads(1); tg=onCleanup(@()maxNumCompThreads(threads));
diary(fullfile(out,'diary.log')); dg=onCleanup(@()diary('off'));
m=k1_pm_io('read',fullfile(master,'manifest.json')); assert(strcmp(m.scope,'K1_PM_RECOVERY_MASTER'));
pm_sensitivity_io('verify',master,m.files);
cfg=k1_pm_io('read',fullfile(master,'config.json'));
if isfield(cfg,'codex_target_mode') && cfg.codex_target_mode
    assert(cfg.inherited_used_wall_s>=818 && cfg.max_processes==4);
else
    assert(cfg.inherited_used_wall_s==724 && cfg.max_processes==1);
end
assert(cfg.inherited_RAM_baseline_bytes>=119021486080);
assert(cfg.new_BIE_invocations_max==25 && cfg.imported_BIE_invocations==5 && cfg.N40_models_limit_s==180);
guards=fullfile(out,'solve_sentinels'); mkdir(guards);
names={'run_k1_recovery_model_job','run_k1_piston_pm_job', ...
    'solve_helmholtz_bie_halfplane_vertical','solve_helmholtz_bie_halfplane_cached_vertical', ...
    'closeout_s3a_pm_fullfk','run_pe_1d_surface_reflection_validation', ...
    'baffled_strip_piston_field_vertical','piston_independent_angular_spectrum_vertical'};
for j=1:numel(names)
    fid=fopen(fullfile(guards,[names{j} '.m']),'w'); assert(fid>=0);
    fprintf(fid,'function varargout=%s(varargin)\nerror(''k1recover:NoSolve'',''Recovery preparation forbids acoustics'');\nend\n',names{j}); fclose(fid);
end
addpath(guards,'-begin'); rehash; for j=1:numel(names), assert(startsWith(which(names{j}),guards)); end
checks=3;
for j=2:6
    c=cfg.cases(j); r=run_k1_pm_recovery_case(fullfile(master,'cases',c.case_id));
    assert(strcmp(r.status,'DRY_RUN_INPUT_CHECKS_PASS') && r.new_acoustic_solves==0); checks=checks+1;
end
old=fileread(fullfile(master,'provenance','original_model_job.m'));
current=fileread(fullfile(master,'code','run_k1_recovery_model_job.m'));
clean=regexprep(current,'(?m)^[ \t]*k1_recovery_event[^\n]*\n','');
expected=strrep(old,'run_k1_piston_pm_job','run_k1_recovery_model_job');
expected=strrep(expected,'    for j=1:numel(PS), PS{j}=local_ps(x,flat,source,profile,recipe(j,:)); end', ...
    sprintf('    for j=1:numel(PS)\n        PS{j}=local_ps(x,flat,source,profile,recipe(j,:));\n    end'));
expected=strrep(expected,'    for j=1:4, FK{j}=closeout_s3a_pm_fullfk(src,profile,x,z,recipes(j,1),recipes(j,2),flat); end', ...
    sprintf('    for j=1:4\n        FK{j}=closeout_s3a_pm_fullfk(src,profile,x,z,recipes(j,1),recipes(j,2),flat);\n    end'));
assert(strcmp(strrep(clean,sprintf('\r'),''),strrep(expected,sprintf('\r'),'')), ...
    'k1recover:PhysicsDrift','Instrumentation changed more than event records/control-line expansion'); checks=checks+1;
assert(strcmp(strrep(extractAfter(old,'function r=local_ps'),sprintf('\r'),''), ...
    strrep(extractAfter(current,'function r=local_ps'),sprintf('\r'),''))); checks=checks+1;
legacy=fullfile(master,'lineage','original_inputs');
assert(strcmp(k1_pm_io('sha',fullfile(master,'cases','U6_N07','manifest.json')), ...
    k1_pm_io('sha',fullfile(legacy,'cases','U6_N07','manifest.json')))); checks=checks+1;
numeric=k1_pm_io('read',fullfile(master,'lineage','N7_local_numeric_receipt.json')); assert(numeric.failed_checks==0); checks=checks+1;
logtest=fullfile(out,'event_fixture'); mkdir(logtest); k1_recovery_event(logtest,'TEST_ONLY','BEGIN'); k1_recovery_event(logtest,'TEST_ONLY','END');
lines=splitlines(strtrim(fileread(fullfile(logtest,'recovery_events.jsonl')))); assert(numel(lines)==2);
event=jsondecode(lines{1}); assert(strcmp(event.state,'BEGIN')); checks=checks+1;
% Reuse original state/gate tests with the immutable original input snapshot.
original=validate_k1_pm_preparation(legacy,fullfile(out,'original_state_gate_tests'));
assert(original.failed_checks==0); checks=checks+original.checks;
pm_sensitivity_io('verify',master,m.files);
receipt=struct('status','K1_RECOVERY_PREPARATION_TESTS_PASS','failed_checks',0,'checks',checks, ...
    'new_acoustic_solves',0,'solve_sentinels',numel(names),'remaining_case_dry_runs',5, ...
    'model_math_unchanged_after_stripping_only_event_lines',true,'elapsed_s',toc(t), ...
    'inherited_used_wall_s',cfg.inherited_used_wall_s,'inherited_RAM_baseline_bytes',cfg.inherited_RAM_baseline_bytes, ...
    'graphics_note','Linux qualified; runtime saved-only preflight is mandatory; Windows export FAIL retained');
k1_pm_io('write',fullfile(out,'receipt.json'),receipt);
fprintf('K1_RECOVERY_PREPARATION_PASS checks=%d elapsed=%.6f new_solves=0\n',checks,toc(t));
end
