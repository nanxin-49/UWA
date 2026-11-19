function receipt=validate_d1_recovery_preparation(root,out)
%VALIDATE_D1_RECOVERY_PREPARATION Real imported arrays and non-solving entry paths.
arguments
    root (1,:) char
    out (1,:) char
end
assert(~isfolder(out)); mkdir(out); diary(fullfile(out,'diary.log')); cleanup=onCleanup(@()diary('off')); t0=tic;
run_d1_depth_comparison(root,'dry-run'); run_d1_depth_comparison(root,'audit-import');
assert(~isfolder(fullfile(root,'output')),'Dry-run must not dispatch');
cfg=jsondecode(fileread(fullfile(root,'config.json'))); q=load(fullfile(root,'imported','depth_0100m','common.mat'));
b=load(fullfile(root,'imported','depth_0100m','jobs','bie_coarse','data.mat'));
T=table(repmat("PS",411,1),q.x,q.z,real(b.field),imag(b.field),true(411,1), ...
    'VariableNames',{'method','x_m','z_m','real_Pa','imag_Pa','phase_valid'});
write_d1_csv(T,fullfile(out,'receiver_mixed.csv')); R=readtable(fullfile(out,'receiver_mixed.csv'));
assert(isequal(R.x_m,q.x) && isequal(R.z_m,q.z) && isequal(R.real_Pa,real(b.field)) && isequal(R.imag_Pa,imag(b.field)));
assert(all(string(R.method)=="PS") && all(logical(R.phase_valid)));
try
    write_d1_csv(T,fullfile(out,'receiver_mixed.csv')); error('d1:MissingOverwriteGuard','Must reject overwrite');
catch ME
    assert(strcmp(ME.identifier,'d1:CSVOverwrite'));
end
assert(ismember(cfg.recovery.remaining_BIE_solves,[5 9]) && cfg.max_wall_s==5400 && cfg.max_depth_wall_s==2400);
if cfg.recovery.remaining_BIE_solves==5
    imported=fullfile(root,'imported','depth_0100m');
    assert(~isfile(fullfile(imported,'independent_reference_resolution.csv')) && ~isfolder(fullfile(imported,'audit')) && ~isfolder(fullfile(imported,'figures')));
    replay=fullfile(out,'finalization_replay'); mkdir(replay);
    for name={'common.mat','comparison.mat','model_metrics.csv','policy.json'}
        copyfile(fullfile(imported,name{1}),fullfile(replay,name{1}));
    end
    copyfile(fullfile(imported,'jobs'),fullfile(replay,'jobs')); audit_d1_final(replay);
    assert(isfile(fullfile(replay,'independent_reference_resolution.csv')),'Finalization must generate a new CSV');
end
paths={'scripts/comparisons/run_d1_depth_comparison.m','scripts/comparisons/run_d1_model_job.m', ...
    'scripts/validation/support/write_d1_csv.m','scripts/validation/support/verify_d1_saved_import.m', ...
    'scripts/validation/audit_d1_recovery_origin.m','scripts/validation/validate_d1_recovery_preparation.m', ...
    'scripts/validation/support/audit_d1_final.m','scripts/reporting/render_d1_depth_results.m','scripts/validation/audit_d1_saved_tx100.m'};
repo=fileparts(fileparts(fileparts(mfilename('fullpath')))); issues=struct;
for j=1:numel(paths)
    msg=checkcode(fullfile(repo,paths{j}),'-id'); issues.(sprintf('file%d',j))=msg;
    assert(isempty(msg),'d1:Analyzer','Code Analyzer messages require review');
end
receipt=struct('status','RECOVERY_MATLAB_PREPARATION_PASS','failed_checks',0,'new_BIE_solves',0, ...
    'elapsed_s',toc(t0),'analyzer_files',numel(paths),'analyzer_messages',0,'imported_receivers',411, ...
    'coordinate_and_complex_roundtrip_exact',true,'remaining_BIE_solves',cfg.recovery.remaining_BIE_solves,'issues',issues, ...
    'unrun_checks',{{'remote native launcher and continuation; user manual only'}},'matlab_version',version);
pm_sensitivity_io('write',fullfile(out,'receipt.json'),receipt); fprintf('%s\n',jsonencode(receipt));
end
