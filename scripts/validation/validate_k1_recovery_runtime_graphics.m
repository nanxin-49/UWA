function receipt=validate_k1_recovery_runtime_graphics(master,out)
%VALIDATE_K1_RECOVERY_RUNTIME_GRAPHICS Linux saved replay before any new solve.
assert(isunix && ~ismac && ~isfolder(out)); mkdir(out); t=tic;
source=fullfile(master,'cases','U6_N07','finalization','attempts','attempt_0001','work');
render_k1_pm_stable_saved({source},fullfile(out,'figures'));
a=audit_r2_final_figures({source},fullfile(out,'figures'),fullfile(out,'audit'));
assert(a.failed_checks==0 && a.checked_error_curves==12);
receipt=struct('status','RECOVERY_RUNTIME_SAVED_GRAPHICS_PASS','failed_checks',0, ...
    'audit',a,'elapsed_s',toc(t),'new_acoustic_solves',0,'platform',computer);
k1_pm_io('write',fullfile(out,'receipt.json'),receipt);
end
