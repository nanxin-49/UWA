function receipt=probe_k1_saved_graphics(root)
%PROBE_K1_SAVED_GRAPHICS Approved Linux-only saved-field plotting diagnostic.
assert(startsWith(root,'/home/misc/cwx/k1_saved_graphics_diag_'));
cfg=jsondecode(fileread(fullfile(root,'source.json'))); t=tic;
source='/home/misc/cwx/k1_pm_20261118_05/cases/U6_N07/finalization/attempts/attempt_0001/work';
assert(strcmp(cfg.source,source));
for j=1:numel(cfg.files)
    assert(strcmp(r1_bie_manual_io('sha',fullfile(source,cfg.files(j).path)),cfg.files(j).sha256));
end
diary(fullfile(root,'diary.log')); guard=onCleanup(@()diary('off'));
render_k1_pm_stable_saved({source},fullfile(root,'figures'));
a=audit_r2_final_figures({source},fullfile(root,'figures'),fullfile(root,'audit'));
assert(a.failed_checks==0 && a.checked_error_curves==12 && toc(t)<180);
receipt=struct('status','LINUX_SAVED_GRAPHICS_ORIGINAL_GATES_PASS','failed_checks',0, ...
    'audit',a,'elapsed_s',toc(t),'new_acoustic_solves',0,'source_mutated',false, ...
    'MATLAB_version',version,'PID',feature('getpid'));
r1_bie_manual_io('write',fullfile(root,'plot_receipt.json'),receipt);
fprintf('K1_LINUX_SAVED_GRAPHICS_PASS curves=12 elapsed=%.6f\n',toc(t));
end
