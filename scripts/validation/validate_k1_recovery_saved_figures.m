function receipt=validate_k1_recovery_saved_figures(work,out)
%VALIDATE_K1_RECOVERY_SAVED_FIGURES Saved-only, scoped graphics defaults.
assert(isfolder(work) && ~isfolder(out)); mkdir(out);
style=get(groot,'DefaultFigureWindowStyle'); pos=get(groot,'DefaultFigurePosition');
g=onCleanup(@()set(groot,'DefaultFigureWindowStyle',style,'DefaultFigurePosition',pos));
set(groot,'DefaultFigureWindowStyle','normal','DefaultFigurePosition',[40 40 1800 1400]);
t=tic; diary(fullfile(out,'diary.log')); dg=onCleanup(@()diary('off'));
try
    render_k1_recovery_saved({work},fullfile(out,'figures'));
    audit=audit_k1_recovery_figures({work},fullfile(out,'figures'),fullfile(out,'audit'));
    assert(toc(t)<90 && audit.failed_checks==0);
    receipt=struct('status','A0_SCOPED_NORMAL_WINDOW_FIGURE_PASS','failed_checks',0, ...
        'elapsed_s',toc(t),'audit',audit,'new_acoustic_solves',0, ...
        'original_default_style',style,'original_default_position',pos,'defaults_restored_on_exit',true);
    k1_pm_io('write',fullfile(out,'receipt.json'),receipt);
    fprintf('A0_NORMAL_WINDOW_PASS elapsed=%.3f curves=%d\n',toc(t),audit.checked_error_curves);
catch ME
    k1_pm_io('write',fullfile(out,'failure.json'),struct('status','FAIL','error',ME.message,'identifier',ME.identifier,'elapsed_s',toc(t)));
    rethrow(ME);
end
end
