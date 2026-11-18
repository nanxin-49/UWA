function receipt=validate_k1_pm_plotting(master,out)
%VALIDATE_K1_PM_PLOTTING Saved R1 arrays as explicit software-only fixtures.
assert(~isfolder(out)); mkdir(out); timer=tic;
oldpath=path; pg=onCleanup(@()path(oldpath)); addpath(fullfile(master,'code'),'-begin');
q=load(fullfile(master,'anchor','common.mat')); a=load(fullfile(master,'anchor','comparison.mat'));
q.cfg.U_PM_mps=6; q.cfg.modes=12; q.cfg.case_id='TEST_ONLY_R1_ANCHOR_REPLAY';
work=fullfile(out,'anchor_replay'); mkdir(work); save(fullfile(work,'common.mat'),'-struct','q','-v7');
copyfile(fullfile(master,'anchor','comparison.mat'),fullfile(work,'comparison.mat'));
render_k1_pm_stable_saved({work},fullfile(out,'case_figures'));
single=audit_r2_final_figures({work},fullfile(out,'case_figures'),fullfile(out,'case_figure_audit'));
ids={'U6_N12_ANCHOR','U6_N07','U6_N40','U5_N12','U8_N12','U5_N40','U8_N40'};
winds=[6 6 6 5 8 5 8]; counts=[12 7 40 12 12 40 40]; records=struct([]);
for j=1:7
    item=struct('case_id',ids{j},'U_PM_mps',winds(j),'N',counts(j),'x',q.x,'U',a.U, ...
        'scale',a.scale,'mask',a.mask,'test_only',true);
    if isempty(records), records=item; else, records(j)=item; end
end
summary=render_k1_pm_summary_fields(records,fullfile(out,'TEST_ONLY_summary_figures'));
assert(single.failed_checks==0 && summary.failed_checks==0 && summary.test_only);
receipt=struct('status','K1_SOFTWARE_PLOTTING_FIXTURE_PASS','failed_checks',0, ...
    'single_case_curves',single.checked_error_curves,'summary_curves',summary.curves_checked, ...
    'zero_pixel_reexports',8,'new_acoustic_solves',0,'new_K1_physical_results',0,'elapsed_s',toc(timer), ...
    'fixture_scope','Historical R1 fields replayed, not six new acoustic cases');
k1_pm_io('write',fullfile(out,'receipt.json'),receipt); disp(jsonencode(receipt));
end
