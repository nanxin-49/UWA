function verify_d1_saved_import(root,work,recovery)
%VERIFY_D1_SAVED_IMPORT Exact same-source, same-grid saved-field admission.
assert(strcmp(pm_sensitivity_io('sha',fullfile(work,'common.mat')),recovery.common_sha256));
assert(strcmp(pm_sensitivity_io('sha',fullfile(work,'jobs','bie_coarse','data.mat')),recovery.C0_data_sha256));
q=load(fullfile(work,'common.mat')); T=readtable(fullfile(root,'inputs','receivers.csv'));
p=readtable(fullfile(root,'inputs','coefficients.csv'));
assert(isequal(q.x,T.x_m) && isequal(q.z,T.z_m) && isequal(q.profile,p));
assert(q.cfg.source.z_source_m==100 && q.cfg.source.total_width_m==.1 && q.cfg.source.frequency_hz==4000);
assert(q.cfg.seed==260001 && q.cfg.modes==12 && q.cfg.U_PM_mps==6 && strcmp(q.cfg.phase_contract_version,'2.0.0'));
assert(isequal(q.cfg.ps_recipe,[128 .05 .1875;256 .05 .1875;512 .05 .1875;1024 .05 .1875;1024 .025 .1875;1024 .025 .09375]));
assert(isequal(q.cfg.fk_recipe,[400 50;800 50;960 60;1440 90]));
assert(all(q.ps_values(1:2)<1e-4) && q.ps_values(3)<1e-10 && all(q.fk_values<1e-5));
a=jsondecode(fileread(fullfile(root,'recovery_audit','receipt.json')));
if strcmp(recovery.kind,'IMPORT_AUDITED_TX100_NO_RESOLVE')
    assert(strcmp(a.status,'SAVED_TX100_DEPTH_AUDIT_PASS') && a.completed_BIE_layers==5 && recovery.remaining_BIE_solves==5);
    ids={'bie_coarse','bie_main','bie_self','bie_space','bie_window'};
    for j=1:5
        assert(strcmp(pm_sensitivity_io('sha',fullfile(work,'jobs',ids{j},'data.mat')),recovery.layer_data_sha256{j}));
    end
    assert(strcmp(pm_sensitivity_io('sha',fullfile(work,'comparison.mat')),recovery.comparison_sha256));
    assert(a.figure_audit.differing_pixel_channels==0 && a.figure_audit.reopened_curves==24);
else
    assert(strcmp(a.status,'SAVED_C0_INDEPENDENT_AUDIT_PASS'));
end
assert(a.failed_checks==0 && a.new_BIE_solves==0);
assert(strcmp(a.origin_manifest_sha256,recovery.origin_manifest_sha256));
end
