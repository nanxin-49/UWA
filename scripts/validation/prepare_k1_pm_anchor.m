function receipt=prepare_k1_pm_anchor(package_root)
%PREPARE_K1_PM_ANCHOR Strict saved R1/K0 anchor match; never evaluates acoustics.
pm_sensitivity_io('absolute',package_root); root=fullfile(package_root,'anchor');
target=fullfile(root,'anchor_proof.json'); assert(~isfile(target));
q=load(fullfile(root,'common.mat')); a=load(fullfile(root,'comparison.mat'));
original=readtable(fullfile(root,'R1_coefficients.csv'));
k0=readtable(fullfile(root,'K0_U6_N12_coefficients.csv'));
receivers=readtable(fullfile(root,'receivers.csv')); prior=readtable(fullfile(root,'R1_receiver_origin.csv'));
assert(height(original)==12 && height(k0)==12 && isequal(q.profile,original));
assert(isequal(k0.k_rad_per_m,original.k_rad_per_m) && isequal(k0.Sk_m3,original.Sk_m3));
assert(isequal(k0.cos_depth_coeff_m,original.cos_coeff_m) && isequal(k0.sin_depth_coeff_m,original.sin_coeff_m));
source=struct('frequency_hz',4000,'c_mps',1500,'rho_kgpm3',1025,'v0_mps',.001,'total_width_m',.1,'z_source_m',100);
assert(isequal(orderfields(q.cfg.source),orderfields(source)) && q.cfg.seed==260001);
assert(q.cfg.period_m==160 && q.cfg.surface_plateau_halfwidth_m==42 && q.cfg.surface_support_halfwidth_m==50);
assert(strcmp(q.cfg.phase_contract_version,'2.0.0') && strcmp(q.cfg.representation,'absolute_physical peak Pa; exp(-i omega t)'));
assert(q.cfg.phase_mask_fraction_of_BIE_RMS==.01);
assert(isequal(q.x,(-205:205)'*.1953125) && isequal(q.z,3*ones(411,1)));
assert(isequal(q.x,receivers.x_m,prior.x_m) && isequal(q.z,receivers.z_m,prior.z_m));
[depth,slope,second]=closeout_s3a_pm_surface(original,q.sx);
geometry_error=max(abs([depth-q.eta;slope-q.slope;second-q.curvature])); assert(geometry_error<1e-12);
assert(isequal(q.cfg.ps_recipe,[128 .05 .1875;256 .05 .1875;512 .05 .1875;1024 .05 .1875;1024 .025 .1875;1024 .025 .09375]));
assert(isequal(q.cfg.fk_recipe,[400 50;800 50;960 60;1440 90]));
assert(all(q.ps_values(1:2)<1e-4) && q.ps_values(3)<1e-10 && all(q.fk_values<1e-5));
review=jsondecode(fileread(fullfile(root,'R1_formal_review.json')));
assert(strcmp(pm_sensitivity_io('sha',fullfile(root,'R1_formal_review.json')), ...
    'eb6d86a8e845ee2785efab4b93470f91d26026ed0d415bdbd808023152b43265'));
assert(strcmp(review.status,'SAVED_RESULT_REVIEW_PASS') && review.failed_checks==0 && review.reference_stage==3);
numeric=jsondecode(fileread(fullfile(root,'audit','receipt.json')));
final=jsondecode(fileread(fullfile(root,'independent_final_metrics_audit.json')));
assert(numeric.failed_checks==0 && final.failed_checks==0 && a.reference_qualified && strcmp(a.reference_job,'bie_space'));
assert(all(a.bie_values<1e-4) && isequal(a.mask,abs(a.U(:,2))>=.01*norm(a.U(:,2))/sqrt(411)));
policy=jsondecode(fileread(fullfile(root,'policy.json'))); r1_bie_continuation_io('policy',policy);
seal=jsondecode(fileread(fullfile(root,'R1_final_seal.json')));
for file={'common.mat','comparison.mat','model_metrics.csv','receiver_fields_and_errors.csv','independent_final_metrics_audit.json','audit/receipt.json'}
    suffix=['work/' file{1}]; matches=arrayfun(@(v)endsWith(v.path,suffix),seal.files);
    assert(nnz(matches)==1,'k1:AnchorSeal','Missing R1 final qualification binding');
    assert(strcmp(seal.files(matches).sha256,pm_sensitivity_io('sha',fullfile(root,file{1}))));
end
receipt=struct('status','STRICT_R1_K0_ANCHOR_MATCH_PASS','failed_checks',0,'seed',260001,'U_PM_mps',6,'N',12, ...
    'source',source,'phase_contract_version','2.0.0','receivers',411,'C2_geometry_max_abs_difference',geometry_error, ...
    'receiver_coordinates_sha256',pm_sensitivity_io('sha',fullfile(root,'receivers.csv')), ...
    'R1_coefficients_sha256',pm_sensitivity_io('sha',fullfile(root,'R1_coefficients.csv')), ...
    'K0_coefficients_sha256',pm_sensitivity_io('sha',fullfile(root,'K0_U6_N12_coefficients.csv')), ...
    'R1_formal_qualification_sha256',pm_sensitivity_io('sha',fullfile(root,'R1_formal_review.json')), ...
    'R1_final_seal_sha256',pm_sensitivity_io('sha',fullfile(root,'R1_final_seal.json')), ...
    'policy_sha256',pm_sensitivity_io('sha',fullfile(root,'policy.json')), ...
    'reused_scope','exact U6/N12 historical anchor only; no density transferred to new surfaces', ...
    'acoustic_calls',0);
pm_sensitivity_io('write',target,receipt); disp(jsonencode(receipt));
end
