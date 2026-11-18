function out=r2_pm_io(action,varargin)
%R2_PM_IO R2 package verification; reuse R1 atomic SHA checkpoint primitives.
if ~strcmp(action,'verify_manifest')
    out=r1_bie_manual_io(action,varargin{:}); return;
end
root=varargin{1}; out=r1_bie_manual_io('read',fullfile(root,'manifest.json'));
assert(out.schema_version==1 && strcmp(out.scope,'R2_PM_TWO_SEED_REVIEW'));
assert(out.max_parallel_seed_processes==2 && out.aggregate_admission_bytes==24*1024^3);
assert(ismember(out.seed,[260002 260003]) && out.single_process && out.single_thread && ...
    ~out.parpool && ~out.automatic_retry && out.RAM_reserve_fraction==.3);
for j=1:numel(out.files)
    f=out.files(j); p=r1_bie_manual_io('path',root,f.path);
    assert(strcmp(r1_bie_manual_io('sha',p),f.sha256),'r2:InputDrift','Changed input/code: %s',f.path);
end
policy=r1_bie_manual_io('read',fullfile(root,'policy.json'));
r1_bie_continuation_io('policy',policy);
assert(strcmp(out.policy_sha256,r1_bie_manual_io('sha',fullfile(root,'policy.json'))));
frozen=r1_bie_manual_io('read',fullfile(root,'input','frozen_config.json'));
expected=struct('frequency_hz',4000,'c_mps',1500,'rho_kgpm3',1025, ...
    'v0_mps',.001,'total_width_m',.1,'z_source_m',100);
assert(frozen.seed==out.seed && isequal(orderfields(frozen.source),orderfields(expected)));
assert(frozen.period_m==160 && frozen.modes==12 && frozen.C2_plateau_halfwidth_m==42 && ...
    frozen.C2_support_halfwidth_m==50 && frozen.receivers==411 && frozen.receiver_z_m==3 && ...
    frozen.phase_mask_fraction_of_BIE_RMS==.01 && strcmp(frozen.phase_contract_version,'2.0.0'));
p=readtable(fullfile(root,'input','coefficients.csv'));
anchor=readtable(fullfile(root,'input','r1_coefficients.csv'));
r=readtable(fullfile(root,'input','receivers.csv'));
assert(height(p)==12 && isequal(p.k_rad_per_m,anchor.k_rad_per_m) && isequal(p.Sk_m3,anchor.Sk_m3));
assert(all(isfinite(p{:,:}),'all') && abs(max(p.k_rad_per_m)-.471238898038469)<1e-14);
assert(height(r)==411 && numel(unique(r.x_m))==411 && issorted(r.x_m) && all(r.z_m==3));
assert(isequal(r.x_m,(-205:205)'*.1953125),'r2:Receivers','Original physical coordinates changed');
% Geometry-only preflight; no source, Green function or acoustic solve.
[eta,slope,curvature]=closeout_s3a_pm_surface(p,linspace(-70,70,14001)');
assert(all(isfinite([eta;slope;curvature])) && min(eta)>-2 && max(eta)<3);
assert(~isequal(p.cos_coeff_m,anchor.cos_coeff_m) || ~isequal(p.sin_coeff_m,anchor.sin_coeff_m));
end
