function audit = audit_two_dimensional_closeout_s2c_b_p1_boundary(run_dir)
%AUDIT_TWO_DIMENSIONAL_CLOSEOUT_S2C_B_P1_BOUNDARY Read-only boundary audit.
% Independently reconstructs frozen K_HIGH C2 geometry; no field evaluation.
arguments
    run_dir (1,:) char
end
root=fileparts(fileparts(fileparts(fileparts(run_dir))));
prior=fullfile(root,'results/comparisons/two_dimensional_channel_closeout/20261023_s2c_a_1642');
q=load(fullfile(run_dir,'jobs','bie_main96.mat'),'boundary','src','receipt');
b=q.boundary;
cfg=jsondecode(fileread(fullfile(prior,'protocol','K_HIGH_configuration.json')));
assert(cfg.A_m==.05 && cfg.K_rad_per_m==.2 && cfg.phi_rad==0);
assert(strcmp(q.receipt.status,'COMPLETE') && b.node_count==3840);
assert(strcmp(b.time_convention,'exp(-i*omega*t)'));
assert(strcmp(b.normal_orientation,'into water domain z>eta(x)'));
assert(b.frequency_hz==cfg.frequency_hz && b.c0_mps==cfg.c_m_per_s);
assert(b.halfplane_h_m==-2 && b.half_width_m==60);
assert(b.panel_order==8 && b.self_quadrature_order==96);
assert(numel(q.src.coeff)==984 && all(isfinite(q.src.coeff)));
source=load(fullfile(prior,'protocol','frozen_source.mat'),'src');
assert(isequal(q.src.coeff,source.src.coeff));
assert(strcmp(closeout_s2b_hash(fullfile(prior,'protocol','frozen_source.mat')),cfg.source_sha256));
x=b.node_x_m(:); r=abs(x); w=zeros(size(x)); wp=w;
inside=r<=42; transition=r>42 & r<50;
w(inside)=1;
t=(r(transition)-42)/8;
w(transition)=1-10*t.^3+15*t.^4-6*t.^5;
wp(transition)=(-30*t.^2+60*t.^3-30*t.^4).*sign(x(transition))/8;
eta=.05*sin(.2*x).*w;
slope=.05*(.2*cos(.2*x).*w+sin(.2*x).*wp);
normal=[-slope,ones(size(slope))]./sqrt(1+slope.^2);
height_max=max(abs(b.node_z_m(:)-eta));
normal_max=max(abs(b.node_normal-normal),[],'all');
assert(height_max<1e-13 && normal_max<1e-13);
outside=r>=50;
assert(any(outside) && all(b.node_z_m(outside)==0));
assert(all(b.node_normal(outside,1)==0) && all(b.node_normal(outside,2)==1));
assert(all(isfinite(real(b.density))) && all(isfinite(imag(b.density))));
assert(b.offgrid_boundary_residual<1e-6 && b.linear_residual<1e-10);
audit=struct('pass',true,'source','saved BIE MAT plus independent frozen C2 polynomial', ...
    'node_count',numel(x),'source_modes',numel(q.src.coeff), ...
    'max_height_difference_m',height_max,'max_normal_difference',normal_max, ...
    'flat_exterior_nodes',nnz(outside),'flat_exterior_exact',true, ...
    'source_spectrum_exact_match',true,'density_finite',true, ...
    'offgrid_boundary_residual',b.offgrid_boundary_residual, ...
    'linear_residual',b.linear_residual,'reference_qualified',false, ...
    'no_new_acoustic_calculation',true);
closeout_s2b_json(fullfile(run_dir,'audit','P1_boundary_independent_audit.json'),audit);
disp(audit);
end
