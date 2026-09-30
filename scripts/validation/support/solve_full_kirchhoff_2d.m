function out = solve_full_kirchhoff_2d(cfg)
%SOLVE_FULL_KIRCHHOFF_2D Validation-only 2-D Kirchhoff line integral.
%   Uses the physical parametric boundary z=eta(x), the free-space outgoing
%   Helmholtz Green function and the Dirichlet Kirchhoff approximation
%       u_ref(X) = - integral [ G d_n u_inc + (d_n G) u_inc ] ds.
%   The normal points into the water domain. This file is intentionally
%   independent of production PE, BIE, and surface-model code.

arguments
    cfg (1,1) struct
end
required = {'frequency_hz','c0_mps','surface_x_m','surface_z_m', ...
    'surface_normal','surface_incident','surface_incident_dn', ...
    'receiver_x_m','receiver_z_m'};
for ii=1:numel(required)
    assert(isfield(cfg,required{ii}),'Missing field %s.',required{ii});
end
k = 2*pi*cfg.frequency_hz/cfg.c0_mps;
xs = cfg.surface_x_m(:); zs = cfg.surface_z_m(:);
assert(numel(xs)>=3 && all(diff(xs)>0),'Surface samples must be ordered.');
assert(numel(zs)==numel(xs),'Surface coordinate sizes differ.');
nv = cfg.surface_normal;
assert(size(nv,1)==numel(xs) && size(nv,2)==2,'surface_normal must be N-by-2.');
ui = cfg.surface_incident(:); dni = cfg.surface_incident_dn(:);
assert(numel(ui)==numel(xs) && numel(dni)==numel(xs),'Incident data size mismatch.');

% Trapezoidal quadrature in the parameter x. The default retains the exact
% local ds. The parameter-x option is a validation-only ablation used to
% isolate the surface-Jacobian contribution.
measure = "arc_length";
if isfield(cfg,'quadrature_measure'), measure=string(cfg.quadrature_measure); end
if measure=="arc_length"
    segment_length = hypot(diff(xs),diff(zs));
elseif measure=="parameter_x"
    segment_length = diff(xs);
else
    error('Unknown quadrature_measure %s.',measure);
end
qw = [0.5*segment_length(1); ...
    0.5*(segment_length(1:end-1)+segment_length(2:end)); ...
    0.5*segment_length(end)];

xr = cfg.receiver_x_m(:); zr = cfg.receiver_z_m(:);
assert(numel(zr)==numel(xr),'Receiver coordinate sizes differ.');
ur = complex(zeros(size(xr)));
local_radius = Inf;
if isfield(cfg,'local_radius_m'), local_radius=cfg.local_radius_m; end
assert(isscalar(local_radius) && local_radius>0,'local_radius_m must be positive.');
block = 256;
for first=1:block:numel(xr)
    rows=first:min(first+block-1,numel(xr));
    dxm = xr(rows)-xs.'; dzm = zr(rows)-zs.';
    rr = hypot(dxm,dzm); rr=max(rr,eps);
    G = 1i/4*besselh(0,1,k*rr);
    % Derivative with respect to the source point in the water normal.
    dGdn = 1i*k/4*besselh(1,1,k*rr).* ...
        (dxm.*nv(:,1).'+dzm.*nv(:,2).')./rr;
    integrand = (G.*dni.' + dGdn.*ui.').*qw.';
    if isfinite(local_radius)
        integrand(abs(dxm)>local_radius)=0;
    end
    ur(rows) = -sum(integrand,2);
end

out = struct('schema_version','1.0.0','formulation', ...
    'Dirichlet Kirchhoff line integral on z=eta(x)', ...
    'time_convention','exp(-i*omega*t)', ...
    'normal_orientation','into water domain z>eta(x)', ...
    'frequency_hz',cfg.frequency_hz,'c0_mps',cfg.c0_mps,'k_radpm',k, ...
    'surface_x_m',xs,'surface_z_m',zs,'surface_normal',nv, ...
    'surface_incident',ui,'surface_incident_dn',dni, ...
    'quadrature_measure',measure,'local_radius_m',local_radius, ...
    'receiver_x_m',xr,'receiver_z_m',zr,'receiver_field',ur, ...
    'quadrature_weights_m',qw,'finite',all(isfinite([real(ur);imag(ur)])));
end
