function out = solve_helmholtz_bie_halfplane_vertical(cfg)
%SOLVE_HELMHOLTZ_BIE_HALFPLANE_VERTICAL Validation-only 2-D sound-soft BIE.
%   Solves the R0-frozen half-plane-Green combined-layer equation on a
%   smooth finite section of z=eta(x). This helper is deliberately isolated
%   from the PE and Bellhop physics implementations.

arguments
    cfg (1,1) struct
end

required = {'frequency_hz','c0_mps','half_width_m','points_per_wavelength', ...
    'panel_order','self_quadrature_order','halfplane_h_m','window_plateau_ratio', ...
    'eta_fn','eta_prime_fn','incident_fn','receiver_x_m','receiver_z_m'};
for ii = 1:numel(required)
    if ~isfield(cfg, required{ii})
        error('Missing BIE configuration field %s.', required{ii});
    end
end

k = 2*pi*cfg.frequency_hz/cfg.c0_mps;
lambda = cfg.c0_mps/cfg.frequency_hz;
if cfg.half_width_m <= 0 || cfg.points_per_wavelength <= 0
    error('half_width_m and points_per_wavelength must be positive.');
end
if cfg.panel_order < 4 || mod(cfg.panel_order, 2) ~= 0
    error('panel_order must be an even integer >=4.');
end
if cfg.self_quadrature_order < 2*cfg.panel_order
    error('self_quadrature_order must be at least twice panel_order.');
end
if cfg.window_plateau_ratio <= 0 || cfg.window_plateau_ratio >= 1
    error('window_plateau_ratio must lie in (0,1).');
end

target_spacing = lambda/cfg.points_per_wavelength;
panel_width_target = target_spacing*cfg.panel_order;
panel_count = max(4, ceil(2*cfg.half_width_m/panel_width_target));
edges = linspace(-cfg.half_width_m, cfg.half_width_m, panel_count+1);
[gx, gw] = local_gauss_legendre(cfg.panel_order);

node_count = panel_count*cfg.panel_order;
x = zeros(node_count,1); z = x; nx = x; nz = x; arc_weight = x;
panel_of_node = zeros(node_count,1); local_coordinate = x;
for pp = 1:panel_count
    jj = (pp-1)*cfg.panel_order + (1:cfg.panel_order);
    half = 0.5*(edges(pp+1)-edges(pp));
    center = 0.5*(edges(pp+1)+edges(pp));
    xp = center + half*gx;
    ep = cfg.eta_fn(xp);
    dep = cfg.eta_prime_fn(xp);
    jac = sqrt(1+dep.^2);
    x(jj) = xp;
    z(jj) = ep;
    nx(jj) = -dep./jac;
    nz(jj) = 1./jac;
    arc_weight(jj) = half*gw.*jac;
    panel_of_node(jj) = pp;
    local_coordinate(jj) = gx;
end

if cfg.halfplane_h_m >= min(z)
    error('halfplane_h_m must be strictly smaller than every wall z value.');
end

window = local_slow_rise_window(x, cfg.half_width_m, cfg.window_plateau_ratio);
weighted_arc = arc_weight.*window;
A = complex(zeros(node_count,node_count));
block_size = min(128,node_count);
for first = 1:block_size:node_count
    rows = first:min(first+block_size-1,node_count);
    [S,K] = local_halfplane_kernels(x(rows),z(rows),x,z,nx,nz,k,cfg.halfplane_h_m);
    A(rows,:) = 2*(K - 1i*k*S).*weighted_arc.';
end

% Replace same- and adjacent-panel blocks by high-order interpolation-based
% quadrature. This includes the logarithmic self singularity without moving
% the boundary or deleting a diagonal term.
for it = 1:node_count
    pt = panel_of_node(it);
    source_panels = max(1,pt-1):min(panel_count,pt+1);
    for ps = source_panels
        cols = (ps-1)*cfg.panel_order + (1:cfg.panel_order);
        singular = ps == pt;
        weights = local_panel_operator_weights(x(it),z(it), ...
            local_coordinate(it),edges(ps),edges(ps+1),gx,k, ...
            cfg.halfplane_h_m,cfg.window_plateau_ratio,cfg.half_width_m, ...
            cfg.eta_fn,cfg.eta_prime_fn,cfg.self_quadrature_order,singular);
        A(it,cols) = 2*weights;
    end
end
A(1:node_count+1:end) = A(1:node_count+1:end) + 1;

u_inc_nodes = cfg.incident_fn(x,z);
rhs = -2*u_inc_nodes;
psi = A\rhs;
linear_residual = norm(A*psi-rhs)/max(norm(rhs),realmin);

receiver_x = cfg.receiver_x_m(:);
receiver_z = cfg.receiver_z_m + zeros(size(receiver_x));
u_receiver = local_evaluate_off_surface(receiver_x,receiver_z,x,z,nx,nz, ...
    weighted_arc,psi,k,cfg.halfplane_h_m);

% Boundary residual is evaluated at panel midpoints, independent of the
% Gauss collocation nodes, and only in the w=1 region.
x_check = 0.5*(edges(1:end-1)+edges(2:end)).';
mask_check = abs(x_check) <= cfg.window_plateau_ratio*cfg.half_width_m;
x_check = x_check(mask_check);
z_check = cfg.eta_fn(x_check);
u_check = cfg.incident_fn(x_check,z_check);
psi_check = zeros(size(x_check));
u_layer_check = complex(zeros(size(x_check)));
for tt = 1:numel(x_check)
    pp = find(x_check(tt) >= edges(1:end-1) & x_check(tt) <= edges(2:end),1,'first');
    if isempty(pp), error('Could not locate boundary check panel.'); end
    for ps = 1:panel_count
        cols = (ps-1)*cfg.panel_order + (1:cfg.panel_order);
        singular = ps == pp;
        weights = local_panel_operator_weights(x_check(tt),z_check(tt),0, ...
            edges(ps),edges(ps+1),gx,k,cfg.halfplane_h_m, ...
            cfg.window_plateau_ratio,cfg.half_width_m,cfg.eta_fn, ...
            cfg.eta_prime_fn,cfg.self_quadrature_order,singular);
        u_layer_check(tt) = u_layer_check(tt) + weights*psi(cols);
    end
    cols = (pp-1)*cfg.panel_order + (1:cfg.panel_order);
    basis0 = local_lagrange_matrix(gx,local_barycentric_weights(gx),0);
    psi_check(tt) = basis0*psi(cols);
end
boundary_total = u_check + 0.5*psi_check + u_layer_check;
boundary_residual = norm(boundary_total)/max(norm(u_check),realmin);

out = struct( ...
    'schema_version','1.0.0', ...
    'formulation','Dirichlet half-plane Green D_h-i*k*S_h; smooth finite section', ...
    'time_convention','exp(-i*omega*t)', ...
    'normal_orientation','into water domain z>eta(x)', ...
    'jump_relation','gamma_D D_h = +1/2 I + K_h', ...
    'frequency_hz',cfg.frequency_hz,'c0_mps',cfg.c0_mps,'k_radpm',k, ...
    'lambda_m',lambda,'halfplane_h_m',cfg.halfplane_h_m, ...
    'half_width_m',cfg.half_width_m,'window_plateau_ratio',cfg.window_plateau_ratio, ...
    'points_per_wavelength_requested',cfg.points_per_wavelength, ...
    'points_per_wavelength_actual',node_count*lambda/(2*cfg.half_width_m), ...
    'panel_order',cfg.panel_order,'panel_count',panel_count,'node_count',node_count, ...
    'self_quadrature_order',cfg.self_quadrature_order, ...
    'node_x_m',x,'node_z_m',z,'node_normal',[nx,nz],'window',window, ...
    'density',psi,'linear_residual',linear_residual, ...
    'boundary_check_x_m',x_check,'boundary_total_residual',boundary_total, ...
    'offgrid_boundary_residual',boundary_residual, ...
    'receiver_x_m',receiver_x,'receiver_z_m',receiver_z,'receiver_field',u_receiver, ...
    'finite',all(isfinite([real(psi);imag(psi);real(u_receiver);imag(u_receiver); ...
        real(boundary_total);imag(boundary_total)])));
end

function weights = local_panel_operator_weights(xt,zt,target_xi,a,b,gx,k,h,c,L,eta_fn,deta_fn,qa,singular)
[ug,uw] = local_gauss_legendre(qa);
if singular
    xi0 = target_xi;
    vl = 0.5*(ug+1); wl = 0.5*uw;
    xi_left = xi0 - (xi0+1)*vl.^2;
    jac_left = 2*(xi0+1)*vl;
    xi_right = xi0 + (1-xi0)*vl.^2;
    jac_right = 2*(1-xi0)*vl;
    xi = [xi_left;xi_right];
    qw = [wl.*jac_left;wl.*jac_right];
else
    xi = ug;
    qw = uw;
end
half = 0.5*(b-a); center = 0.5*(a+b);
xs = center + half*xi;
zs = eta_fn(xs); dzs = deta_fn(xs); jac = sqrt(1+dzs.^2);
nxs = -dzs./jac; nzs = 1./jac;
w = local_slow_rise_window(xs,L,c);
[S,K] = local_halfplane_kernels(xt,zt,xs,zs,nxs,nzs,k,h);
kernel = K - 1i*k*S;
bary = local_barycentric_weights(gx);
B = local_lagrange_matrix(gx,bary,xi);
physical_weights = half*qw.*jac.*w;
weights = (kernel(:).*physical_weights(:)).'*B;
end

function u = local_evaluate_off_surface(xt,zt,xs,zs,nxs,nzs,weights,psi,k,h)
u = complex(zeros(numel(xt),1));
block_size = 256;
source_strength = weights.*psi;
for first = 1:block_size:numel(xt)
    rows = first:min(first+block_size-1,numel(xt));
    [S,K] = local_halfplane_kernels(xt(rows),zt(rows),xs,zs,nxs,nzs,k,h);
    u(rows) = (K-1i*k*S)*source_strength;
end
end

function [S,K] = local_halfplane_kernels(xt,zt,xs,zs,nxs,nzs,k,h)
xt = xt(:); zt = zt(:); xs = xs(:).'; zs = zs(:).';
nxs = nxs(:).'; nzs = nzs(:).';
dx = xt-xs; dz = zt-zs; r = hypot(dx,dz);
singular = r == 0;
r_safe = r; r_safe(singular) = 1;
phi = 1i/4*besselh(0,1,k*r_safe);
dphi = 1i*k/4*besselh(1,1,k*r_safe).*((dx./r_safe).*nxs+(dz./r_safe).*nzs);

z_image = 2*h-zs;
dz_image = zt-z_image;
r_image = hypot(dx,dz_image);
phi_image = 1i/4*besselh(0,1,k*r_image);
dphi_image = 1i*k/4*besselh(1,1,k*r_image).* ...
    ((dx./r_image).*nxs+(dz_image./r_image).*(-nzs));
S = phi-phi_image;
K = dphi-dphi_image;
S(singular) = 0;
K(singular) = 0;
end

function w = local_slow_rise_window(x,L,c)
r = abs(x)/L;
w = ones(size(r));
w(r >= 1) = 0;
mid = r > c & r < 1;
t = (r(mid)-c)/(1-c);
a = exp(-1./t);
b = exp(-1./(1-t));
w(mid) = b./(a+b);
end

function B = local_lagrange_matrix(nodes,bary,x)
nodes = nodes(:).'; bary = bary(:).'; x = x(:);
B = zeros(numel(x),numel(nodes));
for ii = 1:numel(x)
    d = x(ii)-nodes;
    [dm,jj] = min(abs(d));
    if dm <= 32*eps(max(1,abs(x(ii))))
        B(ii,jj) = 1;
    else
        row = bary./d;
        B(ii,:) = row/sum(row);
    end
end
end

function bary = local_barycentric_weights(nodes)
nodes = nodes(:);
n = numel(nodes); bary = ones(n,1);
for jj = 1:n
    bary(jj) = 1/prod(nodes(jj)-nodes([1:jj-1,jj+1:n]));
end
end

function [x,w] = local_gauss_legendre(n)
j = (1:n-1).'; beta = j./sqrt(4*j.^2-1);
[V,D] = eig(diag(beta,1)+diag(beta,-1));
[x,ord] = sort(diag(D)); V = V(:,ord);
w = 2*(V(1,:).^2).';
end
