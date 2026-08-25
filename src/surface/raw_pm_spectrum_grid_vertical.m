function spec = raw_pm_spectrum_grid_vertical(U, nx, ny, xw, yw)
%RAW_PM_SPECTRUM_GRID_VERTICAL Audit-ready raw PM spectrum on a periodic grid.
%   SPEC = RAW_PM_SPECTRUM_GRID_VERTICAL(U,NX,NY,XW,YW) reproduces the
%   raw-PM convention used by pm_surface_boundary_model.m and adds analytic
%   infinite-domain variance and finite-grid coverage diagnostics.

arguments
    U (1,1) double {mustBeFinite, mustBePositive}
    nx (1,1) double {mustBeInteger, mustBePositive}
    ny (1,1) double {mustBeInteger, mustBePositive}
    xw (1,1) double {mustBeFinite, mustBePositive}
    yw (1,1) double {mustBeFinite, mustBePositive}
end
if mod(nx, 2) ~= 0 || mod(ny, 2) ~= 0
    error('raw_pm_spectrum_grid_vertical:EvenGridRequired', ...
        'nx and ny must be even to match the project FFT grid convention.');
end

g = 9.81;
alpha_PM = 8.10e-3;
beta_PM = 0.74;
dx = xw / nx;
dy = yw / ny;
dkx = 2*pi / xw;
dky = 2*pi / yw;
kx = dkx * [0:(nx/2-1), -nx/2:-1];
ky = dky * [0:(ny/2-1), -ny/2:-1];
[KX, KY] = meshgrid(kx, ky);
K = hypot(KX, KY);
mask = K > 0;

E1D = zeros(size(K));
Phi2D = zeros(size(K));
A_pm = beta_PM * g^2 / U^4;
E1D(mask) = alpha_PM ./ (2 .* K(mask).^3) .* exp(-A_pm ./ K(mask).^2);
Phi2D(mask) = E1D(mask) ./ (2*pi*K(mask));

variance_discrete = sum(Phi2D(:)) * dkx * dky;
variance_infinite = alpha_PM * U^4 / (4 * beta_PM * g^2);
k_positive = K(mask);
k_min = min(k_positive);
k_peak = sqrt(2 * beta_PM / 3) * g / U^2;
k_nyquist_axis = min(pi / dx, pi / dy);
k_max_corner = max(K(:));
support_capture_radial = max(0, ...
    exp(-A_pm / k_max_corner^2) - exp(-A_pm / k_min^2));

spec = struct();
spec.model = 'raw_pm_isotropic_k_grid_v1';
spec.U_mps = U;
spec.wind_definition = ['U in m/s; measurement height/averaging convention is not ', ...
    'defined by the legacy project formula and must be supplied by the caller.'];
spec.g_mps2 = g;
spec.alpha_PM = alpha_PM;
spec.beta_PM = beta_PM;
spec.nx = nx;
spec.ny = ny;
spec.xw_m = xw;
spec.yw_m = yw;
spec.dx_m = dx;
spec.dy_m = dy;
spec.dkx_rad_per_m = dkx;
spec.dky_rad_per_m = dky;
spec.kx_rad_per_m = kx;
spec.ky_rad_per_m = ky;
spec.KX_rad_per_m = KX;
spec.KY_rad_per_m = KY;
spec.K_rad_per_m = K;
spec.E1D = E1D;
spec.Phi2D = Phi2D;
spec.W_eta_kstat = (2*pi)^2 .* Phi2D;
spec.K_min_rad_per_m = k_min;
spec.K_peak_rad_per_m = k_peak;
spec.K_nyquist_axis_rad_per_m = k_nyquist_axis;
spec.K_max_corner_rad_per_m = k_max_corner;
spec.sigma_eta_discrete2_m2 = variance_discrete;
spec.sigma_eta_infinite2_m2 = variance_infinite;
spec.sigma_eta_discrete_m = sqrt(max(variance_discrete, 0));
spec.Hs_implied_discrete_m = 4 * spec.sigma_eta_discrete_m;
spec.Hs_implied_infinite_m = 4 * sqrt(variance_infinite);
spec.capture_ratio_discrete_to_infinite = variance_discrete / variance_infinite;
spec.capture_ratio_radial_support_idealized = support_capture_radial;
spec.peak_bins_from_origin = k_peak / max(dkx, dky);
spec.formula = ['E(K)=alpha/(2*K^3)*exp(-beta*g^2/(U^4*K^2)); ', ...
    'Phi2D=E/(2*pi*K); sigma_inf^2=alpha*U^4/(4*beta*g^2).'];
spec.capture_note = ['capture_ratio_discrete_to_infinite includes both finite-support ', ...
    'truncation and rectangular-grid quadrature error. The idealized radial-support ', ...
    'ratio is diagnostic only and does not replace aperture/grid convergence.'];
end
