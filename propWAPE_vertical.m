function [psiout, psifinal_xy, x, y, z_track, Axz, Ayz, ...
          A_center, R_center, fit_slope, fit_err_rms, pass_1_over_R, fit_mask, ...
          surface_elevation, delta_phi, psi_ref, roughness_meta] = ...
          propWAPE_vertical(cfg)
%PROP WAPE VERTICAL
% Upward marching PE in z (z axis is positive downward).
% Starts at z=z_max and marches to z=0 using negative dz_step.

k0 = 2*pi*cfg.f0/cfg.c0;
% z coordinate marches upward with negative dz_step, while ds is the
% positive path-length increment used in the propagator.
ds = abs(cfg.dz_step);

% --- transverse grid (x-y plane) ---
x = (-0.5*cfg.xw) : cfg.dx : (0.5*cfg.xw - cfg.dx);
y = (-0.5*cfg.yw) : cfg.dy : (0.5*cfg.yw - cfg.dy);
[X, Y] = meshgrid(x, y); % size: ny x nx

[~, ix0] = min(abs(x - cfg.xs));
[~, iy0] = min(abs(y - cfg.ys));

% --- initial Gaussian source at z=z_max ---
psi_space = exp(-((X - cfg.xs).^2 + (Y - cfg.ys).^2) / (2*cfg.sigma_src_m^2));
psi_space = complex(psi_space, 0);

% --- lateral sponge/taper window ---
wx = local_edge_taper(cfg.nx, cfg.taper_ratio);
wy = local_edge_taper(cfg.ny, cfg.taper_ratio);
Wxy = wy * wx.'; % size: ny x nx

% --- spectral grid ---
kx = (2*pi/cfg.xw) * [0:(cfg.nx/2-1), -cfg.nx/2:-1];
ky = (2*pi/cfg.yw) * [0:(cfg.ny/2-1), -cfg.ny/2:-1];
[KX, KY] = meshgrid(kx, ky);
kappa2 = KX.^2 + KY.^2;
denom = sqrt(complex(k0^2 - kappa2, 0)) + k0;
fr0 = exp(-1i * 0.5 * ds * kappa2 ./ denom);

% --- output allocation ---
nout = cfg.nout;
nnout = round(cfg.numstep * (1:nout) / nout);
psiout = complex(zeros(nout, cfg.ny, cfg.nx));

z_track = zeros(cfg.numstep + 1, 1);
z_track(1) = cfg.z_max;

A_center = zeros(cfg.numstep + 1, 1);
A_center(1) = abs(psi_space(iy0, ix0));

R_center = zeros(cfg.numstep + 1, 1);
R_center(1) = 0;

Axz = complex(zeros(cfg.nx, cfg.numstep + 1)); % x-z slice at y=ys
Ayz = complex(zeros(cfg.ny, cfg.numstep + 1)); % y-z slice at x=xs
Axz(:, 1) = psi_space(iy0, :).';
Ayz(:, 1) = psi_space(:, ix0);

% --- upward marching ---
psi_k = fft2(psi_space);
for jj = 1:cfg.numstep
    z_curr = cfg.z_max + jj * cfg.dz_step;
    c_local = local_sound_speed(z_curr, cfg);
    if ~isfinite(c_local) || c_local <= 0
        error('Invalid local sound speed at z=%g m.', z_curr);
    end

    U = (c_local - cfg.c0) / c_local;
    screen_scalar = exp(-1i * k0 * ds * U);

    psi_k = fr0 .* fft2(Wxy .* (screen_scalar .* ifft2(fr0 .* psi_k)));
    psi_space = ifft2(psi_k);

    z_track(jj + 1) = z_curr;
    R_center(jj + 1) = cfg.z_max - z_curr;
    A_center(jj + 1) = abs(psi_space(iy0, ix0));

    Axz(:, jj + 1) = psi_space(iy0, :).';
    Ayz(:, jj + 1) = psi_space(:, ix0);

    hit_idx = find(nnout == jj);
    for kk = 1:numel(hit_idx)
        psiout(hit_idx(kk), :, :) = psi_space;
    end
end

psifinal_xy = psi_space;

[surface_elevation, delta_phi, psi_ref, roughness_meta] = ...
    pm_surface_kirchhoff_module(psifinal_xy, KX, KY, x, y, cfg.xw, cfg.yw, cfg.lambda0);

[fit_slope, fit_err_rms, pass_1_over_R, fit_mask] = ...
    local_validate_one_over_R(A_center, R_center, cfg.lambda0, cfg.z_max);

end

function w = local_edge_taper(n, taper_ratio)
n_taper = round(taper_ratio * n);
w = ones(n, 1);
if n_taper <= 1
    return
end

edge = sin(linspace(0, pi/2, n_taper)).^2;
edge = edge(:);
w(1:n_taper) = edge;
w((n - n_taper + 1):n) = flipud(edge);
end

function c_local = local_sound_speed(z_curr, cfg)
if isstring(cfg.env_mode)
    mode = lower(char(cfg.env_mode));
else
    mode = lower(cfg.env_mode);
end
switch mode
    case 'uniform'
        c_local = cfg.c0;
    case 'layered'
        if isfield(cfg, 'cz_func') && ~isempty(cfg.cz_func)
            c_local = cfg.cz_func(z_curr);
            return
        end

        if isfield(cfg, 'cz_table_z') && isfield(cfg, 'cz_table_c') && ...
           ~isempty(cfg.cz_table_z) && ~isempty(cfg.cz_table_c)
            c_local = interp1(cfg.cz_table_z, cfg.cz_table_c, z_curr, 'linear', 'extrap');
            return
        end

        % Default layered profile adapted from existing test profile.
        wid = 8;
        zla = 20;
        c_local = cfg.c0 + 12*erfc((z_curr/wid) - (zla/wid)) - 24;
    otherwise
        error('Unsupported env_mode: %s', cfg.env_mode);
end
end

function [fit_slope, fit_err_rms, pass_flag, fit_mask] = ...
    local_validate_one_over_R(A_center, R_center, lambda0, z_max)
fit_slope = NaN;
fit_err_rms = NaN;
pass_flag = false;

min_R = max(10*lambda0, 5);
max_R = min(90, z_max - eps);
fit_mask = (R_center >= min_R) & (R_center <= max_R) & ...
           isfinite(A_center) & (A_center > 0);

if nnz(fit_mask) < 10
    return
end

R_fit = R_center(fit_mask);
A_fit = A_center(fit_mask);

p = polyfit(log10(R_fit), log10(A_fit), 1);
fit_slope = p(1);

C_fit = mean(A_fit .* R_fit);
A_model = C_fit ./ R_fit;
fit_err_rms = sqrt(mean(((A_fit - A_model) ./ A_model).^2));

pass_flag = (fit_slope >= -1.05) && (fit_slope <= -0.95) && (fit_err_rms <= 0.10);
end
