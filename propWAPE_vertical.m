function [psiout, psifinal_xy, x, y, z_track, Axz, Ayz, ...
          A_center, R_center, fit_slope, fit_err_rms, pass_1_over_R, fit_mask, ...
          surface_elevation, delta_phi, psi_ref, roughness_meta, ...
          h_direct, h_reflect, h_total, rx_state_used, fd_hz_used] = ...
          propWAPE_vertical(cfg)
%PROP WAPE VERTICAL
% Upward marching PE in z (z axis is positive downward).
% Starts at z=z_tx and marches to z=z_rx using negative dz_step.

k0 = 2*pi*cfg.f0/cfg.c0;

% --- transverse grid (x-y plane) ---
x = (-0.5*cfg.xw) : cfg.dx : (0.5*cfg.xw - cfg.dx);
y = (-0.5*cfg.yw) : cfg.dy : (0.5*cfg.yw - cfg.dy);
[X, Y] = meshgrid(x, y); % size: ny x nx

[~, ix_tx] = min(abs(x - cfg.x_tx));
[~, iy_tx] = min(abs(y - cfg.y_tx));

rx_state_used = struct('x_rx', cfg.x_rx, 'y_rx', cfg.y_rx, ...
                       'z_rx', cfg.z_rx, 't_s', 0, 'source', 'fixed');
if ~isempty(cfg.rx_position_fn)
    rx_state = struct('x_tx', cfg.x_tx, 'y_tx', cfg.y_tx, 'z_tx', cfg.z_tx, ...
                      'x_rx_nominal', cfg.x_rx, 'y_rx_nominal', cfg.y_rx, 'z_rx_nominal', cfg.z_rx, ...
                      'x', x, 'y', y);
    rx_xyz = cfg.rx_position_fn(0, rx_state);
    if ~(isnumeric(rx_xyz) && numel(rx_xyz) == 3 && all(isfinite(rx_xyz(:))))
        error('rx_position_fn must return a finite numeric [x_rx, y_rx, z_rx].');
    end
    rx_state_used.x_rx = rx_xyz(1);
    rx_state_used.y_rx = rx_xyz(2);
    rx_state_used.z_rx = rx_xyz(3);
    rx_state_used.source = 'rx_position_fn';
end

if abs(rx_state_used.x_rx) > 0.5*cfg.xw || abs(rx_state_used.y_rx) > 0.5*cfg.yw
    error('Resolved Rx position must lie inside the transverse domain.');
end
if rx_state_used.z_rx < 0 || rx_state_used.z_rx >= cfg.z_tx
    error('Resolved z_rx must satisfy 0 <= z_rx < z_tx.');
end

[~, ix_rx] = min(abs(x - rx_state_used.x_rx));
[~, iy_rx] = min(abs(y - rx_state_used.y_rx));
rx_state_used.ix_rx = ix_rx;
rx_state_used.iy_rx = iy_rx;

path_span_used = cfg.z_tx - rx_state_used.z_rx;
numstep = ceil(path_span_used / cfg.dz_abs);
numstep = max(1, numstep);
dz_step_used = -path_span_used / numstep;
ds = abs(dz_step_used);

fd_hz_used = 0;
if ~isempty(cfg.doppler_fn)
    tx_state = struct('x_tx', cfg.x_tx, 'y_tx', cfg.y_tx, 'z_tx', cfg.z_tx);
    env_state = struct('z_max', cfg.z_max, 'enable_surface_reflection', cfg.enable_surface_reflection);
    fd_hz_used = cfg.doppler_fn(0, tx_state, rx_state_used, env_state);
    if ~(isscalar(fd_hz_used) && isnumeric(fd_hz_used) && isfinite(fd_hz_used))
        error('doppler_fn must return a finite scalar fd_hz.');
    end
end

% --- initial Gaussian source at z=z_tx ---
psi_space = exp(-((X - cfg.x_tx).^2 + (Y - cfg.y_tx).^2) / (2*cfg.sigma_src_m^2));
psi_space = complex(psi_space, 0);
psi_init = psi_space;

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
nout = min(cfg.nout, numstep);
nnout = round(numstep * (1:nout) / nout);
psiout = complex(zeros(nout, cfg.ny, cfg.nx));

z_track = zeros(numstep + 1, 1);
z_track(1) = cfg.z_tx;

A_center = zeros(numstep + 1, 1);
A_center(1) = abs(psi_space(iy_tx, ix_tx));

R_center = zeros(numstep + 1, 1);
R_center(1) = 0;

Axz = complex(zeros(cfg.nx, numstep + 1)); % x-z slice at y=y_tx
Ayz = complex(zeros(cfg.ny, numstep + 1)); % y-z slice at x=x_tx
Axz(:, 1) = psi_space(iy_tx, :).';
Ayz(:, 1) = psi_space(:, ix_tx);

% --- upward marching ---
psi_k = fft2(psi_space);
for jj = 1:numstep
    z_curr = cfg.z_tx + jj * dz_step_used;
    c_local = local_sound_speed(z_curr, cfg);
    if ~isfinite(c_local) || c_local <= 0
        error('Invalid local sound speed at z=%g m.', z_curr);
    end

    U = (c_local - cfg.c0) / c_local;
    screen_scalar = exp(-1i * k0 * ds * U);

    psi_k = fr0 .* fft2(Wxy .* (screen_scalar .* ifft2(fr0 .* psi_k)));
    psi_space = ifft2(psi_k);

    z_track(jj + 1) = z_curr;
    R_center(jj + 1) = cfg.z_tx - z_curr;
    A_center(jj + 1) = abs(psi_space(iy_tx, ix_tx));

    Axz(:, jj + 1) = psi_space(iy_tx, :).';
    Ayz(:, jj + 1) = psi_space(:, ix_tx);

    hit_idx = find(nnout == jj);
    for kk = 1:numel(hit_idx)
        psiout(hit_idx(kk), :, :) = psi_space;
    end
end

psifinal_xy = psi_space;
h_direct = psifinal_xy(iy_rx, ix_rx);
h_reflect = complex(0, 0);
h_total = h_direct;

if cfg.enable_surface_reflection
    pm_cfg = struct();
    pm_cfg.U = cfg.sea_wind_speed;
    pm_cfg.Hs_target = cfg.sea_hs_target;
    pm_cfg.seed = cfg.sea_seed;
    pm_cfg.show_figure = cfg.show_figures;
    pm_cfg.reflect_coeff = cfg.surface_reflect_coeff;
    pm_cfg.phase_mode = cfg.surface_phase_mode;
    pm_cfg.oblique_clip = cfg.surface_oblique_clip;
    pm_cfg.tx_xyz = [cfg.x_tx, cfg.y_tx, cfg.z_tx];
    pm_cfg.rx_xyz = [rx_state_used.x_rx, rx_state_used.y_rx, rx_state_used.z_rx];
    pm_cfg.z_surface = 0;

    psi_surface_inc = local_march_field(psi_init, cfg.z_tx, 0, cfg, Wxy, kappa2, k0);

    [surface_elevation, delta_phi, psi_ref, roughness_meta] = ...
        pm_surface_kirchhoff_module(psi_surface_inc, KX, KY, x, y, cfg.xw, cfg.yw, cfg.lambda0, pm_cfg);

    psi_ref_at_rx = local_march_field(psi_ref, 0, rx_state_used.z_rx, cfg, Wxy, kappa2, k0);
    h_reflect = psi_ref_at_rx(iy_rx, ix_rx);
    h_total = h_direct + h_reflect;

    roughness_meta.reflection_model = 'two_segment_tx_surface_rx';
else
    surface_elevation = zeros(size(psifinal_xy));
    delta_phi = zeros(size(psifinal_xy));
    psi_ref = complex(zeros(size(psifinal_xy)));
    roughness_meta = struct( ...
        'enabled', false, ...
        'reflection_coeff_used', cfg.surface_reflect_coeff, ...
        'phase_mode_used', cfg.surface_phase_mode, ...
        'phase_factor_stats', struct('min', NaN, 'max', NaN, 'mean', NaN), ...
        'phi2d_transform_error', struct('max_rel', NaN, 'mean_rel', NaN));
end

[fit_slope, fit_err_rms, pass_1_over_R, fit_mask] = ...
    local_validate_one_over_R(A_center, R_center, cfg.lambda0, path_span_used);

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

function psi_end = local_march_field(psi_start, z_start, z_end, cfg, Wxy, kappa2, k0)
%LOCAL_MARCH_FIELD March one complex field between two depths.

if abs(z_end - z_start) <= eps
    psi_end = psi_start;
    return
end

n_step = ceil(abs(z_end - z_start) / cfg.dz_abs);
n_step = max(1, n_step);
dz_step_local = (z_end - z_start) / n_step;
ds_local = abs(dz_step_local);

denom = sqrt(complex(k0^2 - kappa2, 0)) + k0;
fr_local = exp(-1i * 0.5 * ds_local * kappa2 ./ denom);

psi_k = fft2(psi_start);
for jj = 1:n_step
    z_curr = z_start + jj * dz_step_local;
    c_local = local_sound_speed(z_curr, cfg);
    if ~isfinite(c_local) || c_local <= 0
        error('Invalid local sound speed at z=%g m.', z_curr);
    end

    U = (c_local - cfg.c0) / c_local;
    screen_scalar = exp(-1i * k0 * ds_local * U);

    psi_k = fr_local .* fft2(Wxy .* (screen_scalar .* ifft2(fr_local .* psi_k)));
end

psi_end = ifft2(psi_k);
end
