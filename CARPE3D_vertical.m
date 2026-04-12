function output = CARPE3D_vertical(paramsV)
%CARPE3D_VERTICAL
% Independent upward vertical PE driver.

cfg = local_prepare_config(paramsV);

disp('--- CARPE3D_vertical ---')
if numel(cfg.f0) == 1
    if cfg.enable_wideband
        disp(['f0 band (Hz): [', num2str(cfg.f_band_hz(1)), ', ', num2str(cfg.f_band_hz(2)), ...
              '], seed f0=', num2str(cfg.f0)])
    else
        disp(['f0 (Hz): ', num2str(cfg.f0)])
    end
else
    disp(['f0 sweep (Hz): [', num2str(cfg.f0(1)), ', ', num2str(cfg.f0(end)), ...
          '], Nf=', num2str(numel(cfg.f0))])
end
disp(['z march: ', num2str(cfg.z_tx), ' -> ', num2str(cfg.z_rx), ' m, dz_step=', num2str(cfg.dz_step), ' m'])
disp(['grid: nx=', num2str(cfg.nx), ', ny=', num2str(cfg.ny), ...
      ', xw=', num2str(cfg.xw), ' m, yw=', num2str(cfg.yw), ' m'])
disp(['sigma_src_m=', num2str(cfg.sigma_src_m), ', sponge_ratio=', num2str(cfg.sponge_ratio), ...
      ', alpha_max_np_per_m=', num2str(cfg.alpha_max_np_per_m)])
disp(['tx(x,y,z)=(', num2str(cfg.x_tx), ',', num2str(cfg.y_tx), ',', num2str(cfg.z_tx), ') m'])
disp(['rx(x,y,z)=(', num2str(cfg.x_rx), ',', num2str(cfg.y_rx), ',', num2str(cfg.z_rx), ') m'])
disp(['surface reflection enabled=', num2str(cfg.enable_surface_reflection)])
disp(['wideband=', num2str(cfg.enable_wideband), ', save_mode=', cfg.save_mode, ...
      ', use_gpu=', num2str(cfg.use_gpu)])
disp(['surface_reflect_coeff=', num2str(cfg.surface_reflect_coeff), ...
      ', surface_phase_mode=', cfg.surface_phase_mode, ...
      ', surface_oblique_clip=[', num2str(cfg.surface_oblique_clip(1)), ',', ...
      num2str(cfg.surface_oblique_clip(2)), ']'])

[psiout, psifinal_xy, x, y, z_track, Axz, Ayz, ...
 A_center, R_center, fit_slope, fit_err_rms, pass_1_over_R, fit_mask, ...
 surface_elevation, delta_phi, psi_ref, roughness_meta, ...
 h_direct, h_reflect, h_total, rx_state_used, fd_hz_used, ...
 f_axis, H_direct_f, H_reflect_f, H_f, idx_f_ref] = ...
    propWAPE_vertical(cfg);

if cfg.enforce_1_over_R
    if strcmpi(cfg.env_mode, 'uniform') && ~pass_1_over_R
        error('1/R validation failed. slope=%.4f, rel_rms=%.4f', fit_slope, fit_err_rms);
    end
end

output = struct();
output.psifinal_xy = psifinal_xy;
output.psiout = psiout;
output.x = x;
output.y = y;
output.z_track = z_track;
output.Axz = Axz;
output.Ayz = Ayz;
output.A_center = A_center;
output.R_center = R_center;
output.fit_slope = fit_slope;
output.fit_err_rms = fit_err_rms;
output.pass_1_over_R = pass_1_over_R;
output.fit_mask = fit_mask;
output.surface_elevation = surface_elevation;
output.delta_phi = delta_phi;
output.psi_ref = psi_ref;
output.roughness_meta = roughness_meta;
output.h_direct = h_direct;
output.h_reflect = h_reflect;
output.h_total = h_total;
output.f_axis = f_axis;
output.H_direct_f = H_direct_f;
output.H_reflect_f = H_reflect_f;
output.H_f = H_f;
output.idx_f_ref = idx_f_ref;
output.H_baseband = [];
output.rx_state_used = rx_state_used;
output.fd_hz_used = fd_hz_used;
output.rx_amplitude = abs(h_total);
output.rx_phase_rad = angle(h_total);
output.path_loss_db = -20*log10(max(abs(h_total), eps));
output.direct_to_reflect_db = 20*log10(max(abs(h_direct), eps) / max(abs(h_reflect), eps));
if abs(h_reflect) > 0
    output.phase_diff_rad = angle(h_direct) - angle(h_reflect);
else
    output.phase_diff_rad = NaN;
end
output.config = cfg;

if cfg.show_figures
    local_plot_results(output);
end

end

function cfg = local_prepare_config(paramsV)
defaults = struct( ...
    'f0', 4000, ...
    'enable_wideband', false, ...
    'f_band_hz', [4000, 8000], ...
    'Nf_min', 32, ...
    'Nf_max', 128, ...
    'f_ref_hz', [], ...
    'c0', 1500, ...
    'z_max', 100, ...
    'stepz_lamb', 0.5, ...
    'xw', 50, ...
    'yw', 50, ...
    'nx', 1024, ...
    'ny', 1024, ...
    'x_tx', 0, ...
    'y_tx', 0, ...
    'z_tx', 100, ...
    'x_rx', 0, ...
    'y_rx', 0, ...
    'z_rx', 3, ...
    'rx_position_fn', [], ...
    'doppler_fn', [], ...
    'nout', 6, ...
    'sigma_src_m', 0.3, ...
    'taper_ratio', 0.12, ...
    'sponge_ratio', 0.12, ...
    'alpha_max_np_per_m', 0.15, ...
    'env_mode', 'uniform', ...
    'show_figures', true, ...
    'enforce_1_over_R', true, ...
    'enable_surface_reflection', true, ...
    'save_mode', 'rx_only', ...
    'use_gpu', false, ...
    'sea_wind_speed', 5.0, ...
    'sea_hs_target', 0.5, ...
    'sea_seed', 12345, ...
    'surface_reflect_coeff', -1, ...
    'surface_phase_mode', 'normal', ...
    'surface_oblique_clip', [0, 1]);

cfg = defaults;
if nargin > 0 && ~isempty(paramsV)
    fields = fieldnames(paramsV);
    for k = 1:numel(fields)
        cfg.(fields{k}) = paramsV.(fields{k});
    end

    % Backward compatibility for legacy field names.
    if ~isfield(paramsV, 'x_tx') && isfield(paramsV, 'xs')
        cfg.x_tx = paramsV.xs;
    end
    if ~isfield(paramsV, 'y_tx') && isfield(paramsV, 'ys')
        cfg.y_tx = paramsV.ys;
    end
    if ~isfield(paramsV, 'z_tx')
        cfg.z_tx = cfg.z_max;
    end
    if ~isfield(paramsV, 'x_rx')
        cfg.x_rx = cfg.x_tx;
    end
    if ~isfield(paramsV, 'y_rx')
        cfg.y_rx = cfg.y_tx;
    end
    if ~isfield(paramsV, 'sponge_ratio') && isfield(paramsV, 'taper_ratio')
        cfg.sponge_ratio = paramsV.taper_ratio;
    end
end

cfg.nx = local_force_int(cfg.nx, 'nx');
cfg.ny = local_force_int(cfg.ny, 'ny');
cfg.nout = local_force_int(cfg.nout, 'nout');
cfg.sea_seed = local_force_int(cfg.sea_seed, 'sea_seed');

if cfg.nx <= 0 || cfg.ny <= 0
    error('nx and ny must be positive.');
end
if mod(cfg.nx,2) ~= 0 || mod(cfg.ny,2) ~= 0
    error('nx and ny must be even for current FFT spectral layout.');
end
if cfg.nx > 2048 || cfg.ny > 2048
    error('nx and ny must be <= 2048 to avoid memory risk.');
end

if cfg.xw <= 0 || cfg.yw <= 0
    error('xw and yw must be positive.');
end
if cfg.xw > 100 || cfg.yw > 100
    error('xw and yw must be <= 100 m (kilometer-scale domains are forbidden in this mode).');
end

if cfg.z_max <= 0
    error('z_max must be positive.');
end
if ~isnumeric(cfg.f0) || isempty(cfg.f0) || any(~isfinite(cfg.f0(:))) || any(cfg.f0(:) <= 0)
    error('f0 must contain positive finite frequency values.');
end
cfg.f0 = cfg.f0(:).';
if cfg.c0 <= 0
    error('c0 must be positive.');
end
if cfg.stepz_lamb <= 0
    error('stepz_lamb must be positive.');
end
if cfg.z_tx <= 0 || cfg.z_tx > cfg.z_max
    error('z_tx must satisfy 0 < z_tx <= z_max.');
end
if cfg.z_rx < 0 || cfg.z_rx >= cfg.z_tx
    error('z_rx must satisfy 0 <= z_rx < z_tx.');
end

if isstring(cfg.env_mode)
    cfg.env_mode = char(cfg.env_mode);
end
if ~ischar(cfg.env_mode)
    error('env_mode must be ''uniform'' or ''layered''.');
end
if ~isscalar(cfg.show_figures)
    error('show_figures must be a scalar logical flag.');
end
cfg.show_figures = logical(cfg.show_figures);
if ~isscalar(cfg.enforce_1_over_R)
    error('enforce_1_over_R must be a scalar logical flag.');
end
cfg.enforce_1_over_R = logical(cfg.enforce_1_over_R);
if ~isscalar(cfg.enable_surface_reflection)
    error('enable_surface_reflection must be a scalar logical flag.');
end
cfg.enable_surface_reflection = logical(cfg.enable_surface_reflection);
if ~isscalar(cfg.enable_wideband)
    error('enable_wideband must be a scalar logical flag.');
end
cfg.enable_wideband = logical(cfg.enable_wideband);
if ~isscalar(cfg.use_gpu)
    error('use_gpu must be a scalar logical flag.');
end
cfg.use_gpu = logical(cfg.use_gpu);
if ~(isnumeric(cfg.f_band_hz) && numel(cfg.f_band_hz) == 2 && all(isfinite(cfg.f_band_hz(:))) && ...
     cfg.f_band_hz(1) > 0 && cfg.f_band_hz(2) > cfg.f_band_hz(1))
    error('f_band_hz must be [fmin,fmax] with 0<fmin<fmax.');
end
cfg.f_band_hz = cfg.f_band_hz(:).';
cfg.Nf_min = local_force_int(cfg.Nf_min, 'Nf_min');
cfg.Nf_max = local_force_int(cfg.Nf_max, 'Nf_max');
if cfg.Nf_min < 2 || cfg.Nf_max < cfg.Nf_min
    error('Require 2 <= Nf_min <= Nf_max.');
end
if ~(isempty(cfg.f_ref_hz) || (isscalar(cfg.f_ref_hz) && isnumeric(cfg.f_ref_hz) && isfinite(cfg.f_ref_hz) && cfg.f_ref_hz > 0))
    error('f_ref_hz must be empty or a positive finite scalar.');
end
if isstring(cfg.save_mode)
    cfg.save_mode = char(cfg.save_mode);
end
if ~ischar(cfg.save_mode)
    error('save_mode must be ''rx_only'' or ''slice''.');
end
cfg.save_mode = lower(cfg.save_mode);
if ~strcmp(cfg.save_mode, 'rx_only') && ~strcmp(cfg.save_mode, 'slice')
    error('save_mode must be ''rx_only'' or ''slice''.');
end
if ~(isempty(cfg.rx_position_fn) || isa(cfg.rx_position_fn, 'function_handle'))
    error('rx_position_fn must be empty or a function handle.');
end
if ~(isempty(cfg.doppler_fn) || isa(cfg.doppler_fn, 'function_handle'))
    error('doppler_fn must be empty or a function handle.');
end

if cfg.sea_wind_speed <= 0
    error('sea_wind_speed must be positive.');
end
if cfg.sea_hs_target < 0
    error('sea_hs_target must be non-negative.');
end
if ~(isscalar(cfg.surface_reflect_coeff) && isnumeric(cfg.surface_reflect_coeff) && isfinite(cfg.surface_reflect_coeff))
    error('surface_reflect_coeff must be a finite scalar (real or complex).');
end
if isstring(cfg.surface_phase_mode)
    cfg.surface_phase_mode = char(cfg.surface_phase_mode);
end
if ~ischar(cfg.surface_phase_mode)
    error('surface_phase_mode must be ''normal'' or ''oblique''.');
end
cfg.surface_phase_mode = lower(cfg.surface_phase_mode);
if ~strcmp(cfg.surface_phase_mode, 'normal') && ~strcmp(cfg.surface_phase_mode, 'oblique')
    error('surface_phase_mode must be ''normal'' or ''oblique''.');
end
if ~(isnumeric(cfg.surface_oblique_clip) && numel(cfg.surface_oblique_clip) == 2 && ...
     all(isfinite(cfg.surface_oblique_clip(:))))
    error('surface_oblique_clip must be a finite 1x2 numeric range.');
end
cfg.surface_oblique_clip = sort(cfg.surface_oblique_clip(:).');
if cfg.surface_oblique_clip(1) < 0 || cfg.surface_oblique_clip(2) > 1
    error('surface_oblique_clip must satisfy 0 <= min <= max <= 1.');
end

if cfg.sigma_src_m > 2
    error('sigma_src_m > 2 m is not allowed: beam divergence becomes too small for this 100 m upward test.');
end
if cfg.sigma_src_m < 0.2 || cfg.sigma_src_m > 0.5
    error('sigma_src_m must be in [0.2, 0.5] m for this v1 implementation.');
end

if cfg.sponge_ratio < 0.10 || cfg.sponge_ratio > 0.15
    error('sponge_ratio must be in [0.10, 0.15].');
end
if ~(isscalar(cfg.alpha_max_np_per_m) && isfinite(cfg.alpha_max_np_per_m) && cfg.alpha_max_np_per_m > 0)
    error('alpha_max_np_per_m must be a positive finite scalar.');
end

f_for_grid = cfg.f0(1);
if numel(cfg.f0) > 1
    if isempty(cfg.f_ref_hz)
        f_target = 0.5 * (min(cfg.f0) + max(cfg.f0));
    else
        f_target = cfg.f_ref_hz;
    end
    [~, idx_tmp] = min(abs(cfg.f0 - f_target));
    f_for_grid = cfg.f0(idx_tmp);
elseif cfg.enable_wideband
    if isempty(cfg.f_ref_hz)
        f_for_grid = 0.5 * (cfg.f_band_hz(1) + cfg.f_band_hz(2));
    else
        f_for_grid = cfg.f_ref_hz;
    end
end
cfg.lambda0 = cfg.c0 / f_for_grid;
cfg.dx = cfg.xw / cfg.nx;
cfg.dy = cfg.yw / cfg.ny;

if abs(cfg.x_tx) > 0.5*cfg.xw || abs(cfg.y_tx) > 0.5*cfg.yw
    error('Tx location (x_tx,y_tx) must lie inside the transverse domain.');
end
if abs(cfg.x_rx) > 0.5*cfg.xw || abs(cfg.y_rx) > 0.5*cfg.yw
    error('Rx location (x_rx,y_rx) must lie inside the transverse domain.');
end

if cfg.sigma_src_m < max(cfg.dx, cfg.dy)
    error('sigma_src_m must be >= max(dx,dy) to avoid spatial aliasing.');
end

cfg.dz_abs = cfg.stepz_lamb * cfg.lambda0;
if cfg.dz_abs <= 0
    error('Computed dz_abs must be positive.');
end

cfg.path_span = cfg.z_tx - cfg.z_rx;
cfg.numstep = ceil(cfg.path_span / cfg.dz_abs);
cfg.dz_step = -cfg.path_span / cfg.numstep;

if cfg.nout < 1
    error('nout must be >= 1.');
end
cfg.nout = min(cfg.nout, cfg.numstep);

end

function out = local_force_int(v, name)
if ~isscalar(v) || ~isfinite(v)
    error('%s must be a finite scalar.', name);
end
out = round(v);
if abs(out - v) > 1e-9
    error('%s must be an integer.', name);
end
end

function local_plot_results(output)
if ~isempty(output.psifinal_xy)
    abs_final = abs(output.psifinal_xy);
    ref_final = max(abs_final(:));
    final_db = 20*log10(abs_final / max(ref_final, eps));

    figure(11); clf
    pcolor(output.x, output.y, final_db); shading flat
    xlabel('x (m)')
    ylabel('y (m)')
    title('|psi| at z=0 (dB)')
    caxis([-40 0])
    colorbar
    axis equal tight
end

if ~isempty(output.Axz) && ~isempty(output.z_track)
    abs_axz = abs(output.Axz);
    ref_axz = max(abs_axz(:));
    axz_db = 20*log10(abs_axz / max(ref_axz, eps));

    figure(12); clf
    pcolor(output.x, output.z_track, axz_db.'); shading flat
    xlabel('x (m)')
    ylabel('z (m, positive downward)')
    title('|psi(x,z)| at y=y_{tx} (dB)')
    caxis([-40 0])
    colorbar
    set(gca, 'YDir', 'reverse')
end

if ~isempty(output.Ayz) && ~isempty(output.z_track)
    abs_ayz = abs(output.Ayz);
    ref_ayz = max(abs_ayz(:));
    ayz_db = 20*log10(abs_ayz / max(ref_ayz, eps));

    figure(13); clf
    pcolor(output.y, output.z_track, ayz_db.'); shading flat
    xlabel('y (m)')
    ylabel('z (m, positive downward)')
    title('|psi(y,z)| at x=x_{tx} (dB)')
    caxis([-40 0])
    colorbar
    set(gca, 'YDir', 'reverse')
end

if ~isempty(output.R_center) && ~isempty(output.A_center)
    valid = output.R_center > 0 & output.A_center > 0 & isfinite(output.A_center);
    figure(14); clf
    h = [];
    labels = {};
    if any(valid)
        h(end+1) = loglog(output.R_center(valid), output.A_center(valid), 'b-', 'LineWidth', 1.2); %#ok<AGROW>
        labels{end+1} = 'Measured center amplitude'; %#ok<AGROW>
        hold on
    end
    if any(output.fit_mask)
        R_fit = output.R_center(output.fit_mask);
        A_fit = output.A_center(output.fit_mask);
        C_fit = mean(A_fit .* R_fit);
        h(end+1) = loglog(R_fit, C_fit ./ R_fit, 'r--', 'LineWidth', 1.2); %#ok<AGROW>
        labels{end+1} = 'C/R fit'; %#ok<AGROW>
    end
    grid on
    xlabel('R (m)')
    ylabel('|psi_{center}|')
    title(sprintf('1/R check: slope=%.4f, relRMS=%.4f, pass=%d', ...
          output.fit_slope, output.fit_err_rms, output.pass_1_over_R))
    if ~isempty(h)
        legend(h, labels, 'Location', 'southwest')
    end
end

end
