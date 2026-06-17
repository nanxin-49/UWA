function [psiout, psifinal_xy, x, y, z_track, Axz, Ayz, ...
          A_center, R_center, fit_slope, fit_err_rms, pass_1_over_R, fit_mask, ...
          surface_elevation, delta_phi, psi_ref, roughness_meta, ...
          h_direct, h_reflect, h_total, rx_state_used, fd_hz_used, ...
          f_axis, H_direct_f, H_reflect_f, H_f, idx_f_ref, bubble_meta] = ...
          vertical_wape_propagator(cfg)
%VERTICAL_WAPE_PROPAGATOR
% Upward marching PE in z (z axis is positive downward).
% Supports scalar or wideband frequency sweep.

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
if path_span_used <= 0
    error('Propagation span must be positive (z_tx > z_rx).');
end

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
psi_init_cpu = exp(-((X - cfg.x_tx).^2 + (Y - cfg.y_tx).^2) / (2*cfg.sigma_src_m^2));
psi_init_cpu = complex(psi_init_cpu, 0);

% --- spectral grid ---
kx = (2*pi/cfg.xw) * [0:(cfg.nx/2-1), -cfg.nx/2:-1];
ky = (2*pi/cfg.yw) * [0:(cfg.ny/2-1), -cfg.ny/2:-1];
[KX, KY] = meshgrid(kx, ky);
kappa2 = KX.^2 + KY.^2;

% --- wideband frequency axis ---
[f_axis, idx_f_ref] = local_resolve_frequency_axis(cfg, rx_state_used);
Nf = numel(f_axis);
H_direct_f = complex(zeros(Nf, 1));
H_reflect_f = complex(zeros(Nf, 1));
H_f = complex(zeros(Nf, 1));

% --- complex absorbing boundary ---
alpha_xy = local_absorption_profile(x, y, cfg.xw, cfg.yw, cfg.sponge_ratio, cfg.alpha_max_np_per_m);

use_gpu = logical(cfg.use_gpu) && local_has_gpu();
if use_gpu
    kappa2_work = gpuArray(kappa2);
    alpha_xy_work = gpuArray(alpha_xy);
else
    kappa2_work = kappa2;
    alpha_xy_work = alpha_xy;
end

save_mode = lower(char(cfg.save_mode));
save_slice = strcmp(save_mode, 'slice');

% --- defaults for outputs (filled by reference frequency) ---
psiout = complex(zeros(0, 0, 0));
psifinal_xy = [];
z_track = [];
Axz = complex(zeros(0, 0));
Ayz = complex(zeros(0, 0));
A_center = [];
R_center = [];
fit_slope = NaN;
fit_err_rms = NaN;
pass_1_over_R = false;
fit_mask = false(0, 1);
surface_elevation = [];
delta_phi = [];
psi_ref = [];
roughness_meta = local_disabled_roughness_meta(cfg);
bubble_meta = struct('enabled', false, 'model', 'off');

for ifq = 1:Nf
    f_hz = f_axis(ifq);
    lambda_f = cfg.c0 / f_hz;
    k0 = 2*pi / lambda_f;
    dz_abs_f = cfg.stepz_lamb * lambda_f;
    numstep_f = max(1, ceil(path_span_used / dz_abs_f));
    dz_step_f = -path_span_used / numstep_f;
    ds = abs(dz_step_f);

    denom = sqrt(complex(k0^2 - kappa2, 0)) + k0;
    fr0_cpu = exp(-1i * 0.5 * ds * kappa2 ./ denom);
    if use_gpu
        fr0 = gpuArray(fr0_cpu);
    else
        fr0 = fr0_cpu;
    end

    capture_ref = (ifq == idx_f_ref);
    psi_k = fft2(psi_init_cpu);
    if use_gpu
        psi_k = gpuArray(psi_k);
    end

    if capture_ref
        z_track_f = zeros(numstep_f + 1, 1);
        z_track_f(1) = cfg.z_tx;
        A_center_f = zeros(numstep_f + 1, 1);
        A_center_f(1) = abs(psi_init_cpu(iy_tx, ix_tx));
        R_center_f = zeros(numstep_f + 1, 1);
        bubble_path_meta = struct([]);
        if save_slice
            nout_f = min(cfg.nout, numstep_f);
            nnout_f = round(numstep_f * (1:nout_f) / nout_f);
            psiout_f = complex(zeros(nout_f, cfg.ny, cfg.nx));
            Axz_f = complex(zeros(cfg.nx, numstep_f + 1));
            Ayz_f = complex(zeros(cfg.ny, numstep_f + 1));
            Axz_f(:, 1) = psi_init_cpu(iy_tx, :).';
            Ayz_f(:, 1) = psi_init_cpu(:, ix_tx);
        else
            nnout_f = [];
            psiout_f = complex(zeros(0, 0, 0));
            Axz_f = complex(zeros(0, 0));
            Ayz_f = complex(zeros(0, 0));
        end
    end

    for jj = 1:numstep_f
        z_curr = cfg.z_tx + jj * dz_step_f;
        [screen, bubble_step_meta] = local_phase_screen( ...
            x, y, z_curr, f_hz, cfg, alpha_xy_work, k0, ds, use_gpu);
        psi_k = fr0 .* fft2(screen .* ifft2(fr0 .* psi_k));

        if capture_ref
            bubble_path_meta = local_update_bubble_path_meta(bubble_path_meta, bubble_step_meta);
            psi_step = ifft2(psi_k);
            if use_gpu
                center_val = gather(psi_step(iy_tx, ix_tx));
            else
                center_val = psi_step(iy_tx, ix_tx);
            end
            z_track_f(jj + 1) = z_curr;
            R_center_f(jj + 1) = cfg.z_tx - z_curr;
            A_center_f(jj + 1) = abs(center_val);

            if save_slice
                if use_gpu
                    psi_step_cpu = gather(psi_step);
                else
                    psi_step_cpu = psi_step;
                end
                Axz_f(:, jj + 1) = psi_step_cpu(iy_tx, :).';
                Ayz_f(:, jj + 1) = psi_step_cpu(:, ix_tx);
                hit_idx = find(nnout_f == jj);
                for kk = 1:numel(hit_idx)
                    psiout_f(hit_idx(kk), :, :) = psi_step_cpu;
                end
            end
        end
    end

    if capture_ref
        bubble_meta = bubble_path_meta;
        bubble_meta.f_axis = f_axis;
        bubble_meta.path = 'direct_tx_rx';
    end

    psi_end = ifft2(psi_k);
    if use_gpu
        h_direct_fi = gather(psi_end(iy_rx, ix_rx));
    else
        h_direct_fi = psi_end(iy_rx, ix_rx);
    end

    h_reflect_fi = complex(0, 0);
    if cfg.enable_surface_reflection
        pm_cfg = struct();
        pm_cfg.U = cfg.sea_wind_speed;
        pm_cfg.Hs_target = cfg.sea_hs_target;
        pm_cfg.seed = cfg.sea_seed;
        pm_cfg.show_figure = cfg.show_figures && capture_ref;
        pm_cfg.reflect_coeff = cfg.surface_reflect_coeff;
        pm_cfg.phase_mode = cfg.surface_phase_mode;
        pm_cfg.oblique_clip = cfg.surface_oblique_clip;
        pm_cfg.boundary_model = cfg.surface_boundary_model;
        pm_cfg.boundary_check_equivalence = cfg.surface_boundary_check_equivalence;
        pm_cfg.boundary_debug = cfg.surface_boundary_debug;
        pm_cfg.boundary_equivalence_tol = cfg.surface_boundary_equivalence_tol;
        pm_cfg.boundary_coupling_diagnostics = cfg.surface_boundary_coupling_diagnostics && capture_ref;
        pm_cfg.boundary_coupling_debug = cfg.surface_boundary_coupling_debug && capture_ref;
        pm_cfg.boundary_redistribution_diagnostics = cfg.surface_boundary_redistribution_diagnostics && capture_ref;
        pm_cfg.boundary_redistribution_debug = cfg.surface_boundary_redistribution_debug && capture_ref;
        pm_cfg.ssa_random_scatter = cfg.surface_ssa_random_scatter;
        pm_cfg.ssa_scatter_scale = cfg.surface_ssa_scatter_scale;
        pm_cfg.ssa_seed_offset = cfg.surface_ssa_seed_offset;
        pm_cfg.ssa_kernel_mode = cfg.surface_ssa_kernel_mode;
        pm_cfg.ssa_geometry_source_id = cfg.surface_ssa_geometry_source_id;
        pm_cfg.ssa_kz_branch = cfg.surface_ssa_kz_branch;
        pm_cfg.ssa_conv_padding = cfg.surface_ssa_conv_padding;
        pm_cfg.frequency_index = ifq;
        pm_cfg.tx_xyz = [cfg.x_tx, cfg.y_tx, cfg.z_tx];
        pm_cfg.rx_xyz = [rx_state_used.x_rx, rx_state_used.y_rx, rx_state_used.z_rx];
        pm_cfg.z_surface = 0;

        [psi_surface_inc, bubble_surface_meta] = local_march_field( ...
            psi_init_cpu, cfg.z_tx, 0, cfg, x, y, alpha_xy_work, kappa2_work, k0, f_hz, use_gpu);
        [surface_f, delta_f, psi_ref_f, rough_meta_f] = ...
            pm_surface_boundary_model(psi_surface_inc, KX, KY, x, y, cfg.xw, cfg.yw, lambda_f, pm_cfg);

        [psi_ref_at_rx, bubble_rx_meta] = local_march_field( ...
            psi_ref_f, 0, rx_state_used.z_rx, cfg, x, y, alpha_xy_work, kappa2_work, k0, f_hz, use_gpu);
        h_reflect_fi = psi_ref_at_rx(iy_rx, ix_rx);

        if capture_ref
            surface_elevation = surface_f;
            delta_phi = delta_f;
            psi_ref = psi_ref_f;
            roughness_meta = rough_meta_f;
            roughness_meta.reflection_model = 'two_segment_tx_surface_rx';
            bubble_meta.reflection_tx_surface = bubble_surface_meta;
            bubble_meta.reflection_surface_rx = bubble_rx_meta;
        end
    end

    H_direct_f(ifq) = h_direct_fi;
    H_reflect_f(ifq) = h_reflect_fi;
    H_f(ifq) = h_direct_fi + h_reflect_fi;

    if capture_ref
        if save_slice || cfg.show_figures
            if use_gpu
                psifinal_xy = gather(psi_end);
            else
                psifinal_xy = psi_end;
            end
        else
            psifinal_xy = [];
        end

        z_track = z_track_f;
        A_center = A_center_f;
        R_center = R_center_f;
        Axz = Axz_f;
        Ayz = Ayz_f;
        psiout = psiout_f;
        [fit_slope, fit_err_rms, pass_1_over_R, fit_mask] = ...
            local_validate_one_over_R(A_center, R_center, lambda_f, path_span_used);
    end
end

if ~isfield(roughness_meta, 'enabled') || ~roughness_meta.enabled
    if isempty(psifinal_xy)
        surface_elevation = [];
        delta_phi = [];
        psi_ref = [];
    else
        surface_elevation = zeros(size(psifinal_xy));
        delta_phi = zeros(size(psifinal_xy));
        psi_ref = complex(zeros(size(psifinal_xy)));
    end
    roughness_meta = local_disabled_roughness_meta(cfg);
end

h_direct = H_direct_f(idx_f_ref);
h_reflect = H_reflect_f(idx_f_ref);
h_total = H_f(idx_f_ref);

end

function meta = local_disabled_roughness_meta(cfg)
meta = struct( ...
    'enabled', false, ...
    'reflection_coeff_used', cfg.surface_reflect_coeff, ...
    'phase_mode_used', cfg.surface_phase_mode, ...
    'phase_factor_stats', struct('min', NaN, 'max', NaN, 'mean', NaN), ...
    'phi2d_transform_error', struct('max_rel', NaN, 'mean_rel', NaN), ...
    'boundary_model', cfg.surface_boundary_model, ...
    'boundary_operator_form', 'not_executed', ...
    'boundary_dense_matrix_used', false, ...
    'boundary_fft_convention', local_surface_boundary_fft_convention(), ...
    'boundary_equivalence_error', struct( ...
        'enabled', false, ...
        'reference_model', '', ...
        'comparison_model', '', ...
        'max_abs', NaN, ...
        'rel_l2', NaN, ...
        'tol', cfg.surface_boundary_equivalence_tol, ...
        'passed', false), ...
    'boundary_coupling_diagnostics', local_disabled_boundary_coupling_diagnostics(), ...
    'boundary_coupling_debug', struct(), ...
    'boundary_redistribution_diagnostics', local_disabled_boundary_redistribution_diagnostics(), ...
    'boundary_redistribution_debug', struct(), ...
    'boundary_debug_stats', struct(), ...
    'ssa_stat_kernel_meta', local_disabled_ssa_stat_kernel_meta(cfg));
end

function meta = local_disabled_ssa_stat_kernel_meta(cfg)
meta = struct( ...
    'enabled', false, ...
    'model', 'off', ...
    'surface_realization_generated', false, ...
    'random_scatter_enabled', false, ...
    'kernel_mode', cfg.surface_ssa_kernel_mode, ...
    'geometry_source_id', cfg.surface_ssa_geometry_source_id, ...
    'kz_branch', cfg.surface_ssa_kz_branch, ...
    'conv_padding', cfg.surface_ssa_conv_padding, ...
    'conv_operator', 'not_executed', ...
    'conv_padding_size', [NaN, NaN], ...
    'conv_crop_start_index', [NaN, NaN], ...
    'conv_crop_end_index', [NaN, NaN], ...
    'conv_crop_rule', 'not_executed', ...
    'formula_source', 'not_executed', ...
    'boundary_condition', 'not_executed', ...
    'evanescent_included', false, ...
    'sigma_eta_m', NaN, ...
    'Hs_target_m', cfg.sea_hs_target, ...
    'R_coh', complex(NaN, NaN), ...
    'R_coh_raw', complex(NaN, NaN), ...
    'coherent_exponent', NaN, ...
    'coherent_gamma_sum_eff_rad_per_m', NaN, ...
    'coherent_vertical_factor_eff', NaN, ...
    'coherent_reflection_formula', 'not_executed', ...
    'phase_factor_eff', NaN, ...
    'phase_factor_eff_legacy_note', 'not_executed', ...
    'W_eta_variance_target', NaN, ...
    'W_eta_variance_discrete', NaN, ...
    'W_eta_variance_rel_error', NaN, ...
    'P_sca', struct('min', NaN, 'max', NaN, 'mean', NaN, 'std', NaN, 'sum', NaN, 'finite_count', 0), ...
    'E_inc', NaN, ...
    'E_coh', NaN, ...
    'E_sca_raw', NaN, ...
    'E_sca_limited', NaN, ...
    'E_sca', NaN, ...
    'E_ref', NaN, ...
    'incident_spectrum_stats', local_empty_spectrum_stats(), ...
    'coherent_spectrum_stats', local_empty_spectrum_stats(), ...
    'scatter_power_spectrum_stats', local_empty_spectrum_stats(), ...
    'reflected_spectrum_stats', local_empty_spectrum_stats(), ...
    'energy_scale_applied', NaN, ...
    'energy_limit_applied', false, ...
    'energy_conservation_error', NaN, ...
    'propagating_bin_fraction', NaN, ...
    'seed_ssa', NaN, ...
    'seed_offset', cfg.surface_ssa_seed_offset, ...
    'kernel_detail', struct(), ...
    'kernel_formula', 'not_executed', ...
    'limitations', 'not_executed');
end

function stats = local_empty_spectrum_stats()
stats = struct( ...
    'total_energy', NaN, ...
    'centroid_kx_rad_per_m', NaN, ...
    'centroid_ky_rad_per_m', NaN, ...
    'rms_delta_k_rad_per_m', NaN, ...
    'energy_radius_90_rad_per_m', NaN);
end

function diag = local_disabled_boundary_coupling_diagnostics()
diag = struct( ...
    'enabled', false, ...
    'diagnostic_target', 'Kirchhoff boundary screen G_xy', ...
    'total_spectral_energy', NaN, ...
    'zero_wavenumber_energy', NaN, ...
    'nonzero_wavenumber_energy', NaN, ...
    'nonzero_power_fraction', NaN, ...
    'rms_delta_k_rad_per_m', NaN, ...
    'rms_delta_k_nonzero_rad_per_m', NaN, ...
    'energy_radius_50_rad_per_m', NaN, ...
    'energy_radius_90_rad_per_m', NaN, ...
    'energy_radius_95_rad_per_m', NaN, ...
    'nonzero_energy_radius_90_rad_per_m', NaN, ...
    'interpretation', 'not_executed');
end

function diag = local_disabled_boundary_redistribution_diagnostics()
diag = struct( ...
    'enabled', false, ...
    'diagnostic_target', 'incident and reflected Kirchhoff boundary spectra', ...
    'incident_total_energy', NaN, ...
    'incident_centroid_kx_rad_per_m', NaN, ...
    'incident_centroid_ky_rad_per_m', NaN, ...
    'incident_rms_delta_k_rad_per_m', NaN, ...
    'incident_energy_radius_90_rad_per_m', NaN, ...
    'reflect_total_energy', NaN, ...
    'reflect_centroid_kx_rad_per_m', NaN, ...
    'reflect_centroid_ky_rad_per_m', NaN, ...
    'reflect_rms_delta_k_rad_per_m', NaN, ...
    'reflect_energy_radius_90_rad_per_m', NaN, ...
    'rms_delta_k_increase_rad_per_m', NaN, ...
    'centroid_shift_kx_rad_per_m', NaN, ...
    'centroid_shift_ky_rad_per_m', NaN, ...
    'centroid_shift_mag_rad_per_m', NaN, ...
    'high_k_threshold_rad_per_m', NaN, ...
    'reflect_high_k_fraction', NaN, ...
    'flat_ref_rms_delta_k_rad_per_m', NaN, ...
    'flat_ref_energy_radius_90_rad_per_m', NaN, ...
    'flat_ref_high_k_fraction', NaN, ...
    'rough_vs_flat_rms_delta_k_increase_rad_per_m', NaN, ...
    'rough_vs_flat_energy_radius_90_increase_rad_per_m', NaN, ...
    'rough_vs_flat_high_k_fraction_increase', NaN, ...
    'interpretation', 'not_executed');
end

function convention = local_surface_boundary_fft_convention()
convention = ['MATLAB fft2 forward transform is unnormalized; ifft2 includes 1/(nx*ny). ', ...
    'Transverse wavenumbers use [0:N/2-1,-N/2:-1]*2*pi/L in x and y. ', ...
    'The implicit boundary operator is represented by circular convolution on the periodic transverse grid.'];
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

function [psi_end, march_meta] = local_march_field( ...
    psi_start, z_start, z_end, cfg, x, y, alpha_xy, kappa2, k0, f_hz, use_gpu)
%LOCAL_MARCH_FIELD March one complex field between two depths.

march_meta = struct('enabled', false, 'model', 'off');
if abs(z_end - z_start) <= eps
    psi_end = psi_start;
    return
end

n_step = ceil(abs(z_end - z_start) / (cfg.stepz_lamb * (2*pi/k0)));
n_step = max(1, n_step);
dz_step_local = (z_end - z_start) / n_step;
ds_local = abs(dz_step_local);

denom = sqrt(complex(k0^2 - kappa2, 0)) + k0;
fr_local = exp(-1i * 0.5 * ds_local * kappa2 ./ denom);

if use_gpu
    if ~isa(alpha_xy, 'gpuArray')
        alpha_xy = gpuArray(alpha_xy);
    end
    if ~isa(kappa2, 'gpuArray')
        kappa2 = gpuArray(kappa2);
    end
    if ~isa(fr_local, 'gpuArray')
        fr_local = gpuArray(fr_local);
    end
    if ~isa(psi_start, 'gpuArray')
        psi_start = gpuArray(psi_start);
    end
end

psi_k = fft2(psi_start);
march_path_meta = struct([]);
for jj = 1:n_step
    z_curr = z_start + jj * dz_step_local;
    [screen, step_meta] = local_phase_screen( ...
        x, y, z_curr, f_hz, cfg, alpha_xy, k0, ds_local, use_gpu);
    march_path_meta = local_update_bubble_path_meta(march_path_meta, step_meta);
    psi_k = fr_local .* fft2(screen .* ifft2(fr_local .* psi_k));
end

psi_end = ifft2(psi_k);
if use_gpu
    psi_end = gather(psi_end);
end
march_meta = march_path_meta;
end

function [screen, bubble_step_meta] = local_phase_screen( ...
    x, y, z_curr, f_hz, cfg, alpha_xy, k0, ds, use_gpu)
c_bg = local_sound_speed(z_curr, cfg);
if ~isfinite(c_bg) || c_bg <= 0
    error('Invalid local sound speed at z=%g m.', z_curr);
end

[c_eff_xy, alpha_bub_xy, bubble_step_meta] = ...
    bubble_environment_vertical(x, y, z_curr, f_hz, c_bg, cfg);

if ~(isnumeric(c_eff_xy) && all(isfinite(c_eff_xy(:))) && all(c_eff_xy(:) > 0))
    error('bubble_environment_vertical returned invalid c_eff_xy at z=%g m.', z_curr);
end
if ~(isnumeric(alpha_bub_xy) && all(isfinite(alpha_bub_xy(:))) && all(alpha_bub_xy(:) >= 0))
    error('bubble_environment_vertical returned invalid alpha_bub_xy at z=%g m.', z_curr);
end
if ~isequal(size(c_eff_xy), [numel(y), numel(x)])
    error('c_eff_xy must have size [numel(y), numel(x)].');
end
if ~isequal(size(alpha_bub_xy), [numel(y), numel(x)])
    error('alpha_bub_xy must have size [numel(y), numel(x)].');
end

if use_gpu
    if ~isa(alpha_xy, 'gpuArray')
        alpha_xy = gpuArray(alpha_xy);
    end
    c_eff_xy = gpuArray(c_eff_xy);
    alpha_bub_xy = gpuArray(alpha_bub_xy);
end

alpha_total_xy = alpha_xy + alpha_bub_xy;
U_real_xy = (c_eff_xy - cfg.c0) ./ c_eff_xy;
screen = exp(-1i * k0 * ds * (U_real_xy - 1i * alpha_total_xy / k0));
end

function meta_accum = local_update_bubble_path_meta(meta_accum, step_meta)
if isempty(meta_accum)
    meta_accum = step_meta;
    meta_accum.n_step = 1;
    meta_accum.z_min_m = step_meta.z_m;
    meta_accum.z_max_m = step_meta.z_m;
    return
end

n_old = meta_accum.n_step;
n_new = n_old + 1;
meta_accum.n_step = n_new;
meta_accum.z_min_m = min(meta_accum.z_min_m, step_meta.z_m);
meta_accum.z_max_m = max(meta_accum.z_max_m, step_meta.z_m);
meta_accum.z_m = step_meta.z_m;
meta_accum.alpha_bub_np_per_m = step_meta.alpha_bub_np_per_m;
meta_accum.delta_c_bub_mps = step_meta.delta_c_bub_mps;
meta_accum.c_eff_stats = local_merge_step_stats(meta_accum.c_eff_stats, step_meta.c_eff_stats, n_old);
meta_accum.alpha_bub_stats = local_merge_step_stats(meta_accum.alpha_bub_stats, step_meta.alpha_bub_stats, n_old);
if isfield(meta_accum, 'beta_stats') && isfield(step_meta, 'beta_stats')
    meta_accum.beta_stats = local_merge_step_stats(meta_accum.beta_stats, step_meta.beta_stats, n_old);
end
if isfield(meta_accum, 'resonance_radius_stats') && isfield(step_meta, 'resonance_radius_stats')
    meta_accum.resonance_radius_stats = local_merge_step_stats( ...
        meta_accum.resonance_radius_stats, step_meta.resonance_radius_stats, n_old);
end
if ~isempty(step_meta.warning_flags)
    meta_accum.warning_flags = unique([meta_accum.warning_flags, step_meta.warning_flags]);
end
end

function stats = local_merge_step_stats(stats, step_stats, n_old)
stats.min = min(stats.min, step_stats.min);
stats.max = max(stats.max, step_stats.max);
stats.mean = (stats.mean * n_old + step_stats.mean) / (n_old + 1);
stats.std = NaN;
stats.finite_count = stats.finite_count + step_stats.finite_count;
end

function [f_axis, idx_f_ref] = local_resolve_frequency_axis(cfg, rx_state_used)
if numel(cfg.f0) > 1
    f_axis = cfg.f0(:).';
else
    if cfg.enable_wideband
        f_min = cfg.f_band_hz(1);
        f_max = cfg.f_band_hz(2);
        delta_tau = local_estimate_delay_spread(cfg, rx_state_used);
        delta_f_target = 0.9 / max(delta_tau, eps);
        Nf = ceil((f_max - f_min) / delta_f_target) + 1;
        Nf = max(cfg.Nf_min, min(cfg.Nf_max, Nf));
        f_axis = linspace(f_min, f_max, Nf);
    else
        f_axis = cfg.f0;
    end
end

if any(~isfinite(f_axis)) || any(f_axis <= 0)
    error('Frequency axis must contain positive finite values.');
end
f_axis = unique(f_axis(:).');
if isempty(f_axis)
    error('Frequency axis is empty.');
end

f_ref = cfg.f_ref_hz;
if isempty(f_ref) || ~isfinite(f_ref) || f_ref <= 0
    f_ref = 0.5 * (f_axis(1) + f_axis(end));
end
[~, idx_f_ref] = min(abs(f_axis - f_ref));
end

function delta_tau = local_estimate_delay_spread(cfg, rx_state_used)
dx = rx_state_used.x_rx - cfg.x_tx;
dy = rx_state_used.y_rx - cfg.y_tx;
dz = cfg.z_tx - rx_state_used.z_rx;
L_direct = sqrt(dx.^2 + dy.^2 + dz.^2);
tau_direct = L_direct / cfg.c0;

if cfg.enable_surface_reflection
    dz_img = cfg.z_tx + rx_state_used.z_rx;
    L_reflect = sqrt(dx.^2 + dy.^2 + dz_img.^2);
    tau_reflect = L_reflect / cfg.c0;
    delta_tau = abs(tau_reflect - tau_direct);
else
    delta_tau = 1e-6;
end

delta_tau = max(delta_tau, 1e-6);
end

function alpha_xy = local_absorption_profile(x, y, xw, yw, sponge_ratio, alpha_max)
Lsponge_x = sponge_ratio * xw;
Lsponge_y = sponge_ratio * yw;
if Lsponge_x <= 0 || Lsponge_y <= 0
    error('Sponge thickness must be positive.');
end

x_start = 0.5 * xw - Lsponge_x;
y_start = 0.5 * yw - Lsponge_y;
x_abs = abs(x);
y_abs = abs(y);

alpha_x = zeros(size(x_abs));
mask_x = (x_abs > x_start);
alpha_x(mask_x) = alpha_max * ((x_abs(mask_x) - x_start) / Lsponge_x).^3;
alpha_x = min(max(alpha_x, 0), alpha_max);

alpha_y = zeros(size(y_abs));
mask_y = (y_abs > y_start);
alpha_y(mask_y) = alpha_max * ((y_abs(mask_y) - y_start) / Lsponge_y).^3;
alpha_y = min(max(alpha_y, 0), alpha_max);

alpha_xy = alpha_y(:) * ones(1, numel(x_abs)) + ones(numel(y_abs), 1) * alpha_x(:).';
end

function tf = local_has_gpu()
tf = false;
try
    tf = (gpuDeviceCount("available") > 0);
catch
    tf = false;
end
end
