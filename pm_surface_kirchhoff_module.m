function [surface_elevation, delta_phi, psi_ref, meta] = ...
    pm_surface_kirchhoff_module(psi_inc, KX, KY, x, y, xw, yw, lambda0, pm_cfg)
%PM_SURFACE_KIRCHHOFF_MODULE
% Task 1: 2D Pierson-Moskowitz rough surface synthesis on (x,y) grid.
% Task 2: Kirchhoff phase distortion on incident field psi_inc.
% Task 3: Sanity-check visualization.

if nargin < 9 || isempty(pm_cfg)
    pm_cfg = struct();
end
if ~isfield(pm_cfg, 'U') || isempty(pm_cfg.U)
    pm_cfg.U = 5.0;
end
if ~isfield(pm_cfg, 'Hs_target') || isempty(pm_cfg.Hs_target)
    pm_cfg.Hs_target = 0.5;
end
if ~isfield(pm_cfg, 'seed') || isempty(pm_cfg.seed)
    pm_cfg.seed = 12345;
end
if ~isfield(pm_cfg, 'show_figure') || isempty(pm_cfg.show_figure)
    pm_cfg.show_figure = true;
end
if ~isfield(pm_cfg, 'reflect_coeff') || isempty(pm_cfg.reflect_coeff)
    pm_cfg.reflect_coeff = -1;
end
if ~isfield(pm_cfg, 'phase_mode') || isempty(pm_cfg.phase_mode)
    pm_cfg.phase_mode = 'normal';
end
if ~isfield(pm_cfg, 'oblique_clip') || isempty(pm_cfg.oblique_clip)
    pm_cfg.oblique_clip = [0, 1];
end
if ~isfield(pm_cfg, 'tx_xyz') || isempty(pm_cfg.tx_xyz)
    pm_cfg.tx_xyz = [0, 0, NaN];
end
if ~isfield(pm_cfg, 'rx_xyz') || isempty(pm_cfg.rx_xyz)
    pm_cfg.rx_xyz = [0, 0, NaN];
end
if ~isfield(pm_cfg, 'z_surface') || isempty(pm_cfg.z_surface)
    pm_cfg.z_surface = 0;
end
if ~isfield(pm_cfg, 'boundary_model') || isempty(pm_cfg.boundary_model)
    pm_cfg.boundary_model = 'kirchhoff_spatial';
end
if ~isfield(pm_cfg, 'boundary_check_equivalence') || isempty(pm_cfg.boundary_check_equivalence)
    pm_cfg.boundary_check_equivalence = false;
end
if ~isfield(pm_cfg, 'boundary_debug') || isempty(pm_cfg.boundary_debug)
    pm_cfg.boundary_debug = false;
end
if ~isfield(pm_cfg, 'boundary_equivalence_tol') || isempty(pm_cfg.boundary_equivalence_tol)
    pm_cfg.boundary_equivalence_tol = 1e-10;
end
if ~isfield(pm_cfg, 'boundary_coupling_diagnostics') || isempty(pm_cfg.boundary_coupling_diagnostics)
    pm_cfg.boundary_coupling_diagnostics = false;
end
if ~isfield(pm_cfg, 'boundary_coupling_debug') || isempty(pm_cfg.boundary_coupling_debug)
    pm_cfg.boundary_coupling_debug = false;
end
if ~isfield(pm_cfg, 'boundary_redistribution_diagnostics') || isempty(pm_cfg.boundary_redistribution_diagnostics)
    pm_cfg.boundary_redistribution_diagnostics = false;
end
if ~isfield(pm_cfg, 'boundary_redistribution_debug') || isempty(pm_cfg.boundary_redistribution_debug)
    pm_cfg.boundary_redistribution_debug = false;
end

if ~(isscalar(pm_cfg.reflect_coeff) && isnumeric(pm_cfg.reflect_coeff) && isfinite(pm_cfg.reflect_coeff))
    error('pm_cfg.reflect_coeff must be a finite scalar (real or complex).');
end
if isstring(pm_cfg.phase_mode)
    pm_cfg.phase_mode = char(pm_cfg.phase_mode);
end
phase_mode = lower(pm_cfg.phase_mode);
if ~ischar(phase_mode) || (~strcmp(phase_mode, 'normal') && ~strcmp(phase_mode, 'oblique'))
    error('pm_cfg.phase_mode must be ''normal'' or ''oblique''.');
end
if ~(isnumeric(pm_cfg.oblique_clip) && numel(pm_cfg.oblique_clip) == 2 && all(isfinite(pm_cfg.oblique_clip(:))))
    error('pm_cfg.oblique_clip must be a finite 1x2 numeric range.');
end
clip_pair = sort(pm_cfg.oblique_clip(:).');
if clip_pair(1) < 0 || clip_pair(2) > 1
    error('pm_cfg.oblique_clip must stay within [0,1].');
end
if ~(isnumeric(pm_cfg.tx_xyz) && numel(pm_cfg.tx_xyz) == 3 && all(isfinite(pm_cfg.tx_xyz(:))))
    error('pm_cfg.tx_xyz must be a finite [x,y,z].');
end
if ~(isnumeric(pm_cfg.rx_xyz) && numel(pm_cfg.rx_xyz) == 3 && all(isfinite(pm_cfg.rx_xyz(:))))
    error('pm_cfg.rx_xyz must be a finite [x,y,z].');
end
boundary_model = local_normalize_choice(pm_cfg.boundary_model, 'pm_cfg.boundary_model');
if ~any(strcmp(boundary_model, {'kirchhoff_spatial', 'kirchhoff_kdomain'}))
    error('pm_cfg.boundary_model must be ''kirchhoff_spatial'' or ''kirchhoff_kdomain''.');
end
if ~isscalar(pm_cfg.boundary_check_equivalence)
    error('pm_cfg.boundary_check_equivalence must be a scalar logical flag.');
end
boundary_check_equivalence = logical(pm_cfg.boundary_check_equivalence);
if ~isscalar(pm_cfg.boundary_debug)
    error('pm_cfg.boundary_debug must be a scalar logical flag.');
end
boundary_debug = logical(pm_cfg.boundary_debug);
if ~(isscalar(pm_cfg.boundary_equivalence_tol) && isnumeric(pm_cfg.boundary_equivalence_tol) && ...
     isfinite(pm_cfg.boundary_equivalence_tol) && pm_cfg.boundary_equivalence_tol > 0)
    error('pm_cfg.boundary_equivalence_tol must be a positive finite scalar.');
end
boundary_equivalence_tol = pm_cfg.boundary_equivalence_tol;
if ~isscalar(pm_cfg.boundary_coupling_diagnostics)
    error('pm_cfg.boundary_coupling_diagnostics must be a scalar logical flag.');
end
boundary_coupling_diagnostics = logical(pm_cfg.boundary_coupling_diagnostics);
if ~isscalar(pm_cfg.boundary_coupling_debug)
    error('pm_cfg.boundary_coupling_debug must be a scalar logical flag.');
end
boundary_coupling_debug = logical(pm_cfg.boundary_coupling_debug);
if ~isscalar(pm_cfg.boundary_redistribution_diagnostics)
    error('pm_cfg.boundary_redistribution_diagnostics must be a scalar logical flag.');
end
boundary_redistribution_diagnostics = logical(pm_cfg.boundary_redistribution_diagnostics);
if ~isscalar(pm_cfg.boundary_redistribution_debug)
    error('pm_cfg.boundary_redistribution_debug must be a scalar logical flag.');
end
boundary_redistribution_debug = logical(pm_cfg.boundary_redistribution_debug);

% -----------------------------
% 1) PM rough surface generation
% -----------------------------
U = pm_cfg.U;          % m/s
g = 9.81;              % m/s^2
alpha_PM = 8.10e-3;
beta_PM = 0.74;

K = sqrt(KX.^2 + KY.^2);
E1D = zeros(size(K));
Phi2D = zeros(size(K));

mask = (K > 0);
K_nonzero = K(mask);
E1D(mask) = (alpha_PM ./ (2 .* K_nonzero.^3)) .* ...
            exp(-beta_PM * (g^2) ./ (U^4 .* K_nonzero.^2));
Phi2D(mask) = E1D(mask) ./ (2*pi*K_nonzero);

dkx = 2*pi / xw;
dky = 2*pi / yw;

A = sqrt(Phi2D .* dkx .* dky);
rng(pm_cfg.seed, 'twister')
N = (randn(size(K)) + 1i*randn(size(K))) / sqrt(2);
Zk = A .* N;

eta_raw = real(ifft2(Zk));
eta_raw = eta_raw * numel(eta_raw); % compensate MATLAB ifft2 normalization

% Calibrate to realistic PM sea-state scale (~0.5 m significant wave height at U=5 m/s).
Hs_target = pm_cfg.Hs_target; % m
Hs_raw = 4 * std(eta_raw(:));
if Hs_raw > 0
    scale_factor = Hs_target / Hs_raw;
else
    scale_factor = 0;
end
surface_elevation = eta_raw * scale_factor;

% -----------------------------
% 2) Kirchhoff phase distortion
% -----------------------------
k0 = 2*pi / lambda0;
[X, Y] = meshgrid(x, y);

switch phase_mode
    case 'normal'
        phase_factor = 2 * ones(size(surface_elevation));
    case 'oblique'
        tx_xyz = pm_cfg.tx_xyz(:).';
        rx_xyz = pm_cfg.rx_xyz(:).';
        z_surface = pm_cfg.z_surface;

        dz_i = tx_xyz(3) - z_surface;
        dz_r = rx_xyz(3) - z_surface;
        if dz_i <= 0 || dz_r < 0
            error('For oblique mode, tx/rx must be at or below the surface: z >= z_surface.');
        end

        dx_i = X - tx_xyz(1);
        dy_i = Y - tx_xyz(2);
        Ri = sqrt(dx_i.^2 + dy_i.^2 + dz_i.^2);
        cos_i = dz_i ./ max(Ri, eps);

        dx_r = X - rx_xyz(1);
        dy_r = Y - rx_xyz(2);
        Rr = sqrt(dx_r.^2 + dy_r.^2 + dz_r.^2);
        cos_r = dz_r ./ max(Rr, eps);

        cos_i = min(max(cos_i, clip_pair(1)), clip_pair(2));
        cos_r = min(max(cos_r, clip_pair(1)), clip_pair(2));
        phase_factor = cos_i + cos_r;
    otherwise
        error('Unsupported phase_mode: %s', phase_mode);
end

delta_phi = 2 * k0 * surface_elevation .* phase_factor;
G_xy = pm_cfg.reflect_coeff .* exp(1i * delta_phi);
[coupling_diag, coupling_debug] = local_boundary_coupling_diagnostics( ...
    G_xy, KX, KY, boundary_coupling_diagnostics, boundary_coupling_debug);
[psi_ref, boundary_meta] = local_apply_kirchhoff_boundary( ...
    psi_inc, G_xy, boundary_model, boundary_check_equivalence, ...
    boundary_debug, boundary_equivalence_tol);
[redistribution_diag, redistribution_debug] = local_boundary_redistribution_diagnostics( ...
    psi_inc, psi_ref, pm_cfg.reflect_coeff .* psi_inc, KX, KY, ...
    boundary_redistribution_diagnostics, boundary_redistribution_debug);

% -----------------------------
% 3) Sanity-check visualization
% -----------------------------
if pm_cfg.show_figure
    figure(15); clf
    subplot(1,2,1)
    surf(X, Y, surface_elevation, 'EdgeColor', 'none')
    view(40, 35)
    axis tight
    xlabel('x (m)')
    ylabel('y (m)')
    zlabel('\xi(x,y) (m)')
    title('PM rough sea surface elevation')
    colorbar

    subplot(1,2,2)
    imagesc(x, y, angle(psi_ref))
    axis xy tight
    xlabel('x (m)')
    ylabel('y (m)')
    title(sprintf('Phase of reflected field angle(\\psi_{ref}), mode=%s', phase_mode))
    colorbar
end

E_from_Phi = zeros(size(K));
E_from_Phi(mask) = Phi2D(mask) .* 2*pi .* K(mask);
rel_err = abs(E_from_Phi(mask) - E1D(mask)) ./ max(abs(E1D(mask)), eps);

meta = struct();
meta.U = U;
meta.g = g;
meta.alpha_PM = alpha_PM;
meta.beta_PM = beta_PM;
meta.dkx = dkx;
meta.dky = dky;
meta.spectrum_definition = 'E1D(K)=alpha/(2K^3)exp(-beta g^2/(U^4 K^2)); Phi2D(Kx,Ky)=E1D(K)/(2*pi*K)';
meta.E1D = local_field_stats(E1D(mask));
meta.Phi2D = local_field_stats(Phi2D(mask));
meta.phi2d_transform_error = struct( ...
    'max_rel', max(rel_err(:), [], 'omitnan'), ...
    'mean_rel', mean(rel_err(:), 'omitnan'));
meta.Hs_target = Hs_target;
meta.Hs_raw = Hs_raw;
meta.Hs_scaled = 4 * std(surface_elevation(:));
meta.scale_factor = scale_factor;
meta.seed = pm_cfg.seed;
meta.reflection_coeff_used = pm_cfg.reflect_coeff;
meta.phase_mode_used = phase_mode;
meta.phase_factor_stats = struct( ...
    'min', min(phase_factor(:)), ...
    'max', max(phase_factor(:)), ...
    'mean', mean(phase_factor(:)));
meta.oblique_clip_used = clip_pair;
meta.show_figure = logical(pm_cfg.show_figure);
meta.enabled = true;
meta.boundary_model = boundary_meta.boundary_model;
meta.boundary_operator_form = boundary_meta.boundary_operator_form;
meta.boundary_dense_matrix_used = boundary_meta.boundary_dense_matrix_used;
meta.boundary_fft_convention = boundary_meta.boundary_fft_convention;
meta.boundary_equivalence_error = boundary_meta.boundary_equivalence_error;
meta.boundary_coupling_diagnostics = coupling_diag;
meta.boundary_coupling_debug = coupling_debug;
meta.boundary_redistribution_diagnostics = redistribution_diag;
meta.boundary_redistribution_debug = redistribution_debug;
meta.boundary_debug_stats = boundary_meta.boundary_debug_stats;

end

function [diag, debug] = local_boundary_redistribution_diagnostics( ...
    psi_inc_xy, psi_ref_xy, psi_flat_ref_xy, KX, KY, enabled, debug_enabled)
diag = local_disabled_redistribution_diagnostics();
debug = struct();
if ~enabled
    return
end

Psi_inc_k = fft2(psi_inc_xy);
Psi_ref_k = fft2(psi_ref_xy);
Psi_flat_ref_k = fft2(psi_flat_ref_xy);
P_inc = abs(Psi_inc_k).^2;
P_ref = abs(Psi_ref_k).^2;
P_flat = abs(Psi_flat_ref_k).^2;
DeltaK = sqrt(KX.^2 + KY.^2);

inc_stats = local_spectrum_redistribution_stats(P_inc, KX, KY, DeltaK);
ref_stats = local_spectrum_redistribution_stats(P_ref, KX, KY, DeltaK);
flat_stats = local_spectrum_redistribution_stats(P_flat, KX, KY, DeltaK);

high_k_threshold = inc_stats.energy_radius_90_rad_per_m;
reflect_high_k_fraction = local_high_k_fraction(P_ref, DeltaK, high_k_threshold);
flat_high_k_fraction = local_high_k_fraction(P_flat, DeltaK, high_k_threshold);

centroid_shift_kx = ref_stats.centroid_kx_rad_per_m - inc_stats.centroid_kx_rad_per_m;
centroid_shift_ky = ref_stats.centroid_ky_rad_per_m - inc_stats.centroid_ky_rad_per_m;

diag = struct( ...
    'enabled', true, ...
    'diagnostic_target', 'incident and reflected Kirchhoff boundary spectra', ...
    'spectrum_quantity', 'P_k = abs(fft2(psi_xy)).^2', ...
    'incident_total_energy', inc_stats.total_energy, ...
    'incident_centroid_kx_rad_per_m', inc_stats.centroid_kx_rad_per_m, ...
    'incident_centroid_ky_rad_per_m', inc_stats.centroid_ky_rad_per_m, ...
    'incident_rms_delta_k_rad_per_m', inc_stats.rms_delta_k_rad_per_m, ...
    'incident_energy_radius_90_rad_per_m', inc_stats.energy_radius_90_rad_per_m, ...
    'reflect_total_energy', ref_stats.total_energy, ...
    'reflect_centroid_kx_rad_per_m', ref_stats.centroid_kx_rad_per_m, ...
    'reflect_centroid_ky_rad_per_m', ref_stats.centroid_ky_rad_per_m, ...
    'reflect_rms_delta_k_rad_per_m', ref_stats.rms_delta_k_rad_per_m, ...
    'reflect_energy_radius_90_rad_per_m', ref_stats.energy_radius_90_rad_per_m, ...
    'rms_delta_k_increase_rad_per_m', ref_stats.rms_delta_k_rad_per_m - inc_stats.rms_delta_k_rad_per_m, ...
    'centroid_shift_kx_rad_per_m', centroid_shift_kx, ...
    'centroid_shift_ky_rad_per_m', centroid_shift_ky, ...
    'centroid_shift_mag_rad_per_m', hypot(centroid_shift_kx, centroid_shift_ky), ...
    'high_k_threshold_rad_per_m', high_k_threshold, ...
    'reflect_high_k_fraction', reflect_high_k_fraction, ...
    'flat_ref_rms_delta_k_rad_per_m', flat_stats.rms_delta_k_rad_per_m, ...
    'flat_ref_energy_radius_90_rad_per_m', flat_stats.energy_radius_90_rad_per_m, ...
    'flat_ref_high_k_fraction', flat_high_k_fraction, ...
    'rough_vs_flat_rms_delta_k_increase_rad_per_m', ref_stats.rms_delta_k_rad_per_m - flat_stats.rms_delta_k_rad_per_m, ...
    'rough_vs_flat_energy_radius_90_increase_rad_per_m', ref_stats.energy_radius_90_rad_per_m - flat_stats.energy_radius_90_rad_per_m, ...
    'rough_vs_flat_high_k_fraction_increase', reflect_high_k_fraction - flat_high_k_fraction, ...
    'interpretation', ['Compares the current incident angular spectrum with the Kirchhoff-reflected spectrum. ', ...
        'Larger reflected rms_delta_k, centroid shift, and high-k fraction indicate stronger spectrum redistribution for this incident field.'], ...
    'limitations', ['Diagnostic only: these are incident-weighted spectrum moments under the Kirchhoff phase screen, ', ...
        'not T-matrix entries, SSA/NLSSA terms, scattering cross sections, or statistical channel generation.']);

if debug_enabled
    debug = struct( ...
        'high_k_threshold_rad_per_m', high_k_threshold, ...
        'incident', local_redistribution_debug(P_inc, DeltaK, KX, KY), ...
        'reflect', local_redistribution_debug(P_ref, DeltaK, KX, KY), ...
        'flat_ref', local_redistribution_debug(P_flat, DeltaK, KX, KY));
end
end

function stats = local_spectrum_redistribution_stats(P_k, KX, KY, DeltaK)
total_energy = sum(P_k(:));
den = max(total_energy, eps);
stats = struct( ...
    'total_energy', total_energy, ...
    'centroid_kx_rad_per_m', sum(KX(:) .* P_k(:)) / den, ...
    'centroid_ky_rad_per_m', sum(KY(:) .* P_k(:)) / den, ...
    'rms_delta_k_rad_per_m', sqrt(sum((DeltaK(:).^2) .* P_k(:)) / den), ...
    'energy_radius_90_rad_per_m', local_weighted_radius(DeltaK(:), P_k(:), 0.90));
end

function frac = local_high_k_fraction(P_k, DeltaK, threshold)
total_energy = sum(P_k(:));
if ~(isfinite(threshold) && total_energy > 0)
    frac = NaN;
    return
end
frac = sum(P_k(DeltaK > threshold)) / max(total_energy, eps);
end

function diag = local_disabled_redistribution_diagnostics()
diag = struct( ...
    'enabled', false, ...
    'diagnostic_target', 'incident and reflected Kirchhoff boundary spectra', ...
    'spectrum_quantity', 'P_k = abs(fft2(psi_xy)).^2', ...
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
    'interpretation', 'not_executed', ...
    'limitations', 'not_executed');
end

function debug = local_redistribution_debug(P_k, DeltaK, KX, KY)
nbin = 16;
max_delta_k = max(DeltaK(:));
if max_delta_k <= 0 || ~isfinite(max_delta_k)
    bin_edges = zeros(1, nbin + 1);
    radial_energy = zeros(1, nbin);
else
    bin_edges = linspace(0, max_delta_k, nbin + 1);
    radial_energy = zeros(1, nbin);
    for ib = 1:nbin
        if ib == nbin
            mask_bin = DeltaK >= bin_edges(ib) & DeltaK <= bin_edges(ib + 1);
        else
            mask_bin = DeltaK >= bin_edges(ib) & DeltaK < bin_edges(ib + 1);
        end
        radial_energy(ib) = sum(P_k(mask_bin));
    end
end

total_energy = sum(P_k(:));
top_n = min(8, numel(P_k));
[top_energy, top_idx] = maxk(P_k(:), top_n);
debug = struct( ...
    'radial_bin_edges_rad_per_m', bin_edges, ...
    'radial_energy_fraction', radial_energy / max(total_energy, eps), ...
    'top_peak_energy', top_energy(:).', ...
    'top_peak_fraction', top_energy(:).' / max(total_energy, eps), ...
    'top_peak_kx_rad_per_m', KX(top_idx).', ...
    'top_peak_ky_rad_per_m', KY(top_idx).', ...
    'top_peak_delta_k_rad_per_m', DeltaK(top_idx).');
end

function [diag, debug] = local_boundary_coupling_diagnostics(G_xy, KX, KY, enabled, debug_enabled)
diag = local_disabled_coupling_diagnostics();
debug = struct();
if ~enabled
    return
end

G_hat_k = fft2(G_xy);
P_k = abs(G_hat_k).^2;
DeltaK = sqrt(KX.^2 + KY.^2);

total_energy = sum(P_k(:));
zero_energy = P_k(1, 1);
nonzero_energy = max(total_energy - zero_energy, 0);
nonzero_mask = true(size(P_k));
nonzero_mask(1, 1) = false;

nonzero_fraction = nonzero_energy / max(total_energy, eps);
rms_delta_k = sqrt(sum((DeltaK(:).^2) .* P_k(:)) / max(total_energy, eps));
if nonzero_energy > 0
    rms_delta_k_nonzero = sqrt(sum((DeltaK(nonzero_mask).^2) .* P_k(nonzero_mask)) / nonzero_energy);
else
    rms_delta_k_nonzero = NaN;
end

diag = struct( ...
    'enabled', true, ...
    'diagnostic_target', 'Kirchhoff boundary screen G_xy', ...
    'spectrum_quantity', 'P_k = abs(fft2(G_xy)).^2', ...
    'total_spectral_energy', total_energy, ...
    'zero_wavenumber_energy', zero_energy, ...
    'nonzero_wavenumber_energy', nonzero_energy, ...
    'nonzero_power_fraction', nonzero_fraction, ...
    'rms_delta_k_rad_per_m', rms_delta_k, ...
    'rms_delta_k_nonzero_rad_per_m', rms_delta_k_nonzero, ...
    'energy_radius_50_rad_per_m', local_weighted_radius(DeltaK(:), P_k(:), 0.50), ...
    'energy_radius_90_rad_per_m', local_weighted_radius(DeltaK(:), P_k(:), 0.90), ...
    'energy_radius_95_rad_per_m', local_weighted_radius(DeltaK(:), P_k(:), 0.95), ...
    'nonzero_energy_radius_90_rad_per_m', local_weighted_radius( ...
        DeltaK(nonzero_mask), P_k(nonzero_mask), 0.90), ...
    'interpretation', ['Larger nonzero_power_fraction and rms_delta_k indicate a less uniform ', ...
        'Kirchhoff screen and stronger potential off-diagonal Kprime-to-K coupling in the implicit convolution.'], ...
    'limitations', ['Diagnostic only: these are spectrum metrics of the Kirchhoff phase screen, ', ...
        'not T-matrix entries, SSA/NLSSA terms, scattering cross sections, or statistical channel generation.']);

if debug_enabled
    debug = local_boundary_coupling_debug(P_k, DeltaK, KX, KY, nonzero_mask);
end
end

function diag = local_disabled_coupling_diagnostics()
diag = struct( ...
    'enabled', false, ...
    'diagnostic_target', 'Kirchhoff boundary screen G_xy', ...
    'spectrum_quantity', 'P_k = abs(fft2(G_xy)).^2', ...
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
    'interpretation', 'not_executed', ...
    'limitations', 'not_executed');
end

function radius = local_weighted_radius(radius_values, weights, fraction)
radius_values = radius_values(:);
weights = weights(:);
valid = isfinite(radius_values) & isfinite(weights) & weights >= 0;
radius_values = radius_values(valid);
weights = weights(valid);
if isempty(weights) || sum(weights) <= 0
    radius = NaN;
    return
end

[radius_sorted, order] = sort(radius_values);
weights_sorted = weights(order);
cum_weights = cumsum(weights_sorted);
target = max(0, min(1, fraction)) * cum_weights(end);
idx = find(cum_weights >= target, 1, 'first');
radius = radius_sorted(idx);
end

function debug = local_boundary_coupling_debug(P_k, DeltaK, KX, KY, nonzero_mask)
nbin = 16;
max_delta_k = max(DeltaK(:));
if max_delta_k <= 0 || ~isfinite(max_delta_k)
    bin_edges = zeros(1, nbin + 1);
    radial_energy = zeros(1, nbin);
else
    bin_edges = linspace(0, max_delta_k, nbin + 1);
    radial_energy = zeros(1, nbin);
    for ib = 1:nbin
        if ib == nbin
            mask_bin = DeltaK >= bin_edges(ib) & DeltaK <= bin_edges(ib + 1);
        else
            mask_bin = DeltaK >= bin_edges(ib) & DeltaK < bin_edges(ib + 1);
        end
        radial_energy(ib) = sum(P_k(mask_bin));
    end
end

total_energy = sum(P_k(:));
radial_energy_fraction = radial_energy / max(total_energy, eps);

nonzero_power = P_k;
nonzero_power(~nonzero_mask) = -Inf;
top_n = min(8, nnz(nonzero_mask));
[top_energy, top_idx] = maxk(nonzero_power(:), top_n);
finite_top = isfinite(top_energy);
top_energy = top_energy(finite_top);
top_idx = top_idx(finite_top);

debug = struct( ...
    'radial_bin_edges_rad_per_m', bin_edges, ...
    'radial_energy_fraction', radial_energy_fraction, ...
    'top_nonzero_peak_energy', top_energy(:).', ...
    'top_nonzero_peak_fraction', (top_energy(:).' / max(total_energy, eps)), ...
    'top_nonzero_peak_kx_rad_per_m', KX(top_idx).', ...
    'top_nonzero_peak_ky_rad_per_m', KY(top_idx).', ...
    'top_nonzero_peak_delta_k_rad_per_m', DeltaK(top_idx).');
end

function [psi_ref_xy, boundary_meta] = local_apply_kirchhoff_boundary( ...
    psi_inc_xy, G_xy, boundary_model, check_equivalence, debug_flag, equivalence_tol)
%LOCAL_APPLY_KIRCHHOFF_BOUNDARY Apply the phase screen through the selected interface.

boundary_meta = local_boundary_meta_base(boundary_model, equivalence_tol);

switch boundary_model
    case 'kirchhoff_spatial'
        psi_ref_xy = G_xy .* psi_inc_xy;
        boundary_meta.boundary_operator_form = ...
            'spatial_pointwise: psi_ref_xy = G_xy .* psi_inc_xy';

        if check_equivalence
            psi_ref_kdomain_xy = local_apply_kdomain_boundary(psi_inc_xy, G_xy);
            boundary_meta.boundary_equivalence_error = local_boundary_equivalence_error( ...
                psi_ref_xy, psi_ref_kdomain_xy, equivalence_tol, true, ...
                'kirchhoff_spatial', 'kirchhoff_kdomain');
        end

    case 'kirchhoff_kdomain'
        psi_ref_xy = local_apply_kdomain_boundary(psi_inc_xy, G_xy);
        boundary_meta.boundary_operator_form = ...
            'implicit_fft_convolution: Psi_inc_k -> fft2(G_xy .* ifft2(Psi_inc_k)) -> Psi_ref_k';

        if check_equivalence
            psi_ref_spatial_xy = G_xy .* psi_inc_xy;
            boundary_meta.boundary_equivalence_error = local_boundary_equivalence_error( ...
                psi_ref_spatial_xy, psi_ref_xy, equivalence_tol, true, ...
                'kirchhoff_spatial', 'kirchhoff_kdomain');
        end

    otherwise
        error('Unsupported boundary_model: %s', boundary_model);
end

if debug_flag
    Psi_inc_k = fft2(psi_inc_xy);
    Psi_ref_k = fft2(psi_ref_xy);
    G_hat_k = fft2(G_xy);
    boundary_meta.boundary_debug_stats = struct( ...
        'ny', size(psi_inc_xy, 1), ...
        'nx', size(psi_inc_xy, 2), ...
        'psi_inc_norm_l2', norm(psi_inc_xy(:)), ...
        'psi_ref_norm_l2', norm(psi_ref_xy(:)), ...
        'Psi_inc_k_norm_l2', norm(Psi_inc_k(:)), ...
        'Psi_ref_k_norm_l2', norm(Psi_ref_k(:)), ...
        'G_abs_min', min(abs(G_xy(:))), ...
        'G_abs_max', max(abs(G_xy(:))), ...
        'G_hat_norm_l2', norm(G_hat_k(:)));
end
end

function psi_ref_xy = local_apply_kdomain_boundary(psi_inc_xy, G_xy)
Psi_inc_k = fft2(psi_inc_xy);
Psi_ref_k = fft2(G_xy .* ifft2(Psi_inc_k));
psi_ref_xy = ifft2(Psi_ref_k);
end

function boundary_meta = local_boundary_meta_base(boundary_model, equivalence_tol)
boundary_meta = struct( ...
    'boundary_model', boundary_model, ...
    'boundary_operator_form', '', ...
    'boundary_dense_matrix_used', false, ...
    'boundary_fft_convention', local_boundary_fft_convention(), ...
    'boundary_equivalence_error', local_boundary_equivalence_error( ...
        [], [], equivalence_tol, false, '', ''), ...
    'boundary_debug_stats', struct());
end

function err = local_boundary_equivalence_error(ref_xy, test_xy, tol, enabled, ref_model, test_model)
if ~enabled
    err = struct( ...
        'enabled', false, ...
        'reference_model', ref_model, ...
        'comparison_model', test_model, ...
        'max_abs', NaN, ...
        'rel_l2', NaN, ...
        'tol', tol, ...
        'passed', false);
    return
end

diff_xy = test_xy - ref_xy;
max_abs = max(abs(diff_xy(:)));
rel_l2 = norm(diff_xy(:)) / max(norm(ref_xy(:)), eps);
err = struct( ...
    'enabled', true, ...
    'reference_model', ref_model, ...
    'comparison_model', test_model, ...
    'max_abs', max_abs, ...
    'rel_l2', rel_l2, ...
    'tol', tol, ...
    'passed', rel_l2 <= tol);
end

function convention = local_boundary_fft_convention()
convention = ['MATLAB fft2 forward transform is unnormalized; ifft2 includes 1/(nx*ny). ', ...
    'Transverse wavenumbers use [0:N/2-1,-N/2:-1]*2*pi/L in x and y. ', ...
    'The implicit boundary operator is represented by circular convolution on the periodic transverse grid.'];
end

function out = local_normalize_choice(v, name)
if isstring(v)
    v = char(v);
end
if ~ischar(v)
    error('%s must be a string or char.', name);
end
out = lower(v);
end

function stats = local_field_stats(v)
if isempty(v)
    stats = struct('min', NaN, 'max', NaN, 'mean', NaN, 'std', NaN);
    return
end
stats = struct( ...
    'min', min(v(:)), ...
    'max', max(v(:)), ...
    'mean', mean(v(:)), ...
    'std', std(v(:)));
end
