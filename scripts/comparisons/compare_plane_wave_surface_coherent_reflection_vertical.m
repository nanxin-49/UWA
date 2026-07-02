run(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'bootstrap_project.m'));
%COMPARE_PLANE_WAVE_SURFACE_COHERENT_REFLECTION_VERTICAL
% Sea-surface-only coherent reflection diagnostic for vertical plane waves.
%
% This script does not call PE/WAPE and does not enter the communication
% chain. It compares the normal-incidence pressure-release SSA1 coherent
% reflection with Kirchhoff phase-screen ensemble coherent reflection under
% raw, non-Hs-normalized PM spectra.

clear; clc;

cfg = struct();
cfg.wind_list_mps = local_env_numeric_list('PLANE_WAVE_WIND_LIST', [3, 5, 8, 12]);
cfg.f_list_hz = local_env_numeric_list('PLANE_WAVE_F_LIST_HZ', [4000, 6000, 8000, 10000]);
cfg.seed_list = local_env_numeric_list('PLANE_WAVE_SEED_LIST', 12345 + (0:31));
seed_count_override = local_env_scalar('PLANE_WAVE_SEED_COUNT', NaN);
if isfinite(seed_count_override)
    cfg.seed_list = cfg.seed_list(1:min(numel(cfg.seed_list), max(1, round(seed_count_override))));
end
cfg.nx = round(local_env_scalar('PLANE_WAVE_NX', 128));
cfg.ny = round(local_env_scalar('PLANE_WAVE_NY', 128));
cfg.xw_m = local_env_scalar('PLANE_WAVE_XW_M', 50);
cfg.yw_m = local_env_scalar('PLANE_WAVE_YW_M', 50);
cfg.c0_mps = local_env_scalar('PLANE_WAVE_C0_MPS', 1500);
cfg.R0 = -1;
cfg.gamma_legacy_4k_eta_normal = 2;
cfg.output_prefix = local_env_string('PLANE_WAVE_OUTPUT_PREFIX', 'plane_wave');

fprintf('Plane-wave surface coherent reflection diagnostic\n');
fprintf('  winds: %s m/s\n', mat2str(cfg.wind_list_mps));
fprintf('  frequencies: %s Hz\n', mat2str(cfg.f_list_hz));
fprintf('  seeds: %d, grid: %d x %d, aperture: %.3g m x %.3g m\n', ...
    numel(cfg.seed_list), cfg.nx, cfg.ny, cfg.xw_m, cfg.yw_m);

[KX, KY, dkx, dky] = local_wavenumber_grid(cfg.nx, cfg.ny, cfg.xw_m, cfg.yw_m);
flat_check = local_flat_check(cfg);

summary_rows = struct([]);
run_rows = struct([]);
row_idx = 0;
run_idx = 0;

for iw = 1:numel(cfg.wind_list_mps)
    U = cfg.wind_list_mps(iw);
    pm_spec = local_pm_spectrum(KX, KY, cfg.xw_m, cfg.yw_m, U);
    sigma_eta_raw_m = sqrt(max(sum(pm_spec.Phi2D(:)) * dkx * dky, 0));
    Hs_raw_m = 4 * sigma_eta_raw_m;
    variance_rel_error = local_relative_error(sum(pm_spec.Phi2D(:)) * dkx * dky, sigma_eta_raw_m^2);

    eta_stack_std = zeros(numel(cfg.seed_list), 1);
    for ifq = 1:numel(cfg.f_list_hz)
        f_hz = cfg.f_list_hz(ifq);
        k0 = 2*pi*f_hz / cfg.c0_mps;

        R_ssa1 = cfg.R0 * exp(-2 * k0^2 * sigma_eta_raw_m^2);
        abs_R_ssa1 = abs(R_ssa1);
        power_R_ssa1 = abs_R_ssa1^2;
        loss_R_ssa1_db = local_loss_db(abs_R_ssa1);

        R_k_corrected = complex(zeros(numel(cfg.seed_list), 1));
        R_k_legacy = complex(zeros(numel(cfg.seed_list), 1));
        P_k_corrected = zeros(numel(cfg.seed_list), 1);
        P_k_legacy = zeros(numel(cfg.seed_list), 1);

        for iseed = 1:numel(cfg.seed_list)
            seed = cfg.seed_list(iseed);
            eta_xy = local_pm_surface_realization(pm_spec.Phi2D, dkx, dky, seed);
            eta_stack_std(iseed) = std(eta_xy(:));

            delta_phi_corrected = 2 * k0 * eta_xy;
            G_corrected = cfg.R0 .* exp(1i * delta_phi_corrected);
            R_k_corrected(iseed) = mean(G_corrected(:));
            P_k_corrected(iseed) = mean(abs(G_corrected(:)).^2);

            delta_phi_legacy = 2 * k0 * cfg.gamma_legacy_4k_eta_normal * eta_xy;
            G_legacy = cfg.R0 .* exp(1i * delta_phi_legacy);
            R_k_legacy(iseed) = mean(G_legacy(:));
            P_k_legacy(iseed) = mean(abs(G_legacy(:)).^2);

            run_idx = run_idx + 1;
            run_rows(run_idx).wind_speed_mps = U; %#ok<SAGROW>
            run_rows(run_idx).f_hz = f_hz;
            run_rows(run_idx).seed = seed;
            run_rows(run_idx).k0_rad_per_m = k0;
            run_rows(run_idx).sigma_eta_raw_m = sigma_eta_raw_m;
            run_rows(run_idx).Hs_raw_m = Hs_raw_m;
            run_rows(run_idx).eta_std_realization_m = eta_stack_std(iseed);
            run_rows(run_idx).R_ssa1 = R_ssa1;
            run_rows(run_idx).abs_R_ssa1 = abs_R_ssa1;
            run_rows(run_idx).power_R_ssa1 = power_R_ssa1;
            run_rows(run_idx).loss_R_ssa1_db = loss_R_ssa1_db;
            run_rows(run_idx).R_kirchhoff_corrected = R_k_corrected(iseed);
            run_rows(run_idx).abs_R_kirchhoff_corrected = abs(R_k_corrected(iseed));
            run_rows(run_idx).power_R_kirchhoff_corrected = abs(R_k_corrected(iseed))^2;
            run_rows(run_idx).kirchhoff_total_reflected_power_corrected = P_k_corrected(iseed);
            run_rows(run_idx).R_kirchhoff_legacy_4k_eta = R_k_legacy(iseed);
            run_rows(run_idx).abs_R_kirchhoff_legacy_4k_eta = abs(R_k_legacy(iseed));
            run_rows(run_idx).power_R_kirchhoff_legacy_4k_eta = abs(R_k_legacy(iseed))^2;
            run_rows(run_idx).kirchhoff_total_reflected_power_legacy_4k_eta = P_k_legacy(iseed);
        end

        abs_R_k_corrected = abs(R_k_corrected);
        abs_R_k_legacy = abs(R_k_legacy);
        R_k_corrected_gaussian_expected = cfg.R0 * exp(-2 * k0^2 * sigma_eta_raw_m^2);
        R_k_legacy_gaussian_expected = cfg.R0 * exp(-8 * k0^2 * sigma_eta_raw_m^2);

        row_idx = row_idx + 1;
        summary_rows(row_idx).wind_speed_mps = U; %#ok<SAGROW>
        summary_rows(row_idx).f_hz = f_hz;
        summary_rows(row_idx).k0_rad_per_m = k0;
        summary_rows(row_idx).sigma_eta_raw_m = sigma_eta_raw_m;
        summary_rows(row_idx).Hs_raw_m = Hs_raw_m;
        summary_rows(row_idx).pm_variance_discrete_m2 = sigma_eta_raw_m^2;
        summary_rows(row_idx).pm_variance_rel_error = variance_rel_error;
        summary_rows(row_idx).eta_std_realization_mean_m = mean(eta_stack_std);
        summary_rows(row_idx).eta_std_realization_std_m = std(eta_stack_std, 0);
        summary_rows(row_idx).eta_variance_realization_mean_m2 = mean(eta_stack_std.^2);
        summary_rows(row_idx).eta_variance_to_pm_variance_ratio = ...
            mean(eta_stack_std.^2) / max(sigma_eta_raw_m^2, eps);
        summary_rows(row_idx).R_ssa1 = R_ssa1;
        summary_rows(row_idx).abs_R_ssa1 = abs_R_ssa1;
        summary_rows(row_idx).power_R_ssa1 = power_R_ssa1;
        summary_rows(row_idx).loss_R_ssa1_db = loss_R_ssa1_db;
        summary_rows(row_idx).R_kirchhoff_corrected_mean = mean(R_k_corrected);
        summary_rows(row_idx).abs_mean_R_kirchhoff_corrected_ensemble = abs(mean(R_k_corrected));
        summary_rows(row_idx).abs_R_kirchhoff_corrected_mean = mean(abs_R_k_corrected);
        summary_rows(row_idx).abs_R_kirchhoff_corrected_std = std(abs_R_k_corrected, 0);
        summary_rows(row_idx).power_R_kirchhoff_corrected_mean = mean(abs_R_k_corrected.^2);
        summary_rows(row_idx).kirchhoff_total_reflected_power_corrected_mean = mean(P_k_corrected);
        summary_rows(row_idx).kirchhoff_total_reflected_power_corrected_std = std(P_k_corrected, 0);
        summary_rows(row_idx).R_kirchhoff_corrected_gaussian_expected = R_k_corrected_gaussian_expected;
        summary_rows(row_idx).R_kirchhoff_legacy_4k_eta_mean = mean(R_k_legacy);
        summary_rows(row_idx).abs_mean_R_kirchhoff_legacy_4k_eta_ensemble = abs(mean(R_k_legacy));
        summary_rows(row_idx).abs_R_kirchhoff_legacy_4k_eta_mean = mean(abs_R_k_legacy);
        summary_rows(row_idx).abs_R_kirchhoff_legacy_4k_eta_std = std(abs_R_k_legacy, 0);
        summary_rows(row_idx).power_R_kirchhoff_legacy_4k_eta_mean = mean(abs_R_k_legacy.^2);
        summary_rows(row_idx).kirchhoff_total_reflected_power_legacy_4k_eta_mean = mean(P_k_legacy);
        summary_rows(row_idx).kirchhoff_total_reflected_power_legacy_4k_eta_std = std(P_k_legacy, 0);
        summary_rows(row_idx).R_kirchhoff_legacy_4k_eta_gaussian_expected = R_k_legacy_gaussian_expected;
        summary_rows(row_idx).ssa_vs_kirchhoff_corrected_delta_db = ...
            20*log10(max(abs_R_ssa1, eps) / max(abs(mean(R_k_corrected)), eps));
        summary_rows(row_idx).ssa_vs_kirchhoff_legacy_4k_eta_delta_db = ...
            20*log10(max(abs_R_ssa1, eps) / max(abs(mean(R_k_legacy)), eps));
        summary_rows(row_idx).phase_convention_note = ...
            "corrected Kirchhoff uses delta_phi=2*k0*eta in normal incidence; legacy_4k_eta uses delta_phi=2*k0*2*eta as the pre-fix diagnostic.";
    end
end

run_table = struct2table(run_rows);
summary_table = struct2table(summary_rows);
summary_csv_table = local_csv_summary_table(summary_table);
controlled_gaussian_table = local_controlled_gaussian_validation(summary_table, cfg);
validation_report = local_validation_report(summary_table, flat_check, controlled_gaussian_table);

mat_file = 'compare_plane_wave_surface_coherent_reflection_vertical_result.mat';
csv_file = 'compare_plane_wave_surface_coherent_reflection_vertical_summary.csv';
controlled_csv_file = 'compare_plane_wave_surface_coherent_reflection_vertical_controlled_gaussian.csv';
save(mat_file, 'cfg', 'run_table', 'summary_table', 'summary_csv_table', ...
    'controlled_gaussian_table', 'validation_report', 'flat_check');
writetable(summary_csv_table, csv_file);
writetable(controlled_gaussian_table, controlled_csv_file);

local_plot_absR_vs_wind(summary_table, cfg);
local_plot_power_vs_wind(summary_table, cfg);
local_plot_loss_vs_wind(summary_table, cfg);
local_plot_heatmap(summary_table, cfg);
local_plot_sigma_vs_wind(summary_table, cfg);

disp(summary_table(:, {'wind_speed_mps', 'f_hz', 'sigma_eta_raw_m', 'Hs_raw_m', ...
    'abs_R_ssa1', 'abs_R_kirchhoff_corrected_mean', ...
    'abs_mean_R_kirchhoff_corrected_ensemble', ...
    'abs_R_kirchhoff_legacy_4k_eta_mean', 'kirchhoff_total_reflected_power_corrected_mean'}));
disp(controlled_gaussian_table(:, {'sigma_eta_m', 'f_hz', 'expected_abs_R', ...
    'abs_mean_R_ensemble', 'abs_error', 'mean_abs_R_sample'}));
disp(validation_report);

fprintf('Saved %s, %s, and %s\n', mat_file, csv_file, controlled_csv_file);

function [KX, KY, dkx, dky] = local_wavenumber_grid(nx, ny, xw, yw)
dkx = 2*pi / xw;
dky = 2*pi / yw;
kx = dkx * [0:(nx/2-1), -nx/2:-1];
ky = dky * [0:(ny/2-1), -ny/2:-1];
[KX, KY] = meshgrid(kx, ky);
end

function spec = local_pm_spectrum(KX, KY, xw, yw, U)
if ~(isscalar(U) && isnumeric(U) && isfinite(U) && U > 0)
    error('Wind speed U must be positive and finite.');
end
g = 9.81;
alpha_PM = 8.10e-3;
beta_PM = 0.74;
K = sqrt(KX.^2 + KY.^2);
E1D = zeros(size(K));
Phi2D = zeros(size(K));
mask = K > 0;
K_nonzero = K(mask);
E1D(mask) = (alpha_PM ./ (2 .* K_nonzero.^3)) .* ...
    exp(-beta_PM * (g^2) ./ (U^4 .* K_nonzero.^2));
Phi2D(mask) = E1D(mask) ./ (2*pi*K_nonzero);
spec = struct('U', U, 'g', g, 'alpha_PM', alpha_PM, 'beta_PM', beta_PM, ...
    'K', K, 'E1D', E1D, 'Phi2D', Phi2D, 'dkx', 2*pi/xw, 'dky', 2*pi/yw);
end

function eta_xy = local_pm_surface_realization(Phi2D, dkx, dky, seed)
A = sqrt(Phi2D .* dkx .* dky);
rng(mod(round(seed), 2^32), 'twister');
N = (randn(size(Phi2D)) + 1i*randn(size(Phi2D))) / sqrt(2);
Zk = A .* N;
eta_xy = real(ifft2(Zk));
eta_xy = eta_xy * numel(eta_xy);
end

function flat_check = local_flat_check(cfg)
k0 = 2*pi*cfg.f_list_hz(1) / cfg.c0_mps;
sigma0 = 0;
eta0 = zeros(cfg.ny, cfg.nx);
R_ssa = cfg.R0 * exp(-2*k0^2*sigma0^2);
R_k = mean(cfg.R0 .* exp(1i * 2*k0*eta0), 'all');
flat_check = struct( ...
    'abs_R_ssa1_flat', abs(R_ssa), ...
    'abs_R_kirchhoff_flat', abs(R_k), ...
    'passed', abs(abs(R_ssa) - 1) <= 1e-14 && abs(abs(R_k) - 1) <= 1e-14);
end

function T = local_controlled_gaussian_validation(summary_table, cfg)
rows = struct([]);
row_idx = 0;
for ir = 1:height(summary_table)
    sigma_eta = summary_table.sigma_eta_raw_m(ir);
    f_hz = summary_table.f_hz(ir);
    k0 = summary_table.k0_rad_per_m(ir);
    expected_abs_R = exp(-2 * k0^2 * sigma_eta^2);
    R_samples = complex(zeros(numel(cfg.seed_list), 1));
    for iseed = 1:numel(cfg.seed_list)
        rng(mod(round(cfg.seed_list(iseed) + 200000 + ir), 2^32), 'twister');
        eta_xy = sigma_eta .* randn(cfg.ny, cfg.nx);
        G_xy = cfg.R0 .* exp(1i * 2 * k0 .* eta_xy);
        R_samples(iseed) = mean(G_xy(:));
    end
    row_idx = row_idx + 1;
    rows(row_idx).sigma_eta_m = sigma_eta; %#ok<AGROW>
    rows(row_idx).f_hz = f_hz;
    rows(row_idx).k0_rad_per_m = k0;
    rows(row_idx).expected_abs_R = expected_abs_R;
    rows(row_idx).abs_mean_R_ensemble = abs(mean(R_samples));
    rows(row_idx).mean_abs_R_sample = mean(abs(R_samples));
    rows(row_idx).std_abs_R_sample = std(abs(R_samples), 0);
    rows(row_idx).abs_error = abs(rows(row_idx).abs_mean_R_ensemble - expected_abs_R);
    rows(row_idx).sample_count = numel(R_samples) * cfg.nx * cfg.ny;
    rows(row_idx).interpretation = ...
        "controlled iid Gaussian eta validates coherent ensemble; mean_abs_R_sample is finite-aperture residual";
end
T = struct2table(rows);
end

function report = local_validation_report(summary_table, flat_check, controlled_gaussian_table)
power_err = max(abs(summary_table.kirchhoff_total_reflected_power_corrected_mean - 1));
pm_var_err = max(abs(summary_table.pm_variance_rel_error));
finite_values = all(isfinite(summary_table.abs_R_ssa1)) && ...
    all(isfinite(summary_table.abs_R_kirchhoff_corrected_mean)) && ...
    all(isfinite(summary_table.abs_R_kirchhoff_legacy_4k_eta_mean)) && ...
    all(isfinite(controlled_gaussian_table.abs_mean_R_ensemble));
if height(controlled_gaussian_table) == 0
    gaussian_tolerance = Inf;
else
    sample_count_min = min(controlled_gaussian_table.sample_count);
    gaussian_tolerance = max(0.03, 8 / sqrt(max(sample_count_min, 1)));
end
gaussian_max_abs_error = max(controlled_gaussian_table.abs_error);
report = struct();
report.flat_check_passed = flat_check.passed;
report.max_pm_variance_rel_error = pm_var_err;
report.max_kirchhoff_total_power_error = power_err;
report.controlled_gaussian_max_abs_error = gaussian_max_abs_error;
report.controlled_gaussian_tolerance = gaussian_tolerance;
report.controlled_gaussian_passed = gaussian_max_abs_error <= gaussian_tolerance;
report.all_finite = finite_values;
report.all_passed = flat_check.passed && pm_var_err <= 1e-14 && power_err <= 1e-12 && ...
    finite_values && report.controlled_gaussian_passed;
report.note = ['This is a surface-boundary-only coherent reflection diagnostic. ', ...
    'It does not run PE/WAPE and does not test H_f channel invariants. ', ...
    'abs(mean(R_sample)) is the coherent ensemble estimate; mean(abs(R_sample)) is a finite-aperture residual.'];
end

function Tcsv = local_csv_summary_table(T)
keep = {'wind_speed_mps', 'f_hz', 'k0_rad_per_m', 'sigma_eta_raw_m', 'Hs_raw_m', ...
    'pm_variance_discrete_m2', 'pm_variance_rel_error', ...
    'eta_std_realization_mean_m', 'eta_std_realization_std_m', ...
    'eta_variance_realization_mean_m2', 'eta_variance_to_pm_variance_ratio', ...
    'abs_R_ssa1', 'power_R_ssa1', 'loss_R_ssa1_db', ...
    'abs_mean_R_kirchhoff_corrected_ensemble', ...
    'abs_R_kirchhoff_corrected_mean', 'abs_R_kirchhoff_corrected_std', ...
    'power_R_kirchhoff_corrected_mean', ...
    'kirchhoff_total_reflected_power_corrected_mean', ...
    'kirchhoff_total_reflected_power_corrected_std', ...
    'abs_mean_R_kirchhoff_legacy_4k_eta_ensemble', ...
    'abs_R_kirchhoff_legacy_4k_eta_mean', 'abs_R_kirchhoff_legacy_4k_eta_std', ...
    'power_R_kirchhoff_legacy_4k_eta_mean', ...
    'kirchhoff_total_reflected_power_legacy_4k_eta_mean', ...
    'kirchhoff_total_reflected_power_legacy_4k_eta_std', ...
    'ssa_vs_kirchhoff_corrected_delta_db', ...
    'ssa_vs_kirchhoff_legacy_4k_eta_delta_db'};
Tcsv = T(:, keep);
end

function local_plot_absR_vs_wind(T, cfg)
fig = figure('Visible', 'off', 'Color', 'w');
tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
for ifq = 1:numel(cfg.f_list_hz)
    nexttile;
    f_hz = cfg.f_list_hz(ifq);
    S = T(T.f_hz == f_hz, :);
    semilogy(S.wind_speed_mps, S.abs_R_ssa1, '-o', 'LineWidth', 1.2); hold on;
    errorbar(S.wind_speed_mps, S.abs_R_kirchhoff_corrected_mean, ...
        S.abs_R_kirchhoff_corrected_std, '-s', 'LineWidth', 1.2);
    errorbar(S.wind_speed_mps, S.abs_R_kirchhoff_legacy_4k_eta_mean, ...
        S.abs_R_kirchhoff_legacy_4k_eta_std, '-^', 'LineWidth', 1.2);
    grid on; xlabel('Wind speed U (m/s)'); ylabel('|R_{coh}|');
    title(sprintf('f = %.1f kHz', f_hz/1000));
    legend('SSA1', 'Kirchhoff corrected', 'Kirchhoff legacy 4k\eta', 'Location', 'southwest');
end
sgtitle('Vertical plane wave coherent reflection amplitude');
exportgraphics(fig, 'plane_wave_coherent_absR_vs_wind.png', 'Resolution', 180);
close(fig);
end

function local_plot_power_vs_wind(T, cfg)
fig = figure('Visible', 'off', 'Color', 'w');
tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
for ifq = 1:numel(cfg.f_list_hz)
    nexttile;
    f_hz = cfg.f_list_hz(ifq);
    S = T(T.f_hz == f_hz, :);
    semilogy(S.wind_speed_mps, S.power_R_ssa1, '-o', 'LineWidth', 1.2); hold on;
    semilogy(S.wind_speed_mps, S.power_R_kirchhoff_corrected_mean, '-s', 'LineWidth', 1.2);
    semilogy(S.wind_speed_mps, S.power_R_kirchhoff_legacy_4k_eta_mean, '-^', 'LineWidth', 1.2);
    grid on; xlabel('Wind speed U (m/s)'); ylabel('|R_{coh}|^2');
    title(sprintf('f = %.1f kHz', f_hz/1000));
    legend('SSA1', 'Kirchhoff corrected', 'Kirchhoff legacy 4k\eta', 'Location', 'southwest');
end
sgtitle('Vertical plane wave coherent reflection power');
exportgraphics(fig, 'plane_wave_coherent_power_vs_wind.png', 'Resolution', 180);
close(fig);
end

function local_plot_loss_vs_wind(T, cfg)
fig = figure('Visible', 'off', 'Color', 'w');
tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
for ifq = 1:numel(cfg.f_list_hz)
    nexttile;
    f_hz = cfg.f_list_hz(ifq);
    S = T(T.f_hz == f_hz, :);
    plot(S.wind_speed_mps, S.loss_R_ssa1_db, '-o', 'LineWidth', 1.2); hold on;
    plot(S.wind_speed_mps, local_loss_db(S.abs_R_kirchhoff_corrected_mean), '-s', 'LineWidth', 1.2);
    plot(S.wind_speed_mps, local_loss_db(S.abs_R_kirchhoff_legacy_4k_eta_mean), '-^', 'LineWidth', 1.2);
    grid on; xlabel('Wind speed U (m/s)'); ylabel('Coherent loss (dB)');
    title(sprintf('f = %.1f kHz', f_hz/1000));
    legend('SSA1', 'Kirchhoff corrected', 'Kirchhoff legacy 4k\eta', 'Location', 'northwest');
end
sgtitle('Vertical plane wave coherent reflection loss');
exportgraphics(fig, 'plane_wave_coherent_loss_vs_wind.png', 'Resolution', 180);
close(fig);
end

function local_plot_heatmap(T, cfg)
[W, F] = ndgrid(cfg.wind_list_mps, cfg.f_list_hz);
Z = nan(size(W));
for idx = 1:numel(W)
    row = T(T.wind_speed_mps == W(idx) & T.f_hz == F(idx), :);
    if ~isempty(row)
        Z(idx) = row.ssa_vs_kirchhoff_legacy_4k_eta_delta_db(1);
    end
end
fig = figure('Visible', 'off', 'Color', 'w');
imagesc(cfg.f_list_hz/1000, cfg.wind_list_mps, Z);
set(gca, 'YDir', 'normal');
colorbar; grid on;
xlabel('Frequency (kHz)'); ylabel('Wind speed U (m/s)');
title('SSA1 minus Kirchhoff legacy 4k\eta coherent amplitude (dB)');
exportgraphics(fig, 'plane_wave_kirchhoff_vs_ssa_absR_heatmap.png', 'Resolution', 180);
close(fig);
end

function local_plot_sigma_vs_wind(T, cfg)
fig = figure('Visible', 'off', 'Color', 'w');
[winds_unique, ia] = unique(T.wind_speed_mps, 'stable');
sigma_unique = T.sigma_eta_raw_m(ia);
Hs_unique = T.Hs_raw_m(ia);
yyaxis left;
plot(winds_unique, sigma_unique, '-o', 'LineWidth', 1.4);
ylabel('\sigma_\eta from raw PM spectrum (m)');
yyaxis right;
plot(winds_unique, Hs_unique, '-s', 'LineWidth', 1.4);
ylabel('H_{s,raw}=4\sigma_\eta (m)');
grid on; xlabel('Wind speed U (m/s)');
title(sprintf('Raw PM roughness without H_s normalization, %d x %d grid', cfg.nx, cfg.ny));
exportgraphics(fig, 'plane_wave_raw_pm_sigma_vs_wind.png', 'Resolution', 180);
close(fig);
end

function values = local_env_numeric_list(name, default_values)
raw = getenv(name);
if isempty(raw)
    values = default_values;
    return
end
raw = strrep(raw, ',', ' ');
values = sscanf(raw, '%f').';
if isempty(values)
    values = default_values;
end
end

function value = local_env_scalar(name, default_value)
raw = getenv(name);
if isempty(raw)
    value = default_value;
    return
end
value = str2double(raw);
if ~isfinite(value)
    value = default_value;
end
end

function value = local_env_string(name, default_value)
raw = getenv(name);
if isempty(raw)
    value = default_value;
else
    value = raw;
end
end

function e = local_relative_error(value, reference)
if abs(reference) <= eps
    e = abs(value - reference);
else
    e = abs(value - reference) / abs(reference);
end
end

function loss_db = local_loss_db(abs_R)
loss_db = -20 * log10(max(abs_R, eps));
end

