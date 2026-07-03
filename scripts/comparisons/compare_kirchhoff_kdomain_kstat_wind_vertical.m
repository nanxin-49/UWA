run(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'bootstrap_project.m'));
%COMPARE_KIRCHHOFF_KDOMAIN_KSTAT_WIND_VERTICAL
% Wind-driven raw-PM comparison between explicit Kirchhoff k-domain and
% Kirchhoff K-Stat phase-screen branches.

clear
clc
format compact

cfg = struct();
cfg.wind_list_mps = local_env_numeric_vector('KIRCH_WIND_LIST', [3, 5, 8, 10, 12, 15]);
cfg.seed_count = local_env_int('KIRCH_SEED_COUNT', 32);
cfg.seed_list = local_env_numeric_vector('KIRCH_SEED_LIST', 12345 + (0:(cfg.seed_count - 1)));
cfg.seed_list = cfg.seed_list(1:min(numel(cfg.seed_list), cfg.seed_count));
cfg.f0_hz = local_env_scalar('KIRCH_F0_HZ', 6000);
cfg.grid_n = local_env_int('KIRCH_GRID_N', 128);
cfg.nx = local_env_int('KIRCH_NX', cfg.grid_n);
cfg.ny = local_env_int('KIRCH_NY', cfg.grid_n);
cfg.xw_m = local_env_scalar('KIRCH_XW_M', 50 * cfg.nx / 128);
cfg.yw_m = local_env_scalar('KIRCH_YW_M', 50 * cfg.ny / 128);
cfg.save_mode = 'rx_only';
cfg.show_figures = false;
cfg.roughness_scale_mode = 'raw_pm';
cfg.result_mat = project_result_file('comparisons', ...
    'compare_kirchhoff_kdomain_kstat_wind_vertical_result.mat');
cfg.summary_csv = project_result_file('comparisons', ...
    'compare_kirchhoff_kdomain_kstat_wind_vertical_summary.csv');
cfg.figure_dir = project_result_dir('comparisons');
cfg.energy_tol = 1e-8;
cfg.invariant_tol = 1e-10;

fprintf('Raw-PM Kirchhoff comparison\n');
fprintf('  winds: %s m/s\n', mat2str(cfg.wind_list_mps));
fprintf('  seeds: %d, grid: %d x %d, f0: %.3f Hz\n', ...
    numel(cfg.seed_list), cfg.nx, cfg.ny, cfg.f0_hz);
fprintf('  aperture: %.3g m x %.3g m\n', cfg.xw_m, cfg.yw_m);

summary_rows = struct([]);
run_rows = struct([]);
row_idx = 0;
run_idx = 0;

for iw = 1:numel(cfg.wind_list_mps)
    U = cfg.wind_list_mps(iw);
    fprintf('Wind %.3g m/s: running coherent K-Stat reference.\n', U);

    kstat_ref_params = local_base_params(cfg, U, cfg.seed_list(1), 'kirchhoff_kstat');
    kstat_ref_params.surface_kstat_random_scatter = false;
    kstat_ref = vertical_channel_model(kstat_ref_params);
    kstat_ref_meta = kstat_ref.roughness_meta.kirchhoff_kstat_meta;

    kd_abs_h_reflect = zeros(numel(cfg.seed_list), 1);
    kd_abs_h_total = zeros(numel(cfg.seed_list), 1);
    kd_screen_R = complex(zeros(numel(cfg.seed_list), 1));
    kd_Hs_raw = zeros(numel(cfg.seed_list), 1);
    kd_sigma_raw = zeros(numel(cfg.seed_list), 1);
    kd_invariant = zeros(numel(cfg.seed_list), 1);

    ks_abs_h_reflect = zeros(numel(cfg.seed_list), 1);
    ks_abs_h_total = zeros(numel(cfg.seed_list), 1);
    ks_invariant = zeros(numel(cfg.seed_list), 1);
    ks_energy_error = zeros(numel(cfg.seed_list), 1);

    for iseed = 1:numel(cfg.seed_list)
        seed = cfg.seed_list(iseed);
        fprintf('  seed %d/%d: %d\n', iseed, numel(cfg.seed_list), seed);

        kd_params = local_base_params(cfg, U, seed, 'kirchhoff_kdomain');
        kd_channel = vertical_channel_model(kd_params);
        kd_meta = kd_channel.roughness_meta;
        kd_abs_h_reflect(iseed) = abs(kd_channel.h_reflect);
        kd_abs_h_total(iseed) = abs(kd_channel.h_total);
        kd_invariant(iseed) = local_invariant_error(kd_channel);
        kd_Hs_raw(iseed) = kd_meta.Hs_raw_m;
        kd_sigma_raw(iseed) = kd_meta.sigma_eta_raw_m;
        kd_screen_R(iseed) = local_screen_coherent_R(kd_channel.surface_elevation, ...
            kd_params.surface_reflect_coeff, 2*pi*cfg.f0_hz/kd_params.c0);

        ks_params = local_base_params(cfg, U, seed, 'kirchhoff_kstat');
        ks_params.surface_kstat_random_scatter = true;
        ks_channel = vertical_channel_model(ks_params);
        ks_meta = ks_channel.roughness_meta.kirchhoff_kstat_meta;
        ks_abs_h_reflect(iseed) = abs(ks_channel.h_reflect);
        ks_abs_h_total(iseed) = abs(ks_channel.h_total);
        ks_invariant(iseed) = local_invariant_error(ks_channel);
        ks_energy_error(iseed) = ks_meta.phase_screen_energy_error;

        run_idx = run_idx + 1;
        run_rows(run_idx).wind_speed_mps = U; %#ok<SAGROW>
        run_rows(run_idx).seed = seed;
        run_rows(run_idx).kdomain_abs_h_reflect = kd_abs_h_reflect(iseed);
        run_rows(run_idx).kdomain_abs_h_total = kd_abs_h_total(iseed);
        run_rows(run_idx).kdomain_screen_R_real = real(kd_screen_R(iseed));
        run_rows(run_idx).kdomain_screen_R_imag = imag(kd_screen_R(iseed));
        run_rows(run_idx).kdomain_screen_abs_R = abs(kd_screen_R(iseed));
        run_rows(run_idx).kdomain_sigma_eta_raw_m = kd_sigma_raw(iseed);
        run_rows(run_idx).kdomain_Hs_raw_m = kd_Hs_raw(iseed);
        run_rows(run_idx).kdomain_invariant_error = kd_invariant(iseed);
        run_rows(run_idx).kstat_abs_h_reflect = ks_abs_h_reflect(iseed);
        run_rows(run_idx).kstat_abs_h_total = ks_abs_h_total(iseed);
        run_rows(run_idx).kstat_invariant_error = ks_invariant(iseed);
        run_rows(run_idx).kstat_phase_screen_energy_error = ks_energy_error(iseed);
    end

    kd_R_mean = mean(kd_screen_R);
    row_idx = row_idx + 1;
    summary_rows(row_idx).wind_speed_mps = U; %#ok<SAGROW>
    summary_rows(row_idx).f0_hz = cfg.f0_hz;
    summary_rows(row_idx).seed_count = numel(cfg.seed_list);
    summary_rows(row_idx).nx = cfg.nx;
    summary_rows(row_idx).ny = cfg.ny;
    summary_rows(row_idx).roughness_scale_mode_code = 2; % 2 means raw_pm.
    summary_rows(row_idx).kdomain_sigma_eta_raw_mean_m = mean(kd_sigma_raw);
    summary_rows(row_idx).kdomain_sigma_eta_raw_std_m = std(kd_sigma_raw, 0);
    summary_rows(row_idx).kdomain_Hs_raw_mean_m = mean(kd_Hs_raw);
    summary_rows(row_idx).kdomain_Hs_raw_std_m = std(kd_Hs_raw, 0);
    summary_rows(row_idx).kstat_sigma_eta_raw_m = kstat_ref_meta.sigma_eta_raw_m;
    summary_rows(row_idx).kstat_Hs_raw_m = kstat_ref_meta.Hs_raw_m;
    summary_rows(row_idx).Hs_raw_rel_error = ...
        abs(summary_rows(row_idx).kdomain_Hs_raw_mean_m - summary_rows(row_idx).kstat_Hs_raw_m) / ...
        max(summary_rows(row_idx).kstat_Hs_raw_m, eps);
    summary_rows(row_idx).pm_variance_raw_discrete_m2 = kstat_ref.roughness_meta.pm_variance_raw_discrete_m2;
    summary_rows(row_idx).pm_variance_raw_continuous_m2 = kstat_ref.roughness_meta.pm_variance_raw_continuous_m2;
    summary_rows(row_idx).kdomain_screen_R_mean_real = real(kd_R_mean);
    summary_rows(row_idx).kdomain_screen_R_mean_imag = imag(kd_R_mean);
    summary_rows(row_idx).kdomain_screen_abs_mean_R = abs(kd_R_mean);
    summary_rows(row_idx).kdomain_screen_abs_R_mean = mean(abs(kd_screen_R));
    summary_rows(row_idx).kdomain_screen_abs_R_std = std(abs(kd_screen_R), 0);
    summary_rows(row_idx).kstat_R_coh_real = real(kstat_ref_meta.R_coh);
    summary_rows(row_idx).kstat_R_coh_imag = imag(kstat_ref_meta.R_coh);
    summary_rows(row_idx).kstat_abs_R_coh = abs(kstat_ref_meta.R_coh);
    summary_rows(row_idx).coherent_abs_delta = abs(abs(kd_R_mean) - abs(kstat_ref_meta.R_coh));
    summary_rows(row_idx).P_sca_sum = kstat_ref_meta.P_sca.sum;
    summary_rows(row_idx).phase_screen_incoherent_energy = kstat_ref_meta.phase_screen_incoherent_energy;
    summary_rows(row_idx).phase_screen_incoherent_energy_propagating = ...
        kstat_ref_meta.phase_screen_incoherent_energy_propagating;
    summary_rows(row_idx).phase_screen_energy_error = kstat_ref_meta.phase_screen_energy_error;
    summary_rows(row_idx).kdomain_abs_h_reflect_mean = mean(kd_abs_h_reflect);
    summary_rows(row_idx).kdomain_abs_h_reflect_std = std(kd_abs_h_reflect, 0);
    summary_rows(row_idx).kdomain_abs_h_reflect_p10 = local_percentile(kd_abs_h_reflect, 10);
    summary_rows(row_idx).kdomain_abs_h_reflect_p50 = local_percentile(kd_abs_h_reflect, 50);
    summary_rows(row_idx).kdomain_abs_h_reflect_p90 = local_percentile(kd_abs_h_reflect, 90);
    summary_rows(row_idx).kdomain_abs_h_total_mean = mean(kd_abs_h_total);
    summary_rows(row_idx).kdomain_abs_h_total_std = std(kd_abs_h_total, 0);
    summary_rows(row_idx).kdomain_abs_h_total_p10 = local_percentile(kd_abs_h_total, 10);
    summary_rows(row_idx).kdomain_abs_h_total_p50 = local_percentile(kd_abs_h_total, 50);
    summary_rows(row_idx).kdomain_abs_h_total_p90 = local_percentile(kd_abs_h_total, 90);
    summary_rows(row_idx).kstat_abs_h_reflect_mean = mean(ks_abs_h_reflect);
    summary_rows(row_idx).kstat_abs_h_reflect_std = std(ks_abs_h_reflect, 0);
    summary_rows(row_idx).kstat_abs_h_reflect_p10 = local_percentile(ks_abs_h_reflect, 10);
    summary_rows(row_idx).kstat_abs_h_reflect_p50 = local_percentile(ks_abs_h_reflect, 50);
    summary_rows(row_idx).kstat_abs_h_reflect_p90 = local_percentile(ks_abs_h_reflect, 90);
    summary_rows(row_idx).kstat_abs_h_total_mean = mean(ks_abs_h_total);
    summary_rows(row_idx).kstat_abs_h_total_std = std(ks_abs_h_total, 0);
    summary_rows(row_idx).kstat_abs_h_total_p10 = local_percentile(ks_abs_h_total, 10);
    summary_rows(row_idx).kstat_abs_h_total_p50 = local_percentile(ks_abs_h_total, 50);
    summary_rows(row_idx).kstat_abs_h_total_p90 = local_percentile(ks_abs_h_total, 90);
    summary_rows(row_idx).max_kdomain_invariant_error = max(kd_invariant);
    summary_rows(row_idx).max_kstat_invariant_error = max(ks_invariant);
    summary_rows(row_idx).max_kstat_seed_energy_error = max(ks_energy_error);
    summary_rows(row_idx).kstat_ref_invariant_error = local_invariant_error(kstat_ref);
end

run_table = struct2table(run_rows);
summary_table = struct2table(summary_rows);
validation_report = local_validation_report(summary_table, cfg);

save(cfg.result_mat, 'cfg', 'run_table', 'summary_table', 'validation_report');
writetable(summary_table, cfg.summary_csv);

local_plot_Hs_vs_wind(summary_table, cfg);
local_plot_R_vs_wind(summary_table, cfg);
local_plot_abs_h_reflect_vs_wind(summary_table, cfg);
local_plot_incoherent_energy_vs_wind(summary_table, cfg);
local_plot_propagating_energy_vs_wind(summary_table, cfg);

disp(summary_table(:, {'wind_speed_mps', 'kdomain_Hs_raw_mean_m', 'kstat_Hs_raw_m', ...
    'Hs_raw_rel_error', 'kdomain_screen_abs_mean_R', 'kstat_abs_R_coh', ...
    'kdomain_abs_h_reflect_mean', 'kstat_abs_h_reflect_mean', ...
    'phase_screen_incoherent_energy', 'phase_screen_energy_error'}));
disp(validation_report);
fprintf('Saved %s\n', cfg.result_mat);
fprintf('Saved %s\n', cfg.summary_csv);

function paramsV = local_base_params(cfg, wind_speed_mps, seed, boundary_model)
paramsV = struct();
paramsV.f0 = cfg.f0_hz;
paramsV.enable_wideband = false;
paramsV.c0 = 1500;
paramsV.z_max = 100;
paramsV.stepz_lamb = 0.5;
paramsV.xw = cfg.xw_m;
paramsV.yw = cfg.yw_m;
paramsV.nx = cfg.nx;
paramsV.ny = cfg.ny;
paramsV.x_tx = 0;
paramsV.y_tx = 0;
paramsV.z_tx = 100;
paramsV.x_rx = 0;
paramsV.y_rx = 0;
paramsV.z_rx = 3;
paramsV.nout = 2;
paramsV.sigma_src_m = 0.4;
paramsV.sponge_ratio = 0.12;
paramsV.alpha_max_np_per_m = 0.15;
paramsV.env_mode = 'uniform';
paramsV.enforce_1_over_R = false;
paramsV.show_figures = cfg.show_figures;
paramsV.save_mode = cfg.save_mode;
paramsV.use_gpu = false;
paramsV.enable_surface_reflection = true;
paramsV.sea_wind_speed = wind_speed_mps;
paramsV.sea_hs_target = 0.5;
paramsV.surface_roughness_scale_mode = cfg.roughness_scale_mode;
paramsV.sea_seed = seed;
paramsV.surface_reflect_coeff = -1;
paramsV.surface_phase_mode = 'normal';
paramsV.surface_boundary_model = boundary_model;
paramsV.surface_boundary_coupling_diagnostics = false;
paramsV.surface_boundary_redistribution_diagnostics = false;
paramsV.surface_kstat_random_scatter = true;
paramsV.surface_kstat_seed_offset = 200000;
paramsV.surface_kstat_conv_padding = 'periodic';
paramsV.surface_kstat_trusted_angle_deg = NaN;
end

function err = local_invariant_error(channel)
err = max(abs(channel.H_f(:) - channel.H_direct_f(:) - channel.H_reflect_f(:)));
end

function R = local_screen_coherent_R(surface_elevation, R0, k0)
if isempty(surface_elevation)
    R = complex(NaN, NaN);
else
    R = mean(R0 .* exp(1i * 2 * k0 .* surface_elevation(:)));
end
end

function value = local_percentile(values, p)
values = sort(values(isfinite(values(:))));
if isempty(values)
    value = NaN;
    return
end
if numel(values) == 1
    value = values(1);
    return
end
idx = 1 + (numel(values) - 1) * p / 100;
lo = floor(idx);
hi = ceil(idx);
if lo == hi
    value = values(lo);
else
    value = values(lo) + (idx - lo) * (values(hi) - values(lo));
end
end

function report = local_validation_report(summary_table, cfg)
report = struct();
report.max_kdomain_invariant_error = max(summary_table.max_kdomain_invariant_error);
report.max_kstat_invariant_error = max(summary_table.max_kstat_invariant_error);
report.max_kstat_energy_error = max(summary_table.phase_screen_energy_error);
report.kdomain_invariant_passed = report.max_kdomain_invariant_error <= cfg.invariant_tol;
report.kstat_invariant_passed = report.max_kstat_invariant_error <= cfg.invariant_tol;
report.kstat_energy_passed = report.max_kstat_energy_error <= cfg.energy_tol;
report.Hs_raw_varies_with_wind = ...
    (max(summary_table.kdomain_Hs_raw_mean_m) - min(summary_table.kdomain_Hs_raw_mean_m)) > 0 && ...
    (max(summary_table.kstat_Hs_raw_m) - min(summary_table.kstat_Hs_raw_m)) > 0;
report.output_directory = cfg.figure_dir;
end

function local_plot_Hs_vs_wind(T, cfg)
figure('Visible', 'off');
plot(T.wind_speed_mps, T.kdomain_Hs_raw_mean_m, '-o', 'LineWidth', 1.5);
hold on
plot(T.wind_speed_mps, T.kstat_Hs_raw_m, '-s', 'LineWidth', 1.5);
grid on
xlabel('Wind speed U (m/s)');
ylabel('Raw H_s (m)');
legend('kirchhoff\_kdomain realization mean', 'kirchhoff\_kstat discrete PM', 'Location', 'northwest');
title('Raw PM H_s vs wind');
exportgraphics(gcf, fullfile(cfg.figure_dir, 'kirchhoff_raw_pm_Hs_vs_wind.png'), 'Resolution', 160);
close(gcf);
end

function local_plot_R_vs_wind(T, cfg)
figure('Visible', 'off');
plot(T.wind_speed_mps, T.kdomain_screen_abs_mean_R, '-o', 'LineWidth', 1.5);
hold on
plot(T.wind_speed_mps, T.kstat_abs_R_coh, '-s', 'LineWidth', 1.5);
grid on
xlabel('Wind speed U (m/s)');
ylabel('|coherent reflection|');
legend('explicit kdomain ensemble', 'kstat R_{coh}', 'Location', 'best');
title('Coherent reflection vs wind');
exportgraphics(gcf, fullfile(cfg.figure_dir, 'kirchhoff_coherent_R_vs_wind.png'), 'Resolution', 160);
close(gcf);
end

function local_plot_abs_h_reflect_vs_wind(T, cfg)
figure('Visible', 'off');
errorbar(T.wind_speed_mps, T.kdomain_abs_h_reflect_mean, T.kdomain_abs_h_reflect_std, ...
    '-o', 'LineWidth', 1.5);
hold on
errorbar(T.wind_speed_mps, T.kstat_abs_h_reflect_mean, T.kstat_abs_h_reflect_std, ...
    '-s', 'LineWidth', 1.5);
grid on
xlabel('Wind speed U (m/s)');
ylabel('mean |h_{reflect}|');
legend('kirchhoff\_kdomain', 'kirchhoff\_kstat random', 'Location', 'best');
title('Reflected channel tap magnitude vs wind');
exportgraphics(gcf, fullfile(cfg.figure_dir, 'kirchhoff_abs_h_reflect_vs_wind.png'), 'Resolution', 160);
close(gcf);
end

function local_plot_incoherent_energy_vs_wind(T, cfg)
figure('Visible', 'off');
semilogy(T.wind_speed_mps, max(T.phase_screen_incoherent_energy, realmin), ...
    '-o', 'LineWidth', 1.5);
grid on
xlabel('Wind speed U (m/s)');
ylabel('phase-screen incoherent energy');
title('K-Stat incoherent phase-screen energy');
exportgraphics(gcf, fullfile(cfg.figure_dir, 'kirchhoff_incoherent_energy_vs_wind.png'), 'Resolution', 160);
close(gcf);
end

function local_plot_propagating_energy_vs_wind(T, cfg)
figure('Visible', 'off');
semilogy(T.wind_speed_mps, max(T.phase_screen_incoherent_energy_propagating, realmin), ...
    '-o', 'LineWidth', 1.5);
grid on
xlabel('Wind speed U (m/s)');
ylabel('propagating-window incoherent energy');
title('K-Stat propagating-window energy');
exportgraphics(gcf, fullfile(cfg.figure_dir, 'kirchhoff_propagating_energy_vs_wind.png'), 'Resolution', 160);
close(gcf);
end

function values = local_env_numeric_vector(name, default_value)
raw = getenv(name);
if isempty(raw)
    values = default_value;
    return
end
values = str2num(raw); %#ok<ST2NM>
if isempty(values)
    error('Environment variable %s must be a numeric vector.', name);
end
values = values(:).';
end

function value = local_env_scalar(name, default_value)
values = local_env_numeric_vector(name, default_value);
value = values(1);
end

function value = local_env_int(name, default_value)
value = round(local_env_scalar(name, default_value));
if value < 1
    error('%s must be a positive integer.', name);
end
end
