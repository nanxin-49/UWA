run(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'bootstrap_project.m'));
% Generate compact SSA1 group-meeting figures and traceable summary tables.
% This script does not modify the propagation, boundary, or communication code.

clear
format compact

report_dir = project_result_dir('reports', 'ssa1_group_meeting_2026-06-24');
if ~exist(report_dir, 'dir')
    mkdir(report_dir);
end

required_files = { ...
    project_result_file('comparisons', 'compare_ssa1_kirchhoff_surface_statistics_vertical_result.mat'), ...
    project_result_file('validation', 'validate_ssa_physical_trends_vertical_result.mat'), ...
    project_result_file('validation', 'validate_ssa_stat_kernel_vertical_result.mat')};
for ii = 1:numel(required_files)
    if ~exist(required_files{ii}, 'file')
        error('generate_ssa1_group_meeting_report_vertical:MissingInput', ...
            'Required result file is missing: %s', required_files{ii});
    end
end

compare_data = load(required_files{1});
trend_data = load(required_files{2});
validation_data = load(required_files{3});

local_copy_existing_figure(project_result_file('validation', 'validate_ssa_abs_R_coh_vs_Hs.png'), ...
    fullfile(report_dir, '01_Rcoh_vs_Hs.png'));
local_copy_existing_figure(project_result_file('validation', 'validate_ssa_abs_R_coh_vs_frequency.png'), ...
    fullfile(report_dir, '02_Rcoh_vs_frequency.png'));
local_copy_existing_figure(project_result_file('validation', 'validate_ssa_W_eta_variance_check.png'), ...
    fullfile(report_dir, '03_PM_variance_normalization.png'));

local_plot_energy_budget(trend_data.summary_table, ...
    fullfile(report_dir, '04_coherent_and_scatter_budget.png'));

coherent_channel_proxy = local_build_coherent_proxy(compare_data.run_table);
writetable(coherent_channel_proxy, fullfile(report_dir, 'coherent_channel_proxy.csv'));
local_plot_coherent_proxy(coherent_channel_proxy, ...
    fullfile(report_dir, '05_coherent_channel_proxy_compare.png'));

broadening_table = local_build_broadening_table(compare_data.summary_table);
local_plot_broadening(broadening_table, ...
    fullfile(report_dir, '06_rms_k_broadening_increment.png'));

rng(24680, 'twister');
[phase_ci_table, seed_stability] = local_bootstrap_statistics(compare_data.run_table);
local_plot_phase_ci(phase_ci_table, ...
    fullfile(report_dir, '07_phase_circular_variance_ci.png'));
writetable(seed_stability, fullfile(report_dir, 'seed_stability.csv'));
local_plot_seed_stability(seed_stability, ...
    fullfile(report_dir, '08_seed_stability.png'));

[angular_spectrum_summary, angular_spectrum_radial_profiles] = ...
    local_run_representative_spectra(report_dir);
writetable(angular_spectrum_radial_profiles, ...
    fullfile(report_dir, 'angular_spectrum_radial_profiles.csv'));

trend_compatibility = local_build_trend_compatibility( ...
    trend_data.summary_table, broadening_table);
local_plot_trend_compatibility(trend_compatibility, ...
    fullfile(report_dir, '11_trend_compatibility.png'));

ssa1_applicability_metrics = local_build_applicability_metrics( ...
    trend_data.summary_table);
writetable(ssa1_applicability_metrics, ...
    fullfile(report_dir, 'ssa1_applicability_metrics.csv'));
local_plot_applicability(ssa1_applicability_metrics, ...
    fullfile(report_dir, '12_ssa1_applicability_parameters.png'));

validation_scorecard = local_build_validation_scorecard( ...
    validation_data, compare_data);
writetable(validation_scorecard, ...
    fullfile(report_dir, 'validation_scorecard.csv'));

report_meta = struct();
report_meta.created_at = char(datetime('now'));
report_meta.report_dir = report_dir;
report_meta.source_files = string(required_files);
report_meta.bootstrap_rng_seed = 24680;
report_meta.bootstrap_count = 1000;
report_meta.representative_frequency_hz = 6000;
report_meta.representative_Hs_list_m = [0, 0.2, 0.5];
report_meta.representative_rough_seed_count = 8;
report_meta.notes = [ ...
    'Coherent reflected-channel loss uses abs(mean(h_reflect)), not mean(abs(h_reflect)). ', ...
    'Bootstrap intervals quantify finite-seed uncertainty. Discrete PM RMS slope is a ', ...
    'finite-grid diagnostic only and is not used as an unreferenced SSA validity threshold.'];

save(fullfile(report_dir, 'ssa1_group_meeting_summary.mat'), ...
    'validation_scorecard', 'coherent_channel_proxy', 'broadening_table', ...
    'phase_ci_table', 'seed_stability', 'angular_spectrum_summary', ...
    'angular_spectrum_radial_profiles', 'trend_compatibility', ...
    'ssa1_applicability_metrics', 'report_meta');

local_write_readme(report_dir, validation_scorecard, trend_compatibility);

fprintf('Generated SSA1 group-meeting report in %s\n', report_dir);
fprintf('Validation scorecard rows: %d\n', height(validation_scorecard));
fprintf('Representative spectrum runs: %d\n', ...
    sum(angular_spectrum_summary.seed_count));

function local_copy_existing_figure(source_file, target_file)
if ~exist(source_file, 'file')
    error('Missing source figure: %s', source_file);
end
copyfile(source_file, target_file);
end

function local_plot_energy_budget(summary_table, file_name)
data = summary_table(summary_table.kernel_mode == "ssa1_geometry", :);
freqs = unique(data.f_hz).';
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1100, 760]);
tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for ii = 1:numel(freqs)
    nexttile
    rows = sortrows(data(data.f_hz == freqs(ii), :), 'Hs_target');
    plot(rows.Hs_target, rows.E_coh_over_E_inc, '-o', 'LineWidth', 1.4, ...
        'DisplayName', 'E_{coh}/E_{inc}');
    hold on
    plot(rows.Hs_target, rows.E_sca_limit_over_E_inc, '-s', 'LineWidth', 1.4, ...
        'DisplayName', 'E_{sca}^{max}/E_{inc}');
    grid on
    ylim([0, 1.05])
    xlabel('H_s (m)')
    ylabel('Energy fraction')
    title(sprintf('f = %g Hz', freqs(ii)))
    legend('Location', 'best')
end
sgtitle('SSA1 coherent energy and available scatter budget')
exportgraphics(fig, file_name, 'Resolution', 180);
close(fig)
end

function table_out = local_build_coherent_proxy(run_table)
models = unique(run_table.model_name, 'stable').';
Hs_values = unique(run_table.Hs_target).';
freqs = unique(run_table.f_hz).';
rows = struct([]);
for im = 1:numel(models)
    for ff = 1:numel(freqs)
        flat = run_table(run_table.model_name == models(im) & ...
            run_table.f_hz == freqs(ff) & run_table.Hs_target == 0, :);
        flat_ref = abs(mean(flat.h_reflect));
        for hh = 1:numel(Hs_values)
            subset = run_table(run_table.model_name == models(im) & ...
                run_table.f_hz == freqs(ff) & ...
                run_table.Hs_target == Hs_values(hh), :);
            if isempty(subset)
                continue
            end
            proxy = abs(mean(subset.h_reflect)) / max(flat_ref, eps);
            row = struct();
            row.model_name = models(im);
            row.Hs_target = Hs_values(hh);
            row.f_hz = freqs(ff);
            row.seed_count = height(subset);
            row.abs_complex_mean_h_reflect = abs(mean(subset.h_reflect));
            row.mean_abs_h_reflect = mean(abs(subset.h_reflect));
            row.coherent_loss_proxy = proxy;
            row.coherent_loss_db = -20*log10(max(proxy, eps));
            if models(im) == "ssa1_geometry"
                row.analytic_abs_R_coh = local_nanmean(subset.abs_R_coh);
            else
                k0 = 2*pi*freqs(ff)/1500;
                sigma_eta = Hs_values(hh)/4;
                row.analytic_abs_R_coh = exp(-2*k0^2*sigma_eta^2);
            end
            rows = local_append_struct(rows, row);
        end
    end
end
table_out = struct2table(rows);
end

function local_plot_coherent_proxy(data, file_name)
freqs = unique(data.f_hz).';
models = unique(data.model_name, 'stable').';
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1100, 760]);
tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for ff = 1:numel(freqs)
    nexttile
    hold on
    for im = 1:numel(models)
        rows = sortrows(data(data.f_hz == freqs(ff) & ...
            data.model_name == models(im), :), 'Hs_target');
        plot(rows.Hs_target, rows.coherent_loss_proxy, '-o', 'LineWidth', 1.3, ...
            'DisplayName', char(models(im)));
    end
    analytic = sortrows(data(data.f_hz == freqs(ff) & ...
        data.model_name == "ssa1_geometry", :), 'Hs_target');
    plot(analytic.Hs_target, analytic.analytic_abs_R_coh, 'k--', ...
        'LineWidth', 1.5, 'DisplayName', 'analytic |R_{coh}|');
    grid on
    ylim([0, 1.05])
    xlabel('H_s (m)')
    ylabel('|E[h_{reflect}]| / |h_{flat}|')
    title(sprintf('f = %g Hz', freqs(ff)))
    legend('Location', 'best')
end
sgtitle('Coherent reflected-channel proxy from complex ensemble mean')
exportgraphics(fig, file_name, 'Resolution', 180);
close(fig)
end

function table_out = local_build_broadening_table(summary_table)
rows = summary_table(:, {'model_name', 'Hs_target', 'f_hz', ...
    'reflected_rms_delta_k_mean', 'reflected_rms_delta_k_std'});
models = unique(rows.model_name, 'stable').';
freqs = unique(rows.f_hz).';
delta = nan(height(rows), 1);
for im = 1:numel(models)
    for ff = 1:numel(freqs)
        mask = rows.model_name == models(im) & rows.f_hz == freqs(ff);
        flat = rows.reflected_rms_delta_k_mean(mask & rows.Hs_target == 0);
        delta(mask) = rows.reflected_rms_delta_k_mean(mask) - flat(1);
    end
end
rows.rms_delta_k_increment = delta;
table_out = rows;
end

function local_plot_broadening(data, file_name)
freqs = unique(data.f_hz).';
models = unique(data.model_name, 'stable').';
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1100, 760]);
tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for ff = 1:numel(freqs)
    nexttile
    hold on
    for im = 1:numel(models)
        rows = sortrows(data(data.f_hz == freqs(ff) & ...
            data.model_name == models(im), :), 'Hs_target');
        errorbar(rows.Hs_target, rows.rms_delta_k_increment, ...
            rows.reflected_rms_delta_k_std, '-o', 'LineWidth', 1.2, ...
            'DisplayName', char(models(im)));
    end
    grid on
    xlabel('H_s (m)')
    ylabel('\Delta K_{rms} from H_s=0 (rad/m)')
    title(sprintf('f = %g Hz', freqs(ff)))
    legend('Location', 'best')
end
sgtitle('Reflected angular-spectrum broadening increment')
exportgraphics(fig, file_name, 'Resolution', 180);
close(fig)
end

function [phase_table, stability_table] = local_bootstrap_statistics(run_table)
bootstrap_count = 1000;
models = unique(run_table.model_name, 'stable').';
Hs_values = unique(run_table.Hs_target).';
freqs = unique(run_table.f_hz).';
phase_rows = struct([]);
for im = 1:numel(models)
    for ff = 1:numel(freqs)
        for hh = 1:numel(Hs_values)
            subset = run_table(run_table.model_name == models(im) & ...
                run_table.f_hz == freqs(ff) & ...
                run_table.Hs_target == Hs_values(hh), :);
            values = subset.phase_h_reflect;
            boots = zeros(bootstrap_count, 1);
            for bb = 1:bootstrap_count
                idx = randi(numel(values), numel(values), 1);
                boots(bb) = local_circular_variance(values(idx));
            end
            row = struct();
            row.model_name = models(im);
            row.Hs_target = Hs_values(hh);
            row.f_hz = freqs(ff);
            row.seed_count = numel(values);
            row.circular_variance = local_circular_variance(values);
            row.ci_low = local_percentile(boots, 2.5);
            row.ci_high = local_percentile(boots, 97.5);
            phase_rows = local_append_struct(phase_rows, row);
        end
    end
end
phase_table = struct2table(phase_rows);

N_values = [2, 4, 8, 16];
target_Hs = 0.2;
target_f = 6000;
stability_rows = struct([]);
for im = 1:numel(models)
    subset = run_table(run_table.model_name == models(im) & ...
        run_table.f_hz == target_f & run_table.Hs_target == target_Hs, :);
    flat = run_table(run_table.model_name == models(im) & ...
        run_table.f_hz == target_f & run_table.Hs_target == 0, :);
    flat_ref = abs(mean(flat.h_reflect));
    for nn = 1:numel(N_values)
        N = N_values(nn);
        proxy_boot = zeros(bootstrap_count, 1);
        rms_boot = zeros(bootstrap_count, 1);
        phase_boot = zeros(bootstrap_count, 1);
        for bb = 1:bootstrap_count
            idx = randi(height(subset), N, 1);
            proxy_boot(bb) = abs(mean(subset.h_reflect(idx))) / max(flat_ref, eps);
            rms_boot(bb) = mean(subset.reflected_rms_delta_k_rad_per_m(idx));
            phase_boot(bb) = local_circular_variance(subset.phase_h_reflect(idx));
        end
        stability_rows = local_append_struct(stability_rows, ...
            local_stability_row(models(im), target_Hs, target_f, N, ...
            'coherent_proxy', proxy_boot));
        stability_rows = local_append_struct(stability_rows, ...
            local_stability_row(models(im), target_Hs, target_f, N, ...
            'rms_delta_k', rms_boot));
        stability_rows = local_append_struct(stability_rows, ...
            local_stability_row(models(im), target_Hs, target_f, N, ...
            'phase_circular_variance', phase_boot));
    end
end
stability_table = struct2table(stability_rows);
end

function row = local_stability_row(model, Hs, f_hz, N, metric, values)
row = struct();
row.model_name = model;
row.Hs_target = Hs;
row.f_hz = f_hz;
row.seed_count = N;
row.metric = string(metric);
row.bootstrap_mean = mean(values);
row.bootstrap_std = std(values);
row.ci_low = local_percentile(values, 2.5);
row.ci_high = local_percentile(values, 97.5);
end

function local_plot_phase_ci(data, file_name)
freqs = unique(data.f_hz).';
models = unique(data.model_name, 'stable').';
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1100, 760]);
tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for ff = 1:numel(freqs)
    nexttile
    hold on
    for im = 1:numel(models)
        rows = sortrows(data(data.f_hz == freqs(ff) & ...
            data.model_name == models(im), :), 'Hs_target');
        low = rows.circular_variance - rows.ci_low;
        high = rows.ci_high - rows.circular_variance;
        errorbar(rows.Hs_target, rows.circular_variance, low, high, ...
            '-o', 'LineWidth', 1.2, 'DisplayName', char(models(im)));
    end
    grid on
    ylim([0, 1.05])
    xlabel('H_s (m)')
    ylabel('Circular variance (95% bootstrap CI)')
    title(sprintf('f = %g Hz', freqs(ff)))
    legend('Location', 'best')
end
sgtitle('Finite-seed uncertainty of reflected-channel phase variance')
exportgraphics(fig, file_name, 'Resolution', 180);
close(fig)
end

function local_plot_seed_stability(data, file_name)
metrics = ["coherent_proxy", "rms_delta_k", "phase_circular_variance"];
models = unique(data.model_name, 'stable').';
labels = {'Coherent proxy', 'Mean RMS \DeltaK (rad/m)', 'Circular variance'};
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1200, 420]);
tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
for mm = 1:numel(metrics)
    nexttile
    hold on
    for im = 1:numel(models)
        rows = sortrows(data(data.metric == metrics(mm) & ...
            data.model_name == models(im), :), 'seed_count');
        low = rows.bootstrap_mean - rows.ci_low;
        high = rows.ci_high - rows.bootstrap_mean;
        errorbar(rows.seed_count, rows.bootstrap_mean, low, high, ...
            '-o', 'LineWidth', 1.2, 'DisplayName', char(models(im)));
    end
    grid on
    xlabel('Number of seeds')
    ylabel(labels{mm})
    title(sprintf('H_s=0.2 m, f=6000 Hz'))
    legend('Location', 'best')
end
sgtitle('Bootstrap seed-count stability')
exportgraphics(fig, file_name, 'Resolution', 180);
close(fig)
end

function [summary_table, radial_table] = local_run_representative_spectra(report_dir)
f_hz = 6000;
Hs_list = [0, 0.2, 0.5];
seed_list = 12345 + (0:7);
models = ["kirchhoff_spatial", "ssa1_geometry"];
grid_n = 64;
nbin = 28;
summary_rows = struct([]);
radial_rows = struct([]);
map_data = cell(numel(models), numel(Hs_list));

for im = 1:numel(models)
    for hh = 1:numel(Hs_list)
        if Hs_list(hh) == 0
            seeds = seed_list(1);
        else
            seeds = seed_list;
        end
        P_sum = zeros(grid_n, grid_n);
        radial_stack = [];
        for ss = 1:numel(seeds)
            paramsV = local_representative_params( ...
                models(im), Hs_list(hh), f_hz, seeds(ss), grid_n);
            channel = vertical_channel_model(paramsV);
            P = abs(fftshift(fft2(channel.psi_ref))).^2;
            P = P / max(sum(P(:)), eps);
            P_sum = P_sum + P;
            [KX, KY] = local_k_grid(channel.x, channel.y);
            [centers, radial] = local_radial_profile(P, KX, KY, nbin);
            radial_stack(:, ss) = radial; %#ok<AGROW>
        end
        P_mean = P_sum / numel(seeds);
        map_data{im, hh} = struct('P_mean', P_mean, 'KX', KX, 'KY', KY);
        radial_mean = mean(radial_stack, 2);
        radial_std = std(radial_stack, 0, 2);
        for bb = 1:numel(centers)
            row = struct();
            row.model_name = models(im);
            row.Hs_target = Hs_list(hh);
            row.f_hz = f_hz;
            row.seed_count = numel(seeds);
            row.k_radius_rad_per_m = centers(bb);
            row.radial_energy_fraction_mean = radial_mean(bb);
            row.radial_energy_fraction_std = radial_std(bb);
            radial_rows = local_append_struct(radial_rows, row);
        end
        srow = struct();
        srow.model_name = models(im);
        srow.Hs_target = Hs_list(hh);
        srow.f_hz = f_hz;
        srow.seed_count = numel(seeds);
        srow.rms_k_rad_per_m = sqrt(sum((KX(:).^2 + KY(:).^2) .* P_mean(:)));
        summary_rows = local_append_struct(summary_rows, srow);
    end
end

summary_table = struct2table(summary_rows);
radial_table = struct2table(radial_rows);
local_plot_spectrum_maps(map_data, models, Hs_list, ...
    fullfile(report_dir, '09_representative_angular_spectrum_maps.png'));
local_plot_radial_profiles(radial_table, models, Hs_list, ...
    fullfile(report_dir, '10_radial_angular_spectrum_profiles.png'));
end

function paramsV = local_representative_params(model, Hs, f_hz, seed, grid_n)
paramsV = struct();
paramsV.f0 = f_hz;
paramsV.c0 = 1500;
paramsV.enable_wideband = false;
paramsV.f_ref_hz = f_hz;
paramsV.z_max = 100;
paramsV.stepz_lamb = 0.5;
paramsV.xw = 32;
paramsV.yw = 32;
paramsV.nx = grid_n;
paramsV.ny = grid_n;
paramsV.x_tx = 0;
paramsV.y_tx = 0;
paramsV.z_tx = 100;
paramsV.x_rx = 0;
paramsV.y_rx = 0;
paramsV.z_rx = 3;
paramsV.nout = 2;
paramsV.sigma_src_m = 0.5;
paramsV.sponge_ratio = 0.12;
paramsV.alpha_max_np_per_m = 0.15;
paramsV.env_mode = 'uniform';
paramsV.enforce_1_over_R = false;
paramsV.show_figures = false;
paramsV.save_mode = 'rx_only';
paramsV.use_gpu = false;
paramsV.enable_surface_reflection = true;
paramsV.sea_wind_speed = 5;
paramsV.sea_hs_target = Hs;
paramsV.sea_seed = seed;
paramsV.surface_reflect_coeff = -1;
paramsV.surface_boundary_redistribution_diagnostics = true;
if model == "kirchhoff_spatial"
    paramsV.surface_boundary_model = 'kirchhoff_spatial';
else
    paramsV.surface_boundary_model = 'ssa_stat_kernel';
    paramsV.surface_ssa_kernel_mode = 'ssa1_geometry';
    paramsV.surface_ssa_random_scatter = true;
    paramsV.surface_ssa_scatter_scale = 1;
    paramsV.surface_ssa_conv_padding = 'periodic';
    paramsV.surface_ssa_geometry_source_id = ...
        'SSA.md; first-order pressure-release / Dirichlet geometry';
end
end

function [KX, KY] = local_k_grid(x, y)
dx = mean(diff(x));
dy = mean(diff(y));
nx = numel(x);
ny = numel(y);
kx = 2*pi*((-floor(nx/2)):(ceil(nx/2)-1))/(nx*dx);
ky = 2*pi*((-floor(ny/2)):(ceil(ny/2)-1))/(ny*dy);
[KX, KY] = meshgrid(kx, ky);
end

function [centers, radial] = local_radial_profile(P, KX, KY, nbin)
radius = sqrt(KX.^2 + KY.^2);
edges = linspace(0, max(radius(:)), nbin + 1);
centers = 0.5*(edges(1:end-1) + edges(2:end));
radial = zeros(nbin, 1);
for ii = 1:nbin
    if ii == nbin
        mask = radius >= edges(ii) & radius <= edges(ii + 1);
    else
        mask = radius >= edges(ii) & radius < edges(ii + 1);
    end
    radial(ii) = sum(P(mask));
end
end

function local_plot_spectrum_maps(map_data, models, Hs_list, file_name)
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80, 80, 1250, 720]);
tiledlayout(fig, numel(models), numel(Hs_list), ...
    'TileSpacing', 'compact', 'Padding', 'compact');
for im = 1:numel(models)
    for hh = 1:numel(Hs_list)
        nexttile
        item = map_data{im, hh};
        P_db = 10*log10(max(item.P_mean, eps) / max(item.P_mean(:)));
        imagesc(item.KX(1, :), item.KY(:, 1), P_db);
        axis image xy
        clim([-45, 0])
        colorbar
        xlabel('K_x (rad/m)')
        ylabel('K_y (rad/m)')
        title(sprintf('%s, H_s=%.1f m', models(im), Hs_list(hh)), ...
            'Interpreter', 'none')
    end
end
sgtitle('Representative ensemble-mean reflected angular spectra, f=6000 Hz')
exportgraphics(fig, file_name, 'Resolution', 180);
close(fig)
end

function local_plot_radial_profiles(data, models, Hs_list, file_name)
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1150, 480]);
tiledlayout(fig, 1, numel(models), 'TileSpacing', 'compact', 'Padding', 'compact');
colors = lines(numel(Hs_list));
for im = 1:numel(models)
    nexttile
    hold on
    for hh = 1:numel(Hs_list)
        rows = data(data.model_name == models(im) & ...
            data.Hs_target == Hs_list(hh), :);
        plot(rows.k_radius_rad_per_m, rows.radial_energy_fraction_mean, ...
            'LineWidth', 1.5, 'Color', colors(hh, :), ...
            'DisplayName', sprintf('H_s=%.1f m', Hs_list(hh)));
    end
    grid on
    xlabel('|K| (rad/m)')
    ylabel('Radial energy fraction per bin')
    title(char(models(im)), 'Interpreter', 'none')
    legend('Location', 'best')
end
sgtitle('Radial reflected angular-spectrum profiles, f=6000 Hz')
exportgraphics(fig, file_name, 'Resolution', 180);
close(fig)
end

function table_out = local_build_trend_compatibility(trend_summary, broadening)
freqs = unique(trend_summary.f_hz).';
rows = struct([]);
for ff = 1:numel(freqs)
    ssa = sortrows(trend_summary(trend_summary.kernel_mode == "ssa1_geometry" & ...
        trend_summary.f_hz == freqs(ff), :), 'Hs_target');
    kir = sortrows(broadening(broadening.model_name == "kirchhoff_spatial" & ...
        broadening.f_hz == freqs(ff), :), 'Hs_target');
    x = ssa.E_sca_limit_over_E_inc;
    y = kir.rms_delta_k_increment;
    rho = local_spearman(x, y);
    row = struct();
    row.f_hz = freqs(ff);
    row.spearman_rho = rho;
    row.ssa_scatter_budget_at_max_Hs = x(end);
    row.kirchhoff_rms_increment_at_max_Hs = y(end);
    row.interpretation = "descriptive trend compatibility, not model equivalence";
    rows = local_append_struct(rows, row);
end
table_out = struct2table(rows);
end

function local_plot_trend_compatibility(data, file_name)
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 900, 430]);
yyaxis left
bar(data.f_hz, data.spearman_rho, 0.45)
ylabel('Spearman \rho')
ylim([-1, 1])
yyaxis right
plot(data.f_hz, data.kirchhoff_rms_increment_at_max_Hs, '-o', ...
    'LineWidth', 1.5)
ylabel('Kirchhoff \DeltaK_{rms} at H_s=1 m (rad/m)')
grid on
xlabel('Frequency (Hz)')
title('SSA1 scatter-budget and Kirchhoff broadening trend compatibility')
exportgraphics(fig, file_name, 'Resolution', 180);
close(fig)
end

function table_out = local_build_applicability_metrics(summary_table)
ssa = summary_table(summary_table.kernel_mode == "ssa1_geometry", :);
Hs_values = unique(ssa.Hs_target).';
freqs = unique(ssa.f_hz).';
[rms_slope_by_Hs, variance_error_by_Hs] = local_discrete_pm_slope(Hs_values);
rows = struct([]);
for ff = 1:numel(freqs)
    for hh = 1:numel(Hs_values)
        subset = ssa(ssa.f_hz == freqs(ff) & ssa.Hs_target == Hs_values(hh), :);
        k0 = 2*pi*freqs(ff)/1500;
        sigma_eta = Hs_values(hh)/4;
        row = struct();
        row.Hs_target = Hs_values(hh);
        row.f_hz = freqs(ff);
        row.sigma_eta_m = sigma_eta;
        row.k0_rad_per_m = k0;
        row.k0_sigma_eta = k0*sigma_eta;
        row.abs_R_coh = subset.abs_R_coh;
        row.coherent_loss_db = -20*log10(max(subset.abs_R_coh, eps));
        row.discrete_PM_rms_slope = rms_slope_by_Hs(hh);
        row.PM_variance_rel_error = variance_error_by_Hs(hh);
        row.applicability_note = string(local_applicability_note(Hs_values(hh)));
        rows = local_append_struct(rows, row);
    end
end
table_out = struct2table(rows);
end

function [rms_slope, variance_error] = local_discrete_pm_slope(Hs_values)
nx = 128;
ny = 128;
xw = 50;
yw = 50;
U = 5;
kx = 2*pi*[0:(nx/2-1), -nx/2:-1]/xw;
ky = 2*pi*[0:(ny/2-1), -ny/2:-1]/yw;
[KX, KY] = meshgrid(kx, ky);
K = sqrt(KX.^2 + KY.^2);
g = 9.81;
alpha_PM = 8.10e-3;
beta_PM = 0.74;
Phi2D = zeros(size(K));
mask = K > 0;
E1D = zeros(size(K));
E1D(mask) = alpha_PM ./ (2*K(mask).^3) .* ...
    exp(-beta_PM*g^2 ./ (U^4*K(mask).^2));
Phi2D(mask) = E1D(mask) ./ (2*pi*K(mask));
dkx = 2*pi/xw;
dky = 2*pi/yw;
base_variance = sum(Phi2D(:))*dkx*dky;
rms_slope = zeros(numel(Hs_values), 1);
variance_error = zeros(numel(Hs_values), 1);
for ii = 1:numel(Hs_values)
    target = (Hs_values(ii)/4)^2;
    if target == 0
        W = zeros(size(Phi2D));
    else
        W = Phi2D * target / base_variance;
    end
    discrete = sum(W(:))*dkx*dky;
    rms_slope(ii) = sqrt(sum((K(:).^2).*W(:))*dkx*dky);
    variance_error(ii) = abs(discrete-target)/max(target, eps);
end
end

function note = local_applicability_note(Hs)
if Hs == 0
    note = 'flat reference';
elseif Hs < 1
    note = 'trend-validation sea state; RMS slope is diagnostic only';
else
    note = 'stress-test sea state; first-order SSA applicability not established';
end
end

function local_plot_applicability(data, file_name)
freqs = unique(data.f_hz).';
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1200, 400]);
tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile
hold on
for ff = 1:numel(freqs)
    rows = sortrows(data(data.f_hz == freqs(ff), :), 'Hs_target');
    plot(rows.Hs_target, rows.k0_sigma_eta, '-o', ...
        'DisplayName', sprintf('%g Hz', freqs(ff)));
end
grid on
xlabel('H_s (m)')
ylabel('k_0 \sigma_\eta')
title('Dimensionless roughness')
legend('Location', 'best')

nexttile
hold on
for ff = 1:numel(freqs)
    rows = sortrows(data(data.f_hz == freqs(ff), :), 'Hs_target');
    plot(rows.Hs_target, rows.coherent_loss_db, '-o', ...
        'DisplayName', sprintf('%g Hz', freqs(ff)));
end
grid on
xlabel('H_s (m)')
ylabel('Coherent loss (dB)')
title('Analytic coherent attenuation')

nexttile
rows = sortrows(data(data.f_hz == freqs(1), :), 'Hs_target');
plot(rows.Hs_target, rows.discrete_PM_rms_slope, '-o', 'LineWidth', 1.5)
grid on
xlabel('H_s (m)')
ylabel('Discrete PM RMS slope')
title('Finite-grid slope diagnostic')
sgtitle('SSA1 applicability diagnostics; no unreferenced hard threshold applied')
exportgraphics(fig, file_name, 'Resolution', 180);
close(fig)
end

function table_out = local_build_validation_scorecard(validation_data, compare_data)
checks = validation_data.summary_table;
rows = struct([]);
rows = local_score_row(rows, 'coherent_formula_max_abs_error', ...
    local_check_value(checks, 'normal_R_coh_matches_broschat_formula'), ...
    '<= 1e-12', 'validate_ssa_stat_kernel_vertical_result.mat');
rows = local_score_row(rows, 'PM_variance_max_relative_error', ...
    local_check_value(checks, 'Hs_W_eta_variance_rel_error_small'), ...
    '<= 1e-10', 'validate_ssa_stat_kernel_vertical_result.mat');
rows = local_score_row(rows, 'dense_vs_periodic_FFT_relative_error', ...
    local_check_value(checks, 'ssa1_dense_fft_P_sca_raw_sum_match'), ...
    '<= 1e-10', 'validate_ssa_stat_kernel_vertical_result.mat');
ssa_rows = compare_data.run_table(compare_data.run_table.model_name == "ssa1_geometry", :);
rows = local_score_row(rows, 'SSA1_max_energy_conservation_error', ...
    local_nanmax(ssa_rows.energy_conservation_error), ...
    '<= 1e-12', 'compare_ssa1_kirchhoff_surface_statistics_vertical_result.mat');
Hs0_checks = compare_data.trend_table(contains(compare_data.trend_table.check_name, ...
    'Hs0_h_'), :);
rows = local_score_row(rows, 'Hs0_SSA1_Kirchhoff_max_channel_difference', ...
    max(Hs0_checks.value), '<= 1e-10', ...
    'compare_ssa1_kirchhoff_surface_statistics_vertical_result.mat');
rows = local_score_row(rows, 'max_channel_invariant_error', ...
    max(compare_data.run_table.invariant_error), '<= 1e-10', ...
    'compare_ssa1_kirchhoff_surface_statistics_vertical_result.mat');
table_out = struct2table(rows);
end

function rows = local_score_row(rows, metric, value, criterion, source)
row = struct();
row.metric = string(metric);
row.value = value;
row.criterion = string(criterion);
row.source_file = string(source);
rows = local_append_struct(rows, row);
end

function value = local_check_value(checks, name)
row = checks(checks.check_name == string(name), :);
if isempty(row)
    value = NaN;
else
    value = row.value(1);
end
end

function local_write_readme(report_dir, scorecard, compatibility)
file_name = fullfile(report_dir, 'README.md');
fid = fopen(file_name, 'w', 'n', 'UTF-8');
if fid < 0
    error('Could not write %s.', file_name);
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '# SSA1 缁勪細姹囨姤鏉愭枡锛?026-06-24锛塡n\n');
fprintf(fid, '鏈洰褰曠敱 `generate_ssa1_group_meeting_report_vertical.m` 鐢熸垚銆?);
fprintf(fid, '鏍稿績妯″瀷鏈慨鏀癸紝鏂板鍐呭浠呯敤浜庢眹鎶ョ粺璁°€佸彲瑙嗗寲鍜岄€傜敤鎬у璁°€俓n\n');
fprintf(fid, '## 寤鸿姹囨姤椤哄簭\n\n');
fprintf(fid, '1. PM 璋卞綊涓€鍖栵細`03_PM_variance_normalization.png`銆俓n');
fprintf(fid, '2. 鐩稿共鍙嶅皠闅忔捣鍐靛拰棰戠巼鍙樺寲锛歚01`銆乣02`銆俓n');
fprintf(fid, '3. 鐩稿共鑳介噺涓庢暎灏勯绠楋細`04`銆俓n');
fprintf(fid, '4. 澶嶆暟闆嗗悎骞冲潎寰楀埌鐨勭浉骞查€氶亾浠ｇ悊锛歚05`銆俓n');
fprintf(fid, '5. SSA1 涓?Kirchhoff 鐨勮璋卞睍瀹藉拰鐩镐綅缁熻锛歚06`銆乣07`銆俓n');
fprintf(fid, '6. 16-seed 绋冲畾鎬э細`08`銆俓n');
fprintf(fid, '7. 浠ｈ〃鎬т簩缁磋璋变笌寰勫悜璋憋細`09`銆乣10`銆俓n');
fprintf(fid, '8. 瓒嬪娍鐩稿鎬у拰閫傜敤鎬ц竟鐣岋細`11`銆乣12`銆俓n\n');
fprintf(fid, '## 鍥捐〃璇存槑\n\n');
fprintf(fid, '- `01_Rcoh_vs_Hs.png`锛氶獙璇佺矖绯欏害澧炲ぇ鏃惰В鏋愮浉骞插弽灏勪笅闄嶃€俓n');
fprintf(fid, '- `02_Rcoh_vs_frequency.png`锛氶獙璇佸浐瀹氭捣鍐典笅楂橀 coherent loss 鏇村己銆俓n');
fprintf(fid, '- `03_PM_variance_normalization.png`锛氶獙璇佺鏁?PM 璋辩Н鍒嗙瓑浜庣洰鏍囨捣闈㈤珮搴︽柟宸€俓n');
fprintf(fid, '- `04_coherent_and_scatter_budget.png`锛氬睍绀虹浉骞茶兘閲忎笅闄嶄笌鍙敤鏁ｅ皠棰勭畻涓婂崌鐨勪簰琛ュ叧绯汇€俓n');
fprintf(fid, '- `05_coherent_channel_proxy_compare.png`锛氫娇鐢ㄥ鏁伴泦鍚堝钩鍧囨瘮杈?SSA1銆並irchhoff 鍜岃В鏋愮浉骞查」銆傚己绮楃硻搴︿笅 SSA1 鐨勯潪闆跺钩鍙颁富瑕佹槸 16-seed 闅忔満鏁ｅ皠娈嬪樊锛屼笉搴旇В閲婁负瑙ｆ瀽鐩稿共椤广€俓n');
fprintf(fid, '- `06_rms_k_broadening_increment.png`锛氭墸闄ゅ钩鏁存捣闈㈠熀绾垮悗鐨勮璋卞睍瀹姐€侹irchhoff 灞曞鏄庢樉寮轰簬褰撳墠 SSA1銆俓n');
fprintf(fid, '- `07_phase_circular_variance_ci.png`锛氬弽灏勯€氶亾 circular variance 鍙?bootstrap 95%% 鍖洪棿锛屽尯闂磋緝瀹借鏄?16 seeds 瀵圭浉浣嶉珮闃剁粺璁′粛鏈夐檺銆俓n');
fprintf(fid, '- `08_seed_stability.png`锛氬睍绀?N=2銆?銆?銆?6 鏃?bootstrap 缁熻銆俁MS 灞曞杈冪ǔ瀹氾紝鐩稿共浠ｇ悊鍜岀浉浣嶆柟宸敹鏁涙洿鎱€俓n');
fprintf(fid, '- `09_representative_angular_spectrum_maps.png`锛?000 Hz 鐨勯泦鍚堝钩鍧囦簩缁村弽灏勮璋便€侹irchhoff 绮楃硻娴烽潰鎺ヨ繎鍏ㄧ獥鍙ｆ墿鏁ｏ紝SSA1 淇濈暀鏇撮泦涓殑浣庢í鍚戞尝鏁板垎甯冦€俓n');
fprintf(fid, '- `10_radial_angular_spectrum_profiles.png`锛氫簩缁磋璋辩殑寰勫悜绉垎鐗堟湰锛屼究浜庡畾閲忔瘮杈冭氨鑳介噺鍚戦珮妯悜娉㈡暟杩佺Щ銆俓n');
fprintf(fid, '- `11_trend_compatibility.png`锛歋SA1 鏁ｅ皠棰勭畻涓?Kirchhoff 灞曞鐨勬弿杩版€хЗ鐩稿叧锛涘彧琛ㄧず瓒嬪娍鍚屽悜锛屼笉琛ㄧず鏁板€肩瓑浠枫€俓n');
fprintf(fid, '- `12_ssa1_applicability_parameters.png`锛氬睍绀?$k_0\\sigma_\\eta$銆乧oherent loss 鍜屾湁闄愮綉鏍?PM RMS slope銆傛湭浣跨敤鏃犳枃鐚緷鎹殑纭槇鍊笺€俓n\n');
fprintf(fid, '## 鏍稿績鍏紡\n\n');
fprintf(fid, '$$\\sum W_\\eta\\Delta K_x\\Delta K_y=\\sigma_\\eta^2,\\quad \\sigma_\\eta=H_s/4.$$\n\n');
fprintf(fid, '$$R_{\\rm coh}=R_0\\exp(-2k_0^2\\sigma_\\eta^2).$$\n\n');
fprintf(fid, '$$G_{\\rm SSA1}=4\\gamma_s\\gamma_i.$$\n\n');
fprintf(fid, '$$E_{\\rm coh}+E_{\\rm sca}\\le E_{\\rm inc}.$$\n\n');
fprintf(fid, 'Kirchhoff 鐨勫叡鍚岀浉骞叉寚鏍囦娇鐢╘n\n');
fprintf(fid, '$$L_{\\rm coh}=|\\mathbb{E}[h_{\\rm reflect}]|/|h_{\\rm flat}|,$$\n\n');
fprintf(fid, '鑰屼笉鏄?$\\mathbb{E}[|h_{\\rm reflect}|]$銆俓n\n');
fprintf(fid, '## 楠岃瘉璁板垎琛╘n\n');
fprintf(fid, '| 鎸囨爣 | 鏁板€?| 楠屾敹鏉′欢 |\n| --- | ---: | --- |\n');
for ii = 1:height(scorecard)
    fprintf(fid, '| %s | %.6g | %s |\n', scorecard.metric(ii), ...
        scorecard.value(ii), scorecard.criterion(ii));
end
fprintf(fid, '\n## 瓒嬪娍鐩稿鎬n\n');
fprintf(fid, '| 棰戠巼 / Hz | Spearman rho | Kirchhoff 鏈€澶ф捣鍐靛睍瀹?/ rad m^-1 |\n');
fprintf(fid, '| ---: | ---: | ---: |\n');
for ii = 1:height(compatibility)
    fprintf(fid, '| %.0f | %.4f | %.4f |\n', compatibility.f_hz(ii), ...
        compatibility.spearman_rho(ii), ...
        compatibility.kirchhoff_rms_increment_at_max_Hs(ii));
end
fprintf(fid, '\n杩欎簺鐩稿叧绯绘暟鍙弿杩版湁闄愭祴璇曠偣涓婄殑瓒嬪娍鐩稿鎬э紝涓嶈瘉鏄庝袱涓ā鍨嬬瓑浠枫€俓n\n');
fprintf(fid, '## 缁撹杈圭晫\n\n');
fprintf(fid, '- 宸查獙璇佸叕寮忛€€鍖栥€丳M 褰掍竴鍖栥€丗FT 瀹炵幇銆佽兘閲忕害鏉熷拰 16-seed 缁熻瓒嬪娍銆俓n');
fprintf(fid, '- Kirchhoff 鐨勮璋卞睍瀹介€氬父寮轰簬褰撳墠 SSA1锛屽己娴峰喌鍙嶅皠骞呭害涔熷瓨鍦ㄦ槑鏄惧樊寮傘€俓n');
fprintf(fid, '- SSA1 鍦ㄤ腑寮烘捣鍐典笅鍙兘鍥犺兘閲忛檺鍒跺嚭鐜版暎灏勮氨鍜屽弽灏勭粺璁￠ケ鍜岋紝涓嶈兘灏嗗钩鍙拌璇讳负鐪熷疄娴烽潰鏁ｅ皠宸茬粡涓嶅啀澧炲己銆俓n');
fprintf(fid, '- `surface_ssa_scatter_scale` 浠嶆槸宸ョ▼褰掍竴鍖栧弬鏁帮紝涓嶆槸瀹為獙鏍囧畾鎴潰甯告暟銆俓n');
fprintf(fid, '- $H_s=1$ m 鏄帇鍔涙祴璇曠偣锛涚鏁?RMS slope 浠呬綔璇婃柇锛屽皻鏈瘉鏄庡叾婊¤冻涓€闃?SSA 閫傜敤鏉′欢銆俓n');
fprintf(fid, '- 灏氭湭楠岃瘉鏂滃叆灏?benchmark銆佺粷瀵规暎灏勬埅闈㈠拰鐗╃悊璺ㄩ鐜囩浉鍏虫ā鍨嬨€俓n');
end

function out = local_append_struct(out, row)
if isempty(out)
    out = row;
else
    out(end + 1) = row; %#ok<AGROW>
end
end

function value = local_circular_variance(theta)
value = 1 - abs(mean(exp(1i*theta)));
value = max(value, 0);
end

function value = local_percentile(x, p)
x = sort(x(:));
if isempty(x)
    value = NaN;
    return
end
position = 1 + (numel(x)-1)*p/100;
lo = floor(position);
hi = ceil(position);
if lo == hi
    value = x(lo);
else
    value = x(lo) + (position-lo)*(x(hi)-x(lo));
end
end

function rho = local_spearman(x, y)
rx = local_rank(x(:));
ry = local_rank(y(:));
C = corrcoef(rx, ry);
rho = C(1, 2);
end

function ranks = local_rank(x)
[sorted, order] = sort(x);
ranks = zeros(size(x));
ii = 1;
while ii <= numel(x)
    jj = ii;
    while jj < numel(x) && sorted(jj + 1) == sorted(ii)
        jj = jj + 1;
    end
    ranks(order(ii:jj)) = mean(ii:jj);
    ii = jj + 1;
end
end

function value = local_nanmean(x)
x = x(isfinite(x));
if isempty(x)
    value = NaN;
else
    value = mean(x);
end
end

function value = local_nanmax(x)
x = x(isfinite(x));
if isempty(x)
    value = NaN;
else
    value = max(x);
end
end

