run(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'bootstrap_project.m'));
% Multi-sea-state Monte Carlo statistics for Kirchhoff rough-surface channels.
% This script scans Hs and wind speed, then varies only sea_seed inside each
% sea-state condition. Outputs are empirical diagnostics for the implemented
% phase-screen channel, not a new scattering model or stochastic generator.

clear
format compact

result_file = getenv('SURFACE_MC_SWEEP_RESULT_FILE');
if isempty(result_file)
    result_file = 'sweep_monte_carlo_surface_channel_vertical_result.mat';
end
figure_prefix = getenv('SURFACE_MC_SWEEP_FIGURE_PREFIX');
if isempty(figure_prefix)
    figure_prefix = 'sweep_monte_carlo_surface_channel_vertical_';
end
if exist(result_file, 'file')
    warning('Result file %s already exists and will be overwritten.', result_file);
end

sea_hs_values = [0.05, 0.5, 1.0];
sea_wind_values = [3, 5, 8, 12];
mc_count = 16;
mc_override = str2double(getenv('SURFACE_MC_SWEEP_COUNT'));
if isfinite(mc_override) && mc_override >= 1
    mc_count = round(mc_override);
end
seed_list = 12345 + (0:(mc_count - 1));

max_conditions = Inf;
max_conditions_override = str2double(getenv('SURFACE_MC_SWEEP_MAX_CONDITIONS'));
if isfinite(max_conditions_override) && max_conditions_override >= 1
    max_conditions = round(max_conditions_override);
end

symbol_rate_hz = 1000;
tap_fft_len = 512;
tap_energy_ratio = 0.999;

base_params = local_base_params();
condition_specs = local_build_condition_specs(base_params, sea_hs_values, sea_wind_values);
if isfinite(max_conditions)
    condition_specs = condition_specs(1:min(max_conditions, numel(condition_specs)));
end

condition_results = struct([]);
run_rows = struct([]);
condition_rows = struct([]);

for cc = 1:numel(condition_specs)
    [condition_result, rows_this_condition, condition_row] = local_run_condition( ...
        condition_specs(cc), cc, numel(condition_specs), seed_list, ...
        symbol_rate_hz, tap_fft_len, tap_energy_ratio);
    condition_results = local_append_struct(condition_results, condition_result);
    run_rows = local_append_struct_array(run_rows, rows_this_condition);
    condition_rows = local_append_struct(condition_rows, condition_row);
end

run_summary_table = struct2table(run_rows);
condition_summary_table = struct2table(condition_rows);
sweep_stats = local_build_sweep_stats(condition_summary_table, sea_hs_values, sea_wind_values);

disp(condition_summary_table(:, {'condition_index', 'sea_hs_target', 'sea_wind_speed', ...
    'mc_count', 'abs_H_ref_mean', 'abs_H_ref_std', ...
    'reflect_rms_delta_k_mean', 'reflect_high_k_fraction_mean', ...
    'tap_rms_delay_symbols_mean', 'max_invariant_error', 'direct_drift_max_abs'}))

local_plot_heatmap(sweep_stats, 'abs_H_ref_mean', ...
    'mean |H(f_ref)|', figure_prefix, 'abs_H_ref_mean_heatmap');
local_plot_heatmap(sweep_stats, 'abs_H_ref_std', ...
    'std |H(f_ref)|', figure_prefix, 'abs_H_ref_std_heatmap');
local_plot_heatmap(sweep_stats, 'reflect_rms_delta_k_mean', ...
    'mean reflected rms delta k', figure_prefix, 'reflect_rms_delta_k_mean_heatmap');
local_plot_heatmap(sweep_stats, 'reflect_high_k_fraction_mean', ...
    'mean reflected high-k fraction', figure_prefix, 'reflect_high_k_fraction_mean_heatmap');
local_plot_heatmap(sweep_stats, 'rms_delta_k_increase_mean', ...
    'mean rms delta k increase', figure_prefix, 'rms_delta_k_increase_mean_heatmap');
local_plot_heatmap(sweep_stats, 'tap_rms_delay_symbols_mean', ...
    'mean tap RMS delay', figure_prefix, 'tap_rms_delay_mean_heatmap');
local_plot_lines_vs_wind(sweep_stats, 'reflect_high_k_fraction_mean', ...
    'mean reflected high-k fraction vs wind', figure_prefix, 'reflect_high_k_fraction_vs_wind');
local_plot_lines_vs_wind(sweep_stats, 'reflect_rms_delta_k_mean', ...
    'mean reflected rms delta k vs wind', figure_prefix, 'reflect_rms_delta_k_vs_wind');

save(result_file, 'base_params', 'sea_hs_values', 'sea_wind_values', ...
    'seed_list', 'symbol_rate_hz', 'tap_fft_len', 'tap_energy_ratio', ...
    'condition_specs', 'condition_results', 'run_summary_table', ...
    'condition_summary_table', 'sweep_stats');

function paramsV = local_base_params()
paramsV = struct();
paramsV.f0 = 4000;
paramsV.enable_wideband = true;
paramsV.f_band_hz = [4000, 8000];
paramsV.Nf_min = 16;
paramsV.Nf_max = 16;
paramsV.f_ref_hz = 6000;
paramsV.c0 = 1500;
paramsV.z_max = 100;
paramsV.stepz_lamb = 0.5;
paramsV.xw = 50;
paramsV.yw = 50;
paramsV.nx = 128;
paramsV.ny = 128;
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
paramsV.show_figures = false;
paramsV.save_mode = 'rx_only';
paramsV.use_gpu = false;
paramsV.enable_surface_reflection = true;
paramsV.sea_wind_speed = 5.0;
paramsV.sea_hs_target = 0.5;
paramsV.sea_seed = 12345;
paramsV.surface_boundary_model = 'kirchhoff_kdomain';
paramsV.surface_boundary_coupling_diagnostics = true;
paramsV.surface_boundary_redistribution_diagnostics = true;
end

function condition_specs = local_build_condition_specs(base_params, sea_hs_values, sea_wind_values)
condition_specs = struct([]);
condition_index = 0;
for ih = 1:numel(sea_hs_values)
    for iw = 1:numel(sea_wind_values)
        condition_index = condition_index + 1;
        paramsV = base_params;
        paramsV.sea_hs_target = sea_hs_values(ih);
        paramsV.sea_wind_speed = sea_wind_values(iw);
        spec = struct();
        spec.condition_index = condition_index;
        spec.name = sprintf('Hs_%g_wind_%g', sea_hs_values(ih), sea_wind_values(iw));
        spec.sea_hs_target = sea_hs_values(ih);
        spec.sea_wind_speed = sea_wind_values(iw);
        spec.paramsV = paramsV;
        condition_specs = local_append_struct(condition_specs, spec);
    end
end
end

function [condition_result, rows, condition_row] = local_run_condition( ...
    spec, condition_number, condition_total, seed_list, symbol_rate_hz, tap_fft_len, tap_energy_ratio)

mc_count = numel(seed_list);
H_samples = [];
H_direct_samples = [];
H_reflect_samples = [];
h_total_samples = complex(zeros(mc_count, 1));
h_direct_samples = complex(zeros(mc_count, 1));
h_reflect_samples = complex(zeros(mc_count, 1));
tap_metric_rows = struct([]);
rows = struct([]);
f_axis_ref = [];
idx_f_ref_ref = NaN;

for ii = 1:mc_count
    paramsV = spec.paramsV;
    paramsV.sea_seed = seed_list(ii);
    fprintf(['Surface MC sweep condition %d/%d (%s), seed %d/%d, ' ...
        'sea_seed=%d\n'], condition_number, condition_total, spec.name, ...
        ii, mc_count, seed_list(ii));
    channel = vertical_channel_model(paramsV);

    invariant_error = max(abs(channel.H_f(:) - channel.H_direct_f(:) - channel.H_reflect_f(:)));
    if invariant_error > 1e-10
        error('sweep_monte_carlo_surface_channel_vertical:Invariant', ...
            'H_f invariant failed for condition %s seed %d: %.3e', ...
            spec.name, seed_list(ii), invariant_error);
    end

    if ii == 1
        f_axis_ref = channel.f_axis(:);
        idx_f_ref_ref = channel.idx_f_ref;
        n_freq = numel(f_axis_ref);
        H_samples = complex(zeros(n_freq, mc_count));
        H_direct_samples = complex(zeros(n_freq, mc_count));
        H_reflect_samples = complex(zeros(n_freq, mc_count));
    else
        if numel(channel.f_axis) ~= numel(f_axis_ref) || any(abs(channel.f_axis(:) - f_axis_ref) > 1e-9)
            error('Frequency axis changed for condition %s seed %d.', spec.name, seed_list(ii));
        end
        if channel.idx_f_ref ~= idx_f_ref_ref
            error('idx_f_ref changed for condition %s seed %d.', spec.name, seed_list(ii));
        end
    end

    H_samples(:, ii) = channel.H_f(:);
    H_direct_samples(:, ii) = channel.H_direct_f(:);
    H_reflect_samples(:, ii) = channel.H_reflect_f(:);
    h_total_samples(ii) = channel.h_total;
    h_direct_samples(ii) = channel.h_direct;
    h_reflect_samples(ii) = channel.h_reflect;

    [~, H_baseband] = local_build_baseband_response( ...
        channel.f_axis, channel.H_f, channel.idx_f_ref, symbol_rate_hz, tap_fft_len);
    [h_bb, tap_count, tap_energy_kept] = local_build_channel_taps(H_baseband, tap_energy_ratio);
    tap_metrics = local_tap_metrics(h_bb, tap_count, tap_energy_kept);
    tap_metric_rows = local_append_struct(tap_metric_rows, tap_metrics);

    coupling = channel.roughness_meta.boundary_coupling_diagnostics;
    redistribution = channel.roughness_meta.boundary_redistribution_diagnostics;
    if ~coupling.enabled || ~redistribution.enabled
        error('Boundary diagnostics were not enabled for condition %s seed %d.', ...
            spec.name, seed_list(ii));
    end

    ratio_ref = local_safe_complex_ratio(channel.h_reflect, channel.h_direct);
    row = struct();
    row.condition_index = spec.condition_index;
    row.sea_hs_target = spec.sea_hs_target;
    row.sea_wind_speed = spec.sea_wind_speed;
    row.seed = seed_list(ii);
    row.invariant_error = invariant_error;
    row.abs_h_ref = abs(channel.h_total);
    row.phase_h_ref_rad = angle(channel.h_total);
    row.abs_h_direct = abs(channel.h_direct);
    row.abs_h_reflect = abs(channel.h_reflect);
    row.reflect_direct_abs = abs(ratio_ref);
    row.reflect_direct_phase_rad = angle(ratio_ref);
    row.coupling_nonzero_power_fraction = coupling.nonzero_power_fraction;
    row.coupling_rms_delta_k_rad_per_m = coupling.rms_delta_k_rad_per_m;
    row.redistribution_reflect_rms_delta_k_rad_per_m = redistribution.reflect_rms_delta_k_rad_per_m;
    row.redistribution_reflect_high_k_fraction = redistribution.reflect_high_k_fraction;
    row.redistribution_rms_delta_k_increase_rad_per_m = redistribution.rms_delta_k_increase_rad_per_m;
    row.tap_count = tap_metrics.tap_count;
    row.tap_energy_kept = tap_metrics.tap_energy_kept;
    row.tap_rms_delay_symbols = tap_metrics.tap_rms_delay_symbols;
    row.tap_peak_index = tap_metrics.tap_peak_index;
    row.tap_peak_fraction = tap_metrics.tap_peak_fraction;
    rows = local_append_struct(rows, row);
end

mc_stats = local_build_mc_stats( ...
    f_axis_ref, idx_f_ref_ref, H_samples, H_direct_samples, H_reflect_samples, ...
    h_total_samples, h_direct_samples, h_reflect_samples, rows, tap_metric_rows);

condition_result = struct();
condition_result.condition_index = spec.condition_index;
condition_result.name = spec.name;
condition_result.sea_hs_target = spec.sea_hs_target;
condition_result.sea_wind_speed = spec.sea_wind_speed;
condition_result.paramsV = spec.paramsV;
condition_result.seed_list = seed_list;
condition_result.f_axis = f_axis_ref;
condition_result.idx_f_ref = idx_f_ref_ref;
condition_result.run_rows = rows;
condition_result.mc_stats = mc_stats;

condition_row = local_build_condition_row(spec, mc_stats, mc_count);
end

function row = local_build_condition_row(spec, mc_stats, mc_count)
row = struct();
row.condition_index = spec.condition_index;
row.sea_hs_target = spec.sea_hs_target;
row.sea_wind_speed = spec.sea_wind_speed;
row.mc_count = mc_count;
row.idx_f_ref = mc_stats.H_f.idx_f_ref;
row.f_ref_hz = mc_stats.H_f.f_axis(mc_stats.H_f.idx_f_ref);
row.max_invariant_error = mc_stats.validation.max_invariant_error;
row.direct_drift_max_abs = mc_stats.validation.direct_drift_max_abs;
row.abs_H_ref_mean = mc_stats.reference_frequency.abs_h_ref_mean;
row.abs_H_ref_std = mc_stats.reference_frequency.abs_h_ref_std;
row.phase_H_ref_circular_mean_rad = mc_stats.reference_frequency.phase_h_ref_circular_mean_rad;
row.phase_H_ref_circular_variance = mc_stats.reference_frequency.phase_h_ref_circular_variance;
row.reflect_direct_abs_mean = mc_stats.reflect_direct_ratio.abs_mean;
row.reflect_direct_abs_std = mc_stats.reflect_direct_ratio.abs_std;
row.reflect_direct_phase_circular_mean_rad = mc_stats.reflect_direct_ratio.phase_circular_mean_rad;
row.reflect_direct_phase_circular_variance = mc_stats.reflect_direct_ratio.phase_circular_variance;
row.coupling_nonzero_power_fraction_mean = mc_stats.diagnostics.coupling_nonzero_power_fraction.mean;
row.coupling_nonzero_power_fraction_std = mc_stats.diagnostics.coupling_nonzero_power_fraction.std;
row.coupling_rms_delta_k_mean = mc_stats.diagnostics.coupling_rms_delta_k_rad_per_m.mean;
row.coupling_rms_delta_k_std = mc_stats.diagnostics.coupling_rms_delta_k_rad_per_m.std;
row.reflect_rms_delta_k_mean = mc_stats.diagnostics.redistribution_reflect_rms_delta_k_rad_per_m.mean;
row.reflect_rms_delta_k_std = mc_stats.diagnostics.redistribution_reflect_rms_delta_k_rad_per_m.std;
row.reflect_high_k_fraction_mean = mc_stats.diagnostics.redistribution_reflect_high_k_fraction.mean;
row.reflect_high_k_fraction_std = mc_stats.diagnostics.redistribution_reflect_high_k_fraction.std;
row.rms_delta_k_increase_mean = mc_stats.diagnostics.redistribution_rms_delta_k_increase_rad_per_m.mean;
row.rms_delta_k_increase_std = mc_stats.diagnostics.redistribution_rms_delta_k_increase_rad_per_m.std;
row.tap_count_mean = mc_stats.tap_metrics.tap_count.mean;
row.tap_count_std = mc_stats.tap_metrics.tap_count.std;
row.tap_energy_kept_mean = mc_stats.tap_metrics.tap_energy_kept.mean;
row.tap_rms_delay_symbols_mean = mc_stats.tap_metrics.tap_rms_delay_symbols.mean;
row.tap_rms_delay_symbols_std = mc_stats.tap_metrics.tap_rms_delay_symbols.std;
row.tap_peak_fraction_mean = mc_stats.tap_metrics.tap_peak_fraction.mean;
row.tap_peak_fraction_std = mc_stats.tap_metrics.tap_peak_fraction.std;
end

function out = local_append_struct(out, row)
if isempty(out)
    out = row;
else
    out(end + 1) = row; %#ok<AGROW>
end
end

function out = local_append_struct_array(out, rows)
for ii = 1:numel(rows)
    out = local_append_struct(out, rows(ii));
end
end

function ratio = local_safe_complex_ratio(num, den)
if abs(den) <= eps
    den = complex(eps, 0);
end
ratio = num ./ den;
end

function [f_bb_axis, H_baseband_shifted] = local_build_baseband_response(f_axis, H_f, idx_f_ref, fs_hz, n_fft)
f_axis = f_axis(:);
H_f = H_f(:);
if numel(f_axis) ~= numel(H_f)
    error('f_axis and H_f length mismatch.');
end
if idx_f_ref < 1 || idx_f_ref > numel(f_axis)
    error('idx_f_ref out of range.');
end

fc = f_axis(idx_f_ref);
f_rel = f_axis - fc;
f_bb_axis = ((0:n_fft-1).' - floor(n_fft/2)) * (fs_hz / n_fft);
H_baseband_shifted = interp1(f_rel, H_f, f_bb_axis, 'linear', 0);
end

function [h_eff, n_tap_eff, energy_kept] = local_build_channel_taps(H_baseband_shifted, energy_ratio)
H_baseband_shifted = H_baseband_shifted(:);
h_full = ifft(ifftshift(H_baseband_shifted));
energy = abs(h_full).^2;
if all(energy == 0)
    h_eff = complex(0, 0);
    n_tap_eff = 1;
    energy_kept = 0;
    return
end

cum_energy = cumsum(energy);
target = max(min(energy_ratio, 1), 0);
n_tap_eff = find(cum_energy >= target * cum_energy(end), 1, 'first');
n_tap_eff = max(1, n_tap_eff);
h_eff = h_full(1:n_tap_eff);
energy_kept = sum(abs(h_eff).^2) / sum(abs(h_full).^2);
end

function metrics = local_tap_metrics(h_bb, tap_count, energy_kept)
h_bb = h_bb(:);
tap_energy = abs(h_bb).^2;
total_energy = sum(tap_energy);
if total_energy <= 0
    mean_delay = NaN;
    rms_delay = NaN;
    peak_index = NaN;
    peak_fraction = NaN;
else
    delay_idx = (0:(numel(h_bb) - 1)).';
    mean_delay = sum(delay_idx .* tap_energy) / total_energy;
    rms_delay = sqrt(sum(((delay_idx - mean_delay).^2) .* tap_energy) / total_energy);
    [peak_energy, peak_index] = max(tap_energy);
    peak_fraction = peak_energy / total_energy;
end

metrics = struct( ...
    'tap_count', tap_count, ...
    'tap_energy_kept', energy_kept, ...
    'tap_mean_delay_symbols', mean_delay, ...
    'tap_rms_delay_symbols', rms_delay, ...
    'tap_peak_index', peak_index, ...
    'tap_peak_fraction', peak_fraction);
end

function mc_stats = local_build_mc_stats( ...
    f_axis, idx_f_ref, H_samples, H_direct_samples, H_reflect_samples, ...
    h_total_samples, h_direct_samples, h_reflect_samples, rows, tap_rows)

mean_f = mean(H_samples, 2);
var_f = mean(abs(H_samples - mean_f).^2, 2);
direct_mean_f = mean(H_direct_samples, 2);
direct_drift_max_abs = max(max(abs(H_direct_samples - direct_mean_f)));

abs_h_ref = abs(h_total_samples);
phase_h_ref = angle(h_total_samples);
ratio_ref = arrayfun(@local_safe_complex_ratio, h_reflect_samples, h_direct_samples);
ratio_abs = abs(ratio_ref);
ratio_phase = angle(ratio_ref);

mc_stats = struct();
mc_stats.H_f = struct( ...
    'f_axis', f_axis, ...
    'idx_f_ref', idx_f_ref, ...
    'samples', H_samples, ...
    'mean_f', mean_f, ...
    'var_f', var_f, ...
    'mean_abs_f', mean(abs(H_samples), 2), ...
    'std_abs_f', std(abs(H_samples), 0, 2));
mc_stats.H_direct_f = struct( ...
    'samples', H_direct_samples, ...
    'mean_f', direct_mean_f, ...
    'direct_drift_max_abs', direct_drift_max_abs);
mc_stats.H_reflect_f = struct( ...
    'samples', H_reflect_samples, ...
    'mean_f', mean(H_reflect_samples, 2), ...
    'var_f', mean(abs(H_reflect_samples - mean(H_reflect_samples, 2)).^2, 2));
mc_stats.reference_frequency = struct( ...
    'h_total_samples', h_total_samples, ...
    'abs_h_ref_mean', mean(abs_h_ref), ...
    'abs_h_ref_std', std(abs_h_ref), ...
    'abs_h_ref_quantiles_5_25_50_75_95', local_quantiles(abs_h_ref, [5, 25, 50, 75, 95]), ...
    'phase_h_ref_samples_rad', phase_h_ref, ...
    'phase_h_ref_circular_mean_rad', local_circular_mean(phase_h_ref), ...
    'phase_h_ref_circular_variance', local_circular_variance(phase_h_ref), ...
    'h_total_mean', mean(h_total_samples), ...
    'h_total_var', mean(abs(h_total_samples - mean(h_total_samples)).^2));
mc_stats.reflect_direct_ratio = struct( ...
    'samples', ratio_ref, ...
    'abs_mean', mean(ratio_abs), ...
    'abs_std', std(ratio_abs), ...
    'abs_quantiles_5_25_50_75_95', local_quantiles(ratio_abs, [5, 25, 50, 75, 95]), ...
    'phase_circular_mean_rad', local_circular_mean(ratio_phase), ...
    'phase_circular_variance', local_circular_variance(ratio_phase));
mc_stats.diagnostics = struct( ...
    'coupling_nonzero_power_fraction', local_distribution_stats([rows.coupling_nonzero_power_fraction].'), ...
    'coupling_rms_delta_k_rad_per_m', local_distribution_stats([rows.coupling_rms_delta_k_rad_per_m].'), ...
    'redistribution_reflect_rms_delta_k_rad_per_m', local_distribution_stats([rows.redistribution_reflect_rms_delta_k_rad_per_m].'), ...
    'redistribution_reflect_high_k_fraction', local_distribution_stats([rows.redistribution_reflect_high_k_fraction].'), ...
    'redistribution_rms_delta_k_increase_rad_per_m', local_distribution_stats([rows.redistribution_rms_delta_k_increase_rad_per_m].'));
mc_stats.tap_metrics = struct( ...
    'tap_count', local_distribution_stats([tap_rows.tap_count].'), ...
    'tap_energy_kept', local_distribution_stats([tap_rows.tap_energy_kept].'), ...
    'tap_rms_delay_symbols', local_distribution_stats([tap_rows.tap_rms_delay_symbols].'), ...
    'tap_peak_index', local_distribution_stats([tap_rows.tap_peak_index].'), ...
    'tap_peak_fraction', local_distribution_stats([tap_rows.tap_peak_fraction].'));
mc_stats.validation = struct( ...
    'max_invariant_error', max([rows.invariant_error]), ...
    'direct_drift_max_abs', direct_drift_max_abs);
end

function stats = local_distribution_stats(x)
x = x(:);
x = x(isfinite(x));
if isempty(x)
    stats = struct('mean', NaN, 'var', NaN, 'std', NaN, ...
        'quantiles_5_25_50_75_95', NaN(1, 5));
    return
end
mu = mean(x);
stats = struct( ...
    'mean', mu, ...
    'var', mean((x - mu).^2), ...
    'std', std(x), ...
    'quantiles_5_25_50_75_95', local_quantiles(x, [5, 25, 50, 75, 95]));
end

function q = local_quantiles(x, pct)
x = sort(x(:));
pct = pct(:).';
if isempty(x)
    q = NaN(size(pct));
    return
end
if numel(x) == 1
    q = repmat(x, size(pct));
    return
end
pos = 1 + (pct / 100) * (numel(x) - 1);
q = interp1(1:numel(x), x, pos, 'linear');
end

function mu = local_circular_mean(theta)
z = mean(exp(1i * theta(:)));
mu = angle(z);
end

function v = local_circular_variance(theta)
z = mean(exp(1i * theta(:)));
v = 1 - abs(z);
end

function sweep_stats = local_build_sweep_stats(condition_summary_table, sea_hs_values, sea_wind_values)
sweep_stats = struct();
sweep_stats.sea_hs_values = sea_hs_values(:);
sweep_stats.sea_wind_values = sea_wind_values(:).';
metric_names = { ...
    'abs_H_ref_mean', ...
    'abs_H_ref_std', ...
    'reflect_rms_delta_k_mean', ...
    'reflect_high_k_fraction_mean', ...
    'rms_delta_k_increase_mean', ...
    'tap_rms_delay_symbols_mean', ...
    'max_invariant_error', ...
    'direct_drift_max_abs'};
for im = 1:numel(metric_names)
    sweep_stats.(metric_names{im}) = NaN(numel(sea_hs_values), numel(sea_wind_values));
end
sweep_stats.condition_mask = false(numel(sea_hs_values), numel(sea_wind_values));

for ir = 1:height(condition_summary_table)
    ih = find(abs(sea_hs_values - condition_summary_table.sea_hs_target(ir)) <= 1e-12, 1, 'first');
    iw = find(abs(sea_wind_values - condition_summary_table.sea_wind_speed(ir)) <= 1e-12, 1, 'first');
    if isempty(ih) || isempty(iw)
        continue
    end
    sweep_stats.condition_mask(ih, iw) = true;
    for im = 1:numel(metric_names)
        sweep_stats.(metric_names{im})(ih, iw) = condition_summary_table.(metric_names{im})(ir);
    end
end
end

function local_plot_heatmap(sweep_stats, field_name, title_text, figure_prefix, suffix)
values = sweep_stats.(field_name);
if ~any(isfinite(values(:)))
    return
end
figure('Visible', 'off');
imagesc(sweep_stats.sea_wind_values, sweep_stats.sea_hs_values, values)
set(gca, 'YDir', 'normal')
colorbar
grid on
xlabel('sea wind speed (m/s)')
ylabel('sea Hs target (m)')
title(title_text, 'Interpreter', 'none')
print(gcf, '-dpng', '-r200', [figure_prefix suffix '.png'])
close(gcf)
end

function local_plot_lines_vs_wind(sweep_stats, field_name, title_text, figure_prefix, suffix)
values = sweep_stats.(field_name);
if ~any(isfinite(values(:)))
    return
end
figure('Visible', 'off');
hold on
legend_labels = {};
for ih = 1:numel(sweep_stats.sea_hs_values)
    y = values(ih, :);
    if any(isfinite(y))
        plot(sweep_stats.sea_wind_values, y, 'o-', 'LineWidth', 1.2)
        legend_labels{end + 1} = sprintf('Hs=%g m', sweep_stats.sea_hs_values(ih)); %#ok<AGROW>
    end
end
grid on
xlabel('sea wind speed (m/s)')
ylabel(field_name, 'Interpreter', 'none')
title(title_text, 'Interpreter', 'none')
if ~isempty(legend_labels)
    legend(legend_labels, 'Location', 'best')
end
print(gcf, '-dpng', '-r200', [figure_prefix suffix '.png'])
close(gcf)
end

