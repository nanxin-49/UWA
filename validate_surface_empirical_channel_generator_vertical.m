% Validate C4 empirical random channel generator without running PE/WAPE.
% This script reads existing C3/C3.5 Monte Carlo result files, builds a
% compact empirical bootstrap model, and compares generated samples against
% the source Monte Carlo samples.

clear
format compact

result_file = local_default_result_file();
output_file = 'surface_empirical_channel_generator_validation_result.mat';
figure_prefix = 'surface_empirical_channel_generator_';

n_fast_samples = 1000;
n_override = str2double(getenv('SURFACE_EMPIRICAL_SAMPLE_COUNT'));
if isfinite(n_override) && n_override >= 1
    n_fast_samples = round(n_override);
end

max_conditions = 3;
max_conditions_override = str2double(getenv('SURFACE_EMPIRICAL_MAX_CONDITIONS'));
if isfinite(max_conditions_override) && max_conditions_override >= 1
    max_conditions = round(max_conditions_override);
end

model = build_surface_empirical_channel_model_vertical(result_file);
condition_positions = local_select_condition_positions(numel(model.conditions), max_conditions);
metric_names = { ...
    'abs_h_ref', ...
    'phase_h_ref_rad', ...
    'reflect_direct_abs', ...
    'redistribution_reflect_rms_delta_k_rad_per_m', ...
    'redistribution_reflect_high_k_fraction', ...
    'tap_rms_delay_symbols'};

comparison_rows = struct([]);
sample_sets = struct([]);
mode_validation = struct([]);
for ii = 1:numel(condition_positions)
    condition = model.conditions(condition_positions(ii));
    query = struct( ...
        'sea_hs_target', condition.sea_hs_target, ...
        'sea_wind_speed', condition.sea_wind_speed, ...
        'f_ref_hz', condition.f_ref_hz);
    [generated_samples, sample_meta] = sample_surface_empirical_channel_vertical( ...
        model, query, n_fast_samples, 50000 + ii);

    for im = 1:numel(metric_names)
        metric_name = metric_names{im};
        row = local_compare_metric(condition, generated_samples, sample_meta, metric_name);
        comparison_rows = local_append_struct(comparison_rows, row);
    end

    sample_set = struct();
    sample_set.query = query;
    sample_set.sample_meta = sample_meta;
    sample_set.generated_sample_table = generated_samples.sample_table;
    sample_sets = local_append_struct(sample_sets, sample_set);

    if ii == 1
        mode_validation = local_validate_wideband_modes(model, query, min(n_fast_samples, 50));
    end
end

comparison_table = struct2table(comparison_rows);
validation_meta = struct();
validation_meta.result_file = result_file;
validation_meta.output_file = output_file;
validation_meta.n_fast_samples = n_fast_samples;
validation_meta.condition_positions = condition_positions;
validation_meta.metric_names = metric_names;
validation_meta.no_pe_wape_propagation_run = true;
validation_meta.full_channel_structs_saved = false;
validation_meta.large_spectral_arrays_saved = false;
validation_meta.wideband_mode_checked = ~isempty(mode_validation);

disp(model.condition_table)
disp(comparison_table(:, {'condition_index', 'sea_hs_target', 'sea_wind_speed', ...
    'metric_name', 'source_mean', 'generated_mean', 'abs_mean_diff', ...
    'source_std', 'generated_std', 'abs_std_diff'}))

local_plot_selected_histograms(model, sample_sets, metric_names, figure_prefix);

save(output_file, 'result_file', 'model', 'comparison_table', ...
    'validation_meta', 'sample_sets', 'mode_validation');
file_info = dir(output_file);
validation_meta.output_file_size_bytes = file_info.bytes;
save(output_file, 'result_file', 'model', 'comparison_table', ...
    'validation_meta', 'sample_sets', 'mode_validation');

fprintf('C4 empirical generator validation saved %s (%d bytes).\n', ...
    output_file, validation_meta.output_file_size_bytes);

function result_file = local_default_result_file()
if exist('sweep_monte_carlo_surface_channel_vertical_result.mat', 'file')
    result_file = 'sweep_monte_carlo_surface_channel_vertical_result.mat';
elseif exist('monte_carlo_surface_channel_vertical_result.mat', 'file')
    result_file = 'monte_carlo_surface_channel_vertical_result.mat';
else
    error('validate_surface_empirical_channel_generator_vertical:MissingInput', ...
        ['No C3/C3.5 Monte Carlo result file found. Expected ' ...
        'sweep_monte_carlo_surface_channel_vertical_result.mat or ' ...
        'monte_carlo_surface_channel_vertical_result.mat.']);
end
end

function positions = local_select_condition_positions(n_conditions, max_conditions)
if n_conditions <= max_conditions
    positions = 1:n_conditions;
    return
end
positions = unique(round(linspace(1, n_conditions, max_conditions)));
end

function row = local_compare_metric(condition, generated_samples, sample_meta, metric_name)
source_x = condition.samples.(metric_name);
generated_x = generated_samples.(metric_name);
source_stats = local_distribution_stats(source_x);
generated_stats = local_distribution_stats(generated_x);

row = struct();
row.condition_index = condition.condition_index;
row.sea_hs_target = condition.sea_hs_target;
row.sea_wind_speed = condition.sea_wind_speed;
row.f_ref_hz = condition.f_ref_hz;
row.metric_name = metric_name;
row.mc_count = condition.mc_count;
row.n_fast_samples = sample_meta.n_samples;
row.source_mean = source_stats.mean;
row.generated_mean = generated_stats.mean;
row.abs_mean_diff = abs(generated_stats.mean - source_stats.mean);
row.source_variance = source_stats.variance;
row.generated_variance = generated_stats.variance;
row.abs_variance_diff = abs(generated_stats.variance - source_stats.variance);
row.source_std = source_stats.std;
row.generated_std = generated_stats.std;
row.abs_std_diff = abs(generated_stats.std - source_stats.std);
row.source_q05 = source_stats.quantiles_5_25_50_75_95(1);
row.source_q25 = source_stats.quantiles_5_25_50_75_95(2);
row.source_q50 = source_stats.quantiles_5_25_50_75_95(3);
row.source_q75 = source_stats.quantiles_5_25_50_75_95(4);
row.source_q95 = source_stats.quantiles_5_25_50_75_95(5);
row.generated_q05 = generated_stats.quantiles_5_25_50_75_95(1);
row.generated_q25 = generated_stats.quantiles_5_25_50_75_95(2);
row.generated_q50 = generated_stats.quantiles_5_25_50_75_95(3);
row.generated_q75 = generated_stats.quantiles_5_25_50_75_95(4);
row.generated_q95 = generated_stats.quantiles_5_25_50_75_95(5);
row.abs_q50_diff = abs(row.generated_q50 - row.source_q50);
end

function mode_validation = local_validate_wideband_modes(model, query, n_samples)
mode_validation = struct();
if ~isfield(model, 'supports_wideband_hf') || ~model.supports_wideband_hf
    mode_validation.available = false;
    mode_validation.reason = 'model_has_no_wideband_hf_samples';
    return
end

[wideband_samples, wideband_meta] = sample_surface_empirical_channel_vertical( ...
    model, query, n_samples, 61001, 'wideband_hf');
[tap_samples, tap_meta] = sample_surface_empirical_channel_vertical( ...
    model, query, n_samples, 61002, 'tap_level');

mode_validation.available = true;
mode_validation.wideband_sample_mode = wideband_meta.sample_mode;
mode_validation.tap_sample_mode = tap_meta.sample_mode;
mode_validation.n_samples = n_samples;
mode_validation.n_freq = numel(wideband_samples.f_axis);
mode_validation.H_f_sample_size = size(wideband_samples.H_f_samples);
mode_validation.idx_f_ref = wideband_samples.idx_f_ref;
mode_validation.tap_count_mean = mean(tap_samples.tap_level_tap_count, 'omitnan');
mode_validation.tap_rms_delay_mean = mean(tap_samples.tap_level_rms_delay_symbols, 'omitnan');
mode_validation.tap_peak_fraction_mean = mean(tap_samples.tap_level_peak_fraction, 'omitnan');
mode_validation.tap_metrics_all_finite = all(isfinite(tap_samples.tap_level_tap_count)) && ...
    all(isfinite(tap_samples.tap_level_rms_delay_symbols)) && ...
    all(isfinite(tap_samples.tap_level_peak_fraction));
mode_validation.no_pe_wape_propagation_run = wideband_meta.no_pe_wape_propagation_run && ...
    tap_meta.no_pe_wape_propagation_run;
end

function stats = local_distribution_stats(x)
x = x(:);
x = x(isfinite(x));
if isempty(x)
    stats = struct('mean', NaN, 'variance', NaN, 'std', NaN, ...
        'quantiles_5_25_50_75_95', NaN(1, 5));
    return
end
mu = mean(x);
stats = struct( ...
    'mean', mu, ...
    'variance', mean((x - mu).^2), ...
    'std', std(x), ...
    'quantiles_5_25_50_75_95', local_quantiles(x, [5, 25, 50, 75, 95]));
end

function q = local_quantiles(x, pct)
x = sort(x(:));
pct = pct(:).';
if isempty(x)
    q = NaN(size(pct));
elseif numel(x) == 1
    q = repmat(x, size(pct));
else
    pos = 1 + (pct / 100) * (numel(x) - 1);
    q = interp1(1:numel(x), x, pos, 'linear');
end
end

function local_plot_selected_histograms(model, sample_sets, metric_names, figure_prefix)
if isempty(sample_sets)
    return
end
condition = model.conditions(sample_sets(1).sample_meta.condition_position);
generated_table = sample_sets(1).generated_sample_table;
for ii = 1:numel(metric_names)
    metric_name = metric_names{ii};
    source_x = condition.samples.(metric_name);
    generated_x = generated_table.(metric_name);

    figure('Visible', 'off');
    hold on
    histogram(source_x, 'Normalization', 'probability')
    histogram(generated_x, 'Normalization', 'probability')
    grid on
    xlabel(metric_name, 'Interpreter', 'none')
    ylabel('probability')
    title(sprintf('C4 empirical generator: %s', metric_name), 'Interpreter', 'none')
    legend('source MC', 'generated bootstrap', 'Location', 'best')
    print(gcf, '-dpng', '-r200', [figure_prefix metric_name '_hist.png'])
    close(gcf)
end
end

function out = local_append_struct(out, row)
if isempty(out)
    out = row;
else
    out(end + 1) = row; %#ok<AGROW>
end
end
