function model = build_surface_empirical_channel_model_vertical(result_file)
%BUILD_SURFACE_EMPIRICAL_CHANNEL_MODEL_VERTICAL Build compact empirical channel model.
%   MODEL = BUILD_SURFACE_EMPIRICAL_CHANNEL_MODEL_VERTICAL(RESULT_FILE) reads a
%   C3/C3.5 Monte Carlo result MAT-file and builds a compact empirical sample
%   table for the C4 bootstrap generator. This function does not run PE/WAPE.

if nargin < 1 || isempty(result_file)
    result_file = local_default_result_file();
end
if isstring(result_file)
    result_file = char(result_file);
end
if ~exist(result_file, 'file')
    error('build_surface_empirical_channel_model_vertical:MissingFile', ...
        'Monte Carlo result file not found: %s', result_file);
end

S = load(result_file);
[source_type, rows, condition_summary, base_params, condition_results] = local_extract_source(S);
conditions = local_build_conditions(rows, condition_summary, base_params, condition_results);

model = struct();
model.kind = 'surface_empirical_channel_model_v1';
model.source_file = result_file;
model.source_type = source_type;
model.created_at = datestr(now, 30);
model.selection_mode = 'nearest';
model.sampling_method = 'empirical_bootstrap_with_replacement';
model.description = ['C4 empirical random channel prototype. Samples are drawn ' ...
    'from existing Monte Carlo summary rows; no PE/WAPE propagation is run.'];
model.limitations = ['Nearest-neighbor sea-state matching only; valid only within ' ...
    'the scanned Monte Carlo parameter range; not a T matrix, SSA/NLSSA, ' ...
    'closed-form statistical channel model, or scattering cross-section model.'];
model.sample_fields = local_sample_field_names();
model.supports_wideband_hf = any([conditions.has_wideband_hf]);
model.supports_tap_level = model.supports_wideband_hf;
model.conditions = conditions;
model.condition_table = local_condition_table(conditions);
model.available_sea_hs_target = unique([conditions.sea_hs_target]).';
model.available_sea_wind_speed = unique([conditions.sea_wind_speed]).';
model.available_f_ref_hz = unique([conditions.f_ref_hz]).';
end

function result_file = local_default_result_file()
if exist('sweep_monte_carlo_surface_channel_vertical_result.mat', 'file')
    result_file = 'sweep_monte_carlo_surface_channel_vertical_result.mat';
elseif exist('monte_carlo_surface_channel_vertical_result.mat', 'file')
    result_file = 'monte_carlo_surface_channel_vertical_result.mat';
else
    result_file = 'sweep_monte_carlo_surface_channel_vertical_result.mat';
end
end

function [source_type, rows, condition_summary, base_params, condition_results] = local_extract_source(S)
condition_summary = table();
base_params = struct();
condition_results = struct([]);
if isfield(S, 'base_params')
    base_params = S.base_params;
elseif isfield(S, 'params_base')
    base_params = S.params_base;
end

if isfield(S, 'run_summary_table')
    source_type = 'C3.5_multi_sea_state_sweep';
    rows = S.run_summary_table;
    if isfield(S, 'condition_summary_table')
        condition_summary = S.condition_summary_table;
    end
    if isfield(S, 'condition_results')
        condition_results = S.condition_results;
    end
elseif isfield(S, 'summary_table')
    source_type = 'C3_fixed_sea_state';
    rows = S.summary_table;
    rows = local_add_fixed_condition_columns(rows, base_params);
    if isfield(S, 'mc_stats')
        condition_results = local_fixed_condition_result(S.mc_stats, base_params, rows);
    end
else
    error('build_surface_empirical_channel_model_vertical:UnsupportedFile', ...
        'Result file does not contain run_summary_table or summary_table.');
end

required = {'abs_h_ref', 'phase_h_ref_rad', 'reflect_direct_abs', ...
    'reflect_direct_phase_rad', 'redistribution_reflect_rms_delta_k_rad_per_m', ...
    'redistribution_reflect_high_k_fraction', 'tap_rms_delay_symbols'};
for ii = 1:numel(required)
    if ~ismember(required{ii}, rows.Properties.VariableNames)
        error('build_surface_empirical_channel_model_vertical:MissingColumn', ...
            'Required Monte Carlo summary column is missing: %s', required{ii});
    end
end
end

function condition_results = local_fixed_condition_result(mc_stats, base_params, rows)
condition_results = struct();
condition_results.condition_index = 1;
condition_results.name = 'fixed_sea_state';
condition_results.sea_hs_target = local_first_finite(rows.sea_hs_target);
condition_results.sea_wind_speed = local_first_finite(rows.sea_wind_speed);
condition_results.paramsV = base_params;
condition_results.seed_list = local_optional_column(rows, 'seed', NaN(height(rows), 1));
if isfield(mc_stats, 'H_f') && isfield(mc_stats.H_f, 'f_axis')
    condition_results.f_axis = mc_stats.H_f.f_axis(:);
    condition_results.idx_f_ref = mc_stats.H_f.idx_f_ref;
end
condition_results.mc_stats = mc_stats;
end

function rows = local_add_fixed_condition_columns(rows, base_params)
n = height(rows);
if ~ismember('condition_index', rows.Properties.VariableNames)
    rows.condition_index = ones(n, 1);
end
if ~ismember('sea_hs_target', rows.Properties.VariableNames)
    rows.sea_hs_target = repmat(local_get_field(base_params, 'sea_hs_target', NaN), n, 1);
end
if ~ismember('sea_wind_speed', rows.Properties.VariableNames)
    rows.sea_wind_speed = repmat(local_get_field(base_params, 'sea_wind_speed', NaN), n, 1);
end
end

function conditions = local_build_conditions(rows, condition_summary, base_params, condition_results)
if ~ismember('condition_index', rows.Properties.VariableNames)
    rows.condition_index = ones(height(rows), 1);
end

condition_ids = unique(rows.condition_index(:)).';
conditions = struct([]);
for ii = 1:numel(condition_ids)
    cid = condition_ids(ii);
    mask = rows.condition_index == cid;
    rows_i = rows(mask, :);
    samples = local_extract_samples(rows_i);

    f_ref_hz = local_condition_f_ref(cid, condition_summary, base_params);
    condition = struct();
    condition.condition_index = cid;
    condition.sea_hs_target = local_first_finite(rows_i.sea_hs_target);
    condition.sea_wind_speed = local_first_finite(rows_i.sea_wind_speed);
    condition.f_ref_hz = f_ref_hz;
    condition.mc_count = height(rows_i);
    condition.seed_list = local_optional_column(rows_i, 'seed', NaN(height(rows_i), 1));
    condition.samples = samples;
    condition.wideband = local_extract_wideband_samples(condition_results, cid);
    condition.has_wideband_hf = condition.wideband.available;
    condition.stats = local_sample_stats(samples);
    conditions = local_append_struct(conditions, condition);
end
end

function wideband = local_extract_wideband_samples(condition_results, cid)
wideband = struct( ...
    'available', false, ...
    'f_axis', [], ...
    'idx_f_ref', NaN, ...
    'H_f_samples', [], ...
    'H_direct_f_samples', [], ...
    'H_reflect_f_samples', []);
if isempty(condition_results)
    return
end

idx = [];
for ii = 1:numel(condition_results)
    if isfield(condition_results(ii), 'condition_index') && condition_results(ii).condition_index == cid
        idx = ii;
        break
    end
end
if isempty(idx)
    return
end

result = condition_results(idx);
if ~isfield(result, 'mc_stats') || ~isfield(result.mc_stats, 'H_f') || ...
        ~isfield(result.mc_stats.H_f, 'samples')
    return
end

wideband.available = true;
if isfield(result.mc_stats.H_f, 'f_axis')
    wideband.f_axis = result.mc_stats.H_f.f_axis(:);
elseif isfield(result, 'f_axis')
    wideband.f_axis = result.f_axis(:);
end
if isfield(result.mc_stats.H_f, 'idx_f_ref')
    wideband.idx_f_ref = result.mc_stats.H_f.idx_f_ref;
elseif isfield(result, 'idx_f_ref')
    wideband.idx_f_ref = result.idx_f_ref;
end
wideband.H_f_samples = result.mc_stats.H_f.samples;
if isfield(result.mc_stats, 'H_direct_f') && isfield(result.mc_stats.H_direct_f, 'samples')
    wideband.H_direct_f_samples = result.mc_stats.H_direct_f.samples;
end
if isfield(result.mc_stats, 'H_reflect_f') && isfield(result.mc_stats.H_reflect_f, 'samples')
    wideband.H_reflect_f_samples = result.mc_stats.H_reflect_f.samples;
end
end

function samples = local_extract_samples(rows_i)
samples = struct();
samples.abs_h_ref = rows_i.abs_h_ref(:);
samples.phase_h_ref_rad = rows_i.phase_h_ref_rad(:);
samples.h_ref_complex = samples.abs_h_ref .* exp(1i * samples.phase_h_ref_rad);
samples.reflect_direct_abs = rows_i.reflect_direct_abs(:);
samples.reflect_direct_phase_rad = rows_i.reflect_direct_phase_rad(:);
samples.reflect_direct_ratio_complex = samples.reflect_direct_abs .* exp(1i * samples.reflect_direct_phase_rad);
samples.coupling_nonzero_power_fraction = local_optional_column(rows_i, ...
    'coupling_nonzero_power_fraction', NaN(height(rows_i), 1));
samples.coupling_rms_delta_k_rad_per_m = local_optional_column(rows_i, ...
    'coupling_rms_delta_k_rad_per_m', NaN(height(rows_i), 1));
samples.redistribution_reflect_rms_delta_k_rad_per_m = local_optional_column(rows_i, ...
    'redistribution_reflect_rms_delta_k_rad_per_m', NaN(height(rows_i), 1));
samples.redistribution_reflect_high_k_fraction = local_optional_column(rows_i, ...
    'redistribution_reflect_high_k_fraction', NaN(height(rows_i), 1));
samples.redistribution_rms_delta_k_increase_rad_per_m = local_optional_column(rows_i, ...
    'redistribution_rms_delta_k_increase_rad_per_m', NaN(height(rows_i), 1));
samples.tap_count = local_optional_column(rows_i, 'tap_count', NaN(height(rows_i), 1));
samples.tap_rms_delay_symbols = local_optional_column(rows_i, ...
    'tap_rms_delay_symbols', NaN(height(rows_i), 1));
samples.tap_peak_fraction = local_optional_column(rows_i, ...
    'tap_peak_fraction', NaN(height(rows_i), 1));
end

function stats = local_sample_stats(samples)
fields = local_sample_field_names();
stats = struct();
for ii = 1:numel(fields)
    name = fields{ii};
    if isfield(samples, name)
        x = samples.(name);
        if ~isreal(x)
            stats.(name) = struct( ...
                'mean', mean(x), ...
                'variance', mean(abs(x - mean(x)).^2), ...
                'mean_abs', mean(abs(x)), ...
                'std_abs', std(abs(x)), ...
                'phase_circular_mean_rad', local_circular_mean(angle(x)), ...
                'phase_circular_variance', local_circular_variance(angle(x)));
        else
            stats.(name) = local_distribution_stats(x);
        end
    end
end
end

function names = local_sample_field_names()
names = { ...
    'abs_h_ref', ...
    'phase_h_ref_rad', ...
    'h_ref_complex', ...
    'reflect_direct_abs', ...
    'reflect_direct_phase_rad', ...
    'reflect_direct_ratio_complex', ...
    'coupling_nonzero_power_fraction', ...
    'coupling_rms_delta_k_rad_per_m', ...
    'redistribution_reflect_rms_delta_k_rad_per_m', ...
    'redistribution_reflect_high_k_fraction', ...
    'redistribution_rms_delta_k_increase_rad_per_m', ...
    'tap_count', ...
    'tap_rms_delay_symbols', ...
    'tap_peak_fraction'};
end

function T = local_condition_table(conditions)
n = numel(conditions);
condition_index = [conditions.condition_index].';
sea_hs_target = [conditions.sea_hs_target].';
sea_wind_speed = [conditions.sea_wind_speed].';
f_ref_hz = [conditions.f_ref_hz].';
mc_count = [conditions.mc_count].';
has_wideband_hf = [conditions.has_wideband_hf].';
T = table(condition_index, sea_hs_target, sea_wind_speed, f_ref_hz, mc_count, has_wideband_hf);
end

function f_ref_hz = local_condition_f_ref(cid, condition_summary, base_params)
f_ref_hz = local_get_field(base_params, 'f_ref_hz', NaN);
if ~isempty(condition_summary) && ismember('condition_index', condition_summary.Properties.VariableNames)
    mask = condition_summary.condition_index == cid;
    if any(mask) && ismember('f_ref_hz', condition_summary.Properties.VariableNames)
        f_ref_hz = condition_summary.f_ref_hz(find(mask, 1, 'first'));
    end
end
end

function value = local_get_field(S, name, default_value)
if isstruct(S) && isfield(S, name)
    value = S.(name);
else
    value = default_value;
end
end

function x = local_optional_column(T, name, default_value)
if ismember(name, T.Properties.VariableNames)
    x = T.(name)(:);
else
    x = default_value(:);
end
end

function value = local_first_finite(x)
x = x(:);
idx = find(isfinite(x), 1, 'first');
if isempty(idx)
    value = NaN;
else
    value = x(idx);
end
end

function out = local_append_struct(out, row)
if isempty(out)
    out = row;
else
    out(end + 1) = row; %#ok<AGROW>
end
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

function mu = local_circular_mean(theta)
z = mean(exp(1i * theta(:)));
mu = angle(z);
end

function v = local_circular_variance(theta)
z = mean(exp(1i * theta(:)));
v = 1 - abs(z);
end
