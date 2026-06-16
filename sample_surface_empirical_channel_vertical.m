function [samples, sample_meta] = sample_surface_empirical_channel_vertical(model_or_file, query, n_samples, rng_seed, sample_mode)
%SAMPLE_SURFACE_EMPIRICAL_CHANNEL_VERTICAL Draw C4 empirical channel samples.
%   [SAMPLES, META] = SAMPLE_SURFACE_EMPIRICAL_CHANNEL_VERTICAL(MODEL, QUERY,
%   N_SAMPLES, RNG_SEED) draws bootstrap samples from the nearest sea-state
%   condition in a model built by BUILD_SURFACE_EMPIRICAL_CHANNEL_MODEL_VERTICAL.
%   MODEL_OR_FILE may also be a C3/C3.5 result MAT-file path.

if nargin < 1 || isempty(model_or_file)
    model = build_surface_empirical_channel_model_vertical();
elseif ischar(model_or_file) || isstring(model_or_file)
    model = build_surface_empirical_channel_model_vertical(model_or_file);
else
    model = model_or_file;
end
if nargin < 2 || isempty(query)
    query = struct();
end
if nargin < 3 || isempty(n_samples)
    n_samples = 1000;
end
if nargin < 4 || isempty(rng_seed)
    rng_seed = 424242;
end
if nargin < 5 || isempty(sample_mode)
    sample_mode = local_get_field(query, 'sample_mode', 'summary');
end
if isstring(sample_mode)
    sample_mode = char(sample_mode);
end
sample_mode = lower(sample_mode);
n_samples = max(1, round(n_samples));

[condition, condition_distance, condition_position] = local_find_nearest_condition(model, query);
rng(rng_seed, 'twister');
sample_indices = randi(condition.mc_count, n_samples, 1);

samples = local_take_samples(condition.samples, sample_indices);
sample_table = local_samples_to_table(samples);
samples.sample_mode = sample_mode;
switch sample_mode
    case 'summary'
        % Summary bootstrap fields are already populated above.
    case 'wideband_hf'
        samples = local_add_wideband_samples(samples, condition, sample_indices);
    case 'tap_level'
        samples = local_add_wideband_samples(samples, condition, sample_indices);
        samples = local_add_tap_level_samples(samples, query);
        sample_table = local_samples_to_table(samples);
    otherwise
        error('sample_surface_empirical_channel_vertical:UnsupportedMode', ...
            'Unsupported sample_mode: %s', sample_mode);
end

sample_meta = struct();
sample_meta.generator = 'sample_surface_empirical_channel_vertical';
sample_meta.model_kind = model.kind;
sample_meta.source_file = model.source_file;
sample_meta.source_type = model.source_type;
sample_meta.selection_mode = model.selection_mode;
sample_meta.sampling_method = model.sampling_method;
sample_meta.sample_mode = sample_mode;
sample_meta.query = query;
sample_meta.condition_index = condition.condition_index;
sample_meta.condition_position = condition_position;
sample_meta.matched_hs = condition.sea_hs_target;
sample_meta.matched_wind = condition.sea_wind_speed;
sample_meta.matched_f_ref_hz = condition.f_ref_hz;
sample_meta.condition_distance = condition_distance;
sample_meta.mc_count = condition.mc_count;
sample_meta.n_samples = n_samples;
sample_meta.rng_seed = rng_seed;
sample_meta.sample_indices = sample_indices;
sample_meta.no_pe_wape_propagation_run = true;
sample_meta.has_wideband_hf = isfield(condition, 'wideband') && condition.wideband.available;
sample_meta.tap_level_derived_from_wideband_hf = strcmp(sample_mode, 'tap_level');
sample_meta.limitations = model.limitations;

samples.sample_table = sample_table;
samples.stats = local_sample_stats(samples);
end

function samples = local_add_wideband_samples(samples, condition, sample_indices)
if ~isfield(condition, 'wideband') || ~condition.wideband.available
    error('sample_surface_empirical_channel_vertical:NoWidebandSamples', ...
        'Matched condition %d has no wideband H_f samples.', condition.condition_index);
end
samples.f_axis = condition.wideband.f_axis(:);
samples.idx_f_ref = condition.wideband.idx_f_ref;
samples.H_f_samples = condition.wideband.H_f_samples(:, sample_indices);
if ~isempty(condition.wideband.H_direct_f_samples)
    samples.H_direct_f_samples = condition.wideband.H_direct_f_samples(:, sample_indices);
end
if ~isempty(condition.wideband.H_reflect_f_samples)
    samples.H_reflect_f_samples = condition.wideband.H_reflect_f_samples(:, sample_indices);
end
end

function samples = local_add_tap_level_samples(samples, query)
symbol_rate_hz = local_get_field(query, 'symbol_rate_hz', 1000);
tap_fft_len = local_get_field(query, 'tap_fft_len', 512);
tap_energy_ratio = local_get_field(query, 'tap_energy_ratio', 0.999);

n_samples = size(samples.H_f_samples, 2);
tap_cells = cell(n_samples, 1);
tap_count = NaN(n_samples, 1);
tap_energy_kept = NaN(n_samples, 1);
tap_rms_delay_symbols = NaN(n_samples, 1);
tap_peak_fraction = NaN(n_samples, 1);
for ii = 1:n_samples
    [~, H_baseband] = local_build_baseband_response( ...
        samples.f_axis, samples.H_f_samples(:, ii), samples.idx_f_ref, ...
        symbol_rate_hz, tap_fft_len);
    [h_taps, n_tap, energy_kept] = local_build_channel_taps(H_baseband, tap_energy_ratio);
    metrics = local_tap_metrics(h_taps, n_tap, energy_kept);
    tap_cells{ii} = h_taps;
    tap_count(ii) = metrics.tap_count;
    tap_energy_kept(ii) = metrics.tap_energy_kept;
    tap_rms_delay_symbols(ii) = metrics.tap_rms_delay_symbols;
    tap_peak_fraction(ii) = metrics.tap_peak_fraction;
end

samples.tap_level_symbol_rate_hz = symbol_rate_hz;
samples.tap_level_fft_len = tap_fft_len;
samples.tap_level_energy_ratio = tap_energy_ratio;
samples.tap_h_taps = tap_cells;
samples.tap_level_tap_count = tap_count;
samples.tap_level_energy_kept = tap_energy_kept;
samples.tap_level_rms_delay_symbols = tap_rms_delay_symbols;
samples.tap_level_peak_fraction = tap_peak_fraction;
end

function [condition, d_min, idx_min] = local_find_nearest_condition(model, query)
conditions = model.conditions;
if isempty(conditions)
    error('sample_surface_empirical_channel_vertical:EmptyModel', ...
        'The empirical model contains no sea-state conditions.');
end

hs = [conditions.sea_hs_target].';
wind = [conditions.sea_wind_speed].';
f_ref = [conditions.f_ref_hz].';

query_hs = local_get_field(query, 'sea_hs_target', hs(1));
query_wind = local_get_field(query, 'sea_wind_speed', wind(1));
query_f_ref = local_get_field(query, 'f_ref_hz', f_ref(1));

hs_scale = max(max(hs) - min(hs), eps);
wind_scale = max(max(wind) - min(wind), eps);
f_scale = max(max(f_ref) - min(f_ref), eps);
d = sqrt(((hs - query_hs) / hs_scale).^2 + ...
    ((wind - query_wind) / wind_scale).^2 + ...
    ((f_ref - query_f_ref) / f_scale).^2);
[d_min, idx_min] = min(d);
condition = conditions(idx_min);
end

function value = local_get_field(S, name, default_value)
if isstruct(S) && isfield(S, name) && ~isempty(S.(name))
    value = S.(name);
else
    value = default_value;
end
end

function samples = local_take_samples(pool, sample_indices)
fields = fieldnames(pool);
samples = struct();
for ii = 1:numel(fields)
    name = fields{ii};
    x = pool.(name);
    if isnumeric(x) || islogical(x)
        samples.(name) = x(sample_indices);
    end
end
end

function T = local_samples_to_table(samples)
fields = fieldnames(samples);
T = table();
n_rows = local_infer_sample_count(samples);
for ii = 1:numel(fields)
    name = fields{ii};
    x = samples.(name);
    if iscell(x) && numel(x) == n_rows
        T.(name) = x(:);
    elseif (isnumeric(x) || islogical(x)) && isvector(x) && numel(x) == n_rows
        T.(name) = x(:);
    end
end
end

function n_rows = local_infer_sample_count(samples)
if isfield(samples, 'abs_h_ref')
    n_rows = numel(samples.abs_h_ref);
else
    n_rows = 0;
end
end

function stats = local_sample_stats(samples)
fields = fieldnames(samples);
stats = struct();
n_rows = local_infer_sample_count(samples);
for ii = 1:numel(fields)
    name = fields{ii};
    if strcmp(name, 'sample_table') || strcmp(name, 'stats')
        continue
    end
    x = samples.(name);
    if isnumeric(x) && isvector(x) && numel(x) == n_rows && ~isempty(x)
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
