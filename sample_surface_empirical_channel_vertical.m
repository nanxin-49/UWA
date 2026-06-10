function [samples, sample_meta] = sample_surface_empirical_channel_vertical(model_or_file, query, n_samples, rng_seed)
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
n_samples = max(1, round(n_samples));

[condition, condition_distance, condition_position] = local_find_nearest_condition(model, query);
rng(rng_seed, 'twister');
sample_indices = randi(condition.mc_count, n_samples, 1);

samples = local_take_samples(condition.samples, sample_indices);
sample_table = local_samples_to_table(samples);

sample_meta = struct();
sample_meta.generator = 'sample_surface_empirical_channel_vertical';
sample_meta.model_kind = model.kind;
sample_meta.source_file = model.source_file;
sample_meta.source_type = model.source_type;
sample_meta.selection_mode = model.selection_mode;
sample_meta.sampling_method = model.sampling_method;
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
sample_meta.limitations = model.limitations;

samples.sample_table = sample_table;
samples.stats = local_sample_stats(samples);
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
for ii = 1:numel(fields)
    name = fields{ii};
    x = samples.(name);
    if isnumeric(x) || islogical(x)
        T.(name) = x(:);
    end
end
end

function stats = local_sample_stats(samples)
fields = fieldnames(samples);
stats = struct();
for ii = 1:numel(fields)
    name = fields{ii};
    if strcmp(name, 'sample_table') || strcmp(name, 'stats')
        continue
    end
    x = samples.(name);
    if isnumeric(x) && ~isempty(x)
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
