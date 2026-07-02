run(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'bootstrap_project.m'));
% Validate optional frequency-correlated random scatter for ssa_stat_kernel.
% This script isolates the random reflected contribution by subtracting the
% random_scatter=false reflected channel from the random_scatter=true result.

clear
format compact

result_file = getenv('SSA_FREQ_CORR_RESULT_FILE');
if isempty(result_file)
    result_file = 'validate_ssa_frequency_correlation_vertical_result.mat';
end

mode_list = {'independent', 'shared_seed_phase', 'ar1_frequency'};
seed_count = local_env_scalar_int('SSA_FREQ_CORR_SEED_COUNT', 4);
grid_n = local_env_scalar_int('SSA_FREQ_CORR_GRID_N', 64);
seed_list = 12345 + (0:(seed_count - 1));
rho = local_env_scalar('SSA_FREQ_CORR_RHO', 0.85);

run_rows = struct([]);
run_index = 0;
coherent_params = local_base_params(grid_n, rho);
coherent_params.surface_ssa_random_scatter = false;
coherent_channel = vertical_channel_model(coherent_params);

for im = 1:numel(mode_list)
    for ss = 1:numel(seed_list)
        paramsV = local_base_params(grid_n, rho);
        paramsV.surface_ssa_frequency_correlation_mode = mode_list{im};
        paramsV.sea_seed = seed_list(ss);

        fprintf('Frequency correlation validation: mode=%s, seed=%d.\n', ...
            mode_list{im}, paramsV.sea_seed);
        channel = vertical_channel_model(paramsV);

        repeat_channel = vertical_channel_model(paramsV);
        random_reflect_f = channel.H_reflect_f - coherent_channel.H_reflect_f;
        repeat_error = max(abs(channel.H_f(:) - repeat_channel.H_f(:)));
        direct_error = max(abs(channel.H_direct_f(:) - coherent_channel.H_direct_f(:)));

        run_index = run_index + 1;
        run_rows = local_append_struct(run_rows, local_run_row( ...
            run_index, paramsV, channel, random_reflect_f, repeat_error, direct_error));
    end
end

run_table = struct2table(run_rows);
summary_table = local_build_summary_table(run_table, string(mode_list));
independent_ref = summary_table.adjacent_random_reflect_coherence_abs_mean( ...
    summary_table.frequency_correlation_mode == "independent");
shared_ref = summary_table.adjacent_random_reflect_coherence_abs_mean( ...
    summary_table.frequency_correlation_mode == "shared_seed_phase");
ar1_ref = summary_table.adjacent_random_reflect_coherence_abs_mean( ...
    summary_table.frequency_correlation_mode == "ar1_frequency");

checks = struct([]);
checks = local_add_check(checks, 'channel_invariant_all', ...
    max(run_table.invariant_error), 1e-10, '<=');
checks = local_add_check(checks, 'repeatability_all', ...
    max(run_table.repeat_error), 1e-12, '<=');
checks = local_add_check(checks, 'direct_path_stable_all', ...
    max(run_table.direct_error_vs_metadata_only), 1e-12, '<=');
checks = local_add_check(checks, 'shared_seed_phase_increases_random_frequency_coherence', ...
    shared_ref - independent_ref, 0, '>=');
checks = local_add_check(checks, 'ar1_frequency_increases_random_frequency_coherence', ...
    ar1_ref - independent_ref, 0, '>=');

trend_table = struct2table(checks);
validation_report = struct();
validation_report.script = mfilename;
validation_report.created_at = char(datetime('now'));
validation_report.mode_list = string(mode_list);
validation_report.seed_list = seed_list;
validation_report.grid_n = grid_n;
validation_report.rho = rho;
validation_report.all_passed = all(trend_table.passed);
validation_report.notes = ['frequency correlation modes are optional and default to independent. ', ...
    'The random reflected contribution is approximated by subtracting the random_scatter=false reflected response.'];

local_plot_summary(summary_table);

save(result_file, 'run_table', 'summary_table', 'trend_table', ...
    'validation_report', 'mode_list', 'seed_list', 'grid_n', 'rho');

disp(summary_table)
disp(trend_table)
fprintf('Saved %s\n', result_file);
fprintf('All frequency-correlation checks passed: %d\n', validation_report.all_passed);
if ~validation_report.all_passed
    error('validate_ssa_frequency_correlation_vertical:Failed', ...
        'One or more frequency-correlation checks failed.');
end

function paramsV = local_base_params(grid_n, rho)
domain_width_m = min(50, 0.5 * grid_n);
paramsV = struct();
paramsV.f0 = 6000;
paramsV.c0 = 1500;
paramsV.enable_wideband = true;
paramsV.f_band_hz = [4000, 8000];
paramsV.Nf_min = 16;
paramsV.Nf_max = 16;
paramsV.f_ref_hz = 6000;
paramsV.z_max = 100;
paramsV.stepz_lamb = 0.5;
paramsV.xw = domain_width_m;
paramsV.yw = domain_width_m;
paramsV.nx = grid_n;
paramsV.ny = grid_n;
paramsV.x_tx = 0;
paramsV.y_tx = 0;
paramsV.z_tx = 100;
paramsV.x_rx = 0;
paramsV.y_rx = 0;
paramsV.z_rx = 3;
paramsV.nout = 2;
paramsV.sigma_src_m = max(0.4, domain_width_m / grid_n);
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
paramsV.surface_reflect_coeff = -1;
paramsV.surface_phase_mode = 'normal';
paramsV.surface_boundary_model = 'ssa_stat_kernel';
paramsV.surface_ssa_random_scatter = true;
paramsV.surface_ssa_scatter_scale = 4.0;
paramsV.surface_ssa_seed_offset = 100000;
paramsV.surface_ssa_kernel_mode = 'ssa1_geometry';
paramsV.surface_ssa_coherent_order = 'ssa1';
paramsV.surface_ssa_frequency_correlation_mode = 'independent';
paramsV.surface_ssa_frequency_correlation_rho = rho;
paramsV.surface_ssa_geometry_source_id = ['SSA.md; Thorsos & Broschat 1995 JASA, ', ...
    'pressure-release / Dirichlet first-order perturbation-limit geometry'];
paramsV.surface_ssa_kz_branch = 'downward_positive_real';
paramsV.surface_ssa_conv_padding = 'periodic';
end

function row = local_run_row(run_index, paramsV, channel, random_reflect_f, repeat_error, direct_error)
ssa = channel.roughness_meta.ssa_stat_kernel_meta;
row = struct();
row.run_index = run_index;
row.frequency_correlation_mode = string(paramsV.surface_ssa_frequency_correlation_mode);
row.rho = paramsV.surface_ssa_frequency_correlation_rho;
row.seed = paramsV.sea_seed;
row.adjacent_random_reflect_coherence_abs = local_adjacent_coherence_abs(random_reflect_f);
row.random_reflect_energy = sum(abs(random_reflect_f(:)).^2);
row.invariant_error = max(abs(channel.H_f(:) - channel.H_direct_f(:) - channel.H_reflect_f(:)));
row.repeat_error = repeat_error;
row.direct_error_vs_metadata_only = direct_error;
row.ref_frequency_correlation_mode = string(local_get_field_or_text(ssa, 'frequency_correlation_mode', 'missing'));
row.ref_frequency_correlation_rho = local_get_field_or_nan(ssa, 'frequency_correlation_rho');
row.ref_random_spectrum_generation_rule = string(local_get_field_or_text(ssa, 'random_spectrum_generation_rule', 'missing'));
end

function value = local_adjacent_coherence_abs(x)
x = x(:);
if numel(x) < 2 || all(abs(x) == 0)
    value = NaN;
    return
end
a = x(1:end-1);
b = x(2:end);
value = abs(sum(conj(a) .* b)) / sqrt(max(sum(abs(a).^2) * sum(abs(b).^2), eps));
end

function summary_table = local_build_summary_table(run_table, mode_list)
rows = struct([]);
for ii = 1:numel(mode_list)
    subset = run_table(run_table.frequency_correlation_mode == mode_list(ii), :);
    row = struct();
    row.frequency_correlation_mode = mode_list(ii);
    row.seed_count = height(subset);
    row.adjacent_random_reflect_coherence_abs_mean = mean(subset.adjacent_random_reflect_coherence_abs);
    row.adjacent_random_reflect_coherence_abs_std = std(subset.adjacent_random_reflect_coherence_abs);
    row.random_reflect_energy_mean = mean(subset.random_reflect_energy);
    row.repeat_error_max = max(subset.repeat_error);
    row.direct_error_max = max(subset.direct_error_vs_metadata_only);
    rows = local_append_struct(rows, row);
end
summary_table = struct2table(rows);
end

function local_plot_summary(summary_table)
fig = figure('Visible', 'off');
bar(categorical(summary_table.frequency_correlation_mode), ...
    summary_table.adjacent_random_reflect_coherence_abs_mean)
hold on
errorbar(categorical(summary_table.frequency_correlation_mode), ...
    summary_table.adjacent_random_reflect_coherence_abs_mean, ...
    summary_table.adjacent_random_reflect_coherence_abs_std, '.k')
grid on
ylabel('Adjacent random-reflection coherence')
title('SSA random scatter frequency correlation')
saveas(fig, 'validate_ssa_frequency_correlation_summary.png')
close(fig)
end

function value = local_get_field_or_nan(s, field_name)
if isstruct(s) && isfield(s, field_name)
    value = s.(field_name);
else
    value = NaN;
end
end

function value = local_get_field_or_text(s, field_name, default_value)
if isstruct(s) && isfield(s, field_name)
    value = s.(field_name);
else
    value = default_value;
end
end

function rows = local_append_struct(rows, row)
if isempty(rows)
    rows = row;
else
    rows(end+1) = row; %#ok<AGROW>
end
end

function checks = local_add_check(checks, name, value, threshold, comparator)
switch comparator
    case '<='
        passed = value <= threshold;
    case '>='
        passed = value >= threshold;
    case '=='
        passed = value == threshold;
    otherwise
        error('Unsupported comparator %s', comparator);
end
row = struct('name', string(name), 'value', value, 'threshold', threshold, ...
    'comparator', string(comparator), 'passed', logical(passed));
checks = local_append_struct(checks, row);
end

function value = local_env_scalar_int(name, default_value)
raw = getenv(name);
if isempty(raw)
    value = default_value;
else
    value = round(str2double(raw));
    if ~isfinite(value) || value <= 0
        error('%s must be a positive integer.', name);
    end
end
end

function value = local_env_scalar(name, default_value)
raw = getenv(name);
if isempty(raw)
    value = default_value;
else
    value = str2double(raw);
    if ~isfinite(value)
        error('%s must be finite.', name);
    end
end
end

