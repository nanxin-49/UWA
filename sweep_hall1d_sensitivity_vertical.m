% Reduced-grid Hall1D bubble-model sensitivity diagnostics.

clear
format compact

result_file = 'sweep_hall1d_sensitivity_vertical_result.mat';
figure_prefix = 'sweep_hall1d_sensitivity_vertical_';
if exist(result_file, 'file')
    warning('Result file %s already exists and will be overwritten.', result_file);
end

base_params = local_base_params();

baseline = local_run_case(base_params, 'baseline', 'no_bubble', ...
    'enable_bubbles', false);

case_specs = struct([]);
case_specs = local_add_cases(case_specs, 'strength_scale', ...
    'bubble_strength_scale', [1, 1e2, 1e4, 1e6], ...
    base_params, @(p, v) local_set(p, 'bubble_strength_scale', v));
case_specs = local_add_cases(case_specs, 'wind_speed', ...
    'sea_wind_speed', [3, 5, 8, 12, 15], ...
    base_params, @(p, v) local_set(p, 'sea_wind_speed', v));
case_specs = local_add_cases(case_specs, 'damping_const', ...
    'bubble_delta_const', [0.03, 0.1, 0.3, 1.0], ...
    local_set(base_params, 'bubble_strength_scale', 1e4), ...
    @(p, v) local_set(p, 'bubble_delta_const', v));

radius_grids = { ...
    logspace(-5, -3, 80), ...
    logspace(-5, log10(2e-3), 100), ...
    logspace(-5, log10(3e-3), 120)};
radius_values = cellfun(@max, radius_grids);
for ii = 1:numel(radius_grids)
    p = local_set(base_params, 'bubble_strength_scale', 1e4);
    p.bubble_radius_grid_m = radius_grids{ii};
    case_specs(end+1).group = 'radius_grid_max'; %#ok<SAGROW>
    case_specs(end).name = sprintf('radius_grid_max_%.4g_m', radius_values(ii));
    case_specs(end).parameter = 'radius_grid_max_m';
    case_specs(end).value = radius_values(ii);
    case_specs(end).paramsV = p;
    case_specs(end).baseline_key = 'zrx_3';
end

z_rx_values = [1, 3, 5];
for ii = 1:numel(z_rx_values)
    p = local_set(base_params, 'bubble_strength_scale', 1e4);
    p.z_rx = z_rx_values(ii);
    case_specs(end+1).group = 'receiver_depth'; %#ok<SAGROW>
    case_specs(end).name = sprintf('z_rx_%g_m', z_rx_values(ii));
    case_specs(end).parameter = 'z_rx_m';
    case_specs(end).value = z_rx_values(ii);
    case_specs(end).paramsV = p;
    case_specs(end).baseline_key = sprintf('zrx_%g', z_rx_values(ii));
end

baselines = struct();
baselines.zrx_3 = baseline;
for z_rx = [1, 5]
    p0 = local_set(base_params, 'z_rx', z_rx);
    key = sprintf('zrx_%g', z_rx);
    baselines.(key) = local_run_case(p0, 'baseline', ...
        sprintf('no_bubble_z_rx_%g_m', z_rx), 'enable_bubbles', false);
end

results = struct([]);
summary = struct([]);
for ii = 1:numel(case_specs)
    spec = case_specs(ii);
    case_result = local_run_case(spec.paramsV, spec.group, spec.name);
    baseline_match = baselines.(spec.baseline_key);
    [case_result, summary_row] = local_attach_metrics( ...
        case_result, baseline_match, spec.parameter, spec.value);
    if ii == 1
        results = case_result;
        summary = summary_row;
    else
        results(ii) = case_result; %#ok<SAGROW>
        summary(ii) = summary_row; %#ok<SAGROW>
    end
end

summary_table = struct2table(summary);
disp(summary_table(:, {'group', 'sweep_value', 'max_delta_TL_dB', ...
    'mean_delta_TL_dB', 'max_abs_phase_diff_rad', 'max_alpha_bub', ...
    'max_beta', 'invariant_error'}))

local_plot_sweep(summary_table, 'strength_scale', 'bubble_strength_scale', true, figure_prefix);
local_plot_sweep(summary_table, 'wind_speed', 'sea_wind_speed', false, figure_prefix);
local_plot_sweep(summary_table, 'damping_const', 'bubble_delta_const', true, figure_prefix);
local_plot_sweep(summary_table, 'radius_grid_max', 'radius_grid_max_m', true, figure_prefix);
local_plot_sweep(summary_table, 'receiver_depth', 'z_rx_m', false, figure_prefix);
local_plot_alpha_summary(summary_table, figure_prefix);

save(result_file, 'base_params', 'baseline', 'baselines', ...
    'case_specs', 'results', 'summary_table');

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
paramsV.enable_bubbles = true;
paramsV.bubble_model = 'hall1d';
paramsV.bubble_spatial_mode = '1d';
paramsV.bubble_strength_scale = 1;
end

function case_specs = local_add_cases(case_specs, group, parameter, values, base_params, make_params)
for ii = 1:numel(values)
    value = values(ii);
    p = make_params(base_params, value);
    case_specs(end+1).group = group; %#ok<AGROW>
    case_specs(end).name = sprintf('%s_%g', parameter, value);
    case_specs(end).parameter = parameter;
    case_specs(end).value = value;
    case_specs(end).paramsV = p;
    case_specs(end).baseline_key = 'zrx_3';
end
end

function paramsV = local_set(paramsV, field_name, value)
paramsV.(field_name) = value;
end

function result = local_run_case(paramsV, group, name, varargin)
if mod(numel(varargin), 2) ~= 0
    error('local_run_case expects name/value overrides.');
end
for ii = 1:2:numel(varargin)
    paramsV.(varargin{ii}) = varargin{ii+1};
end

fprintf('Running %s / %s\n', group, name);
channel = vertical_channel_model(paramsV);
invariant_error = norm(channel.H_f(:) - channel.H_direct_f(:) - channel.H_reflect_f(:));
if invariant_error > 1e-10
    error('sweep_hall1d_sensitivity_vertical:Invariant', ...
        'H_f invariant failed for %s: %.3e', name, invariant_error);
end

result.group = group;
result.name = name;
result.paramsV = paramsV;
result.channel = channel;
result.f_axis = channel.f_axis;
result.H_f = channel.H_f;
result.H_direct_f = channel.H_direct_f;
result.H_reflect_f = channel.H_reflect_f;
result.h_total = channel.h_total;
result.abs_h_total = abs(channel.h_total);
result.phase_h_total = angle(channel.h_total);
result.bubble_meta = channel.bubble_meta;
result.invariant_error = invariant_error;
end

function [result, row] = local_attach_metrics(result, baseline, parameter, value)
H_case = result.H_f(:);
H_base = baseline.H_f(:);
if numel(H_case) ~= numel(H_base) || any(abs(result.f_axis(:) - baseline.f_axis(:)) > 1e-9)
    error('Frequency axis mismatch for case %s.', result.name);
end

delta_TL_dB = -20*log10(abs(H_case) ./ max(abs(H_base), eps));
phase_diff_rad = unwrap(angle(H_case)) - unwrap(angle(H_base));
result.delta_TL_dB = delta_TL_dB;
result.phase_diff_rad = phase_diff_rad;
result.scenario_group = result.group;
result.scenario_name = result.name;
result.sweep_parameter = parameter;
result.sweep_value = value;
result.max_delta_TL_dB = max(delta_TL_dB);
result.mean_delta_TL_dB = mean(delta_TL_dB);
result.max_abs_phase_diff_rad = max(abs(phase_diff_rad));
result.max_alpha_bub = local_get_stat(result.bubble_meta, 'alpha_bub_stats', 'max');
result.mean_alpha_bub = local_get_stat(result.bubble_meta, 'alpha_bub_stats', 'mean');
result.max_beta = local_get_stat(result.bubble_meta, 'beta_stats', 'max');
result.max_resonance_radius_m = local_get_stat(result.bubble_meta, 'resonance_radius_stats', 'max');

row.scenario_group = result.group;
row.group = result.group;
row.scenario_name = result.name;
row.sweep_parameter = parameter;
row.sweep_value = value;
row.abs_h_total = result.abs_h_total;
row.phase_h_total = result.phase_h_total;
row.max_delta_TL_dB = result.max_delta_TL_dB;
row.mean_delta_TL_dB = result.mean_delta_TL_dB;
row.max_abs_phase_diff_rad = result.max_abs_phase_diff_rad;
row.max_alpha_bub = result.max_alpha_bub;
row.mean_alpha_bub = result.mean_alpha_bub;
row.max_beta = result.max_beta;
row.max_resonance_radius_m = result.max_resonance_radius_m;
row.invariant_error = result.invariant_error;
end

function value = local_get_stat(meta, stats_field, stat_name)
value = NaN;
if isstruct(meta) && isfield(meta, stats_field)
    stats = meta.(stats_field);
    if isstruct(stats) && isfield(stats, stat_name)
        value = stats.(stat_name);
    end
end
end

function local_plot_sweep(summary_table, group, parameter_label, use_log_x, figure_prefix)
mask = strcmp(summary_table.group, group);
x = summary_table.sweep_value(mask);
y = summary_table.max_delta_TL_dB(mask);
alpha = summary_table.max_alpha_bub(mask);
[x, order] = sort(x);
y = y(order);
alpha = alpha(order);

fig = figure('Visible', 'off');
yyaxis left
if use_log_x
    semilogx(x, y, '-o', 'LineWidth', 1.3)
else
    plot(x, y, '-o', 'LineWidth', 1.3)
end
grid on
xlabel(strrep(parameter_label, '_', '\_'))
ylabel('max \DeltaTL (dB)')

yyaxis right
if use_log_x
    semilogx(x, alpha, '--s', 'LineWidth', 1.1)
else
    plot(x, alpha, '--s', 'LineWidth', 1.1)
end
ylabel('max \alpha_{bub} (Np/m)')
title(sprintf('Hall1D sensitivity: %s', strrep(group, '_', '\_')))
legend({'max \DeltaTL', 'max \alpha_{bub}'}, 'Location', 'best')
print(fig, '-dpng', '-r200', [figure_prefix group '.png'])
close(fig)
end

function local_plot_alpha_summary(summary_table, figure_prefix)
groups = unique(summary_table.group, 'stable');
fig = figure('Visible', 'off');
tiledlayout(numel(groups), 1)
for ii = 1:numel(groups)
    nexttile
    mask = strcmp(summary_table.group, groups{ii});
    x = summary_table.sweep_value(mask);
    y = summary_table.max_alpha_bub(mask);
    [x, order] = sort(x);
    y = y(order);
    plot(x, y, '-o', 'LineWidth', 1.2)
    grid on
    ylabel('max \alpha')
    title(strrep(groups{ii}, '_', '\_'))
end
xlabel('sweep value')
print(fig, '-dpng', '-r200', [figure_prefix 'max_alpha_summary.png'])
close(fig)
end
