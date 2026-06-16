% Reduced-grid Hall1D calibration with decoupled bubble and PM wind forcing.

clear
format compact

result_file = 'calibrate_hall1d_bubble_vertical_result.mat';
figure_prefix = 'calibrate_hall1d_bubble_vertical_';
if exist(result_file, 'file')
    warning('Result file %s already exists and will be overwritten.', result_file);
end

base_params = local_base_params();

baseline = local_run_case(local_set(base_params, 'enable_bubbles', false), ...
    'baseline', 'no_bubble');
validation_table = local_run_validation(base_params, baseline);

case_specs = struct([]);
bubble_wind_values = [3, 5, 8, 12, 15];
for ii = 1:numel(bubble_wind_values)
    p = base_params;
    p.bubble_wind_speed = bubble_wind_values(ii);
    p.bubble_strength_scale = 1;
    case_specs = local_add_case(case_specs, 'bubble_wind_speed', ...
        sprintf('bubble_wind_speed_%g', bubble_wind_values(ii)), ...
        p, bubble_wind_values(ii));
end

strength_values = [1, 1e2, 1e3, 1e4, 1e5, 1e6];
for ii = 1:numel(strength_values)
    p = base_params;
    p.bubble_wind_speed = 8;
    p.bubble_strength_scale = strength_values(ii);
    case_specs = local_add_case(case_specs, 'strength_scale', ...
        sprintf('strength_scale_%g_at_bubble_wind_8', strength_values(ii)), ...
        p, strength_values(ii));
end

results = struct([]);
summary = struct([]);
for ii = 1:numel(case_specs)
    spec = case_specs(ii);
    case_result = local_run_case(spec.paramsV, spec.group, spec.name);
    [case_result, summary_row] = local_attach_metrics(case_result, baseline, spec.value);
    if ii == 1
        results = case_result;
        summary = summary_row;
    else
        results(ii) = case_result; %#ok<SAGROW>
        summary(ii) = summary_row; %#ok<SAGROW>
    end
end

summary_table = struct2table(summary);
target_table = local_target_matching(summary_table, [0.5, 1, 3]);

disp('Validation table:')
disp(validation_table)
disp('Calibration summary table:')
disp(summary_table(:, {'group', 'sweep_value', 'sea_wind_speed', ...
    'bubble_wind_speed', 'bubble_strength_scale', 'max_delta_TL_dB', ...
    'mean_delta_TL_dB', 'max_abs_phase_diff_rad', 'max_alpha_bub', ...
    'max_beta', 'invariant_error'}))
disp('Target matching table:')
disp(target_table)

local_plot_sweep(summary_table, 'bubble_wind_speed', 'bubble wind speed (m/s)', ...
    false, 'max_delta_TL_dB', 'max \DeltaTL (dB)', ...
    [figure_prefix 'max_delta_TL_vs_bubble_wind_speed.png']);
local_plot_sweep(summary_table, 'bubble_wind_speed', 'bubble wind speed (m/s)', ...
    false, 'max_alpha_bub', 'max \alpha_{bub} (Np/m)', ...
    [figure_prefix 'max_alpha_vs_bubble_wind_speed.png']);
local_plot_sweep(summary_table, 'strength_scale', 'bubble strength scale', ...
    true, 'max_delta_TL_dB', 'max \DeltaTL (dB)', ...
    [figure_prefix 'max_delta_TL_vs_strength_scale.png']);
local_plot_sweep(summary_table, 'strength_scale', 'bubble strength scale', ...
    true, 'max_alpha_bub', 'max \alpha_{bub} (Np/m)', ...
    [figure_prefix 'max_alpha_vs_strength_scale.png']);
local_plot_targets(summary_table, target_table, figure_prefix);

save(result_file, 'base_params', 'baseline', 'validation_table', ...
    'case_specs', 'results', 'summary_table', 'target_table');

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
paramsV.bubble_wind_speed = [];
end

function paramsV = local_set(paramsV, field_name, value)
paramsV.(field_name) = value;
end

function case_specs = local_add_case(case_specs, group, name, paramsV, value)
case_specs(end+1).group = group; %#ok<AGROW>
case_specs(end).name = name;
case_specs(end).paramsV = paramsV;
case_specs(end).value = value;
end

function validation_table = local_run_validation(base_params, baseline)
disabled_override = local_run_case( ...
    local_set(local_set(base_params, 'enable_bubbles', false), 'bubble_wind_speed', 15), ...
    'validation', 'disabled_with_bubble_wind');

off_params = base_params;
off_params.enable_bubbles = true;
off_params.bubble_model = 'off';
off_params.bubble_spatial_mode = 'none';
off_params.bubble_wind_speed = 15;
off_case = local_run_case(off_params, 'validation', 'off_with_bubble_wind');

level0_params = base_params;
level0_params.bubble_model = 'level0_empirical';
level0_params.bubble_spatial_mode = '1d';
level0_params.bubble_alpha0_np_per_m = 0.02;
level0_params.bubble_layer_decay_m = 20;
level0_params.bubble_delta_c0_mps = 0;
level0_params.bubble_wind_speed = [];
level0_ref = local_run_case(level0_params, 'validation', 'level0_no_bubble_wind');
level0_params.bubble_wind_speed = 15;
level0_override = local_run_case(level0_params, 'validation', 'level0_with_bubble_wind');

hall_legacy = local_run_case(base_params, 'validation', 'hall1d_legacy_wind');
hall_explicit = local_run_case(local_set(base_params, 'bubble_wind_speed', base_params.sea_wind_speed), ...
    'validation', 'hall1d_explicit_bubble_wind_5');

rows = [ ...
    local_validation_row('disabled override unchanged', disabled_override, baseline); ...
    local_validation_row('off mode unchanged', off_case, baseline); ...
    local_validation_row('level0 ignores bubble_wind_speed', level0_override, level0_ref); ...
    local_validation_row('hall1d legacy equals explicit sea wind', hall_explicit, hall_legacy)];
validation_table = struct2table(rows);
end

function row = local_validation_row(name, candidate, reference)
row.check = {name};
row.rel_H_f = local_rel_error(candidate.H_f, reference.H_f);
row.rel_H_direct_f = local_rel_error(candidate.H_direct_f, reference.H_direct_f);
row.rel_H_reflect_f = local_rel_error(candidate.H_reflect_f, reference.H_reflect_f);
row.invariant_error = candidate.invariant_error;
row.pass = row.rel_H_f < 1e-10 && row.rel_H_direct_f < 1e-10 && ...
    row.rel_H_reflect_f < 1e-10 && row.invariant_error < 1e-10;
end

function result = local_run_case(paramsV, group, name)
fprintf('Running %s / %s\n', group, name);
channel = vertical_channel_model(paramsV);
invariant_error = norm(channel.H_f(:) - channel.H_direct_f(:) - channel.H_reflect_f(:));
if invariant_error > 1e-10
    error('calibrate_hall1d_bubble_vertical:Invariant', ...
        'H_f invariant failed for %s: %.3e', name, invariant_error);
end
result.group = group;
result.scenario_name = name;
result.paramsV = paramsV;
result.channel = channel;
result.f_axis = channel.f_axis;
result.H_f = channel.H_f;
result.H_direct_f = channel.H_direct_f;
result.H_reflect_f = channel.H_reflect_f;
result.h_total = channel.h_total;
result.bubble_meta = channel.bubble_meta;
result.invariant_error = invariant_error;
end

function [result, row] = local_attach_metrics(result, baseline, sweep_value)
H_case = result.H_f(:);
H_base = baseline.H_f(:);
if numel(H_case) ~= numel(H_base) || any(abs(result.f_axis(:) - baseline.f_axis(:)) > 1e-9)
    error('Frequency axis mismatch for case %s.', result.scenario_name);
end

delta_TL_dB = -20*log10(abs(H_case) ./ max(abs(H_base), eps));
phase_diff_rad = unwrap(angle(H_case)) - unwrap(angle(H_base));

result.delta_TL_dB = delta_TL_dB;
result.phase_diff_rad = phase_diff_rad;
result.max_delta_TL_dB = max(delta_TL_dB);
result.mean_delta_TL_dB = mean(delta_TL_dB);
result.max_abs_phase_diff_rad = max(abs(phase_diff_rad));
result.max_alpha_bub = local_get_stat(result.bubble_meta, 'alpha_bub_stats', 'max');
result.mean_alpha_bub = local_get_stat(result.bubble_meta, 'alpha_bub_stats', 'mean');
result.max_beta = local_get_stat(result.bubble_meta, 'beta_stats', 'max');

row.group = result.group;
row.scenario_name = result.scenario_name;
row.sea_wind_speed = result.paramsV.sea_wind_speed;
row.bubble_wind_speed = local_result_bubble_wind_speed(result);
row.bubble_strength_scale = result.paramsV.bubble_strength_scale;
row.sweep_value = sweep_value;
row.h_total = result.h_total;
row.max_delta_TL_dB = result.max_delta_TL_dB;
row.mean_delta_TL_dB = result.mean_delta_TL_dB;
row.max_abs_phase_diff_rad = result.max_abs_phase_diff_rad;
row.max_alpha_bub = result.max_alpha_bub;
row.mean_alpha_bub = result.mean_alpha_bub;
row.max_beta = result.max_beta;
row.invariant_error = result.invariant_error;
end

function U10 = local_result_bubble_wind_speed(result)
U10 = result.paramsV.sea_wind_speed;
if isfield(result.bubble_meta, 'bubble_wind_speed') && ~isempty(result.bubble_meta.bubble_wind_speed)
    U10 = result.bubble_meta.bubble_wind_speed;
elseif isfield(result.paramsV, 'bubble_wind_speed') && ~isempty(result.paramsV.bubble_wind_speed)
    U10 = result.paramsV.bubble_wind_speed;
end
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

function err = local_rel_error(a, b)
err = norm(a(:) - b(:)) / max(norm(b(:)), eps);
end

function target_table = local_target_matching(summary_table, targets)
mask = strcmp(summary_table.group, 'strength_scale');
strength_rows = summary_table(mask, :);
rows = struct([]);
for ii = 1:numel(targets)
    [~, idx] = min(abs(strength_rows.max_delta_TL_dB - targets(ii)));
    rows(ii).target_delta_TL_dB = targets(ii); %#ok<AGROW>
    rows(ii).matched_strength_scale = strength_rows.bubble_strength_scale(idx);
    rows(ii).matched_bubble_wind_speed = strength_rows.bubble_wind_speed(idx);
    rows(ii).matched_max_delta_TL_dB = strength_rows.max_delta_TL_dB(idx);
    rows(ii).matched_mean_delta_TL_dB = strength_rows.mean_delta_TL_dB(idx);
    rows(ii).matched_max_alpha_bub = strength_rows.max_alpha_bub(idx);
    rows(ii).matched_max_beta = strength_rows.max_beta(idx);
    rows(ii).abs_error_dB = abs(strength_rows.max_delta_TL_dB(idx) - targets(ii));
end
target_table = struct2table(rows);
end

function local_plot_sweep(summary_table, group, x_label, use_log_x, y_field, y_label, file_name)
mask = strcmp(summary_table.group, group);
x = summary_table.sweep_value(mask);
y = summary_table.(y_field)(mask);
[x, order] = sort(x);
y = y(order);

fig = figure('Visible', 'off');
if use_log_x
    semilogx(x, y, '-o', 'LineWidth', 1.3)
else
    plot(x, y, '-o', 'LineWidth', 1.3)
end
grid on
xlabel(x_label)
ylabel(y_label)
title(strrep(sprintf('%s vs %s', y_field, group), '_', '\_'))
print(fig, '-dpng', '-r200', file_name)
close(fig)
end

function local_plot_targets(summary_table, target_table, figure_prefix)
mask = strcmp(summary_table.group, 'strength_scale');
x = summary_table.bubble_strength_scale(mask);
y = summary_table.max_delta_TL_dB(mask);
[x, order] = sort(x);
y = y(order);

fig = figure('Visible', 'off');
semilogx(x, y, '-o', 'LineWidth', 1.3)
hold on
for ii = 1:height(target_table)
    yline(target_table.target_delta_TL_dB(ii), '--', ...
        sprintf('target %.1f dB', target_table.target_delta_TL_dB(ii)));
end
grid on
xlabel('bubble strength scale')
ylabel('max \DeltaTL (dB)')
title('Target matching by Hall1D strength scale')
print(fig, '-dpng', '-r200', [figure_prefix 'target_matching.png'])
close(fig)
end
