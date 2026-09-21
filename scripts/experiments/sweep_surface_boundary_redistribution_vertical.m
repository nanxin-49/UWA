run(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'bootstrap_project.m'));
% Reduced-grid diagnostics for Kirchhoff incident-weighted k-domain redistribution.
% This script varies sea-state parameters and scalar frequency, then records
% how the current incident angular spectrum is redistributed by the surface
% phase screen. The metrics are diagnostics, not scattering cross sections.

clear
format compact

result_file = 'sweep_surface_boundary_redistribution_vertical_result.mat';
figure_prefix = 'sweep_surface_boundary_redistribution_vertical_';
if exist(result_file, 'file')
    warning('Result file %s already exists and will be overwritten.', result_file);
end

base_params = local_base_params();

case_specs = struct([]);
case_specs = local_add_cases(case_specs, 'Hs', 'sea_hs_target', [0.05, 0.5, 1.0], ...
    base_params, @(p, v) local_set(p, 'sea_hs_target', v));
case_specs = local_add_cases(case_specs, 'wind', 'sea_wind_speed', [3, 5, 8, 12], ...
    base_params, @(p, v) local_set(p, 'sea_wind_speed', v));
case_specs = local_add_cases(case_specs, 'frequency', 'f0', [4000, 6000, 8000], ...
    base_params, @(p, v) local_set(p, 'f0', v));

results = struct([]);
summary = struct([]);
for ii = 1:numel(case_specs)
    [result, row] = local_run_case(case_specs(ii));
    if ii == 1
        results = result;
        summary = row;
    else
        results(ii) = result; %#ok<SAGROW>
        summary(ii) = row; %#ok<SAGROW>
    end
end

summary_table = struct2table(summary);
disp(summary_table(:, {'group', 'parameter', 'sweep_value', ...
    'incident_rms_delta_k', 'reflect_rms_delta_k', 'rms_delta_k_increase', ...
    'centroid_shift_mag', 'reflect_high_k_fraction', ...
    'rough_vs_flat_rms_delta_k_increase', 'invariant_error'}))

local_plot_group(summary_table, 'Hs', 'sea_hs_target', 'rms_delta_k_increase', ...
    'RMS broadening vs Hs', figure_prefix, 'Hs_rms_broadening');
local_plot_group(summary_table, 'wind', 'sea_wind_speed', 'rms_delta_k_increase', ...
    'RMS broadening vs wind speed', figure_prefix, 'wind_rms_broadening');
local_plot_group(summary_table, 'frequency', 'f0', 'rms_delta_k_increase', ...
    'RMS broadening vs frequency', figure_prefix, 'frequency_rms_broadening');
local_plot_group(summary_table, 'Hs', 'sea_hs_target', 'reflect_high_k_fraction', ...
    'Reflected high-k fraction vs Hs', figure_prefix, 'Hs_high_k_fraction');

save(result_file, 'base_params', 'case_specs', 'results', 'summary_table');

function paramsV = local_base_params()
paramsV = struct();
paramsV.f0 = 4000;
paramsV.enable_wideband = false;
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
paramsV.surface_boundary_redistribution_diagnostics = true;
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
end
end

function paramsV = local_set(paramsV, field_name, value)
paramsV.(field_name) = value;
end

function [result, row] = local_run_case(case_spec)
fprintf('Running redistribution diagnostic %s / %s\n', case_spec.group, case_spec.name);
channel = vertical_channel_model(case_spec.paramsV);
invariant_error = max(abs(channel.H_f(:) - channel.H_direct_f(:) - channel.H_reflect_f(:)));
if invariant_error > 1e-10
    error('sweep_surface_boundary_redistribution_vertical:Invariant', ...
        'H_f invariant failed for %s: %.3e', case_spec.name, invariant_error);
end

diag = channel.roughness_meta.boundary_redistribution_diagnostics;
if ~diag.enabled
    error('Redistribution diagnostics were not enabled for %s.', case_spec.name);
end

result.group = case_spec.group;
result.name = case_spec.name;
result.parameter = case_spec.parameter;
result.value = case_spec.value;
result.paramsV = case_spec.paramsV;
result.channel = channel;
result.redistribution = diag;
result.invariant_error = invariant_error;

row.group = string(case_spec.group);
row.name = string(case_spec.name);
row.parameter = string(case_spec.parameter);
row.sweep_value = case_spec.value;
row.f0 = case_spec.paramsV.f0;
row.sea_hs_target = case_spec.paramsV.sea_hs_target;
row.sea_wind_speed = case_spec.paramsV.sea_wind_speed;
row.incident_rms_delta_k = diag.incident_rms_delta_k_rad_per_m;
row.reflect_rms_delta_k = diag.reflect_rms_delta_k_rad_per_m;
row.rms_delta_k_increase = diag.rms_delta_k_increase_rad_per_m;
row.centroid_shift_mag = diag.centroid_shift_mag_rad_per_m;
row.reflect_high_k_fraction = diag.reflect_high_k_fraction;
row.flat_ref_rms_delta_k = diag.flat_ref_rms_delta_k_rad_per_m;
row.rough_vs_flat_rms_delta_k_increase = diag.rough_vs_flat_rms_delta_k_increase_rad_per_m;
row.rough_vs_flat_high_k_fraction_increase = diag.rough_vs_flat_high_k_fraction_increase;
row.invariant_error = invariant_error;
end

function local_plot_group(summary_table, group_name, x_name, y_name, title_text, figure_prefix, suffix)
mask = summary_table.group == string(group_name);
if ~any(mask)
    return
end
figure('Visible', 'off'); %#ok<UNRCH>
x = summary_table.sweep_value(mask);
y = summary_table.(y_name)(mask);
[x, order] = sort(x);
y = y(order);
plot(x, y, 'o-', 'LineWidth', 1.2)
grid on
xlabel(x_name, 'Interpreter', 'none')
ylabel(y_name, 'Interpreter', 'none')
title(title_text, 'Interpreter', 'none')
print(gcf, '-dpng', '-r200', [figure_prefix suffix '.png'])
close(gcf)
end

