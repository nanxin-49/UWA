% Scatter-scale sensitivity sweep for ssa_stat_kernel kernels.
% Compares pm_convolution and ssa1_geometry over Hs, scatter scale, and seed.

clear
format compact

result_file = getenv('SSA_SCALE_SWEEP_RESULT_FILE');
if isempty(result_file)
    result_file = 'sweep_ssa_scatter_scale_sensitivity_vertical_result.mat';
end
if exist(result_file, 'file')
    warning('Result file %s already exists and will be overwritten.', result_file);
end

sea_hs_values = [0.05, 0.2, 0.5];
scatter_scale_values = [0, 0.25, 1, 4];
model_specs = local_model_specs();
model_names = {model_specs.model_name};

mc_count = 4;
mc_override = str2double(getenv('SSA_SCALE_SWEEP_MC_COUNT'));
if isfinite(mc_override) && mc_override >= 1
    mc_count = round(mc_override);
end
seed_list = 12345 + (0:(mc_count - 1));

max_conditions = Inf;
max_conditions_override = str2double(getenv('SSA_SCALE_SWEEP_MAX_CONDITIONS'));
if isfinite(max_conditions_override) && max_conditions_override >= 1
    max_conditions = round(max_conditions_override);
end

base_params = local_base_params();
condition_specs = local_build_condition_specs(base_params, model_specs, sea_hs_values, scatter_scale_values);
if isfinite(max_conditions)
    condition_specs = condition_specs(1:min(max_conditions, numel(condition_specs)));
end

run_rows = struct([]);
condition_rows = struct([]);
condition_results = struct([]);

for cc = 1:numel(condition_specs)
    [condition_result, rows_this_condition, condition_row] = local_run_condition( ...
        condition_specs(cc), cc, numel(condition_specs), seed_list);
    condition_results = local_append_struct(condition_results, condition_result);
    run_rows = local_append_struct_array(run_rows, rows_this_condition);
    condition_rows = local_append_struct(condition_rows, condition_row);
end

run_summary_table = struct2table(run_rows);
condition_summary_table = struct2table(condition_rows);
sensitivity_meta = struct();
sensitivity_meta.script = mfilename;
sensitivity_meta.created_at = char(datetime('now'));
sensitivity_meta.notes = ['surface_ssa_scatter_scale is an engineering normalization ', ...
    'factor for sensitivity analysis, not a physical calibration parameter.'];

disp(condition_summary_table(:, {'model_name', 'sea_hs_target', 'surface_ssa_scatter_scale', ...
    'mc_count', 'abs_h_reflect_mean', 'E_sca_raw_over_E_inc_mean', ...
    'E_sca_limited_over_E_inc_mean', 'E_sca_over_E_inc_mean', ...
    'energy_limit_applied_fraction', ...
    'energy_conservation_error_max'}))

save(result_file, 'base_params', 'sea_hs_values', 'scatter_scale_values', ...
    'model_names', 'model_specs', 'seed_list', 'condition_specs', ...
    'condition_results', 'run_summary_table', 'condition_summary_table', ...
    'sensitivity_meta');
fprintf('Saved %s\n', result_file);

function paramsV = local_base_params()
paramsV = struct();
paramsV.f0 = 6000;
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
paramsV.sea_hs_target = 0.2;
paramsV.sea_seed = 12345;
paramsV.surface_reflect_coeff = -1;
paramsV.surface_boundary_model = 'ssa_stat_kernel';
paramsV.surface_boundary_coupling_diagnostics = false;
paramsV.surface_boundary_redistribution_diagnostics = false;
paramsV.surface_ssa_random_scatter = false;
paramsV.surface_ssa_scatter_scale = 1.0;
paramsV.surface_ssa_seed_offset = 100000;
paramsV.surface_ssa_kernel_mode = 'pm_convolution';
paramsV.surface_ssa_geometry_source_id = '';
paramsV.surface_ssa_kz_branch = 'downward_positive_real';
paramsV.surface_ssa_conv_padding = 'periodic';
end

function model_specs = local_model_specs()
source_id = ['SSA.md; Thorsos & Broschat 1995 JASA, ', ...
    'Dirichlet SSA first-order / perturbation-limit geometry'];
model_specs = struct([]);
model_specs(1).model_name = 'ssa_pm_convolution';
model_specs(1).surface_ssa_kernel_mode = 'pm_convolution';
model_specs(1).surface_ssa_geometry_source_id = '';
model_specs(2).model_name = 'ssa1_geometry';
model_specs(2).surface_ssa_kernel_mode = 'ssa1_geometry';
model_specs(2).surface_ssa_geometry_source_id = source_id;
end

function specs = local_build_condition_specs(base_params, model_specs, sea_hs_values, scatter_scale_values)
specs = struct([]);
condition_index = 0;
for im = 1:numel(model_specs)
    for ih = 1:numel(sea_hs_values)
        for iscale = 1:numel(scatter_scale_values)
            condition_index = condition_index + 1;
            paramsV = base_params;
            paramsV.surface_ssa_kernel_mode = model_specs(im).surface_ssa_kernel_mode;
            paramsV.surface_ssa_geometry_source_id = model_specs(im).surface_ssa_geometry_source_id;
            paramsV.sea_hs_target = sea_hs_values(ih);
            paramsV.surface_ssa_scatter_scale = scatter_scale_values(iscale);
            spec = struct();
            spec.condition_index = condition_index;
            spec.model_index = im;
            spec.model_name = model_specs(im).model_name;
            spec.sea_hs_target = sea_hs_values(ih);
            spec.surface_ssa_scatter_scale = scatter_scale_values(iscale);
            spec.paramsV = paramsV;
            spec.name = sprintf('%s_Hs_%g_scale_%g', ...
                spec.model_name, spec.sea_hs_target, spec.surface_ssa_scatter_scale);
            specs = local_append_struct(specs, spec);
        end
    end
end
end

function [condition_result, rows, condition_row] = local_run_condition(spec, condition_number, condition_total, seed_list)
rows = struct([]);
mc_count = numel(seed_list);
for ii = 1:mc_count
    paramsV = spec.paramsV;
    paramsV.sea_seed = seed_list(ii);
    fprintf('SSA scale sweep condition %d/%d (%s), seed %d/%d, sea_seed=%d\n', ...
        condition_number, condition_total, spec.name, ii, mc_count, seed_list(ii));
    channel = vertical_channel_model(paramsV);
    invariant_error = max(abs(channel.H_f(:) - channel.H_direct_f(:) - channel.H_reflect_f(:)));
    if invariant_error > 1e-10
        error('sweep_ssa_scatter_scale_sensitivity_vertical:Invariant', ...
            'H_f invariant failed for %s seed %d: %.3e', spec.name, seed_list(ii), invariant_error);
    end
    rows = local_append_struct(rows, local_build_run_row(spec, ii, seed_list(ii), channel, invariant_error));
end

condition_result = struct();
condition_result.condition_index = spec.condition_index;
condition_result.model_name = string(spec.model_name);
condition_result.sea_hs_target = spec.sea_hs_target;
condition_result.surface_ssa_scatter_scale = spec.surface_ssa_scatter_scale;
condition_result.seed_list = seed_list;
condition_result.run_rows = rows;

condition_row = local_build_condition_row(spec, rows, mc_count);
end

function row = local_build_run_row(spec, seed_index, sea_seed, channel, invariant_error)
ssa = channel.roughness_meta.ssa_stat_kernel_meta;
row = struct();
row.condition_index = spec.condition_index;
row.model_index = spec.model_index;
row.model_name = string(spec.model_name);
row.seed_index = seed_index;
row.sea_seed = sea_seed;
row.sea_hs_target = spec.sea_hs_target;
row.surface_ssa_scatter_scale = spec.surface_ssa_scatter_scale;
row.kernel_mode = string(ssa.kernel_mode);
row.geometry_source_id = string(ssa.geometry_source_id);
row.conv_padding = string(ssa.conv_padding);
row.formula_source = string(ssa.formula_source);
row.boundary_condition = string(ssa.boundary_condition);
row.random_scatter_enabled = logical(ssa.random_scatter_enabled);
row.abs_h_total = abs(channel.h_total);
row.abs_h_reflect = abs(channel.h_reflect);
row.phase_h_total_rad = angle(channel.h_total);
row.invariant_error = invariant_error;
row.R_coh_abs = abs(ssa.R_coh);
row.P_sca_raw_sum = ssa.P_sca_raw.sum;
row.P_sca_sum = ssa.P_sca.sum;
row.E_inc = ssa.E_inc;
row.E_coh = ssa.E_coh;
row.E_sca_raw = ssa.E_sca_raw;
row.E_sca_limited = ssa.E_sca_limited;
row.E_sca = ssa.E_sca;
row.E_ref = ssa.E_ref;
row.E_sca_raw_over_E_inc = ssa.E_sca_raw / max(ssa.E_inc, eps);
row.E_sca_limited_over_E_inc = ssa.E_sca_limited / max(ssa.E_inc, eps);
row.E_sca_over_E_inc = ssa.E_sca / max(ssa.E_inc, eps);
row.E_ref_over_E_inc = ssa.E_ref / max(ssa.E_inc, eps);
row.energy_scale_applied = ssa.energy_scale_applied;
row.energy_limit_applied = logical(ssa.energy_limit_applied);
row.energy_conservation_error = ssa.energy_conservation_error;
row.propagating_bin_fraction = ssa.propagating_bin_fraction;
row.seed_ssa = ssa.seed_ssa;
end

function row = local_build_condition_row(spec, rows, mc_count)
row = struct();
row.condition_index = spec.condition_index;
row.model_index = spec.model_index;
row.model_name = string(spec.model_name);
row.sea_hs_target = spec.sea_hs_target;
row.surface_ssa_scatter_scale = spec.surface_ssa_scatter_scale;
row.mc_count = mc_count;
row.abs_h_total_mean = mean([rows.abs_h_total]);
row.abs_h_total_std = std([rows.abs_h_total]);
row.abs_h_reflect_mean = mean([rows.abs_h_reflect]);
row.abs_h_reflect_std = std([rows.abs_h_reflect]);
row.R_coh_abs_mean = mean([rows.R_coh_abs]);
row.P_sca_raw_sum_mean = mean([rows.P_sca_raw_sum]);
row.P_sca_raw_sum_std = std([rows.P_sca_raw_sum]);
row.E_sca_raw_over_E_inc_mean = mean([rows.E_sca_raw_over_E_inc]);
row.E_sca_raw_over_E_inc_std = std([rows.E_sca_raw_over_E_inc]);
row.E_sca_limited_over_E_inc_mean = mean([rows.E_sca_limited_over_E_inc]);
row.E_sca_limited_over_E_inc_std = std([rows.E_sca_limited_over_E_inc]);
row.E_sca_over_E_inc_mean = mean([rows.E_sca_over_E_inc]);
row.E_sca_over_E_inc_std = std([rows.E_sca_over_E_inc]);
row.E_ref_over_E_inc_mean = mean([rows.E_ref_over_E_inc]);
row.energy_scale_applied_mean = mean([rows.energy_scale_applied]);
row.energy_limit_applied_fraction = mean([rows.energy_limit_applied]);
row.energy_conservation_error_max = max([rows.energy_conservation_error]);
row.propagating_bin_fraction_mean = mean([rows.propagating_bin_fraction]);
row.max_invariant_error = max([rows.invariant_error]);
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
