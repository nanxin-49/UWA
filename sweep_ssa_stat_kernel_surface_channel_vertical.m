% Multi-Hs, multi-seed channel statistics for ssa_stat_kernel comparison.
% The script compares kirchhoff_spatial and the SSA-like statistical kernel
% on a reduced wideband grid. It stores compact frequency-response samples
% and scalar metadata, not spatial fields or full 2-D spectra.

clear
format compact

result_file = getenv('SSA_STAT_SWEEP_RESULT_FILE');
if isempty(result_file)
    result_file = 'sweep_ssa_stat_kernel_surface_channel_vertical_result.mat';
end
if exist(result_file, 'file')
    warning('Result file %s already exists and will be overwritten.', result_file);
end

sea_hs_values = [0, 0.05, 0.2, 0.5];
sea_wind_speed = 5.0;
model_specs = local_model_specs();
model_names = {model_specs.model_name};

mc_count = 8;
mc_override = str2double(getenv('SSA_STAT_MC_COUNT'));
if isfinite(mc_override) && mc_override >= 1
    mc_count = round(mc_override);
end
seed_list = 12345 + (0:(mc_count - 1));

max_conditions = Inf;
max_conditions_override = str2double(getenv('SSA_STAT_SWEEP_MAX_CONDITIONS'));
if isfinite(max_conditions_override) && max_conditions_override >= 1
    max_conditions = round(max_conditions_override);
end

base_params = local_base_params(sea_wind_speed);
condition_specs = local_build_condition_specs(base_params, model_specs, sea_hs_values);
if isfinite(max_conditions)
    condition_specs = condition_specs(1:min(max_conditions, numel(condition_specs)));
end

condition_results = struct([]);
run_rows = struct([]);
condition_rows = struct([]);
frequency_stats = struct([]);

for cc = 1:numel(condition_specs)
    [condition_result, rows_this_condition, condition_row, freq_row] = local_run_condition( ...
        condition_specs(cc), cc, numel(condition_specs), seed_list);
    condition_results = local_append_struct(condition_results, condition_result);
    run_rows = local_append_struct_array(run_rows, rows_this_condition);
    condition_rows = local_append_struct(condition_rows, condition_row);
    frequency_stats = local_append_struct(frequency_stats, freq_row);
end

run_summary_table = struct2table(run_rows);
condition_summary_table = struct2table(condition_rows);
metadata_stats = local_build_metadata_stats(condition_summary_table);
flat_check_table = local_build_flat_check_table(condition_results, model_names);

disp(condition_summary_table(:, {'model_name', 'sea_hs_target', 'mc_count', ...
    'abs_h_total_mean', 'abs_h_total_std', 'abs_h_reflect_mean', ...
    'abs_h_reflect_std', 'phase_h_total_circular_mean_rad', ...
    'E_sca_over_E_inc_mean', 'E_ref_over_E_inc_mean', ...
    'energy_scale_applied_mean', 'max_invariant_error'}))
disp(flat_check_table)

save(result_file, 'base_params', 'sea_hs_values', 'sea_wind_speed', ...
    'model_names', 'seed_list', 'condition_specs', 'condition_results', ...
    'run_summary_table', 'condition_summary_table', 'frequency_stats', ...
    'metadata_stats', 'flat_check_table');
fprintf('Saved %s\n', result_file);

function paramsV = local_base_params(sea_wind_speed)
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
paramsV.sea_wind_speed = sea_wind_speed;
paramsV.sea_hs_target = 0.5;
paramsV.sea_seed = 12345;
paramsV.surface_reflect_coeff = -1;
paramsV.surface_boundary_model = 'kirchhoff_spatial';
paramsV.surface_boundary_coupling_diagnostics = false;
paramsV.surface_boundary_redistribution_diagnostics = false;
paramsV.surface_ssa_random_scatter = true;
paramsV.surface_ssa_scatter_scale = 1.0;
paramsV.surface_ssa_seed_offset = 100000;
paramsV.surface_ssa_kernel_mode = 'pm_convolution';
paramsV.surface_ssa_geometry_source_id = '';
paramsV.surface_ssa_kz_branch = 'downward_positive_real';
paramsV.surface_ssa_conv_padding = 'periodic';
end

function specs = local_build_condition_specs(base_params, model_specs, sea_hs_values)
specs = struct([]);
idx = 0;
for im = 1:numel(model_specs)
    for ih = 1:numel(sea_hs_values)
        idx = idx + 1;
        paramsV = base_params;
        paramsV.surface_boundary_model = model_specs(im).surface_boundary_model;
        paramsV.surface_ssa_kernel_mode = model_specs(im).surface_ssa_kernel_mode;
        paramsV.surface_ssa_geometry_source_id = model_specs(im).surface_ssa_geometry_source_id;
        paramsV.sea_hs_target = sea_hs_values(ih);
        spec = struct();
        spec.condition_index = idx;
        spec.model_index = im;
        spec.model_name = model_specs(im).model_name;
        spec.surface_boundary_model = model_specs(im).surface_boundary_model;
        spec.surface_ssa_kernel_mode = model_specs(im).surface_ssa_kernel_mode;
        spec.sea_hs_target = sea_hs_values(ih);
        spec.name = sprintf('%s_Hs_%g', model_specs(im).model_name, sea_hs_values(ih));
        spec.paramsV = paramsV;
        specs = local_append_struct(specs, spec);
    end
end
end

function model_specs = local_model_specs()
source_id = ['SSA.md; Thorsos & Broschat 1995 JASA, ', ...
    'Dirichlet SSA first-order / perturbation-limit geometry'];
model_specs = struct([]);
model_specs(1).model_name = 'kirchhoff_spatial';
model_specs(1).surface_boundary_model = 'kirchhoff_spatial';
model_specs(1).surface_ssa_kernel_mode = 'pm_convolution';
model_specs(1).surface_ssa_geometry_source_id = '';
model_specs(2).model_name = 'ssa_pm_convolution';
model_specs(2).surface_boundary_model = 'ssa_stat_kernel';
model_specs(2).surface_ssa_kernel_mode = 'pm_convolution';
model_specs(2).surface_ssa_geometry_source_id = '';
model_specs(3).model_name = 'ssa1_geometry';
model_specs(3).surface_boundary_model = 'ssa_stat_kernel';
model_specs(3).surface_ssa_kernel_mode = 'ssa1_geometry';
model_specs(3).surface_ssa_geometry_source_id = source_id;
end

function [condition_result, rows, condition_row, freq_row] = local_run_condition( ...
    spec, condition_number, condition_total, seed_list)

mc_count = numel(seed_list);
H_samples = [];
H_reflect_samples = [];
H_direct_samples = [];
h_total = complex(zeros(mc_count, 1));
h_reflect = complex(zeros(mc_count, 1));
h_direct = complex(zeros(mc_count, 1));
rows = struct([]);
f_axis_ref = [];
idx_f_ref_ref = NaN;

for ii = 1:mc_count
    paramsV = spec.paramsV;
    paramsV.sea_seed = seed_list(ii);
    fprintf('SSA stat sweep condition %d/%d (%s), seed %d/%d, sea_seed=%d\n', ...
        condition_number, condition_total, spec.name, ii, mc_count, seed_list(ii));
    channel = CARPE3D_vertical(paramsV);

    invariant_error = max(abs(channel.H_f(:) - channel.H_direct_f(:) - channel.H_reflect_f(:)));
    if invariant_error > 1e-10
        error('sweep_ssa_stat_kernel_surface_channel_vertical:Invariant', ...
            'H_f invariant failed for %s seed %d: %.3e', ...
            spec.name, seed_list(ii), invariant_error);
    end

    if ii == 1
        f_axis_ref = channel.f_axis(:);
        idx_f_ref_ref = channel.idx_f_ref;
        n_freq = numel(f_axis_ref);
        H_samples = complex(zeros(n_freq, mc_count));
        H_reflect_samples = complex(zeros(n_freq, mc_count));
        H_direct_samples = complex(zeros(n_freq, mc_count));
    else
        if numel(channel.f_axis) ~= numel(f_axis_ref) || any(abs(channel.f_axis(:) - f_axis_ref) > 1e-9)
            error('Frequency axis changed for %s seed %d.', spec.name, seed_list(ii));
        end
        if channel.idx_f_ref ~= idx_f_ref_ref
            error('idx_f_ref changed for %s seed %d.', spec.name, seed_list(ii));
        end
    end

    H_samples(:, ii) = channel.H_f(:);
    H_reflect_samples(:, ii) = channel.H_reflect_f(:);
    H_direct_samples(:, ii) = channel.H_direct_f(:);
    h_total(ii) = channel.h_total;
    h_reflect(ii) = channel.h_reflect;
    h_direct(ii) = channel.h_direct;

    row = local_build_run_row(spec, seed_list(ii), ii, channel, invariant_error);
    rows = local_append_struct(rows, row);
end

freq_row = local_build_frequency_stats(spec, f_axis_ref, idx_f_ref_ref, H_samples, H_reflect_samples);
condition_result = struct();
condition_result.condition_index = spec.condition_index;
condition_result.model_name = string(spec.model_name);
condition_result.sea_hs_target = spec.sea_hs_target;
condition_result.seed_list = seed_list;
condition_result.f_axis = f_axis_ref;
condition_result.idx_f_ref = idx_f_ref_ref;
condition_result.H_f_samples = H_samples;
condition_result.H_reflect_f_samples = H_reflect_samples;
condition_result.H_direct_f_samples = H_direct_samples;
condition_result.h_total_samples = h_total;
condition_result.h_reflect_samples = h_reflect;
condition_result.h_direct_samples = h_direct;
condition_result.run_rows = rows;

condition_row = local_build_condition_row(spec, rows, h_total, h_reflect, h_direct, mc_count, idx_f_ref_ref, f_axis_ref);
end

function row = local_build_run_row(spec, sea_seed, seed_index, channel, invariant_error)
ssa = channel.roughness_meta.ssa_stat_kernel_meta;
row = struct();
row.condition_index = spec.condition_index;
row.model_index = spec.model_index;
row.model_name = string(spec.model_name);
row.seed_index = seed_index;
row.sea_seed = sea_seed;
row.sea_hs_target = spec.sea_hs_target;
row.sea_wind_speed = spec.paramsV.sea_wind_speed;
row.abs_h_total = abs(channel.h_total);
row.abs_h_reflect = abs(channel.h_reflect);
row.abs_h_direct = abs(channel.h_direct);
row.phase_h_total_rad = angle(channel.h_total);
row.phase_h_reflect_rad = angle(channel.h_reflect);
row.max_abs_H_f = max(abs(channel.H_f(:)));
row.max_abs_H_reflect_f = max(abs(channel.H_reflect_f(:)));
row.reflect_frequency_energy_ratio = sum(abs(channel.H_reflect_f(:)).^2) / ...
    max(sum(abs(channel.H_direct_f(:)).^2), eps);
row.invariant_error = invariant_error;
row.h_total_ref_consistency_error = abs(channel.h_total - channel.H_f(channel.idx_f_ref));
row.h_reflect_ref_consistency_error = abs(channel.h_reflect - channel.H_reflect_f(channel.idx_f_ref));
row.ssa_enabled = logical(ssa.enabled);
row.ssa_surface_realization_generated = logical(ssa.surface_realization_generated);
row.ssa_random_scatter_enabled = logical(ssa.random_scatter_enabled);
row.ssa_kernel_mode = string(ssa.kernel_mode);
row.ssa_geometry_source_id = string(ssa.geometry_source_id);
row.sigma_eta_m = ssa.sigma_eta_m;
row.Hs_target_m = ssa.Hs_target_m;
row.R_coh_abs = abs(ssa.R_coh);
row.R_coh_phase_rad = angle(ssa.R_coh);
row.P_sca_sum = ssa.P_sca.sum;
row.E_inc = ssa.E_inc;
row.E_coh = ssa.E_coh;
row.E_sca_raw = ssa.E_sca_raw;
row.E_sca_limited = local_get_field_or_nan(ssa, 'E_sca_limited');
row.E_sca = ssa.E_sca;
row.E_ref = local_get_field_or_nan(ssa, 'E_ref');
row.E_sca_over_E_inc = row.E_sca / max(row.E_inc, eps);
row.E_ref_over_E_inc = row.E_ref / max(row.E_inc, eps);
row.energy_scale_applied = ssa.energy_scale_applied;
row.energy_limit_applied = local_get_field_or_nan(ssa, 'energy_limit_applied');
row.energy_conservation_error = local_get_field_or_nan(ssa, 'energy_conservation_error');
row.propagating_bin_fraction = local_get_field_or_nan(ssa, 'propagating_bin_fraction');
row.seed_ssa = ssa.seed_ssa;
end

function row = local_build_condition_row(spec, rows, h_total, h_reflect, h_direct, mc_count, idx_f_ref, f_axis)
row = struct();
row.condition_index = spec.condition_index;
row.model_index = spec.model_index;
row.model_name = string(spec.model_name);
row.sea_hs_target = spec.sea_hs_target;
row.sea_wind_speed = spec.paramsV.sea_wind_speed;
row.mc_count = mc_count;
row.idx_f_ref = idx_f_ref;
row.f_ref_hz = f_axis(idx_f_ref);
row.abs_h_total_mean = mean(abs(h_total));
row.abs_h_total_std = std(abs(h_total));
row.abs_h_reflect_mean = mean(abs(h_reflect));
row.abs_h_reflect_std = std(abs(h_reflect));
row.abs_h_direct_mean = mean(abs(h_direct));
row.abs_h_direct_std = std(abs(h_direct));
row.phase_h_total_circular_mean_rad = local_circular_mean(angle(h_total));
row.phase_h_total_circular_std_rad = local_circular_std(angle(h_total));
row.reflect_frequency_energy_ratio_mean = mean([rows.reflect_frequency_energy_ratio]);
row.reflect_frequency_energy_ratio_std = std([rows.reflect_frequency_energy_ratio]);
row.E_sca_over_E_inc_mean = local_nanmean([rows.E_sca_over_E_inc]);
row.E_sca_over_E_inc_std = local_nanstd([rows.E_sca_over_E_inc]);
row.E_ref_over_E_inc_mean = local_nanmean([rows.E_ref_over_E_inc]);
row.E_ref_over_E_inc_std = local_nanstd([rows.E_ref_over_E_inc]);
row.energy_scale_applied_mean = local_nanmean([rows.energy_scale_applied]);
row.energy_scale_applied_std = local_nanstd([rows.energy_scale_applied]);
row.energy_scale_applied_fraction = local_fraction_less_than_one([rows.energy_scale_applied]);
row.energy_limit_applied_fraction = local_nanmean([rows.energy_limit_applied]);
row.energy_conservation_error_max = local_nanmax([rows.energy_conservation_error]);
row.propagating_bin_fraction_mean = local_nanmean([rows.propagating_bin_fraction]);
row.sigma_eta_m_mean = local_nanmean([rows.sigma_eta_m]);
row.R_coh_abs_mean = local_nanmean([rows.R_coh_abs]);
row.R_coh_abs_std = local_nanstd([rows.R_coh_abs]);
row.P_sca_sum_mean = local_nanmean([rows.P_sca_sum]);
row.P_sca_sum_std = local_nanstd([rows.P_sca_sum]);
row.max_invariant_error = max([rows.invariant_error]);
row.max_h_total_ref_consistency_error = max([rows.h_total_ref_consistency_error]);
row.max_h_reflect_ref_consistency_error = max([rows.h_reflect_ref_consistency_error]);
end

function freq = local_build_frequency_stats(spec, f_axis, idx_f_ref, H_samples, H_reflect_samples)
freq = struct();
freq.condition_index = spec.condition_index;
freq.model_name = string(spec.model_name);
freq.sea_hs_target = spec.sea_hs_target;
freq.f_axis = f_axis(:);
freq.idx_f_ref = idx_f_ref;
freq.abs_H_f_mean = mean(abs(H_samples), 2);
freq.abs_H_f_std = std(abs(H_samples), 0, 2);
freq.abs_H_reflect_f_mean = mean(abs(H_reflect_samples), 2);
freq.abs_H_reflect_f_std = std(abs(H_reflect_samples), 0, 2);
freq.H_f_mean = mean(H_samples, 2);
freq.H_reflect_f_mean = mean(H_reflect_samples, 2);
freq.phase_H_f_circular_mean = local_circular_mean_dim(angle(H_samples), 2);
freq.phase_H_reflect_f_circular_mean = local_circular_mean_dim(angle(H_reflect_samples), 2);
end

function metadata_stats = local_build_metadata_stats(condition_summary_table)
metadata_stats = struct();
metadata_stats.condition_summary_table = condition_summary_table(:, { ...
    'model_name', 'sea_hs_target', 'sigma_eta_m_mean', 'R_coh_abs_mean', ...
    'P_sca_sum_mean', 'E_sca_over_E_inc_mean', 'E_sca_over_E_inc_std', ...
    'E_ref_over_E_inc_mean', 'E_ref_over_E_inc_std', ...
    'energy_scale_applied_mean', 'energy_scale_applied_fraction', ...
    'energy_limit_applied_fraction', 'energy_conservation_error_max', ...
    'propagating_bin_fraction_mean'});
metadata_stats.description = ['Kirchhoff rows have NaN SSA-only metadata by design. ', ...
    'ssa_stat_kernel rows summarize the energy-limited statistical kernel.'];
end

function flat_table = local_build_flat_check_table(condition_results, model_names)
flat_table = table();
idx_k = local_find_condition(condition_results, 'kirchhoff_spatial', 0);
if isempty(idx_k)
    return
end
rows = struct([]);
for im = 1:numel(model_names)
    if strcmp(model_names{im}, 'kirchhoff_spatial')
        continue
    end
    idx_s = local_find_condition(condition_results, model_names{im}, 0);
    if isempty(idx_s)
        continue
    end
    H_k = condition_results(idx_k).H_f_samples;
    H_s = condition_results(idx_s).H_f_samples;
    H_ref_k = condition_results(idx_k).H_reflect_f_samples;
    H_ref_s = condition_results(idx_s).H_reflect_f_samples;
    row = struct();
    row.compare_name = string(['Hs0_kirchhoff_vs_' model_names{im}]);
    row.max_abs_H_f_diff = max(abs(H_k(:) - H_s(:)));
    row.max_abs_H_reflect_f_diff = max(abs(H_ref_k(:) - H_ref_s(:)));
    row.tolerance = 1e-10;
    row.passed = row.max_abs_H_f_diff <= row.tolerance && ...
        row.max_abs_H_reflect_f_diff <= row.tolerance;
    rows = local_append_struct(rows, row);
end
if ~isempty(rows)
    flat_table = struct2table(rows);
end
end

function idx = local_find_condition(condition_results, model_name, hs_value)
idx = [];
for ii = 1:numel(condition_results)
    if strcmp(char(condition_results(ii).model_name), model_name) && ...
            abs(condition_results(ii).sea_hs_target - hs_value) <= 1e-12
        idx = ii;
        return
    end
end
end

function out = local_get_field_or_nan(s, field_name)
if isfield(s, field_name)
    out = s.(field_name);
else
    out = NaN;
end
end

function mu = local_circular_mean(theta)
theta = theta(:);
theta = theta(isfinite(theta));
if isempty(theta)
    mu = NaN;
else
    mu = angle(mean(exp(1i*theta)));
end
end

function sd = local_circular_std(theta)
theta = theta(:);
theta = theta(isfinite(theta));
if isempty(theta)
    sd = NaN;
else
    r = abs(mean(exp(1i*theta)));
    sd = sqrt(max(-2*log(max(r, eps)), 0));
end
end

function mu = local_circular_mean_dim(theta, dim)
z = mean(exp(1i*theta), dim);
mu = angle(z);
end

function y = local_nanmean(x)
x = x(:);
x = x(isfinite(x));
if isempty(x)
    y = NaN;
else
    y = mean(x);
end
end

function y = local_nanstd(x)
x = x(:);
x = x(isfinite(x));
if numel(x) <= 1
    y = 0;
    if isempty(x)
        y = NaN;
    end
else
    y = std(x);
end
end

function y = local_nanmax(x)
x = x(:);
x = x(isfinite(x));
if isempty(x)
    y = NaN;
else
    y = max(x);
end
end

function y = local_fraction_less_than_one(x)
x = x(:);
x = x(isfinite(x));
if isempty(x)
    y = NaN;
else
    y = mean(x < 1);
end
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
