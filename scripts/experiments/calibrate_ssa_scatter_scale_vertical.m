run(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'bootstrap_project.m'));
% Calibrate the engineering SSA scatter scale against Kirchhoff statistics.
% This script compares reduced-grid realization-based Kirchhoff statistics
% with ssa1_geometry over candidate surface_ssa_scatter_scale values.

clear
format compact

result_file = getenv('SSA_SCALE_CAL_RESULT_FILE');
if isempty(result_file)
    result_file = 'calibrate_ssa_scatter_scale_vertical_result.mat';
end

Hs_list = local_env_numeric_vector('SSA_SCALE_CAL_HS_LIST', [0.05, 0.2, 0.5]);
f_list_hz = local_env_numeric_vector('SSA_SCALE_CAL_F_LIST_HZ', [4000, 6000, 8000]);
scale_values = local_env_numeric_vector('SSA_SCALE_CAL_SCALE_VALUES', [0, 0.1, 0.25, 0.5, 1, 2, 4]);
seed_count = local_env_scalar_int('SSA_SCALE_CAL_SEED_COUNT', 4);
grid_n = local_env_scalar_int('SSA_SCALE_CAL_GRID_N', 64);
seed_list = 12345 + (0:(seed_count - 1));
sea_wind_speed = 5;

run_rows = struct([]);
summary_rows = struct([]);
run_index = 0;

for ff = 1:numel(f_list_hz)
    for hh = 1:numel(Hs_list)
        for ss = 1:numel(seed_list)
            paramsV = local_base_params(grid_n, sea_wind_speed);
            paramsV.f0 = f_list_hz(ff);
            paramsV.f_ref_hz = paramsV.f0;
            paramsV.sea_hs_target = Hs_list(hh);
            paramsV.sea_seed = seed_list(ss);
            paramsV.surface_boundary_model = 'kirchhoff_spatial';

            fprintf('Scale calibration Kirchhoff: f=%g Hz, Hs=%g m, seed=%d.\n', ...
                paramsV.f0, paramsV.sea_hs_target, paramsV.sea_seed);
            channel = vertical_channel_model(paramsV);
            run_index = run_index + 1;
            run_rows = local_append_struct(run_rows, ...
                local_run_row(run_index, "kirchhoff_spatial", paramsV, channel));
        end

        for iscale = 1:numel(scale_values)
            for ss = 1:numel(seed_list)
                paramsV = local_base_params(grid_n, sea_wind_speed);
                paramsV.f0 = f_list_hz(ff);
                paramsV.f_ref_hz = paramsV.f0;
                paramsV.sea_hs_target = Hs_list(hh);
                paramsV.sea_seed = seed_list(ss);
                paramsV.surface_boundary_model = 'ssa_stat_kernel';
                paramsV.surface_ssa_kernel_mode = 'ssa1_geometry';
                paramsV.surface_ssa_random_scatter = true;
                paramsV.surface_ssa_scatter_scale = scale_values(iscale);

                fprintf('Scale calibration SSA1: f=%g Hz, Hs=%g m, scale=%g, seed=%d.\n', ...
                    paramsV.f0, paramsV.sea_hs_target, ...
                    paramsV.surface_ssa_scatter_scale, paramsV.sea_seed);
                channel = vertical_channel_model(paramsV);
                run_index = run_index + 1;
                run_rows = local_append_struct(run_rows, ...
                    local_run_row(run_index, "ssa1_geometry", paramsV, channel));
            end
        end
    end
end

run_table = struct2table(run_rows);
summary_table = local_build_summary_table(run_table);
[calibration_table, candidate_table] = local_build_calibration_table( ...
    summary_table, Hs_list, f_list_hz, scale_values);
global_scale_table = local_build_global_scale_table(candidate_table, scale_values);
[~, best_idx] = min(global_scale_table.objective_mean);
recommended_global_scale = global_scale_table.surface_ssa_scatter_scale(best_idx);

validation_report = struct();
validation_report.script = mfilename;
validation_report.created_at = char(datetime('now'));
validation_report.Hs_list = Hs_list;
validation_report.f_list_hz = f_list_hz;
validation_report.scale_values = scale_values;
validation_report.seed_list = seed_list;
validation_report.grid_n = grid_n;
validation_report.recommended_global_scale = recommended_global_scale;
validation_report.notes = ['surface_ssa_scatter_scale remains an engineering ', ...
    'normalization parameter fitted against Kirchhoff reduced-grid statistics, ', ...
    'not an experimentally calibrated scattering cross section.'];

local_plot_calibration(calibration_table, Hs_list, f_list_hz);

save(result_file, 'run_table', 'summary_table', 'calibration_table', 'candidate_table', ...
    'global_scale_table', 'validation_report', 'Hs_list', 'f_list_hz', ...
    'scale_values', 'seed_list', 'grid_n');

disp(calibration_table(:, {'Hs_target', 'f_hz', 'recommended_scale', ...
    'best_objective', 'best_abs_h_reflect_rel_error', 'best_rms_delta_k_rel_error'}))
disp(global_scale_table)
fprintf('Saved %s\n', result_file);
fprintf('Recommended global surface_ssa_scatter_scale = %.6g\n', recommended_global_scale);

function paramsV = local_base_params(grid_n, sea_wind_speed)
domain_width_m = min(50, 0.5 * grid_n);
paramsV = struct();
paramsV.f0 = 6000;
paramsV.c0 = 1500;
paramsV.enable_wideband = false;
paramsV.f_band_hz = [4000, 8000];
paramsV.Nf_min = 8;
paramsV.Nf_max = 8;
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
paramsV.sea_wind_speed = sea_wind_speed;
paramsV.sea_hs_target = 0.2;
paramsV.sea_seed = 12345;
paramsV.surface_reflect_coeff = -1;
paramsV.surface_phase_mode = 'normal';
paramsV.surface_boundary_model = 'kirchhoff_spatial';
paramsV.surface_boundary_redistribution_diagnostics = true;
paramsV.surface_ssa_random_scatter = true;
paramsV.surface_ssa_scatter_scale = 1.0;
paramsV.surface_ssa_seed_offset = 100000;
paramsV.surface_ssa_kernel_mode = 'ssa1_geometry';
paramsV.surface_ssa_coherent_order = 'ssa1';
paramsV.surface_ssa_frequency_correlation_mode = 'independent';
paramsV.surface_ssa_frequency_correlation_rho = 0.8;
paramsV.surface_ssa_geometry_source_id = ['vertical_comm_guide.md; Thorsos & Broschat 1995 JASA, ', ...
    'pressure-release / Dirichlet first-order perturbation-limit geometry'];
paramsV.surface_ssa_kz_branch = 'downward_positive_real';
paramsV.surface_ssa_conv_padding = 'periodic';
end

function row = local_run_row(run_index, model_name, paramsV, channel)
ssa = channel.roughness_meta.ssa_stat_kernel_meta;
redistribution = channel.roughness_meta.boundary_redistribution_diagnostics;
if abs(channel.h_direct) > 0
    ratio = channel.h_reflect / channel.h_direct;
else
    ratio = complex(NaN, NaN);
end
row = struct();
row.run_index = run_index;
row.model_name = model_name;
row.Hs_target = paramsV.sea_hs_target;
row.f_hz = paramsV.f0;
row.seed = paramsV.sea_seed;
row.surface_ssa_scatter_scale = paramsV.surface_ssa_scatter_scale;
row.abs_h_reflect = abs(channel.h_reflect);
row.abs_h_total = abs(channel.h_total);
row.reflection_direct_ratio_abs = abs(ratio);
row.phase_h_reflect = angle(channel.h_reflect);
row.invariant_error = max(abs(channel.H_f(:) - channel.H_direct_f(:) - channel.H_reflect_f(:)));
row.E_inc = local_get_field_or_nan(ssa, 'E_inc');
row.E_sca_limit = local_get_field_or_nan(ssa, 'E_sca_limit');
row.E_sca_limit_over_E_inc = row.E_sca_limit / max(row.E_inc, eps);
row.energy_conservation_error = local_get_field_or_nan(ssa, 'energy_conservation_error');
if model_name == "ssa1_geometry"
    row.reflected_rms_delta_k_rad_per_m = local_nested_field_or_nan(ssa, {'reflected_spectrum_stats', 'rms_delta_k_rad_per_m'});
else
    row.reflected_rms_delta_k_rad_per_m = local_nested_field_or_nan(redistribution, {'reflect_rms_delta_k_rad_per_m'});
end
end

function summary_table = local_build_summary_table(run_table)
models = unique(run_table.model_name, 'stable');
Hs_values = unique(run_table.Hs_target).';
f_values = unique(run_table.f_hz).';
scales = unique(run_table.surface_ssa_scatter_scale).';
rows = struct([]);
for im = 1:numel(models)
    for ih = 1:numel(Hs_values)
        for ifq = 1:numel(f_values)
            if models(im) == "kirchhoff_spatial"
                scale_iter = NaN;
            else
                scale_iter = scales;
            end
            for is = 1:numel(scale_iter)
                idx = run_table.model_name == models(im) & ...
                    run_table.Hs_target == Hs_values(ih) & run_table.f_hz == f_values(ifq);
                if models(im) ~= "kirchhoff_spatial"
                    idx = idx & run_table.surface_ssa_scatter_scale == scale_iter(is);
                end
                subset = run_table(idx, :);
                if isempty(subset)
                    continue
                end
                row = struct();
                row.model_name = models(im);
                row.Hs_target = Hs_values(ih);
                row.f_hz = f_values(ifq);
                row.surface_ssa_scatter_scale = scale_iter(is);
                row.seed_count = height(subset);
                row.abs_h_reflect_mean = mean(subset.abs_h_reflect);
                row.abs_h_reflect_std = std(subset.abs_h_reflect);
                row.reflection_direct_ratio_abs_mean = mean(subset.reflection_direct_ratio_abs);
                row.reflection_direct_ratio_abs_std = std(subset.reflection_direct_ratio_abs);
                row.reflected_rms_delta_k_mean = local_nanmean(subset.reflected_rms_delta_k_rad_per_m);
                row.reflected_rms_delta_k_std = local_nanstd(subset.reflected_rms_delta_k_rad_per_m);
                row.E_sca_limit_over_E_inc_mean = local_nanmean(subset.E_sca_limit_over_E_inc);
                row.energy_conservation_error_max = max(subset.energy_conservation_error);
                row.invariant_error_max = max(subset.invariant_error);
                rows = local_append_struct(rows, row);
            end
        end
    end
end
summary_table = struct2table(rows);
end

function [calibration_table, candidate_table] = local_build_calibration_table(summary_table, Hs_list, f_list_hz, scale_values)
rows = struct([]);
candidate_rows = struct([]);
for ih = 1:numel(Hs_list)
    for ifq = 1:numel(f_list_hz)
        ref = summary_table(summary_table.model_name == "kirchhoff_spatial" & ...
            summary_table.Hs_target == Hs_list(ih) & summary_table.f_hz == f_list_hz(ifq), :);
        candidates = summary_table(summary_table.model_name == "ssa1_geometry" & ...
            summary_table.Hs_target == Hs_list(ih) & summary_table.f_hz == f_list_hz(ifq), :);
        objective = NaN(height(candidates), 1);
        abs_err = NaN(height(candidates), 1);
        ratio_err = NaN(height(candidates), 1);
        rms_err = NaN(height(candidates), 1);
        for ii = 1:height(candidates)
            abs_err(ii) = local_rel_error(candidates.abs_h_reflect_mean(ii), ref.abs_h_reflect_mean);
            ratio_err(ii) = local_rel_error(candidates.reflection_direct_ratio_abs_mean(ii), ref.reflection_direct_ratio_abs_mean);
            rms_err(ii) = local_rel_error(candidates.reflected_rms_delta_k_mean(ii), ref.reflected_rms_delta_k_mean);
            objective(ii) = local_nanmean([abs_err(ii)^2, ratio_err(ii)^2, rms_err(ii)^2]);
            candidate = struct();
            candidate.Hs_target = Hs_list(ih);
            candidate.f_hz = f_list_hz(ifq);
            candidate.surface_ssa_scatter_scale = candidates.surface_ssa_scatter_scale(ii);
            candidate.objective = objective(ii);
            candidate.abs_h_reflect_rel_error = abs_err(ii);
            candidate.ratio_rel_error = ratio_err(ii);
            candidate.rms_delta_k_rel_error = rms_err(ii);
            candidate.E_sca_limit_over_E_inc_mean = candidates.E_sca_limit_over_E_inc_mean(ii);
            candidate_rows = local_append_struct(candidate_rows, candidate);
        end
        [best_objective, best_idx] = min(objective);
        row = struct();
        row.Hs_target = Hs_list(ih);
        row.f_hz = f_list_hz(ifq);
        row.recommended_scale = candidates.surface_ssa_scatter_scale(best_idx);
        row.best_objective = best_objective;
        row.best_abs_h_reflect_rel_error = abs_err(best_idx);
        row.best_ratio_rel_error = ratio_err(best_idx);
        row.best_rms_delta_k_rel_error = rms_err(best_idx);
        row.kirchhoff_abs_h_reflect_mean = ref.abs_h_reflect_mean;
        row.ssa_abs_h_reflect_mean = candidates.abs_h_reflect_mean(best_idx);
        row.kirchhoff_rms_delta_k_mean = ref.reflected_rms_delta_k_mean;
        row.ssa_rms_delta_k_mean = candidates.reflected_rms_delta_k_mean(best_idx);
        row.E_sca_limit_over_E_inc_mean = candidates.E_sca_limit_over_E_inc_mean(best_idx);
        row.scale_values_tested = string(mat2str(scale_values));
        rows = local_append_struct(rows, row);
    end
end
calibration_table = struct2table(rows);
candidate_table = struct2table(candidate_rows);
end

function global_scale_table = local_build_global_scale_table(candidate_table, scale_values)
rows = struct([]);
for ii = 1:numel(scale_values)
    idx = candidate_table.surface_ssa_scatter_scale == scale_values(ii);
    row = struct();
    row.surface_ssa_scatter_scale = scale_values(ii);
    row.condition_count = nnz(idx);
    if any(idx)
        row.objective_mean = mean(candidate_table.objective(idx));
        row.objective_median = median(candidate_table.objective(idx));
        row.abs_h_reflect_rel_error_mean = mean(candidate_table.abs_h_reflect_rel_error(idx));
        row.rms_delta_k_rel_error_mean = mean(candidate_table.rms_delta_k_rel_error(idx));
    else
        row.objective_mean = Inf;
        row.objective_median = Inf;
        row.abs_h_reflect_rel_error_mean = Inf;
        row.rms_delta_k_rel_error_mean = Inf;
    end
    rows = local_append_struct(rows, row);
end
global_scale_table = struct2table(rows);
end

function local_plot_calibration(calibration_table, Hs_list, f_list_hz)
fig = figure('Visible', 'off');
tiledlayout(1, 2, 'TileSpacing', 'compact')
nexttile
hold on
for ii = 1:numel(f_list_hz)
    subset = sortrows(calibration_table(calibration_table.f_hz == f_list_hz(ii), :), 'Hs_target');
    plot(subset.Hs_target, subset.recommended_scale, '-o', 'DisplayName', sprintf('%g Hz', f_list_hz(ii)))
end
grid on
xlabel('H_s (m)')
ylabel('Recommended scale')
title('SSA scatter scale calibration')
legend('Location', 'best')

nexttile
hold on
for ii = 1:numel(f_list_hz)
    subset = sortrows(calibration_table(calibration_table.f_hz == f_list_hz(ii), :), 'Hs_target');
    plot(subset.Hs_target, sqrt(subset.best_objective), '-o', 'DisplayName', sprintf('%g Hz', f_list_hz(ii)))
end
grid on
xlabel('H_s (m)')
ylabel('RMS normalized residual')
title('Calibration residual')
legend('Location', 'best')
saveas(fig, 'calibrate_ssa_scatter_scale_summary.png')
close(fig)
end

function err = local_rel_error(value, reference)
err = abs(value - reference) / max(abs(reference), eps);
end

function value = local_get_field_or_nan(s, field_name)
if isstruct(s) && isfield(s, field_name)
    value = s.(field_name);
else
    value = NaN;
end
end

function value = local_nested_field_or_nan(s, path)
value = NaN;
cur = s;
for ii = 1:numel(path)
    if isstruct(cur) && isfield(cur, path{ii})
        cur = cur.(path{ii});
    else
        return
    end
end
if isnumeric(cur) && isscalar(cur)
    value = cur;
end
end

function rows = local_append_struct(rows, row)
if isempty(rows)
    rows = row;
else
    rows(end+1) = row; %#ok<AGROW>
end
end

function value = local_nanmean(x)
x = x(isfinite(x));
if isempty(x)
    value = NaN;
else
    value = mean(x);
end
end

function value = local_nanstd(x)
x = x(isfinite(x));
if numel(x) <= 1
    value = 0;
else
    value = std(x);
end
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

function values = local_env_numeric_vector(name, default_values)
raw = getenv(name);
if isempty(raw)
    values = default_values;
    return
end
parts = regexp(strtrim(raw), '[,;\s]+', 'split');
values = str2double(parts);
if isempty(values) || any(~isfinite(values))
    error('%s must be a numeric vector.', name);
end
end

