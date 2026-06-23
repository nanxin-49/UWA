% Reduced-grid statistical comparison between Kirchhoff realizations and SSA1.
% This script validates rough-surface model trends only. It does not run the
% communication chain and does not modify propagation or receiver logic.

clear
format compact

result_file = getenv('SSA1_KIRCHHOFF_COMPARE_RESULT_FILE');
if isempty(result_file)
    result_file = 'compare_ssa1_kirchhoff_surface_statistics_vertical_result.mat';
end
if exist(result_file, 'file')
    warning('Result file %s already exists and will be overwritten.', result_file);
end

Hs_list = local_env_numeric_vector('SSA1_KIRCHHOFF_COMPARE_HS_LIST', [0, 0.05, 0.2, 0.5, 1.0]);
f_list_hz = local_env_numeric_vector('SSA1_KIRCHHOFF_COMPARE_F_LIST_HZ', [4000, 6000, 8000, 10000]);
seed_count = local_env_scalar_int('SSA1_KIRCHHOFF_COMPARE_SEED_COUNT', 16);
seed_list = 12345 + (0:(seed_count - 1));
grid_n = local_env_scalar_int('SSA1_KIRCHHOFF_COMPARE_GRID_N', 64);
include_pm = local_env_flag('SSA1_KIRCHHOFF_COMPARE_INCLUDE_PM', false);

sea_wind_speed = 5.0;
model_specs = local_model_specs(include_pm);
model_names = string({model_specs.model_name});

run_rows = struct([]);
run_index = 0;

for im = 1:numel(model_specs)
    for ff = 1:numel(f_list_hz)
        for hh = 1:numel(Hs_list)
            for ss = 1:numel(seed_list)
                paramsV = local_base_params(sea_wind_speed, grid_n);
                paramsV.f0 = f_list_hz(ff);
                paramsV.f_ref_hz = paramsV.f0;
                paramsV.sea_hs_target = Hs_list(hh);
                paramsV.sea_seed = seed_list(ss);
                paramsV.surface_boundary_model = model_specs(im).surface_boundary_model;
                if strcmp(model_specs(im).surface_boundary_model, 'ssa_stat_kernel')
                    paramsV.surface_ssa_kernel_mode = model_specs(im).surface_ssa_kernel_mode;
                else
                    paramsV.surface_ssa_kernel_mode = 'pm_convolution';
                end
                paramsV.surface_ssa_geometry_source_id = model_specs(im).surface_ssa_geometry_source_id;
                paramsV.surface_ssa_random_scatter = model_specs(im).surface_ssa_random_scatter;

                fprintf('Compare SSA1/Kirchhoff: model=%s, f=%g Hz, Hs=%g m, seed=%d (%d/%d).\n', ...
                    model_specs(im).model_name, paramsV.f0, paramsV.sea_hs_target, ...
                    paramsV.sea_seed, ss, numel(seed_list));
                channel = vertical_channel_model(paramsV);

                run_index = run_index + 1;
                run_rows = local_append_struct(run_rows, ...
                    local_build_run_row(run_index, model_specs(im), paramsV, channel));
            end
        end
    end
end

run_table = struct2table(run_rows);
summary_table = local_build_summary_table(run_table, model_names, Hs_list, f_list_hz);
trend_table = local_build_trend_table(run_table, summary_table, Hs_list, f_list_hz);

hard_rows = trend_table(trend_table.check_type == "hard_pass", :);
compat_rows = trend_table(trend_table.check_type == "compatibility_observed", :);
validation_report = struct();
validation_report.script = mfilename;
validation_report.created_at = char(datetime('now'));
validation_report.Hs_list = Hs_list;
validation_report.f_list_hz = f_list_hz;
validation_report.seed_list = seed_list;
validation_report.seed_count = numel(seed_list);
validation_report.grid_n = grid_n;
validation_report.model_names = model_names;
validation_report.include_pm_convolution = include_pm;
validation_report.hard_checks_passed = all(hard_rows.passed);
validation_report.compatibility_observed_fraction = local_nanmean(double(compat_rows.passed));
validation_report.notes = ['Kirchhoff is a realization-based phase-screen model; ', ...
    'ssa1_geometry is a PM-spectrum-driven first-order pressure-release Dirichlet ', ...
    'statistical reflection/scattering model. The comparison checks statistical ', ...
    'trend compatibility, not sample-wise equivalence.'];

local_plot_all(summary_table, f_list_hz, model_names);

save(result_file, 'run_table', 'summary_table', 'trend_table', 'validation_report', ...
    'Hs_list', 'f_list_hz', 'seed_list', 'model_names', 'grid_n', 'sea_wind_speed');

disp(summary_table(:, {'model_name', 'Hs_target', 'f_hz', 'seed_count', ...
    'abs_h_reflect_mean', 'abs_h_reflect_std', 'phase_h_reflect_circular_variance', ...
    'coherent_h_reflect_loss_proxy', 'reflected_rms_delta_k_mean', ...
    'reflect_high_k_fraction_mean'}))
disp(trend_table)
fprintf('Saved %s\n', result_file);
fprintf('Hard checks passed: %d\n', validation_report.hard_checks_passed);
fprintf('Compatibility observed fraction: %.3f\n', validation_report.compatibility_observed_fraction);

if ~validation_report.hard_checks_passed
    error('compare_ssa1_kirchhoff_surface_statistics_vertical:HardCheckFailed', ...
        'One or more hard validation checks failed.');
end

function paramsV = local_base_params(sea_wind_speed, grid_n)
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
paramsV.sea_hs_target = 0.5;
paramsV.sea_seed = 12345;
paramsV.surface_reflect_coeff = -1;
paramsV.surface_boundary_model = 'kirchhoff_spatial';
paramsV.surface_boundary_coupling_diagnostics = false;
paramsV.surface_boundary_redistribution_diagnostics = true;
paramsV.surface_boundary_redistribution_debug = false;
paramsV.surface_ssa_random_scatter = true;
paramsV.surface_ssa_scatter_scale = 1.0;
paramsV.surface_ssa_seed_offset = 100000;
paramsV.surface_ssa_kernel_mode = 'pm_convolution';
paramsV.surface_ssa_geometry_source_id = '';
paramsV.surface_ssa_kz_branch = 'downward_positive_real';
paramsV.surface_ssa_conv_padding = 'periodic';
end

function specs = local_model_specs(include_pm)
source_id = ['SSA.md; Thorsos & Broschat 1995 JASA, ', ...
    'pressure-release / Dirichlet first-order perturbation-limit geometry'];
specs = struct([]);
specs(1).model_name = 'kirchhoff_spatial';
specs(1).surface_boundary_model = 'kirchhoff_spatial';
specs(1).surface_ssa_kernel_mode = 'not_applicable';
specs(1).surface_ssa_geometry_source_id = '';
specs(1).surface_ssa_random_scatter = false;
specs(2).model_name = 'ssa1_geometry';
specs(2).surface_boundary_model = 'ssa_stat_kernel';
specs(2).surface_ssa_kernel_mode = 'ssa1_geometry';
specs(2).surface_ssa_geometry_source_id = source_id;
specs(2).surface_ssa_random_scatter = true;
if include_pm
    specs(3).model_name = 'ssa_pm_convolution';
    specs(3).surface_boundary_model = 'ssa_stat_kernel';
    specs(3).surface_ssa_kernel_mode = 'pm_convolution';
    specs(3).surface_ssa_geometry_source_id = '';
    specs(3).surface_ssa_random_scatter = true;
end
end

function row = local_build_run_row(run_index, spec, paramsV, channel)
ssa = channel.roughness_meta.ssa_stat_kernel_meta;
redistribution = channel.roughness_meta.boundary_redistribution_diagnostics;
if abs(channel.h_direct) > 0
    ratio = channel.h_reflect / channel.h_direct;
else
    ratio = complex(NaN, NaN);
end
row = struct();
row.run_index = run_index;
row.model_name = string(spec.model_name);
row.Hs_target = paramsV.sea_hs_target;
row.f_hz = paramsV.f0;
row.seed = paramsV.sea_seed;
row.boundary_model = string(spec.surface_boundary_model);
row.kernel_mode = string(spec.surface_ssa_kernel_mode);
row.h_reflect = channel.h_reflect;
row.abs_h_reflect = abs(channel.h_reflect);
row.phase_h_reflect = angle(channel.h_reflect);
row.h_total = channel.h_total;
row.abs_h_total = abs(channel.h_total);
row.h_direct = channel.h_direct;
row.abs_h_direct = abs(channel.h_direct);
row.reflection_direct_ratio = ratio;
row.reflection_direct_ratio_abs = abs(ratio);
row.reflection_direct_ratio_phase = angle(ratio);
row.H_f = channel.H_f(channel.idx_f_ref);
row.H_reflect_f = channel.H_reflect_f(channel.idx_f_ref);
row.invariant_error = max(abs(channel.H_f(:) - channel.H_direct_f(:) - channel.H_reflect_f(:)));
row.R_coh = local_get_field_or_nan(ssa, 'R_coh');
row.abs_R_coh = abs(row.R_coh);
row.E_inc = local_get_field_or_nan(ssa, 'E_inc');
row.E_coh = local_get_field_or_nan(ssa, 'E_coh');
row.E_sca = local_get_field_or_nan(ssa, 'E_sca');
row.E_sca_limit = local_get_field_or_nan(ssa, 'E_sca_limit');
row.E_ref = local_get_field_or_nan(ssa, 'E_ref');
row.energy_conservation_error = local_get_field_or_nan(ssa, 'energy_conservation_error');
row.E_coh_over_E_inc = row.E_coh / max(row.E_inc, eps);
row.E_sca_limit_over_E_inc = row.E_sca_limit / max(row.E_inc, eps);
row.ssa_boundary_condition = string(local_get_field_or_text(ssa, 'boundary_condition', 'not_applicable'));
row.ssa_G_SSA1_formula = string(local_nested_text_or_default(ssa, {'kernel_detail', 'G_SSA1_formula'}, 'not_applicable'));
row.ssa_formula_source = string(local_get_field_or_text(ssa, 'formula_source', 'not_applicable'));

if strcmp(spec.surface_boundary_model, 'ssa_stat_kernel')
    row.reflected_rms_delta_k_rad_per_m = local_nested_field_or_nan(ssa, {'reflected_spectrum_stats', 'rms_delta_k_rad_per_m'});
    row.reflected_energy_radius_90_rad_per_m = local_nested_field_or_nan(ssa, {'reflected_spectrum_stats', 'energy_radius_90_rad_per_m'});
    row.reflected_centroid_shift_mag_rad_per_m = NaN;
    row.reflect_high_k_fraction = NaN;
    row.boundary_reflect_rms_delta_k_rad_per_m = local_nested_field_or_nan(redistribution, {'reflect_rms_delta_k_rad_per_m'});
    row.boundary_reflect_high_k_fraction = local_nested_field_or_nan(redistribution, {'reflect_high_k_fraction'});
else
    row.reflected_rms_delta_k_rad_per_m = local_nested_field_or_nan(redistribution, {'reflect_rms_delta_k_rad_per_m'});
    row.reflected_energy_radius_90_rad_per_m = local_nested_field_or_nan(redistribution, {'reflect_energy_radius_90_rad_per_m'});
    row.reflected_centroid_shift_mag_rad_per_m = local_nested_field_or_nan(redistribution, {'centroid_shift_mag_rad_per_m'});
    row.reflect_high_k_fraction = local_nested_field_or_nan(redistribution, {'reflect_high_k_fraction'});
    row.boundary_reflect_rms_delta_k_rad_per_m = row.reflected_rms_delta_k_rad_per_m;
    row.boundary_reflect_high_k_fraction = row.reflect_high_k_fraction;
end
end

function summary_table = local_build_summary_table(run_table, model_names, Hs_list, f_list_hz)
rows = struct([]);
for im = 1:numel(model_names)
    model_name = model_names(im);
    for ff = 1:numel(f_list_hz)
        f_hz = f_list_hz(ff);
        flat_subset = run_table(run_table.model_name == model_name & ...
            run_table.f_hz == f_hz & run_table.Hs_target == 0, :);
        flat_coherent_ref = abs(mean(flat_subset.h_reflect));
        for hh = 1:numel(Hs_list)
            Hs = Hs_list(hh);
            subset = run_table(run_table.model_name == model_name & ...
                run_table.f_hz == f_hz & run_table.Hs_target == Hs, :);
            if isempty(subset)
                continue
            end
            row = struct();
            row.model_name = model_name;
            row.Hs_target = Hs;
            row.f_hz = f_hz;
            row.seed_count = height(subset);
            row.abs_h_reflect_mean = mean(subset.abs_h_reflect);
            row.abs_h_reflect_std = std(subset.abs_h_reflect);
            row.abs_h_total_mean = mean(subset.abs_h_total);
            row.abs_h_total_std = std(subset.abs_h_total);
            row.phase_h_reflect_circular_mean = local_circular_mean(subset.phase_h_reflect);
            row.phase_h_reflect_circular_variance = local_circular_variance(subset.phase_h_reflect);
            row.reflection_direct_ratio_abs_mean = mean(subset.reflection_direct_ratio_abs);
            row.reflection_direct_ratio_abs_std = std(subset.reflection_direct_ratio_abs);
            row.abs_R_coh_mean = local_nanmean(subset.abs_R_coh);
            row.abs_R_coh_std = local_nanstd(subset.abs_R_coh);
            row.E_coh_over_E_inc_mean = local_nanmean(subset.E_coh_over_E_inc);
            row.E_coh_over_E_inc_std = local_nanstd(subset.E_coh_over_E_inc);
            row.E_sca_limit_over_E_inc_mean = local_nanmean(subset.E_sca_limit_over_E_inc);
            row.E_sca_limit_over_E_inc_std = local_nanstd(subset.E_sca_limit_over_E_inc);
            row.reflected_rms_delta_k_mean = local_nanmean(subset.reflected_rms_delta_k_rad_per_m);
            row.reflected_rms_delta_k_std = local_nanstd(subset.reflected_rms_delta_k_rad_per_m);
            row.reflected_energy_radius_90_mean = local_nanmean(subset.reflected_energy_radius_90_rad_per_m);
            row.reflected_energy_radius_90_std = local_nanstd(subset.reflected_energy_radius_90_rad_per_m);
            row.reflect_high_k_fraction_mean = local_nanmean(subset.reflect_high_k_fraction);
            row.reflect_high_k_fraction_std = local_nanstd(subset.reflect_high_k_fraction);
            row.centroid_shift_mag_mean = local_nanmean(subset.reflected_centroid_shift_mag_rad_per_m);
            row.centroid_shift_mag_std = local_nanstd(subset.reflected_centroid_shift_mag_rad_per_m);
            row.coherent_h_reflect_loss_proxy = abs(mean(subset.h_reflect)) / max(flat_coherent_ref, eps);
            row.invariant_error_max = max(subset.invariant_error);
            row.energy_conservation_error_max = local_nanmax(subset.energy_conservation_error);
            rows = local_append_struct(rows, row);
        end
    end
end
summary_table = struct2table(rows);
end

function trend_table = local_build_trend_table(run_table, summary_table, Hs_list, f_list_hz)
checks = struct([]);
checks = local_add_check(checks, 'all_channel_invariants', 'hard_pass', ...
    max(run_table.invariant_error), 1e-10, '<=', 'H_f = H_direct_f + H_reflect_f for all runs.');

for ff = 1:numel(f_list_hz)
    f_hz = f_list_hz(ff);
    for ss = unique(run_table.seed).'
        k = run_table(run_table.model_name == "kirchhoff_spatial" & run_table.f_hz == f_hz & ...
            run_table.Hs_target == 0 & run_table.seed == ss, :);
        s = run_table(run_table.model_name == "ssa1_geometry" & run_table.f_hz == f_hz & ...
            run_table.Hs_target == 0 & run_table.seed == ss, :);
        if ~isempty(k) && ~isempty(s)
            checks = local_add_check(checks, sprintf('Hs0_h_reflect_match_f%d_seed%d', f_hz, ss), ...
                'hard_pass', abs(k.h_reflect - s.h_reflect), 1e-10, '<=', ...
                'Hs=0 SSA1 and Kirchhoff flat reflected channel should match.');
            checks = local_add_check(checks, sprintf('Hs0_h_total_match_f%d_seed%d', f_hz, ss), ...
                'hard_pass', abs(k.h_total - s.h_total), 1e-10, '<=', ...
                'Hs=0 SSA1 and Kirchhoff flat total channel should match.');
        end
    end

    ssa_rows = sortrows(summary_table(summary_table.model_name == "ssa1_geometry" & ...
        summary_table.f_hz == f_hz, :), 'Hs_target');
    checks = local_add_check(checks, sprintf('ssa1_Rcoh_nonincreasing_Hs_f%d', f_hz), ...
        'hard_pass', max(diff(ssa_rows.abs_R_coh_mean)), 1e-12, '<=', ...
        'SSA1 coherent reflection should decrease with roughness.');
    checks = local_add_check(checks, sprintf('ssa1_Ecoh_nonincreasing_Hs_f%d', f_hz), ...
        'hard_pass', max(diff(ssa_rows.E_coh_over_E_inc_mean)), 1e-12, '<=', ...
        'SSA1 coherent energy fraction should decrease with roughness.');
    checks = local_add_check(checks, sprintf('ssa1_Esca_limit_nondecreasing_Hs_f%d', f_hz), ...
        'hard_pass', min(diff(ssa_rows.E_sca_limit_over_E_inc_mean)), -1e-12, '>=', ...
        'SSA1 scatter budget should increase with roughness.');

    k_rows = sortrows(summary_table(summary_table.model_name == "kirchhoff_spatial" & ...
        summary_table.f_hz == f_hz, :), 'Hs_target');
    high_k_change = k_rows.reflect_high_k_fraction_mean(end) - k_rows.reflect_high_k_fraction_mean(1);
    rms_change = k_rows.reflected_rms_delta_k_mean(end) - k_rows.reflected_rms_delta_k_mean(1);
    checks = local_add_check(checks, sprintf('kirchhoff_spread_increases_Hs_f%d', f_hz), ...
        'compatibility_observed', max(high_k_change, rms_change), 0, '>', ...
        'Kirchhoff reflected spectrum should broaden from Hs=0 to the roughest tested sea.');
end

for hh = 2:numel(Hs_list)
    Hs = Hs_list(hh);
    ssa_rows = sortrows(summary_table(summary_table.model_name == "ssa1_geometry" & ...
        summary_table.Hs_target == Hs, :), 'f_hz');
    checks = local_add_check(checks, sprintf('ssa1_Rcoh_nonincreasing_frequency_Hs%g', Hs), ...
        'hard_pass', max(diff(ssa_rows.abs_R_coh_mean)), 1e-12, '<=', ...
        'SSA1 coherent loss should strengthen with frequency for Hs>0.');

    k_rows = sortrows(summary_table(summary_table.model_name == "kirchhoff_spatial" & ...
        summary_table.Hs_target == Hs, :), 'f_hz');
    phase_var_change = k_rows.phase_h_reflect_circular_variance(end) - ...
        k_rows.phase_h_reflect_circular_variance(1);
    rms_change = k_rows.reflected_rms_delta_k_mean(end) - k_rows.reflected_rms_delta_k_mean(1);
    checks = local_add_check(checks, sprintf('kirchhoff_frequency_variability_observed_Hs%g', Hs), ...
        'compatibility_observed', max(phase_var_change, rms_change), 0, '>', ...
        'Kirchhoff frequency trend is recorded as compatibility evidence, not a formula pass/fail.');
end

ssa_run_rows = run_table(run_table.model_name == "ssa1_geometry", :);
checks = local_add_check(checks, 'ssa1_energy_conservation', 'hard_pass', ...
    local_nanmax(ssa_run_rows.energy_conservation_error), 1e-12, '<=', ...
    'SSA1 energy audit should remain bounded.');
checks = local_add_check(checks, 'ssa1_G_formula_recorded', 'hard_pass', ...
    double(all(contains(ssa_run_rows.ssa_G_SSA1_formula, '4*gamma'))), 1, '==', ...
    'SSA1 metadata must record the first-order Dirichlet geometry factor.');
checks = local_add_check(checks, 'ssa1_dirichlet_boundary_recorded', 'hard_pass', ...
    double(all(ssa_run_rows.ssa_boundary_condition == "pressure-release / Dirichlet")), 1, '==', ...
    'SSA1 metadata must record pressure-release / Dirichlet boundary condition.');

trend_table = struct2table(checks);
end

function checks = local_add_check(checks, check_name, check_type, value, tolerance, comparison, interpretation)
row = struct();
row.check_name = string(check_name);
row.check_type = string(check_type);
row.value = double(value);
row.tolerance = double(tolerance);
row.comparison = string(comparison);
row.interpretation = string(interpretation);
switch comparison
    case '<='
        row.passed = value <= tolerance;
    case '>='
        row.passed = value >= tolerance;
    case '>'
        row.passed = value > tolerance;
    case '=='
        row.passed = value == tolerance;
    otherwise
        error('Unsupported comparison %s.', comparison);
end
checks = local_append_struct(checks, row);
end

function local_plot_all(summary_table, f_list_hz, model_names)
local_plot_metric(summary_table, f_list_hz, model_names, ...
    'abs_h_reflect_mean', 'abs_h_reflect_std', ...
    '|h_{reflect}|', 'Reflected channel amplitude vs H_s', ...
    'compare_ssa1_kirchhoff_abs_h_reflect_vs_Hs.png');
local_plot_metric(summary_table, f_list_hz, model_names, ...
    'abs_R_coh_mean', 'abs_R_coh_std', ...
    '|R_{coh}|', 'SSA1 coherent reflection vs H_s', ...
    'compare_ssa1_kirchhoff_abs_R_coh_vs_Hs.png');
local_plot_metric(summary_table, f_list_hz, model_names, ...
    'E_coh_over_E_inc_mean', 'E_coh_over_E_inc_std', ...
    'E_{coh}/E_{inc}', 'Coherent energy fraction vs H_s', ...
    'compare_ssa1_kirchhoff_coherent_fraction_vs_Hs.png');
local_plot_metric(summary_table, f_list_hz, model_names, ...
    'E_sca_limit_over_E_inc_mean', 'E_sca_limit_over_E_inc_std', ...
    'E_{sca}^{max}/E_{inc}', 'SSA1 scatter budget vs H_s', ...
    'compare_ssa1_kirchhoff_scatter_budget_vs_Hs.png');
local_plot_metric(summary_table, f_list_hz, model_names, ...
    'reflected_rms_delta_k_mean', 'reflected_rms_delta_k_std', ...
    'RMS \DeltaK (rad/m)', 'Reflected angular-spectrum width vs H_s', ...
    'compare_ssa1_kirchhoff_rms_delta_k_vs_Hs.png');
local_plot_metric(summary_table, f_list_hz, model_names, ...
    'phase_h_reflect_circular_variance', '', ...
    'Circular phase variance', 'Reflected-channel phase variance vs H_s', ...
    'compare_ssa1_kirchhoff_phase_variance_vs_Hs.png');
end

function local_plot_metric(summary_table, f_list_hz, model_names, mean_field, std_field, y_label, plot_title, file_name)
fig = figure('Visible', 'off');
for ff = 1:numel(f_list_hz)
    subplot(2, 2, ff)
    hold on
    for im = 1:numel(model_names)
        subset = sortrows(summary_table(summary_table.f_hz == f_list_hz(ff) & ...
            summary_table.model_name == model_names(im), :), 'Hs_target');
        if isempty(subset)
            continue
        end
        y = subset.(mean_field);
        if all(isnan(y))
            continue
        end
        if ~isempty(std_field)
            e = subset.(std_field);
            errorbar(subset.Hs_target, y, e, '-o', 'DisplayName', char(model_names(im)));
        else
            plot(subset.Hs_target, y, '-o', 'DisplayName', char(model_names(im)));
        end
    end
    grid on
    xlabel('H_s (m)')
    ylabel(y_label)
    title(sprintf('%s, f=%g Hz', plot_title, f_list_hz(ff)))
    legend('Location', 'best')
end
saveas(fig, file_name)
close(fig)
end

function out = local_env_numeric_vector(name, default_value)
txt = strtrim(getenv(name));
if isempty(txt)
    out = default_value;
else
    out = str2num(txt); %#ok<ST2NM>
    if isempty(out) || any(~isfinite(out))
        error('Environment variable %s must contain a numeric vector.', name);
    end
end
out = out(:).';
end

function out = local_env_scalar_int(name, default_value)
txt = strtrim(getenv(name));
if isempty(txt)
    out = default_value;
else
    out = round(str2double(txt));
    if ~isfinite(out) || out < 1
        error('Environment variable %s must be a positive integer.', name);
    end
end
end

function out = local_env_flag(name, default_value)
txt = lower(strtrim(getenv(name)));
if isempty(txt)
    out = default_value;
else
    out = any(strcmp(txt, {'1', 'true', 'yes', 'on'}));
end
end

function s = local_append_struct(s, row)
if isempty(s)
    s = row;
else
    s(end + 1) = row; %#ok<AGROW>
end
end

function value = local_get_field_or_nan(s, field_name)
if isstruct(s) && isfield(s, field_name)
    value = s.(field_name);
else
    value = NaN;
end
if ~(isnumeric(value) && isscalar(value))
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

function out = local_nested_field_or_nan(s, field_path)
out = s;
for ii = 1:numel(field_path)
    if isstruct(out) && isfield(out, field_path{ii})
        out = out.(field_path{ii});
    else
        out = NaN;
        return
    end
end
if ~(isnumeric(out) && isscalar(out) && isfinite(out))
    out = NaN;
end
end

function out = local_nested_text_or_default(s, field_path, default_value)
out = s;
for ii = 1:numel(field_path)
    if isstruct(out) && isfield(out, field_path{ii})
        out = out.(field_path{ii});
    else
        out = default_value;
        return
    end
end
end

function m = local_nanmean(x)
x = x(isfinite(x));
if isempty(x)
    m = NaN;
else
    m = mean(x);
end
end

function s = local_nanstd(x)
x = x(isfinite(x));
if isempty(x)
    s = NaN;
elseif numel(x) <= 1
    s = 0;
else
    s = std(x);
end
end

function m = local_nanmax(x)
x = x(isfinite(x));
if isempty(x)
    m = NaN;
else
    m = max(x);
end
end

function mu = local_circular_mean(theta)
theta = theta(isfinite(theta));
if isempty(theta)
    mu = NaN;
else
    mu = angle(mean(exp(1i * theta)));
end
end

function v = local_circular_variance(theta)
theta = theta(isfinite(theta));
if isempty(theta)
    v = NaN;
else
    v = 1 - abs(mean(exp(1i * theta)));
end
end
