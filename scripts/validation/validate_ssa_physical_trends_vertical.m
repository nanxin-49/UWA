run(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'bootstrap_project.m'));
% Reduced-grid physical trend validation for the SSA statistical sea-surface branch.
% This script does not run the communication chain. It checks coherent
% reflection, PM spectrum variance normalization, energy audit metadata, and
% first-order Dirichlet SSA geometry labels over Hs and frequency sweeps.

clear
format compact

result_file = getenv('SSA_PHYSICAL_TRENDS_RESULT_FILE');
if isempty(result_file)
    result_file = 'validate_ssa_physical_trends_vertical_result.mat';
end

Hs_list = [0, 0.05, 0.2, 0.5, 1.0];
f_list_hz = [4000, 6000, 8000, 10000];
kernel_modes = {'pm_convolution', 'ssa1_geometry'};
geometry_source_id = ['vertical_comm_guide.md; Thorsos & Broschat 1995 JASA, ', ...
    'pressure-release / Dirichlet first-order perturbation-limit geometry'];

tol_roundoff = 1e-12;
tol_w_eta_rel = 1e-10;
tol_monotone = 1e-12;
tol_energy = 1e-12;

rows = struct([]);
checks = struct([]);
run_index = 0;

for kk = 1:numel(kernel_modes)
    for ff = 1:numel(f_list_hz)
        for hh = 1:numel(Hs_list)
            paramsV = local_base_params();
            paramsV.f0 = f_list_hz(ff);
            paramsV.f_ref_hz = paramsV.f0;
            paramsV.surface_boundary_model = 'ssa_stat_kernel';
            paramsV.surface_ssa_kernel_mode = kernel_modes{kk};
            paramsV.surface_ssa_random_scatter = false;
            paramsV.surface_ssa_conv_padding = 'periodic';
            paramsV.sea_hs_target = Hs_list(hh);
            paramsV.sea_wind_speed = 5;
            paramsV.sea_seed = 12345;
            if strcmp(kernel_modes{kk}, 'ssa1_geometry')
                paramsV.surface_ssa_geometry_source_id = geometry_source_id;
            end

            fprintf('Running SSA trend case: kernel=%s, f=%g Hz, Hs=%g m.\n', ...
                kernel_modes{kk}, paramsV.f0, paramsV.sea_hs_target);
            channel = vertical_channel_model(paramsV);
            meta = channel.roughness_meta.ssa_stat_kernel_meta;

            run_index = run_index + 1;
            rows = local_append_row(rows, local_summary_row(run_index, paramsV, channel, meta));
        end
    end
end

summary_table = struct2table(rows);

checks = local_add_check(checks, 'channel_invariant_all', ...
    max(summary_table.invariant_error), 1e-10, '<=');
checks = local_add_check(checks, 'energy_conservation_all', ...
    max(summary_table.energy_conservation_error), tol_energy, '<=');
checks = local_add_check(checks, 'W_eta_variance_rel_error_all', ...
    max(summary_table.W_eta_variance_rel_error), tol_w_eta_rel, '<=');
checks = local_add_check(checks, 'normal_coherent_formula_all', ...
    max(summary_table.R_coh_formula_abs_error), tol_roundoff, '<=');

for kk = 1:numel(kernel_modes)
    kernel = string(kernel_modes{kk});
    kernel_name = char(kernel);
    kernel_rows = summary_table(summary_table.kernel_mode == kernel, :);

    hs0_rows = kernel_rows(kernel_rows.Hs_target == 0, :);
    checks = local_add_check(checks, sprintf('%s_Hs0_abs_R_coh_equals_R0', kernel_name), ...
        max(abs(hs0_rows.abs_R_coh - abs(hs0_rows.R0))), tol_roundoff, '<=');
    checks = local_add_check(checks, sprintf('%s_Hs0_E_sca_raw_zero', kernel_name), ...
        max(abs(hs0_rows.E_sca_raw)), tol_roundoff, '<=');
    checks = local_add_check(checks, sprintf('%s_Hs0_E_sca_zero', kernel_name), ...
        max(abs(hs0_rows.E_sca)), tol_roundoff, '<=');
    checks = local_add_check(checks, sprintf('%s_Hs0_W_eta_variance_zero', kernel_name), ...
        max(abs(hs0_rows.W_eta_variance_discrete)), tol_roundoff, '<=');

    for ff = 1:numel(f_list_hz)
        f_hz = f_list_hz(ff);
        subset = sortrows(kernel_rows(kernel_rows.f_hz == f_hz, :), 'Hs_target');
        checks = local_add_check(checks, sprintf('%s_f%d_R_coh_nonincreasing_Hs', kernel_name, f_hz), ...
            max(diff(subset.abs_R_coh)), tol_monotone, '<=');
        checks = local_add_check(checks, sprintf('%s_f%d_E_coh_frac_nonincreasing_Hs', kernel_name, f_hz), ...
            max(diff(subset.E_coh_over_E_inc)), tol_monotone, '<=');
        checks = local_add_check(checks, sprintf('%s_f%d_E_sca_limit_frac_nondecreasing_Hs', kernel_name, f_hz), ...
            min(diff(subset.E_sca_limit_over_E_inc)), -tol_monotone, '>=');
        checks = local_add_check(checks, sprintf('%s_f%d_W_eta_matches_sigma2', kernel_name, f_hz), ...
            max(subset.W_eta_variance_rel_error), tol_w_eta_rel, '<=');
    end

    for hh = 2:numel(Hs_list)
        Hs_target = Hs_list(hh);
        subset = sortrows(kernel_rows(kernel_rows.Hs_target == Hs_target, :), 'f_hz');
        checks = local_add_check(checks, sprintf('%s_Hs%g_R_coh_nonincreasing_frequency', kernel_name, Hs_target), ...
            max(diff(subset.abs_R_coh)), tol_monotone, '<=');
        checks = local_add_check(checks, sprintf('%s_Hs%g_exponent_more_negative_frequency', kernel_name, Hs_target), ...
            max(diff(subset.coherent_exponent)), tol_monotone, '<=');
        checks = local_add_check(checks, sprintf('%s_Hs%g_energy_conservation_frequency', kernel_name, Hs_target), ...
            max(subset.energy_conservation_error), tol_energy, '<=');
    end
end

ssa1_rows = summary_table(summary_table.kernel_mode == "ssa1_geometry", :);
pm_rows = summary_table(summary_table.kernel_mode == "pm_convolution", :);
checks = local_add_check(checks, 'ssa1_G_formula_recorded_all', ...
    double(all(contains(ssa1_rows.G_SSA1_formula, '4*gamma'))), 1, '==');
checks = local_add_check(checks, 'ssa1_boundary_condition_dirichlet_all', ...
    double(all(ssa1_rows.boundary_condition == "pressure-release / Dirichlet")), 1, '==');
checks = local_add_check(checks, 'ssa1_formula_source_recorded_all', ...
    double(all(contains(ssa1_rows.formula_source, 'vertical_comm_guide.md'))), 1, '==');
checks = local_add_check(checks, 'pm_convolution_engineering_baseline_all', ...
    double(all(contains(lower(pm_rows.formula_source), 'engineering') & ...
               contains(lower(pm_rows.limitations), 'engineering'))), 1, '==');
checks = local_add_check(checks, 'pm_convolution_not_strict_ssa_all', ...
    double(all(pm_rows.G_SSA1_formula == "not_applicable")), 1, '==');

bad_boundary = local_base_params();
bad_boundary.surface_boundary_model = 'ssa_stat_kernel';
bad_boundary.surface_ssa_kernel_mode = 'ssa1_geometry';
bad_boundary.surface_ssa_geometry_source_id = geometry_source_id;
bad_boundary.surface_ssa_random_scatter = false;
bad_boundary.surface_reflect_coeff = -0.8;
[bad_boundary_error_id, bad_boundary_error_message] = local_expect_error(@() vertical_channel_model(bad_boundary));
checks = local_add_check(checks, 'ssa1_non_dirichlet_rejected', ...
    double(strcmp(bad_boundary_error_id, 'pm_surface_boundary_model:SsaDirichletReflectCoeffRequired')), 1, '==');

trend_table = struct2table(checks);
validation_report = struct();
validation_report.script = mfilename;
validation_report.created_at = char(datetime('now'));
validation_report.Hs_list = Hs_list;
validation_report.f_list_hz = f_list_hz;
validation_report.kernel_modes = string(kernel_modes);
validation_report.conv_padding = 'periodic';
validation_report.random_scatter_enabled = false;
validation_report.tolerances = struct( ...
    'roundoff', tol_roundoff, ...
    'W_eta_variance_rel_error', tol_w_eta_rel, ...
    'monotone_slack', tol_monotone, ...
    'energy_conservation_error', tol_energy);
validation_report.non_dirichlet_error_id = bad_boundary_error_id;
validation_report.non_dirichlet_error_message = bad_boundary_error_message;
validation_report.all_passed = all(trend_table.passed);
validation_report.notes = ['Coherent reflection uses R0*exp(-2*k0^2*sigma_eta^2) ', ...
    'for normal specular pressure-release reflection. E_sca is zero because ', ...
    'surface_ssa_random_scatter=false; E_sca_limited records the scatter budget.'];

local_plot_trends(summary_table);

save(result_file, 'summary_table', 'trend_table', 'validation_report');

fprintf('\nSaved %s\n', result_file);
fprintf('All physical trend checks passed: %d\n', validation_report.all_passed);
if ~validation_report.all_passed
    disp(trend_table(~trend_table.passed, :))
    error('validate_ssa_physical_trends_vertical:Failed', ...
        'One or more SSA physical trend checks failed.');
end

function paramsV = local_base_params()
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

function row = local_summary_row(run_index, paramsV, channel, meta)
k0 = 2*pi*paramsV.f0/paramsV.c0;
R_expected = paramsV.surface_reflect_coeff * exp(-2*k0^2*meta.sigma_eta_m^2);
row = struct();
row.run_index = run_index;
row.kernel_mode = string(meta.kernel_mode);
row.conv_padding = string(meta.conv_padding);
row.f_hz = paramsV.f0;
row.Hs_target = paramsV.sea_hs_target;
row.sigma_eta = meta.sigma_eta_m;
row.k0 = k0;
row.R0 = paramsV.surface_reflect_coeff;
row.R_coh = meta.R_coh;
row.abs_R_coh = abs(meta.R_coh);
row.R_expected_normal = R_expected;
row.R_coh_formula_abs_error = abs(meta.R_coh - R_expected);
row.coherent_exponent = meta.coherent_exponent;
row.coherent_gamma_sum_eff_rad_per_m = meta.coherent_gamma_sum_eff_rad_per_m;
row.coherent_vertical_factor_eff = meta.coherent_vertical_factor_eff;
row.W_eta_variance_target = meta.W_eta_variance_target;
row.W_eta_variance_discrete = meta.W_eta_variance_discrete;
row.W_eta_variance_rel_error = meta.W_eta_variance_rel_error;
row.E_inc = meta.E_inc;
row.E_coh = meta.E_coh;
row.E_sca_raw = meta.E_sca_raw;
row.E_sca_limited = meta.E_sca_limited;
row.E_sca = meta.E_sca;
row.E_sca_limit = meta.E_sca_limit;
row.E_ref = meta.E_ref;
row.E_coh_over_E_inc = meta.E_coh / max(meta.E_inc, eps);
row.E_sca_raw_over_E_inc = meta.E_sca_raw / max(meta.E_inc, eps);
row.E_sca_limited_over_E_inc = meta.E_sca_limited / max(meta.E_inc, eps);
row.E_sca_limit_over_E_inc = meta.E_sca_limit / max(meta.E_inc, eps);
row.E_ref_over_E_inc = meta.E_ref / max(meta.E_inc, eps);
row.energy_conservation_error = meta.energy_conservation_error;
row.propagating_bin_fraction = meta.propagating_bin_fraction;
row.kernel_formula = string(meta.kernel_formula);
row.G_SSA1_formula = string(meta.kernel_detail.G_SSA1_formula);
row.boundary_condition = string(meta.boundary_condition);
row.formula_source = string(meta.formula_source);
row.limitations = string(meta.limitations);
row.random_scatter_enabled = logical(meta.random_scatter_enabled);
row.invariant_error = max(abs(channel.H_f(:) - channel.H_direct_f(:) - channel.H_reflect_f(:)));
row.abs_h_reflect = abs(channel.h_reflect);
row.abs_h_total = abs(channel.h_total);
end

function rows = local_append_row(rows, row)
if isempty(rows)
    rows = row;
else
    rows(end + 1) = row; %#ok<AGROW>
end
end

function checks = local_add_check(checks, check_name, value, tolerance, comparison)
row = struct();
row.check_name = string(check_name);
row.value = double(value);
row.tolerance = double(tolerance);
row.comparison = string(comparison);
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
if isempty(checks)
    checks = row;
else
    checks(end + 1) = row; %#ok<AGROW>
end
end

function [err_id, err_message] = local_expect_error(fn)
try
    fn();
    err_id = '';
    err_message = '';
catch ME
    err_id = ME.identifier;
    err_message = ME.message;
end
end

function local_plot_trends(summary_table)
primary = summary_table(summary_table.kernel_mode == "ssa1_geometry", :);

fig = figure('Visible', 'off');
hold on
f_values = unique(primary.f_hz).';
for f_hz = f_values
    subset = sortrows(primary(primary.f_hz == f_hz, :), 'Hs_target');
    plot(subset.Hs_target, subset.abs_R_coh, '-o', 'DisplayName', sprintf('%g Hz', f_hz));
end
grid on
xlabel('H_s (m)')
ylabel('|R_{coh}|')
title('|R_{coh}| vs H_s, ssa1\_geometry')
legend('Location', 'best')
saveas(fig, 'validate_ssa_abs_R_coh_vs_Hs.png')
close(fig)

fig = figure('Visible', 'off');
hold on
Hs_values = unique(primary.Hs_target).';
for Hs = Hs_values
    subset = sortrows(primary(primary.Hs_target == Hs, :), 'f_hz');
    plot(subset.f_hz, subset.abs_R_coh, '-o', 'DisplayName', sprintf('H_s=%g m', Hs));
end
grid on
xlabel('Frequency (Hz)')
ylabel('|R_{coh}|')
title('|R_{coh}| vs frequency, ssa1\_geometry')
legend('Location', 'best')
saveas(fig, 'validate_ssa_abs_R_coh_vs_frequency.png')
close(fig)

fig = figure('Visible', 'off');
hold on
for f_hz = f_values
    subset = sortrows(primary(primary.f_hz == f_hz, :), 'Hs_target');
    plot(subset.Hs_target, subset.E_coh_over_E_inc, '-o', 'DisplayName', sprintf('%g Hz', f_hz));
end
grid on
xlabel('H_s (m)')
ylabel('E_{coh}/E_{inc}')
title('Coherent energy fraction vs H_s, ssa1\_geometry')
legend('Location', 'best')
saveas(fig, 'validate_ssa_coherent_energy_fraction_vs_Hs.png')
close(fig)

fig = figure('Visible', 'off');
hold on
for f_hz = f_values
    subset = sortrows(primary(primary.f_hz == f_hz, :), 'Hs_target');
    plot(subset.Hs_target, subset.E_sca_limit_over_E_inc, '-o', 'DisplayName', sprintf('%g Hz', f_hz));
end
grid on
xlabel('H_s (m)')
ylabel('E_{sca}^{max}/E_{inc}')
title('Scatter budget fraction vs H_s, ssa1\_geometry')
legend('Location', 'best')
saveas(fig, 'validate_ssa_scatter_budget_fraction_vs_Hs.png')
close(fig)

fig = figure('Visible', 'off');
plot(primary.W_eta_variance_target, primary.W_eta_variance_discrete, 'o')
hold on
lims = [min([primary.W_eta_variance_target; primary.W_eta_variance_discrete]), ...
        max([primary.W_eta_variance_target; primary.W_eta_variance_discrete])];
if lims(2) <= lims(1)
    lims = [0, 1];
end
plot(lims, lims, 'k--', 'DisplayName', 'target=discrete')
grid on
xlabel('\sigma_\eta^2 target')
ylabel('sum(W_\eta)\DeltaK_x\DeltaK_y')
title('W_\eta variance normalization check')
legend('Location', 'best')
saveas(fig, 'validate_ssa_W_eta_variance_check.png')
close(fig)
end

