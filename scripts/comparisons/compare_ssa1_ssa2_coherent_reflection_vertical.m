run(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'bootstrap_project.m'));
% Compare SSA1 and Broschat-style SSA2 coherent reflection coefficients.
% This reduced-grid script does not run the communication chain. It only
% exercises the ssa_stat_kernel coherent reflection branch and its metadata.

clear
format compact

result_file = getenv('SSA1_SSA2_COHERENT_RESULT_FILE');
if isempty(result_file)
    result_file = 'compare_ssa1_ssa2_coherent_reflection_vertical_result.mat';
end

Hs_list = local_env_numeric_vector('SSA1_SSA2_COHERENT_HS_LIST', [0, 0.05, 0.2, 0.5, 1.0]);
f_list_hz = local_env_numeric_vector('SSA1_SSA2_COHERENT_F_LIST_HZ', [4000, 6000, 8000, 10000]);
grid_n = local_env_scalar_int('SSA1_SSA2_COHERENT_GRID_N', 128);
threshold_db = local_env_scalar('SSA1_SSA2_COHERENT_THRESHOLD_DB', 0.1);

orders = {'ssa1', 'ssa2_broschat_coherent'};
tol_roundoff = 1e-12;
tol_w_eta_rel = 1e-10;
tol_energy = 1e-12;

run_rows = struct([]);
summary_rows = struct([]);
checks = struct([]);
run_index = 0;
summary_index = 0;

for ff = 1:numel(f_list_hz)
    for hh = 1:numel(Hs_list)
        order_rows = struct([]);
        for oo = 1:numel(orders)
            paramsV = local_base_params(grid_n);
            paramsV.f0 = f_list_hz(ff);
            paramsV.f_ref_hz = paramsV.f0;
            paramsV.sea_hs_target = Hs_list(hh);
            paramsV.surface_ssa_coherent_order = orders{oo};

            fprintf('SSA1/SSA2 coherent compare: order=%s, f=%g Hz, Hs=%g m.\n', ...
                orders{oo}, paramsV.f0, paramsV.sea_hs_target);
            channel = vertical_channel_model(paramsV);
            meta = channel.roughness_meta.ssa_stat_kernel_meta;

            run_index = run_index + 1;
            row = local_run_row(run_index, paramsV, channel, meta);
            run_rows = local_append_struct(run_rows, row);
            order_rows = local_append_struct(order_rows, row);
        end

        summary_index = summary_index + 1;
        summary_rows = local_append_struct(summary_rows, ...
            local_summary_row(summary_index, struct2table(order_rows)));
    end
end

run_table = struct2table(run_rows);
summary_table = struct2table(summary_rows);

checks = local_add_check(checks, 'channel_invariant_all', ...
    max(run_table.invariant_error), 1e-10, '<=');
checks = local_add_check(checks, 'energy_conservation_all', ...
    max(run_table.energy_conservation_error), tol_energy, '<=');
checks = local_add_check(checks, 'W_eta_variance_rel_error_all', ...
    max(summary_table.W_eta_variance_rel_error), tol_w_eta_rel, '<=');
checks = local_add_check(checks, 'Hs0_R1_equals_R0', ...
    max(abs(summary_table.R_coh_ssa1(summary_table.Hs_target == 0) - summary_table.R0(summary_table.Hs_target == 0))), ...
    tol_roundoff, '<=');
checks = local_add_check(checks, 'Hs0_R2_equals_R0', ...
    max(abs(summary_table.R_coh_ssa2(summary_table.Hs_target == 0) - summary_table.R0(summary_table.Hs_target == 0))), ...
    tol_roundoff, '<=');
checks = local_add_check(checks, 'Hs0_delta_R_zero', ...
    max(abs(summary_table.delta_R_abs(summary_table.Hs_target == 0))), tol_roundoff, '<=');
checks = local_add_check(checks, 'normal_ssa1_formula_all', ...
    max(summary_table.R_coh_ssa1_formula_abs_error), tol_roundoff, '<=');
checks = local_add_check(checks, 'ssa2_finite_all', ...
    double(all(isfinite(real(summary_table.R_coh_ssa2)) & isfinite(imag(summary_table.R_coh_ssa2)))), 1, '==');

bad_params = local_base_params(grid_n);
bad_params.surface_ssa_coherent_order = 'ssa2_broschat_coherent';
bad_params.surface_reflect_coeff = -0.8;
[bad_error_id, bad_error_message] = local_expect_error(@() vertical_channel_model(bad_params));
checks = local_add_check(checks, 'ssa2_non_dirichlet_rejected', ...
    double(strcmp(bad_error_id, 'pm_surface_boundary_model:Ssa2DirichletReflectCoeffRequired')), 1, '==');

negligible_mask = abs(summary_table.coherent_loss_delta_db) <= threshold_db;
checks = local_add_check(checks, 'ssa2_loss_delta_finite_all', ...
    double(all(isfinite(summary_table.coherent_loss_delta_db))), 1, '==');

trend_table = struct2table(checks);
validation_report = struct();
validation_report.script = mfilename;
validation_report.created_at = char(datetime('now'));
validation_report.Hs_list = Hs_list;
validation_report.f_list_hz = f_list_hz;
validation_report.grid_n = grid_n;
validation_report.threshold_db = threshold_db;
validation_report.all_passed = all(trend_table.passed);
validation_report.ssa2_negligible_all = all(negligible_mask);
validation_report.ssa2_negligible_count = nnz(negligible_mask);
validation_report.ssa2_total_count = numel(negligible_mask);
validation_report.max_abs_coherent_loss_delta_db = max(abs(summary_table.coherent_loss_delta_db));
validation_report.non_negligible_rows = summary_table(~negligible_mask, ...
    {'Hs_target', 'f_hz', 'coherent_loss_delta_db', 'ssa2_evanescent_fraction_abs'});
validation_report.non_dirichlet_error_id = bad_error_id;
validation_report.non_dirichlet_error_message = bad_error_message;
if validation_report.ssa2_negligible_all
    validation_report.interpretation = ...
        'SSA2 coherent correction is negligible under the current near-vertical setting.';
else
    validation_report.interpretation = ...
        'SSA2 coherent correction exceeds the threshold for at least one Hs/frequency point.';
end
validation_report.notes = ['ssa2_broschat_coherent replaces only R_coh. ', ...
    'The incoherent P_sca branch remains the first-order Dirichlet SSA1 geometry kernel.'];

local_plot_all(summary_table, Hs_list, f_list_hz, threshold_db);

save(result_file, 'run_table', 'summary_table', 'trend_table', ...
    'validation_report', 'Hs_list', 'f_list_hz', 'threshold_db', 'grid_n');

disp(summary_table(:, {'Hs_target', 'f_hz', 'abs_R_coh_ssa1', 'abs_R_coh_ssa2', ...
    'coherent_loss_ssa1_db', 'coherent_loss_ssa2_db', ...
    'coherent_loss_delta_db', 'ssa2_evanescent_fraction_abs'}))
disp(trend_table)
fprintf('Saved %s\n', result_file);
fprintf('All checks passed: %d\n', validation_report.all_passed);
fprintf('Max |SSA2-SSA1 coherent loss delta| = %.6g dB\n', ...
    validation_report.max_abs_coherent_loss_delta_db);
if ~validation_report.all_passed
    error('compare_ssa1_ssa2_coherent_reflection_vertical:Failed', ...
        'One or more SSA1/SSA2 coherent reflection checks failed.');
end

function paramsV = local_base_params(grid_n)
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
paramsV.surface_phase_mode = 'normal';
paramsV.surface_boundary_model = 'ssa_stat_kernel';
paramsV.surface_boundary_coupling_diagnostics = false;
paramsV.surface_boundary_redistribution_diagnostics = false;
paramsV.surface_ssa_random_scatter = false;
paramsV.surface_ssa_scatter_scale = 1.0;
paramsV.surface_ssa_seed_offset = 100000;
paramsV.surface_ssa_kernel_mode = 'ssa1_geometry';
paramsV.surface_ssa_coherent_order = 'ssa1';
paramsV.surface_ssa_geometry_source_id = ['SSA.md; Thorsos & Broschat 1995 JASA, ', ...
    'pressure-release / Dirichlet first-order perturbation-limit geometry'];
paramsV.surface_ssa_kz_branch = 'downward_positive_real';
paramsV.surface_ssa_conv_padding = 'periodic';
end

function row = local_run_row(run_index, paramsV, channel, meta)
k0 = 2*pi*paramsV.f0/paramsV.c0;
R_expected = paramsV.surface_reflect_coeff * exp(-2*k0^2*meta.sigma_eta_m^2);
row = struct();
row.run_index = run_index;
row.Hs_target = paramsV.sea_hs_target;
row.f_hz = paramsV.f0;
row.k0 = k0;
row.sigma_eta = meta.sigma_eta_m;
row.R0 = paramsV.surface_reflect_coeff;
row.coherent_order = string(meta.coherent_order);
row.R_coh_selected = meta.R_coh;
row.R_coh_raw_selected = meta.R_coh_raw;
row.R_coh_ssa1 = meta.R_coh_ssa1;
row.R_coh_ssa2 = meta.R_coh_ssa2;
row.abs_R_coh_ssa1 = meta.abs_R_coh_ssa1;
row.abs_R_coh_ssa2 = meta.abs_R_coh_ssa2;
row.R_expected_ssa1_normal = R_expected;
row.R_coh_ssa1_formula_abs_error = abs(meta.R_coh_ssa1 - R_expected);
row.coherent_loss_ssa1_db = meta.coherent_loss_ssa1_db;
row.coherent_loss_ssa2_db = meta.coherent_loss_ssa2_db;
row.coherent_loss_delta_db = meta.coherent_loss_delta_db;
row.delta_R_abs = meta.delta_R_abs;
row.delta_R_rel = meta.delta_R_rel;
row.ssa2_correction_integral = meta.ssa2_correction_integral;
row.ssa2_correction_integral_real = meta.ssa2_correction_integral_real;
row.ssa2_correction_integral_imag = meta.ssa2_correction_integral_imag;
row.ssa2_gamma_i_eff_rad_per_m = meta.ssa2_gamma_i_eff_rad_per_m;
row.ssa2_K_i_eff_kx_rad_per_m = meta.ssa2_K_i_eff_rad_per_m(1);
row.ssa2_K_i_eff_ky_rad_per_m = meta.ssa2_K_i_eff_rad_per_m(2);
row.ssa2_propagating_integral = meta.ssa2_propagating_integral;
row.ssa2_evanescent_integral = meta.ssa2_evanescent_integral;
row.ssa2_evanescent_fraction_abs = meta.ssa2_evanescent_fraction_abs;
row.ssa2_sqrt_branch = string(meta.ssa2_sqrt_branch);
row.ssa2_formula = string(meta.ssa2_formula);
row.ssa2_formula_source = string(meta.ssa2_formula_source);
row.W_eta_variance_target = meta.W_eta_variance_target;
row.W_eta_variance_discrete = meta.W_eta_variance_discrete;
row.W_eta_variance_rel_error = meta.W_eta_variance_rel_error;
row.E_inc = meta.E_inc;
row.E_coh = meta.E_coh;
row.E_sca_raw = meta.E_sca_raw;
row.E_sca_limited = meta.E_sca_limited;
row.E_sca = meta.E_sca;
row.E_ref = meta.E_ref;
row.energy_conservation_error = meta.energy_conservation_error;
row.invariant_error = max(abs(channel.H_f(:) - channel.H_direct_f(:) - channel.H_reflect_f(:)));
end

function row = local_summary_row(summary_index, order_table)
ssa2_rows = order_table(order_table.coherent_order == "ssa2_broschat_coherent", :);
if isempty(ssa2_rows)
    source = order_table(1, :);
else
    source = ssa2_rows(1, :);
end
row = struct();
row.summary_index = summary_index;
row.Hs_target = source.Hs_target;
row.f_hz = source.f_hz;
row.k0 = source.k0;
row.sigma_eta = source.sigma_eta;
row.R0 = source.R0;
row.R_coh_ssa1 = source.R_coh_ssa1;
row.R_coh_ssa2 = source.R_coh_ssa2;
row.abs_R_coh_ssa1 = source.abs_R_coh_ssa1;
row.abs_R_coh_ssa2 = source.abs_R_coh_ssa2;
row.R_expected_ssa1_normal = source.R_expected_ssa1_normal;
row.R_coh_ssa1_formula_abs_error = source.R_coh_ssa1_formula_abs_error;
row.coherent_loss_ssa1_db = source.coherent_loss_ssa1_db;
row.coherent_loss_ssa2_db = source.coherent_loss_ssa2_db;
row.coherent_loss_delta_db = source.coherent_loss_delta_db;
row.delta_R_abs = source.delta_R_abs;
row.delta_R_rel = source.delta_R_rel;
row.ssa2_correction_integral = source.ssa2_correction_integral;
row.ssa2_correction_integral_real = source.ssa2_correction_integral_real;
row.ssa2_correction_integral_imag = source.ssa2_correction_integral_imag;
row.ssa2_evanescent_fraction_abs = source.ssa2_evanescent_fraction_abs;
row.energy_conservation_error = max(order_table.energy_conservation_error);
row.W_eta_variance_target = source.W_eta_variance_target;
row.W_eta_variance_discrete = source.W_eta_variance_discrete;
row.W_eta_variance_rel_error = source.W_eta_variance_rel_error;
end

function local_plot_all(summary_table, Hs_list, f_list_hz, threshold_db)
local_plot_vs_hs(summary_table, f_list_hz, ...
    'abs_R_coh_ssa1', 'abs_R_coh_ssa2', '|R_{coh}|', ...
    'SSA1 vs SSA2 coherent magnitude', 'compare_ssa1_ssa2_abs_R_coh_vs_Hs.png');
local_plot_vs_hs(summary_table, f_list_hz, ...
    'coherent_loss_ssa1_db', 'coherent_loss_ssa2_db', 'Coherent loss (dB)', ...
    'SSA1 vs SSA2 coherent loss', 'compare_ssa1_ssa2_coherent_loss_vs_Hs.png');
local_plot_single_vs_hs(summary_table, f_list_hz, ...
    'coherent_loss_delta_db', 'SSA2-SSA1 loss delta (dB)', ...
    'SSA2 coherent loss delta vs Hs', 'compare_ssa1_ssa2_loss_delta_vs_Hs.png');
local_plot_single_vs_frequency(summary_table, Hs_list, ...
    'coherent_loss_delta_db', 'SSA2-SSA1 loss delta (dB)', ...
    'SSA2 coherent loss delta vs frequency', 'compare_ssa1_ssa2_loss_delta_vs_frequency.png');
local_plot_heatmap(summary_table, Hs_list, f_list_hz, ...
    double(abs(summary_table.coherent_loss_delta_db) <= threshold_db), ...
    '1 means |loss delta| <= threshold', ...
    'SSA2 negligible region', 'compare_ssa1_ssa2_negligible_region_heatmap.png');
local_plot_heatmap(summary_table, Hs_list, f_list_hz, ...
    summary_table.ssa2_evanescent_fraction_abs, ...
    'Evanescent integral fraction', ...
    'SSA2 evanescent fraction', 'compare_ssa1_ssa2_evanescent_fraction_heatmap.png');
end

function local_plot_vs_hs(T, f_list_hz, field1, field2, y_label, plot_title, file_name)
fig = figure('Visible', 'off');
tiledlayout(2, 2, 'TileSpacing', 'compact')
for ii = 1:numel(f_list_hz)
    nexttile
    subset = sortrows(T(T.f_hz == f_list_hz(ii), :), 'Hs_target');
    plot(subset.Hs_target, subset.(field1), '-o', 'LineWidth', 1.2)
    hold on
    plot(subset.Hs_target, subset.(field2), '-s', 'LineWidth', 1.2)
    grid on
    xlabel('H_s (m)')
    ylabel(y_label)
    title(sprintf('%g Hz', f_list_hz(ii)))
    legend('SSA1', 'SSA2', 'Location', 'best')
end
sgtitle(plot_title)
saveas(fig, file_name)
close(fig)
end

function local_plot_single_vs_hs(T, f_list_hz, field_name, y_label, plot_title, file_name)
fig = figure('Visible', 'off');
hold on
for ii = 1:numel(f_list_hz)
    subset = sortrows(T(T.f_hz == f_list_hz(ii), :), 'Hs_target');
    plot(subset.Hs_target, subset.(field_name), '-o', 'LineWidth', 1.2, ...
        'DisplayName', sprintf('%g Hz', f_list_hz(ii)))
end
grid on
xlabel('H_s (m)')
ylabel(y_label)
title(plot_title)
legend('Location', 'best')
saveas(fig, file_name)
close(fig)
end

function local_plot_single_vs_frequency(T, Hs_list, field_name, y_label, plot_title, file_name)
fig = figure('Visible', 'off');
hold on
for ii = 1:numel(Hs_list)
    subset = sortrows(T(T.Hs_target == Hs_list(ii), :), 'f_hz');
    plot(subset.f_hz, subset.(field_name), '-o', 'LineWidth', 1.2, ...
        'DisplayName', sprintf('H_s=%g m', Hs_list(ii)))
end
grid on
xlabel('Frequency (Hz)')
ylabel(y_label)
title(plot_title)
legend('Location', 'best')
saveas(fig, file_name)
close(fig)
end

function local_plot_heatmap(T, Hs_list, f_list_hz, values, color_label, plot_title, file_name)
Z = NaN(numel(Hs_list), numel(f_list_hz));
for ii = 1:numel(Hs_list)
    for jj = 1:numel(f_list_hz)
        idx = T.Hs_target == Hs_list(ii) & T.f_hz == f_list_hz(jj);
        if any(idx)
            Z(ii, jj) = values(find(idx, 1, 'first'));
        end
    end
end
fig = figure('Visible', 'off');
imagesc(f_list_hz, Hs_list, Z)
set(gca, 'YDir', 'normal')
cb = colorbar;
xlabel('Frequency (Hz)')
ylabel('H_s (m)')
title(plot_title)
ylabel(cb, color_label)
saveas(fig, file_name)
close(fig)
end

function rows = local_append_struct(rows, row)
if isempty(rows)
    rows = row;
else
    rows(end+1) = row; %#ok<AGROW>
end
end

function checks = local_add_check(checks, name, value, threshold, comparator)
passed = false;
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

function [err_id, err_msg] = local_expect_error(fn)
err_id = '';
err_msg = '';
try
    fn();
catch ME
    err_id = ME.identifier;
    err_msg = ME.message;
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

function value = local_env_scalar(name, default_value)
raw = getenv(name);
if isempty(raw)
    value = default_value;
else
    value = str2double(raw);
    if ~isfinite(value)
        error('%s must be a finite scalar.', name);
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

