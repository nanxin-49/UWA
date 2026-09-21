run(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'bootstrap_project.m'));
% Compare specular coherent and incoherent SSA1 surface-reflection components.
%
% This diagnostic can either keep Hs_target fixed while sweeping wind speed
% and frequency, or scan an explicit Hs list. At each Hs, wind changes the
% normalized PM spectral shape, not the total roughness variance. It does not
% run the communication chain.

clear
format compact

result_file = getenv('SPEC_INCOH_RESULT_FILE');
if isempty(result_file)
    result_file = 'compare_specular_incoherent_surface_reflection_vertical_result.mat';
end

wind_list = local_env_numeric_vector('SPEC_INCOH_WIND_LIST', [3, 5, 8, 12]);
f_list_hz = local_env_numeric_vector('SPEC_INCOH_F_LIST_HZ', [4000, 6000, 8000, 10000]);
Hs_list = local_env_hs_list();
seed_count = local_env_scalar_int('SPEC_INCOH_SEED_COUNT', 8);
grid_n = local_env_scalar_int('SPEC_INCOH_GRID_N', 128);
scatter_scale = local_env_scalar('SPEC_INCOH_SCATTER_SCALE', 1.0);
conv_padding = getenv('SPEC_INCOH_CONV_PADDING');
if isempty(conv_padding)
    conv_padding = 'periodic';
end
seed_list = 12345 + (0:(seed_count - 1));

tol_invariant = 1e-10;
tol_component = 1e-10;
tol_energy = 1e-12;
tol_w_eta = 1e-10;
threshold_negligible_db = -20;
threshold_not_negligible_db = -10;
threshold_std_large_db = 6;

run_rows = struct([]);
run_index = 0;

fprintf('Running diagnostics-disabled equivalence check.\n');
disabled_params = local_base_params(grid_n, Hs_list(1), wind_list(1), scatter_scale, conv_padding);
disabled_params.f0 = f_list_hz(1);
disabled_params.f_ref_hz = disabled_params.f0;
disabled_params.sea_seed = seed_list(1);
enabled_params = disabled_params;
enabled_params.surface_ssa_component_diagnostics = true;
channel_disabled = vertical_channel_model(disabled_params);
channel_enabled = vertical_channel_model(enabled_params);
disabled_equivalence = struct( ...
    'H_f_max_abs_diff', max(abs(channel_disabled.H_f(:) - channel_enabled.H_f(:))), ...
    'H_direct_f_max_abs_diff', max(abs(channel_disabled.H_direct_f(:) - channel_enabled.H_direct_f(:))), ...
    'H_reflect_f_max_abs_diff', max(abs(channel_disabled.H_reflect_f(:) - channel_enabled.H_reflect_f(:))));

for ih = 1:numel(Hs_list)
    for iw = 1:numel(wind_list)
        for ff = 1:numel(f_list_hz)
            for ss = 1:numel(seed_list)
                paramsV = local_base_params(grid_n, Hs_list(ih), wind_list(iw), scatter_scale, conv_padding);
                paramsV.f0 = f_list_hz(ff);
                paramsV.f_ref_hz = paramsV.f0;
                paramsV.sea_seed = seed_list(ss);

                fprintf('Specular/incoherent case: Hs=%g m, U=%g m/s, f=%g Hz, seed=%d.\n', ...
                    paramsV.sea_hs_target, paramsV.sea_wind_speed, ...
                    paramsV.f0, paramsV.sea_seed);
                channel = vertical_channel_model(paramsV);
                run_index = run_index + 1;
                run_rows = local_append_struct(run_rows, local_run_row(run_index, paramsV, channel));
            end
        end
    end
end

run_table = struct2table(run_rows);
summary_table = local_build_summary_table(run_table, Hs_list, wind_list, f_list_hz);
decision_table = local_build_decision_table( ...
    summary_table, threshold_negligible_db, threshold_not_negligible_db, threshold_std_large_db);
validation_report = local_build_validation_report( ...
    run_table, disabled_equivalence, tol_invariant, tol_component, tol_energy, tol_w_eta);

local_plot_results(summary_table, decision_table, Hs_list, wind_list, f_list_hz);

summary_csv = 'compare_specular_incoherent_surface_reflection_vertical_summary.csv';
writetable(summary_table, summary_csv);
save(result_file, 'run_table', 'summary_table', 'decision_table', 'validation_report', ...
    'wind_list', 'f_list_hz', 'Hs_list', 'seed_list', 'grid_n', 'scatter_scale', ...
    'conv_padding', 'threshold_negligible_db', 'threshold_not_negligible_db', ...
    'threshold_std_large_db', 'disabled_equivalence');

disp(decision_table)
fprintf('Saved %s\n', result_file);
fprintf('Saved %s\n', summary_csv);
fprintf('All checks passed: %d\n', validation_report.all_passed);
if ~validation_report.all_passed
    error('compare_specular_incoherent_surface_reflection_vertical:Failed', ...
        'One or more specular/incoherent diagnostics checks failed.');
end

function paramsV = local_base_params(grid_n, Hs_target, wind_speed, scatter_scale, conv_padding)
domain_width_m = min(50, 0.5 * grid_n);
paramsV = struct();
paramsV.f0 = 6000;
paramsV.c0 = 1500;
paramsV.enable_wideband = false;
paramsV.f_ref_hz = paramsV.f0;
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
paramsV.surface_reflect_coeff = -1;
paramsV.surface_phase_mode = 'normal';
paramsV.surface_boundary_model = 'ssa_stat_kernel';
paramsV.surface_ssa_kernel_mode = 'ssa1_geometry';
paramsV.surface_ssa_coherent_order = 'ssa1';
paramsV.surface_ssa_random_scatter = true;
paramsV.surface_ssa_scatter_scale = scatter_scale;
paramsV.surface_ssa_seed_offset = 100000;
paramsV.surface_ssa_geometry_source_id = ...
    'vertical_comm_guide.md; first-order pressure-release / Dirichlet geometry';
paramsV.surface_ssa_kz_branch = 'downward_positive_real';
paramsV.surface_ssa_conv_padding = char(conv_padding);
paramsV.surface_ssa_frequency_correlation_mode = 'independent';
paramsV.surface_ssa_component_diagnostics = true;
paramsV.sea_wind_speed = wind_speed;
paramsV.sea_hs_target = Hs_target;
paramsV.sea_seed = 12345;
end

function row = local_run_row(run_index, paramsV, channel)
ssa = channel.roughness_meta.ssa_stat_kernel_meta;
component = channel.surface_ssa_component_meta;
idx = find(component.recorded_frequency_mask, 1, 'first');
if isempty(idx)
    idx = 1;
end

row = struct();
row.run_index = run_index;
row.wind_speed_mps = paramsV.sea_wind_speed;
row.f_hz = paramsV.f0;
row.Hs_target_m = paramsV.sea_hs_target;
row.seed = paramsV.sea_seed;
row.grid_n = paramsV.nx;
row.scatter_scale = paramsV.surface_ssa_scatter_scale;
row.kernel_mode = string(ssa.kernel_mode);
row.conv_padding = string(ssa.conv_padding);
row.random_scatter_enabled = logical(ssa.random_scatter_enabled);
row.R_coh = ssa.R_coh;
row.abs_R_coh = abs(ssa.R_coh);
row.R_coh_power = abs(ssa.R_coh)^2;
row.h_direct = channel.h_direct;
row.h_reflect_total = component.h_reflect_total_f(idx);
row.h_reflect_coh = component.h_reflect_coh_f(idx);
row.h_reflect_sca = component.h_reflect_sca_f(idx);
row.h_total = channel.h_total;
row.abs_h_direct = abs(channel.h_direct);
row.abs_h_reflect_total = abs(component.h_reflect_total_f(idx));
row.abs_h_reflect_coh = abs(component.h_reflect_coh_f(idx));
row.abs_h_reflect_sca = abs(component.h_reflect_sca_f(idx));
row.abs_h_total = abs(channel.h_total);
row.phase_h_reflect_total = angle(component.h_reflect_total_f(idx));
row.phase_h_reflect_coh = angle(component.h_reflect_coh_f(idx));
row.phase_h_reflect_sca = angle(component.h_reflect_sca_f(idx));
row.component_phase_diff_rad = wrapToPiLocal(row.phase_h_reflect_sca - row.phase_h_reflect_coh);
row.scatter_to_coherent_rx_db = component.scatter_to_coherent_rx_db_f(idx);
row.scatter_to_total_rx_db = component.scatter_to_total_rx_db_f(idx);
row.coherent_to_total_rx_db = component.coherent_to_total_rx_db_f(idx);
row.scatter_to_coherent_rx_abs = row.abs_h_reflect_sca / max(row.abs_h_reflect_coh, eps);
row.E_inc = ssa.E_inc;
row.E_coh = ssa.E_coh;
row.E_sca_raw = ssa.E_sca_raw;
row.E_sca_limited = ssa.E_sca_limited;
row.E_sca = ssa.E_sca;
row.E_ref = ssa.E_ref;
row.E_coh_over_E_inc = ssa.E_coh / max(ssa.E_inc, eps);
row.E_sca_limited_over_E_inc = ssa.E_sca_limited / max(ssa.E_inc, eps);
row.E_sca_over_E_inc = ssa.E_sca / max(ssa.E_inc, eps);
row.E_ref_over_E_inc = ssa.E_ref / max(ssa.E_inc, eps);
row.scatter_energy_fraction_db = 10 * log10(max(ssa.E_sca_limited, eps) / max(ssa.E_coh, eps));
row.energy_conservation_error = ssa.energy_conservation_error;
row.energy_limit_applied = logical(ssa.energy_limit_applied);
row.energy_scale_applied = ssa.energy_scale_applied;
row.W_eta_variance_rel_error = ssa.W_eta_variance_rel_error;
row.component_sum_error_abs = component.component_sum_error_abs_f(idx);
row.component_sum_error_rel = component.component_sum_error_rel_f(idx);
row.channel_invariant_error = max(abs(channel.H_f(:) - channel.H_direct_f(:) - channel.H_reflect_f(:)));
row.reflected_rms_delta_k_rad_per_m = ssa.reflected_spectrum_stats.rms_delta_k_rad_per_m;
row.scatter_rms_delta_k_rad_per_m = ssa.scatter_power_spectrum_stats.rms_delta_k_rad_per_m;
end

function summary_table = local_build_summary_table(run_table, Hs_list, wind_list, f_list_hz)
summary_rows = struct([]);
idx_out = 0;
for ih = 1:numel(Hs_list)
    for iw = 1:numel(wind_list)
        for ff = 1:numel(f_list_hz)
            mask = run_table.Hs_target_m == Hs_list(ih) & ...
                run_table.wind_speed_mps == wind_list(iw) & ...
                run_table.f_hz == f_list_hz(ff);
            subset = run_table(mask, :);
            idx_out = idx_out + 1;
            row = struct();
            row.summary_index = idx_out;
            row.Hs_target_m = Hs_list(ih);
            row.wind_speed_mps = wind_list(iw);
            row.f_hz = f_list_hz(ff);
            row.seed_count = height(subset);
            row.abs_R_coh_mean = mean(subset.abs_R_coh);
            row.R_coh_power_mean = mean(subset.R_coh_power);
            row.E_coh_over_E_inc_mean = mean(subset.E_coh_over_E_inc);
            row.E_sca_limited_over_E_inc_mean = mean(subset.E_sca_limited_over_E_inc);
            row.E_sca_over_E_inc_mean = mean(subset.E_sca_over_E_inc);
            row.E_ref_over_E_inc_mean = mean(subset.E_ref_over_E_inc);
            row.scatter_energy_fraction_db_mean = mean(subset.scatter_energy_fraction_db);
            row.scatter_energy_fraction_db_std = std(subset.scatter_energy_fraction_db);
            row.abs_h_reflect_coh_mean = mean(subset.abs_h_reflect_coh);
            row.abs_h_reflect_coh_std = std(subset.abs_h_reflect_coh);
            row.abs_h_reflect_sca_mean = mean(subset.abs_h_reflect_sca);
            row.abs_h_reflect_sca_std = std(subset.abs_h_reflect_sca);
            row.abs_h_reflect_total_mean = mean(subset.abs_h_reflect_total);
            row.abs_h_reflect_total_std = std(subset.abs_h_reflect_total);
            row.abs_h_total_mean = mean(subset.abs_h_total);
            row.abs_h_total_std = std(subset.abs_h_total);
            row.scatter_to_coherent_rx_db_mean = mean(subset.scatter_to_coherent_rx_db);
            row.scatter_to_coherent_rx_db_std = std(subset.scatter_to_coherent_rx_db);
            row.scatter_to_total_rx_db_mean = mean(subset.scatter_to_total_rx_db);
            row.coherent_to_total_rx_db_mean = mean(subset.coherent_to_total_rx_db);
            row.phase_diff_circular_mean_rad = angle(mean(exp(1i * subset.component_phase_diff_rad)));
            row.phase_diff_circular_variance = 1 - abs(mean(exp(1i * subset.component_phase_diff_rad)));
            row.energy_limit_trigger_fraction = mean(double(subset.energy_limit_applied));
            row.energy_scale_applied_mean = mean(subset.energy_scale_applied);
            row.energy_conservation_error_max = max(subset.energy_conservation_error);
            row.W_eta_variance_rel_error_max = max(subset.W_eta_variance_rel_error);
            row.component_sum_error_abs_max = max(subset.component_sum_error_abs);
            row.component_sum_error_rel_max = max(subset.component_sum_error_rel);
            row.channel_invariant_error_max = max(subset.channel_invariant_error);
            row.reflected_rms_delta_k_mean = mean(subset.reflected_rms_delta_k_rad_per_m);
            row.scatter_rms_delta_k_mean = mean(subset.scatter_rms_delta_k_rad_per_m);
            summary_rows = local_append_struct(summary_rows, row);
        end
    end
end
summary_table = struct2table(summary_rows);
end

function decision_table = local_build_decision_table(summary_table, negligible_db, not_negligible_db, std_large_db)
decision_rows = struct([]);
for ii = 1:height(summary_table)
    rx_db = summary_table.scatter_to_coherent_rx_db_mean(ii);
    energy_db = summary_table.scatter_energy_fraction_db_mean(ii);
    ratio_std_db = summary_table.scatter_to_coherent_rx_db_std(ii);
    if rx_db <= negligible_db && energy_db <= negligible_db
        decision = "negligible";
    elseif rx_db > not_negligible_db || energy_db > not_negligible_db || ratio_std_db > std_large_db
        decision = "not_negligible";
    else
        decision = "borderline";
    end
    row = struct();
    row.Hs_target_m = summary_table.Hs_target_m(ii);
    row.wind_speed_mps = summary_table.wind_speed_mps(ii);
    row.f_hz = summary_table.f_hz(ii);
    row.scatter_to_coherent_rx_db_mean = rx_db;
    row.scatter_to_coherent_rx_db_std = ratio_std_db;
    row.scatter_energy_fraction_db_mean = energy_db;
    row.energy_limit_trigger_fraction = summary_table.energy_limit_trigger_fraction(ii);
    row.decision = decision;
    decision_rows = local_append_struct(decision_rows, row);
end
decision_table = struct2table(decision_rows);
end

function report = local_build_validation_report( ...
    run_table, disabled_equivalence, tol_invariant, tol_component, tol_energy, tol_w_eta)
report = struct();
report.script = mfilename;
report.created_at = char(datetime('now'));
report.disabled_equivalence = disabled_equivalence;
report.max_channel_invariant_error = max(run_table.channel_invariant_error);
report.max_component_sum_error_rel = max(run_table.component_sum_error_rel);
report.max_energy_conservation_error = max(run_table.energy_conservation_error);
report.max_W_eta_variance_rel_error = max(run_table.W_eta_variance_rel_error);
report.energy_limit_trigger_fraction = mean(double(run_table.energy_limit_applied));
report.tolerances = struct( ...
    'channel_invariant_error', tol_invariant, ...
    'component_sum_error_rel', tol_component, ...
    'energy_conservation_error', tol_energy, ...
    'W_eta_variance_rel_error', tol_w_eta);
report.all_passed = ...
    disabled_equivalence.H_f_max_abs_diff <= tol_invariant && ...
    disabled_equivalence.H_direct_f_max_abs_diff <= tol_invariant && ...
    disabled_equivalence.H_reflect_f_max_abs_diff <= tol_invariant && ...
    report.max_channel_invariant_error <= tol_invariant && ...
    report.max_component_sum_error_rel <= tol_component && ...
    report.max_energy_conservation_error <= tol_energy && ...
    report.max_W_eta_variance_rel_error <= tol_w_eta;
report.notes = ['Fixed-Hs wind sweep: sea_wind_speed changes the normalized PM ', ...
    'spectrum shape only. surface_ssa_scatter_scale is an engineering normalization, ', ...
    'not an experimentally calibrated scattering cross section.'];
end

function local_plot_results(summary_table, decision_table, Hs_list, wind_list, f_list_hz)
colors = lines(numel(wind_list));

for ih = 1:numel(Hs_list)
Hs_target = Hs_list(ih);
hs_mask = summary_table.Hs_target_m == Hs_target;
plot_table = summary_table(hs_mask, :);
suffix = local_hs_suffix(Hs_target, numel(Hs_list) > 1);

fig = figure('Visible', 'off', 'Color', 'w');
hold on
for iw = 1:numel(wind_list)
    subset = sortrows(plot_table(plot_table.wind_speed_mps == wind_list(iw), :), 'f_hz');
    plot(subset.f_hz, subset.E_coh_over_E_inc_mean, '-o', 'Color', colors(iw, :), ...
        'DisplayName', sprintf('U=%g coh', wind_list(iw)));
    plot(subset.f_hz, subset.E_sca_over_E_inc_mean, '--s', 'Color', colors(iw, :), ...
        'DisplayName', sprintf('U=%g sca', wind_list(iw)));
end
xlabel('Frequency (Hz)')
ylabel('Surface spectral energy fraction')
title(sprintf('Specular coherent and incoherent scatter energy, Hs=%g m', Hs_target))
legend('Location', 'bestoutside')
grid on
saveas(fig, ['specular_incoherent_Ecoh_Esca_vs_frequency', suffix, '.png'])
close(fig)

fig = figure('Visible', 'off', 'Color', 'w');
hold on
for iw = 1:numel(wind_list)
    subset = sortrows(plot_table(plot_table.wind_speed_mps == wind_list(iw), :), 'f_hz');
    errorbar(subset.f_hz, subset.scatter_to_coherent_rx_db_mean, ...
        subset.scatter_to_coherent_rx_db_std, '-o', 'Color', colors(iw, :), ...
        'DisplayName', sprintf('U=%g m/s', wind_list(iw)));
end
yline(-20, ':k', 'negligible')
yline(-10, '--k', 'not negligible')
xlabel('Frequency (Hz)')
ylabel('|h_{sca}| / |h_{coh}| (dB)')
title(sprintf('Receiver scatter-to-coherent reflected component, Hs=%g m', Hs_target))
legend('Location', 'best')
grid on
saveas(fig, ['specular_incoherent_rx_ratio_vs_frequency', suffix, '.png'])
close(fig)

fig = figure('Visible', 'off', 'Color', 'w');
hold on
for iw = 1:numel(wind_list)
    subset = sortrows(plot_table(plot_table.wind_speed_mps == wind_list(iw), :), 'f_hz');
    plot(subset.f_hz, subset.abs_h_reflect_coh_mean, '-o', 'Color', colors(iw, :), ...
        'DisplayName', sprintf('U=%g coh', wind_list(iw)));
    plot(subset.f_hz, subset.abs_h_reflect_sca_mean, '--s', 'Color', colors(iw, :), ...
        'DisplayName', sprintf('U=%g sca', wind_list(iw)));
    plot(subset.f_hz, subset.abs_h_reflect_total_mean, ':^', 'Color', colors(iw, :), ...
        'DisplayName', sprintf('U=%g total', wind_list(iw)));
end
xlabel('Frequency (Hz)')
ylabel('Receiver reflected amplitude')
title(sprintf('Receiver reflected component amplitudes, Hs=%g m', Hs_target))
legend('Location', 'bestoutside')
grid on
saveas(fig, ['specular_incoherent_rx_abs_vs_frequency', suffix, '.png'])
close(fig)

rx_matrix = local_matrix_from_table(plot_table, wind_list, f_list_hz, 'scatter_to_coherent_rx_db_mean');
fig = figure('Visible', 'off', 'Color', 'w');
imagesc(f_list_hz, wind_list, rx_matrix)
axis xy
colorbar
xlabel('Frequency (Hz)')
ylabel('Wind speed (m/s)')
title(sprintf('Scatter-to-coherent receiver ratio (dB), Hs=%g m', Hs_target))
saveas(fig, ['specular_incoherent_decision_heatmap', suffix, '.png'])
close(fig)

limit_matrix = local_matrix_from_table(plot_table, wind_list, f_list_hz, 'energy_limit_trigger_fraction');
fig = figure('Visible', 'off', 'Color', 'w');
imagesc(f_list_hz, wind_list, limit_matrix)
axis xy
colorbar
xlabel('Frequency (Hz)')
ylabel('Wind speed (m/s)')
title(sprintf('Energy-limit trigger fraction, Hs=%g m', Hs_target))
saveas(fig, ['specular_incoherent_energy_limit_heatmap', suffix, '.png'])
close(fig)
end

% Touch decision_table so MATLAB's code analyzer knows this input is intentional.
if isempty(decision_table)
    warning('Decision table is empty.')
end
end

function suffix = local_hs_suffix(Hs_target, use_suffix)
if ~use_suffix
    suffix = '';
    return
end
label = strrep(sprintf('%.3g', Hs_target), '.', 'p');
label = strrep(label, '-', 'm');
suffix = ['_Hs', label];
end

function matrix = local_matrix_from_table(table_in, wind_list, f_list_hz, field_name)
matrix = NaN(numel(wind_list), numel(f_list_hz));
for iw = 1:numel(wind_list)
    for ff = 1:numel(f_list_hz)
        mask = table_in.wind_speed_mps == wind_list(iw) & table_in.f_hz == f_list_hz(ff);
        if any(mask)
            matrix(iw, ff) = table_in.(field_name)(find(mask, 1, 'first'));
        end
    end
end
end

function values = local_env_numeric_vector(name, default_values)
raw = getenv(name);
if isempty(raw)
    values = default_values;
    return
end
values = str2num(raw); %#ok<ST2NM>
if isempty(values) || any(~isfinite(values(:)))
    error('%s must be a finite numeric vector.', name);
end
values = values(:).';
end

function values = local_env_hs_list()
raw = getenv('SPEC_INCOH_HS_LIST');
if isempty(raw)
    values = local_env_scalar('SPEC_INCOH_HS_TARGET', 0.5);
else
    values = str2num(raw); %#ok<ST2NM>
    if isempty(values) || any(~isfinite(values(:))) || any(values(:) < 0)
        error('SPEC_INCOH_HS_LIST must be a finite nonnegative numeric vector.');
    end
    values = values(:).';
end
end

function value = local_env_scalar(name, default_value)
values = local_env_numeric_vector(name, default_value);
if numel(values) ~= 1
    error('%s must be a scalar.', name);
end
value = values;
end

function value = local_env_scalar_int(name, default_value)
value = round(local_env_scalar(name, default_value));
if value < 1
    error('%s must be a positive integer.', name);
end
end

function out = local_append_struct(in, row)
if isempty(in)
    out = row;
else
    out = [in; row]; %#ok<AGROW>
end
end

function y = wrapToPiLocal(x)
y = mod(x + pi, 2*pi) - pi;
end

