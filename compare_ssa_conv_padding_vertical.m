% Compare periodic and zero-padded convolution paths for ssa_stat_kernel.
% This reduced scalar-frequency script stores compact tables and light
% figures only; it does not save full fields or 2-D spectra.

clear
format compact

result_file = getenv('SSA_CONV_COMPARE_RESULT_FILE');
if isempty(result_file)
    result_file = 'compare_ssa_conv_padding_vertical_result.mat';
end

figure_prefix = getenv('SSA_CONV_COMPARE_FIGURE_PREFIX');
if isempty(figure_prefix)
    figure_prefix = 'ssa_conv_padding_compare_';
end

sea_hs_values = [0, 0.05, 0.2];
kernel_specs = local_kernel_specs();
padding_modes = {'periodic', 'zero_padded'};

mc_count = 2;
mc_override = str2double(getenv('SSA_CONV_COMPARE_MC_COUNT'));
if isfinite(mc_override) && mc_override >= 1
    mc_count = round(mc_override);
end
seed_list = 12345 + (0:(mc_count - 1));

max_pairs = Inf;
max_pairs_override = str2double(getenv('SSA_CONV_COMPARE_MAX_PAIRS'));
if isfinite(max_pairs_override) && max_pairs_override >= 1
    max_pairs = round(max_pairs_override);
end

base_params = local_base_params();
run_rows = struct([]);
pair_rows = struct([]);
pair_count = 0;

for ik = 1:numel(kernel_specs)
    for ih = 1:numel(sea_hs_values)
        for iseed = 1:numel(seed_list)
            if pair_count >= max_pairs
                break
            end
            pair_count = pair_count + 1;
            hs_value = sea_hs_values(ih);
            sea_seed = seed_list(iseed);
            fprintf('SSA conv compare pair %d: kernel=%s, Hs=%g, seed=%d\n', ...
                pair_count, kernel_specs(ik).kernel_mode, hs_value, sea_seed);

            channel_by_padding = cell(1, numel(padding_modes));
            row_by_padding = cell(1, numel(padding_modes));
            for ip = 1:numel(padding_modes)
                paramsV = base_params;
                paramsV.sea_hs_target = hs_value;
                paramsV.sea_seed = sea_seed;
                paramsV.surface_ssa_kernel_mode = kernel_specs(ik).kernel_mode;
                paramsV.surface_ssa_geometry_source_id = kernel_specs(ik).geometry_source_id;
                paramsV.surface_ssa_conv_padding = padding_modes{ip};

                channel = vertical_channel_model(paramsV);
                flat_diff = NaN;
                if hs_value == 0
                    flat_channel = local_flat_kirchhoff_reference(base_params, sea_seed);
                    flat_diff = max(abs(channel.H_f(:) - flat_channel.H_f(:)));
                end

                row = local_build_run_row(pair_count, kernel_specs(ik), hs_value, ...
                    sea_seed, padding_modes{ip}, channel, flat_diff);
                run_rows = local_append_struct(run_rows, row);
                channel_by_padding{ip} = channel;
                row_by_padding{ip} = row;
            end

            pair_row = local_build_pair_row(pair_count, kernel_specs(ik), hs_value, ...
                sea_seed, row_by_padding{1}, row_by_padding{2}, ...
                channel_by_padding{1}, channel_by_padding{2});
            pair_rows = local_append_struct(pair_rows, pair_row);
        end
        if pair_count >= max_pairs
            break
        end
    end
    if pair_count >= max_pairs
        break
    end
end

run_summary_table = struct2table(run_rows);
pair_summary_table = struct2table(pair_rows);
condition_summary_table = local_build_condition_summary(pair_summary_table);

disp(run_summary_table(:, {'kernel_mode', 'sea_hs_target', 'sea_seed', ...
    'conv_padding', 'E_sca_raw_over_E_inc', 'E_sca_limited_over_E_inc', ...
    'E_ref_over_E_inc', 'abs_h_reflect', 'abs_h_total', ...
    'energy_conservation_error', 'flat_Hs0_diff'}))
disp(pair_summary_table(:, {'kernel_mode', 'sea_hs_target', 'sea_seed', ...
    'rel_diff_H_f', 'rel_diff_H_reflect_f', ...
    'rel_diff_E_sca_raw_over_E_inc', 'rel_diff_scatter_rms_delta_k'}))
disp(condition_summary_table)

figure_files = local_write_figures(figure_prefix, run_summary_table, pair_summary_table);

save(result_file, 'base_params', 'sea_hs_values', 'kernel_specs', ...
    'padding_modes', 'seed_list', 'run_summary_table', 'pair_summary_table', ...
    'condition_summary_table', 'figure_files');
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
paramsV.surface_ssa_random_scatter = true;
paramsV.surface_ssa_scatter_scale = 1.0;
paramsV.surface_ssa_seed_offset = 100000;
paramsV.surface_ssa_kernel_mode = 'pm_convolution';
paramsV.surface_ssa_geometry_source_id = '';
paramsV.surface_ssa_kz_branch = 'downward_positive_real';
paramsV.surface_ssa_conv_padding = 'periodic';
end

function specs = local_kernel_specs()
source_id = ['SSA.md; Thorsos & Broschat 1995 JASA, ', ...
    'Dirichlet SSA first-order / perturbation-limit geometry'];
specs = struct([]);
specs(1).kernel_mode = 'pm_convolution';
specs(1).geometry_source_id = '';
specs(2).kernel_mode = 'ssa1_geometry';
specs(2).geometry_source_id = source_id;
end

function flat_channel = local_flat_kirchhoff_reference(base_params, sea_seed)
paramsV = base_params;
paramsV.surface_boundary_model = 'kirchhoff_spatial';
paramsV.surface_ssa_kernel_mode = 'pm_convolution';
paramsV.surface_ssa_conv_padding = 'periodic';
paramsV.sea_hs_target = 0;
paramsV.sea_seed = sea_seed;
flat_channel = vertical_channel_model(paramsV);
end

function row = local_build_run_row(pair_index, kernel_spec, hs_value, sea_seed, ...
    conv_padding, channel, flat_diff)
ssa = channel.roughness_meta.ssa_stat_kernel_meta;
invariant_error = max(abs(channel.H_f(:) - channel.H_direct_f(:) - channel.H_reflect_f(:)));
row = struct();
row.pair_index = pair_index;
row.kernel_mode = string(kernel_spec.kernel_mode);
row.sea_hs_target = hs_value;
row.sea_seed = sea_seed;
row.conv_padding = string(conv_padding);
row.conv_operator = string(ssa.conv_operator);
row.invariant_error = invariant_error;
row.flat_Hs0_diff = flat_diff;
row.abs_h_total = abs(channel.h_total);
row.abs_h_reflect = abs(channel.h_reflect);
row.abs_h_direct = abs(channel.h_direct);
row.phase_h_total_rad = angle(channel.h_total);
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
row.energy_limit_applied = ssa.energy_limit_applied;
row.energy_conservation_error = ssa.energy_conservation_error;
row.scatter_rms_delta_k_rad_per_m = ssa.scatter_power_spectrum_stats.rms_delta_k_rad_per_m;
row.reflected_rms_delta_k_rad_per_m = ssa.reflected_spectrum_stats.rms_delta_k_rad_per_m;
row.scatter_energy_radius_90_rad_per_m = ssa.scatter_power_spectrum_stats.energy_radius_90_rad_per_m;
end

function row = local_build_pair_row(pair_index, kernel_spec, hs_value, sea_seed, ...
    periodic_row, zero_row, periodic_channel, zero_channel)
row = struct();
row.pair_index = pair_index;
row.kernel_mode = string(kernel_spec.kernel_mode);
row.sea_hs_target = hs_value;
row.sea_seed = sea_seed;
row.periodic_conv_operator = string(periodic_row.conv_operator);
row.zero_padded_conv_operator = string(zero_row.conv_operator);
row.periodic_invariant_error = periodic_row.invariant_error;
row.zero_padded_invariant_error = zero_row.invariant_error;
row.periodic_energy_conservation_error = periodic_row.energy_conservation_error;
row.zero_padded_energy_conservation_error = zero_row.energy_conservation_error;
row.flat_Hs0_diff_periodic = periodic_row.flat_Hs0_diff;
row.flat_Hs0_diff_zero_padded = zero_row.flat_Hs0_diff;
row.periodic_E_sca_raw_over_E_inc = periodic_row.E_sca_raw_over_E_inc;
row.zero_padded_E_sca_raw_over_E_inc = zero_row.E_sca_raw_over_E_inc;
row.periodic_E_sca_limited_over_E_inc = periodic_row.E_sca_limited_over_E_inc;
row.zero_padded_E_sca_limited_over_E_inc = zero_row.E_sca_limited_over_E_inc;
row.periodic_E_ref_over_E_inc = periodic_row.E_ref_over_E_inc;
row.zero_padded_E_ref_over_E_inc = zero_row.E_ref_over_E_inc;
row.periodic_abs_h_reflect = periodic_row.abs_h_reflect;
row.zero_padded_abs_h_reflect = zero_row.abs_h_reflect;
row.periodic_abs_h_total = periodic_row.abs_h_total;
row.zero_padded_abs_h_total = zero_row.abs_h_total;
row.periodic_scatter_rms_delta_k = periodic_row.scatter_rms_delta_k_rad_per_m;
row.zero_padded_scatter_rms_delta_k = zero_row.scatter_rms_delta_k_rad_per_m;
row.periodic_reflected_rms_delta_k = periodic_row.reflected_rms_delta_k_rad_per_m;
row.zero_padded_reflected_rms_delta_k = zero_row.reflected_rms_delta_k_rad_per_m;
row.rel_diff_H_f = local_rel_l2(zero_channel.H_f(:), periodic_channel.H_f(:));
row.rel_diff_H_reflect_f = local_rel_l2(zero_channel.H_reflect_f(:), periodic_channel.H_reflect_f(:));
row.rel_diff_abs_h_reflect = local_rel_abs(zero_row.abs_h_reflect, periodic_row.abs_h_reflect);
row.rel_diff_abs_h_total = local_rel_abs(zero_row.abs_h_total, periodic_row.abs_h_total);
row.rel_diff_E_sca_raw_over_E_inc = local_rel_abs( ...
    zero_row.E_sca_raw_over_E_inc, periodic_row.E_sca_raw_over_E_inc);
row.rel_diff_E_sca_limited_over_E_inc = local_rel_abs( ...
    zero_row.E_sca_limited_over_E_inc, periodic_row.E_sca_limited_over_E_inc);
row.rel_diff_E_ref_over_E_inc = local_rel_abs( ...
    zero_row.E_ref_over_E_inc, periodic_row.E_ref_over_E_inc);
row.rel_diff_scatter_rms_delta_k = local_rel_abs( ...
    zero_row.scatter_rms_delta_k_rad_per_m, periodic_row.scatter_rms_delta_k_rad_per_m);
row.rel_diff_reflected_rms_delta_k = local_rel_abs( ...
    zero_row.reflected_rms_delta_k_rad_per_m, periodic_row.reflected_rms_delta_k_rad_per_m);
end

function summary_table = local_build_condition_summary(pair_table)
if isempty(pair_table)
    summary_table = table();
    return
end
[groups, kernel_modes, hs_values] = findgroups(pair_table.kernel_mode, pair_table.sea_hs_target);
summary_table = table();
summary_table.kernel_mode = kernel_modes;
summary_table.sea_hs_target = hs_values;
summary_table.seed_count = splitapply(@numel, pair_table.sea_seed, groups);
summary_table.rel_diff_H_f_mean = splitapply(@local_nanmean, pair_table.rel_diff_H_f, groups);
summary_table.rel_diff_H_reflect_f_mean = splitapply(@local_nanmean, pair_table.rel_diff_H_reflect_f, groups);
summary_table.rel_diff_E_sca_raw_mean = splitapply(@local_nanmean, pair_table.rel_diff_E_sca_raw_over_E_inc, groups);
summary_table.rel_diff_scatter_rms_delta_k_mean = splitapply(@local_nanmean, pair_table.rel_diff_scatter_rms_delta_k, groups);
summary_table.periodic_E_sca_raw_over_E_inc_mean = splitapply(@local_nanmean, ...
    pair_table.periodic_E_sca_raw_over_E_inc, groups);
summary_table.zero_padded_E_sca_raw_over_E_inc_mean = splitapply(@local_nanmean, ...
    pair_table.zero_padded_E_sca_raw_over_E_inc, groups);
summary_table.periodic_E_ref_over_E_inc_mean = splitapply(@local_nanmean, ...
    pair_table.periodic_E_ref_over_E_inc, groups);
summary_table.zero_padded_E_ref_over_E_inc_mean = splitapply(@local_nanmean, ...
    pair_table.zero_padded_E_ref_over_E_inc, groups);
pair_max_invariant = max([pair_table.periodic_invariant_error, ...
    pair_table.zero_padded_invariant_error], [], 2);
summary_table.max_invariant_error = splitapply(@max, pair_max_invariant, groups);
pair_max_energy_error = max([pair_table.periodic_energy_conservation_error, ...
    pair_table.zero_padded_energy_conservation_error], [], 2);
summary_table.max_energy_conservation_error = splitapply(@max, pair_max_energy_error, groups);
end

function figure_files = local_write_figures(prefix, run_table, pair_table)
figure_files = strings(0, 1);
if isempty(run_table) || isempty(pair_table)
    return
end

label = strcat(cellstr(run_table.kernel_mode), "_Hs", compose('%.2g', run_table.sea_hs_target), ...
    "_", cellstr(run_table.conv_padding));
fig = figure('Visible', 'off', 'Color', 'w');
bar(categorical(label), run_table.E_sca_raw_over_E_inc);
ylabel('E_{sca}^{raw}/E_{inc}');
title('Raw scatter energy by convolution path');
xtickangle(45);
grid on
file1 = [prefix 'raw_scatter_energy.png'];
exportgraphics(fig, file1, 'Resolution', 180);
close(fig);
figure_files(end + 1, 1) = string(file1);

pair_label = strcat(cellstr(pair_table.kernel_mode), "_Hs", compose('%.2g', pair_table.sea_hs_target), ...
    "_seed", compose('%d', pair_table.sea_seed));
fig = figure('Visible', 'off', 'Color', 'w');
bar(categorical(pair_label), [pair_table.rel_diff_H_reflect_f, ...
    pair_table.rel_diff_E_sca_raw_over_E_inc, pair_table.rel_diff_scatter_rms_delta_k]);
ylabel('relative difference');
title('zero-padded vs periodic differences');
legend({'H_{reflect}', 'E_{sca}^{raw}', 'scatter RMS K'}, 'Location', 'best');
xtickangle(45);
grid on
file2 = [prefix 'relative_difference.png'];
exportgraphics(fig, file2, 'Resolution', 180);
close(fig);
figure_files(end + 1, 1) = string(file2);
end

function y = local_rel_l2(a, b)
y = norm(a(:) - b(:)) / max(norm(b(:)), eps);
end

function y = local_rel_abs(a, b)
if ~(isfinite(a) && isfinite(b))
    y = NaN;
else
    y = abs(a - b) / max(abs(b), eps);
end
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

function out = local_append_struct(out, row)
if isempty(out)
    out = row;
else
    out(end + 1) = row; %#ok<AGROW>
end
end
