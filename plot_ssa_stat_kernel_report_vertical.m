% Generate report-ready figures for the ssa_stat_kernel validation sweep.
% This script only reads saved .mat result files and does not run PE/WAPE.

clear
format compact

channel_result_file = getenv('SSA_STAT_CHANNEL_RESULT_FILE');
if isempty(channel_result_file)
    channel_result_file = 'sweep_ssa_stat_kernel_surface_channel_vertical_result.mat';
end
comm_result_file = getenv('SSA_STAT_COMM_RESULT_FILE');
if isempty(comm_result_file)
    comm_result_file = 'monte_carlo_comm_ssa_stat_kernel_psk_vertical_result.mat';
end
validate_result_file = getenv('SSA_STAT_VALIDATE_RESULT_FILE');
if isempty(validate_result_file)
    validate_result_file = 'validate_ssa_stat_kernel_vertical_result.mat';
end
figure_prefix = getenv('SSA_STAT_REPORT_PREFIX');
if isempty(figure_prefix)
    figure_prefix = 'ssa_stat_kernel_report_';
end

if ~exist(channel_result_file, 'file')
    error('Channel result file not found: %s', channel_result_file);
end
if ~exist(comm_result_file, 'file')
    error('Communication result file not found: %s', comm_result_file);
end

channel_data = load(channel_result_file);
comm_data = load(comm_result_file);
if exist(validate_result_file, 'file')
    validate_data = load(validate_result_file);
else
    validate_data = struct();
end

local_plot_frequency_stats(channel_data.frequency_stats, channel_data.model_names, ...
    channel_data.sea_hs_values, 'abs_H_f_mean', 'abs_H_f_std', ...
    '|H(f)| mean/std', [figure_prefix 'abs_H_f_mean_std.png']);
local_plot_frequency_stats(channel_data.frequency_stats, channel_data.model_names, ...
    channel_data.sea_hs_values, 'abs_H_reflect_f_mean', 'abs_H_reflect_f_std', ...
    '|H_reflect(f)| mean/std', [figure_prefix 'abs_H_reflect_f_mean_std.png']);
local_plot_phase_vs_hs(channel_data.condition_summary_table, channel_data.model_names, ...
    channel_data.sea_hs_values, [figure_prefix 'phase_H_f_ref_vs_Hs.png']);
local_plot_scalar_vs_hs(channel_data.condition_summary_table, channel_data.model_names, ...
    'reflect_frequency_energy_ratio_mean', 'reflect_frequency_energy_ratio_std', ...
    'sum |H_reflect(f)|^2 / sum |H_direct(f)|^2', ...
    [figure_prefix 'reflect_energy_vs_Hs.png']);
ssa_kernel_plot_names = local_ssa_kernel_plot_names(channel_data.model_names);
local_plot_scalar_vs_hs(channel_data.condition_summary_table, ssa_kernel_plot_names, ...
    'E_sca_over_E_inc_mean', 'E_sca_over_E_inc_std', ...
    'E_sca / E_inc', [figure_prefix 'Esca_Einc_vs_Hs.png']);
local_plot_scalar_vs_hs(channel_data.condition_summary_table, ssa_kernel_plot_names, ...
    'energy_scale_applied_mean', 'energy_scale_applied_std', ...
    'energy scale applied', [figure_prefix 'energy_scale_vs_Hs.png']);
local_plot_scalar_vs_hs(channel_data.condition_summary_table, channel_data.model_names, ...
    'abs_h_total_mean', 'abs_h_total_std', '|h_total| at f_ref', ...
    [figure_prefix 'abs_h_total_model_compare.png']);
local_plot_scalar_vs_hs(channel_data.condition_summary_table, channel_data.model_names, ...
    'abs_h_reflect_mean', 'abs_h_reflect_std', '|h_reflect| at f_ref', ...
    [figure_prefix 'abs_h_reflect_model_compare.png']);
local_plot_comm_metric(comm_data.curve_stats, comm_data.model_names, ...
    comm_data.sea_hs_values, 'BER', [figure_prefix 'BER_model_compare.png']);
local_plot_comm_metric(comm_data.curve_stats, comm_data.model_names, ...
    comm_data.sea_hs_values, 'SER', [figure_prefix 'SER_model_compare.png']);
local_plot_hs0_flat_check(channel_data.flat_check_table, ...
    [figure_prefix 'Hs0_flat_check.png']);

% Kernel-mode report aliases for the SSA upgrade interface.
local_plot_frequency_stats(channel_data.frequency_stats, channel_data.model_names, ...
    channel_data.sea_hs_values, 'abs_H_f_mean', 'abs_H_f_std', ...
    '|H(f)| mean/std by kernel mode', 'ssa_kernel_mode_abs_H_f_mean_std.png');
local_plot_frequency_stats(channel_data.frequency_stats, channel_data.model_names, ...
    channel_data.sea_hs_values, 'abs_H_reflect_f_mean', 'abs_H_reflect_f_std', ...
    '|H_reflect(f)| mean/std by kernel mode', 'ssa_kernel_mode_abs_H_reflect_f_mean_std.png');
local_plot_scalar_vs_hs(channel_data.condition_summary_table, ssa_kernel_plot_names, ...
    'R_coh_abs_mean', 'R_coh_abs_std', ...
    '|R_coh| for pm_convolution baseline', 'ssa_kernel_mode_Rcoh_vs_Hs.png');
local_plot_scalar_vs_hs(channel_data.condition_summary_table, ssa_kernel_plot_names, ...
    'E_sca_over_E_inc_mean', 'E_sca_over_E_inc_std', ...
    'E_sca / E_inc by kernel mode', 'ssa_kernel_mode_Esca_Einc_vs_Hs.png');
local_plot_scalar_vs_hs(channel_data.condition_summary_table, ssa_kernel_plot_names, ...
    'energy_conservation_error_max', 'energy_scale_applied_std', ...
    'max energy conservation error', 'ssa_kernel_mode_energy_conservation_error.png');
local_plot_scalar_vs_hs(channel_data.condition_summary_table, ssa_kernel_plot_names, ...
    'P_sca_sum_mean', 'P_sca_sum_std', ...
    'P_sca sum frequency trend proxy', 'ssa_kernel_mode_frequency_trend.png');
local_plot_comm_metric_pair(comm_data.curve_stats, comm_data.model_names, ...
    'ssa_kernel_mode_BER_SER_compare.png');
local_plot_dense_vs_fft(validate_data, 'ssa_kernel_mode_dense_vs_fft_error.png');

channel_condition_summary_table = channel_data.condition_summary_table;
channel_metadata_stats = channel_data.metadata_stats;
flat_check_table = channel_data.flat_check_table;
comm_curve_summary_table = comm_data.curve_summary_table;
save([figure_prefix 'metadata_energy_table.mat'], ...
    'channel_condition_summary_table', 'channel_metadata_stats', ...
    'flat_check_table', 'comm_curve_summary_table');

fprintf('Saved report figures with prefix %s\n', figure_prefix);

function local_plot_frequency_stats(freq_stats, model_names, sea_hs_values, mean_field, std_field, title_text, file_name)
figure('Visible', 'off');
n_hs = numel(sea_hs_values);
n_row = ceil(n_hs / 2);
for ih = 1:n_hs
    subplot(n_row, 2, ih)
    hold on
    for im = 1:numel(model_names)
        idx = local_find_freq(freq_stats, model_names{im}, sea_hs_values(ih));
        if isempty(idx)
            continue
        end
        f_khz = freq_stats(idx).f_axis(:) / 1000;
        y = freq_stats(idx).(mean_field)(:);
        e = freq_stats(idx).(std_field)(:);
        errorbar(f_khz, y, e, 'o-', 'LineWidth', 1.0, ...
            'DisplayName', model_names{im});
    end
    grid on
    xlabel('frequency (kHz)')
    ylabel(strrep(mean_field, '_', '\_'))
    title(sprintf('Hs=%g m', sea_hs_values(ih)), 'Interpreter', 'none')
    if ih == 1
        legend('Location', 'best')
    end
end

sgtitle(title_text, 'Interpreter', 'none')
print(gcf, '-dpng', '-r220', file_name)
close(gcf)
end

function names = local_ssa_kernel_plot_names(model_names)
names = {};
for ii = 1:numel(model_names)
    name = char(model_names{ii});
    if startsWith(name, 'ssa')
        names{end + 1} = name; %#ok<AGROW>
    end
end
if isempty(names)
    names = {'ssa_stat_kernel'};
end
end

function local_plot_text_notice(message_text, file_name)
figure('Visible', 'off');
axis off
text(0.05, 0.55, message_text, 'Interpreter', 'none', ...
    'FontSize', 11, 'HorizontalAlignment', 'left')
print(gcf, '-dpng', '-r220', file_name)
close(gcf)
end

function local_plot_dense_vs_fft(validate_data, file_name)
if ~isfield(validate_data, 'case_results') || ...
        ~isfield(validate_data.case_results, 'ssa1_dense_compare')
    local_plot_text_notice('Dense vs FFT SSA1 result not found. Rerun validate_ssa_stat_kernel_vertical.', file_name);
    return
end
cmp = validate_data.case_results.ssa1_dense_compare;
values = [cmp.P_sca_raw_sum_diff, cmp.P_sca_sum_diff];
figure('Visible', 'off');
bar(values)
set(gca, 'XTickLabel', {'P_sca_raw sum diff', 'P_sca sum diff'})
ylabel('absolute difference')
title('ssa1 debug dense vs FFT periodic sum', 'Interpreter', 'none')
grid on
print(gcf, '-dpng', '-r220', file_name)
close(gcf)
end

function local_plot_phase_vs_hs(T, model_names, sea_hs_values, file_name)
figure('Visible', 'off');
hold on
for im = 1:numel(model_names)
    [x, y, e] = local_table_series(T, model_names{im}, sea_hs_values, ...
        'phase_h_total_circular_mean_rad', 'phase_h_total_circular_std_rad');
    errorbar(x, y, e, 'o-', 'LineWidth', 1.1, 'DisplayName', model_names{im});
end
grid on
xlabel('Hs target (m)')
ylabel('circular phase mean at f_ref (rad)')
title('phase(H(f_ref)) vs Hs', 'Interpreter', 'none')
legend('Location', 'best')
print(gcf, '-dpng', '-r220', file_name)
close(gcf)
end

function local_plot_scalar_vs_hs(T, model_names, mean_field, std_field, ylabel_text, file_name)
if ~ismember(mean_field, T.Properties.VariableNames)
    local_plot_text_notice(sprintf('Field %s is not available. Rerun the sweep script to regenerate current metadata.', mean_field), ...
        file_name);
    return
end
figure('Visible', 'off');
hold on
sea_hs_values = unique(T.sea_hs_target).';
for im = 1:numel(model_names)
    [x, y, e] = local_table_series(T, model_names{im}, sea_hs_values, mean_field, std_field);
    if any(isfinite(y))
        errorbar(x, y, e, 'o-', 'LineWidth', 1.1, 'DisplayName', model_names{im});
    end
end
grid on
xlabel('Hs target (m)')
ylabel(ylabel_text, 'Interpreter', 'none')
title(ylabel_text, 'Interpreter', 'none')
legend('Location', 'best')
print(gcf, '-dpng', '-r220', file_name)
close(gcf)
end

function local_plot_comm_metric(curve_stats, model_names, sea_hs_values, metric_name, file_name)
figure('Visible', 'off');
n_hs = numel(sea_hs_values);
n_row = ceil(n_hs / 2);
for ih = 1:n_hs
    subplot(n_row, 2, ih)
    hold on
    for im = 1:numel(model_names)
        local_plot_metric_curve(curve_stats, metric_name, im, ih, model_names{im});
    end
    grid on
    set(gca, 'YScale', 'log')
    xlabel('Eb/N0 (dB)')
    ylabel(metric_name)
    title(sprintf('Hs=%g m', sea_hs_values(ih)), 'Interpreter', 'none')
    if ih == 1
        legend('Location', 'southwest')
    end
end
sgtitle([metric_name ' mean with +/-1 std'], 'Interpreter', 'none')
print(gcf, '-dpng', '-r220', file_name)
close(gcf)
end

function local_plot_comm_metric_pair(curve_stats, model_names, file_name)
figure('Visible', 'off');
metric_names = {'BER', 'SER'};
for imet = 1:numel(metric_names)
    subplot(1, 2, imet)
    hold on
    metric_name = metric_names{imet};
    stats = curve_stats.(metric_name);
    for ih = 1:numel(curve_stats.sea_hs_values)
        for im = 1:numel(model_names)
            y = stats.mean(:, im, ih);
            plot(curve_stats.EbN0_dB_list(:), max(y, 1e-5), 'o-', ...
                'LineWidth', 1.0, ...
                'DisplayName', sprintf('%s Hs=%g', model_names{im}, curve_stats.sea_hs_values(ih)));
        end
    end
    grid on
    set(gca, 'YScale', 'log')
    xlabel('Eb/N0 (dB)')
    ylabel(metric_name)
    title(metric_name, 'Interpreter', 'none')
    if imet == 1
        legend('Location', 'southwest')
    end
end
sgtitle('BER/SER by kernel mode', 'Interpreter', 'none')
print(gcf, '-dpng', '-r220', file_name)
close(gcf)
end

function local_plot_metric_curve(curve_stats, metric_name, model_idx, hs_idx, label_text)
stats = curve_stats.(metric_name);
x = curve_stats.EbN0_dB_list(:);
y = stats.mean(:, model_idx, hs_idx);
e = stats.std(:, model_idx, hs_idx);
plot_floor = 1e-5;
plot(x, max(y, plot_floor), 'o-', 'LineWidth', 1.1, 'DisplayName', label_text)
plot(x, max(y - e, plot_floor), '--', 'LineWidth', 0.8, 'HandleVisibility', 'off')
plot(x, max(y + e, plot_floor), '--', 'LineWidth', 0.8, 'HandleVisibility', 'off')
end

function local_plot_hs0_flat_check(flat_check_table, file_name)
figure('Visible', 'off');
if isempty(flat_check_table)
    text(0.1, 0.5, 'Hs=0 flat check not available', 'Interpreter', 'none')
    axis off
else
    values = [flat_check_table.max_abs_H_f_diff, flat_check_table.max_abs_H_reflect_f_diff];
    bar(values)
    hold on
    yline(flat_check_table.tolerance(1), 'r--', 'Tolerance')
    set(gca, 'XTickLabel', {'H_f diff', 'H_reflect_f diff'})
    ylabel('max absolute difference')
    title('Hs=0 flat-surface degeneration check', 'Interpreter', 'none')
    grid on
end
print(gcf, '-dpng', '-r220', file_name)
close(gcf)
end

function idx = local_find_freq(freq_stats, model_name, hs_value)
idx = [];
for ii = 1:numel(freq_stats)
    if strcmp(char(freq_stats(ii).model_name), model_name) && ...
            abs(freq_stats(ii).sea_hs_target - hs_value) <= 1e-12
        idx = ii;
        return
    end
end
end

function [x, y, e] = local_table_series(T, model_name, sea_hs_values, mean_field, std_field)
x = sea_hs_values(:);
y = NaN(size(x));
e = NaN(size(x));
for ii = 1:numel(x)
    mask = strcmp(string(T.model_name), model_name) & abs(T.sea_hs_target - x(ii)) <= 1e-12;
    idx = find(mask, 1, 'first');
    if isempty(idx)
        continue
    end
    y(ii) = T.(mean_field)(idx);
    if ismember(std_field, T.Properties.VariableNames)
        e(ii) = T.(std_field)(idx);
    else
        e(ii) = 0;
    end
end
end
