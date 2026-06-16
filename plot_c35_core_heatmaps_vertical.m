% Plot the five core C3.5 Monte Carlo heatmaps from a saved sweep result.
% This script reads compact statistics only; it does not run PE/WAPE.

clear
format compact

result_file = getenv('C35_RESULT_FILE');
if isempty(result_file)
    result_file = 'sweep_monte_carlo_surface_channel_vertical_result.mat';
end
figure_prefix = getenv('C35_CORE_FIGURE_PREFIX');
if isempty(figure_prefix)
    figure_prefix = 'c35_core_';
end

if ~exist(result_file, 'file')
    error('plot_c35_core_heatmaps_vertical:MissingFile', ...
        'C3.5 result file not found: %s', result_file);
end

S = load(result_file, 'sweep_stats', 'condition_summary_table', ...
    'sea_hs_values', 'sea_wind_values', 'seed_list');
if ~isfield(S, 'sweep_stats')
    error('plot_c35_core_heatmaps_vertical:MissingSweepStats', ...
        'Result file does not contain sweep_stats: %s', result_file);
end

plot_specs = { ...
    'abs_H_ref_mean', '|H(f_ref)| mean', 'abs_H_ref_mean_heatmap'; ...
    'abs_H_ref_std', '|H(f_ref)| std', 'abs_H_ref_std_heatmap'; ...
    'reflect_rms_delta_k_mean', 'reflect rms delta k mean (rad/m)', 'reflect_rms_delta_k_mean_heatmap'; ...
    'reflect_high_k_fraction_mean', 'reflect high-k fraction mean', 'reflect_high_k_fraction_mean_heatmap'; ...
    'tap_rms_delay_symbols_mean', 'tap RMS delay mean (symbols)', 'tap_rms_delay_mean_heatmap'};

for ii = 1:size(plot_specs, 1)
    local_plot_heatmap(S.sweep_stats, plot_specs{ii, 1}, ...
        plot_specs{ii, 2}, [figure_prefix plot_specs{ii, 3} '.png']);
end

summary = struct();
summary.result_file = result_file;
summary.figure_prefix = figure_prefix;
summary.generated_files = strcat(figure_prefix, plot_specs(:, 3), '.png');
if isfield(S, 'condition_summary_table')
    summary.condition_count = height(S.condition_summary_table);
    summary.max_invariant_error = max(S.condition_summary_table.max_invariant_error);
    summary.max_direct_drift_abs = max(S.condition_summary_table.direct_drift_max_abs);
    summary.mc_count_unique = unique(S.condition_summary_table.mc_count).';
end
if isfield(S, 'seed_list')
    summary.seed_count = numel(S.seed_list);
end

disp(summary)

function local_plot_heatmap(sweep_stats, field_name, title_text, file_name)
if ~isfield(sweep_stats, field_name)
    error('plot_c35_core_heatmaps_vertical:MissingField', ...
        'sweep_stats.%s is missing.', field_name);
end

figure('Visible', 'off');
values = sweep_stats.(field_name);
imagesc(sweep_stats.sea_wind_values, sweep_stats.sea_hs_values, values)
set(gca, 'YDir', 'normal')
colorbar
grid on
xlabel('Wind speed (m/s)')
ylabel('H_s target (m)')
title(title_text, 'Interpreter', 'none')
for ih = 1:numel(sweep_stats.sea_hs_values)
    for iw = 1:numel(sweep_stats.sea_wind_values)
        if isfinite(values(ih, iw))
            text(sweep_stats.sea_wind_values(iw), sweep_stats.sea_hs_values(ih), ...
                sprintf('%.3g', values(ih, iw)), ...
                'HorizontalAlignment', 'center', 'Color', 'w', 'FontWeight', 'bold');
        end
    end
end
print(gcf, '-dpng', '-r220', file_name)
close(gcf)
end
