run(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'bootstrap_project.m'));
%PLOT_LFM_TX_RX_COMPARISON_VERTICAL
% Plot transmitted analytic baseband LFM against representative received
% total/reflected LFM waveforms from the saved K-Stat vs K-Domain result.
%
% This reporting script does not rerun vertical_channel_model. It only reads
% the saved validation MAT file and writes comparison figures.

clear
clc
format compact

result_mat = project_result_file('validation', ...
    'validate_kstat_vs_kdomain_lfm_channel_vertical_result.mat');
out_dir = project_result_dir('validation');

if ~isfile(result_mat)
    error('Result MAT file not found: %s. Run validate_kstat_vs_kdomain_lfm_channel_vertical first.', result_mat);
end

S = load(result_mat, 'lfm', 'case_results', 'cfg');
lfm = S.lfm;
case_results = S.case_results;
cfg = S.cfg;

[kd, ks] = local_get_case_pair(case_results);

figure_total = fullfile(out_dir, 'lfm_tx_vs_rx_total_waveform_compare.png');
local_plot_tx_rx_compare(figure_total, lfm, kd, ks, 'rx_total', ...
    'TX LFM vs total received LFM');

figure_reflect = fullfile(out_dir, 'lfm_tx_vs_rx_reflect_waveform_compare.png');
local_plot_tx_rx_compare(figure_reflect, lfm, kd, ks, 'rx_reflect', ...
    'TX LFM vs reflected received LFM');

fprintf('Saved %s\n', figure_total);
fprintf('Saved %s\n', figure_reflect);

function local_plot_tx_rx_compare(file, lfm, kd, ks, field_name, fig_title)
t_ms = lfm.t_s(:) * 1e3;
tx = lfm.tx_bb(:);
rx_kdomain = kd.(field_name)(:);
rx_kstat = ks.(field_name)(:);

tx_real = local_norm_real(tx);
rx_kdomain_real = local_norm_real(rx_kdomain);
rx_kstat_real = local_norm_real(rx_kstat);

tx_env = local_norm_abs(tx);
rx_kdomain_env = local_norm_abs(rx_kdomain);
rx_kstat_env = local_norm_abs(rx_kstat);

figure('Visible', 'off');
subplot(2, 1, 1);
plot(t_ms, tx_real, 'k-', 'LineWidth', 1.1); hold on
plot(t_ms, rx_kdomain_real, 'b-', 'LineWidth', 1.0);
plot(t_ms, rx_kstat_real, 'r--', 'LineWidth', 1.0);
grid on
xlabel('time (ms)');
ylabel('normalized real part');
legend('TX LFM', 'RX kdomain', 'RX kstat', 'Location', 'best');
title(fig_title);

subplot(2, 1, 2);
plot(t_ms, tx_env, 'k-', 'LineWidth', 1.1); hold on
plot(t_ms, rx_kdomain_env, 'b-', 'LineWidth', 1.0);
plot(t_ms, rx_kstat_env, 'r--', 'LineWidth', 1.0);
grid on
xlabel('time (ms)');
ylabel('normalized envelope');
legend('TX LFM', 'RX kdomain', 'RX kstat', 'Location', 'best');

exportgraphics(gcf, file, 'Resolution', 160);
close(gcf);
end

function y = local_norm_real(x)
x = x(:);
scale = max(abs(real(x)));
if scale <= eps
    y = zeros(size(x));
else
    y = real(x) ./ scale;
end
end

function y = local_norm_abs(x)
x = x(:);
scale = max(abs(x));
if scale <= eps
    y = zeros(size(x));
else
    y = abs(x) ./ scale;
end
end

function [kd, ks] = local_get_case_pair(case_results)
if isempty(case_results)
    error('Representative case was not recorded in the validation result.');
end
kd_idx = find([case_results.branch_code] == 1, 1, 'first');
ks_idx = find([case_results.branch_code] == 2, 1, 'first');
if isempty(kd_idx) || isempty(ks_idx)
    error('Representative case must include both kdomain and kstat branches.');
end
kd = case_results(kd_idx);
ks = case_results(ks_idx);
end
