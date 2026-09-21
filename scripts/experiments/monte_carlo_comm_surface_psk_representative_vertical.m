run(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'bootstrap_project.m'));
% Representative-sea-state QPSK BER/SER Monte Carlo for the D2 link.
% Uses weak/mid/strong sea states and keeps the D2 receive-window and Eb/N0
% policies fixed. This script calls vertical_channel_model; it does not change the
% propagation model or communication reference implementation.

clear
format compact

result_file = getenv('SURFACE_COMM_REP_RESULT_FILE');
if isempty(result_file)
    result_file = 'monte_carlo_comm_surface_psk_representative_result.mat';
end
figure_prefix = getenv('SURFACE_COMM_REP_FIGURE_PREFIX');
if isempty(figure_prefix)
    figure_prefix = 'monte_carlo_comm_surface_psk_representative_';
end
if exist(result_file, 'file')
    warning('Result file %s already exists and will be overwritten.', result_file);
end

params_base = local_base_params();
sea_conditions = local_representative_sea_conditions();
scenarios = local_scenarios();
comm_cfg = local_comm_cfg();

mc_count = 16;
mc_override = str2double(getenv('SURFACE_COMM_REP_MC_COUNT'));
if isfinite(mc_override) && mc_override >= 1
    mc_count = round(mc_override);
end
max_conditions = Inf;
max_conditions_override = str2double(getenv('SURFACE_COMM_REP_MAX_CONDITIONS'));
if isfinite(max_conditions_override) && max_conditions_override >= 1
    max_conditions = round(max_conditions_override);
end
if isfinite(max_conditions)
    sea_conditions = sea_conditions(1:min(max_conditions, numel(sea_conditions)));
end
seed_list = 12345 + (0:(mc_count - 1));

rng(comm_cfg.bits_seed, 'twister');
bits_per_symbol = round(log2(comm_cfg.M));
bits_tx = randi([0, 1], comm_cfg.n_sym * bits_per_symbol, 1);
[tx_symbols, bits_used] = modem_psk('modulate', bits_tx, comm_cfg.M);
bits_tx = bits_used;

n_cond = numel(sea_conditions);
n_scen = numel(scenarios);
n_eb = numel(comm_cfg.EbN0_dB_list);
ber_samples = NaN(n_eb, mc_count, n_scen, n_cond);
ser_samples = NaN(n_eb, mc_count, n_scen, n_cond);
effective_snr_samples = NaN(n_eb, mc_count, n_scen, n_cond);
run_rows = struct([]);

for cc = 1:n_cond
    for ss = 1:n_scen
        for mm = 1:mc_count
            paramsV = params_base;
            paramsV.sea_hs_target = sea_conditions(cc).sea_hs_target;
            paramsV.sea_wind_speed = sea_conditions(cc).sea_wind_speed;
            paramsV.enable_surface_reflection = scenarios(ss).enable_surface_reflection;
            paramsV.sea_seed = seed_list(mm);

            fprintf(['Representative comm MC condition %d/%d (%s), scenario %d/%d (%s), ' ...
                'seed %d/%d, sea_seed=%d\n'], ...
                cc, n_cond, sea_conditions(cc).name, ss, n_scen, scenarios(ss).name, ...
                mm, mc_count, seed_list(mm));
            channel = vertical_channel_model(paramsV);

            invariant_error = max(abs(channel.H_f(:) - channel.H_direct_f(:) - channel.H_reflect_f(:)));
            if invariant_error > 1e-10
                error('monte_carlo_comm_surface_psk_representative_vertical:Invariant', ...
                    'H_f invariant failed for %s/%s seed %d: %.3e', ...
                    sea_conditions(cc).name, scenarios(ss).name, seed_list(mm), invariant_error);
            end
            if ~paramsV.enable_surface_reflection && max(abs(channel.H_reflect_f(:))) > 1e-12
                error('Direct-only reflected response is nonzero for %s seed %d.', ...
                    sea_conditions(cc).name, seed_list(mm));
            end

            [~, H_baseband_shifted] = local_build_baseband_response( ...
                channel.f_axis, channel.H_f, channel.idx_f_ref, ...
                comm_cfg.symbol_rate_hz, comm_cfg.n_sym);
            [h_bb, n_tap_eff, energy_kept] = local_build_channel_taps( ...
                H_baseband_shifted, comm_cfg.tap_energy_ratio);
            tap_metrics = local_tap_metrics(h_bb, n_tap_eff, energy_kept);
            [rx_clean, h_eq, receive_meta] = local_apply_channel_window( ...
                tx_symbols, h_bb, comm_cfg.receive_window_mode);
            noise_signal_ref = local_select_noise_reference( ...
                comm_cfg.ebn0_reference, tx_symbols, rx_clean);

            ber = NaN(n_eb, 1);
            ser = NaN(n_eb, 1);
            effective_snr_db = NaN(n_eb, 1);
            for ee = 1:n_eb
                noise_cfg = struct();
                noise_cfg.enable_noise = true;
                noise_cfg.model = comm_cfg.noise_model;
                noise_cfg.ebn0_db = comm_cfg.EbN0_dB_list(ee);
                noise_cfg.bits_per_symbol = bits_per_symbol;
                noise_cfg.seed = comm_cfg.noise_seed_base + cc*100000 + ss*10000 + mm*100 + ee;
                noise_cfg.custom_noise_fn = [];

                [rx_noisy, ~, noise_meta] = noise_inject_vertical(rx_clean, noise_cfg, noise_signal_ref);
                rx_eq = local_mmse_equalize(rx_noisy, h_eq, comm_cfg.mmse_reg_eps);
                bits_rx = modem_psk('demodulate', rx_eq, comm_cfg.M);
                [ber(ee), ser(ee)] = modem_psk('error_rate', bits_tx, bits_rx, comm_cfg.M);
                effective_snr_db(ee) = noise_meta.effective_snr_db;
            end

            ber_samples(:, mm, ss, cc) = ber;
            ser_samples(:, mm, ss, cc) = ser;
            effective_snr_samples(:, mm, ss, cc) = effective_snr_db;

            row = struct();
            row.condition_index = cc;
            row.condition_name = string(sea_conditions(cc).name);
            row.scenario_index = ss;
            row.scenario_name = string(scenarios(ss).name);
            row.seed_index = mm;
            row.sea_seed = seed_list(mm);
            row.sea_hs_target = sea_conditions(cc).sea_hs_target;
            row.sea_wind_speed = sea_conditions(cc).sea_wind_speed;
            row.enable_surface_reflection = scenarios(ss).enable_surface_reflection;
            row.abs_h_total = abs(channel.h_total);
            row.phase_h_total_rad = angle(channel.h_total);
            row.abs_h_direct = abs(channel.h_direct);
            row.abs_h_reflect = abs(channel.h_reflect);
            row.h_bb_tap_count = n_tap_eff;
            row.h_eq_tap_count = numel(h_eq);
            row.h_bb_energy_kept = energy_kept;
            row.receive_peak_index_original = receive_meta.peak_index_original;
            row.receive_discarded_pre_peak_energy_fraction = receive_meta.discarded_pre_peak_energy_fraction;
            row.tap_rms_delay_symbols = tap_metrics.tap_rms_delay_symbols;
            row.tap_peak_fraction = tap_metrics.tap_peak_fraction;
            row.invariant_error = invariant_error;
            row.max_abs_H_reflect_f = max(abs(channel.H_reflect_f(:)));
            row.BER = ber.';
            row.SER = ser.';
            row.effective_snr_db = effective_snr_db.';
            run_rows = local_append_struct(run_rows, row);
        end
    end
end

run_summary_table = struct2table(run_rows);
curve_stats = local_build_curve_stats( ...
    ber_samples, ser_samples, effective_snr_samples, comm_cfg.EbN0_dB_list, ...
    sea_conditions, scenarios);
curve_summary_table = local_build_curve_summary_table(curve_stats, sea_conditions, scenarios);

disp(curve_summary_table(:, {'condition_name', 'scenario_name', 'EbN0_dB', ...
    'BER_mean', 'BER_std', 'BER_q25', 'BER_q75', ...
    'SER_mean', 'SER_std', 'SER_q25', 'SER_q75'}))

local_plot_condition_compare(curve_stats, sea_conditions, scenarios, 'BER', ...
    'direct_plus_reflect', 'BER: representative direct+reflect', ...
    figure_prefix, 'BER_representative_direct_plus_reflect');
local_plot_condition_compare(curve_stats, sea_conditions, scenarios, 'SER', ...
    'direct_plus_reflect', 'SER: representative direct+reflect', ...
    figure_prefix, 'SER_representative_direct_plus_reflect');
for cc = 1:n_cond
    local_plot_scenario_compare(curve_stats, sea_conditions, scenarios, 'BER', ...
        sea_conditions(cc).name, ['BER: ' sea_conditions(cc).name ' direct-only vs direct+reflect'], ...
        figure_prefix, ['BER_' sea_conditions(cc).name '_direct_vs_reflect']);
    local_plot_scenario_compare(curve_stats, sea_conditions, scenarios, 'SER', ...
        sea_conditions(cc).name, ['SER: ' sea_conditions(cc).name ' direct-only vs direct+reflect'], ...
        figure_prefix, ['SER_' sea_conditions(cc).name '_direct_vs_reflect']);
end

save(result_file, 'params_base', 'comm_cfg', 'sea_conditions', 'scenarios', ...
    'seed_list', 'bits_tx', 'ber_samples', 'ser_samples', ...
    'effective_snr_samples', 'run_summary_table', ...
    'curve_summary_table', 'curve_stats');

function paramsV = local_base_params()
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
paramsV.sea_wind_speed = 5.0;
paramsV.sea_hs_target = 0.5;
paramsV.sea_seed = 12345;
paramsV.surface_boundary_model = 'kirchhoff_kdomain';
paramsV.surface_boundary_coupling_diagnostics = false;
paramsV.surface_boundary_redistribution_diagnostics = false;
end

function comm_cfg = local_comm_cfg()
comm_cfg = struct();
comm_cfg.M = 4;
comm_cfg.n_sym = 1000;
comm_cfg.EbN0_dB_list = 0:2:20;
comm_cfg.symbol_rate_hz = 1000;
comm_cfg.tap_energy_ratio = 0.999;
comm_cfg.bits_seed = 9000;
comm_cfg.noise_seed_base = 7000;
comm_cfg.noise_model = 'awgn';
comm_cfg.mmse_reg_eps = 1e-6;
comm_cfg.receive_window_mode = 'peak_sync';
comm_cfg.ebn0_reference = 'rx_clean';
end

function sea_conditions = local_representative_sea_conditions()
sea_conditions = struct([]);
sea_conditions(1).name = 'weak';
sea_conditions(1).sea_hs_target = 0.05;
sea_conditions(1).sea_wind_speed = 5;
sea_conditions(2).name = 'mid';
sea_conditions(2).sea_hs_target = 0.5;
sea_conditions(2).sea_wind_speed = 8;
sea_conditions(3).name = 'strong';
sea_conditions(3).sea_hs_target = 1.0;
sea_conditions(3).sea_wind_speed = 12;
end

function scenarios = local_scenarios()
scenarios = struct([]);
scenarios(1).name = 'direct_only';
scenarios(1).enable_surface_reflection = false;
scenarios(2).name = 'direct_plus_reflect';
scenarios(2).enable_surface_reflection = true;
end

function out = local_append_struct(out, row)
if isempty(out)
    out = row;
else
    out(end + 1) = row; %#ok<AGROW>
end
end

function [rx_clean, h_eq, meta] = local_apply_channel_window(tx_symbols, h_bb, receive_window_mode)
tx_symbols = tx_symbols(:);
h_bb = h_bb(:);
if nargin < 3 || isempty(receive_window_mode)
    receive_window_mode = 'peak_sync';
end
if isstring(receive_window_mode)
    receive_window_mode = char(receive_window_mode);
end

tap_energy = abs(h_bb).^2;
if all(tap_energy == 0)
    peak_index = 1;
else
    [~, peak_index] = max(tap_energy);
end

switch lower(receive_window_mode)
    case 'same_legacy'
        rx_clean = conv(tx_symbols, h_bb, 'same');
        h_eq = h_bb;
        start_index = NaN;
        discarded_energy_fraction = 0;
    case 'causal_head'
        rx_full = conv(tx_symbols, h_bb, 'full');
        rx_clean = rx_full(1:numel(tx_symbols));
        h_eq = h_bb;
        start_index = 1;
        discarded_energy_fraction = 0;
    case 'peak_sync'
        h_eq = h_bb(peak_index:end);
        rx_full = conv(tx_symbols, h_eq, 'full');
        rx_clean = rx_full(1:numel(tx_symbols));
        start_index = peak_index;
        discarded_energy_fraction = sum(tap_energy(1:max(peak_index - 1, 0))) / max(sum(tap_energy), eps);
    otherwise
        error('Unsupported receive_window_mode: %s', receive_window_mode);
end

meta = struct();
meta.mode = lower(receive_window_mode);
meta.peak_index_original = peak_index;
meta.start_index_original = start_index;
meta.original_tap_count = numel(h_bb);
meta.equalizer_tap_count = numel(h_eq);
meta.discarded_pre_peak_energy_fraction = discarded_energy_fraction;
end

function signal_ref = local_select_noise_reference(ebn0_reference, tx_symbols, rx_clean)
if nargin < 1 || isempty(ebn0_reference)
    ebn0_reference = 'rx_clean';
end
if isstring(ebn0_reference)
    ebn0_reference = char(ebn0_reference);
end
switch lower(ebn0_reference)
    case 'rx_clean'
        signal_ref = rx_clean;
    case 'tx_symbols'
        signal_ref = tx_symbols;
    otherwise
        error('Unsupported Eb/N0 reference: %s', ebn0_reference);
end
end

function [f_bb_axis, H_baseband_shifted] = local_build_baseband_response(f_axis, H_f, idx_f_ref, fs_hz, n_fft)
f_axis = f_axis(:);
H_f = H_f(:);
if numel(f_axis) ~= numel(H_f)
    error('f_axis and H_f length mismatch.');
end
if idx_f_ref < 1 || idx_f_ref > numel(f_axis)
    error('idx_f_ref out of range.');
end

fc = f_axis(idx_f_ref);
f_rel = f_axis - fc;
f_bb_axis = ((0:n_fft-1).' - floor(n_fft/2)) * (fs_hz / n_fft);
H_baseband_shifted = interp1(f_rel, H_f, f_bb_axis, 'linear', 0);
end

function [h_eff, n_tap_eff, energy_kept] = local_build_channel_taps(H_baseband_shifted, energy_ratio)
H_baseband_shifted = H_baseband_shifted(:);
h_full = ifft(ifftshift(H_baseband_shifted));
energy = abs(h_full).^2;
if all(energy == 0)
    h_eff = complex(0, 0);
    n_tap_eff = 1;
    energy_kept = 0;
    return
end

cum_energy = cumsum(energy);
target = max(min(energy_ratio, 1), 0);
n_tap_eff = find(cum_energy >= target * cum_energy(end), 1, 'first');
n_tap_eff = max(1, n_tap_eff);
h_eff = h_full(1:n_tap_eff);
energy_kept = sum(abs(h_eff).^2) / sum(abs(h_full).^2);
end

function metrics = local_tap_metrics(h_bb, tap_count, energy_kept)
h_bb = h_bb(:);
tap_energy = abs(h_bb).^2;
total_energy = sum(tap_energy);
if total_energy <= 0
    mean_delay = NaN;
    rms_delay = NaN;
    peak_index = NaN;
    peak_fraction = NaN;
else
    delay_idx = (0:(numel(h_bb) - 1)).';
    mean_delay = sum(delay_idx .* tap_energy) / total_energy;
    rms_delay = sqrt(sum(((delay_idx - mean_delay).^2) .* tap_energy) / total_energy);
    [peak_energy, peak_index] = max(tap_energy);
    peak_fraction = peak_energy / total_energy;
end

metrics = struct( ...
    'tap_count', tap_count, ...
    'tap_energy_kept', energy_kept, ...
    'tap_mean_delay_symbols', mean_delay, ...
    'tap_rms_delay_symbols', rms_delay, ...
    'tap_peak_index', peak_index, ...
    'tap_peak_fraction', peak_fraction);
end

function rx_eq = local_mmse_equalize(rx_noisy, h_taps, reg_eps)
rx_noisy = rx_noisy(:);
h_taps = h_taps(:);
N = numel(rx_noisy);
L = min(numel(h_taps), N);
h_pad = [h_taps(1:L); zeros(N - L, 1)];
H = fft(h_pad);
W = conj(H) ./ (abs(H).^2 + reg_eps);
rx_eq = ifft(fft(rx_noisy) .* W);
end

function curve_stats = local_build_curve_stats( ...
    ber_samples, ser_samples, effective_snr_samples, EbN0_dB_list, sea_conditions, scenarios)

curve_stats = struct();
curve_stats.EbN0_dB_list = EbN0_dB_list(:);
curve_stats.condition_names = string({sea_conditions.name});
curve_stats.scenario_names = string({scenarios.name});
curve_stats.BER = local_sample_stats(ber_samples);
curve_stats.SER = local_sample_stats(ser_samples);
curve_stats.effective_snr_db = local_sample_stats(effective_snr_samples);
curve_stats.validation = struct();
curve_stats.validation.BER_nonmonotonic_count = local_nonmonotonic_counts(curve_stats.BER.mean);
curve_stats.validation.SER_nonmonotonic_count = local_nonmonotonic_counts(curve_stats.SER.mean);
end

function stats = local_sample_stats(samples)
stats = struct();
mu = mean(samples, 2, 'omitnan');
stats.samples = samples;
stats.mean = squeeze(mu);
stats.var = squeeze(mean((samples - mu).^2, 2, 'omitnan'));
stats.std = squeeze(std(samples, 0, 2, 'omitnan'));
stats.quantiles_5_25_50_75_95 = local_quantiles_dim(samples, [5, 25, 50, 75, 95]);
end

function q = local_quantiles_dim(samples, pct)
[n_eb, ~, n_scen, n_cond] = size(samples);
q = NaN(n_eb, numel(pct), n_scen, n_cond);
for ee = 1:n_eb
    for ss = 1:n_scen
        for cc = 1:n_cond
            x = squeeze(samples(ee, :, ss, cc));
            q(ee, :, ss, cc) = local_quantiles(x, pct);
        end
    end
end
end

function q = local_quantiles(x, pct)
x = sort(x(:));
x = x(isfinite(x));
pct = pct(:).';
if isempty(x)
    q = NaN(size(pct));
elseif numel(x) == 1
    q = repmat(x, size(pct));
else
    pos = 1 + (pct / 100) * (numel(x) - 1);
    q = interp1(1:numel(x), x, pos, 'linear');
end
end

function counts = local_nonmonotonic_counts(mean_curves)
[~, n_scen, n_cond] = size(mean_curves);
counts = zeros(n_scen, n_cond);
for ss = 1:n_scen
    for cc = 1:n_cond
        y = mean_curves(:, ss, cc);
        counts(ss, cc) = sum(diff(y) > 1e-12);
    end
end
end

function curve_summary_table = local_build_curve_summary_table(curve_stats, sea_conditions, scenarios)
rows = struct([]);
pct_idx = struct('q5', 1, 'q25', 2, 'q50', 3, 'q75', 4, 'q95', 5);
for cc = 1:numel(sea_conditions)
    for ss = 1:numel(scenarios)
        for ee = 1:numel(curve_stats.EbN0_dB_list)
            row = struct();
            row.condition_index = cc;
            row.condition_name = string(sea_conditions(cc).name);
            row.sea_hs_target = sea_conditions(cc).sea_hs_target;
            row.sea_wind_speed = sea_conditions(cc).sea_wind_speed;
            row.scenario_index = ss;
            row.scenario_name = string(scenarios(ss).name);
            row.enable_surface_reflection = scenarios(ss).enable_surface_reflection;
            row.EbN0_dB = curve_stats.EbN0_dB_list(ee);
            row.BER_mean = curve_stats.BER.mean(ee, ss, cc);
            row.BER_var = curve_stats.BER.var(ee, ss, cc);
            row.BER_std = curve_stats.BER.std(ee, ss, cc);
            row.BER_q5 = curve_stats.BER.quantiles_5_25_50_75_95(ee, pct_idx.q5, ss, cc);
            row.BER_q25 = curve_stats.BER.quantiles_5_25_50_75_95(ee, pct_idx.q25, ss, cc);
            row.BER_q50 = curve_stats.BER.quantiles_5_25_50_75_95(ee, pct_idx.q50, ss, cc);
            row.BER_q75 = curve_stats.BER.quantiles_5_25_50_75_95(ee, pct_idx.q75, ss, cc);
            row.BER_q95 = curve_stats.BER.quantiles_5_25_50_75_95(ee, pct_idx.q95, ss, cc);
            row.SER_mean = curve_stats.SER.mean(ee, ss, cc);
            row.SER_var = curve_stats.SER.var(ee, ss, cc);
            row.SER_std = curve_stats.SER.std(ee, ss, cc);
            row.SER_q5 = curve_stats.SER.quantiles_5_25_50_75_95(ee, pct_idx.q5, ss, cc);
            row.SER_q25 = curve_stats.SER.quantiles_5_25_50_75_95(ee, pct_idx.q25, ss, cc);
            row.SER_q50 = curve_stats.SER.quantiles_5_25_50_75_95(ee, pct_idx.q50, ss, cc);
            row.SER_q75 = curve_stats.SER.quantiles_5_25_50_75_95(ee, pct_idx.q75, ss, cc);
            row.SER_q95 = curve_stats.SER.quantiles_5_25_50_75_95(ee, pct_idx.q95, ss, cc);
            row.effective_snr_db_mean = curve_stats.effective_snr_db.mean(ee, ss, cc);
            row.effective_snr_db_std = curve_stats.effective_snr_db.std(ee, ss, cc);
            rows = local_append_struct(rows, row);
        end
    end
end
curve_summary_table = struct2table(rows);
end

function local_plot_condition_compare(curve_stats, sea_conditions, scenarios, metric_name, scenario_name, title_text, figure_prefix, suffix)
scenario_idx = local_find_name({scenarios.name}, scenario_name);
if isempty(scenario_idx)
    return
end
figure('Visible', 'off');
hold on
for cc = 1:numel(sea_conditions)
    local_plot_mean_with_iqr(curve_stats, metric_name, scenario_idx, cc, sea_conditions(cc).name);
end
grid on
xlabel('Eb/N0 (dB)')
ylabel(metric_name)
title(title_text, 'Interpreter', 'none')
legend('Location', 'southwest')
set(gca, 'YScale', 'log')
print(gcf, '-dpng', '-r200', [figure_prefix suffix '.png'])
close(gcf)
end

function local_plot_scenario_compare(curve_stats, sea_conditions, scenarios, metric_name, condition_name, title_text, figure_prefix, suffix)
condition_idx = local_find_name({sea_conditions.name}, condition_name);
if isempty(condition_idx)
    return
end
figure('Visible', 'off');
hold on
for ss = 1:numel(scenarios)
    local_plot_mean_with_iqr(curve_stats, metric_name, ss, condition_idx, scenarios(ss).name);
end
grid on
xlabel('Eb/N0 (dB)')
ylabel(metric_name)
title(title_text, 'Interpreter', 'none')
legend('Location', 'southwest')
set(gca, 'YScale', 'log')
print(gcf, '-dpng', '-r200', [figure_prefix suffix '.png'])
close(gcf)
end

function local_plot_mean_with_iqr(curve_stats, metric_name, scenario_idx, condition_idx, label_text)
stats = curve_stats.(metric_name);
mean_y = stats.mean(:, scenario_idx, condition_idx);
q25 = stats.quantiles_5_25_50_75_95(:, 2, scenario_idx, condition_idx);
q75 = stats.quantiles_5_25_50_75_95(:, 4, scenario_idx, condition_idx);
x = curve_stats.EbN0_dB_list;
plot_floor = 1e-5;
plot(x, max(mean_y, plot_floor), 'o-', 'LineWidth', 1.2, ...
    'DisplayName', [char(label_text) ' mean'])
plot(x, max(q25, plot_floor), '--', 'LineWidth', 0.8, ...
    'HandleVisibility', 'off')
plot(x, max(q75, plot_floor), '--', 'LineWidth', 0.8, ...
    'HandleVisibility', 'off')
end

function idx = local_find_name(names, target_name)
idx = [];
for ii = 1:numel(names)
    if strcmpi(names{ii}, target_name)
        idx = ii;
        return
    end
end
end

