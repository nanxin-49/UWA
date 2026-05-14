% Compare QPSK BER/SER across vertical channel bubble scenarios.
% Fairness: all scenarios use the same bits_tx. For each Eb/N0 index, all
% scenarios reuse the same AWGN seed so channel differences dominate.

clear
format compact

result_file = 'comm_compare_bubble_models_vertical_result.mat';
figure_prefix = 'comm_compare_bubble_models_vertical_';
if exist(result_file, 'file')
    warning('Result file %s already exists and will be overwritten.', result_file);
end

M = 4;
n_sym = 2000;
EbN0_dB_list = 0:2:20;
k = log2(M);
symbol_rate_hz = 1000;
isi_mode = 'linear_conv';

noise_control = struct();
noise_control.enable_noise = true;
noise_control.model = 'awgn';
noise_control.seed_base = 7000;

rng(20240514, 'twister');
bits_tx = randi([0, 1], n_sym*k, 1);
[tx_symbols, bits_used] = modem_psk('modulate', bits_tx, M);
bits_tx = bits_used;

params_base = local_base_channel_params();
scenarios = local_scenarios();

results = struct([]);
for ss = 1:numel(scenarios)
    paramsV = local_apply_scenario(params_base, scenarios(ss));
    fprintf('Running channel/comm scenario %s\n', scenarios(ss).name);
    channel = CARPE3D_vertical(paramsV);
    invariant_error = norm(channel.H_f(:) - channel.H_direct_f(:) - channel.H_reflect_f(:));
    if invariant_error > 1e-10
        error('comm_compare_bubble_models_vertical:Invariant', ...
            'H_f invariant failed for %s: %.3e', scenarios(ss).name, invariant_error);
    end

    if ~strcmpi(isi_mode, 'linear_conv')
        error('Unsupported isi_mode: %s', isi_mode);
    end

    [f_bb_axis, H_baseband] = local_build_baseband_response( ...
        channel.f_axis, channel.H_f, channel.idx_f_ref, symbol_rate_hz, n_sym);
    [h_bb, n_tap_eff, energy_kept] = local_build_channel_taps(H_baseband, 0.999);
    rx_clean = conv(tx_symbols, h_bb, 'same');

    BER = zeros(numel(EbN0_dB_list), 1);
    SER = zeros(numel(EbN0_dB_list), 1);
    effective_snr_db = zeros(numel(EbN0_dB_list), 1);
    for ii = 1:numel(EbN0_dB_list)
        noise_cfg = struct();
        noise_cfg.enable_noise = noise_control.enable_noise;
        noise_cfg.model = noise_control.model;
        noise_cfg.ebn0_db = EbN0_dB_list(ii);
        noise_cfg.bits_per_symbol = k;
        noise_cfg.seed = noise_control.seed_base + ii;
        noise_cfg.custom_noise_fn = [];

        [rx_noisy, ~, noise_meta] = noise_inject_vertical(rx_clean, noise_cfg, tx_symbols);
        rx_eq = local_mmse_equalize(rx_noisy, h_bb, 1e-6);
        bits_rx = modem_psk('demodulate', rx_eq, M);
        [BER(ii), SER(ii)] = modem_psk('error_rate', bits_tx, bits_rx, M);
        effective_snr_db(ii) = noise_meta.effective_snr_db;
    end

    if any(~isfinite(BER)) || any(~isfinite(SER)) || any(~isfinite(effective_snr_db))
        error('Non-finite communication metric detected for %s.', scenarios(ss).name);
    end

    results(ss).name = scenarios(ss).name; %#ok<SAGROW>
    results(ss).paramsV = paramsV;
    results(ss).channel = channel;
    results(ss).f_axis = channel.f_axis;
    results(ss).H_f = channel.H_f;
    results(ss).H_direct_f = channel.H_direct_f;
    results(ss).H_reflect_f = channel.H_reflect_f;
    results(ss).h_total = channel.h_total;
    results(ss).f_bb_axis = f_bb_axis;
    results(ss).H_baseband = H_baseband;
    results(ss).h_bb = h_bb;
    results(ss).h_bb_tap_count = n_tap_eff;
    results(ss).h_bb_energy_kept = energy_kept;
    results(ss).EbN0_dB_list = EbN0_dB_list;
    results(ss).BER = BER;
    results(ss).SER = SER;
    results(ss).effective_snr_db = effective_snr_db;
    results(ss).invariant_error = invariant_error;
    results(ss).bubble_meta = channel.bubble_meta;
end

results = local_attach_channel_metrics(results);
summary_table = local_summary_table(results, EbN0_dB_list);
disp(summary_table)

local_plot_ber(results, figure_prefix);
local_plot_ser(results, figure_prefix);
local_plot_channel_magnitude(results, figure_prefix);
local_plot_channel_phase(results, figure_prefix);
local_plot_taps(results, figure_prefix);
local_plot_delta_tl(results, figure_prefix);

save(result_file, 'params_base', 'scenarios', 'results', ...
    'noise_control', 'EbN0_dB_list', 'M', 'n_sym', 'symbol_rate_hz', ...
    'isi_mode', 'bits_tx', 'summary_table');

function paramsV = local_base_channel_params()
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
paramsV.nout = 6;
% Current CARPE3D validation requires sigma_src_m >= max(dx,dy).
% With nx=ny=128 and xw=yw=50, dx=dy=0.390625 m, so use 0.4 m.
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
end

function scenarios = local_scenarios()
scenarios = struct([]);
scenarios(1).name = 'no_bubble';
scenarios(1).enable_bubbles = false;
scenarios(1).bubble_model = 'off';
scenarios(1).bubble_spatial_mode = 'none';

scenarios(2).name = 'level0_empirical';
scenarios(2).enable_bubbles = true;
scenarios(2).bubble_model = 'level0_empirical';
scenarios(2).bubble_spatial_mode = '1d';
scenarios(2).bubble_alpha0_np_per_m = 0.02;
scenarios(2).bubble_layer_decay_m = 20;
scenarios(2).bubble_delta_c0_mps = 0;

scenarios(3).name = 'hall1d_default';
scenarios(3).enable_bubbles = true;
scenarios(3).bubble_model = 'hall1d';
scenarios(3).bubble_spatial_mode = '1d';
scenarios(3).bubble_strength_scale = 1;
scenarios(3).bubble_wind_speed = [];

scenarios(4).name = 'hall1d_calibrated_3dB';
scenarios(4).enable_bubbles = true;
scenarios(4).bubble_model = 'hall1d';
scenarios(4).bubble_spatial_mode = '1d';
scenarios(4).sea_wind_speed = 5.0;
scenarios(4).bubble_wind_speed = 8;
scenarios(4).bubble_strength_scale = 1e2;

scenarios(5).name = 'hall1d_high_bubble_wind';
scenarios(5).enable_bubbles = true;
scenarios(5).bubble_model = 'hall1d';
scenarios(5).bubble_spatial_mode = '1d';
scenarios(5).sea_wind_speed = 5.0;
scenarios(5).bubble_wind_speed = 12;
scenarios(5).bubble_strength_scale = 1;
end

function paramsV = local_apply_scenario(params_base, scenario)
paramsV = params_base;
fields = fieldnames(scenario);
for ii = 1:numel(fields)
    if strcmp(fields{ii}, 'name')
        continue
    end
    if isempty(scenario.(fields{ii}))
        continue
    end
    paramsV.(fields{ii}) = scenario.(fields{ii});
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

function results = local_attach_channel_metrics(results)
H0 = results(1).H_f(:);
for ss = 1:numel(results)
    H = results(ss).H_f(:);
    if numel(H) ~= numel(H0)
        error('H_f length mismatch against no_bubble baseline.');
    end
    results(ss).delta_TL_dB = -20*log10(abs(H) ./ max(abs(H0), eps));
    results(ss).phase_diff_rad = unwrap(angle(H)) - unwrap(angle(H0));
    results(ss).max_delta_TL_dB = max(results(ss).delta_TL_dB);
    results(ss).max_abs_phase_diff_rad = max(abs(results(ss).phase_diff_rad));
end
end

function summary_table = local_summary_table(results, EbN0_dB_list)
idx0 = find(EbN0_dB_list == 0, 1);
idx10 = find(EbN0_dB_list == 10, 1);
idx20 = find(EbN0_dB_list == 20, 1);
rows = struct([]);
for ss = 1:numel(results)
    rows(ss).scenario = {results(ss).name}; %#ok<AGROW>
    rows(ss).h_total_abs = abs(results(ss).h_total);
    rows(ss).h_bb_tap_count = results(ss).h_bb_tap_count;
    rows(ss).max_delta_TL_dB = results(ss).max_delta_TL_dB;
    rows(ss).max_abs_phase_diff_rad = results(ss).max_abs_phase_diff_rad;
    rows(ss).BER_0dB = results(ss).BER(idx0);
    rows(ss).BER_10dB = results(ss).BER(idx10);
    rows(ss).BER_20dB = results(ss).BER(idx20);
    rows(ss).SER_0dB = results(ss).SER(idx0);
    rows(ss).SER_10dB = results(ss).SER(idx10);
    rows(ss).SER_20dB = results(ss).SER(idx20);
    rows(ss).invariant_error = results(ss).invariant_error;
end
summary_table = struct2table(rows);
end

function local_plot_ber(results, figure_prefix)
fig = figure('Visible', 'off');
hold on
for ss = 1:numel(results)
    semilogy(results(ss).EbN0_dB_list, max(results(ss).BER, 1e-6), ...
        'o-', 'LineWidth', 1.2)
end
grid on
xlabel('Eb/N0 (dB)')
ylabel('BER')
title('QPSK BER across bubble channel scenarios')
legend({results.name}, 'Interpreter', 'none', 'Location', 'southwest')
print(fig, '-dpng', '-r200', [figure_prefix 'BER_vs_EbN0.png'])
close(fig)
end

function local_plot_ser(results, figure_prefix)
fig = figure('Visible', 'off');
hold on
for ss = 1:numel(results)
    semilogy(results(ss).EbN0_dB_list, max(results(ss).SER, 1e-6), ...
        's-', 'LineWidth', 1.2)
end
grid on
xlabel('Eb/N0 (dB)')
ylabel('SER')
title('QPSK SER across bubble channel scenarios')
legend({results.name}, 'Interpreter', 'none', 'Location', 'southwest')
print(fig, '-dpng', '-r200', [figure_prefix 'SER_vs_EbN0.png'])
close(fig)
end

function local_plot_channel_magnitude(results, figure_prefix)
fig = figure('Visible', 'off');
hold on
for ss = 1:numel(results)
    plot(results(ss).f_axis, 20*log10(max(abs(results(ss).H_f), eps)), ...
        'LineWidth', 1.2)
end
grid on
xlabel('Frequency (Hz)')
ylabel('|H(f)| (dB)')
title('Channel magnitude')
legend({results.name}, 'Interpreter', 'none', 'Location', 'best')
print(fig, '-dpng', '-r200', [figure_prefix 'H_magnitude.png'])
close(fig)
end

function local_plot_channel_phase(results, figure_prefix)
fig = figure('Visible', 'off');
hold on
for ss = 1:numel(results)
    plot(results(ss).f_axis, unwrap(angle(results(ss).H_f)), 'LineWidth', 1.2)
end
grid on
xlabel('Frequency (Hz)')
ylabel('Unwrapped phase (rad)')
title('Channel phase')
legend({results.name}, 'Interpreter', 'none', 'Location', 'best')
print(fig, '-dpng', '-r200', [figure_prefix 'H_phase.png'])
close(fig)
end

function local_plot_taps(results, figure_prefix)
fig = figure('Visible', 'off');
hold on
for ss = 1:numel(results)
    stem(0:(numel(results(ss).h_bb)-1), abs(results(ss).h_bb), ...
        'DisplayName', results(ss).name, 'LineWidth', 1.0)
end
grid on
xlabel('Tap index')
ylabel('|h_{bb}|')
title('Baseband tap magnitudes')
legend('Interpreter', 'none', 'Location', 'best')
print(fig, '-dpng', '-r200', [figure_prefix 'h_bb_taps.png'])
close(fig)
end

function local_plot_delta_tl(results, figure_prefix)
fig = figure('Visible', 'off');
hold on
for ss = 1:numel(results)
    plot(results(ss).f_axis, results(ss).delta_TL_dB, 'LineWidth', 1.2)
end
grid on
xlabel('Frequency (Hz)')
ylabel('\DeltaTL relative to no bubble (dB)')
title('Bubble excess TL')
legend({results.name}, 'Interpreter', 'none', 'Location', 'best')
print(fig, '-dpng', '-r200', [figure_prefix 'delta_TL.png'])
close(fig)
end
