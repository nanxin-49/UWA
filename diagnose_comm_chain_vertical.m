% Diagnose the vertical MPSK communication chain before using BER/SER results.
% This script does not modify propagation or bubble physics. It isolates:
% 1) ideal flat-channel modem/noise behavior;
% 2) scalar PE h_total behavior;
% 3) wideband H_f -> h_bb timing/tap-spread behavior;
% 4) main-tap circular alignment;
% 5) ideal training-aided delay/phase compensation.

clear
format compact

result_file = 'diagnose_comm_chain_vertical_result.mat';
if exist(result_file, 'file')
    warning('Result file %s already exists and will be overwritten.', result_file);
end

M = 4;
k = log2(M);
n_sym = 2000;
EbN0_dB_list = 0:2:20;
symbol_rate_hz = 1000;
noise_seed_base = 7000;

rng(20240514, 'twister');
bits_tx = randi([0, 1], n_sym*k, 1);
[tx_symbols, bits_tx] = modem_psk('modulate', bits_tx, M);

paramsV = local_base_params();
channel = vertical_channel_model(paramsV);
invariant_error = norm(channel.H_f(:) - channel.H_direct_f(:) - channel.H_reflect_f(:));
if invariant_error > 1e-10
    error('H_f invariant failed: %.3e', invariant_error);
end

[f_bb_axis, H_baseband] = local_build_baseband_response( ...
    channel.f_axis, channel.H_f, channel.idx_f_ref, symbol_rate_hz, n_sym);
[h_original, tap_count_original, energy_original, h_full] = ...
    local_build_channel_taps(H_baseband, 0.999, false);
[h_aligned, tap_count_aligned, energy_aligned, h_aligned_full, main_tap_index] = ...
    local_build_aligned_channel_taps(H_baseband, 0.999);

h_diag = local_hbb_diagnostics(f_bb_axis, H_baseband, h_full, h_aligned_full, ...
    tap_count_original, tap_count_aligned, energy_original, energy_aligned, ...
    main_tap_index, symbol_rate_hz);
disp('Wideband H_f -> h_bb diagnostics:')
disp(struct2table(h_diag, 'AsArray', true))

ideal_table = local_eval_flat(tx_symbols, bits_tx, M, EbN0_dB_list, noise_seed_base);
scalar_table = local_eval_scalar(tx_symbols, bits_tx, M, EbN0_dB_list, noise_seed_base, channel.h_total);
original_table = local_eval_wideband(tx_symbols, bits_tx, M, EbN0_dB_list, ...
    noise_seed_base, h_original, false);
aligned_table = local_eval_wideband(tx_symbols, bits_tx, M, EbN0_dB_list, ...
    noise_seed_base, h_aligned, false);
aligned_comp_table = local_eval_wideband(tx_symbols, bits_tx, M, EbN0_dB_list, ...
    noise_seed_base, h_aligned, true);

disp('Ideal flat channel:')
disp(ideal_table)
disp('Scalar PE channel using h_total:')
disp(scalar_table)
disp('Original wideband h_bb:')
disp(original_table)
disp('Aligned wideband h_bb:')
disp(aligned_table)
disp('Aligned wideband h_bb with ideal training-aided delay/phase compensation:')
disp(aligned_comp_table)

save(result_file, 'paramsV', 'channel', 'invariant_error', 'f_bb_axis', ...
    'H_baseband', 'h_full', 'h_original', 'h_aligned_full', 'h_aligned', ...
    'h_diag', 'ideal_table', 'scalar_table', 'original_table', ...
    'aligned_table', 'aligned_comp_table', 'EbN0_dB_list', 'M', 'n_sym');

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
paramsV.nout = 6;
% Current CARPE3D validation requires sigma_src_m >= max(dx,dy).
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
paramsV.enable_bubbles = false;
end

function [f_bb_axis, H_baseband_shifted] = local_build_baseband_response(f_axis, H_f, idx_f_ref, fs_hz, n_fft)
f_axis = f_axis(:);
H_f = H_f(:);
fc = f_axis(idx_f_ref);
f_rel = f_axis - fc;
f_bb_axis = ((0:n_fft-1).' - floor(n_fft/2)) * (fs_hz / n_fft);
H_baseband_shifted = interp1(f_rel, H_f, f_bb_axis, 'linear', 0);
end

function [h_eff, n_tap_eff, energy_kept, h_full] = local_build_channel_taps(H_baseband_shifted, energy_ratio, align_main)
H_baseband_shifted = H_baseband_shifted(:);
h_full = ifft(ifftshift(H_baseband_shifted));
if align_main
    [~, i0] = max(abs(h_full));
    h_full = circshift(h_full, 1 - i0);
end
[h_eff, n_tap_eff, energy_kept] = local_truncate_taps(h_full, energy_ratio);
end

function [h_eff, n_tap_eff, energy_kept, h_aligned_full, main_tap_index] = local_build_aligned_channel_taps(H_baseband_shifted, energy_ratio)
[~, ~, ~, h_full] = local_build_channel_taps(H_baseband_shifted, energy_ratio, false);
[~, main_tap_index] = max(abs(h_full));
h_aligned_full = circshift(h_full, 1 - main_tap_index);
[h_eff, n_tap_eff, energy_kept] = local_truncate_taps(h_aligned_full, energy_ratio);
end

function [h_eff, n_tap_eff, energy_kept] = local_truncate_taps(h_full, energy_ratio)
energy = abs(h_full(:)).^2;
if all(energy == 0)
    h_eff = complex(0, 0);
    n_tap_eff = 1;
    energy_kept = 0;
    return
end
cum_energy = cumsum(energy);
n_tap_eff = find(cum_energy >= energy_ratio * cum_energy(end), 1, 'first');
n_tap_eff = max(1, n_tap_eff);
h_eff = h_full(1:n_tap_eff);
energy_kept = sum(abs(h_eff).^2) / sum(energy);
end

function diag_row = local_hbb_diagnostics(f_bb_axis, H_baseband, h_full, h_aligned_full, ...
    tap_count_original, tap_count_aligned, energy_original, energy_aligned_ratio, main_tap_index, fs_hz)

H_abs = abs(H_baseband(:));
nonzero_mask = H_abs > 0;
active_mask = H_abs > max(H_abs) * 1e-6;
phase = unwrap(angle(H_baseband(active_mask)));
freq = f_bb_axis(active_mask);
if numel(freq) >= 2
    p = polyfit(freq(:), phase(:), 1);
    group_delay_s = -p(1) / (2*pi);
else
    group_delay_s = NaN;
end

energy_full = abs(h_full(:)).^2;
energy_aligned_vec = abs(h_aligned_full(:)).^2;
diag_row = struct();
diag_row.H_nonzero_bins = nnz(nonzero_mask);
diag_row.H_active_bins = nnz(active_mask);
diag_row.H_abs_min = min(H_abs(nonzero_mask));
diag_row.H_abs_max = max(H_abs);
diag_row.H_abs_mean = mean(H_abs(nonzero_mask));
diag_row.h_full_main_tap_index = main_tap_index;
diag_row.h_full_main_tap_abs = max(abs(h_full));
diag_row.h_full_energy_first_10 = sum(energy_full(1:min(10,end))) / sum(energy_full);
diag_row.h_full_energy_first_100 = sum(energy_full(1:min(100,end))) / sum(energy_full);
diag_row.h_aligned_energy_first_10 = sum(energy_aligned_vec(1:min(10,end))) / sum(energy_aligned_vec);
diag_row.h_aligned_energy_first_100 = sum(energy_aligned_vec(1:min(100,end))) / sum(energy_aligned_vec);
diag_row.tap_count_original = tap_count_original;
diag_row.tap_count_aligned = tap_count_aligned;
diag_row.conv_same_delay_est = floor(tap_count_original / 2);
diag_row.energy_kept_original = energy_original;
diag_row.energy_kept_aligned = energy_aligned_ratio;
diag_row.group_delay_s = group_delay_s;
diag_row.group_delay_samples = group_delay_s * fs_hz;
diag_row.circular_delay_samples = mod(round(diag_row.group_delay_samples), numel(h_full));
end

function result_table = local_eval_flat(tx_symbols, bits_tx, M, EbN0_dB_list, seed_base)
BER = zeros(numel(EbN0_dB_list), 1);
SER = zeros(numel(EbN0_dB_list), 1);
effective_snr_db = zeros(numel(EbN0_dB_list), 1);
for ii = 1:numel(EbN0_dB_list)
    rx_clean = tx_symbols;
    [rx_noisy, meta] = local_add_noise(rx_clean, tx_symbols, M, EbN0_dB_list(ii), seed_base + ii);
    bits_rx = modem_psk('demodulate', rx_noisy, M);
    [BER(ii), SER(ii)] = modem_psk('error_rate', bits_tx, bits_rx, M);
    effective_snr_db(ii) = meta.effective_snr_db;
end
result_table = table(EbN0_dB_list(:), BER, SER, effective_snr_db, ...
    'VariableNames', {'EbN0_dB', 'BER', 'SER', 'EffectiveSNR_dB'});
end

function result_table = local_eval_scalar(tx_symbols, bits_tx, M, EbN0_dB_list, seed_base, h)
BER = zeros(numel(EbN0_dB_list), 1);
SER = zeros(numel(EbN0_dB_list), 1);
effective_snr_db = zeros(numel(EbN0_dB_list), 1);
for ii = 1:numel(EbN0_dB_list)
    rx_clean = tx_symbols .* h;
    [rx_noisy, meta] = local_add_noise(rx_clean, tx_symbols, M, EbN0_dB_list(ii), seed_base + ii);
    rx_eq = rx_noisy ./ h;
    bits_rx = modem_psk('demodulate', rx_eq, M);
    [BER(ii), SER(ii)] = modem_psk('error_rate', bits_tx, bits_rx, M);
    effective_snr_db(ii) = meta.effective_snr_db;
end
result_table = table(EbN0_dB_list(:), BER, SER, effective_snr_db, ...
    'VariableNames', {'EbN0_dB', 'BER', 'SER', 'EffectiveSNR_dB'});
end

function result_table = local_eval_wideband(tx_symbols, bits_tx, M, EbN0_dB_list, seed_base, h_taps, use_training_comp)
BER = zeros(numel(EbN0_dB_list), 1);
SER = zeros(numel(EbN0_dB_list), 1);
effective_snr_db = zeros(numel(EbN0_dB_list), 1);
best_delay = zeros(numel(EbN0_dB_list), 1);
gain_abs = ones(numel(EbN0_dB_list), 1);
for ii = 1:numel(EbN0_dB_list)
    rx_clean = conv(tx_symbols, h_taps, 'same');
    [rx_noisy, meta] = local_add_noise(rx_clean, tx_symbols, M, EbN0_dB_list(ii), seed_base + ii);
    rx_eq = local_mmse_equalize(rx_noisy, h_taps, 1e-6);
    if use_training_comp
        max_delay = min(numel(tx_symbols) - 1, max(100, numel(h_taps)));
        [rx_eq, best_delay(ii), g] = local_training_compensate(rx_eq, tx_symbols, max_delay);
        gain_abs(ii) = abs(g);
    end
    bits_rx = modem_psk('demodulate', rx_eq, M);
    [BER(ii), SER(ii)] = modem_psk('error_rate', bits_tx, bits_rx, M);
    effective_snr_db(ii) = meta.effective_snr_db;
end
result_table = table(EbN0_dB_list(:), BER, SER, effective_snr_db, best_delay, gain_abs, ...
    'VariableNames', {'EbN0_dB', 'BER', 'SER', 'EffectiveSNR_dB', 'BestDelay', 'CompGainAbs'});
end

function [rx_noisy, meta] = local_add_noise(rx_clean, signal_ref, M, EbN0_dB, seed)
noise_cfg = struct();
noise_cfg.enable_noise = true;
noise_cfg.model = 'awgn';
noise_cfg.ebn0_db = EbN0_dB;
noise_cfg.bits_per_symbol = log2(M);
noise_cfg.seed = seed;
noise_cfg.custom_noise_fn = [];
[rx_noisy, ~, meta] = noise_inject_vertical(rx_clean, noise_cfg, signal_ref);
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

function [rx_comp, best_delay, best_gain] = local_training_compensate(rx_eq, tx_symbols, max_delay)
rx_eq = rx_eq(:);
tx_symbols = tx_symbols(:);
best_mse = Inf;
best_delay = 0;
best_gain = 1;
best_rx = rx_eq;
for dd = -max_delay:max_delay
    [rx_seg, tx_seg] = local_delay_segments(rx_eq, tx_symbols, dd);
    if numel(rx_seg) < 10
        continue
    end
    gain = (rx_seg' * tx_seg) / max(rx_seg' * rx_seg, eps);
    err = tx_seg - gain .* rx_seg;
    mse = mean(abs(err).^2);
    if mse < best_mse
        best_mse = mse;
        best_delay = dd;
        best_gain = gain;
        best_rx = local_apply_delay(rx_eq, dd);
    end
end
rx_comp = best_gain .* best_rx;
end

function [rx_seg, tx_seg] = local_delay_segments(rx, tx, delay)
N = min(numel(rx), numel(tx));
if delay >= 0
    rx_idx = (1+delay):N;
    tx_idx = 1:(N-delay);
else
    rx_idx = 1:(N+delay);
    tx_idx = (1-delay):N;
end
rx_seg = rx(rx_idx);
tx_seg = tx(tx_idx);
end

function rx_shifted = local_apply_delay(rx, delay)
rx = rx(:);
rx_shifted = zeros(size(rx));
N = numel(rx);
if delay >= 0
    rx_shifted(1:(N-delay)) = rx((1+delay):N);
else
    d = -delay;
    rx_shifted((1+d):N) = rx(1:(N-d));
end
end
