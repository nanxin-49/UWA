% End-to-end MPSK communication demo over vertical PE channel.
% Scenario: seabed instrument TX -> hydrophone at fixed z_rx=3 m below buoy.

clear
format compact

% --- channel / geometry setup ---
paramsV = struct();
paramsV.f0 = 4000;
paramsV.enable_wideband = true;
paramsV.f_band_hz = [4000, 8000];
paramsV.Nf_min = 32;
paramsV.Nf_max = 64;
paramsV.f_ref_hz = 6000;
paramsV.c0 = 1500;
paramsV.z_max = 100;
paramsV.stepz_lamb = 0.5;

paramsV.xw = 50;
paramsV.yw = 50;
paramsV.nx = 256;
paramsV.ny = 256;

paramsV.x_tx = 0;
paramsV.y_tx = 0;
paramsV.z_tx = 100;
paramsV.x_rx = 0;
paramsV.y_rx = 0;
paramsV.z_rx = 3; % fixed 3 m for current stage

paramsV.rx_position_fn = []; % future interface: rx_position_fn(t_s, state) -> [x_rx,y_rx,z_rx]
paramsV.doppler_fn = [];     % future interface: doppler_fn(t_s, tx_state, rx_state, env_state) -> fd_hz

paramsV.nout = 6;
paramsV.sigma_src_m = 0.3;
paramsV.sponge_ratio = 0.12;
paramsV.alpha_max_np_per_m = 0.15;
paramsV.env_mode = 'uniform';
paramsV.enforce_1_over_R = true;
paramsV.show_figures = false;
paramsV.save_mode = 'rx_only';
paramsV.use_gpu = false;

paramsV.sea_wind_speed = 5.0;
paramsV.sea_hs_target = 0.5;
paramsV.sea_seed = 12345;

% --- communication setup ---
M = 4;                % QPSK (MPSK can be changed here)
n_sym = 2000;
EbN0_dB_list = 0:2:20;
k = log2(M);
symbol_rate_hz = 1000;
isi_mode = 'linear_conv';

noise_control = struct();
noise_control.enable_noise = true;
noise_control.model = 'awgn';
noise_control.seed_base = 7000;

bits_tx = randi([0, 1], n_sym*k, 1);
[tx_symbols, bits_used] = modem_psk('modulate', bits_tx, M);
bits_tx = bits_used;

scenarios = struct([]);
scenarios(1).name = 'direct_only';
scenarios(1).enable_surface_reflection = false;
scenarios(2).name = 'direct_plus_reflect';
scenarios(2).enable_surface_reflection = true;

results = struct([]);
for ss = 1:numel(scenarios)
    paramsV.enable_surface_reflection = scenarios(ss).enable_surface_reflection;
    channel = CARPE3D_vertical(paramsV);

    h = channel.h_total;
    if abs(h) < 1e-12
        error('Channel gain is too small for coherent equalization. scenario=%s', scenarios(ss).name);
    end

    ber = zeros(numel(EbN0_dB_list), 1);
    ser = zeros(numel(EbN0_dB_list), 1);
    effective_snr_db = inf(numel(EbN0_dB_list), 1);

    if ~strcmpi(isi_mode, 'linear_conv')
        error('Unsupported isi_mode: %s', isi_mode);
    end

    [f_bb_axis, H_baseband_shifted] = local_build_baseband_response( ...
        channel.f_axis, channel.H_f, channel.idx_f_ref, symbol_rate_hz, n_sym);
    [h_bb, n_tap_eff, energy_kept] = local_build_channel_taps(H_baseband_shifted, 0.999);
    rx_clean = conv(tx_symbols, h_bb, 'same');

    for ii = 1:numel(EbN0_dB_list)
        noise_cfg = struct();
        noise_cfg.enable_noise = noise_control.enable_noise;
        noise_cfg.model = noise_control.model;
        noise_cfg.ebn0_db = EbN0_dB_list(ii);
        noise_cfg.bits_per_symbol = k;
        noise_cfg.seed = noise_control.seed_base + 1000*ss + ii;
        noise_cfg.custom_noise_fn = [];

        [rx_noisy, ~, noise_meta] = noise_inject_vertical(rx_clean, noise_cfg, tx_symbols);

        % Frequency-domain MMSE equalization with known effective channel taps.
        rx_eq = local_mmse_equalize(rx_noisy, h_bb, 1e-6);
        bits_rx = modem_psk('demodulate', rx_eq, M);
        [ber(ii), ser(ii)] = modem_psk('error_rate', bits_tx, bits_rx, M);
        effective_snr_db(ii) = noise_meta.effective_snr_db;
    end

    results(ss).name = scenarios(ss).name;
    results(ss).enable_surface_reflection = scenarios(ss).enable_surface_reflection;
    results(ss).noise_enabled = noise_control.enable_noise;
    results(ss).noise_model = noise_control.model;
    results(ss).h_direct = channel.h_direct;
    results(ss).h_reflect = channel.h_reflect;
    results(ss).h_total = channel.h_total;
    results(ss).fd_hz_used = channel.fd_hz_used;
    results(ss).rx_state_used = channel.rx_state_used;
    results(ss).EbN0_dB_list = EbN0_dB_list;
    results(ss).BER = ber;
    results(ss).SER = ser;
    results(ss).effective_snr_db = effective_snr_db;
    results(ss).symbol_rate_hz = symbol_rate_hz;
    results(ss).isi_mode = isi_mode;
    results(ss).f_axis = channel.f_axis;
    results(ss).H_f = channel.H_f;
    results(ss).f_bb_axis = f_bb_axis;
    results(ss).H_baseband = H_baseband_shifted;
    results(ss).h_bb = h_bb;
    results(ss).h_bb_tap_count = n_tap_eff;
    results(ss).h_bb_energy_kept = energy_kept;
    results(ss).channel = channel;

    disp(['--- Scenario: ', scenarios(ss).name, ' ---'])
    disp(['h_direct=', num2str(channel.h_direct), ', h_reflect=', num2str(channel.h_reflect), ', h_total=', num2str(channel.h_total)])
    disp(['|h_total|=', num2str(abs(channel.h_total)), ', phase(rad)=', num2str(angle(channel.h_total)), ...
          ', fd_hz_used=', num2str(channel.fd_hz_used)])
    disp(['Nf=', num2str(numel(channel.f_axis)), ', f_ref=', num2str(channel.f_axis(channel.idx_f_ref)), ...
          ' Hz, effective taps=', num2str(n_tap_eff), ', tap_energy=', num2str(energy_kept)])
    T = table(EbN0_dB_list(:), ber, ser, effective_snr_db, ...
              'VariableNames', {'EbN0_dB', 'BER', 'SER', 'EffectiveSNR_dB'});
    disp(T)
end

figure(31); clf
semilogy(EbN0_dB_list, max(results(1).BER, 1e-6), 'o-', 'LineWidth', 1.2)
hold on
semilogy(EbN0_dB_list, max(results(2).BER, 1e-6), 's-', 'LineWidth', 1.2)
grid on
xlabel('Eb/N0 (dB)')
ylabel('BER')
title('QPSK BER: direct-only vs direct+reflect')
legend(results(1).name, results(2).name, 'Location', 'southwest')

figure(32); clf
semilogy(EbN0_dB_list, max(results(1).SER, 1e-6), 'o-', 'LineWidth', 1.2)
hold on
semilogy(EbN0_dB_list, max(results(2).SER, 1e-6), 's-', 'LineWidth', 1.2)
grid on
xlabel('Eb/N0 (dB)')
ylabel('SER')
title('QPSK SER: direct-only vs direct+reflect')
legend(results(1).name, results(2).name, 'Location', 'southwest')

save('psk_comm_result.mat', ...
     'paramsV', 'scenarios', 'results', 'noise_control', 'EbN0_dB_list', 'M', 'n_sym')

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
