% End-to-end MPSK communication demo over vertical PE channel.
% Scenario: seabed instrument TX -> hydrophone at fixed z_rx=3 m below buoy.

clear
format compact

% --- channel / geometry setup ---
paramsV = struct();
paramsV.f0 = 4000;
paramsV.c0 = 1500;
paramsV.z_max = 100;
paramsV.stepz_lamb = 0.5;

paramsV.xw = 50;
paramsV.yw = 50;
paramsV.nx = 1024;
paramsV.ny = 1024;

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
paramsV.taper_ratio = 0.12;
paramsV.env_mode = 'uniform';
paramsV.enforce_1_over_R = true;
paramsV.show_figures = false;

paramsV.sea_wind_speed = 5.0;
paramsV.sea_hs_target = 0.5;
paramsV.sea_seed = 12345;

% --- communication setup ---
M = 4;                % QPSK (MPSK can be changed here)
n_sym = 2000;
EbN0_dB_list = 0:2:20;
k = log2(M);

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

    rx_clean = h * tx_symbols;
    for ii = 1:numel(EbN0_dB_list)
        noise_cfg = struct();
        noise_cfg.enable_noise = noise_control.enable_noise;
        noise_cfg.model = noise_control.model;
        noise_cfg.ebn0_db = EbN0_dB_list(ii);
        noise_cfg.bits_per_symbol = k;
        noise_cfg.seed = noise_control.seed_base + 1000*ss + ii;
        noise_cfg.custom_noise_fn = [];

        [rx_noisy, ~, noise_meta] = noise_inject_vertical(rx_clean, noise_cfg, tx_symbols);

        % Coherent equalization assumes perfect channel estimate.
        rx_eq = rx_noisy / h;
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
    results(ss).channel = channel;

    disp(['--- Scenario: ', scenarios(ss).name, ' ---'])
    disp(['h_direct=', num2str(channel.h_direct), ', h_reflect=', num2str(channel.h_reflect), ', h_total=', num2str(channel.h_total)])
    disp(['|h_total|=', num2str(abs(channel.h_total)), ', phase(rad)=', num2str(angle(channel.h_total)), ...
          ', fd_hz_used=', num2str(channel.fd_hz_used)])
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
