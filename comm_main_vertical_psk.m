run(fullfile(fileparts(mfilename('fullpath')), 'scripts', 'bootstrap_project.m'));
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

% Optional smoke-test overrides. Defaults above remain the reference demo.
paramsV.nx = local_env_scalar('COMM_NX', paramsV.nx);
paramsV.ny = local_env_scalar('COMM_NY', paramsV.ny);
paramsV.Nf_min = local_env_scalar('COMM_NF_MIN', paramsV.Nf_min);
paramsV.Nf_max = local_env_scalar('COMM_NF_MAX', paramsV.Nf_max);
paramsV.enforce_1_over_R = local_env_bool('COMM_ENFORCE_1_OVER_R', paramsV.enforce_1_over_R);

% --- communication setup ---
M = 4;                % QPSK (MPSK can be changed here)
n_sym = 2000;
EbN0_dB_list = 0:2:20;
M = local_env_scalar('COMM_M', M);
n_sym = local_env_scalar('COMM_N_SYM', n_sym);
EbN0_dB_list = local_env_numeric_list('COMM_EBN0_DB_LIST', EbN0_dB_list);
k = log2(M);
symbol_rate_hz = 1000;
isi_mode = 'linear_conv';
receive_window_mode = 'peak_sync';

noise_control = struct();
noise_control.enable_noise = true;
noise_control.model = 'awgn';
noise_control.seed_base = 7000;
noise_control.ebn0_reference = 'rx_clean';

bits_tx = randi([0, 1], n_sym*k, 1);
[tx_symbols, bits_used] = modem_psk('modulate', bits_tx, M);
bits_tx = bits_used;

scenarios = struct([]);
scenarios(1).name = 'direct_only';
scenarios(1).enable_surface_reflection = false;
scenarios(2).name = 'direct_plus_reflect';
scenarios(2).enable_surface_reflection = true;

% Optional external-channel entry. The default empty path preserves the
% original PE scenarios exactly. The MAT file must contain external_channels
% (struct array) or external_channel (scalar), each with H_f+f_axis_hz or h_t.
external_channel_file = strtrim(getenv('COMM_EXTERNAL_CHANNEL_FILE'));
use_external_channels = ~isempty(external_channel_file);
if use_external_channels
    loaded_external = load(external_channel_file);
    if isfield(loaded_external,'external_channels')
        external_channels = loaded_external.external_channels;
    elseif isfield(loaded_external,'external_channel')
        external_channels = loaded_external.external_channel;
    else
        error('COMM_EXTERNAL_CHANNEL_FILE must contain external_channels or external_channel.');
    end
    scenarios = repmat(struct('name','','enable_surface_reflection',true),numel(external_channels),1);
    for jj=1:numel(external_channels)
        if isfield(external_channels(jj),'name'), scenarios(jj).name=external_channels(jj).name;
        else, scenarios(jj).name=sprintf('external_%d',jj); end
    end
end

results = struct([]);
for ss = 1:numel(scenarios)
    if use_external_channels
        [channel,external_h_t]=local_external_channel(external_channels(ss),paramsV.f_ref_hz);
    else
        paramsV.enable_surface_reflection = scenarios(ss).enable_surface_reflection;
        channel = vertical_channel_model(paramsV); external_h_t=[];
    end

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

    if isempty(external_h_t)
        [f_bb_axis, H_baseband_shifted] = local_build_baseband_response( ...
            channel.f_axis, channel.H_f, channel.idx_f_ref, symbol_rate_hz, n_sym);
        [h_bb, n_tap_eff, energy_kept] = local_build_channel_taps(H_baseband_shifted, 0.999);
    else
        h_bb=external_h_t(:); n_tap_eff=numel(h_bb); energy_kept=1;
        f_bb_axis=[]; H_baseband_shifted=[];
    end
    [rx_clean, h_eq, receive_meta] = local_apply_channel_window(tx_symbols, h_bb, receive_window_mode);
    noise_signal_ref = local_select_noise_reference(noise_control.ebn0_reference, tx_symbols, rx_clean);

    for ii = 1:numel(EbN0_dB_list)
        noise_cfg = struct();
        noise_cfg.enable_noise = noise_control.enable_noise;
        noise_cfg.model = noise_control.model;
        noise_cfg.ebn0_db = EbN0_dB_list(ii);
        noise_cfg.bits_per_symbol = k;
        noise_cfg.seed = noise_control.seed_base + 1000*ss + ii;
        noise_cfg.custom_noise_fn = [];

        [rx_noisy, ~, noise_meta] = noise_inject_vertical(rx_clean, noise_cfg, noise_signal_ref);

        % Frequency-domain MMSE equalization with known effective channel taps.
        rx_eq = local_mmse_equalize(rx_noisy, h_eq, 1e-6);
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
    results(ss).receive_window_mode = receive_window_mode;
    results(ss).receive_meta = receive_meta;
    results(ss).ebn0_reference = noise_control.ebn0_reference;
    results(ss).f_axis = channel.f_axis;
    results(ss).H_f = channel.H_f;
    results(ss).f_bb_axis = f_bb_axis;
    results(ss).H_baseband = H_baseband_shifted;
    results(ss).h_bb = h_bb;
    results(ss).h_eq = h_eq;
    results(ss).h_bb_tap_count = n_tap_eff;
    results(ss).h_bb_energy_kept = energy_kept;
    results(ss).channel = channel;

    disp(['--- Scenario: ', scenarios(ss).name, ' ---'])
    disp(['h_direct=', num2str(channel.h_direct), ', h_reflect=', num2str(channel.h_reflect), ', h_total=', num2str(channel.h_total)])
    disp(['|h_total|=', num2str(abs(channel.h_total)), ', phase(rad)=', num2str(angle(channel.h_total)), ...
          ', fd_hz_used=', num2str(channel.fd_hz_used)])
    disp(['Nf=', num2str(numel(channel.f_axis)), ', f_ref=', num2str(channel.f_axis(channel.idx_f_ref)), ...
          ' Hz, effective taps=', num2str(n_tap_eff), ', tap_energy=', num2str(energy_kept), ...
          ', rx_window=', receive_window_mode, ', peak_idx=', num2str(receive_meta.peak_index_original), ...
          ', ebn0_ref=', noise_control.ebn0_reference])
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

function [channel,h_t]=local_external_channel(ext,f_ref_hz)
h_t=[];
if isfield(ext,'H_f') && ~isempty(ext.H_f)
    if ~isfield(ext,'f_axis_hz'), error('External H_f requires f_axis_hz.'); end
    f=ext.f_axis_hz(:); H=ext.H_f(:); if numel(f)~=numel(H), error('External H_f/f_axis_hz mismatch.'); end
    [~,idx]=min(abs(f-f_ref_hz));
elseif isfield(ext,'h_t') && ~isempty(ext.h_t)
    h_t=ext.h_t(:); H=fftshift(fft(h_t)); f=(1:numel(H)).'; idx=ceil(numel(H)/2);
else
    error('Each external channel requires H_f+f_axis_hz or h_t.');
end
channel=struct('H_f',H,'H_direct_f',complex(zeros(size(H))), ...
    'H_reflect_f',H,'f_axis',f,'idx_f_ref',idx,'h_total',H(idx), ...
    'h_direct',0,'h_reflect',H(idx),'fd_hz_used',0,'rx_state_used',struct(), ...
    'external_input',true,'noise_included',false);
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

function value = local_env_scalar(name, default_value)
raw = getenv(name);
if isempty(raw)
    value = default_value;
    return
end
parsed = str2double(raw);
if ~isfinite(parsed)
    warning('Ignoring invalid numeric environment override %s=%s.', name, raw);
    value = default_value;
    return
end
value = parsed;
end

function tf = local_env_bool(name, default_value)
raw = getenv(name);
if isempty(raw)
    tf = default_value;
    return
end
switch lower(strtrim(raw))
    case {'1', 'true', 'yes', 'on'}
        tf = true;
    case {'0', 'false', 'no', 'off'}
        tf = false;
    otherwise
        warning('Ignoring invalid logical environment override %s=%s.', name, raw);
        tf = default_value;
end
end

function values = local_env_numeric_list(name, default_values)
raw = getenv(name);
if isempty(raw)
    values = default_values;
    return
end
parts = regexp(raw, '[,;\s]+', 'split');
parts = parts(~cellfun('isempty', parts));
values = str2double(parts);
if isempty(values) || any(~isfinite(values))
    warning('Ignoring invalid numeric-list environment override %s=%s.', name, raw);
    values = default_values;
end
end

