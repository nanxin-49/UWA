% Minimal closed-loop validation for the vertical QPSK communication chain.
% The first cases use synthetic channels and do not depend on PE propagation.
% The final cases diagnose whether PE-derived baseband taps are misaligned.

clear
format compact

result_file = 'validate_comm_link_minimal_vertical_result.mat';
figure_prefix = 'validate_comm_link_minimal_vertical_';
if exist(result_file, 'file')
    warning('Result file %s already exists and will be overwritten.', result_file);
end

comm_cfg = struct();
comm_cfg.M = 4;
comm_cfg.n_sym = 2000;
comm_cfg.EbN0_dB_list = 0:2:20;
comm_cfg.symbol_rate_hz = 1000;
comm_cfg.tap_energy_ratio = 0.999;
comm_cfg.bits_seed = 9000;
comm_cfg.noise_seed_base = 7000;
comm_cfg.noise_model = 'awgn';
comm_cfg.mmse_reg_eps = 1e-6;
comm_cfg.ebn0_reference = 'rx_clean';
comm_cfg.receive_window_mode = 'peak_sync';

rng(comm_cfg.bits_seed, 'twister');
bits_per_symbol = round(log2(comm_cfg.M));
bits_tx = randi([0, 1], comm_cfg.n_sym * bits_per_symbol, 1);
[tx_symbols, bits_used] = modem_psk('modulate', bits_tx, comm_cfg.M);
bits_tx = bits_used;

test_cases = struct([]);
test_cases = local_append_struct(test_cases, local_case_spec('unit_channel', 1, 'same'));
test_cases = local_append_struct(test_cases, ...
    local_case_spec('single_tap_complex_gain', 0.7 * exp(1i * 0.8), 'same'));
test_cases = local_append_struct(test_cases, ...
    local_case_spec('known_short_multipath_same', ...
    [1; 0.35 * exp(1i * 0.5); 0.15 * exp(-1i * 0.7)], 'same'));
test_cases = local_append_struct(test_cases, ...
    local_case_spec('known_short_multipath_peak_sync', ...
    [1; 0.35 * exp(1i * 0.5); 0.15 * exp(-1i * 0.7)], 'peak_sync'));

case_results = struct([]);
for ii = 1:numel(test_cases)
    result = local_run_comm_case(test_cases(ii), ii, comm_cfg, bits_tx, tx_symbols);
    case_results = local_append_struct(case_results, result);
end

pe_alignment = local_run_pe_alignment_diagnostic(comm_cfg, bits_tx, tx_symbols, numel(case_results));
case_results = local_append_struct(case_results, pe_alignment.current_result);
case_results = local_append_struct(case_results, pe_alignment.aligned_result);
case_results = local_append_struct(case_results, pe_alignment.peak_sync_result);

summary_table = local_build_summary_table(case_results);
disp(summary_table(:, {'case_name', 'rx_mode', 'tap_count', 'peak_index', 'peak_fraction', ...
    'noiseless_BER', 'noiseless_SER', 'BER_at_20dB', 'SER_at_20dB', ...
    'BER_nonmonotonic_count', 'SER_nonmonotonic_count'}))

local_plot_comm_curves(case_results, comm_cfg, {'unit_channel', 'single_tap_complex_gain', ...
    'known_short_multipath_same', 'known_short_multipath_peak_sync'}, ...
    'BER', [figure_prefix 'simple_channel_BER.png']);
local_plot_comm_curves(case_results, comm_cfg, {'unit_channel', 'single_tap_complex_gain', ...
    'known_short_multipath_same', 'known_short_multipath_peak_sync'}, ...
    'SER', [figure_prefix 'simple_channel_SER.png']);
local_plot_taps(pe_alignment.h_current, [figure_prefix 'pe_h_bb_current_taps.png'], ...
    'PE direct-only h\_bb current taps');
local_plot_taps(pe_alignment.h_aligned, [figure_prefix 'pe_h_bb_aligned_taps.png'], ...
    'PE direct-only h\_bb peak-aligned taps');

save(result_file, 'comm_cfg', 'test_cases', 'summary_table', 'case_results', 'pe_alignment');

function spec = local_case_spec(name, h_bb, rx_mode)
if nargin < 3 || isempty(rx_mode)
    rx_mode = 'same';
end
spec = struct();
spec.name = name;
spec.h_bb = h_bb(:);
spec.source = 'synthetic';
spec.rx_mode = rx_mode;
end

function result = local_run_comm_case(spec, case_index, comm_cfg, bits_tx, tx_symbols)
fprintf('Validating communication case %d: %s\n', case_index, spec.name);
[noiseless_ber, noiseless_ser] = local_run_once(spec.h_bb, spec.rx_mode, false, NaN, NaN, ...
    comm_cfg, bits_tx, tx_symbols);

n_eb = numel(comm_cfg.EbN0_dB_list);
ber = NaN(n_eb, 1);
ser = NaN(n_eb, 1);
effective_snr_db = NaN(n_eb, 1);
for ee = 1:n_eb
    noise_seed = comm_cfg.noise_seed_base + case_index * 1000 + ee;
    [ber(ee), ser(ee), effective_snr_db(ee)] = local_run_once(spec.h_bb, spec.rx_mode, true, ...
        comm_cfg.EbN0_dB_list(ee), noise_seed, comm_cfg, bits_tx, tx_symbols);
end

tap_metrics = local_tap_metrics(spec.h_bb);
result = struct();
result.case_name = spec.name;
result.source = spec.source;
result.rx_mode = spec.rx_mode;
result.h_bb = spec.h_bb;
result.tap_metrics = tap_metrics;
result.noiseless_BER = noiseless_ber;
result.noiseless_SER = noiseless_ser;
result.BER = ber;
result.SER = ser;
result.effective_snr_db = effective_snr_db;
result.BER_nonmonotonic_count = sum(diff(ber) > 1e-12);
result.SER_nonmonotonic_count = sum(diff(ser) > 1e-12);
result.BER_at_20dB = ber(end);
result.SER_at_20dB = ser(end);
end

function [ber, ser, effective_snr_db] = local_run_once(h_bb, rx_mode, enable_noise, ebn0_db, noise_seed, ...
    comm_cfg, bits_tx, tx_symbols)
[rx_clean, h_eq] = local_apply_channel(tx_symbols, h_bb(:), rx_mode);
if enable_noise
    noise_cfg = struct();
    noise_cfg.enable_noise = true;
    noise_cfg.model = comm_cfg.noise_model;
    noise_cfg.ebn0_db = ebn0_db;
    noise_cfg.bits_per_symbol = round(log2(comm_cfg.M));
    noise_cfg.seed = noise_seed;
    noise_cfg.custom_noise_fn = [];
    noise_signal_ref = local_select_noise_reference(comm_cfg.ebn0_reference, tx_symbols, rx_clean);
    [rx_noisy, ~, noise_meta] = noise_inject_vertical(rx_clean, noise_cfg, noise_signal_ref);
    effective_snr_db = noise_meta.effective_snr_db;
else
    rx_noisy = rx_clean;
    effective_snr_db = Inf;
end
rx_eq = local_mmse_equalize(rx_noisy, h_eq, comm_cfg.mmse_reg_eps);
bits_rx = modem_psk('demodulate', rx_eq, comm_cfg.M);
[ber, ser] = modem_psk('error_rate', bits_tx, bits_rx, comm_cfg.M);
end

function [rx_clean, h_eq] = local_apply_channel(tx_symbols, h_bb, rx_mode)
switch lower(rx_mode)
    case 'same'
        rx_clean = conv(tx_symbols, h_bb, 'same');
        h_eq = h_bb;
    case 'causal_head'
        rx_full = conv(tx_symbols, h_bb, 'full');
        rx_clean = rx_full(1:numel(tx_symbols));
        h_eq = h_bb;
    case 'peak_sync'
        tap_energy = abs(h_bb).^2;
        if all(tap_energy == 0)
            peak_index = 1;
        else
            [~, peak_index] = max(tap_energy);
        end
        h_eq = h_bb(peak_index:end);
        rx_full = conv(tx_symbols, h_eq, 'full');
        rx_clean = rx_full(1:numel(tx_symbols));
    otherwise
        error('Unsupported rx_mode: %s', rx_mode);
end
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

function pe_alignment = local_run_pe_alignment_diagnostic(comm_cfg, bits_tx, tx_symbols, case_offset)
fprintf('Running PE direct-only baseband tap alignment diagnostic\n');
paramsV = local_pe_params();
channel = CARPE3D_vertical(paramsV);
invariant_error = max(abs(channel.H_f(:) - channel.H_direct_f(:) - channel.H_reflect_f(:)));
if invariant_error > 1e-10
    error('validate_comm_link_minimal_vertical:PEInvariant', ...
        'PE direct-only H_f invariant failed: %.3e', invariant_error);
end
if max(abs(channel.H_reflect_f(:))) > 1e-12
    error('PE direct-only reflected response is nonzero.');
end

[f_bb_axis, H_baseband_shifted] = local_build_baseband_response( ...
    channel.f_axis, channel.H_f, channel.idx_f_ref, comm_cfg.symbol_rate_hz, comm_cfg.n_sym);
[h_current, current_tap_count, current_energy_kept, h_full] = local_build_channel_taps( ...
    H_baseband_shifted, comm_cfg.tap_energy_ratio);
full_metrics = local_tap_metrics(h_full);

[~, peak_index] = max(abs(h_full).^2);
h_full_aligned = circshift(h_full, 1 - peak_index);
[h_aligned, aligned_tap_count, aligned_energy_kept] = local_trim_channel_taps( ...
    h_full_aligned, comm_cfg.tap_energy_ratio);

current_spec = local_case_spec('pe_direct_current_taps', h_current, 'same');
current_spec.source = 'pe_direct_only';
aligned_spec = local_case_spec('pe_direct_peak_aligned_taps', h_aligned, 'same');
aligned_spec.source = 'pe_direct_only_diagnostic_aligned';
peak_sync_spec = local_case_spec('pe_direct_peak_sync', h_current, 'peak_sync');
peak_sync_spec.source = 'pe_direct_only_d2_policy';

current_result = local_run_comm_case(current_spec, case_offset + 1, comm_cfg, bits_tx, tx_symbols);
aligned_result = local_run_comm_case(aligned_spec, case_offset + 2, comm_cfg, bits_tx, tx_symbols);
peak_sync_result = local_run_comm_case(peak_sync_spec, case_offset + 3, comm_cfg, bits_tx, tx_symbols);

pe_alignment = struct();
pe_alignment.paramsV = paramsV;
pe_alignment.channel_summary = struct( ...
    'idx_f_ref', channel.idx_f_ref, ...
    'f_ref_hz', channel.f_axis(channel.idx_f_ref), ...
    'invariant_error', invariant_error, ...
    'max_abs_H_reflect_f', max(abs(channel.H_reflect_f(:))));
pe_alignment.f_bb_axis = f_bb_axis;
pe_alignment.H_baseband_shifted = H_baseband_shifted;
pe_alignment.h_full = h_full;
pe_alignment.h_current = h_current;
pe_alignment.h_aligned = h_aligned;
pe_alignment.full_tap_metrics = full_metrics;
pe_alignment.current_tap_count = current_tap_count;
pe_alignment.current_energy_kept = current_energy_kept;
pe_alignment.aligned_tap_count = aligned_tap_count;
pe_alignment.aligned_energy_kept = aligned_energy_kept;
pe_alignment.peak_index_before_alignment = peak_index;
pe_alignment.current_result = current_result;
pe_alignment.aligned_result = aligned_result;
pe_alignment.peak_sync_result = peak_sync_result;
end

function paramsV = local_pe_params()
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
paramsV.enable_surface_reflection = false;
paramsV.sea_wind_speed = 5.0;
paramsV.sea_hs_target = 0.5;
paramsV.sea_seed = 12345;
paramsV.surface_boundary_model = 'kirchhoff_kdomain';
end

function summary_table = local_build_summary_table(case_results)
rows = struct([]);
for ii = 1:numel(case_results)
    metrics = case_results(ii).tap_metrics;
    row = struct();
    row.case_name = string(case_results(ii).case_name);
    row.source = string(case_results(ii).source);
    row.rx_mode = string(case_results(ii).rx_mode);
    row.tap_count = metrics.tap_count;
    row.peak_index = metrics.peak_index;
    row.peak_fraction = metrics.peak_fraction;
    row.tap_energy_kept = metrics.tap_energy_kept;
    row.tap_rms_delay_symbols = metrics.tap_rms_delay_symbols;
    row.noiseless_BER = case_results(ii).noiseless_BER;
    row.noiseless_SER = case_results(ii).noiseless_SER;
    row.BER_at_20dB = case_results(ii).BER_at_20dB;
    row.SER_at_20dB = case_results(ii).SER_at_20dB;
    row.BER_nonmonotonic_count = case_results(ii).BER_nonmonotonic_count;
    row.SER_nonmonotonic_count = case_results(ii).SER_nonmonotonic_count;
    rows = local_append_struct(rows, row);
end
summary_table = struct2table(rows);
end

function out = local_append_struct(out, row)
if isempty(out)
    out = row;
else
    out(end + 1) = row; %#ok<AGROW>
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

function [h_eff, n_tap_eff, energy_kept, h_full] = local_build_channel_taps(H_baseband_shifted, energy_ratio)
H_baseband_shifted = H_baseband_shifted(:);
h_full = ifft(ifftshift(H_baseband_shifted));
[h_eff, n_tap_eff, energy_kept] = local_trim_channel_taps(h_full, energy_ratio);
end

function [h_eff, n_tap_eff, energy_kept] = local_trim_channel_taps(h_full, energy_ratio)
h_full = h_full(:);
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

function metrics = local_tap_metrics(h_bb)
h_bb = h_bb(:);
tap_energy = abs(h_bb).^2;
total_energy = sum(tap_energy);
if total_energy <= 0
    mean_delay = NaN;
    rms_delay = NaN;
    peak_index = NaN;
    peak_fraction = NaN;
    energy_kept = 0;
else
    delay_idx = (0:(numel(h_bb) - 1)).';
    mean_delay = sum(delay_idx .* tap_energy) / total_energy;
    rms_delay = sqrt(sum(((delay_idx - mean_delay).^2) .* tap_energy) / total_energy);
    [peak_energy, peak_index] = max(tap_energy);
    peak_fraction = peak_energy / total_energy;
    energy_kept = 1;
end

metrics = struct( ...
    'tap_count', numel(h_bb), ...
    'tap_energy_kept', energy_kept, ...
    'tap_mean_delay_symbols', mean_delay, ...
    'tap_rms_delay_symbols', rms_delay, ...
    'peak_index', peak_index, ...
    'peak_fraction', peak_fraction);
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

function local_plot_comm_curves(case_results, comm_cfg, case_names, metric_name, file_name)
figure('Visible', 'off');
hold on
for ii = 1:numel(case_names)
    idx = local_find_case(case_results, case_names{ii});
    if isempty(idx)
        continue
    end
    y = case_results(idx).(metric_name);
    semilogy(comm_cfg.EbN0_dB_list, max(y, 1e-6), 'o-', ...
        'LineWidth', 1.2, 'DisplayName', strrep(case_names{ii}, '_', '\_'))
end
grid on
xlabel('Eb/N0 (dB)')
ylabel(metric_name)
title([metric_name ' simple-channel validation'])
legend('Location', 'southwest')
print(gcf, '-dpng', '-r200', file_name)
close(gcf)
end

function idx = local_find_case(case_results, case_name)
idx = [];
for ii = 1:numel(case_results)
    if strcmpi(case_results(ii).case_name, case_name)
        idx = ii;
        return
    end
end
end

function local_plot_taps(h_bb, file_name, title_text)
figure('Visible', 'off');
stem(0:(numel(h_bb)-1), abs(h_bb), 'filled')
grid on
xlabel('tap index (symbols, zero-based)')
ylabel('|h\_bb|')
title(title_text, 'Interpreter', 'none')
print(gcf, '-dpng', '-r200', file_name)
close(gcf)
end
