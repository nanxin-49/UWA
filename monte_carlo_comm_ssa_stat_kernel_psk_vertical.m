% Reduced QPSK BER/SER Monte Carlo comparison for ssa_stat_kernel.
% This script keeps the existing CARPE3D_vertical -> H_f -> baseband taps
% -> peak_sync/MMSE link policy and compares only the surface model.

clear
format compact

result_file = getenv('SSA_STAT_COMM_RESULT_FILE');
if isempty(result_file)
    result_file = 'monte_carlo_comm_ssa_stat_kernel_psk_vertical_result.mat';
end
if exist(result_file, 'file')
    warning('Result file %s already exists and will be overwritten.', result_file);
end

params_base = local_base_params();
sea_hs_values = [0, 0.05, 0.2, 0.5];
sea_wind_speed = 5.0;
model_specs = local_model_specs();
model_names = {model_specs.model_name};

comm_cfg = struct();
comm_cfg.M = 4;
comm_cfg.n_sym = 1000;
comm_cfg.EbN0_dB_list = [0, 10, 20];
comm_cfg.symbol_rate_hz = 1000;
comm_cfg.tap_energy_ratio = 0.999;
comm_cfg.bits_seed = 9000;
comm_cfg.noise_seed_base = 7000;
comm_cfg.noise_model = 'awgn';
comm_cfg.mmse_reg_eps = 1e-6;
comm_cfg.receive_window_mode = 'peak_sync';
comm_cfg.ebn0_reference = 'rx_clean';

mc_count = 4;
mc_override = str2double(getenv('SSA_STAT_COMM_MC_COUNT'));
if isfinite(mc_override) && mc_override >= 1
    mc_count = round(mc_override);
end
seed_list = 12345 + (0:(mc_count - 1));

max_condition_pairs = Inf;
max_condition_override = str2double(getenv('SSA_STAT_COMM_MAX_CONDITIONS'));
if isfinite(max_condition_override) && max_condition_override >= 1
    max_condition_pairs = round(max_condition_override);
end

rng(comm_cfg.bits_seed, 'twister');
bits_per_symbol = round(log2(comm_cfg.M));
bits_tx = randi([0, 1], comm_cfg.n_sym * bits_per_symbol, 1);
[tx_symbols, bits_used] = modem_psk('modulate', bits_tx, comm_cfg.M);
bits_tx = bits_used;

n_hs = numel(sea_hs_values);
n_model = numel(model_names);
n_eb = numel(comm_cfg.EbN0_dB_list);
ber_samples = NaN(n_eb, mc_count, n_model, n_hs);
ser_samples = NaN(n_eb, mc_count, n_model, n_hs);
effective_snr_samples = NaN(n_eb, mc_count, n_model, n_hs);
run_results = struct([]);
run_rows = struct([]);

for ih = 1:n_hs
    for im = 1:n_model
        condition_pair_index = (ih - 1) * n_model + im;
        if condition_pair_index > max_condition_pairs
            continue
        end
        for mm = 1:mc_count
            paramsV = params_base;
            paramsV.sea_hs_target = sea_hs_values(ih);
            paramsV.sea_wind_speed = sea_wind_speed;
            paramsV.surface_boundary_model = model_specs(im).surface_boundary_model;
            paramsV.surface_ssa_kernel_mode = model_specs(im).surface_ssa_kernel_mode;
            paramsV.surface_ssa_geometry_source_id = model_specs(im).surface_ssa_geometry_source_id;
            paramsV.enable_surface_reflection = true;
            paramsV.sea_seed = seed_list(mm);

            fprintf(['SSA stat comm MC Hs %d/%d (%.3g m), model %d/%d (%s), ' ...
                'seed %d/%d, sea_seed=%d\n'], ...
                ih, n_hs, sea_hs_values(ih), im, n_model, model_names{im}, ...
                mm, mc_count, seed_list(mm));
            channel = CARPE3D_vertical(paramsV);

            invariant_error = max(abs(channel.H_f(:) - channel.H_direct_f(:) - channel.H_reflect_f(:)));
            if invariant_error > 1e-10
                error('monte_carlo_comm_ssa_stat_kernel_psk_vertical:Invariant', ...
                    'H_f invariant failed for %s Hs=%g seed %d: %.3e', ...
                    model_names{im}, sea_hs_values(ih), seed_list(mm), invariant_error);
            end

            [f_bb_axis, H_baseband_shifted] = local_build_baseband_response( ...
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
                noise_cfg.seed = comm_cfg.noise_seed_base + ih*100000 + mm*100 + ee;
                noise_cfg.custom_noise_fn = [];

                [rx_noisy, ~, noise_meta] = noise_inject_vertical(rx_clean, noise_cfg, noise_signal_ref);
                rx_eq = local_mmse_equalize(rx_noisy, h_eq, comm_cfg.mmse_reg_eps);
                bits_rx = modem_psk('demodulate', rx_eq, comm_cfg.M);
                [ber(ee), ser(ee)] = modem_psk('error_rate', bits_tx, bits_rx, comm_cfg.M);
                effective_snr_db(ee) = noise_meta.effective_snr_db;
            end

            ber_samples(:, mm, im, ih) = ber;
            ser_samples(:, mm, im, ih) = ser;
            effective_snr_samples(:, mm, im, ih) = effective_snr_db;

            ssa = channel.roughness_meta.ssa_stat_kernel_meta;
            run_result = struct();
            run_result.h_index = ih;
            run_result.model_index = im;
            run_result.model_name = string(model_names{im});
            run_result.seed_index = mm;
            run_result.sea_seed = seed_list(mm);
            run_result.sea_hs_target = sea_hs_values(ih);
            run_result.sea_wind_speed = sea_wind_speed;
            run_result.h_direct = channel.h_direct;
            run_result.h_reflect = channel.h_reflect;
            run_result.h_total = channel.h_total;
            run_result.idx_f_ref = channel.idx_f_ref;
            run_result.f_ref_hz = channel.f_axis(channel.idx_f_ref);
            run_result.invariant_error = invariant_error;
            run_result.max_abs_H_reflect_f = max(abs(channel.H_reflect_f(:)));
            run_result.f_bb_axis = f_bb_axis;
            run_result.h_bb_tap_count = n_tap_eff;
            run_result.h_bb_energy_kept = energy_kept;
            run_result.h_eq_tap_count = numel(h_eq);
            run_result.receive_meta = receive_meta;
            run_result.tap_metrics = tap_metrics;
            run_result.ssa_meta_compact = local_compact_ssa_meta(ssa);
            run_result.BER = ber;
            run_result.SER = ser;
            run_result.effective_snr_db = effective_snr_db;
            run_results = local_append_struct(run_results, run_result);

            row = struct();
            row.h_index = ih;
            row.model_index = im;
            row.model_name = string(model_names{im});
            row.seed_index = mm;
            row.sea_seed = seed_list(mm);
            row.sea_hs_target = sea_hs_values(ih);
            row.sea_wind_speed = sea_wind_speed;
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
            row.ssa_kernel_mode = string(ssa.kernel_mode);
            row.ssa_E_sca_over_E_inc = ssa.E_sca / max(ssa.E_inc, eps);
            row.ssa_energy_scale_applied = ssa.energy_scale_applied;
            row.ssa_energy_conservation_error = local_get_field_or_nan(ssa, 'energy_conservation_error');
            row.ssa_propagating_bin_fraction = local_get_field_or_nan(ssa, 'propagating_bin_fraction');
            row.BER = ber.';
            row.SER = ser.';
            row.effective_snr_db = effective_snr_db.';
            run_rows = local_append_struct(run_rows, row);
        end
    end
end

run_summary_table = struct2table(run_rows);
curve_stats = local_build_curve_stats( ...
    ber_samples, ser_samples, effective_snr_samples, ...
    comm_cfg.EbN0_dB_list, model_names, sea_hs_values);
curve_summary_table = local_build_curve_summary_table(curve_stats, model_names, sea_hs_values);

disp(curve_summary_table(:, {'model_name', 'sea_hs_target', 'EbN0_dB', ...
    'BER_mean', 'BER_std', 'SER_mean', 'SER_std'}))

save(result_file, 'params_base', 'comm_cfg', 'sea_hs_values', ...
    'sea_wind_speed', 'model_names', 'seed_list', 'bits_tx', ...
    'max_condition_pairs', ...
    'ber_samples', 'ser_samples', 'effective_snr_samples', ...
    'run_results', 'run_summary_table', 'curve_summary_table', 'curve_stats');
fprintf('Saved %s\n', result_file);

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
paramsV.surface_reflect_coeff = -1;
paramsV.surface_boundary_model = 'kirchhoff_spatial';
paramsV.surface_boundary_coupling_diagnostics = false;
paramsV.surface_boundary_redistribution_diagnostics = false;
paramsV.surface_ssa_random_scatter = true;
paramsV.surface_ssa_scatter_scale = 1.0;
paramsV.surface_ssa_seed_offset = 100000;
paramsV.surface_ssa_kernel_mode = 'pm_convolution';
paramsV.surface_ssa_geometry_source_id = '';
paramsV.surface_ssa_kz_branch = 'downward_positive_real';
paramsV.surface_ssa_conv_padding = 'periodic';
end

function model_specs = local_model_specs()
source_id = ['SSA.md; Thorsos & Broschat 1995 JASA, ', ...
    'Dirichlet SSA first-order / perturbation-limit geometry'];
model_specs = struct([]);
model_specs(1).model_name = 'kirchhoff_spatial';
model_specs(1).surface_boundary_model = 'kirchhoff_spatial';
model_specs(1).surface_ssa_kernel_mode = 'pm_convolution';
model_specs(1).surface_ssa_geometry_source_id = '';
model_specs(2).model_name = 'ssa_pm_convolution';
model_specs(2).surface_boundary_model = 'ssa_stat_kernel';
model_specs(2).surface_ssa_kernel_mode = 'pm_convolution';
model_specs(2).surface_ssa_geometry_source_id = '';
model_specs(3).model_name = 'ssa1_geometry';
model_specs(3).surface_boundary_model = 'ssa_stat_kernel';
model_specs(3).surface_ssa_kernel_mode = 'ssa1_geometry';
model_specs(3).surface_ssa_geometry_source_id = source_id;
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
    case 'peak_sync'
        h_eq = h_bb(peak_index:end);
        rx_full = conv(tx_symbols, h_eq, 'full');
        rx_clean = rx_full(1:numel(tx_symbols));
        start_index = peak_index;
        discarded_energy_fraction = sum(tap_energy(1:max(peak_index - 1, 0))) / max(sum(tap_energy), eps);
    case 'causal_head'
        rx_full = conv(tx_symbols, h_bb, 'full');
        rx_clean = rx_full(1:numel(tx_symbols));
        h_eq = h_bb;
        start_index = 1;
        discarded_energy_fraction = 0;
    case 'same_legacy'
        rx_clean = conv(tx_symbols, h_bb, 'same');
        h_eq = h_bb;
        start_index = NaN;
        discarded_energy_fraction = 0;
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

function out = local_compact_ssa_meta(meta)
out = struct();
out.enabled = meta.enabled;
out.random_scatter_enabled = meta.random_scatter_enabled;
out.kernel_mode = meta.kernel_mode;
out.geometry_source_id = meta.geometry_source_id;
out.kz_branch = meta.kz_branch;
out.conv_padding = meta.conv_padding;
out.sigma_eta_m = meta.sigma_eta_m;
out.Hs_target_m = meta.Hs_target_m;
out.R_coh = meta.R_coh;
out.P_sca = meta.P_sca;
out.E_inc = meta.E_inc;
out.E_coh = meta.E_coh;
out.E_sca_raw = meta.E_sca_raw;
out.E_sca_limited = local_get_field_or_nan(meta, 'E_sca_limited');
out.E_sca = meta.E_sca;
out.E_ref = local_get_field_or_nan(meta, 'E_ref');
out.energy_scale_applied = meta.energy_scale_applied;
out.energy_conservation_error = local_get_field_or_nan(meta, 'energy_conservation_error');
out.propagating_bin_fraction = local_get_field_or_nan(meta, 'propagating_bin_fraction');
out.seed_ssa = meta.seed_ssa;
end

function curve_stats = local_build_curve_stats( ...
    ber_samples, ser_samples, effective_snr_samples, EbN0_dB_list, model_names, sea_hs_values)

curve_stats = struct();
curve_stats.EbN0_dB_list = EbN0_dB_list(:);
curve_stats.model_names = string(model_names(:));
curve_stats.sea_hs_values = sea_hs_values(:);
curve_stats.BER = local_sample_stats(ber_samples);
curve_stats.SER = local_sample_stats(ser_samples);
curve_stats.effective_snr_db = local_sample_stats(effective_snr_samples);
curve_stats.validation = struct();
curve_stats.validation.BER_nonmonotonic_count = local_nonmonotonic_counts(curve_stats.BER.mean);
curve_stats.validation.SER_nonmonotonic_count = local_nonmonotonic_counts(curve_stats.SER.mean);
end

function stats = local_sample_stats(samples)
% samples has dimensions [EbN0 x mc_count x model x Hs].
stats = struct();
stats.samples = samples;
stats.mean = squeeze(mean(samples, 2, 'omitnan'));
stats.std = squeeze(std(samples, 0, 2, 'omitnan'));
stats.quantiles_5_25_50_75_95 = local_quantiles_dim(samples, [5, 25, 50, 75, 95]);
end

function q = local_quantiles_dim(samples, pct)
[n_eb, ~, n_model, n_hs] = size(samples);
q = NaN(n_eb, numel(pct), n_model, n_hs);
for ee = 1:n_eb
    for im = 1:n_model
        for ih = 1:n_hs
            x = squeeze(samples(ee, :, im, ih));
            q(ee, :, im, ih) = local_quantiles(x, pct);
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
    return
end
if numel(x) == 1
    q = repmat(x, size(pct));
    return
end
pos = 1 + (pct / 100) * (numel(x) - 1);
q = interp1(1:numel(x), x, pos, 'linear');
end

function counts = local_nonmonotonic_counts(mean_curves)
% mean_curves has dimensions [EbN0 x model x Hs].
[~, n_model, n_hs] = size(mean_curves);
counts = zeros(n_model, n_hs);
for im = 1:n_model
    for ih = 1:n_hs
        y = mean_curves(:, im, ih);
        counts(im, ih) = sum(diff(y) > 1e-12);
    end
end
end

function curve_summary_table = local_build_curve_summary_table(curve_stats, model_names, sea_hs_values)
rows = struct([]);
pct_idx = struct('q25', 2, 'q50', 3, 'q75', 4);
for ih = 1:numel(sea_hs_values)
    for im = 1:numel(model_names)
        for ee = 1:numel(curve_stats.EbN0_dB_list)
            row = struct();
            row.h_index = ih;
            row.model_index = im;
            row.model_name = string(model_names{im});
            row.sea_hs_target = sea_hs_values(ih);
            row.EbN0_dB = curve_stats.EbN0_dB_list(ee);
            row.BER_mean = curve_stats.BER.mean(ee, im, ih);
            row.BER_std = curve_stats.BER.std(ee, im, ih);
            row.BER_q25 = curve_stats.BER.quantiles_5_25_50_75_95(ee, pct_idx.q25, im, ih);
            row.BER_q50 = curve_stats.BER.quantiles_5_25_50_75_95(ee, pct_idx.q50, im, ih);
            row.BER_q75 = curve_stats.BER.quantiles_5_25_50_75_95(ee, pct_idx.q75, im, ih);
            row.SER_mean = curve_stats.SER.mean(ee, im, ih);
            row.SER_std = curve_stats.SER.std(ee, im, ih);
            row.SER_q25 = curve_stats.SER.quantiles_5_25_50_75_95(ee, pct_idx.q25, im, ih);
            row.SER_q50 = curve_stats.SER.quantiles_5_25_50_75_95(ee, pct_idx.q50, im, ih);
            row.SER_q75 = curve_stats.SER.quantiles_5_25_50_75_95(ee, pct_idx.q75, im, ih);
            row.effective_snr_db_mean = curve_stats.effective_snr_db.mean(ee, im, ih);
            row.effective_snr_db_std = curve_stats.effective_snr_db.std(ee, im, ih);
            rows = local_append_struct(rows, row);
        end
    end
end
curve_summary_table = struct2table(rows);
end

function out = local_get_field_or_nan(s, field_name)
if isfield(s, field_name)
    out = s.(field_name);
else
    out = NaN;
end
end

function out = local_append_struct(out, row)
if isempty(out)
    out = row;
else
    out(end + 1) = row; %#ok<AGROW>
end
end
