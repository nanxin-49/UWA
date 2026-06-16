% Monte Carlo statistics for random Kirchhoff rough-surface realizations.
% This script varies only sea_seed while keeping the sea-state and solver
% parameters fixed. Outputs are empirical diagnostics for the implemented
% phase-screen channel, not a new scattering model.

clear
format compact

result_file = 'monte_carlo_surface_channel_vertical_result.mat';
figure_prefix = 'monte_carlo_surface_channel_vertical_';
if exist(result_file, 'file')
    warning('Result file %s already exists and will be overwritten.', result_file);
end

mc_count = 16;
mc_override = str2double(getenv('SURFACE_MC_COUNT'));
if isfinite(mc_override) && mc_override >= 1
    mc_count = round(mc_override);
end
seed_list = 12345 + (0:(mc_count - 1));
symbol_rate_hz = 1000;
tap_fft_len = 512;
tap_energy_ratio = 0.999;

params_base = local_base_params();

H_samples = [];
H_direct_samples = [];
H_reflect_samples = [];
h_total_samples = complex(zeros(mc_count, 1));
h_direct_samples = complex(zeros(mc_count, 1));
h_reflect_samples = complex(zeros(mc_count, 1));
tap_metric_rows = struct([]);
channels_summary = struct([]);
f_axis_ref = [];
idx_f_ref_ref = NaN;

for ii = 1:mc_count
    paramsV = params_base;
    paramsV.sea_seed = seed_list(ii);
    fprintf('Monte Carlo surface case %d/%d, sea_seed=%d\n', ii, mc_count, seed_list(ii));
    channel = vertical_channel_model(paramsV);

    invariant_error = max(abs(channel.H_f(:) - channel.H_direct_f(:) - channel.H_reflect_f(:)));
    if invariant_error > 1e-10
        error('monte_carlo_surface_channel_vertical:Invariant', ...
            'H_f invariant failed for seed %d: %.3e', seed_list(ii), invariant_error);
    end

    if ii == 1
        f_axis_ref = channel.f_axis(:);
        idx_f_ref_ref = channel.idx_f_ref;
        n_freq = numel(f_axis_ref);
        H_samples = complex(zeros(n_freq, mc_count));
        H_direct_samples = complex(zeros(n_freq, mc_count));
        H_reflect_samples = complex(zeros(n_freq, mc_count));
    else
        if numel(channel.f_axis) ~= numel(f_axis_ref) || any(abs(channel.f_axis(:) - f_axis_ref) > 1e-9)
            error('Frequency axis changed at seed %d.', seed_list(ii));
        end
        if channel.idx_f_ref ~= idx_f_ref_ref
            error('idx_f_ref changed at seed %d.', seed_list(ii));
        end
    end

    H_samples(:, ii) = channel.H_f(:);
    H_direct_samples(:, ii) = channel.H_direct_f(:);
    H_reflect_samples(:, ii) = channel.H_reflect_f(:);
    h_total_samples(ii) = channel.h_total;
    h_direct_samples(ii) = channel.h_direct;
    h_reflect_samples(ii) = channel.h_reflect;

    [~, H_baseband] = local_build_baseband_response( ...
        channel.f_axis, channel.H_f, channel.idx_f_ref, symbol_rate_hz, tap_fft_len);
    [h_bb, tap_count, tap_energy_kept] = local_build_channel_taps(H_baseband, tap_energy_ratio);
    tap_metrics = local_tap_metrics(h_bb, tap_count, tap_energy_kept);
    tap_metric_rows = local_append_struct(tap_metric_rows, tap_metrics);

    coupling = channel.roughness_meta.boundary_coupling_diagnostics;
    redistribution = channel.roughness_meta.boundary_redistribution_diagnostics;
    if ~coupling.enabled || ~redistribution.enabled
        error('Boundary diagnostics were not enabled for seed %d.', seed_list(ii));
    end

    row = struct();
    row.seed = seed_list(ii);
    row.invariant_error = invariant_error;
    row.abs_h_ref = abs(channel.h_total);
    row.phase_h_ref_rad = angle(channel.h_total);
    row.abs_h_direct = abs(channel.h_direct);
    row.abs_h_reflect = abs(channel.h_reflect);
    row.reflect_direct_abs = abs(local_safe_complex_ratio(channel.h_reflect, channel.h_direct));
    row.reflect_direct_phase_rad = angle(local_safe_complex_ratio(channel.h_reflect, channel.h_direct));
    row.coupling_nonzero_power_fraction = coupling.nonzero_power_fraction;
    row.coupling_rms_delta_k_rad_per_m = coupling.rms_delta_k_rad_per_m;
    row.redistribution_reflect_rms_delta_k_rad_per_m = redistribution.reflect_rms_delta_k_rad_per_m;
    row.redistribution_reflect_high_k_fraction = redistribution.reflect_high_k_fraction;
    row.redistribution_rms_delta_k_increase_rad_per_m = redistribution.rms_delta_k_increase_rad_per_m;
    row.tap_count = tap_metrics.tap_count;
    row.tap_energy_kept = tap_metrics.tap_energy_kept;
    row.tap_rms_delay_symbols = tap_metrics.tap_rms_delay_symbols;
    row.tap_peak_index = tap_metrics.tap_peak_index;
    row.tap_peak_fraction = tap_metrics.tap_peak_fraction;
    channels_summary = local_append_struct(channels_summary, row);
end

summary_table = struct2table(channels_summary);
mc_stats = local_build_mc_stats( ...
    f_axis_ref, idx_f_ref_ref, H_samples, H_direct_samples, H_reflect_samples, ...
    h_total_samples, h_direct_samples, h_reflect_samples, channels_summary, tap_metric_rows);

disp(summary_table(:, {'seed', 'abs_h_ref', 'phase_h_ref_rad', ...
    'reflect_direct_abs', 'coupling_nonzero_power_fraction', ...
    'redistribution_reflect_rms_delta_k_rad_per_m', ...
    'redistribution_reflect_high_k_fraction', 'tap_rms_delay_symbols', ...
    'invariant_error'}))

local_plot_mean_magnitude(mc_stats, figure_prefix);
local_plot_hist(abs(h_total_samples), '|H(f_ref)|', [figure_prefix 'abs_h_ref_hist.png']);
local_plot_hist(angle(h_total_samples), 'angle(H(f_ref)) rad', [figure_prefix 'phase_h_ref_hist.png']);
local_plot_hist([channels_summary.redistribution_reflect_rms_delta_k_rad_per_m].', ...
    'reflect rms delta k (rad/m)', [figure_prefix 'reflect_rms_delta_k_hist.png']);
local_plot_hist([channels_summary.redistribution_reflect_high_k_fraction].', ...
    'reflect high-k fraction', [figure_prefix 'reflect_high_k_fraction_hist.png']);
local_plot_hist([tap_metric_rows.tap_rms_delay_symbols].', ...
    'tap RMS delay (symbols)', [figure_prefix 'tap_rms_delay_hist.png']);

save(result_file, 'params_base', 'seed_list', 'symbol_rate_hz', 'tap_fft_len', ...
    'tap_energy_ratio', 'channels_summary', 'summary_table', 'mc_stats');

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
paramsV.surface_boundary_coupling_diagnostics = true;
paramsV.surface_boundary_redistribution_diagnostics = true;
end

function out = local_append_struct(out, row)
if isempty(out)
    out = row;
else
    out(end + 1) = row; %#ok<AGROW>
end
end

function ratio = local_safe_complex_ratio(num, den)
if abs(den) <= eps
    den = complex(eps, 0);
end
ratio = num ./ den;
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

function mc_stats = local_build_mc_stats( ...
    f_axis, idx_f_ref, H_samples, H_direct_samples, H_reflect_samples, ...
    h_total_samples, h_direct_samples, h_reflect_samples, rows, tap_rows)

mean_f = mean(H_samples, 2);
var_f = mean(abs(H_samples - mean_f).^2, 2);
direct_mean_f = mean(H_direct_samples, 2);
direct_drift_max_abs = max(max(abs(H_direct_samples - direct_mean_f)));

abs_h_ref = abs(h_total_samples);
phase_h_ref = angle(h_total_samples);
ratio_ref = arrayfun(@local_safe_complex_ratio, h_reflect_samples, h_direct_samples);
ratio_abs = abs(ratio_ref);
ratio_phase = angle(ratio_ref);

mc_stats = struct();
mc_stats.H_f = struct( ...
    'f_axis', f_axis, ...
    'idx_f_ref', idx_f_ref, ...
    'samples', H_samples, ...
    'mean_f', mean_f, ...
    'var_f', var_f, ...
    'mean_abs_f', mean(abs(H_samples), 2), ...
    'std_abs_f', std(abs(H_samples), 0, 2));
mc_stats.H_direct_f = struct( ...
    'samples', H_direct_samples, ...
    'mean_f', direct_mean_f, ...
    'direct_drift_max_abs', direct_drift_max_abs);
mc_stats.H_reflect_f = struct( ...
    'samples', H_reflect_samples, ...
    'mean_f', mean(H_reflect_samples, 2), ...
    'var_f', mean(abs(H_reflect_samples - mean(H_reflect_samples, 2)).^2, 2));
mc_stats.reference_frequency = struct( ...
    'h_total_samples', h_total_samples, ...
    'abs_h_ref_mean', mean(abs_h_ref), ...
    'abs_h_ref_std', std(abs_h_ref), ...
    'abs_h_ref_quantiles_5_25_50_75_95', local_quantiles(abs_h_ref, [5, 25, 50, 75, 95]), ...
    'phase_h_ref_samples_rad', phase_h_ref, ...
    'phase_h_ref_circular_mean_rad', local_circular_mean(phase_h_ref), ...
    'phase_h_ref_circular_variance', local_circular_variance(phase_h_ref), ...
    'h_total_mean', mean(h_total_samples), ...
    'h_total_var', mean(abs(h_total_samples - mean(h_total_samples)).^2));
mc_stats.reflect_direct_ratio = struct( ...
    'samples', ratio_ref, ...
    'abs_mean', mean(ratio_abs), ...
    'abs_std', std(ratio_abs), ...
    'abs_quantiles_5_25_50_75_95', local_quantiles(ratio_abs, [5, 25, 50, 75, 95]), ...
    'phase_circular_mean_rad', local_circular_mean(ratio_phase), ...
    'phase_circular_variance', local_circular_variance(ratio_phase));
mc_stats.diagnostics = struct( ...
    'coupling_nonzero_power_fraction', local_distribution_stats([rows.coupling_nonzero_power_fraction].'), ...
    'coupling_rms_delta_k_rad_per_m', local_distribution_stats([rows.coupling_rms_delta_k_rad_per_m].'), ...
    'redistribution_reflect_rms_delta_k_rad_per_m', local_distribution_stats([rows.redistribution_reflect_rms_delta_k_rad_per_m].'), ...
    'redistribution_reflect_high_k_fraction', local_distribution_stats([rows.redistribution_reflect_high_k_fraction].'), ...
    'redistribution_rms_delta_k_increase_rad_per_m', local_distribution_stats([rows.redistribution_rms_delta_k_increase_rad_per_m].'));
mc_stats.tap_metrics = struct( ...
    'tap_count', local_distribution_stats([tap_rows.tap_count].'), ...
    'tap_energy_kept', local_distribution_stats([tap_rows.tap_energy_kept].'), ...
    'tap_rms_delay_symbols', local_distribution_stats([tap_rows.tap_rms_delay_symbols].'), ...
    'tap_peak_index', local_distribution_stats([tap_rows.tap_peak_index].'), ...
    'tap_peak_fraction', local_distribution_stats([tap_rows.tap_peak_fraction].'));
mc_stats.validation = struct( ...
    'max_invariant_error', max([rows.invariant_error]), ...
    'direct_drift_max_abs', direct_drift_max_abs);
end

function stats = local_distribution_stats(x)
x = x(:);
x = x(isfinite(x));
if isempty(x)
    stats = struct('mean', NaN, 'var', NaN, 'std', NaN, ...
        'quantiles_5_25_50_75_95', NaN(1, 5));
    return
end
mu = mean(x);
stats = struct( ...
    'mean', mu, ...
    'var', mean((x - mu).^2), ...
    'std', std(x), ...
    'quantiles_5_25_50_75_95', local_quantiles(x, [5, 25, 50, 75, 95]));
end

function q = local_quantiles(x, pct)
x = sort(x(:));
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

function mu = local_circular_mean(theta)
z = mean(exp(1i * theta(:)));
mu = angle(z);
end

function v = local_circular_variance(theta)
z = mean(exp(1i * theta(:)));
v = 1 - abs(z);
end

function local_plot_mean_magnitude(mc_stats, figure_prefix)
f_axis = mc_stats.H_f.f_axis;
mean_abs = mc_stats.H_f.mean_abs_f;
std_abs = mc_stats.H_f.std_abs_f;
upper_abs = mean_abs + std_abs;
lower_abs = max(mean_abs - std_abs, eps);

figure('Visible', 'off');
plot(f_axis, 20*log10(max(mean_abs, eps)), 'k-', 'LineWidth', 1.4)
hold on
plot(f_axis, 20*log10(max(upper_abs, eps)), 'r--', 'LineWidth', 1.0)
plot(f_axis, 20*log10(lower_abs), 'b--', 'LineWidth', 1.0)
grid on
xlabel('Frequency (Hz)')
ylabel('|H_f| (dB)')
title('Monte Carlo mean |H_f| with +/- one std envelope')
legend('mean', 'mean+std', 'mean-std', 'Location', 'best')
print(gcf, '-dpng', '-r200', [figure_prefix 'H_f_mean_magnitude.png'])
close(gcf)
end

function local_plot_hist(x, x_label_text, file_name)
x = x(:);
x = x(isfinite(x));
figure('Visible', 'off');
if isempty(x)
    histogram(NaN)
else
    histogram(x)
end
grid on
xlabel(x_label_text, 'Interpreter', 'none')
ylabel('count')
title(x_label_text, 'Interpreter', 'none')
print(gcf, '-dpng', '-r200', file_name)
close(gcf)
end
