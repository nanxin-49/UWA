run(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'bootstrap_project.m'));
%VALIDATE_KSTAT_VS_KDOMAIN_CHANNEL_STATS_VERTICAL
% Channel-level Monte Carlo validation of Kirchhoff K-Stat against explicit
% Kirchhoff k-domain ensemble statistics after PE/WAPE propagation.
%
% This script compares receiver-side statistics only. It does not require
% single-seed complex fields, h_reflect, h_total, or PDPs to match pointwise.

clear
clc
format compact

cfg = struct();
cfg.Hs_list_m = local_env_numeric_vector('KSTAT_CHAN_HS_LIST', [0.05, 0.1, 0.2, 0.5]);
cfg.seed_count = local_env_int('KSTAT_CHAN_SEED_COUNT', 64);
cfg.grid_n = local_env_int('KSTAT_CHAN_GRID_N', 128);
cfg.nx = cfg.grid_n;
cfg.ny = cfg.grid_n;
cfg.xw_m = local_env_scalar('KSTAT_CHAN_XW_M', 50 * cfg.nx / 128);
cfg.yw_m = local_env_scalar('KSTAT_CHAN_YW_M', 50 * cfg.ny / 128);
cfg.Nf = local_env_int('KSTAT_CHAN_NF', 32);
cfg.f_band_hz = local_env_numeric_vector('KSTAT_CHAN_F_BAND', [4000, 8000]);
cfg.wind_mps = local_env_scalar('KSTAT_CHAN_WIND_MPS', 5.0);
cfg.seed_list = 12345 + (0:(cfg.seed_count - 1));
cfg.f_ref_hz = 6000;
cfg.symbol_rate_hz = 1000;
cfg.n_fft_pdp = 2048;
cfg.tap_energy_ratio = 0.999;
cfg.roughness_scale_mode = 'target_hs';
cfg.invariant_tol = 1e-10;
cfg.energy_tol = 1e-8;
cfg.E_href2_rel_tol = 0.25;
cfg.pdp_l2_tol = 0.25;
cfg.pdp_corr_tol = 0.90;
cfg.abs_h_total_quantile_nrmse_tol = 0.20;
cfg.figure_resolution = 160;

if numel(cfg.f_band_hz) ~= 2 || cfg.f_band_hz(2) <= cfg.f_band_hz(1)
    error('KSTAT_CHAN_F_BAND must contain [f_min f_max] with f_max > f_min.');
end
cfg.seed_count = max(1, cfg.seed_count);
cfg.Nf = max(2, cfg.Nf);

cfg.result_mat = project_result_file('validation', ...
    'validate_kstat_vs_kdomain_channel_stats_vertical_result.mat');
cfg.summary_csv = project_result_file('validation', ...
    'validate_kstat_vs_kdomain_channel_stats_vertical_summary.csv');
cfg.runs_csv = project_result_file('validation', ...
    'validate_kstat_vs_kdomain_channel_stats_vertical_runs.csv');
cfg.report_file = project_result_file('validation', ...
    'validate_kstat_vs_kdomain_channel_stats_vertical_report.md');
cfg.figure_dir = project_result_dir('validation');

fprintf('K-Stat vs K-Domain channel-stat validation\n');
fprintf('  Hs: %s m\n', mat2str(cfg.Hs_list_m));
fprintf('  seeds: %d, grid: %d x %d, aperture: %.3g x %.3g m, Nf: %d\n', ...
    cfg.seed_count, cfg.nx, cfg.ny, cfg.xw_m, cfg.yw_m, cfg.Nf);
fprintf('  f_band: %s Hz, f_ref: %.3f Hz, PM wind: %.3g m/s\n', ...
    mat2str(cfg.f_band_hz), cfg.f_ref_hz, cfg.wind_mps);

branch_names = {'kirchhoff_kdomain', 'kirchhoff_kstat'};
n_hs = numel(cfg.Hs_list_m);
n_branch = numel(branch_names);
n_seed = cfg.seed_count;
n_fft = cfg.n_fft_pdp;

h_reflect_samples = complex(NaN(n_seed, n_branch, n_hs), NaN(n_seed, n_branch, n_hs));
h_total_samples = complex(NaN(n_seed, n_branch, n_hs), NaN(n_seed, n_branch, n_hs));
pdp_reflect_raw = NaN(n_fft, n_seed, n_branch, n_hs);
pdp_reflect_sync = NaN(n_fft, n_seed, n_branch, n_hs);
pdp_total_raw = NaN(n_fft, n_seed, n_branch, n_hs);
pdp_total_sync = NaN(n_fft, n_seed, n_branch, n_hs);
tap_metrics_reflect = repmat(local_empty_tap_metrics(), n_seed, n_branch, n_hs);
tap_metrics_total = repmat(local_empty_tap_metrics(), n_seed, n_branch, n_hs);

run_rows = struct([]);
run_idx = 0;
summary_rows = struct([]);
summary_idx = 0;

for ih = 1:n_hs
    Hs_m = cfg.Hs_list_m(ih);
    direct_ref = [];
    fprintf('Hs %.4g m\n', Hs_m);

    for iseed = 1:n_seed
        seed = cfg.seed_list(iseed);
        fprintf('  seed %d/%d: %d\n', iseed, n_seed, seed);

        for ib = 1:n_branch
            branch = branch_names{ib};
            paramsV = local_base_params(cfg, Hs_m, seed, branch);
            channel = vertical_channel_model(paramsV);
            invariant_error = local_invariant_error(channel);
            if invariant_error > cfg.invariant_tol
                error('validate_kstat_vs_kdomain_channel_stats_vertical:Invariant', ...
                    'H_f invariant failed for %s Hs %.4g seed %d: %.3e', ...
                    branch, Hs_m, seed, invariant_error);
            end

            if isempty(direct_ref)
                direct_ref = channel.H_direct_f(:);
            end
            direct_delta = max(abs(channel.H_direct_f(:) - direct_ref));

            kstat_energy_error = 0;
            if strcmp(branch, 'kirchhoff_kstat')
                kstat_energy_error = channel.roughness_meta.kirchhoff_kstat_meta.phase_screen_energy_error;
                if kstat_energy_error > cfg.energy_tol
                    error('validate_kstat_vs_kdomain_channel_stats_vertical:KStatEnergy', ...
                        'K-Stat phase-screen energy failed Hs %.4g seed %d: %.3e', ...
                        Hs_m, seed, kstat_energy_error);
                end
            end

            [~, H_reflect_bb] = local_build_baseband_response( ...
                channel.f_axis, channel.H_reflect_f, channel.idx_f_ref, ...
                cfg.symbol_rate_hz, n_fft);
            [~, H_total_bb] = local_build_baseband_response( ...
                channel.f_axis, channel.H_f, channel.idx_f_ref, ...
                cfg.symbol_rate_hz, n_fft);

            h_reflect_full = ifft(ifftshift(H_reflect_bb));
            h_total_full = ifft(ifftshift(H_total_bb));
            [h_reflect_sync, h_total_sync, peak_index] = local_peak_sync_taps( ...
                h_reflect_full, h_total_full);

            pdp_ref_raw = abs(h_reflect_full).^2;
            pdp_tot_raw = abs(h_total_full).^2;
            pdp_ref_sync = abs(h_reflect_sync).^2;
            pdp_tot_sync = abs(h_total_sync).^2;

            h_reflect_samples(iseed, ib, ih) = channel.h_reflect;
            h_total_samples(iseed, ib, ih) = channel.h_total;
            pdp_reflect_raw(:, iseed, ib, ih) = pdp_ref_raw;
            pdp_reflect_sync(:, iseed, ib, ih) = pdp_ref_sync;
            pdp_total_raw(:, iseed, ib, ih) = pdp_tot_raw;
            pdp_total_sync(:, iseed, ib, ih) = pdp_tot_sync;
            tap_metrics_reflect(iseed, ib, ih) = local_tap_metrics(h_reflect_sync, cfg.tap_energy_ratio);
            tap_metrics_total(iseed, ib, ih) = local_tap_metrics(h_total_sync, cfg.tap_energy_ratio);

            run_idx = run_idx + 1;
            run_rows(run_idx).Hs_m = Hs_m; %#ok<SAGROW>
            run_rows(run_idx).seed = seed;
            run_rows(run_idx).branch_code = ib;
            run_rows(run_idx).abs_h_reflect = abs(channel.h_reflect);
            run_rows(run_idx).abs_h_reflect2 = abs(channel.h_reflect)^2;
            run_rows(run_idx).abs_h_total = abs(channel.h_total);
            run_rows(run_idx).abs_h_total2 = abs(channel.h_total)^2;
            run_rows(run_idx).pdp_reflect_energy_sync = sum(pdp_ref_sync);
            run_rows(run_idx).pdp_total_energy_sync = sum(pdp_tot_sync);
            run_rows(run_idx).pdp_reflect_rms_delay_samples = tap_metrics_reflect(iseed, ib, ih).rms_delay_samples;
            run_rows(run_idx).pdp_total_rms_delay_samples = tap_metrics_total(iseed, ib, ih).rms_delay_samples;
            run_rows(run_idx).pdp_reflect_peak_fraction = tap_metrics_reflect(iseed, ib, ih).peak_fraction;
            run_rows(run_idx).pdp_total_peak_fraction = tap_metrics_total(iseed, ib, ih).peak_fraction;
            run_rows(run_idx).total_peak_index_raw = peak_index;
            run_rows(run_idx).invariant_error = invariant_error;
            run_rows(run_idx).direct_path_max_abs_delta = direct_delta;
            run_rows(run_idx).kstat_phase_screen_energy_error = kstat_energy_error;
        end
    end

    kd = 1;
    ks = 2;
    kd_h_ref = h_reflect_samples(:, kd, ih);
    ks_h_ref = h_reflect_samples(:, ks, ih);
    kd_h_total = h_total_samples(:, kd, ih);
    ks_h_total = h_total_samples(:, ks, ih);
    kd_abs_total = abs(kd_h_total);
    ks_abs_total = abs(ks_h_total);

    mean_pdp_ref_kd = mean(pdp_reflect_sync(:, :, kd, ih), 2);
    mean_pdp_ref_ks = mean(pdp_reflect_sync(:, :, ks, ih), 2);
    mean_pdp_total_kd = mean(pdp_total_sync(:, :, kd, ih), 2);
    mean_pdp_total_ks = mean(pdp_total_sync(:, :, ks, ih), 2);

    [pdp_ref_l2, pdp_ref_corr] = local_compare_pdp(mean_pdp_ref_kd, mean_pdp_ref_ks);
    [pdp_total_l2, pdp_total_corr] = local_compare_pdp(mean_pdp_total_kd, mean_pdp_total_ks);
    q_kd = local_quantiles(kd_abs_total, [10, 50, 90]);
    q_ks = local_quantiles(ks_abs_total, [10, 50, 90]);
    quantile_nrmse = sqrt(mean((q_kd - q_ks).^2)) / max(mean([q_kd(:); q_ks(:)]), eps);

    summary_idx = summary_idx + 1;
    summary_rows(summary_idx).Hs_m = Hs_m; %#ok<SAGROW>
    summary_rows(summary_idx).seed_count = n_seed;
    summary_rows(summary_idx).nx = cfg.nx;
    summary_rows(summary_idx).ny = cfg.ny;
    summary_rows(summary_idx).Nf = cfg.Nf;
    summary_rows(summary_idx).E_abs_h_reflect2_kdomain = mean(abs(kd_h_ref).^2);
    summary_rows(summary_idx).E_abs_h_reflect2_kstat = mean(abs(ks_h_ref).^2);
    summary_rows(summary_idx).E_abs_h_reflect2_rel_error = local_rel_error( ...
        summary_rows(summary_idx).E_abs_h_reflect2_kdomain, ...
        summary_rows(summary_idx).E_abs_h_reflect2_kstat);
    summary_rows(summary_idx).abs_h_total_mean_kdomain = mean(kd_abs_total);
    summary_rows(summary_idx).abs_h_total_std_kdomain = std(kd_abs_total, 0);
    summary_rows(summary_idx).abs_h_total_p10_kdomain = q_kd(1);
    summary_rows(summary_idx).abs_h_total_p50_kdomain = q_kd(2);
    summary_rows(summary_idx).abs_h_total_p90_kdomain = q_kd(3);
    summary_rows(summary_idx).abs_h_total_mean_kstat = mean(ks_abs_total);
    summary_rows(summary_idx).abs_h_total_std_kstat = std(ks_abs_total, 0);
    summary_rows(summary_idx).abs_h_total_p10_kstat = q_ks(1);
    summary_rows(summary_idx).abs_h_total_p50_kstat = q_ks(2);
    summary_rows(summary_idx).abs_h_total_p90_kstat = q_ks(3);
    summary_rows(summary_idx).abs_h_total_quantile_nrmse = quantile_nrmse;
    summary_rows(summary_idx).abs_h_total_ecdf_max_delta = local_ecdf_max_delta(kd_abs_total, ks_abs_total);
    summary_rows(summary_idx).pdp_reflect_energy_kdomain = sum(mean_pdp_ref_kd);
    summary_rows(summary_idx).pdp_reflect_energy_kstat = sum(mean_pdp_ref_ks);
    summary_rows(summary_idx).pdp_reflect_energy_rel_error = local_rel_error( ...
        summary_rows(summary_idx).pdp_reflect_energy_kdomain, ...
        summary_rows(summary_idx).pdp_reflect_energy_kstat);
    summary_rows(summary_idx).pdp_reflect_l2_error = pdp_ref_l2;
    summary_rows(summary_idx).pdp_reflect_corr = pdp_ref_corr;
    summary_rows(summary_idx).pdp_total_energy_kdomain = sum(mean_pdp_total_kd);
    summary_rows(summary_idx).pdp_total_energy_kstat = sum(mean_pdp_total_ks);
    summary_rows(summary_idx).pdp_total_energy_rel_error = local_rel_error( ...
        summary_rows(summary_idx).pdp_total_energy_kdomain, ...
        summary_rows(summary_idx).pdp_total_energy_kstat);
    summary_rows(summary_idx).pdp_total_l2_error = pdp_total_l2;
    summary_rows(summary_idx).pdp_total_corr = pdp_total_corr;
    summary_rows(summary_idx).pdp_reflect_rms_delay_mean_kdomain = ...
        local_metric_mean(tap_metrics_reflect(:, kd, ih), 'rms_delay_samples');
    summary_rows(summary_idx).pdp_reflect_rms_delay_std_kdomain = ...
        local_metric_std(tap_metrics_reflect(:, kd, ih), 'rms_delay_samples');
    summary_rows(summary_idx).pdp_reflect_rms_delay_mean_kstat = ...
        local_metric_mean(tap_metrics_reflect(:, ks, ih), 'rms_delay_samples');
    summary_rows(summary_idx).pdp_reflect_rms_delay_std_kstat = ...
        local_metric_std(tap_metrics_reflect(:, ks, ih), 'rms_delay_samples');
    summary_rows(summary_idx).pdp_total_rms_delay_mean_kdomain = ...
        local_metric_mean(tap_metrics_total(:, kd, ih), 'rms_delay_samples');
    summary_rows(summary_idx).pdp_total_rms_delay_std_kdomain = ...
        local_metric_std(tap_metrics_total(:, kd, ih), 'rms_delay_samples');
    summary_rows(summary_idx).pdp_total_rms_delay_mean_kstat = ...
        local_metric_mean(tap_metrics_total(:, ks, ih), 'rms_delay_samples');
    summary_rows(summary_idx).pdp_total_rms_delay_std_kstat = ...
        local_metric_std(tap_metrics_total(:, ks, ih), 'rms_delay_samples');
    summary_rows(summary_idx).pdp_reflect_peak_fraction_mean_kdomain = ...
        local_metric_mean(tap_metrics_reflect(:, kd, ih), 'peak_fraction');
    summary_rows(summary_idx).pdp_reflect_peak_fraction_mean_kstat = ...
        local_metric_mean(tap_metrics_reflect(:, ks, ih), 'peak_fraction');
    summary_rows(summary_idx).pdp_total_peak_fraction_mean_kdomain = ...
        local_metric_mean(tap_metrics_total(:, kd, ih), 'peak_fraction');
    summary_rows(summary_idx).pdp_total_peak_fraction_mean_kstat = ...
        local_metric_mean(tap_metrics_total(:, ks, ih), 'peak_fraction');
    summary_rows(summary_idx).max_kdomain_invariant_error = ...
        max([run_rows([run_rows.Hs_m] == Hs_m & [run_rows.branch_code] == kd).invariant_error]);
    summary_rows(summary_idx).max_kstat_invariant_error = ...
        max([run_rows([run_rows.Hs_m] == Hs_m & [run_rows.branch_code] == ks).invariant_error]);
    summary_rows(summary_idx).max_kstat_energy_error = ...
        max([run_rows([run_rows.Hs_m] == Hs_m & [run_rows.branch_code] == ks).kstat_phase_screen_energy_error]);
    summary_rows(summary_idx).max_direct_path_abs_delta = ...
        max([run_rows([run_rows.Hs_m] == Hs_m).direct_path_max_abs_delta]);
end

run_table = struct2table(run_rows);
summary_table = struct2table(summary_rows);
validation_report = local_validation_report(summary_table, run_table, cfg);
figure_files = local_write_figures(summary_table, pdp_reflect_sync, pdp_total_sync, ...
    h_total_samples, cfg);
local_write_report(cfg.report_file, validation_report, summary_table, cfg, figure_files);

save(cfg.result_mat, 'cfg', 'run_table', 'summary_table', 'validation_report', ...
    'figure_files', 'h_reflect_samples', 'h_total_samples', ...
    'pdp_reflect_raw', 'pdp_reflect_sync', 'pdp_total_raw', 'pdp_total_sync', ...
    'tap_metrics_reflect', 'tap_metrics_total', '-v7.3');
writetable(summary_table, cfg.summary_csv);
writetable(run_table, cfg.runs_csv);

disp(summary_table(:, {'Hs_m', 'seed_count', 'E_abs_h_reflect2_kdomain', ...
    'E_abs_h_reflect2_kstat', 'E_abs_h_reflect2_rel_error', ...
    'pdp_reflect_l2_error', 'pdp_reflect_corr', ...
    'pdp_total_l2_error', 'pdp_total_corr', 'abs_h_total_quantile_nrmse', ...
    'max_kstat_energy_error'}));
disp(validation_report);
fprintf('Saved %s\n', cfg.result_mat);
fprintf('Saved %s\n', cfg.summary_csv);
fprintf('Saved %s\n', cfg.runs_csv);
fprintf('Saved %s\n', cfg.report_file);

if ~validation_report.all_hard_checks_passed
    error('validate_kstat_vs_kdomain_channel_stats_vertical:HardChecks', ...
        'One or more hard checks failed. See validation_report.');
end

function paramsV = local_base_params(cfg, Hs_m, seed, boundary_model)
paramsV = struct();
paramsV.f0 = cfg.f_band_hz(1);
paramsV.enable_wideband = true;
paramsV.f_band_hz = cfg.f_band_hz;
paramsV.Nf_min = cfg.Nf;
paramsV.Nf_max = cfg.Nf;
paramsV.f_ref_hz = cfg.f_ref_hz;
paramsV.c0 = 1500;
paramsV.z_max = 100;
paramsV.stepz_lamb = 0.5;
paramsV.xw = cfg.xw_m;
paramsV.yw = cfg.yw_m;
paramsV.nx = cfg.nx;
paramsV.ny = cfg.ny;
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
paramsV.sea_wind_speed = cfg.wind_mps;
paramsV.sea_hs_target = Hs_m;
paramsV.surface_roughness_scale_mode = cfg.roughness_scale_mode;
paramsV.sea_seed = seed;
paramsV.surface_reflect_coeff = -1;
paramsV.surface_phase_mode = 'normal';
paramsV.surface_boundary_model = boundary_model;
paramsV.surface_boundary_coupling_diagnostics = false;
paramsV.surface_boundary_redistribution_diagnostics = false;
paramsV.surface_kstat_random_scatter = true;
paramsV.surface_kstat_seed_offset = 200000;
paramsV.surface_kstat_conv_padding = 'periodic';
paramsV.surface_kstat_trusted_angle_deg = NaN;
end

function err = local_invariant_error(channel)
err = max(abs(channel.H_f(:) - channel.H_direct_f(:) - channel.H_reflect_f(:)));
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

function [h_reflect_sync, h_total_sync, peak_index] = local_peak_sync_taps(h_reflect_full, h_total_full)
h_reflect_full = h_reflect_full(:);
h_total_full = h_total_full(:);
tap_energy = abs(h_total_full).^2;
if all(tap_energy == 0)
    peak_index = 1;
else
    [~, peak_index] = max(tap_energy);
end
h_total_sync = [h_total_full(peak_index:end); zeros(peak_index - 1, 1)];
h_reflect_sync = [h_reflect_full(peak_index:end); zeros(peak_index - 1, 1)];
end

function metrics = local_tap_metrics(h_taps, energy_ratio)
h_taps = h_taps(:);
energy = abs(h_taps).^2;
total_energy = sum(energy);
if total_energy <= 0
    metrics = local_empty_tap_metrics();
    return
end
delay_idx = (0:(numel(h_taps) - 1)).';
mean_delay = sum(delay_idx .* energy) / total_energy;
rms_delay = sqrt(sum(((delay_idx - mean_delay).^2) .* energy) / total_energy);
[peak_energy, peak_index] = max(energy);
cum_energy = cumsum(energy);
tap_count = find(cum_energy >= max(min(energy_ratio, 1), 0) * total_energy, 1, 'first');
metrics = struct( ...
    'energy', total_energy, ...
    'mean_delay_samples', mean_delay, ...
    'rms_delay_samples', rms_delay, ...
    'peak_index', peak_index, ...
    'peak_fraction', peak_energy / total_energy, ...
    'tap_count', tap_count);
end

function metrics = local_empty_tap_metrics()
metrics = struct( ...
    'energy', NaN, ...
    'mean_delay_samples', NaN, ...
    'rms_delay_samples', NaN, ...
    'peak_index', NaN, ...
    'peak_fraction', NaN, ...
    'tap_count', NaN);
end

function [l2_error, corr_value] = local_compare_pdp(a, b)
a = a(:);
b = b(:);
a_sum = sum(a);
b_sum = sum(b);
if a_sum > 0
    an = a / a_sum;
else
    an = a;
end
if b_sum > 0
    bn = b / b_sum;
else
    bn = b;
end
l2_error = norm(an - bn) / max(norm(an), eps);
aa = an - mean(an);
bb = bn - mean(bn);
den = norm(aa) * norm(bb);
if den > 0
    corr_value = real((aa' * bb) / den);
else
    corr_value = NaN;
end
end

function rel = local_rel_error(a, b)
rel = abs(a - b) / max([abs(a), abs(b), eps]);
end

function q = local_quantiles(values, p_list)
values = sort(values(isfinite(values(:))));
q = NaN(size(p_list));
if isempty(values)
    return
end
for ii = 1:numel(p_list)
    p = p_list(ii);
    if numel(values) == 1
        q(ii) = values(1);
    else
        idx = 1 + (numel(values) - 1) * p / 100;
        lo = floor(idx);
        hi = ceil(idx);
        if lo == hi
            q(ii) = values(lo);
        else
            q(ii) = values(lo) + (idx - lo) * (values(hi) - values(lo));
        end
    end
end
end

function delta = local_ecdf_max_delta(a, b)
a = sort(a(isfinite(a(:))));
b = sort(b(isfinite(b(:))));
if isempty(a) || isempty(b)
    delta = NaN;
    return
end
grid = unique([a; b]);
Fa = zeros(size(grid));
Fb = zeros(size(grid));
ia = 1;
ib = 1;
for ii = 1:numel(grid)
    while ia <= numel(a) && a(ia) <= grid(ii)
        ia = ia + 1;
    end
    while ib <= numel(b) && b(ib) <= grid(ii)
        ib = ib + 1;
    end
    Fa(ii) = (ia - 1) / numel(a);
    Fb(ii) = (ib - 1) / numel(b);
end
delta = max(abs(Fa - Fb));
end

function value = local_metric_mean(metrics, field_name)
values = local_metric_values(metrics, field_name);
value = mean(values, 'omitnan');
end

function value = local_metric_std(metrics, field_name)
values = local_metric_values(metrics, field_name);
value = std(values, 0, 'omitnan');
end

function values = local_metric_values(metrics, field_name)
values = NaN(numel(metrics), 1);
for ii = 1:numel(metrics)
    values(ii) = metrics(ii).(field_name);
end
end

function report = local_validation_report(T, run_table, cfg)
report = struct();
report.max_kdomain_invariant_error = max(T.max_kdomain_invariant_error);
report.max_kstat_invariant_error = max(T.max_kstat_invariant_error);
report.max_kstat_energy_error = max(T.max_kstat_energy_error);
report.max_direct_path_abs_delta = max(T.max_direct_path_abs_delta);
report.max_E_abs_h_reflect2_rel_error = max(T.E_abs_h_reflect2_rel_error);
report.max_pdp_reflect_l2_error = max(T.pdp_reflect_l2_error);
report.min_pdp_reflect_corr = min(T.pdp_reflect_corr);
report.max_pdp_total_l2_error = max(T.pdp_total_l2_error);
report.min_pdp_total_corr = min(T.pdp_total_corr);
report.max_abs_h_total_quantile_nrmse = max(T.abs_h_total_quantile_nrmse);
report.max_abs_h_total_ecdf_delta = max(T.abs_h_total_ecdf_max_delta);
report.invariants_passed = report.max_kdomain_invariant_error <= cfg.invariant_tol && ...
    report.max_kstat_invariant_error <= cfg.invariant_tol;
report.kstat_energy_passed = report.max_kstat_energy_error <= cfg.energy_tol;
report.finite_summary_passed = local_table_finite(T);
report.finite_runs_passed = local_table_finite(run_table);
report.stat_E_abs_h_reflect2_passed = report.max_E_abs_h_reflect2_rel_error <= cfg.E_href2_rel_tol;
report.stat_pdp_reflect_l2_passed = report.max_pdp_reflect_l2_error <= cfg.pdp_l2_tol;
report.stat_pdp_reflect_corr_passed = report.min_pdp_reflect_corr >= cfg.pdp_corr_tol;
report.stat_abs_h_total_quantile_passed = ...
    report.max_abs_h_total_quantile_nrmse <= cfg.abs_h_total_quantile_nrmse_tol;
report.all_hard_checks_passed = report.invariants_passed && report.kstat_energy_passed && ...
    report.finite_summary_passed && report.finite_runs_passed;
report.statistical_targets_passed = report.stat_E_abs_h_reflect2_passed && ...
    report.stat_pdp_reflect_l2_passed && report.stat_pdp_reflect_corr_passed && ...
    report.stat_abs_h_total_quantile_passed;
report.interpretation = ['Hard checks test implementation invariants. Statistical targets are ', ...
    'Monte Carlo consistency targets for the current implementation, not a proof of ', ...
    'single-realization field equality or calibrated sea-surface frequency correlation.'];
end

function tf = local_table_finite(T)
tf = true;
for ii = 1:width(T)
    values = T{:, ii};
    if isnumeric(values) || islogical(values)
        tf = tf && all(isfinite(double(values(:))));
    end
end
end

function figure_files = local_write_figures(T, pdp_reflect_sync, pdp_total_sync, h_total_samples, cfg)
figure_files = strings(0, 1);

file = fullfile(cfg.figure_dir, 'kstat_kdomain_E_abs_h_reflect2_vs_Hs.png');
figure('Visible', 'off');
plot(T.Hs_m, T.E_abs_h_reflect2_kdomain, '-o', 'LineWidth', 1.4); hold on
plot(T.Hs_m, T.E_abs_h_reflect2_kstat, '--s', 'LineWidth', 1.4);
grid on
xlabel('H_s (m)');
ylabel('E[|h_{ref}|^2]');
legend('kirchhoff\_kdomain', 'kirchhoff\_kstat', 'Location', 'best');
title('Reflected receiver energy statistic');
exportgraphics(gcf, file, 'Resolution', cfg.figure_resolution);
close(gcf);
figure_files(end + 1, 1) = string(file);

file = fullfile(cfg.figure_dir, 'kstat_kdomain_abs_h_total_distribution.png');
figure('Visible', 'off');
for ih = 1:numel(cfg.Hs_list_m)
    subplot(2, ceil(numel(cfg.Hs_list_m)/2), ih);
    kd = sort(abs(h_total_samples(:, 1, ih)));
    ks = sort(abs(h_total_samples(:, 2, ih)));
    p = linspace(0, 1, numel(kd)).';
    plot(kd, p, '-o', 'LineWidth', 1.0); hold on
    plot(ks, p, '--s', 'LineWidth', 1.0);
    grid on
    xlabel('|h_{total}|');
    ylabel('ECDF');
    title(sprintf('H_s=%.3g m', cfg.Hs_list_m(ih)));
end
legend('kdomain', 'kstat', 'Location', 'best');
exportgraphics(gcf, file, 'Resolution', cfg.figure_resolution);
close(gcf);
figure_files(end + 1, 1) = string(file);

file = fullfile(cfg.figure_dir, 'kstat_kdomain_mean_pdp_reflect_compare.png');
local_plot_mean_pdp(file, pdp_reflect_sync, cfg, 'Mean reflected PDP, peak synced');
figure_files(end + 1, 1) = string(file);

file = fullfile(cfg.figure_dir, 'kstat_kdomain_mean_pdp_total_compare.png');
local_plot_mean_pdp(file, pdp_total_sync, cfg, 'Mean total PDP, peak synced');
figure_files(end + 1, 1) = string(file);

file = fullfile(cfg.figure_dir, 'kstat_kdomain_pdp_error_vs_Hs.png');
figure('Visible', 'off');
plot(T.Hs_m, T.pdp_reflect_l2_error, '-o', 'LineWidth', 1.4); hold on
plot(T.Hs_m, T.pdp_total_l2_error, '--s', 'LineWidth', 1.4);
grid on
xlabel('H_s (m)');
ylabel('normalized L2 error');
legend('reflected PDP', 'total PDP', 'Location', 'best');
title('Mean PDP shape error');
exportgraphics(gcf, file, 'Resolution', cfg.figure_resolution);
close(gcf);
figure_files(end + 1, 1) = string(file);
end

function local_plot_mean_pdp(file, pdp_sync, cfg, fig_title)
figure('Visible', 'off');
for ih = 1:numel(cfg.Hs_list_m)
    subplot(2, ceil(numel(cfg.Hs_list_m)/2), ih);
    kd = mean(pdp_sync(:, :, 1, ih), 2);
    ks = mean(pdp_sync(:, :, 2, ih), 2);
    if sum(kd) > 0
        kd = kd / sum(kd);
    end
    if sum(ks) > 0
        ks = ks / sum(ks);
    end
    n_show = min(120, numel(kd));
    semilogy(0:(n_show - 1), max(kd(1:n_show), realmin), '-o', 'LineWidth', 1.0);
    hold on
    semilogy(0:(n_show - 1), max(ks(1:n_show), realmin), '--s', 'LineWidth', 1.0);
    grid on
    xlabel('tap index after peak sync');
    ylabel('normalized PDP');
    title(sprintf('H_s=%.3g m', cfg.Hs_list_m(ih)));
end
legend('kdomain', 'kstat', 'Location', 'best');
sgtitle(fig_title);
exportgraphics(gcf, file, 'Resolution', cfg.figure_resolution);
close(gcf);
end

function local_write_report(report_file, report, T, cfg, figure_files)
fid = fopen(report_file, 'w');
if fid < 0
    error('Cannot open report file %s for writing.', report_file);
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '# K-Stat vs K-Domain Channel Statistics Validation\n\n');
fprintf(fid, 'Generated: %s\n\n', datestr(now, 31));
fprintf(fid, '## Setup\n\n');
fprintf(fid, '- Full `vertical_channel_model` path with PE/WAPE propagation.\n');
fprintf(fid, '- No modulation, noise, demodulation, or BER/SER loop.\n');
fprintf(fid, '- Hs values: `%s` m; seed count: `%d`; grid `%d x %d`; Nf `%d`.\n', ...
    mat2str(cfg.Hs_list_m), cfg.seed_count, cfg.nx, cfg.ny, cfg.Nf);
fprintf(fid, '- Wideband band: `%s` Hz; reference frequency: `%.3f` Hz.\n', ...
    mat2str(cfg.f_band_hz), cfg.f_ref_hz);
fprintf(fid, '- Roughness mode: `%s`; PM shape wind: `%.3g m/s`.\n\n', ...
    cfg.roughness_scale_mode, cfg.wind_mps);

fprintf(fid, '## Key Results\n\n');
fprintf(fid, '- Max H invariant error, kdomain: `%.4g`.\n', report.max_kdomain_invariant_error);
fprintf(fid, '- Max H invariant error, kstat: `%.4g`.\n', report.max_kstat_invariant_error);
fprintf(fid, '- Max K-Stat phase-screen energy error: `%.4g`.\n', report.max_kstat_energy_error);
fprintf(fid, '- Max direct-path branch/seed delta: `%.4g`.\n', report.max_direct_path_abs_delta);
fprintf(fid, '- Max E[|h_ref|^2] relative error: `%.4g`.\n', report.max_E_abs_h_reflect2_rel_error);
fprintf(fid, '- Max reflected mean-PDP L2 error: `%.4g`.\n', report.max_pdp_reflect_l2_error);
fprintf(fid, '- Min reflected mean-PDP correlation: `%.4g`.\n', report.min_pdp_reflect_corr);
fprintf(fid, '- Max |h_total| quantile NRMSE: `%.4g`.\n', report.max_abs_h_total_quantile_nrmse);
fprintf(fid, '- Hard checks passed: `%d`.\n', report.all_hard_checks_passed);
fprintf(fid, '- Statistical targets passed: `%d`.\n\n', report.statistical_targets_passed);

fprintf(fid, '## Interpretation\n\n');
fprintf(fid, '- `kirchhoff_kdomain` is an ensemble of explicit sea-surface phase-screen realizations.\n');
fprintf(fid, '- `kirchhoff_kstat` is an ensemble of statistical phase-screen reflected-field realizations.\n');
fprintf(fid, '- Single-seed `h_ref`, `h_total`, or PDP samples are not expected to match pointwise.\n');
fprintf(fid, '- Valid comparison targets are `E[|h_ref|^2]`, mean PDP shape, `|h_total|` distribution, and implementation invariants.\n');
fprintf(fid, '- PDP is a wideband receiver-side statistic. The current kstat frequency correlation is still an engineering simplification, so PDP agreement is an implementation-level statistical check, not a calibrated sea-surface time-frequency model.\n\n');

fprintf(fid, '## Figures\n\n');
for ii = 1:numel(figure_files)
    fprintf(fid, '- `%s`\n', char(figure_files(ii)));
end

fprintf(fid, '\n## Summary Table\n\n');
preview = T(:, {'Hs_m', 'seed_count', 'E_abs_h_reflect2_rel_error', ...
    'pdp_reflect_l2_error', 'pdp_reflect_corr', 'pdp_total_l2_error', ...
    'pdp_total_corr', 'abs_h_total_quantile_nrmse', 'max_kstat_energy_error'});
fprintf(fid, '%s\n', evalc('disp(preview)'));
end

function values = local_env_numeric_vector(name, default_value)
raw = getenv(name);
if isempty(raw)
    values = default_value;
    return
end
values = str2num(raw); %#ok<ST2NM>
if isempty(values) || ~isnumeric(values) || any(~isfinite(values(:)))
    error('Environment variable %s must be a finite numeric vector.', name);
end
values = values(:).';
end

function value = local_env_scalar(name, default_value)
values = local_env_numeric_vector(name, default_value);
value = values(1);
end

function value = local_env_int(name, default_value)
value = round(local_env_scalar(name, default_value));
end
