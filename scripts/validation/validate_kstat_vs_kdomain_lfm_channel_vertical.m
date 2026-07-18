run(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'bootstrap_project.m'));
%VALIDATE_KSTAT_VS_KDOMAIN_LFM_CHANNEL_VERTICAL
% Channel-level LFM test-signal validation of Kirchhoff K-Stat against
% explicit Kirchhoff k-domain ensemble statistics after PE/WAPE propagation.
%
% This script does not call comm_main_vertical_psk and does not use PSK,
% BER/SER, noise injection, or equalization. It sends an analytic baseband
% LFM through the channel frequency responses H_f and H_reflect_f.

clear
clc
format compact

cfg = struct();
cfg.Hs_list_m = local_env_numeric_vector('KSTAT_LFM_HS_LIST', [0.05, 0.1, 0.2, 0.5]);
cfg.seed_count = local_env_int('KSTAT_LFM_SEED_COUNT', 64);
cfg.grid_n = local_env_int('KSTAT_LFM_GRID_N', 128);
cfg.nx = cfg.grid_n;
cfg.ny = cfg.grid_n;
cfg.xw_m = local_env_scalar('KSTAT_LFM_XW_M', 50 * cfg.nx / 128);
cfg.yw_m = local_env_scalar('KSTAT_LFM_YW_M', 50 * cfg.ny / 128);
cfg.Nf = local_env_int('KSTAT_LFM_NF', 32);
cfg.f_band_hz = local_env_numeric_vector('KSTAT_LFM_F_BAND', [4000, 8000]);
cfg.f_ref_hz = local_env_scalar('KSTAT_LFM_F_REF_HZ', 6000);
cfg.wind_mps = local_env_scalar('KSTAT_LFM_WIND_MPS', 5.0);
cfg.seed_list = 12345 + (0:(cfg.seed_count - 1));

cfg.lfm_fs_hz = local_env_scalar('KSTAT_LFM_FS_HZ', 48000);
cfg.lfm_duration_s = local_env_scalar('KSTAT_LFM_DURATION_S', 0.05);
cfg.lfm_taper_fraction = local_env_scalar('KSTAT_LFM_TAPER_FRACTION', 0.05);
cfg.representative_Hs_m = local_env_scalar('KSTAT_LFM_REP_HS', 0.2);
cfg.representative_seed = local_env_int('KSTAT_LFM_REP_SEED', cfg.seed_list(1));

cfg.roughness_scale_mode = 'target_hs';
cfg.invariant_tol = 1e-10;
cfg.energy_tol = 1e-8;
cfg.E_href2_rel_tol = 0.25;
cfg.pdp_corr_tol = 0.90;
cfg.rx_envelope_corr_tol = 0.90;
cfg.mf_peak_rel_tol = 0.25;
cfg.figure_resolution = 160;

if numel(cfg.f_band_hz) ~= 2 || cfg.f_band_hz(2) <= cfg.f_band_hz(1)
    error('KSTAT_LFM_F_BAND must contain [f_min f_max] with f_max > f_min.');
end
if cfg.f_ref_hz <= cfg.f_band_hz(1) || cfg.f_ref_hz >= cfg.f_band_hz(2)
    error('KSTAT_LFM_F_REF_HZ must lie inside KSTAT_LFM_F_BAND.');
end
cfg.seed_count = max(1, cfg.seed_count);
cfg.Nf = max(2, cfg.Nf);
cfg.lfm_fs_hz = max(cfg.lfm_fs_hz, 2 * diff(cfg.f_band_hz));
cfg.lfm_duration_s = max(cfg.lfm_duration_s, 1 / diff(cfg.f_band_hz));

cfg.result_mat = project_result_file('validation', ...
    'validate_kstat_vs_kdomain_lfm_channel_vertical_result.mat');
cfg.summary_csv = project_result_file('validation', ...
    'validate_kstat_vs_kdomain_lfm_channel_vertical_summary.csv');
cfg.runs_csv = project_result_file('validation', ...
    'validate_kstat_vs_kdomain_lfm_channel_vertical_runs.csv');
cfg.report_file = project_result_file('validation', ...
    'validate_kstat_vs_kdomain_lfm_channel_vertical_report.md');
cfg.figure_dir = project_result_dir('validation');

[lfm, cfg] = local_build_lfm(cfg);
branch_names = {'kirchhoff_kdomain', 'kirchhoff_kstat'};
n_hs = numel(cfg.Hs_list_m);
n_branch = numel(branch_names);
n_seed = cfg.seed_count;
n_sig = numel(lfm.tx_bb);

fprintf('K-Stat vs K-Domain LFM channel validation\n');
fprintf('  Hs: %s m\n', mat2str(cfg.Hs_list_m));
fprintf('  seeds: %d, grid: %d x %d, aperture: %.3g x %.3g m, Nf: %d\n', ...
    n_seed, cfg.nx, cfg.ny, cfg.xw_m, cfg.yw_m, cfg.Nf);
fprintf('  LFM fs=%.3g Hz, duration=%.4g s, samples=%d, band=%s Hz\n', ...
    cfg.lfm_fs_hz, cfg.lfm_duration_s, n_sig, mat2str(cfg.f_band_hz));

h_reflect_samples = complex(NaN(n_seed, n_branch, n_hs), NaN(n_seed, n_branch, n_hs));
h_total_samples = complex(NaN(n_seed, n_branch, n_hs), NaN(n_seed, n_branch, n_hs));
rx_total_samples = complex(NaN(n_sig, n_seed, n_branch, n_hs), NaN(n_sig, n_seed, n_branch, n_hs));
rx_reflect_samples = complex(NaN(n_sig, n_seed, n_branch, n_hs), NaN(n_sig, n_seed, n_branch, n_hs));
mf_total_samples = complex(NaN(n_sig, n_seed, n_branch, n_hs), NaN(n_sig, n_seed, n_branch, n_hs));
mf_reflect_samples = complex(NaN(n_sig, n_seed, n_branch, n_hs), NaN(n_sig, n_seed, n_branch, n_hs));
pdp_total_sync = NaN(n_sig, n_seed, n_branch, n_hs);
pdp_reflect_sync = NaN(n_sig, n_seed, n_branch, n_hs);
mf_total_metrics = repmat(local_empty_mf_metrics(), n_seed, n_branch, n_hs);
mf_reflect_metrics = repmat(local_empty_mf_metrics(), n_seed, n_branch, n_hs);

run_rows = struct([]);
summary_rows = struct([]);
case_results = struct([]);
run_idx = 0;
summary_idx = 0;
case_idx = 0;

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
                error('validate_kstat_vs_kdomain_lfm_channel_vertical:Invariant', ...
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
                    error('validate_kstat_vs_kdomain_lfm_channel_vertical:KStatEnergy', ...
                        'K-Stat phase-screen energy failed Hs %.4g seed %d: %.3e', ...
                        Hs_m, seed, kstat_energy_error);
                end
            end

            H_total_shift = local_interpolate_channel(channel.f_axis, channel.H_f, ...
                channel.idx_f_ref, lfm.f_rel_shifted);
            H_reflect_shift = local_interpolate_channel(channel.f_axis, channel.H_reflect_f, ...
                channel.idx_f_ref, lfm.f_rel_shifted);
            rx_total = local_apply_freq_response(lfm.tx_bb, H_total_shift);
            rx_reflect = local_apply_freq_response(lfm.tx_bb, H_reflect_shift);
            mf_total = local_matched_filter(rx_total, lfm.tx_bb);
            mf_reflect = local_matched_filter(rx_reflect, lfm.tx_bb);
            h_total_taps = ifft(ifftshift(H_total_shift));
            h_reflect_taps = ifft(ifftshift(H_reflect_shift));
            [h_reflect_sync, h_total_sync] = local_peak_sync_taps(h_reflect_taps, h_total_taps);

            mt = local_mf_metrics(mf_total, cfg.lfm_fs_hz);
            mr = local_mf_metrics(mf_reflect, cfg.lfm_fs_hz);
            mf_total_metrics(iseed, ib, ih) = mt;
            mf_reflect_metrics(iseed, ib, ih) = mr;
            h_reflect_samples(iseed, ib, ih) = channel.h_reflect;
            h_total_samples(iseed, ib, ih) = channel.h_total;
            rx_total_samples(:, iseed, ib, ih) = rx_total;
            rx_reflect_samples(:, iseed, ib, ih) = rx_reflect;
            mf_total_samples(:, iseed, ib, ih) = mf_total;
            mf_reflect_samples(:, iseed, ib, ih) = mf_reflect;
            pdp_total_sync(:, iseed, ib, ih) = abs(h_total_sync).^2;
            pdp_reflect_sync(:, iseed, ib, ih) = abs(h_reflect_sync).^2;

            run_idx = run_idx + 1;
            run_rows(run_idx).Hs_m = Hs_m; %#ok<SAGROW>
            run_rows(run_idx).seed = seed;
            run_rows(run_idx).branch_code = ib;
            run_rows(run_idx).abs_h_reflect2 = abs(channel.h_reflect)^2;
            run_rows(run_idx).abs_h_total = abs(channel.h_total);
            run_rows(run_idx).rx_total_energy = sum(abs(rx_total).^2);
            run_rows(run_idx).rx_reflect_energy = sum(abs(rx_reflect).^2);
            run_rows(run_idx).mf_total_peak_abs = mt.peak_abs;
            run_rows(run_idx).mf_reflect_peak_abs = mr.peak_abs;
            run_rows(run_idx).mf_total_peak_delay_samples = mt.peak_delay_samples;
            run_rows(run_idx).mf_reflect_peak_delay_samples = mr.peak_delay_samples;
            run_rows(run_idx).mf_total_mainlobe_width_samples = mt.mainlobe_width_samples;
            run_rows(run_idx).mf_reflect_mainlobe_width_samples = mr.mainlobe_width_samples;
            run_rows(run_idx).mf_total_peak_sidelobe_db = mt.peak_sidelobe_db;
            run_rows(run_idx).mf_reflect_peak_sidelobe_db = mr.peak_sidelobe_db;
            run_rows(run_idx).invariant_error = invariant_error;
            run_rows(run_idx).direct_path_max_abs_delta = direct_delta;
            run_rows(run_idx).kstat_phase_screen_energy_error = kstat_energy_error;

            if local_is_representative_case(cfg, Hs_m, seed)
                case_idx = case_idx + 1;
                case_results(case_idx).Hs_m = Hs_m; %#ok<SAGROW>
                case_results(case_idx).seed = seed;
                case_results(case_idx).branch_code = ib;
                case_results(case_idx).branch_name = branch;
                case_results(case_idx).f_axis = channel.f_axis;
                case_results(case_idx).H_f = channel.H_f;
                case_results(case_idx).H_reflect_f = channel.H_reflect_f;
                case_results(case_idx).H_total_shift = H_total_shift;
                case_results(case_idx).H_reflect_shift = H_reflect_shift;
                case_results(case_idx).rx_total = rx_total;
                case_results(case_idx).rx_reflect = rx_reflect;
                case_results(case_idx).mf_total = mf_total;
                case_results(case_idx).mf_reflect = mf_reflect;
                case_results(case_idx).pdp_total = abs(h_total_sync).^2;
                case_results(case_idx).pdp_reflect = abs(h_reflect_sync).^2;
            end
        end
    end

    kd = 1;
    ks = 2;
    kd_h_ref = h_reflect_samples(:, kd, ih);
    ks_h_ref = h_reflect_samples(:, ks, ih);
    kd_abs_total = abs(h_total_samples(:, kd, ih));
    ks_abs_total = abs(h_total_samples(:, ks, ih));

    mean_rx_total_env_kd = mean(abs(rx_total_samples(:, :, kd, ih)), 2);
    mean_rx_total_env_ks = mean(abs(rx_total_samples(:, :, ks, ih)), 2);
    mean_rx_reflect_env_kd = mean(abs(rx_reflect_samples(:, :, kd, ih)), 2);
    mean_rx_reflect_env_ks = mean(abs(rx_reflect_samples(:, :, ks, ih)), 2);
    [rx_total_env_l2, rx_total_env_corr] = local_compare_shape(mean_rx_total_env_kd, mean_rx_total_env_ks);
    [rx_reflect_env_l2, rx_reflect_env_corr] = local_compare_shape(mean_rx_reflect_env_kd, mean_rx_reflect_env_ks);

    mean_pdp_total_kd = mean(pdp_total_sync(:, :, kd, ih), 2);
    mean_pdp_total_ks = mean(pdp_total_sync(:, :, ks, ih), 2);
    mean_pdp_reflect_kd = mean(pdp_reflect_sync(:, :, kd, ih), 2);
    mean_pdp_reflect_ks = mean(pdp_reflect_sync(:, :, ks, ih), 2);
    [pdp_total_l2, pdp_total_corr] = local_compare_shape(mean_pdp_total_kd, mean_pdp_total_ks);
    [pdp_reflect_l2, pdp_reflect_corr] = local_compare_shape(mean_pdp_reflect_kd, mean_pdp_reflect_ks);

    kd_rx_total_energy = squeeze(sum(abs(rx_total_samples(:, :, kd, ih)).^2, 1));
    ks_rx_total_energy = squeeze(sum(abs(rx_total_samples(:, :, ks, ih)).^2, 1));
    kd_rx_reflect_energy = squeeze(sum(abs(rx_reflect_samples(:, :, kd, ih)).^2, 1));
    ks_rx_reflect_energy = squeeze(sum(abs(rx_reflect_samples(:, :, ks, ih)).^2, 1));
    q_kd = local_quantiles(kd_abs_total, [10, 50, 90]);
    q_ks = local_quantiles(ks_abs_total, [10, 50, 90]);
    quantile_nrmse = sqrt(mean((q_kd - q_ks).^2)) / max(mean([q_kd(:); q_ks(:)]), eps);

    summary_idx = summary_idx + 1;
    summary_rows(summary_idx).Hs_m = Hs_m; %#ok<SAGROW>
    summary_rows(summary_idx).seed_count = n_seed;
    summary_rows(summary_idx).nx = cfg.nx;
    summary_rows(summary_idx).ny = cfg.ny;
    summary_rows(summary_idx).Nf = cfg.Nf;
    summary_rows(summary_idx).lfm_fs_hz = cfg.lfm_fs_hz;
    summary_rows(summary_idx).lfm_duration_s = cfg.lfm_duration_s;
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
    summary_rows(summary_idx).rx_total_energy_mean_kdomain = mean(kd_rx_total_energy);
    summary_rows(summary_idx).rx_total_energy_mean_kstat = mean(ks_rx_total_energy);
    summary_rows(summary_idx).rx_total_energy_rel_error = local_rel_error( ...
        summary_rows(summary_idx).rx_total_energy_mean_kdomain, ...
        summary_rows(summary_idx).rx_total_energy_mean_kstat);
    summary_rows(summary_idx).rx_reflect_energy_mean_kdomain = mean(kd_rx_reflect_energy);
    summary_rows(summary_idx).rx_reflect_energy_mean_kstat = mean(ks_rx_reflect_energy);
    summary_rows(summary_idx).rx_reflect_energy_rel_error = local_rel_error( ...
        summary_rows(summary_idx).rx_reflect_energy_mean_kdomain, ...
        summary_rows(summary_idx).rx_reflect_energy_mean_kstat);
    summary_rows(summary_idx).rx_total_env_l2_error = rx_total_env_l2;
    summary_rows(summary_idx).rx_total_env_corr = rx_total_env_corr;
    summary_rows(summary_idx).rx_reflect_env_l2_error = rx_reflect_env_l2;
    summary_rows(summary_idx).rx_reflect_env_corr = rx_reflect_env_corr;
    summary_rows(summary_idx).pdp_total_l2_error = pdp_total_l2;
    summary_rows(summary_idx).pdp_total_corr = pdp_total_corr;
    summary_rows(summary_idx).pdp_reflect_l2_error = pdp_reflect_l2;
    summary_rows(summary_idx).pdp_reflect_corr = pdp_reflect_corr;
    summary_rows(summary_idx).mf_total_peak_mean_kdomain = local_metric_mean(mf_total_metrics(:, kd, ih), 'peak_abs');
    summary_rows(summary_idx).mf_total_peak_std_kdomain = local_metric_std(mf_total_metrics(:, kd, ih), 'peak_abs');
    summary_rows(summary_idx).mf_total_peak_mean_kstat = local_metric_mean(mf_total_metrics(:, ks, ih), 'peak_abs');
    summary_rows(summary_idx).mf_total_peak_std_kstat = local_metric_std(mf_total_metrics(:, ks, ih), 'peak_abs');
    summary_rows(summary_idx).mf_total_peak_rel_error = local_rel_error( ...
        summary_rows(summary_idx).mf_total_peak_mean_kdomain, ...
        summary_rows(summary_idx).mf_total_peak_mean_kstat);
    summary_rows(summary_idx).mf_reflect_peak_mean_kdomain = local_metric_mean(mf_reflect_metrics(:, kd, ih), 'peak_abs');
    summary_rows(summary_idx).mf_reflect_peak_std_kdomain = local_metric_std(mf_reflect_metrics(:, kd, ih), 'peak_abs');
    summary_rows(summary_idx).mf_reflect_peak_mean_kstat = local_metric_mean(mf_reflect_metrics(:, ks, ih), 'peak_abs');
    summary_rows(summary_idx).mf_reflect_peak_std_kstat = local_metric_std(mf_reflect_metrics(:, ks, ih), 'peak_abs');
    summary_rows(summary_idx).mf_reflect_peak_rel_error = local_rel_error( ...
        summary_rows(summary_idx).mf_reflect_peak_mean_kdomain, ...
        summary_rows(summary_idx).mf_reflect_peak_mean_kstat);
    summary_rows(summary_idx).mf_total_delay_mean_kdomain = local_metric_mean(mf_total_metrics(:, kd, ih), 'peak_delay_samples');
    summary_rows(summary_idx).mf_total_delay_mean_kstat = local_metric_mean(mf_total_metrics(:, ks, ih), 'peak_delay_samples');
    summary_rows(summary_idx).mf_reflect_delay_mean_kdomain = local_metric_mean(mf_reflect_metrics(:, kd, ih), 'peak_delay_samples');
    summary_rows(summary_idx).mf_reflect_delay_mean_kstat = local_metric_mean(mf_reflect_metrics(:, ks, ih), 'peak_delay_samples');
    summary_rows(summary_idx).mf_total_mainlobe_width_mean_kdomain = local_metric_mean(mf_total_metrics(:, kd, ih), 'mainlobe_width_samples');
    summary_rows(summary_idx).mf_total_mainlobe_width_mean_kstat = local_metric_mean(mf_total_metrics(:, ks, ih), 'mainlobe_width_samples');
    summary_rows(summary_idx).mf_reflect_mainlobe_width_mean_kdomain = local_metric_mean(mf_reflect_metrics(:, kd, ih), 'mainlobe_width_samples');
    summary_rows(summary_idx).mf_reflect_mainlobe_width_mean_kstat = local_metric_mean(mf_reflect_metrics(:, ks, ih), 'mainlobe_width_samples');
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
figure_files = local_write_figures(summary_table, lfm, rx_total_samples, ...
    rx_reflect_samples, mf_total_samples, mf_reflect_samples, ...
    pdp_total_sync, pdp_reflect_sync, case_results, cfg);
local_write_report(cfg.report_file, validation_report, summary_table, cfg, figure_files);

save(cfg.result_mat, 'cfg', 'lfm', 'run_table', 'summary_table', ...
    'validation_report', 'figure_files', 'h_reflect_samples', 'h_total_samples', ...
    'rx_total_samples', 'rx_reflect_samples', 'mf_total_samples', ...
    'mf_reflect_samples', 'pdp_total_sync', 'pdp_reflect_sync', ...
    'mf_total_metrics', 'mf_reflect_metrics', 'case_results', '-v7.3');
writetable(summary_table, cfg.summary_csv);
writetable(run_table, cfg.runs_csv);

disp(summary_table(:, {'Hs_m', 'seed_count', 'E_abs_h_reflect2_rel_error', ...
    'rx_total_env_corr', 'rx_reflect_env_corr', 'pdp_total_corr', ...
    'pdp_reflect_corr', 'mf_total_peak_rel_error', ...
    'mf_reflect_peak_rel_error', 'max_kstat_energy_error'}));
disp(validation_report);
fprintf('Saved %s\n', cfg.result_mat);
fprintf('Saved %s\n', cfg.summary_csv);
fprintf('Saved %s\n', cfg.runs_csv);
fprintf('Saved %s\n', cfg.report_file);

if ~validation_report.all_hard_checks_passed
    error('validate_kstat_vs_kdomain_lfm_channel_vertical:HardChecks', ...
        'One or more hard checks failed. See validation_report.');
end

function [lfm, cfg] = local_build_lfm(cfg)
bandwidth_hz = diff(cfg.f_band_hz);
n = max(8, round(cfg.lfm_duration_s * cfg.lfm_fs_hz));
cfg.lfm_duration_s = n / cfg.lfm_fs_hz;
t = (0:(n - 1)).' / cfg.lfm_fs_hz;
tau = t - cfg.lfm_duration_s / 2;
mu = bandwidth_hz / cfg.lfm_duration_s;
window = local_tukey_window(n, cfg.lfm_taper_fraction);
tx_bb = window .* exp(1i * pi * mu .* tau.^2);
tx_bb = tx_bb / sqrt(max(mean(abs(tx_bb).^2), eps));
tx_passband = real(tx_bb .* exp(1i * 2*pi*cfg.f_ref_hz .* t));
f_rel_shifted = ((0:n-1).' - floor(n/2)) * (cfg.lfm_fs_hz / n);
lfm = struct('t_s', t, 'tx_bb', tx_bb, 'tx_passband', tx_passband, ...
    'window', window, 'mu_hz_per_s', mu, 'bandwidth_hz', bandwidth_hz, ...
    'f_rel_shifted', f_rel_shifted, 'fs_hz', cfg.lfm_fs_hz);
end

function w = local_tukey_window(n, alpha)
alpha = max(0, min(alpha, 1));
w = ones(n, 1);
if alpha <= 0
    return
end
if alpha >= 1
    idx = (0:n-1).';
    w = 0.5 - 0.5*cos(2*pi*idx/(n-1));
    return
end
idx = (0:n-1).' / (n-1);
first = idx < alpha/2;
last = idx >= 1 - alpha/2;
w(first) = 0.5 * (1 + cos(2*pi/alpha * (idx(first) - alpha/2)));
w(last) = 0.5 * (1 + cos(2*pi/alpha * (idx(last) - 1 + alpha/2)));
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

function H_shift = local_interpolate_channel(f_axis, H_f, idx_f_ref, f_rel_shifted)
f_axis = f_axis(:);
H_f = H_f(:);
fc = f_axis(idx_f_ref);
f_rel = f_axis - fc;
H_shift = interp1(f_rel, H_f, f_rel_shifted(:), 'linear', 0);
end

function rx = local_apply_freq_response(tx, H_shift)
rx = ifft(fft(tx(:)) .* ifftshift(H_shift(:)));
end

function mf = local_matched_filter(rx, tx)
mf = ifft(fft(rx(:)) .* conj(fft(tx(:))));
end

function [h_reflect_sync, h_total_sync] = local_peak_sync_taps(h_reflect, h_total)
h_reflect = h_reflect(:);
h_total = h_total(:);
energy = abs(h_total).^2;
if all(energy == 0)
    peak_idx = 1;
else
    [~, peak_idx] = max(energy);
end
h_total_sync = [h_total(peak_idx:end); zeros(peak_idx - 1, 1)];
h_reflect_sync = [h_reflect(peak_idx:end); zeros(peak_idx - 1, 1)];
end

function tf = local_is_representative_case(cfg, Hs_m, seed)
[~, ih] = min(abs(cfg.Hs_list_m - cfg.representative_Hs_m));
rep_Hs = cfg.Hs_list_m(ih);
tf = abs(Hs_m - rep_Hs) <= 10*eps(max(1, abs(rep_Hs))) && seed == cfg.representative_seed;
end

function err = local_invariant_error(channel)
err = max(abs(channel.H_f(:) - channel.H_direct_f(:) - channel.H_reflect_f(:)));
end

function metrics = local_mf_metrics(mf, fs_hz)
mag = abs(mf(:));
n = numel(mag);
if all(mag == 0)
    metrics = local_empty_mf_metrics();
    return
end
[peak_abs, peak_idx] = max(mag);
half = 0.5 * peak_abs;
above = mag >= half;
left_count = 0;
idx = peak_idx;
while left_count < n
    idx = mod(idx - 2, n) + 1;
    if ~above(idx)
        break
    end
    left_count = left_count + 1;
end
right_count = 0;
idx = peak_idx;
while right_count < n
    idx = mod(idx, n) + 1;
    if ~above(idx)
        break
    end
    right_count = right_count + 1;
end
width = 1 + left_count + right_count;
main_mask = false(n, 1);
for kk = -left_count:right_count
    main_mask(mod(peak_idx - 1 + kk, n) + 1) = true;
end
side = mag(~main_mask);
if isempty(side)
    peak_sidelobe_db = -300;
else
    peak_sidelobe_db = 20*log10(max(side) / max(peak_abs, eps));
end
metrics = struct('peak_abs', peak_abs, ...
    'peak_index', peak_idx, ...
    'peak_delay_samples', peak_idx - 1, ...
    'peak_delay_s', (peak_idx - 1) / fs_hz, ...
    'mainlobe_width_samples', width, ...
    'mainlobe_width_s', width / fs_hz, ...
    'energy', sum(mag.^2), ...
    'peak_sidelobe_db', peak_sidelobe_db);
end

function metrics = local_empty_mf_metrics()
metrics = struct('peak_abs', NaN, ...
    'peak_index', NaN, ...
    'peak_delay_samples', NaN, ...
    'peak_delay_s', NaN, ...
    'mainlobe_width_samples', NaN, ...
    'mainlobe_width_s', NaN, ...
    'energy', NaN, ...
    'peak_sidelobe_db', NaN);
end

function [l2_error, corr_value] = local_compare_shape(a, b)
a = a(:);
b = b(:);
if sum(abs(a)) > 0
    an = a / sum(abs(a));
else
    an = a;
end
if sum(abs(b)) > 0
    bn = b / sum(abs(b));
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
report.min_rx_total_env_corr = min(T.rx_total_env_corr);
report.min_rx_reflect_env_corr = min(T.rx_reflect_env_corr);
report.min_pdp_total_corr = min(T.pdp_total_corr);
report.min_pdp_reflect_corr = min(T.pdp_reflect_corr);
report.max_mf_total_peak_rel_error = max(T.mf_total_peak_rel_error);
report.max_mf_reflect_peak_rel_error = max(T.mf_reflect_peak_rel_error);
report.invariants_passed = report.max_kdomain_invariant_error <= cfg.invariant_tol && ...
    report.max_kstat_invariant_error <= cfg.invariant_tol;
report.kstat_energy_passed = report.max_kstat_energy_error <= cfg.energy_tol;
report.finite_summary_passed = local_table_finite(T);
report.finite_runs_passed = local_table_finite(run_table);
report.stat_E_abs_h_reflect2_passed = report.max_E_abs_h_reflect2_rel_error <= cfg.E_href2_rel_tol;
report.stat_rx_total_env_corr_passed = report.min_rx_total_env_corr >= cfg.rx_envelope_corr_tol;
report.stat_rx_reflect_env_corr_passed = report.min_rx_reflect_env_corr >= cfg.rx_envelope_corr_tol;
report.stat_pdp_total_corr_passed = report.min_pdp_total_corr >= cfg.pdp_corr_tol;
report.stat_pdp_reflect_corr_passed = report.min_pdp_reflect_corr >= cfg.pdp_corr_tol;
report.stat_mf_peak_passed = report.max_mf_total_peak_rel_error <= cfg.mf_peak_rel_tol && ...
    report.max_mf_reflect_peak_rel_error <= cfg.mf_peak_rel_tol;
report.all_hard_checks_passed = report.invariants_passed && report.kstat_energy_passed && ...
    report.finite_summary_passed && report.finite_runs_passed;
report.statistical_targets_passed = report.stat_E_abs_h_reflect2_passed && ...
    report.stat_rx_total_env_corr_passed && report.stat_rx_reflect_env_corr_passed && ...
    report.stat_pdp_total_corr_passed && report.stat_pdp_reflect_corr_passed && ...
    report.stat_mf_peak_passed;
report.interpretation = ['This LFM validation uses H_f/H_reflect_f as channel ', ...
    'frequency responses. It checks receiver-side waveform and matched-filter ', ...
    'statistics, not pointwise equality between single realizations.'];
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

function figure_files = local_write_figures(T, lfm, rx_total_samples, ...
    rx_reflect_samples, mf_total_samples, mf_reflect_samples, ...
    pdp_total_sync, pdp_reflect_sync, case_results, cfg)
figure_files = strings(0, 1);

file = fullfile(cfg.figure_dir, 'lfm_tx_signal_time_frequency.png');
local_plot_lfm_tx(file, lfm, cfg);
figure_files(end + 1, 1) = string(file);

file = fullfile(cfg.figure_dir, 'lfm_rx_total_waveform_compare.png');
local_plot_representative_waveform(file, case_results, lfm, 'rx_total', 'Total received LFM');
figure_files(end + 1, 1) = string(file);

file = fullfile(cfg.figure_dir, 'lfm_rx_reflect_waveform_compare.png');
local_plot_representative_waveform(file, case_results, lfm, 'rx_reflect', 'Reflected received LFM');
figure_files(end + 1, 1) = string(file);

file = fullfile(cfg.figure_dir, 'lfm_matched_filter_total_compare.png');
local_plot_representative_mf(file, case_results, lfm, 'mf_total', 'Total matched filter');
figure_files(end + 1, 1) = string(file);

file = fullfile(cfg.figure_dir, 'lfm_matched_filter_reflect_compare.png');
local_plot_representative_mf(file, case_results, lfm, 'mf_reflect', 'Reflected matched filter');
figure_files(end + 1, 1) = string(file);

file = fullfile(cfg.figure_dir, 'lfm_ensemble_rx_energy_vs_Hs.png');
figure('Visible', 'off');
plot(T.Hs_m, T.rx_total_energy_mean_kdomain, '-o', 'LineWidth', 1.3); hold on
plot(T.Hs_m, T.rx_total_energy_mean_kstat, '--s', 'LineWidth', 1.3);
plot(T.Hs_m, T.rx_reflect_energy_mean_kdomain, '-.^', 'LineWidth', 1.3);
plot(T.Hs_m, T.rx_reflect_energy_mean_kstat, ':v', 'LineWidth', 1.3);
grid on
xlabel('H_s (m)');
ylabel('LFM received energy');
legend('total kdomain', 'total kstat', 'reflect kdomain', 'reflect kstat', 'Location', 'best');
title('LFM received energy vs H_s');
exportgraphics(gcf, file, 'Resolution', cfg.figure_resolution);
close(gcf);
figure_files(end + 1, 1) = string(file);

file = fullfile(cfg.figure_dir, 'lfm_pdp_and_mf_error_vs_Hs.png');
figure('Visible', 'off');
plot(T.Hs_m, T.pdp_reflect_l2_error, '-o', 'LineWidth', 1.3); hold on
plot(T.Hs_m, T.rx_reflect_env_l2_error, '--s', 'LineWidth', 1.3);
plot(T.Hs_m, T.mf_reflect_peak_rel_error, '-.^', 'LineWidth', 1.3);
grid on
xlabel('H_s (m)');
ylabel('relative / normalized error');
legend('reflected PDP L2', 'reflected envelope L2', 'reflected MF peak rel', 'Location', 'best');
title('LFM/PDP error trends');
exportgraphics(gcf, file, 'Resolution', cfg.figure_resolution);
close(gcf);
figure_files(end + 1, 1) = string(file);

file = fullfile(cfg.figure_dir, 'lfm_Hf_magnitude_phase_representative.png');
local_plot_representative_Hf(file, case_results);
figure_files(end + 1, 1) = string(file);

% Use these arrays so MATLAB does not warn when checking function inputs.
unused = {rx_total_samples, rx_reflect_samples, mf_total_samples, mf_reflect_samples, ...
    pdp_total_sync, pdp_reflect_sync}; %#ok<NASGU>
end

function local_plot_lfm_tx(file, lfm, cfg)
figure('Visible', 'off');
t_ms = lfm.t_s * 1e3;
subplot(3, 1, 1);
plot(t_ms, real(lfm.tx_bb), 'LineWidth', 1.0); hold on
plot(t_ms, imag(lfm.tx_bb), 'LineWidth', 1.0);
grid on
xlabel('time (ms)');
ylabel('baseband');
legend('I', 'Q');
title('Analytic baseband LFM');
subplot(3, 1, 2);
plot(t_ms, lfm.tx_passband, 'LineWidth', 1.0); hold on
plot(t_ms, abs(lfm.tx_bb), 'LineWidth', 1.0);
grid on
xlabel('time (ms)');
ylabel('amplitude');
legend('real passband', 'envelope');
subplot(3, 1, 3);
S = fftshift(fft(lfm.tx_bb));
plot(lfm.f_rel_shifted + cfg.f_ref_hz, 20*log10(abs(S)/max(abs(S)) + eps), 'LineWidth', 1.0);
grid on
xlabel('frequency (Hz)');
ylabel('normalized magnitude (dB)');
title('LFM spectrum');
exportgraphics(gcf, file, 'Resolution', cfg.figure_resolution);
close(gcf);
end

function local_plot_representative_waveform(file, case_results, lfm, field_name, fig_title)
[kd, ks] = local_get_case_pair(case_results);
figure('Visible', 'off');
t_ms = lfm.t_s * 1e3;
subplot(2, 1, 1);
plot(t_ms, real(kd.(field_name)), 'LineWidth', 1.0); hold on
plot(t_ms, real(ks.(field_name)), '--', 'LineWidth', 1.0);
grid on
xlabel('time (ms)');
ylabel('real part');
legend('kdomain', 'kstat');
title(fig_title);
subplot(2, 1, 2);
plot(t_ms, abs(kd.(field_name)), 'LineWidth', 1.0); hold on
plot(t_ms, abs(ks.(field_name)), '--', 'LineWidth', 1.0);
grid on
xlabel('time (ms)');
ylabel('envelope');
legend('kdomain', 'kstat');
exportgraphics(gcf, file, 'Resolution', 160);
close(gcf);
end

function local_plot_representative_mf(file, case_results, lfm, field_name, fig_title)
[kd, ks] = local_get_case_pair(case_results);
figure('Visible', 'off');
t_ms = lfm.t_s * 1e3;
plot(t_ms, 20*log10(abs(kd.(field_name))/max(abs(kd.(field_name))) + eps), 'LineWidth', 1.0);
hold on
plot(t_ms, 20*log10(abs(ks.(field_name))/max(abs(ks.(field_name))) + eps), '--', 'LineWidth', 1.0);
grid on
xlabel('delay (ms)');
ylabel('normalized output (dB)');
legend('kdomain', 'kstat');
title(fig_title);
exportgraphics(gcf, file, 'Resolution', 160);
close(gcf);
end

function local_plot_representative_Hf(file, case_results)
[kd, ks] = local_get_case_pair(case_results);
figure('Visible', 'off');
subplot(2, 1, 1);
plot(kd.f_axis, 20*log10(abs(kd.H_f) + eps), '-o', 'LineWidth', 1.0); hold on
plot(ks.f_axis, 20*log10(abs(ks.H_f) + eps), '--s', 'LineWidth', 1.0);
grid on
xlabel('frequency (Hz)');
ylabel('|H_f| (dB)');
legend('kdomain', 'kstat');
title('Representative H_f magnitude');
subplot(2, 1, 2);
plot(kd.f_axis, unwrap(angle(kd.H_f)), '-o', 'LineWidth', 1.0); hold on
plot(ks.f_axis, unwrap(angle(ks.H_f)), '--s', 'LineWidth', 1.0);
grid on
xlabel('frequency (Hz)');
ylabel('phase (rad)');
legend('kdomain', 'kstat');
title('Representative H_f phase');
exportgraphics(gcf, file, 'Resolution', 160);
close(gcf);
end

function [kd, ks] = local_get_case_pair(case_results)
if isempty(case_results)
    error('Representative case was not recorded.');
end
kd_idx = find([case_results.branch_code] == 1, 1, 'first');
ks_idx = find([case_results.branch_code] == 2, 1, 'first');
if isempty(kd_idx) || isempty(ks_idx)
    error('Representative case must include both kdomain and kstat branches.');
end
kd = case_results(kd_idx);
ks = case_results(ks_idx);
end

function local_write_report(report_file, report, T, cfg, figure_files)
fid = fopen(report_file, 'w');
if fid < 0
    error('Cannot open report file %s for writing.', report_file);
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '# K-Stat vs K-Domain LFM Channel Validation\n\n');
fprintf(fid, 'Generated: %s\n\n', datestr(now, 31));
fprintf(fid, '## Setup\n\n');
fprintf(fid, '- Independent LFM channel validation; no PSK, BER/SER, noise injection, or communication script.\n');
fprintf(fid, '- Hs values: `%s` m; seed count: `%d`; grid `%d x %d`; Nf `%d`.\n', ...
    mat2str(cfg.Hs_list_m), cfg.seed_count, cfg.nx, cfg.ny, cfg.Nf);
fprintf(fid, '- LFM band: `%s` Hz; reference frequency: `%.3f` Hz; fs: `%.3f` Hz; duration: `%.5f` s.\n', ...
    mat2str(cfg.f_band_hz), cfg.f_ref_hz, cfg.lfm_fs_hz, cfg.lfm_duration_s);
fprintf(fid, '- Roughness mode: `%s`; PM shape wind: `%.3g m/s`.\n\n', ...
    cfg.roughness_scale_mode, cfg.wind_mps);

fprintf(fid, '## Key Results\n\n');
fprintf(fid, '- Max H invariant error, kdomain: `%.4g`.\n', report.max_kdomain_invariant_error);
fprintf(fid, '- Max H invariant error, kstat: `%.4g`.\n', report.max_kstat_invariant_error);
fprintf(fid, '- Max K-Stat phase-screen energy error: `%.4g`.\n', report.max_kstat_energy_error);
fprintf(fid, '- Max direct-path branch/seed delta: `%.4g`.\n', report.max_direct_path_abs_delta);
fprintf(fid, '- Max E[|h_ref|^2] relative error: `%.4g`.\n', report.max_E_abs_h_reflect2_rel_error);
fprintf(fid, '- Min total LFM envelope correlation: `%.4g`.\n', report.min_rx_total_env_corr);
fprintf(fid, '- Min reflected LFM envelope correlation: `%.4g`.\n', report.min_rx_reflect_env_corr);
fprintf(fid, '- Min total/reflected PDP correlation: `%.4g` / `%.4g`.\n', ...
    report.min_pdp_total_corr, report.min_pdp_reflect_corr);
fprintf(fid, '- Max total/reflected matched-filter peak relative error: `%.4g` / `%.4g`.\n', ...
    report.max_mf_total_peak_rel_error, report.max_mf_reflect_peak_rel_error);
fprintf(fid, '- Hard checks passed: `%d`.\n', report.all_hard_checks_passed);
fprintf(fid, '- Statistical targets passed: `%d`.\n\n', report.statistical_targets_passed);

fprintf(fid, '## Interpretation\n\n');
fprintf(fid, '- The transmitted test signal is an analytic baseband LFM envelope; figures also show its real passband projection.\n');
fprintf(fid, '- The received LFM is generated by multiplying the LFM spectrum by `H_f` or `H_reflect_f` interpolated onto the LFM FFT grid.\n');
fprintf(fid, '- This validates receiver-side waveform and matched-filter statistics. It does not alter the PE/WAPE source field and does not call the communication script.\n');
fprintf(fid, '- Single-seed waveforms are not expected to match pointwise; ensemble energy, envelope, PDP, and matched-filter trends are the comparison targets.\n\n');

fprintf(fid, '## Figures\n\n');
for ii = 1:numel(figure_files)
    fprintf(fid, '- `%s`\n', char(figure_files(ii)));
end

fprintf(fid, '\n## Summary Table\n\n');
preview = T(:, {'Hs_m', 'seed_count', 'E_abs_h_reflect2_rel_error', ...
    'rx_total_env_corr', 'rx_reflect_env_corr', 'pdp_total_corr', ...
    'pdp_reflect_corr', 'mf_total_peak_rel_error', ...
    'mf_reflect_peak_rel_error', 'max_kstat_energy_error'});
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
