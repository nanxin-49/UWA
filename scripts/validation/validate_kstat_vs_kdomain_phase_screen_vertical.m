run(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'bootstrap_project.m'));
%VALIDATE_KSTAT_VS_KDOMAIN_PHASE_SCREEN_VERTICAL
% Surface-boundary-only validation of Kirchhoff K-Stat against explicit
% Kirchhoff phase-screen ensemble statistics. No PE/WAPE and no communication
% chain are used in this script.

clear
clc
format compact

cfg = struct();
cfg.f_list_hz = local_env_numeric_vector('KSTAT_KDOMAIN_F_LIST', [4000, 6000, 8000]);
cfg.Hs_list_m = local_env_numeric_vector('KSTAT_KDOMAIN_HS_LIST', [0.05, 0.1, 0.2, 0.5]);
cfg.M_list = local_env_int_vector('KSTAT_KDOMAIN_M_LIST', [8, 16, 32, 64]);
cfg.grid_n = local_env_int('KSTAT_KDOMAIN_GRID_N', 128);
cfg.nx = local_env_int('KSTAT_KDOMAIN_NX', cfg.grid_n);
cfg.ny = local_env_int('KSTAT_KDOMAIN_NY', cfg.grid_n);
cfg.xw_m = local_env_scalar('KSTAT_KDOMAIN_XW_M', 50);
cfg.yw_m = local_env_scalar('KSTAT_KDOMAIN_YW_M', 50);
cfg.U_mps = local_env_scalar('KSTAT_KDOMAIN_WIND_MPS', 5);
cfg.c0_mps = local_env_scalar('KSTAT_KDOMAIN_C0_MPS', 1500);
cfg.seed_base = local_env_int('KSTAT_KDOMAIN_SEED_BASE', 12345);
cfg.R0 = -1;
cfg.n_radial_bins = local_env_int('KSTAT_KDOMAIN_RADIAL_BINS', 48);
cfg.energy_tol = 1e-8;
cfg.figure_resolution = 160;

cfg.M_list = unique(max(1, round(cfg.M_list(:).')), 'stable');
cfg.M_max = max(cfg.M_list);
cfg.seed_list = cfg.seed_base + (0:(cfg.M_max - 1));

cfg.result_mat = project_result_file('validation', ...
    'validate_kstat_vs_kdomain_phase_screen_vertical_result.mat');
cfg.summary_csv = project_result_file('validation', ...
    'validate_kstat_vs_kdomain_phase_screen_vertical_summary.csv');
cfg.report_file = project_result_file('validation', ...
    'validate_kstat_vs_kdomain_phase_screen_vertical_report.md');
cfg.figure_dir = project_result_dir('validation');

fprintf('K-Stat vs K-Domain phase-screen validation\n');
fprintf('  f: %s Hz\n', mat2str(cfg.f_list_hz));
fprintf('  Hs: %s m\n', mat2str(cfg.Hs_list_m));
fprintf('  M: %s, grid: %d x %d, aperture: %.3g m x %.3g m\n', ...
    mat2str(cfg.M_list), cfg.nx, cfg.ny, cfg.xw_m, cfg.yw_m);
fprintf('  PM shape wind speed: %.3g m/s\n', cfg.U_mps);

[KX, KY, kx, ky, dkx, dky, dx, dy] = local_k_grid(cfg);
Kmag = sqrt(KX.^2 + KY.^2);
pm_spec = local_pm_spectrum(KX, KY, cfg.U_mps, dkx, dky);
psi_inc_xy = ones(cfg.ny, cfg.nx);
Psi_inc_k = fft2(psi_inc_xy);
P_inc = abs(Psi_inc_k).^2;

summary_rows = struct([]);
row_idx = 0;
case_results = struct([]);
case_idx = 0;
typical = struct();
typical_set = false;

for ih = 1:numel(cfg.Hs_list_m)
    Hs_m = cfg.Hs_list_m(ih);
    sigma_target = Hs_m / 4;
    W_eta = local_target_hs_w_eta(pm_spec.Phi2D, dkx, dky, sigma_target);
    eta_stack = zeros(cfg.ny, cfg.nx, cfg.M_max);
    G_stack = complex(zeros(cfg.ny, cfg.nx, cfg.M_max));
    Psi_ref_stack = complex(zeros(cfg.ny, cfg.nx, cfg.M_max));
    eta_sigma = zeros(cfg.M_max, 1);
    G_power_mean = zeros(cfg.M_max, 1);

    fprintf('Hs %.3g m: generating %d explicit phase screens.\n', Hs_m, cfg.M_max);
    for im = 1:cfg.M_max
        eta_xy = local_pm_surface_realization(pm_spec.Phi2D, dkx, dky, cfg.seed_list(im));
        Hs_raw = 4 * std(eta_xy(:));
        if Hs_raw > 0
            eta_xy = eta_xy * (Hs_m / Hs_raw);
        else
            eta_xy(:) = 0;
        end
        eta_stack(:, :, im) = eta_xy;
        eta_sigma(im) = std(eta_xy(:));
    end

    for ifq = 1:numel(cfg.f_list_hz)
        f_hz = cfg.f_list_hz(ifq);
        k0 = 2*pi*f_hz/cfg.c0_mps;
        propagating_mask = Kmag <= k0;
        kstat = local_kstat_reference(W_eta, Psi_inc_k, P_inc, KX, KY, ...
            dkx, dky, dx, dy, k0, cfg.R0);

        for im = 1:cfg.M_max
            G_xy = cfg.R0 .* exp(1i * 2 * k0 .* eta_stack(:, :, im));
            G_stack(:, :, im) = G_xy;
            Psi_ref_stack(:, :, im) = fft2(G_xy .* psi_inc_xy);
            G_power_mean(im) = mean(abs(G_xy(:)).^2);
        end

        for iM = 1:numel(cfg.M_list)
            M = cfg.M_list(iM);
            Psi_ref_sub = Psi_ref_stack(:, :, 1:M);
            Psi_mean_k = mean(Psi_ref_sub, 3);
            deltaPsi = Psi_ref_sub - Psi_mean_k;
            P_sca_kdomain = mean(abs(deltaPsi).^2, 3);
            Rcoh_kdomain = Psi_mean_k(1, 1) / max(Psi_inc_k(1, 1), eps);

            radial_kd = local_radial_average(Kmag, P_sca_kdomain, cfg.n_radial_bins);
            radial_ks = local_radial_average(Kmag, kstat.P_sca, cfg.n_radial_bins);
            [radial_err, radial_corr] = local_compare_radial(radial_kd.energy, radial_ks.energy);

            incoh_energy_kdomain = sum(P_sca_kdomain(:));
            incoh_energy_kstat = sum(kstat.P_sca(:));
            propagating_energy_kdomain = sum(P_sca_kdomain(propagating_mask));
            propagating_energy_kstat = sum(kstat.P_sca(propagating_mask));
            sigma_kdomain = mean(eta_sigma(1:M));
            energy_closure_kdomain = abs(mean(G_power_mean(1:M)) - 1);

            row_idx = row_idx + 1;
            summary_rows(row_idx).f_hz = f_hz; %#ok<SAGROW>
            summary_rows(row_idx).Hs_m = Hs_m;
            summary_rows(row_idx).M = M;
            summary_rows(row_idx).sigma_eta_kdomain_m = sigma_kdomain;
            summary_rows(row_idx).sigma_eta_kstat_m = kstat.sigma_eta_m;
            summary_rows(row_idx).sigma_eta_rel_error = ...
                abs(sigma_kdomain - kstat.sigma_eta_m) / max(kstat.sigma_eta_m, eps);
            summary_rows(row_idx).Rcoh_kdomain_real = real(Rcoh_kdomain);
            summary_rows(row_idx).Rcoh_kdomain_imag = imag(Rcoh_kdomain);
            summary_rows(row_idx).Rcoh_kdomain_abs = abs(Rcoh_kdomain);
            summary_rows(row_idx).Rcoh_kstat_real = real(kstat.R_coh);
            summary_rows(row_idx).Rcoh_kstat_imag = imag(kstat.R_coh);
            summary_rows(row_idx).Rcoh_kstat_abs = abs(kstat.R_coh);
            summary_rows(row_idx).Rcoh_abs_rel_error = ...
                abs(abs(Rcoh_kdomain) - abs(kstat.R_coh)) / max(abs(kstat.R_coh), eps);
            summary_rows(row_idx).Rcoh_complex_abs_error = abs(Rcoh_kdomain - kstat.R_coh);
            summary_rows(row_idx).incoh_energy_kdomain = incoh_energy_kdomain;
            summary_rows(row_idx).incoh_energy_kstat = incoh_energy_kstat;
            summary_rows(row_idx).incoh_energy_rel_error = ...
                abs(incoh_energy_kdomain - incoh_energy_kstat) / max(incoh_energy_kstat, eps);
            summary_rows(row_idx).propagating_energy_kdomain = propagating_energy_kdomain;
            summary_rows(row_idx).propagating_energy_kstat = propagating_energy_kstat;
            summary_rows(row_idx).propagating_energy_rel_error = ...
                abs(propagating_energy_kdomain - propagating_energy_kstat) / max(propagating_energy_kstat, eps);
            summary_rows(row_idx).radial_spectrum_error = radial_err;
            summary_rows(row_idx).radial_spectrum_corr = radial_corr;
            summary_rows(row_idx).energy_closure_error_kdomain = energy_closure_kdomain;
            summary_rows(row_idx).energy_closure_error_kstat = kstat.energy_closure_error;
            summary_rows(row_idx).phase_screen_incoherent_energy_kstat = ...
                kstat.phase_screen_incoherent_energy;
            summary_rows(row_idx).k0_rad_per_m = k0;
            summary_rows(row_idx).propagating_bin_fraction = nnz(propagating_mask) / numel(propagating_mask);

            if M == cfg.M_max
                case_idx = case_idx + 1;
                case_results(case_idx).f_hz = f_hz; %#ok<SAGROW>
                case_results(case_idx).Hs_m = Hs_m;
                case_results(case_idx).M = M;
                case_results(case_idx).radial_kdomain = radial_kd;
                case_results(case_idx).radial_kstat = radial_ks;
                case_results(case_idx).P_sca_kdomain = P_sca_kdomain;
                case_results(case_idx).P_sca_kstat = kstat.P_sca;
                case_results(case_idx).Psi_mean_k = Psi_mean_k;
            end

            if ~typical_set && M == cfg.M_max && abs(f_hz - 6000) == min(abs(cfg.f_list_hz - 6000)) && ...
                    abs(Hs_m - 0.2) == min(abs(cfg.Hs_list_m - 0.2))
                typical = struct('f_hz', f_hz, 'Hs_m', Hs_m, 'M', M, ...
                    'KX', KX, 'KY', KY, 'P_sca_kdomain', P_sca_kdomain, ...
                    'P_sca_kstat', kstat.P_sca, ...
                    'P_diff', P_sca_kdomain - kstat.P_sca, ...
                    'radial_kdomain', radial_kd, 'radial_kstat', radial_ks);
                typical_set = true;
            end
        end
    end
end

summary_table = struct2table(summary_rows);
validation_report = local_validation_report(summary_table, cfg);
figure_files = local_write_figures(summary_table, case_results, typical, cfg);
local_write_report(cfg.report_file, validation_report, summary_table, cfg, figure_files);

save(cfg.result_mat, 'cfg', 'summary_table', 'validation_report', ...
    'case_results', 'typical', 'figure_files');
writetable(summary_table, cfg.summary_csv);

disp(summary_table(:, {'f_hz', 'Hs_m', 'M', 'sigma_eta_kdomain_m', ...
    'sigma_eta_kstat_m', 'Rcoh_kdomain_abs', 'Rcoh_kstat_abs', ...
    'incoh_energy_kdomain', 'incoh_energy_kstat', ...
    'radial_spectrum_error', 'radial_spectrum_corr', ...
    'energy_closure_error_kdomain', 'energy_closure_error_kstat'}));
disp(validation_report);
fprintf('Saved %s\n', cfg.result_mat);
fprintf('Saved %s\n', cfg.summary_csv);
fprintf('Saved %s\n', cfg.report_file);

if ~validation_report.all_hard_checks_passed
    error('validate_kstat_vs_kdomain_phase_screen_vertical:FailedChecks', ...
        'One or more hard validation checks failed. See validation_report.');
end

function [KX, KY, kx, ky, dkx, dky, dx, dy] = local_k_grid(cfg)
dx = cfg.xw_m / cfg.nx;
dy = cfg.yw_m / cfg.ny;
kx = (2*pi/cfg.xw_m) * [0:(cfg.nx/2-1), -cfg.nx/2:-1];
ky = (2*pi/cfg.yw_m) * [0:(cfg.ny/2-1), -cfg.ny/2:-1];
[KX, KY] = meshgrid(kx, ky);
dkx = 2*pi / cfg.xw_m;
dky = 2*pi / cfg.yw_m;
end

function spec = local_pm_spectrum(KX, KY, U, dkx, dky)
g = 9.81;
alpha_PM = 8.10e-3;
beta_PM = 0.74;
K = sqrt(KX.^2 + KY.^2);
E1D = zeros(size(K));
Phi2D = zeros(size(K));
mask = K > 0;
K_nonzero = K(mask);
E1D(mask) = (alpha_PM ./ (2 .* K_nonzero.^3)) .* ...
    exp(-beta_PM * (g^2) ./ (U^4 .* K_nonzero.^2));
Phi2D(mask) = E1D(mask) ./ (2*pi*K_nonzero);
variance_discrete_m2 = sum(Phi2D(:)) * dkx * dky;
spec = struct('U', U, 'g', g, 'alpha_PM', alpha_PM, 'beta_PM', beta_PM, ...
    'K', K, 'E1D', E1D, 'Phi2D', Phi2D, 'mask', mask, ...
    'variance_discrete_m2', variance_discrete_m2);
end

function eta_xy = local_pm_surface_realization(Phi2D, dkx, dky, seed)
A = sqrt(Phi2D .* dkx .* dky);
rng(mod(round(seed), 2^32), 'twister');
N = (randn(size(Phi2D)) + 1i*randn(size(Phi2D))) / sqrt(2);
Zk = A .* N;
eta_xy = real(ifft2(Zk));
eta_xy = eta_xy * numel(eta_xy);
eta_xy = sqrt(2) * eta_xy;
end

function W_eta = local_target_hs_w_eta(Phi2D, dkx, dky, sigma_target)
spectral_variance = sum(Phi2D(:)) * dkx * dky;
if spectral_variance > 0 && sigma_target > 0
    scale_factor = sigma_target / sqrt(spectral_variance / (2*pi)^2);
else
    scale_factor = 0;
end
W_eta = Phi2D .* scale_factor.^2;
end

function kstat = local_kstat_reference(W_eta, Psi_inc_k, P_inc, KX, KY, ...
    dkx, dky, dx, dy, k0, R0)
grid_count = numel(W_eta);
fourier_area_scale = dkx * dky / (2*pi)^2;
C_eta_xy = real(ifft2(W_eta)) * grid_count * fourier_area_scale;
sigma_eta2 = max(real(C_eta_xy(1, 1)), 0);
sigma_eta = sqrt(sigma_eta2);
alpha = 2 * k0;
alpha2 = alpha^2;
G_mean = exp(-0.5 * alpha2 * sigma_eta2);
R_coh = R0 .* G_mean;
C_deltaG = exp(alpha2 .* (real(C_eta_xy) - sigma_eta2)) - exp(-alpha2 * sigma_eta2);
C_deltaG = real(C_deltaG);
S_deltaG_raw = real(fft2(C_deltaG) * dx * dy);
S_deltaG = max(S_deltaG_raw, 0);
phase_screen_incoherent_energy = sum(S_deltaG(:)) * fourier_area_scale;
energy_closure = abs(G_mean)^2 + phase_screen_incoherent_energy;
energy_closure_error = abs(energy_closure - 1);
conv_power = real(ifft2(fft2(S_deltaG) .* fft2(P_inc)));
P_sca = max(real(abs(R0)^2 .* conv_power .* fourier_area_scale), 0);
Kmag = sqrt(KX.^2 + KY.^2);
kstat = struct('sigma_eta_m', sigma_eta, 'sigma_eta2_m2', sigma_eta2, ...
    'G_mean', G_mean, 'R_coh', R_coh, 'C_eta_xy', C_eta_xy, ...
    'C_deltaG', C_deltaG, 'S_deltaG', S_deltaG, 'S_deltaG_raw', S_deltaG_raw, ...
    'P_sca', P_sca, 'energy_closure', energy_closure, ...
    'energy_closure_error', energy_closure_error, ...
    'phase_screen_incoherent_energy', phase_screen_incoherent_energy, ...
    'E_inc', sum(Psi_inc_k(:) .* conj(Psi_inc_k(:))), ...
    'Kmag', Kmag);
end

function radial = local_radial_average(Kmag, P, nbin)
max_k = max(Kmag(:));
edges = linspace(0, max_k, nbin + 1);
centers = 0.5 * (edges(1:end-1) + edges(2:end));
energy = zeros(1, nbin);
mean_power = zeros(1, nbin);
counts = zeros(1, nbin);
for ib = 1:nbin
    if ib == nbin
        mask = Kmag >= edges(ib) & Kmag <= edges(ib+1);
    else
        mask = Kmag >= edges(ib) & Kmag < edges(ib+1);
    end
    counts(ib) = nnz(mask);
    if counts(ib) > 0
        energy(ib) = sum(P(mask));
        mean_power(ib) = energy(ib) / counts(ib);
    else
        energy(ib) = 0;
        mean_power(ib) = NaN;
    end
end
radial = struct('edges', edges, 'centers', centers, 'energy', energy, ...
    'mean_power', mean_power, 'counts', counts, ...
    'energy_fraction', energy ./ max(sum(energy), eps));
end

function [err, corr_value] = local_compare_radial(a, b)
a = a(:);
b = b(:);
valid = isfinite(a) & isfinite(b) & (a >= 0) & (b >= 0);
a = a(valid);
b = b(valid);
if isempty(a) || sum(a) <= 0 || sum(b) <= 0
    err = NaN;
    corr_value = NaN;
    return
end
an = a ./ sum(a);
bn = b ./ sum(b);
err = norm(an - bn, 2) / max(norm(bn, 2), eps);
aa = an - mean(an);
bb = bn - mean(bn);
corr_value = sum(aa .* bb) / max(sqrt(sum(aa.^2) * sum(bb.^2)), eps);
end

function report = local_validation_report(T, cfg)
is_final_M = T.M == max(cfg.M_list);
report = struct();
report.script = mfilename;
report.created_at = char(datetime('now'));
report.f_list_hz = cfg.f_list_hz;
report.Hs_list_m = cfg.Hs_list_m;
report.M_list = cfg.M_list;
report.grid = [cfg.ny, cfg.nx];
report.wind_shape_mps = cfg.U_mps;
report.max_kstat_energy_closure_error = max(T.energy_closure_error_kstat);
report.max_kdomain_energy_closure_error = max(T.energy_closure_error_kdomain);
report.max_sigma_eta_rel_error = max(T.sigma_eta_rel_error);
coh_rel_threshold = 1e-3;
coh_rel_mask = is_final_M & T.Rcoh_kstat_abs > coh_rel_threshold;
if any(coh_rel_mask)
    report.final_M_mean_Rcoh_abs_rel_error = mean(T.Rcoh_abs_rel_error(coh_rel_mask), 'omitnan');
else
    report.final_M_mean_Rcoh_abs_rel_error = NaN;
end
report.final_M_mean_Rcoh_abs_error = mean(abs(T.Rcoh_kdomain_abs(is_final_M) - ...
    T.Rcoh_kstat_abs(is_final_M)), 'omitnan');
report.Rcoh_relative_error_threshold = coh_rel_threshold;
report.Rcoh_relative_error_note = sprintf(['Relative |Rcoh| error is averaged only ', ...
    'where kstat |Rcoh| > %.1e; use absolute error when coherent reflection is nearly zero.'], ...
    coh_rel_threshold);
report.final_M_mean_incoh_energy_rel_error = mean(T.incoh_energy_rel_error(is_final_M), 'omitnan');
report.final_M_mean_radial_spectrum_error = mean(T.radial_spectrum_error(is_final_M), 'omitnan');
report.final_M_mean_radial_spectrum_corr = mean(T.radial_spectrum_corr(is_final_M), 'omitnan');
report.kstat_energy_passed = report.max_kstat_energy_closure_error <= cfg.energy_tol;
report.kdomain_energy_passed = report.max_kdomain_energy_closure_error <= 1e-12;
report.finite_summary_passed = all(all(isfinite(T{:, :})));
report.all_hard_checks_passed = report.kstat_energy_passed && ...
    report.kdomain_energy_passed && report.finite_summary_passed;
report.interpretation = ['This surface-only validation compares explicit ', ...
    'Kirchhoff phase-screen ensemble statistics against the Kirchhoff K-Stat ', ...
    'statistical phase-screen formulas. Single-seed complex fields are not ', ...
    'expected to match pointwise.'];
report.pe_note = ['If PE/WAPE propagation is added later, compare receiver-side ', ...
    'statistics such as E[|h_ref|^2], std(|h_ref|), or average PDP, not ', ...
    'single-realization h equality.'];
end

function figure_files = local_write_figures(T, case_results, typical, cfg)
figure_files = strings(0, 1);
figure_files(end+1, 1) = local_plot_Rcoh(T, cfg);
figure_files(end+1, 1) = local_plot_incoh_energy(T, cfg);
figure_files(end+1, 1) = local_plot_convergence(T, cfg);
figure_files(end+1, 1) = local_plot_radial(case_results, typical, cfg);
figure_files(end+1, 1) = local_plot_2d_spectrum(typical, cfg);
end

function file = local_plot_Rcoh(T, cfg)
file = fullfile(cfg.figure_dir, 'kstat_kdomain_Rcoh_vs_Hs_f.png');
T = T(T.M == max(cfg.M_list), :);
fig = figure('Visible', 'off', 'Color', 'w');
tiledlayout(1, numel(cfg.f_list_hz), 'Padding', 'compact', 'TileSpacing', 'compact');
for ifq = 1:numel(cfg.f_list_hz)
    nexttile;
    S = T(T.f_hz == cfg.f_list_hz(ifq), :);
    semilogy(S.Hs_m, max(S.Rcoh_kdomain_abs, eps), '-o', 'LineWidth', 1.2); hold on;
    semilogy(S.Hs_m, max(S.Rcoh_kstat_abs, eps), '--s', 'LineWidth', 1.2);
    grid on; xlabel('H_s (m)'); ylabel('|R_{coh}|');
    title(sprintf('%.1f kHz', cfg.f_list_hz(ifq)/1000));
    legend('kdomain ensemble', 'kstat', 'Location', 'southwest');
end
sgtitle(sprintf('Coherent reflection comparison, M=%d', max(cfg.M_list)));
exportgraphics(fig, file, 'Resolution', cfg.figure_resolution);
close(fig);
end

function file = local_plot_incoh_energy(T, cfg)
file = fullfile(cfg.figure_dir, 'kstat_kdomain_incoh_energy_vs_Hs_f.png');
T = T(T.M == max(cfg.M_list), :);
fig = figure('Visible', 'off', 'Color', 'w');
tiledlayout(1, numel(cfg.f_list_hz), 'Padding', 'compact', 'TileSpacing', 'compact');
for ifq = 1:numel(cfg.f_list_hz)
    nexttile;
    S = T(T.f_hz == cfg.f_list_hz(ifq), :);
    semilogy(S.Hs_m, max(S.incoh_energy_kdomain, eps), '-o', 'LineWidth', 1.2); hold on;
    semilogy(S.Hs_m, max(S.incoh_energy_kstat, eps), '--s', 'LineWidth', 1.2);
    grid on; xlabel('H_s (m)'); ylabel('sum P_{sca}(K)');
    title(sprintf('%.1f kHz', cfg.f_list_hz(ifq)/1000));
    legend('kdomain ensemble', 'kstat', 'Location', 'northwest');
end
sgtitle(sprintf('Incoherent spectral energy comparison, M=%d', max(cfg.M_list)));
exportgraphics(fig, file, 'Resolution', cfg.figure_resolution);
close(fig);
end

function file = local_plot_convergence(T, cfg)
file = fullfile(cfg.figure_dir, 'kstat_kdomain_convergence_vs_M.png');
fig = figure('Visible', 'off', 'Color', 'w');
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
nexttile;
for ih = 1:numel(cfg.Hs_list_m)
    S = T(T.f_hz == cfg.f_list_hz(min(2, numel(cfg.f_list_hz))) & T.Hs_m == cfg.Hs_list_m(ih), :);
    plot(S.M, S.Rcoh_abs_rel_error, '-o', 'LineWidth', 1.1); hold on;
end
grid on; xlabel('M'); ylabel('|R_{coh}| relative error');
legend(compose('H_s=%.3g m', cfg.Hs_list_m), 'Location', 'best');
title('Coherent convergence');
nexttile;
for ih = 1:numel(cfg.Hs_list_m)
    S = T(T.f_hz == cfg.f_list_hz(min(2, numel(cfg.f_list_hz))) & T.Hs_m == cfg.Hs_list_m(ih), :);
    plot(S.M, S.radial_spectrum_error, '-s', 'LineWidth', 1.1); hold on;
end
grid on; xlabel('M'); ylabel('normalized radial L2 error');
title('Radial spectrum convergence');
sgtitle(sprintf('Convergence at %.1f kHz', cfg.f_list_hz(min(2, numel(cfg.f_list_hz)))/1000));
exportgraphics(fig, file, 'Resolution', cfg.figure_resolution);
close(fig);
end

function file = local_plot_radial(case_results, typical, cfg)
file = fullfile(cfg.figure_dir, 'kstat_kdomain_radial_spectrum_compare.png');
if isempty(case_results)
    return
end
fig = figure('Visible', 'off', 'Color', 'w');
plot_count = min(numel(case_results), 6);
tiledlayout(2, ceil(plot_count/2), 'Padding', 'compact', 'TileSpacing', 'compact');
selected = round(linspace(1, numel(case_results), plot_count));
for ii = 1:plot_count
    C = case_results(selected(ii));
    nexttile;
    kd = C.radial_kdomain.energy_fraction;
    ks = C.radial_kstat.energy_fraction;
    semilogy(C.radial_kdomain.centers, max(kd, eps), '-o', 'LineWidth', 1.1); hold on;
    semilogy(C.radial_kstat.centers, max(ks, eps), '--s', 'LineWidth', 1.1);
    grid on; xlabel('K_h (rad/m)'); ylabel('radial energy fraction');
    title(sprintf('%.1f kHz, H_s=%.2g m', C.f_hz/1000, C.Hs_m));
end
if isfield(typical, 'f_hz')
    sgtitle(sprintf('Radial P_{sca} comparison, typical %.1f kHz H_s=%.2g m', ...
        typical.f_hz/1000, typical.Hs_m));
else
    sgtitle('Radial P_{sca} comparison');
end
legend('kdomain ensemble', 'kstat', 'Location', 'best');
exportgraphics(fig, file, 'Resolution', cfg.figure_resolution);
close(fig);
end

function file = local_plot_2d_spectrum(typical, cfg)
file = fullfile(cfg.figure_dir, 'kstat_kdomain_2d_spectrum_typical.png');
if ~isfield(typical, 'P_sca_kdomain')
    return
end
kd = fftshift(log10(max(typical.P_sca_kdomain, realmin)));
ks = fftshift(log10(max(typical.P_sca_kstat, realmin)));
df = fftshift(typical.P_diff);
lim = max(abs(df(:)));
fig = figure('Visible', 'off', 'Color', 'w');
tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
nexttile;
imagesc(kd); axis image; colorbar; title('log10 kdomain P_{sca}');
nexttile;
imagesc(ks); axis image; colorbar; title('log10 kstat P_{sca}');
nexttile;
imagesc(df); axis image; colorbar; title('kdomain - kstat');
if lim > 0
    clim([-lim, lim]);
end
sgtitle(sprintf('Typical 2D spectrum, %.1f kHz, H_s=%.2g m, M=%d', ...
    typical.f_hz/1000, typical.Hs_m, typical.M));
exportgraphics(fig, file, 'Resolution', cfg.figure_resolution);
close(fig);
end

function local_write_report(report_file, report, T, cfg, figure_files)
fid = fopen(report_file, 'w');
if fid < 0
    error('Could not open report file: %s', report_file);
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '# K-Stat vs K-Domain Phase-Screen Validation\n\n');
fprintf(fid, 'Generated: %s\n\n', report.created_at);
fprintf(fid, '## Setup\n\n');
fprintf(fid, '- Surface-boundary-only validation; no PE/WAPE and no communication chain.\n');
fprintf(fid, '- Frequencies: `%s` Hz.\n', mat2str(cfg.f_list_hz));
fprintf(fid, '- Hs values: `%s` m.\n', mat2str(cfg.Hs_list_m));
fprintf(fid, '- M values: `%s`; grid `%d x %d`; PM shape wind `%.3g m/s`.\n\n', ...
    mat2str(cfg.M_list), cfg.ny, cfg.nx, cfg.U_mps);
fprintf(fid, '## Key Results\n\n');
fprintf(fid, '- Max K-Stat energy closure error: `%.4g`.\n', report.max_kstat_energy_closure_error);
fprintf(fid, '- Max K-Domain phase-screen power closure error: `%.4g`.\n', report.max_kdomain_energy_closure_error);
fprintf(fid, '- Max sigma_eta relative error: `%.4g`.\n', report.max_sigma_eta_rel_error);
fprintf(fid, '- Final-M mean |Rcoh| relative error where kstat |Rcoh| > %.1e: `%.4g`.\n', ...
    report.Rcoh_relative_error_threshold, report.final_M_mean_Rcoh_abs_rel_error);
fprintf(fid, '- Final-M mean |Rcoh| absolute error: `%.4g`.\n', report.final_M_mean_Rcoh_abs_error);
fprintf(fid, '- Final-M mean radial spectrum L2 error: `%.4g`.\n', report.final_M_mean_radial_spectrum_error);
fprintf(fid, '- Final-M mean radial spectrum correlation: `%.4g`.\n\n', report.final_M_mean_radial_spectrum_corr);
fprintf(fid, '## Interpretation\n\n');
fprintf(fid, '- kdomain single realizations are concrete phase screens; kstat is a statistical phase-screen model.\n');
fprintf(fid, '- Single-seed complex fields are not expected to match pointwise.\n');
fprintf(fid, '- Valid comparisons are ensemble coherent mean, incoherent power spectrum, radial spectral shape, energy closure, and statistical trends.\n');
fprintf(fid, '- If PE/WAPE propagation is added later, compare receiver-side statistics such as `E[|h_ref|^2]`, `std(|h_ref|)`, or average PDP, not single-realization `h` equality.\n\n');
fprintf(fid, '## Figures\n\n');
for ii = 1:numel(figure_files)
    fprintf(fid, '- `%s`\n', figure_files(ii));
end
fprintf(fid, '\n## Summary Table Preview\n\n');
preview = T(T.M == max(cfg.M_list), :);
preview = preview(:, {'f_hz', 'Hs_m', 'M', 'Rcoh_kdomain_abs', ...
    'Rcoh_kstat_abs', 'incoh_energy_kdomain', 'incoh_energy_kstat', ...
    'radial_spectrum_error', 'energy_closure_error_kstat'});
fprintf(fid, '%s\n', evalc('disp(preview)'));
end

function values = local_env_numeric_vector(name, default_values)
raw = strtrim(getenv(name));
if isempty(raw)
    values = default_values;
    return
end
raw = strrep(raw, '[', ' ');
raw = strrep(raw, ']', ' ');
raw = strrep(raw, ',', ' ');
values = sscanf(raw, '%f').';
if isempty(values)
    values = default_values;
end
end

function values = local_env_int_vector(name, default_values)
values = round(local_env_numeric_vector(name, default_values));
values = values(isfinite(values) & values > 0);
if isempty(values)
    values = default_values;
end
end

function value = local_env_scalar(name, default_value)
raw = strtrim(getenv(name));
if isempty(raw)
    value = default_value;
    return
end
value = str2double(raw);
if ~isfinite(value)
    value = default_value;
end
end

function value = local_env_int(name, default_value)
value = round(local_env_scalar(name, default_value));
if ~isfinite(value) || value <= 0
    value = default_value;
end
end
