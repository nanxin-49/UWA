run(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'bootstrap_project.m'));
%VALIDATE_KSTAT_JOINT_FREQUENCY_BOUNDARY_VERTICAL Boundary-only F=32 audit.

clear
format compact

cfg = struct();
cfg.U_mps = 5;
cfg.f_axis_hz = linspace(4000, 8000, 32).';
cfg.c0_mps = 1500;
cfg.reflect_coeff = -1;
cfg.grid_n = local_env_int('KSTAT_JOINT_GRID_N', 32);
cfg.aperture_m = local_env_scalar('KSTAT_JOINT_APERTURE_M', 50);
cfg.realization_count = local_env_int('KSTAT_JOINT_REALIZATIONS', 64);
cfg.seed_list = 31001 + (0:(cfg.realization_count-1));
cfg.joint_rng_seed = 41001;
cfg.independent_rng_seed = 51001;
cfg.lfm_fs_hz = 12000;
cfg.lfm_duration_s = 0.020;
cfg.lfm_fft_len = 512;

spec = raw_pm_spectrum_grid_vertical(cfg.U_mps, cfg.grid_n, cfg.grid_n, ...
    cfg.aperture_m, cfg.aperture_m);
F = numel(cfg.f_axis_hz);
M = cfg.realization_count;
alpha_f = 4*pi*cfg.f_axis_hz/cfg.c0_mps;
R_coh_f = cfg.reflect_coeff .* exp(-0.5 .* alpha_f.^2 .* ...
    spec.sigma_eta_discrete2_m2);

fprintf('Generating exact same-surface Kirchhoff ensemble.\n');
tic_explicit = tic;
explicit = local_empty_ensemble(spec.ny, spec.nx, F, M);
for mm = 1:M
    eta_xy = sample_raw_pm_surface_vertical(spec, cfg.seed_list(mm));
    delta_f = complex(zeros(spec.ny, spec.nx, F));
    for ii = 1:F
        G_xy = cfg.reflect_coeff .* exp(1i * alpha_f(ii) .* eta_xy);
        delta_f(:,:,ii) = G_xy - R_coh_f(ii);
    end
    explicit = local_accumulate_ensemble(explicit, delta_f, R_coh_f, mm);
end
explicit.generation_time_s = toc(tic_explicit);
explicit = local_finalize_ensemble(explicit, M);

fprintf('Generating independent-frequency K-stat ensemble.\n');
[delta_independent, independent_meta] = ...
    sample_kirchhoff_kstat_joint_frequency_vertical(spec, cfg.f_axis_hz, ...
    cfg.c0_mps, cfg.reflect_coeff, M, cfg.independent_rng_seed, 'independent');
independent = local_ensemble_from_array(delta_independent, R_coh_f);
clear delta_independent

fprintf('Generating covariance+pseudo-covariance joint-frequency K-stat ensemble.\n');
[delta_joint, joint_meta] = ...
    sample_kirchhoff_kstat_joint_frequency_vertical(spec, cfg.f_axis_hz, ...
    cfg.c0_mps, cfg.reflect_coeff, M, cfg.joint_rng_seed, 'joint');
joint = local_ensemble_from_array(delta_joint, R_coh_f);
clear delta_joint

lfm = local_lfm_setup(cfg);
explicit.waveform = local_waveform_stats(explicit.H_proxy_fm, cfg.f_axis_hz, lfm);
independent.waveform = local_waveform_stats(independent.H_proxy_fm, cfg.f_axis_hz, lfm);
joint.waveform = local_waveform_stats(joint.H_proxy_fm, cfg.f_axis_hz, lfm);

metrics = struct();
metrics.joint_cov_rel_error = local_rel_fro(joint.C_f, explicit.C_f);
metrics.independent_cov_rel_error = local_rel_fro(independent.C_f, explicit.C_f);
metrics.joint_pseudo_rel_error = local_rel_fro(joint.P_f, explicit.P_f);
metrics.independent_pseudo_rel_error = local_rel_fro(independent.P_f, explicit.P_f);
metrics.joint_adjacent_corr_rmse = local_rmse(joint.adjacent_corr, explicit.adjacent_corr);
metrics.independent_adjacent_corr_rmse = local_rmse( ...
    independent.adjacent_corr, explicit.adjacent_corr);
metrics.joint_pdp_corr = local_corr(joint.waveform.mean_pdp, explicit.waveform.mean_pdp);
metrics.independent_pdp_corr = local_corr( ...
    independent.waveform.mean_pdp, explicit.waveform.mean_pdp);
metrics.joint_matched_filter_corr = local_corr( ...
    joint.waveform.mean_matched_filter_power, ...
    explicit.waveform.mean_matched_filter_power);
metrics.independent_matched_filter_corr = local_corr( ...
    independent.waveform.mean_matched_filter_power, ...
    explicit.waveform.mean_matched_filter_power);
metrics.joint_incoherent_spectrum_rel_error = local_rel_fro( ...
    joint.mean_power_spectrum_kf, explicit.mean_power_spectrum_kf);
metrics.independent_incoherent_spectrum_rel_error = local_rel_fro( ...
    independent.mean_power_spectrum_kf, explicit.mean_power_spectrum_kf);
metrics.explicit_coherent_max_abs_error = max(abs(mean(explicit.H_proxy_fm,2)-R_coh_f));
metrics.joint_coherent_max_abs_error = max(abs(mean(joint.H_proxy_fm,2)-R_coh_f));
metrics.independent_coherent_max_abs_error = max(abs(mean(independent.H_proxy_fm,2)-R_coh_f));
metrics.explicit_pseudo_to_cov_ratio = norm(explicit.P_f,'fro') / max(norm(explicit.C_f,'fro'),eps);
metrics.joint_pseudo_to_cov_ratio = norm(joint.P_f,'fro') / max(norm(joint.C_f,'fro'),eps);
metrics.theory_pseudo_to_cov_ratio = joint_meta.pseudo_to_cov_fro_ratio_zero_lag;
metrics.joint_negative_to_positive_eigenvalue_ratio = ...
    joint_meta.negative_to_positive_eigenvalue_ratio;

figure_files = local_make_figures(cfg, explicit, independent, joint, metrics, R_coh_f);
result_file = project_result_file('validation', ...
    'validate_kstat_joint_frequency_boundary_vertical_result.mat');
save(result_file, 'cfg', 'spec', 'R_coh_f', 'explicit', 'independent', ...
    'joint', 'independent_meta', 'joint_meta', 'lfm', 'metrics', 'figure_files', '-v7.3');

disp(metrics)
fprintf('Saved %s\n', result_file);

function out = local_empty_ensemble(ny, nx, F, M)
out = struct();
out.power_spectrum_sum = zeros(ny,nx,F);
out.C_sum = complex(zeros(F,F));
out.P_sum = complex(zeros(F,F));
out.mean_delta_sum = complex(zeros(F,1));
out.H_proxy_fm = complex(zeros(F,M));
out.sample_count = 0;
end

function out = local_accumulate_ensemble(out, delta_f, R_coh_f, sample_index)
[ny,nx,F] = size(delta_f);
D = reshape(delta_f, ny*nx, F);
out.C_sum = out.C_sum + D.' * conj(D);
out.P_sum = out.P_sum + D.' * D;
out.mean_delta_sum = out.mean_delta_sum + sum(D,1).';
for ii = 1:F
    out.power_spectrum_sum(:,:,ii) = out.power_spectrum_sum(:,:,ii) + ...
        abs(fft2(delta_f(:,:,ii))).^2;
end
out.H_proxy_fm(:,sample_index) = R_coh_f + mean(D,1).';
out.sample_count = out.sample_count + ny*nx;
end

function out = local_finalize_ensemble(out, M)
out.C_f = out.C_sum / out.sample_count;
out.P_f = out.P_sum / out.sample_count;
out.mean_delta_f = out.mean_delta_sum / out.sample_count;
out.mean_power_spectrum_kf = out.power_spectrum_sum / M;
out.corr_f = local_cov_to_corr(out.C_f);
out.adjacent_corr = diag(out.corr_f,1);
out.pseudo_to_cov_ratio = norm(out.P_f,'fro') / max(norm(out.C_f,'fro'),eps);
out = rmfield(out, {'C_sum','P_sum','mean_delta_sum','power_spectrum_sum','sample_count'});
end

function out = local_ensemble_from_array(delta_fm, R_coh_f)
[ny,nx,~,M] = size(delta_fm);
out = local_empty_ensemble(ny,nx,numel(R_coh_f),M);
for mm = 1:M
    out = local_accumulate_ensemble(out, delta_fm(:,:,:,mm), R_coh_f, mm);
end
out.generation_time_s = NaN;
out = local_finalize_ensemble(out, M);
end

function corr_f = local_cov_to_corr(C_f)
d = sqrt(max(real(diag(C_f)),0));
corr_f = C_f ./ max(d*d.', eps);
end

function lfm = local_lfm_setup(cfg)
N = cfg.lfm_fft_len;
fs = cfg.lfm_fs_hz;
t = (0:N-1).' / fs;
active = t < cfg.lfm_duration_s;
T = cfg.lfm_duration_s;
B = cfg.f_axis_hz(end) - cfg.f_axis_hz(1);
tx = complex(zeros(N,1));
ta = t(active);
tx(active) = exp(1i*2*pi*((-B/2).*ta + 0.5*(B/T).*ta.^2));
tx = tx / max(norm(tx),eps);
lfm = struct('t_s',t,'tx',tx,'fs_hz',fs,'duration_s',T, ...
    'bandwidth_hz',B,'fft_len',N);
end

function stats = local_waveform_stats(H_fm, f_axis, lfm)
M = size(H_fm,2);
N = lfm.fft_len;
f_rel = f_axis(:) - mean(f_axis);
f_fft_shift = ((0:N-1).' - floor(N/2)) * (lfm.fs_hz/N);
tx_fft = fft(lfm.tx);
pdp = zeros(N,M);
mf_power = zeros(N,M);
for mm = 1:M
    H_shift = interp1(f_rel,H_fm(:,mm),f_fft_shift,'linear',0);
    h = ifft(ifftshift(H_shift));
    rx = ifft(tx_fft .* ifftshift(H_shift));
    mf = conv(rx,flipud(conj(lfm.tx)),'same');
    pdp(:,mm) = abs(h).^2;
    mf_power(:,mm) = abs(mf).^2;
end
stats = struct();
stats.mean_pdp = mean(pdp,2);
stats.mean_matched_filter_power = mean(mf_power,2);
stats.delay_axis_s = (0:N-1).' / lfm.fs_hz;
stats.matched_filter_axis_s = lfm.t_s;
end

function files = local_make_figures(cfg, kd, ind, joint, metrics, R_coh_f)
out_dir = project_result_dir('validation');
files = strings(0,1);

file = fullfile(out_dir,'kstat_joint_frequency_correlation_compare.png');
figure('Visible','off');
models = {kd,ind,joint}; names = {'same-surface kdomain','independent kstat','joint kstat'};
for ii=1:3
    subplot(1,3,ii); imagesc(cfg.f_axis_hz,cfg.f_axis_hz,abs(models{ii}.corr_f));
    axis xy; clim([0 1]); colorbar; title(names{ii}); xlabel('Hz'); ylabel('Hz');
end
exportgraphics(gcf,file,'Resolution',160); close(gcf); files(end+1)=string(file);

file = fullfile(out_dir,'kstat_joint_frequency_adjacent_correlation.png');
figure('Visible','off');
plot(cfg.f_axis_hz(1:end-1),abs(kd.adjacent_corr),'-o'); hold on
plot(cfg.f_axis_hz(1:end-1),abs(ind.adjacent_corr),'--s');
plot(cfg.f_axis_hz(1:end-1),abs(joint.adjacent_corr),'-.d'); grid on
xlabel('frequency (Hz)'); ylabel('|rho(f_i,f_{i+1})|'); legend(names,'Location','best');
exportgraphics(gcf,file,'Resolution',160); close(gcf); files(end+1)=string(file);

file = fullfile(out_dir,'kstat_joint_frequency_pseudo_covariance.png');
figure('Visible','off');
for ii=1:3
    subplot(1,3,ii); imagesc(cfg.f_axis_hz,cfg.f_axis_hz,abs(models{ii}.P_f));
    axis xy; colorbar; title(names{ii}); xlabel('Hz'); ylabel('Hz');
end
exportgraphics(gcf,file,'Resolution',160); close(gcf); files(end+1)=string(file);

file = fullfile(out_dir,'kstat_joint_frequency_pdp_compare.png');
figure('Visible','off');
plot(kd.waveform.delay_axis_s*1e3,kd.waveform.mean_pdp/max(kd.waveform.mean_pdp)); hold on
plot(ind.waveform.delay_axis_s*1e3,ind.waveform.mean_pdp/max(ind.waveform.mean_pdp),'--');
plot(joint.waveform.delay_axis_s*1e3,joint.waveform.mean_pdp/max(joint.waveform.mean_pdp),'-.');
xlim([0,20]); grid on; xlabel('delay (ms)'); ylabel('normalized mean PDP'); legend(names);
title(sprintf('PDP corr: independent %.3f, joint %.3f', ...
    metrics.independent_pdp_corr,metrics.joint_pdp_corr));
exportgraphics(gcf,file,'Resolution',160); close(gcf); files(end+1)=string(file);

file = fullfile(out_dir,'kstat_joint_frequency_matched_filter_compare.png');
figure('Visible','off');
plot(kd.waveform.matched_filter_axis_s*1e3,kd.waveform.mean_matched_filter_power/max(kd.waveform.mean_matched_filter_power)); hold on
plot(ind.waveform.matched_filter_axis_s*1e3,ind.waveform.mean_matched_filter_power/max(ind.waveform.mean_matched_filter_power),'--');
plot(joint.waveform.matched_filter_axis_s*1e3,joint.waveform.mean_matched_filter_power/max(joint.waveform.mean_matched_filter_power),'-.');
grid on; xlabel('time (ms)'); ylabel('normalized matched-filter power'); legend(names);
title(sprintf('MF corr: independent %.3f, joint %.3f', ...
    metrics.independent_matched_filter_corr,metrics.joint_matched_filter_corr));
exportgraphics(gcf,file,'Resolution',160); close(gcf); files(end+1)=string(file);

file = fullfile(out_dir,'kstat_joint_frequency_coherent_reflection.png');
figure('Visible','off');
plot(cfg.f_axis_hz,abs(R_coh_f),'-','LineWidth',1.4); hold on
plot(cfg.f_axis_hz,abs(mean(kd.H_proxy_fm,2)),'--');
plot(cfg.f_axis_hz,abs(mean(joint.H_proxy_fm,2)),':'); grid on
xlabel('frequency (Hz)'); ylabel('|coherent/aperture-mean reflection|');
legend('analytic R_{coh}','explicit ensemble','joint ensemble');
exportgraphics(gcf,file,'Resolution',160); close(gcf); files(end+1)=string(file);
end

function value = local_rel_fro(a,b)
value = norm(a-b,'fro')/max(norm(b,'fro'),eps);
end
function value = local_rmse(a,b)
value = sqrt(mean(abs(a(:)-b(:)).^2));
end
function value = local_corr(a,b)
a=a(:); b=b(:); a=a-mean(a); b=b-mean(b);
value = real((a'*b)/(max(norm(a)*norm(b),eps)));
end
function value = local_env_scalar(name,default_value)
raw=getenv(name); if isempty(raw), value=default_value; else, value=str2double(raw); end
if ~isfinite(value), error('%s must be finite.',name); end
end
function value = local_env_int(name,default_value)
value=round(local_env_scalar(name,default_value));
if value<1, error('%s must be positive.',name); end
end
