%VALIDATE_2K_OCEAN_SPECTRA Independent PM/JONSWAP/TMA 2k-phase benchmark.
% This script is deliberately separate from the project's PE, KStat, and
% SSA implementations. It uses omnidirectional radial spectra only.

clear
clc

script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(fileparts(script_dir));
setup_vertical_project();
result_dir = fullfile(project_root,'results','validation','ssa_2k_phase');
report_file = fullfile(project_root,'reports','validation_2k_ocean_spectra_report.md');
if ~isfolder(result_dir), mkdir(result_dir); end

%% Representative ocean spectra and grids
Hs = 0.20;                 % m
Tp = 4.0;                  % s
gamma_j = 3.3;             % JONSWAP peak enhancement
water_depth = 20;          % m, TMA
g = 9.81;                  % m/s^2
c = 1500;                  % m/s
f = linspace(0.03,2.0,2000).';
fp = 1/Tp;
f_acoustic = [4000 6000 8000];
theta_s_deg = linspace(0,80,801).';
theta_s = deg2rad(theta_s_deg);
target_m0 = (Hs/4)^2;

% PM-type and JONSWAP shapes are normalized only after their requested
% shape factors are formed. TMA uses the JONSWAP raw shape plus Phi.
S_pm_raw = f.^(-5).*exp(-1.25*(fp./f).^4);
sigma_j = 0.07*ones(size(f));
sigma_j(f > fp) = 0.09;
r_j = exp(-(f-fp).^2./(2*sigma_j.^2*fp^2));
S_jonswap_raw = S_pm_raw.*gamma_j.^r_j;
omega_h = 2*pi*f*sqrt(water_depth/g);
phi_tma = local_tma_factor(omega_h);
S_tma_raw = S_jonswap_raw.*phi_tma;

spectrum_names = {'PM-type','JONSWAP','TMA'};
branch_names = {'full','2k'};
S_f = {local_normalize_spectrum(S_pm_raw,f,target_m0), ...
    local_normalize_spectrum(S_jonswap_raw,f,target_m0), ...
    local_normalize_spectrum(S_tma_raw,f,target_m0)};

%% Map each frequency spectrum to a radial wavenumber spectrum
spectra = repmat(struct(),1,numel(spectrum_names));
mapping_relative_error = zeros(1,numel(spectrum_names));
covariance_tail_ratio = zeros(1,numel(spectrum_names));
for ss = 1:numel(spectrum_names)
    Sf = S_f{ss};
    if ss < 3
        K = (2*pi*f).^2/g;
    else
        K = local_tma_wavenumber(f,g,water_depth);
    end
    dK_df = gradient(K,f);
    SK = Sf./dK_df;
    m0_f = trapz(f,Sf);
    m0_K = trapz(K,SK);
    mapping_relative_error(ss) = abs(m0_K-m0_f)/m0_f;

    m2 = trapz(K,K.^2.*SK);
    sigma_eta = sqrt(m0_f);
    s_rms = sqrt(m2);
    L_eff = sqrt(m0_f/m2);
    rho = linspace(0,6*L_eff,6001).';
    C_eta = local_covariance_from_radial_spectrum(K,SK,rho);
    covariance_tail_ratio(ss) = C_eta(end)/C_eta(1);
    covariance_zero_error = abs(C_eta(1)-m0_f)/m0_f;

    spectra(ss).name = spectrum_names{ss};
    spectra(ss).K = K;
    spectra(ss).S_K = SK;
    spectra(ss).m0 = m0_f;
    spectra(ss).sigma_eta = sigma_eta;
    spectra(ss).s_rms = s_rms;
    spectra(ss).L_eff = L_eff;
    spectra(ss).rho = rho;
    spectra(ss).C_eta = C_eta;
    spectra(ss).mapping_error = mapping_relative_error(ss);
    spectra(ss).covariance_zero_error = covariance_zero_error;
end

surface_statistics = table(string(spectrum_names(:)), ...
    repmat(Hs,numel(spectrum_names),1),repmat(Tp,numel(spectrum_names),1), ...
    [spectra.sigma_eta].',[spectra.s_rms].',[spectra.L_eff].', ...
    [spectra.m0].',repmat(water_depth,numel(spectrum_names),1), ...
    'VariableNames',{'spectrum','Hs_m','Tp_s','sigma_eta_m','s_rms', ...
    'L_eff_m','m0','water_depth_m'});
writetable(surface_statistics,fullfile(result_dir,'ocean_surface_statistics.csv'))

%% Figure 1: radial wavenumber spectra
fig = figure('Visible','off','Color','w');
hold on
for ss = 1:numel(spectra)
    loglog(spectra(ss).K,spectra(ss).S_K,'LineWidth',1.4)
end
set(gca,'XScale','log','YScale','log')
all_sk = vertcat(spectra.S_K);
all_sk = all_sk(all_sk > 0);
peak_sk = max(all_sk);
ylim([peak_sk*1e-8,peak_sk*1.2])
grid on
xlabel('K (rad/m)')
ylabel('S_K(K)')
title('Omnidirectional radial ocean-wave wavenumber spectra')
legend(spectrum_names,'Location','best')
exportgraphics(fig,fullfile(result_dir,'ocean_wave_spectra.png'),'Resolution',180)
close(fig)

%% Acoustic SSA comparisons: 3 spectra x 3 acoustic frequencies
n_cases = numel(spectra)*numel(f_acoustic);
summary_values = nan(n_cases,24);
p_full_8khz = cell(numel(spectra),1);
p_2k_8khz = cell(numel(spectra),1);
case_row = 0;
largest_negative_relative = 0;
largest_negative_case = struct('spectrum','','frequency_hz',NaN,'branch','','value',NaN);
for ss = 1:numel(spectra)
    for ff = 1:numel(f_acoustic)
        case_row = case_row + 1;
        fa = f_acoustic(ff);
        k = 2*pi*fa/c;
        h = spectra(ss).sigma_eta;
        [p_full,p_2k,metrics] = local_ssa_case(theta_s,spectra(ss),k,h);
        if fa == 8000
            p_full_8khz{ss} = p_full;
            p_2k_8khz{ss} = p_2k;
        end

        kh = k*h;
        kL_eff = k*spectra(ss).L_eff;
        inside_kh = kh >= 0.25 && kh <= 2;
        inside_slope = spectra(ss).s_rms <= 0.2;
        summary_values(case_row,:) = [ss fa Hs Tp h spectra(ss).s_rms ...
            spectra(ss).L_eff spectra(ss).m0 kh kL_eff inside_kh inside_slope ...
            metrics.theta90_full metrics.theta95_full metrics.theta99_full ...
            metrics.theta90_2k metrics.theta95_2k metrics.theta99_2k ...
            metrics.TV_distance metrics.weighted_phase_rms ...
            metrics.F_phase_lt_0p1rad metrics.F_phase_lt_0p3rad ...
            metrics.max_negative_relative_full metrics.max_negative_relative_2k];

        neg_values = [metrics.max_negative_relative_full,metrics.max_negative_relative_2k];
        [case_negative,branch_idx] = max(neg_values);
        if case_negative > largest_negative_relative
            largest_negative_relative = case_negative;
            largest_negative_case.spectrum = spectrum_names{ss};
            largest_negative_case.frequency_hz = fa;
            largest_negative_case.branch = branch_names{branch_idx};
            largest_negative_case.value = case_negative;
        end
    end
end

summary_table = array2table(summary_values,'VariableNames',{ ...
    'spectrum_index','acoustic_frequency_hz','Hs_m','Tp_s','sigma_eta_m', ...
    's_rms','L_eff_m','m0','kh','kL_eff', ...
    'inside_previous_kh_range','inside_previous_slope_range', ...
    'theta90_full_deg','theta95_full_deg','theta99_full_deg', ...
    'theta90_2k_deg','theta95_2k_deg','theta99_2k_deg', ...
    'TV_distance','weighted_phase_rms_rad','F_phase_lt_0p1rad', ...
    'F_phase_lt_0p3rad','max_negative_relative_full', ...
    'max_negative_relative_2k'});
summary_table.spectrum = reshape(string(spectrum_names(summary_table.spectrum_index)),[],1);
summary_table = movevars(summary_table,'spectrum','Before','spectrum_index');
writetable(summary_table,fullfile(result_dir,'ocean_spectra_2k_summary.csv'))

%% Figure 2: 8 kHz angular-spectrum examples
fig = figure('Visible','off','Color','w');
layout = tiledlayout(fig,3,1,'TileSpacing','compact','Padding','compact');
for ss = 1:numel(spectra)
    nexttile
    plot(theta_s_deg,p_full_8khz{ss},'LineWidth',1.3)
    hold on
    plot(theta_s_deg,p_2k_8khz{ss},'--','LineWidth',1.3)
    grid on
    ylabel('p(theta) (rad^{-1})')
    title(sprintf('%s, 8 kHz',spectrum_names{ss}))
    if ss == 1
        legend('full phase: a = v_z','2k phase: a = 2k', ...
            'Location','best')
    end
end
xlabel(layout,'theta_s (deg)')
title(layout,'Normalized diffuse angular energy density')
exportgraphics(fig,fullfile(result_dir,'ocean_spectra_2k_examples_8khz.png'), ...
    'Resolution',180)
close(fig)

%% Figure 3: TV distance versus acoustic frequency
fig = figure('Visible','off','Color','w');
hold on
for ss = 1:numel(spectra)
    rows = summary_table.spectrum_index == ss;
    plot(summary_table.acoustic_frequency_hz(rows)/1000, ...
        summary_table.TV_distance(rows),'-o','LineWidth',1.3)
end
grid on
xlabel('acoustic frequency (kHz)')
ylabel('TV distance')
title('2k-phase angular-spectrum difference for ocean spectra')
legend(spectrum_names,'Location','best')
exportgraphics(fig,fullfile(result_dir,'ocean_spectra_2k_summary.png'), ...
    'Resolution',180)
close(fig)

%% Report and concise command-line summary
local_write_report(report_file, ...
    surface_statistics,summary_table,spectra,mapping_relative_error, ...
    covariance_tail_ratio,largest_negative_relative,largest_negative_case)

[worst_TV,worst_idx] = max(summary_table.TV_distance);
worst = summary_table(worst_idx,:);
largest_variance_error = max(mapping_relative_error);
largest_C0_error = max([spectra.covariance_zero_error]);

fprintf('\nOcean-spectrum 2k validation completed.\n\n')
fprintf('Surface statistics:\n')
for ss = 1:height(surface_statistics)
    fprintf('%s: sigma_eta = %.5g m, s_rms = %.5g, L_eff = %.5g m\n', ...
        surface_statistics.spectrum(ss),surface_statistics.sigma_eta_m(ss), ...
        surface_statistics.s_rms(ss),surface_statistics.L_eff_m(ss))
end
fprintf('\nWorst TV case:\n')
fprintf('spectrum = %s\n',worst.spectrum)
fprintf('frequency = %.3g kHz\n',worst.acoustic_frequency_hz/1000)
fprintf('kh = %.6g\n',worst.kh)
fprintf('s_rms = %.6g\n',worst.s_rms)
fprintf('theta95_full = %.6g deg\n',worst.theta95_full_deg)
fprintf('theta99_full = %.6g deg\n',worst.theta99_full_deg)
fprintf('TV distance = %.6g\n',worst_TV)
fprintf('weighted phase RMS = %.6g rad\n',worst.weighted_phase_rms_rad)
fprintf('F(<0.1 rad) = %.6g\n',worst.F_phase_lt_0p1rad)
fprintf('F(<0.3 rad) = %.6g\n',worst.F_phase_lt_0p3rad)
fprintf('\nLargest spectrum-mapping variance error = %.6g\n',largest_variance_error)
fprintf('Largest C(0) error = %.6g\n',largest_C0_error)
fprintf('Covariance tail ratios C(end)/C(0): PM-type = %.6g, JONSWAP = %.6g, TMA = %.6g\n', ...
    covariance_tail_ratio(1),covariance_tail_ratio(2),covariance_tail_ratio(3))
fprintf('Largest negative Hankel / positive peak = %.6g (%s, %.3g kHz, %s)\n', ...
    largest_negative_relative,largest_negative_case.spectrum, ...
    largest_negative_case.frequency_hz/1000,largest_negative_case.branch)

function S = local_normalize_spectrum(S_raw,f,target_m0)
S = S_raw*target_m0/trapz(f,S_raw);
end

function phi = local_tma_factor(omega_h)
phi = ones(size(omega_h));
low = omega_h <= 1;
mid = omega_h > 1 & omega_h < 2;
phi(low) = 0.5*omega_h(low).^2;
phi(mid) = 1-0.5*(2-omega_h(mid)).^2;
end

function K = local_tma_wavenumber(f,g,d)
omega = 2*pi*f;
K = zeros(size(f));
for ii = 1:numel(f)
    target = omega(ii)^2;
    fun = @(x) g*x*tanh(x*d)-target;
    K_deep = target/g;
    upper = max(1,2*K_deep+1);
    while fun(upper) < 0
        upper = 2*upper;
    end
    K(ii) = fzero(fun,[0 upper]);
end
end

function C_eta = local_covariance_from_radial_spectrum(K,S_K,rho)
C_eta = zeros(size(rho));
for ii = 1:numel(rho)
    C_eta(ii) = trapz(K,S_K.*besselj(0,K*rho(ii)));
end
end

function [p_full,p_2k,metrics] = local_ssa_case(theta_s,spectrum,k,h)
gamma_i = k;
gamma_s = k*cos(theta_s);
v_z = gamma_i+gamma_s;
q = k*sin(theta_s);
geometry = gamma_i^2*gamma_s.^2./v_z.^2;
rho = spectrum.rho;
C_eta = spectrum.C_eta;
F_full = zeros(size(theta_s));
F_2k = zeros(size(theta_s));

for jj = 1:numel(theta_s)
    J0 = besselj(0,q(jj)*rho);
    rough_full = exp(-(v_z(jj)*h)^2).*expm1(v_z(jj)^2*C_eta);
    rough_2k = exp(-(2*k*h)^2).*expm1((2*k)^2*C_eta);
    F_full(jj) = 2*pi*trapz(rho,rho.*J0.*rough_full);
    F_2k(jj) = 2*pi*trapz(rho,rho.*J0.*rough_2k);
end

[F_full,negative_full] = local_clip_hankel(F_full);
[F_2k,negative_2k] = local_clip_hankel(F_2k);
sigma_full = geometry.*F_full;
sigma_2k = geometry.*F_2k;
w_full = sigma_full.*sin(theta_s);
w_2k = sigma_2k.*sin(theta_s);
p_full = w_full/trapz(theta_s,w_full);
p_2k = w_2k/trapz(theta_s,w_2k);

metrics = local_distribution_metrics(theta_s,rad2deg(theta_s),p_full,p_2k);
sigma_delta_phi = k*h*(1-cos(theta_s));
rho_phi = exp(-0.5*sigma_delta_phi.^2); %#ok<NASGU>
metrics.weighted_phase_rms = sqrt(trapz(theta_s, ...
    p_full.*sigma_delta_phi.^2));
metrics.F_phase_lt_0p1rad = trapz(theta_s,p_full.*double(sigma_delta_phi < 0.1));
metrics.F_phase_lt_0p3rad = trapz(theta_s,p_full.*double(sigma_delta_phi < 0.3));
metrics.max_negative_relative_full = negative_full;
metrics.max_negative_relative_2k = negative_2k;
end

function [values,max_negative_relative] = local_clip_hankel(values)
scale = max(values);
negative_values = values(values < 0);
if isempty(negative_values)
    max_negative_relative = 0;
else
    max_negative_relative = max(abs(negative_values))/scale;
end
tiny = values < 0 & abs(values) < 1e-8*scale;
values(tiny) = 0;
if any(values < 0)
    warning(['Negative Hankel values clipped to zero; maximum ', ...
        'relative magnitude = %.3g.'],max_negative_relative)
    values(values < 0) = 0;
end
end

function metrics = local_distribution_metrics(theta,theta_deg,p_full,p_2k)
cdf_full = cumtrapz(theta,p_full);
cdf_2k = cumtrapz(theta,p_2k);
metrics.theta90_full = local_quantile(theta_deg,cdf_full,0.90);
metrics.theta95_full = local_quantile(theta_deg,cdf_full,0.95);
metrics.theta99_full = local_quantile(theta_deg,cdf_full,0.99);
metrics.theta90_2k = local_quantile(theta_deg,cdf_2k,0.90);
metrics.theta95_2k = local_quantile(theta_deg,cdf_2k,0.95);
metrics.theta99_2k = local_quantile(theta_deg,cdf_2k,0.99);
metrics.TV_distance = 0.5*trapz(theta,abs(p_2k-p_full));
end

function angle_deg = local_quantile(theta_deg,cdf_values,fraction)
idx = find(cdf_values >= fraction,1,'first');
if idx == 1
    angle_deg = theta_deg(1);
else
    angle_deg = interp1(cdf_values(idx-1:idx),theta_deg(idx-1:idx),fraction);
end
end

function local_write_report(file_name,surface_statistics,summary_table, ...
        spectra,mapping_error,covariance_tail_ratio,largest_negative_relative, ...
        largest_negative_case)
[~,worst_idx] = max(summary_table.TV_distance);
worst = summary_table(worst_idx,:);
largest_C0_error = max([spectra.covariance_zero_error]);
largest_kh = max(summary_table.kh);
smallest_L = min(surface_statistics.L_eff_m);
largest_slope = max(surface_statistics.s_rms);
count_kh = sum(summary_table.inside_previous_kh_range);
count_slope = sum(summary_table.inside_previous_slope_range);

fid = fopen(file_name,'w');
fprintf(fid,'# Validation of 2k phase with PM-type, JONSWAP, and TMA spectra\n\n');
fprintf(fid,['## 1. Goal\n\nExtend the Gaussian 2k-phase benchmark to PM-type, ', ...
    'JONSWAP, and TMA ocean-wave spectral shapes. This report studies only ', ...
    'the replacement of the roughness height-phase argument `(gamma_i + gamma_s) ', ...
    'eta` by `2k eta`. The lowest-order SSA calculation is a diagnostic, not an ', ...
    'exact solution. No PE, KStat, directional spreading, random realization, or ', ...
    'Monte Carlo calculation is used. The PM/JONSWAP/TMA spectra are treated as ', ...
    'omnidirectional radial spectra in this first comparison; directional spreading ', ...
    'is not included.\n\n']);

fprintf(fid,'## 2. Surface statistics\n\n');
fprintf(fid,'| spectrum | Hs (m) | Tp (s) | sigma_eta (m) | s_rms | L_eff (m) | m0 | water depth (m) |\n');
fprintf(fid,'|---|---:|---:|---:|---:|---:|---:|---:|\n');
for ii = 1:height(surface_statistics)
    fprintf(fid,'| %s | %.3f | %.2f | %.5f | %.5f | %.5f | %.6g | %.1f |\n', ...
        surface_statistics.spectrum(ii),surface_statistics.Hs_m(ii), ...
        surface_statistics.Tp_s(ii),surface_statistics.sigma_eta_m(ii), ...
        surface_statistics.s_rms(ii),surface_statistics.L_eff_m(ii), ...
        surface_statistics.m0(ii),surface_statistics.water_depth_m(ii));
end
fprintf(fid,['\nThe largest combined slope is %.5f for %s; the smallest effective ', ...
    'scale is %.5f m for %s. `L_eff` is only a convenient statistic ', ...
    'defined by sigma_eta/s_rms, not a unique correlation length.\n\n'], ...
    largest_slope,surface_statistics.spectrum(surface_statistics.s_rms == largest_slope), ...
    smallest_L,surface_statistics.spectrum(surface_statistics.L_eff_m == smallest_L));

fprintf(fid,'## 3. 2k validation results\n\n');
fprintf(fid,'| spectrum | frequency (kHz) | kh | s_rms | theta95 full | theta99 full | theta95 2k | TV | weighted phase RMS | F(<0.1 rad) | F(<0.3 rad) |\n');
fprintf(fid,'|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n');
for ii = 1:height(summary_table)
    fprintf(fid,'| %s | %.1f | %.4f | %.5f | %.3f | %.3f | %.3f | %.6f | %.6f | %.5f | %.5f |\n', ...
        summary_table.spectrum(ii),summary_table.acoustic_frequency_hz(ii)/1000, ...
        summary_table.kh(ii),summary_table.s_rms(ii), ...
        summary_table.theta95_full_deg(ii),summary_table.theta99_full_deg(ii), ...
        summary_table.theta95_2k_deg(ii),summary_table.TV_distance(ii), ...
        summary_table.weighted_phase_rms_rad(ii), ...
        summary_table.F_phase_lt_0p1rad(ii),summary_table.F_phase_lt_0p3rad(ii));
end

fprintf(fid,'\n## 4. Comparison with the Gaussian benchmark\n\n');
fprintf(fid,['The previous Gaussian 80-deg reference covered approximately `kh = ', ...
    '0.25--2` and `s_rms <= 0.2`. In the nine ocean-spectrum cases, %d/%d ', ...
    'are inside the kh interval and %d/%d are inside the slope interval. ', ...
    'These flags are positional references only, not validity labels. The ', ...
    'largest ocean kh is %.4f, and the largest ocean slope is %.5f.\n\n'], ...
    count_kh,height(summary_table),count_slope,height(summary_table), ...
    largest_kh,largest_slope);
fprintf(fid,['For PM/JONSWAP, the deep-water mapping is `K=(2*pi*f)^2/g`. For TMA, ', ...
    '`(2*pi*f)^2 = g K tanh(K d)` is solved independently at each frequency. ', ...
    'All three frequency spectra are normalized to `m0=(Hs/4)^2`; the radial ', ...
    'mapping and covariance checks are reported below.\n\n']);
fprintf(fid,['The frequency grid is 0.03--2.0 Hz with 2000 points; acoustic ', ...
    'frequencies are 4, 6, and 8 kHz; each covariance uses ', ...
    '`rho = linspace(0,6*L_eff,6001)` and the angular grid is 0--80 deg ', ...
    'with 801 points.\n\n']);
fprintf(fid,['The PM-type shape is `f^(-5)*exp(-1.25*(fp/f)^4)` and is explicitly ', ...
    'normalized to the prescribed Hs and Tp rather than tied to a wind-speed ', ...
    'convention. JONSWAP multiplies that shape by `gamma_j^r(f)` with ', ...
    '`gamma_j=3.3`, `sigma=0.07/0.09`. TMA multiplies the JONSWAP raw shape ', ...
    'by `Phi(omega_h)=0.5*omega_h^2` for `omega_h<=1`, ', ...
    '`1-0.5*(2-omega_h)^2` for `1<omega_h<2`, and 1 otherwise.\n\n']);
fprintf(fid,'Largest frequency-to-wavenumber variance error: %.3e\n\n',max(mapping_error));
fprintf(fid,'Largest covariance zero-lag relative error: %.3e\n\n',largest_C0_error);
fprintf(fid,'Covariance tail ratios `C_eta(end)/C_eta(1)`: PM-type %.3e, JONSWAP %.3e, TMA %.3e.\n\n', ...
    covariance_tail_ratio(1),covariance_tail_ratio(2),covariance_tail_ratio(3));
fprintf(fid,['Largest pre-clipping negative Hankel value relative to the positive ', ...
    'peak: %.3e (%s, %.1f kHz, %s branch).\n\n'],largest_negative_relative, ...
    largest_negative_case.spectrum,largest_negative_case.frequency_hz/1000, ...
    largest_negative_case.branch);
if largest_negative_relative > 1e-5
    fprintf(fid,['This exceeds 1e-5 and is retained as a finite-rho-window ', ...
        'numerical limitation; the negative tail is clipped for the normalized ', ...
        'diagnostic spectra and is not interpreted as an exact positive Hankel ', ...
        'transform.\n\n']);
end

fprintf(fid,'## 5. Figures\n\n');
fprintf(fid,'![Radial ocean-wave spectra](../results/validation/ssa_2k_phase/ocean_wave_spectra.png)\n\n');
fprintf(fid,'![8 kHz angular-spectrum examples](../results/validation/ssa_2k_phase/ocean_spectra_2k_examples_8khz.png)\n\n');
fprintf(fid,'![TV distance summary](../results/validation/ssa_2k_phase/ocean_spectra_2k_summary.png)\n\n');

fprintf(fid,'## 6. Short conclusion\n\n');
fprintf(fid,['Across the nine cases, full-phase 95%% and 99%% energy angles range ', ...
    'from %.2f--%.2f deg and %.2f--%.2f deg, respectively, over the 0--80 deg ', ...
    'normal-incidence scattering window. The maximum TV distance is %.6f for ', ...
    '%s at %.1f kHz (`kh=%.4f`).\n\n'], ...
    min(summary_table.theta95_full_deg),max(summary_table.theta95_full_deg), ...
    min(summary_table.theta99_full_deg),max(summary_table.theta99_full_deg), ...
    worst.TV_distance,worst.spectrum,worst.acoustic_frequency_hz/1000,worst.kh);
fprintf(fid,['Within this deliberately small comparison, the spectra do not produce ', ...
    'a dramatically broader full-phase tail than the previous Gaussian reference. ', ...
    'The 2k replacement changes the normalized angular spectra by TV distances ', ...
    'of %.6f--%.6f, with the largest energy-weighted phase RMS %.6f rad. ', ...
    'This is consistent with the earlier small-impact trend at the level of these ', ...
    'diagnostic cases, while the result is not a validation of real-ocean scattering ', ...
    'or a universal statement over wind speed, Hs, or Tp.\n\n'], ...
    min(summary_table.TV_distance),max(summary_table.TV_distance), ...
    max(summary_table.weighted_phase_rms_rad));
fprintf(fid,['The PM/JONSWAP/TMA shapes here are one normalized Hs/Tp condition and ', ...
    'an omnidirectional first comparison. Directional spreading, other sea states, ', ...
    'and exact higher-order scattering remain outside the scope.\n']);
fclose(fid);
end
