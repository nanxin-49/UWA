%VALIDATE_2K_PHASE_APPROX Independent check of the near-normal 2k phase approximation.
% Angles are measured from the surface normal. This script does not call
% the project PE, KStat, or SSA implementations.

clear
clc

script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(fileparts(script_dir));
setup_vertical_project();
result_dir = fullfile(project_root,'results','validation','ssa_2k_phase');
report_file = fullfile(project_root,'reports','validation_2k_phase_report.md');
if ~isfolder(result_dir), mkdir(result_dir); end

%% Part A: analytic geometry error
k = 1;
theta_i_deg = 0:0.5:30;
theta_s_deg = 0:0.5:60;
[theta_s_grid_deg, theta_i_grid_deg] = meshgrid(theta_s_deg, theta_i_deg);
epsilon_alpha = 1 - (cosd(theta_i_grid_deg) + cosd(theta_s_grid_deg))/2;

fig = figure('Visible','off','Color','w');
imagesc(theta_s_deg,theta_i_deg,100*epsilon_alpha)
axis xy
xlabel('theta_s (deg)')
ylabel('theta_i (deg)')
title({'2k phase-coefficient relative error', ...
    'Angles measured from surface normal'})
cb = colorbar;
cb.Label.String = 'epsilon_alpha (%)';
hold on
[contours,handle] = contour(theta_s_deg,theta_i_deg,100*epsilon_alpha, ...
    [1 3 5],'k','LineWidth',1);
clabel(contours,handle,'Color','k')
exportgraphics(fig,fullfile(result_dir,'phase_coefficient_error_map.png'), ...
    'Resolution',180)
close(fig)

theta_i_check_deg = [0 5 10 15 20].';
thresholds = [0.01 0.03 0.05];
theta_s_limits = nan(numel(theta_i_check_deg),numel(thresholds));
for ii = 1:numel(theta_i_check_deg)
    error_line = 1 - (cosd(theta_i_check_deg(ii)) + cosd(theta_s_deg))/2;
    for jj = 1:numel(thresholds)
        allowed = theta_s_deg(error_line <= thresholds(jj));
        if ~isempty(allowed)
            theta_s_limits(ii,jj) = max(allowed);
        end
    end
end

geometry_table = table(theta_i_check_deg,theta_s_limits(:,1), ...
    theta_s_limits(:,2),theta_s_limits(:,3), ...
    'VariableNames',{'theta_i_deg','theta_s_max_1pct_deg', ...
    'theta_s_max_3pct_deg','theta_s_max_5pct_deg'});
writetable(geometry_table,fullfile(result_dir,'phase_geometry_limits.csv'))

%% Parts B and C: Gaussian roughness and lowest-order two-dimensional SSA
theta_s_deg = linspace(0,80,801).';
theta_s = deg2rad(theta_s_deg);
idx_60 = theta_s_deg <= 60;
kh_values = [0.25 0.5 1.0 2.0];
kl_values = [10 20 40];
max_s_rms = 0.25;

n_valid = sum(reshape(2*kh_values(:)./kl_values,[],1) <= max_s_rms);
values = nan(n_valid,16);
comparison_values = nan(n_valid,17);
p_full_all = cell(n_valid,1);
p_2k_all = cell(n_valid,1);
row = 0;

for kh = kh_values
    for kl = kl_values
        s_rms = 2*kh/kl;
        if s_rms > max_s_rms
            fprintf('Skipping kh = %.2g, kl = %.2g because s_rms = %.3g > %.2f.\n', ...
                kh,kl,s_rms,max_s_rms)
            continue
        end

        row = row + 1;
        [p_full,p_2k,w_full,w_2k,negative_rel] = ...
            local_ssa_spectra(theta_s,kh,kl,k);
        [theta_full,theta_2k,TV_distance] = local_distribution_metrics( ...
            theta_s,theta_s_deg,p_full,p_2k);

        theta_60 = theta_s(idx_60);
        theta_60_deg = theta_s_deg(idx_60);
        p_full_60 = w_full(idx_60)/trapz(theta_60,w_full(idx_60));
        p_2k_60 = w_2k(idx_60)/trapz(theta_60,w_2k(idx_60));
        [theta_full_60,theta_2k_60,TV_60] = local_distribution_metrics( ...
            theta_60,theta_60_deg,p_full_60,p_2k_60);

        sigma_delta_phi = kh*(1-cos(theta_s));
        rho_phi = exp(-0.5*sigma_delta_phi.^2);
        weighted_phase_rms = sqrt(trapz(theta_s, ...
            p_full.*sigma_delta_phi.^2));
        phase_fractions = arrayfun(@(delta) trapz(theta_s, ...
            p_full.*double(sigma_delta_phi < delta)),[0.1 0.3]);

        values(row,:) = [kh kl s_rms theta_full theta_2k ...
            theta_2k(2)-theta_full(2) TV_distance weighted_phase_rms ...
            phase_fractions negative_rel];
        comparison_values(row,:) = [kh kl ...
            theta_full_60(2) theta_full(2) theta_full(2)-theta_full_60(2) ...
            theta_full_60(3) theta_full(3) theta_full(3)-theta_full_60(3) ...
            theta_2k_60(2) theta_2k(2) theta_2k(2)-theta_2k_60(2) ...
            theta_2k_60(3) theta_2k(3) theta_2k(3)-theta_2k_60(3) ...
            TV_60 TV_distance TV_distance-TV_60];
        p_full_all{row} = p_full;
        p_2k_all{row} = p_2k;
    end
end

summary_table = array2table(values,'VariableNames',{ ...
    'kh','kl','s_rms', ...
    'theta90_full_deg','theta95_full_deg','theta99_full_deg', ...
    'theta90_2k_deg','theta95_2k_deg','theta99_2k_deg', ...
    'theta95_difference_deg','TV_distance','weighted_phase_rms_rad', ...
    'F_phase_lt_0p1rad','F_phase_lt_0p3rad', ...
    'max_negative_relative_full','max_negative_relative_2k'});
writetable(summary_table,fullfile(result_dir,'ssa_2k_summary.csv'))

comparison_table = array2table(comparison_values,'VariableNames',{ ...
    'kh','kl', ...
    'theta95_full_60deg','theta95_full_80deg','delta_theta95_full_deg', ...
    'theta99_full_60deg','theta99_full_80deg','delta_theta99_full_deg', ...
    'theta95_2k_60deg','theta95_2k_80deg','delta_theta95_2k_deg', ...
    'theta99_2k_60deg','theta99_2k_80deg','delta_theta99_2k_deg', ...
    'TV_60deg','TV_80deg','delta_TV'});
writetable(comparison_table, ...
    fullfile(result_dir,'ssa_angle_limit_comparison.csv'))

%% Figure 2: representative normalized diffuse angular spectra
targets = [0.5 10; 1.0 20; 2.0 40];
fig = figure('Visible','off','Color','w');
layout = tiledlayout(fig,3,1,'TileSpacing','compact','Padding','compact');
for ii = 1:size(targets,1)
    distance = (summary_table.kh-targets(ii,1)).^2 + ...
        ((summary_table.kl-targets(ii,2))/10).^2;
    [~,idx] = min(distance);
    nexttile
    plot(theta_s_deg,p_full_all{idx},'LineWidth',1.3)
    hold on
    plot(theta_s_deg,p_2k_all{idx},'--','LineWidth',1.3)
    grid on
    ylabel('p(theta) (rad^{-1})')
    title(sprintf('kh = %.2g, kl = %.2g, s_{rms} = %.3g', ...
        summary_table.kh(idx),summary_table.kl(idx),summary_table.s_rms(idx)))
    if ii == 1
        legend('full phase: a = gamma_i + gamma_s', ...
            '2k phase approximation','Location','best')
    end
end
xlabel(layout,'theta_s (deg)')
title(layout,'Normalized diffuse angular energy density')
exportgraphics(fig,fullfile(result_dir,'ssa_spectrum_examples.png'), ...
    'Resolution',180)
close(fig)

%% Figure 3: summary of the normalized-spectrum difference
fig = figure('Visible','off','Color','w','Position',[100 100 850 600]);
scatter(summary_table.kl,summary_table.kh,90,summary_table.TV_distance,'filled')
xlabel('kl')
ylabel('kh')
title('2k-phase angular-spectrum difference')
grid on
xlim([8 45])
ylim([0.15 2.1])
xticks(kl_values)
cb = colorbar;
cb.Label.String = 'TV distance';
for ii = 1:height(summary_table)
    text(summary_table.kl(ii),summary_table.kh(ii), ...
        sprintf('  %.1f^\\circ',summary_table.theta95_full_deg(ii)), ...
        'FontSize',8,'Interpreter','tex')
end
exportgraphics(fig,fullfile(result_dir,'ssa_2k_error_summary.png'), ...
    'Resolution',180)
close(fig)

%% Sanity checks
epsilon_normal = 1-(cos(0)+cos(0))/2;
sigma_delta_phi_normal = 1*(2-cos(0)-cos(0));
fprintf('\nSanity check A at theta_i = theta_s = 0 deg:\n')
fprintf('epsilon_alpha = %.3g\n',epsilon_normal)
fprintf('sigma_delta_phi = %.3g rad\n',sigma_delta_phi_normal)

[p_full_test,p_2k_test,~,~,negative_rel_test] = ...
    local_ssa_spectra(theta_s,0.01,20,k);
TV_test = 0.5*trapz(theta_s,abs(p_2k_test-p_full_test));
fprintf('\nSanity check B:\n')
fprintf('TV distance for kh = 0.01, kl = 20: %.6g\n',TV_test)

%% Markdown report and command-line summary
local_write_report(report_file, ...
    geometry_table,summary_table,comparison_table,negative_rel_test,max_s_rms)

[~,worst_idx] = max(summary_table.TV_distance);
worst = summary_table(worst_idx,:);
worst_comparison = comparison_table(worst_idx,:);
max_negative_relative = max([summary_table.max_negative_relative_full; ...
    summary_table.max_negative_relative_2k]);
fprintf('\n2k phase approximation validation completed.\n\n')
fprintf('Analytic table:\nphase_geometry_limits.csv\n\n')
fprintf('SSA summary:\nssa_2k_summary.csv\n\n')
fprintf('Angle-limit comparison:\nssa_angle_limit_comparison.csv\n\n')
fprintf('Report:\n%s\n\n',report_file)
fprintf('Worst tested case by 80-deg TV distance:\n')
fprintf('kh = %.3g\n',worst.kh)
fprintf('kl = %.3g\n',worst.kl)
fprintf('s_rms = %.3g\n',worst.s_rms)
fprintf('theta95_full: 60 deg = %.3f, 80 deg = %.3f, delta = %.3g deg\n', ...
    worst_comparison.theta95_full_60deg,worst_comparison.theta95_full_80deg, ...
    worst_comparison.delta_theta95_full_deg)
fprintf('theta99_full: 60 deg = %.3f, 80 deg = %.3f, delta = %.3g deg\n', ...
    worst_comparison.theta99_full_60deg,worst_comparison.theta99_full_80deg, ...
    worst_comparison.delta_theta99_full_deg)
fprintf('TV: 60 deg = %.6g, 80 deg = %.6g, delta = %.3g\n', ...
    worst_comparison.TV_60deg,worst_comparison.TV_80deg, ...
    worst_comparison.delta_TV)
fprintf('weighted_phase_rms = %.6g rad\n',worst.weighted_phase_rms_rad)
fprintf('F_phase_lt_0.1rad = %.6g\n',worst.F_phase_lt_0p1rad)
fprintf('F_phase_lt_0.3rad = %.6g\n',worst.F_phase_lt_0p3rad)
fprintf('Maximum clipped Hankel negative / positive peak = %.6g\n', ...
    max_negative_relative)

function [p_full,p_2k,w_full,w_2k,negative_rel] = ...
        local_ssa_spectra(theta_s,kh,kl,k)
% Lowest-order isotropic two-dimensional SSA angular shape.
h = kh/k;
l = kl/k;
rho = linspace(0,6*l,4001);
C_eta = h^2*exp(-(rho/l).^2);
gamma_i = k;
gamma_s = k*cos(theta_s);
v_z = gamma_i+gamma_s;
q = k*sin(theta_s);
geometry = gamma_i^2*gamma_s.^2./v_z.^2;
F_full = zeros(size(theta_s));
F_2k = zeros(size(theta_s));

for jj = 1:numel(theta_s)
    J0 = besselj(0,q(jj)*rho);
    % exp(-a^2*h^2)*expm1(a^2*C_eta) is the stable form of the
    % difference of exponentials in the requested Hankel integral.
    rough_full = exp(-(v_z(jj)*h)^2)*expm1(v_z(jj)^2*C_eta);
    rough_2k = exp(-(2*k*h)^2)*expm1((2*k)^2*C_eta);
    F_full(jj) = 2*pi*trapz(rho,rho.*J0.*rough_full);
    F_2k(jj) = 2*pi*trapz(rho,rho.*J0.*rough_2k);
end

[F_full,negative_rel_full] = ...
    local_check_hankel_values(F_full,kh,kl,'full');
[F_2k,negative_rel_2k] = ...
    local_check_hankel_values(F_2k,kh,kl,'2k');
negative_rel = [negative_rel_full negative_rel_2k];
sigma_full = geometry.*F_full;
sigma_2k = geometry.*F_2k;
w_full = sigma_full.*sin(theta_s);
w_2k = sigma_2k.*sin(theta_s);
p_full = w_full/trapz(theta_s,w_full);
p_2k = w_2k/trapz(theta_s,w_2k);
end

function [values,max_negative_relative] = ...
        local_check_hankel_values(values,kh,kl,label)
scale = max(values);
negative_values = values(values < 0);
if isempty(negative_values)
    max_negative_relative = 0;
else
    max_negative_relative = max(abs(negative_values))/scale;
end
tiny_negative = values < 0 & abs(values) < 1e-8*scale;
values(tiny_negative) = 0;
if any(values < 0)
    warning(['Negative Hankel values for kh=%g, kl=%g, %s phase; ', ...
        'clipping them to zero (max relative magnitude %.3g).'], ...
        kh,kl,label,max_negative_relative)
    values(values < 0) = 0;
end
end

function [theta_full,theta_2k,TV_distance] = ...
        local_distribution_metrics(theta,theta_deg,p_full,p_2k)
cdf_full = cumtrapz(theta,p_full);
cdf_2k = cumtrapz(theta,p_2k);
theta_full = arrayfun(@(fraction) local_quantile( ...
    theta_deg,cdf_full,fraction),[0.90 0.95 0.99]);
theta_2k = arrayfun(@(fraction) local_quantile( ...
    theta_deg,cdf_2k,fraction),[0.90 0.95 0.99]);
TV_distance = 0.5*trapz(theta,abs(p_2k-p_full));
end

function angle_deg = local_quantile(theta_deg,cdf_values,fraction)
idx = find(cdf_values >= fraction,1,'first');
if idx == 1
    angle_deg = theta_deg(1);
else
    angle_deg = interp1(cdf_values(idx-1:idx),theta_deg(idx-1:idx),fraction);
end
end

function local_write_report(file_name,geometry_table,summary_table, ...
        comparison_table,negative_rel_test,max_s_rms)
[max_TV,max_TV_idx] = max(summary_table.TV_distance);
[max_phase,max_phase_idx] = max(summary_table.weighted_phase_rms_rad);
[max_negative_full,max_negative_full_idx] = ...
    max(summary_table.max_negative_relative_full);
[max_negative_2k,max_negative_2k_idx] = ...
    max(summary_table.max_negative_relative_2k);
max_delta95_full = max(abs(comparison_table.delta_theta95_full_deg));
[max_delta99_full,max_delta99_full_idx] = ...
    max(abs(comparison_table.delta_theta99_full_deg));
max_delta95_2k = max(abs(comparison_table.delta_theta95_2k_deg));
[max_delta99_2k,max_delta99_2k_idx] = ...
    max(abs(comparison_table.delta_theta99_2k_deg));
[max_delta_TV,max_delta_TV_idx] = max(abs(comparison_table.delta_TV));
[~,TV_order] = sort(summary_table.TV_distance,'descend');
worst_rows = TV_order(1:min(3,numel(TV_order)));
n_01 = sum(summary_table.F_phase_lt_0p1rad >= 0.9);
n_03 = sum(summary_table.F_phase_lt_0p3rad >= 0.9);

fid = fopen(file_name,'w');
fprintf(fid,'# Validation of the near-normal 2k height-phase approximation\n\n');
fprintf(fid,'## 1. Validation goal\n\n');
fprintf(fid,['This independent calculation studies only the replacement ', ...
    '`(gamma_i + gamma_s) eta -> 2k eta` in the roughness height phase. ', ...
    'The `2k-phase` branch is a diagnostic construction, not a new SSA theory. ', ...
    'It does not validate the Kirchhoff approximation, SSA as an ocean truth model, ', ...
    'PM/TMA/JONSWAP, PE propagation, or KStat. Angles are measured from the ', ...
    'surface normal. The primary SSA angular distributions are normalized over ', ...
    '0--80 deg and compared with the same raw spectra normalized over 0--60 deg.\n\n']);

fprintf(fid,'## 2. Analytic geometry result\n\n');
fprintf(fid,'| theta_i (deg) | max theta_s at 1%% (deg) | max theta_s at 3%% (deg) | max theta_s at 5%% (deg) |\n');
fprintf(fid,'|---:|---:|---:|---:|\n');
for ii = 1:height(geometry_table)
    fprintf(fid,'| %.1f | %s | %s | %s |\n',geometry_table.theta_i_deg(ii), ...
        local_number(geometry_table.theta_s_max_1pct_deg(ii),'%.1f'), ...
        local_number(geometry_table.theta_s_max_3pct_deg(ii),'%.1f'), ...
        local_number(geometry_table.theta_s_max_5pct_deg(ii),'%.1f'));
end

fprintf(fid,'\n## 3. Gaussian SSA results\n\n');
fprintf(fid,['The lowest-order two-dimensional isotropic Gaussian SSA shape retains ', ...
    'all angle-dependent geometry factors. Only the roughness-statistics argument ', ...
    '`a` changes between the full and 2k-phase branches. Cases with ', ...
    '`s_rms > %.2f` are skipped. Each Hankel integral uses 4001 points over ', ...
    '`rho = 0--6l`; negative finite-quadrature tail values are warned about and ', ...
    'clipped to zero before angular normalization.\n\n'],max_s_rms);
fprintf(fid,'| kh | kl | s_rms | theta95 full | theta99 full | theta95 2k | TV | weighted phase RMS | F(<0.1 rad) | F(<0.3 rad) |\n');
fprintf(fid,'|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n');
for ii = 1:height(summary_table)
    fprintf(fid,'| %.2g | %.2g | %.3f | %.2f | %.2f | %.2f | %.5f | %.5f | %.4f | %.4f |\n', ...
        summary_table.kh(ii),summary_table.kl(ii),summary_table.s_rms(ii), ...
        summary_table.theta95_full_deg(ii),summary_table.theta99_full_deg(ii), ...
        summary_table.theta95_2k_deg(ii),summary_table.TV_distance(ii), ...
        summary_table.weighted_phase_rms_rad(ii), ...
        summary_table.F_phase_lt_0p1rad(ii),summary_table.F_phase_lt_0p3rad(ii));
end

fprintf(fid,['\nBefore clipping, the largest negative Hankel value relative to the ', ...
    'positive peak is %.3e for the full branch at `(kh, kl) = (%.2g, %.2g)` ', ...
    'and %.3e for the 2k branch at `(kh, kl) = (%.2g, %.2g)`. The weak ', ...
    '`kh=0.01, kl=20` sanity case gives %.3e and %.3e, respectively.\n'], ...
    max_negative_full,summary_table.kh(max_negative_full_idx), ...
    summary_table.kl(max_negative_full_idx),max_negative_2k, ...
    summary_table.kh(max_negative_2k_idx), ...
    summary_table.kl(max_negative_2k_idx),negative_rel_test(1), ...
    negative_rel_test(2));

fprintf(fid,'\n## 4. Scattering-angle limit check\n\n');
fprintf(fid,['The table below shows the three largest 80-deg TV cases. Deltas are ', ...
    '`80-deg result - 60-deg result`, computed from the same raw Hankel spectra.\n\n']);
fprintf(fid,'| kh | kl | theta95 full: 60 | theta95 full: 80 | delta95 | theta99 full: 60 | theta99 full: 80 | delta99 | TV60 | TV80 | delta TV |\n');
fprintf(fid,'|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n');
for jj = 1:numel(worst_rows)
    ii = worst_rows(jj);
    fprintf(fid,'| %.2g | %.2g | %.4f | %.4f | %.3g | %.4f | %.4f | %.3g | %.6f | %.6f | %.3g |\n', ...
        comparison_table.kh(ii),comparison_table.kl(ii), ...
        comparison_table.theta95_full_60deg(ii), ...
        comparison_table.theta95_full_80deg(ii), ...
        comparison_table.delta_theta95_full_deg(ii), ...
        comparison_table.theta99_full_60deg(ii), ...
        comparison_table.theta99_full_80deg(ii), ...
        comparison_table.delta_theta99_full_deg(ii), ...
        comparison_table.TV_60deg(ii),comparison_table.TV_80deg(ii), ...
        comparison_table.delta_TV(ii));
end
fprintf(fid,['\nAcross all retained cases, the largest absolute full-phase changes are ', ...
    '%.3g deg in theta95 and %.3g deg in theta99; the largest 2k-phase ', ...
    'changes are %.3g deg and %.3g deg. The largest absolute TV change is ', ...
    '%.3g. Full results are in `ssa_angle_limit_comparison.csv`.\n'], ...
    max_delta95_full,max_delta99_full,max_delta95_2k,max_delta99_2k, ...
    max_delta_TV);
fprintf(fid,['\nThe worst-TV case remains `(kh, kl) = (%.2g, %.2g)`. The full-phase ', ...
    'percentiles are essentially stable: their largest theta95 and theta99 ', ...
    'changes occur at `(%.2g, %.2g)` and are %.3g deg and %.3g deg. ', ...
    'The 2k-phase high-angle tail is less converged at 60 deg: `(%.2g, %.2g)` ', ...
    'changes by %.3g deg in theta95 and %.3g deg in theta99. TV changes by at ', ...
    'most %.3g at `(%.2g, %.2g)`, so the qualitative TV ranking and scale are ', ...
    'stable, but the worst-case TV values are not numerically identical.\n'], ...
    summary_table.kh(max_TV_idx),summary_table.kl(max_TV_idx), ...
    comparison_table.kh(max_delta99_full_idx), ...
    comparison_table.kl(max_delta99_full_idx),max_delta95_full, ...
    max_delta99_full,comparison_table.kh(max_delta99_2k_idx), ...
    comparison_table.kl(max_delta99_2k_idx),max_delta95_2k, ...
    max_delta99_2k,max_delta_TV,comparison_table.kh(max_delta_TV_idx), ...
    comparison_table.kl(max_delta_TV_idx));

fprintf(fid,'\n## 5. Figures\n\n');
fprintf(fid,'![Analytic phase-coefficient error](../results/validation/ssa_2k_phase/phase_coefficient_error_map.png)\n\n');
fprintf(fid,'![Representative normalized diffuse spectra](../results/validation/ssa_2k_phase/ssa_spectrum_examples.png)\n\n');
fprintf(fid,'![SSA 2k error summary](../results/validation/ssa_2k_phase/ssa_2k_error_summary.png)\n\n');

fprintf(fid,'## 6. Short interpretation\n\n');
fprintf(fid,['Across the retained Gaussian cases, the full-phase 90%%, 95%%, and 99%% ', ...
    'energy angles span %.2f--%.2f deg, %.2f--%.2f deg, and %.2f--%.2f deg, ', ...
    'respectively, within the normalized 0--80 deg window.\n\n'], ...
    min(summary_table.theta90_full_deg),max(summary_table.theta90_full_deg), ...
    min(summary_table.theta95_full_deg),max(summary_table.theta95_full_deg), ...
    min(summary_table.theta99_full_deg),max(summary_table.theta99_full_deg));
fprintf(fid,['Replacing only the height-phase argument by 2k gives TV distances ', ...
    'from %.5f to %.5f. The largest TV distance occurs at `(kh, kl) = ', ...
    '(%.2g, %.2g)`; its full-phase theta95 is %.2f deg.\n\n'], ...
    min(summary_table.TV_distance),max_TV,summary_table.kh(max_TV_idx), ...
    summary_table.kl(max_TV_idx),summary_table.theta95_full_deg(max_TV_idx));
fprintf(fid,['The largest energy-weighted phase RMS is %.5f rad at `(kh, kl) = ', ...
    '(%.2g, %.2g)`. This and TV identify the strongest tested diagnostic ', ...
    'differences without defining a valid/invalid threshold.\n\n'],max_phase, ...
    summary_table.kh(max_phase_idx),summary_table.kl(max_phase_idx));
fprintf(fid,['%d of %d retained cases place at least 90%% of their normalized ', ...
    'diffuse energy below 0.1 rad phase mismatch; %d of %d do so below ', ...
    '0.3 rad. These are diagnostic fractions, not universal applicability limits.\n\n'], ...
    n_01,height(summary_table),n_03,height(summary_table));
fprintf(fid,['The result addresses only the tested Gaussian covariance and angular ', ...
    'window. It does not establish applicability for other sea spectra, and the ', ...
    'lowest-order SSA calculation is not an exact scattering result.\n']);
fclose(fid);
end

function text_value = local_number(value,format_spec)
if isnan(value)
    text_value = 'NaN';
else
    text_value = sprintf(format_spec,value);
end
end
