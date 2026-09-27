function summary = validate_kachoyan_macaskill_1987(overrides)
%VALIDATE_KACHOYAN_MACASKILL_1987 Reproduce the 1987 rough-surface method.
%   This is a validation-only workflow. It writes auditable MAT/CSV/PNG
%   artifacts under results/validation/kachoyan_macaskill_1987/ and a
%   report under reports/. The paper's original random arrays are not
%   published; all random surfaces here therefore record a fixed seed.

if nargin < 1 || isempty(overrides), overrides = struct(); end
if ~isstruct(overrides) || ~isscalar(overrides), error('overrides must be a scalar struct.'); end
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
setup_vertical_project();
out_root = fullfile(root,'results','validation','kachoyan_macaskill_1987');
if ~exist(out_root,'dir'), mkdir(out_root); end
report_file = fullfile(root,'reports','kachoyan_macaskill_1987_reproduction_report.md');

cfg = struct('run_mode','paper','N_gaussian',256,'points_per_correlation',32, ...
    'gaussian_seed',198705,'pm_seed',198706,'surface_window','none', ...
    'z_count',96,'z_margin',0.002,'make_figures',true,'fail_on_check',false, ...
    'output_dir',out_root,'report_file',report_file);
names = fieldnames(overrides);
for ii = 1:numel(names)
    if ~isfield(cfg,names{ii}), error('Unknown override %s.',names{ii}); end
    cfg.(names{ii}) = overrides.(names{ii});
end
if strcmpi(cfg.run_mode,'smoke')
    cfg.N_gaussian = min(cfg.N_gaussian,64);
    cfg.points_per_correlation = min(cfg.points_per_correlation,16);
    cfg.z_count = min(cfg.z_count,24);
    cfg.make_figures = false;
end

% A fixed Gaussian-correlation realization, with 32 points per correlation
% length in the paper's convention. Its seed is a replacement provenance,
% not the unavailable Fig. 1 array.
g = local_gaussian_surface(cfg.N_gaussian,cfg.points_per_correlation, ...
    cfg.gaussian_seed);
pm = local_pm_surface(cfg.N_gaussian, cfg.N_gaussian/cfg.points_per_correlation, ...
    cfg.pm_seed);

cases = struct('id',{},'family',{},'description',{},'x',{},'xi',{}, ...
    'xi_prime',{},'xi_second',{},'dx',{},'khat',{},'hhat',{}, ...
    'theta_deg',{},'periodic_replicas',{},'z_samples',{},'seed',{}, ...
    'surface_source',{},'mode',{},'window',{});
ci = 0;
% Gaussian Fig. 2/3/4/5 analogs and one PM Fig. 8 analog.
ci=ci+1; cases(ci)=local_case('gaussian_fig2','gaussian', ...
    'Gaussian substitute; khat=50 hhat=0.1 theta=90',g,50,0.1,90,0, ...
    cfg.gaussian_seed,cfg,'gaussian_seed_substitute');
ci=ci+1; cases(ci)=local_case('gaussian_center128','gaussian', ...
    'Center 128-point truncation of Fig. 2 substitute', ...
    local_center_surface(g,128),50,0.1,90,0,cfg.gaussian_seed,cfg, ...
    'center_128_of_gaussian_seed_substitute');
ci=ci+1; cases(ci)=local_case('gaussian_fig4','gaussian', ...
    'Gaussian substitute; khat=15 hhat=1/3 theta=90',g,15,1/3,90,0, ...
    cfg.gaussian_seed,cfg,'gaussian_seed_substitute');
ci=ci+1; cases(ci)=local_case('gaussian_fig5','gaussian', ...
    'Gaussian substitute; khat=15 hhat=1/3 theta=45',g,15,1/3,45,0, ...
    cfg.gaussian_seed,cfg,'gaussian_seed_substitute');
ci=ci+1; cases(ci)=local_case('pm_fig8','modified_pm', ...
    'Modified Pierson-Moskowitz substitute; khat=50 hhat=0.36 theta=90', ...
    pm,50,0.36,90,0,cfg.pm_seed,cfg,'pm_spectrum_substitute');

% Table II controlled periodic sine cases: N=32, NQ=4, one period.
sine = local_sine_surface(32,1);
ci=ci+1; cases(ci)=local_case('sine_tableII_50_01_90','sinusoid', ...
    'Table II sine; khat=50 hhat=0.1 theta=90',sine,50,0.1,90,4, ...
    NaN,cfg,'analytic_sine');
ci=ci+1; cases(ci)=local_case('sine_tableII_15_13_45','sinusoid', ...
    'Table II sine; khat=15 hhat=1/3 theta=45',sine,15,1/3,45,4, ...
    NaN,cfg,'analytic_sine');
ci=ci+1; cases(ci)=local_case('sine_tableII_15_13_90','sinusoid', ...
    'Table II sine; khat=15 hhat=1/3 theta=90',sine,15,1/3,90,4, ...
    NaN,cfg,'analytic_sine');

records = repmat(struct('case_id','','family','','description','','N',0,'dx',0, ...
    'khat',0,'hhat',0,'theta_deg',0,'periodic_replicas',0,'seed',NaN, ...
    'surface_source','','propagating_count',0,'linear_residual_full',NaN, ...
    'linear_residual_ka',NaN,'dft_system_error',NaN,'energy_full',NaN,'energy_ka',NaN, ...
    'finite',false,'z_min',NaN,'surface_max',NaN),numel(cases),1);
results = struct();
for ii=1:numel(cases)
    c = cases(ii);
    res = solve_kachoyan_macaskill_1987(c);
    results.(c.id) = res;
    case_dir = fullfile(cfg.output_dir,c.id);
    if ~exist(case_dir,'dir'), mkdir(case_dir); end
    save(fullfile(case_dir,[c.id '.mat']),'res','c','-v7');
    if cfg.make_figures
        local_write_case_figures(case_dir,c,res);
    end
    records(ii).case_id = c.id; records(ii).family = c.family;
    records(ii).description = c.description; records(ii).N = numel(c.x);
    records(ii).dx = c.dx; records(ii).khat = c.khat; records(ii).hhat = c.hhat;
    records(ii).theta_deg = c.theta_deg; records(ii).periodic_replicas = c.periodic_replicas;
    records(ii).seed = c.seed; records(ii).surface_source = c.surface_source;
    records(ii).propagating_count = res.propagating_count;
    records(ii).linear_residual_full = res.linear_residual_full;
    records(ii).linear_residual_ka = res.linear_residual_ka; records(ii).dft_system_error = res.dft_system_error;
    records(ii).energy_full = res.energy_full; records(ii).energy_ka = res.energy_ka;
    records(ii).finite = res.finite; records(ii).z_min = min(c.z_samples);
    records(ii).surface_max = max(c.hhat*c.xi);
end
summary_table = struct2table(records);
writetable(summary_table,fullfile(cfg.output_dir,'kachoyan_macaskill_case_summary.csv'));
save(fullfile(cfg.output_dir,'kachoyan_macaskill_1987_results.mat'), ...
    'results','summary_table','cfg','g','pm','-v7');
if cfg.make_figures
    local_write_overview_figure(cfg.output_dir,cases,results);
end
local_write_report(report_file,cfg,summary_table,g);
summary = struct('output_dir',cfg.output_dir,'report_file',report_file, ...
    'summary_table',summary_table,'results',results, ...
    'all_finite',all(summary_table.finite), ...
    'max_full_residual',max(summary_table.linear_residual_full), ...
    'case_count',height(summary_table));
if cfg.fail_on_check && (~summary.all_finite || summary.max_full_residual > 1e-8)
    error('Kachoyan validation hard check failed.');
end
end

function c = local_case(id,family,description,surface,khat,hhat,theta,reps,seed,cfg,source)
c = struct('id',id,'family',family,'description',description,'x',surface.x, ...
    'xi',surface.xi,'xi_prime',surface.xi_prime,'xi_second',surface.xi_second, ...
    'dx',surface.dx,'khat',khat,'hhat',hhat,'theta_deg',theta, ...
    'periodic_replicas',reps,'z_samples',local_z_grid(surface.xi,hhat,cfg.z_count,cfg.z_margin), ...
    'seed',seed,'surface_source',source,'mode','full_and_ka','window',local_surface_window(numel(surface.x),cfg.surface_window));
end

function z = local_z_grid(xi,hhat,n,margin)
z = linspace(max(hhat*xi)+margin,max(hhat*xi)+0.30,n);
end

function w = local_surface_window(N,kind)
switch lower(char(kind))
    case 'none'
        w = ones(N,1);
    case 'hann'
        w = 0.5 - 0.5*cos(2*pi*(0:N-1)'/N);
    otherwise
        error('Unsupported surface_window %s; use none or hann.',char(kind));
end
end

function s = local_gaussian_surface(N,ppc,seed)
if mod(N,2)~=0 || ppc<2, error('Gaussian N must be even and ppc >=2.'); end
T = N/ppc; dx = 1/ppc; x = ((0:N-1)-N/2).'*dx;
q = 2*pi/T*[0:(N/2-1),-N/2:-1].';
rng(mod(round(seed),2^32),'twister');
amp = exp(-0.25*(q*T/(2*pi)).^2); amp(1)=0;
coef = amp .* (randn(N,1)+1i*randn(N,1));
coef(1)=0; coef(N/2+1)=real(coef(N/2+1));
xi = real(ifft(coef)); xi = xi-mean(xi); xi = xi/max(std(xi,1),realmin);
xi_prime = real(ifft(1i*q.*fft(xi))); xi_second = real(ifft(-(q.^2).*fft(xi)));
s = struct('x',x,'xi',xi,'xi_prime',xi_prime,'xi_second',xi_second,'dx',dx, ...
    'period_length',T,'points_per_correlation',ppc);
end

function s = local_pm_surface(N,T,seed)
if mod(N,2)~=0, error('PM N must be even.'); end
dx=T/N; x=((0:N-1)-N/2).'*dx;
q=2*pi/T*[0:(N/2-1),-N/2:-1].'; qabs=abs(q);
g=9.81; U=6; alpha=8.10e-3; beta=0.74;
S=zeros(N,1); m=qabs>0; S(m)=alpha./(2*qabs(m).^3).*exp(-beta*g^2./(U^4*qabs(m).^2));
rng(mod(round(seed),2^32),'twister');
coef=zeros(N,1); pos=(2:N/2).'; amp=sqrt(max(S(pos),0)*2/T);
phase=2*pi*rand(size(pos)); coef(pos)=amp.*exp(1i*phase); coef(N-pos+2)=conj(coef(pos));
xi=real(ifft(coef))*N; xi=xi-mean(xi); xi=xi/max(std(xi,1),realmin);
xi_prime=real(ifft(1i*q.*fft(xi))); xi_second=real(ifft(-(q.^2).*fft(xi)));
s=struct('x',x,'xi',xi,'xi_prime',xi_prime,'xi_second',xi_second,'dx',dx, ...
    'period_length',T,'wind_speed_mps',U,'spectrum_q',S);
end

function s = local_sine_surface(N,T)
dx=T/N; x=(0:N-1).'*dx; xi=sin(2*pi*x/T);
q=2*pi/T*[0:(N/2-1),-N/2:-1].';
xi_prime=real(ifft(1i*q.*fft(xi))); xi_second=real(ifft(-(q.^2).*fft(xi)));
s=struct('x',x,'xi',xi,'xi_prime',xi_prime,'xi_second',xi_second,'dx',dx);
end

function s = local_center_surface(s0,Nc)
if Nc>=numel(s0.x), s=s0; return; end
mid=floor(numel(s0.x)/2)+1; idx=(mid-Nc/2):(mid+Nc/2-1);
s=s0; s.x=s0.x(idx)-mean(s0.x(idx)); s.xi=s0.xi(idx);
s.dx=s0.dx; q=2*pi/(Nc*s.dx)*[0:(Nc/2-1),-Nc/2:-1].';
s.xi=s.xi-mean(s.xi); s.xi=s.xi/max(std(s.xi,1),realmin);
s.xi_prime=real(ifft(1i*q.*fft(s.xi))); s.xi_second=real(ifft(-(q.^2).*fft(s.xi)));
end

function local_write_case_figures(case_dir,c,res)
fig=figure('Visible','off'); plot(c.x,c.hhat*c.xi,'k-','LineWidth',1.2); grid on;
xlabel('x/L'); ylabel('h/L \\xi(x)'); title(c.description,'Interpreter','none');
exportgraphics(fig,fullfile(case_dir,'surface.png'),'Resolution',150); close(fig);
fig=figure('Visible','off'); imagesc(c.x,res.z_samples,10*log10(max(res.nearfield_intensity_full,1e-14))); axis xy; colorbar;
xlabel('x/L'); ylabel('z'); title('Full-solution scattered intensity (dB)','Interpreter','none');
exportgraphics(fig,fullfile(case_dir,'nearfield_full.png'),'Resolution',150); close(fig);
fig=figure('Visible','off');
mask=res.propagating_mask; polarplot(deg2rad(atan2(res.gamma_mode(mask),res.q_mode(mask))), ...
    abs(res.R_full(mask)).^2.*res.gamma_mode(mask)/res.gamma0,'k.');
title(sprintf('Full reflected modal energy, E=%.6g',res.energy_full));
exportgraphics(fig,fullfile(case_dir,'polar_energy_full.png'),'Resolution',150); close(fig);
fig=figure('Visible','off'); imagesc(c.x,res.z_samples,10*log10(max(res.nearfield_intensity_ka,1e-14))); axis xy; colorbar;
xlabel('x/L'); ylabel('z'); title('Kirchhoff approximation intensity (dB)','Interpreter','none');
exportgraphics(fig,fullfile(case_dir,'nearfield_ka.png'),'Resolution',150); close(fig);
end

function local_write_overview_figure(outdir,cases,results)
fig=figure('Visible','off'); hold on; grid on;
for ii=1:numel(cases)
    r=results.(cases(ii).id); plot(cases(ii).theta_deg,r.energy_full,'ko','MarkerSize',5);
end
xlabel('incident angle from horizontal (deg)'); ylabel('E'); title('Kachoyan-Macaskill energy checks');
exportgraphics(fig,fullfile(outdir,'energy_overview.png'),'Resolution',150); close(fig);
end

function local_write_report(path,cfg,tab,g)
fid=fopen(path,'w'); if fid<0, error('Cannot open report %s.',path); end
cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'# Kachoyan--Macaskill (1987) reproduction audit\n\n');
fprintf(fid,'Generated by `scripts/validation/validate_kachoyan_macaskill_1987.m` on %s.\n\n',char(datetime("now","Format","yyyy-MM-dd HH:mm:ss")));
fprintf(fid,'## Scope and reproduction level\n\n');
fprintf(fid,'This is a validation-only MATLAB implementation of the dimensionless Fredholm collocation method and the reflection-mode reconstruction in B. J. Kachoyan and C. Macaskill, JASA 82(5), 1720--1726 (1987), DOI 10.1121/1.395165. The paper does not publish the exact 256-point Gaussian or modified Pierson--Moskowitz arrays, random seeds, or all kernel implementation details needed for pixel-identical figures. The Gaussian and PM figures below therefore use fixed, auditable replacement realizations. Equation/discretization/metric checks are strict at the implemented formula level; figure-level numerical agreement is approximate and is not claimed as a successful pixel-level reproduction.\n\n');
fprintf(fid,'## Formula mapping\n\n');
fprintf(fid,'- Surface: `z=hhat*xi(x)` with `khat=kL`, `hhat=h/L`; incident angle is from the horizontal, `alpha0=cos(theta)`, `gamma0=sin(theta)`.\n');
fprintf(fid,"- Boundary forcing: psi_inc uses gamma0+alpha0*hhat*xi-prime and exp(i*khat*(alpha0*x-gamma0*hhat*xi)).\n");
fprintf(fid,"- Fredholm equation: (I+dx*K)*psi=2*psi_inc; Eq.(5) uses the H1_1 kernel, rho^2=(xprime-x)^2+hhat^2*(xiprime-xi)^2, and the curvature diagonal limit.\n");
fprintf(fid,"- Kernel implementation: K=-(i*pi/2)*khat*hhat^2*[xi_source-xi_target-dx*xi_prime_target]*H1_1(khat*rho)/rho; diagonal=0.5*hhat^2*xi_second/(1+hhat^2*xi_prime^2).\n");
fprintf(fid,"- Reflection modes: q_r=alpha0+2*pi*r/(khat*N*dx), gamma_r=sqrt(1-q_r^2); carrier-demodulated psi_hat and C use the paper DFT convention and only real propagating modes.\n");
fprintf(fid,"- Reflection coefficient evaluation: psi_hat uses carrier-demodulated ifft; C uses fft(exp(-i*khat*gamma*hhat*xi))/N; the circular lookup is C at r-minus-m.\n");
fprintf(fid,'- Energy: `E=Re(sum(abs(R_r).^2.*gamma_r/gamma0))`. Near-field reconstruction is only evaluated at `z>hhat*max(xi)` as required by Appendix A.\n\n');
fprintf(fid,'## Provenance and settings\n\n');
fprintf(fid,'- MATLAB run mode: `%s`; Gaussian N=%d, points per correlation length=%d, seed=%d; PM seed=%d.\n',cfg.run_mode,numel(g.x),g.points_per_correlation,cfg.gaussian_seed,cfg.pm_seed);
fprintf(fid,'- Gaussian substitute is periodic, zero-mean and RMS-normalized to one before applying each requested hhat. PM substitute uses the repository PM spectrum shape with U=6 m/s and is RMS-normalized to one before applying hhat. These normalizations are explicit replacements, not recovered paper data.\n');
fprintf(fid,'- Table-II sine cases use N=32 and NQ=4 over one analytic period; window mode is `%s`.\n\n',cfg.surface_window);
fprintf(fid,'## Results from this MATLAB run\n\n');
fprintf(fid,'| case | N | khat | hhat | theta (deg) | Nprop | residual full | DFT-system err | E full | E KA | finite |\n|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|\n');
for ii=1:height(tab)
    fprintf(fid,'| %s | %d | %.6g | %.6g | %.6g | %d | %.3e | %.6g | %.6g | %.6g | %d |\n', ...
        tab.case_id{ii},tab.N(ii),tab.khat(ii),tab.hhat(ii),tab.theta_deg(ii), ...
        tab.propagating_count(ii),tab.linear_residual_full(ii),tab.dft_system_error(ii),tab.energy_full(ii), ...
        tab.energy_ka(ii),tab.finite(ii));
end
fprintf(fid,'\n');
fprintf(fid,'## Paper readouts and comparison limits\n\n');
fprintf(fid,'Paper readouts used for orientation are Fig. 2 `E approximately 1.000`, Fig. 4 `E=0.993`, Fig. 5 `E=0.927`, Fig. 6 KA `E=0.813`, Fig. 8 PM `E=0.954`, and Table II FS/KA values `(0.9969,0.9852)`, `(1.108,3.435)`, `(1.190,0.231)`. They are not substituted into the run table and are not treated as validation evidence for the replacement realizations.\n\n');
fprintf(fid,'## Validation gates\n\n');
fprintf(fid,'- Finite values: `%d/%d` cases passed.\n',sum(tab.finite),height(tab));
fprintf(fid,'- Full linear-system residual: maximum `%.3e`; independent DFT-system solve error: maximum `%.3e`.\n',max(tab.linear_residual_full),max(tab.dft_system_error));
fprintf(fid,'- All near-field z samples satisfy the surface-clearance inequality by construction; the case table records the minimum z and maximum surface height.\n');
fprintf(fid,'- The full and Kirchhoff paths are both computed from the same surface and forcing; KA is obtained by setting the integral operator to zero.\n\n');
fprintf(fid,'## Known limitations and next steps\n\n');
fprintf(fid,'1. Exact random surface arrays and the original seed are unavailable, so Fig. 2/4/5/8 and the central-128 comparison can only be approximate.\n');
fprintf(fid,'2. Eq.(5) is implemented with the hhat-squared graph-normal form and its curvature limit documented above; the paper does not expose the authors'' quadrature/roundoff choices, so the numerical implementation remains an auditable interpretation.\n');
fprintf(fid,'3. The modified PM substitute is a one-dimensional spectral realization based on the repository PM convention, not a recovered 1987 sea realization.\n');
fprintf(fid,'4. If closer numerical agreement is required, obtain the authors'' original realization arrays or digitize the plotted surface and rerun with that exact input.\n');
end
