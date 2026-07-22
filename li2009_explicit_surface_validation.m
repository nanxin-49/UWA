function result = li2009_explicit_surface_validation(paramsV)
%LI2009_EXPLICIT_SURFACE_VALIDATION Independent Li et al. (2009) acoustic case.
%   This validation-only path implements a monostatic bottom-to-surface-to-
%   bottom geometry without calling or changing the public communication path.
%   One explicit raw-PM surface is shared by every frequency of each pulse.

if nargin < 1 || isempty(paramsV)
    paramsV = struct();
end
cfg = local_config(paramsV);
if cfg.verbose
    fprintf('Li2009 independent case: H=%g m, PE=%dx%d over %gx%g m, PM=%dx%d over %gx%g m\n', ...
        cfg.H_m, cfg.pe_nx, cfg.pe_ny, cfg.pe_xw_m, cfg.pe_yw_m, ...
        cfg.pm_nx, cfg.pm_ny, cfg.pm_xw_m, cfg.pm_yw_m);
end

cache_timer = tic;
cache = local_build_pe_cache(cfg);
cache_elapsed_s = toc(cache_timer);
[tx_waveform, tx_metrics] = local_synthesize_waveform(ones(cfg.Nf,1), cfg);

flat_h_f = local_receiver_response(cache, zeros(cfg.pe_ny,cfg.pe_nx));
[flat_waveform, flat_relative] = local_synthesize_waveform(flat_h_f, cfg);
tau0_s = 2*cfg.H_m/cfg.c0_mps;
flat_uncalibrated_metrics = local_absolute_metrics(flat_relative, tx_metrics, tau0_s);
% A single deterministic flat-surface reference removes source/finite-grid
% group delay. No rough realization is shifted or peak-aligned individually.
flat_metrics = local_absolute_metrics(flat_relative, flat_relative, tau0_s);

U10_list = cfg.U10_list_mps(:).';
Nu = numel(U10_list);
M = numel(cfg.sea_seed_list);
peak_time_s = zeros(M,Nu);
first_time_s = zeros(M,Nu);
surface_audit = repmat(local_empty_surface_audit(),Nu,1);
example_surface = cell(Nu,1);
example_waveform = cell(Nu,1);
radial = cell(Nu,1);
reflection_checks = repmat(local_empty_reflection_check(),Nu,1);
U19p5_list = zeros(size(U10_list));

mc_timer = tic;
for iu = 1:Nu
    [U19p5,ustar,z0] = li2009_u10_to_u195(U10_list(iu));
    U19p5_list(iu) = U19p5;
    pm_spec = raw_pm_spectrum_grid_vertical( ...
        U19p5,cfg.pm_nx,cfg.pm_ny,cfg.pm_xw_m,cfg.pm_yw_m);
    eta_std = zeros(M,1);
    hs = zeros(M,1);
    variance_rel_error = zeros(M,1);
    radial_accum = [];
    radial_count = min(M,cfg.radial_ensemble_count);
    for im = 1:M
        [eta_pm,sample_meta] = sample_raw_pm_surface_vertical( ...
            pm_spec,cfg.sea_seed_list(im));
        eta_pe = eta_pm(cache.pm_iy,cache.pm_ix);
        h_f = local_receiver_response(cache,eta_pe);
        [waveform,relative_metrics] = local_synthesize_waveform(h_f,cfg);
        absolute_metrics = local_absolute_metrics(relative_metrics,flat_relative,tau0_s);
        peak_time_s(im,iu) = absolute_metrics.peak_time_s;
        first_time_s(im,iu) = absolute_metrics.first_threshold_time_s;
        eta_std(im) = std(eta_pm(:),1);
        hs(im) = 4*eta_std(im);
        variance_rel_error(im) = sample_meta.energy_rel_error_about_zero;
        if im <= radial_count
            sample_radial = local_radial_spectrum(eta_pm,pm_spec,cfg.radial_bin_count);
            if isempty(radial_accum)
                radial_accum = sample_radial;
                radial_accum.E1D_sample_mean = zeros(size(sample_radial.E1D_sample));
            end
            radial_accum.E1D_sample_mean = radial_accum.E1D_sample_mean + ...
                sample_radial.E1D_sample/radial_count;
        end
        if im == 1
            example_surface{iu} = eta_pm;
            example_waveform{iu} = waveform;
            reflection_checks(iu) = local_reflection_checks( ...
                cache.incident_surface_xy_f(:,:,cache.idx_f_ref),eta_pe, ...
                cache.k0_rad_per_m(cache.idx_f_ref));
        end
    end
    radial{iu} = radial_accum;
    surface_audit(iu) = struct( ...
        'U10_mps',U10_list(iu),'U19p5_mps',U19p5,'ustar_mps',ustar, ...
        'z0_m',z0,'sigma_eta_sample_mean_m',mean(eta_std), ...
        'sigma_eta_sample_std_m',std(eta_std,1), ...
        'Hs_sample_mean_m',mean(hs),'Hs_sample_std_m',std(hs,1), ...
        'sigma_eta_discrete_target_m',pm_spec.sigma_eta_discrete_m, ...
        'Hs_discrete_target_m',pm_spec.Hs_implied_discrete_m, ...
        'sigma_eta_infinite_target_m',sqrt(pm_spec.sigma_eta_infinite2_m2), ...
        'Hs_infinite_target_m',pm_spec.Hs_implied_infinite_m, ...
        'discrete_spectrum_variance_m2',pm_spec.sigma_eta_discrete2_m2, ...
        'spectrum_capture_ratio',pm_spec.capture_ratio_discrete_to_infinite, ...
        'radial_support_capture_ratio',pm_spec.capture_ratio_radial_support_idealized, ...
        'mean_sample_variance_rel_error',mean(variance_rel_error), ...
        'max_sample_variance_rel_error',max(variance_rel_error));
end
mc_elapsed_s = toc(mc_timer);

[summary_table,fit] = local_summary_table( ...
    U10_list,peak_time_s,first_time_s,flat_metrics);
sample_convergence_table = local_sample_convergence( ...
    U10_list,peak_time_s,first_time_s,cfg.sample_convergence_counts);
surface_table = struct2table(surface_audit);
reflection_table = struct2table(reflection_checks);
cost = local_cost_estimate(cfg,cache,cache_elapsed_s,mc_elapsed_s);

result = struct();
result.schema_version = '1.0.0';
result.kind = 'li2009_explicit_raw_pm_kirchhoff_pe_validation';
result.created_at = char(datetime('now','Format','yyyy-MM-dd HH:mm:ss Z'));
result.config = cfg;
result.source_model = struct( ...
    'name','circular_piston_kspace_directivity', ...
    'paper_formula','Eq. (7): 2*J1(k*a*sin(theta))/(k*a*sin(theta))', ...
    'implementation',['The same factor is applied as 2*J1(a*K)/(a*K) on the ', ...
        'transverse PE spectrum; amplitude includes f/fc. Absolute source level is normalized.'], ...
    'limitation',['This is a band-limited PE initialization of the paper piston ', ...
        'directivity, not a resolved 5-cm aperture on the reduced grid.']);
result.geometry = struct('type','monostatic_bottom_surface_bottom', ...
    'H_m',cfg.H_m,'tau0_s',tau0_s,'c0_mps',cfg.c0_mps);
result.shared_surface_rule = ['For each seed eta_pm is generated exactly once, cropped once ', ...
    'to eta_pe, and the same eta_pe array is passed to every frequency.'];
result.cache_audit = cache.audit;
result.cost_estimate = cost;
result.flat = struct('H_reduced_f',flat_h_f,'waveform',flat_waveform, ...
    'tx_waveform',tx_waveform,'metrics',flat_metrics, ...
    'uncalibrated_metrics',flat_uncalibrated_metrics, ...
    'calibration_policy',['One fixed flat-surface source/grid group-delay reference; ', ...
        'no per-realization or per-wind alignment.'], ...
    'uncalibrated_peak_bias_s',flat_uncalibrated_metrics.peak_time_s-tau0_s, ...
    'uncalibrated_first_bias_s',flat_uncalibrated_metrics.first_threshold_time_s-tau0_s, ...
    'tau0_error_peak_s',flat_metrics.peak_time_s-tau0_s, ...
    'tau0_error_first_s',flat_metrics.first_threshold_time_s-tau0_s);
result.U10_list_mps = U10_list;
result.U19p5_list_mps = U19p5_list;
result.sea_seed_list = cfg.sea_seed_list;
result.peak_time_s = peak_time_s;
result.first_threshold_time_s = first_time_s;
result.surface_audit = surface_audit;
result.surface_table = surface_table;
result.reflection_checks = reflection_checks;
result.reflection_table = reflection_table;
result.radial_spectra = radial;
result.example_surface_pm = example_surface;
result.example_waveform = example_waveform;
result.summary_table = summary_table;
result.rayleigh_fit = fit;
result.sample_convergence_table = sample_convergence_table;
result.acceptance = local_acceptance(result);

if cfg.save_artifacts
    local_write_artifacts(result);
end
if cfg.make_figures
    local_make_figures(result);
end
end

function cfg = local_config(paramsV)
defaults = struct( ...
    'run_label','baseline', ...
    'output_dir',fullfile('results','validation','li2009_explicit_surface'), ...
    'H_m',256, ...
    'c0_mps',1500, ...
    'fc_hz',12000, ...
    'pulse_duration_s',6e-3, ...
    'piston_radius_m',0.05, ...
    'U10_list_mps',[5,10], ...
    'sea_seed_list',41001+(0:63), ...
    'pe_nx',128,'pe_ny',128,'pe_xw_m',100,'pe_yw_m',100, ...
    'pm_nx',256,'pm_ny',256,'pm_xw_m',200,'pm_yw_m',200, ...
    'dz_m',4, ...
    'sponge_ratio',0.12, ...
    'alpha_max_np_per_m',0.15, ...
    'frequency_half_span_hz',500, ...
    'Nf',65, ...
    'waveform_nfft',4096, ...
    'threshold_fraction',0.20, ...
    'waveform_window_s',12e-3, ...
    'radial_bin_count',40, ...
    'radial_ensemble_count',16, ...
    'sample_convergence_counts',[8,16,32,64], ...
    'make_figures',true, ...
    'save_artifacts',true, ...
    'verbose',true);
cfg = defaults;
names = fieldnames(paramsV);
for ii = 1:numel(names)
    if ~isfield(defaults,names{ii})
        error('li2009_explicit_surface_validation:UnknownParameter', ...
            'Unknown paramsV field: %s',names{ii});
    end
    cfg.(names{ii}) = paramsV.(names{ii});
end
positive_scalars = {'H_m','c0_mps','fc_hz','pulse_duration_s','piston_radius_m', ...
    'pe_xw_m','pe_yw_m','pm_xw_m','pm_yw_m','dz_m','alpha_max_np_per_m', ...
    'frequency_half_span_hz','waveform_window_s'};
for ii = 1:numel(positive_scalars)
    validateattributes(cfg.(positive_scalars{ii}),{'numeric'},{'scalar','finite','positive'});
end
integer_fields = {'pe_nx','pe_ny','pm_nx','pm_ny','Nf','waveform_nfft', ...
    'radial_bin_count','radial_ensemble_count'};
for ii = 1:numel(integer_fields)
    validateattributes(cfg.(integer_fields{ii}),{'numeric'},{'scalar','integer','positive'});
end
if any(mod([cfg.pe_nx,cfg.pe_ny,cfg.pm_nx,cfg.pm_ny],2)~=0)
    error('All PE and PM grid sizes must be even.');
end
if mod(cfg.Nf,2)~=1 || cfg.Nf<3
    error('Nf must be an odd integer >= 3.');
end
if mod(cfg.waveform_nfft,2)~=0 || cfg.waveform_nfft<cfg.Nf
    error('waveform_nfft must be even and at least Nf.');
end
if cfg.pm_nx<cfg.pe_nx || cfg.pm_ny<cfg.pe_ny
    error('PM grid must not be smaller than the PE grid.');
end
dx_pe = cfg.pe_xw_m/cfg.pe_nx; dy_pe = cfg.pe_yw_m/cfg.pe_ny;
dx_pm = cfg.pm_xw_m/cfg.pm_nx; dy_pm = cfg.pm_yw_m/cfg.pm_ny;
if abs(dx_pe-dx_pm)>1e-12 || abs(dy_pe-dy_pm)>1e-12
    error('PE and PM grids require identical dx/dy for the explicit central crop.');
end
if ~(cfg.sponge_ratio>0 && cfg.sponge_ratio<0.5)
    error('sponge_ratio must lie in (0,0.5).');
end
if ~(cfg.threshold_fraction>0 && cfg.threshold_fraction<1)
    error('threshold_fraction must lie in (0,1).');
end
validateattributes(cfg.U10_list_mps,{'numeric'},{'vector','finite','positive'});
validateattributes(cfg.sea_seed_list,{'numeric'},{'vector','finite'});
cfg.U10_list_mps = cfg.U10_list_mps(:).';
cfg.sea_seed_list = round(cfg.sea_seed_list(:).');
cfg.sample_convergence_counts = unique(round(cfg.sample_convergence_counts(:).'));
cfg.sample_convergence_counts = cfg.sample_convergence_counts( ...
    cfg.sample_convergence_counts>=2 & cfg.sample_convergence_counts<=numel(cfg.sea_seed_list));
cfg.frequency_axis_hz = linspace(cfg.fc_hz-cfg.frequency_half_span_hz, ...
    cfg.fc_hz+cfg.frequency_half_span_hz,cfg.Nf).';
cfg.frequency_offset_hz = cfg.frequency_axis_hz-cfg.fc_hz;
cfg.frequency_spacing_hz = cfg.frequency_axis_hz(2)-cfg.frequency_axis_hz(1);
cfg.maximum_unambiguous_relative_time_s = 1/cfg.frequency_spacing_hz;
cfg.physical_time_resolution_s = 1/(2*cfg.frequency_half_span_hz);
cfg.dx_m = dx_pe; cfg.dy_m = dy_pe;
cfg.output_dir = char(cfg.output_dir);
cfg.run_label = char(cfg.run_label);
end

function cache = local_build_pe_cache(cfg)
x = (-cfg.pe_xw_m/2):cfg.dx_m:(cfg.pe_xw_m/2-cfg.dx_m);
y = (-cfg.pe_yw_m/2):cfg.dy_m:(cfg.pe_yw_m/2-cfg.dy_m);
kx = (2*pi/cfg.pe_xw_m)*[0:(cfg.pe_nx/2-1),-cfg.pe_nx/2:-1];
ky = (2*pi/cfg.pe_yw_m)*[0:(cfg.pe_ny/2-1),-cfg.pe_ny/2:-1];
[KX,KY] = meshgrid(kx,ky); K = hypot(KX,KY); kappa2 = K.^2;
[~,ix0] = min(abs(x)); [~,iy0] = min(abs(y));
alpha_xy = local_absorption_profile(x,y,cfg);
nstep = max(1,ceil(cfg.H_m/cfg.dz_m));
ds = cfg.H_m/nstep;
F = cfg.Nf;
incident = complex(zeros(cfg.pe_ny,cfg.pe_nx,F));
q = complex(zeros(cfg.pe_ny,cfg.pe_nx,F));
k0_list = 2*pi*cfg.frequency_axis_hz/cfg.c0_mps;
piston = ones(size(K));
mask = K>0;
piston(mask) = 2*besselj(1,cfg.piston_radius_m*K(mask)) ./ ...
    (cfg.piston_radius_m*K(mask));
receiver_source = complex(zeros(cfg.pe_ny,cfg.pe_nx));
receiver_source(iy0,ix0) = 1;
for ifq = 1:F
    k0 = k0_list(ifq);
    denom = sqrt(complex(k0^2-kappa2,0))+k0;
    fr = exp(-1i*0.5*ds*kappa2./denom);
    screen = exp(-alpha_xy*ds);
    psi0 = fftshift(ifft2(piston));
    psi0 = psi0/max(abs(psi0(iy0,ix0)),eps) * ...
        (cfg.frequency_axis_hz(ifq)/cfg.fc_hz);
    psi_k = fft2(psi0);
    for jj = 1:nstep
        psi_k = fr.*fft2(screen.*ifft2(fr.*psi_k));
    end
    incident(:,:,ifq) = ifft2(psi_k);
    q_k = fft2(receiver_source);
    fr_h = conj(fr); screen_h = conj(screen);
    for jj = nstep:-1:1 %#ok<NASGU>
        q_k = fr_h.*fft2(screen_h.*ifft2(fr_h.*q_k));
    end
    q(:,:,ifq) = ifft2(q_k);
end
ix0_pm = floor((cfg.pm_nx-cfg.pe_nx)/2)+1;
iy0_pm = floor((cfg.pm_ny-cfg.pe_ny)/2)+1;
adjoint_error = local_adjoint_audit(cfg,kappa2,alpha_xy,nstep,ds,ix0,iy0,k0_list(ceil(F/2)));
cache = struct('x_m',x,'y_m',y,'KX',KX,'KY',KY,'K',K, ...
    'incident_surface_xy_f',incident,'q_surface_xy_f',q, ...
    'k0_rad_per_m',k0_list,'idx_f_ref',ceil(F/2), ...
    'pm_ix',ix0_pm:(ix0_pm+cfg.pe_nx-1), ...
    'pm_iy',iy0_pm:(iy0_pm+cfg.pe_ny-1), ...
    'audit',struct('nstep_each_segment',nstep,'actual_dz_m',ds, ...
        'adjoint_inner_product_relative_error',adjoint_error, ...
        'source_directivity_min',min(piston(:)), ...
        'source_directivity_max',max(piston(:)), ...
        'pm_to_pe_rule','same-dx central crop; no interpolation; no de-meaning'));
end

function error_rel = local_adjoint_audit(cfg,kappa2,alpha_xy,nstep,ds,ix0,iy0,k0)
denom = sqrt(complex(k0^2-kappa2,0))+k0;
fr = exp(-1i*0.5*ds*kappa2./denom); screen = exp(-alpha_xy*ds);
rng(91901,'twister');
u = randn(cfg.pe_ny,cfg.pe_nx)+1i*randn(cfg.pe_ny,cfg.pe_nx);
v = complex(zeros(cfg.pe_ny,cfg.pe_nx)); v(iy0,ix0)=1;
Au = local_march_uniform(u,fr,screen,nstep);
AHv = local_march_adjoint(v,fr,screen,nstep);
lhs = sum(conj(v(:)).*Au(:)); rhs = sum(conj(AHv(:)).*u(:));
error_rel = abs(lhs-rhs)/max([abs(lhs),abs(rhs),eps]);
end

function out = local_march_uniform(in,fr,screen,nstep)
field_k = fft2(in);
for jj = 1:nstep
    field_k = fr.*fft2(screen.*ifft2(fr.*field_k));
end
out = ifft2(field_k);
end

function out = local_march_adjoint(in,fr,screen,nstep)
field_k = fft2(in); fr_h=conj(fr); screen_h=conj(screen);
for jj = nstep:-1:1 %#ok<NASGU>
    field_k = fr_h.*fft2(screen_h.*ifft2(fr_h.*field_k));
end
out = ifft2(field_k);
end

function alpha_xy = local_absorption_profile(x,y,cfg)
Ls_x = cfg.sponge_ratio*cfg.pe_xw_m; Ls_y = cfg.sponge_ratio*cfg.pe_yw_m;
x_start = cfg.pe_xw_m/2-Ls_x; y_start = cfg.pe_yw_m/2-Ls_y;
alpha_x = zeros(size(x)); alpha_y = zeros(size(y));
mx = abs(x)>x_start; my = abs(y)>y_start;
alpha_x(mx) = cfg.alpha_max_np_per_m*((abs(x(mx))-x_start)/Ls_x).^2;
alpha_y(my) = cfg.alpha_max_np_per_m*((abs(y(my))-y_start)/Ls_y).^2;
[AX,AY] = meshgrid(alpha_x,alpha_y); alpha_xy = AX+AY;
end

function h_f = local_receiver_response(cache,eta_pe)
F = numel(cache.k0_rad_per_m); h_f = complex(zeros(F,1));
for ifq = 1:F
    delta_phi = 2*cache.k0_rad_per_m(ifq)*eta_pe;
    G_xy = -exp(1i*delta_phi);
    inc = cache.incident_surface_xy_f(:,:,ifq);
    q = cache.q_surface_xy_f(:,:,ifq);
    h_f(ifq) = sum(conj(q(:)).*(G_xy(:).*inc(:)));
end
end

function [waveform,metrics] = local_synthesize_waveform(h_f,cfg)
source_f = cfg.pulse_duration_s*sinc(cfg.frequency_offset_hz*cfg.pulse_duration_s);
Y = source_f(:).*h_f(:);
N = cfg.waveform_nfft; F = cfg.Nf; mid = N/2+1; half=(F-1)/2;
Ypad = complex(zeros(N,1)); Ypad((mid-half):(mid+half)) = Y;
y = fftshift(fft(ifftshift(Ypad)))/N;
t = ((-N/2):(N/2-1)).'/(N*cfg.frequency_spacing_hz);
env = abs(y);
mask = abs(t)<=cfg.waveform_window_s;
idx = find(mask); envw=env(mask); tw=t(mask);
[peak_value,iloc] = max(envw); ipeak=idx(iloc); peak_time=tw(iloc);
threshold = cfg.threshold_fraction*peak_value;
icross_local = find(envw>=threshold,1,'first');
if isempty(icross_local)
    first_time = NaN;
elseif icross_local==1
    first_time = tw(1);
else
    t1=tw(icross_local-1); t2=tw(icross_local);
    y1=envw(icross_local-1); y2=envw(icross_local);
    first_time=t1+(threshold-y1)*(t2-t1)/max(y2-y1,eps);
end
waveform = struct('relative_time_s',t,'complex_envelope',y, ...
    'envelope',env,'source_spectrum',source_f,'received_spectrum',Y);
metrics = struct('peak_time_relative_s',peak_time, ...
    'first_threshold_time_relative_s',first_time,'peak_envelope',peak_value, ...
    'threshold_envelope',threshold,'threshold_fraction',cfg.threshold_fraction, ...
    'peak_index',ipeak);
end

function absolute = local_absolute_metrics(relative,tx,tau0)
absolute = struct( ...
    'peak_time_s',tau0+relative.peak_time_relative_s-tx.peak_time_relative_s, ...
    'first_threshold_time_s',tau0+relative.first_threshold_time_relative_s- ...
        tx.first_threshold_time_relative_s, ...
    'relative_peak_shift_s',relative.peak_time_relative_s-tx.peak_time_relative_s, ...
    'relative_first_shift_s',relative.first_threshold_time_relative_s- ...
        tx.first_threshold_time_relative_s, ...
    'threshold_fraction',relative.threshold_fraction, ...
    'threshold_reference','fraction of the peak envelope of that waveform');
end

function check = local_reflection_checks(incident,eta,k0)
flat_ref = -incident;
flat_error = max(abs((-incident)-flat_ref),[],'all');
delta = 2*k0*eta; reflected = -incident.*exp(1i*delta);
amplitude_error = max(abs(abs(reflected)-abs(incident)),[],'all');
mask = abs(incident)>max(abs(incident(:)))*1e-12;
factor_measured = reflected(mask)./incident(mask);
phase_error = angle((-factor_measured).*exp(-1i*delta(mask)));
operator_error = max(abs(factor_measured+exp(1i*delta(mask))));
check = struct('flat_pressure_release_max_abs_error',flat_error, ...
    'local_amplitude_max_abs_error',amplitude_error, ...
    'phase_2keta_max_wrapped_error_rad',max(abs(phase_error)), ...
    'boundary_operator_max_abs_error',operator_error, ...
    'formula','psi_ref=-psi_inc*exp(1i*2*k*eta); no 4*k*eta branch');
end

function radial = local_radial_spectrum(eta,spec,nbin)
Feta = fft2(eta); Phi = abs(Feta).^2/(numel(eta)^2* ...
    spec.dkx_rad_per_m*spec.dky_rad_per_m);
K = spec.K_rad_per_m; edges=linspace(0,max(K(:)),nbin+1);
kc=zeros(nbin,1); sample=NaN(nbin,1); target=NaN(nbin,1); count=zeros(nbin,1);
for ib=1:nbin
    m=K>=edges(ib)&K<edges(ib+1);
    kc(ib)=0.5*(edges(ib)+edges(ib+1)); count(ib)=nnz(m);
    if any(m(:))
        sample(ib)=mean(Phi(m))*2*pi*kc(ib);
        target(ib)=mean(spec.Phi2D(m))*2*pi*kc(ib);
    end
end
radial=struct('K_rad_per_m',kc,'E1D_sample',sample, ...
    'E1D_target',target,'bin_count',count, ...
    'sample_variance_from_spectrum_m2',sum(Phi(:))* ...
        spec.dkx_rad_per_m*spec.dky_rad_per_m, ...
    'sample_variance_spatial_m2',mean(eta(:).^2));
end

function [tbl,fit] = local_summary_table(U10,peak,first,flat)
metric_names={'peak','first_threshold'}; arrays={peak,first};
flat_values=[flat.peak_time_s,flat.first_threshold_time_s];
rows=cell(0,1); empty_fit=struct('t0_s',NaN,'b_s',NaN,'cdf_rmse',NaN,'method','');
fit=repmat(empty_fit,numel(U10),2); ir=0;
for imetric=1:2
    for iu=1:numel(U10)
        x=arrays{imetric}(:,iu); stats=local_stats(x); ray=local_fit_shifted_rayleigh(x);
        fit(iu,imetric)=ray; ir=ir+1;
        rows{ir,1}=struct('U10_mps',U10(iu),'metric',metric_names{imetric}, ...
            'sample_count',numel(x),'mean_ms',1e3*stats.mean, ...
            'std_ms',1e3*stats.std,'q05_ms',1e3*stats.q05, ...
            'q25_ms',1e3*stats.q25,'median_ms',1e3*stats.q50, ...
            'q75_ms',1e3*stats.q75,'q95_ms',1e3*stats.q95, ...
            'minimum_ms',1e3*min(x),'skewness',stats.skewness, ...
            'mean_shift_from_flat_ms',1e3*(stats.mean-flat_values(imetric)), ...
            'rayleigh_t0_ms',1e3*ray.t0_s,'rayleigh_b_ms',1e3*ray.b_s, ...
            'rayleigh_cdf_rmse',ray.cdf_rmse);
    end
end
tbl=struct2table(vertcat(rows{:}));
end

function stats=local_stats(x)
x=x(isfinite(x)); mu=mean(x); sd=sqrt(mean((x-mu).^2));
if sd>0, sk=mean(((x-mu)/sd).^3); else, sk=0; end
q=local_quantile(x,[0.05,0.25,0.5,0.75,0.95]);
stats=struct('mean',mu,'std',sd,'q05',q(1),'q25',q(2), ...
    'q50',q(3),'q75',q(4),'q95',q(5),'skewness',sk);
end

function q=local_quantile(x,p)
x=sort(x(:)); n=numel(x); q=zeros(size(p));
for ii=1:numel(p)
    pos=1+(n-1)*p(ii); lo=floor(pos); hi=ceil(pos);
    q(ii)=x(lo)+(pos-lo)*(x(hi)-x(lo));
end
end

function fit=local_fit_shifted_rayleigh(x)
x=sort(x(:)); n=numel(x); Femp=((1:n).'-0.5)/n;
s=max(sqrt(mean((x-mean(x)).^2)),eps);
b0=max(s/sqrt((4-pi)/2),eps); t00=mean(x)-b0*sqrt(pi/2);
gap0=max(min(x)-t00,0.05*s+eps);
objective=@(p) local_rayleigh_objective(p,x,Femp);
options=optimset('Display','off','MaxIter',1000,'MaxFunEvals',2000);
p=fminsearch(objective,log([gap0,b0]),options);
t0=min(x)-exp(p(1)); b=exp(p(2));
Fmodel=local_rayleigh_cdf(x,t0,b);
fit=struct('t0_s',t0,'b_s',b,'cdf_rmse',sqrt(mean((Fmodel-Femp).^2)), ...
    'method','two-parameter empirical-CDF least squares with t0<min(sample)');
end

function value=local_rayleigh_objective(p,x,Femp)
t0=min(x)-exp(p(1)); b=exp(p(2)); F=local_rayleigh_cdf(x,t0,b);
value=mean((F-Femp).^2);
if ~isfinite(value), value=realmax; end
end

function F=local_rayleigh_cdf(x,t0,b)
u=max(x-t0,0); F=1-exp(-u.^2/(2*b^2));
end

function tbl=local_sample_convergence(U10,peak,first,counts)
rows=cell(0,1); ir=0; metric_names={'peak','first_threshold'}; arrays={peak,first};
for imetric=1:2
    for iu=1:numel(U10)
        previous=[];
        for in=1:numel(counts)
            n=counts(in); s=local_stats(arrays{imetric}(1:n,iu));
            if isempty(previous)
                dmean=NaN; dstd=NaN; dq=NaN;
            else
                dmean=1e3*abs(s.mean-previous.mean);
                dstd=abs(s.std-previous.std)/max(previous.std,eps);
                dq=1e3*max(abs([s.q05,s.q50,s.q95]- ...
                    [previous.q05,previous.q50,previous.q95]));
            end
            ir=ir+1; rows{ir,1}=struct('U10_mps',U10(iu), ...
                'metric',metric_names{imetric},'sample_count',n, ...
                'mean_ms',1e3*s.mean,'std_ms',1e3*s.std, ...
                'q05_ms',1e3*s.q05,'median_ms',1e3*s.q50,'q95_ms',1e3*s.q95, ...
                'delta_mean_from_previous_ms',dmean, ...
                'relative_std_change_from_previous',dstd, ...
                'max_q_change_from_previous_ms',dq, ...
                'converged_vs_previous',isfinite(dmean)&&dmean<=0.10&&dstd<=0.10&&dq<=0.20);
            previous=s;
        end
    end
end
tbl=struct2table(vertcat(rows{:}));
end

function cost=local_cost_estimate(cfg,cache,cache_s,mc_s)
nstep=cache.audit.nstep_each_segment; F=cfg.Nf; nxy=cfg.pe_nx*cfg.pe_ny;
cost=struct('measured_cache_s',cache_s,'measured_monte_carlo_s',mc_s, ...
    'pe_fft2_count_estimate',4*F*nstep+4*F, ...
    'cache_array_bytes_estimate',2*nxy*F*16, ...
    'monte_carlo_boundary_complex_exponential_count', ...
        numel(cfg.U10_list_mps)*numel(cfg.sea_seed_list)*F*nxy, ...
    'paper_2048_cache_bytes_same_Nf',2*2048^2*F*16, ...
    'paper_2048_fft_grid_ratio_vs_current',(2048^2)/nxy, ...
    'note',['FFT count covers one upward incident march and one exact adjoint ', ...
        'receiver projection per frequency; MC uses the cached projection.']);
end

function acceptance=local_acceptance(result)
peak_std=zeros(numel(result.U10_list_mps),1); first_std=peak_std;
for iu=1:numel(peak_std)
    peak_std(iu)=std(result.peak_time_s(:,iu),1);
    first_std(iu)=std(result.first_threshold_time_s(:,iu),1);
end
checks=result.reflection_checks;
phase_err=max([checks.phase_2keta_max_wrapped_error_rad]);
amp_err=max([checks.local_amplitude_max_abs_error]);
flat_err=max(abs([result.flat.tau0_error_peak_s,result.flat.tau0_error_first_s]));
acceptance=struct( ...
    'flat_time_within_one_band_resolution',flat_err<=result.config.physical_time_resolution_s, ...
    'phase_uses_2keta_machine_precision',phase_err<=1e-12, ...
    'local_amplitude_preserved',amp_err<=1e-12, ...
    'same_surface_all_frequencies_by_construction',true, ...
    'pm_capture_all_cases_above_95pct',all([result.surface_audit.spectrum_capture_ratio]>=0.95), ...
    'broadening_gate_relative_increase',0.05, ...
    'peak_distribution_broadens',peak_std(end)>1.05*peak_std(1), ...
    'first_distribution_broadens',first_std(end)>1.05*first_std(1), ...
    'any_travel_time_metric_broadens', ...
        peak_std(end)>1.05*peak_std(1) || first_std(end)>1.05*first_std(1), ...
    'peak_std_ms',1e3*peak_std,'first_std_ms',1e3*first_std, ...
    'interpretation',['Trend gates are reported, never forced. A failed widening gate ', ...
        'identifies a model/numerical difference for investigation.']);
end

function local_write_artifacts(result)
out=result.config.output_dir; if ~exist(out,'dir'), mkdir(out); end
label=result.config.run_label;
save(fullfile(out,[label,'_result.mat']),'result','-v7.3');
writetable(result.summary_table,fullfile(out,[label,'_travel_time_summary.csv']));
writetable(result.surface_table,fullfile(out,[label,'_surface_audit.csv']));
writetable(result.reflection_table,fullfile(out,[label,'_reflection_checks.csv']));
writetable(result.sample_convergence_table, ...
    fullfile(out,[label,'_sample_convergence.csv']));
end

function local_make_figures(result)
out=result.config.output_dir; if ~exist(out,'dir'), mkdir(out); end
label=result.config.run_label; colors=lines(numel(result.U10_list_mps));
f=figure('Visible','off','Color','w'); tiledlayout(1,numel(result.U10_list_mps));
for iu=1:numel(result.U10_list_mps)
    nexttile; r=result.radial_spectra{iu};
    loglog(r.K_rad_per_m,r.E1D_target,'k-','LineWidth',1.5); hold on
    loglog(r.K_rad_per_m,r.E1D_sample_mean,'Color',colors(iu,:),'LineWidth',1.2);
    grid on; xlabel('K (rad m^{-1})'); ylabel('W(K) (m^3)');
    title(sprintf('U_{10}=%g m/s',result.U10_list_mps(iu)));
    legend('target PM','sample radial mean','Location','best');
end
exportgraphics(f,fullfile(out,[label,'_pm_target_vs_radial_sample.png']),'Resolution',180); close(f)

f=figure('Visible','off','Color','w'); tiledlayout(1,numel(result.U10_list_mps));
xpm=(-result.config.pm_xw_m/2):result.config.dx_m: ...
    (result.config.pm_xw_m/2-result.config.dx_m);
ypm=(-result.config.pm_yw_m/2):result.config.dy_m: ...
    (result.config.pm_yw_m/2-result.config.dy_m);
for iu=1:numel(result.U10_list_mps)
    nexttile; imagesc(xpm,ypm,result.example_surface_pm{iu}); axis xy image; c=colorbar;
    ylabel(c,'eta (m)');
    title(sprintf('U_{10}=%g m/s, explicit eta',result.U10_list_mps(iu)));
    xlabel('x (m)'); ylabel('y (m)');
end
exportgraphics(f,fullfile(out,[label,'_example_surfaces.png']),'Resolution',180); close(f)

f=figure('Visible','off','Color','w');
t0=result.flat.tx_waveform.relative_time_s*1e3;
plot(t0,result.flat.tx_waveform.envelope/max(result.flat.tx_waveform.envelope),'k--','LineWidth',1.2); hold on
plot(t0,result.flat.waveform.envelope/max(result.flat.waveform.envelope),'b-','LineWidth',1.2);
xlim(1e3*[-result.config.waveform_window_s,result.config.waveform_window_s]); grid on
xlabel('Time relative to nominal propagation (ms)'); ylabel('Normalized envelope');
legend('transmit pulse','flat-surface receive','Location','best');
title(sprintf('Flat surface: tau_0=%.3f ms, peak error=%.3f ms', ...
    1e3*result.geometry.tau0_s,1e3*result.flat.tau0_error_peak_s));
exportgraphics(f,fullfile(out,[label,'_flat_tx_rx_waveform.png']),'Resolution',180); close(f)

f=figure('Visible','off','Color','w');
for iu=1:numel(result.U10_list_mps)
    w=result.example_waveform{iu}; env=w.envelope/max(w.envelope);
    plot(1e3*w.relative_time_s,env,'Color',colors(iu,:),'LineWidth',1.2); hold on
end
xlim(1e3*[-result.config.waveform_window_s,result.config.waveform_window_s]); grid on
xlabel('Time relative to nominal propagation (ms)'); ylabel('Normalized envelope');
legend(compose('U_{10}=%g m/s',result.U10_list_mps),'Location','best');
title('Example rough-surface received pulses');
exportgraphics(f,fullfile(out,[label,'_rough_received_waveforms.png']),'Resolution',180); close(f)

f=figure('Visible','off','Color','w'); tiledlayout(1,2);
arrays={result.peak_time_s,result.first_threshold_time_s}; names={'Peak time','First 20% threshold'};
for imetric=1:2
    nexttile; hold on
    for iu=1:numel(result.U10_list_mps)
        histogram(1e3*(arrays{imetric}(:,iu)-result.geometry.tau0_s),12, ...
            'Normalization','probability','DisplayStyle','stairs','LineWidth',1.5, ...
            'EdgeColor',colors(iu,:));
    end
    grid on; xlabel('Offset from tau_0 (ms)'); ylabel('Probability'); title(names{imetric});
    legend(compose('U_{10}=%g',result.U10_list_mps),'Location','best');
end
exportgraphics(f,fullfile(out,[label,'_travel_time_distributions.png']),'Resolution',180); close(f)

f=figure('Visible','off','Color','w'); tiledlayout(1,3);
nexttile; plot(result.U10_list_mps,1e3*std(result.peak_time_s,1,1),'o-','LineWidth',1.3); hold on
plot(result.U10_list_mps,1e3*std(result.first_threshold_time_s,1,1),'s-','LineWidth',1.3);
grid on; xlabel('U_{10} (m/s)'); ylabel('Standard deviation (ms)'); title('Distribution width');
legend('peak','first threshold','Location','best');
nexttile; plot(result.U10_list_mps,1e3*min(result.first_threshold_time_s,[],1),'o-','LineWidth',1.3);
grid on; xlabel('U_{10} (m/s)'); ylabel('Earliest time (ms)'); title('First arrival');
nexttile; plot(result.U10_list_mps,1e3*mean(result.peak_time_s,1),'o-','LineWidth',1.3);
grid on; xlabel('U_{10} (m/s)'); ylabel('Mean peak time (ms)'); title('Peak time');
exportgraphics(f,fullfile(out,[label,'_width_first_peak_vs_wind.png']),'Resolution',180); close(f)
end

function a=local_empty_surface_audit()
a=struct('U10_mps',NaN,'U19p5_mps',NaN,'ustar_mps',NaN,'z0_m',NaN, ...
    'sigma_eta_sample_mean_m',NaN,'sigma_eta_sample_std_m',NaN, ...
    'Hs_sample_mean_m',NaN,'Hs_sample_std_m',NaN, ...
    'sigma_eta_discrete_target_m',NaN,'Hs_discrete_target_m',NaN, ...
    'sigma_eta_infinite_target_m',NaN,'Hs_infinite_target_m',NaN, ...
    'discrete_spectrum_variance_m2',NaN,'spectrum_capture_ratio',NaN, ...
    'radial_support_capture_ratio',NaN,'mean_sample_variance_rel_error',NaN, ...
    'max_sample_variance_rel_error',NaN);
end

function c=local_empty_reflection_check()
c=struct('flat_pressure_release_max_abs_error',NaN, ...
    'local_amplitude_max_abs_error',NaN, ...
    'phase_2keta_max_wrapped_error_rad',NaN, ...
    'boundary_operator_max_abs_error',NaN,'formula','');
end
