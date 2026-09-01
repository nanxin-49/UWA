function [s_m,eta_m,meta] = sample_fixed_pm_profile_vertical(options)
%SAMPLE_FIXED_PM_PROFILE_VERTICAL Generate one deterministic 1-D PM profile.
% The Fourier coefficients are sampled once and the returned (s,eta) pair
% is the sole realization used by both native ATI and rotated-wall cases.
arguments
    options.wind_speed_mps (1,1) double {mustBePositive} = 6
    options.span_m (1,1) double {mustBePositive} = 160
    options.sample_count (1,1) double {mustBeInteger,mustBeGreaterThan(options.sample_count,2)} = 129
    options.seed (1,1) double {mustBeFinite} = 260001
    options.origin_m (1,1) double = NaN
end
N=options.sample_count;
L=options.span_m;
if mod(N,2)==0
    error('sample_count must be odd so s=0 is represented exactly.');
end
ds=L/N;
if isnan(options.origin_m)
    % Use one non-duplicated periodic interval and center the sample grid at
    % s=0 for every profile density.
    origin_m=-0.5*L*(N-1)/N;
else
    if ~isfinite(options.origin_m)
        error('origin_m must be finite or NaN (auto-centered).');
    end
    origin_m=options.origin_m;
end
s_m=origin_m+(0:N-1).' * ds;
if abs(s_m((N+1)/2))>1e-12
    error('Profile origin must place s=0 on a sample.');
end
k=2*pi/L*[0:(N-1)/2 -(N-1)/2:-1].';
g=9.81; alpha_PM=8.10e-3; beta_PM=0.74;
E=zeros(N,1);
mask=abs(k)>0;
E(mask)=alpha_PM./(2*abs(k(mask)).^3).*exp(-beta_PM*g^2./(options.wind_speed_mps^4*abs(k(mask)).^2));
rng(mod(round(options.seed),2^32),'twister');
coeff=sqrt(E*2*pi/L).*((randn(N,1)+1i*randn(N,1))/sqrt(2));
eta_m=sqrt(2)*N*real(ifft(coeff));
deta=sqrt(2)*N*real(ifft(1i*k.*coeff));
ddeta=sqrt(2)*N*real(ifft(-k.^2.*coeff));
meta=struct;
meta.model='raw_1d_pm_fourier_v1';
meta.wind_speed_mps=options.wind_speed_mps;
meta.span_m=L;
meta.sample_count=N;
meta.sample_spacing_m=ds;
meta.seed=mod(round(options.seed),2^32);
meta.origin_m=origin_m;
meta.alpha_PM=alpha_PM;
meta.beta_PM=beta_PM;
meta.g_mps2=g;
meta.mean_height_m=mean(eta_m);
meta.rms_height_m=sqrt(mean(eta_m.^2));
meta.rms_slope=sqrt(mean(deta.^2));
meta.max_slope=max(abs(deta));
meta.rms_curvature_per_m=sqrt(mean(ddeta.^2));
meta.max_curvature_per_m=max(abs(ddeta));
meta.curvature_radius_min_m=1/max(meta.max_curvature_per_m,realmin);
meta.height_samples_m=eta_m;
meta.slope_samples= deta;
meta.curvature_samples_per_m=ddeta;
meta.formula=['E1D(k)=alpha/(2*|k|^3)*exp(-beta*g^2/(U^4*|k|^2)); ', ...
    'eta=sqrt(2)*N*real(ifft(sqrt(E1D*dk).*CN(0,1)))'];
end
