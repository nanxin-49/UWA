function receipt = audit_pm_wind_cutoff_k0_case(case_dir)
%AUDIT_PM_WIND_CUTOFF_K0_CASE Independent saved Fourier/CSV/FFT audit, no solve.
q=load(fullfile(case_dir,'surface.mat'),'data'); d=q.data;
c=readtable(fullfile(case_dir,'coefficients.csv'));
s=readtable(fullfile(case_dir,'surface.csv'));
p=readtable(fullfile(case_dir,'spectra.csv'));
m=d.metrics; n=height(s); x=s.x_m; dk=2*pi/160;
assert(max(abs(table2array(c)-table2array(d.coefficients)),[],'all')<1e-12);
assert(max(abs(table2array(s)-table2array(d.samples)),[],'all')<1e-11);
assert(max(abs(table2array(p)-table2array(d.spectra)),[],'all')<1e-11);
assert(max(abs(c.k_rad_per_m-(1:m.N)'*dk))<1e-14);
expected=.0081./(2*c.k_rad_per_m.^3).*exp(-.74*9.81^2./(m.U_PM_mps^4*c.k_rad_per_m.^2));
assert(max(abs(c.Sk_m3./expected-1))<1e-11);
assert(max(abs(c.cos_depth_coeff_m./sqrt(c.Sk_m3*dk)-c.gaussian_cos))<1e-11);
assert(max(abs(c.sin_depth_coeff_m./sqrt(c.Sk_m3*dk)-c.gaussian_sin))<1e-11);
legacy=RandStream('mt19937ar','Seed',m.seed);
za=randn(legacy,12,1); zb=randn(legacy,12,1);
extension=RandStream('mt19937ar','Seed',m.seed+1000000);
za=[za;randn(extension,28,1)]; zb=[zb;randn(extension,28,1)];
assert(max(abs(c.gaussian_cos-za(1:m.N)))<1e-11);
assert(max(abs(c.gaussian_sin-zb(1:m.N)))<1e-11);
% Separate scalar loop reconstruction, including independent analytic taper.
ix=unique([1:31:n n/2+1]); xx=x(ix); raw=zeros(numel(ix),3);
for j=1:height(c)
    k=c.k_rad_per_m(j); a=-c.cos_depth_coeff_m(j); b=-c.sin_depth_coeff_m(j);
    raw(:,1)=raw(:,1)+a*cos(k*xx)+b*sin(k*xx);
    raw(:,2)=raw(:,2)-k*a*sin(k*xx)+k*b*cos(k*xx);
    raw(:,3)=raw(:,3)-k^2*(a*cos(k*xx)+b*sin(k*xx));
end
w=double(abs(xx)<=42); wp=zeros(size(xx)); wpp=wp;
for j=1:numel(xx)
    if abs(xx(j))>42 && abs(xx(j))<50
        t=(abs(xx(j))-42)/8;
        w(j)=1-(6*t^5-15*t^4+10*t^3);
        wp(j)=-30*t^2*(t-1)^2*sign(xx(j))/8;
        wpp(j)=(-120*t^3+180*t^2-60*t)/64;
    end
end
emb=[w.*raw(:,1),w.*raw(:,2)+wp.*raw(:,1),w.*raw(:,3)+2*wp.*raw(:,2)+wpp.*raw(:,1)];
reconstruction=max(abs([raw emb]-table2array(s(ix,2:7))),[],'all');
assert(reconstruction<1e-11);
raw_variance=sum(c.cos_depth_coeff_m.^2+c.sin_depth_coeff_m.^2)/2;
assert(abs(raw_variance-m.raw_rms_height_m^2)<1e-12);
assert(abs(sum(p.raw_bin_power_m2)-mean(s.raw_elevation_m.^2))<1e-12);
assert(abs(sum(p.embedded_bin_power_m2)-mean(s.embedded_elevation_m.^2))<1e-12);
assert(max(abs(p.raw_bin_power_m2(2:m.N+1)-(c.cos_depth_coeff_m.^2+c.sin_depth_coeff_m.^2)/2))<1e-12);
assert(sum(p.raw_bin_power_m2(m.N+2:end))<1e-20);
assert(max(abs(p.PSD_difference_m3-(p.embedded_bin_power_m2-p.raw_bin_power_m2)/dk))<1e-12);
probes=unique([0 1 m.N m.N+1 53 127]);
for kk=probes
    % x-origin adds only phase; squared Fourier amplitudes are invariant.
    value=abs(sum(s.embedded_elevation_m.*exp(-2i*pi*kk*(0:n-1)'/n))/n)^2;
    if kk>0, value=2*value; end
    assert(abs(value-p.embedded_bin_power_m2(kk+1))<1e-12);
end
metric_names={'raw_rms_height_m','raw_rms_slope','raw_rms_second_per_m', ...
    'embedded_rms_height_m','embedded_rms_slope','embedded_rms_second_per_m'};
for kk=1:6
    assert(abs(m.(metric_names{kk})-sqrt(mean(s{:,kk+1}.^2)))<1e-12);
end
assert(abs(m.embedded_max_abs_slope-max(abs(s.embedded_slope)))<1e-12);
assert(abs(m.raw_max_abs_slope-max(abs(s.raw_slope)))<1e-12);
A=.74*9.81^2/m.U_PM_mps^4;
density=@(k).0081./(2*k.^3).*exp(-A./k.^2);
numerical_inside=integral(density,0,m.Kmax_rad_per_m,'AbsTol',1e-13,'RelTol',1e-10);
numerical_outside=integral(density,m.Kmax_rad_per_m,Inf,'AbsTol',1e-13,'RelTol',1e-10);
assert(abs(numerical_inside/m.variance_infinite_m2-m.continuous_coverage_0_Kmax)<1e-9);
assert(abs(numerical_outside-m.variance_above_Kmax_m2)<1e-10);
assert(all(s.embedded_elevation_m(abs(x)>=50)==0));
assert(isequal(s.raw_elevation_m(abs(x)<=42),s.embedded_elevation_m(abs(x)<=42)));
% Sampling audit changes dx only; physical coefficients/cutoff stay frozen.
coarse=s(1:2:end,:);
cols={'raw_elevation_m','raw_slope','raw_second_per_m','embedded_elevation_m','embedded_slope','embedded_second_per_m'};
grid_rms_relative=zeros(1,6);
for j=1:6
    fine_rms=sqrt(mean(s.(cols{j}).^2)); coarse_rms=sqrt(mean(coarse.(cols{j}).^2));
    grid_rms_relative(j)=abs(fine_rms-coarse_rms)/max(fine_rms,eps);
end
grid_max_slope_absolute=abs(max(abs(coarse.embedded_slope))-m.embedded_max_abs_slope);
assert(max(grid_rms_relative)<1e-4 && grid_max_slope_absolute<1e-4,'Spatial sampling audit failed');
receipt=struct('status','AUDITED_PASS','failed_checks',0,'reconstruction_max_abs',reconstruction, ...
    'grid_rms_relative',grid_rms_relative,'grid_max_slope_absolute',grid_max_slope_absolute, ...
    'physical_cutoff_fixed',true,'reference_grid_is_sampling_only',true);
pm_sensitivity_io('write',fullfile(case_dir,'audit.json'),receipt);
end
