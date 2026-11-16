function [eta,slope,curvature,raw] = closeout_s3a_pm_surface(profile,x)
%CLOSEOUT_S3A_PM_SURFACE Saved Fourier realization with frozen C2 embedding.
x=x(:); k=profile.k_rad_per_m(:).';
a=profile.cos_coeff_m(:); b=profile.sin_coeff_m(:);
phase=x*k;
raw=cos(phase)*a+sin(phase)*b;
hp=sin(phase)*(-k(:).*a)+cos(phase)*(k(:).*b);
hpp=cos(phase)*(-k(:).^2.*a)+sin(phase)*(-k(:).^2.*b);
r=abs(x); chi=double(r<=42); cp=zeros(size(x)); cpp=cp;
mid=r>42 & r<50; t=(r(mid)-42)/8;
chi(mid)=1-10*t.^3+15*t.^4-6*t.^5;
cp(mid)=(-30*t.^2+60*t.^3-30*t.^4).*sign(x(mid))/8;
cpp(mid)=(-60*t+180*t.^2-120*t.^3)/64;
eta=raw.*chi;
slope=hp.*chi+raw.*cp;
curvature=hpp.*chi+2*hp.*cp+raw.*cpp;
end
