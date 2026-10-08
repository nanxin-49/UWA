function [eta,slope,second] = evaluate_pm_pilot_surface(profile,x,embed)
%EVALUATE_PM_PILOT_SURFACE Fixed Fourier PM with the frozen C2 embedding.
% Raw coefficients are unchanged on |x|<=42; the surface is flat beyond 50.
x=x(:); k=profile.k_rad_per_m(:).';
if nargin<3, embed=true; end
phase=x*k; a=profile.cos_coeff_m(:); b=profile.sin_coeff_m(:);
h=cos(phase)*a+sin(phase)*b;
hp=sin(phase)*(-k(:).*a)+cos(phase)*(k(:).*b);
hpp=cos(phase)*(-k(:).^2.*a)+sin(phase)*(-k(:).^2.*b);
if ~embed, eta=h; slope=hp; second=hpp; return; end
v=double(abs(x)<=42); vp=zeros(size(x)); vpp=vp;
mid=abs(x)>42 & abs(x)<50; t=(abs(x(mid))-42)/8;
v(mid)=1-10*t.^3+15*t.^4-6*t.^5;
vp(mid)=(-30*t.^2+60*t.^3-30*t.^4).*sign(x(mid))/8;
vpp(mid)=(-60*t+180*t.^2-120*t.^3)/64;
eta=h.*v; slope=hp.*v+h.*vp; second=hpp.*v+2*hp.*vp+h.*vpp;
end
