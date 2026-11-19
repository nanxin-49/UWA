function [z,zp,zpp,raw,rp,rpp] = d1b_pm_surface(p,x)
%D1B_PM_SURFACE Physical depth and analytic derivatives, positive down.
x = x(:); r = abs(x); theta = x*p.k.';
C = cos(theta); S = sin(theta);
raw = C*p.a+S*p.b;
rp = C*(p.k.*p.b)-S*(p.k.*p.a);
rpp = -C*(p.k.^2.*p.a)-S*(p.k.^2.*p.b);
width = p.support-p.plateau;
w = double(r<=p.plateau); wp = zeros(size(x)); wpp = wp;
mask = r>p.plateau & r<p.support; t = (r(mask)-p.plateau)/width;
w(mask) = (1-t).^3.*(1+3*t+6*t.^2);
wp(mask) = -30*t.^2.*(1-t).^2.*sign(x(mask))/width;
wpp(mask) = -60*t.*(1-t).*(1-2*t)/width^2;
z = w.*raw; zp = wp.*raw+w.*rp; zpp = wpp.*raw+2*wp.*rp+w.*rpp;
end
