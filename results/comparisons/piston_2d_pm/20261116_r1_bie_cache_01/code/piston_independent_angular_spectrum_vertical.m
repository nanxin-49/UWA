function [p,px,pz] = piston_independent_angular_spectrum_vertical(x,z,s)
% Independent continuous spectrum including propagating and evanescent parts.
x=x(:); z=z(:)+zeros(size(x)); p=complex(zeros(size(x))); px=p; pz=p;
k=2*pi*s.frequency_hz/s.c_mps; a=k*s.total_width_m/2;
C=s.rho_kgpm3*2*pi*s.frequency_hz*s.v0_mps*s.total_width_m/pi;
for j=1:numel(x)
    d=s.z_source_m-z(j); assert(d>0); tmax=asinh(45/(k*d));
    fun=@(b)local_sinc(a*sin(b)).*exp(1i*k*d*cos(b));
    ev=@(t)local_sinc(a*cosh(t)).*exp(-k*d*sinh(t));
    opts={'RelTol',1e-11,'AbsTol',1e-12};
    pr=integral(@(b)fun(b).*cos(k*x(j)*sin(b)),0,pi/2,opts{:});
    pe=integral(@(t)ev(t).*cos(k*x(j)*cosh(t)),0,tmax,opts{:});
    p(j)=C*(pr-1i*pe);
    px(j)=C*(integral(@(b)-k*sin(b).*fun(b).*sin(k*x(j)*sin(b)),0,pi/2,opts{:}) ...
        +1i*integral(@(t)k*cosh(t).*ev(t).*sin(k*x(j)*cosh(t)),0,tmax,opts{:}));
    pz(j)=C*(integral(@(b)-1i*k*cos(b).*fun(b).*cos(k*x(j)*sin(b)),0,pi/2,opts{:}) ...
        -1i*integral(@(t)k*sinh(t).*ev(t).*cos(k*x(j)*cosh(t)),0,tmax,opts{:}));
end
end
function v=local_sinc(q)
v=ones(size(q)); ix=q~=0; v(ix)=sin(q(ix))./q(ix);
end
