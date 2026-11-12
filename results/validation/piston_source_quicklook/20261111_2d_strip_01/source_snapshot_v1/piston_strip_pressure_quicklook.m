function p_f = piston_strip_pressure_quicklook(x_m,s_m,width_m,f_hz,cfg,nquad)
%PISTON_STRIP_PRESSURE_QUICKLOOK Validation-only 2-D Neumann aperture integral.
% e^(-i*omega*t); s=100-z is upward distance; v0 is a peak velocity.
% p = rho*omega*v0/2 * integral_{-W/2}^{W/2} H0^(1)(k*R) dxi.
% Infinite y-invariant strip in an infinite rigid source baffle, s>0.
% No surface, bottom, PE, periodic window, or circular-piston approximation.
assert(isequal(size(x_m),size(s_m)) && all(s_m(:)>0));
assert(width_m>0 && f_hz>0 && nquad>=2);
beta=(1:nquad-1)./sqrt(4*(1:nquad-1).^2-1);
[Q,L]=eig(diag(beta,1)+diag(beta,-1));
[nodes,idx]=sort(diag(L));
weights=2*Q(1,idx).^2;
nodes=nodes*width_m/2; weights=weights*width_m/2;
p_f=complex(zeros(size(x_m)));
k=2*pi*f_hz/cfg.c_mps;
for jj=1:nquad
    R_m=hypot(x_m-nodes(jj),s_m);
    p_f=p_f+weights(jj)*besselh(0,1,k*R_m);
end
p_f=p_f*(cfg.rho_kgpm3*2*pi*f_hz*cfg.v0_mps/2);
end
