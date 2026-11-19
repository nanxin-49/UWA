function [field,owner,setup_s,receiver_s] = d1b_ps_prototype(src,p,x,z,flat,workers)
%D1B_PS_PROTOTYPE D1 phase-screen correction with coarse frozen FFT sampling.
t0 = tic; W = 640; dx = .125; N = W/dx;
xx = -W/2+(0:N-1)'*dx;
k = 2*pi*src.frequency_hz/src.c_mps;
wall = d1b_pm_surface(p,xx);
inc = baffled_strip_piston_field_vertical(xx,0,src);
delta = (1-exp(-2i*k*wall)).*inc;
kx = (2*pi/W)*[0:N/2-1 -N/2:-1]';
kz = sqrt(complex(k^2-kx.^2,0));
assert(all(z==6));
spec = (fft(delta)/N).*exp(1i*kz*z(1));
setup_s = toc(t0);
t0 = tic;
[correction,owner] = d1b_tiles(numel(x),32,workers,@(ii) exp(1i*(x(ii)+W/2)*kx.')*spec);
field = flat+correction;
receiver_s = toc(t0);
end
