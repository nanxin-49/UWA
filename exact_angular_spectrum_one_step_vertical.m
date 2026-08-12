function [psi_end_xy,meta] = exact_angular_spectrum_one_step_vertical( ...
    psi_init_xy,xw_m,yw_m,f_hz,c0_mps,distance_m)
%EXACT_ANGULAR_SPECTRUM_ONE_STEP_VERTICAL Independent reduced-field reference.
% This validation helper deliberately performs one total-distance FFT step;
% it does not call or reuse the production PE marching implementation.

arguments
    psi_init_xy double
    xw_m (1,1) double {mustBeFinite,mustBePositive}
    yw_m (1,1) double {mustBeFinite,mustBePositive}
    f_hz (1,1) double {mustBeFinite,mustBePositive}
    c0_mps (1,1) double {mustBeFinite,mustBePositive}
    distance_m (1,1) double {mustBeFinite,mustBeNonnegative}
end
[ny,nx]=size(psi_init_xy);
if mod(nx,2)~=0 || mod(ny,2)~=0
    error('The production FFT layout requires even nx and ny.');
end
if any(~isfinite(psi_init_xy(:)))
    error('psi_init_xy must be finite.');
end
k=2*pi*f_hz/c0_mps;
kx=(2*pi/xw_m)*[0:(nx/2-1),-nx/2:-1];
ky=(2*pi/yw_m)*[0:(ny/2-1),-ny/2:-1];
[KX,KY]=meshgrid(kx,ky);
kz=sqrt(complex(k^2-KX.^2-KY.^2,0));
factor=exp(1i*distance_m*(kz-k));
psi_end_xy=ifft2(fft2(psi_init_xy).*factor);
meta=struct('operator','one-step exact discrete angular spectrum', ...
    'reduced_factor','exp(i*L*(kz-k))','evanescent_branch','imag(kz)>=0', ...
    'frequency_hz',f_hz,'distance_m',distance_m,'nx',nx,'ny',ny, ...
    'max_negative_imag_kz',max(max(-imag(kz))));
end
