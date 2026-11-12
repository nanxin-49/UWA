function [p_f,meta] = piston_continuous_pe_march(x_m,z_m,source,n_modes,step_m,offset_m)
%PISTON_CONTINUOUS_PE_MARCH Continuous square-root PE modes for a strip.
% Integrate propagating modes in theta and evanescent modes in sinh(t).
% This uses the same exp(i*(kz-k)*ds) reduced PE march. The continuum source
% avoids periodic replicas of its slowly decaying grazing spatial tails.
arguments
    x_m (:,1) double
    z_m (1,1) double
    source (1,1) struct
    n_modes (1,1) double {mustBeInteger,mustBePositive}
    step_m (1,1) double {mustBePositive}
    offset_m (1,1) double {mustBePositive}
end
d=source.z_source_m-z_m; assert(d>=offset_m);
k=2*pi*source.frequency_hz/source.c_mps;
[g,w]=local_gauss(n_modes);
theta=(g+1)*pi/4; wtheta=w*pi/4;
tmax=asinh(45/(k*d)); t=(g+1)*tmax/2; wt=w*tmax/2;
kx=[k*sin(theta);k*cosh(t)]; kz=[k*cos(theta);1i*k*sinh(t)];
coeff=[wtheta.*local_sinc(kx(1:n_modes)*source.total_width_m/2); ...
    -1i*wt.*local_sinc(kx(n_modes+1:end)*source.total_width_m/2)];
coeff=coeff.*exp(1i*(kz-k)*offset_m);
nsteps=max(1,ceil((d-offset_m)/step_m)); ds=(d-offset_m)/nsteps;
half=exp(1i*(kz-k)*ds/2);
for j=1:nsteps
    coeff=half.*(half.*coeff);
end
p_f=complex(zeros(size(x_m)));
for first=1:128:numel(x_m)
    ii=first:min(first+127,numel(x_m));
    p_f(ii)=cos(x_m(ii)*kx.')*coeff;
end
p_f=p_f*(source.rho_kgpm3*2*pi*source.frequency_hz*source.v0_mps*source.total_width_m/pi)*exp(1i*k*d);
meta=struct('operator','same reduced square-root PE, continuous kx quadrature', ...
    'modes_per_branch',n_modes,'nsteps',nsteps,'step_m',ds,'offset_m',offset_m, ...
    'physical_source_z_m',source.z_source_m,'numerical_plane_z_m',source.z_source_m-offset_m, ...
    'evanescent_tail_exponent',45,'representation','absolute_physical Pa, phase contract v2', ...
    'source_support','full aperture spectrum, propagating and evanescent branches');
end
function [g,w]=local_gauss(n)
j=(1:n-1)'; b=j./sqrt(4*j.^2-1); [V,D]=eig(diag(b,1)+diag(b,-1));
[g,ix]=sort(diag(D)); w=2*(V(1,ix).^2)';
end
function v=local_sinc(q)
v=ones(size(q)); nz=q~=0; v(nz)=sin(q(nz))./q(nz);
end
