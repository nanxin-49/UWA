function [pressure,meta] = weyl_point_source_reference_vertical(rho_m,z_m,f_hz,c0_mps)
%WEYL_POINT_SOURCE_REFERENCE_VERTICAL Continuous free-space Weyl integral.
% The propagating substitution kz=sqrt(k^2-kappa^2) and evanescent
% substitution q=sqrt(kappa^2-k^2) remove the 1/kz grazing singularity:
% G = i/(4*pi) int_0^k J0(rho*sqrt(k^2-kz^2))*exp(i*kz*z) dkz
%   + 1/(4*pi) int_0^Inf J0(rho*sqrt(k^2+q^2))*exp(-q*z) dq.

arguments
    rho_m double {mustBeFinite,mustBeNonnegative}
    z_m (1,1) double {mustBeFinite,mustBePositive}
    f_hz (1,1) double {mustBeFinite,mustBePositive}
    c0_mps (1,1) double {mustBeFinite,mustBePositive}
end
k=2*pi*f_hz/c0_mps;
pressure=complex(zeros(size(rho_m)));
for ii=1:numel(rho_m)
    rho=rho_m(ii);
    prop=@(kz)besselj(0,rho*sqrt(max(k^2-kz.^2,0))).*exp(1i*kz*z_m);
    evan=@(q)besselj(0,rho*sqrt(k^2+q.^2)).*exp(-q*z_m);
    Iprop=integral(prop,0,k,'RelTol',1e-11,'AbsTol',1e-13);
    Ievan=integral(evan,0,Inf,'RelTol',1e-11,'AbsTol',1e-13);
    pressure(ii)=1i*Iprop/(4*pi)+Ievan/(4*pi);
end
R=hypot(z_m,rho_m);
analytic=exp(1i*k*R)./(4*pi*R);
meta=struct('representation','continuous Weyl integral with grazing-singularity substitutions', ...
    'time_convention','exp(-i*omega*t)','relative_error_to_green', ...
    norm(pressure(:)-analytic(:))/norm(analytic(:)), ...
    'max_tl_error_db',max(abs(20*log10(abs(pressure(:)./analytic(:))))), ...
    'max_phase_error_rad',max(abs(angle(pressure(:).*conj(analytic(:))))));
end
