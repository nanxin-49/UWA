function [p_f,px_f,pz_f,meta] = baffled_strip_piston_field_vertical(x_m,z_m,source,nquad)
%BAFFLED_STRIP_PISTON_FIELD_VERTICAL Physical 2-D uniform velocity aperture.
% e^(-i omega t), z down, source radiates upward; all pressures are peak Pa.
% No production defaults or Gaussian implementation are changed.
arguments
    x_m double
    z_m double
    source (1,1) struct
    nquad (1,1) double {mustBeInteger,mustBePositive} = 64
end
fields={'frequency_hz','c_mps','rho_kgpm3','v0_mps','total_width_m','z_source_m'};
assert(all(isfield(source,fields)) && nquad>=4);
assert(source.frequency_hz>0 && source.c_mps>0 && source.rho_kgpm3>0 && source.total_width_m>0);
if isscalar(z_m), z_m=z_m+zeros(size(x_m)); end
assert(isequal(size(x_m),size(z_m)) && all(isfinite(x_m(:))) && all(isfinite(z_m(:))));
s_m=source.z_source_m-z_m; assert(all(s_m(:)>0),'Points must lie above the physical source baffle.');
j=(1:nquad-1)'; beta=j./sqrt(4*j.^2-1);
[Q,L]=eig(diag(beta,1)+diag(beta,-1)); [g,ix]=sort(diag(L));
weights=source.total_width_m*(Q(1,ix).^2); xi=g*source.total_width_m/2;
k=2*pi*source.frequency_hz/source.c_mps;
scale=source.rho_kgpm3*2*pi*source.frequency_hz*source.v0_mps/2;
p_f=complex(zeros(size(x_m))); px_f=p_f; pz_f=p_f;
for q=1:nquad
    dx=x_m-xi(q); R=hypot(dx,s_m); h0=besselh(0,1,k*R);
    p_f=p_f+weights(q)*h0;
    if nargout>=2
        h1=besselh(1,1,k*R);
        px_f=px_f-weights(q)*k*h1.*dx./R;
        pz_f=pz_f+weights(q)*k*h1.*s_m./R;
    end
end
p_f=p_f*scale; px_f=px_f*scale; pz_f=pz_f*scale;
meta=struct('kind','2D infinite strip in rigid baffle','units','Pa; Pa/m gradients', ...
    'time_convention','exp(-i omega t)','velocity_direction',[0 -1], ...
    'phase_contract_version','2.0.0','representation','absolute_physical', ...
    'physical_source_z_m',source.z_source_m,'total_width_m',source.total_width_m, ...
    'quadrature_order',nquad,'velocity_peak_mps',source.v0_mps);
end
