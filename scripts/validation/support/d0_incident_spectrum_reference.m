function branches = d0_incident_spectrum_reference(x_m,z_m,source,order,phase_span)
%D0_INCIDENT_SPECTRUM_REFERENCE Segmented quadrature of qualified P0 integrands.
% Same Weyl integrands/normalization as piston_independent_angular_spectrum_vertical.
% Segmentation changes only integration, resolving deep-water oscillations.
arguments
    x_m (:,1) double
    z_m double
    source (1,1) struct
    order (1,1) double = 16
    phase_span (1,1) double = 8
end
z_m=z_m(:)+zeros(size(x_m)); d=source.z_source_m-z_m;
assert(all(d>0)); k=2*pi*source.frequency_hz/source.c_mps;
C=source.rho_kgpm3*2*pi*source.frequency_hz*source.v0_mps*source.total_width_m/pi;
j=(1:order-1)'; b=j./sqrt(4*j.^2-1);
[V,D]=eig(diag(b,1)+diag(b,-1)); [g,ix]=sort(diag(D)); w=2*(V(1,ix).^2)';
prop=complex(zeros(numel(x_m),3)); ev=prop;
segments=ceil(k*(max(d)+max(abs(x_m)))*(pi/2)/phase_span);
edges=linspace(0,pi/2,segments+1);
theta=reshape((g+1)/2.*diff(edges)+edges(1:end-1),[],1);
weights=reshape(w/2.*diff(edges),[],1);
qx=k*sin(theta); qz=k*cos(theta);
ap=local_sinc(qx*source.total_width_m/2);
for ii=1:numel(x_m)
    f=weights.*ap.*exp(1i*qz*d(ii));
    prop(ii,:)=C*[sum(f.*cos(qx*x_m(ii))), ...
        sum(-qx.*f.*sin(qx*x_m(ii))),sum(-1i*qz.*f.*cos(qx*x_m(ii)))];
    % sinh transform with exponent cutoff 45, matching P0. No evanescent deletion.
    top=asinh(45/(k*d(ii))); te=linspace(0,top,9);
    t=reshape((g+1)/2.*diff(te)+te(1:end-1),[],1);
    wt=reshape(w/2.*diff(te),[],1);
    qxe=k*cosh(t); alpha=k*sinh(t);
    f=wt.*local_sinc(qxe*source.total_width_m/2).*exp(-alpha*d(ii));
    ev(ii,:)=C*[sum(-1i*f.*cos(qxe*x_m(ii))), ...
        sum(1i*qxe.*f.*sin(qxe*x_m(ii))), ...
        sum(-1i*alpha.*f.*cos(qxe*x_m(ii)))];
end
branches=struct('propagating',prop,'evanescent',ev,'total',prop+ev, ...
    'columns',{{'p_Pa','px_Pa_per_m','pz_Pa_per_m'}},'segments',segments, ...
    'order',order,'phase_span',phase_span,'evanescent_segments',8,'tail_exponent',45, ...
    'pressure_tail_bound_Pa',C*exp(-45)./(k*d), ...
    'representation','absolute_physical','phase_contract_version','2.0.0');
end
function y=local_sinc(t)
y=ones(size(t)); ii=t~=0; y(ii)=sin(t(ii))./t(ii);
end
