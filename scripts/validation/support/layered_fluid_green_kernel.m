function [g,meta] = layered_fluid_green_kernel(cfg,kappa,z)
%LAYERED_FLUID_GREEN_KERNEL Exact per-layer 1-D outgoing Green kernel.
% exp(-i omega t), g''+(k(z)^2-kappa^2)g=-delta(z-zs).
% Constant water density; pressure and p_z/rho continuous at fluid interfaces.
% Bounded reflection recursion uses only decaying propagation factors.
arguments
    cfg (1,1) struct
    kappa (1,:) double
    z (:,1) double {mustBeFinite}
end
e=cfg.layer_edges_m(:); c=cfg.layer_c_mps(:); n=numel(c);
assert(numel(e)==n+1 && all(diff(e)>0) && all(c>0));
assert(all(z>=e(1) & z<=e(end)) && cfg.source_depth_m>e(1) && cfg.source_depth_m<e(end));
assert(all(real(kappa)>=0) && all(imag(kappa)<=1e-14),'Use lower-half-plane contour.');
w=2*pi*cfg.frequency_hz; h=diff(e); kw=w./c;
q=sqrt(complex(kw.^2-kappa.^2)); q(imag(q)<0)=-q(imag(q)<0);
assert(all(abs(q)>1e-14,'all'),'Contour crosses vertical branch point.');
E=exp(2i*q.*h); down_bottom=complex(zeros(size(q))); down_top=down_bottom;
up_top=down_bottom; up_bottom=down_bottom;
if strcmp(cfg.bottom_type,'outgoing')
    down_bottom(n,:)=0;
else
    attenuation=cfg.bottom_atten_db_wavelength*log(10)/20*cfg.frequency_hz/cfg.bottom_cp_mps;
    kb=w/cfg.bottom_cp_mps+1i*attenuation;
    qb=sqrt(complex(kb^2-kappa.^2)); qb(imag(qb)<0)=-qb(imag(qb)<0);
    yw=q(n,:)/cfg.water_rho_gcm3; yb=qb/cfg.bottom_rho_gcm3;
    down_bottom(n,:)=(yw-yb)./(yw+yb);
end
down_top(n,:)=down_bottom(n,:).*E(n,:);
for j=n-1:-1:1
    r=(q(j,:)-q(j+1,:))./(q(j,:)+q(j+1,:));
    down_bottom(j,:)=(r+down_top(j+1,:))./(1+r.*down_top(j+1,:));
    down_top(j,:)=down_bottom(j,:).*E(j,:);
end
if strcmp(cfg.top_type,'pressure_release'), up_top(1,:)=-1;
else, assert(strcmp(cfg.top_type,'outgoing')); end
up_bottom(1,:)=up_top(1,:).*E(1,:);
for j=2:n
    r=(q(j,:)-q(j-1,:))./(q(j,:)+q(j-1,:));
    up_top(j,:)=(r+up_bottom(j-1,:))./(1+r.*up_bottom(j-1,:));
    up_bottom(j,:)=up_top(j,:).*E(j,:);
end
zs=cfg.source_depth_m; s=find(zs>=e(1:end-1) & zs<e(2:end),1,'last');
qs=q(s,:); ht=zs-e(s); hb=e(s+1)-zs;
ru=up_top(s,:).*exp(2i*qs*ht); rd=down_bottom(s,:).*exp(2i*qs*hb);
den=1-ru.*rd; emission=1i./(2*qs);
U=emission.*(1+rd)./den; B=emission.*(1+ru)./den;
g=complex(zeros(numel(z),numel(kappa))); iz=discretize(z,e);
rows=find(iz==s & z<=zs); delta=z(rows)-e(s);
g(rows,:)=U.*(exp(1i*(zs-z(rows)).*qs)+up_top(s,:).*exp(1i*(delta+ht).*qs));
rows=find(iz==s & z>zs); delta=e(s+1)-z(rows);
g(rows,:)=B.*(exp(1i*(z(rows)-zs).*qs)+down_bottom(s,:).*exp(1i*(delta+hb).*qs));
A=B.*exp(1i*qs*hb);
for j=s+1:n
    r=(q(j-1,:)-q(j,:))./(q(j-1,:)+q(j,:));
    A=A.*(1+r)./(1+r.*down_top(j,:));
    rows=find(iz==j); dz=z(rows)-e(j);
    g(rows,:)=A.*(exp(1i*dz.*q(j,:))+down_bottom(j,:).*exp(1i*(2*h(j)-dz).*q(j,:)));
    A=A.*exp(1i*q(j,:)*h(j));
end
A=U.*exp(1i*qs*ht);
for j=s-1:-1:1
    r=(q(j+1,:)-q(j,:))./(q(j+1,:)+q(j,:));
    A=A.*(1+r)./(1+r.*up_bottom(j,:));
    rows=find(iz==j); dz=e(j+1)-z(rows);
    g(rows,:)=A.*(exp(1i*dz.*q(j,:))+up_top(j,:).*exp(1i*(2*h(j)-dz).*q(j,:)));
    A=A.*exp(1i*q(j,:)*h(j));
end
assert(all(isfinite(g),'all'),'Nonfinite layered kernel.');
meta=struct('source_layer',s,'source_layer_c_mps',c(s), ...
    'minimum_source_denominator',min(abs(den)), ...
    'max_propagation_factor',max(abs(E),[],'all'), ...
    'representation','exact piecewise-constant layered Green; density constant in water', ...
    'time_convention','exp(-i omega t)');
end
