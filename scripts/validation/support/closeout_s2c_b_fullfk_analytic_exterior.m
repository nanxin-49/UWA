function r = closeout_s2c_b_fullfk_analytic_exterior(src,x,z,panels,extent_m,flat_exact)
%CLOSEOUT_S2C_B_FULLFK_ANALYTIC_EXTERIOR Full double-term K_HIGH integral.
% Same complete Kirchhoff algebra as the qualified S2b GL8 formulation.
arguments
    src (1,1) struct
    x (:,1) double
    z (:,1) double
    panels (1,1) double {mustBeInteger,mustBePositive}
    extent_m (1,1) double {mustBePositive}
    flat_exact (:,1) double
end
assert(isequal(size(x),size(z),size(flat_exact)) && extent_m>=50);
j=(1:7)'; beta=j./sqrt(4*j.^2-1); [v,d]=eig(diag(beta,1)+diag(beta,-1));
[gx,order]=sort(diag(d)); v=v(:,order); gw=2*(v(1,:).^2)';
edges=linspace(-extent_m,extent_m,panels+1)';
centers=(edges(1:end-1)+edges(2:end))/2; half=extent_m/panels;
xs=reshape((centers+half*gx.').',[],1); qw=repmat(half*gw,panels,1);
[eta,ep]=closeout_s2c_b_surface(xs,.05,.2,0);
jac=hypot(1,ep); nx=-ep./jac; nz=1./jac;
[ui,ux,uz]=fk_single_receiver_incident(src,xs,eta);
[u0,~,uz0]=fk_single_receiver_incident(src,xs,zeros(size(xs)));
dni=nx.*ux+nz.*uz; k=2*pi*4000/1500;
t1=complex(zeros(size(x))); t2=t1;
for first=1:32:numel(x)
    ii=first:min(first+31,numel(x)); dx=x(ii)-xs.'; dz=z(ii)-eta.'; rr=hypot(dx,dz);
    g=1i/4*besselh(0,1,k*rr);
    gn=1i*k/4*besselh(1,1,k*rr).*(dx.*nx.'+dz.*nz.')./rr;
    r0=hypot(dx,z(ii)); g0=1i/4*besselh(0,1,k*r0);
    gn0=1i*k/4*besselh(1,1,k*r0).*z(ii)./r0;
    t1(ii)=-sum((g.*(dni.*jac).'-g0.*uz0.').*qw.',2);
    t2(ii)=-sum((gn.*(ui.*jac).'-gn0.*u0.').*qw.',2);
end
r=struct('receiver_field',flat_exact+t1+t2, ...
    'term_incident_dn',flat_exact/2+t1,'term_green_dn',flat_exact/2+t2, ...
    'rough_flat_term1_difference',t1,'rough_flat_term2_difference',t2, ...
    'surface_x_m',xs,'surface_z_m',eta,'surface_normal',[nx nz], ...
    'arc_jacobian',jac,'quadrature_dx_weights',qw, ...
    'surface_incident',ui,'surface_incident_dn',dni, ...
    'flat_surface_incident',u0,'flat_surface_incident_dn',uz0, ...
    'panels',panels,'gauss_order',8,'extent_m',extent_m,'local_radius_m',Inf, ...
    'quadrature_measure','exact analytic arc_length', ...
    'formulation','complete double-term integral: analytic infinite flat exterior plus compact rough-minus-flat integrand', ...
    'version','FullFK_complete_analytic_flat_exterior_GL8_v1', ...
    'representation','native exp(-i omega t); outgoing exp(+ikR)');
assert(max(abs(r.receiver_field-r.term_incident_dn-r.term_green_dn))<1e-13);
end
