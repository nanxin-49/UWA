function receipt = audit_piston_pm_r1(root)
%AUDIT_PISTON_PM_R1 Independent saved-field, kernel, source and CSV checks.
arguments
    root (1,:) char
end
dest=fullfile(root,'audit'); assert(~isfolder(dest)); mkdir(dest); t0=tic;
diary(fullfile(dest,'diary.log')); guard=onCleanup(@()diary('off')); %#ok<NASGU>
q=load(fullfile(root,'common.mat')); a=load(fullfile(root,'comparison.mat'));
b=load(fullfile(root,'jobs',a.reference_job,'data.mat')); s=q.cfg.source;
checks=struct; values=struct; k=2*pi*4000/1500;
T=readtable(fullfile(root,'receiver_fields_and_errors.csv'));
saved=reshape(T.real_Pa+1i*T.imag_Pa,411,3);
values.csv_complex_max_abs_Pa=max(abs(saved-a.U),[],'all');
checks.csv_complex=values.csv_complex_max_abs_Pa<1e-11;
names=["PS" "BIE" "FullFK"];
for m=1:3
    rows=(m-1)*411+(1:411); err=a.U(:,m)-a.U(:,2);
    assert(all(T.x_m(rows)==q.x) && all(T.z_m(rows)==q.z));
    assert(max(abs(T.residual_abs_Pa(rows)-abs(err)))<1e-11);
    assert(max(abs(T.amplitude_difference_Pa(rows)-(abs(a.U(:,m))-abs(a.U(:,2)))))<1e-11);
    assert(max(abs(T.residual_over_BIE_RMS(rows)-abs(err)/a.scale))<1e-11);
    phase=angle(a.U(:,m).*conj(a.U(:,2)));
    assert(max(abs(T.phase_difference_rad(rows(a.mask))-phase(a.mask)))<1e-11);
    assert(all(isnan(T.phase_difference_rad(rows(~a.mask)))));
    assert(abs(norm(err)/norm(a.U(:,2))-a.metrics.complex_relative_L2(m))<1e-12);
    assert(all(string(T.method(rows))==names(m)));
end
checks.points_phase_metrics=true;
% Independent adaptive angular spectrum at selected actual surface points.
ii=unique(round(linspace(1,b.boundary.node_count,7)));
xx=b.boundary.node_x_m(ii); zz=b.boundary.node_z_m(ii);
[pi0,px0,pz0]=piston_independent_angular_spectrum_vertical(xx,zz,s);
[p,px,pz]=baffled_strip_piston_field_vertical(xx,zz,s,128);
values.source_pressure_relative=norm(p-pi0)/norm(pi0);
values.source_gradient_relative=norm([px-px0;pz-pz0])/norm([px0;pz0]);
checks.source=values.source_pressure_relative<1e-7 && values.source_gradient_relative<1e-7;
inc=baffled_strip_piston_field_vertical(b.boundary.node_x_m,b.boundary.node_z_m,s,128) ...
    -baffled_strip_piston_field_vertical(b.boundary.node_x_m,-b.boundary.node_z_m,s,128);
values.rhs_relative=norm(inc-b.incident)/norm(inc); checks.rhs=values.rhs_relative<1e-11;
% Reconstruct GL8 arc/window weights and half-plane Green field independently.
bn=b.boundary; [gx,gw]=local_gauss(8); edges=linspace(-bn.half_width_m,bn.half_width_m,bn.panel_count+1)';
half=bn.half_width_m/bn.panel_count; centers=(edges(1:end-1)+edges(2:end))/2;
xs=reshape((centers+half*gx.').',[],1); dxw=repmat(half*gw,bn.panel_count,1);
[zs,sl]=closeout_s3a_pm_surface(q.profile,xs); jac=hypot(1,sl); nx=-sl./jac; nz=1./jac;
win=local_window(xs,bn.half_width_m,.85); weights=dxw.*jac.*win;
checks.geometry=max(abs(xs-bn.node_x_m))<1e-12 && max(abs(zs-bn.node_z_m))<1e-12 ...
    && max(abs([nx nz]-bn.node_normal),[],'all')<1e-12;
reconstructed=complex(zeros(411,1));
for first=1:32:411
    jj=first:min(first+31,411); kernel=local_kernel(q.x(jj),q.z(jj),xs,zs,nx,nz,k,-2);
    reconstructed(jj)=kernel*(weights.*bn.density);
end
reconstructed=reconstructed+q.flat;
values.BIE_density_reconstruction_relative=norm(reconstructed-a.U(:,2))/norm(a.U(:,2));
checks.BIE_density=values.BIE_density_reconstruction_relative<1e-11;
% Independent selected off-grid boundary residual using GL112 and explicit
% polynomial interpolation. This does not solve for density again.
mid=(edges(1:end-1)+edges(2:end))/2; mid=mid(abs(mid)<=.85*bn.half_width_m);
targets=unique(round(linspace(1,numel(mid),19))); xt=mid(targets); zt=closeout_s3a_pm_surface(q.profile,xt);
[g,w]=local_gauss(112); layer=complex(zeros(numel(xt),1)); density0=layer;
for it=1:numel(xt)
    panel=find(xt(it)>=edges(1:end-1) & xt(it)<=edges(2:end),1);
    density0(it)=local_basis(gx,0)*bn.density((panel-1)*8+(1:8));
    for j=1:bn.panel_count
        if j==panel
            v=(g+1)/2; xi=[-v.^2;v.^2]; qw=[w.*v;w.*v];
        else
            xi=g; qw=w;
        end
        sx=centers(j)+half*xi; [sz,sp]=closeout_s3a_pm_surface(q.profile,sx); sj=hypot(1,sp);
        ker=local_kernel(xt(it),zt(it),sx,sz,-sp./sj,1./sj,k,-2);
        den=local_basis(gx,xi)*bn.density((j-1)*8+(1:8));
        layer(it)=layer(it)+ker*(half*qw.*sj.*local_window(sx,bn.half_width_m,.85).*den);
    end
end
rhs=baffled_strip_piston_field_vertical(xt,zt,s,128)-baffled_strip_piston_field_vertical(xt,-zt,s,128);
values.independent_boundary_relative=norm(rhs+.5*density0+layer)/norm(rhs);
checks.independent_boundary=values.independent_boundary_relative<1e-4;
% Full two-term FK reconstruction with saved geometry and independent source.
f=q.FK{2}; [ui,ux,uz]=baffled_strip_piston_field_vertical(f.surface_x_m,f.surface_z_m,s,128);
[u0,~,uz0]=baffled_strip_piston_field_vertical(f.surface_x_m,0,s,128);
dni=f.surface_normal(:,1).*ux+f.surface_normal(:,2).*uz;
rf=q.flat; term1=rf/2; term2=rf/2;
for first=1:32:411
    jj=first:min(first+31,411); dx=q.x(jj)-f.surface_x_m.'; dz=q.z(jj)-f.surface_z_m.'; rr=hypot(dx,dz);
    G=1i/4*besselh(0,1,k*rr); Gn=1i*k/4*besselh(1,1,k*rr).*(dx.*f.surface_normal(:,1).'+dz.*f.surface_normal(:,2).')./rr;
    rr0=hypot(dx,q.z(jj)); G0=1i/4*besselh(0,1,k*rr0); Gn0=1i*k/4*besselh(1,1,k*rr0).*q.z(jj)./rr0;
    term1(jj)=term1(jj)-sum((G.*(dni.*f.arc_jacobian).'-G0.*uz0.').*f.quadrature_dx_weights.',2);
    term2(jj)=term2(jj)-sum((Gn.*(ui.*f.arc_jacobian).'-Gn0.*u0.').*f.quadrature_dx_weights.',2);
end
values.FullFK_reconstruction_relative=norm(term1+term2-a.U(:,3))/norm(a.U(:,3));
checks.FullFK=values.FullFK_reconstruction_relative<1e-11 && norm(term1-f.term_incident_dn)/norm(term1)<1e-11 ...
    && norm(term2-f.term_green_dn)/norm(term2)<1e-11;
% Independent direct Green normal-derivative propagation of the compact PS
% screen: (i k/2) H1(kR) z/R. This checks Fourier origin, carrier and sign.
screen=q.PS{6}; xp=screen.surface_x_m; delta=screen.delta_Pa; dps=q.flat;
for first=1:32:411
    jj=first:min(first+31,411); rr=hypot(q.x(jj)-xp.',3);
    dps(jj)=dps(jj)+(.5i*k*besselh(1,1,k*rr).*3./rr)*(delta*screen.dx_m);
end
values.PS_independent_Green_relative=norm(dps-a.U(:,1))/norm(a.U(:,1));
checks.PS_green=values.PS_independent_Green_relative<1e-4;
checks.executed_model_numerical_gates=all(q.ps_values(1:2)<q.cfg.ps_gate) ...
    && q.ps_values(3)<q.cfg.ps_step_gate && all(q.fk_values<q.cfg.fk_gate);
unrun_checks={};
if a.reference_qualified
    checks.BIE_convergence=all(a.bie_values<q.cfg.bie_gate);
else
    unrun_checks={'BIE grid convergence','BIE self-quadrature convergence','BIE window convergence','reference-qualified model resolution'};
end
save(fullfile(dest,'reconstruction.mat'),'reconstructed','term1','term2','dps','xt','zt','rhs','layer','density0','-v7');
receipt=struct('status','INDEPENDENT_AUDIT_COMPLETE','failed_checks',nnz(~cell2mat(struct2cell(checks))), ...
    'elapsed_s',toc(t0),'checks',checks,'metrics',values,'boundary_audit_points',numel(xt), ...
    'unrun_checks',{unrun_checks},'reference_qualified',a.reference_qualified);
if ~a.reference_qualified, receipt.status='EXECUTED_SAVED_DATA_CHECKS_ONLY_REFERENCE_UNQUALIFIED'; end
fid=fopen(fullfile(dest,'receipt.json'),'w'); fprintf(fid,'%s',jsonencode(receipt,PrettyPrint=true)); fclose(fid);
disp(jsonencode(receipt)); assert(receipt.failed_checks==0,'Independent R1 audit failed');
end
function [g,w]=local_gauss(n)
j=(1:n-1)'; b=j./sqrt(4*j.^2-1); [V,D]=eig(diag(b,1)+diag(b,-1)); [g,ix]=sort(diag(D)); w=2*(V(1,ix).^2)';
end
function B=local_basis(g,x)
x=x(:); B=ones(numel(x),numel(g));
for j=1:numel(g), for l=1:numel(g), if l~=j, B(:,j)=B(:,j).*(x-g(l))/(g(j)-g(l)); end, end, end
end
function w=local_window(x,L,c)
r=abs(x)/L; w=ones(size(r)); w(r>=1)=0; ii=r>c & r<1; t=(r(ii)-c)/(1-c);
left=exp(-1./t); right=exp(-1./(1-t)); w(ii)=right./(left+right);
end
function K=local_kernel(xt,zt,xs,zs,nx,nz,k,h)
dx=xt(:)-xs(:).'; dz=zt(:)-zs(:).'; rr=hypot(dx,dz); dzI=zt(:)-(2*h-zs(:).'); rI=hypot(dx,dzI);
S=1i/4*(besselh(0,1,k*rr)-besselh(0,1,k*rI));
D=1i*k/4*(besselh(1,1,k*rr).*(dx.*nx(:).'+dz.*nz(:).')./rr ...
    -besselh(1,1,k*rI).*(dx.*nx(:).'-dzI.*nz(:).')./rI);
K=D-1i*k*S;
end
