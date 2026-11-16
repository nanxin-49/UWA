function receipt=audit_r1_bie_manual_stage(stage_dir,common_path,recipe)
%AUDIT_R1_BIE_MANUAL_STAGE Independent saved-density/CSV/source audit.
% No BIE solve. All solver platform midpoints remain mandatory; the 19
% GL112 boundary reconstruction points are an additional independent check.
q=load(common_path); b=load(fullfile(stage_dir,'data.mat')); bn=b.boundary;
s=q.cfg.source; k=2*pi*4000/1500; n=8;
assert(bn.half_width_m==recipe(1) && bn.points_per_wavelength_requested==recipe(2) && ...
    bn.self_quadrature_order==recipe(3) && bn.panel_order==8 && bn.halfplane_h_m==-2);
expected_nodes=8*max(4,ceil(2*recipe(1)/((1500/4000)/recipe(2)*8)));
assert(bn.node_count==expected_nodes && bn.frequency_hz==4000 && bn.c0_mps==1500 && bn.window_plateau_ratio==.85);
assert(isfield(b.numeric_config,'enable_static_cache') && b.numeric_config.enable_static_cache);
assert(isequal(bn.receiver_x_m,q.x) && isequal(bn.receiver_z_m,q.z) && numel(b.field)==411);
assert(bn.finite && all(bn.linear_residual<1e-10) && all(bn.offgrid_boundary_residual<1e-4));
[gx,gw]=local_gauss(n); L=recipe(1); P=bn.panel_count;
edges=linspace(-L,L,P+1)'; centers=(edges(1:end-1)+edges(2:end))/2; half=diff(edges)/2;
xs=reshape((centers+half.*gx.').',[],1); dxw=reshape((half.*gw.').',[],1);
[zs,sp]=closeout_s3a_pm_surface(q.profile,xs); jac=hypot(1,sp); nx=-sp./jac; nz=1./jac;
assert(max(abs(xs-bn.node_x_m))<1e-12 && max(abs(zs-bn.node_z_m))<1e-12 && ...
    max(abs([nx nz]-bn.node_normal),[],'all')<1e-12);
weights=dxw.*jac.*local_window(xs,L,.85); field=q.flat;
for first=1:32:411
    ii=first:min(first+31,411);
    field(ii)=field(ii)+local_kernel(q.x(ii),q.z(ii),xs,zs,nx,nz,k,-2)*(weights.*bn.density);
end
T=readtable(fullfile(stage_dir,'points.csv'));
assert(isequal(T.x_m,q.x) && isequal(T.z_m,q.z));
csv=T.real_Pa+1i*T.imag_Pa;
values=struct('density_receiver_relative',norm(field-b.field)/norm(b.field), ...
    'csv_max_abs_Pa',max(abs(csv-b.field)),'component_closure_abs_Pa',max(abs(b.field-q.flat-bn.receiver_field)));
inc=baffled_strip_piston_field_vertical(xs,zs,s,128)-baffled_strip_piston_field_vertical(xs,-zs,s,128);
values.incident_rhs_relative=norm(inc-b.incident)/max(norm(inc),realmin);
mid=centers(abs(centers)<=.85*L); zmid=closeout_s3a_pm_surface(q.profile,mid);
drive=baffled_strip_piston_field_vertical(mid,zmid,s,128)-baffled_strip_piston_field_vertical(mid,-zmid,s,128);
assert(isequal(mid,bn.boundary_check_x_m),'Missing full-platform midpoint support');
values.saved_all_midpoint_residual=norm(bn.boundary_total_residual)/max(norm(drive),realmin);
ii=unique(round(linspace(1,numel(mid),19))); xt=mid(ii); zt=zmid(ii);
[g,w]=local_gauss(112); layer=complex(zeros(numel(xt),1)); density0=layer;
for it=1:numel(xt)
    pt=find(xt(it)>=edges(1:end-1) & xt(it)<=edges(2:end),1);
    density0(it)=local_basis(gx,0)*bn.density((pt-1)*n+(1:n));
    for ps=1:P
        if ps==pt
            v=(g+1)/2; xi=[-v.^2;v.^2]; qw=[w.*v;w.*v];
        else
            xi=g; qw=w;
        end
        x=centers(ps)+half(ps)*xi; [z,slope]=closeout_s3a_pm_surface(q.profile,x); j=hypot(1,slope);
        den=local_basis(gx,xi)*bn.density((ps-1)*n+(1:n));
        layer(it)=layer(it)+local_kernel(xt(it),zt(it),x,z,-slope./j,1./j,k,-2)* ...
            (half(ps)*qw.*j.*local_window(x,L,.85).*den);
    end
end
values.independent_19_boundary_relative=norm(drive(ii)+.5*density0+layer)/max(norm(drive(ii)),realmin);
checks=struct('density_receiver',values.density_receiver_relative<1e-11, ...
    'MAT_CSV',values.csv_max_abs_Pa<1e-11,'component_closure',values.component_closure_abs_Pa<1e-11, ...
    'incident_rhs',values.incident_rhs_relative<1e-11, ...
    'all_midpoints',values.saved_all_midpoint_residual<1e-4, ...
    'independent_boundary',values.independent_19_boundary_relative<1e-4);
save(fullfile(stage_dir,'audit_reconstruction.mat'),'field','xt','zt','layer','density0','-v7');
receipt=struct('status','PASS','failed_checks',nnz(~cell2mat(struct2cell(checks))), ...
    'checks',checks,'metrics',values,'platform_midpoints',numel(mid),'independent_boundary_points',numel(ii));
if receipt.failed_checks>0, receipt.status='FAIL'; end
r1_bie_manual_io('write',fullfile(stage_dir,'audit.json'),receipt);
assert(receipt.failed_checks==0,'r1manual:StageAudit','Independent stage audit failed');
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
a=exp(-1./t); b=exp(-1./(1-t)); w(ii)=b./(a+b);
end
function K=local_kernel(xt,zt,xs,zs,nx,nz,k,h)
dx=xt(:)-xs(:).'; dz=zt(:)-zs(:).'; rr=hypot(dx,dz); dzI=zt(:)-(2*h-zs(:).'); rI=hypot(dx,dzI);
S=1i/4*(besselh(0,1,k*rr)-besselh(0,1,k*rI));
D=1i*k/4*(besselh(1,1,k*rr).*(dx.*nx(:).'+dz.*nz(:).')./rr- ...
    besselh(1,1,k*rI).*(dx.*nx(:).'-dzI.*nz(:).')./rI);
K=D-1i*k*S;
end
