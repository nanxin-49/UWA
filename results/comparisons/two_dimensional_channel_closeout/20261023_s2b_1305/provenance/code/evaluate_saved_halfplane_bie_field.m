function u = evaluate_saved_halfplane_bie_field(boundary, x_m, z_m)
%EVALUATE_SAVED_HALFPLANE_BIE_FIELD Evaluate saved density; never solve BIE.
% Same off-surface D_h-i*k*S_h kernel as frozen halfplane BIE validation.
% Quadrature is reconstructed from its saved uniform panels and GL nodes.
assert(isequal(size(x_m),size(z_m)) && all(isfinite([x_m(:);z_m(:)])));
b=boundary; n=numel(b.node_x_m);
assert(n==b.panel_order*b.panel_count && all(b.node_normal(:,2)>0));
j=(1:b.panel_order-1)'; beta=j./sqrt(4*j.^2-1);
[v,d]=eig(diag(beta,1)+diag(beta,-1));
[~,order]=sort(diag(d)); v=v(:,order); gw=2*(v(1,:).^2)';
weights=(b.half_width_m/b.panel_count)*repmat(gw,b.panel_count,1) ...
    ./b.node_normal(:,2).*b.window(:);
strength=weights.*b.density(:); k=b.k_radpm; h=b.halfplane_h_m;
xs=b.node_x_m(:).'; zs=b.node_z_m(:).';
nx=b.node_normal(:,1).'; nz=b.node_normal(:,2).';
shape=size(x_m); xt=x_m(:); zt=z_m(:); value=complex(zeros(numel(xt),1));
assert(all(zt>max(zs)),'This evaluator is restricted to planes below the rough surface.');
for first=1:128:numel(xt)
    rows=first:min(first+127,numel(xt));
    dx=xt(rows)-xs; dz=zt(rows)-zs; r=hypot(dx,dz);
    assert(all(r>0,'all'));
    phi=1i/4*besselh(0,1,k*r);
    dp=1i*k/4*besselh(1,1,k*r).*((dx./r).*nx+(dz./r).*nz);
    dzi=zt(rows)-(2*h-zs); ri=hypot(dx,dzi);
    phii=1i/4*besselh(0,1,k*ri);
    dpi=1i*k/4*besselh(1,1,k*ri).*((dx./ri).*nx-(dzi./ri).*nz);
    value(rows)=((dp-dpi)-1i*k*(phi-phii))*strength;
end
u=reshape(value,shape); assert(all(isfinite(u),'all'));
end
