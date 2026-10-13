function out=observe_deepsea3d_fg_residual(W,p,cfg,multiplier)
%OBSERVE_DEEPSEA3D_FG_RESIDUAL Independent oversampled physical operator.
% Analytic continuous C2 geometry is evaluated afresh. W alone is the
% retained Fourier/Chebyshev field. No projection before measuring residual.
timer=tic;n=p.n;m=multiplier*n;M=numel(p.s);
x=-p.L+(0:m-1)*2*p.L/m;[xx,yy]=meshgrid(x,x);
g=deepsea3d_tfe_geometry(xx,yy,cfg.profile);
Ws=W*p.D.';Wss=W*p.D2.';num=0;den=0;low=0;poly_num=0;
alias_num=0;retained_mapped=0;poly_retained_mapped=0;poly_full_mapped=0;
geometry_diff=0;mapped_geometry_diff=0;
% Also evaluate the discrete coefficient polynomial exactly at this grid,
% separating geometry interpolation from retained-product aliasing.
gp=struct();for key=["eta","ex","ey","lap"],gp.(key)=real(deepsea3d_fg_resize(p.geo.(key),m));end
gp.slope2=gp.ex.^2+gp.ey.^2;
for j=2:M-1
    F=fft2(reshape(W(:,j),n,n));Fs=fft2(reshape(Ws(:,j),n,n));
    H=deepsea3d_fg_resize(ifft2((p.k^2-p.kv.^2-p.kv.'.^2).*F),m);
    X=deepsea3d_fg_resize(ifft2(1i*p.kv.*Fs),m);Y=deepsea3d_fg_resize(ifft2(1i*p.kv.'.*Fs),m);
    S=deepsea3d_fg_resize(reshape(Ws(:,j),n,n),m);SS=deepsea3d_fg_resize(reshape(Wss(:,j),n,n),m);
    v=deepsea3d_fg_resize(reshape(W(:,j),n,n),m);r=1-p.s(j)/p.b;
    raw=local_operator(H,X,Y,S,SS,g,r,p.b);
    rawpoly=local_operator(H,X,Y,S,SS,gp,r,p.b);
    R=raw./(1-g.eta/p.b).^2;RP=rawpoly./(1-gp.eta/p.b).^2;
    qpoly=deepsea3d_fg_resize(rawpoly,n);qtrue=deepsea3d_fg_resize(raw,n);
    hb=ifft2((p.k^2-p.kv.^2-p.kv.'.^2).*F);
    xb=ifft2(1i*p.kv.*Fs);yb=ifft2(1i*p.kv.'.*Fs);
    colloc=local_operator(hb,xb,yb,reshape(Ws(:,j),n,n),reshape(Wss(:,j),n,n),p.geo,r,p.b);
    colloc=deepsea3d_fg_resize(colloc,n);
    alias_num=alias_num+norm(colloc-qpoly,'fro')^2*(m/n)^2;
    retained_mapped=retained_mapped+norm(qtrue,'fro')^2*(m/n)^2;
    poly_retained_mapped=poly_retained_mapped+norm(qpoly,'fro')^2*(m/n)^2;
    poly_full_mapped=poly_full_mapped+norm(rawpoly,'fro')^2;
    num=num+sum(abs(R).^2,'all');den=den+p.k^4*sum(abs(v).^2,'all');
    geometry_diff=geometry_diff+norm(R-RP,'fro')^2;
    mapped_geometry_diff=mapped_geometry_diff+norm(raw-rawpoly,'fro')^2;
    pr=deepsea3d_fg_resize(R,n);low=low+sum(abs(pr).^2,'all')*(m/n)^2;
    poly_num=poly_num+sum(abs(RP).^2,'all');
end
bottom=deepsea3d_fg_resize(reshape(W(:,1),n,n),m);
points=[xx(:)+cfg.center_m(1),yy(:)+cfg.center_m(2),g.eta(:)];
image=cfg.tx_m;image(3)=-image(3);
flat=reshape(deepsea3d_point_field(points,cfg.tx_m,p.k)-deepsea3d_point_field(points,image,p.k),m,m);
mask=hypot(xx,yy)<cfg.profile.support_m;
boundary=norm(flat(mask)+bottom(mask))/max(norm(flat(mask)),realmin);
top=reshape(W(:,end),n,n);
T=deepsea3d_fg_resize(ifft2(1i*p.gamma.*fft2(top)),m);
S=deepsea3d_fg_resize(reshape(W*p.D(end,:).',n,n),m);
dt=S./(1-g.eta/p.b)-T;
out=struct('oversampling',multiplier,'physical_PDE',sqrt(num/max(den,realmin)), ...
    'projected_physical_PDE',sqrt(low/max(den,realmin)), ...
    'coefficient_interpolant_PDE',sqrt(poly_num/max(den,realmin)), ...
    'collocation_alias_operator',sqrt(alias_num/max(den,realmin)), ...
    'retained_true_mapped_PDE',sqrt(retained_mapped/max(den,realmin)), ...
    'retained_polynomial_mapped_PDE',sqrt(poly_retained_mapped/max(den,realmin)), ...
    'full_polynomial_mapped_PDE',sqrt(poly_full_mapped/max(den,realmin)), ...
    'geometry_operator_difference',sqrt(geometry_diff/max(den,realmin)), ...
    'mapped_geometry_operator_difference',sqrt(mapped_geometry_diff/max(den,realmin)), ...
    'true_boundary',boundary,'true_DtN',norm(dt,'fro')/max(norm(T,'fro'),realmin), ...
    'wall_s',toc(timer),'scope','unprojected analytic geometry, same retained field; vertical Chebyshev nodes');
end
function R=local_operator(H,X,Y,S,SS,g,r,b)
e=g.eta/b;mix=g.ex.*X+g.ey.*Y;
R=(1-e).^2.*H+(1+g.slope2*r^2).*SS ...
    -2*(1-e)*r.*mix-((1-e).*g.lap+2*g.slope2/b)*r.*S;
end
