function out=deepsea3d_dirichlet_fe(cfg,opts)
%DEEPSEA3D_DIRICHLET_FE Exploratory independent Dirichlet field expansion.
% Homotopy t*eta; delta field has outgoing half-space Fourier coefficients.
% b_n=-[t^n]u_flat(t eta)-sum_j eta^j/j! D_z^j b_(n-j).
% Retains all sampled evanescent modes. No filter, fitting or KA boundary.
arguments
    cfg (1,1) struct
    opts (1,1) struct
end
assert(cfg.profile.height_m==0,'deepsea3dFE:Plane','This version uses z=0 background.');
assert(opts.order>=1 && opts.order<=20 && opts.residual_order>=opts.order, ...
    'deepsea3dFE:Order','Unsupported expansion order.');
timer=tic; cpu0=cputime;
L=opts.domain_halfwidth_m; n=2*ceil(L/opts.spacing_m); h=2*L/n;
x=-L+(0:n-1)*h; [xx,yy]=meshgrid(x,x);
eta=deepsea3d_surface(xx,yy,cfg.profile);
k=2*pi*cfg.frequency_hz/cfg.c0_mps;
kv=2*pi/(2*L)*[0:n/2-1 -n/2:-1];
kz=sqrt(complex(k^2-kv.^2-kv.'.^2)); dz=1i*kz;
points=[xx(:)+cfg.center_m(1),yy(:)+cfg.center_m(2),zeros(n^2,1)];
R0=reshape(vecnorm(points-cfg.tx_m,2,2),n,n);
zt=cfg.tx_m(3);
N=opts.order; R=cell(N+1,1); invR=R; E=R; forcing=R;
R{1}=R0; invR{1}=1./R0; E{1}=exp(1i*k*R0); forcing{1}=zeros(n,n);
for m=1:N
    a=zeros(n,n);
    if m==1, a=-2*zt*eta; elseif m==2, a=eta.^2; end
    for j=1:m-1, a=a-R{j+1}.*R{m-j+1}; end
    R{m+1}=a./(2*R0);
    a=zeros(n,n); b=complex(a);
    for j=1:m
        a=a+R{j+1}.*invR{m-j+1};
        b=b+j*R{j+1}.*E{m-j+1};
    end
    invR{m+1}=-a./R0; E{m+1}=1i*k*b/m;
    g=complex(zeros(n,n));
    for j=0:m, g=g+E{j+1}.*invR{m-j+1}; end
    if mod(m,2)==1, forcing{m+1}=2*g/(4*pi); else, forcing{m+1}=zeros(n,n); end
end
clear R invR E points xx yy R0;
offset=cfg.rx_m(1:2)-cfg.center_m;
rx_kernel=exp(1i*kz*cfg.rx_m(3)).*exp(1i*kv*(offset(1)+L)).*exp(1i*kv.'*(offset(2)+L));
[B,terms,bnorm]=deepsea3d_fe_recursion(eta,dz,forcing,rx_kernel);
partials=cumsum(terms);
totalB=complex(zeros(n,n));
for m=1:N, totalB=totalB+B{m}; end
% True rough-boundary evaluation: Taylor expansion of outgoing exp(i*kz*eta)
% uses a separately chosen larger M, not the N-degree enforced identity.
surface=complex(zeros(n,n)); power=ones(n,n); factorial=1;
for j=0:opts.residual_order
    if j>0, power=power.*eta; factorial=factorial*j; end
    surface=surface+power/factorial.*ifft2(dz.^j.*totalB);
end
[xx,yy]=meshgrid(x,x);
pts=[xx(:)+cfg.center_m(1),yy(:)+cfg.center_m(2),eta(:)];
ui=deepsea3d_point_field(pts,cfg.tx_m,k);
image=cfg.tx_m; image(3)=-image(3);
flat=reshape(ui-deepsea3d_point_field(pts,image,k),n,n);
mask=hypot(xx,yy)<cfg.profile.support_m;
res=flat+surface;
den=norm(flat(mask));
forcing_sum=zeros(n,n);
for m=1:N, forcing_sum=forcing_sum+forcing{m+1}; end
bg=deepsea3d_flat(cfg,0);
out=struct('rough_correction',partials(end),'reflect_native',bg.reflect_native+partials(end), ...
    'total_native',bg.direct_native+bg.reflect_native+partials(end), ...
    'terms',terms,'partials',partials,'boundary_coefficient_norms',bnorm, ...
    'boundary_residual_relative',norm(res(mask))/max(den,realmin), ...
    'boundary_residual_abs_max',max(abs(res(mask))), ...
    'forcing_series_relative',norm(forcing_sum(mask)-flat(mask))/max(den,realmin), ...
    'boundary_residual_all_to_incident',norm(res(:))/norm(ui), ...
    'spacing_m',h,'nodes',n^2,'order',N,'residual_order',opts.residual_order, ...
    'wall_s',toc(timer),'cpu_s',cputime-cpu0, ...
    'stored_spectral_MiB',16*n^2*N/2^20,'single_array_MiB',16*n^2/2^20, ...
    'max_eta_kz',max(abs(eta(:)))*max(abs(kz(:))), ...
    'max_height_m',max(abs(eta(:))),'time_convention','exp(-i omega t)');
if isfield(opts,'offgrid') && opts.offgrid
    % Half-cell shift evaluates the same Fourier field on independently
    % sampled true geometry; not the collocation grid used by the recurrence.
    [xx,yy]=meshgrid(x+h/2,x+h/2);
    eta_shift=deepsea3d_surface(xx,yy,cfg.profile);
    phase=exp(1i*kv*h/2).*exp(1i*kv.'*h/2);
    pressure=complex(zeros(n,n)); power=ones(n,n); factorial=1;
    for j=0:opts.residual_order
        if j>0, power=power.*eta_shift; factorial=factorial*j; end
        pressure=pressure+power/factorial.*ifft2(dz.^j.*totalB.*phase);
    end
    pts=[xx(:)+cfg.center_m(1),yy(:)+cfg.center_m(2),eta_shift(:)];
    exact_flat=reshape(deepsea3d_point_field(pts,cfg.tx_m,k)-deepsea3d_point_field(pts,image,k),n,n);
    mask=hypot(xx,yy)<cfg.profile.support_m;
    residual=pressure+exact_flat;
    out.offgrid_boundary_residual_relative=norm(residual(mask))/norm(exact_flat(mask));
    out.offgrid_boundary_abs_max=max(abs(residual(mask)));
    out.wall_s=toc(timer); out.cpu_s=cputime-cpu0;
end
end
