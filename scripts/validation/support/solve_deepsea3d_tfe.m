function out=solve_deepsea3d_tfe(cfg,opts)
%SOLVE_DEEPSEA3D_TFE Independent validation-only mapped Helmholtz expansion.
% Exact infinite-flat background + periodic compact correction; not FullFK.
arguments
    cfg (1,1) struct
    opts (1,1) struct
end
timer=tic;cpu0=cputime;L=opts.domain_halfwidth_m;N=opts.order;b=opts.top_depth_m;
n=2*ceil(L/opts.spacing_m);h=2*L/n;M=opts.vertical_degree;
assert(N>=1 && N<=40 && b>0 && b<cfg.rx_m(3) && b<cfg.tx_m(3));
estimated_workspace_MiB=12*16*n^2*(M+1)/2^20;
assert(estimated_workspace_MiB<=opts.max_workspace_MiB,'TFE:MemoryBudget','Named-array budget exceeded.');
x=-L+(0:n-1)*h;[xx,yy]=meshgrid(x,x);geo=deepsea3d_tfe_geometry(xx,yy,cfg.profile);
assert(cfg.profile.height_m==0 && max(geo.eta(:))<b,'TFE:Map','Invalid map or plane background.');
assert(L>cfg.profile.support_m,'TFE:PeriodicBuffer','Surface must vanish before lateral cell edge.');
[s,D,D2]=deepsea3d_tfe_chebyshev(M,b);k=2*pi*cfg.frequency_hz/cfg.c0_mps;
kv=2*pi/(2*L)*[0:n/2-1 -n/2:-1];gamma=sqrt(complex(k^2-kv.^2-kv.'.^2));
plan=deepsea3d_tfe_modal_plan(D,D2,gamma);
forcing=deepsea3d_tfe_forcing(cfg,xx,yy,geo.eta,N);
W1=complex(zeros(n^2,M+1));W2=W1;total=W1;terms=complex(zeros(N,1));
norms=zeros(N,1);conditions=zeros(N,1);linear_residual=zeros(N,1);
offset=cfg.rx_m(1:2)-cfg.center_m;
rx_kernel=exp(1i*gamma*(cfg.rx_m(3)-b)).*exp(1i*kv*(offset(1)+L)).*exp(1i*kv.'*(offset(2)+L));
for m=1:N
    [~,A,~]=deepsea3d_tfe_operators(W1,geo,s,D,D2,kv,k,b);
    [~,~,B]=deepsea3d_tfe_operators(W2,geo,s,D,D2,kv,k,b);
    F=-A-B;conditions(m)=(norm(A(:))+norm(B(:)))/max(norm(F(:)),realmin);
    clear A B;
    forcing_f=reshape(fft2(reshape(F,n,n,M+1)),n^2,M+1).';
    bottom=fft2(forcing{m+1});
    prevtop=reshape(W1(:,end),n,n);
    Tprev=ifft2(1i*gamma.*fft2(prevtop));
    top=fft2(-geo.eta/b.*Tprev);
    W_f=deepsea3d_tfe_modal_solve(plan,forcing_f,bottom(:),top(:));
    interior=plan.interior;
    residual=D2(interior,:)*W_f+gamma(:).'.^2.*W_f(interior,:)-forcing_f(interior,:);
    linear_residual(m)=norm(residual(:))/max(norm(forcing_f(interior,:),'fro')+norm(W_f,'fro')*k^2,realmin);
    clear forcing_f residual F;
    terms(m)=sum(reshape(W_f(end,:),n,n).*rx_kernel,'all')/n^2;
    W=reshape(ifft2(reshape(W_f.',n,n,M+1)),n^2,M+1);clear W_f;
    norms(m)=norm(W(:))*h/sqrt(n^2*(M+1));total=total+W;W2=W1;W1=W;clear W;
end
clear forcing;
% Independent true boundary data on the physical graph and half-cell graph.
bottom=reshape(total(:,1),n,n);pts=[xx(:)+cfg.center_m(1),yy(:)+cfg.center_m(2),geo.eta(:)];
image=cfg.tx_m;image(3)=-image(3);
flat=deepsea3d_point_field(pts,cfg.tx_m,k)-deepsea3d_point_field(pts,image,k);
mask=hypot(xx,yy)<cfg.profile.support_m;
boundary=norm(flat(mask(:))+bottom(mask))/max(norm(flat(mask(:))),realmin);
[xs,ys]=meshgrid(x+h/2,x+h/2);gs=deepsea3d_tfe_geometry(xs,ys,cfg.profile);
shift=exp(1i*kv*h/2).*exp(1i*kv.'*h/2);
bottomshift=ifft2(fft2(bottom).*shift);
pts=[xs(:)+cfg.center_m(1),ys(:)+cfg.center_m(2),gs.eta(:)];
flatshift=deepsea3d_point_field(pts,cfg.tx_m,k)-deepsea3d_point_field(pts,image,k);
maskshift=hypot(xs,ys)<cfg.profile.support_m;
offgrid=norm(flatshift(maskshift(:))+bottomshift(maskshift))/max(norm(flatshift(maskshift(:))),realmin);
% Physical Helmholtz residual recomputed at t=1, including half-grid geometry.
[P,A,B]=deepsea3d_tfe_operators(total,geo,s,D,D2,kv,k,b);
interior=2:M;physical=(P(:,interior)+A(:,interior)+B(:,interior))./(1-geo.eta(:)/b).^2;
pde=norm(physical(:))/max(k^2*norm(total(:,interior),'fro'),realmin);clear P A B physical;
shifted=reshape(ifft2(fft2(reshape(total,n,n,M+1)).*shift),n^2,M+1);
[P,A,B]=deepsea3d_tfe_operators(shifted,gs,s,D,D2,kv,k,b);
physical=(P(:,interior)+A(:,interior)+B(:,interior))./(1-gs.eta(:)/b).^2;
offgrid_pde=norm(physical(:))/max(k^2*norm(shifted(:,interior),'fro'),realmin);clear P A B physical shifted;
% True transformed DtN at t=1: physical dz=ds/(1-eta/b).
Ws=total*D(end,:).';Ttop=ifft2(1i*gamma.*fft2(reshape(total(:,end),n,n)));
transparent=norm(Ws./(1-geo.eta(:)/b)-Ttop(:))/max(norm(Ttop(:)),realmin);
bg=deepsea3d_flat(cfg,0);delta=sum(terms);
out=struct('rough_correction',delta,'reflect_native',bg.reflect_native+delta, ...
    'total_native',bg.total_native+delta,'terms',terms,'partials',cumsum(terms), ...
    'volume_coefficient_norms',norms,'recursion_cancellation',conditions, ...
    'linear_residual',linear_residual,'boundary_residual_relative',boundary, ...
    'offgrid_boundary_residual_relative',offgrid,'PDE_residual_relative',pde, ...
    'offgrid_PDE_residual_relative',offgrid_pde,'transparent_residual_relative',transparent, ...
    'order',N,'spacing_m',h,'domain_halfwidth_m',L,'vertical_degree',M,'top_depth_m',b, ...
    'min_map_jacobian',min(1-geo.eta(:)/b),'modal_basis_rcond',plan.basis_rcond, ...
    'near_pole_direct_columns',numel(plan.fallback),'volume_MiB',16*n^2*(M+1)/2^20, ...
    'estimated_workspace_MiB',estimated_workspace_MiB,'wall_s',toc(timer),'cpu_s',cputime-cpu0, ...
    'time_convention','exp(-i omega t)','lateral_boundary','periodic compact correction', ...
    'radiation','unfiltered outgoing/decaying Fourier DtN at z=b', ...
    'normalization','physical PDE residual/(k^2 volume norm); boundary residual/flat true-boundary norm');
end
