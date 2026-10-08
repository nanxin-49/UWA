function plan = prepare_spectral_fast_kirchhoff_2d(cfg)
%PREPARE_SPECTRAL_FAST_KIRCHHOFF_2D Global low-rank FFT kernel plan.
% No receiver-by-surface kernel is assembled. Tables are lag-by-fixed-height.
% Required cfg: frequency_hz,c0_mps,surface_x_m,height_bound_m,
% receiver_depth_m,spectral_source. Ranks use one global basis on [-A,A].
k=2*pi*cfg.frequency_hz/cfg.c0_mps;
x=cfg.surface_x_m(:); n=numel(x); dx=mean(diff(x));
assert(n>=3 && dx>0 && max(abs(diff(x)-dx))<1e-10*dx);
A=cfg.height_bound_m; z=cfg.receiver_depth_m;
assert(A>0 && z>A && isfinite(A));
source=cfg.spectral_source;
kx=source.kx(:); kz=source.kz(:); co=source.coeff(:);
assert(numel(kx)==numel(co) && numel(kz)==numel(co));
assert(all(abs(kx)<=k),'This first incident basis requires propagating saved bins.');
assert(max(abs(imag(kz)))<1e-12 && max(abs(kz.^2+kx.^2-k^2))<1e-8);
max_rank=4; if isfield(cfg,'max_rank'), max_rank=cfg.max_rank; end
q=33; if isfield(cfg,'height_nodes'), q=cfg.height_nodes; end
assert(q>=max_rank+1 && max_rank>=1);
clock=tic; cpu0=cputime;
period=source.period_m; n_period=round(period/dx);
assert(abs(n_period*dx-period)<1e-9*period,'Source period must align with dx.');
harmonic=round(kx*period/(2*pi));
assert(max(abs(kx-harmonic*2*pi/period))<1e-10);
source_index=round((x-source.x0_m)/dx);
assert(max(abs(x-source.x0_m-source_index*dx))<1e-8*dx);
% Same weighted right subspace as the preceding audit, not a new sweep.
weights=abs(co).^2; positive=weights>0;
kx=kx(positive); kz=real(kz(positive)); co=co(positive);
harmonic=harmonic(positive); weights=weights(positive)/sum(weights);
eta_grid=linspace(-A,A,257).';
C=exp(1i*eta_grid*(k-kz.')).*sqrt(weights.');
[~,S,V]=svd(C,'econ'); V=V(:,1:max_rank);
nodes=A*cos(pi*(0:q-1)'/(q-1));
height_inc=exp(1i*nodes*(k-kz.')).*sqrt(weights.')*V;
% u_l(x), ux_l(x), uz_l(x) are evaluated by IFFT on the original period.
% Accumulating aliased harmonic indices is an EXACT evaluation on this grid;
% no saved spectral bin is cut, even on the coarse scaling-control grid.
coeff=(co.*exp(1i*kz*source.z_tx_m)).*conj(V)./sqrt(weights);
period_fields=complex(zeros(n_period,max_rank,3));
for r=1:max_rank
 for derivative=1:3
  multiplier=ones(size(kx));
  if derivative==2, multiplier=1i*kx; end
  if derivative==3, multiplier=-1i*kz; end
  bins=mod(harmonic,n_period)+1;
  a=coeff(:,r).*multiplier;
  embedded=complex(accumarray(bins,real(a),[n_period 1]),accumarray(bins,imag(a),[n_period 1]));
  period_fields(:,r,derivative)=ifft(embedded)*n_period/source.n_original;
 end
end
incoming_modes=period_fields(mod(source_index,n_period)+1,:,:);
incoming_height=fast_fk_height_interpolate(nodes,height_inc,cfg.surface_z_m(:));
% Exact sampled free-space Green kernel and SOURCE-coordinate derivatives.
% Residualize the height carrier; never Taylor-expand exp(-i*k*eta).
lag=(-(n-1):(n-1))'*dx; separation=z-nodes.';
R=hypot(lag,separation);
H0=1i/4*besselh(0,1,k*R);
H1=1i*k/4*besselh(1,1,k*R)./R;
carrier=exp(1i*k*nodes.');
B=[H0;H1.*separation/k;H1.*lag/k].*carrier;
[U,So,Vo]=svd(B,'econ');
U=U(:,1:max_rank)*So(1:max_rank,1:max_rank);
m=numel(lag);
kernels=cat(3,U(1:m,:),k*U(m+(1:m),:),k*U(2*m+(1:m),:));
outgoing_height=fast_fk_height_interpolate(nodes,conj(Vo(:,1:max_rank)),cfg.surface_z_m(:));
plan=struct('schema_version','1.0','n',n,'dx_m',dx,'x_m',x,'eta_m',cfg.surface_z_m(:), ...
 'receiver_depth_m',z,'k_radpm',k,'height_bound_m',A,'height_nodes',q, ...
 'max_rank',max_rank,'n_source_period',n_period,'incoming_modes',incoming_modes, ...
 'incoming_height',incoming_height,'outgoing_height',outgoing_height,'lag_kernels',kernels, ...
 'incoming_singular_values',diag(S),'outgoing_singular_values',diag(So), ...
 'setup_wall_seconds',toc(clock),'setup_cpu_seconds',cputime-cpu0, ...
 'incoming_evanescent','none in the saved incident bins', ...
 'outgoing_evanescent','retained implicitly by sampled outgoing Hankel kernels; no spectral cutoff');
payload=whos('C','B','U','V','Vo','kernels','period_fields');
plan.setup_named_arrays_bytes=sum([payload.bytes]);
payload=whos('plan'); plan.plan_bytes=payload.bytes;
end
