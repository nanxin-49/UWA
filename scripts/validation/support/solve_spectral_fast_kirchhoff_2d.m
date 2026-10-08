function out = solve_spectral_fast_kirchhoff_2d(cfg, plan)
%SOLVE_SPECTRAL_FAST_KIRCHHOFF_2D Validation-only low-rank linear FFT FK.
% cfg.rank, padding_factor (1/2/4), geometry ('A'/'B'), surface_normal.
% Receivers are the aligned uniform surface x grid at a single depth.
% Padding 1 is a deliberately aliased control; >=2 is linear convolution.
if nargin<2, plan=prepare_spectral_fast_kirchhoff_2d(cfg); end
r=cfg.rank; ro=r; if isfield(cfg,'outgoing_rank'), ro=cfg.outgoing_rank; end
assert(r>=1 && r<=plan.max_rank && ro>=1 && ro<=plan.max_rank);
n=plan.n; p=cfg.padding_factor; assert(ismember(p,[1 2 4]));
clock=tic; cpu0=cputime;
carrier=exp(-1i*plan.k_radpm*plan.eta_m);
u=carrier.*sum(plan.incoming_height(:,1:r).*plan.incoming_modes(:,1:r,1),2);
ux=carrier.*sum(plan.incoming_height(:,1:r).*plan.incoming_modes(:,1:r,2),2);
uz=carrier.*sum(plan.incoming_height(:,1:r).*plan.incoming_modes(:,1:r,3),2);
normal=cfg.surface_normal;
assert(isequal(size(normal),[n 2]));
if strcmp(cfg.geometry,'A')
 normal=[zeros(n,1),ones(n,1)];
 q=plan.dx_m*ones(n,1); q([1 end])=q([1 end])/2;
elseif strcmp(cfg.geometry,'B')
 ds=hypot(diff(plan.x_m),diff(plan.eta_m));
 q=[ds(1)/2;(ds(1:end-1)+ds(2:end))/2;ds(end)/2];
else
 error('solve_spectral_fast_kirchhoff_2d:Geometry','geometry must be A or B.');
end
dn=normal(:,1).*ux+normal(:,2).*uz;
nfft=p*n;
indices=mod((-(n-1):(n-1))',nfft)+1;
field=complex(zeros(n,1));
% Each term has FFT(source)*FFT(lag kernel), followed by one inverse FFT.
% Source density already includes trapezoidal dx/ds: no extra FFT dx factor.
for ell=1:ro
 height=carrier.*plan.outgoing_height(:,ell).*q;
 densities=[height.*dn,height.*u.*normal(:,2),height.*u.*normal(:,1)];
 sum_spectrum=complex(zeros(nfft,1));
 for derivative=1:3
  if strcmp(cfg.geometry,'A') && derivative==3, continue; end
  kernel=plan.lag_kernels(:,ell,derivative);
  embedded=complex(accumarray(indices,real(kernel),[nfft 1]),accumarray(indices,imag(kernel),[nfft 1]));
  sum_spectrum=sum_spectrum+fft(embedded).*fft(densities(:,derivative),nfft);
 end
 contribution=ifft(sum_spectrum);
 field=field-contribution(1:n);
end
payload=whos('u','ux','uz','dn','field','densities','sum_spectrum','embedded','contribution');
out=struct('receiver_x_m',plan.x_m,'receiver_z_m',plan.receiver_depth_m, ...
 'receiver_field',field,'surface_incident',u,'surface_incident_dn',dn, ...
 'rank',r,'outgoing_rank',ro,'geometry',cfg.geometry,'padding_factor',p,'fft_size',nfft, ...
 'wall_seconds',toc(clock),'cpu_seconds',cputime-cpu0, ...
 'named_apply_arrays_bytes',sum([payload.bytes]),'plan_bytes',plan.plan_bytes, ...
 'linear_convolution',nfft>=2*n-1,'dense_receiver_surface_entries',0, ...
 'fft_normalization','MATLAB forward unnormalized, inverse 1/L; quadrature in source density', ...
 'physical_support_m',[plan.x_m(1) plan.x_m(end)],'sponge',false);
assert(all(isfinite(field)));
end
