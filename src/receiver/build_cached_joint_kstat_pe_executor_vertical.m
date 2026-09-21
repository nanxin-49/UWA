function cache = build_cached_joint_kstat_pe_executor_vertical(cfg,pm_spec)
%BUILD_CACHED_JOINT_KSTAT_PE_EXECUTOR_VERTICAL Validation-only cached PE path.
%   The implementation intentionally supports the fixed uniform, no-bubble,
%   no-Doppler prototype only. It does not replace the public propagator.

arguments
    cfg (1,1) struct
    pm_spec (1,1) struct
end
required={'f_axis_hz','c0','xw','yw','nx','ny','x_tx','y_tx','z_tx', ...
    'x_rx','y_rx','z_rx','sigma_src_m','stepz_lamb','sponge_ratio', ...
    'alpha_max_np_per_m','reflect_coeff'};
for ii=1:numel(required)
    if ~isfield(cfg,required{ii}), error('Missing cfg.%s.',required{ii}); end
end
if isfield(cfg,'env_mode') && ~strcmpi(cfg.env_mode,'uniform')
    error('Cached prototype currently supports only env_mode=uniform.');
end
if isfield(cfg,'enable_bubbles') && cfg.enable_bubbles
    error('Cached prototype requires bubbles disabled.');
end

f_axis=cfg.f_axis_hz(:);
F=numel(f_axis);
dx=cfg.xw/cfg.nx; dy=cfg.yw/cfg.ny;
if abs(pm_spec.dx_m-dx)>1e-12 || abs(pm_spec.dy_m-dy)>1e-12
    error('PM and PE grids must have identical dx/dy for central-crop mapping.');
end
if pm_spec.nx<cfg.nx || pm_spec.ny<cfg.ny
    error('PM grid must not be smaller than the PE grid.');
end
x=(-0.5*cfg.xw):dx:(0.5*cfg.xw-dx);
y=(-0.5*cfg.yw):dy:(0.5*cfg.yw-dy);
[X,Y]=meshgrid(x,y);
[~,ix_rx]=min(abs(x-cfg.x_rx)); [~,iy_rx]=min(abs(y-cfg.y_rx));
psi_init=complex(exp(-((X-cfg.x_tx).^2+(Y-cfg.y_tx).^2)/(2*cfg.sigma_src_m^2)),0);
kx=(2*pi/cfg.xw)*[0:(cfg.nx/2-1),-cfg.nx/2:-1];
ky=(2*pi/cfg.yw)*[0:(cfg.ny/2-1),-cfg.ny/2:-1];
[KX,KY]=meshgrid(kx,ky); kappa2=KX.^2+KY.^2;
alpha_xy=local_absorption_profile(x,y,cfg.xw,cfg.yw,cfg.sponge_ratio,cfg.alpha_max_np_per_m);

incident_surface_xy_f=complex(zeros(cfg.ny,cfg.nx,F));
H_direct_f=complex(zeros(F,1));
H_ref_coh_f=complex(zeros(F,1));
R_coh_f=complex(zeros(F,1));
surface_rx_fr=complex(zeros(cfg.ny,cfg.nx,F));
surface_rx_screen=zeros(cfg.ny,cfg.nx,F);
surface_rx_nstep=zeros(F,1);
build_timer=tic;
for ii=1:F
    f_hz=f_axis(ii); k0=2*pi*f_hz/cfg.c0;
    psi_direct=local_march_uniform(psi_init,cfg.z_tx,cfg.z_rx,k0,f_hz,cfg,kappa2,alpha_xy);
    H_direct_f(ii)=psi_direct(iy_rx,ix_rx);
    incident_surface_xy_f(:,:,ii)=local_march_uniform( ...
        psi_init,cfg.z_tx,0,k0,f_hz,cfg,kappa2,alpha_xy);
    kernel=local_surface_rx_kernel(cfg.z_rx,k0,f_hz,cfg,kappa2,alpha_xy);
    surface_rx_fr(:,:,ii)=kernel.fr;
    surface_rx_screen(:,:,ii)=kernel.screen;
    surface_rx_nstep(ii)=kernel.n_step;
    R_coh_f(ii)=cfg.reflect_coeff*exp(-0.5*(2*k0)^2*pm_spec.sigma_eta_discrete2_m2);
    psi_coh=local_apply_kernel(R_coh_f(ii).*incident_surface_xy_f(:,:,ii),kernel);
    H_ref_coh_f(ii)=psi_coh(iy_rx,ix_rx);
end
build_time_s=toc(build_timer);

ix0=floor((pm_spec.nx-cfg.nx)/2)+1;
iy0=floor((pm_spec.ny-cfg.ny)/2)+1;
mapping=struct('ix',ix0:(ix0+cfg.nx-1),'iy',iy0:(iy0+cfg.ny-1), ...
    'pm_grid',[pm_spec.ny,pm_spec.nx],'pe_grid',[cfg.ny,cfg.nx], ...
    'dx_m',dx,'rule','same-dx central crop; no interpolation; no de-meaning');

cache=struct();
cache.kind='cached_joint_kstat_pe_executor_uniform_v1';
cache.cfg=cfg; cache.f_axis_hz=f_axis; cache.x=x; cache.y=y;
cache.ix_rx=ix_rx; cache.iy_rx=iy_rx;
cache.incident_surface_xy_f=incident_surface_xy_f;
cache.H_direct_reduced_f=H_direct_f;
cache.H_ref_coh_reduced_f=H_ref_coh_f;
% Legacy cache fields remain reduced so older validation scripts can load
% existing caches without changing the stored PE operator semantics.
cache.H_direct_f=H_direct_f;
cache.H_ref_coh_f=H_ref_coh_f;
cache.R_coh_f=R_coh_f;
cache.surface_rx_fr=surface_rx_fr;
cache.surface_rx_screen=surface_rx_screen;
cache.surface_rx_nstep=surface_rx_nstep;
cache.pm_mapping=mapping;
cache.phase_geometry=struct('z_tx',cfg.z_tx,'z_rx',cfg.z_rx, ...
    'z_surface',0,'c0',cfg.c0);
[deterministic_dsp,phase_meta]=apply_pe_channel_phase_reference_vertical( ...
    f_axis,struct('direct_f',H_direct_f,'reflect_coh_f',H_ref_coh_f), ...
    cache.phase_geometry,'direct_dsp');
cache.H_direct_dsp_f=deterministic_dsp.direct_f;
cache.H_ref_coh_dsp_f=deterministic_dsp.reflect_coh_f;
cache.phase_reference_meta=phase_meta;
cache.build_time_s=build_time_s;
cache.cache_array_bytes=numel(incident_surface_xy_f)*16+ ...
    numel(surface_rx_fr)*16+numel(surface_rx_screen)*8+ ...
    numel(H_direct_f)*16+numel(H_ref_coh_f)*16;
cache.memory_snapshot_bytes=local_memory_used_bytes();
cache.limitations=['Uniform c, no bubbles, no Doppler, CPU double PE only. ', ...
    'Validation interface; cached PE arrays remain reduced-envelope operators.'];
end

function psi_end=local_march_uniform(psi_start,z_start,z_end,k0,f_hz,cfg,kappa2,alpha_xy) %#ok<INUSD>
kernel=local_surface_rx_kernel(abs(z_end-z_start),k0,f_hz,cfg,kappa2,alpha_xy);
psi_end=local_apply_kernel(psi_start,kernel);
end

function kernel=local_surface_rx_kernel(distance,k0,f_hz,cfg,kappa2,alpha_xy) %#ok<INUSD>
n_step=max(1,ceil(abs(distance)/(cfg.stepz_lamb*(2*pi/k0))));
ds=abs(distance)/n_step;
denom=sqrt(complex(k0^2-kappa2,0))+k0;
kernel=struct('fr',exp(-1i*0.5*ds*kappa2./denom), ...
    'screen',exp(-ds*alpha_xy),'n_step',n_step,'ds_m',ds);
end

function psi_end=local_apply_kernel(psi_start,kernel)
psi_k=fft2(psi_start);
for jj=1:kernel.n_step
    psi_k=kernel.fr.*fft2(kernel.screen.*ifft2(kernel.fr.*psi_k));
end
psi_end=ifft2(psi_k);
end

function alpha_xy=local_absorption_profile(x,y,xw,yw,sponge_ratio,alpha_max)
Lx=sponge_ratio*xw; Ly=sponge_ratio*yw;
xs=0.5*xw-Lx; ys=0.5*yw-Ly;
ax=zeros(size(x)); ay=zeros(size(y));
mx=abs(x)>xs; my=abs(y)>ys;
ax(mx)=alpha_max*((abs(x(mx))-xs)/Lx).^3;
ay(my)=alpha_max*((abs(y(my))-ys)/Ly).^3;
ax=min(max(ax,0),alpha_max); ay=min(max(ay,0),alpha_max);
alpha_xy=ay(:)*ones(1,numel(x))+ones(numel(y),1)*ax(:).';
end

function value=local_memory_used_bytes()
value=NaN; try, info=memory; value=info.MemUsedMATLAB; catch, end
end
