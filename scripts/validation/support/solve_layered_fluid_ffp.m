function out=solve_layered_fluid_ffp(cfg)
%SOLVE_LAYERED_FLUID_FFP Exact-J0 inverse Hankel with a closed contour detour.
% Reuses project Green and Gauss references; no Bellhop or PE calls.
arguments
    cfg (1,1) struct
end
t=tic; x=cfg.range_m(:).'; z=cfg.depth_m(:); zs=cfg.source_depth_m;
assert(cfg.contour_offset_radpm>0 && cfg.kappa_max_radpm>0 && cfg.panel_width_radpm>0);
assert(numel(x)>0 && all(isfinite(x)) && all(isfinite(z)));
[nodes,weights]=deepsea_bie_gauss(cfg.gauss_order);
[kappa,wk]=local_contour(cfg,nodes,weights);
s=find(zs>=cfg.layer_edges_m(1:end-1) & zs<cfg.layer_edges_m(2:end),1,'last');
c_ref=cfg.layer_c_mps(s); k_ref=2*pi*cfg.frequency_hz/c_ref;
source_mask=(z==zs) & (x==0); green=complex(zeros(numel(z),numel(x)));
use_ref=true(size(z));
if ~cfg.subtract_reference, use_ref(:)=false; end
% Surface rows integrate the exact kernel directly; avoid subtraction cancellation.
if strcmp(cfg.top_type,'pressure_release'), use_ref(z==0)=false; end
if any(use_ref)
    ref_rows=find(use_ref);
    for first_row=1:128:numel(ref_rows)
        rows=ref_rows(first_row:min(first_row+127,numel(ref_rows)));
        [zz,xx]=ndgrid(z(rows),x); points=[xx(:),zeros(numel(xx),1),zz(:)];
        valid=any(points~=[0 0 zs],2); p=complex(NaN(numel(xx),1));
        p(valid)=deepsea3d_point_field(points(valid,:),[0 0 zs],k_ref);
        green(rows,:)=reshape(p,size(xx));
    end
end
min_den=Inf; max_factor=0;
for first=1:cfg.spectral_block:numel(kappa)
    if toc(t)>cfg.max_wall_s, error('munk_ffp:TimeGate','Reached per-job wall-time gate.'); end
    if isfield(cfg,'minimum_available_bytes')
        [~,mem]=memory;
        assert(mem.PhysicalMemory.Available>=cfg.minimum_available_bytes,'munk_ffp:RAMGate','Available RAM below saved gate.');
    end
    ind=first:min(first+cfg.spectral_block-1,numel(kappa)); kv=kappa(ind).';
    [g,m]=layered_fluid_green_kernel(cfg,kv,z);
    if any(use_ref)
        q=sqrt(complex(k_ref^2-kv.^2)); q(imag(q)<0)=-q(imag(q)<0);
        ref=(1i./(2*q)).*exp(1i*abs(z(use_ref)-zs).*q);
        g(use_ref,:)=g(use_ref,:)-ref;
    end
    transform=besselj(0,kappa(ind)*abs(x));
    scaled=g.*(wk(ind).'.*kv/(2*pi));
    for first_row=1:256:numel(z)
        rows=first_row:min(first_row+255,numel(z));
        green(rows,:)=green(rows,:)+scaled(rows,:)*transform;
    end
    min_den=min(min_den,m.minimum_source_denominator); max_factor=max(max_factor,m.max_propagation_factor);
end
green(source_mask)=NaN;
assert(all(isfinite(green(~source_mask))));
out=struct('cfg',cfg,'green',green,'pressure_1m',4*pi*green,'source_singularity_mask',source_mask, ...
    'range_m',x,'depth_m',z,'kappa',kappa,'kappa_weights',wk,'elapsed_s',toc(t), ...
    'meta',struct('normalization','Green=e^(ikR)/(4piR); pressure_1m=4pi*Green, not calibrated Pa', ...
    'inverse_transform','(1/(2pi))*integral g(kappa,z)*J0(kappa*abs(x))*kappa dkappa', ...
    'contour','0 -> -i*a -> K-i*a -> K; all connectors retained', ...
    'time_convention','exp(-i omega t)','minimum_source_denominator',min_den, ...
    'maximum_propagation_factor',max_factor,'source_layer_c_mps',c_ref,'spectral_node_count',numel(kappa)));
end
function [q,w]=local_contour(c,xg,wg)
a=c.contour_offset_radpm; K=c.kappa_max_radpm;
ends=[0;-1i*a;K-1i*a;K]; q=[]; w=[];
for s=1:3
    count=max(1,ceil(abs(ends(s+1)-ends(s))/c.panel_width_radpm));
    edges=linspace(ends(s),ends(s+1),count+1);
    half=(edges(2:end)-edges(1:end-1))/2; mid=(edges(2:end)+edges(1:end-1))/2;
    points=xg*half+mid; weights=wg*half;
    q=[q;points(:)]; w=[w;weights(:)]; %#ok<AGROW>
end
end
