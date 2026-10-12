function out = fk_single_receiver_prototype(action,cfg,plan)
%FK_SINGLE_RECEIVER_PROTOTYPE Isolated fixed-Tx/single-Rx research operators.
% No production path, reference solver, outgoing rank or convolution changes.
% FastIn plans bind incident dependencies; composite plans also bind receiver.
if nargin<3, plan=[]; end
x=cfg.surface_x_m(:); z=cfg.surface_z_m(:); n=numel(x);
k=2*pi*cfg.frequency_hz/cfg.c0_mps;
assert(n>=3 && all(diff(x)>0) && numel(z)==n && all(isfinite(z)));
assert(isscalar(cfg.receiver_x_m) && isscalar(cfg.receiver_z_m));
assert(cfg.receiver_z_m>cfg.height_bound_m && cfg.height_bound_m>0);
assert(max(abs(z))<=cfg.height_bound_m,'fk_single_receiver:HeightDomain','Height exceeds fixed plan domain.');
assert(isequal(size(cfg.surface_normal),[n 2]));
assert(max(abs(sum(cfg.surface_normal.^2,2)-1))<1e-10);
assert(all(cfg.surface_normal(:,2)>0));
src=cfg.spectral_source;
assert(all(abs(src.kx)<=k) && max(abs(imag(src.kz)))<1e-12, ...
    'fk_single_receiver:InvalidSpectrum','Propagating spectrum required.');
assert(max(abs(src.kz.^2+src.kx.^2-k^2))<1e-8, ...
    'fk_single_receiver:InvalidSpectrum','Spectrum dispersion mismatch.');
% Bind only fields read by precomputation; descriptive source metadata is inert.
source_fields={'kx','kz','coeff','period_m','n_original','x0_m','z_tx_m'};
source_signature=struct;
for j=1:numel(source_fields), source_signature.(source_fields{j})=src.(source_fields{j}); end
signature=struct('x',x,'source',source_signature,'f',cfg.frequency_hz,'c',cfg.c0_mps, ...
    'A',cfg.height_bound_m,'q',cfg.height_nodes);
is_fastin=strcmp(action,'prepare_fastin') || strcmp(action,'fastin');
if is_fastin
    signature.rank=cfg.rank; kind='fastin';
else
    signature.rx=cfg.receiver_x_m; signature.rz=cfg.receiver_z_m; kind='composite';
end
if startsWith(action,'prepare')
    nodes=cfg.height_bound_m*cos(pi*(0:cfg.height_nodes-1)'/(cfg.height_nodes-1));
    assert(cfg.height_nodes>=2 && (~is_fastin || (cfg.height_nodes>=cfg.rank+1 && cfg.rank>=1)));
    plan=struct('signature',signature,'nodes',nodes,'k',k,'kind',kind);
    if strcmp(action,'prepare_fastin')
        % Same energy-weighted residual-height SVD as existing FastFK.
        w=abs(src.coeff(:)).^2; keep=w>0; w=w(keep)/sum(w);
        kz=real(src.kz(keep)); co=src.coeff(keep); kx=src.kx(keep);
        eta=linspace(-cfg.height_bound_m,cfg.height_bound_m,257)';
        C=exp(1i*eta*(k-kz.')).*sqrt(w.');
        [~,S,V]=svd(C,'econ'); V=V(:,1:cfg.rank);
        plan.height=exp(1i*nodes*(k-kz.')).*sqrt(w.')*V;
        coefficients=co.*exp(1i*kz*src.z_tx_m).*conj(V)./sqrt(w);
        plan.modes=local_fourier(src,x,kx,kz,coefficients);
        plan.singular_values=diag(S);
    elseif strcmp(action,'prepare_composite')
        % Exact incident field and PHYSICAL derivatives at every height node.
        co=src.coeff(:).*exp(-1i*(nodes.'-src.z_tx_m).*src.kz(:));
        values=local_fourier(src,x,src.kx,src.kz,co);
        [G,Gx,Gz]=local_green(k,x,nodes.',cfg);
        carrier=exp(2i*k*nodes.');
        plan.Px=(G.*values(:,:,2)+Gx.*values(:,:,1)).*carrier;
        plan.Pz=(G.*values(:,:,3)+Gz.*values(:,:,1)).*carrier;
    else
        error('Unknown preparation action.');
    end
    out=plan; return
end
if strcmp(action,'fastin') || strcmp(action,'composite')
    assert(~isempty(plan),'fk_single_receiver:MissingPlan','Prepared plan required.');
    assert(isfield(plan,'kind') && strcmp(plan.kind,kind), ...
        'fk_single_receiver:PlanKind','Wrong plan kind.');
    assert(isequaln(plan.signature,signature), ...
        'fk_single_receiver:PlanMismatch','Stale plan/config mismatch.');
end
if strcmp(action,'full') || strcmp(action,'exact_composite')
    [u,ux,uz]=fk_single_receiver_incident(src,x,z);
elseif strcmp(action,'fastin')
    assert(isfield(plan,'modes'));
    h=fast_fk_height_interpolate(plan.nodes,plan.height,z).*exp(-1i*k*z);
    u=sum(h.*plan.modes(:,:,1),2);
    ux=sum(h.*plan.modes(:,:,2),2); uz=sum(h.*plan.modes(:,:,3),2);
elseif strcmp(action,'composite')
    assert(isfield(plan,'Px'));
    % Each x row has its own height table: no spatial interpolation.
    px=local_row_interpolate(plan.nodes,plan.Px,z).*exp(-2i*k*z);
    pz=local_row_interpolate(plan.nodes,plan.Pz,z).*exp(-2i*k*z);
    out=-sum((cfg.surface_normal(:,1).*px+cfg.surface_normal(:,2).*pz).*local_weights(x,z));
    assert(isfinite(out)); return
else
    error('Unknown application action.');
end
if strcmp(action,'exact_composite')
    [G,Gx,Gz]=local_green(k,x,z,cfg);
    out=-sum((cfg.surface_normal(:,1).*(G.*ux+Gx.*u)+ ...
        cfg.surface_normal(:,2).*(G.*uz+Gz.*u)).*local_weights(x,z));
else
    refcfg=struct('frequency_hz',cfg.frequency_hz,'c0_mps',cfg.c0_mps, ...
        'surface_x_m',x,'surface_z_m',z,'surface_normal',cfg.surface_normal, ...
        'surface_incident',u,'surface_incident_dn',cfg.surface_normal(:,1).*ux+cfg.surface_normal(:,2).*uz, ...
        'receiver_x_m',cfg.receiver_x_m,'receiver_z_m',cfg.receiver_z_m);
    ref=solve_full_kirchhoff_2d(refcfg); out=ref.receiver_field;
end
assert(isscalar(out) && isfinite(out));
end

function values=local_fourier(src,x,kx,kz,co)
% Exact IFFT synthesis on the original source period; no discarded bins.
dx=mean(diff(x)); np=round(src.period_m/dx);
assert(max(abs(diff(x)-dx))<1e-10*dx && abs(np*dx-src.period_m)<1e-9*src.period_m);
harm=round(kx(:)*src.period_m/(2*pi)); ix=round((x-src.x0_m)/dx);
assert(max(abs(kx(:)-harm*2*pi/src.period_m))<1e-10);
assert(max(abs(x-src.x0_m-ix*dx))<1e-8*dx);
values=complex(zeros(numel(x),size(co,2),3));
for d=1:3
    multiplier=ones(numel(kx),1);
    if d==2, multiplier=1i*kx(:); elseif d==3, multiplier=-1i*kz(:); end
    for j=1:size(co,2)
        a=co(:,j).*multiplier; bins=mod(harm,np)+1;
        embedded=complex(accumarray(bins,real(a),[np 1]),accumarray(bins,imag(a),[np 1]));
        field=ifft(embedded)*np/src.n_original;
        values(:,j,d)=field(mod(ix,np)+1);
    end
end
end

function [G,Gx,Gz]=local_green(k,x,z,cfg)
dx=cfg.receiver_x_m-x; dz=cfg.receiver_z_m-z; r=hypot(dx,dz);
assert(all(r(:)>0));
G=1i/4*besselh(0,1,k*r); a=1i*k/4*besselh(1,1,k*r)./r;
Gx=a.*dx; Gz=a.*dz;
end

function q=local_weights(x,z)
ds=hypot(diff(x),diff(z)); q=[ds(1)/2;(ds(1:end-1)+ds(2:end))/2;ds(end)/2];
end

function value=local_row_interpolate(nodes,table,z)
q=numel(nodes); w=(-1).^(0:q-1); w([1 end])=w([1 end])/2;
delta=z-nodes.'; [distance,nearest]=min(abs(delta),[],2);
exact=distance<=16*eps(max(1,max(abs(nodes))));
ratio=w./delta; ratio(exact,:)=0; den=sum(ratio,2); den(exact)=1;
value=sum(ratio.*table,2)./den;
rows=find(exact); value(rows)=table(sub2ind(size(table),rows,nearest(rows)));
end
