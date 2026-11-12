function model=build_kirchhoff_kstat_joint_model_streaming_vertical( ...
    spec,f_axis_hz,c0_mps,reflect_coeff,options)
%BUILD_KIRCHHOFF_KSTAT_JOINT_MODEL_STREAMING_VERTICAL Compressed exact-series builder.
%   Re-expresses exp(alpha_i*alpha_j*C_eta)-1 as a converged power series.
%   It stores N_K-by-N_series scalar spectra, never F-by-F-by-N_K tensors,
%   uses one global augmented-frequency basis, and performs blockwise small EVDs.

arguments
    spec (1,1) struct
    f_axis_hz (:,1) double {mustBeFinite,mustBePositive}
    c0_mps (1,1) double {mustBeFinite,mustBePositive}
    reflect_coeff (1,1) double {mustBeFinite}
    options (1,1) struct = struct()
end
if ~isfield(options,'series_tail_tolerance'), options.series_tail_tolerance=1e-10; end
if ~isfield(options,'basis_variance_keep'), options.basis_variance_keep=1-1e-12; end
if ~isfield(options,'factor_variance_keep'), options.factor_variance_keep=1-1e-7; end
if ~isfield(options,'k_block_size'), options.k_block_size=512; end
f_axis_hz=f_axis_hz(:); F=numel(f_axis_hz); ny=spec.ny; nx=spec.nx; Nxy=nx*ny;
dxdy=spec.dx_m*spec.dy_m; fourier_scale=spec.dkx_rad_per_m*spec.dky_rad_per_m/(2*pi)^2;
C_eta=real(ifft2(spec.W_eta_kstat))*Nxy*fourier_scale; sigma2=max(C_eta(1,1),0);
alpha=4*pi*f_axis_hz/c0_mps; lambda=alpha.^2*sigma2;
R0=repmat(complex(reflect_coeff),F,1); R_coh=R0.*exp(-0.5*lambda);
n_series=local_poisson_order(max(lambda),options.series_tail_tolerance);

total_timer=tic; spectrum_timer=tic;
t=C_eta/max(sigma2,eps); t_power=ones(ny,nx);
s_modes=zeros(Nxy,n_series,'single'); q_scale=Nxy/dxdy;
for nn=1:n_series
    t_power=t_power.*t;
    sn=q_scale.*real(fft2(t_power))*dxdy;
    s_modes(:,nn)=single(sn(:));
end
spectrum_time=toc(spectrum_timer);

% Stable Poisson-amplitude modes: exp(-lambda/2)*lambda^(n/2)/sqrt(n!).
B=zeros(F,n_series); B(:,1)=real(R0).*exp(-0.5*lambda).*sqrt(lambda);
for nn=2:n_series, B(:,nn)=B(:,nn-1).*sqrt(lambda/nn); end
sign_n=(-1).^(1:n_series);
Baug=[B;B.*sign_n];
[Q,S,V]=svd(Baug,'econ'); singular=diag(S); energy=cumsum(singular.^2)/sum(singular.^2);
r_basis=find(energy>=options.basis_variance_keep,1,'first');
Q=Q(:,1:r_basis); T=S(1:r_basis,1:r_basis)*V(:,1:r_basis)';
basis_discard=1-sum(singular(1:r_basis).^2)/sum(singular.^2);

[pair_iy,pair_ix,pair_jy,pair_jx,pair_self]=local_pairs(nx,ny);
pair_count=numel(pair_iy); factor_cells=cell(pair_count,1); rank_list=zeros(pair_count,1,'uint8');
negative_sum=0; positive_sum=0; discarded_sum=0; evd_timer=tic;
for block_start=1:options.k_block_size:pair_count
    block_end=min(pair_count,block_start+options.k_block_size-1);
    for pp=block_start:block_end
        lin=sub2ind([ny,nx],double(pair_iy(pp)),double(pair_ix(pp)));
        s=double(s_modes(lin,:));
        if pair_self(pp)
            Ck=(B.*s)*B'; Pk=(B.*(s.*sign_n))*B.';
            RR=0.5*real(Ck+Pk); II=0.5*real(Ck-Pk);
            RI=0.5*(imag(Pk)-imag(Ck)); A=real([RR,RI;RI.',II]);
        else
            A=(T.*s)*T';
        end
        A=0.5*(A+A');
        [factor,em]=local_factor(A,options.factor_variance_keep);
        factor_cells{pp}=single(factor); rank_list(pp)=size(factor,2);
        negative_sum=negative_sum+em.negative; positive_sum=positive_sum+em.positive;
        discarded_sum=discarded_sum+em.discarded;
    end
end
evd_time=toc(evd_timer); peak_memory=local_memory_used_bytes();

% Retain the exact diagonal-spectrum independent comparison path.
sqrt_q_independent=zeros(ny,nx,F,'single'); independent_timer=tic;
for ii=1:F
    exponent=-lambda(ii); Cii=abs(R0(ii))^2.*(exp(exponent+alpha(ii)^2.*C_eta)-exp(exponent));
    Sii=real(fft2(Cii))*dxdy;
    sqrt_q_independent(:,:,ii)=sqrt(single(max(q_scale.*Sii,0)));
end
independent_time=toc(independent_timer);

factor_bytes=sum(cellfun(@(x)numel(x)*4,factor_cells));
model=struct('kind','kirchhoff_kstat_joint_streaming_series_v1', ...
    'f_axis_hz',f_axis_hz,'F',F,'nx',nx,'ny',ny, ...
    'alpha_f_rad_per_m',alpha,'R_coh_f',R_coh,'sigma_eta2_m2',sigma2, ...
    'reflect_coeff',reflect_coeff,'variance_keep',options.factor_variance_keep, ...
    'factor_cells',{factor_cells},'factor_basis',single(Q), ...
    'factor_basis_mode','nonself factors are compressed in global augmented-frequency basis; self factors are full real-augmented', ...
    'pair_iy',pair_iy,'pair_ix',pair_ix,'pair_jy',pair_jy,'pair_jx',pair_jx, ...
    'pair_is_self',pair_self,'rank_list',rank_list, ...
    'sqrt_q_independent',sqrt_q_independent,'n_series',n_series, ...
    'series_tail_tolerance',options.series_tail_tolerance, ...
    'basis_rank',r_basis,'basis_variance_keep',options.basis_variance_keep, ...
    'basis_discarded_energy_fraction',basis_discard, ...
    'k_block_size',options.k_block_size,'active_k_count',Nxy, ...
    'active_k_note','No K crop: spectral audits show broad K support; compression is frequency-basis based.', ...
    'spectrum_fft_count',n_series+F,'spectrum_build_time_s',spectrum_time, ...
    'factor_build_time_s',evd_time,'independent_build_time_s',independent_time, ...
    'total_build_time_s',toc(total_timer), ...
    'spectrum_peak_storage_bytes',numel(s_modes)*4, ...
    'factor_storage_bytes',factor_bytes+numel(Q)*4, ...
    'independent_storage_bytes',numel(sqrt_q_independent)*4, ...
    'max_rank',double(max(rank_list)),'mean_rank',mean(double(rank_list)), ...
    'negative_to_positive_eigenvalue_ratio',negative_sum/max(positive_sum,eps), ...
    'discarded_positive_variance_fraction',discarded_sum/max(positive_sum,eps), ...
    'peak_memory_snapshot_bytes',peak_memory, ...
    'formula_note',['Uses exact converged series in t=C_eta/sigma_eta^2: ', ...
        'C_G=sum_n b_n*b_n^H*t^n and P_G=sum_n(-1)^n*b_n*b_n^T*t^n.']);
end

function n=local_poisson_order(lambda,tol)
p=exp(-lambda); cumulative=p; n=0;
while 1-cumulative>tol
    n=n+1; p=p*lambda/n; cumulative=cumulative+p;
    if n>1000, error('Series order failed to converge.'); end
end
n=max(n,1);
end

function [iy,ix,jy,jx,self]=local_pairs(nx,ny)
count=(nx*ny+(1+double(mod(nx,2)==0))*(1+double(mod(ny,2)==0)))/2;
iy=zeros(count,1,'uint16'); ix=iy; jy=iy; jx=iy; self=false(count,1); visited=false(ny,nx); p=0;
for yy=1:ny
    yy2=mod(-(yy-1),ny)+1;
    for xx=1:nx
        if visited(yy,xx), continue; end
        xx2=mod(-(xx-1),nx)+1; visited(yy,xx)=true; visited(yy2,xx2)=true; p=p+1;
        iy(p)=yy; ix(p)=xx; jy(p)=yy2; jx(p)=xx2; self(p)=yy==yy2&&xx==xx2;
    end
end
end

function [factor,meta]=local_factor(A,keep_fraction)
[V,d0]=eig(A,'vector'); d0=real(d0); scale=max(max(abs(d0)),1); tol=1e-11*scale;
negative=d0<-tol; dpos=max(d0,0); [d,order]=sort(dpos,'descend'); total=sum(d);
if total<=tol, keep=false(size(d)); else, n=find(cumsum(d)>=keep_fraction*total,1); keep=false(size(d)); keep(1:n)=d(1:n)>tol; end
idx=order(keep); factor=V(:,idx).*sqrt(d(keep)).';
meta=struct('negative',sum(abs(d0(negative))), ...
    'positive',total,'discarded',sum(d(~keep)));
end

function value=local_memory_used_bytes()
value=NaN; try, info=memory; value=info.MemUsedMATLAB; catch, end
end
