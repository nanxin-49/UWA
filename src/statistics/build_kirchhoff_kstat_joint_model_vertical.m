function model = build_kirchhoff_kstat_joint_model_vertical( ...
    spec, f_axis_hz, c0_mps, reflect_coeff, variance_keep)
%BUILD_KIRCHHOFF_KSTAT_JOINT_MODEL_VERTICAL Cache joint K/-K factors.
%   Validation-only model. It stores truncated single-precision augmented
%   spectral factors so many realization batches can be drawn without
%   rebuilding F-by-F covariance and pseudo-covariance spectra.

arguments
    spec (1,1) struct
    f_axis_hz (:,1) double {mustBeFinite, mustBePositive}
    c0_mps (1,1) double {mustBeFinite, mustBePositive}
    reflect_coeff (1,1) double {mustBeFinite}
    variance_keep (1,1) double {mustBeGreaterThan(variance_keep,0),mustBeLessThanOrEqual(variance_keep,1)} = 1-1e-7
end

f_axis_hz = f_axis_hz(:);
F = numel(f_axis_hz);
ny = spec.ny;
nx = spec.nx;
Nxy = nx*ny;
dxdy = spec.dx_m*spec.dy_m;
fourier_area_scale = spec.dkx_rad_per_m*spec.dky_rad_per_m/(2*pi)^2;
C_eta = real(ifft2(spec.W_eta_kstat))*Nxy*fourier_area_scale;
sigma2 = max(C_eta(1,1),0);
alpha_f = 4*pi*f_axis_hz/c0_mps;
R0_f = repmat(complex(reflect_coeff),F,1);
R_coh_f = R0_f.*exp(-0.5*alpha_f.^2*sigma2);

build_timer = tic;
S_cov = complex(zeros(ny,nx,F,F,'single'));
S_pseudo = complex(zeros(ny,nx,F,F,'single'));
for ii=1:F
    for jj=ii:F
        exponent0 = -0.5*(alpha_f(ii)^2+alpha_f(jj)^2)*sigma2;
        product_alpha = alpha_f(ii)*alpha_f(jj);
        C_xy = R0_f(ii)*conj(R0_f(jj)).*( ...
            exp(exponent0+product_alpha.*C_eta)-exp(exponent0));
        P_xy = R0_f(ii)*R0_f(jj).*( ...
            exp(exponent0-product_alpha.*C_eta)-exp(exponent0));
        Sc = fft2(C_xy)*dxdy;
        Sp = fft2(P_xy)*dxdy;
        S_cov(:,:,ii,jj)=single(Sc);
        S_cov(:,:,jj,ii)=single(conj(Sc));
        S_pseudo(:,:,ii,jj)=single(Sp);
        S_pseudo(:,:,jj,ii)=single(Sp);
    end
end
spectrum_build_time_s = toc(build_timer);
spectrum_storage_bytes = numel(S_cov)*8+numel(S_pseudo)*8;

pair_count = (Nxy+local_self_bin_count(nx,ny))/2;
factor_cells = cell(pair_count,1);
pair_iy = zeros(pair_count,1,'uint16');
pair_ix = zeros(pair_count,1,'uint16');
pair_jy = zeros(pair_count,1,'uint16');
pair_jx = zeros(pair_count,1,'uint16');
pair_is_self = false(pair_count,1);
rank_list = zeros(pair_count,1,'uint8');
visited = false(ny,nx);
q_scale = Nxy/dxdy;
negative_abs_sum=0;
positive_sum=0;
discarded_positive_sum=0;
pair_index=0;
factor_timer=tic;
for iy=1:ny
    jy=mod(-(iy-1),ny)+1;
    for ix=1:nx
        if visited(iy,ix), continue; end
        jx=mod(-(ix-1),nx)+1;
        visited(iy,ix)=true;
        visited(jy,jx)=true;
        pair_index=pair_index+1;
        pair_iy(pair_index)=iy; pair_ix(pair_index)=ix;
        pair_jy(pair_index)=jy; pair_jx(pair_index)=jx;
        pair_is_self(pair_index)=(iy==jy && ix==jx);

        Ck=q_scale.*squeeze(double(S_cov(iy,ix,:,:)));
        Cminus=q_scale.*squeeze(double(S_cov(jy,jx,:,:)));
        Pk=q_scale.*squeeze(double(S_pseudo(iy,ix,:,:)));
        Ck=0.5*(Ck+Ck'); Cminus=0.5*(Cminus+Cminus');
        if pair_is_self(pair_index)
            RR=0.5*real(Ck+Pk);
            II=0.5*real(Ck-Pk);
            RI=0.5*(imag(Pk)-imag(Ck));
            A=real([RR,RI;RI.',II]);
            A=0.5*(A+A.');
        else
            A=[Ck,Pk;Pk',conj(Cminus)];
            A=0.5*(A+A');
        end
        [factor,eig_meta]=local_truncated_factor(A,variance_keep);
        factor_cells{pair_index}=single(factor);
        rank_list(pair_index)=size(factor,2);
        negative_abs_sum=negative_abs_sum+eig_meta.negative_abs_sum;
        positive_sum=positive_sum+eig_meta.positive_sum;
        discarded_positive_sum=discarded_positive_sum+eig_meta.discarded_positive_sum;
    end
end
factor_build_time_s=toc(factor_timer);
peak_memory_snapshot_bytes=local_memory_used_bytes();
clear S_pseudo

sqrt_q_independent=zeros(ny,nx,F,'single');
for ii=1:F
    sqrt_q_independent(:,:,ii)=sqrt(single(max(q_scale.*real(double(S_cov(:,:,ii,ii))),0)));
end
clear S_cov

factor_storage_bytes=0;
for pp=1:pair_count
    factor_storage_bytes=factor_storage_bytes+numel(factor_cells{pp})*8;
end

model=struct();
model.kind='kirchhoff_kstat_joint_factor_model_v1';
model.f_axis_hz=f_axis_hz;
model.F=F; model.nx=nx; model.ny=ny;
model.alpha_f_rad_per_m=alpha_f;
model.R_coh_f=R_coh_f;
model.sigma_eta2_m2=sigma2;
model.reflect_coeff=reflect_coeff;
model.variance_keep=variance_keep;
model.factor_cells=factor_cells;
model.pair_iy=pair_iy; model.pair_ix=pair_ix;
model.pair_jy=pair_jy; model.pair_jx=pair_jx;
model.pair_is_self=pair_is_self;
model.rank_list=rank_list;
model.sqrt_q_independent=sqrt_q_independent;
model.spectrum_fft_count=F*(F+1);
model.spectrum_build_time_s=spectrum_build_time_s;
model.factor_build_time_s=factor_build_time_s;
model.total_build_time_s=toc(build_timer);
model.spectrum_peak_storage_bytes=spectrum_storage_bytes;
model.factor_storage_bytes=factor_storage_bytes;
model.independent_storage_bytes=numel(sqrt_q_independent)*4;
model.max_rank=double(max(rank_list));
model.mean_rank=mean(double(rank_list));
model.negative_to_positive_eigenvalue_ratio=negative_abs_sum/max(positive_sum,eps);
model.discarded_positive_variance_fraction=discarded_positive_sum/max(positive_sum,eps);
model.peak_memory_snapshot_bytes=peak_memory_snapshot_bytes;
end

function count=local_self_bin_count(nx,ny)
count=(1+double(mod(nx,2)==0))*(1+double(mod(ny,2)==0));
end

function [factor,meta]=local_truncated_factor(A,variance_keep)
[V,D]=eig(A,'vector');
D=real(D);
scale=max(max(abs(D)),1);
tol=1e-11*scale;
negative=D<-tol;
Dpos=max(D,0);
[Dsort,order]=sort(Dpos,'descend');
total=sum(Dsort);
if total<=tol
    keep=false(size(Dsort));
else
    nkeep=find(cumsum(Dsort)>=variance_keep*total,1,'first');
    nkeep=max(1,nkeep);
    keep=false(size(Dsort)); keep(1:nkeep)=Dsort(1:nkeep)>tol;
end
order_keep=order(keep);
factor=V(:,order_keep).*reshape(sqrt(D(order_keep)).',1,[]);
discarded_mask=true(size(Dsort)); discarded_mask(keep)=false;
meta=struct('negative_abs_sum',sum(abs(D(negative))), ...
    'positive_sum',total, ...
    'discarded_positive_sum',sum(Dsort(discarded_mask)));
end

function value=local_memory_used_bytes()
value=NaN;
try
    info=memory;
    value=info.MemUsedMATLAB;
catch
end
end
