%% Validate streaming-series F=64 builder against the original full builder
clear; close all; clc;
root=fileparts(fileparts(fileparts(mfilename('fullpath')))); addpath(root);
folder=fullfile(root,'results','validation','u5_conditional_channel_f64');
load(fullfile(folder,'f64_joint_pe_cache.mat'),'cache','joint_model'); original=joint_model; clear joint_model
load(fullfile(folder,'f64_joint_streaming_model_u5.mat'),'joint_model'); optimized=joint_model; clear joint_model

F=original.F; num2=0; den2=0;
timer=tic;
audit_index=unique([round(linspace(1,numel(original.factor_cells),256)),find(original.pair_is_self).']);
for pp=audit_index
    A=double(original.factor_cells{pp});
    if optimized.pair_is_self(pp), B=double(optimized.factor_cells{pp});
    else, B=double(optimized.factor_basis)*double(optimized.factor_cells{pp}); end
    aa=A'*A; bb=B'*B; cross=A'*B;
    num2=num2+norm(aa,'fro')^2+norm(bb,'fro')^2-2*norm(cross,'fro')^2;
    den2=den2+norm(aa,'fro')^2;
end
exact_time=toc(timer);
[power_error,adjacent_error,pseudo_error]=local_series_spectrum_errors(original,optimized,cache);
exact=struct('augmented_covariance_relative_error',sqrt(max(num2,0)/den2), ...
    'pseudo_diagonal_relative_error',pseudo_error, ...
    'per_frequency_power_relative_error',power_error, ...
    'adjacent_C_relative_error',adjacent_error, ...
    'audited_k_pair_count',numel(audit_index),'comparison_time_s',exact_time);
checkpoint=fullfile(folder,'validate_streaming_joint_builder_u5_exact_checkpoint.mat');
save(checkpoint,'exact');
original_meta=struct('build_s',original.total_build_time_s, ...
    'peak_memory_bytes',original.peak_memory_snapshot_bytes, ...
    'factor_bytes',original.factor_storage_bytes+original.independent_storage_bytes);
timer=tic; sample_kirchhoff_kstat_factor_model_vertical(original,8,2100001,'joint'); old_sample_s=toc(timer);
clear original

M=128; batch=4;
timer=tic; sample_kirchhoff_kstat_factor_model_vertical(optimized,8,2100001,'joint'); opt_sample_s=toc(timer);
[Hopt,topt,component_error]=generate_cached_kstat_ensemble_vertical(cache,optimized,M,2200000,'joint',batch);
E=load(fullfile(folder,'f64_full_ensembles.mat'),'ensemble'); Hold=E.ensemble.Hj_train;
if ~isfield(E.ensemble,'phase_reference_meta')
    geometry=struct('z_tx',cache.cfg.z_tx,'z_rx',cache.cfg.z_rx, ...
        'z_surface',0,'c0',cache.cfg.c0);
    [~,phase_meta]=apply_pe_channel_phase_reference_vertical(cache.f_axis_hz(:), ...
        struct('direct_f',cache.H_direct_f,'reflect_coh_f',cache.H_ref_coh_f), ...
        geometry,'direct_dsp');
    Hold=phase_meta.reflect_dsp_factor_f.*Hold;
end
S0=local_stats(Hold); S1=local_stats(Hopt);
t0=local_temporal(Hold,cache.f_axis_hz,0); t1=local_temporal(Hopt,cache.f_axis_hz,0);
receiver=struct('covariance_relative_error',norm(S1.C-S0.C,'fro')/norm(S0.C,'fro'), ...
    'correlation_relative_error',norm(S1.R-S0.R,'fro')/norm(S0.R,'fro'), ...
    'pdp_correlation',local_corr(t1.pdp,t0.pdp), ...
    'lfm_correlation',local_corr(local_lfm(Hopt,cache.f_axis_hz),local_lfm(Hold,cache.f_axis_hz)), ...
    'eigenvalue_spectrum_correlation',local_corr(S1.eig,S0.eig), ...
    'component_sum_error',component_error);
performance=struct('original_build_s',original_meta.build_s, ...
    'optimized_build_s',optimized.total_build_time_s,'build_speedup',original_meta.build_s/optimized.total_build_time_s, ...
    'original_peak_memory_bytes',original_meta.peak_memory_bytes, ...
    'optimized_peak_memory_bytes',optimized.peak_memory_snapshot_bytes, ...
    'original_factor_bytes',original_meta.factor_bytes, ...
    'optimized_factor_bytes',optimized.factor_storage_bytes+optimized.independent_storage_bytes, ...
    'old_sample_8_s',old_sample_s,'optimized_sample_8_s',opt_sample_s, ...
    'optimized_receiver_128_s',topt.total_s);
result=struct('exact',exact,'receiver',receiver,'performance',performance, ...
    'optimized_meta',rmfield(optimized,{'factor_cells','factor_basis','sqrt_q_independent'}));
save(fullfile(folder,'validate_streaming_joint_builder_u5_result.mat'),'result','Hopt','-v7.3');
fprintf('exact A %.3g C %.3g P %.3g power %.3g adj %.3g; PDP %.6f LFM %.6f; build %.3fx sample %.3fx\n', ...
    exact.augmented_covariance_relative_error,exact.augmented_covariance_relative_error, ...
    exact.pseudo_diagonal_relative_error,exact.per_frequency_power_relative_error,exact.adjacent_C_relative_error, ...
    receiver.pdp_correlation,receiver.lfm_correlation,performance.build_speedup,old_sample_s/opt_sample_s);

function s=local_stats(H)
X=H-mean(H,2); C=(X*X')/(size(H,2)-1); d=sqrt(max(real(diag(C)),0));
s=struct('C',C,'R',C./max(d*d.',eps),'eig',sort(max(real(eig(0.5*(C+C'))),0),'descend'));
end
function t=local_temporal(H,f,tref)
c=build_channel_cir_vertical(H,f,struct('input_reference','direct_dsp', ...
    'time_origin_shift_s',tref));
p=abs(c.h_tau).^2; p=mean(p,2); t=struct('pdp',p/sum(p));
end
function env=local_lfm(H,f)
fs=12000; N=512; tt=(0:N-1).'/fs; D=.02; active=tt<D; tx=zeros(N,1);
tx(active)=exp(1i*pi*((f(end)-f(1))/D)*(tt(active)-D/2).^2); fb=(-N/2:N/2-1).'*fs/N;
Hb=complex(zeros(N,size(H,2))); for m=1:size(H,2), Hb(:,m)=interp1(f-mean(f),H(:,m),fb,'linear',0); end
TX=fftshift(fft(tx)); rx=ifft(ifftshift(TX.*Hb),[],1); mf=ifft(fft(rx).*conj(fft(tx)),[],1);
env=mean(abs(mf),2); env=env/max(env);
end
function r=local_corr(a,b), q=corrcoef(real(a(:)),real(b(:))); r=q(1,2); end

function [power_error,adjacent_error,pseudo_error]=local_series_spectrum_errors(original,opt,cache)
F=original.F; spec=raw_pm_spectrum_grid_vertical(5,opt.nx,opt.ny,100,100); N=opt.nx*opt.ny;
Ceta=real(ifft2(spec.W_eta_kstat))*N*spec.dkx_rad_per_m*spec.dky_rad_per_m/(2*pi)^2;
sig=Ceta(1,1); t=Ceta/sig; alpha=original.alpha_f_rad_per_m; lambda=alpha.^2*sig;
B=zeros(F,opt.n_series); B(:,1)=-exp(-.5*lambda).*sqrt(lambda);
for n=2:opt.n_series, B(:,n)=B(:,n-1).*sqrt(lambda/n); end
tp=ones(size(t)); Sm=cell(opt.n_series,1);
for n=1:opt.n_series, tp=tp.*t; Sm{n}=real(fft2(tp)); end
pow_num=0; pow_den=0; adj_num=0; adj_den=0; p_num=0; p_den=0; dxdy=spec.dx_m*spec.dy_m;
for i=1:F
    e=-lambda(i); direct=real(fft2(exp(e+alpha(i)^2*Ceta)-exp(e)))*dxdy;
    series=zeros(size(t)); for n=1:opt.n_series, series=series+B(i,n)^2.*Sm{n}*dxdy; end
    pow_num=pow_num+norm(series-direct,'fro')^2; pow_den=pow_den+norm(direct,'fro')^2;
    pdirect=real(fft2(exp(e-alpha(i)^2*Ceta)-exp(e)))*dxdy;
    pseries=zeros(size(t)); for n=1:opt.n_series, pseries=pseries+(-1)^n*B(i,n)^2.*Sm{n}*dxdy; end
    p_num=p_num+norm(pseries-pdirect,'fro')^2; p_den=p_den+norm(pdirect,'fro')^2;
    if i<F
        e2=-.5*(lambda(i)+lambda(i+1)); prod=alpha(i)*alpha(i+1);
        d=real(fft2(exp(e2+prod*Ceta)-exp(e2)))*dxdy;
        s=zeros(size(t)); for n=1:opt.n_series, s=s+B(i,n)*B(i+1,n).*Sm{n}*dxdy; end
        adj_num=adj_num+norm(s-d,'fro')^2; adj_den=adj_den+norm(d,'fro')^2;
    end
end
power_error=sqrt(pow_num/pow_den); adjacent_error=sqrt(adj_num/adj_den); pseudo_error=sqrt(p_num/max(p_den,eps));
end
