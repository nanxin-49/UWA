%% Receiver validation for the cached joint-frequency kstat prototype
% Validation-only entry point. Public defaults and vertical_wape_propagator
% are intentionally not modified.
clear; close all; clc;
root_dir=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(root_dir);
setup_vertical_project();
addpath(fileparts(mfilename('fullpath')));
out_dir=fullfile(root_dir,'results','validation','cached_joint_kstat_pe_receiver');
if ~exist(out_dir,'dir'), mkdir(out_dir); end

cfg=struct('f_axis_hz',linspace(4000,8000,32).','c0',1500, ...
    'xw',50,'yw',50,'nx',128,'ny',128,'x_tx',0,'y_tx',0, ...
    'z_tx',100,'x_rx',0,'y_rx',0,'z_rx',3,'sigma_src_m',0.4, ...
    'stepz_lamb',0.5,'sponge_ratio',0.12,'alpha_max_np_per_m',0.15, ...
    'reflect_coeff',-1,'env_mode','uniform','enable_bubbles',false);
U_mps=5; L_train=128; L_test=64; main_batch_size=4;
seed=struct('kdomain_train',100001+(0:L_train-1), ...
    'kdomain_test',200001+(0:L_test-1),'joint_train_base',300000, ...
    'joint_test_base',400000,'independent_train_base',500000, ...
    'independent_test_base',600000,'public_consistency',700001, ...
    'batch_benchmark_base',800000);
assert(numel(unique([seed.kdomain_train,seed.kdomain_test, ...
    seed.joint_train_base+(1:ceil(L_train/main_batch_size)), ...
    seed.joint_test_base+(1:ceil(L_test/main_batch_size)), ...
    seed.independent_train_base+(1:ceil(L_train/main_batch_size)), ...
    seed.independent_test_base+(1:ceil(L_test/main_batch_size))])) == ...
    L_train+L_test+2*ceil(L_train/main_batch_size)+2*ceil(L_test/main_batch_size));

fprintf('Building raw-PM grid, cached PE executor, and joint model...\n');
pm_spec=raw_pm_spectrum_grid_vertical(U_mps,256,256,100,100);
cache=build_cached_joint_kstat_pe_executor_vertical(cfg,pm_spec);
joint_model=build_kirchhoff_kstat_joint_model_vertical( ...
    pm_spec,cfg.f_axis_hz,cfg.c0,cfg.reflect_coeff,1-1e-7);

fprintf('Generating explicit same-surface kdomain train/test ensembles...\n');
[Hkd_train,t_kd_train,map_train,err_kd_train]=local_explicit_ensemble( ...
    cache,pm_spec,seed.kdomain_train,main_batch_size);
[Hkd_test,t_kd_test,map_test,err_kd_test]=local_explicit_ensemble( ...
    cache,pm_spec,seed.kdomain_test,main_batch_size);
fprintf('Generating joint and independent kstat train/test ensembles...\n');
[Hj_train,t_j_train,err_j_train]=local_stat_ensemble( ...
    cache,joint_model,L_train,seed.joint_train_base,'joint',main_batch_size);
[Hj_test,t_j_test,err_j_test]=local_stat_ensemble( ...
    cache,joint_model,L_test,seed.joint_test_base,'joint',main_batch_size);
[Hi_train,t_i_train,err_i_train]=local_stat_ensemble( ...
    cache,joint_model,L_train,seed.independent_train_base,'independent',main_batch_size);
[Hi_test,t_i_test,err_i_test]=local_stat_ensemble( ...
    cache,joint_model,L_test,seed.independent_test_base,'independent',main_batch_size);

stats=struct();
stats.kdomain_train=local_frequency_stats(Hkd_train);
stats.kdomain_test=local_frequency_stats(Hkd_test);
stats.joint_train=local_frequency_stats(Hj_train);
stats.joint_test=local_frequency_stats(Hj_test);
stats.independent_train=local_frequency_stats(Hi_train);
stats.independent_test=local_frequency_stats(Hi_test);
metrics=struct();
metrics.joint=local_compare_stats(stats.joint_test,stats.kdomain_test);
metrics.independent=local_compare_stats(stats.independent_test,stats.kdomain_test);
metrics.train_test=struct( ...
    'kdomain_C_rel',local_rel_fro(stats.kdomain_train.C,stats.kdomain_test.C), ...
    'joint_C_rel',local_rel_fro(stats.joint_train.C,stats.joint_test.C), ...
    'independent_C_rel',local_rel_fro(stats.independent_train.C,stats.independent_test.C));

tau_ref_signed_s=0;
cir=struct();
cir.kdomain=local_cir_stats(Hkd_test,cfg.f_axis_hz,tau_ref_signed_s);
cir.joint=local_cir_stats(Hj_test,cfg.f_axis_hz,tau_ref_signed_s);
cir.independent=local_cir_stats(Hi_test,cfg.f_axis_hz,tau_ref_signed_s);
metrics.joint.pdp_correlation=local_vector_corr(cir.joint.pdp,cir.kdomain.pdp);
metrics.independent.pdp_correlation=local_vector_corr(cir.independent.pdp,cir.kdomain.pdp);

lfm=struct();
lfm.kdomain=local_lfm_stats(Hkd_test,cfg.f_axis_hz);
lfm.joint=local_lfm_stats(Hj_test,cfg.f_axis_hz);
lfm.independent=local_lfm_stats(Hi_test,cfg.f_axis_hz);
metrics.joint.lfm_correlation=local_vector_corr(lfm.joint.mean_envelope,lfm.kdomain.mean_envelope);
metrics.independent.lfm_correlation=local_vector_corr(lfm.independent.mean_envelope,lfm.kdomain.mean_envelope);

distribution=struct('kdomain',local_distribution_stats(Hkd_test), ...
    'joint',local_distribution_stats(Hj_test), ...
    'independent',local_distribution_stats(Hi_test));
mapping=local_mapping_summary([map_train;map_test]);
component_sum_max_abs=max([err_kd_train,err_kd_test,err_j_train,err_j_test,err_i_train,err_i_test]);

fprintf('Running same-input public/cached consistency and wall-time benchmark...\n');
performance=local_public_cached_benchmark(cfg,seed.public_consistency);
performance.cache_build_time_s=cache.build_time_s;
performance.cache_bytes=cache.cache_array_bytes;
performance.joint_model_build_time_s=joint_model.total_build_time_s;
performance.joint_model_bytes=joint_model.factor_storage_bytes+ ...
    joint_model.independent_storage_bytes;
performance.train_timing=struct('kdomain',t_kd_train,'joint',t_j_train, ...
    'independent',t_i_train);
performance.cached_L=local_cumulative_timing(t_j_train,[16,64,128]);
performance.batch=local_batch_benchmark(cache,joint_model,seed.batch_benchmark_base);
performance.peak_memory_snapshot_bytes=max([cache.memory_snapshot_bytes, ...
    joint_model.peak_memory_snapshot_bytes,performance.batch.peak_memory_snapshot_bytes],[],'omitnan');

properness_receiver=max([stats.kdomain_test.pseudo_ratio,stats.joint_test.pseudo_ratio]);
pass=struct();
pass.joint_covariance_better=metrics.joint.epsilon_C < metrics.independent.epsilon_C;
pass.pdp_correlation=metrics.joint.pdp_correlation>=0.9;
pass.lfm_correlation=metrics.joint.lfm_correlation>=0.9;
pass.component_sum=component_sum_max_abs<=1e-10;
pass.mapping_ci_stable=mapping.energy_ratio_ci95_halfwidth<=0.02;
pass.cached_public_consistency=performance.max_abs_total_error<=1e-10;
pass.cached_actual_speedup=performance.speedup_excluding_build>1;
pass.pseudo_covariance_negligible=properness_receiver<=0.05;
pass.ready_for_receiver_statistics=all(structfun(@(x)logical(x),pass));

validation=struct();
validation.generated_at=char(datetime('now','Format','yyyy-MM-dd HH:mm:ss Z'));
validation.config=cfg; validation.U_mps=U_mps; validation.pm_spec_meta= ...
    rmfield(pm_spec,{'Phi2D','W_eta_kstat','E1D','KX_rad_per_m','KY_rad_per_m','K_rad_per_m'});
validation.seeds=seed; validation.L_train=L_train; validation.L_test=L_test;
validation.main_batch_size=main_batch_size; validation.cache_meta=rmfield(cache, ...
    {'incident_surface_xy_f','surface_rx_fr','surface_rx_screen'});
validation.joint_model_meta=rmfield(joint_model,{'factor_cells','sqrt_q_independent'});
validation.H_ref_sca=struct('kdomain_train',Hkd_train,'kdomain_test',Hkd_test, ...
    'joint_train',Hj_train,'joint_test',Hj_test, ...
    'independent_train',Hi_train,'independent_test',Hi_test);
deterministic_f=cache.H_direct_dsp_f+cache.H_ref_coh_dsp_f;
validation.H_total=struct('kdomain_train',deterministic_f+Hkd_train, ...
    'kdomain_test',deterministic_f+Hkd_test,'joint_train',deterministic_f+Hj_train, ...
    'joint_test',deterministic_f+Hj_test,'independent_train',deterministic_f+Hi_train, ...
    'independent_test',deterministic_f+Hi_test);
validation.components=struct('H_direct_f',cache.H_direct_dsp_f, ...
    'H_ref_coh_f',cache.H_ref_coh_dsp_f, ...
    'H_direct_reduced_f',cache.H_direct_f, ...
    'H_ref_coh_reduced_f',cache.H_ref_coh_f, ...
    'phase_reference_meta',cache.phase_reference_meta, ...
    'component_sum_max_abs',component_sum_max_abs);
validation.stats=stats; validation.metrics=metrics; validation.cir=cir;
validation.lfm=lfm; validation.distribution=distribution;
validation.mapping=mapping; validation.performance=performance; validation.pass=pass;
save(fullfile(out_dir,'cached_joint_kstat_pe_receiver_validation.mat'),'validation','-v7.3');
local_make_figures(validation,out_dir);
local_write_text_summary(validation,fullfile(out_dir,'summary.txt'));
analyze_receiver_properness_null_vertical( ...
    fullfile(out_dir,'cached_joint_kstat_pe_receiver_validation.mat'),2000,900001);
fprintf('Saved validation outputs to %s\n',out_dir);

function [H,timing,mapping_error,max_component_error]=local_explicit_ensemble(cache,spec,seeds,batch_size)
F=numel(cache.f_axis_hz); L=numel(seeds); H=complex(zeros(F,L));
mapping_error=zeros(L,2); max_component_error=0; timing=local_empty_timing(L,batch_size);
nb=ceil(L/batch_size);
for bb=1:nb
    idx=(bb-1)*batch_size+1:min(bb*batch_size,L); B=numel(idx);
    timer=tic; delta=complex(zeros(cache.cfg.ny,cache.cfg.nx,F,B,'single'));
    for mm=1:B
        [eta,meta]=sample_raw_pm_surface_vertical(spec,seeds(idx(mm)));
        eta_pe=eta(cache.pm_mapping.iy,cache.pm_mapping.ix);
        mapping_error(idx(mm),:)=[meta.variance_about_zero_m2,mean(eta_pe(:).^2)];
        for ff=1:F
            alpha=4*pi*cache.f_axis_hz(ff)/cache.cfg.c0;
            delta(:,:,ff,mm)=single(cache.cfg.reflect_coeff.*exp(1i*alpha.*eta_pe)-cache.R_coh_f(ff));
        end
    end
    timing.generator_s(bb)=toc(timer);
    [r,m]=run_cached_joint_kstat_pe_executor_vertical(cache,delta);
    H(:,idx)=r.H_ref_sca_fm; timing.pe_s(bb)=m.elapsed_s;
    timing.memory_bytes(bb)=m.memory_snapshot_bytes;
    max_component_error=max(max_component_error,r.component_sum_max_abs_error);
end
timing=local_finish_timing(timing);
end

function [H,timing,max_component_error]=local_stat_ensemble(cache,model,L,seed_base,mode,batch_size)
F=model.F; H=complex(zeros(F,L)); max_component_error=0;
timing=local_empty_timing(L,batch_size); nb=ceil(L/batch_size);
for bb=1:nb
    idx=(bb-1)*batch_size+1:min(bb*batch_size,L);
    [delta,gm]=sample_kirchhoff_kstat_factor_model_vertical(model,numel(idx),seed_base+bb,mode);
    [r,pm]=run_cached_joint_kstat_pe_executor_vertical(cache,delta);
    H(:,idx)=r.H_ref_sca_fm; timing.generator_s(bb)=gm.elapsed_s;
    timing.pe_s(bb)=pm.elapsed_s; timing.memory_bytes(bb)=max(gm.peak_memory_snapshot_bytes,pm.memory_snapshot_bytes);
    max_component_error=max(max_component_error,r.component_sum_max_abs_error);
end
timing=local_finish_timing(timing);
end

function t=local_empty_timing(L,batch)
nb=ceil(L/batch); t=struct('L',L,'batch_size',batch,'generator_s',zeros(nb,1), ...
    'pe_s',zeros(nb,1),'memory_bytes',nan(nb,1));
end
function t=local_finish_timing(t)
t.total_s=sum(t.generator_s+t.pe_s); t.per_realization_s=t.total_s/t.L;
t.peak_memory_snapshot_bytes=max(t.memory_bytes,[],'omitnan');
end

function s=local_frequency_stats(H)
mu=mean(H,2); X=H-mu; L=size(H,2);
C=(X*X')/max(L-1,1); P=(X*X.')/max(L-1,1);
d=sqrt(max(real(diag(C)),0)); denom=d*d.'; R=C./max(denom,eps);
e=sort(max(real(eig(0.5*(C+C'))),0),'descend');
phase=unwrap(angle(H),[],1);
s=struct('mu',mu,'C',C,'P',P,'R',R,'adjacent_corr',diag(R,1), ...
    'phase_step_rms_rad',sqrt(mean(diff(phase,1,1).^2,2)), ...
    'eigvals',e,'pseudo_ratio',norm(P,'fro')/max(norm(C,'fro'),eps));
end

function m=local_compare_stats(a,b)
m=struct('epsilon_mu',norm(a.mu-b.mu)/(norm(b.mu)+eps), ...
    'epsilon_C',local_rel_fro(a.C,b.C),'epsilon_P',local_rel_fro(a.P,b.P), ...
    'correlation_matrix_error',local_rel_fro(a.R,b.R), ...
    'adjacent_corr_rmse',sqrt(mean(abs(a.adjacent_corr-b.adjacent_corr).^2)), ...
    'pseudo_ratio',a.pseudo_ratio);
end
function v=local_rel_fro(a,b), v=norm(a-b,'fro')/max(norm(b,'fro'),eps); end

function c=local_cir_stats(H,f,tau_ref_signed)
F=numel(f); df=mean(diff(f)); T=1/df; tau=(0:F-1).'*T/F;
phase=exp(-1i*2*pi*f(:)*tau_ref_signed);
h=ifft(H.*phase,[],1); p=abs(h).^2; pdp=mean(p,2); pdp=pdp/max(sum(pdp),eps);
L=size(H,2); rms_delay=zeros(L,1); tail99=zeros(L,1);
for mm=1:L
    [~,ip]=max(p(:,mm)); q=circshift(p(:,mm),1-ip); q=q/max(sum(q),eps);
    mt=sum(tau.*q); rms_delay(mm)=sqrt(sum((tau-mt).^2.*q));
    jj=find(cumsum(q)>=0.99,1); tail99(mm)=tau(jj);
end
[~,peak]=max(pdp);
c=struct('h_t',h,'pdp',pdp,'delay_axis_s',tau,'delta_tau_s',1/(f(end)-f(1)), ...
    'Tmax_s',T,'tau_ref_signed_s',tau_ref_signed,'peak_delay_s',tau(peak), ...
    'mean_rms_delay_s',mean(rms_delay),'mean_tail99_s',mean(tail99));
end

function s=local_lfm_stats(H,f)
F=numel(f); fs=12000; N=512; t=(0:N-1).'/fs; duration=0.02;
active=t<duration; k=(f(end)-f(1))/duration;
tx=zeros(N,1); tx(active)=exp(1i*pi*k*(t(active)-duration/2).^2);
fbb=(-N/2:N/2-1).'*fs/N; Hbb=complex(zeros(N,size(H,2)));
for mm=1:size(H,2)
    Hbb(:,mm)=interp1(f-mean(f),H(:,mm),fbb,'linear',0);
end
TX=fftshift(fft(tx)); rx=ifft(ifftshift(TX.*Hbb),[],1);
mf=ifft(fft(rx).*conj(fft(tx)),[],1); env=abs(mf);
s=struct('time_s',t,'mean_envelope',mean(env,2)/max(mean(env,2)), ...
    'mean_power',mean(env.^2,2)/max(mean(env.^2,2)),'fs_hz',fs,'duration_s',duration);
end

function d=local_distribution_stats(H)
idx=unique([1,ceil(size(H,1)/2),size(H,1)]); d=struct('frequency_indices',idx);
for kk=1:numel(idx)
    z=H(idx(kk),:); xr=real(z); xi=imag(z); a=abs(z); ph=angle(z);
    d.points(kk)=struct('real_skewness',local_skew(xr),'imag_skewness',local_skew(xi), ...
        'magnitude_skewness',local_skew(a),'real_kurtosis',local_kurt(xr), ...
        'imag_kurtosis',local_kurt(xi),'magnitude_kurtosis',local_kurt(a), ...
        'phase_resultant',abs(mean(exp(1i*ph)))); %#ok<AGROW>
end
end
function v=local_skew(x), q=x-mean(x); v=mean(q.^3)/max(mean(q.^2)^(3/2),eps); end
function v=local_kurt(x), q=x-mean(x); v=mean(q.^4)/max(mean(q.^2)^2,eps); end
function r=local_vector_corr(a,b), R=corrcoef(real(a(:)),real(b(:))); r=R(1,2); end

function s=local_mapping_summary(v)
ratio=v(:,2)./v(:,1); n=numel(ratio); hw=1.96*std(ratio)/sqrt(n);
s=struct('n',n,'full_variance_m2',v(:,1),'cropped_variance_m2',v(:,2), ...
    'energy_ratio_mean',mean(ratio),'energy_ratio_std',std(ratio), ...
    'energy_ratio_ci95',[mean(ratio)-hw,mean(ratio)+hw], ...
    'energy_ratio_ci95_halfwidth',hw,'relative_energy_error_mean',mean(abs(ratio-1)));
end

function p=local_public_cached_benchmark(cfg,sea_seed)
spec=raw_pm_spectrum_grid_vertical(5,cfg.nx,cfg.ny,cfg.xw,cfg.yw);
cache=build_cached_joint_kstat_pe_executor_vertical(cfg,spec);
[eta,~]=sample_raw_pm_surface_vertical(spec,sea_seed); F=numel(cfg.f_axis_hz);
delta=complex(zeros(cfg.ny,cfg.nx,F));
for ff=1:F
    alpha=4*pi*cfg.f_axis_hz(ff)/cfg.c0;
    delta(:,:,ff)=cfg.reflect_coeff*exp(1i*alpha*eta)-cache.R_coh_f(ff);
end
timer=tic; [r,m]=run_cached_joint_kstat_pe_executor_vertical(cache,delta); cached_s=toc(timer);
params=struct('f0',cfg.f_axis_hz.','enable_wideband',false,'c0',cfg.c0, ...
    'z_max',cfg.z_tx,'z_tx',cfg.z_tx,'z_rx',cfg.z_rx,'xw',cfg.xw,'yw',cfg.yw, ...
    'nx',cfg.nx,'ny',cfg.ny,'x_tx',cfg.x_tx,'y_tx',cfg.y_tx, ...
    'x_rx',cfg.x_rx,'y_rx',cfg.y_rx,'sigma_src_m',cfg.sigma_src_m, ...
    'stepz_lamb',cfg.stepz_lamb,'sponge_ratio',cfg.sponge_ratio, ...
    'alpha_max_np_per_m',cfg.alpha_max_np_per_m,'env_mode','uniform', ...
    'enable_bubbles',false,'enable_surface_reflection',true,'surface_reflect_coeff',cfg.reflect_coeff, ...
    'surface_boundary_model','kirchhoff_kdomain','surface_roughness_scale_mode','raw_pm', ...
    'sea_wind_speed',5,'sea_seed',sea_seed,'show_figures',false,'enforce_1_over_R',false, ...
    'save_mode','rx_only','use_gpu',false);
timer=tic; public=vertical_channel_model(params); public_s=toc(timer);
err_direct=max(abs(r.H_direct_f-public.H_direct_f));
err_ref=max(abs(r.H_ref_coh_f+r.H_ref_sca_fm-public.H_reflect_f));
err_total=max(abs(r.H_total_fm-public.H_f));
p=struct('public_single_s',public_s,'cached_single_s',cached_s, ...
    'cached_executor_reported_s',m.elapsed_s,'speedup_excluding_build',public_s/cached_s, ...
    'same_grid_cache_build_s',cache.build_time_s,'max_abs_direct_error',err_direct, ...
    'max_abs_reflect_error',err_ref,'max_abs_total_error',err_total);
end

function s=local_cumulative_timing(t,Lvals)
c=cumsum(t.generator_s+t.pe_s); B=t.batch_size;
s=struct('L',num2cell(Lvals),'wall_s',num2cell(c(Lvals/B).'), ...
    'per_realization_s',num2cell(c(Lvals/B).'./Lvals));
end

function p=local_batch_benchmark(cache,model,seed_base)
bs=[1,2,4,8]; L=16; total=zeros(size(bs)); peak=nan(size(bs));
for kk=1:numel(bs)
    [~,t,~]=local_stat_ensemble(cache,model,L,seed_base+100*kk,'joint',bs(kk));
    total(kk)=t.total_s; peak(kk)=t.peak_memory_snapshot_bytes;
end
p=struct('batch_size',bs,'L',L,'wall_s',total,'per_realization_s',total/L, ...
    'peak_memory_by_batch_bytes',peak,'peak_memory_snapshot_bytes',max(peak,[],'omitnan'));
end

function local_make_figures(v,out_dir)
names={'kdomain','independent','joint'};
fig=figure('Visible','off','Color','w','Position',[80 80 1200 380]);
for k=1:3, subplot(1,3,k); imagesc(abs(v.stats.([names{k} '_test']).R)); axis image; colorbar;
    title([names{k} ' |R|']); xlabel('f_j'); ylabel('f_i'); end
exportgraphics(fig,fullfile(out_dir,'receiver_frequency_correlation.png'),'Resolution',180); close(fig);
fig=figure('Visible','off','Color','w'); hold on;
for k=1:3, semilogy(v.stats.([names{k} '_test']).eigvals/max(v.stats.([names{k} '_test']).eigvals),'-o'); end
grid on; legend(names,'Location','best'); xlabel('mode'); ylabel('normalized eigenvalue');
exportgraphics(fig,fullfile(out_dir,'receiver_covariance_eigenvalues.png'),'Resolution',180); close(fig);
fig=figure('Visible','off','Color','w'); hold on;
for k=1:3, plot(v.cir.(names{k}).delay_axis_s*1e3,v.cir.(names{k}).pdp,'LineWidth',1.2); end
grid on; legend(names); xlabel('relative delay (ms)'); ylabel('normalized PDP');
exportgraphics(fig,fullfile(out_dir,'receiver_reflected_only_pdp.png'),'Resolution',180); close(fig);
fig=figure('Visible','off','Color','w'); hold on;
for k=1:3, plot(v.lfm.(names{k}).time_s*1e3,v.lfm.(names{k}).mean_envelope,'LineWidth',1.1); end
grid on; legend(names); xlabel('matched-filter delay (ms)'); ylabel('normalized envelope');
exportgraphics(fig,fullfile(out_dir,'receiver_lfm_matched_filter.png'),'Resolution',180); close(fig);
fig=figure('Visible','off','Color','w'); hold on;
for k=1:3
    z=real(v.H_ref_sca.([names{k} '_test'])(ceil(end/2),:)); z=(z-mean(z))/max(std(z),eps);
    q=sort(z); n=numel(q); qn=sqrt(2)*erfinv(2*((1:n)-0.5)/n-1); plot(qn,q,'.-');
end
plot([-3 3],[-3 3],'k--'); grid on; axis equal; legend(names); xlabel('normal quantile'); ylabel('sample quantile');
exportgraphics(fig,fullfile(out_dir,'receiver_center_frequency_qq.png'),'Resolution',180); close(fig);
end

function local_write_text_summary(v,path)
fid=fopen(path,'w'); cleaner=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'cached joint-kstat + PE receiver validation\n');
fprintf(fid,'U=%.1f m/s, F=%d, train=%d, test=%d\n',v.U_mps,numel(v.config.f_axis_hz),v.L_train,v.L_test);
fprintf(fid,'joint epsilon_C=%.6g, independent epsilon_C=%.6g\n',v.metrics.joint.epsilon_C,v.metrics.independent.epsilon_C);
fprintf(fid,'joint/independent PDP corr=%.6g / %.6g\n',v.metrics.joint.pdp_correlation,v.metrics.independent.pdp_correlation);
fprintf(fid,'joint/independent LFM corr=%.6g / %.6g\n',v.metrics.joint.lfm_correlation,v.metrics.independent.lfm_correlation);
fprintf(fid,'component max error=%.6g\n',v.components.component_sum_max_abs);
fprintf(fid,'mapping energy ratio mean=%.6g, CI95=[%.6g %.6g]\n',v.mapping.energy_ratio_mean,v.mapping.energy_ratio_ci95);
fprintf(fid,'public/cached total max error=%.6g, speedup=%.6gx\n',v.performance.max_abs_total_error,v.performance.speedup_excluding_build);
fprintf(fid,'kdomain/joint P-to-C=%.6g / %.6g\n',v.stats.kdomain_test.pseudo_ratio,v.stats.joint_test.pseudo_ratio);
fprintf(fid,'ready=%d\n',v.pass.ready_for_receiver_statistics);
end
