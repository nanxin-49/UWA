%% U=5 m/s, F=64 conditional receiver-channel generator validation
% Set U5_CONDITIONAL_MODE=smoke or full. Public defaults remain unchanged.
clear; close all; clc;
root_dir=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(root_dir); addpath(fileparts(mfilename('fullpath')));
mode=lower(strtrim(getenv('U5_CONDITIONAL_MODE'))); if isempty(mode), mode='smoke'; end
if ~ismember(mode,{'smoke','full'}), error('U5_CONDITIONAL_MODE must be smoke or full.'); end
out_dir=fullfile(root_dir,'results','validation','u5_conditional_channel_f64');
if ~exist(out_dir,'dir'), mkdir(out_dir); end

F=64; f_axis=linspace(4000,8000,F).'; U=5; batch_size=4;
if strcmp(mode,'smoke'), L_train=32; L_test=16; M_gen=1000; else, L_train=128; L_test=64; M_gen=10000; end
cfg=struct('f_axis_hz',f_axis,'c0',1500,'xw',50,'yw',50,'nx',128,'ny',128, ...
    'x_tx',0,'y_tx',0,'z_tx',100,'x_rx',0,'y_rx',0,'z_rx',3,'sigma_src_m',0.4, ...
    'stepz_lamb',0.5,'sponge_ratio',0.12,'alpha_max_np_per_m',0.15, ...
    'reflect_coeff',-1,'env_mode','uniform','enable_bubbles',false);
seed=struct('kdomain_train',1100001+(0:L_train-1),'kdomain_test',1200001+(0:L_test-1), ...
    'joint_train_base',1300000,'joint_test_base',1400000, ...
    'independent_train_base',1500000,'independent_test_base',1600000, ...
    'generator_base',1700000,'properness_base',1800000,'public',1900001);
cir_common_shift_s=0;
validation_run_meta=pe_phase_release_run_meta_vertical(root_dir,struct( ...
    'frequency_axis_hz',f_axis,'seed_definition',seed));

cache_file=fullfile(out_dir,'f64_joint_pe_cache.mat');
if exist(cache_file,'file')
    fprintf('Loading reusable F=64 PM/joint/PE cache...\n'); load(cache_file,'pm_spec','cache','joint_model');
else
    fprintf('Building reusable F=64 PM/joint/PE cache...\n');
    pm_spec=raw_pm_spectrum_grid_vertical(U,256,256,100,100);
    cache=build_cached_joint_kstat_pe_executor_vertical(cfg,pm_spec);
    joint_model=build_kirchhoff_kstat_joint_model_vertical(pm_spec,f_axis,cfg.c0,cfg.reflect_coeff,1-1e-7);
    schema_version='2.0.0'; phase_reference_meta=struct('target_reference','phase_neutral');
    save(cache_file,'pm_spec','cache','joint_model','schema_version', ...
        'phase_reference_meta','validation_run_meta','-v7.3');
end
phase_neutral_cache_audit=validate_phase_neutral_pe_cache_vertical(cache,pm_spec,joint_model, ...
    struct('f_axis_hz',f_axis,'nx',128,'ny',128,'pm_nx',256,'pm_ny',256, ...
    'z_tx',100,'z_rx',3,'c0',1500));
phase_geometry=struct('z_tx',cfg.z_tx,'z_rx',cfg.z_rx,'z_surface',0,'c0',cfg.c0);
[deterministic,phase_meta]=apply_pe_channel_phase_reference_vertical(f_axis, ...
    struct('direct_f',cache.H_direct_f,'reflect_coh_f',cache.H_ref_coh_f), ...
    phase_geometry,'direct_dsp');

ensemble_file=fullfile(out_dir,['f64_' mode '_ensembles.mat']);
if exist(ensemble_file,'file')
    fprintf('Loading checkpointed %s ensembles...\n',mode); load(ensemble_file,'ensemble');
    assert_pe_phase_release_artifact_vertical(ensemble,validation_run_meta.run_id,'U=5 ensemble checkpoint');
else
    fprintf('Generating %s kdomain/joint/independent train/test ensembles...\n',mode);
    [Hkd_train,t_kd_train,map_train,e1]=generate_cached_kdomain_ensemble_vertical(cache,pm_spec,seed.kdomain_train,batch_size);
    [Hkd_test,t_kd_test,map_test,e2]=generate_cached_kdomain_ensemble_vertical(cache,pm_spec,seed.kdomain_test,batch_size);
    [Hj_train,t_j_train,e3]=generate_cached_kstat_ensemble_vertical(cache,joint_model,L_train,seed.joint_train_base,'joint',batch_size);
    [Hj_test,t_j_test,e4]=generate_cached_kstat_ensemble_vertical(cache,joint_model,L_test,seed.joint_test_base,'joint',batch_size);
    [Hi_train,t_i_train,e5]=generate_cached_kstat_ensemble_vertical(cache,joint_model,L_train,seed.independent_train_base,'independent',batch_size);
    [Hi_test,t_i_test,e6]=generate_cached_kstat_ensemble_vertical(cache,joint_model,L_test,seed.independent_test_base,'independent',batch_size);
    ensemble=struct('Hkd_train',Hkd_train,'Hkd_test',Hkd_test,'Hj_train',Hj_train,'Hj_test',Hj_test, ...
        'Hi_train',Hi_train,'Hi_test',Hi_test,'timing',struct('kdomain_train',t_kd_train, ...
        'kdomain_test',t_kd_test,'joint_train',t_j_train,'joint_test',t_j_test, ...
        'independent_train',t_i_train,'independent_test',t_i_test), ...
        'mapping_train',map_train,'mapping_test',map_test,'component_error',max([e1,e2,e3,e4,e5,e6]), ...
        'phase_reference_meta',phase_meta,'schema_version','2.0.0', ...
        'validation_run_meta',validation_run_meta);
    schema_version=ensemble.schema_version; phase_reference_meta=phase_meta;
    save(ensemble_file,'ensemble','schema_version','phase_reference_meta', ...
        'validation_run_meta','-v7.3');
end
if ~isfield(ensemble,'phase_reference_meta')
    names={'Hkd_train','Hkd_test','Hj_train','Hj_test','Hi_train','Hi_test'};
    for kk=1:numel(names), ensemble.(names{kk})=phase_meta.reflect_dsp_factor_f.*ensemble.(names{kk}); end
    ensemble.phase_reference_meta=phase_meta;
end
Hkd_train=ensemble.Hkd_train; Hkd_test=ensemble.Hkd_test;
Hj_train=ensemble.Hj_train; Hj_test=ensemble.Hj_test; Hi_train=ensemble.Hi_train; Hi_test=ensemble.Hi_test;

stats=struct('kdomain_train',local_stats(Hkd_train),'kdomain_test',local_stats(Hkd_test), ...
    'joint_train',local_stats(Hj_train),'joint_test',local_stats(Hj_test), ...
    'independent_train',local_stats(Hi_train),'independent_test',local_stats(Hi_test));
properness=struct(); properness.kdomain=properness_null_test_vertical(stats.kdomain_train.C,Hkd_test,2000,seed.properness_base+1);
properness.joint=properness_null_test_vertical(stats.joint_train.C,Hj_test,2000,seed.properness_base+2);
properness.independent=properness_null_test_vertical(stats.independent_train.C,Hi_test,2000,seed.properness_base+3);

receiver=struct();
receiver.joint=local_compare(stats.joint_test,stats.kdomain_test,Hj_test,Hkd_test,f_axis,cir_common_shift_s);
receiver.independent=local_compare(stats.independent_test,stats.kdomain_test,Hi_test,Hkd_test,f_axis,cir_common_shift_s);
receiver.kdomain_temporal=local_temporal(Hkd_test,f_axis,cir_common_shift_s);
receiver.joint_temporal=local_temporal(Hj_test,f_axis,cir_common_shift_s);
receiver.independent_temporal=local_temporal(Hi_test,f_axis,cir_common_shift_s);

condition=struct('wind_speed_mps',U,'wind_convention',pm_spec.wind_definition, ...
    'Hs_implied_m',pm_spec.Hs_implied_discrete_m, ...
    'raw_pm_capture_ratio',pm_spec.capture_ratio_discrete_to_infinite, ...
    'surface_mode','raw_pm','surface_model','kirchhoff_kstat_joint_frequency', ...
    'bubbles','off','doppler','off','sound_speed','uniform', ...
    'pm_grid','100 m / 256 x 256','pe_grid','50 m / 128 x 128');
options=struct('frequency_axis_hz',f_axis,'H_direct_f',deterministic.direct_f, ...
    'H_ref_coh_f',deterministic.reflect_coh_f,'shrinkage_parameter',0, ...
    'properness_result',properness.joint,'train_seed_list',seed.joint_train_base+(1:ceil(L_train/batch_size)), ...
    'condition',condition,'reference_delay_s',0,'phase_reference_meta',phase_meta,'rank_selection','full');
model=estimate_conditional_channel_stats_vertical(Hj_train,options);

rank_names={'full','99.9','99'}; generator=struct();
for rr=1:numel(rank_names)
    name=rank_names{rr}; field=local_rank_field(name);
    [draw,tm]=sample_conditional_channel_vertical(model,M_gen,seed.generator_base+rr,struct('path','proper','rank',name));
    generator.(field)=local_generator_validation(draw.H_ref_sca_f,Hkd_test,Hj_test,f_axis,cir_common_shift_s);
    generator.(field).timing=tm; generator.(field).rank_name=name; generator.(field).rank_used=tm.rank_used;
end
[~,improper_meta]=sample_conditional_channel_vertical(model,32,seed.generator_base+10,struct('path','improper','rank','full'));
generator.improper_smoke=improper_meta;
selected=local_select_rank(generator,rank_names);
model.stats.rank_selected=generator.(local_rank_field(selected)).rank_used;
model.validation=struct('rank_candidates',generator,'rank_selected_by_heldout',selected, ...
    'heldout_reference','independent explicit kdomain + cached PE test set');

performance=struct();
counts=unique([100,1000,M_gen]); performance.sample_counts=counts; performance.sample_wall_s=zeros(size(counts));
performance.cir_wall_s=zeros(size(counts));
for kk=1:numel(counts)
    [draw,tm]=sample_conditional_channel_vertical(model,counts(kk),seed.generator_base+100+kk, ...
        struct('path','proper','rank',selected));
    performance.sample_wall_s(kk)=tm.elapsed_s;
    cc=build_channel_cir_vertical(draw.H_total_f,f_axis,struct('input_reference','direct_dsp'));
    performance.cir_wall_s(kk)=cc.ifft_elapsed_s;
end
performance.sample_per_channel_s=performance.sample_wall_s./counts;
performance.cached_joint_single_s=ensemble.timing.joint_train.per_realization_s;
performance.cache_build_s=cache.build_time_s;
performance.joint_factor_build_s=joint_model.total_build_time_s;
performance.training_ensemble_s=ensemble.timing.joint_train.total_s;
performance.stats_estimate_s=model.timing.estimate_s;
performance.build_total_s=cache.build_time_s+joint_model.total_build_time_s+ ...
    ensemble.timing.joint_train.total_s+model.timing.estimate_s;
if strcmp(mode,'full')
    performance.public=local_public_baseline(cfg,seed.public);
    sample_single=performance.sample_per_channel_s(counts==M_gen);
    performance.break_even_channels=performance.build_total_s/ ...
        max(performance.public.wall_s-sample_single,eps);
else
    performance.public=struct('wall_s',NaN); performance.break_even_channels=NaN;
end
performance.peak_memory_snapshot_bytes=max([joint_model.peak_memory_snapshot_bytes, ...
    ensemble.timing.joint_train.peak_memory_snapshot_bytes],[],'omitnan');

ratio=[ensemble.mapping_train.crop_variance_m2;ensemble.mapping_test.crop_variance_m2]./ ...
    [ensemble.mapping_train.full_variance_m2;ensemble.mapping_test.full_variance_m2];
mapping=struct('mean_ratio',mean(ratio),'ci95',mean(ratio)+[-1,1]*1.96*std(ratio)/sqrt(numel(ratio)));
pass=struct('proper_joint_not_rejected',~properness.joint.reject_proper_at_5pct_one_sided, ...
    'joint_pdp',receiver.joint.pdp_correlation>=0.9,'joint_lfm',receiver.joint.lfm_correlation>=0.9, ...
    'tail_window',receiver.kdomain_temporal.tail_to_window_ratio<0.8, ...
    'generator_pdp',generator.(local_rank_field(selected)).pdp_correlation>=0.9, ...
    'proper_path',true,'improper_path',improper_meta.rank_used>0, ...
    'component_sum',ensemble.component_error<=1e-10,'no_pe_for_generation',true);
pass.all=all(structfun(@(x)logical(x),pass));

model.timing=performance; model.validation.receiver=receiver; model.validation.properness=properness;
model.validation_run_meta=validation_run_meta;
model_file=fullfile(out_dir,['u5_conditional_channel_model_f64_' mode '.mat']);
schema_version=model.schema_version; phase_reference_meta=model.phase_reference_meta;
save(model_file,'model','schema_version','phase_reference_meta','validation_run_meta','-v7.3'); info=dir(model_file); performance.model_file_bytes=info.bytes;
result=struct('mode',mode,'config',cfg,'condition',condition,'L_train',L_train,'L_test',L_test, ...
    'M_gen',M_gen,'seeds',seed,'stats',stats,'properness',properness,'receiver',receiver, ...
    'generator',generator,'rank_selected',selected,'mapping',mapping,'performance',performance, ...
    'component_sum_max_abs_error',ensemble.component_error,'pass',pass, ...
    'delay_resolution_s',1/(f_axis(end)-f_axis(1)),'Tmax_s',1/mean(diff(f_axis)));
result.schema_version='2.0.0'; result.phase_reference_meta=phase_meta;
result.validation_run_meta=validation_run_meta;
result.phase_neutral_cache_audit=phase_neutral_cache_audit;
result_file=fullfile(out_dir,['u5_conditional_channel_validation_f64_' mode '.mat']);
schema_version=result.schema_version; phase_reference_meta=result.phase_reference_meta;
save(result_file,'result','schema_version','phase_reference_meta','validation_run_meta','-v7.3');
local_plots(result,out_dir); local_summary(result,fullfile(out_dir,['summary_f64_' mode '.txt']));
fprintf('F=64 %s complete. joint PDP %.4f, LFM %.4f, tail/T %.4f, selected %s, pass %d\n', ...
    mode,receiver.joint.pdp_correlation,receiver.joint.lfm_correlation, ...
    receiver.kdomain_temporal.tail_to_window_ratio,selected,pass.all);

function s=local_stats(H)
mu=mean(H,2); X=H-mu; L=size(H,2); C=(X*X')/(L-1); C=0.5*(C+C'); P=(X*X.')/(L-1);
d=sqrt(max(real(diag(C)),0)); R=C./max(d*d.',eps);
s=struct('mu',mu,'C',C,'P',P,'R',R,'adjacent_corr',diag(R,1), ...
    'eigenvalues',sort(max(real(eig(C)),0),'descend'), ...
    'properness_ratio',norm(P,'fro')/max(norm(C,'fro'),eps));
end

function m=local_compare(a,b,Ha,Hb,f,tref)
ta=local_temporal(Ha,f,tref); tb=local_temporal(Hb,f,tref);
m=struct('covariance_relative_error',norm(a.C-b.C,'fro')/max(norm(b.C,'fro'),eps), ...
    'correlation_matrix_relative_error',norm(a.R-b.R,'fro')/max(norm(b.R,'fro'),eps), ...
    'adjacent_correlation_rmse',sqrt(mean(abs(a.adjacent_corr-b.adjacent_corr).^2)), ...
    'mean_absolute_error',mean(abs(a.mu-b.mu)), ...
    'pdp_correlation',local_corr(ta.pdp,tb.pdp), ...
    'lfm_correlation',local_corr(local_lfm(Ha,f),local_lfm(Hb,f)));
end

function t=local_temporal(H,f,tref)
c=build_channel_cir_vertical(H,f,struct('input_reference','direct_dsp', ...
    'time_origin_shift_s',tref));
p=abs(c.h_tau).^2; pdp=mean(p,2); pdp=pdp/max(sum(pdp),eps);
tau=c.delay_axis_s; L=size(H,2); N=numel(tau); dt=c.delay_axis_spacing_s;
signed_tau=((-floor(N/2)):(ceil(N/2)-1)).'*dt;
mean_delay=zeros(L,1); rms_delay=zeros(L,1); tail=zeros(L,1);
for mm=1:L
    [~,ip]=max(p(:,mm)); center=floor(N/2)+1;
    q=circshift(p(:,mm),center-ip); q=q/max(sum(q),eps);
    mean_delay(mm)=sum(signed_tau.*q); rms_delay(mm)=sqrt(sum((signed_tau-mean_delay(mm)).^2.*q));
    tail(mm)=local_min_circular_energy_span(p(:,mm),0.99,dt);
end
[~,ip]=max(pdp);
t=struct('pdp',pdp,'delay_axis_s',tau,'peak_delay_s',tau(ip), ...
    'mean_delay_s',mean(mean_delay),'rms_delay_s',mean(rms_delay), ...
    'tail99_s',mean(tail),'tail99_definition','shortest circular contiguous interval containing 99% energy', ...
    'tail_to_window_ratio',mean(tail)/c.maximum_unambiguous_delay_s, ...
    'Tmax_s',c.maximum_unambiguous_delay_s,'physical_resolution_s',c.physical_delay_resolution_s);
end

function span=local_min_circular_energy_span(power,target,dt)
q=real(power(:)); q=q/max(sum(q),eps); N=numel(q); q2=[q;q]; best=N;
right=0; running=0;
for left=1:N
    while right<left+N-1 && running<target
        right=right+1; running=running+q2(right);
    end
    if running>=target, best=min(best,right-left); end
    running=running-q2(left);
    if right<left, right=left; running=0; end
end
span=best*dt;
end

function env=local_lfm(H,f)
fs=12000; N=512; t=(0:N-1).'/fs; duration=0.02; active=t<duration;
tx=zeros(N,1); tx(active)=exp(1i*pi*((f(end)-f(1))/duration)*(t(active)-duration/2).^2);
fbb=(-N/2:N/2-1).'*fs/N; Hbb=complex(zeros(N,size(H,2)));
for mm=1:size(H,2), Hbb(:,mm)=interp1(f-mean(f),H(:,mm),fbb,'linear',0); end
TX=fftshift(fft(tx)); rx=ifft(ifftshift(TX.*Hbb),[],1); mf=ifft(fft(rx).*conj(fft(tx)),[],1);
env=mean(abs(mf),2); env=env/max(env);
end

function v=local_generator_validation(H,Hkd,Hj,f,tref)
s=local_stats(H); sk=local_stats(Hkd); sj=local_stats(Hj); tk=local_temporal(Hkd,f,tref); t=local_temporal(H,f,tref);
se=sqrt(max(real(diag(sk.C)),0)/size(Hkd,2)); mean_z=abs(s.mu-sk.mu)./max(se,eps);
v=struct('C_vs_kdomain',norm(s.C-sk.C,'fro')/max(norm(sk.C,'fro'),eps), ...
    'C_vs_joint',norm(s.C-sj.C,'fro')/max(norm(sj.C,'fro'),eps), ...
    'P_error_normalized',norm(s.P-sk.P,'fro')/max(norm(sk.C,'fro'),eps), ...
    'absolute_mean_error',mean(abs(s.mu-sk.mu)), ...
    'mean_over_standard_error_median',median(mean_z),'mean_over_standard_error_max',max(mean_z), ...
    'pdp_correlation',local_corr(t.pdp,tk.pdp), ...
    'lfm_correlation',local_corr(local_lfm(H,f),local_lfm(Hkd,f)), ...
    'temporal',t,'coherence_bandwidth_hz',local_coherence_bandwidth(s.C,f), ...
    'distribution',local_distribution(H,Hkd));
end

function d=local_distribution(H,Hr)
ii=ceil(size(H,1)/2); z=H(ii,:); r=Hr(ii,:);
d=struct('ks_real',local_ks(real(z),real(r)),'ks_imag',local_ks(imag(z),imag(r)), ...
    'ks_magnitude',local_ks(abs(z),abs(r)),'skew_real',local_skew(real(z)), ...
    'skew_imag',local_skew(imag(z)),'kurtosis_real',local_kurt(real(z)), ...
    'kurtosis_imag',local_kurt(imag(z)),'K_factor',abs(mean(z))^2/max(var(z,1),eps), ...
    'reference_K_factor',abs(mean(r))^2/max(var(r,1),eps));
end

function bw=local_coherence_bandwidth(C,f)
d=sqrt(max(real(diag(C)),0)); R=abs(C./max(d*d.',eps)); F=numel(f); rho=zeros(F,1);
for lag=0:F-1, rho(lag+1)=mean(diag(R,lag)); end
idx=find(rho<0.5,1); if isempty(idx), bw=f(end)-f(1); else, bw=(idx-1)*mean(diff(f)); end
end
function r=local_corr(a,b), q=corrcoef(real(a(:)),real(b(:))); r=q(1,2); end
function v=local_ks(a,b)
a=sort(a(:)); b=sort(b(:)); x=sort([a;b]); Fa=arrayfun(@(q)mean(a<=q),x); Fb=arrayfun(@(q)mean(b<=q),x); v=max(abs(Fa-Fb));
end
function v=local_skew(x), q=x-mean(x); v=mean(q.^3)/max(mean(q.^2)^(3/2),eps); end
function v=local_kurt(x), q=x-mean(x); v=mean(q.^4)/max(mean(q.^2)^2,eps); end
function field=local_rank_field(name), field=['rank_' strrep(name,'.','_')]; end

function selected=local_select_rank(g,names)
valid=false(size(names)); err=inf(size(names));
for ii=1:numel(names), q=g.(local_rank_field(names{ii})); valid(ii)=q.pdp_correlation>=0.9&&q.lfm_correlation>=0.9; err(ii)=q.C_vs_kdomain; end
if any(valid), idx=find(valid); [~,jj]=min(err(valid)); selected=names{idx(jj)}; else, selected='full'; end
end

function p=local_public_baseline(cfg,sea_seed)
params=struct('f0',cfg.f_axis_hz.','enable_wideband',false,'c0',cfg.c0,'z_max',cfg.z_tx, ...
    'z_tx',cfg.z_tx,'z_rx',cfg.z_rx,'xw',cfg.xw,'yw',cfg.yw,'nx',cfg.nx,'ny',cfg.ny, ...
    'x_tx',0,'y_tx',0,'x_rx',0,'y_rx',0,'sigma_src_m',cfg.sigma_src_m, ...
    'stepz_lamb',cfg.stepz_lamb,'sponge_ratio',cfg.sponge_ratio, ...
    'alpha_max_np_per_m',cfg.alpha_max_np_per_m,'env_mode','uniform','enable_bubbles',false, ...
    'enable_surface_reflection',true,'surface_reflect_coeff',-1, ...
    'surface_boundary_model','kirchhoff_kdomain','surface_roughness_scale_mode','raw_pm', ...
    'sea_wind_speed',5,'sea_seed',sea_seed,'show_figures',false,'enforce_1_over_R',false, ...
    'save_mode','rx_only','use_gpu',false);
timer=tic; out=vertical_channel_model(params); p=struct('wall_s',toc(timer), ...
    'component_sum_error',max(abs(out.H_f-out.H_direct_f-out.H_reflect_f)));
end

function local_plots(r,out_dir)
fig=figure('Visible','off','Color','w'); hold on;
plot(r.receiver.kdomain_temporal.delay_axis_s*1e3,r.receiver.kdomain_temporal.pdp,'LineWidth',1.3);
plot(r.receiver.joint_temporal.delay_axis_s*1e3,r.receiver.joint_temporal.pdp,'LineWidth',1.3);
plot(r.receiver.independent_temporal.delay_axis_s*1e3,r.receiver.independent_temporal.pdp,'LineWidth',1.3);
for name={'full','99.9','99'}, q=r.generator.(local_rank_field(name{1})); plot(q.temporal.delay_axis_s*1e3,q.temporal.pdp,'--'); end
grid on; legend('kdomain','joint','independent','gen full','gen 99.9','gen 99'); xlabel('relative delay (ms)'); ylabel('PDP');
exportgraphics(fig,fullfile(out_dir,['f64_pdp_' r.mode '.png']),'Resolution',180); close(fig);
fig=figure('Visible','off','Color','w');
bar([r.generator.rank_full.C_vs_kdomain,r.generator.rank_99_9.C_vs_kdomain,r.generator.rank_99.C_vs_kdomain; ...
    r.generator.rank_full.pdp_correlation,r.generator.rank_99_9.pdp_correlation,r.generator.rank_99.pdp_correlation; ...
    r.generator.rank_full.lfm_correlation,r.generator.rank_99_9.lfm_correlation,r.generator.rank_99.lfm_correlation]);
set(gca,'XTickLabel',{'C error','PDP corr','LFM corr'}); legend('full','99.9%','99%'); grid on;
exportgraphics(fig,fullfile(out_dir,['f64_rank_comparison_' r.mode '.png']),'Resolution',180); close(fig);
end

function local_summary(r,path)
fid=fopen(path,'w'); c=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'U5 conditional receiver generator F64 %s\n',r.mode);
fprintf(fid,'F=%d df=%.9g Hz dt_phys=%.9g s Tmax=%.9g s train=%d test=%d Mgen=%d\n', ...
    numel(r.config.f_axis_hz),mean(diff(r.config.f_axis_hz)),r.delay_resolution_s,r.Tmax_s,r.L_train,r.L_test,r.M_gen);
fprintf(fid,'proper joint observed=%.6g null95=[%.6g %.6g] p=%.6g decision=%s\n', ...
    r.properness.joint.observed_ratio,r.properness.joint.null_interval_95, ...
    r.properness.joint.upper_tail_p_value,r.properness.joint.decision);
fprintf(fid,'joint C=%.6g R=%.6g adj=%.6g PDP=%.6g LFM=%.6g\n', ...
    r.receiver.joint.covariance_relative_error,r.receiver.joint.correlation_matrix_relative_error, ...
    r.receiver.joint.adjacent_correlation_rmse,r.receiver.joint.pdp_correlation,r.receiver.joint.lfm_correlation);
fprintf(fid,'tail99=%.9g s tail/T=%.6g selected=%s pass=%d\n', ...
    r.receiver.kdomain_temporal.tail99_s,r.receiver.kdomain_temporal.tail_to_window_ratio,r.rank_selected,r.pass.all);
for name={'full','99.9','99'}, q=r.generator.(local_rank_field(name{1})); fprintf(fid,'rank %s r=%d C=%.6g PDP=%.6g LFM=%.6g KSabs=%.6g\n', ...
    name{1},q.rank_used,q.C_vs_kdomain,q.pdp_correlation,q.lfm_correlation,q.distribution.ks_magnitude); end
fprintf(fid,'sample counts/times'); fprintf(fid,' %d:%.9g', [r.performance.sample_counts;r.performance.sample_wall_s]); fprintf(fid,'\n');
fprintf(fid,'public=%.9g cached_joint=%.9g build=%.9g break_even=%.6g model_bytes=%d\n', ...
    r.performance.public.wall_s,r.performance.cached_joint_single_s,r.performance.build_total_s, ...
    r.performance.break_even_channels,r.performance.model_file_bytes);
end
