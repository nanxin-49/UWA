%% Exact discrete-adjoint PE receiver-projection feasibility validation
clear; close all; clc;
root_dir=fileparts(fileparts(fileparts(mfilename('fullpath')))); addpath(root_dir); setup_vertical_project();
addpath(fullfile(root_dir,'scripts','reporting'));
out_dir=fullfile(root_dir,'results','validation','adjoint_pe_receiver_projection');
if ~exist(out_dir,'dir'), mkdir(out_dir); end
mode=lower(strtrim(getenv('ADJOINT_PE_VALIDATION_MODE')));
if isempty(mode), mode='full'; end
if ~ismember(mode,{'smoke','full'}), error('Mode must be smoke or full.'); end
validation_run_meta=pe_phase_release_run_meta_vertical(root_dir,struct( ...
    'frequency_axis_hz',linspace(4000,8000,64).', ...
    'seed_definition',struct('adjoint',710001,'projection',720001, ...
    'f9_ensemble',730001,'f64_ensemble',740001)));

fprintf('Adjoint PE receiver validation mode=%s\n',mode);
validation=struct('mode',mode,'created_at',char(datetime('now')), ...
    'scope','uniform CPU double; one nearest-grid receiver; no bubbles or Doppler');

%% Layer 1: exact adjoint and receiver projection at 4/6/8 kHz.
f3=[4000;6000;8000];
[cfg3,spec3]=local_case(f3,64,128);
cache3=build_cached_joint_kstat_pe_executor_vertical(cfg3,spec3);
[projection3,projection3_meta]=build_adjoint_receiver_projection_vertical(cache3,spec3);
adjoint=local_adjoint_tests(cache3,5,710001);
projection_cases=local_projection_cases(cache3,projection3,spec3,720001);
validation.adjoint=adjoint;
validation.projection_cases=projection_cases;
validation.projection_build=projection3_meta;
pass_exact=adjoint.max_relative_error<=1e-10 && ...
    projection_cases.max_relative_error<=1e-10 && ...
    projection_cases.legacy_runner_relative_error<=1e-12;
if ~pass_exact
    error('Exact adjoint/projection gate failed; statistical validation was not run.');
end
fprintf('Exact gate: adjoint %.3g, projection %.3g, runner %.3g\n', ...
    adjoint.max_relative_error,projection_cases.max_relative_error, ...
    projection_cases.legacy_runner_relative_error);

%% Layer 2: nontrivial PM 16^2 -> PE 8^2 dense/FFT baseline.
[cfg_dense,spec_dense]=local_case(f3,8,16);
cfg_dense.z_tx=10; cfg_dense.z_rx=1;
cache_dense=build_cached_joint_kstat_pe_executor_vertical(cfg_dense,spec_dense);
[projection_dense,~]=build_adjoint_receiver_projection_vertical(cache_dense,spec_dense);
[dense_stats,dense_meta]=contract_kstat_receiver_stats_vertical( ...
    spec_dense,projection_dense,struct('method','dense'));
[fft_dense_stats,fft_dense_meta]=contract_kstat_receiver_stats_vertical( ...
    spec_dense,projection_dense,struct('method','fft'));
dense_fft=struct( ...
    'C_relative_error',norm(dense_stats.C_H-fft_dense_stats.C_H,'fro')/max(norm(dense_stats.C_H,'fro'),eps), ...
    'P_relative_error',norm(dense_stats.P_H-fft_dense_stats.P_H,'fro')/max(norm(dense_stats.P_H,'fro'),eps), ...
    'dense_meta',dense_meta,'fft_meta',fft_dense_meta, ...
    'dense_stats',dense_stats,'fft_stats',fft_dense_stats);
validation.dense_fft=dense_fft;
if dense_fft.C_relative_error>1e-10 || dense_fft.P_relative_error>1e-10 || ...
        fft_dense_meta.embedding_inside_max_abs_error~=0 || ...
        fft_dense_meta.embedding_outside_max_abs~=0
    error('Dense/FFT or PM zero-embedding gate failed.');
end
fprintf('Dense/FFT gate: C %.3g, P %.3g\n', ...
    dense_fft.C_relative_error,dense_fft.P_relative_error);

%% Layer 3: F=9 wideband ensemble and analytic statistics.
f9=(4000:500:8000).';
[cfg9,spec9]=local_case(f9,64,128);
cache9=build_cached_joint_kstat_pe_executor_vertical(cfg9,spec9);
[projection9,projection9_meta]=build_adjoint_receiver_projection_vertical(cache9,spec9);
joint9=build_kirchhoff_kstat_joint_model_streaming_vertical( ...
    spec9,f9,cfg9.c0,cfg9.reflect_coeff);
if strcmp(mode,'full'), M9=4096; else, M9=64; end
[ensemble9,timing9]=local_joint_ensemble(cache9,projection9,joint9,M9,16,730000);
[analytic9,analytic9_meta]=contract_kstat_receiver_stats_vertical( ...
    spec9,projection9,struct('method','fft'));
statistics9=local_statistical_comparison( ...
    ensemble9.H_projection_fm,analytic9,cfg9,0);
statistics9.forward_projection_max_relative_error=ensemble9.max_relative_error;
statistics9.pass_forward_projection=ensemble9.max_relative_error<=1e-10;
statistics9.pass_split_floor=local_split_floor_pass(statistics9,M9>=512);
validation.f9=struct('sample_count',M9,'batch_size',16, ...
    'projection_build',projection9_meta,'joint_model_meta',local_model_meta(joint9), ...
    'ensemble_timing',timing9,'analytic_meta',analytic9_meta, ...
    'analytic_stats',analytic9,'comparison',statistics9);
if ~statistics9.pass_forward_projection || ~statistics9.pass_split_floor
    error('F=9 forward/projection or analytic split-floor gate failed.');
end
fprintf('F=9: projection %.3g, C/floor %.3f, P(C-norm)/floor %.3f\n', ...
    ensemble9.max_relative_error,statistics9.C_to_split_floor_ratio, ...
    statistics9.P_Cnorm_to_split_floor_ratio);

%% Layer 4: production F=64 scale and performance matrix.
if strcmp(mode,'full')
    [cache64,joint64,spec64]=local_load_or_build_f64(root_dir);
    [projection64,projection64_meta]=build_adjoint_receiver_projection_vertical(cache64,spec64);
    [ensemble64,timing64]=local_joint_ensemble(cache64,projection64,joint64,512,8,740000);
    [analytic64,analytic64_meta]=contract_kstat_receiver_stats_vertical( ...
        spec64,projection64,struct('method','fft'));
    statistics64=local_statistical_comparison(ensemble64.H_projection_fm,analytic64, ...
        cache64.cfg,0);
    statistics64.forward_projection_max_relative_error=ensemble64.max_relative_error;
    statistics64.pass_forward_projection=ensemble64.max_relative_error<=1e-10;
    statistics64.pass_split_floor=local_split_floor_pass(statistics64,true);
    validation.f64=struct('sample_count',512,'batch_size',8, ...
        'projection_build',projection64_meta,'joint_model_meta',local_model_meta(joint64), ...
        'ensemble_timing',timing64,'analytic_meta',analytic64_meta, ...
        'analytic_stats',analytic64,'comparison',statistics64);
    validation.performance=local_performance_matrix(cfg9,spec9,cache9,projection9, ...
        joint9,cache64,projection64,joint64,spec64);
    if ~statistics64.pass_forward_projection || ~statistics64.pass_split_floor
        error('F=64 forward/projection or analytic split-floor gate failed.');
    end
else
    validation.f64=struct('executed',false,'reason','smoke mode');
    validation.performance=local_smoke_performance(cache9,projection9,joint9);
end

validation.pass=struct('exact_adjoint_projection',pass_exact, ...
    'dense_fft',dense_fft.C_relative_error<=1e-10&&dense_fft.P_relative_error<=1e-10, ...
    'f9',statistics9.pass_forward_projection&&statistics9.pass_split_floor, ...
    'f64',strcmp(mode,'smoke')||(...
        validation.f64.comparison.pass_forward_projection&&validation.f64.comparison.pass_split_floor));
validation.pass.all=all(structfun(@(x)logical(x),validation.pass));
validation.schema_version='2.0.0';
validation.phase_reference_meta=projection3.phase_reference_meta;
validation.validation_run_meta=validation_run_meta;
schema_version=validation.schema_version;
phase_reference_meta=validation.phase_reference_meta;
save(fullfile(out_dir,['adjoint_pe_receiver_projection_' mode '.mat']),'validation', ...
    'schema_version','phase_reference_meta','validation_run_meta','-v7.3');
local_make_figures(validation,out_dir,mode);
local_write_summary(validation,fullfile(out_dir,['summary_' mode '.txt']));
fprintf('Saved adjoint PE validation outputs to %s\n',out_dir);

function [cfg,spec]=local_case(f_axis,pe_n,pm_n)
dx=0.390625; pe_w=pe_n*dx; pm_w=pm_n*dx;
cfg=struct('f_axis_hz',f_axis(:),'c0',1500,'xw',pe_w,'yw',pe_w, ...
    'nx',pe_n,'ny',pe_n,'x_tx',0,'y_tx',0,'z_tx',100, ...
    'x_rx',0,'y_rx',0,'z_rx',3,'sigma_src_m',0.4, ...
    'stepz_lamb',0.5,'sponge_ratio',0.12,'alpha_max_np_per_m',0.15, ...
    'reflect_coeff',-1,'env_mode','uniform','enable_bubbles',false);
spec=raw_pm_spectrum_grid_vertical(5,pm_n,pm_n,pm_w,pm_w);
end

function result=local_adjoint_tests(cache,n_seed,seed0)
F=numel(cache.f_axis_hz); rel=zeros(F,n_seed); normalized=rel; absolute=rel;
for ii=1:F
    for ss=1:n_seed
        rng(seed0+100*ii+ss,'twister');
        x=complex(randn(cache.cfg.ny,cache.cfg.nx),randn(cache.cfg.ny,cache.cfg.nx));
        y=complex(randn(cache.cfg.ny,cache.cfg.nx),randn(cache.cfg.ny,cache.cfg.nx));
        Ax=apply_forward_surface_to_receiver_vertical(cache,x,ii);
        AHy=apply_adjoint_receiver_to_surface_vertical(cache,y,ii);
        lhs=sum(conj(Ax(:)).*y(:)); rhs=sum(conj(x(:)).*AHy(:));
        absolute(ii,ss)=abs(lhs-rhs);
        rel(ii,ss)=absolute(ii,ss)/max([abs(lhs),abs(rhs),eps]);
        normalized(ii,ss)=absolute(ii,ss)/max(norm(Ax(:))*norm(y(:))+norm(x(:))*norm(AHy(:)),eps);
    end
end
result=struct('f_axis_hz',cache.f_axis_hz(:),'seed_count',n_seed, ...
    'relative_error',rel,'normalized_error',normalized,'absolute_error',absolute, ...
    'max_relative_error',max(rel(:)),'max_normalized_error',max(normalized(:)));
end

function result=local_projection_cases(cache,projection,spec,seed)
F=numel(cache.f_axis_hz); errors=zeros(4,F); labels={'random','flat','explicit','joint'};
rng(seed,'twister'); random_fields=complex(randn(cache.cfg.ny,cache.cfg.nx,F),randn(cache.cfg.ny,cache.cfg.nx,F));
[eta_pm,~]=sample_raw_pm_surface_vertical(spec,seed+1);
eta_pe=eta_pm(cache.pm_mapping.iy,cache.pm_mapping.ix);
joint=build_kirchhoff_kstat_joint_model_streaming_vertical( ...
    spec,cache.f_axis_hz,cache.cfg.c0,cache.cfg.reflect_coeff);
[delta_joint,~]=sample_kirchhoff_kstat_factor_model_vertical(joint,1,seed+2,'joint');
for ii=1:F
    alpha=4*pi*cache.f_axis_hz(ii)/cache.cfg.c0;
    fields=cell(4,1);
    fields{1}=random_fields(:,:,ii);
    fields{2}=cache.R_coh_f(ii).*cache.incident_surface_xy_f(:,:,ii);
    fields{3}=cache.cfg.reflect_coeff.*exp(1i*alpha.*eta_pe).*cache.incident_surface_xy_f(:,:,ii);
    fields{4}=double(delta_joint(cache.pm_mapping.iy,cache.pm_mapping.ix,ii,1)).* ...
        cache.incident_surface_xy_f(:,:,ii);
    for cc=1:4
        errors(cc,ii)=local_receiver_projection_error(cache,projection,fields{cc},ii);
    end
end
[legacy,new_result]=local_legacy_runner_comparison(cache,delta_joint);
legacy_error=norm(legacy.H_ref_sca_fm-new_result.H_ref_sca_reduced_fm,'fro')/max(norm(legacy.H_ref_sca_fm,'fro'),eps);
result=struct('labels',{labels},'relative_error',errors, ...
    'max_relative_error',max(errors(:)),'legacy_runner_relative_error',legacy_error);
end

function error_value=local_receiver_projection_error(cache,projection,field,ii)
out=apply_forward_surface_to_receiver_vertical(cache,double(field),ii);
h_forward=out(cache.iy_rx,cache.ix_rx);
q=projection.q_surface_xy_f(:,:,ii);
h_projection=sum(conj(q(:)).*double(field(:)));
error_value=abs(h_forward-h_projection)/max(abs(h_forward),eps);
end

function [legacy,new_result]=local_legacy_runner_comparison(cache,delta_pm)
F=numel(cache.f_axis_hz); M=size(delta_pm,4); delta=delta_pm(cache.pm_mapping.iy,cache.pm_mapping.ix,:,:);
H=complex(zeros(F,M));
for ii=1:F
    pages=double(reshape(delta(:,:,ii,:),cache.cfg.ny,cache.cfg.nx,M)).*cache.incident_surface_xy_f(:,:,ii);
    fr=cache.surface_rx_fr(:,:,ii); screen=cache.surface_rx_screen(:,:,ii); psi_k=fft2(pages);
    for jj=1:cache.surface_rx_nstep(ii), psi_k=fr.*fft2(screen.*ifft2(fr.*psi_k)); end
    psi=ifft2(psi_k); H(ii,:)=reshape(psi(cache.iy_rx,cache.ix_rx,:),1,M);
end
legacy=struct('H_ref_sca_fm',H);
new_result=run_cached_joint_kstat_pe_executor_vertical(cache,delta_pm);
end

function [ensemble,timing]=local_joint_ensemble(cache,projection,model,M,batch,seed0)
F=model.F; Hf=complex(zeros(F,M)); Hp=Hf; generation_s=0; forward_s=0; projection_s=0;
max_error=0; offset=0; batch_index=0;
while offset<M
    count=min(batch,M-offset); batch_index=batch_index+1;
    [delta,gm]=sample_kirchhoff_kstat_factor_model_vertical(model,count,seed0+batch_index,'joint');
    [forward,fm]=run_cached_joint_kstat_pe_executor_vertical(cache,delta);
    [projected,pm]=run_adjoint_receiver_projection_vertical(projection,delta);
    idx=offset+(1:count); Hf(:,idx)=forward.H_ref_sca_fm; Hp(:,idx)=projected.H_ref_sca_fm;
    difference=forward.H_ref_sca_fm-projected.H_ref_sca_fm;
    max_error=max(max_error,max(abs(difference(:))./max(abs(forward.H_ref_sca_fm(:)),eps)));
    generation_s=generation_s+gm.elapsed_s; forward_s=forward_s+fm.elapsed_s; projection_s=projection_s+pm.elapsed_s;
    offset=offset+count;
end
ensemble=struct('H_forward_fm',Hf,'H_projection_fm',Hp,'max_relative_error',max_error);
timing=struct('generation_s',generation_s,'cached_forward_s',forward_s, ...
    'projection_s',projection_s,'per_realization_generation_s',generation_s/M, ...
    'per_realization_cached_forward_s',forward_s/M, ...
    'per_realization_projection_s',projection_s/M, ...
    'online_projection_speedup',forward_s/max(projection_s,eps));
end

function result=local_statistical_comparison(H,analytic,cfg,reference_delay)
M=size(H,2); split=floor(M/2); full=local_sample_stats(H);
one=local_sample_stats(H(:,1:split)); two=local_sample_stats(H(:,split+1:2*split));
Cden=max(norm(analytic.C_H,'fro'),eps); Pden=max(norm(analytic.P_H,'fro'),eps);
C_error=norm(full.C-analytic.C_H,'fro')/Cden;
P_error=norm(full.P-analytic.P_H,'fro')/Pden;
P_Cnorm=norm(full.P-analytic.P_H,'fro')/Cden;
C_floor=norm(one.C-two.C,'fro')/Cden;
P_floor=norm(one.P-two.P,'fro')/Pden;
P_Cfloor=norm(one.P-two.P,'fro')/Cden;
sample_power=mean(abs(H).^2,2); power_error=norm(sample_power-analytic.E_abs_H2_f)/max(norm(analytic.E_abs_H2_f),eps);
cir_basis=build_channel_cir_vertical(eye(size(H,1)),analytic.f_axis_hz, ...
    struct('input_reference','direct_dsp','time_origin_shift_s',reference_delay));
B=cir_basis.h_tau;
pdp_analytic=real(diag(B*analytic.C_H*B')); pdp_sample=mean(abs(B*H).^2,2);
[R,G]=local_lfm_operators(analytic.f_axis_hz);
lfm_analytic=real(diag(R*analytic.C_H*R')); lfm_sample=mean(abs(R*H).^2,2);
mf_analytic=real(diag(G*analytic.C_H*G')); mf_sample=mean(abs(G*H).^2,2);
standard_error=sqrt(max(real(diag(full.C)),0)/M);
mean_z=max(abs(full.mu)./max(standard_error,eps));
result=struct('sample_count',M,'split_count',split, ...
    'C_relative_error',C_error,'P_relative_error',P_error,'P_C_normalized_error',P_Cnorm, ...
    'C_split_floor',C_floor,'P_split_floor',P_floor,'P_C_normalized_split_floor',P_Cfloor, ...
    'C_to_split_floor_ratio',C_error/max(C_floor,eps), ...
    'P_to_split_floor_ratio',P_error/max(P_floor,eps), ...
    'P_Cnorm_to_split_floor_ratio',P_Cnorm/max(P_Cfloor,eps), ...
    'power_relative_error',power_error,'mean_z_max',mean_z, ...
    'pdp_correlation',local_corr(pdp_analytic,pdp_sample), ...
    'lfm_correlation',local_corr(lfm_analytic,lfm_sample), ...
    'matched_filter_correlation',local_corr(mf_analytic,mf_sample), ...
    'pdp_analytic',pdp_analytic,'pdp_sample',pdp_sample, ...
    'lfm_analytic',lfm_analytic,'lfm_sample',lfm_sample, ...
    'matched_filter_analytic',mf_analytic,'matched_filter_sample',mf_sample, ...
    'sample_stats',full,'cfg',cfg);
end

function tf=local_split_floor_pass(result,enforce)
if ~enforce, tf=true; return; end
tf=result.C_to_split_floor_ratio<=1.25 && ...
    result.P_Cnorm_to_split_floor_ratio<=1.25 && result.mean_z_max<=4;
end

function stats=local_sample_stats(H)
mu=mean(H,2); X=H-mu; L=size(H,2);
C=(X*X')/(L-1); C=0.5*(C+C');
P=(X*X.')/(L-1); P=0.5*(P+P.');
stats=struct('mu',mu,'C',C,'P',P,'pseudo_to_cov_ratio',norm(P,'fro')/max(norm(C,'fro'),eps));
end

function [R,G]=local_lfm_operators(f)
fs=12000; N=512; t=(0:N-1).'/fs; duration=.02; active=t<duration; tx=zeros(N,1);
tx(active)=exp(1i*pi*((f(end)-f(1))/duration)*(t(active)-duration/2).^2);
fb=(-N/2:N/2-1).'*fs/N; F=numel(f); R=complex(zeros(N,F)); G=R; TX=fftshift(fft(tx));
for ii=1:F
    basis=zeros(F,1); basis(ii)=1;
    Hb=interp1(f-mean(f),basis,fb,'linear',0);
    R(:,ii)=ifft(ifftshift(TX.*Hb));
    G(:,ii)=ifft(fft(R(:,ii)).*conj(fft(tx)));
end
end

function r=local_corr(a,b)
cc=corrcoef(real(a(:)),real(b(:))); r=cc(1,2);
end

function meta=local_model_meta(model)
drop={'factor_cells','factor_basis','sqrt_q_independent'}; present=intersect(drop,fieldnames(model));
meta=rmfield(model,present);
end

function [cache,model,spec]=local_load_or_build_f64(root_dir)
folder=fullfile(root_dir,'results','validation','u5_conditional_channel_f64');
cache_file=fullfile(folder,'f64_joint_pe_cache.mat'); model_file=fullfile(folder,'f64_joint_streaming_model_u5.mat');
if exist(cache_file,'file')&&exist(model_file,'file')
    c=load(cache_file,'cache'); m=load(model_file,'joint_model'); cache=c.cache; model=m.joint_model;
    spec=raw_pm_spectrum_grid_vertical(5,cache.pm_mapping.pm_grid(2), ...
        cache.pm_mapping.pm_grid(1),cache.pm_mapping.pm_grid(2)*cache.pm_mapping.dx_m, ...
        cache.pm_mapping.pm_grid(1)*cache.pm_mapping.dx_m);
else
    f=linspace(4000,8000,64).'; [cfg,spec]=local_case(f,128,256);
    cache=build_cached_joint_kstat_pe_executor_vertical(cfg,spec);
    model=build_kirchhoff_kstat_joint_model_streaming_vertical(spec,f,cfg.c0,cfg.reflect_coeff);
end
end

function result=local_performance_matrix(cfg9,spec9,cache9,projection9,model9,cache64,projection64,model64,spec64)
one9=sample_kirchhoff_kstat_factor_model_vertical(model9,1,750001,'joint');
one64=sample_kirchhoff_kstat_factor_model_vertical(model64,1,750002,'joint');
result=struct();
result.f9=local_timing_repeats(cache9,projection9,one9,5);
result.f64=local_timing_repeats(cache64,projection64,one64,5);
f32=linspace(4000,8000,32).'; [cfg32,spec32]=local_case(f32,64,128);
cache32=build_cached_joint_kstat_pe_executor_vertical(cfg32,spec32);
[projection32,pmeta32]=build_adjoint_receiver_projection_vertical(cache32,spec32);
[~,smeta32]=contract_kstat_receiver_stats_vertical(spec32,projection32,struct('method','fft'));
result.f32_small=struct('projection_build',pmeta32,'analytic_stats',smeta32);
result.f9.projection_array_bytes=projection9.build_meta.projection_array_bytes;
result.f64.projection_array_bytes=projection64.build_meta.projection_array_bytes;
result.grid_note='F9/F32 PE64-PM128; production F64 PE128-PM256';
result.cfg9=rmfield(cfg9,intersect({'f_axis_hz'},fieldnames(cfg9))); %#ok<STRNU>
result.pm9_grid=[spec9.ny,spec9.nx]; result.pm64_grid=[spec64.ny,spec64.nx];
end

function result=local_smoke_performance(cache,projection,model)
one=sample_kirchhoff_kstat_factor_model_vertical(model,1,760001,'joint');
result=struct('f9',local_timing_repeats(cache,projection,one,3), ...
    'note','Smoke mode omits F32/F64 scaling.');
end

function result=local_timing_repeats(cache,projection,delta,repeats)
run_cached_joint_kstat_pe_executor_vertical(cache,delta);
run_adjoint_receiver_projection_vertical(projection,delta);
tf=zeros(repeats,1); tp=tf;
for ii=1:repeats
    t=tic; run_cached_joint_kstat_pe_executor_vertical(cache,delta); tf(ii)=toc(t);
    t=tic; run_adjoint_receiver_projection_vertical(projection,delta); tp(ii)=toc(t);
end
result=struct('cached_forward_median_s',median(tf),'projection_median_s',median(tp), ...
    'online_speedup',median(tf)/max(median(tp),eps),'repeats',repeats);
end

function local_make_figures(v,out_dir,mode)
plot_adjoint_pe_receiver_projection_validation_vertical(v,out_dir,mode);
end

function local_write_summary(v,path)
fid=fopen(path,'w'); cleaner=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'Adjoint PE receiver projection validation (%s)\n',v.mode);
fprintf(fid,'adjoint max relative error %.12g\n',v.adjoint.max_relative_error);
fprintf(fid,'projection max relative error %.12g\n',v.projection_cases.max_relative_error);
fprintf(fid,'legacy/shared runner relative error %.12g\n',v.projection_cases.legacy_runner_relative_error);
fprintf(fid,'dense/FFT C relative error %.12g\n',v.dense_fft.C_relative_error);
fprintf(fid,'dense/FFT P relative error %.12g\n',v.dense_fft.P_relative_error);
fprintf(fid,'F9 samples %d, projection %.12g, C/floor %.6g, P(Cnorm)/floor %.6g\n', ...
    v.f9.sample_count,v.f9.comparison.forward_projection_max_relative_error, ...
    v.f9.comparison.C_to_split_floor_ratio,v.f9.comparison.P_Cnorm_to_split_floor_ratio);
if isfield(v.f64,'comparison')
    fprintf(fid,'F64 samples %d, projection %.12g, C/floor %.6g, P(Cnorm)/floor %.6g\n', ...
        v.f64.sample_count,v.f64.comparison.forward_projection_max_relative_error, ...
        v.f64.comparison.C_to_split_floor_ratio,v.f64.comparison.P_Cnorm_to_split_floor_ratio);
end
fprintf(fid,'all pass %d\n',v.pass.all);
end
