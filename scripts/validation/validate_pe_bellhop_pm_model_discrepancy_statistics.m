function study = validate_pe_bellhop_pm_model_discrepancy_statistics(overrides)
%VALIDATE_PE_BELLHOP_PM_MODEL_DISCREPANCY_STATISTICS Stage 4 PM study.
% Validation-only seed statistics for the frozen PE/Bellhop PM comparison.
if nargin < 1 || isempty(overrides), overrides = struct(); end
if ~isstruct(overrides) || ~isscalar(overrides), error('overrides must be a scalar struct.'); end
script_dir = fileparts(mfilename('fullpath')); root = fileparts(fileparts(script_dir));
addpath(root); addpath(script_dir); addpath(fullfile(script_dir,'support')); setup_vertical_project();
cfg = local_config(root,overrides); if ~exist(cfg.output_dir,'dir'), mkdir(cfg.output_dir); end
common = cfg.common_overrides; common.seeds=cfg.seeds32; common.minimum_seed_count=32;
common.ensemble_mode='amplitude'; common.output_dir=cfg.output_dir; common.case_cache_dir=cfg.case_cache_dir;
common.report_path=fullfile(cfg.output_dir,'stage3_engine_intermediate_report.md');
m32 = validate_pe_bellhop_pm_ensemble(common); rows=m32.rows(:);
running32=local_running_summary(rows,[8 16 24 32]); boot32=local_bootstrap(rows,cfg.bootstrap_repetitions,cfg.bootstrap_seed);
gates32=local_convergence_gates(running32,boot32); extended=false; m64=struct();
if cfg.auto_extend_to_64 && ~gates32.all
    ext=common; ext.seeds=cfg.seeds64(33:end); ext.minimum_seed_count=1;
    m64=validate_pe_bellhop_pm_ensemble(ext); rows=[rows;m64.rows(:)]; extended=true;
end
counts=[8 16 24 32]; if extended, counts=[counts 48 64]; end
running=local_running_summary(rows,counts); boot=local_bootstrap(rows,cfg.bootstrap_repetitions,cfg.bootstrap_seed);
gates=local_convergence_gates(running,boot); correlations=local_correlations(rows,cfg.bootstrap_repetitions,cfg.bootstrap_seed+17);
bins=local_roughness_bins(rows); outliers=local_outlier_audit(rows,cfg.outlier_count); guards=local_guard_summary(rows,cfg);
dim=local_dimensionality_reference(cfg.stage1_dimensionality_file); classification=local_classification(guards,gates);
study=struct('schema_version','1.0.0','stage','4_pm_model_discrepancy_statistics','classification',classification, ...
    'config',cfg,'stage3b_32',m32,'stage3b_64',m64,'rows',rows,'running',running,'convergence_gates',gates, ...
    'bootstrap',boot,'correlations',correlations,'roughness_bins',bins,'outliers',outliers,'guards',guards, ...
    'dimensionality_reference',dim,'extended_to_m64',extended);
local_write_outputs(study); local_make_figures(study); local_write_report(study,cfg.report_path);
save(fullfile(cfg.output_dir,'result.mat'),'study','-v7.3');
if cfg.fail_on_guard && strcmp(classification,'NUMERICAL_OR_APPLICABILITY_LIMIT')
    error('Stage 4 PM study hit a numerical/applicability limit; see %s.',cfg.report_path);
end
end

function cfg=local_config(root,o)
cfg=struct('output_dir',fullfile(root,'results','validation','pe_bellhop_pm_model_discrepancy_statistics'), ...
 'report_path',fullfile(root,'reports','pe_bellhop_pm_model_discrepancy_statistics_report.md'), ...
 'case_cache_dir',fullfile(root,'results','validation','pe_bellhop_pm_amplitude_ensemble','cases'), ...
 'stage1_dimensionality_file',fullfile(root,'results','validation','pe_bellhop_pm_stage1_dimensionality','stage1b_dimensionality_summary.csv'), ...
 'seeds32',260001:260032,'seeds64',260001:260064,'bootstrap_repetitions',2000,'bootstrap_seed',42032, ...
 'outlier_count',5,'auto_extend_to_64',true,'fail_on_guard',true,'pe_edge5_limit',1e-4, ...
 'wall_residual_limit_m',1e-9,'phase_jump_limit_rad',1e-10,'state_limit',1e-12,'grazing_mu_threshold',.1, ...
 'common_overrides',struct('frequency_hz',4000,'c0_mps',1500,'z_tx_m',100,'z_rx_m',3,'sigma_src_m',.3, ...
  'xw_m',192.1875,'nx',984,'step_m',.05,'profile_count',4097,'beam_count',5001,'angle_limits_deg',[-15 15], ...
  'domain_half_depth_m',1000,'source_pattern_clip_db',-120,'source_geometry','X', ...
  'mapped_receiver_range_m',103,'selector_tolerance_m',1e-4, ...
  'wall_seed',260001,'reference_coeff_file',fullfile(root,'results','validation','bellhop_internal_pm_fixed_realization','fixed_pm_fourier_coefficients.csv'), ...
  'stage2_mat',fullfile(root,'results','validation','pe_bellhop_pm_frequency_extension','frequency_extension.mat'), ...
  'flat_validation_exe',fullfile(root,'results','validation','bellhop_internal_flat_wall_poc','bin','bellhop_iwall_flat_2020.exe'), ...
  'pm_validation_exe',fullfile(root,'results','validation','bellhop_internal_pm_wall_poc','bin','bellhop_iwall_pm_2020.exe'),'fail_on_check',true));
names=fieldnames(o);for ii=1:numel(names),if ~isfield(cfg,names{ii}),error('Unknown Stage 4 override: %s.',names{ii});end;cfg.(names{ii})=o.(names{ii});end
end

function s=local_running_summary(rows,counts)
s=repmat(struct('sample_count',NaN,'delta_tl_mean_db',NaN,'delta_tl_std_db',NaN,'delta_tl_median_db',NaN,'delta_tl_p05_db',NaN,'delta_tl_p95_db',NaN,'pe_power_mean',NaN,'bh_power_mean',NaN,'power_difference_mean',NaN,'phase_circular_mean_rad',NaN,'phase_circular_std_rad',NaN,'phase_resultant_length',NaN),numel(counts),1);
for ii=1:numel(counts)
 q=rows(1:counts(ii));dt=[q.delta_tl_db];ph=[q.delta_phase_rad];z=mean(exp(1i*ph));s(ii).sample_count=counts(ii);s(ii).delta_tl_mean_db=mean(dt);s(ii).delta_tl_std_db=std(dt,0,2);s(ii).delta_tl_median_db=median(dt);s(ii).delta_tl_p05_db=local_quantile(dt,.05);s(ii).delta_tl_p95_db=local_quantile(dt,.95);s(ii).pe_power_mean=mean([q.G_pe_power]);s(ii).bh_power_mean=mean([q.G_bellhop_power]);s(ii).power_difference_mean=mean([q.G_pe_power]-[q.G_bellhop_power]);s(ii).phase_circular_mean_rad=angle(z);s(ii).phase_circular_std_rad=sqrt(max(0,-2*log(max(abs(z),realmin))));s(ii).phase_resultant_length=abs(z);
end
end

function g=local_convergence_gates(running,boot)
g=struct('mean_delta_tl_change_24_to_32_db',NaN,'mean_delta_tl_pass',false,'pe_power_relative_change_24_to_32',NaN,'bh_power_relative_change_24_to_32',NaN,'power_pass',false,'phase_change_24_to_32_rad',NaN,'phase_pass',false,'resultant_change_24_to_32',NaN,'bootstrap_delta_tl_half_width_db',boot.mean_delta_tl_db.half_width,'bootstrap_pass',boot.mean_delta_tl_db.half_width<=.25,'all',false);
i24=find([running.sample_count]==24,1);i32=find([running.sample_count]==32,1);if isempty(i24)||isempty(i32),return;end;a=running(i24);b=running(i32);g.mean_delta_tl_change_24_to_32_db=abs(b.delta_tl_mean_db-a.delta_tl_mean_db);g.mean_delta_tl_pass=g.mean_delta_tl_change_24_to_32_db<=.10;g.pe_power_relative_change_24_to_32=abs(b.pe_power_mean-a.pe_power_mean)/max(abs(a.pe_power_mean),realmin);g.bh_power_relative_change_24_to_32=abs(b.bh_power_mean-a.bh_power_mean)/max(abs(a.bh_power_mean),realmin);g.power_pass=g.pe_power_relative_change_24_to_32<=.05&&g.bh_power_relative_change_24_to_32<=.05;g.phase_change_24_to_32_rad=abs(angle(exp(1i*(b.phase_circular_mean_rad-a.phase_circular_mean_rad))));g.phase_pass=g.phase_change_24_to_32_rad<=.20;g.resultant_change_24_to_32=abs(b.phase_resultant_length-a.phase_resultant_length);g.all=g.mean_delta_tl_pass&&g.power_pass&&g.phase_pass&&g.bootstrap_pass;
end

function b=local_bootstrap(rows,B,seed)
rng_state=rng;cleanup=onCleanup(@()rng(rng_state));rng(seed,'twister');n=numel(rows);ix=randi(n,B,n);dt=[rows.delta_tl_db];pe=[rows.G_pe_power];bh=[rows.G_bellhop_power];pd=pe-bh;ph=[rows.delta_phase_rad];v=struct('mean_delta_tl_db',zeros(B,1),'median_delta_tl_db',zeros(B,1),'mean_pe_power',zeros(B,1),'mean_bh_power',zeros(B,1),'mean_power_difference',zeros(B,1),'circular_mean_phase_rad',zeros(B,1),'resultant_length',zeros(B,1));
for ii=1:B,j=ix(ii,:);z=mean(exp(1i*ph(j)));v.mean_delta_tl_db(ii)=mean(dt(j));v.median_delta_tl_db(ii)=median(dt(j));v.mean_pe_power(ii)=mean(pe(j));v.mean_bh_power(ii)=mean(bh(j));v.mean_power_difference(ii)=mean(pd(j));v.circular_mean_phase_rad(ii)=angle(z);v.resultant_length(ii)=abs(z);end
names=fieldnames(v);b=struct('repetitions',B,'seed',seed);for ii=1:numel(names),x=v.(names{ii});b.(names{ii})=struct('estimate',local_boot_estimate(names{ii},dt,pe,bh,pd,ph),'ci95_low',local_quantile(x,.025),'ci95_high',local_quantile(x,.975),'half_width',.5*(local_quantile(x,.975)-local_quantile(x,.025)));end
end

function x=local_boot_estimate(key,dt,pe,bh,pd,ph)
switch key
 case 'mean_delta_tl_db',x=mean(dt); case 'median_delta_tl_db',x=median(dt); case 'mean_pe_power',x=mean(pe); case 'mean_bh_power',x=mean(bh); case 'mean_power_difference',x=mean(pd); case 'circular_mean_phase_rad',x=angle(mean(exp(1i*ph))); case 'resultant_length',x=abs(mean(exp(1i*ph)));
end
end

function q=local_quantile(x,p)
x=sort(double(x(:)));if isempty(x),q=NaN;return;end;pos=1+(numel(x)-1)*p;lo=floor(pos);hi=ceil(pos);if lo==hi,q=x(lo);else,q=x(lo)+(pos-lo)*(x(hi)-x(lo));end
end

function c=local_correlations(rows,B,seed)
predictors={'surface_rms_eta_m','profile_rms_slope','profile_max_slope','profile_rms_curvature_per_m','profile_max_curvature_per_m','profile_min_radius_m','two_k_sigma_eta','kmax_over_k','min_mu','mean_mu','hit_curvature_rms','hit_curvature_p95'};
responses={'delta_tl_db','complex_relative_error','power_difference','delta_power_db'};
rng_state=rng;cleanup=onCleanup(@()rng(rng_state));rng(seed,'twister');c=repmat(struct('predictor','','response','','pearson_r',NaN,'pearson_p',NaN,'pearson_ci95_low',NaN,'pearson_ci95_high',NaN,'spearman_rho',NaN,'spearman_p',NaN,'spearman_ci95_low',NaN,'spearman_ci95_high',NaN),0,1);
for ip=1:numel(predictors)
 x=[rows.(predictors{ip})].';
 for ir=1:numel(responses)
  if strcmp(responses{ir},'power_difference'),y=[rows.G_pe_power].'-[rows.G_bellhop_power].';else,y=[rows.(responses{ir})].';end
  [pr,pp]=local_corr(x,y,false);[sr,sp]=local_corr(x,y,true);[pci,sci]=local_corr_boot(x,y,B);
  c(end+1)=struct('predictor',predictors{ip},'response',responses{ir},'pearson_r',pr,'pearson_p',pp,'pearson_ci95_low',pci(1),'pearson_ci95_high',pci(2),'spearman_rho',sr,'spearman_p',sp,'spearman_ci95_low',sci(1),'spearman_ci95_high',sci(2)); %#ok<AGROW>
 end
 ph=[rows.delta_phase_rad].'; phase_labels={'phase_cos','phase_sin'}; phase_values={cos(ph),sin(ph)};
 for ik=1:2
  label=phase_labels{ik};y=phase_values{ik};
  [pr,pp]=local_corr(x,y,false);[sr,sp]=local_corr(x,y,true);[pci,sci]=local_corr_boot(x,y,B);
  c(end+1)=struct('predictor',predictors{ip},'response',label,'pearson_r',pr,'pearson_p',pp,'pearson_ci95_low',pci(1),'pearson_ci95_high',pci(2),'spearman_rho',sr,'spearman_p',sp,'spearman_ci95_low',sci(1),'spearman_ci95_high',sci(2)); %#ok<AGROW>
 end
end
end

function [r,p]=local_corr(x,y,ranked)
x=x(:);y=y(:);if ranked,x=local_rank(x);y=local_rank(y);end;x=x-mean(x);y=y-mean(y);den=sqrt(sum(x.^2)*sum(y.^2));if den<=eps,r=NaN;p=NaN;return;end;r=sum(x.*y)/den;n=numel(x);if n<=2||abs(r)>=1,p=0;return;end;t=abs(r)*sqrt((n-2)/max(1-r^2));p=betainc((n-2)/((n-2)+t^2),(n-2)/2,.5);
end

function ranks=local_rank(x)
[s,ord]=sort(x(:));ranks=zeros(size(s));ii=1;while ii<=numel(s),jj=ii;while jj<numel(s)&&s(jj+1)==s(ii),jj=jj+1;end;ranks(ii:jj)=.5*(ii+jj);ii=jj+1;end;tmp=zeros(size(ranks));tmp(ord)=ranks;ranks=tmp;
end

function [pci,sci]=local_corr_boot(x,y,B)
n=numel(x);pr=zeros(B,1);sr=zeros(B,1);for ii=1:B,j=randi(n,n,1);[pr(ii),~]=local_corr(x(j),y(j),false);[sr(ii),~]=local_corr(x(j),y(j),true);end;pci=[local_quantile(pr,.025) local_quantile(pr,.975)];sci=[local_quantile(sr,.025) local_quantile(sr,.975)];
end

function b=local_roughness_bins(rows)
x=[rows.profile_rms_slope].';e=[local_quantile(x,1/3) local_quantile(x,2/3)];labels={'low','middle','high'};b=repmat(struct('bin','','count',0,'slope_low',NaN,'slope_high',NaN,'mean_delta_tl_db',NaN,'mean_pe_power',NaN,'mean_bh_power',NaN,'phase_circular_mean_rad',NaN,'phase_resultant_length',NaN),3,1);
for k=1:3
 if k==1,m=x<=e(1);lo=-Inf;hi=e(1);elseif k==2,m=x>e(1)&x<=e(2);lo=e(1);hi=e(2);else,m=x>e(2);lo=e(2);hi=Inf;end
 q=rows(m);ph=[q.delta_phase_rad];z=mean(exp(1i*ph));b(k).bin=labels{k};b(k).count=numel(q);b(k).slope_low=lo;b(k).slope_high=hi;
 if ~isempty(q),b(k).mean_delta_tl_db=mean([q.delta_tl_db]);b(k).mean_pe_power=mean([q.G_pe_power]);b(k).mean_bh_power=mean([q.G_bellhop_power]);b(k).phase_circular_mean_rad=angle(z);b(k).phase_resultant_length=abs(z);end
end
end

function o=local_outlier_audit(rows,K)
n=numel(rows);ph=[rows.delta_phase_rad];mu=angle(mean(exp(1i*ph)));vals={abs([rows.delta_tl_db]),[rows.complex_relative_error],abs(angle(exp(1i*(ph-mu)))), [rows.min_mu],[rows.profile_max_curvature_per_m]};labels={'abs_delta_tl','complex_error','phase_distance','smallest_mu','largest_curvature'};o=repmat(local_outlier_row('',0,rows(1)),0,1);
for k=1:numel(labels),if k==4,[~,ix]=sort(vals{k},'ascend');else,[~,ix]=sort(vals{k},'descend');end;ix=ix(1:min(K,n));for j=1:numel(ix),o(end+1)=local_outlier_row(labels{k},j,rows(ix(j)));end,end
end

function r=local_outlier_row(metric,rank,q)
bad=~isfinite(q.complex_relative_error)||q.grazing_fraction>0||q.wall_residual_max_m>1e-9||q.phase_jump_error_rad>1e-10||q.q_reflect_error>1e-12||q.p_rotation_error>1e-12||q.q_rotation_error>1e-12||~q.all_post_range_positive;
if bad&&q.grazing_fraction>0,cat='C';reason='Bellhop grazing/beam issue';elseif bad&&q.pe_outer5_receiver_reflected_energy_fraction>1e-4,cat='D';reason='PE boundary/seam issue';elseif bad,cat='B';reason='numerical guard failure';else,cat='A';reason='legitimate physical realization';end
r=struct('metric',metric,'rank',rank,'seed',q.seed,'G_pe',q.G_pe,'G_bellhop',q.G_bellhop,'delta_tl_db',q.delta_tl_db,'delta_phase_rad',q.delta_phase_rad,'complex_relative_error',q.complex_relative_error,'profile_rms_slope',q.profile_rms_slope,'profile_rms_curvature_per_m',q.profile_rms_curvature_per_m,'profile_max_curvature_per_m',q.profile_max_curvature_per_m,'min_mu',q.min_mu,'grazing_fraction',q.grazing_fraction,'wall_residual_max_m',q.wall_residual_max_m,'phase_jump_error_rad',q.phase_jump_error_rad,'q_reflect_error',q.q_reflect_error,'p_rotation_error',q.p_rotation_error,'q_rotation_error',q.q_rotation_error,'classification',cat,'reason',reason);
end

function g=local_guard_summary(rows,cfg)
g=struct('all_finite',all(arrayfun(@(q)isfinite(real(q.G_pe))&&isfinite(imag(q.G_pe))&&isfinite(real(q.G_bellhop))&&isfinite(imag(q.G_bellhop)),rows)), ...
 'all_profile_provenance',all(~cellfun(@isempty,{rows.coeff_file_sha256})),'all_beams_hit',all([rows.wall_hit_count]==[rows.expected_beam_count]), ...
 'all_failed_rays_zero',all([rows.failed_ray_count]==0),'all_rejected_rays_zero',all([rows.rejected_ray_count]==0), ...
 'all_non_grazing',all([rows.grazing_fraction]==0)&&all([rows.min_mu]>=cfg.grazing_mu_threshold),'all_wall_residual',all([rows.wall_residual_max_m]<=cfg.wall_residual_limit_m), ...
 'all_phase',all([rows.phase_jump_error_rad]<=cfg.phase_jump_limit_rad),'all_beam_state',all([rows.q_reflect_error]<=cfg.state_limit)&&all([rows.p_rotation_error]<=cfg.state_limit)&&all([rows.q_rotation_error]<=cfg.state_limit), ...
 'all_positive_post_range',all([rows.all_post_range_positive]),'all_pe_edge',all([rows.pe_outer5_incident_energy_fraction]<=cfg.pe_edge5_limit)&&all([rows.pe_outer5_surface_reflected_energy_fraction]<=cfg.pe_edge5_limit)&&all([rows.pe_outer5_receiver_reflected_energy_fraction]<=cfg.pe_edge5_limit));g.all=all(structfun(@(v)logical(v),g));
end

function d=local_dimensionality_reference(path)
d=struct('file',path,'delta_tl_db',NaN,'delta_phase_rad',NaN,'complex_relative_error',NaN);if exist(path,'file')~=2,return;end;t=readtable(path);d.delta_tl_db=t.delta_tl_db(1);d.delta_phase_rad=t.delta_phase_rad(1);d.complex_relative_error=t.complex_relative_error(1);
end

function c=local_classification(guards,gates)
if ~guards.all,c='NUMERICAL_OR_APPLICABILITY_LIMIT';elseif gates.all,c='STATISTICALLY_ESTABLISHED_MODEL_DISCREPANCY';else,c='PRELIMINARY_MODEL_DISCREPANCY';end
end

function local_write_outputs(study)
out=study.config.output_dir;
writetable(struct2table(study.rows),fullfile(out,'per_seed_results.csv'));
writetable(struct2table(study.rows),fullfile(out,'surface_geometry_statistics.csv'));
writetable(struct2table(study.running),fullfile(out,'convergence_by_sample_count.csv'));
writetable(local_final_statistics_table(study.running(end)),fullfile(out,'ensemble_statistics.csv'));
writetable(local_boot_table(study.bootstrap),fullfile(out,'bootstrap_confidence_intervals.csv'));
writetable(struct2table(study.correlations),fullfile(out,'discrepancy_correlations.csv'));
writetable(struct2table(study.outliers),fullfile(out,'outlier_audit.csv'));
writetable(struct2table(study.roughness_bins),fullfile(out,'roughness_bins.csv'));
end

function t=local_final_statistics_table(r)
t=table(r.sample_count,r.delta_tl_mean_db,r.delta_tl_std_db,r.delta_tl_median_db,r.delta_tl_p05_db,r.delta_tl_p95_db,r.pe_power_mean,r.bh_power_mean,r.power_difference_mean,r.phase_circular_mean_rad,r.phase_circular_std_rad,r.phase_resultant_length,'VariableNames',{'sample_count','delta_tl_mean_db','delta_tl_std_db','delta_tl_median_db','delta_tl_p05_db','delta_tl_p95_db','pe_power_mean','bh_power_mean','power_difference_mean','phase_circular_mean_rad','phase_circular_std_rad','phase_resultant_length'});
end

function t=local_boot_table(b)
names=fieldnames(b);names=setdiff(names,{'repetitions','seed'});
t=table('Size',[numel(names) 5],'VariableTypes',{'string','double','double','double','double'},'VariableNames',{'metric','estimate','ci95_low','ci95_high','half_width'});
for ii=1:numel(names),q=b.(names{ii});t.metric(ii)=string(names{ii});t.estimate(ii)=q.estimate;t.ci95_low(ii)=q.ci95_low;t.ci95_high(ii)=q.ci95_high;t.half_width(ii)=q.half_width;end
end

function local_make_figures(study)
out=fullfile(study.config.output_dir,'figures');if ~exist(out,'dir'),mkdir(out);end;n=[study.running.sample_count];
local_line(n,[study.running.delta_tl_mean_db],'sample count','running mean delta TL (dB)',fullfile(out,'running_mean_delta_tl.png'));
local_line(n,[study.running.delta_tl_std_db],'sample count','running std delta TL (dB)',fullfile(out,'running_std_delta_tl.png'));
local_ci(study.bootstrap.mean_delta_tl_db,fullfile(out,'bootstrap_ci_mean_delta_tl.png'));
local_hist([study.rows.G_pe_power],[study.rows.G_bellhop_power],'|G|^2','reflected power',fullfile(out,'power_distribution.png'));
local_hist([study.rows.delta_tl_db],[],'delta TL (dB)','count',fullfile(out,'delta_tl_histogram.png'));
local_phase([study.rows.delta_phase_rad],fullfile(out,'delta_phase_circular_histogram.png'));
local_scatter([study.rows.profile_rms_slope],[study.rows.delta_tl_db],'RMS slope','delta TL (dB)',fullfile(out,'delta_tl_vs_rms_slope.png'));
local_scatter([study.rows.profile_rms_curvature_per_m],[study.rows.delta_tl_db],'RMS curvature (1/m)','delta TL (dB)',fullfile(out,'delta_tl_vs_rms_curvature.png'));
local_scatter([study.rows.profile_rms_curvature_per_m],[study.rows.complex_relative_error],'RMS curvature (1/m)','complex error',fullfile(out,'complex_error_vs_rms_curvature.png'));
local_outlier_profiles(study,fullfile(out,'selected_outlier_surface_profiles.png'));
end

function local_line(x,y,xlab,ylab,path)
f=figure('Visible','off');plot(x,y,'o-','LineWidth',1.2);grid on;xlabel(xlab);ylabel(ylab);print(f,path,'-dpng','-r120');close(f);
end
function local_ci(q,path)
f=figure('Visible','off');errorbar(1,q.estimate,q.half_width,'o','LineWidth',1.2);grid on;xlim([0 2]);set(gca,'XTick',1,'XTickLabel',{'mean delta TL'});ylabel('95% CI half-width (dB)');print(f,path,'-dpng','-r120');close(f);
end
function local_hist(a,b,xlab,ylab,path)
f=figure('Visible','off');if isempty(b),histogram(a,10);else,hold on;histogram(a,10,'DisplayName','PE');histogram(b,10,'DisplayName','Bellhop');legend('Location','best');end;grid on;xlabel(xlab);ylabel(ylab);print(f,path,'-dpng','-r120');close(f);
end
function local_phase(p,path)
f=figure('Visible','off');polarhistogram(p,12);title('delta phase circular histogram');print(f,path,'-dpng','-r120');close(f);
end
function local_scatter(x,y,xlab,ylab,path)
f=figure('Visible','off');scatter(x,y,18,'filled');grid on;xlabel(xlab);ylabel(ylab);print(f,path,'-dpng','-r120');close(f);
end
function local_outlier_profiles(study,path)
ids=unique([study.outliers.seed]);ids=ids(1:min(3,numel(ids)));f=figure('Visible','off');hold on;
for ii=1:numel(ids)
 file=fullfile(study.config.output_dir,sprintf('seed_%d_coefficients.csv',ids(ii)));if exist(file,'file')~=2,continue;end
 p=load_fixed_pm_profile_for_pe_bellhop_validation(file);s=linspace(-p.span_m/2,p.span_m/2,4097);e=evaluate_fixed_pm_fourier_profile(p,s);plot(s,e.eta_m,'DisplayName',sprintf('seed %d',ids(ii)));
end
grid on;xlabel('profile coordinate (m)');ylabel('eta (m)');legend('Location','best');print(f,path,'-dpng','-r120');close(f);
end

function local_write_report(study,path)
fid=fopen(path,'w','n','UTF-8');if fid<0,error('Cannot write %s.',path);end;cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
c=study.config;r=study.running(end);fprintf(fid,'# PE--Bellhop PM model-discrepancy statistical study\n\n状态：**%s**\n\n',study.classification);
fprintf(fid,'冻结条件：4 kHz、uniform c=1500 m/s、z_tx=100 m、z_rx=3 m、sigma=0.3 m、PM span=160 m、requested Kmax=0.5 rad/m、realized Kmax≈0.471238898 rad/m；PE W=192.1875 m/nx=984/step=0.05 m/sponge off；Bellhop source=%s、profile N=4097、5001 beams、step=0.05 m、sector ±15°。\n\n',c.common_overrides.source_geometry);
fprintf(fid,'The primary ensemble is Stage 3B independent Gaussian cosine/sine coefficient amplitude sampling with variance S(k)*Delta-k, shared identically by PE and Bellhop. No Hs normalization, recentering, smoothing, taper, source fit or Kmax tuning is used. Seed-level bootstrap repetitions: %d. Stage 3A phase-only results remain a separate diagnostic and are not pooled with this ensemble.\n\n',study.bootstrap.repetitions);
fprintf(fid,'## Ensemble result\n\n| quantity | value |\n|---|---:|\n| final seed count | %d |\n| mean delta TL (dB) | %.8g |\n| std delta TL (dB) | %.8g |\n| median delta TL (dB) | %.8g |\n| p5 / p95 delta TL (dB) | %.8g / %.8g |\n| mean PE reflected power | %.8g |\n| mean Bellhop reflected power | %.8g |\n| mean power difference | %.8g |\n| circular mean delta phase (rad) | %.8g |\n| circular std (rad) | %.8g |\n| mean resultant length | %.8g |\n',r.sample_count,r.delta_tl_mean_db,r.delta_tl_std_db,r.delta_tl_median_db,r.delta_tl_p05_db,r.delta_tl_p95_db,r.pe_power_mean,r.bh_power_mean,r.power_difference_mean,r.phase_circular_mean_rad,r.phase_circular_std_rad,r.phase_resultant_length);
fprintf(fid,'\n## Running-prefix convergence\n\n| M | mean delta TL | std | median | p5 | p95 | PE power | BH power | circular mean phase | circular std | R |\n|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n');for ii=1:numel(study.running),q=study.running(ii);fprintf(fid,'| %d | %.8g | %.8g | %.8g | %.8g | %.8g | %.8g | %.8g | %.8g | %.8g | %.8g |\n',q.sample_count,q.delta_tl_mean_db,q.delta_tl_std_db,q.delta_tl_median_db,q.delta_tl_p05_db,q.delta_tl_p95_db,q.pe_power_mean,q.bh_power_mean,q.phase_circular_mean_rad,q.phase_circular_std_rad,q.phase_resultant_length);end
g=study.convergence_gates;fprintf(fid,'\nGates 24→32: mean delta-TL change %.8g dB (%s), PE/Bellhop power relative changes %.8g / %.8g (%s), circular-mean phase change %.8g rad (%s), mean delta-TL bootstrap half-width %.8g dB (%s). These are engineering gates, not physical laws.\n\n',g.mean_delta_tl_change_24_to_32_db,local_pass(g.mean_delta_tl_pass),g.pe_power_relative_change_24_to_32,g.bh_power_relative_change_24_to_32,local_pass(g.power_pass),g.phase_change_24_to_32_rad,local_pass(g.phase_pass),g.bootstrap_delta_tl_half_width_db,local_pass(g.bootstrap_pass));
fprintf(fid,'## Bootstrap confidence intervals\n\n| metric | estimate | low | high | half-width |\n|---|---:|---:|---:|---:|\n');bt=local_boot_table(study.bootstrap);for ii=1:height(bt),fprintf(fid,'| %s | %.8g | %.8g | %.8g | %.8g |\n',char(bt.metric(ii)),bt.estimate(ii),bt.ci95_low(ii),bt.ci95_high(ii),bt.half_width(ii));end
fprintf(fid,'\n## Numerical and applicability guards\n\n');n=fieldnames(study.guards);for ii=1:numel(n),fprintf(fid,'- %s: %s\n',n{ii},local_pass(study.guards.(n{ii})));end
fprintf(fid,'\n## Dimensionality reference\n\nThe existing 1T→2T sensitivity is %.8g dB and %.8g rad from `%s`; it is reported alongside, not subtracted from, the model discrepancy.\n\n',study.dimensionality_reference.delta_tl_db,study.dimensionality_reference.delta_phase_rad,study.dimensionality_reference.file);
fprintf(fid,'## Geometry association\n\nPearson and Spearman statistics are exploratory associations with seed-level bootstrap CIs. Wrapped phase is not correlated directly; phase_cos and phase_sin rows provide circular-linear diagnostics.\n\n');
fprintf(fid,'## Outlier audit\n\nTop %d cases per metric are retained; no seed is deleted. Classes are A legitimate realization, B numerical contamination, C Bellhop grazing/beam issue, and D PE boundary/seam issue.\n\n',c.outlier_count);fprintf(fid,'| metric | rank | seed | delta TL | phase | complex error | RMS slope | RMS curvature | min |u.n| | grazing | wall residual | class |\n|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n');for ii=1:numel(study.outliers),o=study.outliers(ii);fprintf(fid,'| %s | %d | %d | %.8g | %.8g | %.8g | %.8g | %.8g | %.8g | %.8g | %.8g | %s |\n',o.metric,o.rank,o.seed,o.delta_tl_db,o.delta_phase_rad,o.complex_relative_error,o.profile_rms_slope,o.profile_rms_curvature_per_m,o.min_mu,o.grazing_fraction,o.wall_residual_max_m,o.classification);end
fprintf(fid,'\n## Classification\n\nThis study asks whether the discrepancy between two independently validated approximate models is statistically stable; it does not claim that either model is physically correct or force their agreement. Final classification: **%s**.\n\nArtifacts are under `results/validation/pe_bellhop_pm_model_discrepancy_statistics/`, including `per_seed_results.csv`, `ensemble_statistics.csv`, `convergence_by_sample_count.csv`, `bootstrap_confidence_intervals.csv`, `surface_geometry_statistics.csv`, `discrepancy_correlations.csv`, `outlier_audit.csv`, `roughness_bins.csv`, `result.mat`, and `figures/`.\n',study.classification);
end

function s=local_pass(x),if x,s='PASS';else,s='FAIL';end,end
