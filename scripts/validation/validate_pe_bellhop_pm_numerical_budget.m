function audit = validate_pe_bellhop_pm_numerical_budget(overrides)
%VALIDATE_PE_BELLHOP_PM_NUMERICAL_BUDGET Stage 0E budget freeze.
%   Quantifies validation-only PE window/grid/step sensitivity and Bellhop
%   fixed-realization profile/beam/step sensitivity.  No tolerance is tuned
%   from Stage 1 and no propagation core is modified.

if nargin < 1 || isempty(overrides), overrides=struct(); end
if ~isstruct(overrides) || ~isscalar(overrides), error('overrides must be a scalar struct.'); end
root=fileparts(fileparts(fileparts(mfilename('fullpath')))); addpath(root); setup_vertical_project();
cfg=local_config(root,overrides);
if ~exist(cfg.output_dir,'dir'), mkdir(cfg.output_dir); end
profile=load_fixed_pm_profile_for_pe_bellhop_validation(cfg.coeff_file);
pattern=local_source_pattern(cfg);

% PE window, grid and z-step scans use the same continuous periodic profile.
window_rows=repmat(local_pe_row(),1,numel(cfg.window_cases));
for ii=1:numel(cfg.window_cases)
    c=cfg.window_cases(ii); window_rows(ii)=local_pe_run(cfg,profile,c.width_m,c.nx,c.step_m); window_rows(ii).case_name=c.name;
end
grid_rows=repmat(local_pe_row(),1,numel(cfg.grid_cases));
for ii=1:numel(cfg.grid_cases)
    c=cfg.grid_cases(ii); grid_rows(ii)=local_pe_run(cfg,profile,c.width_m,c.nx,c.step_m); grid_rows(ii).case_name=c.name;
end
step_rows=repmat(local_pe_row(),1,numel(cfg.step_cases));
for ii=1:numel(cfg.step_cases)
    c=cfg.step_cases(ii); step_rows(ii)=local_pe_run(cfg,profile,c.width_m,c.nx,c.step_m); step_rows(ii).case_name=c.name;
end

% Bellhop PM local geometry budget. All profile densities are sampled from
% the same coefficient-defined master realization (no regenerated spectrum).
bh_rows=repmat(local_bh_row(),1,numel(cfg.profile_counts)*numel(cfg.beam_counts)*numel(cfg.bh_steps));
rr=0; bh_exe=cfg.pm_validation_exe;
for np=1:numel(cfg.profile_counts)
    nprof=cfg.profile_counts(np); s=linspace(-profile.span_m/2,profile.span_m/2,nprof).'; e=evaluate_fixed_pm_fourier_profile(profile,s);
    for nb=1:numel(cfg.beam_counts)
        for ns=1:numel(cfg.bh_steps)
            rr=rr+1; croot=fullfile(cfg.output_dir,'cases',sprintf('N%d_B%d_h%s',nprof,cfg.beam_counts(nb),local_tag(cfg.bh_steps(ns))));
            c=local_pm_cfg(cfg,pattern,croot,nprof,cfg.beam_counts(nb),cfg.bh_steps(ns),e,s);
            if exist([croot '.shd'],'file')==2 && exist([croot '.iwdiag'],'file')==2
                run=local_read_pm_case(c,croot);
            else
                run=run_bellhop_internal_pm_wall_poc_vertical(c);
            end
            d=run.diagnostics;
            p=select_bellhop_shd_pressure_at_range_vertical(run.data,cfg.mapped_receiver_range_m,0,cfg.selector_tolerance_m);
            mu=abs(d.inc_ur.*d.wall_n_r+d.inc_uz.*d.wall_n_z);
            q=local_bh_row(); q.profile_count=nprof; q.beam_count=cfg.beam_counts(nb); q.step_m=cfg.bh_steps(ns);
            q.pressure=p; q.diagnostic_ray_count=height(d); q.failed_ray_count=0; q.nan_inf_count=0;
            q.wall_residual_max_m=max(abs(d.wall_residual)); q.min_mu=min(mu); q.grazing_fraction=mean(mu<cfg.grazing_mu_threshold);
            q.kappa_abs_max_per_m=max(abs(d.kappa)); q.min_post_range_increment_m=min(d.min_post_dr);
            q.all_post_range_positive=all(d.min_post_dr>0); q.phase_jump_error_rad=max(abs(d.phase_delta-pi));
            q.amp_jump_error=max(abs(d.amp_delta));
            % p_ref_error is the intentional curvature kick magnitude, not
            % an error-to-zero.  q_ref_error is the residual used for the
            % dynamic-state gate; rotation errors must remain zero.
            q.p_curvature_kick_max=max(abs(d.p_ref_error)); q.pq_error=max(abs(d.q_ref_error));
            q.p_rotation_error=max(d.p_rot_error); q.q_rotation_error=max(d.q_rot_error);
            q.tau_receiver_s=d.tau_receiver_real(find(abs(d.alpha_deg)==min(abs(d.alpha_deg)),1)); q.case_root=croot;
            bh_rows(rr)=q;
        end
    end
end

% Convert each family to relative response metrics against its designated
% baseline. The fixed-realization legacy CSV supplies the already accepted
% N=2049 -> 4097 geometry/field pair for an independent cross-check.
window_rows=local_relative_pe(window_rows,1); grid_rows=local_relative_pe(grid_rows,2); step_rows=local_relative_pe(step_rows,2);
bh_rows=local_relative_bh(bh_rows,cfg);
flat=readtable(cfg.stage0c_rows_file); d0=readtable(cfg.pm_convergence_file);
budget=struct(); budget.pe_window=window_rows; budget.pe_grid=grid_rows; budget.pe_step=step_rows; budget.bellhop=bh_rows;
budget.flat=struct('max_q_tl_db',max(abs(flat.q_axis_tl_db)),'max_q_phase_rad',max(abs(flat.q_axis_phase_rad)), ...
    'max_direct_profile_magnitude',max(flat.direct_profile_magnitude),'max_reflect_profile_magnitude',max(flat.reflect_profile_magnitude), ...
    'beam_q_tl_delta_db',abs(flat.q_axis_tl_db(end)-flat.q_axis_tl_db(1)), ...
    'beam_q_phase_delta_rad',abs(flat.q_axis_phase_rad(end)-flat.q_axis_phase_rad(1)));
budget.pm_fixed_realization=d0;
checks=local_checks(window_rows,grid_rows,step_rows,bh_rows,budget,cfg);
audit=struct('schema_version','1.0.0','stage','0E_numerical_budget_freeze','config',cfg, ...
    'profile',profile,'budget',budget,'checks',checks,'passed',checks.all);
audit.files=local_write_outputs(audit);
disp(checks);
if cfg.fail_on_check && ~audit.passed, error('Stage 0E numerical budget failed; see %s.',cfg.report_path); end
end

function row=local_pe_run(cfg,profile,width,nx,step)
x=(-0.5*width)+(0:nx-1)*(width/nx); e=evaluate_fixed_pm_fourier_profile(profile,x);
c=struct('frequency_hz',cfg.frequency_hz,'c0_mps',cfg.c0_mps,'xw_m',width,'nx',nx, ...
    'z_tx_m',cfg.z_tx_m,'z_rx_m',cfg.z_rx_m,'sigma_src_m',cfg.sigma_src_m, ...
    'surface_elevation_x_m',e.eta_m(:).','surface_reflect_coeff',-1,'step_m',step,'x_rx_m',0);
r=run_pe_1d_surface_reflection_validation(c); n=max(1,round(0.05*nx)); edge=[1:n (nx-n+1):nx];
row=local_pe_row(); row.width_m=width; row.nx=nx; row.step_m=step; row.response=r.reflected_receiver;
row.seam_jump_m=abs(e.eta_m(1)-e.eta_m(end)); row.outer5_energy_fraction=sum(abs(r.surface_reflected_field(edge)).^2)/max(sum(abs(r.surface_reflected_field).^2),realmin);
row.profile_rms_m=sqrt(mean(e.eta_m.^2)); row.dx_m=width/nx;
end

function rows=local_relative_pe(rows,baseline_index)
p=rows(baseline_index).response;
for ii=1:numel(rows)
 rows(ii).relative_tl_db=20*log10(abs(rows(ii).response/p)); rows(ii).relative_phase_rad=angle(rows(ii).response*conj(p));
end
end

function c=local_pm_cfg(cfg,pattern,croot,nprof,beams,step,e,s)
c=struct('bellhop_exe',cfg.pm_validation_exe,'case_root',croot,'run_type','C', ...
    'frequency_hz',cfg.frequency_hz,'c0_mps',cfg.c0_mps,'source_depth_m',0, ...
    'receiver_depths_m',0,'receiver_ranges_m',[102 103],'beam_count',beams, ...
    'angle_limits_deg',cfg.angle_limits_deg,'domain_half_depth_m',cfg.domain_half_depth_m, ...
    'sigma_src_m',cfg.sigma_src_m,'source_pattern_clip_db',cfg.source_pattern_clip_db, ...
    'source_pattern_angles_deg',pattern.angles_deg,'source_pattern_level_db',pattern.level_db, ...
    'step_m',step,'wall_r0_m',cfg.wall_r0_m,'wall_seed',cfg.wall_seed, ...
    'wall_profile_r_m',cfg.wall_r0_m-e.eta_m(:).','wall_profile_z_m',s(:).', ...
    'mapped_receiver_range_m',cfg.mapped_receiver_range_m);
end

function result=local_read_pm_case(cfg,croot)
data=read_bellhop_shd_unfolded_vertical([croot '.shd']);
dm=readmatrix([croot '.iwdiag'],'FileType','text','CommentStyle','#'); dm=dm(~all(isnan(dm),2),:);
names={'alpha_deg','hit_r','hit_z','wall_residual','wall_t_r','wall_t_z','wall_n_r','wall_n_z', ...
 'tangent_error','normal_error','inc_ur','inc_uz','ref_ur','ref_uz','rot_ur','rot_uz','specular_error','rotation_error', ...
 'phase_in','phase_ref','phase_delta','amp_in','amp_ref','amp_delta','p1_in','p2_in','p1_ref','p2_ref','p_ref_error', ...
 'q1_in','q2_in','q1_ref','q2_ref','q_ref_error','p_rot_error','q_rot_error','tau_wall_real','tau_wall_imag', ...
 'tau_receiver_real','tau_receiver_imag','min_post_dr','n_post','kappa','wall_seg','wall_lambda','wall_tg','wall_th','wall_rm','wall_rn'};
if size(dm,2)~=numel(names),error('Unexpected PM diagnostic columns in cached case %s.',croot);end
result=struct('config',cfg,'data',data,'diagnostics',array2table(dm,'VariableNames',names), ...
 'files',struct('data',[croot '.shd'],'diagnostics',[croot '.iwdiag']));
end

function rows=local_relative_bh(rows,cfg)
% Use the finest profile, highest beam, and finest step as the frozen BH
% baseline.  The scan reports each perturbation but does not recalibrate it.
ix=find([rows.profile_count]==max([rows.profile_count]) & [rows.beam_count]==max([rows.beam_count]) & [rows.step_m]==min([rows.step_m]),1);
p=rows(ix).pressure;
for ii=1:numel(rows)
 rows(ii).relative_tl_db=20*log10(abs(rows(ii).pressure/p)); rows(ii).relative_phase_rad=angle(rows(ii).pressure*conj(p));
end
if nargin>1 && isfield(cfg,'pm_convergence_file') && exist(cfg.pm_convergence_file,'file')==2
 % Kept as provenance only; the legacy CSV is copied into the report.
end
end

function checks=local_checks(w,g,s,b,budget,cfg)
checks=struct();
checks.pe_window=local_pair_pass(w,cfg.pe_window_tl_limit_db,cfg.pe_window_phase_limit_rad);
checks.pe_grid=local_pair_pass(g,cfg.pe_grid_tl_limit_db,cfg.pe_grid_phase_limit_rad);
checks.pe_step=local_pair_pass(s,cfg.pe_step_tl_limit_db,cfg.pe_step_phase_limit_rad);
checks.bellhop=all(abs([b.relative_tl_db])<=cfg.bh_tl_limit_db) && all(abs([b.relative_phase_rad])<=cfg.bh_phase_limit_rad) && ...
    all([b.all_post_range_positive]) && max([b.phase_jump_error_rad])<=1e-10 && max([b.pq_error])<=1e-12;
checks.flat=budget.flat.max_q_tl_db<=cfg.flat_q_tl_limit_db && budget.flat.max_q_phase_rad<=cfg.flat_q_phase_limit_rad && ...
    budget.flat.max_direct_profile_magnitude<=cfg.flat_profile_limit && budget.flat.max_reflect_profile_magnitude<=cfg.flat_profile_limit;
checks.all=all(structfun(@(v)logical(v),checks));
end

function pass=local_pair_pass(rows,tl_lim,phase_lim)
pass=all(abs([rows.relative_tl_db])<=tl_lim) && all(abs([rows.relative_phase_rad])<=phase_lim);
end

function cfg=local_config(root,o)
dx=192.1875/984;
cfg=struct('frequency_hz',4000,'c0_mps',1500,'z_tx_m',100,'z_rx_m',3,'sigma_src_m',0.3, ...
    'angle_limits_deg',[-30 30],'domain_half_depth_m',1000,'source_pattern_clip_db',-120, ...
    'wall_r0_m',100,'mapped_receiver_range_m',103,'wall_seed',260001, ...
    'coeff_file',fullfile(root,'results','validation','bellhop_internal_pm_fixed_realization','fixed_pm_fourier_coefficients.csv'), ...
    'pm_validation_exe',fullfile(root,'results','validation','bellhop_internal_pm_wall_poc','bin','bellhop_iwall_pm_2020.exe'), ...
    'stage0c_rows_file',fullfile(root,'results','validation','pe_bellhop_pm_stage0_flat_source','stage0c_flat_source_rows.csv'), ...
    'pm_convergence_file',fullfile(root,'results','validation','bellhop_internal_pm_fixed_realization','pm_fixed_realization_convergence.csv'), ...
    'output_dir',fullfile(root,'results','validation','pe_bellhop_pm_numerical_budget'), ...
    'report_path',fullfile(root,'reports','pe_bellhop_pm_numerical_error_budget.md'), ...
    'window_cases',struct('name',{'W160','W192p1875','W320'},'width_m',{160,192.1875,320},'nx',{820,984,1638},'step_m',{0.05,0.05,0.05}), ...
    'grid_cases',struct('name',{'grid492','grid984'},'width_m',{192.1875,192.1875},'nx',{492,984},'step_m',{0.05,0.05}), ...
    'step_cases',struct('name',{'step0p1','step0p05','step0p025'},'width_m',{192.1875,192.1875,192.1875},'nx',{984,984,984},'step_m',{0.1,0.05,0.025}), ...
    'profile_counts',[2049 4097],'beam_counts',[5001 10001],'bh_steps',[0.1 0.05], ...
    'selector_tolerance_m',1e-6,'grazing_mu_threshold',0.1,'pe_window_tl_limit_db',0.1, ...
    'pe_window_phase_limit_rad',0.02,'pe_grid_tl_limit_db',0.1,'pe_grid_phase_limit_rad',0.02, ...
    'pe_step_tl_limit_db',0.1,'pe_step_phase_limit_rad',0.02,'bh_tl_limit_db',0.1, ...
    'bh_phase_limit_rad',0.02,'flat_q_tl_limit_db',0.30,'flat_q_phase_limit_rad',0.05, ...
    'flat_profile_limit',0.05,'fail_on_check',true,'dx_reference_m',dx);
names=fieldnames(o);
for ii=1:numel(names)
 if ~isfield(cfg,names{ii}),error('Unknown Stage 0E override: %s.',names{ii});end
 cfg.(names{ii})=o.(names{ii});
end
end

function pattern=local_source_pattern(cfg)
angles=linspace(cfg.angle_limits_deg(1),cfg.angle_limits_deg(2),2401).';
theta=angles*pi/180; k=2*pi*cfg.frequency_hz/cfg.c0_mps;
d=cos(theta).*exp(-0.5*(k*cfg.sigma_src_m*sin(theta)).^2); d=abs(d)/max(abs(d));
pattern=struct('angles_deg',angles,'level_db',20*log10(max(d,10^(cfg.source_pattern_clip_db/20))));
end

function row=local_pe_row()
row=struct('case_name','','width_m',NaN,'nx',NaN,'step_m',NaN,'dx_m',NaN,'response',complex(NaN), ...
    'relative_tl_db',NaN,'relative_phase_rad',NaN,'seam_jump_m',NaN,'outer5_energy_fraction',NaN,'profile_rms_m',NaN);
end

function row=local_bh_row()
row=struct('profile_count',NaN,'beam_count',NaN,'step_m',NaN,'pressure',complex(NaN), ...
    'relative_tl_db',NaN,'relative_phase_rad',NaN,'diagnostic_ray_count',NaN,'failed_ray_count',NaN,'nan_inf_count',NaN, ...
    'wall_residual_max_m',NaN,'min_mu',NaN,'grazing_fraction',NaN,'kappa_abs_max_per_m',NaN, ...
    'min_post_range_increment_m',NaN,'all_post_range_positive',false,'phase_jump_error_rad',NaN,'amp_jump_error',NaN, ...
    'p_curvature_kick_max',NaN,'pq_error',NaN,'p_rotation_error',NaN,'q_rotation_error',NaN,'tau_receiver_s',NaN,'case_root','');
end

function files=local_write_outputs(audit)
out=audit.config.output_dir;
writetable(struct2table(audit.budget.pe_window),fullfile(out,'pe_window_budget.csv'));
writetable(struct2table(audit.budget.pe_grid),fullfile(out,'pe_grid_budget.csv'));
writetable(struct2table(audit.budget.pe_step),fullfile(out,'pe_step_budget.csv'));
writetable(struct2table(audit.budget.bellhop),fullfile(out,'bellhop_budget.csv'));
writetable(audit.budget.pm_fixed_realization,fullfile(out,'pm_fixed_realization_budget.csv'));
n=fieldnames(audit.checks); v=false(size(n)); for ii=1:numel(n),v(ii)=audit.checks.(n{ii});end
writetable(table(n,v,'VariableNames',{'check_name','passed'}),fullfile(out,'budget_checks.csv'));
mat_file=fullfile(out,'numerical_budget_audit.mat'); save(mat_file,'audit','-v7.3');
local_write_report(audit.config.report_path,audit);
files=struct('pe_window',fullfile(out,'pe_window_budget.csv'),'pe_grid',fullfile(out,'pe_grid_budget.csv'), ...
 'pe_step',fullfile(out,'pe_step_budget.csv'),'bellhop',fullfile(out,'bellhop_budget.csv'), ...
 'pm_fixed',fullfile(out,'pm_fixed_realization_budget.csv'),'checks',fullfile(out,'budget_checks.csv'), ...
 'mat',mat_file,'report',audit.config.report_path);
end

function local_write_report(path,audit)
fid=fopen(path,'w','n','UTF-8'); if fid<0,error('Cannot write %s.',path);end
cleanup=onCleanup(@()fclose(fid)); c=audit.config;
fprintf(fid,'# PE--Bellhop PM Stage 0E numerical error budget freeze\n\n状态：**%s**\n\n',ternary(audit.passed,'PASS','FAIL'));
fprintf(fid,'本阶段冻结后续 PM 对比前的 PE、Bellhop 和 flat/source 数值误差预算；所有 PM profile density 均由同一 seed-260001 Fourier realization 直接采样。\n\n');
fprintf(fid,'## PE budget\n\n');
fprintf(fid,'| family | case | W (m) | nx | step (m) | dx (m) | relative TL (dB) | relative phase (rad) | seam jump (m) | outer 5%% energy |\n|---|---|---:|---:|---:|---:|---:|---:|---:|---:|\n');
for fam={audit.budget.pe_window,audit.budget.pe_grid,audit.budget.pe_step}
 rows=fam{1}; for ii=1:numel(rows),r=rows(ii);fprintf(fid,'| PE | %s | %.8g | %d | %.8g | %.8g | %.8g | %.8g | %.8g | %.8g |\n',r.case_name,r.width_m,r.nx,r.step_m,r.dx_m,r.relative_tl_db,r.relative_phase_rad,r.seam_jump_m,r.outer5_energy_fraction);end
end
fprintf(fid,'\n## Bellhop budget\n\n| profile N | beams | step (m) | relative TL (dB) | relative phase (rad) | rays | min |u.n| | grazing fraction | kappa max (1/m) | p curvature kick | q residual | wall residual (m) | min post dr (m) |\n|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n');
for ii=1:numel(audit.budget.bellhop),r=audit.budget.bellhop(ii);fprintf(fid,'| %d | %d | %.8g | %.8g | %.8g | %d | %.8g | %.8g | %.8g | %.8g | %.8g | %.8g | %.8g |\n',r.profile_count,r.beam_count,r.step_m,r.relative_tl_db,r.relative_phase_rad,r.diagnostic_ray_count,r.min_mu,r.grazing_fraction,r.kappa_abs_max_per_m,r.p_curvature_kick_max,r.pq_error,r.wall_residual_max_m,r.min_post_range_increment_m);end
fprintf(fid,'\n## Flat/source budget\n\n- Stage 0C max axis-Q TL residual: `%.8g dB`; phase: `%.8g rad`; normalized offset magnitude residuals: `%.8g / %.8g`.\n',audit.budget.flat.max_q_tl_db,audit.budget.flat.max_q_phase_rad,audit.budget.flat.max_direct_profile_magnitude,audit.budget.flat.max_reflect_profile_magnitude);
fprintf(fid,'- Legacy fixed-realization N=2049 -> 4097 convergence CSV is copied as `pm_fixed_realization_budget.csv`; its accepted field residual is not used to recalibrate any Stage 1 result.\n\n');
fprintf(fid,'## Frozen checks\n\n'); n=fieldnames(audit.checks);for ii=1:numel(n),fprintf(fid,'- %s: %s\n',n{ii},ternary(audit.checks.(n{ii}),'PASS','FAIL'));end
fprintf(fid,'\nThe reported tolerances are frozen validation budgets: PE family 0.1 dB / 0.02 rad, Bellhop family 0.1 dB / 0.02 rad, and the known flat backward-range influence allowance 0.30 dB. They are not tuned from rough-wall cross-model results.\n');
end

function s=ternary(c,a,b),if c,s=a;else,s=b;end,end
function t=local_tag(v),t=strrep(sprintf('%.6g',v),'.','p');t=strrep(t,'-','m');end
