function comparison = validate_pe_bellhop_pm_stage1_tier1(overrides)
%VALIDATE_PE_BELLHOP_PM_STAGE1_TIER1 Fixed-PM reflected-ratio comparison.
if nargin<1 || isempty(overrides),overrides=struct();end
if ~isstruct(overrides)||~isscalar(overrides),error('overrides must be a scalar struct.');end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));addpath(root);setup_vertical_project();
cfg=local_config(root,overrides);if ~exist(cfg.output_dir,'dir'),mkdir(cfg.output_dir);end
profile=load_fixed_pm_profile_for_pe_bellhop_validation(cfg.coeff_file);pattern=local_source_pattern(cfg);
x=(-0.5*cfg.xw_m)+(0:cfg.nx-1)*(cfg.xw_m/cfg.nx);
base=struct('frequency_hz',cfg.frequency_hz,'c0_mps',cfg.c0_mps,'xw_m',cfg.xw_m,'nx',cfg.nx,'z_tx_m',cfg.z_tx_m,'z_rx_m',cfg.z_rx_m,'sigma_src_m',cfg.sigma_src_m,'surface_reflect_coeff',-1,'step_m',cfg.step_m,'x_rx_m',0);
flat_pe=run_pe_1d_surface_reflection_validation(setfield(base,'surface_elevation_x_m',zeros(1,cfg.nx))); %#ok<SFLD>
e=evaluate_fixed_pm_fourier_profile(profile,x);
rough_pe=run_pe_1d_surface_reflection_validation(setfield(base,'surface_elevation_x_m',e.eta_m(:).')); %#ok<SFLD>
flat_bh=local_flat_bellhop(cfg,pattern);rough_bh=local_rough_bellhop(cfg,pattern,profile);
flat_p=select_bellhop_shd_pressure_at_range_vertical(flat_bh.data,cfg.mapped_receiver_range_m,0,cfg.selector_tolerance_m);
rough_p=select_bellhop_shd_pressure_at_range_vertical(rough_bh.data,cfg.mapped_receiver_range_m,0,cfg.selector_tolerance_m);
Gpe=rough_pe.reflected_receiver/flat_pe.reflected_receiver;Gbh=rough_p/flat_p;delta=Gpe/Gbh;
d=rough_bh.diagnostics;mu=abs(d.inc_ur.*d.wall_n_r+d.inc_uz.*d.wall_n_z);
comparison=struct('schema_version','1.0.0','stage','1A_fixed_pm_tier1','config',cfg,'profile',profile,'flat_pe',flat_pe,'rough_pe',rough_pe,'flat_bellhop',flat_bh,'rough_bellhop',rough_bh,'G_pe',Gpe,'G_bellhop',Gbh,'delta_G',delta,'delta_tl_db',20*log10(abs(delta)),'delta_phase_rad',angle(delta),'delta_complex_relative_error',abs(Gpe-Gbh)/max(abs(Gbh),realmin),'geometry',struct('wall_residual_max_m',max(abs(d.wall_residual)),'min_mu',min(mu),'grazing_fraction',mean(mu<cfg.grazing_mu_threshold),'kappa_abs_max_per_m',max(abs(d.kappa)),'min_post_range_increment_m',min(d.min_post_dr),'all_post_range_positive',all(d.min_post_dr>0),'phase_jump_error_rad',max(abs(d.phase_delta-pi)),'q_reflect_error',max(abs(d.q_ref_error)),'p_rotation_error',max(d.p_rot_error),'q_rotation_error',max(d.q_rot_error),'tau_receiver_s',d.tau_receiver_real(find(abs(d.alpha_deg)==min(abs(d.alpha_deg)),1))));
comparison.checks=struct('profile_provenance',strcmp(profile.coeff_file_sha256,local_file_sha256(cfg.coeff_file)),'bellhop_geometry',comparison.geometry.wall_residual_max_m<=cfg.wall_residual_limit_m&&comparison.geometry.all_post_range_positive&&comparison.geometry.grazing_fraction==0,'reflection_phase',comparison.geometry.phase_jump_error_rad<=cfg.phase_jump_limit_rad,'beam_state',comparison.geometry.q_reflect_error<=cfg.state_limit&&comparison.geometry.p_rotation_error<=cfg.state_limit&&comparison.geometry.q_rotation_error<=cfg.state_limit,'field_finite',all(isfinite([real(Gpe) imag(Gpe) real(Gbh) imag(Gbh)])));
comparison.checks.all=all(structfun(@(v)logical(v),comparison.checks));comparison.passed=comparison.checks.all;comparison.files=local_write_outputs(comparison);disp(struct('G_pe',Gpe,'G_bellhop',Gbh,'delta_tl_db',comparison.delta_tl_db,'delta_phase_rad',comparison.delta_phase_rad,'delta_complex_relative_error',comparison.delta_complex_relative_error));disp(comparison.checks);
if cfg.fail_on_check&&~comparison.passed,error('Stage 1A Tier-1 audit failed; see %s.',cfg.report_path);end
end

function result=local_flat_bellhop(cfg,pattern)
case_root=cfg.reuse_flat_case_root;c=struct('bellhop_exe',cfg.flat_validation_exe,'case_root',case_root,'run_type','C','source_geometry',cfg.source_geometry,'frequency_hz',cfg.frequency_hz,'c0_mps',cfg.c0_mps,'source_depth_m',0,'receiver_depths_m',0,'receiver_ranges_m',[102 103],'beam_count',cfg.beam_count,'angle_limits_deg',cfg.angle_limits_deg,'domain_half_depth_m',cfg.domain_half_depth_m,'sigma_src_m',cfg.sigma_src_m,'source_pattern_clip_db',cfg.source_pattern_clip_db,'source_pattern_angles_deg',pattern.angles_deg,'source_pattern_level_db',pattern.level_db,'step_m',cfg.step_m,'wall_range_m',cfg.wall_r0_m,'mapped_receiver_range_m',cfg.mapped_receiver_range_m);
if exist([case_root '.shd'],'file')==2&&exist([case_root '.iwdiag'],'file')==2,result=local_read_flat_case(c,case_root);else,result=run_bellhop_internal_flat_wall_poc_vertical(c);end
end

function result=local_read_flat_case(c,case_root)
result.data=read_bellhop_shd_unfolded_vertical([case_root '.shd']);dm=readmatrix([case_root '.iwdiag'],'FileType','text','CommentStyle','#');dm=dm(~all(isnan(dm),2),:);
names={'alpha_deg','hit_r','hit_z','residual','inc_ur','inc_uz','ref_ur','ref_uz','rot_ur','rot_uz','phase_in','phase_ref','phase_delta','amp_in','amp_ref','amp_delta','p1_in','p2_in','p1_ref','p2_ref','p_ref_error','q1_in','q2_in','q1_ref','q2_ref','q_ref_error','p_rot_error','q_rot_error','tau_wall_real','tau_wall_imag','tau_receiver_real','tau_receiver_imag','min_post_dr','n_post','kappa'};
if size(dm,2)~=numel(names),error('Unexpected flat diagnostics in %s.',case_root);end
result.config=c;result.diagnostics=array2table(dm,'VariableNames',names);result.files=struct('data',[case_root '.shd'],'diagnostics',[case_root '.iwdiag']);
end

function result=local_rough_bellhop(cfg,pattern,profile)
n=cfg.profile_count;s=linspace(-profile.span_m/2,profile.span_m/2,n).';e=evaluate_fixed_pm_fourier_profile(profile,s);case_root=cfg.reuse_rough_case_root;
c=struct('bellhop_exe',cfg.pm_validation_exe,'case_root',case_root,'run_type','C','source_geometry',cfg.source_geometry,'frequency_hz',cfg.frequency_hz,'c0_mps',cfg.c0_mps,'source_depth_m',0,'receiver_depths_m',0,'receiver_ranges_m',[102 103],'beam_count',cfg.beam_count,'angle_limits_deg',cfg.angle_limits_deg,'domain_half_depth_m',cfg.domain_half_depth_m,'sigma_src_m',cfg.sigma_src_m,'source_pattern_clip_db',cfg.source_pattern_clip_db,'source_pattern_angles_deg',pattern.angles_deg,'source_pattern_level_db',pattern.level_db,'step_m',cfg.step_m,'wall_r0_m',cfg.wall_r0_m,'wall_seed',cfg.wall_seed,'wall_profile_r_m',cfg.wall_r0_m-e.eta_m(:).','wall_profile_z_m',s(:).','mapped_receiver_range_m',cfg.mapped_receiver_range_m);
if exist([case_root '.shd'],'file')==2&&exist([case_root '.iwdiag'],'file')==2,result=local_read_pm_case(c,case_root);else,result=run_bellhop_internal_pm_wall_poc_vertical(c);end
end

function result=local_read_pm_case(c,case_root)
result.data=read_bellhop_shd_unfolded_vertical([case_root '.shd']);dm=readmatrix([case_root '.iwdiag'],'FileType','text','CommentStyle','#');dm=dm(~all(isnan(dm),2),:);
names={'alpha_deg','hit_r','hit_z','wall_residual','wall_t_r','wall_t_z','wall_n_r','wall_n_z','tangent_error','normal_error','inc_ur','inc_uz','ref_ur','ref_uz','rot_ur','rot_uz','specular_error','rotation_error','phase_in','phase_ref','phase_delta','amp_in','amp_ref','amp_delta','p1_in','p2_in','p1_ref','p2_ref','p_ref_error','q1_in','q2_in','q1_ref','q2_ref','q_ref_error','p_rot_error','q_rot_error','tau_wall_real','tau_wall_imag','tau_receiver_real','tau_receiver_imag','min_post_dr','n_post','kappa','wall_seg','wall_lambda','wall_tg','wall_th','wall_rm','wall_rn'};
if size(dm,2)~=numel(names),error('Unexpected PM diagnostic columns in %s.',case_root);end
result.config=c;result.diagnostics=array2table(dm,'VariableNames',names);result.files=struct('data',[case_root '.shd'],'diagnostics',[case_root '.iwdiag']);
end

function cfg=local_config(root,o)
cfg=struct('frequency_hz',4000,'c0_mps',1500,'z_tx_m',100,'z_rx_m',3,'sigma_src_m',0.3,'xw_m',192.1875,'nx',984,'step_m',0.05,'wall_r0_m',100,'mapped_receiver_range_m',103,'wall_seed',260001,'profile_count',4097,'beam_count',10001,'angle_limits_deg',[-30 30],'domain_half_depth_m',1000,'source_pattern_clip_db',-120,'source_geometry','X','coeff_file',fullfile(root,'results','validation','bellhop_internal_pm_fixed_realization','fixed_pm_fourier_coefficients.csv'),'flat_validation_exe',fullfile(root,'results','validation','bellhop_internal_flat_wall_poc','bin','bellhop_iwall_flat_2020.exe'),'pm_validation_exe',fullfile(root,'results','validation','bellhop_internal_pm_wall_poc','bin','bellhop_iwall_pm_2020.exe'),'output_dir',fullfile(root,'results','validation','pe_bellhop_pm_stage1_tier1'),'report_path',fullfile(root,'reports','pe_bellhop_pm_stage1_tier1_report.md'),'reuse_flat_case_root',fullfile(root,'results','validation','pe_bellhop_pm_stage1_tier1','cases_X','flat'),'reuse_rough_case_root',fullfile(root,'results','validation','pe_bellhop_pm_stage1_tier1','cases_X','rough'),'selector_tolerance_m',1e-6,'grazing_mu_threshold',0.1,'wall_residual_limit_m',1e-9,'phase_jump_limit_rad',1e-10,'state_limit',1e-12,'fail_on_check',true);
names=fieldnames(o);for ii=1:numel(names),if ~isfield(cfg,names{ii}),error('Unknown Stage 1A override: %s.',names{ii});end;cfg.(names{ii})=o.(names{ii});end
cfg.source_geometry=upper(char(cfg.source_geometry));
if ~isscalar(cfg.source_geometry)||~ismember(cfg.source_geometry,['R','X']),error('source_geometry must be R (point) or X (line).');end
end

function p=local_source_pattern(cfg)
a=linspace(cfg.angle_limits_deg(1),cfg.angle_limits_deg(2),2401).';th=a*pi/180;k=2*pi*cfg.frequency_hz/cfg.c0_mps;d=cos(th).*exp(-0.5*(k*cfg.sigma_src_m*sin(th)).^2);d=abs(d)/max(abs(d));p=struct('angles_deg',a,'level_db',20*log10(max(d,10^(cfg.source_pattern_clip_db/20))));
end

function files=local_write_outputs(c)
out=c.config.output_dir;writetable(table(c.G_pe,c.G_bellhop,c.delta_G,c.delta_tl_db,c.delta_phase_rad,c.delta_complex_relative_error,'VariableNames',{'G_pe','G_bellhop','delta_G','delta_tl_db','delta_phase_rad','delta_complex_relative_error'}),fullfile(out,'stage1a_tier1_metrics.csv'));save(fullfile(out,'stage1a_tier1_comparison.mat'),'c','-v7.3');n=fieldnames(c.checks);v=false(size(n));for ii=1:numel(n),v(ii)=c.checks.(n{ii});end;writetable(table(n,v,'VariableNames',{'check_name','passed'}),fullfile(out,'stage1a_tier1_checks.csv'));local_write_report(c.config.report_path,c);files=struct('metrics',fullfile(out,'stage1a_tier1_metrics.csv'),'checks',fullfile(out,'stage1a_tier1_checks.csv'),'mat',fullfile(out,'stage1a_tier1_comparison.mat'),'report',c.config.report_path);
end

function local_write_report(path,c)
fid=fopen(path,'w','n','UTF-8');if fid<0,error('Cannot write %s.',path);end;cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'# PE--Bellhop PM Stage 1A Tier-1 reflected-ratio comparison\n\n状态：**%s**\n\n',ternary(c.passed,'PASS_WITH_LIMITS','FAIL'));
fprintf(fid,'固定条件：4 kHz、uniform c=1500 m/s、z_tx=100 m、z_rx=3 m、sigma=0.3 m、xw=192.1875 m、nx=984、step=0.05 m、seed=260001、profile N=4097、beam=10001、Bellhop source=%s。PE 使用 Kirchhoff phase screen，Bellhop 使用 native Reflect2D local-specular internal wall；只比较 reflected-only rough/flat ratios。\n\n',c.config.source_geometry);
fprintf(fid,'## Primary metrics\n\n| quantity | value |\n|---|---:|\n| G_PE | %.16g%+.16gi |\n| G_Bellhop | %.16g%+.16gi |\n| delta TL (dB) | %.8g |\n| delta phase (rad) | %.8g |\n| complex relative error | %.8g |\n',real(c.G_pe),imag(c.G_pe),real(c.G_bellhop),imag(c.G_bellhop),c.delta_tl_db,c.delta_phase_rad,c.delta_complex_relative_error);
fprintf(fid,'\n## Bellhop rough-wall diagnostics\n\nwall residual `%.8g m`; min `|u.n|` `%.8g`; grazing fraction `%.8g`; max `|kappa|` `%.8g 1/m`; min post-wall `dr` `%.8g m`; pressure-release phase error `%.8g rad`; q residual `%.8g`; p/q rotation errors `%.8g / %.8g`; tau `%.12g s`.\n\n',c.geometry.wall_residual_max_m,c.geometry.min_mu,c.geometry.grazing_fraction,c.geometry.kappa_abs_max_per_m,c.geometry.min_post_range_increment_m,c.geometry.phase_jump_error_rad,c.geometry.q_reflect_error,c.geometry.p_rotation_error,c.geometry.q_rotation_error,c.geometry.tau_receiver_s);
fprintf(fid,'## Checks\n\n');n=fieldnames(c.checks);for ii=1:numel(n),fprintf(fid,'- %s: %s\n',n{ii},ternary(c.checks.(n{ii}),'PASS','FAIL'));end
fprintf(fid,'\nThe cross-model delta is intentionally not a closeness gate: it quantifies the expected Kirchhoff phase-screen versus local-specular Gaussian-beam model discrepancy after each flat denominator removes source/absolute-amplitude factors.\n');
end
function s=ternary(c,a,b),if c,s=a;else,s=b;end,end
function t=local_tag(v),t=strrep(sprintf('%.6g',v),'.','p');t=strrep(t,'-','m');end
function d=local_file_sha256(path),md=java.security.MessageDigest.getInstance('SHA-256');fid=fopen(path,'rb');b=fread(fid,Inf,'*uint8');fclose(fid);md.update(typecast(b,'int8'));d=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));end
