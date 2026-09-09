function extension = validate_pe_bellhop_pm_frequency_extension(overrides)
%VALIDATE_PE_BELLHOP_PM_FREQUENCY_EXTENSION Stage 2 frequency extension.
if nargin<1 || isempty(overrides), overrides=struct(); end
if ~isstruct(overrides)||~isscalar(overrides), error('overrides must be a scalar struct.'); end
root=fileparts(fileparts(fileparts(mfilename('fullpath')))); addpath(root); setup_vertical_project();
cfg=local_config(root,overrides); if ~exist(cfg.output_dir,'dir'), mkdir(cfg.output_dir); end
profile=load_fixed_pm_profile_for_pe_bellhop_validation(cfg.coeff_file);
rows=repmat(local_empty_row(),1,numel(cfg.frequencies_hz));
for ii=1:numel(cfg.frequencies_hz)
    f=cfg.frequencies_hz(ii);
    x=(-0.5*cfg.xw_m)+(0:cfg.nx-1)*(cfg.xw_m/cfg.nx);
    e=evaluate_fixed_pm_fourier_profile(profile,x);
    pattern=local_source_pattern(cfg,f);
    base=struct('frequency_hz',f,'c0_mps',cfg.c0_mps,'xw_m',cfg.xw_m,'nx',cfg.nx, ...
        'z_tx_m',cfg.z_tx_m,'z_rx_m',cfg.z_rx_m,'sigma_src_m',cfg.sigma_src_m, ...
        'surface_reflect_coeff',-1,'step_m',cfg.step_m,'x_rx_m',cfg.x_rx_m);
    pe_flat=run_pe_1d_surface_reflection_validation(setfield(base,'surface_elevation_x_m',zeros(1,cfg.nx))); %#ok<SFLD>
    pe_rough=run_pe_1d_surface_reflection_validation(setfield(base,'surface_elevation_x_m',e.eta_m(:).')); %#ok<SFLD>
    bh_flat=local_run_flat_bh(cfg,f,pattern);
    bh_rough=local_run_rough_bh(cfg,f,pattern,profile);
    pflat=select_bellhop_shd_pressure_at_range_vertical(bh_flat.data,cfg.mapped_receiver_range_m,0,cfg.selector_tolerance_m);
    prough=select_bellhop_shd_pressure_at_range_vertical(bh_rough.data,cfg.mapped_receiver_range_m,0,cfg.selector_tolerance_m);
    Gpe=pe_rough.reflected_receiver/pe_flat.reflected_receiver;
    Gbh=prough/pflat; d=bh_rough.diagnostics;
    mu=abs(d.inc_ur.*d.wall_n_r+d.inc_uz.*d.wall_n_z);
    r=local_empty_row(); r.frequency_hz=f; r.G_pe=Gpe; r.G_bellhop=Gbh;
    r.delta_tl_db=20*log10(abs(Gpe/Gbh)); r.delta_phase_rad=angle(Gpe*conj(Gbh));
    r.complex_relative_error=abs(Gpe-Gbh)/max(abs(Gbh),realmin);
    r.pe_seam_jump_m=abs(e.eta_m(1)-e.eta_m(end)); nedge=max(1,round(.05*cfg.nx));
    edge=[1:nedge,(cfg.nx-nedge+1):cfg.nx];
    r.pe_outer5_energy_fraction=sum(abs(pe_rough.surface_reflected_field(edge)).^2)/max(sum(abs(pe_rough.surface_reflected_field).^2),realmin);
    r.wall_residual_max_m=max(abs(d.wall_residual)); r.min_mu=min(mu); r.grazing_fraction=mean(mu<cfg.grazing_mu_threshold);
    r.kappa_abs_max_per_m=max(abs(d.kappa)); r.phase_jump_error_rad=max(abs(d.phase_delta-pi));
    r.q_reflect_error=max(abs(d.q_ref_error)); r.p_rotation_error=max(d.p_rot_error); r.q_rotation_error=max(d.q_rot_error);
    r.min_post_range_increment_m=min(d.min_post_dr); r.all_post_range_positive=all(d.min_post_dr>0);
    r.tau_receiver_s=d.tau_receiver_real(find(abs(d.alpha_deg)==min(abs(d.alpha_deg)),1));
    r.flat=bh_flat; r.rough=bh_rough; r.pe_flat=pe_flat; r.pe_rough=pe_rough; rows(ii)=r;
end
checks=struct('profile_provenance',strcmp(profile.coeff_file_sha256,local_file_sha256(cfg.coeff_file)), ...
    'pe_finite',all(arrayfun(@(r)isfinite(real(r.G_pe))&&isfinite(imag(r.G_pe)),rows)), ...
    'bellhop_finite',all(arrayfun(@(r)isfinite(real(r.G_bellhop))&&isfinite(imag(r.G_bellhop)),rows)), ...
    'bellhop_geometry',all([rows.wall_residual_max_m]<=cfg.wall_residual_limit_m)&&all([rows.all_post_range_positive])&&all([rows.grazing_fraction]==0), ...
    'reflection_phase',max([rows.phase_jump_error_rad])<=cfg.phase_jump_limit_rad, ...
    'beam_state',max([rows.q_reflect_error])<=cfg.state_limit&&max([rows.p_rotation_error])<=cfg.state_limit&&max([rows.q_rotation_error])<=cfg.state_limit);
checks.all=all(structfun(@(v)logical(v),checks));
extension=struct('schema_version','1.0.0','stage','2_frequency_extension','config',cfg,'profile',profile,'rows',rows,'checks',checks,'passed',checks.all);
extension.files=local_write_outputs(extension);
if cfg.fail_on_check&&~extension.passed, error('Stage 2 frequency extension failed; see %s.',cfg.report_path); end
end

function row=local_empty_row()
row=struct('frequency_hz',NaN,'G_pe',complex(NaN),'G_bellhop',complex(NaN),'delta_tl_db',NaN,'delta_phase_rad',NaN,'complex_relative_error',NaN,'pe_seam_jump_m',NaN,'pe_outer5_energy_fraction',NaN,'wall_residual_max_m',NaN,'min_mu',NaN,'grazing_fraction',NaN,'kappa_abs_max_per_m',NaN,'phase_jump_error_rad',NaN,'q_reflect_error',NaN,'p_rotation_error',NaN,'q_rotation_error',NaN,'min_post_range_increment_m',NaN,'all_post_range_positive',false,'tau_receiver_s',NaN,'flat',struct(),'rough',struct(),'pe_flat',struct(),'pe_rough',struct());
end

function files=local_write_outputs(a)
out=a.config.output_dir; summary=rmfield(a.rows,{'flat','rough','pe_flat','pe_rough'});
writetable(struct2table(summary),fullfile(out,'frequency_summary.csv'));
n=fieldnames(a.checks); v=false(size(n)); for ii=1:numel(n),v(ii)=a.checks.(n{ii});end
writetable(table(n,v,'VariableNames',{'check_name','passed'}),fullfile(out,'frequency_checks.csv'));
save(fullfile(out,'frequency_extension.mat'),'a','-v7.3'); local_write_report(a.config.report_path,a);
files=struct('summary',fullfile(out,'frequency_summary.csv'),'checks',fullfile(out,'frequency_checks.csv'),'mat',fullfile(out,'frequency_extension.mat'),'report',a.config.report_path);
end

function local_write_report(path,a)
fid=fopen(path,'w','n','UTF-8'); if fid<0,error('Cannot write %s.',path);end; cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'# PE--Bellhop fixed-PM frequency extension\n\n状态：**%s**\n\n',ternary(a.passed,'PASS_WITH_MODEL_DISCREPANCY','FAIL'));
fprintf(fid,'固定 seed=260001、U=6 m/s、span=160 m、master N=4097、realized Kmax=0.471238898 rad/m；uniform c=1500 m/s、z_tx=100 m、z_rx=3 m、sigma=0.3 m、W=192.1875 m、nx=984、step=0.05 m、Bellhop source=%s、beam=%d、sector=[%.0f,%.0f] deg。PE 与 Bellhop 核心均未修改。\n\n',a.config.source_geometry,a.config.beam_count,a.config.angle_limits_deg(1),a.config.angle_limits_deg(2));
fprintf(fid,'## Frequency results\n\n| f (kHz) | G_PE | G_Bellhop | delta TL (dB) | delta phase (rad) | complex error | wall residual (m) | min |u.n| | max |kappa| (1/m) | tau (s) |\n|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n');
for ii=1:numel(a.rows),r=a.rows(ii); fprintf(fid,'| %.0f | %.12g%+.12gi | %.12g%+.12gi | %.8g | %.8g | %.8g | %.8g | %.8g | %.8g | %.12g |\n',r.frequency_hz/1000,real(r.G_pe),imag(r.G_pe),real(r.G_bellhop),imag(r.G_bellhop),r.delta_tl_db,r.delta_phase_rad,r.complex_relative_error,r.wall_residual_max_m,r.min_mu,r.kappa_abs_max_per_m,r.tau_receiver_s); end
fprintf(fid,'\n## Native numerical diagnostics\n\n'); for ii=1:numel(a.rows),r=a.rows(ii); fprintf(fid,'- %.0f kHz: PE seam jump %.8g m, outer-5%% reflected energy %.8g; Bellhop grazing fraction %.8g, phase-jump error %.8g rad, q residual %.8g, p/q rotation errors %.8g / %.8g, min post-wall dr %.8g m.\n',r.frequency_hz/1000,r.pe_seam_jump_m,r.pe_outer5_energy_fraction,r.grazing_fraction,r.phase_jump_error_rad,r.q_reflect_error,r.p_rotation_error,r.q_rotation_error,r.min_post_range_increment_m); end
fprintf(fid,'\n## Checks\n\n'); n=fieldnames(a.checks); for ii=1:numel(n),fprintf(fid,'- %s: %s\n',n{ii},ternary(a.checks.(n{ii}),'PASS','FAIL'));end
fprintf(fid,'\nThe cross-model residual is intentionally not forced toward zero. This stage tests finite, geometrically valid operation at 4/6/8 kHz with the same PM realization; no group delay is inferred from these three frequencies.\n');
end
function s=ternary(c,a,b),if c,s=a;else,s=b;end,end
function d=local_file_sha256(path),md=java.security.MessageDigest.getInstance('SHA-256');fid=fopen(path,'rb');b=fread(fid,Inf,'*uint8');fclose(fid);md.update(typecast(b,'int8'));d=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));end

function result=local_run_flat_bh(cfg,f,pattern)
case_root=fullfile(cfg.output_dir,'cases',sprintf('f%dk_%s_B%d_flat',round(f),cfg.source_geometry,cfg.beam_count));
c=struct('bellhop_exe',cfg.flat_validation_exe,'case_root',case_root,'run_type','C','source_geometry',cfg.source_geometry,'frequency_hz',f,'c0_mps',cfg.c0_mps,'source_depth_m',0,'receiver_depths_m',0,'receiver_ranges_m',[102 103],'beam_count',cfg.beam_count,'angle_limits_deg',cfg.angle_limits_deg,'domain_half_depth_m',cfg.domain_half_depth_m,'sigma_src_m',cfg.sigma_src_m,'source_pattern_clip_db',cfg.source_pattern_clip_db,'source_pattern_angles_deg',pattern.angles_deg,'source_pattern_level_db',pattern.level_db,'step_m',cfg.step_m,'wall_range_m',cfg.wall_r0_m,'mapped_receiver_range_m',cfg.mapped_receiver_range_m);
if exist([case_root '.shd'],'file')==2&&exist([case_root '.iwdiag'],'file')==2
    try, result=local_read_flat(c,case_root); catch, result=run_bellhop_internal_flat_wall_poc_vertical(c); end
else, result=run_bellhop_internal_flat_wall_poc_vertical(c); end
end

function result=local_run_rough_bh(cfg,f,pattern,profile)
s=linspace(-profile.span_m/2,profile.span_m/2,cfg.profile_count).'; e=evaluate_fixed_pm_fourier_profile(profile,s);
case_root=fullfile(cfg.output_dir,'cases',sprintf('f%dk_%s_B%d_rough',round(f),cfg.source_geometry,cfg.beam_count));
c=struct('bellhop_exe',cfg.pm_validation_exe,'case_root',case_root,'run_type','C','source_geometry',cfg.source_geometry,'frequency_hz',f,'c0_mps',cfg.c0_mps,'source_depth_m',0,'receiver_depths_m',0,'receiver_ranges_m',[102 103],'beam_count',cfg.beam_count,'angle_limits_deg',cfg.angle_limits_deg,'domain_half_depth_m',cfg.domain_half_depth_m,'sigma_src_m',cfg.sigma_src_m,'source_pattern_clip_db',cfg.source_pattern_clip_db,'source_pattern_angles_deg',pattern.angles_deg,'source_pattern_level_db',pattern.level_db,'step_m',cfg.step_m,'wall_r0_m',cfg.wall_r0_m,'wall_seed',cfg.wall_seed,'wall_profile_r_m',cfg.wall_r0_m-e.eta_m(:).','wall_profile_z_m',s(:).','mapped_receiver_range_m',cfg.mapped_receiver_range_m);
if exist([case_root '.shd'],'file')==2&&exist([case_root '.iwdiag'],'file')==2
    try, result=local_read_pm(c,case_root); catch, result=run_bellhop_internal_pm_wall_poc_vertical(c); end
else, result=run_bellhop_internal_pm_wall_poc_vertical(c); end
end

function result=local_read_flat(c,root)
result.data=read_bellhop_shd_unfolded_vertical([root '.shd']); dm=readmatrix([root '.iwdiag'],'FileType','text','CommentStyle','#'); dm=dm(~all(isnan(dm),2),:);
names={'alpha_deg','hit_r','hit_z','residual','inc_ur','inc_uz','ref_ur','ref_uz','rot_ur','rot_uz','phase_in','phase_ref','phase_delta','amp_in','amp_ref','amp_delta','p1_in','p2_in','p1_ref','p2_ref','p_ref_error','q1_in','q2_in','q1_ref','q2_ref','q_ref_error','p_rot_error','q_rot_error','tau_wall_real','tau_wall_imag','tau_receiver_real','tau_receiver_imag','min_post_dr','n_post','kappa'};
if size(dm,2)~=numel(names), error('Unexpected flat diagnostics in %s.',root); end; result.config=c; result.diagnostics=array2table(dm,'VariableNames',names); result.files=struct('data',[root '.shd'],'diagnostics',[root '.iwdiag']);
end

function result=local_read_pm(c,root)
result.data=read_bellhop_shd_unfolded_vertical([root '.shd']); dm=readmatrix([root '.iwdiag'],'FileType','text','CommentStyle','#'); dm=dm(~all(isnan(dm),2),:);
names={'alpha_deg','hit_r','hit_z','wall_residual','wall_t_r','wall_t_z','wall_n_r','wall_n_z','tangent_error','normal_error','inc_ur','inc_uz','ref_ur','ref_uz','rot_ur','rot_uz','specular_error','rotation_error','phase_in','phase_ref','phase_delta','amp_in','amp_ref','amp_delta','p1_in','p2_in','p1_ref','p2_ref','p_ref_error','q1_in','q2_in','q1_ref','q2_ref','q_ref_error','p_rot_error','q_rot_error','tau_wall_real','tau_wall_imag','tau_receiver_real','tau_receiver_imag','min_post_dr','n_post','kappa','wall_seg','wall_lambda','wall_tg','wall_th','wall_rm','wall_rn'};
if size(dm,2)~=numel(names), error('Unexpected PM diagnostics in %s.',root); end; result.config=c; result.diagnostics=array2table(dm,'VariableNames',names); result.files=struct('data',[root '.shd'],'diagnostics',[root '.iwdiag']);
end

function pattern=local_source_pattern(cfg,f)
a=linspace(cfg.angle_limits_deg(1),cfg.angle_limits_deg(2),2401).'; th=a*pi/180; k=2*pi*f/cfg.c0_mps;
d=cos(th).*exp(-.5*(k*cfg.sigma_src_m*sin(th)).^2); d=abs(d)/max(abs(d)); pattern=struct('angles_deg',a,'level_db',20*log10(max(d,10^(cfg.source_pattern_clip_db/20))));
end

function cfg=local_config(root,o)
cfg=struct('frequencies_hz',[4000 6000 8000],'c0_mps',1500,'z_tx_m',100,'z_rx_m',3,'sigma_src_m',.3,'xw_m',192.1875,'nx',984,'x_rx_m',0,'step_m',.05,'wall_r0_m',100,'mapped_receiver_range_m',103,'wall_seed',260001,'profile_count',4097,'beam_count',5001,'angle_limits_deg',[-15 15],'domain_half_depth_m',1000,'source_pattern_clip_db',-120,'source_geometry','X','selector_tolerance_m',1e-4,'grazing_mu_threshold',.1,'wall_residual_limit_m',1e-9,'phase_jump_limit_rad',1e-10,'state_limit',1e-12,'coeff_file',fullfile(root,'results','validation','bellhop_internal_pm_fixed_realization','fixed_pm_fourier_coefficients.csv'),'flat_validation_exe',fullfile(root,'results','validation','bellhop_internal_flat_wall_poc','bin','bellhop_iwall_flat_2020.exe'),'pm_validation_exe',fullfile(root,'results','validation','bellhop_internal_pm_wall_poc','bin','bellhop_iwall_pm_2020.exe'),'output_dir',fullfile(root,'results','validation','pe_bellhop_pm_frequency_extension'),'report_path',fullfile(root,'reports','pe_bellhop_fixed_pm_frequency_extension_report.md'),'fail_on_check',true);
names=fieldnames(o); for ii=1:numel(names), if ~isfield(cfg,names{ii}), error('Unknown Stage 2 override: %s.',names{ii}); end; cfg.(names{ii})=o.(names{ii}); end
cfg.source_geometry=upper(char(cfg.source_geometry));
if ~isscalar(cfg.source_geometry)||~ismember(cfg.source_geometry,['R','X'])
    error('source_geometry must be R (point) or X (line).');
end
end
