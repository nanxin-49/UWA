function ensemble = validate_pe_bellhop_pm_ensemble(overrides)
%VALIDATE_PE_BELLHOP_PM_ENSEMBLE Stage 3 paired PM statistics.
%   Uses the canonical fixed spectral band.  The default phase mode preserves
%   the canonical per-mode amplitudes; ensemble_mode='amplitude' instead
%   draws independent Gaussian cosine/sine coefficients from the canonical
%   PM spectral density. Each seed is sent to unchanged 1-D PE and
%   validation-only Bellhop internal-wall propagation at 4 kHz.
if nargin<1||isempty(overrides),overrides=struct();end
if ~isstruct(overrides)||~isscalar(overrides),error('overrides must be a scalar struct.');end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));addpath(root);setup_vertical_project();
cfg=local_config(root,overrides);if ~exist(cfg.output_dir,'dir'),mkdir(cfg.output_dir);end
reference=load_fixed_pm_profile_for_pe_bellhop_validation(cfg.reference_coeff_file);
stage2=local_load_stage2(cfg.stage2_mat);
pattern=local_source_pattern(cfg);
flat_pe=local_run_flat_pe(cfg);
[flat_bh,flat_case_meta]=local_resolve_flat_bellhop(cfg,stage2,pattern);
pflat=select_bellhop_shd_pressure_at_range_vertical(flat_bh.data,cfg.mapped_receiver_range_m,0,cfg.selector_tolerance_m);
seeds=cfg.seeds(:).'; rows=repmat(local_empty_row(),1,numel(seeds));
for ii=1:numel(seeds)
    seed=seeds(ii); coeff_file=fullfile(cfg.output_dir,sprintf('seed_%d_coefficients.csv',seed));
    if strcmpi(cfg.ensemble_mode,'phase')
        profile=build_paired_pm_phase_realization(reference,seed,coeff_file);
    elseif strcmpi(cfg.ensemble_mode,'amplitude')
        profile=build_paired_pm_gaussian_realization(reference,seed,coeff_file);
    else
        error('Unsupported ensemble_mode: %s.',cfg.ensemble_mode);
    end
    band_match=local_band_match(reference,profile,coeff_file);
    x=(-0.5*cfg.xw_m)+(0:cfg.nx-1)*(cfg.xw_m/cfg.nx); ex=evaluate_fixed_pm_fourier_profile(profile,x);
    base=struct('frequency_hz',cfg.frequency_hz,'c0_mps',cfg.c0_mps,'xw_m',cfg.xw_m,'nx',cfg.nx,'z_tx_m',cfg.z_tx_m,'z_rx_m',cfg.z_rx_m,'sigma_src_m',cfg.sigma_src_m,'surface_elevation_x_m',ex.eta_m(:).','surface_reflect_coeff',-1,'step_m',cfg.step_m,'x_rx_m',0);
    pe=run_pe_1d_surface_reflection_validation(base); Gpe=pe.reflected_receiver/flat_pe.reflected_receiver;
    pe_edge=local_pe_edge_metrics(pe,ex.eta_m,cfg);
    s=linspace(-profile.span_m/2,profile.span_m/2,cfg.profile_count).'; es=evaluate_fixed_pm_fourier_profile(profile,s);
    c=local_rough_bellhop_config(cfg,pattern,profile,es,seed,s);
    [bh,case_meta]=local_resolve_rough_bellhop(cfg,stage2,c,profile,seed,reference.seed);
    p=select_bellhop_shd_pressure_at_range_vertical(bh.data,cfg.mapped_receiver_range_m,0,cfg.selector_tolerance_m); Gbh=p/pflat; d=bh.diagnostics; mu=abs(d.inc_ur.*d.wall_n_r+d.inc_uz.*d.wall_n_z);
    kappa=es.eta_second_m_per_m2./(1+es.eta_prime_m_per_m.^2).^(3/2); r=local_empty_row();r.seed=seed;r.coeff_file_sha256=profile.coeff_file_sha256;r.profile_band_match=band_match;r.surface_mean_eta_m=mean(es.eta_m);r.surface_rms_eta_m=sqrt(mean(es.eta_m.^2));r.surface_sigma_eta_m=std(es.eta_m,1);r.surface_min_eta_m=min(es.eta_m);r.surface_max_eta_m=max(es.eta_m);r.profile_rms_height_m=r.surface_rms_eta_m;r.profile_rms_slope=sqrt(mean(es.eta_prime_m_per_m.^2));r.profile_max_slope=max(abs(es.eta_prime_m_per_m));r.profile_rms_curvature_per_m=sqrt(mean(kappa.^2));r.profile_max_curvature_per_m=max(abs(kappa));r.profile_min_radius_m=1/max(max(abs(kappa)),eps);r.two_k_sigma_eta=2*(2*pi*cfg.frequency_hz/cfg.c0_mps)*r.surface_sigma_eta_m;r.kmax_over_k=max(profile.k_rad_per_m)/(2*pi*cfg.frequency_hz/cfg.c0_mps);r.G_pe=Gpe;r.G_bellhop=Gbh;r.G_pe_abs=abs(Gpe);r.G_bellhop_abs=abs(Gbh);r.G_pe_power=abs(Gpe)^2;r.G_bellhop_power=abs(Gbh)^2;r.delta_tl_db=20*log10(abs(Gpe/Gbh));r.delta_power_db=10*log10(max(r.G_pe_power,realmin)/max(r.G_bellhop_power,realmin));r.delta_phase_rad=angle(Gpe*conj(Gbh));r.complex_relative_error=abs(Gpe-Gbh)/max(abs(Gbh),realmin);r.pe_seam_jump_m=pe_edge.seam_jump_m;r.pe_outer5_incident_energy_fraction=pe_edge.outer5_incident_energy_fraction;r.pe_outer5_surface_reflected_energy_fraction=pe_edge.outer5_surface_reflected_energy_fraction;r.pe_outer5_receiver_reflected_energy_fraction=pe_edge.outer5_receiver_reflected_energy_fraction;r.wall_residual_max_m=max(abs(d.wall_residual));r.min_mu=min(mu);r.mean_mu=mean(mu);r.mu_p05=local_quantile(mu,0.05);r.grazing_fraction=mean(mu<cfg.grazing_mu_threshold);r.wall_hit_count=height(d);r.expected_beam_count=cfg.beam_count;r.failed_ray_count=max(0,cfg.beam_count-r.wall_hit_count);r.rejected_ray_count=0;r.successful_wall_reflection_fraction=r.wall_hit_count/r.expected_beam_count;r.hit_slope_mean=mean(d.wall_t_r./max(abs(d.wall_t_z),realmin));r.hit_curvature_mean=mean(abs(d.kappa));r.hit_curvature_rms=sqrt(mean(d.kappa.^2));r.hit_curvature_p90=local_quantile(abs(d.kappa),0.9);r.hit_curvature_p95=local_quantile(abs(d.kappa),0.95);r.incidence_angle_mean_deg=mean(acosd(min(max(mu,-1),1)));r.incidence_angle_std_deg=std(acosd(min(max(mu,-1),1)),1);r.bellhop_case_source=case_meta.source;r.bellhop_case_fingerprint=case_meta.fingerprint_sha256;r.bellhop_source_geometry=cfg.source_geometry;r.bellhop_beam_count=bh.config.beam_count;r.phase_jump_error_rad=max(abs(d.phase_delta-pi));r.q_reflect_error=max(abs(d.q_ref_error));r.p_rotation_error=max(d.p_rot_error);r.q_rotation_error=max(d.q_rot_error);r.min_post_range_increment_m=min(d.min_post_dr);r.all_post_range_positive=all(d.min_post_dr>0);r.tau_receiver_s=d.tau_receiver_real(find(abs(d.alpha_deg)==min(abs(d.alpha_deg)),1));rows(ii)=r;
    writetable(struct2table(r),fullfile(cfg.output_dir,sprintf('seed_%d_metrics.csv',seed)));
end
stats=local_statistics(rows);
checks=struct('reference_band',all([rows.profile_band_match])&&all(isfinite([rows.profile_max_curvature_per_m])), ...
    'profile_provenance',all(~cellfun(@isempty,{rows.coeff_file_sha256})), ...
    'case_fingerprints',all(~cellfun(@isempty,{rows.bellhop_case_fingerprint})), ...
    'uniform_bellhop_numerics',flat_bh.config.beam_count==cfg.beam_count&&all([rows.bellhop_beam_count]==cfg.beam_count), ...
    'reflection_success',all([rows.wall_hit_count]==cfg.beam_count)&&all(abs([rows.successful_wall_reflection_fraction]-1)<=eps), ...
    'fields_finite',all(arrayfun(@(r)isfinite(real(r.G_pe))&&isfinite(imag(r.G_pe))&&isfinite(real(r.G_bellhop))&&isfinite(imag(r.G_bellhop)),rows)), ...
    'bellhop_geometry',all([rows.wall_residual_max_m]<=cfg.wall_residual_limit_m)&&all([rows.all_post_range_positive])&&all([rows.grazing_fraction]==0)&&all([rows.failed_ray_count]==0), ...
    'reflection_phase',max([rows.phase_jump_error_rad])<=cfg.phase_jump_limit_rad, ...
    'beam_state',max([rows.q_reflect_error])<=cfg.state_limit&&max([rows.p_rotation_error])<=cfg.state_limit&&max([rows.q_rotation_error])<=cfg.state_limit, ...
    'sample_count',numel(rows)>=cfg.minimum_seed_count);
checks.all=all(structfun(@(v)logical(v),checks)); classification=ternary(checks.all,'PASS_WITH_LIMITS','FAIL');
ensemble=struct('schema_version','1.1.0','stage','3_pm_ensemble','config',cfg,'reference_profile',reference,'flat_case_meta',flat_case_meta,'rows',rows,'statistics',stats,'checks',checks,'classification',classification,'passed',checks.all);
ensemble.files=local_write_outputs(ensemble);
if cfg.fail_on_check&&~ensemble.passed,error('Stage 3 PM ensemble validation failed; see %s.',cfg.report_path);end
end

function pe=local_run_flat_pe(cfg)
base=struct('frequency_hz',cfg.frequency_hz,'c0_mps',cfg.c0_mps, ...
    'xw_m',cfg.xw_m,'nx',cfg.nx,'z_tx_m',cfg.z_tx_m,'z_rx_m',cfg.z_rx_m, ...
    'sigma_src_m',cfg.sigma_src_m,'surface_elevation_x_m',zeros(1,cfg.nx), ...
    'surface_reflect_coeff',-1,'step_m',cfg.step_m,'x_rx_m',0);
pe=run_pe_1d_surface_reflection_validation(base);
end

function m=local_pe_edge_metrics(pe,eta,cfg)
n=max(1,round(0.05*cfg.nx)); edge=[1:n,(cfg.nx-n+1):cfg.nx];
m=struct('seam_jump_m',abs(eta(1)-eta(end)), ...
    'outer5_incident_energy_fraction',sum(abs(pe.surface_incident_field(edge)).^2)/max(sum(abs(pe.surface_incident_field).^2),realmin), ...
    'outer5_surface_reflected_energy_fraction',sum(abs(pe.surface_reflected_field(edge)).^2)/max(sum(abs(pe.surface_reflected_field).^2),realmin), ...
    'outer5_receiver_reflected_energy_fraction',sum(abs(pe.reflected_field(edge)).^2)/max(sum(abs(pe.reflected_field).^2),realmin));
end

function c=local_flat_bellhop_config(cfg,pattern)
c=struct('bellhop_exe',cfg.flat_validation_exe, ...
    'case_root',fullfile(cfg.case_cache_dir,sprintf('flat_%s_B%d',cfg.source_geometry,cfg.beam_count)), ...
    'run_type','C','source_geometry',cfg.source_geometry, ...
    'frequency_hz',cfg.frequency_hz,'c0_mps',cfg.c0_mps, ...
    'source_depth_m',0,'receiver_depths_m',0,'receiver_ranges_m',[102 103], ...
    'beam_count',cfg.beam_count,'angle_limits_deg',cfg.angle_limits_deg, ...
    'domain_half_depth_m',cfg.domain_half_depth_m,'sigma_src_m',cfg.sigma_src_m, ...
    'source_pattern_clip_db',cfg.source_pattern_clip_db, ...
    'source_pattern_angles_deg',pattern.angles_deg,'source_pattern_level_db',pattern.level_db, ...
    'step_m',cfg.step_m,'wall_range_m',100, ...
    'mapped_receiver_range_m',cfg.mapped_receiver_range_m);
end

function c=local_rough_bellhop_config(cfg,pattern,profile,es,seed,s)
c=struct('bellhop_exe',cfg.pm_validation_exe, ...
    'case_root',fullfile(cfg.case_cache_dir,sprintf('seed_%d_%s_B%d',seed,cfg.source_geometry,cfg.beam_count)), ...
    'run_type','C','source_geometry',cfg.source_geometry, ...
    'frequency_hz',cfg.frequency_hz,'c0_mps',cfg.c0_mps, ...
    'source_depth_m',0,'receiver_depths_m',0,'receiver_ranges_m',[102 103], ...
    'beam_count',cfg.beam_count,'angle_limits_deg',cfg.angle_limits_deg, ...
    'domain_half_depth_m',cfg.domain_half_depth_m,'sigma_src_m',cfg.sigma_src_m, ...
    'source_pattern_clip_db',cfg.source_pattern_clip_db, ...
    'source_pattern_angles_deg',pattern.angles_deg,'source_pattern_level_db',pattern.level_db, ...
    'step_m',cfg.step_m,'wall_r0_m',100,'wall_seed',seed, ...
    'wall_profile_r_m',100-es.eta_m(:).','wall_profile_z_m',s(:).', ...
    'mapped_receiver_range_m',cfg.mapped_receiver_range_m, ...
    'profile_coeff_sha256',profile.coeff_file_sha256);
end

function [bh,meta]=local_resolve_flat_bellhop(cfg,stage2,pattern)
c=local_flat_bellhop_config(cfg,pattern);
fingerprint=local_case_fingerprint(c,'flat','');
if numel(stage2.rows)>=1 && local_stage2_result_matches( ...
        stage2.rows(1).flat,c,'flat','',cfg.stage2_mat,{c.bellhop_exe})
    bh=stage2.rows(1).flat;
    source='stage2_verified';
elseif local_cached_case_valid(c.case_root,fingerprint)
    bh=local_read_flat(c,c.case_root);
    source='cache_verified';
else
    bh=run_bellhop_internal_flat_wall_poc_vertical(c);
    local_write_case_manifest(c.case_root,fingerprint);
    source='fresh_run';
end
meta=struct('source',source,'fingerprint_sha256',fingerprint.fingerprint_sha256, ...
    'beam_count',bh.config.beam_count,'case_root',c.case_root);
end

function [bh,meta]=local_resolve_rough_bellhop(cfg,stage2,c,profile,seed,reference_seed)
fingerprint=local_case_fingerprint(c,'rough',profile.coeff_file_sha256);
can_reuse_stage2=seed==reference_seed && numel(stage2.rows)>=1 && ...
    local_stage2_result_matches(stage2.rows(1).rough,c,'rough', ...
        profile.coeff_file_sha256,cfg.stage2_mat,{c.bellhop_exe,profile.coeff_file});
if can_reuse_stage2
    bh=stage2.rows(1).rough;
    source='stage2_verified';
elseif local_cached_case_valid(c.case_root,fingerprint)
    bh=local_read_pm(c,c.case_root);
    source='cache_verified';
else
    bh=run_bellhop_internal_pm_wall_poc_vertical(c);
    local_write_case_manifest(c.case_root,fingerprint);
    source='fresh_run';
end
meta=struct('source',source,'fingerprint_sha256',fingerprint.fingerprint_sha256, ...
    'beam_count',bh.config.beam_count,'case_root',c.case_root);
end

function tf=local_stage2_result_matches(result,c,kind,profile_hash,stage2_mat,input_files)
tf=isstruct(result)&&isfield(result,'config')&&isfield(result,'data')&& ...
    isfield(result,'diagnostics')&&local_artifact_newer_than_inputs(stage2_mat,input_files);
if ~tf,return;end
try
    requested=local_case_fingerprint(c,kind,profile_hash);
    stored=local_case_fingerprint(result.config,kind,profile_hash);
    tf=strcmp(requested.fingerprint_sha256,stored.fingerprint_sha256);
catch
    tf=false;
end
end

function tf=local_artifact_newer_than_inputs(artifact,input_files)
if exist(artifact,'file')~=2,tf=false;return;end
a=dir(artifact);tf=true;
for ii=1:numel(input_files)
    if exist(input_files{ii},'file')~=2,tf=false;return;end
    d=dir(input_files{ii});
    if a.datenum<d.datenum,tf=false;return;end %#ok<DATNM>
end
end

function fingerprint=local_case_fingerprint(c,kind,profile_hash)
summary=struct('schema_version','1.0.0','case_kind',kind, ...
    'run_type',upper(char(c.run_type(1))),'source_geometry',upper(char(c.source_geometry)), ...
    'frequency_hz',double(c.frequency_hz), ...
    'c0_mps',double(c.c0_mps),'source_depth_m',double(c.source_depth_m), ...
    'receiver_depths_sha256',local_numeric_sha256(c.receiver_depths_m), ...
    'receiver_ranges_sha256',local_numeric_sha256(c.receiver_ranges_m), ...
    'beam_count',double(c.beam_count),'angle_limits_sha256',local_numeric_sha256(c.angle_limits_deg), ...
    'domain_half_depth_m',double(c.domain_half_depth_m),'sigma_src_m',double(c.sigma_src_m), ...
    'source_pattern_clip_db',double(c.source_pattern_clip_db),'step_m',double(c.step_m), ...
    'mapped_receiver_range_m',double(c.mapped_receiver_range_m), ...
    'executable_sha256',local_file_sha256(c.bellhop_exe), ...
    'source_angles_sha256',local_numeric_sha256(c.source_pattern_angles_deg), ...
    'source_levels_sha256',local_numeric_sha256(c.source_pattern_level_db), ...
    'profile_coeff_sha256',char(profile_hash));
if strcmp(kind,'flat')
    summary.wall_range_m=double(c.wall_range_m);
else
    summary.wall_r0_m=double(c.wall_r0_m);
    summary.wall_seed=double(c.wall_seed);
    summary.wall_r_sha256=local_numeric_sha256(c.wall_profile_r_m);
    summary.wall_z_sha256=local_numeric_sha256(c.wall_profile_z_m);
end
fingerprint=summary;
fingerprint.fingerprint_sha256=local_text_sha256(jsonencode(summary));
end

function tf=local_cached_case_valid(case_root,fingerprint)
manifest_file=[case_root '.case_manifest.json'];
data_file=[case_root '.shd'];diag_file=[case_root '.iwdiag'];
tf=exist(manifest_file,'file')==2&&exist(data_file,'file')==2&&exist(diag_file,'file')==2;
if ~tf,return;end
try
    manifest=jsondecode(fileread(manifest_file));
    tf=isfield(manifest,'fingerprint_sha256')&& ...
        strcmp(manifest.fingerprint_sha256,fingerprint.fingerprint_sha256)&& ...
        isfield(manifest,'output_shd_sha256')&& ...
        strcmp(manifest.output_shd_sha256,local_file_sha256(data_file))&& ...
        isfield(manifest,'output_iwdiag_sha256')&& ...
        strcmp(manifest.output_iwdiag_sha256,local_file_sha256(diag_file));
catch
    tf=false;
end
end

function local_write_case_manifest(case_root,fingerprint)
manifest=fingerprint;
manifest.output_shd_sha256=local_file_sha256([case_root '.shd']);
manifest.output_iwdiag_sha256=local_file_sha256([case_root '.iwdiag']);
manifest.created_utc=char(datetime('now','TimeZone','UTC', ...
    'Format','yyyy-MM-dd''T''HH:mm:ss''Z'''));
manifest_file=[case_root '.case_manifest.json'];
fid=fopen(manifest_file,'w','n','UTF-8');
if fid<0,error('Cannot write case manifest: %s.',manifest_file);end
cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(manifest,'PrettyPrint',true));
end

function digest=local_numeric_sha256(value)
value=double(value);
header=unicode2native(sprintf('double|%s|',mat2str(size(value))),'UTF-8');
body=typecast(value(:),'uint8');
payload=[uint8(header(:));body(:)];
digest=local_bytes_sha256(payload);
end

function digest=local_text_sha256(value)
digest=local_bytes_sha256(unicode2native(char(value),'UTF-8'));
end

function digest=local_file_sha256(path)
fid=fopen(path,'rb');if fid<0,error('Cannot hash file: %s.',path);end
cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
digest=local_bytes_sha256(fread(fid,Inf,'*uint8'));
end

function digest=local_bytes_sha256(bytes)
md=java.security.MessageDigest.getInstance('SHA-256');
md.update(typecast(uint8(bytes(:)),'int8'));
digest=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end

function stage2=local_load_stage2(path)
if exist(path,'file')~=2,error('Stage 2 artifact missing: %s',path);end
s=load(path,'a');stage2=s.a;
end
function pattern=local_source_pattern(cfg)
a=linspace(cfg.angle_limits_deg(1),cfg.angle_limits_deg(2),2401).';th=a*pi/180;k=2*pi*cfg.frequency_hz/cfg.c0_mps;d=cos(th).*exp(-.5*(k*cfg.sigma_src_m*sin(th)).^2);d=abs(d)/max(abs(d));pattern=struct('angles_deg',a,'level_db',20*log10(max(d,10^(cfg.source_pattern_clip_db/20))));
end
function result=local_read_pm(c,root)
result.data=read_bellhop_shd_unfolded_vertical([root '.shd']);dm=readmatrix([root '.iwdiag'],'FileType','text','CommentStyle','#');dm=dm(~all(isnan(dm),2),:);names={'alpha_deg','hit_r','hit_z','wall_residual','wall_t_r','wall_t_z','wall_n_r','wall_n_z','tangent_error','normal_error','inc_ur','inc_uz','ref_ur','ref_uz','rot_ur','rot_uz','specular_error','rotation_error','phase_in','phase_ref','phase_delta','amp_in','amp_ref','amp_delta','p1_in','p2_in','p1_ref','p2_ref','p_ref_error','q1_in','q2_in','q1_ref','q2_ref','q_ref_error','p_rot_error','q_rot_error','tau_wall_real','tau_wall_imag','tau_receiver_real','tau_receiver_imag','min_post_dr','n_post','kappa','wall_seg','wall_lambda','wall_tg','wall_th','wall_rm','wall_rn'};if size(dm,2)~=numel(names),error('Unexpected PM diagnostics in %s.',root);end;result.config=c;result.diagnostics=array2table(dm,'VariableNames',names);result.files=struct('data',[root '.shd'],'diagnostics',[root '.iwdiag']);
end
function result=local_read_flat(c,root)
result.data=read_bellhop_shd_unfolded_vertical([root '.shd']);dm=readmatrix([root '.iwdiag'],'FileType','text','CommentStyle','#');dm=dm(~all(isnan(dm),2),:);names={'alpha_deg','hit_r','hit_z','residual','inc_ur','inc_uz','ref_ur','ref_uz','rot_ur','rot_uz','phase_in','phase_ref','phase_delta','amp_in','amp_ref','amp_delta','p1_in','p2_in','p1_ref','p2_ref','p_ref_error','q1_in','q2_in','q1_ref','q2_ref','q_ref_error','p_rot_error','q_rot_error','tau_wall_real','tau_wall_imag','tau_receiver_real','tau_receiver_imag','min_post_dr','n_post','kappa'};if size(dm,2)~=numel(names),error('Unexpected flat diagnostics in %s.',root);end;result.config=c;result.diagnostics=array2table(dm,'VariableNames',names);result.files=struct('data',[root '.shd'],'diagnostics',[root '.iwdiag']);
end
function tf=local_band_match(reference,profile,coeff_file)
tf=numel(profile.k_rad_per_m)==numel(reference.k_rad_per_m) && ...
    max(abs(double(profile.k_rad_per_m(:))-double(reference.k_rad_per_m(:)))) <= 1e-12;
if ~tf,return;end
if exist(coeff_file,'file')~=2
    coeff_file=reference.coeff_file;
end
tbl=readtable(coeff_file);ref=readtable(reference.coeff_file);
if ismember('Sk_m3',tbl.Properties.VariableNames) && ismember('Sk_m3',ref.Properties.VariableNames)
    tf=tf && max(abs(double(tbl.Sk_m3(:))-double(ref.Sk_m3(:)))) <= 1e-12*max(1,max(abs(double(ref.Sk_m3(:)))));
end
end
function stats=local_statistics(rows)
pe=20*log10(abs([rows.G_pe]));bh=20*log10(abs([rows.G_bellhop]));dt=[rows.delta_tl_db];dp=[rows.delta_phase_rad];
q=[5 25 50 75 95];
stats=struct('n',numel(rows),'percentile_levels_pct',q, ...
    'pe_tl_mean_db',mean(pe),'pe_tl_std_db',std(pe,0,2),'pe_tl_percentiles_db',local_percentiles(pe,q), ...
    'bh_tl_mean_db',mean(bh),'bh_tl_std_db',std(bh,0,2),'bh_tl_percentiles_db',local_percentiles(bh,q), ...
    'delta_tl_mean_db',mean(dt),'delta_tl_std_db',std(dt,0,2),'delta_tl_percentiles_db',local_percentiles(dt,q), ...
    'delta_phase_circular_mean_rad',angle(mean(exp(1i*dp))),'delta_phase_circular_std_rad',sqrt(max(0,-2*log(max(abs(mean(exp(1i*dp))),realmin)))),'delta_phase_percentiles_rad',local_percentiles(dp,q), ...
    'pe_power_mean',mean(abs([rows.G_pe]).^2),'bh_power_mean',mean(abs([rows.G_bellhop]).^2));
end
function v=local_percentiles(x,q)
x=sort(x(:));n=numel(x);v=zeros(size(q));
for ii=1:numel(q)
    pos=1+(n-1)*q(ii)/100;lo=floor(pos);hi=ceil(pos);
    if lo==hi,v(ii)=x(lo);else,v(ii)=x(lo)+(pos-lo)*(x(hi)-x(lo));end
end
end
function cfg=local_config(root,o)
cfg=struct('frequency_hz',4000,'c0_mps',1500,'z_tx_m',100,'z_rx_m',3,'sigma_src_m',.3,'xw_m',192.1875,'nx',984,'step_m',.05,'profile_count',4097,'beam_count',5001,'angle_limits_deg',[-15 15],'domain_half_depth_m',1000,'source_pattern_clip_db',-120,'source_geometry','X','mapped_receiver_range_m',103,'selector_tolerance_m',1e-4,'wall_seed',260001,'grazing_mu_threshold',.1,'wall_residual_limit_m',1e-9,'phase_jump_limit_rad',1e-10,'state_limit',1e-12,'minimum_seed_count',8,'seeds',260001:260008,'ensemble_mode','phase','reference_coeff_file',fullfile(root,'results','validation','bellhop_internal_pm_fixed_realization','fixed_pm_fourier_coefficients.csv'),'stage2_mat',fullfile(root,'results','validation','pe_bellhop_pm_frequency_extension','frequency_extension.mat'),'flat_validation_exe',fullfile(root,'results','validation','bellhop_internal_flat_wall_poc','bin','bellhop_iwall_flat_2020.exe'),'pm_validation_exe',fullfile(root,'results','validation','bellhop_internal_pm_wall_poc','bin','bellhop_iwall_pm_2020.exe'),'output_dir',fullfile(root,'results','validation','pe_bellhop_pm_ensemble'),'case_cache_dir','','report_path',fullfile(root,'reports','pe_bellhop_pm_ensemble_comparison_report.md'),'fail_on_check',true);
names=fieldnames(o);for ii=1:numel(names),if ~isfield(cfg,names{ii}),error('Unknown Stage 3 override: %s.',names{ii});end;cfg.(names{ii})=o.(names{ii});end
if isempty(cfg.case_cache_dir),cfg.case_cache_dir=fullfile(cfg.output_dir,'cases');end
cfg.source_geometry=upper(char(cfg.source_geometry));
if ~isscalar(cfg.source_geometry)||~ismember(cfg.source_geometry,['R','X'])
    error('source_geometry must be R (point) or X (line).');
end
end
function row=local_empty_row()
row=struct('seed',NaN,'coeff_file_sha256','','profile_band_match',false,'surface_mean_eta_m',NaN,'surface_rms_eta_m',NaN,'surface_sigma_eta_m',NaN,'surface_min_eta_m',NaN,'surface_max_eta_m',NaN,'profile_rms_height_m',NaN,'profile_rms_slope',NaN,'profile_max_slope',NaN,'profile_rms_curvature_per_m',NaN,'profile_max_curvature_per_m',NaN,'profile_min_radius_m',NaN,'two_k_sigma_eta',NaN,'kmax_over_k',NaN,'G_pe',complex(NaN),'G_bellhop',complex(NaN),'G_pe_abs',NaN,'G_bellhop_abs',NaN,'G_pe_power',NaN,'G_bellhop_power',NaN,'delta_tl_db',NaN,'delta_power_db',NaN,'delta_phase_rad',NaN,'complex_relative_error',NaN,'pe_seam_jump_m',NaN,'pe_outer5_incident_energy_fraction',NaN,'pe_outer5_surface_reflected_energy_fraction',NaN,'pe_outer5_receiver_reflected_energy_fraction',NaN,'wall_residual_max_m',NaN,'min_mu',NaN,'mean_mu',NaN,'mu_p05',NaN,'grazing_fraction',NaN,'wall_hit_count',NaN,'expected_beam_count',NaN,'failed_ray_count',NaN,'rejected_ray_count',NaN,'successful_wall_reflection_fraction',NaN,'hit_slope_mean',NaN,'hit_curvature_mean',NaN,'hit_curvature_rms',NaN,'hit_curvature_p90',NaN,'hit_curvature_p95',NaN,'incidence_angle_mean_deg',NaN,'incidence_angle_std_deg',NaN,'bellhop_case_source','','bellhop_case_fingerprint','','bellhop_source_geometry','','bellhop_beam_count',NaN,'phase_jump_error_rad',NaN,'q_reflect_error',NaN,'p_rotation_error',NaN,'q_rotation_error',NaN,'min_post_range_increment_m',NaN,'all_post_range_positive',false,'tau_receiver_s',NaN);
end
function files=local_write_outputs(a)
out=a.config.output_dir;summary=a.rows;writetable(struct2table(summary),fullfile(out,'ensemble_seed_summary.csv'));s=a.statistics;
scalar_names={'n','pe_tl_mean_db','pe_tl_std_db','bh_tl_mean_db','bh_tl_std_db','delta_tl_mean_db','delta_tl_std_db','delta_phase_circular_mean_rad','delta_phase_circular_std_rad','pe_power_mean','bh_power_mean'};
scalar_struct=rmfield(s,setdiff(fieldnames(s),scalar_names));writetable(struct2table(scalar_struct),fullfile(out,'ensemble_statistics.csv'));
q=s.percentile_levels_pct(:);pt=table(q,s.pe_tl_percentiles_db(:),s.bh_tl_percentiles_db(:),s.delta_tl_percentiles_db(:),s.delta_phase_percentiles_rad(:),'VariableNames',{'percentile_pct','pe_tl_db','bellhop_tl_db','delta_tl_db','delta_phase_rad'});writetable(pt,fullfile(out,'ensemble_percentiles.csv'));
n=fieldnames(a.checks);v=false(size(n));for ii=1:numel(n),v(ii)=a.checks.(n{ii});end;writetable(table(n,v,'VariableNames',{'check_name','passed'}),fullfile(out,'ensemble_checks.csv'));save(fullfile(out,'ensemble_comparison.mat'),'a','-v7.3');local_write_report(a.config.report_path,a);files=struct('summary',fullfile(out,'ensemble_seed_summary.csv'),'statistics',fullfile(out,'ensemble_statistics.csv'),'percentiles',fullfile(out,'ensemble_percentiles.csv'),'checks',fullfile(out,'ensemble_checks.csv'),'mat',fullfile(out,'ensemble_comparison.mat'),'report',a.config.report_path);
end
function local_write_report(path,a)
fid=fopen(path,'w','n','UTF-8');if fid<0,error('Cannot write %s.',path);end;cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
if strcmpi(a.config.ensemble_mode,'amplitude')
    report_title='PE--Bellhop PM coefficient-amplitude ensemble comparison';
    mode_text='The canonical PM spectral-density column is fixed and each non-reference seed independently draws Gaussian cosine/sine coefficients with variance S(k)*Delta-k. No Hs renormalization, recentering, smoothing, tapering, bandwidth change, source fit, or core-physics modification is applied.';
    limit_text='This is a first reduced independent PM coefficient-amplitude ensemble, not a converged ocean Monte Carlo. Bellhop uses the explicitly fingerprinted line-source X convention; no empirical amplitude correction is applied. Geometry, pressure-release phase, p/q state and delay remain hard checks.';
else
    report_title='PE--Bellhop PM paired multi-realization comparison';
    mode_text='Seed 260001 is the canonical fixed realization; the remaining paired profiles preserve the canonical per-mode spectral amplitudes and assign deterministic random phases. No Hs renormalization, recentering, smoothing, tapering, bandwidth change, source fit, or core-physics modification is applied.';
    limit_text='This is a paired fixed-band phase-ensemble smoke/first statistical comparison, not a claim of an independently sampled PM amplitude ensemble. Bellhop uses the explicitly fingerprinted line-source X convention; no empirical amplitude correction is applied. Geometry, pressure-release phase, p/q state and delay remain hard checks.';
end
fprintf(fid,'# %s\n\n状态：**%s**\n\n',report_title,a.classification);fprintf(fid,'固定 4 kHz、uniform c=1500 m/s、W=192.1875 m、nx=984、step=0.05 m、Bellhop source=%s、%d beams、sector ±15°。%s\n\n',a.config.source_geometry,a.config.beam_count,mode_text);
fprintf(fid,'Canonical coefficient source: `%s`; SHA-256: `%s`; requested Kmax `%.12g rad/m`, realized Kmax `%.12g rad/m`.\n\n',a.reference_profile.coeff_file,a.reference_profile.coeff_file_sha256,a.reference_profile.requested_kmax_rad_per_m,a.reference_profile.realized_kmax_rad_per_m);
fprintf(fid,'Flat Bellhop baseline: source `%s`, beams `%d`, fingerprint `%s`. Every rough case must match this beam count; cached outputs are accepted only when their request fingerprint and output hashes match.\n\n',a.flat_case_meta.source,a.flat_case_meta.beam_count,a.flat_case_meta.fingerprint_sha256);
fprintf(fid,'## Per-seed results\n\n| seed | beams | case source | hits/beams | RMS height (m) | RMS slope | max slope | RMS curvature (1/m) | min radius (m) | G_PE | G_Bellhop | delta TL (dB) | delta phase (rad) | complex error | min |u.n| | wall residual (m) | tau (s) |\n|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n');for ii=1:numel(a.rows),r=a.rows(ii);fprintf(fid,'| %d | %d | %s | %d/%d | %.8g | %.8g | %.8g | %.8g | %.8g | %.10g%+.10gi | %.10g%+.10gi | %.8g | %.8g | %.8g | %.8g | %.8g | %.12g |\n',r.seed,r.bellhop_beam_count,r.bellhop_case_source,r.wall_hit_count,r.expected_beam_count,r.profile_rms_height_m,r.profile_rms_slope,r.profile_max_slope,r.profile_rms_curvature_per_m,r.profile_min_radius_m,real(r.G_pe),imag(r.G_pe),real(r.G_bellhop),imag(r.G_bellhop),r.delta_tl_db,r.delta_phase_rad,r.complex_relative_error,r.min_mu,r.wall_residual_max_m,r.tau_receiver_s);end
s=a.statistics;fprintf(fid,'\n## Statistics\n\n| statistic | PE TL (dB) | Bellhop TL (dB) | model delta TL (dB) | model delta phase (rad) |\n|---|---:|---:|---:|---:|\n| mean / circular mean | %.8g | %.8g | %.8g | %.8g |\n| std / circular std | %.8g | %.8g | %.8g | %.8g |\n',s.pe_tl_mean_db,s.bh_tl_mean_db,s.delta_tl_mean_db,s.delta_phase_circular_mean_rad,s.pe_tl_std_db,s.bh_tl_std_db,s.delta_tl_std_db,s.delta_phase_circular_std_rad);
fprintf(fid,'| median / p50 | %.8g | %.8g | %.8g | %.8g |\n',s.pe_tl_percentiles_db(3),s.bh_tl_percentiles_db(3),s.delta_tl_percentiles_db(3),s.delta_phase_percentiles_rad(3));
for ii=1:numel(s.percentile_levels_pct)
    if s.percentile_levels_pct(ii)==50,continue;end
    fprintf(fid,'| p%d | %.8g | %.8g | %.8g | %.8g |\n',s.percentile_levels_pct(ii),s.pe_tl_percentiles_db(ii),s.bh_tl_percentiles_db(ii),s.delta_tl_percentiles_db(ii),s.delta_phase_percentiles_rad(ii));
end
fprintf(fid,'\nMean reflected roughness power: PE %.8g, Bellhop %.8g.\n\n## Checks\n\n',s.pe_power_mean,s.bh_power_mean);n=fieldnames(a.checks);for ii=1:numel(n),fprintf(fid,'- %s: %s\n',n{ii},ternary(a.checks.(n{ii}),'PASS','FAIL'));end;fprintf(fid,'\n## Limits\n\n%s\n',limit_text);
end
function s=ternary(c,a,b),if c,s=a;else,s=b;end,end
function q=local_quantile(x,p)
x=sort(double(x(:)));if isempty(x),q=NaN;return;end;pos=1+(numel(x)-1)*p;lo=floor(pos);hi=ceil(pos);if lo==hi,q=x(lo);else,q=x(lo)+(pos-lo)*(x(hi)-x(lo));end
end
