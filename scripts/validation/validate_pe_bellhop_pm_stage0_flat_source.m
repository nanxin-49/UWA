function audit = validate_pe_bellhop_pm_stage0_flat_source(overrides)
%VALIDATE_PE_BELLHOP_PM_STAGE0_FLAT_SOURCE Stage 0C flat source audit.
%   Compares the same flat reflected/direct geometry between the 1-D PE
%   bridge and Bellhop 2020.  No rough profile is used in this stage and no
%   source fitting or absolute-amplitude calibration is performed.

if nargin < 1 || isempty(overrides), overrides = struct(); end
if ~isstruct(overrides) || ~isscalar(overrides)
    error('overrides must be a scalar struct.');
end
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
setup_vertical_project();
cfg = local_config(root, overrides);
if isempty(cfg.official_exe), error('Set BELLHOP_EXE or overrides.official_exe to the AcousticsToolbox 2020 bellhop.exe.'); end
if exist(cfg.validation_exe,'file') ~= 2
    error('Missing Bellhop internal-flat validation executable: %s', cfg.validation_exe);
end
if exist(cfg.official_exe,'file') ~= 2
    error('Missing official Bellhop 2020 executable: %s', cfg.official_exe);
end
if ~exist(cfg.output_dir,'dir'), mkdir(cfg.output_dir); end

x = (-0.5*cfg.xw_m) + (0:cfg.nx-1)*(cfg.xw_m/cfg.nx);
one_cfg = struct('frequency_hz',cfg.frequency_hz,'c0_mps',cfg.c0_mps, ...
    'xw_m',cfg.xw_m,'nx',cfg.nx,'z_tx_m',cfg.z_tx_m,'z_rx_m',cfg.z_rx_m, ...
    'sigma_src_m',cfg.sigma_src_m,'surface_elevation_x_m',zeros(1,cfg.nx), ...
    'surface_reflect_coeff',-1,'step_m',cfg.step_m,'x_rx_m',0);
one = run_pe_1d_surface_reflection_validation(one_cfg);
offsets = cfg.receiver_offsets_m(:).';
ix_offsets = zeros(size(offsets));
for ii=1:numel(offsets), [~,ix_offsets(ii)] = min(abs(x-offsets(ii))); end

pattern = local_source_pattern(cfg);
rows = repmat(local_empty_row(),1,numel(cfg.beam_counts));
for ib=1:numel(cfg.beam_counts)
    beams = cfg.beam_counts(ib);
    tag = sprintf('beams_%d',beams);
    ref_cfg = local_bellhop_config(cfg, pattern, ...
        fullfile(cfg.output_dir,'cases','reference',tag), cfg.official_exe, ...
        [cfg.direct_range_m cfg.reflect_range_m], beams);
    wall_cfg = local_bellhop_config(cfg, pattern, ...
        fullfile(cfg.output_dir,'cases','wall',tag), cfg.validation_exe, ...
        [cfg.reflect_range_m-1 cfg.reflect_range_m], beams);
    ref = run_bellhop_unfolded_gaussian_vertical(ref_cfg);
    wall = run_bellhop_internal_flat_wall_poc_vertical(wall_cfg);
    bh_direct = local_profile(ref.data,cfg.direct_range_m,offsets);
    bh_image = local_profile(ref.data,cfg.reflect_range_m,offsets);
    bh_reflect_reference = -bh_image;
    bh_wall = local_profile(wall.data,cfg.reflect_range_m,offsets);
    q = local_profile_metrics(one, bh_direct, bh_reflect_reference, offsets, ix_offsets, x);
    q.wall_vs_image_complex = max(abs(bh_wall-bh_reflect_reference)) / max(max(abs(bh_reflect_reference)),realmin);
    q.wall_vs_image_phase_rad = max(abs(angle(bh_wall.*conj(bh_reflect_reference))));
    q.wall_vs_image_tl_db = max(abs(20*log10(max(abs(bh_wall),realmin)./max(abs(bh_reflect_reference),realmin))));
    d = wall.diagnostics;
    q.beam_count = beams;
    q.intersection_residual_m = max(abs(d.residual));
    q.phase_jump_error_rad = max(abs(d.phase_delta-pi));
    q.amp_jump_error = max(abs(d.amp_delta));
    q.p_rotation_error = max(d.p_rot_error);
    q.q_rotation_error = max(d.q_rot_error);
    q.kappa_abs_max = max(abs(d.kappa));
    q.travel_time_error_s = d.tau_receiver_real(find(abs(d.alpha_deg)==min(abs(d.alpha_deg)),1)) - cfg.reflect_range_m/cfg.c0_mps;
    q.min_post_range_increment_m = min(d.min_post_dr);
    q.all_post_range_positive = all(d.min_post_dr>0);
    rows(ib)=q;
end

table_rows = struct2table(rows);
checks = struct();
checks.path_and_phase = max(abs(table_rows.wall_vs_image_phase_rad)) <= cfg.phase_limit_rad && ...
    max(abs(table_rows.wall_vs_image_tl_db)) <= cfg.tl_limit_db;
checks.wall_geometry = max(table_rows.intersection_residual_m) <= 1e-9 && ...
    max(table_rows.kappa_abs_max) == 0;
checks.pressure_release = max(table_rows.phase_jump_error_rad) <= 1e-12 && ...
    max(table_rows.amp_jump_error) <= 1e-12;
checks.rotation_state = max(table_rows.p_rotation_error) <= 1e-12 && ...
    max(table_rows.q_rotation_error) <= 1e-12;
checks.positive_range = all(table_rows.all_post_range_positive);
checks.source_axis = max(abs(table_rows.q_axis_phase_rad)) <= cfg.phase_limit_rad && ...
    max(abs(table_rows.q_axis_tl_db)) <= cfg.tl_limit_db;
checks.source_profiles = max(table_rows.direct_profile_magnitude) <= cfg.profile_magnitude_limit && ...
    max(table_rows.reflect_profile_magnitude) <= cfg.profile_magnitude_limit;
checks.beam_convergence = abs(table_rows.q_axis_tl_db(end)-table_rows.q_axis_tl_db(1)) <= cfg.beam_tl_limit_db && ...
    abs(table_rows.q_axis_phase_rad(end)-table_rows.q_axis_phase_rad(1)) <= cfg.beam_phase_limit_rad && ...
    abs(table_rows.direct_profile_complex(end)-table_rows.direct_profile_complex(1)) <= cfg.profile_complex_limit && ...
    abs(table_rows.reflect_profile_complex(end)-table_rows.reflect_profile_complex(1)) <= cfg.profile_complex_limit;
checks.all = all(structfun(@(v) logical(v),checks));

audit = struct('schema_version','1.0.0','stage','0C_flat_source_normalization', ...
    'config',cfg,'one_d',one,'rows',rows,'table',table_rows, ...
    'checks',checks,'passed',checks.all);
audit.files = local_write_outputs(audit);
if cfg.fail_on_check && ~audit.passed
    error('Stage 0C flat source audit failed; see %s.',cfg.report_path);
end
end

function cfg=local_config(root,o)
% The prior internal-wall audit froze a ~0.26 dB backward-range Cartesian
% influence offset. Keep that known diagnostic outside the source/phase hard
% gate with a predeclared 0.30 dB allowance.
cfg=struct('frequency_hz',4000,'c0_mps',1500,'z_tx_m',100,'z_rx_m',3, ...
    'sigma_src_m',0.3,'xw_m',192.1875,'nx',984,'step_m',0.05, ...
    'direct_range_m',97,'reflect_range_m',103,'r0_m',100, ...
    'receiver_offsets_m',[0 1.953125 4.8828125 9.765625 14.6484375 19.53125 ...
        -1.953125 -4.8828125 -9.765625 -14.6484375 -19.53125], ...
    'angle_limits_deg',[-30 30],'domain_half_depth_m',1000, ...
    'beam_counts',[5001 10001],'source_pattern_clip_db',-120, ...
    'validation_exe',fullfile(root,'results','validation','bellhop_internal_flat_wall_poc','bin','bellhop_iwall_flat_2020.exe'), ...
    'official_exe',getenv('BELLHOP_EXE'), ...
    'output_dir',fullfile(root,'results','validation','pe_bellhop_pm_stage0_flat_source'), ...
    'report_path',fullfile(root,'reports','pe_bellhop_pm_stage0_flat_source_report.md'), ...
    'tl_limit_db',0.30,'phase_limit_rad',0.05,'profile_complex_limit',0.05, ...
    'profile_magnitude_limit',0.05, ...
    'beam_tl_limit_db',0.1,'beam_phase_limit_rad',0.02,'fail_on_check',true);
names=fieldnames(o);
for ii=1:numel(names)
    if ~isfield(cfg,names{ii}), error('Unknown Stage 0C override: %s.',names{ii}); end
    cfg.(names{ii})=o.(names{ii});
end
end

function pattern=local_source_pattern(cfg)
angles=linspace(cfg.angle_limits_deg(1),cfg.angle_limits_deg(2),2401).';
theta=angles*pi/180;
k=2*pi*cfg.frequency_hz/cfg.c0_mps;
d=cos(theta).*exp(-0.5*(k*cfg.sigma_src_m*sin(theta)).^2);
d=abs(d)/max(abs(d));
pattern=struct('angles_deg',angles,'level_db',20*log10(max(d,10^(cfg.source_pattern_clip_db/20))));
end

function out=local_bellhop_config(cfg,pattern,case_root,exe,ranges,beams)
out=struct('bellhop_exe',exe,'case_root',case_root,'run_type','C', ...
    'frequency_hz',cfg.frequency_hz,'c0_mps',cfg.c0_mps,'source_depth_m',0, ...
    'receiver_depths_m',cfg.receiver_offsets_m,'receiver_ranges_m',ranges, ...
    'beam_count',beams,'angle_limits_deg',cfg.angle_limits_deg, ...
    'domain_half_depth_m',cfg.domain_half_depth_m,'sigma_src_m',cfg.sigma_src_m, ...
    'source_pattern_clip_db',cfg.source_pattern_clip_db, ...
    'source_pattern_angles_deg',pattern.angles_deg,'source_pattern_level_db',pattern.level_db, ...
    'step_m',cfg.step_m,'wall_range_m',cfg.r0_m,'mapped_receiver_range_m',cfg.reflect_range_m);
end

function p=local_profile(shade,target_range,depths)
p=zeros(1,numel(depths));
for ii=1:numel(depths)
    p(ii)=select_bellhop_shd_pressure_at_range_vertical(shade,target_range,depths(ii));
end
end

function row=local_profile_metrics(one,bh_direct,bh_reflect,offsets,ix_offsets,x)
one_direct=one.direct_field(ix_offsets);
one_reflect=one.reflected_field(ix_offsets);
axis_index=find(abs(offsets)<eps,1);
row=local_empty_row();
one_direct_n=one_direct/one_direct(axis_index);
bh_direct_n=bh_direct/bh_direct(axis_index);
one_reflect_n=one_reflect/one_reflect(axis_index);
bh_reflect_n=bh_reflect/bh_reflect(axis_index);
% The PE bridge and Bellhop use opposite transverse phasor orientation in
% this unfolded chart.  Keep the complex discrepancy diagnostic, but gate
% the source mapping on the invariant normalized magnitude profile.
row.direct_profile_complex=max(abs(one_direct_n-bh_direct_n));
row.reflect_profile_complex=max(abs(one_reflect_n-bh_reflect_n));
row.direct_profile_magnitude=max(abs(abs(one_direct_n)-abs(bh_direct_n)));
row.reflect_profile_magnitude=max(abs(abs(one_reflect_n)-abs(bh_reflect_n)));
q1=one.reflected_receiver/one.direct_receiver;
q2=bh_reflect(axis_index)/bh_direct(axis_index);
row.q_axis_tl_db=20*log10(abs(q1)/max(abs(q2),realmin));
row.q_axis_phase_rad=angle(q1*conj(q2));
row.q_axis_complex=abs(q1-q2)/max(abs(q2),realmin);
row.direct_axis_scale_ratio=one.direct_receiver/bh_direct(axis_index);
row.reflect_axis_scale_ratio=one.reflected_receiver/bh_reflect(axis_index);
row.receiver_offset_grid_error_m=max(abs(x(ix_offsets)-offsets));
end

function row=local_empty_row()
row=struct('beam_count',NaN,'direct_profile_complex',NaN,'reflect_profile_complex',NaN, ...
    'direct_profile_magnitude',NaN,'reflect_profile_magnitude',NaN, ...
    'q_axis_tl_db',NaN,'q_axis_phase_rad',NaN,'q_axis_complex',NaN, ...
    'direct_axis_scale_ratio',complex(NaN),'reflect_axis_scale_ratio',complex(NaN), ...
    'receiver_offset_grid_error_m',NaN,'wall_vs_image_complex',NaN, ...
    'wall_vs_image_phase_rad',NaN,'wall_vs_image_tl_db',NaN, ...
    'intersection_residual_m',NaN,'phase_jump_error_rad',NaN,'amp_jump_error',NaN, ...
    'p_rotation_error',NaN,'q_rotation_error',NaN,'kappa_abs_max',NaN, ...
    'travel_time_error_s',NaN,'min_post_range_increment_m',NaN, ...
    'all_post_range_positive',false);
end

function files=local_write_outputs(audit)
out=audit.config.output_dir;
writetable(audit.table,fullfile(out,'stage0c_flat_source_rows.csv'));
check_names=fieldnames(audit.checks);
check_values=false(size(check_names));
for ii=1:numel(check_names), check_values(ii)=audit.checks.(check_names{ii}); end
check_table=table(check_names,check_values,'VariableNames',{'check_name','passed'});
writetable(check_table,fullfile(out,'stage0c_flat_source_checks.csv'));
mat_file=fullfile(out,'stage0c_flat_source_audit.mat');
save(mat_file,'audit','-v7.3');
local_write_report(audit.config.report_path,audit);
files=struct('rows',fullfile(out,'stage0c_flat_source_rows.csv'), ...
    'checks',fullfile(out,'stage0c_flat_source_checks.csv'),'mat',mat_file, ...
    'report',audit.config.report_path);
end

function local_write_report(path,audit)
fid=fopen(path,'w','n','UTF-8'); if fid<0,error('Cannot write %s.',path);end
cleanup=onCleanup(@()fclose(fid)); c=audit.config;
fprintf(fid,'# PE--Bellhop PM Stage 0C flat source/normalization audit\n\n');
fprintf(fid,'状态：**%s**\n\n',ternary(audit.passed,'PASS','FAIL'));
fprintf(fid,'本阶段只运行 flat reflected/direct geometry；不使用随机 PM、不拟合 SBP、 不比较绝对源强。PE 是 Stage 0B 的一横向维 bridge，Bellhop 使用当前 2020 internal-flat validation binary。\n\n');
fprintf(fid,'## 配置\n\n');
fprintf(fid,'- f/c：%.0f Hz / %.0f m/s；Tx/Rx：%.3g / %.3g m；direct/image：%.3g / %.3g m\n',c.frequency_hz,c.c0_mps,c.z_tx_m,c.z_rx_m,c.direct_range_m,c.reflect_range_m);
fprintf(fid,'- PE window/grid：%.9g m / %d；step：%.4g m；sponge off\n',c.xw_m,c.nx,c.step_m);
fprintf(fid,'- Bellhop beams：[%s]；step：%.4g m；SBP：2401 points; receiver offsets：%d\n',sprintf('%d ',c.beam_counts),c.step_m,numel(c.receiver_offsets_m));
fprintf(fid,'\n## Results\n\n| beams | axis Q TL dB | axis Q phase rad | axis Q complex | direct profile mag | reflect profile mag | direct profile complex (diag) | reflect profile complex (diag) | wall/image complex |\n|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n');
for ii=1:numel(audit.rows)
    r=audit.rows(ii);
    fprintf(fid,'| %d | %.8g | %.8g | %.8g | %.8g | %.8g | %.8g | %.8g | %.8g |\n',r.beam_count,r.q_axis_tl_db,r.q_axis_phase_rad,r.q_axis_complex,r.direct_profile_magnitude,r.reflect_profile_magnitude,r.direct_profile_complex,r.reflect_profile_complex,r.wall_vs_image_complex);
end
fprintf(fid,'\n## Native wall invariants\n\n');
r=audit.rows(end);
fprintf(fid,'intersection residual max: %.8g m; kappa max: %.8g 1/m; pressure-release phase error: %.8g rad; Amp jump: %.8g; p/q rotation errors: %.8g / %.8g; travel-time error: %.8g s; minimum transformed range increment: %.8g m.\n\n',r.intersection_residual_m,r.kappa_abs_max,r.phase_jump_error_rad,r.amp_jump_error,r.p_rotation_error,r.q_rotation_error,r.travel_time_error_s,r.min_post_range_increment_m);
fprintf(fid,'## Checks\n\n');
for ii=1:numel(fieldnames(audit.checks))
    n=fieldnames(audit.checks); name=n{ii};
    fprintf(fid,'- %s: %s\n',name,ternary(audit.checks.(name),'PASS','FAIL'));
end
fprintf(fid,'\n## Interpretation\n\n');
fprintf(fid,'The flat axis Q and normalized offset magnitude profiles are the source/directivity audit for the later PM ratio. The complex offset-profile residual is retained as a diagnostic because the PE bridge and Bellhop unfolded chart use opposite transverse phasor orientation; no conjugation, source fitting, or calibration is applied. The reported absolute axis scale ratios are diagnostic only. If this stage passes, Stage 0D constant-height sign audit is the only next step.\n');
end

function out=ternary(condition,yes_value,no_value)
if condition,out=yes_value;else,out=no_value;end
end
