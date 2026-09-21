function audit = validate_pe_bellhop_pm_constant_height_sign(overrides)
%VALIDATE_PE_BELLHOP_PM_CONSTANT_HEIGHT_SIGN Stage 0D sign audit.
%   Validation-only constant-height surface shifts verify that the PE
%   +2*k*eta phase convention agrees with Bellhop's moved internal wall.

if nargin < 1 || isempty(overrides), overrides = struct(); end
if ~isstruct(overrides) || ~isscalar(overrides)
    error('overrides must be a scalar struct.');
end
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(root);
setup_vertical_project();
cfg = local_config(root,overrides);
if exist(cfg.validation_exe,'file') ~= 2, error('Missing validation executable: %s',cfg.validation_exe); end
if exist(cfg.official_exe,'file') ~= 2, error('Missing official Bellhop executable: %s',cfg.official_exe); end
if ~exist(cfg.output_dir,'dir'), mkdir(cfg.output_dir); end

pattern = local_source_pattern(cfg);
flat = local_run_case(cfg,pattern,0,1,[]);
rows = repmat(local_empty_row(),1,numel(cfg.eta0_m));
for ii = 1:numel(cfg.eta0_m)
    rows(ii) = local_run_case(cfg,pattern,cfg.eta0_m(ii),ii+1,flat);
end

k = 2*pi*cfg.frequency_hz/cfg.c0_mps;
for ii = 1:numel(rows)
    rows(ii).eta0_m = cfg.eta0_m(ii);
    rows(ii).expected_phase_rad = angle(exp(1i*2*k*cfg.eta0_m(ii)));
    rows(ii).pe_ratio_phase_rad = angle(rows(ii).pe_ratio*conj(rows(ii).expected_ratio));
    rows(ii).bh_ratio_phase_rad = angle(rows(ii).bh_ratio*conj(rows(ii).expected_ratio));
    rows(ii).pe_bh_ratio_phase_rad = angle(rows(ii).pe_ratio*conj(rows(ii).bh_ratio));
    rows(ii).pe_ratio_tl_db = 20*log10(abs(rows(ii).pe_ratio));
    rows(ii).bh_ratio_tl_db = 20*log10(abs(rows(ii).bh_ratio));
    rows(ii).expected_ratio_tl_db = 20*log10(abs(rows(ii).expected_ratio));
    rows(ii).wall_expected_phase_rad = angle(rows(ii).wall_pressure*conj(rows(ii).expected_pressure));
    rows(ii).wall_expected_tl_db = 20*log10(abs(rows(ii).wall_pressure)/max(abs(rows(ii).expected_pressure),realmin));
end

tbl = struct2table(rows);
nonzero = abs(tbl.eta0_m)>0;
checks = struct();
checks.pe_sign = all(abs(tbl.pe_ratio_phase_rad(nonzero)) <= cfg.phase_limit_rad);
checks.bellhop_sign = all(abs(tbl.bh_ratio_phase_rad(nonzero)) <= cfg.phase_limit_rad);
checks.cross_model_sign = all(abs(tbl.pe_bh_ratio_phase_rad(nonzero)) <= cfg.phase_limit_rad);
checks.phase_magnitude = checks.pe_sign && checks.bellhop_sign;
checks.pressure_release_once = max(tbl.phase_jump_error_rad) <= cfg.phase_jump_limit_rad;
checks.no_extra_pi = max(abs(tbl.wall_expected_phase_rad)) <= cfg.phase_limit_rad;
checks.path_time = max(abs(tbl.travel_time_error_s)) <= cfg.travel_time_limit_s;
checks.geometry = max(abs(tbl.intersection_residual_m)) <= cfg.intersection_limit_m && max(tbl.kappa_abs_max)==0;
checks.rotation_state = max([tbl.p_rotation_error tbl.q_rotation_error],[],'all') <= cfg.state_limit;
checks.all = all(structfun(@(v) logical(v),checks));

audit = struct('schema_version','1.0.0','stage','0D_constant_height_sign', ...
    'config',cfg,'rows',rows,'table',tbl,'checks',checks,'passed',checks.all);
audit.files = local_write_outputs(audit);
disp(tbl(:,{'eta0_m','expected_phase_rad','pe_ratio_phase_rad','bh_ratio_phase_rad', ...
    'pe_bh_ratio_phase_rad','travel_time_error_s','wall_expected_phase_rad'}));
disp(checks);
if cfg.fail_on_check && ~audit.passed
    error('Stage 0D constant-height sign audit failed; see %s.',cfg.report_path);
end
end

function row = local_run_case(cfg,pattern,eta0,case_index,flat)
% Both legs move when a flat surface is translated: L=103-2*eta0.
L = cfg.reflect_range_m - 2*eta0;
pe_cfg = struct('frequency_hz',cfg.frequency_hz,'c0_mps',cfg.c0_mps, ...
    'xw_m',cfg.xw_m,'nx',cfg.nx,'z_tx_m',cfg.z_tx_m,'z_rx_m',cfg.z_rx_m, ...
    'sigma_src_m',cfg.sigma_src_m,'surface_elevation_x_m',eta0*ones(1,cfg.nx), ...
    'surface_reflect_coeff',-1,'step_m',cfg.step_m,'x_rx_m',0);
pe = run_pe_1d_surface_reflection_validation(pe_cfg);
pair = local_bellhop_cfg(cfg,pattern,L,eta0,case_index);
ref = run_bellhop_unfolded_gaussian_vertical(pair.reference);
wall = run_bellhop_internal_flat_wall_poc_vertical(pair.wall);
p_ref = select_bellhop_shd_pressure_at_range_vertical(ref.data,L,0,cfg.selector_tolerance_m);
p_wall = select_bellhop_shd_pressure_at_range_vertical(wall.data,L,0,cfg.selector_tolerance_m);
d = wall.diagnostics; [~,i0] = min(abs(d.alpha_deg));
row = local_empty_row();
row.eta0_m=eta0; row.pe_pressure=pe.reflected_receiver;
row.reference_pressure=p_ref; row.wall_pressure=p_wall; row.expected_pressure=-p_ref;
row.expected_ratio=exp(1i*2*(2*pi*cfg.frequency_hz/cfg.c0_mps)*eta0);
if isempty(flat)
    row.flat_pe_pressure=pe.reflected_receiver;
    row.flat_bh_pressure=p_wall;
else
    row.flat_pe_pressure=flat.flat_pe_pressure;
    row.flat_bh_pressure=flat.flat_bh_pressure;
end
row.pe_ratio=row.pe_pressure/row.flat_pe_pressure;
row.bh_ratio=row.wall_pressure/row.flat_bh_pressure;
row.intersection_residual_m=max(abs(d.residual)); row.kappa_abs_max=max(abs(d.kappa));
row.phase_jump_error_rad=max(abs(d.phase_delta-pi)); row.amp_jump_error=max(abs(d.amp_delta));
row.p_rotation_error=max(d.p_rot_error); row.q_rotation_error=max(d.q_rot_error);
row.travel_time_s=d.tau_receiver_real(i0); row.travel_time_error_s=row.travel_time_s-L/cfg.c0_mps;
row.min_post_range_increment_m=min(d.min_post_dr); row.all_post_range_positive=all(d.min_post_dr>0);
row.wall_expected_phase_rad=angle(row.wall_pressure*conj(row.expected_pressure));
row.wall_expected_tl_db=20*log10(abs(row.wall_pressure)/max(abs(row.expected_pressure),realmin));
row.case_root=pair.wall.case_root;
end

function pair = local_bellhop_cfg(cfg,pattern,L,eta0,case_index)
if case_index==1, tag='flat'; else, tag=['eta_' local_tag(eta0)]; end
base=fullfile(cfg.output_dir,'cases',tag);
common=struct('run_type','C','frequency_hz',cfg.frequency_hz,'c0_mps',cfg.c0_mps, ...
    'source_depth_m',0,'receiver_depths_m',0,'receiver_ranges_m',[L-1 L], ...
    'beam_count',cfg.beam_count,'angle_limits_deg',cfg.angle_limits_deg, ...
    'domain_half_depth_m',cfg.domain_half_depth_m,'sigma_src_m',cfg.sigma_src_m, ...
    'source_pattern_clip_db',cfg.source_pattern_clip_db,'source_pattern_angles_deg',pattern.angles_deg, ...
    'source_pattern_level_db',pattern.level_db,'step_m',cfg.step_m, ...
    'wall_range_m',cfg.r0_m-eta0,'mapped_receiver_range_m',L);
pair.reference=common; pair.reference.case_root=[base '_reference']; pair.reference.bellhop_exe=cfg.official_exe;
pair.wall=common; pair.wall.case_root=[base '_wall']; pair.wall.bellhop_exe=cfg.validation_exe;
end

function cfg=local_config(root,o)
cfg=struct('frequency_hz',4000,'c0_mps',1500,'z_tx_m',100,'z_rx_m',3, ...
    'sigma_src_m',0.3,'xw_m',192.1875,'nx',984,'step_m',0.05, ...
    'reflect_range_m',103,'r0_m',100,'eta0_m',[0.05 0 -0.05], ...
    'beam_count',10001,'angle_limits_deg',[-30 30],'domain_half_depth_m',1000, ...
    'source_pattern_clip_db',-120, ...
    'validation_exe',fullfile(root,'results','validation','bellhop_internal_flat_wall_poc','bin','bellhop_iwall_flat_2020.exe'), ...
    'official_exe','AcousticsToolbox_2020/windows-bin-20201102/bellhop.exe', ...
    'output_dir',fullfile(root,'results','validation','pe_bellhop_pm_constant_height_sign'), ...
    'report_path',fullfile(root,'reports','pe_bellhop_pm_constant_height_sign_audit.md'), ...
    'phase_limit_rad',0.02,'phase_jump_limit_rad',1e-10,'travel_time_limit_s',1e-10, ...
    'intersection_limit_m',1e-9,'selector_tolerance_m',1e-4,'state_limit',1e-12,'fail_on_check',true);
names=fieldnames(o);
for ii=1:numel(names)
    if ~isfield(cfg,names{ii}), error('Unknown Stage 0D override: %s.',names{ii}); end
    cfg.(names{ii})=o.(names{ii});
end
end

function pattern=local_source_pattern(cfg)
angles=linspace(cfg.angle_limits_deg(1),cfg.angle_limits_deg(2),2401).';
theta=angles*pi/180; k=2*pi*cfg.frequency_hz/cfg.c0_mps;
d=cos(theta).*exp(-0.5*(k*cfg.sigma_src_m*sin(theta)).^2); d=abs(d)/max(abs(d));
pattern=struct('angles_deg',angles,'level_db',20*log10(max(d,10^(cfg.source_pattern_clip_db/20))));
end

function row=local_empty_row()
row=struct('eta0_m',NaN,'expected_phase_rad',NaN,'pe_pressure',complex(NaN), ...
    'reference_pressure',complex(NaN),'wall_pressure',complex(NaN),'expected_pressure',complex(NaN), ...
    'flat_pe_pressure',complex(NaN),'flat_bh_pressure',complex(NaN),'pe_ratio',complex(NaN), ...
    'bh_ratio',complex(NaN),'expected_ratio',complex(NaN),'pe_ratio_phase_rad',NaN, ...
    'bh_ratio_phase_rad',NaN,'pe_bh_ratio_phase_rad',NaN,'pe_ratio_tl_db',NaN,'bh_ratio_tl_db',NaN, ...
    'expected_ratio_tl_db',NaN,'wall_expected_phase_rad',NaN,'wall_expected_tl_db',NaN, ...
    'intersection_residual_m',NaN,'kappa_abs_max',NaN,'phase_jump_error_rad',NaN, ...
    'amp_jump_error',NaN,'p_rotation_error',NaN,'q_rotation_error',NaN,'travel_time_s',NaN, ...
    'travel_time_error_s',NaN,'min_post_range_increment_m',NaN,'all_post_range_positive',false,'case_root','');
end

function files=local_write_outputs(audit)
out=audit.config.output_dir;
writetable(audit.table,fullfile(out,'stage0d_constant_height_rows.csv'));
n=fieldnames(audit.checks); v=false(size(n)); for ii=1:numel(n),v(ii)=audit.checks.(n{ii});end
writetable(table(n,v,'VariableNames',{'check_name','passed'}),fullfile(out,'stage0d_constant_height_checks.csv'));
mat_file=fullfile(out,'stage0d_constant_height_audit.mat'); save(mat_file,'audit','-v7.3');
local_write_report(audit.config.report_path,audit);
files=struct('rows',fullfile(out,'stage0d_constant_height_rows.csv'),'checks',fullfile(out,'stage0d_constant_height_checks.csv'), ...
    'mat',mat_file,'report',audit.config.report_path);
end

function local_write_report(path,audit)
fid=fopen(path,'w','n','UTF-8'); if fid<0,error('Cannot write %s.',path);end
cleanup=onCleanup(@()fclose(fid)); c=audit.config;
fprintf(fid,'# PE--Bellhop PM Stage 0D constant-height sign audit\n\n状态：**%s**\n\n',ternary(audit.passed,'PASS','FAIL'));
fprintf(fid,'本阶段使用 eta0 = +0.05, 0, -0.05 m，保持 Gaussian `.sbp`、uniform c=1500 m/s、pressure-release 和既有 internal-flat POC；不拟合源强或修改 PE/Bellhop 核心。\n\n');
fprintf(fid,'## Results\n\n| eta0 (m) | expected phase | PE ratio phase err | Bellhop ratio phase err | PE--BH phase err | PE TL ratio (dB) | BH TL ratio (dB) | path time err (s) | wall vs -reference phase |\n|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n');
for ii=1:numel(audit.rows)
 r=audit.rows(ii); fprintf(fid,'| %.8g | %.8g | %.8g | %.8g | %.8g | %.8g | %.8g | %.8g | %.8g |\n',r.eta0_m,r.expected_phase_rad,r.pe_ratio_phase_rad,r.bh_ratio_phase_rad,r.pe_bh_ratio_phase_rad,r.pe_ratio_tl_db,r.bh_ratio_tl_db,r.travel_time_error_s,r.wall_expected_phase_rad);
end
fprintf(fid,'\n## Checks\n\n'); n=fieldnames(audit.checks); for ii=1:numel(n),fprintf(fid,'- %s: %s\n',n{ii},ternary(audit.checks.(n{ii}),'PASS','FAIL'));end
fprintf(fid,'\nThe physical image span is L=103-2*eta0 m. The PE screen predicts exp(+i 2 k eta0); Bellhop moves the wall to R0-eta0 and reads the transformed branch at the same L. The known Cartesian backward-range amplitude offset is not used as a gate.\n');
end

function s=ternary(c,a,b),if c,s=a;else,s=b;end,end
function t=local_tag(v),t=strrep(sprintf('%.6g',v),'.','p');t=strrep(t,'-','m');end
