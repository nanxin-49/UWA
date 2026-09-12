function result = validate_pe_bellhop_controlled_comparison(stage, overrides)
%VALIDATE_PE_BELLHOP_CONTROLLED_COMPARISON
% Sequential, validation-only driver for the frozen PE--Bellhop Goal.
%
% P0 and Stages 0--4 are validation-only entry points. Stage 2 is unlocked
% only by the preserved Stage-1 PASS and Stage-1Y theoretical convention
% closure; Stage 3 is unlocked only by a completed Stage-2 height sweep.

if nargin < 1 || isempty(stage), stage = 'p0'; end
if nargin < 2 || isempty(overrides), overrides = struct(); end
if ~(ischar(stage) || isstring(stage)) || (isstring(stage) && ~isscalar(stage)) || isempty(stage)
    error('stage must be a scalar text value.');
end
if ~isstruct(overrides) || ~isscalar(overrides)
    error('overrides must be a scalar struct.');
end

% Stage 0 helpers are subfunctions below.

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
setup_vertical_project();
cfg = local_config(root, overrides);
if ~exist(cfg.output_dir, 'dir'), mkdir(cfg.output_dir); end

switch lower(char(stage))
    case 'p0'
        result = local_p0(cfg, root);
    case 'stage0'
        result = local_stage0(cfg, root);
    case 'stage1'
        result = local_stage1(cfg, root);
    case 'stage2'
        result = local_stage2(cfg, root);
    case 'stage3'
        result = local_stage3(cfg, root);
    case 'stage4'
        result = local_stage4(cfg, root);
    otherwise
        error('Stage %s is not unlocked or implemented.', char(stage));
end
end

function cfg = local_config(root, overrides)
% User execution override (2026-09-10): subsequent runs are limited to
% 10,001 beams; do not launch the 20,001-beam endpoint.
cfg = struct( ...
    'official_exe', 'E:\MISC\BELLHOP\AcousticsToolbox_2020\windows-bin-20201102\bellhop.exe', ...
    'flat_validation_exe', fullfile(root, 'results', 'validation', ...
        'bellhop_internal_flat_wall_poc', 'bin', 'bellhop_iwall_flat_2020.exe'), ...
    'parametric_validation_exe', fullfile(root, 'results', 'validation', ...
        'bellhop_internal_pm_wall_poc', 'bin', 'bellhop_iwall_pm_2020.exe'), ...
    'flat_overlay_bellhop', fullfile(root, 'scripts', 'validation', 'support', ...
        'bellhop_internal_flat_wall_poc', 'bellhop.f90'), ...
    'flat_overlay_step', fullfile(root, 'scripts', 'validation', 'support', ...
        'bellhop_internal_flat_wall_poc', 'Step.f90'), ...
    'parametric_overlay_bellhop', fullfile(root, 'scripts', 'validation', 'support', ...
        'bellhop_internal_pm_wall_poc', 'bellhop.f90'), ...
    'parametric_overlay_step', fullfile(root, 'scripts', 'validation', 'support', ...
        'bellhop_internal_pm_wall_poc', 'Step.f90'), ...
    'output_dir', fullfile(root, 'results', 'validation', ...
        'pe_bellhop_controlled_comparison'), ...
    'frequency_hz', 4000, 'c0_mps', 1500, 'z_tx_m', 100, 'z_rx_m', 3, ...
    'sigma_src_m', 0.3, 'xw_m', 192.1875, 'nx', 984, 'pe_step_m', 0.05, ...
    'bellhop_step_m', 0.05, 'beam_counts', [5001 10001], ...
    'angle_limits_deg', [-30 30], 'source_geometry', 'X', 'run_type', 'C', ...
    'wall_r0_m', 100, 'wall_support_m', [-80 80], 'wall_profile_count', 4097, ...
    'mapped_receiver_range_m', 103, 'guard_receiver_range_m', 102, ...
    'receiver_tolerance_m', 1e-6, 'domain_half_depth_m', 1000, ...
    'source_pattern_samples', 2401, 'source_pattern_clip_db', -120, ...
    'phase_sign', -1, 'phase_floor_db', -40, 'wall_seed', 0, ...
    'stage0_flat_l2_limit', 0.02, ...
    'stage0_flat_phase_limit_rad', 0.02, 'stage0_flat_tl_limit_db', 0.10, ...
    'stage0_flat_phase_p95_limit_rad', 0.05, 'stage0_flat_rho_limit', 0.9995, ...
    'stage0_flat_global_phase_limit_rad', 0.02, 'stage0_flat_aligned_l2_limit', 0.02, ...
    'stage0_chain_l2_limit', 0.002, 'stage0_chain_phase_limit_rad', 0.005, ...
    'stage0_chain_tl_limit_db', 0.02, 'stage0_beam_l2_limit', 0.002, ...
    'stage0_beam_phase_limit_rad', 0.005, 'stage0_beam_tl_limit_db', 0.02, ...
    'stage0_halfdx_l2_limit', 0.005, 'fail_on_check', true, ...
    'stage1_amplitude_m', 0.01, 'stage1_wavenumber_radpm', 0.10, ...
    'stage1_beam_counts', [10001 10001], 'stage1_profile_counts', [2049 4097], ...
    'stage1_output_subdir', 'stage1_convention_fixed', ...
    'stage2_amplitudes_m', [0.01 0.02 0.05 0.10 0.20], ...
    'stage2_wavenumber_radpm', 0.10, ...
    'stage2_output_subdir', 'stage2_height_sweep', ...
    'stage3_wavenumbers_radpm', [0.10 0.20 0.35 0.47], ...
    'stage3_run_wavenumbers_radpm', [], ...
    'stage3_output_subdir', 'stage3_slope_curvature_sweep', ...
    'stage4_output_subdir', 'stage4_validity_map', ...
    'canonical_coeff_source', fullfile('E:', filesep, 'MISC', 'CARPE3D_matlab', ...
        'Explain', 'results', 'validation', 'bellhop_internal_pm_fixed_realization', ...
        'fixed_pm_fourier_coefficients.csv'), ...
    'canonical_coeff_sha256', ...
        '1f7eda465e4ae85b8ac038310edf053f2d062563a3b943bfd83bd59d07015f67');

names = fieldnames(overrides);
for ii = 1:numel(names)
    if ~isfield(cfg, names{ii}), error('Unknown P0 override: %s.', names{ii}); end
    cfg.(names{ii}) = overrides.(names{ii});
end
end

function result = local_p0(cfg, root)
cfg.source_geometry = upper(char(cfg.source_geometry));
cfg.run_type = upper(char(cfg.run_type));
if ~strcmp(cfg.source_geometry, 'X'), error('P0 requires source geometry X.'); end
if ~strcmp(cfg.run_type, 'C'), error('P0 requires coherent run type C.'); end
if cfg.frequency_hz ~= 4000 || cfg.c0_mps ~= 1500
    error('P0 is frozen to 4 kHz and c=1500 m/s.');
end
if mod(cfg.nx, 2) ~= 0 || cfg.nx < 4, error('nx must be even and >=4.'); end

x = (-0.5 * cfg.xw_m) + (0:cfg.nx-1) * (cfg.xw_m / cfg.nx);
z_unsorted = -x;
[z_sorted, order] = sort(z_unsorted, 'ascend');
inverse_order = zeros(size(order)); inverse_order(order) = 1:numel(order);
receiver = struct('x_pe_m', x, 'z_bh_unsorted_m', z_unsorted, ...
    'z_bh_sorted_m', z_sorted, 'sort_order', order, ...
    'inverse_order', inverse_order, ...
    'map_formula', 'z_BH_prime = -x_PE after T(r,z)=(2R0-r,-z)');

pattern = local_pattern(cfg);
files = {cfg.official_exe, cfg.flat_validation_exe, cfg.parametric_validation_exe, ...
    cfg.flat_overlay_bellhop, cfg.flat_overlay_step, cfg.parametric_overlay_bellhop, ...
    cfg.parametric_overlay_step};
labels = {'official_exe','flat_validation_exe','parametric_validation_exe', ...
    'flat_overlay_bellhop','flat_overlay_step','parametric_overlay_bellhop', ...
    'parametric_overlay_step'};
checks = struct();
for ii = 1:numel(files)
    checks.(labels{ii}) = exist(files{ii}, 'file') == 2;
end
checks.receiver_map_finite = all(isfinite([x z_sorted])) && ...
    numel(unique(z_sorted)) == cfg.nx;
checks.receiver_map_within_domain = max(abs(z_sorted)) <= cfg.domain_half_depth_m;
checks.profile_count_valid = cfg.wall_profile_count >= 3 && ...
    cfg.wall_support_m(2) > cfg.wall_support_m(1);
checks.source_x = strcmp(cfg.source_geometry, 'X');
checks.run_c = strcmp(cfg.run_type, 'C');
checks.official_2020_path = contains(lower(cfg.official_exe), 'acousticstoolbox_2020');
checks.all = all(structfun(@(v) logical(v), checks));

hashes = struct();
for ii = 1:numel(files)
    if checks.(labels{ii})
        hashes.(labels{ii}) = local_sha256_file(files{ii});
    else
        hashes.(labels{ii}) = '';
    end
end
manifest = struct('schema_version','1.0.0','stage','P0', ...
    'project_root',root,'config',cfg,'source_pattern',pattern, ...
    'receiver_mapping',receiver,'file_sha256',hashes, ...
    'checks',checks,'passed',checks.all, ...
    'generated_without_solver_calls',true);

manifest_file = fullfile(cfg.output_dir, 'p0_manifest.mat');
save(manifest_file, 'manifest', '-v7');
json_file = fullfile(cfg.output_dir, 'p0_manifest.json');
fid = fopen(json_file, 'w', 'n', 'UTF-8');
if fid < 0, error('Cannot write %s.', json_file); end
cleanup = onCleanup(@() fclose(fid));
fwrite(fid, jsonencode(manifest), 'char');
clear cleanup

names = string(fieldnames(checks)); values = false(size(names));
for ii = 1:numel(names), values(ii) = checks.(char(names(ii))); end
writetable(table(names, values, 'VariableNames', {'check_name','passed'}), ...
    fullfile(cfg.output_dir, 'p0_checks.csv'));
report_file = fullfile(cfg.output_dir, 'p0_report.md');
local_write_report(report_file, manifest, pattern, hashes);

result = struct('manifest',manifest,'files',struct('mat',manifest_file, ...
    'json',json_file,'checks',fullfile(cfg.output_dir,'p0_checks.csv'), ...
    'report',report_file),'passed',manifest.passed);
if ~result.passed, error('P0 preparation failed; see %s.', report_file); end
end

function result = local_stage0(cfg, root)
% Stage 0: flat reflected-line closure.  All fields are axis-normalized;
% the AS field defines the 95/99-percent energy footprints.
out_dir = fullfile(cfg.output_dir,'stage0');
if ~exist(out_dir,'dir'), mkdir(out_dir); end
pat = local_pattern(cfg);
x = (-0.5*cfg.xw_m) + (0:cfg.nx-1)*(cfg.xw_m/cfg.nx);
pe_cfg = struct('frequency_hz',cfg.frequency_hz,'c0_mps',cfg.c0_mps, ...
    'xw_m',cfg.xw_m,'nx',cfg.nx,'z_tx_m',cfg.z_tx_m,'z_rx_m',cfg.z_rx_m, ...
    'sigma_src_m',cfg.sigma_src_m,'surface_elevation_x_m',zeros(1,cfg.nx), ...
    'surface_reflect_coeff',-1,'step_m',cfg.pe_step_m,'x_rx_m',0);
pe = run_pe_1d_surface_reflection_validation(pe_cfg);
as = local_stage0_as(cfg,x);
fp = local_stage0_footprint(x,as.field,cfg);
pe_as = local_stage0_metrics(pe.reflected_field,as.field,fp,pe.ix_rx);
bh_low = local_stage0_internal(cfg,x,cfg.beam_counts(1),'internal_flat_b5001');
bh_high = local_stage0_internal(cfg,x,cfg.beam_counts(2),'internal_flat_b10001');
bh_low.field = local_stage0_convert(bh_low.pressure_raw,cfg.phase_sign);
bh_high.field = local_stage0_convert(bh_high.pressure_raw,cfg.phase_sign);
bh_free = local_stage0_free(cfg,x,cfg.beam_counts(2),'official_free_b10001');
bh_free.field = -local_stage0_convert(bh_free.pressure_raw,cfg.phase_sign);
bh_ded = local_stage0_dedicated_flat(cfg,x,cfg.beam_counts(2),'dedicated_flat_b10001');
bh_ded.field = local_stage0_convert(bh_ded.pressure_raw,cfg.phase_sign);
metrics = struct('pe_as',pe_as, ...
    'flat_internal',local_stage0_metrics(pe.reflected_field,bh_high.field,fp,pe.ix_rx), ...
    'chain_internal_free',local_stage0_metrics(bh_high.field,bh_free.field,fp,pe.ix_rx), ...
    'dedicated_vs_internal',local_stage0_metrics(bh_high.field,bh_ded.field,fp,pe.ix_rx), ...
    'beam_convergence',local_stage0_metrics(bh_low.field,bh_high.field,fp,pe.ix_rx));
xh = (-0.5*cfg.xw_m)+(0:(2*cfg.nx)-1)*(cfg.xw_m/(2*cfg.nx));
bh_half = local_stage0_internal(cfg,xh,cfg.beam_counts(2),'internal_flat_b10001_half_dx');
bh_half.field = local_stage0_convert(bh_half.pressure_raw,cfg.phase_sign);
[mi,hi,err] = local_stage0_shared_grid(x,xh,cfg.receiver_tolerance_m);
sampling = struct('executed',true,'metrics',local_stage0_metrics(bh_high.field(mi),bh_half.field(hi), ...
    local_stage0_footprint_subset(fp,mi),find(mi==pe.ix_rx,1)), ...
    'shared_count',numel(mi),'max_coordinate_error_m',max(err));
finite_count = sum(~isfinite([pe.reflected_field(:);as.field(:);bh_low.field(:);bh_high.field(:);bh_free.field(:);bh_ded.field(:);bh_half.field(:)]));
checks = struct();
checks.pe_as_complex = metrics.pe_as.max_complex_full <= 1e-10;
checks.pe_as_l2 = metrics.pe_as.l2_m99 <= 1e-10;
checks.pe_outer5 = local_stage0_outer5(pe.reflected_field) <= 1e-5;
checks.receiver_coordinates = max([bh_low.receiver_coordinate_error_m bh_high.receiver_coordinate_error_m bh_free.receiver_coordinate_error_m bh_ded.receiver_coordinate_error_m bh_half.receiver_coordinate_error_m]) <= cfg.receiver_tolerance_m;
checks.flat_l2 = metrics.flat_internal.l2_m99 <= cfg.stage0_flat_l2_limit;
checks.flat_phase = metrics.flat_internal.phase_rms_m99 <= cfg.stage0_flat_phase_limit_rad;
checks.flat_phase_p95 = metrics.flat_internal.phase_p95_m95 <= cfg.stage0_flat_phase_p95_limit_rad;
checks.flat_tl_p95 = metrics.flat_internal.tl_p95_m95 <= cfg.stage0_flat_tl_limit_db;
checks.flat_shape = metrics.flat_internal.rho_shape >= cfg.stage0_flat_rho_limit;
checks.flat_global_phase = abs(metrics.flat_internal.global_phase_rad) <= cfg.stage0_flat_global_phase_limit_rad;
checks.flat_aligned = metrics.flat_internal.aligned_l2_m99 <= cfg.stage0_flat_aligned_l2_limit;
checks.chain_l2 = metrics.chain_internal_free.l2_m99 <= cfg.stage0_chain_l2_limit;
checks.chain_phase = metrics.chain_internal_free.phase_rms_m99 <= cfg.stage0_chain_phase_limit_rad;
checks.chain_tl = metrics.chain_internal_free.tl_rms_m99 <= cfg.stage0_chain_tl_limit_db;
checks.beam_l2 = metrics.beam_convergence.l2_m99 <= cfg.stage0_beam_l2_limit;
checks.beam_phase = metrics.beam_convergence.phase_rms_m99 <= cfg.stage0_beam_phase_limit_rad;
checks.beam_tl = metrics.beam_convergence.tl_rms_m99 <= cfg.stage0_beam_tl_limit_db;
checks.half_dx = sampling.metrics.l2_m99 <= cfg.stage0_halfdx_l2_limit;
checks.geometry = local_stage0_geometry_check(bh_high,cfg);
checks.finite = finite_count==0;
checks.all = all(structfun(@(v)logical(v),checks));
validation = struct('schema_version','1.0.0','stage','stage0_flat_reflected_line', ...
    'config',cfg,'source_pattern',pat,'x_m',x,'pe',pe,'as',as,'footprint',fp, ...
    'bellhop_internal_low',bh_low,'bellhop_internal_high',bh_high, ...
    'bellhop_official_free',bh_free,'bellhop_dedicated_flat',bh_ded, ...
    'receiver_sampling',sampling,'metrics',metrics,'checks',checks, ...
    'finite_count',finite_count,'passed',checks.all);
mat_file = fullfile(out_dir,'stage0_validation.mat'); save(mat_file,'validation','-v7');
local_stage0_write_csv(fullfile(out_dir,'stage0_checks.csv'),checks);
report_file = fullfile(out_dir,'stage0_report.md'); local_stage0_write_report(report_file,validation);
validation.files = struct('mat',mat_file,'checks',fullfile(out_dir,'stage0_checks.csv'),'report',report_file);
result=validation;
if isfield(cfg,'fail_on_check') && cfg.fail_on_check && ~result.passed
    error('Stage 0 failed; see %s.',report_file);
end
end

function result = local_stage1(cfg, root)
% Stage 1: weak, smooth sinusoid.  The same flat denominator and AS-defined
% masks from Stage 0 are used for both solvers.
stage0_file=fullfile(cfg.output_dir,'stage0','stage0_validation.mat');
if exist(stage0_file,'file')~=2, error('Stage 0 result is required before Stage 1.'); end
s0=load(stage0_file,'validation'); s0=s0.validation;
if ~s0.passed, error('Stage 0 did not pass; Stage 1 remains locked.'); end
A=cfg.stage1_amplitude_m; K=cfg.stage1_wavenumber_radpm; x=s0.x_m; fp=s0.footprint;
pe_cfg=struct('frequency_hz',cfg.frequency_hz,'c0_mps',cfg.c0_mps,'xw_m',cfg.xw_m,'nx',cfg.nx, ...
    'z_tx_m',cfg.z_tx_m,'z_rx_m',cfg.z_rx_m,'sigma_src_m',cfg.sigma_src_m, ...
    'surface_elevation_x_m',A*cos(K*x),'surface_reflect_coeff',-1, ...
    'step_m',cfg.pe_step_m,'x_rx_m',0);
pe=run_pe_1d_surface_reflection_validation(pe_cfg);
flat_pe=s0.pe.reflected_field; flat_bh=s0.bellhop_internal_high.field;
N=cfg.stage1_profile_counts; bhN=cell(1,numel(N)); GbhN=cell(1,numel(N));
try
    for ii=1:numel(N)
        bhN{ii}=local_stage1_internal(cfg,x,cfg.stage1_beam_counts(1),local_stage1_tag(A,K,N(ii)),A,K,N(ii));
        bhN{ii}.field=local_stage0_convert(bhN{ii}.pressure_raw,cfg.phase_sign);
        GbhN{ii}=zeros(size(flat_bh)); GbhN{ii}(fp.m99.mask)=bhN{ii}.field(fp.m99.mask)./flat_bh(fp.m99.mask);
    end
    bh_low=bhN{end}; bh_low.field=GbhN{end};
    if cfg.stage1_beam_counts(2)==cfg.stage1_beam_counts(1)
        % The user-requested 10,001-beam policy makes the endpoint identical
        % to the profile-convergence beam case; reuse it instead of rerunning
        % an identical Bellhop job.
        bh_hi=bhN{end};
    else
        bh_hi=local_stage1_internal(cfg,x,cfg.stage1_beam_counts(2),local_stage1_tag(A,K,N(end)),A,K,N(end));
        bh_hi.field=local_stage0_convert(bh_hi.pressure_raw,cfg.phase_sign);
    end
catch ME
    out_dir=fullfile(cfg.output_dir,cfg.stage1_output_subdir); if ~exist(out_dir,'dir'),mkdir(out_dir);end
    validation=struct('schema_version','1.0.0','stage','stage1_weak_sinusoid', ...
        'config',cfg,'A_m',A,'K_radpm',K,'passed',false,'blocked',true, ...
        'blocked_reason',ME.message,'error_identifier',ME.identifier, ...
        'profile_counts',N,'beam_counts',cfg.stage1_beam_counts);
    mat_file=fullfile(out_dir,'stage1_validation.mat'); save(mat_file,'validation','-v7');
    report_file=fullfile(out_dir,'stage1_report.md'); local_stage1_blocked_report(report_file,validation);
    validation.files=struct('mat',mat_file,'report',report_file); result=validation; return
end
Gpe=zeros(size(flat_pe)); Gbh=zeros(size(flat_bh));
valid_ratio=fp.m99.mask;
Gpe(valid_ratio)=pe.reflected_field(valid_ratio)./flat_pe(valid_ratio);
Gbh(valid_ratio)=bh_hi.field(valid_ratio)./flat_bh(valid_ratio);
profile=local_stage1_ratio_metrics(GbhN{1},GbhN{2},fp);
beam=local_stage1_ratio_metrics(GbhN{2},Gbh,fp);
% Stage 1X provenance audit identified a fixed Bellhop ratio conjugation:
% the solver field convention remains unchanged, while the rough/flat
% comparison uses G_BH_cmp = conj(G_BH).  Keep the unmodified metric as a
% diagnostic and never feed this transform into PE/Bellhop physics.
Gbh_cmp=conj(Gbh);
model_raw=local_stage1_ratio_metrics(Gpe,Gbh,fp);
model=local_stage1_ratio_metrics(Gpe,Gbh_cmp,fp);
geom=local_stage1_geometry(bh_hi,cfg);
floorE=max([s0.metrics.flat_internal.l2_m99,s0.metrics.chain_internal_free.l2_m99,s0.metrics.beam_convergence.l2_m99,s0.receiver_sampling.metrics.l2_m99]);
floorPhi=max([s0.metrics.flat_internal.phase_rms_m99,s0.metrics.beam_convergence.phase_rms_m99,s0.receiver_sampling.metrics.phase_rms_m99]);
floorTL=max([s0.metrics.flat_internal.tl_rms_m99,s0.metrics.beam_convergence.tl_rms_m99,s0.receiver_sampling.metrics.tl_rms_m99]);
T_E=floorE+max(.005,.5*floorE); T_phi=floorPhi+max(.005,.5*floorPhi); T_TL=floorTL+max(.01,.5*floorTL);
checks=struct('profile_l2',profile.l2_m99<=cfg.stage0_beam_l2_limit,'profile_phase',profile.phase_rms_m99<=cfg.stage0_beam_phase_limit_rad, ...
 'profile_tl',profile.tl_rms_m99<=cfg.stage0_beam_tl_limit_db,'beam_l2',beam.l2_m99<=cfg.stage0_beam_l2_limit, ...
 'beam_phase',beam.phase_rms_m99<=cfg.stage0_beam_phase_limit_rad,'beam_tl',beam.tl_rms_m99<=cfg.stage0_beam_tl_limit_db, ...
 'geometry',geom.all,'weak_E',model.l2_m99<=T_E,'weak_phase',model.phase_rms_m99<=T_phi, ...
 'weak_tl',model.tl_rms_m99<=T_TL,'weak_shape',model.rho_shape>=.9995,'finite',all(isfinite([Gpe(:);Gbh(:)])));
checks.all=all(structfun(@(v)logical(v),checks));
validation=struct('schema_version','1.0.0','stage','stage1_weak_sinusoid','config',cfg,'A_m',A,'K_radpm',K, ...
 'pe',pe,'flat_pe',flat_pe,'flat_bellhop',flat_bh,'bellhop_profile_cases',{bhN}, ...
 'bellhop_high',bh_hi,'G_PE',Gpe,'G_BH',Gbh,'G_BH_comparison',Gbh_cmp, ...
 'comparison_convention','G_BH_comparison = conj(G_BH); fixed by Stage 1X provenance audit', ...
 'metrics',struct('profile',profile,'beam',beam,'model_raw',model_raw,'model',model), ...
 'geometry',geom,'floor',struct('F_E',floorE,'F_phi',floorPhi,'F_TL',floorTL,'T_E',T_E,'T_phi',T_phi,'T_TL',T_TL), ...
 'checks',checks,'passed',checks.all);
out_dir=fullfile(cfg.output_dir,cfg.stage1_output_subdir); if ~exist(out_dir,'dir'),mkdir(out_dir);end
mat_file=fullfile(out_dir,'stage1_validation.mat'); save(mat_file,'validation','-v7');
report_file=fullfile(out_dir,'stage1_report.md'); local_stage1_write_report(report_file,validation);
local_stage1_append_convention_report(report_file,validation);
validation.files=struct('mat',mat_file,'report',report_file); result=validation;
if cfg.fail_on_check && ~result.passed, error('Stage 1 failed; see %s.',report_file); end
end

function run=local_stage1_internal(cfg,x,beams,tag,A,K,N)
[zsort,ord]=sort(-x(:).','ascend'); invord=zeros(size(ord)); invord(ord)=1:numel(ord);
s=linspace(cfg.wall_support_m(1),cfg.wall_support_m(2),N).'; r=cfg.wall_r0_m-A*cos(K*s);
c=local_stage0_env_cfg(cfg,cfg.parametric_validation_exe,beams,zsort,tag);
c.wall_r0_m=cfg.wall_r0_m; c.wall_seed=cfg.wall_seed; c.wall_profile_r_m=r; c.wall_profile_z_m=s;
if isfield(cfg,'reuse_existing'), c.reuse_existing=cfg.reuse_existing; end
raw=run_bellhop_internal_pm_wall_poc_vertical(c); run=raw; run.pressure_raw=zeros(1,numel(x));
for ii=1:numel(x),run.pressure_raw(ii)=select_bellhop_shd_pressure_at_range_vertical(raw.data,cfg.mapped_receiver_range_m,zsort(ii),cfg.receiver_tolerance_m);end
run.pressure_raw=run.pressure_raw(invord); run.requested_receiver_depth_m=(-x(:)).';
run.receiver_coordinate_error_m=max(abs(double(raw.data.receiver_depth_m(:))-zsort(:)))*ones(1,numel(x)); run.profile_A_m=A;run.profile_K_radpm=K;run.profile_N=N;
end

function result = local_stage2(cfg, root)
% Stage 2: fixed-K height sweep.  Stage 1Y has fixed the comparison
% convention, so this stage may vary only the sinusoid height.  Each new
% case is independent and uses the frozen 10,001-beam endpoint.
stage1_file = fullfile(cfg.output_dir, cfg.stage1_output_subdir, ...
    'stage1_validation.mat');
if exist(stage1_file, 'file') ~= 2
    error('Stage 1 convention-fixed result is required before Stage 2.');
end
tmp = load(stage1_file, 'validation'); s1 = tmp.validation;
if ~isfield(s1, 'passed') || ~s1.passed
    error('Stage 1 convention-fixed regression did not pass; Stage 2 remains locked.');
end
stage0_file = fullfile(cfg.output_dir, 'stage0', 'stage0_validation.mat');
if exist(stage0_file, 'file') ~= 2
    error('Stage 0 result is required for the frozen receiver footprint.');
end
tmp0 = load(stage0_file, 'validation'); s0 = tmp0.validation;
if ~isfield(s0, 'passed') || ~s0.passed
    error('Stage 0 did not pass; Stage 2 remains locked.');
end

Avals = double(cfg.stage2_amplitudes_m(:).');
K = double(cfg.stage2_wavenumber_radpm);
if isempty(Avals) || any(~isfinite(Avals)) || any(Avals <= 0) || ...
        ~isfinite(K) || K <= 0
    error('Stage 2 requires positive finite amplitudes and wavenumber.');
end
if any(diff(Avals) < 0)
    error('Stage 2 amplitudes must be supplied in nondecreasing order.');
end

out_dir = fullfile(cfg.output_dir, cfg.stage2_output_subdir);
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
x = s0.x_m; fp = s0.footprint;
cases = cell(1, numel(Avals));
for ii = 1:numel(Avals)
    A = Avals(ii);
    case_dir = fullfile(out_dir, local_stage2_tag(A, K, cfg.wall_profile_count));
    case_file = fullfile(case_dir, 'stage2_case.mat');
    if exist(case_file, 'file') == 2
        cached = load(case_file, 'case_validation');
        if isfield(cached, 'case_validation') && ...
                abs(cached.case_validation.A_m - A) <= 10 * eps(max(1, abs(A))) && ...
                abs(cached.case_validation.K_radpm - K) <= 10 * eps(max(1, abs(K)))
            % Resume is case-granular: completed cases are never rerun.
            cc = cached.case_validation;
            cc.metrics = local_stage1_ratio_metrics(cc.G_PE, cc.G_BH_comparison, fp);
            cc.surface = local_sinusoid_surface_stats(A, K, cfg.wall_support_m, ...
                cfg.wall_profile_count);
            cc.checks.metrics_finite = local_metrics_finite(cc.metrics);
            cc.checks.all = all(structfun(@(v) logical(v), cc.checks));
            cases{ii} = cc;
            case_validation = cc; %#ok<NASGU>
            save(case_file, 'case_validation', '-v7');
            local_stage2_case_report(fullfile(case_dir, 'stage2_case_report.md'), cc);
            continue
        end
    end
    reuse = abs(A - s1.A_m) <= 10 * eps(max(1, abs(A)));
    if reuse && abs(K - s1.K_radpm) <= 10 * eps(max(1, abs(K)))
        % Reuse the already executed, convention-fixed A=0.01 case.  No
        % solver is called and the original Stage-1 MAT remains untouched.
        Gpe = s1.G_PE;
        Gbh_raw = s1.G_BH;
        Gbh_cmp = s1.G_BH_comparison;
        metrics = s1.metrics.model;
        geom = s1.geometry;
        source_tag = 'reused_stage1_convention_fixed';
    else
        pe_cfg = struct('frequency_hz', cfg.frequency_hz, 'c0_mps', cfg.c0_mps, ...
            'xw_m', cfg.xw_m, 'nx', cfg.nx, 'z_tx_m', cfg.z_tx_m, ...
            'z_rx_m', cfg.z_rx_m, 'sigma_src_m', cfg.sigma_src_m, ...
            'surface_elevation_x_m', A * cos(K * x), ...
            'surface_reflect_coeff', -1, 'step_m', cfg.pe_step_m, 'x_rx_m', 0);
        pe = run_pe_1d_surface_reflection_validation(pe_cfg);
        cfg_case = cfg;
        cfg_case.reuse_existing = true;
        bh = local_stage1_internal(cfg_case, x, cfg.stage1_beam_counts(2), ...
            local_stage2_tag(A, K, cfg.wall_profile_count), A, K, cfg.wall_profile_count);
        bh.field = local_stage0_convert(bh.pressure_raw, cfg.phase_sign);
        Gpe = zeros(size(s1.flat_pe)); Gbh_raw = zeros(size(s1.flat_bellhop));
        m = fp.m99.mask;
        Gpe(m) = pe.reflected_field(m) ./ s1.flat_pe(m);
        Gbh_raw(m) = bh.field(m) ./ s1.flat_bellhop(m);
        Gbh_cmp = conj(Gbh_raw);
        metrics = local_stage1_ratio_metrics(Gpe, Gbh_cmp, fp);
        geom = local_stage1_geometry(bh, cfg);
        source_tag = local_stage2_tag(A, K, cfg.wall_profile_count);
    end
    finite = all(isfinite([Gpe(:); Gbh_raw(:); Gbh_cmp(:)]));
    metric_finite = local_metrics_finite(metrics);
    checks = struct('geometry', logical(geom.all), 'finite', finite, ...
        'metrics_finite', metric_finite);
    checks.all = all(structfun(@(v) logical(v), checks));
    cases{ii} = struct('A_m', A, 'K_radpm', K, 'two_k_A', ...
        2 * (2*pi*cfg.frequency_hz/cfg.c0_mps) * A, ...
        'source_tag', source_tag, 'G_PE', Gpe, 'G_BH_raw', Gbh_raw, ...
        'G_BH_comparison', Gbh_cmp, 'metrics', metrics, ...
        'geometry', geom, 'surface', local_sinusoid_surface_stats(A, K, ...
        cfg.wall_support_m, cfg.wall_profile_count), 'checks', checks);
    if ~exist(case_dir, 'dir'), mkdir(case_dir); end
    case_validation = cases{ii}; %#ok<NASGU>
    save(fullfile(case_dir, 'stage2_case.mat'), 'case_validation', '-v7');
    local_stage2_case_report(fullfile(case_dir, 'stage2_case_report.md'), case_validation);
end

rows = repmat(struct('A_m',0,'K_radpm',0,'two_k_A',0,'E_G_M99',0, ...
    'max_slope',0,'max_abs_curvature_per_m',0,'rms_curvature_per_m',0, ...
    'phase_RMS_rad',0,'TL_RMS_dB',0,'phi0_rad',0,'rho_raw',0,'rho_shape',0, ...
    'E_aligned',0,'wall_residual_m',0,'phase_jump_error_rad',0, ...
    'p_rotation_error',0,'q_rotation_error',0,'min_post_dr_m',0, ...
    'center_tau_error_s',0,'geometry_passed',false,'finite_passed',false), ...
    1, numel(cases));
for ii = 1:numel(cases)
    cc = cases{ii}; mm = cc.metrics;
    rows(ii) = struct('A_m',cc.A_m,'K_radpm',cc.K_radpm,'two_k_A',cc.two_k_A, ...
        'max_slope',cc.surface.max_slope,'max_abs_curvature_per_m',cc.surface.max_abs_curvature_per_m, ...
        'rms_curvature_per_m',cc.surface.rms_curvature_per_m, ...
        'E_G_M99',mm.l2_m99,'phase_RMS_rad',mm.phase_rms_m99, ...
        'TL_RMS_dB',mm.tl_rms_m99,'phi0_rad',mm.global_phase_rad, ...
        'rho_raw',mm.rho_raw,'rho_shape',mm.rho_shape,'E_aligned',mm.aligned_l2_m99, ...
        'wall_residual_m',cc.geometry.wall_residual_max_m, ...
        'phase_jump_error_rad',cc.geometry.phase_jump_error_max_rad, ...
        'p_rotation_error',cc.geometry.p_rotation_error_max, ...
        'q_rotation_error',cc.geometry.q_rotation_error_max, ...
        'min_post_dr_m',cc.geometry.min_post_dr_m, ...
        'center_tau_error_s',cc.geometry.center_tau_error_s, ...
        'geometry_passed',cc.checks.geometry,'finite_passed',cc.checks.finite);
end
case_table = struct2table(rows);
checks = struct('stage1_gate', true, ...
    'all_case_geometry', all(case_table.geometry_passed), ...
    'all_case_finite', all(case_table.finite_passed), ...
    'all_case_metrics_finite', all(cellfun(@(c) c.checks.metrics_finite, cases)));
checks.all = all(structfun(@(v) logical(v), checks));
validation = struct('schema_version','1.0.0','stage','stage2_height_sweep', ...
    'config',cfg,'root',root,'K_radpm',K,'amplitudes_m',Avals, ...
    'comparison_convention','G_BH_comparison = conj(G_BH_abs) = G_BH_raw', ...
    'footprint',fp,'cases',{cases},'case_table',case_table, ...
    'checks',checks,'passed',checks.all,'solver_calls_allowed',true);
mat_file = fullfile(out_dir, 'stage2_validation.mat');
csv_file = fullfile(out_dir, 'stage2_metrics.csv');
report_file = fullfile(out_dir, 'stage2_report.md');
save(mat_file, 'validation', '-v7'); writetable(case_table, csv_file);
local_stage2_write_report(report_file, validation);
validation.files = struct('mat',mat_file,'csv',csv_file,'report',report_file);
result = validation;
if cfg.fail_on_check && ~result.passed
    error('Stage 2 failed geometry/finite checks; see %s.', report_file);
end
end

function result = local_stage3(cfg, root)
% Stage 3: fixed-height slope/curvature sweep. The height is frozen once
% from the completed Stage-2 validity classification and cannot be changed
% in response to the K-sweep results.
stage2_file = fullfile(cfg.output_dir, cfg.stage2_output_subdir, ...
    'stage2_validation.mat');
if exist(stage2_file, 'file') ~= 2
    error('Completed Stage 2 result is required before Stage 3.');
end
tmp = load(stage2_file, 'validation'); s2 = tmp.validation;
if ~isfield(s2, 'passed') || ~s2.passed
    error('Stage 2 did not pass its numerical/mapping guards; Stage 3 remains locked.');
end
stage1_file = fullfile(cfg.output_dir, cfg.stage1_output_subdir, ...
    'stage1_validation.mat');
tmp = load(stage1_file, 'validation'); s1 = tmp.validation;
stage0_file = fullfile(cfg.output_dir, 'stage0', 'stage0_validation.mat');
tmp = load(stage0_file, 'validation'); s0 = tmp.validation;

[fixed_A, selection_reason, stage2_regions] = local_stage3_select_amplitude(s2, s1.floor);
planned_K = double(cfg.stage3_wavenumbers_radpm(:).');
run_K = double(cfg.stage3_run_wavenumbers_radpm(:).');
if isempty(run_K), run_K = planned_K; end
if isempty(planned_K) || any(~isfinite(planned_K)) || any(planned_K <= 0) || ...
        any(diff(planned_K) <= 0)
    error('Stage 3 planned wavenumbers must be positive, finite, and strictly increasing.');
end
for ii = 1:numel(run_K)
    if ~any(abs(planned_K-run_K(ii)) <= 10*eps(max(1,abs(run_K(ii)))))
        error('Stage 3 run wavenumber %.9g is not in the frozen planned sweep.', run_K(ii));
    end
end

out_dir = fullfile(cfg.output_dir, cfg.stage3_output_subdir);
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
manifest_file = fullfile(out_dir, 'stage3_manifest.mat');
manifest = struct('schema_version','1.0.0','stage','stage3_manifest', ...
    'fixed_A_m',fixed_A,'selection_reason',selection_reason, ...
    'planned_K_radpm',planned_K,'beam_count',10001, ...
    'profile_count',cfg.wall_profile_count,'stage2_file',stage2_file, ...
    'stage2_sha256',local_sha256_file(stage2_file),'stage2_regions',{stage2_regions}, ...
    'comparison_convention','G_BH_comparison = conj(G_BH_abs) = G_BH_raw');
if exist(manifest_file, 'file') == 2
    old = load(manifest_file, 'manifest'); old = old.manifest;
    if abs(old.fixed_A_m-fixed_A) > 10*eps(max(1,abs(fixed_A))) || ...
            ~isequal(old.planned_K_radpm, planned_K)
        error('Existing Stage-3 manifest conflicts with the current frozen selection.');
    end
    manifest = old;
else
    save(manifest_file, 'manifest', '-v7');
    json_file = fullfile(out_dir, 'stage3_manifest.json');
    fid = fopen(json_file, 'w', 'n', 'UTF-8');
    if fid < 0, error('Cannot write %s.', json_file); end
    cl = onCleanup(@() fclose(fid)); fwrite(fid, jsonencode(manifest), 'char'); clear cl
end

x = s0.x_m; fp = s0.footprint;
for ii = 1:numel(run_K)
    K = run_K(ii);
    tag = local_stage3_tag(fixed_A, K, cfg.wall_profile_count);
    case_dir = fullfile(out_dir, tag);
    case_file = fullfile(case_dir, 'stage3_case.mat');
    if exist(case_file, 'file') == 2
        cached = load(case_file, 'case_validation'); cc = cached.case_validation;
        if abs(cc.A_m-fixed_A) <= 10*eps(max(1,abs(fixed_A))) && ...
                abs(cc.K_radpm-K) <= 10*eps(max(1,abs(K)))
            cc.metrics = local_stage1_ratio_metrics(cc.G_PE, cc.G_BH_comparison, fp);
            cc.surface = local_sinusoid_surface_stats(fixed_A, K, ...
                cfg.wall_support_m, cfg.wall_profile_count);
            cc.checks.metrics_finite = local_metrics_finite(cc.metrics);
            cc.checks.all = all(structfun(@(v) logical(v), cc.checks));
            case_validation = cc; %#ok<NASGU>
            save(case_file, 'case_validation', '-v7');
            local_stage3_case_report(fullfile(case_dir, 'stage3_case_report.md'), cc);
            continue
        end
    end

    stage2_case_file = fullfile(cfg.output_dir, cfg.stage2_output_subdir, ...
        local_stage2_tag(fixed_A, K, cfg.wall_profile_count), 'stage2_case.mat');
    if exist(stage2_case_file, 'file') == 2
        old = load(stage2_case_file, 'case_validation'); old = old.case_validation;
        Gpe = old.G_PE; Gbh_raw = old.G_BH_raw; Gbh_cmp = old.G_BH_comparison;
        geom = old.geometry; source_tag = old.source_tag;
        diag_file = fullfile(cfg.output_dir, 'bellhop', [source_tag '.iwdiag']);
        diagnostics = local_internal_wall_diagnostics(diag_file);
        reuse_note = 'reused matching Stage-2 case and raw Bellhop diagnostics';
    else
        pe_cfg = struct('frequency_hz',cfg.frequency_hz,'c0_mps',cfg.c0_mps, ...
            'xw_m',cfg.xw_m,'nx',cfg.nx,'z_tx_m',cfg.z_tx_m, ...
            'z_rx_m',cfg.z_rx_m,'sigma_src_m',cfg.sigma_src_m, ...
            'surface_elevation_x_m',fixed_A*cos(K*x), ...
            'surface_reflect_coeff',-1,'step_m',cfg.pe_step_m,'x_rx_m',0);
        pe = run_pe_1d_surface_reflection_validation(pe_cfg);
        cfg_case = cfg; cfg_case.reuse_existing = true;
        source_tag = tag;
        bh = local_stage1_internal(cfg_case, x, 10001, source_tag, ...
            fixed_A, K, cfg.wall_profile_count);
        bh.field = local_stage0_convert(bh.pressure_raw, cfg.phase_sign);
        m99 = fp.m99.mask;
        Gpe = zeros(size(s1.flat_pe)); Gbh_raw = zeros(size(s1.flat_bellhop));
        Gpe(m99) = pe.reflected_field(m99)./s1.flat_pe(m99);
        Gbh_raw(m99) = bh.field(m99)./s1.flat_bellhop(m99);
        Gbh_cmp = conj(Gbh_raw);
        geom = local_stage1_geometry(bh, cfg);
        diagnostics = bh.diagnostics;
        reuse_note = 'new Stage-3 PE and Bellhop case (raw Bellhop output is resumable)';
    end
    metrics = local_stage1_ratio_metrics(Gpe, Gbh_cmp, fp);
    surface = local_sinusoid_surface_stats(fixed_A, K, cfg.wall_support_m, ...
        cfg.wall_profile_count);
    beam_diagnostics = local_stage3_select_diagnostics(diagnostics);
    diagnostic_summary = local_stage3_diagnostic_summary(diagnostics, fixed_A, K);
    checks = struct('geometry',logical(geom.all), ...
        'field_finite',all(isfinite([Gpe(:);Gbh_raw(:);Gbh_cmp(:)])), ...
        'diagnostics_finite',diagnostic_summary.finite, ...
        'metrics_finite',local_metrics_finite(metrics));
    checks.all = all(structfun(@(v) logical(v), checks));
    case_validation = struct('schema_version','1.0.0','stage','stage3_case', ...
        'A_m',fixed_A,'K_radpm',K,'source_tag',source_tag,'reuse_note',reuse_note, ...
        'G_PE',Gpe,'G_BH_raw',Gbh_raw,'G_BH_comparison',Gbh_cmp, ...
        'metrics',metrics,'surface',surface,'geometry',geom, ...
        'diagnostic_summary',diagnostic_summary, ...
        'beam_diagnostics',beam_diagnostics,'diagnostic_file',diag_file_if_needed(cfg,source_tag), ...
        'checks',checks);
    if ~exist(case_dir, 'dir'), mkdir(case_dir); end
    save(case_file, 'case_validation', '-v7');
    local_stage3_case_report(fullfile(case_dir, 'stage3_case_report.md'), case_validation);
    if ~checks.all
        error('Stage 3 numerical/mapping guard failed for K=%.9g.', K);
    end
end

cases = cell(1,numel(planned_K)); complete = true(1,numel(planned_K));
for ii = 1:numel(planned_K)
    f = fullfile(out_dir, local_stage3_tag(fixed_A,planned_K(ii),cfg.wall_profile_count), ...
        'stage3_case.mat');
    if exist(f,'file')==2
        q=load(f,'case_validation'); cases{ii}=q.case_validation;
    else
        complete(ii)=false;
    end
end
done_cases = cases(complete);
case_table = local_stage3_table(done_cases, s1.floor);
checks = struct('stage2_gate',true,'all_cases_complete',all(complete), ...
    'all_case_geometry',all(cellfun(@(c)c.checks.geometry,done_cases)), ...
    'all_case_finite',all(cellfun(@(c)c.checks.field_finite && ...
    c.checks.diagnostics_finite && c.checks.metrics_finite,done_cases)));
checks.all = all(structfun(@(v)logical(v),checks));
validation = struct('schema_version','1.0.0','stage','stage3_slope_curvature_sweep', ...
    'config',cfg,'root',root,'manifest',manifest,'fixed_A_m',fixed_A, ...
    'planned_K_radpm',planned_K,'completed_K_radpm',planned_K(complete), ...
    'cases',{done_cases},'case_table',case_table,'checks',checks,'passed',checks.all);
mat_file=fullfile(out_dir,'stage3_validation.mat'); csv_file=fullfile(out_dir,'stage3_metrics.csv');
report_file=fullfile(out_dir,'stage3_report.md'); save(mat_file,'validation','-v7');
writetable(case_table,csv_file); local_stage3_write_report(report_file,validation);
validation.files=struct('manifest',manifest_file,'mat',mat_file,'csv',csv_file,'report',report_file);
result=validation;
if all(complete) && cfg.fail_on_check && ~result.passed
    error('Stage 3 failed numerical/mapping guards; see %s.', report_file);
end
end

function [A, reason, regions] = local_stage3_select_amplitude(s2, floor)
t=s2.case_table; regions=strings(height(t),1); is_I=false(height(t),1);
for ii=1:height(t)
    is_I(ii)=t.E_G_M99(ii)<=floor.T_E && t.E_aligned(ii)<=floor.T_E && ...
        t.TL_RMS_dB(ii)<=floor.T_TL && t.phase_RMS_rad(ii)<=floor.T_phi && ...
        t.rho_shape(ii)>=0.9995 && t.geometry_passed(ii) && t.finite_passed(ii);
    if is_I(ii)
        regions(ii)="I";
    elseif t.E_aligned(ii)<=0.10 && t.TL_RMS_dB(ii)<=0.5 && ...
            t.phase_RMS_rad(ii)<=0.20 && t.rho_shape(ii)>=0.99
        regions(ii)="II";
    else
        regions(ii)="III";
    end
end
if any(is_I)
    A=max(t.A_m(is_I));
    reason='largest Stage-2 amplitude satisfying all frozen Region-I gates';
elseif any(abs(t.A_m-0.05)<1e-12) && all(~is_I(t.A_m>=0.05))
    A=0.05;
    reason='Goal fallback: all sampled A >= 0.05 m are outside Region I';
else
    error('Stage-3 fixed-height selection rule is not resolvable from Stage 2.');
end
end

function tag=local_stage3_tag(A,K,N)
tag=strrep(sprintf('slope_A%.6g_K%.6g_N%d_B10001',A,K,N),'.','p');
end

function path=diag_file_if_needed(cfg,tag)
path=fullfile(cfg.output_dir,'bellhop',[tag '.iwdiag']);
end

function d=local_internal_wall_diagnostics(path)
if exist(path,'file')~=2, error('Missing internal-wall diagnostics: %s',path); end
m=readmatrix(path,'FileType','text','CommentStyle','#'); m=m(~all(isnan(m),2),:);
names={'alpha_deg','hit_r','hit_z','wall_residual','wall_t_r','wall_t_z','wall_n_r','wall_n_z', ...
    'tangent_error','normal_error','inc_ur','inc_uz','ref_ur','ref_uz','rot_ur','rot_uz', ...
    'specular_error','rotation_error','phase_in','phase_ref','phase_delta','amp_in','amp_ref','amp_delta', ...
    'p1_in','p2_in','p1_ref','p2_ref','p_ref_error','q1_in','q2_in','q1_ref','q2_ref', ...
    'q_ref_error','p_rot_error','q_rot_error','tau_wall_real','tau_wall_imag', ...
    'tau_receiver_real','tau_receiver_imag','min_post_dr','n_post','kappa', ...
    'wall_seg','wall_lambda','wall_tg','wall_th','wall_rm','wall_rn'};
if size(m,2)~=numel(names), error('Unexpected diagnostic width in %s.',path); end
d=array2table(m,'VariableNames',names);
end

function out=local_stage3_select_diagnostics(d)
names={'alpha_deg','hit_r','hit_z','wall_residual','wall_t_r','wall_t_z', ...
    'wall_n_r','wall_n_z','inc_ur','inc_uz','ref_ur','ref_uz','rot_ur','rot_uz', ...
    'phase_delta','amp_delta','p1_in','p2_in','p1_ref','p2_ref','q1_in','q2_in', ...
    'q1_ref','q2_ref','tau_wall_real','tau_receiver_real','min_post_dr','kappa', ...
    'wall_tg','wall_th','wall_rm','wall_rn'};
out=d(:,names);
end

function s=local_stage3_diagnostic_summary(d,A,K)
analytic=-A*K^2*cos(K*d.hit_z)./(1+(A*K*sin(K*d.hit_z)).^2).^(3/2);
e=d.kappa-analytic;
s=struct('ray_count',height(d),'finite',all(isfinite(d{:,:}),'all'), ...
    'hit_kappa_max_abs_per_m',max(abs(d.kappa)), ...
    'hit_kappa_rms_per_m',sqrt(mean(d.kappa.^2)), ...
    'kappa_error_max_abs_per_m',max(abs(e)), ...
    'kappa_error_rms_per_m',sqrt(mean(e.^2)), ...
    'min_abs_incidence_Th',min(abs(d.wall_th)), ...
    'max_abs_RN',max(abs(d.wall_rn)),'max_abs_RM',max(abs(d.wall_rm)), ...
    'wall_residual_max_m',max(abs(d.wall_residual)), ...
    'phase_jump_error_max_rad',max(abs(d.phase_delta-pi)), ...
    'p_ref_error_max',max(abs(d.p_ref_error)), ...
    'q_ref_error_max',max(abs(d.q_ref_error)), ...
    'p_rotation_error_max',max(abs(d.p_rot_error)), ...
    'q_rotation_error_max',max(abs(d.q_rot_error)), ...
    'min_post_dr_m',min(d.min_post_dr), ...
    'tau_min_s',min(d.tau_receiver_real),'tau_max_s',max(d.tau_receiver_real));
end

function t=local_stage3_table(cases,floor)
rows=repmat(struct('A_m',0,'K_radpm',0,'max_slope',0, ...
    'max_abs_curvature_per_m',0,'rms_curvature_per_m',0,'minimum_radius_m',0, ...
    'E_G_M99',0,'phase_RMS_rad',0,'TL_RMS_dB',0,'phi0_rad',0, ...
    'rho_raw',0,'rho_shape',0,'E_aligned',0,'hit_kappa_max_abs_per_m',0, ...
    'kappa_error_max_abs_per_m',0,'min_abs_incidence_Th',0,'max_abs_RN',0, ...
    'max_abs_RM',0,'wall_residual_m',0,'min_post_dr_m',0,'region',"", ...
    'geometry_passed',false,'finite_passed',false),1,numel(cases));
for ii=1:numel(cases)
    c=cases{ii};m=c.metrics;s=c.surface;d=c.diagnostic_summary;
    region="III";
    if m.l2_m99<=floor.T_E && m.aligned_l2_m99<=floor.T_E && ...
            m.tl_rms_m99<=floor.T_TL && m.phase_rms_m99<=floor.T_phi && m.rho_shape>=.9995
        region="I";
    elseif m.aligned_l2_m99<=.10 && m.tl_rms_m99<=.5 && ...
            m.phase_rms_m99<=.20 && m.rho_shape>=.99
        region="II";
    end
    rows(ii)=struct('A_m',c.A_m,'K_radpm',c.K_radpm,'max_slope',s.max_slope, ...
        'max_abs_curvature_per_m',s.max_abs_curvature_per_m, ...
        'rms_curvature_per_m',s.rms_curvature_per_m,'minimum_radius_m',s.minimum_radius_m, ...
        'E_G_M99',m.l2_m99,'phase_RMS_rad',m.phase_rms_m99,'TL_RMS_dB',m.tl_rms_m99, ...
        'phi0_rad',m.global_phase_rad,'rho_raw',m.rho_raw,'rho_shape',m.rho_shape, ...
        'E_aligned',m.aligned_l2_m99,'hit_kappa_max_abs_per_m',d.hit_kappa_max_abs_per_m, ...
        'kappa_error_max_abs_per_m',d.kappa_error_max_abs_per_m, ...
        'min_abs_incidence_Th',d.min_abs_incidence_Th,'max_abs_RN',d.max_abs_RN, ...
        'max_abs_RM',d.max_abs_RM,'wall_residual_m',d.wall_residual_max_m, ...
        'min_post_dr_m',d.min_post_dr_m,'region',region, ...
        'geometry_passed',c.checks.geometry, ...
        'finite_passed',c.checks.field_finite && c.checks.diagnostics_finite && c.checks.metrics_finite);
end
t=struct2table(rows);
end

function local_stage3_case_report(path,v)
fid=fopen(path,'w','n','UTF-8');if fid<0,error('Cannot write %s.',path);end
cl=onCleanup(@()fclose(fid));m=v.metrics;s=v.surface;d=v.diagnostic_summary;
fprintf(fid,'# Stage 3 slope/curvature case\n\n- A=%.9g m; K=%.9g rad/m; max slope %.9g.\n',v.A_m,v.K_radpm,s.max_slope);
fprintf(fid,'- max/RMS curvature: %.9g / %.9g 1/m; minimum radius %.9g m.\n',s.max_abs_curvature_per_m,s.rms_curvature_per_m,s.minimum_radius_m);
fprintf(fid,'- %s.\n\n',v.reuse_note);
fprintf(fid,'| E_G | phase RMS | TL RMS (dB) | phi0 | rho_raw | rho_shape | E_aligned |\n|---:|---:|---:|---:|---:|---:|---:|\n');
fprintf(fid,'| %.9g | %.9g | %.9g | %.9g | %.9g | %.9g | %.9g |\n\n',m.l2_m99,m.phase_rms_m99,m.tl_rms_m99,m.global_phase_rad,m.rho_raw,m.rho_shape,m.aligned_l2_m99);
fprintf(fid,'Bellhop diagnostics: hit max |kappa| %.9g 1/m; max curvature error %.3g 1/m; min |Th| %.9g; max |RN/RM| %.9g / %.9g; wall residual %.3g m; min post dr %.9g m.\n', ...
    d.hit_kappa_max_abs_per_m,d.kappa_error_max_abs_per_m,d.min_abs_incidence_Th,d.max_abs_RN,d.max_abs_RM,d.wall_residual_max_m,d.min_post_dr_m);
end

function local_stage3_write_report(path,v)
fid=fopen(path,'w','n','UTF-8');if fid<0,error('Cannot write %s.',path);end
cl=onCleanup(@()fclose(fid));
fprintf(fid,'# Stage 3 fixed-height slope/curvature sweep\n\n');
fprintf(fid,'状态：**%s**；fixed A=`%.9g m`。\n\n',ternary(v.passed,'PASS',ternary(v.checks.all_cases_complete,'FAIL','INCOMPLETE')),v.fixed_A_m);
fprintf(fid,'选择依据：%s。Bellhop 固定 `X/C`, 10,001 beams, profile N=%d。\n\n',v.manifest.selection_reason,v.manifest.profile_count);
fprintf(fid,'| K (rad/m) | max slope | max/RMS curvature (1/m) | min radius (m) | E_G | phase RMS | TL RMS (dB) | phi0 | rho_raw | rho_shape | E_aligned | hit max |kappa| | kappa max error | min |Th| | max |RN/RM| | wall residual | min post dr | region | guards |\n');
fprintf(fid,'|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---:|:---:|\n');
for ii=1:height(v.case_table)
    r=v.case_table(ii,:);
    fprintf(fid,'| %.9g | %.9g | %.9g/%.9g | %.9g | %.9g | %.9g | %.9g | %.9g | %.9g | %.9g | %.9g | %.9g | %.3g | %.9g | %.9g/%.9g | %.3g | %.9g | %s | %s |\n', ...
        r.K_radpm,r.max_slope,r.max_abs_curvature_per_m,r.rms_curvature_per_m,r.minimum_radius_m, ...
        r.E_G_M99,r.phase_RMS_rad,r.TL_RMS_dB,r.phi0_rad,r.rho_raw,r.rho_shape,r.E_aligned, ...
        r.hit_kappa_max_abs_per_m,r.kappa_error_max_abs_per_m,r.min_abs_incidence_Th,r.max_abs_RN,r.max_abs_RM, ...
        r.wall_residual_m,r.min_post_dr_m,char(r.region),ternary(r.geometry_passed&&r.finite_passed,'PASS','FAIL'));
end
fprintf(fid,'\nModel discrepancy and Region transition do not stop Stage 3. Only solver, mapping, geometry, or finite-state guards are hard gates. Stage 4 is not executed here.\n');
end

function result=local_stage4(cfg,root)
% Stage 4: read-only sampled validity map. No solver call is permitted.
s1f=fullfile(cfg.output_dir,cfg.stage1_output_subdir,'stage1_validation.mat');
s2f=fullfile(cfg.output_dir,cfg.stage2_output_subdir,'stage2_validation.mat');
s3f=fullfile(cfg.output_dir,cfg.stage3_output_subdir,'stage3_validation.mat');
for f={s1f,s2f,s3f}
    if exist(f{1},'file')~=2,error('Stage 4 requires %s.',f{1});end
end
q=load(s1f,'validation');s1=q.validation;q=load(s2f,'validation');s2=q.validation;
q=load(s3f,'validation');s3=q.validation;
if ~s1.passed || ~s2.passed || ~s3.passed
    error('Stage 4 requires PASS results from Stages 1--3.');
end

cases=s2.cases;
duplicate_error=NaN;
for ii=1:numel(s3.cases)
    c=s3.cases{ii};
    duplicate=find(cellfun(@(v)abs(v.A_m-c.A_m)<1e-12 && abs(v.K_radpm-c.K_radpm)<1e-12,cases),1);
    if isempty(duplicate)
        cases{end+1}=c; %#ok<AGROW>
    else
        duplicate_error=max(abs(c.G_PE-cases{duplicate}.G_PE),[],'all');
        duplicate_error=max(duplicate_error,max(abs(c.G_BH_comparison-cases{duplicate}.G_BH_comparison),[],'all'));
    end
end
rows=repmat(struct('A_m',0,'K_radpm',0,'E_G_M99',0,'E_aligned',0, ...
    'TL_RMS_dB',0,'phase_RMS_rad',0,'phi0_rad',0,'rho_raw',0,'rho_shape',0, ...
    'region',"",'numerical_mapping_guards',false),1,numel(cases));
for ii=1:numel(cases)
    c=cases{ii};m=c.metrics;
    region=local_stage4_region(m,s1.floor);
    if isfield(c.checks,'field_finite')
        guard=c.checks.geometry && c.checks.field_finite && ...
            c.checks.diagnostics_finite && c.checks.metrics_finite;
    else
        guard=c.checks.geometry && c.checks.finite && c.checks.metrics_finite;
    end
    rows(ii)=struct('A_m',c.A_m,'K_radpm',c.K_radpm,'E_G_M99',m.l2_m99, ...
        'E_aligned',m.aligned_l2_m99,'TL_RMS_dB',m.tl_rms_m99, ...
        'phase_RMS_rad',m.phase_rms_m99,'phi0_rad',m.global_phase_rad, ...
        'rho_raw',m.rho_raw,'rho_shape',m.rho_shape,'region',region, ...
        'numerical_mapping_guards',guard);
end
t=struct2table(rows);t=sortrows(t,{'K_radpm','A_m'});
height_slice=t(abs(t.K_radpm-cfg.stage2_wavenumber_radpm)<1e-12,:);
slope_slice=t(abs(t.A_m-s3.fixed_A_m)<1e-12,:);
intervals=struct( ...
    'height_region_I_sampled_m',height_slice.A_m(height_slice.region=="I").', ...
    'height_nonI_sampled_m',height_slice.A_m(height_slice.region~="I").', ...
    'height_transition_bracket_m',[max(height_slice.A_m(height_slice.region=="I")), ...
        min(height_slice.A_m(height_slice.region~="I"))], ...
    'wavenumber_region_I_sampled_radpm',slope_slice.K_radpm(slope_slice.region=="I").', ...
    'wavenumber_nonI_sampled_radpm',slope_slice.K_radpm(slope_slice.region~="I").', ...
    'wavenumber_transition_bracket_radpm',[max(slope_slice.K_radpm(slope_slice.region=="I")), ...
        min(slope_slice.K_radpm(slope_slice.region~="I"))]);
checks=struct('stage1_pass',s1.passed,'stage2_pass',s2.passed,'stage3_pass',s3.passed, ...
    'all_source_guards',all(t.numerical_mapping_guards), ...
    'duplicate_consistent',isfinite(duplicate_error)&&duplicate_error<=1e-14, ...
    'regions_assigned',all(ismember(t.region,["I","II","III"])));
checks.all=all(structfun(@(v)logical(v),checks));
out_dir=fullfile(cfg.output_dir,cfg.stage4_output_subdir);if ~exist(out_dir,'dir'),mkdir(out_dir);end
validation=struct('schema_version','1.0.0','stage','stage4_validity_map', ...
    'root',root,'config',cfg,'thresholds',s1.floor,'fixed_A_m',s3.fixed_A_m, ...
    'case_table',t, ...
    'sampled_intervals',intervals,'duplicate_case_error',duplicate_error, ...
    'checks',checks,'passed',checks.all,'solver_calls',0, ...
    'interpretation','sampled intervals only; no boundary interpolation or extrapolation');
mat_file=fullfile(out_dir,'stage4_validation.mat');csv_file=fullfile(out_dir,'stage4_validity_map.csv');
report_file=fullfile(out_dir,'stage4_report.md');save(mat_file,'validation','-v7');
writetable(t,csv_file);local_stage4_write_report(report_file,validation);
validation.files=struct('mat',mat_file,'csv',csv_file,'report',report_file);result=validation;
if cfg.fail_on_check && ~result.passed,error('Stage 4 failed; see %s.',report_file);end
end

function region=local_stage4_region(m,floor)
if m.l2_m99<=floor.T_E && m.aligned_l2_m99<=floor.T_E && ...
        m.tl_rms_m99<=floor.T_TL && m.phase_rms_m99<=floor.T_phi && m.rho_shape>=.9995
    region="I";
elseif m.aligned_l2_m99<=.10 && m.tl_rms_m99<=.5 && ...
        m.phase_rms_m99<=.20 && m.rho_shape>=.99
    region="II";
else
    region="III";
end
end

function local_stage4_write_report(path,v)
fid=fopen(path,'w','n','UTF-8');if fid<0,error('Cannot write %s.',path);end
cl=onCleanup(@()fclose(fid));
fprintf(fid,'# Stage 4 sampled validity map\n\n状态：**%s**；solver calls: `0`。\n\n',ternary(v.passed,'PASS','FAIL'));
fprintf(fid,'| A (m) | K (rad/m) | E_G | E_aligned | TL RMS (dB) | phase RMS (rad) | phi0 | rho_raw | rho_shape | region | guards |\n');
fprintf(fid,'|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---:|:---:|\n');
for ii=1:height(v.case_table)
    r=v.case_table(ii,:);
    fprintf(fid,'| %.9g | %.9g | %.9g | %.9g | %.9g | %.9g | %.9g | %.9g | %.9g | %s | %s |\n', ...
        r.A_m,r.K_radpm,r.E_G_M99,r.E_aligned,r.TL_RMS_dB,r.phase_RMS_rad, ...
        r.phi0_rad,r.rho_raw,r.rho_shape,char(r.region),ternary(r.numerical_mapping_guards,'PASS','FAIL'));
end
h=v.sampled_intervals.height_transition_bracket_m;k=v.sampled_intervals.wavenumber_transition_bracket_radpm;
fprintf(fid,'\n## Sampled intervals only\n\n');
fprintf(fid,'- At K=0.10 rad/m, the sampled Region-I to non-I transition is bracketed by A=(%.9g, %.9g] m.\n',h(1),h(2));
fprintf(fid,'- At A=%.9g m, the sampled Region-I to non-I transition is bracketed by K=(%.9g, %.9g] rad/m.\n',v.fixed_A_m,k(1),k(2));
fprintf(fid,'- These are sampled brackets, not interpolated or extrapolated validity boundaries.\n');
end

function tag = local_stage2_tag(A, K, N)
tag = strrep(sprintf('height_A%.6g_K%.6g_N%d_B10001', A, K, N), '.', 'p');
end

function local_stage2_case_report(path, v)
fid = fopen(path, 'w', 'n', 'UTF-8'); if fid < 0, error('Cannot write %s.', path); end
cl = onCleanup(@() fclose(fid)); m = v.metrics;
fprintf(fid, '# Stage 2 height-sweep case\n\n');
fprintf(fid, '- A=%.9g m, K=%.9g rad/m, `2kA=%.9g`.\n', v.A_m, v.K_radpm, v.two_k_A);
fprintf(fid, '- source tag: `%s`\n\n', v.source_tag);
fprintf(fid, '| max slope | max/RMS curvature (1/m) | E_G (M99) | phase RMS (rad) | TL RMS (dB) | phi0 (rad) | rho_raw | rho_shape | E_aligned | geometry | finite |\n');
fprintf(fid, '|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---:|:---:|\n');
fprintf(fid, '| %.9g | %.9g / %.9g | %.9g | %.9g | %.9g | %.9g | %.9g | %.9g | %.9g | %s | %s |\n\n', ...
    v.surface.max_slope, v.surface.max_abs_curvature_per_m, v.surface.rms_curvature_per_m, ...
    m.l2_m99, m.phase_rms_m99, m.tl_rms_m99, m.global_phase_rad, ...
    m.rho_raw, m.rho_shape, m.aligned_l2_m99, ternary(v.checks.geometry,'PASS','FAIL'), ...
    ternary(v.checks.finite,'PASS','FAIL'));
fprintf(fid, 'Geometry guards: wall residual `%.3g m`; pressure-release phase-jump error `%.3g rad`; p/q rotation errors `%.3g/%.3g`; min post-wall dr `%.6g m`; center tau error `%.3g s`.\n\n', ...
    v.geometry.wall_residual_max_m, v.geometry.phase_jump_error_max_rad, ...
    v.geometry.p_rotation_error_max, v.geometry.q_rotation_error_max, ...
    v.geometry.min_post_dr_m, v.geometry.center_tau_error_s);
fprintf(fid, 'Comparison convention: `G_BH_comparison = conj(G_BH_abs) = G_BH_raw`.\n');
end

function local_stage2_write_report(path, v)
fid = fopen(path, 'w', 'n', 'UTF-8'); if fid < 0, error('Cannot write %s.', path); end
cl = onCleanup(@() fclose(fid));
fprintf(fid, '# Stage 2 fixed-K height-phase sweep\n\n');
fprintf(fid, '状态：**%s**\n\n', ternary(v.passed,'PASS','FAIL'));
fprintf(fid, '本阶段在 Stage 1Y `CONVENTION_THEORETICALLY_CLOSED` 后执行，固定 `K=%.9g rad/m`、source `%s`、run `%s`、Bellhop beam count `%d`。A=%.9g m 复用已通过的 Stage-1 convention-fixed MAT；其余点各运行一次相同配置。\n\n', ...
    v.K_radpm, v.config.source_geometry, v.config.run_type, v.config.stage1_beam_counts(2), v.amplitudes_m(1));
fprintf(fid, '| A (m) | 2kA | max slope | max/RMS curvature (1/m) | E_G M99 | phase RMS (rad) | TL RMS (dB) | phi0 (rad) | rho_raw | rho_shape | E_aligned | wall residual (m) | phase jump error (rad) | p/q rotation error | min post dr (m) | tau error (s) | geometry | finite |\n');
fprintf(fid, '|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---:|:---:|\n');
for ii = 1:height(v.case_table)
    r = v.case_table(ii,:);
    fprintf(fid, '| %.9g | %.9g | %.9g | %.9g/%.9g | %.9g | %.9g | %.9g | %.9g | %.9g | %.9g | %.9g | %.3g | %.3g | %.3g/%.3g | %.6g | %.3g | %s | %s |\n', ...
        r.A_m, r.two_k_A, r.max_slope, r.max_abs_curvature_per_m, r.rms_curvature_per_m, ...
        r.E_G_M99, r.phase_RMS_rad, r.TL_RMS_dB, r.phi0_rad, r.rho_raw, ...
        r.rho_shape, r.E_aligned, r.wall_residual_m, ...
        r.phase_jump_error_rad, r.p_rotation_error, r.q_rotation_error, ...
        r.min_post_dr_m, r.center_tau_error_s, ...
        ternary(r.geometry_passed,'PASS','FAIL'), ternary(r.finite_passed,'PASS','FAIL'));
end
fprintf(fid, '\n## Checks\n\n- Stage 1 gate: PASS\n- all case geometry: %s\n- all case finite: %s\n- all case metrics finite: %s\n\n', ...
    ternary(v.checks.all_case_geometry,'PASS','FAIL'), ternary(v.checks.all_case_finite,'PASS','FAIL'), ...
    ternary(v.checks.all_case_metrics_finite,'PASS','FAIL'));
fprintf(fid, '模型差异不触发 Stage 2 停止；本阶段仅以 solver/mapping/finite guards 作为 hard gate。\n');
end

function tag=local_stage1_tag(A,K,N)
tag=strrep(sprintf('sinusoid_A%.6g_K%.6g_N%d',A,K,N),'.','p');
end

function m=local_stage1_ratio_metrics(a,b,fp)
a=a(:).';b=b(:).';w=fp.energy_weights(:).';m99=fp.m99.mask;valid=fp.m95.mask & abs(a)>=10^(fp.phase_floor_db/20) & abs(b)>=10^(fp.phase_floor_db/20);ph=angle(a.*conj(b));tl=20*log10(max(abs(a),realmin)./max(abs(b),realmin));ww=w(m99);S=sum(ww.*a(m99).*conj(b(m99)));den=sqrt(max(sum(ww.*abs(a(m99)).^2)*sum(ww.*abs(b(m99)).^2),realmin));c=S/den;rho=abs(c);g=angle(S);m=struct('l2_m99',sqrt(sum(ww.*abs(a(m99)-b(m99)).^2)/max(sum(ww.*abs(b(m99)).^2),realmin)),'phase_rms_m99',sqrt(sum(ww.*ph(m99).^2)/max(sum(ww),realmin)),'tl_rms_m99',sqrt(sum(ww.*tl(m99).^2)/max(sum(ww),realmin)),'phase_p95_m95',local_stage0_percentile(abs(ph(valid)),.95),'tl_p95_m95',local_stage0_percentile(abs(tl(valid)),.95),'rho_raw',real(c),'rho_shape',rho,'global_phase_rad',g,'aligned_l2_m99',sqrt(sum(ww.*abs(a(m99)-exp(1i*g)*b(m99)).^2)/max(sum(ww.*abs(a(m99)).^2),realmin)),'phase_difference',ph,'tl_difference_db',tl);
end

function tf = local_metrics_finite(m)
tf = all(isfinite([m.l2_m99, m.phase_rms_m99, m.tl_rms_m99, ...
    m.global_phase_rad, m.rho_raw, m.rho_shape, m.aligned_l2_m99]));
end

function stats = local_sinusoid_surface_stats(A, K, support_m, count)
s = linspace(support_m(1), support_m(2), count);
eta_prime = -A * K * sin(K * s);
eta_second = -A * K^2 * cos(K * s);
kappa = eta_second ./ (1 + eta_prime.^2).^(3/2);
stats = struct('sample_count',count,'support_m',support_m, ...
    'max_slope',max(abs(eta_prime)), ...
    'max_abs_curvature_per_m',max(abs(kappa)), ...
    'rms_curvature_per_m',sqrt(trapz(s, kappa.^2) / (s(end)-s(1))), ...
    'minimum_radius_m',1/max(abs(kappa)), ...
    'analytic_max_abs_curvature_per_m',A*K^2, ...
    'sample_s_m',s,'sample_kappa_per_m',kappa);
end

function g=local_stage1_geometry(run,cfg)
d=run.diagnostics; [~,ic]=min(abs(d.alpha_deg));
g=struct('wall_residual_max_m',max(abs(d.wall_residual)),'phase_jump_error_max_rad',max(abs(d.phase_delta-pi)), ...
 'p_ref_error_max',max(abs(d.p_ref_error)),'q_ref_error_max',max(abs(d.q_ref_error)), ...
 'p_rotation_error_max',max(abs(d.p_rot_error)),'q_rotation_error_max',max(abs(d.q_rot_error)), ...
 'min_post_dr_m',min(d.min_post_dr),'center_tau_error_s',abs(d.tau_receiver_real(ic)-(cfg.mapped_receiver_range_m-2*run.profile_A_m)/cfg.c0_mps), ...
 'all',all(isfinite(d{:,:}),'all') && max(abs(d.wall_residual))<=1e-9 && max(abs(d.phase_delta-pi))<=1e-10 && ...
 max(abs(d.p_rot_error))<=1e-12 && max(abs(d.q_rot_error))<=1e-12 && min(d.min_post_dr)>0 && abs(d.tau_receiver_real(ic)-(cfg.mapped_receiver_range_m-2*run.profile_A_m)/cfg.c0_mps)<=1e-6);
end

function local_stage1_write_report(path,v)
fid=fopen(path,'w','n','UTF-8');if fid<0,error('Cannot write %s.',path);end;cl=onCleanup(@()fclose(fid));fprintf(fid,'# Stage 1 weak sinusoid\n\n状态：**%s**\n\n',ternary(v.passed,'PASS','FAIL'));fprintf(fid,'- A=%.9g m, K=%.9g rad/m; source `%s`, run `%s`.\n',v.A_m,v.K_radpm,v.config.source_geometry,v.config.run_type);fprintf(fid,'- Floors: F_E=%.8g, F_phi=%.8g, F_TL=%.8g; thresholds T_E=%.8g, T_phi=%.8g, T_TL=%.8g.\n\n',v.floor.F_E,v.floor.F_phi,v.floor.F_TL,v.floor.T_E,v.floor.T_phi,v.floor.T_TL);fprintf(fid,'| comparison | L2(M99) | phase RMS | TL RMS (dB) | rho |\n|---|---:|---:|---:|---:|\n');for n={'profile','beam','model'},m=v.metrics.(n{1});fprintf(fid,'| %s | %.8g | %.8g | %.8g | %.8g |\n',n{1},m.l2_m99,m.phase_rms_m99,m.tl_rms_m99,m.rho_shape);end;fprintf(fid,'\n## Geometry\n\n- wall residual max: %.4g m; phase jump error max: %.4g rad; min post dr: %.6g m; center tau error: %.4g s.\n\n',v.geometry.wall_residual_max_m,v.geometry.phase_jump_error_max_rad,v.geometry.min_post_dr_m,v.geometry.center_tau_error_s);fprintf(fid,'## Checks\n\n');n=fieldnames(v.checks);for ii=1:numel(n),fprintf(fid,'- %s: %s\n',n{ii},ternary(v.checks.(n{ii}),'PASS','FAIL'));end;fprintf(fid,'\nStage 2 remains locked unless this weak-limit result passes.\n');
end

function local_stage1_append_convention_report(path,v)
fid=fopen(path,'a','n','UTF-8');if fid<0,error('Cannot append %s.',path);end;cl=onCleanup(@()fclose(fid));
if isfield(v,'comparison_convention')
    fprintf(fid,'\n## Comparison convention\n\n- %s\n',v.comparison_convention);
    if isfield(v.metrics,'model_raw')
        m=v.metrics.model_raw;
        fprintf(fid,'- Raw diagnostic: E_G=%.8g; phase RMS=%.8g rad; TL RMS=%.8g dB; rho=%.8g.\n',m.l2_m99,m.phase_rms_m99,m.tl_rms_m99,m.rho_shape);
    end
end
end

function local_stage1_blocked_report(path,v)
fid=fopen(path,'w','n','UTF-8');if fid<0,error('Cannot write %s.',path);end
cl=onCleanup(@()fclose(fid));fprintf(fid,'# Stage 1 weak sinusoid\n\n状态：**BLOCKED / NOT COMPLETED**\n\n');
fprintf(fid,'本阶段未形成物理 PASS/FAIL 结论。首个高 beam 粗糙 internal-wall case 未在本次受控运行窗口内完成，已停止，未启动后续 case。\n\n');
fprintf(fid,'- A=%.9g m, K=%.9g rad/m; source `%s`, run `%s`; requested beams [%d,%d], profile N [%d,%d].\n',v.A_m,v.K_radpm,v.config.source_geometry,v.config.run_type,v.beam_counts,v.profile_counts);
fprintf(fid,'- 阻塞信息：`%s`\n\n',v.blocked_reason);
fprintf(fid,'该状态表示计算成本/运行完成性阻塞，不把未完成结果解释为 reflection-model discrepancy；Stage 2--7 保持锁定。\n');
end

function as = local_stage0_as(cfg,x)
kx=(2*pi/cfg.xw_m)*[0:(cfg.nx/2-1),-cfg.nx/2:-1];
k0=2*pi*cfg.frequency_hz/cfg.c0_mps; kz=sqrt(complex(k0^2-kx.^2,0));
source=exp(-0.5*(x/cfg.sigma_src_m).^2);
field=-ifft(fft(source).*exp(1i*(cfg.z_tx_m+cfg.z_rx_m)*(kz-k0)));
as=struct('x_m',x,'field',field,'k0_rad_per_m',k0, ...
    'operator','one-step flat image angular spectrum');
end

function run = local_stage0_internal(cfg,x,beams,tag)
[zsort,ord] = sort(-x(:).','ascend'); invord=zeros(size(ord)); invord(ord)=1:numel(ord);
% For an exactly flat wall, three support points are mathematically
% identical to a dense profile and avoid turning the Stage-0 smoke into a
% profile-discretization benchmark.
s=[cfg.wall_support_m(1) 0 cfg.wall_support_m(2)].';
r=cfg.wall_r0_m*ones(size(s));
c=local_stage0_env_cfg(cfg,cfg.parametric_validation_exe,beams,zsort,tag);
c.wall_r0_m=cfg.wall_r0_m; c.wall_seed=cfg.wall_seed;
c.wall_profile_r_m=r; c.wall_profile_z_m=s;
raw=run_bellhop_internal_pm_wall_poc_vertical(c); run=raw;
run.pressure_raw=zeros(1,numel(x));
for ii=1:numel(x)
    run.pressure_raw(ii)=select_bellhop_shd_pressure_at_range_vertical(raw.data, ...
        cfg.mapped_receiver_range_m,zsort(ii),cfg.receiver_tolerance_m);
end
run.pressure_raw=run.pressure_raw(invord);
run.requested_receiver_depth_m=(-x(:)).';
run.receiver_coordinate_error_m=max(abs(double(raw.data.receiver_depth_m(:))-zsort(:)))*ones(1,numel(x));
end

function run = local_stage0_dedicated_flat(cfg,x,beams,tag)
[zsort,ord]=sort(-x(:).','ascend'); invord=zeros(size(ord)); invord(ord)=1:numel(ord);
c=local_stage0_env_cfg(cfg,cfg.flat_validation_exe,beams,zsort,tag);
c.wall_range_m=cfg.wall_r0_m;
raw=run_bellhop_internal_flat_wall_poc_vertical(c); run=raw;
run.pressure_raw=zeros(1,numel(x));
for ii=1:numel(x)
    run.pressure_raw(ii)=select_bellhop_shd_pressure_at_range_vertical(raw.data, ...
        cfg.mapped_receiver_range_m,zsort(ii),cfg.receiver_tolerance_m);
end
run.pressure_raw=run.pressure_raw(invord);
run.requested_receiver_depth_m=(-x(:)).';
run.receiver_coordinate_error_m=max(abs(double(raw.data.receiver_depth_m(:))-zsort(:)))*ones(1,numel(x));
end

function run = local_stage0_free(cfg,x,beams,tag)
[zsort,ord]=sort(-x(:).','ascend'); invord=zeros(size(ord)); invord(ord)=1:numel(ord);
c=local_stage0_env_cfg(cfg,cfg.official_exe,beams,zsort,tag);
raw=run_bellhop_unfolded_gaussian_vertical(c); run=raw;
run.pressure_raw=zeros(1,numel(x));
for ii=1:numel(x)
    run.pressure_raw(ii)=select_bellhop_shd_pressure_at_range_vertical(raw.data, ...
        cfg.mapped_receiver_range_m,zsort(ii),cfg.receiver_tolerance_m);
end
run.pressure_raw=run.pressure_raw(invord);
run.requested_receiver_depth_m=(-x(:)).';
run.receiver_coordinate_error_m=max(abs(double(raw.data.receiver_depth_m(:))-zsort(:)))*ones(1,numel(x));
end

function c = local_stage0_env_cfg(cfg,exe,beams,zsort,tag)
p=local_pattern(cfg);
c=struct('bellhop_exe',exe,'case_root',fullfile(cfg.output_dir,'bellhop',tag), ...
 'run_type',cfg.run_type,'source_geometry',cfg.source_geometry, ...
 'frequency_hz',cfg.frequency_hz,'c0_mps',cfg.c0_mps,'source_depth_m',0, ...
 'receiver_depths_m',zsort(:).','receiver_ranges_m',[cfg.guard_receiver_range_m cfg.mapped_receiver_range_m], ...
 'mapped_receiver_range_m',cfg.mapped_receiver_range_m, ...
 'beam_count',beams,'angle_limits_deg',cfg.angle_limits_deg,'step_m',cfg.bellhop_step_m, ...
 'domain_half_depth_m',cfg.domain_half_depth_m,'source_pattern_angles_deg',p.angles_deg, ...
 'source_pattern_level_db',p.level_db);
end

function p = local_stage0_convert(raw,phase_sign)
if phase_sign<0, p=conj(raw); else, p=raw; end
end

function fp=local_stage0_footprint(x,field,cfg)
w=abs(field(:)).^2; radii=sort(unique(abs(x(:))));
fp=struct('x_m',x(:).','energy_weights',w(:).','phase_floor_db',cfg.phase_floor_db, ...
    'axis_index',find(abs(x)==min(abs(x)),1));
for pp=[.95 .99]
    for ii=1:numel(radii)
        m=abs(x(:))<=radii(ii)+10*eps(max(1,radii(ii)));
        if sum(w(m))/max(sum(w),realmin)>=pp
            nm=sprintf('m%d',round(100*pp));
            fp.(nm)=struct('radius_m',radii(ii),'mask',m(:).', ...
                'energy_fraction',sum(w(m))/max(sum(w),realmin),'count',nnz(m));
            break
        end
    end
end
end

function sub=local_stage0_footprint_subset(fp,ix)
sub=fp; sub.x_m=fp.x_m(ix); sub.energy_weights=fp.energy_weights(ix);
sub.axis_index=find(ix==fp.axis_index,1);
for nm={'m95','m99'}
    n=nm{1}; sub.(n).mask=fp.(n).mask(ix); sub.(n).count=nnz(sub.(n).mask);
    sub.(n).energy_fraction=sum(sub.energy_weights.*sub.(n).mask)/max(sum(sub.energy_weights),realmin);
end
end

function m=local_stage0_metrics(ref,cand,fp,axis)
ref=ref(:).'; cand=cand(:).'; rn=ref/ref(axis); cn=cand/cand(axis);
d=rn-cn; ph=angle(rn.*conj(cn));
tl=20*log10(max(abs(rn),realmin)./max(abs(cn),realmin));
w=fp.energy_weights(:).'; m95=fp.m95.mask; m99=fp.m99.mask;
valid=m95 & abs(rn)>=10^(fp.phase_floor_db/20) & abs(cn)>=10^(fp.phase_floor_db/20);
ww=w(m99); S=sum(ww.*rn(m99).*conj(cn(m99)));
rho=abs(S)/sqrt(max(sum(ww.*abs(rn(m99)).^2)*sum(ww.*abs(cn(m99)).^2),realmin));
g=angle(S); alpha=sum(ww.*conj(cn(m99)).*rn(m99))/max(sum(ww.*abs(cn(m99)).^2),realmin);
m=struct('max_complex_full',max(abs(d)),'l2_m99',sqrt(sum(abs(d(m99)).^2)/max(sum(abs(rn(m99)).^2),realmin)), ...
 'phase_rms_m99',sqrt(sum(ww.*ph(m99).^2)/max(sum(ww),realmin)), ...
 'tl_rms_m99',sqrt(sum(ww.*tl(m99).^2)/max(sum(ww),realmin)), ...
 'phase_p95_m95',local_stage0_percentile(abs(ph(valid)),.95), ...
 'tl_p95_m95',local_stage0_percentile(abs(tl(valid)),.95),'rho_shape',rho, ...
 'global_phase_rad',g,'aligned_l2_m99',sqrt(sum(ww.*abs(rn(m99)-exp(1i*g)*cn(m99)).^2)/max(sum(ww.*abs(rn(m99)).^2),realmin)), ...
 'reference_axis',ref(axis),'candidate_axis',cand(axis),'phase_difference',ph, ...
 'tl_difference_db',tl,'alpha_ls',alpha);
end

function m=local_stage0_empty_metrics()
m=struct('l2_m99',NaN,'phase_rms_m99',NaN,'tl_rms_m99',NaN,'phase_p95_m95',NaN,'tl_p95_m95',NaN);
end

function v=local_stage0_percentile(a,p)
a=sort(a(isfinite(a))); if isempty(a),v=Inf;else,v=a(max(1,min(numel(a),ceil(p*numel(a)))));end
end

function v=local_stage0_outer5(a)
n=max(1,round(.05*numel(a))); ix=[1:n numel(a)-n+1:numel(a)];
v=sum(abs(a(ix)).^2)/max(sum(abs(a(:)).^2),realmin);
end

function [mi,hi,e]=local_stage0_shared_grid(x,h,tol)
mi=1:numel(x); hi=zeros(size(mi)); e=zeros(size(mi));
for ii=mi, [e(ii),hi(ii)]=min(abs(h-x(ii))); end
if any(e>tol), error('Half-dx grid mismatch.'); end
end

function ok=local_stage0_geometry_check(run,cfg)
d=run.diagnostics;
[~,ic]=min(abs(d.alpha_deg));
ok=all(isfinite(d{:,:}),'all') && max(abs(d.wall_residual))<=1e-6 && ...
    max(abs(d.kappa))<=1e-12 && max(abs(d.phase_delta-pi))<=1e-6 && ...
    min(d.min_post_dr)>0 && abs(d.tau_receiver_real(ic) - 103/cfg.c0_mps)<1e-6;
end

function local_stage0_write_csv(path,checks)
n=string(fieldnames(checks)); v=false(size(n));
for ii=1:numel(n),v(ii)=logical(checks.(char(n(ii))));end
writetable(table(n,v,'VariableNames',{'check','passed'}),path);
end

function local_stage0_write_report(path,v)
fid=fopen(path,'w','n','UTF-8'); if fid<0,error('Cannot write %s.',path);end
cl=onCleanup(@()fclose(fid)); c=v.config;
fprintf(fid,'# Stage 0 flat reflected-line closure\n\n状态：**%s**\n\n',ternary(v.passed,'PASS','FAIL'));
fprintf(fid,'- f/c: %.0f Hz / %.0f m/s; source `%s`, run `%s`; beams [%d,%d].\n', ...
    c.frequency_hz,c.c0_mps,c.source_geometry,c.run_type,c.beam_counts);
fprintf(fid,'- AS footprints: M95 %.6g m (%d), M99 %.6g m (%d).\n\n', ...
    v.footprint.m95.radius_m,v.footprint.m95.count,v.footprint.m99.radius_m,v.footprint.m99.count);
fprintf(fid,'| metric | L2(M99) | phase RMS(M99) |\n|---|---:|---:|\n');
names={'pe_as','flat_internal','chain_internal_free','dedicated_vs_internal','beam_convergence'};
for ii=1:numel(names)
    m=v.metrics.(names{ii}); fprintf(fid,'| %s | %.8g | %.8g |\n',names{ii},m.l2_m99,m.phase_rms_m99);
end
fprintf(fid,'\n## Checks\n\n'); n=fieldnames(v.checks);
for ii=1:numel(n),fprintf(fid,'- %s: %s\n',n{ii},ternary(v.checks.(n{ii}),'PASS','FAIL'));end
fprintf(fid,'\nNo subsequent stage is unlocked unless this report is reviewed and passes.\n');
end

function pattern = local_pattern(cfg)
angles = linspace(cfg.angle_limits_deg(1), cfg.angle_limits_deg(2), ...
    cfg.source_pattern_samples).';
k = 2*pi*cfg.frequency_hz/cfg.c0_mps; theta = angles*pi/180;
d = cos(theta) .* exp(-0.5*(k*cfg.sigma_src_m*sin(theta)).^2);
d = abs(d) / max(abs(d));
pattern = struct('source_geometry',cfg.source_geometry,'angles_deg',angles, ...
    'level_db',20*log10(max(d,10^(cfg.source_pattern_clip_db/20))), ...
    'formula','cos(theta)*exp(-0.5*(k*sigma*sin(theta))^2)', ...
    'fingerprint',sprintf('%s|cos(theta)*exp(-0.5*(k*sigma*sin(theta))^2)|N=%d|angle=[%.12g,%.12g]deg|f=%.12gHz|sigma=%.12gm', ...
        cfg.source_geometry,cfg.source_pattern_samples,cfg.angle_limits_deg(1), ...
        cfg.angle_limits_deg(2),cfg.frequency_hz,cfg.sigma_src_m));
end

function local_write_report(path, manifest, pattern, hashes)
fid = fopen(path, 'w', 'n', 'UTF-8');
if fid < 0, error('Cannot write %s.', path); end
cleanup = onCleanup(@() fclose(fid));
c = manifest.config;
fprintf(fid, '# PE--Bellhop controlled comparison P0 preparation\n\n');
fprintf(fid, '状态：**%s**\n\n', ternary(manifest.passed,'PASS','FAIL'));
fprintf(fid, '本阶段只建立 source-X/fingerprint-safe validation tree并核对依赖；没有调用 PE 或 Bellhop。\n\n');
fprintf(fid, '## Frozen configuration\n\n');
fprintf(fid, '- Bellhop：official AcousticsToolbox 2020 executable；run type `%s`；source geometry `%s`。\n', c.run_type, c.source_geometry);
fprintf(fid, '- f/c：%.0f Hz / %.0f m/s；PE `W=%.9g m`, `nx=%d`, `step=%.6g m`。\n', c.frequency_hz,c.c0_mps,c.xw_m,c.nx,c.pe_step_m);
fprintf(fid, '- receiver map：`z_BH''=-x_PE`; explicit ranges `[%g,%g] m`; sorted write/read-back permutation saved in MAT/JSON.\n', c.guard_receiver_range_m,c.mapped_receiver_range_m);
fprintf(fid, '- source-pattern fingerprint：`%s`。\n', pattern.fingerprint);
fprintf(fid, '- canonical PM coefficient source is deferred until Stage 6; expected SHA-256 `%s`.\n\n', c.canonical_coeff_sha256);
fprintf(fid, '## File hashes\n\n| file | SHA-256 |\n|---|---|\n');
names = fieldnames(hashes);
for ii = 1:numel(names), fprintf(fid, '| %s | `%s` |\n', names{ii}, hashes.(names{ii})); end
fprintf(fid, '\n## Checks\n\n');
names = fieldnames(manifest.checks);
for ii = 1:numel(names), fprintf(fid, '- %s: %s\n', names{ii}, ternary(manifest.checks.(names{ii}),'PASS','FAIL')); end
fprintf(fid, '\nP0 result is **%s**. Only after this result is reviewed may Stage 0 be unlocked explicitly.\n', ternary(manifest.passed,'READY','BLOCKED'));
end

function digest = local_sha256_file(path)
md = java.security.MessageDigest.getInstance('SHA-256');
fid = fopen(path, 'rb');
if fid < 0, error('Cannot open file for hashing: %s', path); end
cleanup = onCleanup(@() fclose(fid));
bytes = fread(fid, Inf, '*uint8');
clear cleanup
md.update(typecast(bytes, 'int8'));
digest = lower(reshape(dec2hex(typecast(md.digest(), 'uint8'), 2).', 1, []));
end

function out = ternary(condition, yes_value, no_value)
if condition, out = yes_value; else, out = no_value; end
end
