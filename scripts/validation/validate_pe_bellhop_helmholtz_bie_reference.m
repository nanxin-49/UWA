function result = validate_pe_bellhop_helmholtz_bie_reference(stage, overrides)
%VALIDATE_PE_BELLHOP_HELMHOLTZ_BIE_REFERENCE Independent Helmholtz BIE Goal.
%   Implements the gated R1--R5 independent-reference workflow. Each stage
%   requires the preceding authoritative MAT artifact to have passed.

if nargin < 1 || isempty(stage), stage = 'r1'; end
if nargin < 2 || isempty(overrides), overrides = struct(); end
if ~isstruct(overrides) || ~isscalar(overrides)
    error('overrides must be a scalar struct.');
end

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
setup_vertical_project();
switch lower(char(stage))
    case 'r1'
        cfg = local_r1_config(root,overrides);
        result = local_r1(cfg,root);
    case 'r2'
        cfg = local_r2_config(root,overrides);
        local_require_r1_pass(root);
        result = local_r2(cfg,root);
    case 'r3'
        cfg = local_r3_config(root,overrides);
        local_require_r2_pass(root);
        result = local_r3(cfg,root);
    case 'r4'
        cfg = local_r4_config(root,overrides);
        local_require_r3_pass(root);
        result = local_r4(cfg,root);
    case 'r5'
        cfg = local_r5_config(root,overrides);
        local_require_r4_pass(root);
        result = local_r4(cfg,root);
    otherwise
        error('Stage %s is not implemented in the R1--R5 workflow.',char(stage));
end
end

function cfg = local_r1_config(root,overrides)
cfg = struct( ...
    'frequency_hz',4000,'c0_mps',1500,'z_tx_m',100,'z_rx_m',3, ...
    'sigma_src_m',0.3,'xw_m',192.1875,'nx',984, ...
    'halfplane_h_m',-2.0,'window_plateau_ratio',0.7, ...
    'spatial_levels_ppw',[8 10 12],'spatial_half_width_m',60, ...
    'window_levels_m',[60 70 80],'window_sweep_ppw',8, ...
    'panel_order',8,'self_quadrature_order',64, ...
    'phase_floor_db',-40,'fail_on_check',false, ...
    'output_dir',fullfile(root,'results','validation', ...
        'pe_bellhop_helmholtz_bie_reference','R1_flat'), ...
    'report_file',fullfile(root,'reports', ...
        'pe_bellhop_helmholtz_bie_R1_flat_validation_report.md'));
names = fieldnames(overrides);
for ii = 1:numel(names)
    if ~isfield(cfg,names{ii}), error('Unknown R1 override %s.',names{ii}); end
    cfg.(names{ii}) = overrides.(names{ii});
end
if cfg.frequency_hz ~= 4000 || cfg.c0_mps ~= 1500 || ...
        cfg.z_tx_m ~= 100 || cfg.z_rx_m ~= 3 || cfg.sigma_src_m ~= 0.3
    error('R1 physical parameters are frozen by the Goal.');
end
if mod(cfg.nx,2) ~= 0, error('nx must be even.'); end
end

function cfg = local_r2_config(root,overrides)
cfg = struct( ...
    'frequency_hz',4000,'c0_mps',1500,'z_tx_m',100,'z_rx_m',3, ...
    'sigma_src_m',0.3,'xw_m',192.1875,'nx',984, ...
    'amplitude_m',0.01,'wavenumber_radpm',0.10, ...
    'surface_taper_inner_m',42,'surface_support_m',50, ...
    'halfplane_h_m',-2.0,'window_plateau_ratio',0.85, ...
    'spatial_levels_ppw',[8 10 12],'spatial_half_width_m',60, ...
    'quadrature_levels',[64 96],'quadrature_ppw',12, ...
    'window_levels_m',[60 70 80],'window_sweep_ppw',10, ...
    'panel_order',8,'self_quadrature_order',96, ...
    'phase_floor_db',-40,'fail_on_check',false, ...
    'output_dir',fullfile(root,'results','validation', ...
        'pe_bellhop_helmholtz_bie_reference','R2_weak_rough'), ...
    'report_file',fullfile(root,'reports', ...
        'pe_bellhop_helmholtz_bie_R2_weak_rough_convergence_report.md'));
names = fieldnames(overrides);
for ii = 1:numel(names)
    if ~isfield(cfg,names{ii}), error('Unknown R2 override %s.',names{ii}); end
    cfg.(names{ii}) = overrides.(names{ii});
end
if cfg.frequency_hz ~= 4000 || cfg.c0_mps ~= 1500 || ...
        cfg.z_tx_m ~= 100 || cfg.z_rx_m ~= 3 || cfg.sigma_src_m ~= 0.3 || ...
        cfg.amplitude_m ~= 0.01 || cfg.wavenumber_radpm ~= 0.10
    error('R2 physical parameters are frozen by the Goal.');
end
if numel(cfg.spatial_levels_ppw) < 3 || numel(cfg.window_levels_m) < 3 || ...
        numel(cfg.quadrature_levels) < 2
    error('R2 requires 3 spatial, 2 quadrature, and 3 window levels.');
end
if cfg.surface_taper_inner_m <= 0 || ...
        cfg.surface_support_m <= cfg.surface_taper_inner_m
    error('R2 physical surface taper/support is invalid.');
end
if any(cfg.window_plateau_ratio*cfg.window_levels_m <= cfg.surface_support_m) || ...
        cfg.window_plateau_ratio*cfg.spatial_half_width_m <= cfg.surface_support_m
    error('Every BIE w=1 plateau must strictly contain the physical surface support.');
end
end

function local_require_r1_pass(root)
path = fullfile(root,'results','validation','pe_bellhop_helmholtz_bie_reference', ...
    'R1_flat','R1_flat_validation.mat');
if ~exist(path,'file')
    error('R2 locked: authoritative R1 artifact is missing.');
end
s = load(path,'validation');
if ~isfield(s,'validation') || ~isfield(s.validation,'passed') || ~s.validation.passed
    error('R2 locked: authoritative R1 did not pass all hard gates.');
end
end

function cfg = local_r3_config(root,overrides)
cfg = struct( ...
    'frequency_hz',4000,'c0_mps',1500,'z_tx_m',100,'z_rx_m',3, ...
    'sigma_src_m',0.3,'xw_m',192.1875,'nx',984, ...
    'amplitude_m',0.01,'wavenumber_radpm',0.10, ...
    'surface_taper_inner_m',42,'surface_support_m',50, ...
    'wall_profile_support_m',[-80 80],'wall_profile_count',4097, ...
    'beam_count',10001,'bellhop_step_m',0.05, ...
    'receiver_tolerance_m',1e-6,'reuse_existing',true, ...
    'output_dir',fullfile(root,'results','validation', ...
        'pe_bellhop_helmholtz_bie_reference','R3_weak_three_way'), ...
    'report_file',fullfile(root,'reports', ...
        'pe_bellhop_helmholtz_bie_R3_weak_three_way_closure_report.md'));
names = fieldnames(overrides);
for ii = 1:numel(names)
    if ~isfield(cfg,names{ii}), error('Unknown R3 override %s.',names{ii}); end
    cfg.(names{ii}) = overrides.(names{ii});
end
if cfg.frequency_hz ~= 4000 || cfg.c0_mps ~= 1500 || ...
        cfg.z_tx_m ~= 100 || cfg.z_rx_m ~= 3 || cfg.sigma_src_m ~= 0.3 || ...
        cfg.amplitude_m ~= 0.01 || cfg.wavenumber_radpm ~= 0.10 || ...
        cfg.beam_count ~= 10001
    error('R3 physics and the 10,001-beam policy are frozen by the Goal.');
end
end

function local_require_r2_pass(root)
path = fullfile(root,'results','validation','pe_bellhop_helmholtz_bie_reference', ...
    'R2_weak_rough','R2_weak_rough_convergence.mat');
if ~exist(path,'file'), error('R3 locked: authoritative R2 artifact is missing.'); end
s = load(path,'validation');
if ~isfield(s,'validation') || ~isfield(s.validation,'passed') || ~s.validation.passed
    error('R3 locked: authoritative R2 did not pass all gates.');
end
end

function cfg = local_r4_config(root,overrides)
cfg = struct( ...
    'stage_id','R4_region_II','artifact_stem','R4_region_II', ...
    'frequency_hz',4000,'c0_mps',1500,'z_tx_m',100,'z_rx_m',3, ...
    'sigma_src_m',0.3,'xw_m',192.1875,'nx',984, ...
    'amplitude_m',0.05,'wavenumber_radpm',0.10, ...
    'surface_taper_inner_m',42,'surface_support_m',50, ...
    'halfplane_h_m',-2.0,'window_plateau_ratio',0.85, ...
    'spatial_levels_ppw',[8 10 12],'spatial_half_width_m',60, ...
    'window_levels_m',[60 70 80],'window_sweep_ppw',10, ...
    'panel_order',8,'self_quadrature_order',96, ...
    'phase_floor_db',-40,'wall_profile_support_m',[-80 80], ...
    'wall_profile_count',4097,'beam_count',10001,'bellhop_step_m',0.05, ...
    'receiver_tolerance_m',1e-6,'reuse_existing',true, ...
    'output_dir',fullfile(root,'results','validation', ...
        'pe_bellhop_helmholtz_bie_reference','R4_region_II'), ...
    'report_file',fullfile(root,'reports', ...
        'pe_bellhop_helmholtz_bie_R4_region_II_adjudication_report.md'));
names = fieldnames(overrides);
for ii = 1:numel(names)
    if ~isfield(cfg,names{ii}), error('Unknown R4 override %s.',names{ii}); end
    cfg.(names{ii}) = overrides.(names{ii});
end
if cfg.frequency_hz~=4000 || cfg.c0_mps~=1500 || cfg.amplitude_m~=0.05 || ...
        cfg.wavenumber_radpm~=0.10 || cfg.beam_count~=10001
    error('R4 physics and 10,001-beam policy are frozen by the Goal.');
end
if numel(cfg.spatial_levels_ppw)<3 || numel(cfg.window_levels_m)<3
    error('R4 requires three spatial and three window levels.');
end
end

function cfg = local_r5_config(root,overrides)
cfg = struct( ...
    'stage_id','R5_stronger_height','artifact_stem','R5_stronger_height', ...
    'frequency_hz',4000,'c0_mps',1500,'z_tx_m',100,'z_rx_m',3, ...
    'sigma_src_m',0.3,'xw_m',192.1875,'nx',984, ...
    'amplitude_m',0.20,'wavenumber_radpm',0.10, ...
    'surface_taper_inner_m',42,'surface_support_m',50, ...
    'halfplane_h_m',-2.0,'window_plateau_ratio',0.85, ...
    'spatial_levels_ppw',[8 10 12 14],'spatial_half_width_m',60, ...
    'window_levels_m',[60 70 80],'window_sweep_ppw',10, ...
    'panel_order',8,'self_quadrature_order',96, ...
    'boundary_check_panel_order',10,'boundary_check_self_quadrature_order',120, ...
    'boundary_check_ppw',12, ...
    'phase_floor_db',-40,'wall_profile_support_m',[-80 80], ...
    'wall_profile_count',4097,'beam_count',10001,'bellhop_step_m',0.05, ...
    'receiver_tolerance_m',1e-6,'reuse_existing',true, ...
    'output_dir',fullfile(root,'results','validation', ...
        'pe_bellhop_helmholtz_bie_reference','R5_stronger_height'), ...
    'report_file',fullfile(root,'reports', ...
        'pe_bellhop_helmholtz_bie_R5_stronger_height_report.md'));
names=fieldnames(overrides);
for ii=1:numel(names)
    if ~isfield(cfg,names{ii}),error('Unknown R5 override %s.',names{ii});end
    cfg.(names{ii})=overrides.(names{ii});
end
if cfg.frequency_hz~=4000 || cfg.c0_mps~=1500 || cfg.amplitude_m~=0.20 || ...
        cfg.wavenumber_radpm~=0.10 || cfg.beam_count~=10001
    error('R5 physics and 10,001-beam policy are frozen by the Goal.');
end
if numel(cfg.spatial_levels_ppw)<3 || numel(cfg.window_levels_m)<3
    error('R5 requires three spatial and three window levels.');
end
end

function local_require_r3_pass(root)
path = fullfile(root,'results','validation','pe_bellhop_helmholtz_bie_reference', ...
    'R3_weak_three_way','R3_weak_three_way_validation.mat');
if ~exist(path,'file'), error('R4 locked: authoritative R3 artifact is missing.'); end
s = load(path,'validation');
if ~s.validation.passed, error('R4 locked: R3 weak closure did not pass.'); end
end

function local_require_r4_pass(root)
path=fullfile(root,'results','validation','pe_bellhop_helmholtz_bie_reference', ...
    'R4_region_II','R4_region_II_validation.mat');
if ~exist(path,'file'),error('R5 locked: authoritative R4 artifact is missing.');end
s=load(path,'validation');
if ~s.validation.passed,error('R5 locked: R4 did not pass.');end
end

function result = local_r1(cfg,root)
if ~exist(cfg.output_dir,'dir'), mkdir(cfg.output_dir); end
[source,source_check] = local_discrete_gaussian_source(cfg);
x_rx = source.x_grid_m(:);
flat_ref = -source.evaluate(x_rx,cfg.z_rx_m+zeros(size(x_rx)),cfg.z_tx_m+cfg.z_rx_m,'path');
footprint = local_footprint(x_rx,flat_ref,cfg.phase_floor_db);
if numel(cfg.spatial_levels_ppw) < 3 || numel(cfg.window_levels_m) < 3
    error('R1 requires at least three spatial and three window levels.');
end
if any(cfg.window_plateau_ratio*cfg.window_levels_m <= footprint.radius99_m) || ...
        cfg.window_plateau_ratio*cfg.spatial_half_width_m <= footprint.radius99_m
    error('Every R1 w=1 plateau must strictly contain the AS M99 footprint.');
end

spatial = cell(numel(cfg.spatial_levels_ppw),1);
for ii = 1:numel(spatial)
    bcfg = local_bie_config(cfg,source,x_rx,cfg.spatial_half_width_m, ...
        cfg.spatial_levels_ppw(ii));
    spatial{ii} = solve_helmholtz_bie_halfplane_vertical(bcfg);
    spatial{ii}.metrics = local_metrics(spatial{ii}.receiver_field,flat_ref,footprint);
end

window = cell(numel(cfg.window_levels_m),1);
for ii = 1:numel(window)
    duplicate = find(cfg.window_levels_m(ii)==cfg.spatial_half_width_m & ...
        cfg.window_sweep_ppw==cfg.spatial_levels_ppw,1);
    if ~isempty(duplicate)
        window{ii} = spatial{duplicate};
    else
        bcfg = local_bie_config(cfg,source,x_rx,cfg.window_levels_m(ii), ...
            cfg.window_sweep_ppw);
        window{ii} = solve_helmholtz_bie_halfplane_vertical(bcfg);
        window{ii}.metrics = local_metrics(window{ii}.receiver_field,flat_ref,footprint);
    end
end

main = spatial{end};
spatial_convergence = local_convergence_table(spatial,footprint,'points_per_wavelength_actual');
window_convergence = local_convergence_table(window,footprint,'half_width_m');
spatial_last_change = spatial_convergence.successive_complex_l2_M99(end);
window_last_change = window_convergence.successive_complex_l2_M99(end);
checks = struct( ...
    'source_dft',source_check <= 1e-12, ...
    'linear_residual',main.linear_residual <= 1e-10, ...
    'boundary_residual',main.offgrid_boundary_residual <= 1e-8, ...
    'complex_l2',main.metrics.complex_l2_m99 <= 1e-6, ...
    'phase_rms',main.metrics.phase_rms_m99 <= 1e-6, ...
    'tl_rms_reported',isfinite(main.metrics.tl_rms_m95_db), ...
    'rho_shape_reported',isfinite(main.metrics.rho_shape), ...
    'spatial_convergence',spatial_last_change <= 1e-5, ...
    'window_convergence',window_last_change <= 1e-5, ...
    'finite',main.finite && all(isfinite(cell2mat(struct2cell(main.metrics)))));
checks.all = all(cell2mat(struct2cell(checks)));

cfg_saved = cfg;
validation = struct('schema_version','1.0.0','stage','R1_flat', ...
    'root',root,'config',cfg_saved,'source_fingerprint',source.fingerprint, ...
    'source_dft_grid_error',source_check,'receiver_x_m',x_rx, ...
    'flat_reference',flat_ref,'footprint',footprint,'spatial_runs',{spatial}, ...
    'window_runs',{window},'spatial_convergence',spatial_convergence, ...
    'window_convergence',window_convergence,'main',main,'checks',checks, ...
    'spatial_last_change_M99',spatial_last_change, ...
    'window_last_change_M99',window_last_change, ...
    'passed',checks.all);
mat_file = fullfile(cfg.output_dir,'R1_flat_validation.mat');
csv_file = fullfile(cfg.output_dir,'R1_flat_convergence.csv');
save(mat_file,'validation','-v7.3');
writetable([spatial_convergence;window_convergence],csv_file);
local_write_report(cfg.report_file,validation,mat_file,csv_file);
result = validation;
if cfg.fail_on_check && ~result.passed
    error('R1 flat BIE gate failed; rough stages remain locked.');
end
end

function result = local_r2(cfg,root)
if ~exist(cfg.output_dir,'dir'), mkdir(cfg.output_dir); end
[source,source_check] = local_discrete_gaussian_source(cfg);
x_rx = source.x_grid_m(:);
flat_ref = -source.evaluate(x_rx,cfg.z_rx_m+zeros(size(x_rx)), ...
    cfg.z_tx_m+cfg.z_rx_m,'path');
footprint = local_footprint(x_rx,flat_ref,cfg.phase_floor_db);
if cfg.window_plateau_ratio*min(cfg.window_levels_m) <= footprint.radius99_m
    error('Every R2 w=1 plateau must contain both surface support and AS M99.');
end

[eta_fn,eta_prime_fn,surface_meta] = local_benchmark_surface(cfg);
run_cache = struct('key',{},'value',{});

spatial = cell(numel(cfg.spatial_levels_ppw),1);
for ii = 1:numel(spatial)
    bcfg = local_bie_config(cfg,source,x_rx,cfg.spatial_half_width_m, ...
        cfg.spatial_levels_ppw(ii));
    bcfg.eta_fn = eta_fn;
    bcfg.eta_prime_fn = eta_prime_fn;
    bcfg.self_quadrature_order = cfg.self_quadrature_order;
    [spatial{ii},run_cache] = local_cached_bie_run(bcfg,run_cache);
end

quadrature = cell(numel(cfg.quadrature_levels),1);
for ii = 1:numel(quadrature)
    bcfg = local_bie_config(cfg,source,x_rx,cfg.spatial_half_width_m, ...
        cfg.quadrature_ppw);
    bcfg.eta_fn = eta_fn;
    bcfg.eta_prime_fn = eta_prime_fn;
    bcfg.self_quadrature_order = cfg.quadrature_levels(ii);
    [quadrature{ii},run_cache] = local_cached_bie_run(bcfg,run_cache);
end

window = cell(numel(cfg.window_levels_m),1);
for ii = 1:numel(window)
    bcfg = local_bie_config(cfg,source,x_rx,cfg.window_levels_m(ii), ...
        cfg.window_sweep_ppw);
    bcfg.eta_fn = eta_fn;
    bcfg.eta_prime_fn = eta_prime_fn;
    bcfg.self_quadrature_order = cfg.self_quadrature_order;
    [window{ii},run_cache] = local_cached_bie_run(bcfg,run_cache);
end

main = quadrature{end};
spatial_convergence = local_rough_convergence_table(spatial,main.receiver_field, ...
    footprint,'points_per_wavelength_actual');
quadrature_convergence = local_rough_convergence_table(quadrature,main.receiver_field, ...
    footprint,'self_quadrature_order');
window_convergence = local_rough_convergence_table(window,window{end}.receiver_field, ...
    footprint,'half_width_m');
u_spatial = spatial_convergence.successive_complex_l2_M99(end);
u_quadrature = quadrature_convergence.successive_complex_l2_M99(end);
u_window = window_convergence.successive_complex_l2_M99(end);
u_bie = max([u_spatial,u_quadrature,u_window]);
u_phase = max([spatial_convergence.successive_phase_rms_M99(end), ...
    quadrature_convergence.successive_phase_rms_M99(end), ...
    window_convergence.successive_phase_rms_M99(end)]);
u_tl = max([spatial_convergence.successive_tl_rms_M95_db(end), ...
    quadrature_convergence.successive_tl_rms_M95_db(end), ...
    window_convergence.successive_tl_rms_M95_db(end)]);
all_runs = [spatial;quadrature;window];
finite_runs = all(cellfun(@(q) q.finite,all_runs));
checks = struct( ...
    'source_dft',source_check <= 1e-12, ...
    'linear_residual',main.linear_residual <= 1e-10, ...
    'boundary_residual',main.offgrid_boundary_residual <= 1e-8, ...
    'spatial_convergence',u_spatial <= 1e-5, ...
    'quadrature_convergence',u_quadrature <= 1e-5, ...
    'window_convergence',u_window <= 1e-5, ...
    'finite',finite_runs && all(isfinite([u_bie,u_phase,u_tl])));
checks.all = all(cell2mat(struct2cell(checks)));

validation = struct('schema_version','1.0.0','stage','R2_weak_rough', ...
    'root',root,'config',cfg,'source_fingerprint',source.fingerprint, ...
    'source_dft_grid_error',source_check,'receiver_x_m',x_rx, ...
    'flat_reference',flat_ref,'footprint',footprint,'surface_meta',surface_meta, ...
    'spatial_runs',{spatial},'quadrature_runs',{quadrature}, ...
    'window_runs',{window},'spatial_convergence',spatial_convergence, ...
    'quadrature_convergence',quadrature_convergence, ...
    'window_convergence',window_convergence,'main',main,'checks',checks, ...
    'U_BIE_complex_l2_M99',u_bie,'U_BIE_phase_rms_M99_rad',u_phase, ...
    'U_BIE_tl_rms_M95_db',u_tl,'passed',checks.all);
mat_file = fullfile(cfg.output_dir,'R2_weak_rough_convergence.mat');
csv_file = fullfile(cfg.output_dir,'R2_weak_rough_convergence.csv');
save(mat_file,'validation','-v7.3');
writetable([spatial_convergence;quadrature_convergence;window_convergence],csv_file);
local_write_r2_report(cfg.report_file,validation,mat_file,csv_file);
result = validation;
if cfg.fail_on_check && ~result.passed
    error('R2 weak-rough BIE convergence gate failed; R3 remains locked.');
end
end

function result = local_r3(cfg,root)
if ~exist(cfg.output_dir,'dir'), mkdir(cfg.output_dir); end
r2_file = fullfile(root,'results','validation','pe_bellhop_helmholtz_bie_reference', ...
    'R2_weak_rough','R2_weak_rough_convergence.mat');
t = load(r2_file,'validation'); r2 = t.validation;
stage0_file = fullfile(root,'results','validation','pe_bellhop_controlled_comparison', ...
    'stage0','stage0_validation.mat');
stage1_file = fullfile(root,'results','validation','pe_bellhop_controlled_comparison', ...
    'stage1_convention_fixed','stage1_validation.mat');
if ~exist(stage0_file,'file') || ~exist(stage1_file,'file')
    error('R3 requires the authoritative controlled-comparison Stage 0/1 artifacts.');
end
t = load(stage0_file,'validation'); s0 = t.validation;
t = load(stage1_file,'validation'); s1 = t.validation;
if ~s0.passed || ~s1.passed, error('Controlled-comparison prerequisites did not pass.'); end

x = s0.x_m(:);
if max(abs(x-r2.receiver_x_m(:))) > cfg.receiver_tolerance_m
    error('R3 PE/Bellhop/BIE receiver grids are not identical.');
end
[eta_fn,~,surface_meta] = local_benchmark_surface(cfg);
eta_x = eta_fn(x);
pe_cfg = struct('frequency_hz',cfg.frequency_hz,'c0_mps',cfg.c0_mps, ...
    'xw_m',cfg.xw_m,'nx',cfg.nx,'z_tx_m',cfg.z_tx_m,'z_rx_m',cfg.z_rx_m, ...
    'sigma_src_m',cfg.sigma_src_m,'surface_elevation_x_m',eta_x(:).', ...
    'surface_reflect_coeff',-1,'step_m',s0.config.pe_step_m,'x_rx_m',0);
pe = run_pe_1d_surface_reflection_validation(pe_cfg);

s = linspace(cfg.wall_profile_support_m(1),cfg.wall_profile_support_m(2), ...
    cfg.wall_profile_count).';
r_wall = s0.config.wall_r0_m-eta_fn(s);
[zsort,ord] = sort(-x(:).','ascend');
invord = zeros(size(ord)); invord(ord) = 1:numel(ord);
c = struct('bellhop_exe',s0.config.parametric_validation_exe, ...
    'case_root',fullfile(cfg.output_dir,'bellhop','weak_tapered_benchmark'), ...
    'run_type',s0.config.run_type,'source_geometry',s0.config.source_geometry, ...
    'frequency_hz',cfg.frequency_hz,'c0_mps',cfg.c0_mps,'source_depth_m',0, ...
    'receiver_depths_m',zsort(:).','receiver_ranges_m', ...
    [s0.config.guard_receiver_range_m s0.config.mapped_receiver_range_m], ...
    'mapped_receiver_range_m',s0.config.mapped_receiver_range_m, ...
    'beam_count',cfg.beam_count,'angle_limits_deg',s0.config.angle_limits_deg, ...
    'step_m',cfg.bellhop_step_m,'domain_half_depth_m',s0.config.domain_half_depth_m, ...
    'source_pattern_angles_deg',s0.source_pattern.angles_deg, ...
    'source_pattern_level_db',s0.source_pattern.level_db, ...
    'wall_r0_m',s0.config.wall_r0_m,'wall_seed',s0.config.wall_seed, ...
    'wall_profile_r_m',r_wall,'wall_profile_z_m',s);
if cfg.reuse_existing, c.reuse_existing = true; end
bh = run_bellhop_internal_pm_wall_poc_vertical(c);
pressure_raw = complex(zeros(1,numel(x)));
for ii = 1:numel(x)
    pressure_raw(ii) = select_bellhop_shd_pressure_at_range_vertical( ...
        bh.data,s0.config.mapped_receiver_range_m,zsort(ii),cfg.receiver_tolerance_m);
end
pressure_raw = pressure_raw(invord).';
if s0.config.phase_sign < 0, bh_field = conj(pressure_raw); else, bh_field = pressure_raw; end

fp = r2.footprint;
ratio_mask = fp.m99;
G_PE = complex(zeros(size(x)));
G_BH_native = complex(zeros(size(x)));
G_BIE_native = complex(zeros(size(x)));
pe_rough = pe.reflected_field(:); pe_flat = s0.pe.reflected_field(:);
bh_rough = bh_field(:); bh_flat = s0.bellhop_internal_high.field(:);
bie_rough = r2.main.receiver_field(:); bie_flat = r2.flat_reference(:);
G_PE(ratio_mask) = pe_rough(ratio_mask)./pe_flat(ratio_mask);
G_BH_native(ratio_mask) = bh_rough(ratio_mask)./bh_flat(ratio_mask);
G_BIE_native(ratio_mask) = bie_rough(ratio_mask)./bie_flat(ratio_mask);
G_BH = conj(G_BH_native);
G_BIE = conj(G_BIE_native);
metrics = struct( ...
    'PE_BIE',local_three_way_metrics(G_PE,G_BIE,fp), ...
    'BH_BIE',local_three_way_metrics(G_BH,G_BIE,fp), ...
    'PE_BH',local_three_way_metrics(G_PE,G_BH,fp), ...
    'BH_native_BIE_native',local_three_way_metrics(G_BH_native,G_BIE_native,fp), ...
    'PE_BIE_native_diagnostic',local_three_way_metrics(G_PE,G_BIE_native,fp));
thresholds = struct('E',s1.floor.T_E+r2.U_BIE_complex_l2_M99, ...
    'phase_rad',s1.floor.T_phi+r2.U_BIE_phase_rms_M99_rad, ...
    'tl_db',s1.floor.T_TL+r2.U_BIE_tl_rms_M95_db,'rho_shape',0.9995);
geometry = local_r3_geometry(bh,s0.config);
checks = struct( ...
    'receiver_grid',max(abs(x-r2.receiver_x_m(:)))<=cfg.receiver_tolerance_m, ...
    'bellhop_geometry',geometry.all, ...
    'PE_BIE_E',metrics.PE_BIE.E_G<=thresholds.E, ...
    'PE_BIE_phase',metrics.PE_BIE.phase_rms_rad<=thresholds.phase_rad, ...
    'PE_BIE_TL',metrics.PE_BIE.tl_rms_db<=thresholds.tl_db, ...
    'PE_BIE_shape',metrics.PE_BIE.rho_shape>=thresholds.rho_shape, ...
    'BH_BIE_E',metrics.BH_BIE.E_G<=thresholds.E, ...
    'BH_BIE_phase',metrics.BH_BIE.phase_rms_rad<=thresholds.phase_rad, ...
    'BH_BIE_TL',metrics.BH_BIE.tl_rms_db<=thresholds.tl_db, ...
    'BH_BIE_shape',metrics.BH_BIE.rho_shape>=thresholds.rho_shape, ...
    'finite',all(isfinite([G_PE(ratio_mask);G_BH(ratio_mask);G_BIE(ratio_mask)])));
checks.all = all(cell2mat(struct2cell(checks)));
validation = struct('schema_version','1.0.0','stage','R3_weak_three_way', ...
    'config',cfg,'surface_meta',surface_meta,'x_m',x,'eta_x_m',eta_x, ...
    'PE',pe,'Bellhop',bh,'BIE_R2_file',r2_file,'G_PE',G_PE, ...
    'G_BH_native',G_BH_native,'G_BH',G_BH, ...
    'G_BIE_native',G_BIE_native,'G_BIE',G_BIE, ...
    'comparison_convention',['G_BH=conj(G_BH_native) and ', ...
    'G_BIE=conj(G_BIE_native), fixed mapping from exp(-iwt) native ', ...
    'Helmholtz convention to the PE comparison convention; no fitted scalar'], ...
    'metrics',metrics,'thresholds',thresholds,'geometry',geometry, ...
    'checks',checks,'passed',checks.all);
mat_file = fullfile(cfg.output_dir,'R3_weak_three_way_validation.mat');
save(mat_file,'validation','-v7.3');
local_write_r3_report(cfg.report_file,validation,mat_file);
result = validation;
if ~result.passed
    error('R3 weak three-way closure failed; R4/R5 remain locked.');
end
end

function result = local_r4(cfg,root)
if ~exist(cfg.output_dir,'dir'), mkdir(cfg.output_dir); end
[source,source_check] = local_discrete_gaussian_source(cfg);
x = source.x_grid_m(:);
flat_ref = -source.evaluate(x,cfg.z_rx_m+zeros(size(x)), ...
    cfg.z_tx_m+cfg.z_rx_m,'path');
fp = local_footprint(x,flat_ref,cfg.phase_floor_db);
[eta_fn,eta_prime_fn,surface_meta] = local_benchmark_surface(cfg);
cache = struct('key',{},'value',{});
existing_file=fullfile(cfg.output_dir,[cfg.artifact_stem '_validation.mat']);
if exist(existing_file,'file')
    previous=load(existing_file,'validation');
    if isfield(previous,'validation') && ...
            previous.validation.config.amplitude_m==cfg.amplitude_m && ...
            previous.validation.config.wavenumber_radpm==cfg.wavenumber_radpm
        cache=local_seed_run_cache(cache,previous.validation);
    end
end
spatial = cell(numel(cfg.spatial_levels_ppw),1);
for ii = 1:numel(spatial)
    bcfg = local_bie_config(cfg,source,x,cfg.spatial_half_width_m, ...
        cfg.spatial_levels_ppw(ii));
    bcfg.eta_fn=eta_fn; bcfg.eta_prime_fn=eta_prime_fn;
    [spatial{ii},cache] = local_cached_bie_run(bcfg,cache);
end
window = cell(numel(cfg.window_levels_m),1);
for ii = 1:numel(window)
    bcfg = local_bie_config(cfg,source,x,cfg.window_levels_m(ii), ...
        cfg.window_sweep_ppw);
    bcfg.eta_fn=eta_fn; bcfg.eta_prime_fn=eta_prime_fn;
    [window{ii},cache] = local_cached_bie_run(bcfg,cache);
end
main = spatial{end};
boundary_verification=[];
boundary_pair_l2=0;
if isfield(cfg,'boundary_check_panel_order')
    bcfg=local_bie_config(cfg,source,x,cfg.spatial_half_width_m, ...
        cfg.boundary_check_ppw);
    bcfg.panel_order=cfg.boundary_check_panel_order;
    bcfg.self_quadrature_order=cfg.boundary_check_self_quadrature_order;
    bcfg.eta_fn=eta_fn;bcfg.eta_prime_fn=eta_prime_fn;
    [boundary_verification,cache]=local_cached_bie_run(bcfg,cache); %#ok<ASGLU>
    boundary_pair_l2=local_metrics(boundary_verification.receiver_field, ...
        main.receiver_field,fp).complex_l2_m99;
end
spatial_convergence = local_rough_convergence_table(spatial,main.receiver_field, ...
    fp,'points_per_wavelength_actual');
window_convergence = local_rough_convergence_table(window,window{end}.receiver_field, ...
    fp,'half_width_m');
t = load(fullfile(root,'results','validation','pe_bellhop_helmholtz_bie_reference', ...
    'R2_weak_rough','R2_weak_rough_convergence.mat'),'validation');
r2 = t.validation;
u_case = max([r2.U_BIE_complex_l2_M99, ...
    spatial_convergence.successive_complex_l2_M99(end), ...
    window_convergence.successive_complex_l2_M99(end),boundary_pair_l2]);
[fields,pe,bh,geometry] = local_run_pe_bh_case(cfg,root,eta_fn,x);
mask = fp.m99;
G_BIE_native = complex(zeros(size(x)));
G_BIE_native(mask) = main.receiver_field(mask)./flat_ref(mask);
G_BIE = conj(G_BIE_native);
metrics = struct('PE_BIE',local_three_way_metrics(fields.G_PE,G_BIE,fp), ...
    'BH_BIE',local_three_way_metrics(fields.G_BH,G_BIE,fp), ...
    'PE_BH',local_three_way_metrics(fields.G_PE,fields.G_BH,fp), ...
    'BH_native_BIE_native',local_three_way_metrics( ...
    fields.G_BH_native,G_BIE_native,fp));

all_bie = [spatial;window];
kind = [repmat("spatial",numel(spatial),1);repmat("window",numel(window),1)];
parameter = [cfg.spatial_levels_ppw(:);cfg.window_levels_m(:)];
E_PE = zeros(numel(all_bie),1); E_BH=E_PE; separation=E_PE;
winner = strings(numel(all_bie),1);
for ii = 1:numel(all_bie)
    g_native = complex(zeros(size(x)));
    g_native(mask) = all_bie{ii}.receiver_field(mask)./flat_ref(mask);
    g = conj(g_native);
    E_PE(ii) = local_three_way_metrics(fields.G_PE,g,fp).E_G;
    E_BH(ii) = local_three_way_metrics(fields.G_BH,g,fp).E_G;
    separation(ii) = abs(E_PE(ii)-E_BH(ii));
    if E_PE(ii)<E_BH(ii), winner(ii)="PE"; else, winner(ii)="BELLHOP"; end
end
stability = table(kind,parameter,E_PE,E_BH,separation,winner);
main_separation = abs(metrics.PE_BIE.E_G-metrics.BH_BIE.E_G);
main_winner = winner(numel(spatial));
ranking_stable = all(winner==main_winner);
discriminating = main_separation>=5*u_case;
if isempty(boundary_verification)
    boundary_limit=1e-8;
    boundary_ok=main.offgrid_boundary_residual<=1e-8;
else
    boundary_limit=max(1e-8,u_case);
    boundary_ok=boundary_verification.offgrid_boundary_residual<=boundary_limit && ...
        boundary_pair_l2<=u_case;
end
checks = struct('source_dft',source_check<=1e-12, ...
    'boundary_residual',boundary_ok, ...
    'spatial_convergence',spatial_convergence.successive_complex_l2_M99(end)<=1e-5, ...
    'window_convergence',window_convergence.successive_complex_l2_M99(end)<=1e-5, ...
    'bellhop_geometry',geometry.all,'ranking_stable',ranking_stable, ...
    'five_U_separation',discriminating, ...
    'finite',all(isfinite([fields.G_PE(mask);fields.G_BH(mask);G_BIE(mask)])));
checks.all = all(cell2mat(struct2cell(checks)));
if checks.all && main_winner=="PE"
    conclusion = 'PE_CLOSER_TO_HELMHOLTZ_REFERENCE';
elseif checks.all && main_winner=="BELLHOP"
    conclusion = 'BELLHOP_CLOSER_TO_HELMHOLTZ_REFERENCE';
else
    conclusion = 'REFERENCE_NOT_YET_DISCRIMINATING';
end
validation = struct('schema_version','1.0.0','stage',cfg.stage_id, ...
    'config',cfg,'surface_meta',surface_meta,'x_m',x,'PE',pe,'Bellhop',bh, ...
    'G_PE',fields.G_PE,'G_BH_native',fields.G_BH_native,'G_BH',fields.G_BH, ...
    'G_BIE_native',G_BIE_native,'G_BIE',G_BIE,'metrics',metrics, ...
    'geometry',geometry,'spatial_runs',{spatial},'window_runs',{window}, ...
    'boundary_verification',boundary_verification, ...
    'boundary_verification_field_l2',boundary_pair_l2, ...
    'boundary_acceptance_limit',boundary_limit, ...
    'spatial_convergence',spatial_convergence,'window_convergence',window_convergence, ...
    'stability',stability,'U_BIE',u_case,'separation',main_separation, ...
    'five_U_threshold',5*u_case,'winner',main_winner,'checks',checks, ...
    'conclusion',conclusion,'passed',checks.all);
mat_file=fullfile(cfg.output_dir,[cfg.artifact_stem '_validation.mat']);
csv_file=fullfile(cfg.output_dir,[cfg.artifact_stem '_refinement_stability.csv']);
save(mat_file,'validation','-v7.3'); writetable(stability,csv_file);
local_write_r4_report(cfg.report_file,validation,mat_file,csv_file);
result=validation;
if ~result.passed, error('%s is not numerically stable/discriminating.',cfg.stage_id); end
end

function [fields,pe,bh,geometry] = local_run_pe_bh_case(cfg,root,eta_fn,x)
t=load(fullfile(root,'results','validation','pe_bellhop_controlled_comparison', ...
    'stage0','stage0_validation.mat'),'validation'); s0=t.validation;
if ~s0.passed, error('Controlled Stage 0 prerequisite failed.'); end
eta_x=eta_fn(x);
pe_cfg=struct('frequency_hz',cfg.frequency_hz,'c0_mps',cfg.c0_mps, ...
    'xw_m',cfg.xw_m,'nx',cfg.nx,'z_tx_m',cfg.z_tx_m,'z_rx_m',cfg.z_rx_m, ...
    'sigma_src_m',cfg.sigma_src_m,'surface_elevation_x_m',eta_x(:).', ...
    'surface_reflect_coeff',-1,'step_m',s0.config.pe_step_m,'x_rx_m',0);
pe=run_pe_1d_surface_reflection_validation(pe_cfg);
s=linspace(cfg.wall_profile_support_m(1),cfg.wall_profile_support_m(2), ...
    cfg.wall_profile_count).'; r_wall=s0.config.wall_r0_m-eta_fn(s);
[zsort,ord]=sort(-x(:).','ascend'); invord=zeros(size(ord));
invord(ord)=1:numel(ord);
case_tag=strrep(sprintf('A%.4g_K%.4g',cfg.amplitude_m, ...
    cfg.wavenumber_radpm),'.','p');
c=struct('bellhop_exe',s0.config.parametric_validation_exe, ...
    'case_root',fullfile(cfg.output_dir,'bellhop',case_tag), ...
    'run_type',s0.config.run_type,'source_geometry',s0.config.source_geometry, ...
    'frequency_hz',cfg.frequency_hz,'c0_mps',cfg.c0_mps,'source_depth_m',0, ...
    'receiver_depths_m',zsort(:).','receiver_ranges_m', ...
    [s0.config.guard_receiver_range_m s0.config.mapped_receiver_range_m], ...
    'mapped_receiver_range_m',s0.config.mapped_receiver_range_m, ...
    'beam_count',cfg.beam_count,'angle_limits_deg',s0.config.angle_limits_deg, ...
    'step_m',cfg.bellhop_step_m,'domain_half_depth_m',s0.config.domain_half_depth_m, ...
    'source_pattern_angles_deg',s0.source_pattern.angles_deg, ...
    'source_pattern_level_db',s0.source_pattern.level_db, ...
    'wall_r0_m',s0.config.wall_r0_m,'wall_seed',s0.config.wall_seed, ...
    'wall_profile_r_m',r_wall,'wall_profile_z_m',s);
if cfg.reuse_existing,c.reuse_existing=true;end
bh=run_bellhop_internal_pm_wall_poc_vertical(c);
raw=complex(zeros(1,numel(x)));
for ii=1:numel(x)
    raw(ii)=select_bellhop_shd_pressure_at_range_vertical(bh.data, ...
        s0.config.mapped_receiver_range_m,zsort(ii),cfg.receiver_tolerance_m);
end
raw=raw(invord).'; if s0.config.phase_sign<0,bh_field=conj(raw);else,bh_field=raw;end
mask=abs(s0.pe.reflected_field(:))>0 & abs(s0.bellhop_internal_high.field(:))>0;
G_PE=complex(zeros(size(x)));G_BH_native=G_PE;
pe_rough=pe.reflected_field(:);pe_flat=s0.pe.reflected_field(:);
bh_flat=s0.bellhop_internal_high.field(:);
G_PE(mask)=pe_rough(mask)./pe_flat(mask);
G_BH_native(mask)=bh_field(mask)./bh_flat(mask);
fields=struct('G_PE',G_PE,'G_BH_native',G_BH_native,'G_BH',conj(G_BH_native));
geometry=local_r3_geometry(bh,s0.config);
end

function m = local_three_way_metrics(a,b,fp)
a = a(:); b = b(:); mask = fp.m99; w = fp.energy_weights(mask);
w = w/max(sum(w),realmin);
phase = angle(a.*conj(b));
tl = 20*log10(max(abs(a),realmin)./max(abs(b),realmin));
S = sum(w.*a(mask).*conj(b(mask)));
den = sqrt(max(sum(w.*abs(a(mask)).^2)*sum(w.*abs(b(mask)).^2),realmin));
c = S/den; phi0 = angle(S);
m = struct('E_G',sqrt(sum(w.*abs(a(mask)-b(mask)).^2)/ ...
    max(sum(w.*abs(b(mask)).^2),realmin)), ...
    'phase_rms_rad',sqrt(sum(w.*phase(mask).^2)), ...
    'tl_rms_db',sqrt(sum(w.*tl(mask).^2)), ...
    'rho_raw',real(c),'rho_shape',abs(c),'phi0_rad',phi0, ...
    'E_aligned',sqrt(sum(w.*abs(a(mask)-exp(1i*phi0)*b(mask)).^2)/ ...
    max(sum(w.*abs(a(mask)).^2),realmin)));
end

function g = local_r3_geometry(run,cfg)
d = run.diagnostics;
g = struct('wall_residual_max_m',max(abs(d.wall_residual)), ...
    'phase_jump_error_max_rad',max(abs(d.phase_delta-pi)), ...
    'p_rotation_error_max',max(abs(d.p_rot_error)), ...
    'q_rotation_error_max',max(abs(d.q_rot_error)), ...
    'min_post_dr_m',min(d.min_post_dr));
g.all = all(isfinite(d{:,:}),'all') && ...
    g.wall_residual_max_m<=1e-9 && g.phase_jump_error_max_rad<=1e-10 && ...
    g.p_rotation_error_max<=1e-12 && g.q_rotation_error_max<=1e-12 && ...
    g.min_post_dr_m>0 && cfg.mapped_receiver_range_m==103;
end

function [eta_fn,eta_prime_fn,meta] = local_benchmark_surface(cfg)
A = cfg.amplitude_m;
K = cfg.wavenumber_radpm;
inner = cfg.surface_taper_inner_m;
outer = cfg.surface_support_m;
eta_fn = @(x) A*sin(K*x).*local_surface_taper(x,inner,outer,false);
eta_prime_fn = @(x) A*(K*cos(K*x).*local_surface_taper(x,inner,outer,false) + ...
    sin(K*x).*local_surface_taper(x,inner,outer,true));
probe = linspace(-outer,outer,20001).';
eta = eta_fn(probe);
deta = eta_prime_fn(probe);
meta = struct('definition','A*sin(K*x)*chi(x); fixed C2 quintic taper', ...
    'amplitude_m',A,'wavenumber_radpm',K,'taper_inner_m',inner, ...
    'support_m',outer,'max_abs_height_m',max(abs(eta)), ...
    'max_abs_slope',max(abs(deta)));
end

function y = local_surface_taper(x,inner,outer,derivative)
r = abs(x);
y = zeros(size(x));
inside = r <= inner;
middle = r > inner & r < outer;
if ~derivative
    y(inside) = 1;
    t = (r(middle)-inner)/(outer-inner);
    y(middle) = 1-10*t.^3+15*t.^4-6*t.^5;
else
    t = (r(middle)-inner)/(outer-inner);
    dydr = (-30*t.^2+60*t.^3-30*t.^4)/(outer-inner);
    y(middle) = dydr.*sign(x(middle));
end
end

function [out,cache] = local_cached_bie_run(cfg,cache)
key = sprintf('L=%.12g|ppw=%.12g|p=%d|q=%d',cfg.half_width_m, ...
    cfg.points_per_wavelength,cfg.panel_order,cfg.self_quadrature_order);
for ii = 1:numel(cache)
    if strcmp(cache(ii).key,key)
        out = cache(ii).value;
        return
    end
end
out = solve_helmholtz_bie_halfplane_vertical(cfg);
cache(end+1) = struct('key',key,'value',out); %#ok<AGROW>
end

function cache=local_seed_run_cache(cache,validation)
lists={};
if isfield(validation,'spatial_runs'),lists{end+1}=validation.spatial_runs;end %#ok<AGROW>
if isfield(validation,'window_runs'),lists{end+1}=validation.window_runs;end %#ok<AGROW>
if isfield(validation,'boundary_verification') && ...
        ~isempty(validation.boundary_verification)
    lists{end+1}={validation.boundary_verification}; %#ok<AGROW>
end
for ll=1:numel(lists)
    runs=lists{ll};
    for ii=1:numel(runs)
        q=runs{ii};
        key=sprintf('L=%.12g|ppw=%.12g|p=%d|q=%d',q.half_width_m, ...
            q.points_per_wavelength_requested,q.panel_order,q.self_quadrature_order);
        if ~any(arrayfun(@(entry)strcmp(entry.key,key),cache))
            cache(end+1)=struct('key',key,'value',q); %#ok<AGROW>
        end
    end
end
end

function bcfg = local_bie_config(cfg,source,x_rx,L,ppw)
bcfg = struct('frequency_hz',cfg.frequency_hz,'c0_mps',cfg.c0_mps, ...
    'half_width_m',L,'points_per_wavelength',ppw, ...
    'panel_order',cfg.panel_order,'self_quadrature_order',cfg.self_quadrature_order, ...
    'halfplane_h_m',cfg.halfplane_h_m, ...
    'window_plateau_ratio',cfg.window_plateau_ratio, ...
    'eta_fn',@(x) zeros(size(x)),'eta_prime_fn',@(x) zeros(size(x)), ...
    'incident_fn',@(x,z) source.evaluate(x,z,cfg.z_tx_m,'coordinate'), ...
    'receiver_x_m',x_rx,'receiver_z_m',cfg.z_rx_m);
end

function [source,grid_error] = local_discrete_gaussian_source(cfg)
dx = cfg.xw_m/cfg.nx;
x = (-cfg.nx/2:cfg.nx/2-1).'*dx;
kx = (2*pi/cfg.xw_m)*[0:(cfg.nx/2-1),-cfg.nx/2:-1];
k = 2*pi*cfg.frequency_hz/cfg.c0_mps;
kz = sqrt(complex(k^2-kx.^2,0));
initial = exp(-0.5*(x/cfg.sigma_src_m).^2);
coeff = fft(initial).';
x0 = x(1);
evaluate = @(xq,zq,path,mode) local_evaluate_source(xq,zq,path,mode, ...
    coeff,kx,kz,x0,cfg.nx);
grid_reconstructed = evaluate(x,zeros(size(x)),0,'path');
grid_error = norm(grid_reconstructed-initial)/norm(initial);
source = struct('x_grid_m',x,'kx_radpm',kx,'kz_radpm',kz, ...
    'coefficients',coeff,'initial_field',initial,'evaluate',evaluate, ...
    'fingerprint',sprintf(['discrete-periodic-Gaussian|N=%d|W=%.12g|', ...
        'sigma=%.12g|f=%.12g|c=%.12g|exp(-iwt)'],cfg.nx,cfg.xw_m, ...
        cfg.sigma_src_m,cfg.frequency_hz,cfg.c0_mps));
end

function u = local_evaluate_source(xq,zq,path,mode,coeff,kx,kz,x0,n)
xq = xq(:); zq = zq(:);
if isscalar(zq), zq = zq+zeros(size(xq)); end
u = complex(zeros(size(xq)));
block = 256;
for first = 1:block:numel(xq)
    rows = first:min(first+block-1,numel(xq));
    phase_x = exp(1i*(xq(rows)-x0)*kx);
    if strcmp(mode,'path')
        phase_z = exp(1i*path*kz);
        u(rows) = phase_x*(coeff.*phase_z).'/n;
    elseif strcmp(mode,'coordinate')
        phase_z = exp(-1i*(zq(rows)-path)*kz);
        u(rows) = sum(phase_x.*phase_z.*coeff,2)/n;
    else
        error('Unknown source evaluation mode %s.',mode);
    end
end
end

function fp = local_footprint(x,field,phase_floor_db)
energy = abs(field).^2;
[~,axis_index] = min(abs(x));
radius = abs(x-x(axis_index));
[radius_sorted,ord] = sort(radius);
cumulative = cumsum(energy(ord))/sum(energy);
r95 = radius_sorted(find(cumulative>=0.95,1));
r99 = radius_sorted(find(cumulative>=0.99,1));
m95 = radius<=r95;
m99 = radius<=r99;
phase_mask = m95 & abs(field)>=abs(field(axis_index))*10^(phase_floor_db/20);
weights = energy/max(sum(energy(m99)),realmin);
fp = struct('axis_index',axis_index,'radius95_m',r95,'radius99_m',r99, ...
    'm95',m95,'m99',m99,'phase_mask',phase_mask,'energy_weights',weights);
end

function m = local_metrics(test,ref,fp)
test = test(:); ref = ref(:); mask = fp.m99;
w = fp.energy_weights(mask); w = w/max(sum(w),realmin);
delta = test(mask)-ref(mask);
complex_l2 = norm(delta)/max(norm(ref(mask)),realmin);
phase = angle(test.*conj(ref));
phase_rms = sqrt(sum(w.*phase(mask).^2));
pm = fp.phase_mask;
tl = 20*log10(max(abs(test(pm)),realmin)./max(abs(ref(pm)),realmin));
tl_rms = sqrt(mean(tl.^2));
rho = abs(sum(test(mask).*conj(ref(mask))))/ ...
    sqrt(max(sum(abs(test(mask)).^2)*sum(abs(ref(mask)).^2),realmin));
m = struct('complex_l2_m99',complex_l2,'phase_rms_m99',phase_rms, ...
    'tl_rms_m95_db',tl_rms,'rho_shape',rho);
end

function tab = local_convergence_table(runs,fp,parameter_field)
n = numel(runs);
stage = strings(n,1); parameter = zeros(n,1); nodes = zeros(n,1);
boundary = zeros(n,1); linear = zeros(n,1); l2 = zeros(n,1);
phase = zeros(n,1); tl = zeros(n,1); rho = zeros(n,1);
successive_l2 = nan(n,1);
for ii = 1:n
    r = runs{ii}; stage(ii)=string(parameter_field); parameter(ii)=r.(parameter_field);
    nodes(ii)=r.node_count; boundary(ii)=r.offgrid_boundary_residual;
    linear(ii)=r.linear_residual;l2(ii)=r.metrics.complex_l2_m99;
    phase(ii)=r.metrics.phase_rms_m99;tl(ii)=r.metrics.tl_rms_m95_db;
    rho(ii)=r.metrics.rho_shape;
    if ii>1
        successive_l2(ii)=local_metrics(r.receiver_field,runs{ii-1}.receiver_field,fp).complex_l2_m99;
    end
end
tab=table(stage,parameter,nodes,boundary,linear,l2,phase,tl,rho,successive_l2, ...
    'VariableNames',{'sweep','parameter','node_count','boundary_residual', ...
    'linear_residual','image_complex_l2_M99','image_phase_rms_M99', ...
    'image_tl_rms_M95_db','rho_shape','successive_complex_l2_M99'});
end

function tab = local_rough_convergence_table(runs,reference,fp,parameter_field)
n = numel(runs);
sweep = strings(n,1); parameter = zeros(n,1); nodes = zeros(n,1);
boundary = zeros(n,1); linear = zeros(n,1); to_reference = zeros(n,1);
successive_l2 = nan(n,1); successive_phase = nan(n,1); successive_tl = nan(n,1);
for ii = 1:n
    q = runs{ii};
    sweep(ii) = string(parameter_field);
    parameter(ii) = q.(parameter_field);
    nodes(ii) = q.node_count;
    boundary(ii) = q.offgrid_boundary_residual;
    linear(ii) = q.linear_residual;
    to_reference(ii) = local_metrics(q.receiver_field,reference,fp).complex_l2_m99;
    if ii > 1
        pair = local_metrics(q.receiver_field,runs{ii-1}.receiver_field,fp);
        successive_l2(ii) = pair.complex_l2_m99;
        successive_phase(ii) = pair.phase_rms_m99;
        successive_tl(ii) = pair.tl_rms_m95_db;
    end
end
tab = table(sweep,parameter,nodes,boundary,linear,to_reference, ...
    successive_l2,successive_phase,successive_tl, ...
    'VariableNames',{'sweep','parameter','node_count','boundary_residual', ...
    'linear_residual','complex_l2_to_reference_M99', ...
    'successive_complex_l2_M99','successive_phase_rms_M99', ...
    'successive_tl_rms_M95_db'});
end

function local_write_r2_report(path,v,mat_file,csv_file)
fid = fopen(path,'w','n','UTF-8');
if fid < 0, error('Cannot write %s.',path); end
cl = onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'# PE--Bellhop--Helmholtz BIE R2 weak-rough convergence\n\n');
fprintf(fid,'日期：2026-09-14  \n状态：**%s**\n\n', ...
    ternary(v.passed,'R2 PASS','R2 FAIL / R3 LOCKED'));
fprintf(fid,'## Frozen benchmark\n\n');
fprintf(fid,'- `eta_bench(x)=A sin(Kx) chi(x)`, `A=%.9g m`, `K=%.9g rad/m`.\n', ...
    v.surface_meta.amplitude_m,v.surface_meta.wavenumber_radpm);
fprintf(fid,'- `chi=1` for `|x|<=%.9g m`; fixed physical support ends at `|x|=%.9g m`.\n', ...
    v.surface_meta.taper_inner_m,v.surface_meta.support_m);
fprintf(fid,'- This physical taper is fixed in every run. The BIE window is a separate numerical window outside that support.\n');
fprintf(fid,'- source: `%s`\n\n',v.source_fingerprint);
fprintf(fid,'## Numerical uncertainty\n\n');
fprintf(fid,'- `U_BIE` complex L2 M99: `%.9g`\n',v.U_BIE_complex_l2_M99);
fprintf(fid,'- phase RMS M99 uncertainty: `%.9g rad`\n',v.U_BIE_phase_rms_M99_rad);
fprintf(fid,'- TL RMS M95 uncertainty: `%.9g dB`\n\n',v.U_BIE_tl_rms_M95_db);
fprintf(fid,'| sweep | parameter | nodes | boundary residual | linear residual | L2 to reference | successive L2 | successive phase | successive TL dB |\n');
fprintf(fid,'|---|---:|---:|---:|---:|---:|---:|---:|---:|\n');
t = [v.spatial_convergence;v.quadrature_convergence;v.window_convergence];
for ii = 1:height(t)
    fprintf(fid,'| %s | %.9g | %d | %.3g | %.3g | %.3g | %.3g | %.3g | %.3g |\n', ...
        t.sweep(ii),t.parameter(ii),t.node_count(ii),t.boundary_residual(ii), ...
        t.linear_residual(ii),t.complex_l2_to_reference_M99(ii), ...
        t.successive_complex_l2_M99(ii),t.successive_phase_rms_M99(ii), ...
        t.successive_tl_rms_M95_db(ii));
end
fprintf(fid,'\n## Gates\n\n');
names = fieldnames(v.checks);
for ii = 1:numel(names)
    fprintf(fid,'- `%s`: %s\n',names{ii}, ...
        ternary(v.checks.(names{ii}),'PASS','FAIL'));
end
fprintf(fid,'\nArtifacts: `%s`, `%s`.\n\n',mat_file,csv_file);
fprintf(fid,'R3 remains locked unless every R2 gate passes.\n');
end

function local_write_r4_report(path,v,mat_file,csv_file)
fid=fopen(path,'w','n','UTF-8');if fid<0,error('Cannot write %s.',path);end
cl=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'# PE--Bellhop--Helmholtz BIE %s\n\n',v.stage);
fprintf(fid,'日期：2026-09-15  \n状态：**%s**\n\n',v.conclusion);
fprintf(fid,'- Surface: `A=%.9g m`, `K=%.9g rad/m`, same fixed C2 taper.\n', ...
    v.surface_meta.amplitude_m,v.surface_meta.wavenumber_radpm);
fprintf(fid,'- U_BIE: `%.9g`; required 5U separation: `%.9g`; observed: `%.9g`.\n\n', ...
    v.U_BIE,v.five_U_threshold,v.separation);
if ~isempty(v.boundary_verification)
    fprintf(fid,['- Strong-case boundary cross-check: residual `%.9g`, ', ...
        'main/high-order receiver-field L2 `%.9g`, acceptance limit `%.9g`.\n', ...
        '- This is an explicit strong-height numerical limit; the R1 analytic ', ...
        'boundary gate remains unchanged at `1e-8`.\n\n'], ...
        v.boundary_verification.offgrid_boundary_residual, ...
        v.boundary_verification_field_l2,v.boundary_acceptance_limit);
end
fprintf(fid,'| pair | E_G | phase RMS | TL RMS dB | rho | phi0 | E_aligned |\n');
fprintf(fid,'|---|---:|---:|---:|---:|---:|---:|\n');
for name={'PE_BIE','BH_BIE','PE_BH','BH_native_BIE_native'}
    m=v.metrics.(name{1});fprintf(fid,'| %s | %.9g | %.9g | %.9g | %.9g | %.9g | %.9g |\n', ...
        name{1},m.E_G,m.phase_rms_rad,m.tl_rms_db,m.rho_shape,m.phi0_rad,m.E_aligned);
end
fprintf(fid,'\n## Refinement stability\n\n');
fprintf(fid,'| kind | parameter | E_PE,BIE | E_BH,BIE | separation | winner |\n');
fprintf(fid,'|---|---:|---:|---:|---:|---|\n');
for ii=1:height(v.stability)
    q=v.stability(ii,:);fprintf(fid,'| %s | %.9g | %.9g | %.9g | %.9g | %s |\n', ...
        q.kind,q.parameter,q.E_PE,q.E_BH,q.separation,q.winner);
end
fprintf(fid,'\n## Gates\n\n');n=fieldnames(v.checks);
for ii=1:numel(n),fprintf(fid,'- `%s`: %s\n',n{ii},ternary(v.checks.(n{ii}),'PASS','FAIL'));end
fprintf(fid,'\nArtifacts: `%s`, `%s`.\n',mat_file,csv_file);
end

function local_write_r3_report(path,v,mat_file)
fid = fopen(path,'w','n','UTF-8');
if fid < 0, error('Cannot write %s.',path); end
cl = onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'# PE--Bellhop--Helmholtz BIE R3 weak three-way closure\n\n');
fprintf(fid,'日期：2026-09-14  \n状态：**%s**\n\n', ...
    ternary(v.passed,'R3 PASS','R3 FAIL / STOP'));
fprintf(fid,'- Surface: `A=%.9g m`, `K=%.9g rad/m`, fixed C2 physical taper `%g -> %g m`.\n', ...
    v.surface_meta.amplitude_m,v.surface_meta.wavenumber_radpm, ...
    v.surface_meta.taper_inner_m,v.surface_meta.support_m);
fprintf(fid,'- Convention: `%s`.\n',v.comparison_convention);
fprintf(fid,'- Thresholds (existing weak-limit floor + U_BIE): E `%.9g`, phase `%.9g rad`, TL `%.9g dB`.\n\n', ...
    v.thresholds.E,v.thresholds.phase_rad,v.thresholds.tl_db);
fprintf(fid,'| pair | E_G | phase RMS rad | TL RMS dB | rho_raw | rho_shape | phi0 rad | E_aligned |\n');
fprintf(fid,'|---|---:|---:|---:|---:|---:|---:|---:|\n');
names = {'PE_BIE','BH_BIE','PE_BH','BH_native_BIE_native', ...
    'PE_BIE_native_diagnostic'};
for ii = 1:numel(names)
    m = v.metrics.(names{ii});
    fprintf(fid,'| %s | %.9g | %.9g | %.9g | %.9g | %.9g | %.9g | %.9g |\n', ...
        names{ii},m.E_G,m.phase_rms_rad,m.tl_rms_db,m.rho_raw,m.rho_shape, ...
        m.phi0_rad,m.E_aligned);
end
fprintf(fid,'\n## Bellhop geometry\n\n');
fprintf(fid,'- wall residual max: `%.9g m`; phase-jump error max: `%.9g rad`.\n', ...
    v.geometry.wall_residual_max_m,v.geometry.phase_jump_error_max_rad);
fprintf(fid,'- p/q rotation errors: `%.9g / %.9g`; minimum post-wall dr: `%.9g m`.\n\n', ...
    v.geometry.p_rotation_error_max,v.geometry.q_rotation_error_max, ...
    v.geometry.min_post_dr_m);
fprintf(fid,'## Gates\n\n');
checks = fieldnames(v.checks);
for ii = 1:numel(checks)
    fprintf(fid,'- `%s`: %s\n',checks{ii}, ...
        ternary(v.checks.(checks{ii}),'PASS','FAIL'));
end
fprintf(fid,'\nArtifact: `%s`.\n',mat_file);
end

function local_write_report(path,v,mat_file,csv_file)
fid=fopen(path,'w','n','UTF-8');if fid<0,error('Cannot write %s.',path);end
cl=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'# PE--Bellhop--Helmholtz BIE R1 flat validation\n\n');
fprintf(fid,'日期：2026-09-14  \n状态：**%s**\n\n',ternary(v.passed,'R1 PASS','R1 FAIL / ROUGH STAGES LOCKED'));
fprintf(fid,'- source: `%s`\n',v.source_fingerprint);
fprintf(fid,'- formulation: `%s`\n',v.main.formulation);
fprintf(fid,'- convention/normal/jump: `%s`; `%s`; `%s`\n\n', ...
    v.main.time_convention,v.main.normal_orientation,v.main.jump_relation);
fprintf(fid,'- M95/M99 radii: `%.9g / %.9g m`\n\n',v.footprint.radius95_m,v.footprint.radius99_m);
fprintf(fid,'| sweep | parameter | nodes | boundary residual | linear residual | image L2 M99 | phase RMS M99 | TL RMS M95 (dB) | rho | successive L2 |\n');
fprintf(fid,'|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n');
t=[v.spatial_convergence;v.window_convergence];
for ii=1:height(t)
    fprintf(fid,'| %s | %.9g | %d | %.3g | %.3g | %.3g | %.3g | %.3g | %.9g | %.3g |\n', ...
        t.sweep(ii),t.parameter(ii),t.node_count(ii),t.boundary_residual(ii), ...
        t.linear_residual(ii),t.image_complex_l2_M99(ii),t.image_phase_rms_M99(ii), ...
        t.image_tl_rms_M95_db(ii),t.rho_shape(ii),t.successive_complex_l2_M99(ii));
end
fprintf(fid,'\n## Hard gates\n\n');
names=fieldnames(v.checks);
for ii=1:numel(names)
    fprintf(fid,'- `%s`: %s\n',names{ii},ternary(v.checks.(names{ii}),'PASS','FAIL'));
end
fprintf(fid,'\nArtifacts: `%s`, `%s`.\n\n',mat_file,csv_file);
fprintf(fid,'Rough stages remain locked unless every R1 hard gate passes.\n');
end

function y=ternary(c,a,b)
if c,y=a;else,y=b;end
end
