function validation = validate_pe_bellhop_incident_field_vertical(overrides)
%VALIDATE_PE_BELLHOP_INCIDENT_FIELD_VERTICAL
% Compare the pre-reflection transverse complex field at r=100 m.
%
% This is a validation-only 1-transverse-dimensional PE bridge versus
% Bellhop 2-D coherent free-field comparison.  No PM profile, wall,
% Reflect2D call, PE production march, or receiver reflected branch is used.

if nargin < 1 || isempty(overrides), overrides = struct(); end
if ~isstruct(overrides) || ~isscalar(overrides)
    error('overrides must be a scalar struct.');
end

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
run(fullfile(root,'scripts','bootstrap_project.m'));
cfg = local_config(root,overrides);
local_validate_config(cfg);
if ~exist(cfg.output_dir,'dir'), mkdir(cfg.output_dir); end

bh_norm = local_load_normalization(cfg);
pe = local_run_pe(cfg);
as = local_exact_1d(cfg,pe.x_m);
footprint = local_footprint(pe.x_m,as.field,cfg);
source_pattern = local_pattern(cfg);

pe_as = local_pair_metrics(pe.surface_incident_field,as.field,footprint,pe.axis_index);

bh5001 = local_run_bellhop(cfg,cfg.bellhop_beam_counts(1),pe.x_m,'f4000_b5001');
bh10001 = local_run_bellhop(cfg,cfg.bellhop_beam_counts(2),pe.x_m,'f4000_b10001');

bh5001_converted = local_convert(bh5001.pressure_raw,bh_norm,cfg.bellhop_source_geometry);
bh10001_converted = local_convert(bh10001.pressure_raw,bh_norm,cfg.bellhop_source_geometry);
bh5001.pressure_converted = bh5001_converted;
bh10001.pressure_converted = bh10001_converted;
pe.receiver_coordinate_error_m = max(bh10001.receiver_coordinate_error_m);

bh_as_5001 = local_pair_metrics(bh5001_converted,as.field,footprint,pe.axis_index);
bh_as_10001 = local_pair_metrics(bh10001_converted,as.field,footprint,pe.axis_index);
beam = local_pair_metrics(bh5001_converted,bh10001_converted,footprint,pe.axis_index);
pe_bh = local_pair_metrics(pe.surface_incident_field,bh10001_converted,footprint,pe.axis_index);
symmetry = struct('pe',local_symmetry_metric(pe.x_m,pe.surface_incident_field,pe.axis_index), ...
    'as',local_symmetry_metric(as.x_m,as.field,pe.axis_index), ...
    'bellhop',local_symmetry_metric(pe.x_m,bh10001_converted,pe.axis_index));

sampling = struct('executed',false,'metrics',local_empty_metrics(), ...
    'shared_count',0,'max_coordinate_error_m',NaN,'case',struct());
if cfg.run_half_dx_sampling
    x_half = (-0.5*cfg.xw_m) + (0:(2*cfg.nx)-1)*(cfg.xw_m/(2*cfg.nx));
    bh_half = local_run_bellhop(cfg,cfg.bellhop_beam_counts(2),x_half,'f4000_b10001_halfdx');
    bh_half.pressure_converted = local_convert(bh_half.pressure_raw,bh_norm,cfg.bellhop_source_geometry);
    [shared_ix,half_ix,coord_error] = local_shared_grid(pe.x_m,x_half,cfg.receiver_tolerance_m);
    shared_foot = local_footprint_subset(footprint,shared_ix);
    half_shared = bh_half.pressure_converted(half_ix);
    main_shared = bh10001_converted(shared_ix);
    sampling.executed = true;
    sampling.metrics = local_pair_metrics(main_shared,half_shared,shared_foot,find(shared_ix==pe.axis_index,1));
    sampling.shared_count = numel(shared_ix);
    sampling.max_coordinate_error_m = max(coord_error);
    sampling.case = bh_half;
end

finite_count = local_nonfinite_count(pe.surface_incident_field,as.field,bh5001_converted, ...
    bh10001_converted);
checks = local_checks(pe,pe_as,beam,pe_bh,sampling,finite_count,cfg);

validation = struct('schema_version','1.0.0', ...
    'stage','incident_field_flat_freefield_4khz', ...
    'config',cfg,'bellhop_normalization',bh_norm,'source_pattern',source_pattern, ...
    'pe',pe,'as',as, ...
    'footprint',footprint,'bellhop_5001',bh5001,'bellhop_10001',bh10001, ...
    'metrics',struct('pe_as',pe_as,'bellhop_as_5001',bh_as_5001, ...
        'bellhop_as_10001',bh_as_10001,'beam_convergence',beam,'pe_bellhop',pe_bh), ...
    'receiver_sampling',sampling,'symmetry',symmetry,'nonfinite_count',finite_count, ...
    'checks',checks,'passed',all(checks.passed));
validation.files = local_outputs(validation);
if cfg.fail_on_check && ~validation.passed
    error('Incident-field comparison failed; see %s.',validation.files.report);
end
end

function cfg = local_config(root,overrides)
exe = getenv('BELLHOP_EXE');
if isempty(exe)
    candidate = 'E:\MISC\BELLHOP\AcousticsToolbox_2020\windows-bin-20201102\bellhop.exe';
    if exist(candidate,'file') == 2, exe = candidate; end
end
cfg = struct( ...
    'output_dir',fullfile(root,'results','validation','pe_bellhop_incident_field'), ...
    'report_path',fullfile(root,'reports','pe_bellhop_incident_field_comparison_report.md'), ...
    'bellhop_exe',exe,'bellhop_toolbox_version','OALIB AcousticsToolbox 2020_11_4 (2020-11-02 binary)', ...
    'normalization_audit_file',fullfile(root,'results','validation','pe_bellhop_freefield', ...
        'bellhop_normalization','bellhop_freefield_normalization.mat'), ...
    'frequency_hz',4000,'c0_mps',1500,'z_tx_m',100,'sigma_src_m',0.3, ...
    'xw_m',192.1875,'nx',984,'step_m',0.05,'incident_range_m',100, ...
    'guard_range_m',99,'bellhop_beam_counts',[5001 10001], ...
    'bellhop_angle_limits_deg',[-30 30],'bellhop_step_m',0.05, ...
    'bellhop_source_geometry','X','bellhop_run_type','C', ...
    'bellhop_domain_half_depth_m',1000,'source_pattern_samples',2401, ...
    'source_pattern_clip_db',-120,'receiver_tolerance_m',1e-6, ...
    'run_half_dx_sampling',true,'phase_floor_db',-40, ...
    'pe_as_complex_limit',1e-10,'pe_as_l2_limit',1e-10, ...
    'pe_outer5_energy_limit',1e-5,'beam_l2_limit',2e-3, ...
    'beam_phase_rms_limit',5e-3,'beam_tl_rms_limit_db',0.02, ...
    'pe_bh_l2_m99_limit',0.02,'pe_bh_phase_rms_m99_limit',0.02, ...
    'pe_bh_phase_p95_m95_limit',0.05,'pe_bh_tl_p95_m95_limit_db',0.10, ...
    'receiver_sampling_l2_limit',5e-3,'fail_on_check',true);
names = fieldnames(overrides);
for ii = 1:numel(names)
    if ~isfield(cfg,names{ii}), error('Unknown override: %s.',names{ii}); end
    cfg.(names{ii}) = overrides.(names{ii});
end
if isfield(overrides,'output_dir') && ~isfield(overrides,'report_path')
    cfg.report_path = fullfile(cfg.output_dir,'pe_bellhop_incident_field_comparison_report.md');
end
end

function local_validate_config(cfg)
if isempty(cfg.bellhop_exe) || exist(cfg.bellhop_exe,'file') ~= 2
    error('Bellhop 2020 executable is missing; set BELLHOP_EXE or cfg.bellhop_exe.');
end
if exist(cfg.normalization_audit_file,'file') ~= 2
    error('Saved Bellhop normalization audit is missing: %s',cfg.normalization_audit_file);
end
if cfg.frequency_hz <= 0 || cfg.c0_mps <= 0 || cfg.sigma_src_m <= 0 || ...
        cfg.xw_m <= 0 || cfg.nx < 4 || mod(cfg.nx,2) ~= 0 || cfg.step_m <= 0
    error('Invalid incident-field configuration.');
end
if cfg.incident_range_m <= 0 || cfg.guard_range_m <= 0 || ...
        cfg.guard_range_m >= cfg.incident_range_m
    error('guard_range_m must be positive and smaller than incident_range_m.');
end
if numel(cfg.bellhop_beam_counts) ~= 2 || any(cfg.bellhop_beam_counts < 2)
    error('bellhop_beam_counts must contain [low high].');
end
source_geometry=upper(char(cfg.bellhop_source_geometry));
if ~isscalar(source_geometry) || ~ismember(source_geometry,['R','X'])
    error('bellhop_source_geometry must be R (point) or X (line).');
end
if ~strcmpi(char(cfg.bellhop_run_type),'C')
    error('This validation requires coherent Bellhop run type C.');
end
end

function norm = local_load_normalization(cfg)
loaded = load(cfg.normalization_audit_file,'validation');
if ~isfield(loaded,'validation') || ~isfield(loaded.validation,'passed') || ...
        ~loaded.validation.passed
    error('Saved Bellhop normalization audit is missing or failed.');
end
norm = loaded.validation;
required = {'selected_spatial_sign','bellhop_source_constant'};
for ii = 1:numel(required)
    if ~isfield(norm,required{ii}), error('Normalization audit lacks %s.',required{ii}); end
end
if ~isscalar(norm.selected_spatial_sign) || ~ismember(norm.selected_spatial_sign,[-1 1])
    error('Normalization selected_spatial_sign must be +1 or -1.');
end
end

function pe = local_run_pe(cfg)
pe_cfg = struct('frequency_hz',cfg.frequency_hz,'c0_mps',cfg.c0_mps, ...
    'xw_m',cfg.xw_m,'nx',cfg.nx,'z_tx_m',cfg.z_tx_m,'z_rx_m',0, ...
    'sigma_src_m',cfg.sigma_src_m,'surface_elevation_x_m',zeros(1,cfg.nx), ...
    'surface_reflect_coeff',-1,'step_m',cfg.step_m,'x_rx_m',0);
out = run_pe_1d_surface_reflection_validation(pe_cfg);
[~,axis_index] = min(abs(out.x_m));
pe = struct('config',pe_cfg,'x_m',out.x_m(:).','axis_index',axis_index, ...
    'surface_incident_field',out.surface_incident_field(:).', ...
    'edge_energy_fraction',local_outer5_energy(out.surface_incident_field));
end

function as = local_exact_1d(cfg,x)
kx = (2*pi/cfg.xw_m)*[0:(cfg.nx/2-1),-cfg.nx/2:-1];
k0 = 2*pi*cfg.frequency_hz/cfg.c0_mps;
kz = sqrt(complex(k0^2-kx.^2,0));
source = exp(-0.5*(x/cfg.sigma_src_m).^2);
field = ifft(fft(source).*exp(1i*cfg.incident_range_m*(kz-k0)));
as = struct('x_m',x(:).','field',field(:).','k0_rad_per_m',k0, ...
    'operator','independent one-step exact 1-D angular spectrum');
end

function run = local_run_bellhop(cfg,beams,receiver_depths,tag)
pat = local_pattern(cfg);
case_root = fullfile(cfg.output_dir,'bellhop',sprintf('%s_%s',tag,lower(cfg.bellhop_source_geometry)));
c = struct('bellhop_exe',cfg.bellhop_exe,'case_root',case_root, ...
    'run_type',cfg.bellhop_run_type,'source_geometry',cfg.bellhop_source_geometry, ...
    'frequency_hz',cfg.frequency_hz,'c0_mps',cfg.c0_mps, ...
    'source_depth_m',0,'receiver_depths_m',receiver_depths(:).', ...
    'receiver_ranges_m',[cfg.guard_range_m cfg.incident_range_m], ...
    'beam_count',beams,'angle_limits_deg',cfg.bellhop_angle_limits_deg, ...
    'step_m',cfg.bellhop_step_m,'domain_half_depth_m',cfg.bellhop_domain_half_depth_m, ...
    'source_pattern_angles_deg',pat.angles_deg,'source_pattern_level_db',pat.level_db);
raw = run_bellhop_unfolded_gaussian_vertical(c);
data = raw.data;
ranges = double(data.receiver_range_m(:));
range_ix = find(abs(ranges-cfg.incident_range_m)<=cfg.receiver_tolerance_m);
if numel(range_ix) ~= 1
    error('Expected one explicit %.12g m SHD range in %s.',cfg.incident_range_m,raw.files.data);
end
guard_ix = find(abs(ranges-cfg.guard_range_m)<=cfg.receiver_tolerance_m);
if numel(guard_ix) ~= 1
    error('Expected one explicit %.12g m guard range in %s.',cfg.guard_range_m,raw.files.data);
end
depths = double(data.receiver_depth_m(:)); req = double(receiver_depths(:));
depth_ix = zeros(size(req)); coordinate_error = zeros(size(req));
for ii = 1:numel(req)
    [coordinate_error(ii),depth_ix(ii)] = min(abs(depths-req(ii)));
end
if any(coordinate_error>cfg.receiver_tolerance_m)
    error('Bellhop receiver-depth coordinates do not match requested grid.');
end
run = raw;
run.requested_receiver_depth_m = req(:).';
run.receiver_coordinate_error_m = coordinate_error(:).';
run.range_index = range_ix;
run.guard_range_index = guard_ix;
run.pressure_raw = zeros(1,numel(req));
run.guard_pressure_raw = zeros(1,numel(req));
for ii = 1:numel(req)
    run.pressure_raw(ii) = select_bellhop_shd_pressure_at_range_vertical( ...
        data,cfg.incident_range_m,req(ii),cfg.receiver_tolerance_m);
    run.guard_pressure_raw(ii) = select_bellhop_shd_pressure_at_range_vertical( ...
        data,cfg.guard_range_m,req(ii),cfg.receiver_tolerance_m);
end
run.actual_receiver_depth_m = depths(depth_ix).';
run.actual_range_m = ranges(range_ix);
run.actual_guard_range_m = ranges(guard_ix);
run.beam_count = beams;
run.case_tag = tag;
end

function pat = local_pattern(cfg)
angles = linspace(cfg.bellhop_angle_limits_deg(1),cfg.bellhop_angle_limits_deg(2), ...
    cfg.source_pattern_samples).';
k = 2*pi*cfg.frequency_hz/cfg.c0_mps; theta = angles*pi/180;
d = cos(theta).*exp(-0.5*(k*cfg.sigma_src_m*sin(theta)).^2);
d = abs(d)/max(abs(d));
pat = struct('angles_deg',angles,'level_db',20*log10(max(d,10^(cfg.source_pattern_clip_db/20))), ...
    'formula','cos(theta)*exp(-0.5*(k*sigma*sin(theta))^2)', ...
    'fingerprint',sprintf('%s|cos(theta)*exp(-0.5*(k*sigma*sin(theta))^2)|N=%d|angle=[%.12g,%.12g]deg|f=%.12gHz|sigma=%.12gm', ...
        upper(char(cfg.bellhop_source_geometry)),cfg.source_pattern_samples,cfg.bellhop_angle_limits_deg(1), ...
        cfg.bellhop_angle_limits_deg(2),cfg.frequency_hz,cfg.sigma_src_m));
end

function p = local_convert(raw,norm,source_geometry)
if norm.selected_spatial_sign < 0
    p = conj(raw);
else
    p = raw;
end
% The saved source constant belongs to the historical point-source audit.
% It is only applicable to R. The active X comparison is axis-normalized,
% so no empirical or mismatched global amplitude factor is introduced.
if upper(char(source_geometry))=='R'
    if norm.selected_spatial_sign < 0
        p=p./conj(norm.bellhop_source_constant);
    else
        p=p./norm.bellhop_source_constant;
    end
end
end

function fp = local_footprint(x,as,cfg)
weights = abs(as(:)).^2;
radii = sort(unique(abs(x(:))));
fp = struct('x_m',x(:).','energy_weights',weights(:).', ...
    'axis_index',find(abs(x)==min(abs(x)),1),'dx_m',x(2)-x(1));
for pp = [0.95 0.99]
    found = false;
    for ii = 1:numel(radii)
        mask = abs(x(:)) <= radii(ii)+10*eps(max(1,radii(ii)));
        fraction = sum(weights(mask))/max(sum(weights),realmin);
        if fraction >= pp
            name = sprintf('m%d',round(100*pp));
            fp.(name) = struct('radius_m',radii(ii),'mask',mask(:).', ...
                'energy_fraction',fraction,'count',nnz(mask));
            found = true; break;
        end
    end
    if ~found, error('Unable to construct %.0f%% energy footprint.',100*pp); end
end
fp.phase_floor_db = cfg.phase_floor_db;
end

function sub = local_footprint_subset(fp,indices)
indices = indices(:).';
sub = fp;
sub.x_m = fp.x_m(indices);
sub.energy_weights = fp.energy_weights(indices);
sub.axis_index = find(indices==fp.axis_index,1);
sub.m95 = struct('mask',fp.m95.mask(indices),'count',nnz(fp.m95.mask(indices)), ...
    'energy_fraction',sum(fp.energy_weights(indices).*fp.m95.mask(indices))/max(sum(fp.energy_weights(indices)),realmin));
sub.m99 = struct('mask',fp.m99.mask(indices),'count',nnz(fp.m99.mask(indices)), ...
    'energy_fraction',sum(fp.energy_weights(indices).*fp.m99.mask(indices))/max(sum(fp.energy_weights(indices)),realmin));
end

function [main_ix,half_ix,err] = local_shared_grid(x,half,tol)
main_ix = 1:numel(x); half_ix = zeros(size(main_ix)); err = zeros(size(main_ix));
for ii = 1:numel(main_ix)
    [err(ii),half_ix(ii)] = min(abs(half-x(main_ix(ii))));
end
if any(err>tol), error('Half-dx receiver grid does not contain the PE grid.'); end
end

function m = local_pair_metrics(reference,candidate,fp,axis_index)
reference = reference(:).'; candidate = candidate(:).';
if numel(reference) ~= numel(candidate) || numel(reference) ~= numel(fp.x_m)
    error('Pair metric arrays and footprint must have identical lengths.');
end
ref_n = reference/reference(axis_index); cand_n = candidate/candidate(axis_index);
phase = angle(ref_n.*conj(cand_n));
tl = 20*log10(max(abs(ref_n),realmin)./max(abs(cand_n),realmin));
complex_error = abs(ref_n-cand_n);
mask95 = fp.m95.mask; mask99 = fp.m99.mask;
valid95 = mask95 & abs(ref_n)>=10^(fp.phase_floor_db/20) & abs(cand_n)>=10^(fp.phase_floor_db/20);
w = fp.energy_weights(:).';
m = struct('max_complex_full',max(complex_error), ...
    'l2_full',sqrt(sum(complex_error.^2)/max(sum(abs(ref_n).^2),realmin)), ...
    'l2_m99',sqrt(sum(complex_error(mask99).^2)/max(sum(abs(ref_n(mask99)).^2),realmin)), ...
    'phase_rms_m99',local_weighted_rms(phase(mask99),w(mask99)), ...
    'tl_rms_m99',local_weighted_rms(tl(mask99),w(mask99)), ...
    'phase_rms_m95',sqrt(mean(phase(valid95).^2)), ...
    'tl_rms_m95',sqrt(mean(tl(valid95).^2)), ...
    'phase_p95_m95',local_percentile(abs(phase(valid95)),0.95), ...
    'tl_p95_m95',local_percentile(abs(tl(valid95)),0.95), ...
    'alpha_ls',sum(w(mask99).*conj(cand_n(mask99)).*ref_n(mask99))/ ...
        max(sum(w(mask99).*abs(cand_n(mask99)).^2),realmin), ...
    'valid_count_m95',nnz(valid95),'phase_floor_count_m95',nnz(mask95 & ~valid95), ...
    'reference_axis',reference(axis_index),'candidate_axis',candidate(axis_index), ...
    'reference_normalized',ref_n,'candidate_normalized',cand_n, ...
    'phase_difference',phase,'tl_difference_db',tl,'complex_difference',complex_error);
m.aligned_l2_m99 = sqrt(sum(w(mask99).*abs(ref_n(mask99)-m.alpha_ls*cand_n(mask99)).^2)/ ...
    max(sum(w(mask99).*abs(ref_n(mask99)).^2),realmin));
end

function value = local_weighted_rms(values,weights)
weights = weights(:); values = values(:);
value = sqrt(sum(weights.*values.^2)/max(sum(weights),realmin));
end

function value = local_percentile(values,p)
values = sort(values(isfinite(values)));
if isempty(values), value = Inf; return; end
ix = max(1,min(numel(values),ceil(p*numel(values)))); value = values(ix);
end

function value = local_outer5_energy(field)
n = max(1,round(0.05*numel(field))); edge = [1:n (numel(field)-n+1):numel(field)];
value = sum(abs(field(edge)).^2)/max(sum(abs(field(:)).^2),realmin);
end

function s = local_symmetry_metric(x,field,axis_index)
% Diagnostic only: pair each coordinate with its closest negative partner.
x = x(:); field = field(:)/field(axis_index);
pair_error = [];
for ii = 1:numel(x)
    [distance,jj] = min(abs(x+x(ii)));
    if distance <= 1e-9*max(1,max(abs(x)))
        pair_error(end+1,1) = abs(field(ii)-field(jj)); %#ok<AGROW>
    end
end
if isempty(pair_error)
    s = struct('paired_count',0,'max_complex_error',NaN,'rms_complex_error',NaN);
else
    s = struct('paired_count',numel(pair_error),'max_complex_error',max(pair_error), ...
        'rms_complex_error',sqrt(mean(pair_error.^2)));
end
end

function n = local_nonfinite_count(varargin)
n = 0;
for ii = 1:nargin
    n = n + nnz(~isfinite(varargin{ii}));
end
end

function checks = local_checks(pe,pe_as,beam,pe_bh,sampling,finite_count,cfg)
names = ["pe_as_max_complex_full";"pe_as_l2_m99";"pe_outer5_energy"; ...
    "receiver_coordinate_error";"bellhop_beam_l2_m99";"bellhop_beam_phase_rms_m95"; ...
    "bellhop_beam_tl_rms_m95";"pe_bellhop_l2_m99";"pe_bellhop_phase_rms_m99"; ...
    "pe_bellhop_phase_p95_m95";"pe_bellhop_tl_p95_m95"; ...
    "receiver_sampling_l2_m99";"nonfinite_count"];
values = [pe_as.max_complex_full;pe_as.l2_m99;pe.edge_energy_fraction; ...
    max([pe.receiver_coordinate_error_m]);beam.l2_m99;beam.phase_rms_m95; ...
    beam.tl_rms_m95;pe_bh.l2_m99;pe_bh.phase_rms_m99;pe_bh.phase_p95_m95; ...
    pe_bh.tl_p95_m95;local_sampling_value(sampling);finite_count];
limits = [cfg.pe_as_complex_limit;cfg.pe_as_l2_limit;cfg.pe_outer5_energy_limit; ...
    cfg.receiver_tolerance_m;cfg.beam_l2_limit;cfg.beam_phase_rms_limit; ...
    cfg.beam_tl_rms_limit_db;cfg.pe_bh_l2_m99_limit;cfg.pe_bh_phase_rms_m99_limit; ...
    cfg.pe_bh_phase_p95_m95_limit;cfg.pe_bh_tl_p95_m95_limit_db; ...
    cfg.receiver_sampling_l2_limit;0];
passed = values <= limits;
if ~sampling.executed
    passed(12) = true; values(12) = NaN; limits(12) = NaN;
end
passed(13) = finite_count == 0;
checks = table(names,values,limits,passed,'VariableNames',{'check_name','value','limit','passed'});
end

function value = local_sampling_value(sampling)
if sampling.executed, value = sampling.metrics.l2_m99; else, value = NaN; end
end

function files = local_outputs(v)
out = v.config.output_dir; if ~exist(out,'dir'), mkdir(out); end
x = v.pe.x_m(:);
pe = v.pe.surface_incident_field(:); as = v.as.field(:);
bh = v.bellhop_10001.pressure_converted(:); raw = v.bellhop_10001.pressure_raw(:);
fp = v.footprint;
sample_table = table(x,real(pe),imag(pe),real(as),imag(as),real(raw),imag(raw), ...
    real(bh),imag(bh),fp.m95.mask(:),fp.m99.mask(:), ...
    'VariableNames',{'x_m','pe_real','pe_imag','as_real','as_imag', ...
    'bellhop_raw_real','bellhop_raw_imag','bellhop_converted_real', ...
    'bellhop_converted_imag','mask_m95','mask_m99'});
writetable(sample_table,fullfile(out,'incident_field_samples.csv'));
metric_names = ["pe_as";"bellhop_as_5001";"bellhop_as_10001";"beam_convergence";"pe_bellhop"];
metric_values = [v.metrics.pe_as.l2_m99;v.metrics.bellhop_as_5001.l2_m99; ...
    v.metrics.bellhop_as_10001.l2_m99;v.metrics.beam_convergence.l2_m99; ...
    v.metrics.pe_bellhop.l2_m99];
metric_phase = [v.metrics.pe_as.phase_rms_m99;v.metrics.bellhop_as_5001.phase_rms_m99; ...
    v.metrics.bellhop_as_10001.phase_rms_m99;v.metrics.beam_convergence.phase_rms_m99; ...
    v.metrics.pe_bellhop.phase_rms_m99];
alpha_values = [v.metrics.pe_as.alpha_ls;v.metrics.bellhop_as_5001.alpha_ls; ...
    v.metrics.bellhop_as_10001.alpha_ls;v.metrics.beam_convergence.alpha_ls; ...
    v.metrics.pe_bellhop.alpha_ls];
summary_table = table(metric_names,metric_values,metric_phase,abs(alpha_values),angle(alpha_values), ...
    'VariableNames',{'comparison','l2_m99','phase_rms_m99','alpha_ls_abs','alpha_ls_phase_rad'});
writetable(summary_table,fullfile(out,'incident_field_summary.csv'));
writetable(v.checks,fullfile(out,'incident_field_checks.csv'));
mat_file = fullfile(out,'pe_bellhop_incident_field_audit.mat'); validation = v; %#ok<NASGU>
save(mat_file,'validation','-v7.3');
fig_file = fullfile(out,'pe_bellhop_incident_field.png'); local_plot(v,fig_file);
local_report(v);
files = struct('mat',mat_file,'samples',fullfile(out,'incident_field_samples.csv'), ...
    'summary',fullfile(out,'incident_field_summary.csv'), ...
    'checks',fullfile(out,'incident_field_checks.csv'), ...
    'figure',fig_file,'report',v.config.report_path);
end

function local_plot(v,file)
x = v.pe.x_m; m = v.metrics.pe_bellhop; p = v.pe.surface_incident_field/v.pe.surface_incident_field(v.pe.axis_index);
b = v.bellhop_10001.pressure_converted/v.bellhop_10001.pressure_converted(v.pe.axis_index);
a = v.as.field/v.as.field(v.pe.axis_index);
fig = figure('Visible','off','Color','w','Position',[100 100 1200 800]); cleanup = onCleanup(@()close(fig));
subplot(2,2,1); plot(x,20*log10(max(abs(p),realmin)),x,20*log10(max(abs(b),realmin)),'--',x,20*log10(max(abs(a),realmin)),':'); grid on; xlim([-50 50]); xlabel('transverse coordinate (m)'); ylabel('normalized amplitude (dB)'); legend('PE','Bellhop','AS','Location','southwest');
subplot(2,2,2); plot(x,m.phase_difference); grid on; xlim([-50 50]); xlabel('transverse coordinate (m)'); ylabel('PE--Bellhop phase (rad)');
subplot(2,2,3); plot(x,m.tl_difference_db); grid on; xlim([-50 50]); xlabel('transverse coordinate (m)'); ylabel('PE/Bellhop amplitude difference (dB)');
subplot(2,2,4); semilogy(x,max(m.complex_difference,realmin)); grid on; xlim([-50 50]); xlabel('transverse coordinate (m)'); ylabel('normalized complex difference');
exportgraphics(fig,file,'Resolution',180); clear cleanup
end

function local_report(v)
fid = fopen(v.config.report_path,'w','n','UTF-8'); if fid < 0, error('Cannot write %s.',v.config.report_path); end
cleanup = onCleanup(@()fclose(fid)); c = v.config; fp = v.footprint;
fprintf(fid,'# PE--Bellhop 反射前入射复声场比较\n\n');
fprintf(fid,'状态：**%s**\n\n',ternary(v.passed,'PASS','FAIL'));
fprintf(fid,'本报告只比较均匀介质中 `r=%.6g m` 的反射前 incident plane；不启用 PM、海面反射、internal wall、`Reflect2D` 或 PE 生产二维横向场。\n\n',c.incident_range_m);
fprintf(fid,'## 配置\n\n');
fprintf(fid,'- f/c：%.0f Hz / %.0f m/s；Gaussian sigma：%.6g m；PE window/grid/step：%.9g m / %d / %.6g m。\n',c.frequency_hz,c.c0_mps,c.sigma_src_m,c.xw_m,c.nx,c.step_m);
fprintf(fid,'- Bellhop：`%s`；executable：`%s`；run type：`%s`；source：`%s`；beams：%s；step：%.6g m；angle fan：[%g,%g] deg；SHD ranges：[%g,%g] m。\n',c.bellhop_toolbox_version,c.bellhop_exe,c.bellhop_run_type,c.bellhop_source_geometry,sprintf('%d ',c.bellhop_beam_counts),c.bellhop_step_m,c.bellhop_angle_limits_deg,c.guard_range_m,c.incident_range_m);
fprintf(fid,'- source-pattern fingerprint：`%s`。\n',v.source_pattern.fingerprint);
fprintf(fid,'- footprint：M95 radius %.9g m (%d samples, energy %.9g); M99 radius %.9g m (%d samples, energy %.9g)。\n\n',fp.m95.radius_m,fp.m95.count,fp.m95.energy_fraction,fp.m99.radius_m,fp.m99.count,fp.m99.energy_fraction);
fprintf(fid,'## 复相位与归一化\n\n');
fprintf(fid,'PE uses the reduced-envelope one-step convention `exp(i*R*(kz-k0))` for the incident field. Bellhop raw pressure is converted with the saved normalization-audit spatial sign `%d` (the corresponding fixed conjugation is applied); no data-dependent conjugation, carrier re-addition, phase fitting, or empirical amplitude correction is used. Main comparisons are axis-normalized at x=0. For line-source X, the historical point-source global source constant is not applied.\n\n',v.bellhop_normalization.selected_spatial_sign);
fprintf(fid,'## 结果摘要\n\n');
fprintf(fid,'| comparison | L2 M99 | phase RMS M99 (rad) | phase P95 M95 (rad) | TL P95 M95 (dB) | LS scalar abs | LS scalar phase (rad) |\n|---|---:|---:|---:|---:|---:|---:|\n');
names = {'pe_as','bellhop_as_5001','bellhop_as_10001','beam_convergence','pe_bellhop'};
for ii=1:numel(names)
    m=v.metrics.(names{ii}); fprintf(fid,'| %s | %.8g | %.8g | %.8g | %.8g | %.8g | %.8g |\n',names{ii},m.l2_m99,m.phase_rms_m99,m.phase_p95_m95,m.tl_p95_m95,abs(m.alpha_ls),angle(m.alpha_ls));
end
fprintf(fid,'\nPE outer-5%% incident energy: `%.8g`; receiver coordinate error max: `%.8g m`; non-finite count: `%d`.\n\n',v.pe.edge_energy_fraction,max(v.bellhop_10001.receiver_coordinate_error_m),v.nonfinite_count);
fprintf(fid,'Symmetry diagnostic (not a hard gate): paired samples PE/Bellhop/AS = %d/%d/%d; maximum normalized complex pair error = %.8g / %.8g / %.8g.\n\n', ...
    v.symmetry.pe.paired_count,v.symmetry.bellhop.paired_count,v.symmetry.as.paired_count, ...
    v.symmetry.pe.max_complex_error,v.symmetry.bellhop.max_complex_error,v.symmetry.as.max_complex_error);
if v.receiver_sampling.executed
    fprintf(fid,'Half-dx receiver diagnostic: `%d` shared nodes; maximum coordinate error `%.8g m`; shared-node normalized L2 M99 `%.8g`; phase RMS M99 `%.8g rad`.\n\n', ...
        v.receiver_sampling.shared_count,v.receiver_sampling.max_coordinate_error_m, ...
        v.receiver_sampling.metrics.l2_m99,v.receiver_sampling.metrics.phase_rms_m99);
end
fprintf(fid,'## Hard gates\n\n| check | value | limit | pass |\n|---|---:|---:|:---:|\n');
for ii=1:height(v.checks), fprintf(fid,'| %s | %.8g | %.8g | %d |\n',v.checks.check_name(ii),v.checks.value(ii),v.checks.limit(ii),v.checks.passed(ii)); end
fprintf(fid,'\n');
if v.passed
    fprintf(fid,'4 kHz incident-plane gates pass. Conditional next step: repeat the same frozen convention and masks at 6/8 kHz.\n');
else
    fprintf(fid,'At least one predeclared gate failed. Do not extend to 6/8 kHz or rough-PM comparison until the failure is assigned to PE--AS, Bellhop convergence, receiver sampling, phase convention, or source/beam mapping.\n');
end
fprintf(fid,'\n![incident field comparison](../results/validation/pe_bellhop_incident_field/pe_bellhop_incident_field.png)\n');
end

function value = ternary(condition,yes_value,no_value)
if condition, value = yes_value; else, value = no_value; end
end

function m = local_empty_metrics()
m = struct('l2_m99',NaN,'phase_rms_m99',NaN,'tl_rms_m99',NaN, ...
    'phase_p95_m95',NaN,'tl_p95_m95',NaN,'alpha_ls',complex(NaN), ...
    'aligned_l2_m99',NaN);
end
