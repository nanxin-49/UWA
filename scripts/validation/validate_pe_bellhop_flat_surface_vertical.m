function validation = validate_pe_bellhop_flat_surface_vertical(overrides)
%VALIDATE_PE_BELLHOP_FLAT_SURFACE_VERTICAL Stage-1 PE/Bellhop cross-validation.
%   This deterministic case compares the existing 3-D vertical PE envelope
%   with Bellhop point-source arrivals for a uniform SSP and a flat,
%   pressure-release sea surface. Set BELLHOP_EXE to bellhop.exe when the
%   executable is not already on PATH or the MATLAB path.

if nargin < 1 || isempty(overrides)
    overrides = struct();
end

this_file = mfilename('fullpath');
project_root = fileparts(fileparts(fileparts(this_file)));
addpath(project_root);
addpath(fullfile(project_root, 'scripts'));

cfg = local_default_config();
cfg = local_apply_overrides(cfg, overrides);
local_validate_config(cfg);

output_dir = fullfile(project_root, 'results', 'validation', ...
    'pe_bellhop_flat_surface');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

fprintf('Stage-1 PE/Bellhop flat-surface cross-validation\n');
fprintf('Output directory: %s\n', output_dir);

bellhop_exe = local_find_bellhop_executable(cfg.bellhop_exe);
fprintf('Bellhop executable: %s\n', bellhop_exe);

geometry = local_geometry(cfg);
case_root = fullfile(output_dir, 'flat_surface_case');
env_file = [case_root '.env'];
local_write_bellhop_env(env_file, cfg, geometry);
local_run_bellhop(bellhop_exe, case_root);
bellhop_all = local_read_bellhop_arrivals([case_root '.arr']);
bellhop = local_select_stage1_arrivals(bellhop_all);

if height(bellhop) ~= 2
    error('validate_pe_bellhop_flat_surface_vertical:BellhopPathCount', ...
        'Expected direct and one surface arrival, but selected %d paths.', ...
        height(bellhop));
end

params_scalar = local_pe_params(cfg, cfg.center_frequency_hz);
params_scalar.enable_surface_reflection = false;
pe_scalar_direct_only = vertical_channel_model(params_scalar);

params_scalar.enable_surface_reflection = true;
pe_scalar_flat = vertical_channel_model(params_scalar);

f_axis_hz = linspace(cfg.frequency_band_hz(1), ...
    cfg.frequency_band_hz(2), cfg.frequency_count).';
params_wideband = local_pe_params(cfg, f_axis_hz.');
params_wideband.enable_surface_reflection = true;
params_wideband.f_ref_hz = cfg.center_frequency_hz;
pe_wideband = vertical_channel_model(params_wideband);

invariants = local_check_pe_invariants( ...
    pe_scalar_direct_only, pe_scalar_flat, pe_wideband, cfg);

pe_paths = local_build_pe_path_table( ...
    pe_scalar_flat, pe_wideband, geometry);
bellhop_paths = local_build_bellhop_path_table(bellhop);
direct_scale = bellhop_paths.amplitude(1) / pe_paths.amplitude(1);
[pdp, pdp_metrics] = local_build_pdp_comparison( ...
    pe_wideband, bellhop_paths, geometry, direct_scale, cfg);
pe_paths.arrival_time_s = pdp_metrics.pe_peak_time_s;
comparison_table = local_build_comparison_table(pe_paths, bellhop_paths);
comparison_table.pe_direct_calibrated_amplitude = ...
    comparison_table.pe_amplitude .* direct_scale;
comparison_table.pe_direct_calibrated_tl_db = -20 * log10(max( ...
    comparison_table.pe_direct_calibrated_amplitude, realmin));
comparison_table.calibrated_tl_difference_db = ...
    comparison_table.pe_direct_calibrated_tl_db - comparison_table.bellhop_tl_db;

checks = local_build_checks(comparison_table, invariants, pdp_metrics, cfg);
if ~all(checks.passed)
    failed = strjoin(cellstr(checks.check_name(~checks.passed)), ', ');
    error('validate_pe_bellhop_flat_surface_vertical:ValidationFailed', ...
        'One or more validation checks failed: %s', failed);
end

comparison_csv = fullfile(output_dir, 'arrival_time_comparison.csv');
figure_file = fullfile(output_dir, 'pe_bellhop_pdp_comparison.png');
mat_file = fullfile(output_dir, 'pe_bellhop_flat_surface_validation.mat');
report_file = fullfile(output_dir, 'pe_bellhop_flat_surface_report.md');
writetable(comparison_table, comparison_csv);
local_plot_comparison(pdp, comparison_table, figure_file);

validation = struct();
validation.config = cfg;
validation.geometry = geometry;
validation.bellhop_executable = bellhop_exe;
validation.bellhop_all_arrivals = bellhop_all;
validation.bellhop_selected_arrivals = bellhop;
validation.pe_scalar_direct_only = local_compact_channel(pe_scalar_direct_only);
validation.pe_scalar_flat = local_compact_channel(pe_scalar_flat);
validation.pe_wideband = local_compact_channel(pe_wideband);
validation.invariants = invariants;
validation.comparison_table = comparison_table;
validation.pdp = pdp;
validation.pdp_metrics = pdp_metrics;
validation.checks = checks;
validation.direct_amplitude_calibration = direct_scale;
validation.files = struct('env', env_file, 'arr', [case_root '.arr'], ...
    'prt', [case_root '.prt'], 'comparison_csv', comparison_csv, ...
    'figure', figure_file, 'mat', mat_file, 'report', report_file);
save(mat_file, 'validation');
local_write_report(report_file, validation);

disp(comparison_table(:, {'path_name', 'pe_arrival_time_ms', ...
    'bellhop_arrival_time_ms', 'arrival_time_difference_ms', ...
    'pe_tl_db', 'bellhop_tl_db', 'calibrated_tl_difference_db'}));
disp(checks);
fprintf('Validation passed. Report: %s\n', report_file);
end

function cfg = local_default_config()
cfg = struct();
cfg.water_depth_m = 100;
cfg.sound_speed_mps = 1500;
cfg.tx_xyz_m = [0, 0, 80];
cfg.rx_xyz_m = [5, 0, 10];
cfg.center_frequency_hz = 4000;
cfg.frequency_band_hz = [3000, 5000];
cfg.frequency_count = 65;
cfg.pe_nx = 128;
cfg.pe_ny = 128;
cfg.pe_x_width_m = 32;
cfg.pe_y_width_m = 32;
cfg.pe_stepz_lamb = 0.5;
cfg.pe_sigma_source_m = 0.3;
cfg.pe_sponge_ratio = 0.12;
cfg.pe_alpha_max_np_per_m = 0.15;
cfg.bellhop_beam_count = 5001;
cfg.bellhop_angle_margin_deg = 1.0;
cfg.bellhop_exe = '';
cfg.cir_window = 'hann';
cfg.cir_zero_padding_factor = 8;
cfg.arrival_time_tolerance_ms = 0.75;
cfg.relative_path_tl_tolerance_db = 3.0;
cfg.raw_tl_tolerance_db = 6.0;
cfg.calibrated_reflection_tl_tolerance_db = 3.0;
cfg.pdp_peak_tolerance_ms = 0.75;
cfg.invariant_tolerance = 1e-10;
end

function cfg = local_apply_overrides(cfg, overrides)
if ~isstruct(overrides)
    error('overrides must be a struct.');
end
names = fieldnames(overrides);
for ii = 1:numel(names)
    if ~isfield(cfg, names{ii})
        error('Unknown validation override: %s', names{ii});
    end
    cfg.(names{ii}) = overrides.(names{ii});
end
end

function local_validate_config(cfg)
if cfg.water_depth_m <= 0 || cfg.sound_speed_mps <= 0
    error('water_depth_m and sound_speed_mps must be positive.');
end
if numel(cfg.tx_xyz_m) ~= 3 || numel(cfg.rx_xyz_m) ~= 3
    error('tx_xyz_m and rx_xyz_m must each contain [x,y,z].');
end
if cfg.tx_xyz_m(3) <= cfg.rx_xyz_m(3) || cfg.rx_xyz_m(3) < 0 || ...
        cfg.tx_xyz_m(3) > cfg.water_depth_m
    error('Require 0 <= z_rx < z_tx <= water_depth_m.');
end
if hypot(cfg.rx_xyz_m(1) - cfg.tx_xyz_m(1), ...
        cfg.rx_xyz_m(2) - cfg.tx_xyz_m(2)) <= 0
    error(['Bellhop stage-1 validation requires a nonzero horizontal offset; ', ...
        'use a small near-vertical offset instead of the singular zero-range geometry.']);
end
if cfg.frequency_count < 3 || mod(cfg.frequency_count, 2) ~= 1
    error('frequency_count must be an odd integer >= 3 so f_ref is sampled.');
end
if abs(mean(cfg.frequency_band_hz) - cfg.center_frequency_hz) > 1e-9
    error('center_frequency_hz must equal the midpoint of frequency_band_hz.');
end
end

function geometry = local_geometry(cfg)
horizontal_range_m = hypot(cfg.rx_xyz_m(1) - cfg.tx_xyz_m(1), ...
    cfg.rx_xyz_m(2) - cfg.tx_xyz_m(2));
direct_length_m = hypot(horizontal_range_m, ...
    cfg.tx_xyz_m(3) - cfg.rx_xyz_m(3));
surface_length_m = hypot(horizontal_range_m, ...
    cfg.tx_xyz_m(3) + cfg.rx_xyz_m(3));
geometry = struct();
geometry.horizontal_range_m = horizontal_range_m;
geometry.direct_length_m = direct_length_m;
geometry.surface_length_m = surface_length_m;
geometry.direct_time_s = direct_length_m / cfg.sound_speed_mps;
geometry.surface_time_s = surface_length_m / cfg.sound_speed_mps;
geometry.excess_surface_delay_s = ...
    geometry.surface_time_s - geometry.direct_time_s;
geometry.direct_launch_angle_deg = -atan2d( ...
    cfg.tx_xyz_m(3) - cfg.rx_xyz_m(3), horizontal_range_m);
geometry.surface_launch_angle_deg = -atan2d( ...
    cfg.tx_xyz_m(3) + cfg.rx_xyz_m(3), horizontal_range_m);
end

function paramsV = local_pe_params(cfg, frequency_hz)
paramsV = struct();
paramsV.f0 = frequency_hz;
paramsV.enable_wideband = false;
paramsV.c0 = cfg.sound_speed_mps;
paramsV.z_max = cfg.water_depth_m;
paramsV.stepz_lamb = cfg.pe_stepz_lamb;
paramsV.xw = cfg.pe_x_width_m;
paramsV.yw = cfg.pe_y_width_m;
paramsV.nx = cfg.pe_nx;
paramsV.ny = cfg.pe_ny;
paramsV.x_tx = cfg.tx_xyz_m(1);
paramsV.y_tx = cfg.tx_xyz_m(2);
paramsV.z_tx = cfg.tx_xyz_m(3);
paramsV.x_rx = cfg.rx_xyz_m(1);
paramsV.y_rx = cfg.rx_xyz_m(2);
paramsV.z_rx = cfg.rx_xyz_m(3);
paramsV.sigma_src_m = cfg.pe_sigma_source_m;
paramsV.sponge_ratio = cfg.pe_sponge_ratio;
paramsV.alpha_max_np_per_m = cfg.pe_alpha_max_np_per_m;
paramsV.env_mode = 'uniform';
paramsV.show_figures = false;
paramsV.enforce_1_over_R = true;
paramsV.enable_surface_reflection = true;
paramsV.surface_reflect_coeff = -1;
paramsV.surface_phase_mode = 'normal';
paramsV.surface_boundary_model = 'kirchhoff_spatial';
paramsV.channel_phase_reference = 'direct_dsp';
paramsV.sea_hs_target = 0;
paramsV.surface_roughness_scale_mode = 'target_hs';
paramsV.enable_bubbles = false;
paramsV.save_mode = 'rx_only';
paramsV.use_gpu = false;
end

function exe = local_find_bellhop_executable(configured_exe)
candidates = {};
if ~isempty(configured_exe)
    candidates{end + 1} = char(configured_exe);
end
env_exe = getenv('BELLHOP_EXE');
if ~isempty(env_exe)
    candidates{end + 1} = env_exe;
end
matlab_exe = which('bellhop.exe');
if ~isempty(matlab_exe)
    candidates{end + 1} = matlab_exe;
end
for ii = 1:numel(candidates)
    if exist(candidates{ii}, 'file') == 2
        exe = candidates{ii};
        return
    end
end
if ispc
    [status, found] = system('where bellhop.exe');
else
    [status, found] = system('which bellhop');
end
if status == 0
    lines = regexp(strtrim(found), '\r?\n', 'split');
    if ~isempty(lines) && exist(strtrim(lines{1}), 'file') == 2
        exe = strtrim(lines{1});
        return
    end
end
error('validate_pe_bellhop_flat_surface_vertical:BellhopMissing', ...
    ['Bellhop executable not found. Install the Acoustics Toolbox and put ', ...
    'bellhop.exe on PATH, or set the BELLHOP_EXE environment variable.']);
end

function local_write_bellhop_env(env_file, cfg, geometry)
fid = fopen(env_file, 'w');
if fid < 0
    error('Cannot create Bellhop environment file: %s', env_file);
end
cleanup = onCleanup(@() fclose(fid));
angles = sort([geometry.direct_launch_angle_deg, ...
    geometry.surface_launch_angle_deg]);
angle_min = max(-89.9, angles(1) - cfg.bellhop_angle_margin_deg);
angle_max = min(89.9, angles(2) + cfg.bellhop_angle_margin_deg);
fprintf(fid, '''PE-Bellhop flat pressure-release surface validation''\n');
fprintf(fid, '%.12g\n', cfg.center_frequency_hz);
fprintf(fid, '1\n');
fprintf(fid, '''CVW''\n');
fprintf(fid, '2 0.0 %.12g\n', cfg.water_depth_m);
fprintf(fid, '0.0 %.12g /\n', cfg.sound_speed_mps);
fprintf(fid, '%.12g %.12g /\n', cfg.water_depth_m, cfg.sound_speed_mps);
fprintf(fid, '''A'' 0.0\n');
fprintf(fid, '%.12g 1800.0 0.0 2.0 0.0 /\n', cfg.water_depth_m);
fprintf(fid, '1\n%.12g /\n', cfg.tx_xyz_m(3));
fprintf(fid, '1\n%.12g /\n', cfg.rx_xyz_m(3));
fprintf(fid, '1\n%.12g /\n', geometry.horizontal_range_m / 1000);
fprintf(fid, '''A''\n');
fprintf(fid, '%d\n', cfg.bellhop_beam_count);
fprintf(fid, '%.12g %.12g /\n', angle_min, angle_max);
fprintf(fid, '0.0 %.12g %.12g\n', 1.1 * cfg.water_depth_m, ...
    2 * geometry.horizontal_range_m / 1000);
clear cleanup
end

function local_run_bellhop(exe, case_root)
case_dir = fileparts(case_root);
[~, case_name] = fileparts(case_root);
old_dir = pwd;
cleanup = onCleanup(@() cd(old_dir));
cd(case_dir);
if exist([case_name '.arr'], 'file') == 2
    delete([case_name '.arr']);
end
if exist([case_name '.prt'], 'file') == 2
    delete([case_name '.prt']);
end
command = sprintf('"%s" "%s"', exe, case_name);
[status, output] = system(command);
if status ~= 0
    error('Bellhop failed with status %d:\n%s', status, output);
end
required = {[case_root '.arr'], [case_root '.prt']};
for ii = 1:numel(required)
    if exist(required{ii}, 'file') ~= 2
        error('Bellhop did not create expected output: %s', required{ii});
    end
end
clear cleanup
end

function arrivals = local_read_bellhop_arrivals(arr_file)
fid = fopen(arr_file, 'r');
if fid < 0
    error('Cannot open Bellhop arrivals file: %s', arr_file);
end
cleanup = onCleanup(@() fclose(fid));
frequency_hz = fscanf(fid, '%f', 1);
nsd = fscanf(fid, '%d', 1);
nrd = fscanf(fid, '%d', 1);
nrr = fscanf(fid, '%d', 1);
source_depth_m = fscanf(fid, '%f', nsd);
receiver_depth_m = fscanf(fid, '%f', nrd);
receiver_range_m = 1000 * fscanf(fid, '%f', nrr);
if nsd ~= 1 || nrd ~= 1 || nrr ~= 1
    error('Stage-1 parser expects one source depth, receiver depth, and range.');
end
fscanf(fid, '%d', 1); % maximum number of arrivals for this source
narr = fscanf(fid, '%d', 1);
raw = fscanf(fid, '%f', [8, narr]);
if size(raw, 2) ~= narr
    error('Bellhop arrivals file ended before all arrivals were read.');
end
amplitude_complex = raw(1, :).' .* exp(1i * deg2rad(raw(2, :).'));
arrivals = table((1:narr).', repmat(frequency_hz, narr, 1), ...
    repmat(source_depth_m, narr, 1), repmat(receiver_depth_m, narr, 1), ...
    repmat(receiver_range_m, narr, 1), amplitude_complex, abs(amplitude_complex), ...
    raw(2, :).', raw(3, :).', raw(4, :).', raw(5, :).', raw(6, :).', ...
    raw(7, :).', raw(8, :).', ...
    'VariableNames', {'arrival_index', 'frequency_hz', 'source_depth_m', ...
    'receiver_depth_m', 'receiver_range_m', 'amplitude_complex', 'amplitude', ...
    'phase_deg', 'delay_s', 'delay_imag_s', 'source_angle_deg', ...
    'receiver_angle_deg', 'top_bounce_count', 'bottom_bounce_count'});
clear cleanup
end

function selected = local_select_stage1_arrivals(arrivals)
mask = arrivals.bottom_bounce_count == 0 & ...
    arrivals.top_bounce_count <= 1;
selected = sortrows(arrivals(mask, :), 'delay_s');
selected.path_name = strings(height(selected), 1);
selected.path_name(selected.top_bounce_count == 0) = "direct";
selected.path_name(selected.top_bounce_count == 1) = "surface_reflection";
selected = movevars(selected, 'path_name', 'Before', 1);
end

function invariants = local_check_pe_invariants(direct_only, flat, wideband, cfg)
invariants = struct();
invariants.direct_only_reflection_max_abs = max(abs(direct_only.H_reflect_f(:)));
invariants.scalar_direct_disabled_difference_abs = ...
    abs(direct_only.h_direct - flat.h_direct);
invariants.scalar_sum_error_abs = max(abs( ...
    flat.H_f(:) - flat.H_direct_f(:) - flat.H_reflect_f(:)));
invariants.wideband_sum_error_abs = max(abs( ...
    wideband.H_f(:) - wideband.H_direct_f(:) - wideband.H_reflect_f(:)));
invariants.scalar_1_over_r_passed = logical(flat.pass_1_over_R);
invariants.wideband_1_over_r_passed = logical(wideband.pass_1_over_R);
if max([invariants.direct_only_reflection_max_abs, ...
        invariants.scalar_direct_disabled_difference_abs, ...
        invariants.scalar_sum_error_abs, ...
        invariants.wideband_sum_error_abs]) > cfg.invariant_tolerance
    error('PE public-interface invariant failed.');
end
if ~invariants.scalar_1_over_r_passed || ~invariants.wideband_1_over_r_passed
    error('PE 1/R validation failed.');
end
end

function paths = local_build_pe_path_table(scalar, wideband, geometry)
H_direct_physical_f = wideband.H_direct_physical_f(:);
H_surface_physical_f = wideband.H_reflect_physical_f(:);
paths = table(["direct"; "surface_reflection"], ...
    [geometry.direct_time_s; geometry.surface_time_s], ...
    [abs(scalar.h_direct); abs(scalar.h_reflect)], ...
    [-20 * log10(max(abs(scalar.h_direct), realmin)); ...
     -20 * log10(max(abs(scalar.h_reflect), realmin))], ...
    {H_direct_physical_f; H_surface_physical_f}, ...
    'VariableNames', {'path_name', 'arrival_time_s', 'amplitude', 'tl_db', ...
    'physical_frequency_response'});
end

function paths = local_build_bellhop_path_table(bellhop)
paths = table(bellhop.path_name, bellhop.delay_s, bellhop.amplitude, ...
    -20 * log10(max(bellhop.amplitude, realmin)), bellhop.amplitude_complex, ...
    bellhop.top_bounce_count, bellhop.bottom_bounce_count, ...
    'VariableNames', {'path_name', 'arrival_time_s', 'amplitude', 'tl_db', ...
    'amplitude_complex', 'top_bounce_count', 'bottom_bounce_count'});
end

function comparison = local_build_comparison_table(pe, bellhop)
if ~isequal(pe.path_name, bellhop.path_name)
    error('PE and Bellhop path labels do not match.');
end
comparison = table(pe.path_name, 1000 * pe.arrival_time_s, ...
    1000 * bellhop.arrival_time_s, ...
    1000 * (pe.arrival_time_s - bellhop.arrival_time_s), ...
    pe.amplitude, bellhop.amplitude, pe.tl_db, bellhop.tl_db, ...
    pe.tl_db - bellhop.tl_db, ...
    'VariableNames', {'path_name', 'pe_arrival_time_ms', ...
    'bellhop_arrival_time_ms', 'arrival_time_difference_ms', ...
    'pe_amplitude', 'bellhop_amplitude', 'pe_tl_db', 'bellhop_tl_db', ...
    'raw_tl_difference_db'});
end

function [pdp, metrics] = local_build_pdp_comparison( ...
    pe, bellhop_paths, geometry, direct_scale, cfg)
f = pe.f_axis(:);
H_pe_f = direct_scale * pe.H_physical_f(:);
H_bellhop_f = complex(zeros(size(f)));
for ii = 1:height(bellhop_paths)
    H_bellhop_f = H_bellhop_f + bellhop_paths.amplitude_complex(ii) .* ...
        exp(1i * 2 * pi * f * bellhop_paths.arrival_time_s(ii));
end
reference_delay_s = 0;
cir_options=struct('input_reference','absolute_physical', ...
    'window_type',cfg.cir_window,'zero_padding_factor',cfg.cir_zero_padding_factor);
pe_cir = build_channel_cir_vertical(H_pe_f,f,cir_options);
bellhop_cir = build_channel_cir_vertical(H_bellhop_f,f,cir_options);
delay_s = pe_cir.delay_axis_s;
pe_power = abs(pe_cir.h_physical_tau).^2;
bellhop_power = abs(bellhop_cir.h_physical_tau).^2;
pdp = struct();
pdp.delay_s = delay_s;
pdp.pe_h = pe_cir.h_physical_tau;
pdp.bellhop_h = bellhop_cir.h_physical_tau;
pdp.pe_power = pe_power;
pdp.bellhop_power = bellhop_power;
pdp.pe_power_normalized = pe_power / max(pe_power);
pdp.bellhop_power_normalized = bellhop_power / max(bellhop_power);
pdp.frequency_axis_hz = f;
pdp.reference_delay_s = reference_delay_s;
pdp.window = cfg.cir_window;
pdp.zero_padding_factor = cfg.cir_zero_padding_factor;
pdp.physical_delay_resolution_s = pe_cir.physical_delay_resolution_s;

metrics = struct();
metrics.pe_peak_time_s = local_two_peak_times(delay_s, pe_power, geometry);
metrics.bellhop_peak_time_s = local_two_peak_times(delay_s, bellhop_power, geometry);
metrics.peak_time_difference_ms = 1000 * ( ...
    metrics.pe_peak_time_s - metrics.bellhop_peak_time_s);
end

function peak_times = local_two_peak_times(delay_s, power, geometry)
targets = [geometry.direct_time_s; geometry.surface_time_s];
half_window_s = 0.35 * geometry.excess_surface_delay_s;
peak_times = NaN(2, 1);
for ii = 1:2
    mask = abs(delay_s - targets(ii)) <= half_window_s;
    indices = find(mask);
    [~, local_index] = max(power(mask));
    peak_times(ii) = delay_s(indices(local_index));
end
end

function checks = local_build_checks(comparison, invariants, pdp_metrics, cfg)
relative_pe_db = comparison.pe_tl_db(2) - comparison.pe_tl_db(1);
relative_bh_db = comparison.bellhop_tl_db(2) - comparison.bellhop_tl_db(1);
names = ["path_count"; "arrival_time"; "raw_tl"; "relative_path_tl"; ...
    "direct_calibrated_reflection_tl"; "pdp_peak_time"; ...
    "pe_frequency_sum_invariant"; "pe_direct_only_invariant"; "pe_1_over_r"];
values = [height(comparison); max(abs(comparison.arrival_time_difference_ms)); ...
    max(abs(comparison.raw_tl_difference_db)); abs(relative_pe_db - relative_bh_db); ...
    abs(comparison.calibrated_tl_difference_db(2)); ...
    max(abs(pdp_metrics.peak_time_difference_ms)); ...
    invariants.wideband_sum_error_abs; ...
    max(invariants.direct_only_reflection_max_abs, ...
        invariants.scalar_direct_disabled_difference_abs); ...
    double(~(invariants.scalar_1_over_r_passed && invariants.wideband_1_over_r_passed))];
limits = [2; cfg.arrival_time_tolerance_ms; cfg.raw_tl_tolerance_db; ...
    cfg.relative_path_tl_tolerance_db; cfg.calibrated_reflection_tl_tolerance_db; ...
    cfg.pdp_peak_tolerance_ms; cfg.invariant_tolerance; ...
    cfg.invariant_tolerance; 0];
relations = ["=="; "<="; "<="; "<="; "<="; "<="; "<="; "<="; "=="];
passed = (relations == "==" & values == limits) | ...
    (relations == "<=" & values <= limits);
checks = table(names, values, relations, limits, passed, ...
    'VariableNames', {'check_name', 'value', 'relation', 'limit', 'passed'});
end

function compact = local_compact_channel(channel)
compact = struct();
compact.h_direct = channel.h_direct;
compact.h_reflect = channel.h_reflect;
compact.h_total = channel.h_total;
compact.H_direct_f = channel.H_direct_f;
compact.H_reflect_f = channel.H_reflect_f;
compact.H_f = channel.H_f;
compact.H_direct_reduced_f = channel.H_direct_reduced_f;
compact.H_reflect_reduced_f = channel.H_reflect_reduced_f;
compact.H_direct_physical_f = channel.H_direct_physical_f;
compact.H_reflect_physical_f = channel.H_reflect_physical_f;
compact.H_physical_f = channel.H_physical_f;
compact.phase_reference_meta = channel.phase_reference_meta;
compact.f_axis = channel.f_axis;
compact.idx_f_ref = channel.idx_f_ref;
compact.path_loss_db = channel.path_loss_db;
compact.pass_1_over_R = channel.pass_1_over_R;
compact.fit_slope = channel.fit_slope;
compact.fit_err_rms = channel.fit_err_rms;
compact.config = channel.config;
compact.roughness_meta = channel.roughness_meta;
end

function local_plot_comparison(pdp, comparison, figure_file)
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 900 720]);
cleanup = onCleanup(@() close(fig));
subplot(2, 1, 1);
pe_db = 10 * log10(max(pdp.pe_power_normalized, 1e-8));
bh_db = 10 * log10(max(pdp.bellhop_power_normalized, 1e-8));
plot(1000 * pdp.delay_s, pe_db, 'b-', 'LineWidth', 1.5);
hold on;
plot(1000 * pdp.delay_s, bh_db, 'r--', 'LineWidth', 1.5);
grid on;
xlim([min(comparison.bellhop_arrival_time_ms) - 2, ...
    max(comparison.bellhop_arrival_time_ms) + 2]);
ylim([-60, 1]);
xlabel('Absolute delay (ms)');
ylabel('Normalized PDP (dB)');
title('PE and Bellhop flat-surface PDP');
legend('PE (carrier-restored)', 'Bellhop arrivals', 'Location', 'best');

subplot(2, 1, 2);
bar(categorical(comparison.path_name), ...
    [comparison.pe_tl_db, comparison.bellhop_tl_db]);
grid on;
ylabel('Raw path transmission loss (dB)');
title('Path amplitude / TL comparison at center frequency');
legend('PE', 'Bellhop', 'Location', 'best');
exportgraphics(fig, figure_file, 'Resolution', 180);
clear cleanup
end

function local_write_report(report_file, validation)
cfg = validation.config;
g = validation.geometry;
t = validation.comparison_table;
c = validation.checks;
fid = fopen(report_file, 'w', 'n', 'UTF-8');
if fid < 0
    error('Cannot create validation report: %s', report_file);
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '# PE 与 Bellhop 第一阶段平面海面交叉验证报告\n\n');
fprintf(fid, '## 结论\n\n');
fprintf(fid, ['该精简案例通过全部自动检查。PE 与 Bellhop 均识别出直达和一次海面反射两条主要路径；', ...
    '到达时间、相对路径强度及 PDP 峰位一致。绝对幅度存在近似常数偏差，主要来自当前 PE 高斯初场', ...
    '与 Bellhop 单位点源归一化不同，而不是新增的传播路径或粗糙面效应。\n\n']);
fprintf(fid, '## 共同环境与运行设置\n\n');
fprintf(fid, '- 水深：`%.6g m`；均匀声速：`%.6g m/s`。\n', ...
    cfg.water_depth_m, cfg.sound_speed_mps);
fprintf(fid, '- 发射位置：`[%.6g, %.6g, %.6g] m`；接收位置：`[%.6g, %.6g, %.6g] m`。\n', ...
    cfg.tx_xyz_m, cfg.rx_xyz_m);
fprintf(fid, '- 中心频率：`%.6g Hz`；PE 宽带：`%.6g--%.6g Hz`，`%d` 个等间隔频点。\n', ...
    cfg.center_frequency_hz, cfg.frequency_band_hz, cfg.frequency_count);
fprintf(fid, ['- 海面：平面压力释放边界，`surface_reflect_coeff=-1`、`sea_hs_target=0`；', ...
    '未启用 SSA、随机海面、Kirchhoff 粗糙相位或气泡。\n']);
fprintf(fid, '- PE 网格：`%d x %d`，横向窗口 `%.6g m x %.6g m`，`stepz_lamb=%.6g`。\n', ...
    cfg.pe_nx, cfg.pe_ny, cfg.pe_x_width_m, cfg.pe_y_width_m, cfg.pe_stepz_lamb);
fprintf(fid, '- Bellhop：标准 ASCII arrivals 模式，`%d` 条射线，只覆盖直达/一次海面反射所需角扇区。\n', ...
    cfg.bellhop_beam_count);
fprintf(fid, ['- Bellhop 依赖：外部 Acoustic Toolbox `bellhop.exe`；本次使用官方托管的 ', ...
    '[Windows atWin.zip](https://oalib-acoustics.org/website_resources/AcousticsToolbox/versions/atWin.zip) ', ...
    '分发包，第三方二进制未复制进项目。\n']);
fprintf(fid, '- 几何路径长度：直达 `%.6g m`，海面镜像 `%.6g m`，理论超时延 `%.6g ms`。\n\n', ...
    g.direct_length_m, g.surface_length_m, 1000 * g.excess_surface_delay_s);

fprintf(fid, '## 到达与幅度结果\n\n');
fprintf(fid, '| 路径 | PE 到达 (ms) | Bellhop 到达 (ms) | 差值 (ms) | PE TL (dB) | Bellhop TL (dB) | 直达校准后 TL 差 (dB) |\n');
fprintf(fid, '|---|---:|---:|---:|---:|---:|---:|\n');
for ii = 1:height(t)
    fprintf(fid, '| %s | %.6f | %.6f | %+.6f | %.4f | %.4f | %+.4f |\n', ...
        char(t.path_name(ii)), t.pe_arrival_time_ms(ii), ...
        t.bellhop_arrival_time_ms(ii), t.arrival_time_difference_ms(ii), ...
        t.pe_tl_db(ii), t.bellhop_tl_db(ii), ...
        t.calibrated_tl_difference_db(ii));
end
fprintf(fid, '\nPE 原始幅度相对 Bellhop 的直达标定系数为 `%.8g`。', ...
    validation.direct_amplitude_calibration);
fprintf(fid, ['标定只用于区分源归一化常数与路径相对衰减，不回写 PE 公共输出，也不用于', ...
    '掩盖路径间差异。\n\n']);

fprintf(fid, '## 自动检查\n\n');
fprintf(fid, '| 检查 | 数值 | 条件 | 阈值 | 通过 |\n');
fprintf(fid, '|---|---:|:---:|---:|:---:|\n');
for ii = 1:height(c)
    fprintf(fid, '| %s | %.8g | %s | %.8g | %d |\n', ...
        char(c.check_name(ii)), c.value(ii), char(c.relation(ii)), ...
        c.limit(ii), c.passed(ii));
end
fprintf(fid, '\n![PE 与 Bellhop PDP/TL 对比](pe_bellhop_pdp_comparison.png)\n\n');

fprintf(fid, '## 差异解释与限制\n\n');
fprintf(fid, ['1. 当前 PE 以有限宽度高斯场启动，Bellhop arrivals 使用单位点源的几何扩展归一化；', ...
    '因此原始 TL 可出现近似常数偏置。更有辨识力的是直达校准后的反射残差及反射/直达相对 TL。\n']);
fprintf(fid, ['2. PE 内部推进量是去除纵向载波的 reduced envelope；当前验证直接读取公共 ', ...
    '`H_*_physical_f` 字段。在 `exp(-i*omega*t)` 约定下 physical 字段恢复 ', ...
    '`exp(+i 2 pi f tau_0)`，再用 FFT 提取正时延；不得在验证层重复补载波。\n']);
fprintf(fid, ['3. Bellhop 为射线/高斯波束模型，PE 为波动模型；有限频带、有限网格、波束插值和', ...
    '绕射会产生小的幅相差异，不要求逐点完全一致。\n']);
fprintf(fid, ['4. 为避开 Bellhop 在零水平距离的退化几何，接收机设置 5 m 水平偏移；', ...
    '该传播仍为近垂直上行。阶段一不包含海底反射、粗糙海面、SSA、随机信道或通信调制。\n\n']);
fprintf(fid, ['5. 总 PDP 峰间的 PE 旁瓣来自有限带宽重构和 PE 包络随频率的幅相变化；', ...
    '主要路径数由分离的直达/反射分量和 Bellhop bounce count 判定，不能把旁瓣计为额外本征声线。\n\n']);

fprintf(fid, '## 生成文件\n\n');
fprintf(fid, '- `flat_surface_case.env` / `.arr` / `.prt`：Bellhop 输入及原始输出。\n');
fprintf(fid, '- `arrival_time_comparison.csv`：路径对比表。\n');
fprintf(fid, '- `pe_bellhop_pdp_comparison.png`：PDP 和路径 TL 图。\n');
fprintf(fid, '- `pe_bellhop_flat_surface_validation.mat`：完整可复查结果。\n\n');
fprintf(fid, '运行入口：`validate_pe_bellhop_flat_surface_vertical`。\n');
clear cleanup
end
