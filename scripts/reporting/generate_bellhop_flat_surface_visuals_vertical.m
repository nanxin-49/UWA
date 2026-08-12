function output = generate_bellhop_flat_surface_visuals_vertical(overrides)
%GENERATE_BELLHOP_FLAT_SURFACE_VISUALS_VERTICAL Render Bellhop R/C/I results.
% Validation/reporting-only entrypoint. It reuses the accepted environment
% and receiver results from the PE/Bellhop matrix without changing the PE
% propagator, public channel outputs, or communication chain.

if nargin < 1 || isempty(overrides), overrides = struct(); end
this_file = mfilename('fullpath');
project_root = fileparts(fileparts(fileparts(this_file)));
addpath(project_root);

default_source = fullfile(project_root, 'results', 'validation', ...
    'pe_bellhop_flat_surface_matrix', 'pe_bellhop_flat_surface_matrix.mat');
source_file = default_source;
if isfield(overrides, 'source_validation_mat')
    source_file = char(overrides.source_validation_mat);
end
if exist(source_file, 'file') ~= 2
    error('Bellhop matrix result is missing: %s', source_file);
end
loaded = load(source_file, 'validation');
if ~isfield(loaded, 'validation')
    error('source_validation_mat must contain validation.');
end
if isfield(loaded.validation,'matrix_validation')
    matrix_validation = loaded.validation.matrix_validation;
else
    matrix_validation = loaded.validation;
end

cfg = local_defaults(project_root, source_file, matrix_validation);
cfg = local_overrides(cfg, overrides);
local_validate_cfg(cfg, matrix_validation);
if ~exist(cfg.output_dir, 'dir'), mkdir(cfg.output_dir); end

bellhop_exe = local_find_bellhop(cfg.bellhop_exe, matrix_validation);
fprintf('Bellhop flat-surface R/C/I visualization\nBellhop: %s\n', bellhop_exe);
fprintf('Output: %s\n', cfg.output_dir);

% A sparse ray fan is intentionally separate from the dense field runs.
ray_root = fullfile(cfg.output_dir, 'bellhop_visual_ray');
local_write_env([ray_root '.env'], cfg, 'R', cfg.ray_beam_count, false);
local_run_bellhop(bellhop_exe, ray_root, '.ray');
ray_data = local_read_ray([ray_root '.ray']);

modes = ["C", "I"];
field_runs = repmat(struct('mode', "", 'beam_count', 0, 'root', '', ...
    'shade', []), numel(modes), numel(cfg.field_beam_counts));
for mm = 1:numel(modes)
    for bb = 1:numel(cfg.field_beam_counts)
        beam_count = cfg.field_beam_counts(bb);
        label = local_mode_label(modes(mm));
        root = fullfile(cfg.output_dir, sprintf('bellhop_visual_%s_%05d', ...
            label, beam_count));
        fprintf('Bellhop %s field: %d beams\n', label, beam_count);
        local_write_env([root '.env'], cfg, char(modes(mm)), beam_count, true);
        local_run_bellhop(bellhop_exe, root, '.shd');
        field_runs(mm, bb) = struct('mode', modes(mm), ...
            'beam_count', beam_count, 'root', root, ...
            'shade', local_read_shd_2d([root '.shd']));
    end
end

low_index = 1;
high_index = numel(cfg.field_beam_counts);
coherent_low = field_runs(1, low_index).shade;
coherent_high = field_runs(1, high_index).shade;
incoherent_low = field_runs(2, low_index).shade;
incoherent_high = field_runs(2, high_index).shade;

[coherent_tl_db, source_mask] = local_tl_field(coherent_high, cfg);
[incoherent_tl_db, ~] = local_tl_field(incoherent_high, cfg);
[coherent_low_tl_db, ~] = local_tl_field(coherent_low, cfg);
[incoherent_low_tl_db, ~] = local_tl_field(incoherent_low, cfg);

convergence_table = local_convergence_table(coherent_low_tl_db, ...
    coherent_tl_db, incoherent_low_tl_db, incoherent_tl_db, source_mask, cfg);
receiver_table = local_receiver_table(coherent_high, incoherent_high, ...
    matrix_validation, cfg);
saved_invariant_error = local_saved_invariant_error(matrix_validation);
checks = local_checks(ray_data, coherent_high, incoherent_high, ...
    convergence_table, receiver_table, saved_invariant_error, cfg);
passed = all(checks.passed);

ray_figure = fullfile(cfg.output_dir, 'bellhop_ray_geometry.png');
field_figure = fullfile(cfg.output_dir, 'bellhop_tl_fields.png');
slice_figure = fullfile(cfg.output_dir, 'bellhop_receiver_depth_tl.png');
convergence_csv = fullfile(cfg.output_dir, 'bellhop_tl_convergence.csv');
receiver_csv = fullfile(cfg.output_dir, 'bellhop_receiver_tl_comparison.csv');
mat_file = fullfile(cfg.output_dir, 'bellhop_flat_surface_visualization.mat');
report_file = fullfile(cfg.output_dir, 'bellhop_flat_surface_visual_report.md');

local_plot_rays(ray_data, cfg, ray_figure);
local_plot_fields(coherent_high, coherent_tl_db, incoherent_tl_db, ...
    cfg, field_figure);
local_plot_slice(coherent_high, coherent_tl_db, incoherent_tl_db, ...
    receiver_table, cfg, slice_figure);
writetable(convergence_table, convergence_csv);
writetable(receiver_table, receiver_csv);

output = struct();
output.config = cfg;
output.source_matrix_passed = logical(matrix_validation.passed);
output.source_validation_mat = source_file;
output.bellhop_executable = bellhop_exe;
output.ray_data = ray_data;
output.field_runs = field_runs;
output.coherent_tl_db = coherent_tl_db;
output.incoherent_tl_db = incoherent_tl_db;
output.source_mask = source_mask;
output.convergence_table = convergence_table;
output.receiver_table = receiver_table;
output.saved_pe_invariant_error = saved_invariant_error;
output.checks = checks;
output.passed = passed;
output.files = struct('ray_figure', ray_figure, ...
    'field_figure', field_figure, 'slice_figure', slice_figure, ...
    'convergence_csv', convergence_csv, 'receiver_csv', receiver_csv, ...
    'mat', mat_file, 'report', report_file);
save(mat_file, 'output');
local_write_report_ascii(report_file, output);

disp(convergence_table);
disp(receiver_table);
disp(checks);
if ~passed
    failed = strjoin(cellstr(checks.check_name(~checks.passed)), ', ');
    error('generate_bellhop_flat_surface_visuals_vertical:Failed', ...
        'Bellhop visualization checks failed: %s', failed);
end
fprintf('Bellhop visualization passed: %s\n', report_file);
end

function cfg = local_defaults(project_root, source_file, validation)
v = validation.config;
cfg = struct();
cfg.source_validation_mat = source_file;
cfg.output_dir = fullfile(project_root, 'results', 'visualization', ...
    'bellhop_flat_surface');
cfg.bellhop_exe = '';
cfg.water_depth_m = v.water_depth_m;
cfg.c0_mps = v.c0_mps;
cfg.z_tx_m = v.z_tx_m;
cfg.z_rx_m = v.z_rx_m;
cfg.receiver_offsets_m = v.receiver_offsets_m;
cfg.representative_offset_m = v.convergence_offset_m;
cfg.frequency_hz = v.f_ref_hz;
cfg.ray_beam_count = 51;
cfg.ray_step_m = 0.01;
cfg.field_beam_counts = [5001, 10001];
cfg.angle_limits_deg = [-89.5, -60];
cfg.range_limits_m = [0.25, 12];
cfg.range_count = 241;
cfg.depth_limits_m = [0.5, 99.5];
cfg.depth_count = 199;
cfg.z_box_m = 99.9;
cfg.source_mask_radius_m = 1;
cfg.tl_limits_db = [20, 80];
cfg.valid_tl_max_db = 100;
cfg.coherent_rms_tolerance_db = 1;
cfg.coherent_p95_tolerance_db = 3;
cfg.incoherent_rms_tolerance_db = 0.5;
cfg.incoherent_p95_tolerance_db = 1.5;
cfg.receiver_tl_tolerance_db = 0.5;
cfg.invariant_tolerance = 1e-10;
cfg.show_figures = false;
end

function cfg = local_overrides(cfg, overrides)
if ~isstruct(overrides), error('overrides must be a struct.'); end
names = fieldnames(overrides);
for ii = 1:numel(names)
    if ~isfield(cfg, names{ii}), error('Unknown override: %s', names{ii}); end
    cfg.(names{ii}) = overrides.(names{ii});
end
end

function local_validate_cfg(cfg, validation)
if cfg.water_depth_m <= 0 || cfg.c0_mps <= 0
    error('water_depth_m and c0_mps must be positive.');
end
if cfg.z_rx_m < 0 || cfg.z_rx_m >= cfg.z_tx_m || ...
        cfg.z_tx_m > cfg.water_depth_m
    error('Require 0 <= z_rx_m < z_tx_m <= water_depth_m.');
end
if cfg.z_box_m <= cfg.z_tx_m || cfg.z_box_m >= cfg.water_depth_m
    error('z_box_m must be below the source and above the seabed.');
end
if numel(cfg.field_beam_counts) ~= 2 || ...
        any(diff(cfg.field_beam_counts) <= 0)
    error('field_beam_counts must contain two increasing values.');
end
if cfg.ray_beam_count < 2 || any(cfg.angle_limits_deg < -90) || ...
        any(cfg.angle_limits_deg > 90) || diff(cfg.angle_limits_deg) <= 0
    error('Invalid ray beam count or angle limits.');
end
if cfg.ray_step_m <= 0
    error('ray_step_m must be positive so ZBOX terminates rays before the seabed.');
end
if cfg.range_count < 2 || cfg.depth_count < 2 || ...
        cfg.range_count * cfg.depth_count > 50000
    error('Field grid must contain between 2 and 50000 points.');
end
if cfg.range_limits_m(1) <= 0 || diff(cfg.range_limits_m) <= 0 || ...
        cfg.depth_limits_m(1) < 0 || diff(cfg.depth_limits_m) <= 0
    error('Invalid field range/depth limits.');
end
if cfg.depth_limits_m(2) >= cfg.z_box_m
    error('Maximum receiver depth must be shallower than z_box_m.');
end
if ~isequal(cfg.receiver_offsets_m(:).', [3, 6, 9]) || ...
        cfg.representative_offset_m ~= 6
    error('This display is defined for receiver offsets 3/6/9 m and x=6 m overview.');
end
if abs(validation.config.f_ref_hz - cfg.frequency_hz) > 1e-9
    error('Display frequency must match the saved PE/Bellhop reference frequency.');
end
saved_environment = [validation.config.water_depth_m, ...
    validation.config.c0_mps, validation.config.z_tx_m, ...
    validation.config.z_rx_m, validation.config.receiver_offsets_m(:).'];
display_environment = [cfg.water_depth_m, cfg.c0_mps, cfg.z_tx_m, ...
    cfg.z_rx_m, cfg.receiver_offsets_m(:).'];
if max(abs(saved_environment - display_environment)) > 1e-9
    error('Display environment must match the saved PE/Bellhop matrix.');
end
end

function exe = local_find_bellhop(configured, validation)
candidates = {};
if ~isempty(configured), candidates{end+1} = char(configured); end
env_exe = getenv('BELLHOP_EXE');
if ~isempty(env_exe), candidates{end+1} = env_exe; end
if isfield(validation, 'bellhop_executable') && ...
        ~isempty(validation.bellhop_executable)
    candidates{end+1} = char(validation.bellhop_executable);
end
matlab_exe = which('bellhop.exe');
if ~isempty(matlab_exe), candidates{end+1} = matlab_exe; end
for ii = 1:numel(candidates)
    if exist(candidates{ii}, 'file') == 2
        exe = candidates{ii};
        return
    end
end
if ispc, [status, found] = system('where bellhop.exe');
else, [status, found] = system('which bellhop');
end
if status == 0
    lines = regexp(strtrim(found), '\r?\n', 'split');
    if exist(strtrim(lines{1}), 'file') == 2
        exe = strtrim(lines{1});
        return
    end
end
error('Bellhop executable not found; set BELLHOP_EXE or bellhop_exe.');
end

function label = local_mode_label(mode)
if mode == "C", label = 'coherent';
elseif mode == "I", label = 'incoherent';
else, error('Unsupported field mode: %s', mode);
end
end

function local_write_env(file, cfg, run_type, beam_count, field_grid)
fid = fopen(file, 'w');
if fid < 0, error('Cannot create Bellhop environment: %s', file); end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '''Bellhop flat-surface R-C-I visualization''\n');
fprintf(fid, '%.12g\n1\n''CVW''\n', cfg.frequency_hz);
fprintf(fid, '2 0 %.12g\n0 %.12g /\n%.12g %.12g /\n', ...
    cfg.water_depth_m, cfg.c0_mps, cfg.water_depth_m, cfg.c0_mps);
fprintf(fid, '''A'' 0\n%.12g 1800 0 2 0 /\n', cfg.water_depth_m);
fprintf(fid, '1\n%.12g /\n', cfg.z_tx_m);
if field_grid
    fprintf(fid, '%d\n%.12g %.12g /\n', cfg.depth_count, cfg.depth_limits_m);
    fprintf(fid, '%d\n%.12g %.12g /\n', cfg.range_count, ...
        cfg.range_limits_m / 1000);
else
    fprintf(fid, '1\n%.12g /\n', cfg.z_rx_m);
    fprintf(fid, '1\n%.12g /\n', cfg.representative_offset_m / 1000);
end
fprintf(fid, '''%s''\n%d\n%.12g %.12g /\n', run_type, beam_count, ...
    cfg.angle_limits_deg);
if field_grid, step_m = 0; else, step_m = cfg.ray_step_m; end
fprintf(fid, '%.12g %.12g %.12g\n', step_m, cfg.z_box_m, ...
    cfg.range_limits_m(2) / 1000);
clear cleanup
end

function local_run_bellhop(exe, root, expected_extension)
run_dir = fileparts(root);
[~, name] = fileparts(root);
old_dir = pwd;
cleanup = onCleanup(@() cd(old_dir));
cd(run_dir);
generated = {'.ray', '.shd', '.prt'};
for ii = 1:numel(generated)
    target = [name generated{ii}];
    if exist(target, 'file') == 2, delete(target); end
end
[status, command_output] = system(sprintf('"%s" "%s"', exe, name));
if status ~= 0
    error('Bellhop failed for %s:\n%s', name, command_output);
end
expected = [root expected_extension];
if exist(expected, 'file') ~= 2 || exist([root '.prt'], 'file') ~= 2
    error('Bellhop did not create expected output for %s.', name);
end
clear cleanup
end

function data = local_read_ray(file)
fid = fopen(file, 'r');
if fid < 0, error('Cannot open ray file: %s', file); end
cleanup = onCleanup(@() fclose(fid));
title_text = strtrim(fgetl(fid));
frequency_hz = fscanf(fid, '%f', 1);
source_counts = fscanf(fid, '%d', 3);
beam_counts = fscanf(fid, '%d', 2);
top_depth_m = fscanf(fid, '%f', 1);
bottom_depth_m = fscanf(fid, '%f', 1);
fgetl(fid);
ray_type = strtrim(fgetl(fid));
expected = prod(source_counts) * prod(beam_counts);
template = struct('launch_angle_deg', NaN, 'n_steps', 0, ...
    'top_bounce_count', 0, 'bottom_bounce_count', 0, ...
    'range_m', [], 'depth_m', []);
rays = repmat(template, expected, 1);
count = 0;
while count < expected
    angle = fscanf(fid, '%f', 1);
    if isempty(angle), break; end
    n_steps = fscanf(fid, '%d', 1);
    top_count = fscanf(fid, '%d', 1);
    bottom_count = fscanf(fid, '%d', 1);
    coordinates = fscanf(fid, '%f', [2, n_steps]);
    if size(coordinates, 2) ~= n_steps
        error('Ray file ended inside ray %d.', count + 1);
    end
    count = count + 1;
    rays(count) = struct('launch_angle_deg', angle, 'n_steps', n_steps, ...
        'top_bounce_count', top_count, 'bottom_bounce_count', bottom_count, ...
        'range_m', coordinates(1, :), 'depth_m', coordinates(2, :));
end
rays = rays(1:count);
data = struct('title', title_text, 'frequency_hz', frequency_hz, ...
    'source_counts', source_counts(:).', 'beam_counts', beam_counts(:).', ...
    'top_depth_m', top_depth_m, 'bottom_depth_m', bottom_depth_m, ...
    'ray_type', ray_type, 'rays', rays, 'ray_count', count, ...
    'expected_ray_count', expected, 'file', file);
clear cleanup
end

function shade = local_read_shd_2d(file)
fid = fopen(file, 'rb');
if fid < 0, error('Cannot open shade file: %s', file); end
cleanup = onCleanup(@() fclose(fid));
record_words = fread(fid, 1, 'int32');
if isempty(record_words) || record_words <= 0
    error('Invalid shade-file record length.');
end
record_bytes = 4 * record_words;
title_text = deblank(fread(fid, 80, '*char').');
fseek(fid, record_bytes, 'bof');
plot_type = deblank(fread(fid, 10, '*char').');
fseek(fid, 2 * record_bytes, 'bof');
frequency_hz = fread(fid, 1, 'float32');
n_theta = fread(fid, 1, 'int32');
n_sx = fread(fid, 1, 'int32');
n_sy = fread(fid, 1, 'int32');
n_sd = fread(fid, 1, 'int32');
n_rd = fread(fid, 1, 'int32');
n_rr = fread(fid, 1, 'int32');
attenuation = fread(fid, 1, 'float32');
counts = [n_theta, n_sx, n_sy, n_sd, n_rd, n_rr];
if any(counts < 1) || n_theta ~= 1 || n_sd ~= 1
    error('Expected a single-bearing, single-source-depth 2-D shade file.');
end
fseek(fid, 3 * record_bytes, 'bof');
theta_deg = fread(fid, n_theta, 'float32');
fseek(fid, 4 * record_bytes, 'bof'); source_x_m = fread(fid, n_sx, 'float32');
fseek(fid, 5 * record_bytes, 'bof'); source_y_m = fread(fid, n_sy, 'float32');
fseek(fid, 6 * record_bytes, 'bof'); source_depth_m = fread(fid, n_sd, 'float32');
fseek(fid, 7 * record_bytes, 'bof'); receiver_depth_m = fread(fid, n_rd, 'float32');
fseek(fid, 8 * record_bytes, 'bof'); receiver_range_m = fread(fid, n_rr, 'float32');
pressure = complex(zeros(n_rd, n_rr));
for dd = 1:n_rd
    record_number = 9 + dd - 1;
    if fseek(fid, record_number * record_bytes, 'bof') ~= 0
        error('Failed seeking to shade pressure record %d.', record_number);
    end
    raw = fread(fid, 2 * n_rr, 'float32');
    if numel(raw) ~= 2 * n_rr
        error('Shade file ended inside receiver-depth record %d.', dd);
    end
    pressure(dd, :) = raw(1:2:end) + 1i * raw(2:2:end);
end
shade = struct('title', title_text, 'plot_type', plot_type, ...
    'frequency_hz', frequency_hz, 'attenuation', attenuation, ...
    'theta_deg', theta_deg, 'source_x_m', source_x_m, ...
    'source_y_m', source_y_m, 'source_depth_m', source_depth_m, ...
    'receiver_depth_m', receiver_depth_m, ...
    'receiver_range_m', receiver_range_m, 'pressure', pressure, ...
    'counts', counts, 'record_bytes', record_bytes, 'file', file);
clear cleanup
end

function [tl_db, source_mask] = local_tl_field(shade, cfg)
[range_grid_m, depth_grid_m] = meshgrid(shade.receiver_range_m, ...
    shade.receiver_depth_m);
source_mask = hypot(range_grid_m, depth_grid_m - cfg.z_tx_m) < ...
    cfg.source_mask_radius_m;
amplitude = abs(shade.pressure);
tl_db = -20 * log10(max(amplitude, 1e-37));
tl_db(source_mask) = NaN;
end

function t = local_convergence_table(c_low, c_high, i_low, i_high, mask, cfg)
names = ["coherent"; "incoherent"];
low_fields = {c_low; i_low};
high_fields = {c_high; i_high};
rms_values = zeros(2, 1);
p95_values = zeros(2, 1);
point_counts = zeros(2, 1);
for ii = 1:2
    valid = ~mask & isfinite(low_fields{ii}) & isfinite(high_fields{ii}) & ...
        high_fields{ii} <= cfg.valid_tl_max_db;
    delta = low_fields{ii}(valid) - high_fields{ii}(valid);
    if isempty(delta), error('No valid TL convergence points for %s.', names(ii)); end
    rms_values(ii) = sqrt(mean(delta.^2));
    p95_values(ii) = local_percentile(abs(delta), 0.95);
    point_counts(ii) = numel(delta);
end
t = table(names, repmat(cfg.field_beam_counts(1), 2, 1), ...
    repmat(cfg.field_beam_counts(2), 2, 1), rms_values, p95_values, ...
    point_counts, 'VariableNames', {'mode', 'low_beam_count', ...
    'high_beam_count', 'tl_rms_difference_db', ...
    'tl_p95_abs_difference_db', 'valid_point_count'});
end

function value = local_percentile(x, probability)
x = sort(x(:));
position = 1 + (numel(x) - 1) * probability;
lo = floor(position); hi = ceil(position);
if lo == hi, value = x(lo);
else, value = x(lo) + (position - lo) * (x(hi) - x(lo));
end
end

function t = local_receiver_table(coherent, incoherent, validation, cfg)
n = numel(cfg.receiver_offsets_m);
rows = repmat(struct('offset_m', NaN, 'field_range_m', NaN, ...
    'field_depth_m', NaN, 'coherent_field_tl_db', NaN, ...
    'incoherent_field_tl_db', NaN, 'arrival_synthesis_tl_db', NaN, ...
    'field_minus_arrival_tl_db', NaN, 'pe_scaled_total_tl_db', NaN), n, 1);
f = cfg.frequency_hz;
for ii = 1:n
    x = cfg.receiver_offsets_m(ii);
    [~, ir] = min(abs(coherent.receiver_range_m - x));
    [~, iz] = min(abs(coherent.receiver_depth_m - cfg.z_rx_m));
    pressure_c = coherent.pressure(iz, ir);
    pressure_i = incoherent.pressure(iz, ir);
    paths = validation.cases(ii).bellhop_paths;
    pressure_arrival = sum(paths.amplitude_complex .* ...
        exp(1i * 2 * pi * f * paths.delay_s));
    pe = validation.cases(ii).pe;
    [~, ifref] = min(abs(pe.f_axis(:) - f));
    pressure_pe = validation.global_direct_scale * ( ...
        pe.H_direct_physical_f(ifref) + pe.H_reflect_physical_f(ifref));
    field_tl = -20*log10(max(abs(pressure_c), realmin));
    arrival_tl = -20*log10(max(abs(pressure_arrival), realmin));
    rows(ii) = struct('offset_m', x, ...
        'field_range_m', coherent.receiver_range_m(ir), ...
        'field_depth_m', coherent.receiver_depth_m(iz), ...
        'coherent_field_tl_db', field_tl, ...
        'incoherent_field_tl_db', -20*log10(max(abs(pressure_i), realmin)), ...
        'arrival_synthesis_tl_db', arrival_tl, ...
        'field_minus_arrival_tl_db', field_tl - arrival_tl, ...
        'pe_scaled_total_tl_db', -20*log10(max(abs(pressure_pe), realmin)));
end
t = struct2table(rows);
end

function error_value = local_saved_invariant_error(validation)
error_value = 0;
for ii = 1:numel(validation.cases)
    pe = validation.cases(ii).pe;
    error_value = max(error_value, max(abs(pe.H_f(:) - ...
        pe.H_direct_f(:) - pe.H_reflect_f(:))));
end
end

function checks = local_checks(ray, coherent, incoherent, convergence, ...
    receiver, invariant_error, cfg)
rays = ray.rays;
top_counts = [rays.top_bounce_count];
bottom_counts = [rays.bottom_bounce_count];
metadata_error = max([abs(coherent.frequency_hz-cfg.frequency_hz), ...
    abs(incoherent.frequency_hz-cfg.frequency_hz), ...
    abs(coherent.source_depth_m(1)-cfg.z_tx_m), ...
    abs(incoherent.source_depth_m(1)-cfg.z_tx_m), ...
    abs(numel(coherent.receiver_range_m)-cfg.range_count), ...
    abs(numel(coherent.receiver_depth_m)-cfg.depth_count)]);
finite_error = double(any(~isfinite(coherent.pressure(:))) || ...
    any(~isfinite(incoherent.pressure(:))));
axis_error = double(any(diff(coherent.receiver_range_m) <= 0) || ...
    any(diff(coherent.receiver_depth_m) <= 0) || ...
    any(diff(incoherent.receiver_range_m) <= 0) || ...
    any(diff(incoherent.receiver_depth_m) <= 0) || ...
    max(abs(coherent.receiver_range_m([1 end]).' - cfg.range_limits_m)) > 1e-6 || ...
    max(abs(coherent.receiver_depth_m([1 end]).' - cfg.depth_limits_m)) > 1e-6);
color_limit_error = double(numel(cfg.tl_limits_db) ~= 2 || ...
    any(~isfinite(cfg.tl_limits_db)) || diff(cfg.tl_limits_db) <= 0);
c = convergence(convergence.mode == "coherent", :);
i = convergence(convergence.mode == "incoherent", :);
names = ["ray_count"; "ray_direct_present"; "ray_surface_present"; ...
    "ray_bottom_absent"; "shade_metadata"; "shade_finite"; ...
    "shade_axes_positive_down"; "shared_color_limits"; ...
    "coherent_tl_rms"; "coherent_tl_p95"; "incoherent_tl_rms"; ...
    "incoherent_tl_p95"; "receiver_field_arrival_tl"; "saved_pe_invariant"];
values = [ray.ray_count-ray.expected_ray_count; ...
    double(~any(top_counts==0 & bottom_counts==0)); ...
    double(~any(top_counts>=1 & bottom_counts==0)); ...
    max(bottom_counts); metadata_error; finite_error; axis_error; ...
    color_limit_error; ...
    c.tl_rms_difference_db; c.tl_p95_abs_difference_db; ...
    i.tl_rms_difference_db; i.tl_p95_abs_difference_db; ...
    max(abs(receiver.field_minus_arrival_tl_db)); invariant_error];
limits = [0;0;0;0;1e-6;0;0;0;cfg.coherent_rms_tolerance_db; ...
    cfg.coherent_p95_tolerance_db;cfg.incoherent_rms_tolerance_db; ...
    cfg.incoherent_p95_tolerance_db;cfg.receiver_tl_tolerance_db; ...
    cfg.invariant_tolerance];
relations = [repmat("==",4,1); "<="; repmat("==",3,1); repmat("<=",6,1)];
passed = (relations=="==" & values==limits) | ...
    (relations=="<=" & values<=limits);
checks = table(names, values, relations, limits, passed, ...
    'VariableNames', {'check_name','value','relation','limit','passed'});
end

function local_plot_rays(ray_data, cfg, file)
visibility = local_visibility(cfg.show_figures);
fig = figure('Visible', visibility, 'Color', 'w', 'Position', [100 100 980 760]);
cleanup = onCleanup(@() close(fig));
hold on;
rays = ray_data.rays;
for ii = 1:numel(rays)
    if rays(ii).bottom_bounce_count > 0, color = [0.75 0.2 0.2];
    elseif rays(ii).top_bounce_count > 0, color = [0.25 0.55 0.85];
    else, color = [0.65 0.65 0.65];
    end
    plot(rays(ii).range_m, rays(ii).depth_m, '-', 'Color', color, ...
        'LineWidth', 0.65, 'HandleVisibility', 'off');
end
display_range_m = [0, cfg.range_limits_m(2)];
plot(display_range_m, [0 0], 'k-', 'LineWidth', 1.5, ...
    'DisplayName', 'pressure-release surface');
plot(display_range_m, cfg.water_depth_m*[1 1], 'k--', 'LineWidth', 1.2, ...
    'DisplayName', 'seabed (excluded by ZBOX)');
plot(0, cfg.z_tx_m, 'p', 'MarkerSize', 13, 'MarkerFaceColor', [0.9 0.25 0.1], ...
    'MarkerEdgeColor', 'k', 'DisplayName', 'Tx');
plot(cfg.receiver_offsets_m, cfg.z_rx_m*ones(size(cfg.receiver_offsets_m)), ...
    'ks', 'MarkerFaceColor', [1 0.85 0.1], 'MarkerSize', 8, ...
    'DisplayName', 'Rx: x=3/6/9 m');
x = cfg.representative_offset_m;
plot([0 x], [cfg.z_tx_m cfg.z_rx_m], 'r-', 'LineWidth', 2.6, ...
    'DisplayName', 'analytic direct, x=6 m');
bounce_x = x * cfg.z_tx_m / (cfg.z_tx_m + cfg.z_rx_m);
plot([0 bounce_x x], [cfg.z_tx_m 0 cfg.z_rx_m], 'm-', 'LineWidth', 2.6, ...
    'DisplayName', 'analytic surface path, x=6 m');
set(gca, 'YDir', 'reverse');
xlim(display_range_m); ylim([0 cfg.water_depth_m]);
grid on; xlabel('Horizontal range (m)'); ylabel('Depth (m)');
title(sprintf('Bellhop ray fan, %.0f Hz, %d rays', ...
    cfg.frequency_hz, cfg.ray_beam_count));
legend('Location', 'eastoutside');
exportgraphics(fig, file, 'Resolution', 190);
clear cleanup
end

function local_plot_fields(shade, coherent_tl, incoherent_tl, cfg, file)
visibility = local_visibility(cfg.show_figures);
fig = figure('Visible', visibility, 'Color', 'w', 'Position', [100 100 1180 720]);
cleanup = onCleanup(@() close(fig));
fields = {coherent_tl, incoherent_tl};
titles = {'Coherent TL: interference field', 'Incoherent TL: energy envelope'};
for ii = 1:2
    subplot(1,2,ii);
    h = imagesc(shade.receiver_range_m, shade.receiver_depth_m, fields{ii});
    set(h, 'AlphaData', isfinite(fields{ii}));
    set(gca, 'YDir', 'reverse');
    clim(cfg.tl_limits_db); axis tight;
    xlim([0, cfg.range_limits_m(2)]);
    hold on;
    plot(0, cfg.z_tx_m, 'wp', 'MarkerFaceColor', [0.9 0.2 0.1], ...
        'MarkerEdgeColor', 'k', 'MarkerSize', 11);
    plot(cfg.receiver_offsets_m, cfg.z_rx_m*ones(size(cfg.receiver_offsets_m)), ...
        'ws', 'MarkerFaceColor', [1 0.85 0.1], 'MarkerEdgeColor', 'k', ...
        'MarkerSize', 6);
    xlabel('Horizontal range (m)'); ylabel('Depth (m)'); title(titles{ii});
    colorbar;
end
colormap(fig, flipud(turbo(256)));
sgtitle(sprintf('Bellhop flat-surface TL at %.0f Hz, %d beams; color scale %.0f--%.0f dB', ...
    cfg.frequency_hz, cfg.field_beam_counts(end), cfg.tl_limits_db));
exportgraphics(fig, file, 'Resolution', 190);
clear cleanup
end

function local_plot_slice(shade, coherent_tl, incoherent_tl, receiver, cfg, file)
visibility = local_visibility(cfg.show_figures);
fig = figure('Visible', visibility, 'Color', 'w', 'Position', [100 100 980 660]);
cleanup = onCleanup(@() close(fig));
[~, iz] = min(abs(shade.receiver_depth_m - cfg.z_rx_m));
plot(shade.receiver_range_m, coherent_tl(iz,:), 'b-', 'LineWidth', 1.5, ...
    'DisplayName', 'Bellhop coherent TL');
hold on;
plot(shade.receiver_range_m, incoherent_tl(iz,:), 'Color', [0.85 0.35 0.1], ...
    'LineWidth', 1.5, 'DisplayName', 'Bellhop incoherent TL');
plot(receiver.offset_m, receiver.arrival_synthesis_tl_db, 'rd', ...
    'MarkerFaceColor', 'r', 'MarkerSize', 7, ...
    'DisplayName', 'Bellhop selected-arrival synthesis');
plot(receiver.offset_m, receiver.pe_scaled_total_tl_db, 'ko', ...
    'MarkerFaceColor', [1 0.85 0.1], 'MarkerSize', 8, ...
    'DisplayName', 'PE total H, one global scale');
grid on; xlim(cfg.range_limits_m); ylim(cfg.tl_limits_db);
set(gca, 'YDir', 'reverse');
xlabel('Horizontal range (m)'); ylabel('Transmission loss (dB)');
title(sprintf('Receiver-depth TL slice at z = %.2f m', ...
    shade.receiver_depth_m(iz)));
legend('Location', 'best');
exportgraphics(fig, file, 'Resolution', 190);
clear cleanup
end

function value = local_visibility(show_figures)
if show_figures, value = 'on'; else, value = 'off'; end
end

function local_write_report_ascii(file, output)
cfg = output.config;
conv = output.convergence_table;
rx = output.receiver_table;
checks = output.checks;
fid = fopen(file, 'w', 'n', 'UTF-8');
if fid < 0, error('Cannot write report: %s', file); end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '# Bellhop Flat-Surface Ray and TL Visualization\n\n');
fprintf(fid, ['Visualization checks passed: `%s`. The source PE/Bellhop ', ...
    'strict matrix remains `%s`; this display does not change that ', ...
    'validation result and does not modify the PE or communication chain.\n\n'], ...
    string(output.passed), string(output.source_matrix_passed));
if ~output.passed
    failed = strjoin(cellstr(checks.check_name(~checks.passed)), ', ');
    fprintf(fid, ['Failed checks: %s. Figures and raw outputs were retained ', ...
        'for diagnosis.\n\n'], failed);
end
fprintf(fid, '## Environment and Bellhop modes\n\n');
fprintf(fid, ['- Water depth %.1f m, uniform sound speed %.1f m/s, ', ...
    'Tx depth %.1f m, frequency %.1f Hz.\n'], ...
    cfg.water_depth_m, cfg.c0_mps, cfg.z_tx_m, cfg.frequency_hz);
fprintf(fid, ['- `R`: %d central rays; `C`: coherent TL; ', ...
    '`I`: incoherent TL.\n'], cfg.ray_beam_count);
fprintf(fid, ['- TL grid %d x %d, range %.2f--%.2f m, ', ...
    'depth %.2f--%.2f m.\n'], cfg.depth_count, cfg.range_count, ...
    cfg.range_limits_m, cfg.depth_limits_m);
fprintf(fid, ['- `ZBOX=%.1f m < %.1f m`; the explicit %.3f m ray step ', ...
    'terminates displayed rays before the seabed.\n'], ...
    cfg.z_box_m, cfg.water_depth_m, cfg.ray_step_m);
fprintf(fid, ['- Both field plots use %.0f--%.0f dB; the %.1f m ', ...
    'source neighborhood is masked.\n\n'], ...
    cfg.tl_limits_db, cfg.source_mask_radius_m);
fprintf(fid, '## Beam-count convergence\n\n');
fprintf(fid, ['| Mode | Beams low/high | TL RMS difference (dB) | ', ...
    '|difference| 95th percentile (dB) | Valid points |\n']);
fprintf(fid, '|---|---:|---:|---:|---:|\n');
for ii = 1:height(conv)
    fprintf(fid, '| %s | %d/%d | %.6f | %.6f | %d |\n', conv.mode(ii), ...
        conv.low_beam_count(ii), conv.high_beam_count(ii), ...
        conv.tl_rms_difference_db(ii), conv.tl_p95_abs_difference_db(ii), ...
        conv.valid_point_count(ii));
end
fprintf(fid, '\n## Receiver-point TL\n\n');
fprintf(fid, ['| x (m) | C field TL | I field TL | Arrival synthesis TL | ', ...
    'C field - arrival | Globally scaled PE total TL |\n']);
fprintf(fid, '|---:|---:|---:|---:|---:|---:|\n');
for ii = 1:height(rx)
    fprintf(fid, '| %.1f | %.5f | %.5f | %.5f | %+.5f | %.5f |\n', ...
        rx.offset_m(ii), rx.coherent_field_tl_db(ii), ...
        rx.incoherent_field_tl_db(ii), rx.arrival_synthesis_tl_db(ii), ...
        rx.field_minus_arrival_tl_db(ii), rx.pe_scaled_total_tl_db(ii));
end
fprintf(fid, '\n## Automated checks\n\n');
fprintf(fid, '| Check | Value | Relation | Limit | Passed |\n');
fprintf(fid, '|---|---:|:---:|---:|:---:|\n');
for ii = 1:height(checks)
    fprintf(fid, '| %s | %.8g | %s | %.8g | %d |\n', checks.check_name(ii), ...
        checks.value(ii), checks.relation(ii), checks.limit(ii), checks.passed(ii));
end
fprintf(fid, '\n## Figures\n\n');
fprintf(fid, '![Bellhop ray geometry](bellhop_ray_geometry.png)\n\n');
fprintf(fid, '![Bellhop coherent and incoherent TL](bellhop_tl_fields.png)\n\n');
fprintf(fid, '![Receiver-depth TL slice](bellhop_receiver_depth_tl.png)\n\n');
fprintf(fid, '## Limitations\n\n');
fprintf(fid, ['Only the three saved PE receiver responses are overlaid on ', ...
    'the TL slice; no nonexistent PE 2-D field is synthesized. Bellhop TL ', ...
    'is a single-frequency 4 kHz result, not a 3--5 kHz wideband PDP. ', ...
    'Coherent TL contains phase interference, while incoherent TL is an ', ...
    'energy envelope; they are not interchangeable.\n']);
clear cleanup
end
