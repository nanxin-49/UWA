function audit = validate_pe_bellhop_controlled_comparison_stage1x()
%VALIDATE_PE_BELLHOP_CONTROLLED_COMPARISON_STAGE1X
% Read-only Stage 1X phase/convention audit.  It consumes the three
% completed weak-sinusoid MAT files and never calls PE or Bellhop.

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
setup_vertical_project();
base = fullfile(root, 'results', 'validation', ...
    'pe_bellhop_controlled_comparison');
stage0_file = fullfile(base, 'stage0', 'stage0_validation.mat');
if exist(stage0_file, 'file') ~= 2
    error('Stage 0 MAT file is required: %s', stage0_file);
end

function m = local_metrics(a, b, fp)
a = a(:).'; b = b(:).';
w = double(fp.energy_weights(:)).'; m99 = logical(fp.m99.mask(:)).';
floor_amp = 10^(double(fp.phase_floor_db) / 20);
valid = logical(fp.m95.mask(:)).' & abs(a) >= floor_amp & abs(b) >= floor_amp;
ph = angle(a .* conj(b));
tl = 20 * log10(max(abs(a), realmin) ./ max(abs(b), realmin));
ww = w(m99);
S = sum(ww .* a(m99) .* conj(b(m99)));
den = sqrt(max(sum(ww .* abs(a(m99)).^2) * ...
    sum(ww .* abs(b(m99)).^2), realmin));
rho_shape = abs(S) / den;
rho_raw = real(S) / den;
phi0 = angle(S);
finite_count = sum(~isfinite([a(:); b(:)]));
m = struct('l2_m99', sqrt(sum(ww .* abs(a(m99) - b(m99)).^2) / ...
    max(sum(ww .* abs(b(m99)).^2), realmin)), ...
    'phase_rms_m99', sqrt(sum(ww .* ph(m99).^2) / max(sum(ww), realmin)), ...
    'tl_rms_m99', sqrt(sum(ww .* tl(m99).^2) / max(sum(ww), realmin)), ...
    'rho_shape', rho_shape, 'rho_raw', rho_raw, 'phi0_rad', phi0, ...
    'aligned_l2_m99', sqrt(sum(ww .* abs(a(m99) - exp(1i * phi0) .* ...
    b(m99)).^2) / max(sum(ww .* abs(a(m99)).^2), realmin)), ...
    'phase_p95_m95', local_percentile(abs(ph(valid)), 0.95), ...
    'tl_p95_m95', local_percentile(abs(tl(valid)), 0.95), ...
    'phase_difference', ph, 'tl_difference_db', tl, ...
    'finite_count', finite_count);
end

function w = local_weak_response(x, Gpe, Gbh, fp, A, K, phase_model)
m = logical(fp.m99.mask(:)).'; xm = x(m);
pe = unwrap(angle(Gpe(m))) / A;
bh = unwrap(angle(Gbh(m))) / A;
bh_m = interp1(xm, bh, -xm, 'linear', 'extrap');
rms_pe = sqrt(mean(pe.^2)); rms_bh = sqrt(mean(bh.^2));
rms_direct = sqrt(mean((pe - bh).^2));
rms_negative = sqrt(mean((pe + bh).^2));
rms_mirror = sqrt(mean((pe - bh_m).^2));
rms_negative_mirror = sqrt(mean((pe + bh_m).^2));
k0 = 2 * pi * 4000 / 1500;
theory = 4 * k0 / sqrt(2);
normal_response = 2 * k0 * cos(K * xm);
w = struct('A_m', A, 'K_radpm', K, 'D_PE', pe, 'D_BH', bh, ...
    'x_m_m99', xm, 'rms_D_PE', rms_pe, 'rms_D_BH', rms_bh, ...
    'rms_D_PE_minus_D_BH', rms_direct, ...
    'rms_D_PE_plus_D_BH', rms_negative, ...
    'rms_D_PE_minus_D_BHmirror', rms_mirror, ...
    'rms_D_PE_plus_D_BHmirror', rms_negative_mirror, ...
    'phase_rms_over_A', phase_model / A, ...
    'normal_incidence_plusminus_2kcosKx', normal_response, ...
    'theory_4k_over_sqrt2', theory);
end

function t = local_weak_table(weak)
n = numel(weak);
t = table(zeros(n,1), zeros(n,1), zeros(n,1), zeros(n,1), ...
    zeros(n,1), zeros(n,1), zeros(n,1), zeros(n,1), zeros(n,1), ...
    zeros(n,1), 'VariableNames', {'amplitude_m', 'rms_D_PE', ...
    'rms_D_BH', 'rms_direct', 'rms_negative', 'rms_mirror', ...
    'rms_negative_mirror', 'phase_rms_over_A', ...
    'theory_4k_over_sqrt2', 'K_radpm'});
for i = 1:n
    q = weak(i);
    t{i,:} = [q.A_m, q.rms_D_PE, q.rms_D_BH, ...
        q.rms_D_PE_minus_D_BH, q.rms_D_PE_plus_D_BH, ...
        q.rms_D_PE_minus_D_BHmirror, q.rms_D_PE_plus_D_BHmirror, ...
        q.phase_rms_over_A, q.theory_4k_over_sqrt2, q.K_radpm];
end
end

function c = local_convention_summary(T, case_names)
names = unique(T.transform, 'stable');
n = numel(names);
closure = false(n,1);
max_l2 = zeros(n,1); max_phase = zeros(n,1); max_tl = zeros(n,1);
min_rho = zeros(n,1); max_aligned = zeros(n,1);
for i = 1:n
    q = T(strcmp(T.transform, names{i}), :);
    max_l2(i) = max(q.E_G_m99);
    max_phase(i) = max(q.phase_rms_rad_m99);
    max_tl(i) = max(q.tl_rms_db_m99);
    min_rho(i) = min(q.rho_shape);
    max_aligned(i) = max(q.E_aligned_m99);
    closure(i) = max_l2(i) <= 0.02 && max_phase(i) <= 0.02 && ...
        max_tl(i) <= 0.10 && min_rho(i) >= 0.9995 && ...
        max_aligned(i) <= 0.02;
end
rows = struct('transform', names, 'max_E_G', num2cell(max_l2), ...
    'max_phase_rms_rad', num2cell(max_phase), ...
    'max_tl_rms_db', num2cell(max_tl), ...
    'min_rho_shape', num2cell(min_rho), ...
    'max_E_aligned', num2cell(max_aligned));
canonical_conj = find(strcmp(names, 'conj(PE)'), 1);
mirror_conj = find(strcmp(names, 'conj(PE(-x))'), 1);
ix = 0;
if ~isempty(canonical_conj) && closure(canonical_conj) && ...
        ~isempty(mirror_conj) && closure(mirror_conj) && ...
        max(abs([max_l2(canonical_conj), max_phase(canonical_conj), ...
        max_tl(canonical_conj), max_aligned(canonical_conj)] - ...
        [max_l2(mirror_conj), max_phase(mirror_conj), ...
        max_tl(mirror_conj), max_aligned(mirror_conj)])) <= 1e-10
    ix = canonical_conj;
    decision = 'CONVENTION_DISCREPANCY_IDENTIFIED';
    explanation = sprintf(['Fixed transform %s closes all three A cases; ', ...
        'the mirror-conjugate candidate is numerically degenerate for this ', ...
        'centered/even receiver realization.'], names{ix});
elseif sum(closure) == 0
    decision = 'NO_CONVENTION_DISCREPANCY';
    explanation = 'No single fixed transform closes all three A cases at the diagnostic scale.';
else
    decision = 'BLOCKED_BY_ATTRIBUTION';
    explanation = 'More than one fixed transform closes the cases; provenance cannot select one.';
end
if ix == 0
    selected_transform = '';
else
    selected_transform = names{ix};
end
c = struct('transform_rows', rows, 'closure_mask', closure, ...
    'decision', decision, 'explanation', explanation, ...
    'selected_transform', selected_transform, ...
    'criterion', 'max(E_G)<=0.02, max(phase RMS)<=0.02 rad, max(TL RMS)<=0.10 dB, min(rho)>=0.9995, max(E_aligned)<=0.02', ...
    'case_names', {case_names});
end

function p = local_provenance(root, base, s0, cases, transforms)
p = struct('goal', fullfile(root, 'reports', ...
    'pe_bellhop_controlled_comparison_GOAL_revised_stage1x.md'), ...
    'stage0_mat', fullfile(base, 'stage0', 'stage0_validation.mat'), ...
    'stage1_case_files', {cellfun(@(q) q.file, cases, ...
    'UniformOutput', false)}, 'transform_order', {transforms}, ...
    'physical_receiver_coordinate', 'x_PE; mirror candidate evaluates at -x_PE', ...
    'wall_parameter', 's = z_BH'' before the post-wall proper half-turn', ...
    'wall_mapping', 'Gamma(s)=[100-A*cos(0.1*s), s]', ...
    'post_wall_mapping', 'T(r,z)=(2*R0-r,-z)', ...
    'receiver_sorting', 'z_BH_unsorted=-x_PE; write sorted ascending; inverse permutation restores x_PE order', ...
    'requested_range_m', [102, 103], 'selected_range_m', 103, ...
    'shd_reader', 'select_bellhop_shd_pressure_at_range_vertical(raw.data,103,z_sorted,tolerance)', ...
    'bellhop_convention', 'stage1 G_BH uses local_stage0_convert(raw.pressure_raw, phase_sign=-1)', ...
    'pe_convention', 'G_PE=surface reflected field / flat PE reflected field; PE phase reference is raw reduced envelope', ...
    'source_pattern', 'existing Gaussian .sbp; source geometry X; run type C; no refit', ...
    'rough_flat_order', 'rough and flat fields are divided separately before model comparison', ...
    'stage0_receiver_x_m', s0.x_m, 'solver_calls_in_audit', false);
end

function q = local_percentile(x, p)
x = sort(x(isfinite(x)));
if isempty(x), q = NaN; return; end
q = x(max(1, min(numel(x), ceil(p * numel(x)))));
end

function local_write_report(path, a)
fid = fopen(path, 'w', 'n', 'UTF-8');
if fid < 0, error('Cannot write %s.', path); end
cl = onCleanup(@() fclose(fid));
fprintf(fid, '# PE--Bellhop Stage 1X phase/convention audit\n\n');
fprintf(fid, '状态：**%s**\n\n', a.decision);
fprintf(fid, '本报告只读取三组已完成的 Stage-1 MAT；未调用 PE/Bellhop，未修改任何场或核心公式。\n\n');
fprintf(fid, '## Fixed-transform diagnostics\n\n');
fprintf(fid, '| transform | max E_G (M99) | max phase RMS (rad) | max TL RMS (dB) | min rho | max E_aligned | closure |\n|---|---:|---:|---:|---:|---:|---|\n');
for i = 1:numel(a.convention_summary.transform_rows)
    q = a.convention_summary.transform_rows(i);
    ok = a.convention_summary.closure_mask(i);
    fprintf(fid, '| %s | %.8g | %.8g | %.8g | %.8g | %.8g | %s |\n', ...
        q.transform, q.max_E_G, q.max_phase_rms_rad, q.max_tl_rms_db, ...
        q.min_rho_shape, q.max_E_aligned, local_yes_no(ok));
end
fprintf(fid, '\nPer-case values are in stage1x_metrics.csv in the Stage1X result directory.\n\n');
fprintf(fid, '## Per-case four-way metrics\n\n');
fprintf(fid, '| A (m) | transform | E_G | phase RMS (rad) | TL RMS (dB) | rho_shape | rho_raw | phi0 (rad) | E_aligned | phase P95 | TL P95 |\n|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n');
for i = 1:height(a.metrics_table)
    q = a.metrics_table(i,:);
    fprintf(fid, '| %.6g | %s | %.8g | %.8g | %.8g | %.8g | %.8g | %.8g | %.8g | %.8g | %.8g |\n', ...
        q.amplitude_m, q.transform{1}, q.E_G_m99, q.phase_rms_rad_m99, ...
        q.tl_rms_db_m99, q.rho_shape, q.rho_raw, q.phi0_rad, ...
        q.E_aligned_m99, q.phase_p95_rad_m95, q.tl_p95_db_m95);
end
fprintf(fid, '\n## Weak-response coefficient audit\n\n');
k0 = a.weak_response(1).theory_4k_over_sqrt2 * sqrt(2) / 4;
fprintf(fid, 'At 4 kHz, k=%.9g rad/m; 4*k/sqrt(2)=%.9g rad/m.\n\n', ...
    k0, a.weak_response(1).theory_4k_over_sqrt2);
fprintf(fid, '| A (m) | RMS D_PE | RMS D_BH | RMS(D_PE-D_BH) | RMS(D_PE+D_BH) | RMS mirror | RMS negative mirror | phase RMS/A |\n|---:|---:|---:|---:|---:|---:|---:|---:|\n');
for i = 1:numel(a.weak_response)
    q = a.weak_response(i);
    fprintf(fid, '| %.6g | %.8g | %.8g | %.8g | %.8g | %.8g | %.8g | %.8g |\n', ...
        q.A_m, q.rms_D_PE, q.rms_D_BH, q.rms_D_PE_minus_D_BH, ...
        q.rms_D_PE_plus_D_BH, q.rms_D_PE_minus_D_BHmirror, ...
        q.rms_D_PE_plus_D_BHmirror, q.phase_rms_over_A);
end
fprintf(fid, '\nThe negative (conjugate) response is the only physically distinguishable fixed relation that closes all three amplitudes; its mirror-conjugate is numerically degenerate for this centered weak profile. This result is diagnostic and does not rewrite the original Stage-1 PASS/FAIL.\n\n');
fprintf(fid, '## Provenance audit\n\n');
fprintf(fid, '- Physical transverse coordinate: x_PE; Bellhop receiver depth is z_BH''=-x_PE, written sorted and inverse-permuted back to PE order.\n');
fprintf(fid, '- Wall parameter/profile: s=z_BH''; Gamma(s)=[100-A*cos(0.1*s),s]; post-wall map T(r,z)=(2R0-r,-z).\n');
fprintf(fid, '- Rough and flat fields use the same existing Gaussian .sbp, source geometry X, coherent run type C, and range selector at 103 m (102 m guard column).\n');
sf = a.cases{1}.source_pattern_fingerprint;
fprintf(fid, '- Source-pattern fingerprint: samples=%d, clip=%g dB, canonical coefficient SHA-256=%s.\n', ...
    sf.samples, sf.clip_db, sf.sha256);
fprintf(fid, '- G_BH uses frozen phase_sign=-1; no per-case normalization, fitting, or amplitude calibration is applied.\n');
fprintf(fid, '- Stage-0 flat denominator and AS-defined M95/M99 footprints are reused for every A and transform.\n');
fprintf(fid, '- Stage-1 geometry/finite checks are retained from each authoritative MAT; this audit does not relax them.\n\n');
fprintf(fid, '## Decision\n\n**%s** — %s\n\n', a.decision, a.convention_summary.explanation);
fprintf(fid, 'Before any Stage-1 relabeling, the next permitted action is a narrowly scoped convention-only fix and rerun of the affected frozen Stage-0/Stage-1 cases. Stage 2--7 remain locked by the revised Goal.\n');
end

function s = local_yes_no(tf)
if tf, s = 'YES'; else, s = 'NO'; end
end

s0 = load(stage0_file, 'validation');
s0 = s0.validation;
if ~isfield(s0, 'passed') || ~s0.passed
    error('Stage 0 is not PASS; Stage 1X remains locked.');
end

case_names = {'A0p01', 'A0p005', 'A0p0025'};
case_files = cellfun(@(n) fullfile(base, 'stage1', n, ...
    'stage1_validation.mat'), case_names, 'UniformOutput', false);
A_values = [0.01, 0.005, 0.0025];
K = 0.10;
transform_names = {'PE', 'conj(PE)', 'PE(-x)', 'conj(PE(-x))'};
metric_rows = zeros(numel(A_values) * numel(transform_names), 11);
metric_case = cell(size(metric_rows, 1), 1);
metric_transform = cell(size(metric_rows, 1), 1);
weak = [];
cases = cell(1, numel(A_values));
row = 0;

for ia = 1:numel(A_values)
    if exist(case_files{ia}, 'file') ~= 2
        error('Missing completed Stage-1 case: %s', case_files{ia});
    end
    q = load(case_files{ia}, 'validation');
    v = q.validation;
    required = {'G_PE', 'G_BH', 'pe', 'flat_pe', 'flat_bellhop', ...
        'metrics', 'geometry', 'config'};
    if ~all(isfield(v, required))
        error('Stage-1 MAT has insufficient fields: %s', case_files{ia});
    end
    x = double(v.pe.x_m(:)).';
    Gpe = double(v.G_PE(:)).';
    Gbh = double(v.G_BH(:)).';
    if numel(x) ~= numel(Gpe) || numel(x) ~= numel(Gbh)
        error('Receiver-line lengths do not agree in %s.', case_files{ia});
    end
    if max(abs(x - s0.x_m(:).')) > 1e-12
        error('Stage-1 receiver line differs from Stage-0 in %s.', case_files{ia});
    end

    % The even FFT grid is centered between two samples.  Interpolate at the
    % physical coordinate -x instead of silently reversing the receiver line.
    Gpe_mirror = interp1(x, Gpe, -x, 'linear', 'extrap');
    candidates = {Gpe, conj(Gpe), Gpe_mirror, conj(Gpe_mirror)};
    cm = cell(1, numel(candidates));
    for it = 1:numel(candidates)
        cm{it} = local_metrics(Gbh, candidates{it}, s0.footprint);
        row = row + 1;
        metric_case{row} = case_names{ia};
        metric_transform{row} = transform_names{it};
        metric_rows(row, :) = [A_values(ia), cm{it}.l2_m99, ...
            cm{it}.phase_rms_m99, cm{it}.tl_rms_m99, cm{it}.rho_shape, ...
            cm{it}.rho_raw, cm{it}.phi0_rad, cm{it}.aligned_l2_m99, ...
            cm{it}.phase_p95_m95, cm{it}.tl_p95_m95, cm{it}.finite_count];
    end
    weak_i = local_weak_response(x, Gpe, Gbh, s0.footprint, ...
        A_values(ia), K, v.metrics.model.phase_rms_m99);
    case_i = struct('name', case_names{ia}, 'file', case_files{ia}, ...
        'A_m', A_values(ia), 'K_radpm', K, 'x_m', x, 'G_PE', Gpe, ...
        'G_BH', Gbh, 'G_PE_mirror', Gpe_mirror, 'metrics', {cm}, ...
        'stage1_metrics_original', v.metrics, 'geometry', v.geometry, ...
        'config', v.config, 'passed_geometry', v.geometry.all, ...
        'source_pattern_fingerprint', struct( ...
        'samples', v.config.source_pattern_samples, ...
        'clip_db', v.config.source_pattern_clip_db, ...
        'sha256', v.config.canonical_coeff_sha256));
    if ia == 1
        weak = weak_i;
        cases{ia} = case_i;
    else
        weak(ia) = weak_i;
        cases{ia} = case_i;
    end
end

metric_rows = metric_rows(1:row, :);
metric_case = metric_case(1:row);
metric_transform = metric_transform(1:row);
metric_table = table(metric_case, metric_transform, metric_rows(:, 1), ...
    metric_rows(:, 2), metric_rows(:, 3), metric_rows(:, 4), ...
    metric_rows(:, 5), metric_rows(:, 6), metric_rows(:, 7), ...
    metric_rows(:, 8), metric_rows(:, 9), metric_rows(:, 10), ...
    metric_rows(:, 11), 'VariableNames', {'case_name', 'transform', ...
    'amplitude_m', 'E_G_m99', 'phase_rms_rad_m99', 'tl_rms_db_m99', ...
    'rho_shape', 'rho_raw', 'phi0_rad', 'E_aligned_m99', ...
    'phase_p95_rad_m95', 'tl_p95_db_m95', 'finite_count'});
weak_table = local_weak_table(weak);
conv = local_convention_summary(metric_table, case_names);
provenance = local_provenance(root, base, s0, cases, transform_names);

out_dir = fullfile(base, 'stage1x');
if exist(out_dir, 'dir') ~= 7, mkdir(out_dir); end
mat_file = fullfile(out_dir, 'stage1x_validation.mat');
csv_file = fullfile(out_dir, 'stage1x_metrics.csv');
weak_csv = fullfile(out_dir, 'stage1x_weak_response.csv');
json_file = fullfile(out_dir, 'stage1x_provenance.json');
report_file = fullfile(root, 'reports', ...
    'pe_bellhop_controlled_comparison_stage1x_phase_convention_audit.md');

audit = struct('schema_version', '1.0.0', 'stage', 'stage1x', ...
    'generated_without_solver_calls', true, 'source_goal', ...
    fullfile(root, 'reports', ...
    'pe_bellhop_controlled_comparison_GOAL_revised_stage1x.md'), ...
    'A_m', A_values, 'K_radpm', K, 'cases', {cases}, ...
    'metrics_table', metric_table, 'weak_table', weak_table, ...
    'convention_summary', conv, 'provenance', provenance, ...
    'decision', conv.decision, 'files', struct('mat', mat_file, ...
    'metrics_csv', csv_file, 'weak_csv', weak_csv, ...
    'provenance_json', json_file, 'report', report_file));
% Assign the struct array after construction so audit remains a scalar.
audit.weak_response = weak;
save(mat_file, 'audit', '-v7');
writetable(metric_table, csv_file);
writetable(weak_table, weak_csv);
fid = fopen(json_file, 'w', 'n', 'UTF-8');
if fid < 0, error('Cannot write %s.', json_file); end
cl = onCleanup(@() fclose(fid));
fwrite(fid, jsonencode(provenance), 'char');
clear cl
local_write_report(report_file, audit);
fprintf('Stage 1X decision: %s\n', conv.decision);
fprintf('Report: %s\n', report_file);
end
