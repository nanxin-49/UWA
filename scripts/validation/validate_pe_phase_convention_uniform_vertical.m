function audit = validate_pe_phase_convention_uniform_vertical(overrides)
%VALIDATE_PE_PHASE_CONVENTION_UNIFORM_VERTICAL Audit PE carrier/time signs.
% This validation compares the public reduced and physical fields with
% an independent one-step angular-spectrum reference for the same Gaussian
% source and uses the exp(-i*omega*t) synthesis convention.

if nargin < 1 || isempty(overrides)
    overrides = struct();
end
this_file = mfilename('fullpath');
project_root = fileparts(fileparts(fileparts(this_file)));
addpath(project_root);
setup_vertical_project();

cfg = local_defaults();
cfg = local_overrides(cfg, overrides);
local_validate_cfg(cfg);
out_dir = cfg.output_dir;
if isempty(out_dir)
    out_dir = fullfile(project_root, 'results', 'validation', ...
        'pe_phase_convention_uniform');
end
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

fprintf('PE uniform-medium phase convention audit\n');
f_axis_hz = linspace(cfg.frequency_band_hz(1), ...
    cfg.frequency_band_hz(2), cfg.frequency_count).';
validation_run_meta=cfg.validation_run_meta;
if isempty(validation_run_meta)
    validation_run_meta=pe_phase_release_run_meta_vertical(project_root,struct( ...
        'frequency_axis_hz',f_axis_hz,'seed_definition',struct('deterministic',true)));
end
case_count = numel(cfg.receiver_offsets_m);
case_rows = repmat(local_empty_case_row(), case_count, 1);
phase_records = repmat(struct(), case_count, 1);
channels = cell(case_count, 1);

operator_identity_error = local_operator_identity_error(cfg);
for ii = 1:case_count
    offset_m = cfg.receiver_offsets_m(ii);
    fprintf('Phase audit geometry %d/%d: receiver offset %.3f m\n', ...
        ii, case_count, offset_m);
    paramsV = local_pe_params(cfg, f_axis_hz.', offset_m);
    channel = vertical_channel_model(paramsV);
    channels{ii} = local_compact_channel(channel);

    [H_reduced_reference_f, reference_meta] = ...
        local_angular_spectrum_reference(channel, cfg, offset_m);
    vertical_distance_m = cfg.z_tx_m - cfg.z_rx_m;
    omega = 2 * pi * f_axis_hz;
    carrier_positive = exp(1i * omega * vertical_distance_m / cfg.c0_mps);
    carrier_negative = conj(carrier_positive);
    H_reference_physical_f = H_reduced_reference_f .* carrier_positive;
    candidate = struct();
    candidate.positive_f = channel.H_direct_physical_f(:);
    candidate.none_f = channel.H_direct_reduced_f(:);
    candidate.negative_f = channel.H_direct_reduced_f(:) .* carrier_negative;

    positive_metrics = local_frequency_metrics( ...
        candidate.positive_f, H_reference_physical_f, f_axis_hz, ...
        reference_meta.point_delay_s);
    none_metrics = local_frequency_metrics( ...
        candidate.none_f, H_reference_physical_f, f_axis_hz, ...
        reference_meta.point_delay_s);
    negative_metrics = local_frequency_metrics( ...
        candidate.negative_f, H_reference_physical_f, f_axis_hz, ...
        reference_meta.point_delay_s);

    [pe_peak_time_s, pe_cir] = local_positive_slope_peak( ...
        candidate.positive_f, f_axis_hz, vertical_distance_m / cfg.c0_mps, ...
        reference_meta.point_delay_s, cfg);
    [reference_peak_time_s, reference_cir] = local_positive_slope_peak( ...
        H_reference_physical_f, f_axis_hz, vertical_distance_m / cfg.c0_mps, ...
        reference_meta.point_delay_s, cfg);

    invariant_error = max(abs(channel.H_f(:) - channel.H_direct_f(:) ...
        - channel.H_reflect_f(:)));
    row = local_empty_case_row();
    row.offset_m = offset_m;
    row.rx_grid_x_m = reference_meta.rx_grid_x_m;
    row.point_delay_ms = 1000 * reference_meta.point_delay_s;
    row.positive_phase_rms_rad = positive_metrics.phase_rms_rad;
    row.none_phase_rms_rad = none_metrics.phase_rms_rad;
    row.negative_phase_rms_rad = negative_metrics.phase_rms_rad;
    row.positive_complex_rel_l2 = positive_metrics.complex_rel_l2;
    row.pe_group_delay_ms = 1000 * positive_metrics.group_delay_s;
    row.reference_group_delay_ms = 1000 * positive_metrics.reference_group_delay_s;
    row.group_delay_difference_ms = 1000 * positive_metrics.group_delay_difference_s;
    row.pe_peak_time_ms = 1000 * pe_peak_time_s;
    row.reference_peak_time_ms = 1000 * reference_peak_time_s;
    row.peak_time_difference_ms = 1000 * (pe_peak_time_s - reference_peak_time_s);
    row.channel_invariant_error = invariant_error;
    row.max_abs_reflection = max(abs(channel.H_reflect_f(:)));
    row.pass_1_over_r = logical(channel.pass_1_over_R);
    case_rows(ii) = row;

    phase_records(ii).offset_m = offset_m;
    phase_records(ii).f_axis_hz = f_axis_hz;
    phase_records(ii).positive_residual_phase_rad = positive_metrics.residual_phase_rad;
    phase_records(ii).none_residual_phase_rad = none_metrics.residual_phase_rad;
    phase_records(ii).negative_residual_phase_rad = negative_metrics.residual_phase_rad;
    phase_records(ii).H_pe_reduced_f = channel.H_direct_reduced_f(:);
    phase_records(ii).H_reference_reduced_f = H_reduced_reference_f;
    phase_records(ii).H_pe_physical_f = candidate.positive_f;
    phase_records(ii).H_reference_physical_f = H_reference_physical_f;
    phase_records(ii).pe_cir = pe_cir;
    phase_records(ii).reference_cir = reference_cir;
    phase_records(ii).reference_meta = reference_meta;
end

case_table = struct2table(case_rows);
checks = local_checks(case_table, operator_identity_error, cfg);
passed = all(checks.passed);

csv_file = fullfile(out_dir, 'pe_phase_convention_case_summary.csv');
mat_file = fullfile(out_dir, 'pe_phase_convention_audit.mat');
figure_file = fullfile(out_dir, 'pe_phase_convention_audit.png');
report_file = fullfile(out_dir, 'pe_phase_convention_audit_report.md');
writetable(case_table, csv_file);
local_plot(phase_records, case_table, figure_file);

audit = struct();
audit.config = cfg;
audit.operator_identity_error = operator_identity_error;
audit.selected_carrier_sign = +1;
audit.time_synthesis_convention = 'exp(-i*omega*t); positive phase slope is positive delay';
audit.validation_cir_transform = 'fft after removing positive reference-delay phase';
audit.case_table = case_table;
audit.phase_records = phase_records;
audit.channels = channels;
audit.checks = checks;
audit.passed = passed;
audit.schema_version='2.0.0';
audit.validation_run_meta=validation_run_meta;
phase_geometry=struct('z_tx',cfg.z_tx_m,'z_rx',cfg.z_rx_m,'z_surface',0,'c0',cfg.c0_mps);
[~,audit.phase_reference_meta]=apply_pe_channel_phase_reference_vertical( ...
    f_axis_hz,struct('direct_f',ones(size(f_axis_hz)), ...
    'reflect_fm',ones(size(f_axis_hz))),phase_geometry,'direct_dsp');
audit.files = struct('csv', csv_file, 'mat', mat_file, ...
    'figure', figure_file, 'report', report_file);
schema_version=audit.schema_version;
phase_reference_meta=audit.phase_reference_meta;
save(mat_file,'audit','schema_version','phase_reference_meta','validation_run_meta');
local_write_report(report_file, audit);

disp(case_table);
disp(checks);
if ~passed
    failed = strjoin(cellstr(checks.check_name(~checks.passed)), ', ');
    error('validate_pe_phase_convention_uniform_vertical:Failed', ...
        'Phase convention audit failed: %s', failed);
end
fprintf('Phase convention audit passed: %s\n', report_file);
end

function cfg = local_defaults()
cfg = struct();
cfg.c0_mps = 1500;
cfg.water_depth_m = 100;
cfg.z_tx_m = 80;
cfg.z_rx_m = 10;
cfg.receiver_offsets_m = [3, 6, 9];
cfg.frequency_band_hz = [3000, 5000];
cfg.frequency_count = 65;
cfg.f_ref_hz = 4000;
cfg.nx = 128;
cfg.ny = 128;
cfg.xw_m = 32;
cfg.yw_m = 32;
cfg.stepz_lamb = 0.5;
cfg.sigma_src_m = 0.3;
cfg.sponge_ratio = 0.12;
% The operator audit uses the periodic transverse angular-spectrum problem.
% Absorbing layers are intentionally disabled here because their repeated
% real-space multiplication is not part of the one-step reference operator.
% The public API requires a strictly positive value, so use a numerically
% negligible coefficient instead of changing its validation contract.
cfg.alpha_max_np_per_m = 1e-14;
cfg.cir_window = 'hann';
cfg.cir_zero_padding_factor = 8;
cfg.operator_identity_tolerance = 1e-12;
cfg.phase_rms_tolerance_rad = 0.15;
cfg.convention_margin_factor = 3;
cfg.group_delay_tolerance_ms = 0.25;
cfg.peak_time_tolerance_ms = 0.25;
cfg.invariant_tolerance = 1e-10;
cfg.output_dir = '';
cfg.validation_run_meta = [];
end

function cfg = local_overrides(cfg, overrides)
if ~isstruct(overrides), error('overrides must be a struct.'); end
names = fieldnames(overrides);
for ii = 1:numel(names)
    if ~isfield(cfg, names{ii}), error('Unknown override: %s', names{ii}); end
    cfg.(names{ii}) = overrides.(names{ii});
end
end

function local_validate_cfg(cfg)
if cfg.c0_mps <= 0 || cfg.water_depth_m <= 0
    error('c0_mps and water_depth_m must be positive.');
end
if cfg.z_rx_m < 0 || cfg.z_rx_m >= cfg.z_tx_m || cfg.z_tx_m > cfg.water_depth_m
    error('Require 0 <= z_rx_m < z_tx_m <= water_depth_m.');
end
if any(cfg.receiver_offsets_m <= 0) || any(cfg.receiver_offsets_m >= cfg.xw_m / 2)
    error('Receiver offsets must be positive and inside the transverse window.');
end
if mod(cfg.frequency_count, 2) ~= 1 || cfg.frequency_count < 3
    error('frequency_count must be odd and at least 3.');
end
if abs(mean(cfg.frequency_band_hz) - cfg.f_ref_hz) > 1e-9
    error('f_ref_hz must be the frequency-band midpoint.');
end
end

function paramsV = local_pe_params(cfg, f_axis_hz, offset_m)
paramsV = struct();
paramsV.f0 = f_axis_hz;
paramsV.enable_wideband = false;
paramsV.f_ref_hz = cfg.f_ref_hz;
paramsV.c0 = cfg.c0_mps;
paramsV.z_max = cfg.water_depth_m;
paramsV.stepz_lamb = cfg.stepz_lamb;
paramsV.xw = cfg.xw_m;
paramsV.yw = cfg.yw_m;
paramsV.nx = cfg.nx;
paramsV.ny = cfg.ny;
paramsV.x_tx = 0;
paramsV.y_tx = 0;
paramsV.z_tx = cfg.z_tx_m;
paramsV.x_rx = offset_m;
paramsV.y_rx = 0;
paramsV.z_rx = cfg.z_rx_m;
paramsV.sigma_src_m = cfg.sigma_src_m;
paramsV.sponge_ratio = cfg.sponge_ratio;
paramsV.alpha_max_np_per_m = cfg.alpha_max_np_per_m;
paramsV.env_mode = 'uniform';
paramsV.show_figures = false;
paramsV.enforce_1_over_R = true;
paramsV.enable_surface_reflection = false;
paramsV.sea_hs_target = 0;
paramsV.surface_reflect_coeff = -1;
paramsV.enable_bubbles = false;
paramsV.save_mode = 'rx_only';
paramsV.use_gpu = false;
end

function error_value = local_operator_identity_error(cfg)
k0 = 2 * pi * cfg.f_ref_hz / cfg.c0_mps;
kappa = linspace(0, 1.5 * k0, 4097).';
kz = sqrt(complex(k0^2 - kappa.^2, 0));
d = cfg.z_tx_m - cfg.z_rx_m;
lhs = exp(-1i * d * kappa.^2 ./ (kz + k0));
rhs = exp(1i * d * (kz - k0));
error_value = max(abs(lhs - rhs));
end

function [H_reduced_f, meta] = local_angular_spectrum_reference(channel, cfg, offset_m)
x = channel.x(:).';
y = channel.y(:).';
[X, Y] = meshgrid(x, y);
psi0 = exp(-(X.^2 + Y.^2) / (2 * cfg.sigma_src_m^2));
psi0_k = fft2(psi0);
kx = (2 * pi / cfg.xw_m) * [0:(cfg.nx/2-1), -cfg.nx/2:-1];
ky = (2 * pi / cfg.yw_m) * [0:(cfg.ny/2-1), -cfg.ny/2:-1];
[KX, KY] = meshgrid(kx, ky);
kappa2 = KX.^2 + KY.^2;
[~, ix_rx] = min(abs(x - offset_m));
[~, iy_rx] = min(abs(y));
d = cfg.z_tx_m - cfg.z_rx_m;
f = channel.f_axis(:);
H_reduced_f = complex(zeros(size(f)));
for jj = 1:numel(f)
    k0 = 2 * pi * f(jj) / cfg.c0_mps;
    kz = sqrt(complex(k0^2 - kappa2, 0));
    reduced_factor = exp(1i * d * (kz - k0));
    psi_reference = ifft2(psi0_k .* reduced_factor);
    H_reduced_f(jj) = psi_reference(iy_rx, ix_rx);
end
rx_grid_x_m = x(ix_rx);
point_range_m = hypot(d, rx_grid_x_m);
meta = struct('rx_grid_x_m', rx_grid_x_m, 'rx_grid_y_m', y(iy_rx), ...
    'vertical_distance_m', d, 'point_range_m', point_range_m, ...
    'point_delay_s', point_range_m / cfg.c0_mps, ...
    'reference_model', 'one-step exact discrete angular spectrum without sponge');
end

function metrics = local_frequency_metrics(H_candidate, H_reference, f_hz, expected_delay_s)
[~, idx_ref] = min(abs(f_hz - mean([f_hz(1), f_hz(end)])));
candidate_norm = H_candidate / H_candidate(idx_ref);
reference_norm = H_reference / H_reference(idx_ref);
residual_phase = unwrap(angle(candidate_norm .* conj(reference_norm)));
residual_phase = residual_phase - residual_phase(idx_ref);
metrics = struct();
metrics.residual_phase_rad = residual_phase;
metrics.phase_rms_rad = sqrt(mean(residual_phase.^2));
metrics.complex_rel_l2 = norm(candidate_norm - reference_norm) / max(norm(reference_norm), eps);
metrics.group_delay_s = local_unwrapped_group_delay(f_hz, H_candidate, expected_delay_s);
metrics.reference_group_delay_s = local_unwrapped_group_delay( ...
    f_hz, H_reference, expected_delay_s);
metrics.group_delay_difference_s = metrics.group_delay_s - metrics.reference_group_delay_s;
end

function delay_s = local_unwrapped_group_delay(f_hz, H_f, expected_delay_s)
phase = unwrap(angle(H_f));
p = polyfit(f_hz, phase, 1);
raw_delay_s = p(1) / (2 * pi);
unambiguous_s = 1 / mean(diff(f_hz));
delay_s = raw_delay_s + round((expected_delay_s - raw_delay_s) / unambiguous_s) * unambiguous_s;
end

function [peak_time_s, cir] = local_positive_slope_peak( ...
    H_f, f_hz, reference_delay_s, expected_time_s, cfg)
F = numel(f_hz);
df = mean(diff(f_hz));
Nfft = cfg.cir_zero_padding_factor * F;
n = (0:F-1).';
if strcmpi(cfg.cir_window, 'hann')
    window = 0.5 - 0.5 * cos(2 * pi * n / (F - 1));
else
    window = ones(F, 1);
end
shifted = H_f(:) .* exp(-1i * 2 * pi * f_hz(:) * reference_delay_s) .* window;
h = fft(shifted, Nfft, 1);
delay_residual_s = (0:Nfft-1).' / (df * Nfft);
delay_absolute_s = reference_delay_s + delay_residual_s;
search_half_width_s = 0.002;
mask = abs(delay_absolute_s - expected_time_s) <= search_half_width_s;
indices = find(mask);
[~, local_idx] = max(abs(h(mask)).^2);
peak_time_s = delay_absolute_s(indices(local_idx));
cir = struct('h', h, 'delay_absolute_s', delay_absolute_s, ...
    'delay_residual_s', delay_residual_s, 'reference_delay_s', reference_delay_s, ...
    'transform', 'fft', 'time_convention', 'exp(-i*omega*t)', ...
    'window', cfg.cir_window, 'zero_padding_factor', cfg.cir_zero_padding_factor);
end

function row = local_empty_case_row()
row = struct('offset_m', NaN, 'rx_grid_x_m', NaN, 'point_delay_ms', NaN, ...
    'positive_phase_rms_rad', NaN, 'none_phase_rms_rad', NaN, ...
    'negative_phase_rms_rad', NaN, 'positive_complex_rel_l2', NaN, ...
    'pe_group_delay_ms', NaN, 'reference_group_delay_ms', NaN, ...
    'group_delay_difference_ms', NaN, 'pe_peak_time_ms', NaN, ...
    'reference_peak_time_ms', NaN, 'peak_time_difference_ms', NaN, ...
    'channel_invariant_error', NaN, 'max_abs_reflection', NaN, ...
    'pass_1_over_r', false);
end

function checks = local_checks(t, operator_error, cfg)
dominance_none = max(3 * t.positive_phase_rms_rad ./ max(t.none_phase_rms_rad, eps));
dominance_negative = max(3 * t.positive_phase_rms_rad ./ max(t.negative_phase_rms_rad, eps));
names = ["operator_identity"; "positive_phase_rms"; "positive_vs_none_margin"; ...
    "positive_vs_negative_margin"; "group_delay"; "cir_peak_time"; ...
    "channel_invariant"; "direct_only_reflection"; "one_over_r"];
values = [operator_error; max(t.positive_phase_rms_rad); dominance_none; ...
    dominance_negative; max(abs(t.group_delay_difference_ms)); ...
    max(abs(t.peak_time_difference_ms)); max(t.channel_invariant_error); ...
    max(t.max_abs_reflection); double(~all(t.pass_1_over_r))];
limits = [cfg.operator_identity_tolerance; cfg.phase_rms_tolerance_rad; 1; 1; ...
    cfg.group_delay_tolerance_ms; cfg.peak_time_tolerance_ms; ...
    cfg.invariant_tolerance; cfg.invariant_tolerance; 0];
relations = ["<="; "<="; "<="; "<="; "<="; "<="; "<="; "<="; "=="];
passed = (relations == "<=" & values <= limits) | (relations == "==" & values == limits);
checks = table(names, values, relations, limits, passed, ...
    'VariableNames', {'check_name', 'value', 'relation', 'limit', 'passed'});
end

function compact = local_compact_channel(channel)
compact = struct('H_direct_f', channel.H_direct_f, ...
    'H_reflect_f', channel.H_reflect_f, 'H_f', channel.H_f, ...
    'f_axis', channel.f_axis, 'idx_f_ref', channel.idx_f_ref, ...
    'fit_slope', channel.fit_slope, 'fit_err_rms', channel.fit_err_rms, ...
    'pass_1_over_R', channel.pass_1_over_R, 'rx_state_used', channel.rx_state_used, ...
    'config', channel.config);
end

function local_plot(records, case_table, figure_file)
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1050 760]);
cleanup = onCleanup(@() close(fig));
subplot(2, 2, 1);
hold on;
for ii = 1:numel(records)
    plot(records(ii).f_axis_hz, records(ii).positive_residual_phase_rad, 'LineWidth', 1.2);
end
grid on; xlabel('Frequency (Hz)'); ylabel('Residual phase (rad)');
title('Positive-carrier residual');
legend(compose('x=%.0f m', case_table.offset_m), 'Location', 'best');

subplot(2, 2, 2);
bar(categorical(compose('x=%.0f', case_table.offset_m)), ...
    [case_table.positive_phase_rms_rad, case_table.none_phase_rms_rad, ...
     case_table.negative_phase_rms_rad]);
grid on; ylabel('Phase RMS (rad)'); title('Carrier-sign candidates');
legend('+ carrier', 'none', '- carrier', 'Location', 'best');

subplot(2, 2, 3);
bar(categorical(compose('x=%.0f', case_table.offset_m)), ...
    [case_table.pe_group_delay_ms, case_table.reference_group_delay_ms]);
grid on; ylabel('Group delay (ms)'); title('Frequency-slope delay');
legend('PE', 'Angular reference', 'Location', 'best');

subplot(2, 2, 4);
bar(categorical(compose('x=%.0f', case_table.offset_m)), ...
    [case_table.pe_peak_time_ms, case_table.reference_peak_time_ms]);
grid on; ylabel('CIR peak time (ms)'); title('Validation-local FFT CIR');
legend('PE', 'Angular reference', 'Location', 'best');
exportgraphics(fig, figure_file, 'Resolution', 180);
clear cleanup
end

function local_write_report(report_file, audit)
t = audit.case_table;
c = audit.checks;
fid = fopen(report_file, 'w', 'n', 'UTF-8');
if fid < 0, error('Cannot create report: %s', report_file); end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '# PE 均匀介质相位与载波约定审计\n\n');
fprintf(fid, '## 结论\n\n');
fprintf(fid, ['审计状态：`%s`。该脚本读取公共 PE 的显式 reduced/physical 字段，', ...
    '并独立验证载波恢复符号；载波相位统一由公共相位参考层完成。\n\n'], ...
    string(audit.passed));
fprintf(fid, ['PE 单步约化传播算子对应 `exp(i*d*(kz-k0))`；在 `exp(-i*omega*t)` ', ...
    '时间约定下，验证层选择正号参考载波 `exp(+i*k0*d)`，并用 FFT 将正相位斜率映射到正时延。\n\n']);
fprintf(fid, '算子恒等式最大误差：`%.6g`。\n\n', audit.operator_identity_error);
fprintf(fid, '## 几何结果\n\n');
fprintf(fid, '| x (m) | +载波相位RMS | 无载波相位RMS | -载波相位RMS | PE群时延(ms) | 参考群时延(ms) | CIR峰差(ms) |\n');
fprintf(fid, '|---:|---:|---:|---:|---:|---:|---:|\n');
for ii = 1:height(t)
    fprintf(fid, '| %.3f | %.6g | %.6g | %.6g | %.6f | %.6f | %+.6f |\n', ...
        t.offset_m(ii), t.positive_phase_rms_rad(ii), t.none_phase_rms_rad(ii), ...
        t.negative_phase_rms_rad(ii), t.pe_group_delay_ms(ii), ...
        t.reference_group_delay_ms(ii), t.peak_time_difference_ms(ii));
end
fprintf(fid, '\n## 自动检查\n\n');
fprintf(fid, '| 检查 | 数值 | 条件 | 阈值 | 通过 |\n|---|---:|:---:|---:|:---:|\n');
for ii = 1:height(c)
    fprintf(fid, '| %s | %.8g | %s | %.8g | %d |\n', char(c.check_name(ii)), ...
        c.value(ii), char(c.relation(ii)), c.limit(ii), c.passed(ii));
end
fprintf(fid, '\n![相位约定审计](pe_phase_convention_audit.png)\n\n');
fprintf(fid, ['## 限制\n\n角谱参考使用相同离散高斯初场和无限周期横向域；本审计把 PE 吸收系数设为 ', ...
    '`1e-14`，使海绵层在数值上可忽略而仍满足公共 API 的正数约束。', ...
    '表中的频率斜率群时延包含有限宽高斯波束的频率相关衍射，不应单独当作点源几何时延；', ...
    '本审计判断的是 PE 与同初场角谱参考的闭合，点路径位置另由 CIR 峰与解析解检查。', ...
    '本报告内部的 FFT 只用于 `exp(-i*omega*t)` 物理相量审计；', ...
    '面向 MATLAB IFFT 的公共信道应使用 `direct_dsp` 字段。\n']);
clear cleanup
end
