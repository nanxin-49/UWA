function validation = validate_pe_bellhop_pm_canonical_mapper(overrides)
%VALIDATE_PE_BELLHOP_PM_CANONICAL_MAPPER Stage 0A fixed-profile audit.
%   This validation-only entrypoint verifies that one coefficient-defined PM
%   function can be evaluated identically for PE and Bellhop.  It does not
%   run either propagation model.

if nargin < 1 || isempty(overrides), overrides = struct(); end
if ~isstruct(overrides) || ~isscalar(overrides)
    error('overrides must be a scalar struct.');
end
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
setup_vertical_project();

cfg = local_config(root, overrides);
if ~exist(cfg.output_dir, 'dir'), mkdir(cfg.output_dir); end
profile = load_fixed_pm_profile_for_pe_bellhop_validation(cfg.coeff_file, cfg.profile_options);

% Master samples are an independent audit of the Fourier evaluator; they are
% not the canonical source for new solver inputs.
master = readtable(cfg.master_profile_file);
required_master = {'s_m','eta_m','eta_prime_m_per_m','eta_second_m_per_m2'};
for ii = 1:numel(required_master)
    if ~ismember(required_master{ii}, master.Properties.VariableNames)
        error('Master profile is missing column %s.', required_master{ii});
    end
end
master_samples = evaluate_fixed_pm_fourier_profile(profile, double(master.s_m));
master_eta_err = max(abs(master_samples.eta_m - double(master.eta_m)));
master_slope_err = max(abs(master_samples.eta_prime_m_per_m - double(master.eta_prime_m_per_m)));
master_second_derivative_err = max(abs(master_samples.eta_second_m_per_m2 - ...
    double(master.eta_second_m_per_m2)));

% PE and Bellhop receive different solver grids, but both are evaluated from
% the same immutable continuous Fourier series.
s_pe = linspace(cfg.pe_interval_m(1), cfg.pe_interval_m(2), cfg.pe_sample_count).';
s_bh = linspace(cfg.bellhop_interval_m(1), cfg.bellhop_interval_m(2), cfg.bellhop_sample_count).';
pe = evaluate_fixed_pm_fourier_profile(profile, s_pe);
bh = evaluate_fixed_pm_fourier_profile(profile, s_bh);
pe_roundtrip = evaluate_fixed_pm_fourier_profile(profile, pe.s_m);
bh_roundtrip = evaluate_fixed_pm_fourier_profile(profile, bh.s_m);

% Exact-90-degree chart convention: Gamma_B(s)=[R0-eta(s),s].
wall_r_m = cfg.r0_m - bh.eta_m;
wall_z_m = bh.s_m;
eta_recovered_from_wall = cfg.r0_m - wall_r_m;
wall_roundtrip_err = max(abs(eta_recovered_from_wall - bh.eta_m));

stats = local_statistics(profile, master_samples);
checks = struct();
checks.master_height = master_eta_err <= cfg.master_value_tolerance;
checks.master_slope = master_slope_err <= cfg.master_derivative_tolerance;
checks.master_second_derivative = master_second_derivative_err <= cfg.master_second_derivative_tolerance;
checks.pe_finite = all(isfinite(pe.eta_m(:))) && all(isfinite(pe.eta_prime_m_per_m(:))) && ...
    all(isfinite(pe.eta_second_m_per_m2(:)));
checks.bellhop_finite = all(isfinite(bh.eta_m(:))) && all(isfinite(bh.eta_prime_m_per_m(:))) && ...
    all(isfinite(bh.eta_second_m_per_m2(:)));
checks.pe_inverse = max(abs(pe_roundtrip.eta_m - pe.eta_m)) <= cfg.evaluator_tolerance;
checks.bellhop_inverse = max(abs(bh_roundtrip.eta_m - bh.eta_m)) <= cfg.evaluator_tolerance;
checks.wall_mapping = wall_roundtrip_err <= cfg.evaluator_tolerance;
checks.all = local_all_checks(checks);

validation = struct();
validation.schema_version = '1.0.0';
validation.stage = '0A_canonical_fixed_pm_mapper';
validation.config = cfg;
validation.profile = profile;
validation.master_profile = struct('file',cfg.master_profile_file, ...
    'sha256',local_sha256_file(cfg.master_profile_file), ...
    'sample_count',height(master), 'samples',master_samples);
validation.pe_mapping = struct('s_m',pe.s_m,'eta_m',pe.eta_m, ...
    'eta_prime_m_per_m',pe.eta_prime_m_per_m, ...
    'eta_second_m_per_m2',pe.eta_second_m_per_m2, ...
    'profile_hash_sha256',pe.profile_hash_sha256);
validation.bellhop_mapping = struct('s_m',bh.s_m,'eta_m',bh.eta_m, ...
    'wall_r_m',wall_r_m,'wall_z_m',wall_z_m, ...
    'eta_prime_m_per_m',bh.eta_prime_m_per_m, ...
    'eta_second_m_per_m2',bh.eta_second_m_per_m2, ...
    'profile_hash_sha256',bh.profile_hash_sha256);
validation.statistics = stats;
validation.errors = struct('master_height_max_m',master_eta_err, ...
    'master_slope_max_m_per_m',master_slope_err, ...
    'master_second_derivative_max_m_per_m2',master_second_derivative_err, ...
    'bh_roundtrip_height_max_m',wall_roundtrip_err);
validation.checks = checks;
validation.passed = checks.all;

save(fullfile(cfg.output_dir, 'canonical_mapper_validation.mat'), 'validation', '-v7');
local_write_csv(cfg.output_dir, validation);
local_write_report(cfg.report_path, validation);
if cfg.fail_on_check && ~validation.passed
    error('Canonical PM mapper validation failed; see %s.', cfg.report_path);
end
end

function cfg = local_config(root, overrides)
cfg = struct();
cfg.coeff_file = fullfile(root, 'results', 'validation', ...
    'bellhop_internal_pm_fixed_realization', 'fixed_pm_fourier_coefficients.csv');
cfg.master_profile_file = fullfile(root, 'results', 'validation', ...
    'bellhop_internal_pm_fixed_realization', 'fixed_pm_master_profile.csv');
cfg.output_dir = fullfile(root, 'results', 'validation', 'pe_bellhop_pm_canonical_mapper');
cfg.report_path = fullfile(root, 'reports', 'pe_bellhop_pm_canonical_mapper_report.md');
cfg.r0_m = 100;
cfg.pe_interval_m = [-96.09375, 96.09375];
cfg.bellhop_interval_m = [-80, 80];
cfg.pe_sample_count = 984;
cfg.bellhop_sample_count = 4097;
cfg.master_value_tolerance = 5e-12;
cfg.master_derivative_tolerance = 5e-11;
cfg.master_second_derivative_tolerance = 5e-10;
cfg.evaluator_tolerance = 5e-14;
cfg.fail_on_check = true;
cfg.profile_options = struct('seed',260001,'wind_speed_mps',6, ...
    'requested_kmax_rad_per_m',0.5,'datum_m',0);
names = fieldnames(overrides);
for ii = 1:numel(names)
    name = names{ii};
    if ~isfield(cfg, name), error('Unknown mapper override: %s', name); end
    cfg.(name) = overrides.(name);
end
end

function stats = local_statistics(profile, samples)
eta = samples.eta_m(:);
eta_prime = samples.eta_prime_m_per_m(:);
eta_second = samples.eta_second_m_per_m2(:);
kappa = eta_second ./ (1 + eta_prime.^2).^(3/2);
stats = struct('mean_height_m',mean(eta), ...
    'rms_height_m',sqrt(mean(eta.^2)), ...
    'min_height_m',min(eta), 'max_height_m',max(eta), ...
    'rms_slope',sqrt(mean(eta_prime.^2)), ...
    'max_abs_slope',max(abs(eta_prime)), ...
    'rms_geometric_curvature_per_m',sqrt(mean(kappa.^2)), ...
    'max_abs_geometric_curvature_per_m',max(abs(kappa)), ...
    'minimum_radius_of_curvature_m',1/max(max(abs(kappa)), eps), ...
    'endpoint_height_difference_m',eta(end)-eta(1), ...
    'endpoint_slope_difference',eta_prime(end)-eta_prime(1), ...
    'endpoint_second_derivative_difference',eta_second(end)-eta_second(1), ...
    'source_span_m',profile.span_m, 'coefficient_count',profile.coefficient_count);
end

function tf = local_all_checks(checks)
names = fieldnames(checks);
tf = true;
for ii = 1:numel(names)
    if ~strcmp(names{ii}, 'all'), tf = tf && checks.(names{ii}); end
end
end

function local_write_csv(output_dir, validation)
rows = struct('quantity',{},'value',{},'unit',{});
rows(end+1) = struct('quantity','master_height_max_error','value',validation.errors.master_height_max_m,'unit','m');
rows(end+1) = struct('quantity','master_slope_max_error','value',validation.errors.master_slope_max_m_per_m,'unit','1');
rows(end+1) = struct('quantity','master_second_derivative_max_error','value',validation.errors.master_second_derivative_max_m_per_m2,'unit','1/m');
rows(end+1) = struct('quantity','wall_mapping_roundtrip_max_error','value',validation.errors.bh_roundtrip_height_max_m,'unit','m');
rows(end+1) = struct('quantity','rms_height','value',validation.statistics.rms_height_m,'unit','m');
rows(end+1) = struct('quantity','rms_slope','value',validation.statistics.rms_slope,'unit','1');
rows(end+1) = struct('quantity','rms_curvature','value',validation.statistics.rms_geometric_curvature_per_m,'unit','1/m');
rows(end+1) = struct('quantity','max_curvature','value',validation.statistics.max_abs_geometric_curvature_per_m,'unit','1/m');
writetable(struct2table(rows), fullfile(output_dir,'canonical_mapper_metrics.csv'));
end

function local_write_report(path, v)
fid = fopen(path, 'w');
if fid < 0, error('Cannot write report: %s', path); end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '# PE--Bellhop PM canonical mapper Stage 0A 报告\n\n');
fprintf(fid, '状态：**%s**\n\n', ternary(v.passed,'PASS','FAIL'));
fprintf(fid, '本阶段只验证固定 Fourier PM realization 的跨 solver 数学映射，没有运行 PE/Bellhop propagation，也没有修改任何核心物理。\n\n');
fprintf(fid, '## 输入与映射\n\n');
fprintf(fid, '- coefficients：%s\n', v.profile.coeff_file);
fprintf(fid, '- coefficients SHA-256：%s\n', v.profile.coeff_file_sha256);
fprintf(fid, '- seed/U/span：%d / %.6g m/s / %.6g m\n', v.profile.seed, v.profile.wind_speed_mps, v.profile.span_m);
fprintf(fid, '- coefficient count / requested Kmax / realized Kmax：%d / %.12g / %.12g rad/m\n', v.profile.coefficient_count, v.profile.requested_kmax_rad_per_m, v.profile.realized_kmax_rad_per_m);
fprintf(fid, '- series：%s\n', v.profile.series_formula);
fprintf(fid, '- master profile SHA-256：%s\n', v.master_profile.sha256);
fprintf(fid, '- PE evaluation interval/count：[%g, %g] m / %d\n', v.config.pe_interval_m(1), v.config.pe_interval_m(2), v.config.pe_sample_count);
fprintf(fid, '- Bellhop evaluation interval/count：[%g, %g] m / %d\n', v.config.bellhop_interval_m(1), v.config.bellhop_interval_m(2), v.config.bellhop_sample_count);
fprintf(fid, '- PE/Bellhop sampled profile hashes：%s / %s\n', v.pe_mapping.profile_hash_sha256, v.bellhop_mapping.profile_hash_sha256);
fprintf(fid, '- sampled hashes intentionally differ because the solver grids differ; the common canonical identity is the coefficient-file SHA-256 above.\n');
fprintf(fid, '- PE mapping：eta_PE(x,y)=eta_1D(x) on the configured x samples.\n');
fprintf(fid, '- Bellhop mapping：Gamma_B(s)=[R0-eta(s),s], R0=%.6g m; no sign change or renormalization.\n\n', v.config.r0_m);
fprintf(fid, '## 误差与统计\n\n');
fprintf(fid, '| quantity | value |\n|---|---:|\n');
fprintf(fid, '| master height max error | %.6g m |\n', v.errors.master_height_max_m);
fprintf(fid, '| master slope max error | %.6g |\n', v.errors.master_slope_max_m_per_m);
fprintf(fid, '| master second-derivative max error | %.6g 1/m |\n', v.errors.master_second_derivative_max_m_per_m2);
fprintf(fid, '| Bellhop wall inverse-map error | %.6g m |\n', v.errors.bh_roundtrip_height_max_m);
fprintf(fid, '| RMS height | %.9g m |\n', v.statistics.rms_height_m);
fprintf(fid, '| RMS/max slope | %.9g / %.9g |\n', v.statistics.rms_slope, v.statistics.max_abs_slope);
fprintf(fid, '| RMS/max curvature | %.9g / %.9g 1/m |\n', v.statistics.rms_geometric_curvature_per_m, v.statistics.max_abs_geometric_curvature_per_m);
fprintf(fid, '| minimum radius | %.9g m |\n', v.statistics.minimum_radius_of_curvature_m);
fprintf(fid, '| endpoint height/slope/second-derivative differences | %.6g / %.6g / %.6g |\n', v.statistics.endpoint_height_difference_m, v.statistics.endpoint_slope_difference, v.statistics.endpoint_second_derivative_difference);
fprintf(fid, '\n## Hard checks\n\n');
check_names = fieldnames(v.checks);
for ii = 1:numel(check_names)
    if strcmp(check_names{ii}, 'all'), continue; end
    fprintf(fid, '- %s: %s\n', check_names{ii}, ternary(v.checks.(check_names{ii}),'PASS','FAIL'));
end
fprintf(fid, '\n## 结论\n\n');
fprintf(fid, '同一组 Fourier coefficients 在独立 PE/Bellhop 采样点上直接求值；Bellhop wall 的 r=R0-eta 反变换误差为 %.6g m。本阶段不证明两种传播模型的场一致性，只冻结后续 Stage 0B--1 的 canonical 输入。\n', v.errors.bh_roundtrip_height_max_m);
clear cleanup
end

function out = ternary(condition, yes_value, no_value)
if condition, out = yes_value; else, out = no_value; end
end

function digest = local_sha256_file(path)
md = java.security.MessageDigest.getInstance('SHA-256');
fid = fopen(path, 'rb');
if fid < 0, error('Cannot open file for hashing: %s', path); end
cleanup = onCleanup(@() fclose(fid));
bytes = fread(fid, Inf, '*uint8');
clear cleanup
md.update(typecast(bytes, 'int8'));
digest_bytes = typecast(md.digest(), 'uint8');
digest = lower(reshape(dec2hex(digest_bytes, 2).', 1, []));
end
