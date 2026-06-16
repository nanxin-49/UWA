% Quick reduced-grid validation for the SSA-like statistical surface kernel.
% This script checks interface invariants, Hs=0 flat-surface degeneration,
% energy limiting, deterministic seeding, seed sensitivity, and the
% metadata-only random_scatter=false path.

clear
format compact

result_file = getenv('SSA_STAT_VALIDATE_RESULT_FILE');
if isempty(result_file)
    result_file = 'validate_ssa_stat_kernel_vertical_result.mat';
end
if exist(result_file, 'file')
    warning('Result file %s already exists and will be overwritten.', result_file);
end

tol_invariant = 1e-10;
tol_roundoff = 1e-12;

case_results = struct();
checks = struct([]);

params_base = local_base_params();

flat_ssa = params_base;
flat_ssa.surface_boundary_model = 'ssa_stat_kernel';
flat_ssa.sea_hs_target = 0;
flat_ssa.sea_seed = 12345;

flat_kirchhoff = params_base;
flat_kirchhoff.surface_boundary_model = 'kirchhoff_spatial';
flat_kirchhoff.sea_hs_target = 0;
flat_kirchhoff.sea_seed = 12345;

fprintf('Running Hs=0 ssa_stat_kernel flat-degeneration case.\n');
channel_flat_ssa = CARPE3D_vertical(flat_ssa);
fprintf('Running Hs=0 kirchhoff_spatial flat-reference case.\n');
channel_flat_kirchhoff = CARPE3D_vertical(flat_kirchhoff);
meta_flat = channel_flat_ssa.roughness_meta.ssa_stat_kernel_meta;
case_results.flat_ssa = local_compact_channel_result(channel_flat_ssa);
case_results.flat_kirchhoff = local_compact_channel_result(channel_flat_kirchhoff);
case_results.flat_diff_max_abs = max(abs(channel_flat_ssa.H_f(:) - channel_flat_kirchhoff.H_f(:)));

checks = local_add_check(checks, 'Hs0_invariant', ...
    local_invariant_error(channel_flat_ssa), tol_invariant, '<=');
checks = local_add_check(checks, 'Hs0_sigma_eta_zero', ...
    abs(meta_flat.sigma_eta_m), tol_roundoff, '<=');
checks = local_add_check(checks, 'Hs0_R_coh_equals_R0', ...
    abs(meta_flat.R_coh - flat_ssa.surface_reflect_coeff), tol_roundoff, '<=');
checks = local_add_check(checks, 'Hs0_kernel_mode_pm_convolution', ...
    double(strcmp(meta_flat.kernel_mode, 'pm_convolution')), 1, '==');
checks = local_add_check(checks, 'Hs0_P_sca_zero', ...
    abs(meta_flat.P_sca.sum), tol_roundoff, '<=');
checks = local_add_check(checks, 'Hs0_E_sca_zero', ...
    abs(meta_flat.E_sca), tol_roundoff, '<=');
checks = local_add_check(checks, 'Hs0_flat_model_diff', ...
    case_results.flat_diff_max_abs, tol_invariant, '<=');

flat_ssa1 = flat_ssa;
flat_ssa1.surface_ssa_kernel_mode = 'ssa1_geometry';
flat_ssa1.surface_ssa_geometry_source_id = ['SSA.md; Thorsos & Broschat 1995 JASA, ', ...
    'Dirichlet SSA first-order / perturbation-limit geometry'];
fprintf('Running Hs=0 ssa1_geometry flat-degeneration case.\n');
channel_flat_ssa1 = CARPE3D_vertical(flat_ssa1);
meta_flat_ssa1 = channel_flat_ssa1.roughness_meta.ssa_stat_kernel_meta;
case_results.flat_ssa1 = local_compact_channel_result(channel_flat_ssa1);
case_results.flat_ssa1_diff_max_abs = max(abs(channel_flat_ssa1.H_f(:) - channel_flat_kirchhoff.H_f(:)));
checks = local_add_check(checks, 'Hs0_ssa1_kernel_mode', ...
    double(strcmp(meta_flat_ssa1.kernel_mode, 'ssa1_geometry')), 1, '==');
checks = local_add_check(checks, 'Hs0_ssa1_P_sca_zero', ...
    abs(meta_flat_ssa1.P_sca.sum), tol_roundoff, '<=');
checks = local_add_check(checks, 'Hs0_ssa1_flat_model_diff', ...
    case_results.flat_ssa1_diff_max_abs, tol_invariant, '<=');

direct_only = flat_ssa;
direct_only.enable_surface_reflection = false;
fprintf('Running ssa_stat_kernel direct-only interface case.\n');
channel_direct_only = CARPE3D_vertical(direct_only);
case_results.direct_only = local_compact_channel_result(channel_direct_only);
checks = local_add_check(checks, 'direct_only_reflect_zero', ...
    max(abs(channel_direct_only.H_reflect_f(:))), eps, '<=');
checks = local_add_check(checks, 'direct_only_invariant', ...
    local_invariant_error(channel_direct_only), tol_invariant, '<=');
checks = local_add_check(checks, 'direct_only_core_fields_present', ...
    double(local_has_core_fields(channel_direct_only)), 1, '==');

energy_case = params_base;
energy_case.surface_boundary_model = 'ssa_stat_kernel';
energy_case.sea_hs_target = 0.5;
energy_case.sea_seed = 12345;
fprintf('Running ssa_stat_kernel energy-limit case.\n');
channel_energy = CARPE3D_vertical(energy_case);
meta_energy = channel_energy.roughness_meta.ssa_stat_kernel_meta;
case_results.energy_limit = local_compact_channel_result(channel_energy);
case_results.energy_limit.ssa_meta = local_compact_ssa_meta(meta_energy);
energy_margin = meta_energy.E_inc + 1e-12*max(meta_energy.E_inc, 1) - ...
    (meta_energy.E_coh + meta_energy.E_sca);
ref_energy_margin = meta_energy.E_inc + 1e-12*max(meta_energy.E_inc, 1) - meta_energy.E_ref;
checks = local_add_check(checks, 'energy_coh_plus_sca_limited', energy_margin, 0, '>=');
checks = local_add_check(checks, 'energy_ref_limited', ref_energy_margin, 0, '>=');
checks = local_add_check(checks, 'energy_metadata_complete', ...
    double(local_has_ssa_metadata(meta_energy)), 1, '==');
checks = local_add_check(checks, 'energy_metadata_core_stable', ...
    double(local_has_stable_metadata_values(meta_energy)), 1, '==');
checks = local_add_check(checks, 'energy_conservation_error_small', ...
    meta_energy.energy_conservation_error, tol_roundoff, '<=');
checks = local_add_check(checks, 'propagating_bin_fraction_valid', ...
    double(meta_energy.propagating_bin_fraction >= 0 && meta_energy.propagating_bin_fraction <= 1), 1, '==');

fprintf('Running deterministic same-seed repeat.\n');
channel_repeat = CARPE3D_vertical(energy_case);
case_results.same_seed_repeat_max_abs_H_f = max(abs(channel_energy.H_f(:) - channel_repeat.H_f(:)));
case_results.same_seed_repeat_max_abs_H_reflect_f = max(abs(channel_energy.H_reflect_f(:) - channel_repeat.H_reflect_f(:)));
checks = local_add_check(checks, 'same_seed_H_f_repeat', ...
    case_results.same_seed_repeat_max_abs_H_f, tol_roundoff, '<=');
checks = local_add_check(checks, 'same_seed_H_reflect_repeat', ...
    case_results.same_seed_repeat_max_abs_H_reflect_f, tol_roundoff, '<=');

seed_changed = energy_case;
seed_changed.sea_seed = 12346;
fprintf('Running changed-seed sensitivity case.\n');
channel_seed_changed = CARPE3D_vertical(seed_changed);
case_results.changed_seed = local_compact_channel_result(channel_seed_changed);
case_results.changed_seed_direct_drift = max(abs(channel_energy.H_direct_f(:) - channel_seed_changed.H_direct_f(:)));
case_results.changed_seed_reflect_change = max(abs(channel_energy.H_reflect_f(:) - channel_seed_changed.H_reflect_f(:)));
checks = local_add_check(checks, 'changed_seed_direct_stable', ...
    case_results.changed_seed_direct_drift, tol_roundoff, '<=');
checks = local_add_check(checks, 'changed_seed_reflect_changes', ...
    case_results.changed_seed_reflect_change, 0, '>');

metadata_only = energy_case;
metadata_only.surface_ssa_random_scatter = false;
fprintf('Running random_scatter=false metadata-only case.\n');
channel_metadata_only = CARPE3D_vertical(metadata_only);
meta_metadata_only = channel_metadata_only.roughness_meta.ssa_stat_kernel_meta;
case_results.metadata_only = local_compact_channel_result(channel_metadata_only);
case_results.metadata_only.ssa_meta = local_compact_ssa_meta(meta_metadata_only);
checks = local_add_check(checks, 'metadata_only_P_sca_finite', ...
    double(isfinite(meta_metadata_only.P_sca.sum)), 1, '==');
checks = local_add_check(checks, 'metadata_only_E_sca_zero', ...
    abs(meta_metadata_only.E_sca), tol_roundoff, '<=');
checks = local_add_check(checks, 'metadata_only_random_disabled', ...
    double(~meta_metadata_only.random_scatter_enabled), 1, '==');

scale_zero = energy_case;
scale_zero.surface_ssa_scatter_scale = 0;
fprintf('Running scatter_scale=0 degeneration case.\n');
channel_scale_zero = CARPE3D_vertical(scale_zero);
meta_scale_zero = channel_scale_zero.roughness_meta.ssa_stat_kernel_meta;
case_results.scale_zero = local_compact_channel_result(channel_scale_zero);
case_results.scale_zero.ssa_meta = local_compact_ssa_meta(meta_scale_zero);
checks = local_add_check(checks, 'scale0_P_sca_raw_zero', ...
    abs(meta_scale_zero.P_sca_raw.sum), tol_roundoff, '<=');
checks = local_add_check(checks, 'scale0_E_sca_raw_zero', ...
    abs(meta_scale_zero.E_sca_raw), tol_roundoff, '<=');
checks = local_add_check(checks, 'scale0_E_sca_zero', ...
    abs(meta_scale_zero.E_sca), tol_roundoff, '<=');

scale_values_check = [0, 0.25, 1, 4];
[pm_scale_raw, ssa1_scale_raw] = local_run_scale_monotonic_check( ...
    energy_case, scale_values_check, flat_ssa1.surface_ssa_geometry_source_id);
case_results.scale_monotonic = struct( ...
    'scale_values', scale_values_check, ...
    'pm_convolution_E_sca_raw', pm_scale_raw, ...
    'ssa1_geometry_E_sca_raw', ssa1_scale_raw);
checks = local_add_check(checks, 'pm_scale_E_sca_raw_nondecreasing', ...
    double(all(diff(pm_scale_raw) >= -1e-9*max(max(pm_scale_raw), 1))), 1, '==');
checks = local_add_check(checks, 'ssa1_scale_E_sca_raw_nondecreasing', ...
    double(all(diff(ssa1_scale_raw) >= -1e-9*max(max(ssa1_scale_raw), 1))), 1, '==');

ssa1_energy = energy_case;
ssa1_energy.surface_ssa_kernel_mode = 'ssa1_geometry';
ssa1_energy.surface_ssa_geometry_source_id = flat_ssa1.surface_ssa_geometry_source_id;
fprintf('Running ssa1_geometry energy-limit case.\n');
channel_ssa1_energy = CARPE3D_vertical(ssa1_energy);
meta_ssa1_energy = channel_ssa1_energy.roughness_meta.ssa_stat_kernel_meta;
case_results.ssa1_energy = local_compact_channel_result(channel_ssa1_energy);
case_results.ssa1_energy.ssa_meta = local_compact_ssa_meta(meta_ssa1_energy);
ssa1_energy_margin = meta_ssa1_energy.E_inc + 1e-12*max(meta_ssa1_energy.E_inc, 1) - ...
    (meta_ssa1_energy.E_coh + meta_ssa1_energy.E_sca);
checks = local_add_check(checks, 'ssa1_energy_limited', ssa1_energy_margin, 0, '>=');
checks = local_add_check(checks, 'ssa1_energy_conservation_error_small', ...
    meta_ssa1_energy.energy_conservation_error, tol_roundoff, '<=');
checks = local_add_check(checks, 'ssa1_formula_source_recorded', ...
    double(contains(meta_ssa1_energy.kernel_detail.formula_source, 'SSA.md')), 1, '==');
checks = local_add_check(checks, 'ssa1_metadata_core_stable', ...
    double(local_has_stable_metadata_values(meta_ssa1_energy)), 1, '==');

bad_boundary = ssa1_energy;
bad_boundary.surface_reflect_coeff = -0.8;
fprintf('Running ssa1_geometry non-Dirichlet rejection case.\n');
[bad_boundary_error_id, bad_boundary_error_message] = local_expect_error(@() CARPE3D_vertical(bad_boundary));
case_results.bad_boundary_error_id = string(bad_boundary_error_id);
case_results.bad_boundary_error_message = string(bad_boundary_error_message);
checks = local_add_check(checks, 'ssa1_requires_dirichlet_R0', ...
    double(strcmp(bad_boundary_error_id, 'pm_surface_kirchhoff_module:SsaDirichletReflectCoeffRequired')), 1, '==');

zero_padded_case = energy_case;
zero_padded_case.surface_ssa_conv_padding = 'zero_padded';
fprintf('Running zero_padded convolution rejection case.\n');
[zero_padded_error_id, zero_padded_error_message] = local_expect_error(@() CARPE3D_vertical(zero_padded_case));
case_results.zero_padded_error_id = string(zero_padded_error_id);
case_results.zero_padded_error_message = string(zero_padded_error_message);
checks = local_add_check(checks, 'zero_padded_rejected_explicitly', ...
    double(strcmp(zero_padded_error_id, 'pm_surface_kirchhoff_module:SsaConvPaddingNotImplemented')), 1, '==');

dense_fft = params_base;
dense_fft.nx = 32;
dense_fft.ny = 32;
dense_fft.xw = 8;
dense_fft.yw = 8;
dense_fft.surface_boundary_model = 'ssa_stat_kernel';
dense_fft.sea_hs_target = 0.05;
dense_fft.sea_seed = 12345;
dense_fft.surface_ssa_random_scatter = false;
dense_fft.surface_ssa_kernel_mode = 'ssa1_geometry';
dense_fft.surface_ssa_geometry_source_id = flat_ssa1.surface_ssa_geometry_source_id;
dense_dense = dense_fft;
dense_dense.surface_ssa_kernel_mode = 'ssa1_debug_dense';
fprintf('Running small-grid ssa1_geometry FFT case.\n');
channel_dense_fft = CARPE3D_vertical(dense_fft);
fprintf('Running small-grid ssa1_debug_dense case.\n');
channel_dense_dense = CARPE3D_vertical(dense_dense);
meta_dense_fft = channel_dense_fft.roughness_meta.ssa_stat_kernel_meta;
meta_dense_dense = channel_dense_dense.roughness_meta.ssa_stat_kernel_meta;
case_results.ssa1_dense_compare = struct();
case_results.ssa1_dense_compare.P_sca_sum_diff = abs(meta_dense_fft.P_sca.sum - meta_dense_dense.P_sca.sum);
case_results.ssa1_dense_compare.P_sca_raw_sum_diff = abs(meta_dense_fft.P_sca_raw.sum - meta_dense_dense.P_sca_raw.sum);
dense_scale = max(abs(meta_dense_fft.P_sca_raw.sum), 1);
checks = local_add_check(checks, 'ssa1_dense_fft_P_sca_raw_sum_match', ...
    case_results.ssa1_dense_compare.P_sca_raw_sum_diff / dense_scale, 1e-10, '<=');

summary_table = struct2table(checks);
validation_meta = struct();
validation_meta.script = mfilename;
validation_meta.result_file = result_file;
validation_meta.created_at = char(datetime('now'));
validation_meta.params_base = params_base;
validation_meta.tol_invariant = tol_invariant;
validation_meta.tol_roundoff = tol_roundoff;
validation_meta.all_passed = all(summary_table.passed);
validation_meta.notes = ['SSA statistical kernel validation. pm_convolution remains ', ...
    'the engineering baseline; ssa1_geometry implements the Dirichlet first-order ', ...
    'geometry factor from SSA.md but still does not cover NLSSA or calibration.'];

disp(summary_table(:, {'check_name', 'value', 'tolerance', 'comparison', 'passed'}))
if ~validation_meta.all_passed
    error('validate_ssa_stat_kernel_vertical:Failed', ...
        'One or more ssa_stat_kernel validation checks failed.');
end

save(result_file, 'summary_table', 'case_results', 'validation_meta');
fprintf('Saved %s\n', result_file);

function paramsV = local_base_params()
paramsV = struct();
paramsV.f0 = 6000;
paramsV.enable_wideband = false;
paramsV.c0 = 1500;
paramsV.z_max = 100;
paramsV.stepz_lamb = 0.5;
paramsV.xw = 50;
paramsV.yw = 50;
paramsV.nx = 128;
paramsV.ny = 128;
paramsV.x_tx = 0;
paramsV.y_tx = 0;
paramsV.z_tx = 100;
paramsV.x_rx = 0;
paramsV.y_rx = 0;
paramsV.z_rx = 3;
paramsV.nout = 2;
paramsV.sigma_src_m = 0.4;
paramsV.sponge_ratio = 0.12;
paramsV.alpha_max_np_per_m = 0.15;
paramsV.env_mode = 'uniform';
paramsV.enforce_1_over_R = false;
paramsV.show_figures = false;
paramsV.save_mode = 'rx_only';
paramsV.use_gpu = false;
paramsV.enable_surface_reflection = true;
paramsV.sea_wind_speed = 5.0;
paramsV.sea_hs_target = 0.5;
paramsV.sea_seed = 12345;
paramsV.surface_reflect_coeff = -1;
paramsV.surface_boundary_model = 'kirchhoff_spatial';
paramsV.surface_boundary_coupling_diagnostics = false;
paramsV.surface_boundary_redistribution_diagnostics = false;
paramsV.surface_ssa_random_scatter = true;
paramsV.surface_ssa_scatter_scale = 1.0;
paramsV.surface_ssa_seed_offset = 100000;
paramsV.surface_ssa_kernel_mode = 'pm_convolution';
paramsV.surface_ssa_geometry_source_id = '';
paramsV.surface_ssa_kz_branch = 'downward_positive_real';
paramsV.surface_ssa_conv_padding = 'periodic';
end

function err = local_invariant_error(channel)
err = max(abs(channel.H_f(:) - channel.H_direct_f(:) - channel.H_reflect_f(:)));
end

function tf = local_has_core_fields(channel)
names = {'h_direct', 'h_reflect', 'h_total', 'H_f', 'H_direct_f', ...
    'H_reflect_f', 'f_axis', 'idx_f_ref'};
tf = true;
for ii = 1:numel(names)
    tf = tf && isfield(channel, names{ii});
end
end

function tf = local_has_ssa_metadata(meta)
names = {'sigma_eta_m', 'Hs_target_m', 'R_coh', 'P_sca', 'P_sca_raw', 'E_inc', ...
    'E_coh', 'E_sca_raw', 'E_sca_limited', 'E_sca', 'E_ref', ...
    'energy_scale_applied', 'energy_limit_applied', ...
    'energy_conservation_error', 'propagating_bin_fraction', ...
    'kernel_mode', 'geometry_source_id', 'formula_source', ...
    'boundary_condition', 'evanescent_included', 'conv_padding', ...
    'kernel_detail', 'seed_ssa'};
tf = true;
for ii = 1:numel(names)
    tf = tf && isfield(meta, names{ii});
end
tf = tf && isfield(meta.P_sca, 'sum');
end

function tf = local_has_stable_metadata_values(meta)
tf = isfield(meta, 'kernel_mode') && ~isempty(meta.kernel_mode) && ...
    isfield(meta, 'conv_padding') && ~isempty(meta.conv_padding) && ...
    isfield(meta, 'formula_source') && ~isempty(meta.formula_source) && ...
    isfield(meta, 'boundary_condition') && ~isempty(meta.boundary_condition) && ...
    isfield(meta, 'seed_ssa') && isfinite(meta.seed_ssa) && ...
    isfield(meta, 'P_sca') && isfield(meta.P_sca, 'sum') && ...
    isfield(meta, 'P_sca_raw') && isfield(meta.P_sca_raw, 'sum') && ...
    isfield(meta, 'energy_conservation_error') && isfinite(meta.energy_conservation_error);
end

function out = local_compact_channel_result(channel)
out = struct();
out.h_direct = channel.h_direct;
out.h_reflect = channel.h_reflect;
out.h_total = channel.h_total;
out.H_f = channel.H_f(:);
out.H_direct_f = channel.H_direct_f(:);
out.H_reflect_f = channel.H_reflect_f(:);
out.f_axis = channel.f_axis(:);
out.idx_f_ref = channel.idx_f_ref;
out.invariant_error = local_invariant_error(channel);
out.max_abs_H_reflect_f = max(abs(channel.H_reflect_f(:)));
out.roughness_enabled = channel.roughness_meta.enabled;
out.surface_realization_generated = local_get_field_or_nan( ...
    channel.roughness_meta, 'surface_realization_generated');
end

function out = local_compact_ssa_meta(meta)
out = struct();
out.enabled = meta.enabled;
out.random_scatter_enabled = meta.random_scatter_enabled;
out.kernel_mode = meta.kernel_mode;
out.geometry_source_id = meta.geometry_source_id;
out.kz_branch = meta.kz_branch;
out.conv_padding = meta.conv_padding;
out.formula_source = meta.formula_source;
out.boundary_condition = meta.boundary_condition;
out.evanescent_included = meta.evanescent_included;
out.sigma_eta_m = meta.sigma_eta_m;
out.Hs_target_m = meta.Hs_target_m;
out.R_coh = meta.R_coh;
out.P_sca = meta.P_sca;
out.P_sca_raw = meta.P_sca_raw;
out.E_inc = meta.E_inc;
out.E_coh = meta.E_coh;
out.E_sca_raw = meta.E_sca_raw;
out.E_sca_limited = meta.E_sca_limited;
out.E_sca = meta.E_sca;
out.E_ref = meta.E_ref;
out.energy_scale_applied = meta.energy_scale_applied;
out.energy_limit_applied = meta.energy_limit_applied;
out.energy_conservation_error = meta.energy_conservation_error;
out.propagating_bin_fraction = meta.propagating_bin_fraction;
out.seed_ssa = meta.seed_ssa;
out.seed_offset = meta.seed_offset;
out.limitations = meta.limitations;
end

function [err_id, err_message] = local_expect_error(fn)
try
    fn();
    err_id = '';
    err_message = '';
catch ME
    err_id = ME.identifier;
    err_message = ME.message;
end
end

function [pm_raw, ssa1_raw] = local_run_scale_monotonic_check(base_case, scale_values, source_id)
pm_raw = NaN(size(scale_values));
ssa1_raw = NaN(size(scale_values));
for ii = 1:numel(scale_values)
    params_pm = base_case;
    params_pm.surface_ssa_scatter_scale = scale_values(ii);
    params_pm.surface_ssa_random_scatter = false;
    params_pm.surface_ssa_kernel_mode = 'pm_convolution';
    channel_pm = CARPE3D_vertical(params_pm);
    pm_raw(ii) = channel_pm.roughness_meta.ssa_stat_kernel_meta.E_sca_raw;

    params_ssa1 = params_pm;
    params_ssa1.surface_ssa_kernel_mode = 'ssa1_geometry';
    params_ssa1.surface_ssa_geometry_source_id = source_id;
    channel_ssa1 = CARPE3D_vertical(params_ssa1);
    ssa1_raw(ii) = channel_ssa1.roughness_meta.ssa_stat_kernel_meta.E_sca_raw;
end
end

function checks = local_add_check(checks, check_name, value, tolerance, comparison)
row = struct();
row.check_name = string(check_name);
row.value = value;
row.tolerance = tolerance;
row.comparison = string(comparison);
switch comparison
    case '<='
        row.passed = value <= tolerance;
    case '>='
        row.passed = value >= tolerance;
    case '>'
        row.passed = value > tolerance;
    case '=='
        row.passed = value == tolerance;
    otherwise
        error('Unsupported comparison %s.', comparison);
end
if isempty(checks)
    checks = row;
else
    checks(end + 1) = row; %#ok<AGROW>
end
end

function out = local_get_field_or_nan(s, field_name)
if isfield(s, field_name)
    out = s.(field_name);
else
    out = NaN;
end
end
