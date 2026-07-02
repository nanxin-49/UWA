run(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'bootstrap_project.m'));
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
channel_flat_ssa = vertical_channel_model(flat_ssa);
fprintf('Running Hs=0 kirchhoff_spatial flat-reference case.\n');
channel_flat_kirchhoff = vertical_channel_model(flat_kirchhoff);
meta_flat = channel_flat_ssa.roughness_meta.ssa_stat_kernel_meta;
case_results.flat_ssa = local_compact_channel_result(channel_flat_ssa);
case_results.flat_kirchhoff = local_compact_channel_result(channel_flat_kirchhoff);
case_results.flat_diff_max_abs = max(abs(channel_flat_ssa.H_f(:) - channel_flat_kirchhoff.H_f(:)));

checks = local_add_check(checks, 'Hs0_invariant', ...
    local_invariant_error(channel_flat_ssa), tol_invariant, '<=');
checks = local_add_check(checks, 'Hs0_sigma_eta_zero', ...
    abs(meta_flat.sigma_eta_m), tol_roundoff, '<=');
checks = local_add_check(checks, 'Hs0_W_eta_variance_zero', ...
    abs(meta_flat.W_eta_variance_discrete), tol_roundoff, '<=');
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

kirchhoff_phase = params_base;
kirchhoff_phase.surface_boundary_model = 'kirchhoff_spatial';
kirchhoff_phase.sea_hs_target = 0.2;
fprintf('Running kirchhoff_spatial phase-screen formula case.\n');
channel_kirchhoff_phase = vertical_channel_model(kirchhoff_phase);
k0_kirchhoff_phase = 2*pi*kirchhoff_phase.f0/kirchhoff_phase.c0;
phase_std_expected = k0_kirchhoff_phase * ...
    mean(channel_kirchhoff_phase.roughness_meta.phase_factor_stats.mean) * ...
    std(channel_kirchhoff_phase.surface_elevation(:));
phase_std_observed = std(channel_kirchhoff_phase.delta_phi(:));
case_results.kirchhoff_phase_screen = local_compact_channel_result(channel_kirchhoff_phase);
case_results.kirchhoff_phase_screen.delta_phi_std_observed = phase_std_observed;
case_results.kirchhoff_phase_screen.delta_phi_std_expected = phase_std_expected;
checks = local_add_check(checks, 'kirchhoff_phase_screen_formula_recorded', ...
    double(contains(channel_kirchhoff_phase.roughness_meta.kirchhoff_phase_screen_formula, ...
    'delta_phi=k0*(cos_i+cos_r)*eta')), 1, '==');
checks = local_add_check(checks, 'kirchhoff_normal_delta_phi_std_matches_2k_eta', ...
    abs(phase_std_observed - phase_std_expected) / max(abs(phase_std_expected), eps), ...
    1e-12, '<=');

flat_ssa1 = flat_ssa;
flat_ssa1.surface_ssa_kernel_mode = 'ssa1_geometry';
flat_ssa1.surface_ssa_geometry_source_id = ['SSA.md; Thorsos & Broschat 1995 JASA, ', ...
    'Dirichlet SSA first-order / perturbation-limit geometry'];
fprintf('Running Hs=0 ssa1_geometry flat-degeneration case.\n');
channel_flat_ssa1 = vertical_channel_model(flat_ssa1);
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
channel_direct_only = vertical_channel_model(direct_only);
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
channel_energy = vertical_channel_model(energy_case);
meta_energy = channel_energy.roughness_meta.ssa_stat_kernel_meta;
case_results.energy_limit = local_compact_channel_result(channel_energy);
case_results.energy_limit.ssa_meta = local_compact_ssa_meta(meta_energy);
k0_ref = 2*pi*energy_case.f0/energy_case.c0;
expected_R_coh_normal = energy_case.surface_reflect_coeff * ...
    exp(-2 * k0_ref^2 * meta_energy.sigma_eta_m^2);
old_wrong_R_coh_normal = energy_case.surface_reflect_coeff * ...
    exp(-8 * k0_ref^2 * meta_energy.sigma_eta_m^2);
case_results.normal_coherent_formula = struct( ...
    'k0_rad_per_m', k0_ref, ...
    'sigma_eta_m', meta_energy.sigma_eta_m, ...
    'R_coh', meta_energy.R_coh, ...
    'expected_R_coh_normal', expected_R_coh_normal, ...
    'old_wrong_R_coh_normal', old_wrong_R_coh_normal, ...
    'coherent_exponent', meta_energy.coherent_exponent);
energy_margin = meta_energy.E_inc + 1e-12*max(meta_energy.E_inc, 1) - ...
    (meta_energy.E_coh + meta_energy.E_sca);
ref_energy_margin = meta_energy.E_inc + 1e-12*max(meta_energy.E_inc, 1) - meta_energy.E_ref;
checks = local_add_check(checks, 'energy_coh_plus_sca_limited', energy_margin, 0, '>=');
checks = local_add_check(checks, 'energy_ref_limited', ref_energy_margin, 0, '>=');
checks = local_add_check(checks, 'energy_metadata_complete', ...
    double(local_has_ssa_metadata(meta_energy)), 1, '==');
checks = local_add_check(checks, 'energy_metadata_core_stable', ...
    double(local_has_stable_metadata_values(meta_energy)), 1, '==');
checks = local_add_check(checks, 'pm_W_eta_variance_rel_error_small', ...
    meta_energy.W_eta_variance_rel_error, 1e-10, '<=');
checks = local_add_check(checks, 'pm_convolution_phase_scale_corrected', ...
    abs(meta_energy.phase_scale_corrected_rad_per_m - k0_ref * meta_energy.phase_factor_eff), ...
    tol_roundoff, '<=');
checks = local_add_check(checks, 'pm_convolution_phase_scale_legacy_recorded', ...
    abs(meta_energy.phase_scale_legacy_rad_per_m - 2*k0_ref * meta_energy.phase_factor_eff), ...
    tol_roundoff, '<=');
checks = local_add_check(checks, 'pm_convolution_kernel_formula_corrected', ...
    double(contains(meta_energy.kernel_formula, '(k0*phase_factor_eff)^2')), 1, '==');
checks = local_add_check(checks, 'normal_R_coh_matches_broschat_formula', ...
    abs(meta_energy.R_coh - expected_R_coh_normal), tol_roundoff, '<=');
checks = local_add_check(checks, 'normal_R_coh_not_old_exp8_formula', ...
    abs(meta_energy.R_coh - old_wrong_R_coh_normal), 1e-12, '>');
checks = local_add_check(checks, 'coherent_gamma_sum_normal_is_2k0', ...
    abs(meta_energy.coherent_gamma_sum_eff_rad_per_m - 2*k0_ref), tol_roundoff, '<=');
checks = local_add_check(checks, 'energy_conservation_error_small', ...
    meta_energy.energy_conservation_error, tol_roundoff, '<=');
checks = local_add_check(checks, 'propagating_bin_fraction_valid', ...
    double(meta_energy.propagating_bin_fraction >= 0 && meta_energy.propagating_bin_fraction <= 1), 1, '==');

fprintf('Running deterministic same-seed repeat.\n');
channel_repeat = vertical_channel_model(energy_case);
case_results.same_seed_repeat_max_abs_H_f = max(abs(channel_energy.H_f(:) - channel_repeat.H_f(:)));
case_results.same_seed_repeat_max_abs_H_reflect_f = max(abs(channel_energy.H_reflect_f(:) - channel_repeat.H_reflect_f(:)));
checks = local_add_check(checks, 'same_seed_H_f_repeat', ...
    case_results.same_seed_repeat_max_abs_H_f, tol_roundoff, '<=');
checks = local_add_check(checks, 'same_seed_H_reflect_repeat', ...
    case_results.same_seed_repeat_max_abs_H_reflect_f, tol_roundoff, '<=');

seed_changed = energy_case;
seed_changed.sea_seed = 12346;
fprintf('Running changed-seed sensitivity case.\n');
channel_seed_changed = vertical_channel_model(seed_changed);
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
channel_metadata_only = vertical_channel_model(metadata_only);
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
channel_scale_zero = vertical_channel_model(scale_zero);
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

hs_trend_values = [0, 0.05, 0.2, 0.5];
hs_trend = local_run_hs_coherent_trend(energy_case, hs_trend_values);
case_results.hs_coherent_trend = hs_trend;
checks = local_add_check(checks, 'Hs_R_coh_abs_monotone_decreasing', ...
    double(all(diff(hs_trend.R_coh_abs) <= 1e-12)), 1, '==');
checks = local_add_check(checks, 'Hs_E_sca_limit_nondecreasing', ...
    double(all(diff(hs_trend.E_sca_limit) >= -1e-9*max(max(hs_trend.E_sca_limit), 1))), 1, '==');
checks = local_add_check(checks, 'Hs_W_eta_variance_rel_error_small', ...
    max(hs_trend.W_eta_variance_rel_error), 1e-10, '<=');

ssa1_energy = energy_case;
ssa1_energy.surface_ssa_kernel_mode = 'ssa1_geometry';
ssa1_energy.surface_ssa_geometry_source_id = flat_ssa1.surface_ssa_geometry_source_id;
fprintf('Running ssa1_geometry energy-limit case.\n');
channel_ssa1_energy = vertical_channel_model(ssa1_energy);
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
checks = local_add_check(checks, 'ssa1_G_formula_recorded', ...
    double(contains(meta_ssa1_energy.kernel_detail.G_SSA1_formula, '4*gamma')), 1, '==');
checks = local_add_check(checks, 'ssa1_boundary_condition_dirichlet', ...
    double(strcmp(meta_ssa1_energy.boundary_condition, 'pressure-release / Dirichlet')), 1, '==');
checks = local_add_check(checks, 'ssa1_metadata_core_stable', ...
    double(local_has_stable_metadata_values(meta_ssa1_energy)), 1, '==');

bad_boundary = ssa1_energy;
bad_boundary.surface_reflect_coeff = -0.8;
fprintf('Running ssa1_geometry non-Dirichlet rejection case.\n');
[bad_boundary_error_id, bad_boundary_error_message] = local_expect_error(@() vertical_channel_model(bad_boundary));
case_results.bad_boundary_error_id = string(bad_boundary_error_id);
case_results.bad_boundary_error_message = string(bad_boundary_error_message);
checks = local_add_check(checks, 'ssa1_requires_dirichlet_R0', ...
    double(strcmp(bad_boundary_error_id, 'pm_surface_boundary_model:SsaDirichletReflectCoeffRequired')), 1, '==');

zero_padded_pm = metadata_only;
zero_padded_pm.surface_ssa_conv_padding = 'zero_padded';
fprintf('Running zero_padded pm_convolution case.\n');
channel_zero_padded_pm = vertical_channel_model(zero_padded_pm);
meta_zero_padded_pm = channel_zero_padded_pm.roughness_meta.ssa_stat_kernel_meta;
case_results.zero_padded_pm = local_compact_channel_result(channel_zero_padded_pm);
case_results.zero_padded_pm.ssa_meta = local_compact_ssa_meta(meta_zero_padded_pm);
case_results.zero_padded_pm.periodic_raw_sum = meta_metadata_only.E_sca_raw;
case_results.zero_padded_pm.zero_padded_raw_sum = meta_zero_padded_pm.E_sca_raw;
case_results.zero_padded_pm.raw_sum_relative_delta = ...
    abs(meta_zero_padded_pm.E_sca_raw - meta_metadata_only.E_sca_raw) / max(abs(meta_metadata_only.E_sca_raw), 1);
checks = local_add_check(checks, 'zero_padded_pm_runs', ...
    double(strcmp(meta_zero_padded_pm.conv_padding, 'zero_padded')), 1, '==');
checks = local_add_check(checks, 'zero_padded_pm_operator_recorded', ...
    double(strcmp(meta_zero_padded_pm.conv_operator, 'linear_fft_zero_padded')), 1, '==');
checks = local_add_check(checks, 'zero_padded_pm_energy_conservation_small', ...
    meta_zero_padded_pm.energy_conservation_error, tol_roundoff, '<=');
checks = local_add_check(checks, 'zero_padded_pm_raw_difference_finite', ...
    double(isfinite(case_results.zero_padded_pm.raw_sum_relative_delta)), 1, '==');

zero_padded_ssa1 = ssa1_energy;
zero_padded_ssa1.surface_ssa_conv_padding = 'zero_padded';
zero_padded_ssa1.surface_ssa_random_scatter = false;
fprintf('Running zero_padded ssa1_geometry case.\n');
channel_zero_padded_ssa1 = vertical_channel_model(zero_padded_ssa1);
meta_zero_padded_ssa1 = channel_zero_padded_ssa1.roughness_meta.ssa_stat_kernel_meta;
case_results.zero_padded_ssa1 = local_compact_channel_result(channel_zero_padded_ssa1);
case_results.zero_padded_ssa1.ssa_meta = local_compact_ssa_meta(meta_zero_padded_ssa1);
case_results.zero_padded_ssa1.periodic_raw_sum = meta_ssa1_energy.E_sca_raw;
case_results.zero_padded_ssa1.zero_padded_raw_sum = meta_zero_padded_ssa1.E_sca_raw;
case_results.zero_padded_ssa1.raw_sum_relative_delta = ...
    abs(meta_zero_padded_ssa1.E_sca_raw - meta_ssa1_energy.E_sca_raw) / max(abs(meta_ssa1_energy.E_sca_raw), 1);
checks = local_add_check(checks, 'zero_padded_ssa1_runs', ...
    double(strcmp(meta_zero_padded_ssa1.conv_padding, 'zero_padded')), 1, '==');
checks = local_add_check(checks, 'zero_padded_ssa1_operator_recorded', ...
    double(strcmp(meta_zero_padded_ssa1.conv_operator, 'linear_fft_zero_padded')), 1, '==');
checks = local_add_check(checks, 'zero_padded_ssa1_energy_conservation_small', ...
    meta_zero_padded_ssa1.energy_conservation_error, tol_roundoff, '<=');
checks = local_add_check(checks, 'zero_padded_ssa1_raw_difference_finite', ...
    double(isfinite(case_results.zero_padded_ssa1.raw_sum_relative_delta)), 1, '==');

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
channel_dense_fft = vertical_channel_model(dense_fft);
fprintf('Running small-grid ssa1_debug_dense case.\n');
channel_dense_dense = vertical_channel_model(dense_dense);
meta_dense_fft = channel_dense_fft.roughness_meta.ssa_stat_kernel_meta;
meta_dense_dense = channel_dense_dense.roughness_meta.ssa_stat_kernel_meta;
case_results.ssa1_dense_compare = struct();
case_results.ssa1_dense_compare.P_sca_sum_diff = abs(meta_dense_fft.P_sca.sum - meta_dense_dense.P_sca.sum);
case_results.ssa1_dense_compare.P_sca_raw_sum_diff = abs(meta_dense_fft.P_sca_raw.sum - meta_dense_dense.P_sca_raw.sum);
dense_scale = max(abs(meta_dense_fft.P_sca_raw.sum), 1);
checks = local_add_check(checks, 'ssa1_dense_fft_P_sca_raw_sum_match', ...
    case_results.ssa1_dense_compare.P_sca_raw_sum_diff / dense_scale, 1e-10, '<=');

wideband_zp = params_base;
wideband_zp.enable_wideband = true;
wideband_zp.f0 = 4000;
wideband_zp.f_band_hz = [4000, 8000];
wideband_zp.Nf_min = 4;
wideband_zp.Nf_max = 4;
wideband_zp.f_ref_hz = 6000;
wideband_zp.surface_boundary_model = 'ssa_stat_kernel';
wideband_zp.surface_ssa_kernel_mode = 'ssa1_geometry';
wideband_zp.surface_ssa_conv_padding = 'zero_padded';
wideband_zp.surface_ssa_random_scatter = false;
wideband_zp.surface_ssa_geometry_source_id = flat_ssa1.surface_ssa_geometry_source_id;
fprintf('Running reduced wideband zero_padded ssa1_geometry invariant case.\n');
channel_wideband_zp = vertical_channel_model(wideband_zp);
case_results.wideband_zero_padded_ssa1 = local_compact_channel_result(channel_wideband_zp);
case_results.wideband_zero_padded_ssa1.abs_H_f = abs(channel_wideband_zp.H_f(:));
case_results.wideband_zero_padded_ssa1.abs_H_reflect_f = abs(channel_wideband_zp.H_reflect_f(:));
checks = local_add_check(checks, 'wideband_zero_padded_ssa1_invariant', ...
    local_invariant_error(channel_wideband_zp), tol_invariant, '<=');
checks = local_add_check(checks, 'wideband_zero_padded_ssa1_frequency_count', ...
    numel(channel_wideband_zp.f_axis), 4, '==');
checks = local_add_check(checks, 'wideband_zero_padded_ssa1_ref_consistency', ...
    abs(channel_wideband_zp.h_total - channel_wideband_zp.H_f(channel_wideband_zp.idx_f_ref)), ...
    tol_roundoff, '<=');

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
paramsV.surface_ssa_coherent_order = 'ssa1';
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
    'R_coh_raw', 'R_coh_ssa1', 'R_coh_ssa2', 'abs_R_coh_ssa1', ...
    'abs_R_coh_ssa2', 'coherent_loss_ssa1_db', 'coherent_loss_ssa2_db', ...
    'coherent_loss_delta_db', 'delta_R_abs', 'delta_R_rel', ...
    'coherent_order', 'frequency_correlation_mode', 'frequency_correlation_rho', ...
    'coherent_exponent', 'coherent_gamma_sum_eff_rad_per_m', ...
    'coherent_vertical_factor_eff', 'coherent_reflection_formula', ...
    'phase_factor_eff', 'phase_factor_eff_legacy_note', ...
    'ssa2_correction_integral', 'ssa2_correction_integral_real', ...
    'ssa2_correction_integral_imag', 'ssa2_gamma_i_eff_rad_per_m', ...
    'ssa2_K_i_eff_rad_per_m', 'ssa2_sqrt_branch', 'ssa2_formula', ...
    'ssa2_formula_source', 'ssa2_limitations', 'ssa2_evanescent_handling', ...
    'ssa2_propagating_integral', 'ssa2_evanescent_integral', ...
    'ssa2_evanescent_fraction_abs', 'random_spectrum_seed_base', ...
    'random_spectrum_seed_innovation', 'random_spectrum_frequency_index', ...
    'random_spectrum_generation_rule', ...
    'W_eta_variance_target', 'W_eta_variance_discrete', 'W_eta_variance_rel_error', ...
    'kernel_mode', 'geometry_source_id', 'formula_source', ...
    'boundary_condition', 'evanescent_included', 'conv_padding', ...
    'conv_operator', 'conv_padding_size', 'conv_crop_start_index', ...
    'conv_crop_end_index', 'conv_crop_rule', ...
    'phase_scale_corrected_rad_per_m', 'phase_scale_legacy_rad_per_m', ...
    'phase_scale_formula', 'phase_scale_legacy_formula', ...
    'incident_spectrum_stats', 'scatter_power_spectrum_stats', ...
    'reflected_spectrum_stats', ...
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
    isfield(meta, 'conv_operator') && ~isempty(meta.conv_operator) && ...
    isfield(meta, 'phase_scale_formula') && ~isempty(meta.phase_scale_formula) && ...
    isfield(meta, 'formula_source') && ~isempty(meta.formula_source) && ...
    isfield(meta, 'boundary_condition') && ~isempty(meta.boundary_condition) && ...
    isfield(meta, 'seed_ssa') && isfinite(meta.seed_ssa) && ...
    isfield(meta, 'P_sca') && isfield(meta.P_sca, 'sum') && ...
    isfield(meta, 'P_sca_raw') && isfield(meta.P_sca_raw, 'sum') && ...
    isfield(meta, 'coherent_reflection_formula') && ~isempty(meta.coherent_reflection_formula) && ...
    isfield(meta, 'W_eta_variance_rel_error') && isfinite(meta.W_eta_variance_rel_error) && ...
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
out.conv_operator = meta.conv_operator;
out.conv_padding_size = meta.conv_padding_size;
out.conv_crop_start_index = meta.conv_crop_start_index;
out.conv_crop_end_index = meta.conv_crop_end_index;
out.conv_crop_rule = meta.conv_crop_rule;
out.phase_scale_corrected_rad_per_m = meta.phase_scale_corrected_rad_per_m;
out.phase_scale_legacy_rad_per_m = meta.phase_scale_legacy_rad_per_m;
out.phase_scale_formula = meta.phase_scale_formula;
out.phase_scale_legacy_formula = meta.phase_scale_legacy_formula;
out.formula_source = meta.formula_source;
out.boundary_condition = meta.boundary_condition;
out.evanescent_included = meta.evanescent_included;
out.sigma_eta_m = meta.sigma_eta_m;
out.Hs_target_m = meta.Hs_target_m;
out.R_coh = meta.R_coh;
out.R_coh_raw = meta.R_coh_raw;
out.coherent_exponent = meta.coherent_exponent;
out.coherent_gamma_sum_eff_rad_per_m = meta.coherent_gamma_sum_eff_rad_per_m;
out.coherent_vertical_factor_eff = meta.coherent_vertical_factor_eff;
out.coherent_reflection_formula = meta.coherent_reflection_formula;
out.phase_factor_eff = meta.phase_factor_eff;
out.W_eta_variance_target = meta.W_eta_variance_target;
out.W_eta_variance_discrete = meta.W_eta_variance_discrete;
out.W_eta_variance_rel_error = meta.W_eta_variance_rel_error;
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
out.incident_rms_delta_k_rad_per_m = meta.incident_spectrum_stats.rms_delta_k_rad_per_m;
out.scatter_rms_delta_k_rad_per_m = meta.scatter_power_spectrum_stats.rms_delta_k_rad_per_m;
out.reflected_rms_delta_k_rad_per_m = meta.reflected_spectrum_stats.rms_delta_k_rad_per_m;
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
    channel_pm = vertical_channel_model(params_pm);
    pm_raw(ii) = channel_pm.roughness_meta.ssa_stat_kernel_meta.E_sca_raw;

    params_ssa1 = params_pm;
    params_ssa1.surface_ssa_kernel_mode = 'ssa1_geometry';
    params_ssa1.surface_ssa_geometry_source_id = source_id;
    channel_ssa1 = vertical_channel_model(params_ssa1);
    ssa1_raw(ii) = channel_ssa1.roughness_meta.ssa_stat_kernel_meta.E_sca_raw;
end
end

function trend = local_run_hs_coherent_trend(base_case, hs_values)
R_coh_abs = NaN(size(hs_values));
E_sca_limit = NaN(size(hs_values));
W_eta_rel_error = NaN(size(hs_values));
coherent_exponent = NaN(size(hs_values));
for ii = 1:numel(hs_values)
    params = base_case;
    params.sea_hs_target = hs_values(ii);
    params.surface_ssa_random_scatter = false;
    params.surface_ssa_kernel_mode = 'pm_convolution';
    channel = vertical_channel_model(params);
    meta = channel.roughness_meta.ssa_stat_kernel_meta;
    R_coh_abs(ii) = abs(meta.R_coh);
    E_sca_limit(ii) = meta.E_sca_limit;
    W_eta_rel_error(ii) = meta.W_eta_variance_rel_error;
    coherent_exponent(ii) = meta.coherent_exponent;
end
trend = struct( ...
    'Hs_values', hs_values, ...
    'R_coh_abs', R_coh_abs, ...
    'E_sca_limit', E_sca_limit, ...
    'W_eta_variance_rel_error', W_eta_rel_error, ...
    'coherent_exponent', coherent_exponent);
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

