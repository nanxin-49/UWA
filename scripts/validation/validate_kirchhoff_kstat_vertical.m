run(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'bootstrap_project.m'));
% Reduced-grid validation for the Kirchhoff statistical phase-screen branch.

clear
format compact

result_file = getenv('KSTAT_VALIDATE_RESULT_FILE');
if isempty(result_file)
    result_file = project_result_file('validation', ...
        'validate_kirchhoff_kstat_vertical_result.mat');
end

tol_invariant = 1e-10;
tol_roundoff = 1e-12;
tol_energy = 1e-8;

checks = struct([]);
case_results = struct();
params_base = local_base_params();

flat_kstat = params_base;
flat_kstat.surface_boundary_model = 'kirchhoff_kstat';
flat_kstat.sea_hs_target = 0;

flat_kirchhoff = params_base;
flat_kirchhoff.surface_boundary_model = 'kirchhoff_spatial';
flat_kirchhoff.sea_hs_target = 0;

fprintf('Running Hs=0 kirchhoff_kstat flat-degeneration case.\n');
channel_flat_kstat = vertical_channel_model(flat_kstat);
fprintf('Running Hs=0 kirchhoff_spatial flat-reference case.\n');
channel_flat_kirchhoff = vertical_channel_model(flat_kirchhoff);
meta_flat = channel_flat_kstat.roughness_meta.kirchhoff_kstat_meta;
case_results.flat_kstat = local_compact_channel_result(channel_flat_kstat);
case_results.flat_kirchhoff = local_compact_channel_result(channel_flat_kirchhoff);
case_results.flat_diff_max_abs = max(abs(channel_flat_kstat.H_f(:) - channel_flat_kirchhoff.H_f(:)));

checks = local_add_check(checks, 'Hs0_invariant', ...
    local_invariant_error(channel_flat_kstat), tol_invariant, '<=');
checks = local_add_check(checks, 'Hs0_sigma_eta2_zero', ...
    abs(meta_flat.sigma_eta2_m2), tol_roundoff, '<=');
checks = local_add_check(checks, 'Hs0_R_coh_equals_R0', ...
    abs(meta_flat.R_coh - flat_kstat.surface_reflect_coeff), tol_roundoff, '<=');
checks = local_add_check(checks, 'Hs0_S_deltaG_zero', ...
    abs(meta_flat.S_deltaG.sum), tol_roundoff, '<=');
checks = local_add_check(checks, 'Hs0_P_sca_zero', ...
    abs(meta_flat.P_sca.sum), tol_roundoff, '<=');
checks = local_add_check(checks, 'Hs0_flat_model_diff', ...
    case_results.flat_diff_max_abs, tol_invariant, '<=');

energy_case = params_base;
energy_case.surface_boundary_model = 'kirchhoff_kstat';
energy_case.sea_hs_target = 0.2;
energy_case.surface_kstat_random_scatter = false;
energy_case.surface_kstat_trusted_angle_deg = 20;
fprintf('Running kirchhoff_kstat coherent-only metadata/energy case.\n');
channel_energy = vertical_channel_model(energy_case);
meta_energy = channel_energy.roughness_meta.kirchhoff_kstat_meta;
k0 = 2*pi*energy_case.f0/energy_case.c0;
expected_R_coh = energy_case.surface_reflect_coeff * ...
    exp(-2*k0^2*meta_energy.sigma_eta2_m2);
case_results.energy = local_compact_channel_result(channel_energy);
case_results.energy.kstat_meta = local_compact_kstat_meta(meta_energy);

checks = local_add_check(checks, 'energy_invariant', ...
    local_invariant_error(channel_energy), tol_invariant, '<=');
checks = local_add_check(checks, 'R_coh_normal_formula', ...
    abs(meta_energy.R_coh - expected_R_coh), tol_roundoff, '<=');
checks = local_add_check(checks, 'phase_screen_energy_closure', ...
    meta_energy.phase_screen_energy_error, tol_energy, '<=');
checks = local_add_check(checks, 'propagating_energy_recorded', ...
    double(isfinite(meta_energy.phase_screen_incoherent_energy_propagating)), 1, '==');
checks = local_add_check(checks, 'trusted_angle_energy_recorded', ...
    double(isfinite(meta_energy.phase_screen_incoherent_energy_trusted_angle)), 1, '==');
checks = local_add_check(checks, 'trusted_energy_not_larger_than_propagating', ...
    meta_energy.phase_screen_incoherent_energy_propagating - ...
    meta_energy.phase_screen_incoherent_energy_trusted_angle, -tol_roundoff, '>=');
checks = local_add_check(checks, 'metadata_only_keeps_P_sca', ...
    double(meta_energy.P_sca.sum > 0 && meta_energy.E_sca_realized == 0), 1, '==');
checks = local_add_check(checks, 'metadata_fields_present', ...
    double(local_has_kstat_metadata(meta_energy)), 1, '==');

random_case = energy_case;
random_case.surface_kstat_random_scatter = true;
fprintf('Running deterministic random realization repeat case.\n');
channel_random_a = vertical_channel_model(random_case);
channel_random_b = vertical_channel_model(random_case);
case_results.repeat_max_abs_H_f = max(abs(channel_random_a.H_f(:) - channel_random_b.H_f(:)));
checks = local_add_check(checks, 'same_seed_repeat_H_f', ...
    case_results.repeat_max_abs_H_f, tol_roundoff, '<=');

seed_changed = random_case;
seed_changed.sea_seed = random_case.sea_seed + 1;
fprintf('Running seed sensitivity case.\n');
channel_seed_changed = vertical_channel_model(seed_changed);
case_results.changed_seed_direct_drift = max(abs(channel_random_a.H_direct_f(:) - channel_seed_changed.H_direct_f(:)));
case_results.changed_seed_reflect_delta = max(abs(channel_random_a.H_reflect_f(:) - channel_seed_changed.H_reflect_f(:)));
checks = local_add_check(checks, 'changed_seed_direct_stable', ...
    case_results.changed_seed_direct_drift, tol_roundoff, '<=');
checks = local_add_check(checks, 'changed_seed_reflect_changes', ...
    case_results.changed_seed_reflect_delta, 1e-14, '>');

zero_padded = energy_case;
zero_padded.surface_kstat_conv_padding = 'zero_padded';
fprintf('Running zero-padded convolution metadata case.\n');
channel_zero_padded = vertical_channel_model(zero_padded);
meta_zero_padded = channel_zero_padded.roughness_meta.kirchhoff_kstat_meta;
case_results.zero_padded = local_compact_kstat_meta(meta_zero_padded);
checks = local_add_check(checks, 'zero_padded_energy_closure', ...
    meta_zero_padded.phase_screen_energy_error, tol_energy, '<=');
checks = local_add_check(checks, 'zero_padded_conv_mode', ...
    double(strcmp(meta_zero_padded.conv_padding, 'zero_padded')), 1, '==');

weak_kstat = params_base;
weak_kstat.surface_boundary_model = 'kirchhoff_kstat';
weak_kstat.surface_kstat_random_scatter = false;
weak_kstat.sea_hs_target = 0.05;

weak_ssa1 = params_base;
weak_ssa1.surface_boundary_model = 'ssa_stat_kernel';
weak_ssa1.surface_ssa_kernel_mode = 'ssa1_geometry';
weak_ssa1.surface_ssa_random_scatter = false;
weak_ssa1.sea_hs_target = weak_kstat.sea_hs_target;
weak_ssa1.surface_ssa_geometry_source_id = ['SSA.md; Thorsos & Broschat 1995 JASA, ', ...
    'Dirichlet SSA first-order / perturbation-limit geometry'];
fprintf('Running weak-roughness kstat and SSA1 coherent-reference cases.\n');
channel_weak_kstat = vertical_channel_model(weak_kstat);
channel_weak_ssa1 = vertical_channel_model(weak_ssa1);
meta_weak_kstat = channel_weak_kstat.roughness_meta.kirchhoff_kstat_meta;
meta_weak_ssa1 = channel_weak_ssa1.roughness_meta.ssa_stat_kernel_meta;
case_results.weak_reference = struct( ...
    'R_coh_kstat', meta_weak_kstat.R_coh, ...
    'R_coh_ssa1', meta_weak_ssa1.R_coh, ...
    'abs_diff', abs(meta_weak_kstat.R_coh - meta_weak_ssa1.R_coh));
checks = local_add_check(checks, 'weak_kstat_ssa1_coherent_reference', ...
    case_results.weak_reference.abs_diff, 1e-12, '<=');

mc_count = local_env_int('KSTAT_VALIDATE_ENSEMBLE_COUNT', 6);
fprintf('Running explicit Kirchhoff ensemble comparison with %d seeds.\n', mc_count);
[ensemble_result, ensemble_check_value] = local_explicit_kirchhoff_ensemble( ...
    params_base, energy_case, channel_energy, mc_count);
case_results.explicit_kirchhoff_ensemble = ensemble_result;
checks = local_add_check(checks, 'explicit_kirchhoff_ensemble_same_order', ...
    ensemble_check_value, 1.0, '<=');

save(result_file, 'case_results', 'checks');
fprintf('Saved validation result: %s\n', result_file);

failed = checks(~[checks.passed]);
if ~isempty(failed)
    disp(struct2table(checks))
    error('validate_kirchhoff_kstat_vertical:FailedChecks', ...
        '%d Kirchhoff kstat validation checks failed.', numel(failed));
end

disp(struct2table(checks))
fprintf('Kirchhoff kstat validation completed successfully.\n');

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
paramsV.surface_phase_mode = 'normal';
paramsV.surface_boundary_model = 'kirchhoff_spatial';
paramsV.surface_boundary_coupling_diagnostics = false;
paramsV.surface_boundary_redistribution_diagnostics = false;
paramsV.surface_kstat_random_scatter = true;
paramsV.surface_kstat_seed_offset = 200000;
paramsV.surface_kstat_conv_padding = 'periodic';
paramsV.surface_kstat_trusted_angle_deg = NaN;
end

function err = local_invariant_error(channel)
err = max(abs(channel.H_f(:) - channel.H_direct_f(:) - channel.H_reflect_f(:)));
end

function out = local_compact_channel_result(channel)
out = struct( ...
    'h_direct', channel.h_direct, ...
    'h_reflect', channel.h_reflect, ...
    'h_total', channel.h_total, ...
    'H_f', channel.H_f, ...
    'H_direct_f', channel.H_direct_f, ...
    'H_reflect_f', channel.H_reflect_f, ...
    'roughness_enabled', channel.roughness_meta.enabled, ...
    'surface_realization_generated', channel.roughness_meta.surface_realization_generated);
end

function out = local_compact_kstat_meta(meta)
out = struct( ...
    'enabled', meta.enabled, ...
    'model', meta.model, ...
    'random_scatter_enabled', meta.random_scatter_enabled, ...
    'conv_padding', meta.conv_padding, ...
    'sigma_eta2_m2', meta.sigma_eta2_m2, ...
    'alpha_rad_per_m', meta.alpha_rad_per_m, ...
    'G_mean', meta.G_mean, ...
    'R_coh', meta.R_coh, ...
    'S_deltaG_sum', meta.S_deltaG.sum, ...
    'P_sca_sum', meta.P_sca.sum, ...
    'phase_screen_energy_closure', meta.phase_screen_energy_closure, ...
    'phase_screen_energy_error', meta.phase_screen_energy_error, ...
    'phase_screen_incoherent_energy', meta.phase_screen_incoherent_energy, ...
    'phase_screen_incoherent_energy_propagating', meta.phase_screen_incoherent_energy_propagating, ...
    'phase_screen_incoherent_energy_trusted_angle', meta.phase_screen_incoherent_energy_trusted_angle, ...
    'E_inc', meta.E_inc, ...
    'E_coh', meta.E_coh, ...
    'E_sca_expected', meta.E_sca_expected, ...
    'E_sca_realized', meta.E_sca_realized, ...
    'E_ref', meta.E_ref, ...
    'seed_kstat', meta.seed_kstat);
end

function tf = local_has_kstat_metadata(meta)
names = {'sigma_eta2_m2', 'alpha_rad_per_m', 'G_mean', 'R_coh', ...
    'S_deltaG', 'P_sca', 'phase_screen_energy_closure', ...
    'phase_screen_energy_error', 'phase_screen_incoherent_energy_propagating', ...
    'P_sca_propagating_sum', 'seed_kstat', 'fft_normalization_note', ...
    'incident_spectrum_stats', 'coherent_spectrum_stats', ...
    'scatter_power_spectrum_stats', 'reflected_spectrum_stats'};
tf = true;
for ii = 1:numel(names)
    tf = tf && isfield(meta, names{ii});
end
tf = tf && isfield(meta.P_sca, 'sum') && isfield(meta.S_deltaG, 'sum');
end

function [result, rel_error] = local_explicit_kirchhoff_ensemble(params_base, kstat_case, channel_kstat_coh, mc_count)
mc_count = max(1, mc_count);
h_reflect = complex(zeros(mc_count, 1));
R_screen = complex(zeros(mc_count, 1));
k0 = 2*pi*kstat_case.f0/kstat_case.c0;
for ii = 1:mc_count
    params_i = params_base;
    params_i.surface_boundary_model = 'kirchhoff_spatial';
    params_i.sea_hs_target = kstat_case.sea_hs_target;
    params_i.sea_seed = params_base.sea_seed + ii - 1;
    ch_i = vertical_channel_model(params_i);
    h_reflect(ii) = ch_i.h_reflect;
    R_screen(ii) = params_i.surface_reflect_coeff * ...
        mean(exp(1i * 2 * k0 * ch_i.surface_elevation(:)));
end
ensemble_mean = mean(h_reflect);
kstat_coherent_ref = channel_kstat_coh.h_reflect;
kstat_R_coh = channel_kstat_coh.roughness_meta.kirchhoff_kstat_meta.R_coh;
R_screen_mean = mean(R_screen);
rel_error = abs(R_screen_mean - kstat_R_coh) / max(abs(kstat_R_coh), eps);
result = struct( ...
    'mc_count', mc_count, ...
    'explicit_mean_h_reflect', ensemble_mean, ...
    'explicit_std_abs_h_reflect', std(abs(h_reflect)), ...
    'kstat_coherent_h_reflect', kstat_coherent_ref, ...
    'explicit_screen_R_mean', R_screen_mean, ...
    'explicit_screen_R_std_abs', std(abs(R_screen)), ...
    'kstat_R_coh', kstat_R_coh, ...
    'relative_error', rel_error, ...
    'interpretation', 'Finite-seed explicit Kirchhoff phase-screen mean comparison; receiver h_reflect is recorded but not used as the pass metric.');
end

function checks = local_add_check(checks, name, value, threshold, op)
switch op
    case '<='
        passed = value <= threshold;
    case '>='
        passed = value >= threshold;
    case '>'
        passed = value > threshold;
    case '=='
        passed = value == threshold;
    otherwise
        error('Unsupported check operator: %s', op);
end
entry = struct( ...
    'name', string(name), ...
    'value', value, ...
    'threshold', threshold, ...
    'operator', string(op), ...
    'passed', logical(passed));
checks = [checks; entry]; %#ok<AGROW>
end

function value = local_env_int(name, default_value)
raw = getenv(name);
if isempty(raw)
    value = default_value;
    return
end
parsed = str2double(raw);
if isfinite(parsed) && parsed >= 1
    value = round(parsed);
else
    value = default_value;
end
end
