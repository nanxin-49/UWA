run(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'bootstrap_project.m'));
% Reduced-grid validation for optional incident/reflected wavefield diagnostics.

clear
format compact

tol = 1e-10;
base = local_base_params();

off_params = base;
off_params.surface_wavefield_diagnostics = false;
on_params = base;
on_params.surface_wavefield_diagnostics = true;

fprintf('Running diagnostics-disabled baseline.\n');
channel_off = vertical_channel_model(off_params);
fprintf('Running diagnostics-enabled equivalent case.\n');
channel_on = vertical_channel_model(on_params);

flat_params = on_params;
flat_params.sea_hs_target = 0;
flat_params.surface_ssa_random_scatter = false;
fprintf('Running Hs=0 pressure-release phase-reversal case.\n');
channel_flat = vertical_channel_model(flat_params);
flat_meta = channel_flat.surface_wavefield_meta;

kirchhoff_params = on_params;
kirchhoff_params.surface_boundary_model = 'kirchhoff_spatial';
fprintf('Running rough Kirchhoff wavefield-diagnostic case.\n');
channel_kirchhoff = vertical_channel_model(kirchhoff_params);

wideband_params = on_params;
wideband_params.f0 = [4000, 6000, 8000];
wideband_params.f_ref_hz = 6000;
fprintf('Running reduced wideband invariant case.\n');
channel_wideband = vertical_channel_model(wideband_params);

checks = struct();
checks.disabled_meta_is_disabled = ~channel_off.surface_wavefield_meta.enabled;
checks.enabled_meta_is_enabled = channel_on.surface_wavefield_meta.enabled;
checks.diagnostics_H_f_max_abs_diff = max(abs(channel_off.H_f - channel_on.H_f));
checks.diagnostics_H_direct_max_abs_diff = max(abs(channel_off.H_direct_f - channel_on.H_direct_f));
checks.diagnostics_H_reflect_max_abs_diff = max(abs(channel_off.H_reflect_f - channel_on.H_reflect_f));
checks.flat_surface_phase_reversal_rel_l2 = norm( ...
    flat_meta.surface_reflected_xy(:) + flat_meta.surface_incident_xy(:)) / ...
    max(norm(flat_meta.surface_incident_xy(:)), eps);
checks.flat_surface_magnitude_rel_l2 = norm( ...
    abs(flat_meta.surface_reflected_xy(:)) - abs(flat_meta.surface_incident_xy(:))) / ...
    max(norm(abs(flat_meta.surface_incident_xy(:))), eps);
checks.flat_invariant_error = local_invariant_error(channel_flat);
checks.rough_ssa_invariant_error = local_invariant_error(channel_on);
checks.rough_kirchhoff_invariant_error = local_invariant_error(channel_kirchhoff);
checks.wideband_invariant_error = local_invariant_error(channel_wideband);
checks.ssa_endpoint_error = local_endpoint_error(channel_on.surface_wavefield_meta);
checks.kirchhoff_endpoint_error = local_endpoint_error(channel_kirchhoff.surface_wavefield_meta);
checks.wideband_endpoint_error = local_endpoint_error(channel_wideband.surface_wavefield_meta);
checks.reference_frequency_hz = channel_wideband.surface_wavefield_meta.reference_frequency_hz;

validation_report = struct();
validation_report.tolerance = tol;
validation_report.all_passed = ...
    checks.disabled_meta_is_disabled && checks.enabled_meta_is_enabled && ...
    checks.diagnostics_H_f_max_abs_diff <= tol && ...
    checks.diagnostics_H_direct_max_abs_diff <= tol && ...
    checks.diagnostics_H_reflect_max_abs_diff <= tol && ...
    checks.flat_surface_phase_reversal_rel_l2 <= tol && ...
    checks.flat_surface_magnitude_rel_l2 <= tol && ...
    checks.flat_invariant_error <= tol && ...
    checks.rough_ssa_invariant_error <= tol && ...
    checks.rough_kirchhoff_invariant_error <= tol && ...
    checks.wideband_invariant_error <= tol && ...
    checks.ssa_endpoint_error <= tol && ...
    checks.kirchhoff_endpoint_error <= tol && ...
    checks.wideband_endpoint_error <= tol && ...
    abs(checks.reference_frequency_hz - 6000) <= eps;

save('validate_surface_wavefield_visualization_vertical_result.mat', ...
    'checks', 'validation_report', 'base');
disp(checks)
disp(validation_report)
if ~validation_report.all_passed
    error('Surface wavefield visualization validation failed.');
end

function paramsV = local_base_params()
paramsV = struct();
paramsV.f0 = 6000;
paramsV.c0 = 1500;
paramsV.z_max = 20;
paramsV.z_tx = 20;
paramsV.z_rx = 2;
paramsV.xw = 15;
paramsV.yw = 15;
paramsV.nx = 64;
paramsV.ny = 64;
paramsV.sigma_src_m = 0.3;
paramsV.stepz_lamb = 1;
paramsV.sponge_ratio = 0.12;
paramsV.alpha_max_np_per_m = 0.15;
paramsV.env_mode = 'uniform';
paramsV.show_figures = false;
paramsV.enforce_1_over_R = false;
paramsV.save_mode = 'rx_only';
paramsV.enable_surface_reflection = true;
paramsV.surface_reflect_coeff = -1;
paramsV.surface_phase_mode = 'normal';
paramsV.surface_boundary_model = 'ssa_stat_kernel';
paramsV.surface_ssa_kernel_mode = 'ssa1_geometry';
paramsV.surface_ssa_random_scatter = true;
paramsV.surface_ssa_conv_padding = 'periodic';
paramsV.surface_ssa_geometry_source_id = ...
    'SSA.md; first-order pressure-release / Dirichlet geometry';
paramsV.sea_wind_speed = 5;
paramsV.sea_hs_target = 0.2;
paramsV.sea_seed = 12345;
paramsV.surface_wavefield_slice_axis = 'x';
paramsV.surface_wavefield_max_z_samples = 96;
end

function err = local_invariant_error(channel)
err = max(abs(channel.H_f(:) - channel.H_direct_f(:) - channel.H_reflect_f(:)));
end

function err = local_endpoint_error(meta)
err = max(struct2array(meta.endpoint_consistency));
end

