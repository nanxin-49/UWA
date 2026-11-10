function out = run_pe_physical_surface_model0_validation(cfg)
%RUN_PE_PHYSICAL_SURFACE_MODEL0_VALIDATION Depth-positive Model-0 adapter.
% cfg.eta_z_m is a physical boundary depth. The legacy validation helper
% accepts upward elevation, so convert once here without changing it.
arguments
    cfg (1,1) struct
end
assert(isfield(cfg,'eta_z_m'),'Physical eta_z_m is required.');
legacy = cfg;
legacy.surface_elevation_x_m = -cfg.eta_z_m;
legacy = rmfield(legacy,'eta_z_m');
out = run_pe_1d_surface_reflection_validation(legacy);
out.convention_version = 'S2a-depth-positive-native-v1';
out.physical_eta_z_m = cfg.eta_z_m(:).';
out.adapter_eta_up_m = legacy.surface_elevation_x_m(:).';
out.physical_screen_formula = 'R0*exp(-2i*k*eta_z)';
out.native_carrier = exp(1i*out.k0_rad_per_m*(cfg.z_tx_m+cfg.z_rx_m));
out.reflected_native = out.reflected_field*out.native_carrier;
out.phase_reference = 'native exp(-i omega t), reflected component';
end
