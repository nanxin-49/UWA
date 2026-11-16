function meta = write_bellhop_rotated_REAL8_A_inputs(case_root,cfg)
%WRITE_BELLHOP_ROTATED_MULTIFIELD_Q Explicit requested-order wall sidecars.
arguments
    case_root (1,:) char
    cfg (1,1) struct
end
assert(isfield(cfg,'reflection_order') && cfg.reflection_order>=1);
meta=write_bellhop_rotated_munk_q(case_root,cfg);
fid=fopen([case_root '.iwcov'],'w'); assert(fid>=0); closer=onCleanup(@()fclose(fid));
fprintf(fid,'5\n%.17g\n%.17g\n0\n3\n',cfg.wall_range_m,max(cfg.receiver_range_m));
fprintf(fid,'%.17g %.17g\n',cfg.wall_range_m,cfg.transverse_origin_m-cfg.artificial_halfwidth_m);
fprintf(fid,'%.17g %.17g\n',cfg.wall_range_m,cfg.transverse_origin_m);
fprintf(fid,'%.17g %.17g\n',cfg.wall_range_m,cfg.transverse_origin_m+cfg.artificial_halfwidth_m); clear closer
omega=2*pi*cfg.frequency_hz;
alpha=cfg.bottom_atten_db_wavelength*log(10)/20*cfg.frequency_hz/cfg.bottom_cp_mps;
cp_native=conj(omega/(omega/cfg.bottom_cp_mps+1i*alpha));
kind='V'; if strcmp(cfg.reflection_wall,'bottom'), kind='A'; end
fid=fopen([case_root '.rotq'],'w'); assert(fid>=0); closer=onCleanup(@()fclose(fid));
fprintf(fid,'%.17g %s %.17g %.17g %.17g\n%d %.17g\n',cfg.transverse_origin_m,kind, ...
    real(cp_native),imag(cp_native),cfg.bottom_rho_gcm3,cfg.reflection_order,cfg.water_depth_m); clear closer
write_bellhop_rotated_field_options(case_root,cfg.receiver_range_m,meta.q_receiver_m,point3d=true,fast_q=true,native_range_metadata_repair=false);
meta.reflection_order=cfg.reflection_order;
meta.expected_surface_bounces=cfg.expected_surface_bounces;
meta.expected_bottom_bounces=cfg.expected_bottom_bounces;
meta.last_wall_range_m=cfg.last_wall_range_m;
meta.component=sprintf('%s launch, exactly %d physical reflections; not total field',cfg.direction,cfg.reflection_order);
meta.bottom_cp_native_complex=cp_native;
end
