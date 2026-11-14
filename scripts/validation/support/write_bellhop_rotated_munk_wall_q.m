function meta = write_bellhop_rotated_munk_wall_q(case_root,cfg)
%WRITE_BELLHOP_ROTATED_MUNK_WALL_Q Single physical flat wall in unfolded Q.
% Bottom cp is chosen to match the FFP complex wavenumber exactly after
% phasor conjugation; native CRCI's small-loss c+i*ci approximation is not used.
arguments
    case_root (1,:) char
    cfg (1,1) struct
end
assert(isfield(cfg,'reflection_wall') && isfield(cfg,'wall_range_m'));
meta=write_bellhop_rotated_munk_q(case_root,cfg);
fid=fopen([case_root '.iwcov'],'w'); assert(fid>=0); closer=onCleanup(@()fclose(fid));
fprintf(fid,'3\n%.17g\n%.17g\n0\n3\n',cfg.wall_range_m,max(cfg.receiver_range_m));
fprintf(fid,'%.17g %.17g\n',cfg.wall_range_m,cfg.transverse_origin_m-cfg.artificial_halfwidth_m);
fprintf(fid,'%.17g %.17g\n',cfg.wall_range_m,cfg.transverse_origin_m);
fprintf(fid,'%.17g %.17g\n',cfg.wall_range_m,cfg.transverse_origin_m+cfg.artificial_halfwidth_m);
clear closer
boundary='V'; cp_native=complex(cfg.bottom_cp_mps);
if strcmp(cfg.reflection_wall,'bottom')
    boundary='A'; omega=2*pi*cfg.frequency_hz;
    alpha=cfg.bottom_atten_db_wavelength*log(10)/20*cfg.frequency_hz/cfg.bottom_cp_mps;
    k_ffp=omega/cfg.bottom_cp_mps+1i*alpha;
    cp_native=conj(omega/k_ffp);
end
fid=fopen([case_root '.rotq'],'w'); assert(fid>=0); closer=onCleanup(@()fclose(fid));
fprintf(fid,'%.17g %s %.17g %.17g %.17g\n',cfg.transverse_origin_m,boundary, ...
    real(cp_native),imag(cp_native),cfg.bottom_rho_gcm3);
meta.wall_boundary=boundary; meta.bottom_cp_native_complex=cp_native;
meta.reflection_order=1; meta.component='single reflection only, not total field';
end