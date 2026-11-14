function factor = bellhop_rotated_multifield_reference_factor(cfg,meta,ranges_m,q_m)
%BELLHOP_ROTATED_MULTIFIELD_REFERENCE_FACTOR Independent exact-Q Snell/R product.
% Returns native-phasor reflection factors only, no pressure or amplitude fit.
arguments
    cfg (1,1) struct
    meta (1,1) struct
    ranges_m (1,:) double
    q_m (:,1) double
end
factor=complex(zeros(numel(q_m),numel(ranges_m)));
omega=2*pi*cfg.frequency_hz; cs=meta.c_mps(1);
cw=munk_profile_sound_speed(cfg.water_depth_m,cfg.munk_parameters);
alpha=cfg.bottom_atten_db_wavelength*log(10)/20*cfg.frequency_hz/cfg.bottom_cp_mps;
kb=omega/cfg.bottom_cp_mps+1i*alpha;
bound=sin(deg2rad(max(abs(cfg.angle_limits_deg))))/cs;
for ir=1:numel(ranges_m)
    R=ranges_m(ir); nodes=[meta.range_m(meta.range_m<R);R];
    cc=interp1(meta.range_m,meta.c_mps,nodes);
    h=diff(nodes); c0=cc(1:end-1); c1=cc(2:end);
    for iz=1:numel(q_m)
        target=q_m(iz)-cfg.transverse_origin_m;
        fun=@(s)sum(s*(c0+c1).*h./(sqrt(1-(s*c0).^2)+sqrt(1-(s*c1).^2)))-target;
        if target==0, slow=0; else, slow=fzero(fun,[-bound bound]); end
        qw=sqrt((omega/cw)^2-(omega*slow)^2);
        qb=sqrt(kb^2-(omega*slow)^2); if imag(qb)<0, qb=-qb; end
        Rb=(qw-qb/cfg.bottom_rho_gcm3)/(qw+qb/cfg.bottom_rho_gcm3);
        factor(iz,ir)=(-1)^cfg.expected_surface_bounces*conj(Rb)^cfg.expected_bottom_bounces;
    end
end
end
