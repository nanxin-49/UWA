function cfg = make_bellhop_rotated_multifield_config(direction,order,depths_m,x_m)
%MAKE_BELLHOP_ROTATED_MULTIFIELD_CONFIG Explicit physical receiver permutation.
arguments
    direction (1,:) char {mustBeMember(direction,{'up','down'})}
    order (1,1) double {mustBeInteger,mustBePositive}
    depths_m (1,:) double {mustBeFinite,mustBeNonnegative}
    x_m (1,:) double {mustBeFinite}
end
assert(order<=50);
cfg=make_bellhop_rotated_munk_config(direction);
assert(all(depths_m<=cfg.water_depth_m));
cfg.run_type='C'; cfg.source_geometry='R'; cfg.reflection_order=order;
sg=1-2*strcmp(direction,'up');
if strcmp(direction,'up')
    cfg.reflection_wall='surface'; cfg.wall_range_m=cfg.source_z_m;
    cfg.expected_surface_bounces=ceil(order/2);
else
    cfg.reflection_wall='bottom'; cfg.wall_range_m=cfg.water_depth_m-cfg.source_z_m;
    cfg.expected_surface_bounces=floor(order/2);
end
cfg.expected_bottom_bounces=order-cfg.expected_surface_bounces;
cfg.last_wall_range_m=cfg.wall_range_m+(order-1)*cfg.water_depth_m;
cfg.last_wall_surface=(strcmp(direction,'up') && mod(order,2)==1) || ...
    (strcmp(direction,'down') && mod(order,2)==0);
if cfg.last_wall_surface
    r=cfg.last_wall_range_m+depths_m;
else
    r=cfg.last_wall_range_m+cfg.water_depth_m-depths_m;
end
[cfg.receiver_range_m,cfg.receiver_depth_order]=sort(r);
cfg.receiver_depths_m=depths_m(cfg.receiver_depth_order);
coefficient=-sg*(-1)^order;
[~,cfg.receiver_x_order]=sort(coefficient*x_m);
cfg.receiver_x_m=x_m(cfg.receiver_x_order);
cfg.nominal_receiver_x_m=x_m; cfg.nominal_receiver_depths_m=depths_m;
cfg.transverse_origin_m=0; cfg.artificial_halfwidth_m=1000;
cfg.range_limit_m=cfg.last_wall_range_m+cfg.water_depth_m;
cfg.ssp_spacing_m=5; cfg.step_m=2; cfg.beam_count=301;
if order==1 && strcmp(direction,'down'), cfg.angle_limits_deg=[-2 2];
else, cfg.angle_limits_deg=[-.2 .2]; end
cfg.source_model_note='Omni 3D point with validation geometric outplane transport; full FFP qualification pending.';
end
