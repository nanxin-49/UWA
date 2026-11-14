function cfg = make_bellhop_rotated_direct_field_config(direction,depths_m,x_m,options)
%MAKE_BELLHOP_ROTATED_DIRECT_FIELD_CONFIG Exact direct chart and permutation.
arguments
    direction (1,:) char {mustBeMember(direction,{'up','down','left','right'})}
    depths_m (1,:) double {mustBeFinite,mustBeNonnegative}
    x_m (1,:) double {mustBeFinite}
    options.Q_spacing_m (1,1) double {mustBePositive}=.625
    options.step_m (1,1) double {mustBePositive}=.5
    options.beam_count (1,1) double {mustBeInteger,mustBePositive}=2001
end
cfg=make_bellhop_rotated_munk_config('up'); cfg.run_type='C';
cfg.transverse_origin_m=0; cfg.artificial_halfwidth_m=1000;
cfg.ssp_spacing_m=options.Q_spacing_m; cfg.step_m=options.step_m;
cfg.beam_count=options.beam_count; cfg.point3d=true;
cfg.nominal_receiver_depths_m=depths_m; cfg.nominal_receiver_x_m=x_m;
assert(all(depths_m<=cfg.water_depth_m));
if ismember(direction,{'up','down'})
    cfg.chart_kind='direct'; cfg.direction=direction;
    sg=1-2*strcmp(direction,'up'); rr=sg*(depths_m-cfg.source_z_m);
    assert(all(rr>0));
    [cfg.receiver_range_m,cfg.receiver_depth_order]=sort(rr);
    [~,cfg.receiver_x_order]=sort(-sg*x_m);
    cfg.receiver_x_m=x_m(cfg.receiver_x_order);
    cfg.receiver_depths_m=depths_m(cfg.receiver_depth_order);
    cfg.range_limit_m=max(rr);
    ratio=max(abs(x_m))/min(rr);
else
    cfg.horizontal_ssp_spacing_m=options.Q_spacing_m;
    cfg.chart_kind='horizontal'; cfg.horizontal_side=1-2*strcmp(direction,'left');
    side=cfg.horizontal_side; rr=side*x_m; assert(all(rr>0));
    [cfg.receiver_range_m,cfg.receiver_x_order]=sort(rr);
    [cfg.receiver_q_m,cfg.receiver_depth_order]=sort(side*(depths_m-cfg.source_z_m));
    cfg.range_limit_m=max(rr); ratio=max(abs(cfg.receiver_q_m))/min(rr);
    cfg.step_m=min(cfg.step_m,.01);
end
% Sector coverage is assigned by the planner; diagonal overlap remains raw.
fan=min(46,1.1*rad2deg(atan(1.05*ratio))+.001);
cfg.angle_limits_deg=[-fan fan];
cfg.source_model_note='Derived 3D point-source chart, native complex pressure converted by fixed -conj; no fit.';
end