function plan = make_bellhop_rotated_field_plan(depths_m,x_m,options)
%MAKE_BELLHOP_ROTATED_FIELD_PLAN Exact chart partitions and explicit orders.
% Planning alone runs no solver. Full-grid execution requires separate approval.
arguments
    depths_m (1,:) double {mustBeFinite,mustBeNonnegative}
    x_m (1,:) double {mustBeFinite}
    options.Q_spacing_m (1,1) double {mustBePositive}=.625
    options.step_m (1,1) double {mustBePositive}=.5
    options.reflection_beams (1,1) double {mustBeInteger,mustBePositive}=401
    options.direct_beams (1,1) double {mustBeInteger,mustBePositive}=2001
    options.max_order (1,1) double {mustBeInteger,mustBePositive}=16
end
assert(all(diff(depths_m)>0) && all(diff(x_m)>0));
base=make_bellhop_rotated_munk_config('up'); zs=base.source_z_m;
assert(max(depths_m)<=base.water_depth_m);
jobs=struct('tag',{},'cfg',{},'z_indices',{},'x_indices',{},'coverage',{},'order',{});
for direction={'up','down'}
    sg=1-2*strcmp(direction{1},'up'); distance=sg*(depths_m-zs);
    edges=[0 5 50 inf]; step_limits=[.01 .1 inf]*options.step_m;
    for b=1:3
        iz=find(distance>edges(b) & distance<=edges(b+1)); if isempty(iz), continue; end
        cfg=make_bellhop_rotated_direct_field_config(direction{1},depths_m(iz),x_m, ...
            Q_spacing_m=options.Q_spacing_m,step_m=min(options.step_m,step_limits(b)), ...
            beam_count=options.direct_beams);
        coverage=abs(depths_m(iz).'-zs)>=abs(x_m);
        jobs(end+1)=struct('tag',sprintf('direct_%s_part%d',direction{1},b), ...
            'cfg',cfg,'z_indices',iz,'x_indices',1:numel(x_m),'coverage',coverage,'order',0); %#ok<AGROW>
    end
end
iz=find(abs(depths_m-zs)<=max(abs(x_m)));
for direction={'left','right'}
    side=1-2*strcmp(direction{1},'left'); ix=find(side*x_m>0);
    if isempty(iz) || isempty(ix), continue; end
    cfg=make_bellhop_rotated_direct_field_config(direction{1},depths_m(iz),x_m(ix), ...
        Q_spacing_m=options.Q_spacing_m,step_m=.01*options.step_m,beam_count=options.direct_beams);
    coverage=abs(depths_m(iz).'-zs)<abs(x_m(ix));
    jobs(end+1)=struct('tag',['direct_' direction{1}], ...
        'cfg',cfg,'z_indices',iz,'x_indices',ix,'coverage',coverage,'order',0); %#ok<AGROW>
end
for n=1:options.max_order
    for direction={'up','down'}
        cfg=make_bellhop_rotated_multifield_config(direction{1},n,depths_m,x_m);
        cfg.ssp_spacing_m=options.Q_spacing_m; cfg.step_m=options.step_m;
        cfg.beam_count=options.reflection_beams;
        fan=1.1*rad2deg(atan(1.05*max(abs(x_m))/min(cfg.receiver_range_m)))+.001;
        cfg.angle_limits_deg=[-fan fan];
        jobs(end+1)=struct('tag',sprintf('%s_order%d',direction{1},n), ...
            'cfg',cfg,'z_indices',1:numel(depths_m),'x_indices',1:numel(x_m), ...
            'coverage',[],'order',n); %#ok<AGROW>
    end
end
plan=struct('status','PLANNED_NOT_EXECUTED','depths_m',depths_m,'x_m',x_m, ...
    'settings',options,'jobs',jobs,'physical_contract',base, ...
    'source_singularity_mask',depths_m(:)==zs & x_m==0, ...
    'phase_bridge','-conj(native), no fitted phase/amplitude', ...
    'source_geometry','3D omnidirectional point, 1m pressure reference', ...
    'qualification','Finite explicit-order field; convergence and independent FFP still required');
end