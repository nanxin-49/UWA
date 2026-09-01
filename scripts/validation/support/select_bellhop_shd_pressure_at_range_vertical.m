function [pressure, meta] = select_bellhop_shd_pressure_at_range_vertical(shade, target_range_m, target_depth_m, tolerance_m)
%SELECT_BELLHOP_SHD_PRESSURE_AT_RANGE_VERTICAL Select one explicit SHD sample.
% This validation helper deliberately rejects nearest-column selection.  A
% requested receiver range and depth must each match exactly one stored SHD
% coordinate within the supplied tolerance.
arguments
    shade (1,1) struct
    target_range_m (1,1) double {mustBeFinite}
    target_depth_m (1,1) double {mustBeFinite}
    tolerance_m (1,1) double {mustBePositive} = 1e-6
end
required={'receiver_range_m','receiver_depth_m','pressure'};
for ii=1:numel(required)
    if ~isfield(shade,required{ii})
        error('BellhopRangeSelection:MissingField','SHD data is missing field %s.',required{ii});
    end
end
ranges=double(shade.receiver_range_m(:));
depths=double(shade.receiver_depth_m(:));
range_ix=find(abs(ranges-target_range_m)<=tolerance_m);
depth_ix=find(abs(depths-target_depth_m)<=tolerance_m);
if numel(range_ix)~=1
    error('BellhopRangeSelection:RangeNotUnique', ...
        'Expected exactly one SHD range at %.12g m within %.3g m; found %d in [%s].', ...
        target_range_m,tolerance_m,numel(range_ix),strjoin(compose('%.12g',ranges.'),', '));
end
if numel(depth_ix)~=1
    error('BellhopRangeSelection:DepthNotUnique', ...
        'Expected exactly one SHD depth at %.12g m within %.3g m; found %d in [%s].', ...
        target_depth_m,tolerance_m,numel(depth_ix),strjoin(compose('%.12g',depths.'),', '));
end
if ~isequal(size(shade.pressure),[numel(depths),numel(ranges)])
    error('BellhopRangeSelection:PressureShape', ...
        'SHD pressure size [%s] does not match %d depths by %d ranges.', ...
        num2str(size(shade.pressure)),numel(depths),numel(ranges));
end
pressure=shade.pressure(depth_ix,range_ix);
meta=struct('requested_range_m',target_range_m,'actual_range_m',ranges(range_ix), ...
    'range_index',range_ix,'requested_depth_m',target_depth_m, ...
    'actual_depth_m',depths(depth_ix),'depth_index',depth_ix, ...
    'tolerance_m',tolerance_m);
end
