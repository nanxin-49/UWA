function [sample,meta]=sample_conditional_channel_library_vertical(library,wind_speed_mps,n_samples,seed,options)
%SAMPLE_CONDITIONAL_CHANNEL_LIBRARY_VERTICAL Sample one exact discrete node.
%   No interpolation or implicit nearest-node selection is performed.

arguments
    library (1,1) struct
    wind_speed_mps (1,1) double {mustBeFinite}
    n_samples (1,1) double {mustBeInteger,mustBePositive}
    seed (1,1) double {mustBeFinite}
    options (1,1) struct = struct()
end
[library,library_migration]=upgrade_conditional_channel_phase_vertical(library,struct());
if ~isfield(library,'conditions') || isempty(library.conditions)
    error('The conditional channel library has no condition nodes.');
end
wind=[library.conditions.wind_speed_mps]; idx=find(abs(wind-wind_speed_mps)<=1e-12,1);
if isempty(idx)
    error('Unsupported wind speed %.12g m/s. Available exact nodes: %s. Interpolation is disabled.', ...
        wind_speed_mps,strjoin(string(wind),', '));
end
if ~isfield(library.conditions(idx),'model')
    error('Condition node %.12g m/s does not contain its sampling model.',wind_speed_mps);
end
[sample,meta]=sample_conditional_channel_vertical( ...
    library.conditions(idx).model,n_samples,seed,options);
meta.library_schema_version=library.schema_version;
meta.condition_index=idx;
meta.wind_speed_mps=wind_speed_mps;
meta.interpolation_used=false;
meta.library_phase_migration=library_migration;
end
