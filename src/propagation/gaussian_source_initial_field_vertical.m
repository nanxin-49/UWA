function [psi0_xy,meta]=gaussian_source_initial_field_vertical(X,Y,cfg)
%GAUSSIAN_SOURCE_INITIAL_FIELD_VERTICAL Production Gaussian initial field.
% This shared helper preserves the original production expression exactly.

arguments
    X double
    Y double
    cfg (1,1) struct
end
if ~isequal(size(X),size(Y))
    error('X and Y must have identical sizes.');
end
required={'x_tx','y_tx','sigma_src_m'};
for ii=1:numel(required)
    if ~isfield(cfg,required{ii})
        error('cfg.%s is required.',required{ii});
    end
end
psi0_xy=exp(-((X-cfg.x_tx).^2+(Y-cfg.y_tx).^2)/(2*cfg.sigma_src_m^2));
psi0_xy=complex(psi0_xy,0);
meta=struct('mode','gaussian','expression', ...
    'exp(-((x-x_tx)^2+(y-y_tx)^2)/(2*sigma_src_m^2))', ...
    'sigma_src_m',cfg.sigma_src_m,'x_tx_m',cfg.x_tx,'y_tx_m',cfg.y_tx, ...
    'peak_normalization',1,'frequency_dependent',false);
end
