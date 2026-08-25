function [psi_init_xy,meta] = virtual_point_source_initial_field_vertical(X,Y,f_hz,cfg)
%VIRTUAL_POINT_SOURCE_INITIAL_FIELD_VERTICAL Point-source-referenced PE input.
% The virtual source is s0 behind the PE input plane. The axial phase k*s0
% is removed, so physical source-referenced results require exp(i*k*s0)
% in addition to the production input-plane carrier restoration.

arguments
    X double
    Y double
    f_hz (1,1) double {mustBeFinite,mustBePositive}
    cfg (1,1) struct
end
if ~isequal(size(X),size(Y))
    error('X and Y must have identical sizes.');
end
if ~isfield(cfg,'virtual_source_distance_m')
    error('cfg.virtual_source_distance_m is required.');
end
s0=cfg.virtual_source_distance_m;
validateattributes(s0,{'numeric'},{'scalar','finite','positive'});
if isfield(cfg,'source_green_amplitude')
    amplitude=cfg.source_green_amplitude;
else
    amplitude=1/(4*pi);
end
validateattributes(amplitude,{'numeric'},{'scalar','finite'});
rho2=(X-cfg.x_tx).^2+(Y-cfg.y_tx).^2;
R0=sqrt(s0^2+rho2);
k=2*pi*f_hz/cfg.c0;
psi_init_xy=amplitude.*exp(1i*k*(R0-s0))./R0;
meta=struct('reference','virtual_source_axial_phase_removed', ...
    'virtual_source_distance_m',s0,'green_amplitude',amplitude, ...
    'carrier_to_source_reference',exp(1i*k*s0), ...
    'formula','A0*exp(i*k*(R0-s0))/R0');
end
