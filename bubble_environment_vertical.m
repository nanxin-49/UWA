function [c_eff_xy, alpha_bub_xy, meta] = ...
    bubble_environment_vertical(x, y, z_curr, f_hz, c_bg, cfg)
%BUBBLE_ENVIRONMENT_VERTICAL Disabled-by-default bubble environment hook.
% Priority 1 only wires the propagation framework. Physical bubble models
% are intentionally left unimplemented until later stages.

if nargin < 6
    error('bubble_environment_vertical requires x, y, z_curr, f_hz, c_bg, and cfg.');
end
if ~(isnumeric(x) && isvector(x) && all(isfinite(x(:))))
    error('x must be a finite numeric vector.');
end
if ~(isnumeric(y) && isvector(y) && all(isfinite(y(:))))
    error('y must be a finite numeric vector.');
end
if ~(isscalar(z_curr) && isnumeric(z_curr) && isfinite(z_curr))
    error('z_curr must be a finite numeric scalar.');
end
if ~(isscalar(f_hz) && isnumeric(f_hz) && isfinite(f_hz) && f_hz > 0)
    error('f_hz must be a positive finite numeric scalar.');
end
if ~(isscalar(c_bg) && isnumeric(c_bg) && isfinite(c_bg) && c_bg > 0)
    error('c_bg must be a positive finite numeric scalar.');
end

model = 'off';
if isfield(cfg, 'bubble_model') && ~isempty(cfg.bubble_model)
    model = cfg.bubble_model;
end
if isstring(model)
    model = char(model);
end
if ~ischar(model)
    error('cfg.bubble_model must be a string or char.');
end
model = lower(model);

enabled = false;
if isfield(cfg, 'enable_bubbles') && ~isempty(cfg.enable_bubbles)
    enabled = logical(cfg.enable_bubbles);
end

if enabled && ~strcmp(model, 'off')
    error('bubble_model=%s is not implemented in Priority 1.', model);
end

ny = numel(y);
nx = numel(x);
c_eff_xy = c_bg * ones(ny, nx);
alpha_bub_xy = zeros(ny, nx);

meta = struct();
meta.enabled = false;
meta.model = 'off';
meta.spatial_mode = 'none';
if isfield(cfg, 'bubble_spatial_mode') && ~isempty(cfg.bubble_spatial_mode)
    spatial_mode = cfg.bubble_spatial_mode;
    if isstring(spatial_mode)
        spatial_mode = char(spatial_mode);
    end
    if ischar(spatial_mode)
        meta.spatial_mode = lower(spatial_mode);
    end
end
meta.apply_sound_speed = false;
meta.apply_attenuation = false;
if isfield(cfg, 'bubble_apply_sound_speed') && ~isempty(cfg.bubble_apply_sound_speed)
    meta.apply_sound_speed = logical(cfg.bubble_apply_sound_speed);
end
if isfield(cfg, 'bubble_apply_attenuation') && ~isempty(cfg.bubble_apply_attenuation)
    meta.apply_attenuation = logical(cfg.bubble_apply_attenuation);
end
meta.z_m = z_curr;
meta.f_hz = f_hz;
meta.c_bg_mps = c_bg;
meta.c_eff_stats = local_field_stats(c_eff_xy);
meta.alpha_bub_stats = local_field_stats(alpha_bub_xy);
meta.warning_flags = {};

end

function stats = local_field_stats(v)
if isempty(v)
    stats = struct('min', NaN, 'max', NaN, 'mean', NaN, 'std', NaN, 'finite_count', 0);
    return
end
finite_mask = isfinite(v);
if ~any(finite_mask(:))
    stats = struct('min', NaN, 'max', NaN, 'mean', NaN, 'std', NaN, 'finite_count', 0);
    return
end
vf = v(finite_mask);
stats = struct( ...
    'min', min(vf(:)), ...
    'max', max(vf(:)), ...
    'mean', mean(vf(:)), ...
    'std', std(vf(:)), ...
    'finite_count', nnz(finite_mask));
end
