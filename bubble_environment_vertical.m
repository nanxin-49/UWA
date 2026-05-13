function [c_eff_xy, alpha_bub_xy, meta] = ...
    bubble_environment_vertical(x, y, z_curr, f_hz, c_bg, cfg)
%BUBBLE_ENVIRONMENT_VERTICAL Bubble environment hook for vertical WAPE.
% Level 0 implements a horizontally uniform empirical layer. More advanced
% spectrum and plume models are intentionally left for later priorities.

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

ny = numel(y);
nx = numel(x);
c_eff_xy = c_bg * ones(ny, nx);
alpha_bub_xy = zeros(ny, nx);

spatial_mode = 'none';
if isfield(cfg, 'bubble_spatial_mode') && ~isempty(cfg.bubble_spatial_mode)
    spatial_mode = cfg.bubble_spatial_mode;
    if isstring(spatial_mode)
        spatial_mode = char(spatial_mode);
    end
    if ~ischar(spatial_mode)
        error('cfg.bubble_spatial_mode must be a string or char.');
    end
    spatial_mode = lower(spatial_mode);
end

apply_sound_speed = false;
if isfield(cfg, 'bubble_apply_sound_speed') && ~isempty(cfg.bubble_apply_sound_speed)
    apply_sound_speed = logical(cfg.bubble_apply_sound_speed);
end
apply_attenuation = false;
if isfield(cfg, 'bubble_apply_attenuation') && ~isempty(cfg.bubble_apply_attenuation)
    apply_attenuation = logical(cfg.bubble_apply_attenuation);
end

alpha_bub = 0;
delta_c_bub = 0;
warning_flags = {};

if enabled
    switch model
        case 'off'
            enabled = false;
            spatial_mode = 'none';
        case 'level0_empirical'
            spatial_mode = '1d';
            if apply_attenuation
                alpha_bub = cfg.bubble_alpha0_np_per_m * ...
                    exp(-z_curr / cfg.bubble_layer_decay_m) * ...
                    (f_hz / cfg.bubble_f_ref_hz)^cfg.bubble_alpha_freq_exp;
            end
            if apply_sound_speed
                delta_c_bub = cfg.bubble_delta_c0_mps * ...
                    exp(-z_curr / cfg.bubble_sound_speed_decay_m) * ...
                    (f_hz / cfg.bubble_f_ref_hz)^cfg.bubble_sound_speed_freq_exp;
            end
            if ~(isscalar(alpha_bub) && isnumeric(alpha_bub) && isfinite(alpha_bub) && alpha_bub >= 0)
                error('Level 0 bubble attenuation must be finite and nonnegative.');
            end
            c_eff = c_bg + delta_c_bub;
            if ~(isscalar(c_eff) && isnumeric(c_eff) && isfinite(c_eff) && c_eff > 0)
                error('Level 0 bubble effective sound speed must be finite and positive.');
            end
            c_eff_xy = c_eff * ones(ny, nx);
            alpha_bub_xy = alpha_bub * ones(ny, nx);
        otherwise
            error('bubble_model=%s is not implemented in Priority 2.', model);
    end
end

meta.z_m = z_curr;
meta.f_hz = f_hz;
meta.c_bg_mps = c_bg;
meta.enabled = enabled;
meta.model = model;
if ~enabled
    meta.model = 'off';
end
meta.spatial_mode = spatial_mode;
meta.apply_sound_speed = apply_sound_speed;
meta.apply_attenuation = apply_attenuation;
meta.layer_decay_m = cfg.bubble_layer_decay_m;
meta.sound_speed_decay_m = cfg.bubble_sound_speed_decay_m;
meta.f_ref_hz = cfg.bubble_f_ref_hz;
meta.alpha_freq_exp = cfg.bubble_alpha_freq_exp;
meta.sound_speed_freq_exp = cfg.bubble_sound_speed_freq_exp;
meta.alpha_bub_np_per_m = alpha_bub;
meta.delta_c_bub_mps = delta_c_bub;
meta.c_eff_stats = local_field_stats(c_eff_xy);
meta.alpha_bub_stats = local_field_stats(alpha_bub_xy);
meta.warning_flags = warning_flags;

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
