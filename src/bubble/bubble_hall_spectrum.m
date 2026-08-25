function [n_a, beta, spec_meta] = bubble_hall_spectrum(a_grid_m, z_curr, U10, cfg)
%BUBBLE_HALL_SPECTRUM Hall-type 1D bubble number density spectrum.

if nargin < 4
    error('bubble_hall_spectrum requires a_grid_m, z_curr, U10, and cfg.');
end
if ~(isnumeric(a_grid_m) && isvector(a_grid_m) && numel(a_grid_m) >= 2 && ...
     all(isfinite(a_grid_m(:))) && all(a_grid_m(:) > 0) && all(diff(a_grid_m(:)) > 0))
    error('a_grid_m must be a positive finite increasing vector.');
end
if ~(isscalar(z_curr) && isnumeric(z_curr) && isfinite(z_curr) && z_curr >= 0)
    error('z_curr must be a finite nonnegative scalar depth.');
end
if ~(isscalar(U10) && isnumeric(U10) && isfinite(U10) && U10 > 0)
    error('U10 must be a positive finite scalar wind speed.');
end

a_grid_m = a_grid_m(:).';
p0 = 1.6e10;
if isfield(cfg, 'bubble_p0_m4') && ~isempty(cfg.bubble_p0_m4)
    p0 = cfg.bubble_p0_m4;
end
if ~(isscalar(p0) && isnumeric(p0) && isfinite(p0) && p0 > 0)
    error('bubble_p0_m4 must be a positive finite scalar.');
end

if U10 <= 7.5
    layer_depth_m = 0.4;
else
    layer_depth_m = 0.4 + 0.115 * (U10 - 7.5);
end
depth_factor = exp(-z_curr / layer_depth_m);

a_min_m = 10e-6;
a_max_m = 1000e-6;
a_ref_m = (54.4 + 1.984e-6 * z_curr) * 1e-6;
chi = 4.37 + (z_curr / 2.55)^2;

G = zeros(size(a_grid_m));
mask_low = (a_grid_m >= a_min_m) & (a_grid_m <= a_ref_m);
mask_high = (a_grid_m > a_ref_m) & (a_grid_m <= a_max_m);
G(mask_low) = (a_ref_m ./ a_grid_m(mask_low)).^4;
G(mask_high) = (a_ref_m ./ a_grid_m(mask_high)).^chi;

n_a = p0 .* depth_factor .* G .* (U10 / 13)^3 .* cfg.bubble_strength_scale;
if any(~isfinite(n_a(:))) || any(n_a(:) < 0)
    error('Hall bubble spectrum produced invalid number density.');
end

beta_raw = trapz(a_grid_m, (4*pi/3) .* a_grid_m.^3 .* n_a);
if ~(isscalar(beta_raw) && isfinite(beta_raw) && beta_raw >= 0)
    error('Hall bubble spectrum produced invalid void fraction.');
end

scaled_to_beta_max = false;
scale_factor = 1;
beta = beta_raw;
if beta_raw > cfg.bubble_beta_max
    scale_factor = cfg.bubble_beta_max / beta_raw;
    n_a = n_a .* scale_factor;
    beta = cfg.bubble_beta_max;
    scaled_to_beta_max = true;
end

spec_meta = struct();
spec_meta.model = 'hall1d';
spec_meta.p0_m4 = p0;
spec_meta.U10_mps = U10;
spec_meta.layer_depth_m = layer_depth_m;
spec_meta.depth_factor = depth_factor;
spec_meta.a_ref_m = a_ref_m;
spec_meta.chi = chi;
spec_meta.beta_raw = beta_raw;
spec_meta.beta = beta;
spec_meta.beta_max = cfg.bubble_beta_max;
spec_meta.scaled_to_beta_max = scaled_to_beta_max;
spec_meta.scale_factor = scale_factor;
spec_meta.n_a_stats = local_field_stats(n_a);
spec_meta.radius_grid_m = a_grid_m;
spec_meta.warning_flags = {};
if scaled_to_beta_max
    spec_meta.warning_flags{end+1} = 'beta_clipped_to_max';
end

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
