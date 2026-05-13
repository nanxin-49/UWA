function [c_eff, alpha_bub, em_meta] = ...
    bubble_effective_medium(a_grid_m, n_a, z_curr, f_hz, c_bg, cfg)
%BUBBLE_EFFECTIVE_MEDIUM Convert 1D bubble spectrum to medium properties.

if nargin < 6
    error('bubble_effective_medium requires a_grid_m, n_a, z_curr, f_hz, c_bg, and cfg.');
end
if ~(isnumeric(a_grid_m) && isvector(a_grid_m) && numel(a_grid_m) >= 2 && ...
     all(isfinite(a_grid_m(:))) && all(a_grid_m(:) > 0) && all(diff(a_grid_m(:)) > 0))
    error('a_grid_m must be a positive finite increasing vector.');
end
if ~(isnumeric(n_a) && isvector(n_a) && numel(n_a) == numel(a_grid_m) && ...
     all(isfinite(n_a(:))) && all(n_a(:) >= 0))
    error('n_a must be a finite nonnegative vector matching a_grid_m.');
end
if ~(isscalar(z_curr) && isnumeric(z_curr) && isfinite(z_curr) && z_curr >= 0)
    error('z_curr must be a finite nonnegative scalar depth.');
end
if ~(isscalar(f_hz) && isnumeric(f_hz) && isfinite(f_hz) && f_hz > 0)
    error('f_hz must be a positive finite scalar frequency.');
end
if ~(isscalar(c_bg) && isnumeric(c_bg) && isfinite(c_bg) && c_bg > 0)
    error('c_bg must be a positive finite scalar sound speed.');
end

a_grid_m = a_grid_m(:).';
n_a = n_a(:).';
omega = 2*pi*f_hz;
warning_flags = {};

P0 = cfg.bubble_P_atm_pa + cfg.bubble_rho_w_kg_m3 * cfg.bubble_g_m_s2 * z_curr;
a_res = (1 / (2*pi*f_hz)) * sqrt(3 * cfg.bubble_gamma * P0 / cfg.bubble_rho_w_kg_m3);
if ~(isscalar(a_res) && isfinite(a_res) && a_res > 0)
    error('Computed resonance radius must be positive and finite.');
end

if all(n_a == 0)
    c_eff = c_bg;
    alpha_bub = 0;
    integral_term = complex(0, 0);
    q = 1 / c_bg;
else
    d = cfg.bubble_delta_const;
    if ~(isscalar(d) && isnumeric(d) && isfinite(d) && d > 0)
        error('bubble_delta_const must be positive and finite.');
    end
    % The WAPE screen uses exp(-alpha*ds) for positive alpha. With the
    % exp(-i*omega*t) propagation convention used here, passive damping is
    % represented by -1i*d so omega*imag(q) is positive in the normal path.
    integrand = (a_grid_m .* n_a) ./ (((a_res ./ a_grid_m).^2) - 1 - 1i*d);
    integral_term = trapz(a_grid_m, integrand);
    inv_c_tilde_sq = (1 / c_bg^2) + (1 / (pi * f_hz^2)) * integral_term;
    q = sqrt(inv_c_tilde_sq);
    if real(q) < 0
        q = -q;
        warning_flags{end+1} = 'complex_slowness_sign_flipped';
    end
    c_eff = 1 / real(q);
    alpha_bub = omega * imag(q);
    if alpha_bub < 0
        alpha_bub = max(alpha_bub, 0);
        warning_flags{end+1} = 'negative_attenuation_clipped';
    end
end

if ~(isscalar(c_eff) && isnumeric(c_eff) && isfinite(c_eff) && isreal(c_eff) && c_eff > 0)
    error('Effective sound speed must be finite, real, and positive.');
end
if ~(isscalar(alpha_bub) && isnumeric(alpha_bub) && isfinite(alpha_bub) && isreal(alpha_bub) && alpha_bub >= 0)
    error('Bubble attenuation must be finite, real, and nonnegative.');
end

em_meta = struct();
em_meta.model = 'complex_sound_speed';
em_meta.damping_model = cfg.bubble_damping_model;
em_meta.delta_const = cfg.bubble_delta_const;
em_meta.denominator_damping_sign = '-1i*d';
em_meta.alpha_extraction = 'omega*imag(q)';
em_meta.P0_pa = P0;
em_meta.resonance_radius_m = a_res;
em_meta.integral_term = integral_term;
em_meta.complex_slowness = q;
em_meta.c_eff_mps = c_eff;
em_meta.alpha_bub_np_per_m = alpha_bub;
em_meta.warning_flags = warning_flags;

end
