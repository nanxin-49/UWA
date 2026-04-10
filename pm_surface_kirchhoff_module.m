function [surface_elevation, delta_phi, psi_ref, meta] = ...
    pm_surface_kirchhoff_module(psi_inc, KX, KY, x, y, xw, yw, lambda0, pm_cfg)
%PM_SURFACE_KIRCHHOFF_MODULE
% Task 1: 2D Pierson-Moskowitz rough surface synthesis on (x,y) grid.
% Task 2: Kirchhoff phase distortion on incident field psi_inc.
% Task 3: Sanity-check visualization.

if nargin < 9 || isempty(pm_cfg)
    pm_cfg = struct();
end
if ~isfield(pm_cfg, 'U') || isempty(pm_cfg.U)
    pm_cfg.U = 5.0;
end
if ~isfield(pm_cfg, 'Hs_target') || isempty(pm_cfg.Hs_target)
    pm_cfg.Hs_target = 0.5;
end
if ~isfield(pm_cfg, 'seed') || isempty(pm_cfg.seed)
    pm_cfg.seed = 12345;
end
if ~isfield(pm_cfg, 'show_figure') || isempty(pm_cfg.show_figure)
    pm_cfg.show_figure = true;
end
if ~isfield(pm_cfg, 'reflect_coeff') || isempty(pm_cfg.reflect_coeff)
    pm_cfg.reflect_coeff = -1;
end
if ~isfield(pm_cfg, 'phase_mode') || isempty(pm_cfg.phase_mode)
    pm_cfg.phase_mode = 'normal';
end
if ~isfield(pm_cfg, 'oblique_clip') || isempty(pm_cfg.oblique_clip)
    pm_cfg.oblique_clip = [0, 1];
end
if ~isfield(pm_cfg, 'tx_xyz') || isempty(pm_cfg.tx_xyz)
    pm_cfg.tx_xyz = [0, 0, NaN];
end
if ~isfield(pm_cfg, 'rx_xyz') || isempty(pm_cfg.rx_xyz)
    pm_cfg.rx_xyz = [0, 0, NaN];
end
if ~isfield(pm_cfg, 'z_surface') || isempty(pm_cfg.z_surface)
    pm_cfg.z_surface = 0;
end

if ~(isscalar(pm_cfg.reflect_coeff) && isnumeric(pm_cfg.reflect_coeff) && isfinite(pm_cfg.reflect_coeff))
    error('pm_cfg.reflect_coeff must be a finite scalar (real or complex).');
end
if isstring(pm_cfg.phase_mode)
    pm_cfg.phase_mode = char(pm_cfg.phase_mode);
end
phase_mode = lower(pm_cfg.phase_mode);
if ~ischar(phase_mode) || (~strcmp(phase_mode, 'normal') && ~strcmp(phase_mode, 'oblique'))
    error('pm_cfg.phase_mode must be ''normal'' or ''oblique''.');
end
if ~(isnumeric(pm_cfg.oblique_clip) && numel(pm_cfg.oblique_clip) == 2 && all(isfinite(pm_cfg.oblique_clip(:))))
    error('pm_cfg.oblique_clip must be a finite 1x2 numeric range.');
end
clip_pair = sort(pm_cfg.oblique_clip(:).');
if clip_pair(1) < 0 || clip_pair(2) > 1
    error('pm_cfg.oblique_clip must stay within [0,1].');
end
if ~(isnumeric(pm_cfg.tx_xyz) && numel(pm_cfg.tx_xyz) == 3 && all(isfinite(pm_cfg.tx_xyz(:))))
    error('pm_cfg.tx_xyz must be a finite [x,y,z].');
end
if ~(isnumeric(pm_cfg.rx_xyz) && numel(pm_cfg.rx_xyz) == 3 && all(isfinite(pm_cfg.rx_xyz(:))))
    error('pm_cfg.rx_xyz must be a finite [x,y,z].');
end

% -----------------------------
% 1) PM rough surface generation
% -----------------------------
U = pm_cfg.U;          % m/s
g = 9.81;              % m/s^2
alpha_PM = 8.10e-3;
beta_PM = 0.74;

K = sqrt(KX.^2 + KY.^2);
E1D = zeros(size(K));
Phi2D = zeros(size(K));

mask = (K > 0);
K_nonzero = K(mask);
E1D(mask) = (alpha_PM ./ (2 .* K_nonzero.^3)) .* ...
            exp(-beta_PM * (g^2) ./ (U^4 .* K_nonzero.^2));
Phi2D(mask) = E1D(mask) ./ (2*pi*K_nonzero);

dkx = 2*pi / xw;
dky = 2*pi / yw;

A = sqrt(Phi2D .* dkx .* dky);
rng(pm_cfg.seed, 'twister')
N = (randn(size(K)) + 1i*randn(size(K))) / sqrt(2);
Zk = A .* N;

eta_raw = real(ifft2(Zk));
eta_raw = eta_raw * numel(eta_raw); % compensate MATLAB ifft2 normalization

% Calibrate to realistic PM sea-state scale (~0.5 m significant wave height at U=5 m/s).
Hs_target = pm_cfg.Hs_target; % m
Hs_raw = 4 * std(eta_raw(:));
if Hs_raw > 0
    scale_factor = Hs_target / Hs_raw;
else
    scale_factor = 0;
end
surface_elevation = eta_raw * scale_factor;

% -----------------------------
% 2) Kirchhoff phase distortion
% -----------------------------
k0 = 2*pi / lambda0;
[X, Y] = meshgrid(x, y);

switch phase_mode
    case 'normal'
        phase_factor = 2 * ones(size(surface_elevation));
    case 'oblique'
        tx_xyz = pm_cfg.tx_xyz(:).';
        rx_xyz = pm_cfg.rx_xyz(:).';
        z_surface = pm_cfg.z_surface;

        dz_i = tx_xyz(3) - z_surface;
        dz_r = rx_xyz(3) - z_surface;
        if dz_i <= 0 || dz_r < 0
            error('For oblique mode, tx/rx must be at or below the surface: z >= z_surface.');
        end

        dx_i = X - tx_xyz(1);
        dy_i = Y - tx_xyz(2);
        Ri = sqrt(dx_i.^2 + dy_i.^2 + dz_i.^2);
        cos_i = dz_i ./ max(Ri, eps);

        dx_r = X - rx_xyz(1);
        dy_r = Y - rx_xyz(2);
        Rr = sqrt(dx_r.^2 + dy_r.^2 + dz_r.^2);
        cos_r = dz_r ./ max(Rr, eps);

        cos_i = min(max(cos_i, clip_pair(1)), clip_pair(2));
        cos_r = min(max(cos_r, clip_pair(1)), clip_pair(2));
        phase_factor = cos_i + cos_r;
    otherwise
        error('Unsupported phase_mode: %s', phase_mode);
end

delta_phi = 2 * k0 * surface_elevation .* phase_factor;
psi_ref = pm_cfg.reflect_coeff .* psi_inc .* exp(1i * delta_phi);

% -----------------------------
% 3) Sanity-check visualization
% -----------------------------
if pm_cfg.show_figure
    figure(15); clf
    subplot(1,2,1)
    surf(X, Y, surface_elevation, 'EdgeColor', 'none')
    view(40, 35)
    axis tight
    xlabel('x (m)')
    ylabel('y (m)')
    zlabel('\xi(x,y) (m)')
    title('PM rough sea surface elevation')
    colorbar

    subplot(1,2,2)
    imagesc(x, y, angle(psi_ref))
    axis xy tight
    xlabel('x (m)')
    ylabel('y (m)')
    title(sprintf('Phase of reflected field angle(\\psi_{ref}), mode=%s', phase_mode))
    colorbar
end

E_from_Phi = zeros(size(K));
E_from_Phi(mask) = Phi2D(mask) .* 2*pi .* K(mask);
rel_err = abs(E_from_Phi(mask) - E1D(mask)) ./ max(abs(E1D(mask)), eps);

meta = struct();
meta.U = U;
meta.g = g;
meta.alpha_PM = alpha_PM;
meta.beta_PM = beta_PM;
meta.dkx = dkx;
meta.dky = dky;
meta.spectrum_definition = 'E1D(K)=alpha/(2K^3)exp(-beta g^2/(U^4 K^2)); Phi2D(Kx,Ky)=E1D(K)/(2*pi*K)';
meta.E1D = local_field_stats(E1D(mask));
meta.Phi2D = local_field_stats(Phi2D(mask));
meta.phi2d_transform_error = struct( ...
    'max_rel', max(rel_err(:), [], 'omitnan'), ...
    'mean_rel', mean(rel_err(:), 'omitnan'));
meta.Hs_target = Hs_target;
meta.Hs_raw = Hs_raw;
meta.Hs_scaled = 4 * std(surface_elevation(:));
meta.scale_factor = scale_factor;
meta.seed = pm_cfg.seed;
meta.reflection_coeff_used = pm_cfg.reflect_coeff;
meta.phase_mode_used = phase_mode;
meta.phase_factor_stats = struct( ...
    'min', min(phase_factor(:)), ...
    'max', max(phase_factor(:)), ...
    'mean', mean(phase_factor(:)));
meta.oblique_clip_used = clip_pair;
meta.show_figure = logical(pm_cfg.show_figure);
meta.enabled = true;

end

function stats = local_field_stats(v)
if isempty(v)
    stats = struct('min', NaN, 'max', NaN, 'mean', NaN, 'std', NaN);
    return
end
stats = struct( ...
    'min', min(v(:)), ...
    'max', max(v(:)), ...
    'mean', mean(v(:)), ...
    'std', std(v(:)));
end
