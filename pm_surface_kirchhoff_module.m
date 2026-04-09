function [surface_elevation, delta_phi, psi_ref, meta] = ...
    pm_surface_kirchhoff_module(psi_inc, KX, KY, x, y, xw, yw, lambda0)
%PM_SURFACE_KIRCHHOFF_MODULE
% Task 1: 2D Pierson-Moskowitz rough surface synthesis on (x,y) grid.
% Task 2: Kirchhoff phase distortion on incident field psi_inc.
% Task 3: Sanity-check visualization.

% -----------------------------
% 1) PM rough surface generation
% -----------------------------
U = 5.0;               % m/s
g = 9.81;              % m/s^2
alpha_PM = 8.10e-3;
beta_PM = 0.74;

K = sqrt(KX.^2 + KY.^2);
W = zeros(size(K));

mask = (K > 0);
K_nonzero = K(mask);
W(mask) = (alpha_PM ./ (2 .* K_nonzero.^3)) .* ...
          exp(-beta_PM * (g^2) ./ (U^4 .* K_nonzero.^2));

dkx = 2*pi / xw;
dky = 2*pi / yw;

A = sqrt(W .* dkx .* dky);
N = (randn(size(K)) + 1i*randn(size(K))) / sqrt(2);
Zk = A .* N;

eta_raw = real(ifft2(Zk));
eta_raw = eta_raw * numel(eta_raw); % compensate MATLAB ifft2 normalization

% Calibrate to realistic PM sea-state scale (~0.5 m significant wave height at U=5 m/s).
Hs_target = 0.5; % m
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
delta_phi = 2*pi * (2 * surface_elevation / lambda0);
psi_ref = psi_inc .* exp(1i * delta_phi);

% -----------------------------
% 3) Sanity-check visualization
% -----------------------------
[X, Y] = meshgrid(x, y);

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
title('Phase of reflected field angle(\psi_{ref})')
colorbar

meta = struct();
meta.U = U;
meta.g = g;
meta.alpha_PM = alpha_PM;
meta.beta_PM = beta_PM;
meta.dkx = dkx;
meta.dky = dky;
meta.Hs_target = Hs_target;
meta.Hs_raw = Hs_raw;
meta.Hs_scaled = 4 * std(surface_elevation(:));
meta.scale_factor = scale_factor;

end

