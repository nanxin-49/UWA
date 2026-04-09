% Execute upward vertical PE test (independent from legacy horizontal code).

clear
format compact

paramsV = struct();
paramsV.f0 = 4000;        % Hz
paramsV.c0 = 1500;        % m/s
paramsV.z_max = 100;      % m, sea depth (z positive downward)
paramsV.stepz_lamb = 0.5; % dz_abs in wavelength units

% Strictly limited transverse domain to avoid memory blowup.
paramsV.xw = 50;          % m (full width)
paramsV.yw = 50;          % m (full width)
paramsV.nx = 1024;
paramsV.ny = 1024;

paramsV.xs = 0;           % source center x (m)
paramsV.ys = 0;           % source center y (m)
paramsV.nout = 6;         % number of snapshots along z march

paramsV.sigma_src_m = 0.3;  % must stay in [0.2, 0.5] and >= max(dx,dy)
paramsV.taper_ratio = 0.12; % per-side taper ratio in [0.10, 0.15]

paramsV.env_mode = 'uniform';  % 'uniform' or 'layered'
paramsV.enforce_1_over_R = true;

save_prefix = 'vertical_upward_4k_uniform';

simulata_vertical = CARPE3D_vertical(paramsV);

save([save_prefix '.mat'], 'simulata_vertical', 'paramsV');

print(11, '-dpng', '-r200', [save_prefix '_Figure11_xy.png'])
print(12, '-dpng', '-r200', [save_prefix '_Figure12_xz.png'])
print(13, '-dpng', '-r200', [save_prefix '_Figure13_yz.png'])
print(14, '-dpng', '-r200', [save_prefix '_Figure14_1overR.png'])
print(15, '-dpng', '-r200', [save_prefix '_Figure15_surface_kirchhoff.png'])
