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

paramsV.x_tx = 0;         % tx center x (m)
paramsV.y_tx = 0;         % tx center y (m)
paramsV.z_tx = 100;       % tx depth (m)
paramsV.x_rx = 0;         % rx center x (m)
paramsV.y_rx = 0;         % rx center y (m)
paramsV.z_rx = 4;         % rx depth under buoy (m), typically 3-5 m
paramsV.nout = 6;         % number of snapshots along z march

paramsV.sigma_src_m = 0.3;  % must stay in [0.2, 0.5] and >= max(dx,dy)
paramsV.sponge_ratio = 0.12; % per-side sponge ratio in [0.10, 0.15]
paramsV.alpha_max_np_per_m = 0.15;

paramsV.env_mode = 'uniform';  % 'uniform' or 'layered'
paramsV.enforce_1_over_R = true;
paramsV.enable_surface_reflection = true;
paramsV.save_mode = 'slice';
paramsV.use_gpu = false;
paramsV.sea_wind_speed = 5.0;
paramsV.sea_hs_target = 0.5;
paramsV.sea_seed = 12345;

save_prefix = 'vertical_upward_4k_uniform';

simulata_vertical = CARPE3D_vertical(paramsV);

save([save_prefix '.mat'], 'simulata_vertical', 'paramsV');

print(11, '-dpng', '-r200', [save_prefix '_Figure11_xy.png'])
print(12, '-dpng', '-r200', [save_prefix '_Figure12_xz.png'])
print(13, '-dpng', '-r200', [save_prefix '_Figure13_yz.png'])
print(14, '-dpng', '-r200', [save_prefix '_Figure14_1overR.png'])
print(15, '-dpng', '-r200', [save_prefix '_Figure15_surface_kirchhoff.png'])
