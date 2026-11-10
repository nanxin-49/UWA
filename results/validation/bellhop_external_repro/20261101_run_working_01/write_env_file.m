function write_env_file(filename, title, freq, z_ssp, c_z, water_depth, ...
    source_depth, receiver_depth, receiver_range_km, ...
    bottom_cp, bottom_rho, bottom_atten, ...
    run_type, num_beams, angle_min, angle_max, rmax_km)
% write_env_file - Write BELLHOP environment file
%
% Inputs:
%   filename - full path to .env file
%   title - case title
%   freq - frequency (Hz)
%   z_ssp, c_z - sound speed profile vectors (m, m/s)
%   water_depth - water depth (m)
%   source_depth, receiver_depth - depths (m)
%   receiver_range_km - receiver range (km)
%   bottom_cp, bottom_rho, bottom_atten - bottom properties
%   run_type - 'R' (rays), 'E' (eigenrays), 'A' (arrivals), 'C' (coherent TL)
%   num_beams - number of beams
%   angle_min, angle_max - beam angles (degrees)
%   rmax_km - maximum range for computation (km)

fid = fopen(filename, 'w');

% Title and frequency
fprintf(fid, '''%s''\n', title);
fprintf(fid, '%.1f\n', freq);

% Number of media
fprintf(fid, '1\n');

% SSP interpolation option (C=C-linear, V=vacuum compatible, W=water)
fprintf(fid, '''CVW''\n');

% SSP data
fprintf(fid, '%d 0.0 %.1f\n', length(z_ssp), water_depth);
for i = 1:length(z_ssp)
    fprintf(fid, '%.2f %.4f /\n', z_ssp(i), c_z(i));
end

% Top boundary (vacuum / pressure release)
fprintf(fid, '''V'' 0.0\n');

% Bottom boundary (acousto-elastic half-space)
fprintf(fid, '''A'' %.1f\n', water_depth);
fprintf(fid, '%.1f %.4f %.2f %.2f 0.0 0.0 /\n', ...
    water_depth, bottom_cp, bottom_rho, bottom_atten);

% Source depth
fprintf(fid, '1\n');
fprintf(fid, '%.2f /\n', source_depth);

% Receiver depth
fprintf(fid, '1\n');
fprintf(fid, '%.2f /\n', receiver_depth);

% Receiver range (km)
fprintf(fid, '1\n');
fprintf(fid, '%.6f /\n', receiver_range_km);

% Run type
fprintf(fid, '''%s''\n', run_type);

% Number of beams
fprintf(fid, '%d\n', num_beams);

% Beam angles
fprintf(fid, '%.2f %.2f /\n', angle_min, angle_max);

% Box limits (z_min z_max range_max in km)
fprintf(fid, '0.0 %.1f %.6f\n', water_depth + 5, rmax_km);

fclose(fid);
end
