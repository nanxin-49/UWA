% run_bellhop_working.m
% Working BELLHOP vertical channel simulation
% Uses flat surface for stability, demonstrates the complete workflow

clear; close all; clc;

%% PARAMETERS
freq = 6000;  % Hz
water_depth = 1500;  % m
source_depth = 1495;  % m (near bottom)
receiver_depth = 5;  % m (near surface)
receiver_range_m = 200;  % m horizontal offset
receiver_range_km = receiver_range_m / 1000;  % km

% Munk SSP parameters
c0 = 1500; z_ax = 1300; B = 1300; epsilon = 0.00737;

% Surface (for documentation)
surface_amplitude = 1.0;  % m
surface_wavelength = 50.0;  % m

% Bottom properties
bottom_cp = 1600; bottom_rho = 1.8; bottom_atten = 0.8;

% Beam parameters
num_beams_rays = 500;  % For ray tracing (fewer for clarity)
num_beams_arrivals = 5000;  % For arrivals (more for coverage)

%% SETUP
bellhop_exe = getenv('BELLHOP_EXE');
if ~isfile(bellhop_exe)
    error('BELLHOP executable not found');
end
fprintf('BELLHOP: %s\n', bellhop_exe);

output_dir = fullfile(pwd, 'bellhop_output');
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

%% GENERATE MUNK SSP
z_ssp = (0:10:water_depth)';
eta = 2 * (z_ssp - z_ax) / B;
c_z = c0 * (1 + epsilon * (eta + exp(-eta) - 1));

ssp_data.z = z_ssp;
ssp_data.c = c_z;

fprintf('Munk SSP: %d points\n', length(z_ssp));
fprintf('  c_min = %.2f m/s @ z = %.0f m\n', min(c_z), z_ssp(c_z==min(c_z)));
fprintf('  c_max = %.2f m/s @ z = 0 m\n', max(c_z));

%% GENERATE SURFACE DATA (for plotting)
r_surface = linspace(-100, 500, 500);
z_surface = surface_amplitude * sin(2*pi*r_surface / surface_wavelength);
surface_data.r = r_surface;
surface_data.z = z_surface;

%% RUN 1: RAY TRACE
fprintf('\n=== Ray Trace ===\n');
write_env_file(fullfile(output_dir, 'rays.env'), 'Ray Trace', ...
    freq, z_ssp, c_z, water_depth, source_depth, receiver_depth, ...
    receiver_range_km, bottom_cp, bottom_rho, bottom_atten, ...
    'R', num_beams_rays, -89, 89, receiver_range_km*2.5);

run_bellhop(bellhop_exe, output_dir, 'rays');

%% RUN 2: EIGENRAYS
fprintf('\n=== Eigenrays ===\n');
write_env_file(fullfile(output_dir, 'eigen.env'), 'Eigenrays', ...
    freq, z_ssp, c_z, water_depth, source_depth, receiver_depth, ...
    receiver_range_km, bottom_cp, bottom_rho, bottom_atten, ...
    'E', num_beams_arrivals, -89, 89, receiver_range_km*2.5);

run_bellhop(bellhop_exe, output_dir, 'eigen');

%% RUN 3: ARRIVALS
fprintf('\n=== Arrivals ===\n');
write_env_file(fullfile(output_dir, 'arrivals.env'), 'Arrivals', ...
    freq, z_ssp, c_z, water_depth, source_depth, receiver_depth, ...
    receiver_range_km, bottom_cp, bottom_rho, bottom_atten, ...
    'A', num_beams_arrivals, -89, 89, receiver_range_km*2.5);

run_bellhop(bellhop_exe, output_dir, 'arrivals');

%% PLOTTING
fprintf('\n=== Generating Figures ===\n');

% Figure 1: Munk SSP
fig1 = figure('Position', [100 100 600 800]);
plot(c_z, z_ssp, 'b-', 'LineWidth', 2);
set(gca, 'YDir', 'reverse');
grid on;
xlabel('Sound Speed (m/s)', 'FontSize', 12);
ylabel('Depth (m)', 'FontSize', 12);
title('Munk Deep-Ocean Sound Speed Profile', 'FontSize', 14, 'FontWeight', 'bold');
ylim([0 water_depth]);
saveas(fig1, fullfile(output_dir, 'fig1_munk_ssp.png'));
savefig(fig1, fullfile(output_dir, 'fig1_munk_ssp.fig'));
fprintf('Saved: fig1_munk_ssp\n');

% Figure 2: Surface profile
fig2 = figure('Position', [150 150 900 400]);
plot(r_surface, z_surface, 'b-', 'LineWidth', 2);
hold on;
plot(xlim, [0 0], 'k--', 'LineWidth', 1);
xline(0, 'r-', 'Tx', 'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom');
xline(receiver_range_m, 'g-', 'Rx', 'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom');
grid on;
xlabel('Range (m)', 'FontSize', 12);
ylabel('Surface Elevation (m)', 'FontSize', 12);
title(sprintf('Sinusoidal Surface (A=%.1fm, \\lambda=%.1fm) - For Reference', ...
    surface_amplitude, surface_wavelength), 'FontSize', 13, 'FontWeight', 'bold');
xlim([-100 500]);
text(250, surface_amplitude*0.5, 'Note: Simulation uses flat surface', ...
    'FontSize', 10, 'Color', 'r', 'FontWeight', 'bold');
saveas(fig2, fullfile(output_dir, 'fig2_surface.png'));
savefig(fig2, fullfile(output_dir, 'fig2_surface.fig'));
fprintf('Saved: fig2_surface\n');

% Figure 3: Ray trace
fig3 = plot_bellhop_rays(fullfile(output_dir, 'rays.ray'), ...
    source_depth, receiver_depth, receiver_range_m, water_depth);
if ~isempty(fig3)
    saveas(fig3, fullfile(output_dir, 'fig3_rays.png'));
    savefig(fig3, fullfile(output_dir, 'fig3_rays.fig'));
    fprintf('Saved: fig3_rays\n');
end

% Figure 4: Eigenrays
fig4 = plot_bellhop_rays(fullfile(output_dir, 'eigen.ray'), ...
    source_depth, receiver_depth, receiver_range_m, water_depth);
if ~isempty(fig4)
    saveas(fig4, fullfile(output_dir, 'fig4_eigenrays.png'));
    savefig(fig4, fullfile(output_dir, 'fig4_eigenrays.fig'));
    fprintf('Saved: fig4_eigenrays\n');
end

% Figure 5: Arrivals / Channel Impulse Response
[arrivals, summary] = read_bellhop_arrivals(fullfile(output_dir, 'arrivals.arr'));
if ~isempty(arrivals) && isfield(arrivals, 'delay')
    fig5 = figure('Position', [200 200 1000 600]);
    delays_rel = (arrivals.delay - min(arrivals.delay)) * 1000;  % ms
    stem(delays_rel, abs(arrivals.amplitude), 'b-', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
    grid on;
    xlabel('Relative Delay (ms)', 'FontSize', 12);
    ylabel('Amplitude', 'FontSize', 12);
    title('Channel Impulse Response', 'FontSize', 14, 'FontWeight', 'bold');
    text(0.98, 0.98, sprintf('%d arrivals', summary.num_arrivals), ...
        'Units', 'normalized', 'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'right', 'FontSize', 11, 'BackgroundColor', 'white');
    saveas(fig5, fullfile(output_dir, 'fig5_impulse_response.png'));
    savefig(fig5, fullfile(output_dir, 'fig5_impulse_response.fig'));
    fprintf('Saved: fig5_impulse_response\n');

    % Print summary
    fprintf('\n=== SUMMARY ===\n');
    fprintf('Frequency: %d Hz\n', freq);
    fprintf('Geometry: Tx @ %.0fm depth, Rx @ %.0fm depth, %.0fm range\n', ...
        source_depth, receiver_depth, receiver_range_m);
    fprintf('Arrivals: %d eigenrays found\n', summary.num_arrivals);
    if summary.num_arrivals > 0
        fprintf('First arrival: %.6f s (%.3f ms)\n', summary.first_delay, summary.first_delay*1000);
        fprintf('Delay spread: %.6f s (%.3f ms)\n', summary.delay_spread, summary.delay_spread*1000);
        fprintf('Strongest: amplitude=%.6f @ delay=%.6f s\n', summary.max_amplitude, summary.max_delay);

        % Sanity check
        direct_distance = sqrt(receiver_range_m^2 + (source_depth - receiver_depth)^2);
        c_mean = mean(c_z);
        direct_time_est = direct_distance / c_mean;
        fprintf('\nSanity check:\n');
        fprintf('  Direct path distance: %.1f m\n', direct_distance);
        fprintf('  Estimated travel time: %.6f s (c_mean=%.1f m/s)\n', direct_time_est, c_mean);
        fprintf('  Actual first arrival: %.6f s\n', summary.first_delay);
        fprintf('  Difference: %.3f ms\n', (summary.first_delay - direct_time_est)*1000);
    end

    % Save data
    save(fullfile(output_dir, 'simulation_results.mat'), ...
        'freq', 'source_depth', 'receiver_depth', 'receiver_range_m', ...
        'ssp_data', 'surface_data', 'arrivals', 'summary');
    fprintf('\nSaved: simulation_results.mat\n');
else
    fprintf('WARNING: No arrivals data available\n');
end

fprintf('\n=== COMPLETE ===\n');
fprintf('Output directory: %s\n', output_dir);
