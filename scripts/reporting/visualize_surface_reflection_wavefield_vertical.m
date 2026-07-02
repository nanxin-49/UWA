run(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'bootstrap_project.m'));
% Visualize incident and surface-reflected wavefields at the reference frequency.
% The static plots show PE complex-envelope magnitude/phase. The animation
% restores a nominal-c0 carrier for visualization; it is not a time-domain solve.

clear
format compact

model_name = getenv('SURFACE_WAVEFIELD_MODEL');
if isempty(model_name)
    model_name = 'ssa1_geometry';
end
grid_n = str2double(getenv('SURFACE_WAVEFIELD_GRID_N'));
if ~isfinite(grid_n)
    grid_n = 128;
end
Hs_target = str2double(getenv('SURFACE_WAVEFIELD_HS'));
if ~isfinite(Hs_target)
    Hs_target = 0.2;
end

paramsV = local_base_params(round(grid_n), Hs_target, model_name);
channel = vertical_channel_model(paramsV);
meta = channel.surface_wavefield_meta;
if ~meta.enabled
    error('Surface wavefield diagnostics were not generated.');
end

local_plot_xz(meta, paramsV, 'surface_wavefield_xz_comparison.png');
local_plot_xy(meta, 'surface_wavefield_xy_comparison.png');
spectrum_summary = local_plot_spectra(meta, 'surface_wavefield_spectrum_comparison.png');
animation_written = local_write_animation( ...
    meta, paramsV, 'surface_wavefield_instantaneous_pressure.mp4');

visualization_summary = struct( ...
    'surface_boundary_model', paramsV.surface_boundary_model, ...
    'surface_ssa_kernel_mode', paramsV.surface_ssa_kernel_mode, ...
    'Hs_target_m', paramsV.sea_hs_target, ...
    'reference_frequency_hz', meta.reference_frequency_hz, ...
    'slice_axis', meta.slice_axis, ...
    'endpoint_consistency', meta.endpoint_consistency, ...
    'spectrum_summary', spectrum_summary, ...
    'animation_written', animation_written, ...
    'carrier_reconstruction', meta.carrier_reconstruction);
save('surface_wavefield_visualization_result.mat', 'visualization_summary', 'paramsV');

fprintf('Generated surface wavefield visualization for %s / %s.\n', ...
    paramsV.surface_boundary_model, paramsV.surface_ssa_kernel_mode);
fprintf('Max endpoint consistency error: %.3e\n', max(struct2array(meta.endpoint_consistency)));

function paramsV = local_base_params(grid_n, Hs_target, model_name)
paramsV = struct();
paramsV.f0 = 6000;
paramsV.c0 = 1500;
paramsV.z_max = 30;
paramsV.z_tx = 30;
paramsV.z_rx = 3;
paramsV.xw = 30;
paramsV.yw = 30;
paramsV.nx = grid_n;
paramsV.ny = grid_n;
paramsV.x_tx = 0;
paramsV.y_tx = 0;
paramsV.x_rx = 0;
paramsV.y_rx = 0;
paramsV.sigma_src_m = 0.3;
paramsV.stepz_lamb = 0.75;
paramsV.sponge_ratio = 0.12;
paramsV.alpha_max_np_per_m = 0.15;
paramsV.env_mode = 'uniform';
paramsV.enable_surface_reflection = true;
paramsV.surface_reflect_coeff = -1;
paramsV.surface_phase_mode = 'normal';
paramsV.sea_wind_speed = 5;
paramsV.sea_hs_target = Hs_target;
paramsV.sea_seed = 12345;
paramsV.show_figures = false;
paramsV.enforce_1_over_R = false;
paramsV.save_mode = 'rx_only';
paramsV.use_gpu = false;
paramsV.surface_wavefield_diagnostics = true;
paramsV.surface_wavefield_slice_axis = 'x';
paramsV.surface_wavefield_max_z_samples = 192;
paramsV.surface_boundary_redistribution_diagnostics = true;

if any(strcmpi(model_name, {'kirchhoff_spatial', 'kirchhoff_kdomain'}))
    paramsV.surface_boundary_model = lower(model_name);
else
    paramsV.surface_boundary_model = 'ssa_stat_kernel';
    paramsV.surface_ssa_kernel_mode = lower(model_name);
    paramsV.surface_ssa_random_scatter = true;
    paramsV.surface_ssa_conv_padding = 'zero_padded';
    paramsV.surface_ssa_geometry_source_id = ...
        'SSA.md; first-order pressure-release / Dirichlet geometry';
end
end

function local_plot_xz(meta, paramsV, file_name)
ref_amp = max(meta.normalization_amplitude, eps);
inc_db = 20*log10(max(abs(meta.incident_field_slice), eps) / ref_amp);
ref_db = 20*log10(max(abs(meta.reflected_field_slice), eps) / ref_amp);
axis_name = meta.slice_axis;

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1200, 480]);
tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile
imagesc(meta.transverse_coordinate_m, meta.incident_z_m, inc_db.');
set(gca, 'YDir', 'reverse');
axis tight
clim([-50, 0]);
colorbar
xlabel(sprintf('%s (m)', axis_name));
ylabel('z (m, positive downward)');
title('Incident field envelope, |\Psi_{inc}| (dB)', 'Interpreter', 'tex');

nexttile
imagesc(meta.transverse_coordinate_m, meta.reflected_z_m, ref_db.');
set(gca, 'YDir', 'reverse');
axis tight
clim([-50, 0]);
colorbar
xlabel(sprintf('%s (m)', axis_name));
ylabel('z (m, positive downward)');
title('Surface-reflected field envelope, |\Psi_{ref}| (dB)', 'Interpreter', 'tex');
sgtitle(sprintf('%s, H_s=%.3g m, f=%.0f Hz; common amplitude reference', ...
    local_model_label(meta), paramsV.sea_hs_target, meta.reference_frequency_hz), ...
    'Interpreter', 'none');
exportgraphics(fig, file_name, 'Resolution', 180);
close(fig);
end

function local_plot_xy(meta, file_name)
ref_amp = max(max(abs(meta.surface_incident_xy(:))), eps);
inc_db = 20*log10(max(abs(meta.surface_incident_xy), eps) / ref_amp);
ref_db = 20*log10(max(abs(meta.surface_reflected_xy), eps) / ref_amp);

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80, 80, 1200, 850]);
tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile
imagesc(meta.x_m, meta.y_m, inc_db); axis image xy; clim([-50, 0]); colorbar
xlabel('x (m)'); ylabel('y (m)'); title('Surface incident magnitude (dB)');
nexttile
imagesc(meta.x_m, meta.y_m, ref_db); axis image xy; clim([-50, 0]); colorbar
xlabel('x (m)'); ylabel('y (m)'); title('Surface reflected magnitude (dB)');
nexttile
imagesc(meta.x_m, meta.y_m, angle(meta.surface_incident_xy)); axis image xy
clim([-pi, pi]); colorbar
xlabel('x (m)'); ylabel('y (m)'); title('Surface incident phase (rad)');
nexttile
imagesc(meta.x_m, meta.y_m, angle(meta.surface_reflected_xy)); axis image xy
clim([-pi, pi]); colorbar
xlabel('x (m)'); ylabel('y (m)'); title('Surface reflected phase (rad)');
sgtitle(sprintf('Sea-surface complex envelopes: %s', local_model_label(meta)), ...
    'Interpreter', 'none');
exportgraphics(fig, file_name, 'Resolution', 180);
close(fig);
end

function summary = local_plot_spectra(meta, file_name)
[KX, KY] = local_spectral_grid(meta.x_m, meta.y_m);
P_inc = abs(fftshift(fft2(meta.surface_incident_xy))).^2;
P_ref = abs(fftshift(fft2(meta.surface_reflected_xy))).^2;
stats_inc = local_spectrum_stats(P_inc, KX, KY);
stats_ref = local_spectrum_stats(P_ref, KX, KY);
common_ref = max([P_inc(:); P_ref(:); eps]);
P_inc_db = 10*log10(max(P_inc, eps) / common_ref);
P_ref_db = 10*log10(max(P_ref, eps) / common_ref);

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1200, 500]);
tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile
imagesc(KX(1, :), KY(:, 1), P_inc_db); axis image xy; clim([-60, 0]); colorbar
xlabel('K_x (rad/m)'); ylabel('K_y (rad/m)');
title(sprintf('Incident spectrum: RMS %.3g, r_{90} %.3g rad/m', ...
    stats_inc.rms_k_rad_per_m, stats_inc.energy_radius_90_rad_per_m));
nexttile
imagesc(KX(1, :), KY(:, 1), P_ref_db); axis image xy; clim([-60, 0]); colorbar
xlabel('K_x (rad/m)'); ylabel('K_y (rad/m)');
title(sprintf('Reflected spectrum: RMS %.3g, r_{90} %.3g rad/m', ...
    stats_ref.rms_k_rad_per_m, stats_ref.energy_radius_90_rad_per_m));
sgtitle(sprintf('Angular-spectrum comparison: %s', local_model_label(meta)), ...
    'Interpreter', 'none');
exportgraphics(fig, file_name, 'Resolution', 180);
close(fig);

summary = struct( ...
    'incident', stats_inc, ...
    'reflected', stats_ref, ...
    'rms_k_increase_rad_per_m', stats_ref.rms_k_rad_per_m - stats_inc.rms_k_rad_per_m, ...
    'energy_radius_90_increase_rad_per_m', ...
        stats_ref.energy_radius_90_rad_per_m - stats_inc.energy_radius_90_rad_per_m);
end

function written = local_write_animation(meta, paramsV, file_name)
written = false;
try
    writer = VideoWriter(file_name, 'MPEG-4');
    writer.FrameRate = 12;
    writer.Quality = 95;
    open(writer);

    inc_carrier = exp(1i * meta.k0_rad_per_m * ...
        (paramsV.z_tx - meta.incident_z_m(:).'));
    ref_carrier = exp(1i * meta.k0_rad_per_m * meta.reflected_z_m(:).');
    inc_phasor = meta.incident_field_slice .* inc_carrier;
    ref_phasor = meta.reflected_field_slice .* ref_carrier;
    pressure_ref = max([abs(inc_phasor(:)); abs(ref_phasor(:)); eps]);
    n_frame = 48;

    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80, 80, 1200, 500]);
    for frame_idx = 1:n_frame
        phase_t = 2*pi*(frame_idx - 1)/n_frame;
        p_inc = real(inc_phasor .* exp(-1i*phase_t)) / pressure_ref;
        p_ref = real(ref_phasor .* exp(-1i*phase_t)) / pressure_ref;
        clf(fig);
        layout = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
        nexttile(layout)
        imagesc(meta.transverse_coordinate_m, meta.incident_z_m, p_inc.');
        set(gca, 'YDir', 'reverse'); axis tight; clim([-1, 1]); colorbar
        xlabel(sprintf('%s (m)', meta.slice_axis)); ylabel('z (m)');
        title('Incident instantaneous pressure');
        nexttile(layout)
        imagesc(meta.transverse_coordinate_m, meta.reflected_z_m, p_ref.');
        set(gca, 'YDir', 'reverse'); axis tight; clim([-1, 1]); colorbar
        xlabel(sprintf('%s (m)', meta.slice_axis)); ylabel('z (m)');
        title('Reflected instantaneous pressure');
        colormap(fig, redblue_colormap(256));
        sgtitle(layout, sprintf('Nominal carrier reconstruction, t/T = %.3f', ...
            (frame_idx - 1)/n_frame));
        drawnow
        writeVideo(writer, getframe(fig));
    end
    close(writer);
    close(fig);
    written = true;
catch ME
    warning('Could not write %s: %s', file_name, ME.message);
    if exist('writer', 'var')
        try
            close(writer);
        catch
        end
    end
    if exist('fig', 'var') && isgraphics(fig)
        close(fig);
    end
end
end

function [KX, KY] = local_spectral_grid(x, y)
nx = numel(x);
ny = numel(y);
dx = mean(diff(x));
dy = mean(diff(y));
kx = (2*pi/(nx*dx)) * (-nx/2:(nx/2-1));
ky = (2*pi/(ny*dy)) * (-ny/2:(ny/2-1));
[KX, KY] = meshgrid(kx, ky);
end

function stats = local_spectrum_stats(power_k, KX, KY)
weight = max(real(power_k), 0);
total = sum(weight(:));
if total <= 0
    stats = struct('rms_k_rad_per_m', NaN, 'energy_radius_90_rad_per_m', NaN);
    return
end
k_radius = hypot(KX, KY);
stats.rms_k_rad_per_m = sqrt(sum(weight(:).*k_radius(:).^2) / total);
[radius_sorted, order] = sort(k_radius(:));
cumulative = cumsum(weight(order)) / total;
idx90 = find(cumulative >= 0.9, 1, 'first');
stats.energy_radius_90_rad_per_m = radius_sorted(idx90);
end

function label = local_model_label(meta)
if strcmp(meta.surface_boundary_model, 'ssa_stat_kernel')
    label = sprintf('%s / %s', meta.surface_boundary_model, meta.ssa_kernel_mode);
else
    label = meta.surface_boundary_model;
end
end

function cmap = redblue_colormap(n)
half_n = ceil(n/2);
blue_to_white = [linspace(0, 1, half_n).', linspace(0, 1, half_n).', ones(half_n, 1)];
red_n = n - half_n;
white_to_red = [ones(red_n, 1), linspace(1, 0, red_n).', linspace(1, 0, red_n).'];
cmap = [blue_to_white; white_to_red];
end

