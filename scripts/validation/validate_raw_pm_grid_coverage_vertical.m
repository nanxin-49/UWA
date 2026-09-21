run(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'bootstrap_project.m'));
%VALIDATE_RAW_PM_GRID_COVERAGE_VERTICAL Reproducible PM aperture/grid audit.

clear
format compact

wind_list = [3, 5, 8, 10, 12, 15];
aperture_list = [50, 100, 200, 300, 400];
dx_target = 50 / 128;
grid_for_aperture = round(aperture_list / dx_target);
grid_for_aperture = 2 * round(grid_for_aperture / 2);
resolution_grid_list = [64, 128, 256, 512];
mapping_seed_list = 12345:12376;

rows = struct([]);
idx = 0;
for iu = 1:numel(wind_list)
    for ia = 1:numel(aperture_list)
        idx = idx + 1;
        row = local_row(raw_pm_spectrum_grid_vertical( ...
            wind_list(iu), grid_for_aperture(ia), grid_for_aperture(ia), ...
            aperture_list(ia), aperture_list(ia)), 'fixed_dx_aperture_sweep');
        if isempty(rows), rows = row; else, rows(end+1) = row; end %#ok<SAGROW>
    end
    for in = 1:numel(resolution_grid_list)
        idx = idx + 1;
        row = local_row(raw_pm_spectrum_grid_vertical( ...
            wind_list(iu), resolution_grid_list(in), resolution_grid_list(in), ...
            50, 50), 'fixed_aperture_resolution_sweep');
        rows(end+1) = row; %#ok<SAGROW>
    end
end
coverage_table = struct2table(rows);

% Decoupled mapping audit: synthesize on 100 m / 256^2, crop the central
% 50 m / 128^2 PE window without interpolation, and retain the absolute
% low-frequency offset rather than de-meaning the crop.
spec_stat = raw_pm_spectrum_grid_vertical(5, 256, 256, 100, 100);
n_pe = 128;
i0 = floor((spec_stat.nx - n_pe) / 2) + 1;
crop_idx = i0:(i0 + n_pe - 1);
mapped_energy = zeros(numel(mapping_seed_list), 1);
mapped_variance_demeaned = zeros(numel(mapping_seed_list), 1);
full_energy = zeros(numel(mapping_seed_list), 1);
for iseed = 1:numel(mapping_seed_list)
    eta_big = sample_raw_pm_surface_vertical(spec_stat, mapping_seed_list(iseed));
    eta_crop = eta_big(crop_idx, crop_idx);
    full_energy(iseed) = mean(eta_big(:).^2);
    mapped_energy(iseed) = mean(eta_crop(:).^2);
    mapped_variance_demeaned(iseed) = var(eta_crop(:), 1);
end
mapping_audit = struct();
mapping_audit.U_mps = 5;
mapping_audit.pm_grid = [256, 256];
mapping_audit.pm_aperture_m = [100, 100];
mapping_audit.pe_grid = [128, 128];
mapping_audit.pe_aperture_m = [50, 50];
mapping_audit.seed_count = numel(mapping_seed_list);
mapping_audit.target_energy_m2 = spec_stat.sigma_eta_discrete2_m2;
mapping_audit.full_energy_mean_m2 = mean(full_energy);
mapping_audit.mapped_energy_mean_m2 = mean(mapped_energy);
mapping_audit.mapped_energy_std_m2 = std(mapped_energy);
mapping_audit.mapped_demeaned_variance_mean_m2 = mean(mapped_variance_demeaned);
mapping_audit.mapping_energy_rel_error = abs(mean(mapped_energy) - ...
    spec_stat.sigma_eta_discrete2_m2) / spec_stat.sigma_eta_discrete2_m2;
mapping_audit.demeaning_energy_loss_fraction = 1 - ...
    mean(mapped_variance_demeaned) / max(mean(mapped_energy), eps);
mapping_audit.mapping_rule = ['central crop at identical dx; no interpolation; ', ...
    'absolute elevation offset retained'];

result_file = project_result_file('validation', ...
    'validate_raw_pm_grid_coverage_vertical_result.mat');
csv_file = project_result_file('validation', ...
    'validate_raw_pm_grid_coverage_vertical_table.csv');
save(result_file, 'wind_list', 'aperture_list', 'dx_target', ...
    'grid_for_aperture', 'resolution_grid_list', 'coverage_table', ...
    'mapping_seed_list', 'mapping_audit');
writetable(coverage_table, csv_file);

disp(coverage_table(strcmp(coverage_table.sweep_mode, 'fixed_dx_aperture_sweep'), ...
    {'U_mps','xw_m','nx','K_min_rad_per_m','K_peak_rad_per_m', ...
    'K_nyquist_axis_rad_per_m','capture_ratio_discrete_to_infinite', ...
    'capture_ratio_radial_support_idealized','Hs_implied_discrete_m'}));
disp(mapping_audit)
fprintf('Saved %s\n', result_file);
fprintf('Saved %s\n', csv_file);

function row = local_row(spec, sweep_mode)
row = struct( ...
    'sweep_mode', string(sweep_mode), ...
    'U_mps', spec.U_mps, ...
    'xw_m', spec.xw_m, ...
    'yw_m', spec.yw_m, ...
    'nx', spec.nx, ...
    'ny', spec.ny, ...
    'dx_m', spec.dx_m, ...
    'dkx_rad_per_m', spec.dkx_rad_per_m, ...
    'K_min_rad_per_m', spec.K_min_rad_per_m, ...
    'K_peak_rad_per_m', spec.K_peak_rad_per_m, ...
    'K_nyquist_axis_rad_per_m', spec.K_nyquist_axis_rad_per_m, ...
    'sigma_eta_discrete2_m2', spec.sigma_eta_discrete2_m2, ...
    'sigma_eta_infinite2_m2', spec.sigma_eta_infinite2_m2, ...
    'Hs_implied_discrete_m', spec.Hs_implied_discrete_m, ...
    'Hs_implied_infinite_m', spec.Hs_implied_infinite_m, ...
    'capture_ratio_discrete_to_infinite', spec.capture_ratio_discrete_to_infinite, ...
    'capture_ratio_radial_support_idealized', spec.capture_ratio_radial_support_idealized, ...
    'peak_bins_from_origin', spec.peak_bins_from_origin);
end
