function [eta_xy, meta] = sample_raw_pm_surface_vertical(spec, seed)
%SAMPLE_RAW_PM_SURFACE_VERTICAL Reproduce the explicit raw-PM surface convention.

arguments
    spec (1,1) struct
    seed (1,1) double {mustBeFinite}
end
required = {'Phi2D', 'dkx_rad_per_m', 'dky_rad_per_m', 'nx', 'ny'};
for ii = 1:numel(required)
    if ~isfield(spec, required{ii})
        error('sample_raw_pm_surface_vertical:MissingField', ...
            'spec.%s is required.', required{ii});
    end
end

seed_used = mod(round(seed), 2^32);
rng(seed_used, 'twister');
amplitude = sqrt(spec.Phi2D .* spec.dkx_rad_per_m .* spec.dky_rad_per_m);
z_k = amplitude .* ((randn(size(amplitude)) + 1i*randn(size(amplitude))) / sqrt(2));
eta_xy = sqrt(2) * numel(amplitude) * real(ifft2(z_k));

sample_variance_about_zero = mean(eta_xy(:).^2);
sample_variance_about_mean = var(eta_xy(:), 1);
target_variance = spec.sigma_eta_discrete2_m2;
meta = struct( ...
    'seed', seed_used, ...
    'generation_rule', ['sqrt(2)*real(ifft2(sqrt(Phi2D*dkx*dky).*CN(0,1))) ', ...
        'with MATLAB ifft2 normalization compensated'], ...
    'mean_eta_m', mean(eta_xy(:)), ...
    'variance_about_zero_m2', sample_variance_about_zero, ...
    'variance_about_mean_m2', sample_variance_about_mean, ...
    'target_discrete_variance_m2', target_variance, ...
    'energy_rel_error_about_zero', ...
        abs(sample_variance_about_zero - target_variance) / max(target_variance, eps));
end
