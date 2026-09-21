function profile = build_paired_pm_gaussian_realization(reference_profile, seed, output_file)
%BUILD_PAIRED_PM_GAUSSIAN_REALIZATION Draw one paired PM Fourier realization.
%   The one-sided spectral density and discrete wavenumber band are copied
%   from the canonical coefficient file.  Cosine and sine coefficients are
%   independent zero-mean Gaussian draws with variance S(k)*Delta-k.  No
%   height renormalization, recentering, smoothing, tapering, or bandwidth
%   change is applied.  This helper is validation-only.

if nargin < 2 || isempty(seed), error('seed is required.'); end
if nargin < 3 || isempty(output_file), error('output_file is required.'); end
if ~isstruct(reference_profile) || ~all(isfield(reference_profile, ...
        {'coeff_file','k_rad_per_m','wind_speed_mps','requested_kmax_rad_per_m','datum_m'}))
    error('reference_profile is invalid.');
end
if isfield(reference_profile, 'seed') && seed == reference_profile.seed
    % Preserve the canonical reference byte-for-byte for the paired anchor.
    copyfile(reference_profile.coeff_file, output_file, 'f');
    profile = reference_profile;
    profile.realization_seed = seed;
    profile.ensemble_mode = 'canonical reference realization';
    return
end

tbl = readtable(reference_profile.coeff_file);
required = {'k_rad_per_m','Sk_m3'};
for ii = 1:numel(required)
    if ~ismember(required{ii}, tbl.Properties.VariableNames)
        error('Canonical coefficient file is missing %s.', required{ii});
    end
end
k = double(tbl.k_rad_per_m(:));
sk = double(tbl.Sk_m3(:));
if numel(k) < 2 || any(~isfinite([k; sk])) || any(k <= 0) || any(sk < 0) || any(diff(k) <= 0)
    error('Canonical wavenumbers/spectrum are invalid.');
end
dk = median(diff(k));
if max(abs(diff(k)-dk)) > 1e-10 * max(dk, eps)
    error('Canonical wavenumber grid is not uniform.');
end

stream = RandStream('mt19937ar', 'Seed', double(seed));
sigma_coeff = sqrt(max(sk, 0) * dk);
a = sigma_coeff .* randn(stream, numel(k), 1);
b = sigma_coeff .* randn(stream, numel(k), 1);
out_tbl = table(k, a, b, sk, 'VariableNames', ...
    {'k_rad_per_m','cos_coeff_m','sin_coeff_m','Sk_m3'});
writetable(out_tbl, output_file);

profile = load_fixed_pm_profile_for_pe_bellhop_validation(output_file, ...
    struct('seed', seed, 'wind_speed_mps', reference_profile.wind_speed_mps, ...
    'requested_kmax_rad_per_m', reference_profile.requested_kmax_rad_per_m, ...
    'datum_m', reference_profile.datum_m));
profile.realization_seed = seed;
profile.ensemble_mode = 'independent Gaussian cosine/sine coefficients';
profile.spectral_density_column = 'Sk_m3';
profile.delta_k_rad_per_m = dk;
profile.coefficient_std_m = sigma_coeff;
end
