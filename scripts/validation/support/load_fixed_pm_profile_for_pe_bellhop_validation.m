function profile = load_fixed_pm_profile_for_pe_bellhop_validation(coeff_file, options)
%LOAD_FIXED_PM_PROFILE_FOR_PE_BELLHOP_VALIDATION Load one canonical PM series.
%   PROFILE contains the immutable Fourier coefficients and provenance for
%   the fixed-seed Bellhop/PE comparison.  No random numbers, recentering,
%   variance scaling, smoothing or tapering are performed here.

if nargin < 1 || isempty(coeff_file)
    root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
    coeff_file = fullfile(root, 'results', 'validation', ...
        'bellhop_internal_pm_fixed_realization', ...
        'fixed_pm_fourier_coefficients.csv');
end
if nargin < 2 || isempty(options), options = struct(); end
if ~isstruct(options) || ~isscalar(options)
    error('options must be a scalar struct.');
end
if exist(coeff_file, 'file') ~= 2
    error('Canonical PM coefficient file does not exist: %s', coeff_file);
end

coeff = readtable(coeff_file);
required = {'k_rad_per_m','cos_coeff_m','sin_coeff_m'};
for ii = 1:numel(required)
    if ~ismember(required{ii}, coeff.Properties.VariableNames)
        error('Coefficient file is missing column %s.', required{ii});
    end
end
k = double(coeff.k_rad_per_m(:));
a = double(coeff.cos_coeff_m(:));
b = double(coeff.sin_coeff_m(:));
if isempty(k) || any(~isfinite([k; a; b])) || any(k <= 0) || any(diff(k) <= 0)
    error('Fourier coefficients must be finite with strictly increasing positive k.');
end

% The stored series is eta(s)=sum(a_n cos(k_n s)+b_n sin(k_n s)).
k0 = k(1);
fundamental_tol = max(1e-12, 1e-10 * k0);
harmonic_index = k / k0;
if any(abs(harmonic_index - round(harmonic_index)) > fundamental_tol / max(k0, eps))
    error('Coefficient wavenumbers are not harmonics of the first wavenumber.');
end
span_m = 2*pi/k0;

profile = struct();
profile.schema_version = '1.0.0';
profile.coeff_file = char(coeff_file);
profile.coeff_file_sha256 = local_sha256_file(coeff_file);
profile.k_rad_per_m = k;
profile.cos_coeff_m = a;
profile.sin_coeff_m = b;
profile.harmonic_index = round(harmonic_index);
profile.fundamental_k_rad_per_m = k0;
profile.span_m = span_m;
profile.coefficient_count = numel(k);
profile.seed = local_option(options, 'seed', 260001);
profile.wind_speed_mps = local_option(options, 'wind_speed_mps', 6);
profile.requested_kmax_rad_per_m = local_option(options, 'requested_kmax_rad_per_m', 0.5);
profile.realized_kmax_rad_per_m = max(k);
profile.datum_m = local_option(options, 'datum_m', 0);
profile.sign_convention = local_option(options, 'sign_convention', ...
    'project signed eta; use unchanged in PE and Gamma_BH=[R0-eta,s]');
profile.series_formula = 'eta(s)=sum_n(cos_coeff_m(n)*cos(k_n*s)+sin_coeff_m(n)*sin(k_n*s))';
profile.derivative_formula = 'eta_prime=-sum(k_n*cos_coeff_m(n)*sin(k_n*s))+sum(k_n*sin_coeff_m(n)*cos(k_n*s))';
profile.second_derivative_formula = 'eta_second=-sum(k_n^2*(cos_coeff_m(n)*cos(k_n*s)+sin_coeff_m(n)*sin(k_n*s)))';

if isfield(options, 's_eval_m') && ~isempty(options.s_eval_m)
    s_eval = double(options.s_eval_m(:));
    if any(~isfinite(s_eval))
        error('options.s_eval_m must contain finite values.');
    end
    profile.samples = evaluate_fixed_pm_fourier_profile(profile, s_eval);
else
    profile.samples = struct('s_m', zeros(0,1), 'eta_m', zeros(0,1), ...
        'eta_prime_m_per_m', zeros(0,1), 'eta_second_m_per_m2', zeros(0,1), ...
        'profile_hash_sha256', '');
end

function value = local_option(s, name, default_value)
if isfield(s, name) && ~isempty(s.(name))
    value = s.(name);
else
    value = default_value;
end
end

function digest = local_sha256_file(path)
md = java.security.MessageDigest.getInstance('SHA-256');
fid = fopen(path, 'rb');
if fid < 0, error('Cannot open coefficient file for hashing: %s', path); end
cleanup = onCleanup(@() fclose(fid));
bytes = fread(fid, Inf, '*uint8');
clear cleanup
md.update(typecast(bytes, 'int8'));
digest_bytes = typecast(md.digest(), 'uint8');
digest = lower(reshape(dec2hex(digest_bytes, 2).', 1, []));
end
end
