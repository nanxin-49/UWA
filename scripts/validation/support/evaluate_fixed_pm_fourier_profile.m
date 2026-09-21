function samples = evaluate_fixed_pm_fourier_profile(profile, s_m)
%EVALUATE_FIXED_PM_FOURIER_PROFILE Evaluate the canonical PM series.
%   The evaluator is deterministic and preserves the coefficient-defined
%   periodic function, including its signed datum and nonzero mean.

if ~isstruct(profile) || ~all(isfield(profile, {'k_rad_per_m','cos_coeff_m','sin_coeff_m'}))
    error('profile must be returned by load_fixed_pm_profile_for_pe_bellhop_validation.');
end
s = double(s_m(:));
if any(~isfinite(s))
    error('s_m must contain finite values.');
end
k = double(profile.k_rad_per_m(:));
a = double(profile.cos_coeff_m(:));
b = double(profile.sin_coeff_m(:));
phase = s * k.';
eta = cos(phase) * a + sin(phase) * b;
eta_prime = sin(phase) * (-k .* a) + cos(phase) * (k .* b);
eta_second = cos(phase) * (-(k.^2) .* a) + sin(phase) * (-(k.^2) .* b);

samples = struct();
samples.s_m = s;
samples.eta_m = eta;
samples.eta_prime_m_per_m = eta_prime;
samples.eta_second_m_per_m2 = eta_second;
samples.profile_hash_sha256 = local_hash_samples(s, eta, eta_prime, eta_second);
if isfield(profile, 'span_m'), samples.periodic_span_m = profile.span_m; end
end

function digest = local_hash_samples(varargin)
md = java.security.MessageDigest.getInstance('SHA-256');
for ii = 1:nargin
    values = double(varargin{ii}(:));
    bytes = typecast(values, 'uint8');
    md.update(typecast(bytes, 'int8'));
end
digest_bytes = typecast(md.digest(), 'uint8');
digest = lower(reshape(dec2hex(digest_bytes, 2).', 1, []));
end
