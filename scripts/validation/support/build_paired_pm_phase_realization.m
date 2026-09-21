function profile = build_paired_pm_phase_realization(reference_profile, seed, output_file)
%BUILD_PAIRED_PM_PHASE_REALIZATION Make a deterministic paired PM realization.
%   The canonical per-mode spectral amplitudes are preserved exactly while
%   deterministic random phases are assigned from SEED.  This validation-only
%   ensemble keeps the same band, datum, and Kmax; it performs no Hs scaling,
%   recentering, smoothing, tapering, or bandwidth change.
if nargin<2||isempty(seed),error('seed is required.');end
if ~isstruct(reference_profile)||~all(isfield(reference_profile,{'k_rad_per_m','cos_coeff_m','sin_coeff_m'})),error('reference_profile is invalid.');end
if seed==reference_profile.seed
    profile=reference_profile; profile.realization_seed=seed; return
end
if nargin<3||isempty(output_file),error('output_file is required for non-reference seeds.');end
k=double(reference_profile.k_rad_per_m(:)); a0=double(reference_profile.cos_coeff_m(:)); b0=double(reference_profile.sin_coeff_m(:));
amp=hypot(a0,b0); stream=RandStream('mt19937ar','Seed',double(seed)); phi=2*pi*rand(stream,numel(k),1);
a=amp.*cos(phi); b=amp.*sin(phi);
tbl=table(k,a,b,'VariableNames',{'k_rad_per_m','cos_coeff_m','sin_coeff_m'}); writetable(tbl,output_file);
profile=load_fixed_pm_profile_for_pe_bellhop_validation(output_file,struct('seed',seed,'wind_speed_mps',reference_profile.wind_speed_mps,'requested_kmax_rad_per_m',reference_profile.requested_kmax_rad_per_m,'datum_m',reference_profile.datum_m));
profile.realization_seed=seed; profile.phase_ensemble='canonical per-mode amplitudes with deterministic random phases';
end
