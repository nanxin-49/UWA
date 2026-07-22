function pair=sample_conditional_channel_rank_pair_vertical(model,n_samples,seed,low_rank)
%SAMPLE_CONDITIONAL_CHANNEL_RANK_PAIR_VERTICAL Paired full/low proper draws.
%   Both outputs use the same latent Gaussian coefficients. The difference
%   therefore isolates discarded covariance modes rather than Monte Carlo
%   seed variation. Intended for compression validation only.
arguments
    model (1,1) struct
    n_samples (1,1) double {mustBeInteger,mustBePositive}
    seed (1,1) double {mustBeFinite}
    low_rank (1,:) char {mustBeMember(low_rank,{'99.9','99'})} = '99.9'
end
model=upgrade_conditional_channel_phase_vertical(model,struct());
s=model.stats; rfull=s.rank_candidates.full;
if strcmp(low_rank,'99.9'), rlow=s.rank_candidates.variance_999;
else, rlow=s.rank_candidates.variance_99; end
rng(mod(round(seed),2^32),'twister');
z=(randn(rfull,n_samples)+1i*randn(rfull,n_samples))/sqrt(2);
A=s.eigenvectors(:,1:rfull).*sqrt(max(s.eigenvalues(1:rfull).',0));
Hfull=s.mu_scatter_f+A*z;
Hlow=s.mu_scatter_f+A(:,1:rlow)*z(1:rlow,:);
pair=struct('H_full_f',model.H_direct_f+model.H_ref_coh_f+Hfull, ...
    'H_low_f',model.H_direct_f+model.H_ref_coh_f+Hlow, ...
    'rank_full',rfull,'rank_low',rlow,'low_rank_name',low_rank, ...
    'seed',seed,'n_samples',n_samples,'shared_latent_coefficients',true);
end
