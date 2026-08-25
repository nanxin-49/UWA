function model=estimate_conditional_channel_stats_vertical(H_scatter_fm,options)
%ESTIMATE_CONDITIONAL_CHANNEL_STATS_VERTICAL Estimate one-condition H statistics.
%   H_SCATTER_FM is F-by-L. OPTIONS carries deterministic components and
%   metadata; this independent interface does not alter the public channel API.

arguments
    H_scatter_fm {mustBeNumeric}
    options (1,1) struct = struct()
end

[F,L]=size(H_scatter_fm);
if F<1 || L<2, error('H_scatter_fm must be F-by-L with L >= 2.'); end
options=local_defaults(options,F);
if numel(options.frequency_axis_hz)~=F || numel(options.H_direct_f)~=F || numel(options.H_ref_coh_f)~=F
    error('Frequency axis and deterministic components must have F entries.');
end
timer=tic; H=double(H_scatter_fm);
mu=mean(H,2); X=H-mu;
C_raw=(X*X')/(L-1); C_raw=0.5*(C_raw+C_raw');
P=(X*X.')/(L-1); P=0.5*(P+P.');
gamma=options.shrinkage_parameter;
C=(1-gamma).*C_raw+gamma.*diag(real(diag(C_raw)));
C=0.5*(C+C');
[U,d]=eig(C,'vector'); d=real(d); [d,order]=sort(d,'descend'); U=U(:,order);
scale=max(max(abs(d)),eps); neg=d<0;
negative_energy=sum(abs(d(neg)))/max(sum(max(d,0)),eps);
d=max(d,0); positive=d>scale*1e-12;
full_rank=sum(positive); cumulative=cumsum(d)/max(sum(d),eps);
rank99=local_energy_rank(cumulative,0.99,full_rank);
rank999=local_energy_rank(cumulative,0.999,full_rank);

Xaug=[real(H);imag(H)]; mu_aug=mean(Xaug,2); Xaug=Xaug-mu_aug;
Sigma_aug=(Xaug*Xaug.')/(L-1); Sigma_aug=0.5*(Sigma_aug+Sigma_aug.');
[Uaug,daug]=eig(Sigma_aug,'vector'); daug=real(daug);
[daug,order]=sort(daug,'descend'); Uaug=Uaug(:,order);
aug_scale=max(max(abs(daug)),eps); aug_neg=daug<0;
aug_negative_energy=sum(abs(daug(aug_neg)))/max(sum(max(daug,0)),eps);
daug=max(daug,0); aug_full=sum(daug>aug_scale*1e-12);
aug_cumulative=cumsum(daug)/max(sum(daug),eps);

rank_candidates=struct('variance_99',rank99,'variance_999',rank999,'full',full_rank, ...
    'augmented_variance_99',local_energy_rank(aug_cumulative,0.99,aug_full), ...
    'augmented_variance_999',local_energy_rank(aug_cumulative,0.999,aug_full), ...
    'augmented_full',aug_full);
rank_selected=local_select_rank(options.rank_selection,rank_candidates,false);
augmented_rank_selected=local_select_rank(options.rank_selection,rank_candidates,true);

stats=struct('mu_scatter_f',mu,'C_scatter_f',C,'C_scatter_raw_f',C_raw, ...
    'P_scatter_f',P,'eigenvectors',U,'eigenvalues',d, ...
    'cumulative_variance',cumulative,'rank_candidates',rank_candidates, ...
    'rank_selected',rank_selected,'augmented_rank_selected',augmented_rank_selected, ...
    'variance_retained',sum(d(1:rank_selected))/max(sum(d),eps), ...
    'properness_ratio',norm(P,'fro')/max(norm(C_raw,'fro'),eps), ...
    'properness_test',options.properness_result, ...
    'covariance_estimator','unbiased sample covariance; optional diagonal-target linear shrinkage', ...
    'shrinkage_parameter',gamma,'negative_eigenvalue_energy',negative_energy, ...
    'augmented_mean',mu_aug,'augmented_covariance',Sigma_aug, ...
    'augmented_eigenvectors',Uaug,'augmented_eigenvalues',daug, ...
    'augmented_cumulative_variance',aug_cumulative, ...
    'augmented_negative_eigenvalue_energy',aug_negative_energy);

phase_meta=local_phase_meta(options.phase_reference_meta,options.frequency_axis_hz);
model=struct('kind','conditional_receiver_channel_vertical', ...
    'schema_version','2.0.0','code_revision',options.code_revision, ...
    'condition',options.condition,'frequency_axis',options.frequency_axis_hz(:), ...
    'reference_delay_s',0, ...
    'nominal_reflection_delay_s',local_field(phase_meta,'relative_delay_s',NaN), ...
    'phase_reference_meta',phase_meta, ...
    'H_direct_f',options.H_direct_f(:),'H_ref_coh_f',options.H_ref_coh_f(:), ...
    'stats',stats,'train_seed_list',options.train_seed_list(:).', ...
    'training_sample_count',L,'validation',struct(),'timing',struct('estimate_s',toc(timer)), ...
    'limitations',['One fixed environmental condition; noise excluded; ', ...
        'no interpolation, bubbles, Doppler, or communication-chain integration.']);
end

function o=local_defaults(o,F)
d=struct('frequency_axis_hz',(1:F).','H_direct_f',complex(zeros(F,1)), ...
    'H_ref_coh_f',complex(zeros(F,1)),'shrinkage_parameter',0, ...
    'properness_result',struct(),'train_seed_list',[],'condition',struct(), ...
    'reference_delay_s',0,'phase_reference_meta',struct(), ...
    'rank_selection','full','code_revision','uncommitted');
names=fieldnames(d);
for ii=1:numel(names), if ~isfield(o,names{ii}), o.(names{ii})=d.(names{ii}); end, end
validateattributes(o.shrinkage_parameter,{'numeric'},{'scalar','>=',0,'<=',1});
validateattributes(o.reference_delay_s,{'numeric'},{'scalar','finite'});
end

function meta=local_phase_meta(meta,f_axis)
if isempty(fieldnames(meta))
    meta=struct('schema_version','1.0.0','target_reference','direct_dsp', ...
        'source','assumed_dsp_ready_input_without_project_phase_metadata', ...
        'relative_delay_s',NaN,'frequency_axis_hz',f_axis(:));
elseif ~isfield(meta,'target_reference') || ~strcmp(meta.target_reference,'direct_dsp')
    error('estimate_conditional_channel_stats_vertical:PhaseReference', ...
        'Training samples must use target_reference=direct_dsp.');
end
end

function value=local_field(s,name,default)
if isfield(s,name), value=s.(name); else, value=default; end
end

function rank=local_energy_rank(cumulative,target,full_rank)
if full_rank==0, rank=0; return; end
rank=find(cumulative>=target,1,'first'); rank=min(max(rank,1),full_rank);
end

function rank=local_select_rank(choice,candidates,augmented)
choice=lower(strrep(strtrim(choice),'%',''));
if augmented, prefix='augmented_'; else, prefix=''; end
switch choice
    case {'99','variance_99'}, field=[prefix 'variance_99'];
    case {'99.9','999','variance_999'}, field=[prefix 'variance_999'];
    case {'full','full_rank'}, field=[prefix 'full'];
    otherwise, error('rank_selection must be full, 99.9, or 99.');
end
rank=candidates.(field);
end
