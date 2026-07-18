function [sample,meta]=sample_conditional_channel_vertical(model,n_samples,seed,options)
%SAMPLE_CONDITIONAL_CHANNEL_VERTICAL Draw receiver channels without PE.

arguments
    model (1,1) struct
    n_samples (1,1) double {mustBeInteger,mustBePositive}
    seed (1,1) double {mustBeFinite}
    options (1,1) struct = struct()
end
if ~isfield(options,'path'), options.path='auto'; end
if ~isfield(options,'rank'), options.rank='selected'; end
if ~isfield(options,'return_total'), options.return_total=true; end
required={'frequency_axis','H_direct_f','H_ref_coh_f','stats'};
for ii=1:numel(required), if ~isfield(model,required{ii}), error('Missing model.%s.',required{ii}); end, end
rng(mod(round(seed),2^32),'twister'); timer=tic;
path=local_resolve_path(options.path,model.stats.properness_test);
F=numel(model.frequency_axis);
switch path
    case 'proper'
        r=local_rank(options.rank,model.stats.rank_selected,model.stats.rank_candidates,false);
        U=model.stats.eigenvectors(:,1:r); d=model.stats.eigenvalues(1:r);
        z=(randn(r,n_samples)+1i*randn(r,n_samples))/sqrt(2);
        Hsca=model.stats.mu_scatter_f+U*(sqrt(max(d,0)).*z);
        retained=sum(d)/max(sum(model.stats.eigenvalues),eps);
    case 'improper'
        r=local_rank(options.rank,model.stats.augmented_rank_selected,model.stats.rank_candidates,true);
        U=model.stats.augmented_eigenvectors(:,1:r); d=model.stats.augmented_eigenvalues(1:r);
        x=model.stats.augmented_mean+U*(sqrt(max(d,0)).*randn(r,n_samples));
        Hsca=x(1:F,:)+1i*x(F+1:end,:);
        retained=sum(d)/max(sum(model.stats.augmented_eigenvalues),eps);
    otherwise
        error('Internal sampling path error.');
end
if options.return_total, Htotal=model.H_direct_f+model.H_ref_coh_f+Hsca; else, Htotal=[]; end
sample=struct('sample_id',(1:n_samples).','seed',mod(round(seed),2^32), ...
    'wind_speed_mps',local_condition_value(model.condition,'wind_speed_mps',NaN), ...
    'Hs_implied_m',local_condition_value(model.condition,'Hs_implied_m',NaN), ...
    'H_ref_sca_f',Hsca,'H_total_f',Htotal,'model_version',model.schema_version);
elapsed=toc(timer);
meta=struct('path',path,'rank_used',r,'variance_retained',retained, ...
    'n_samples',n_samples,'elapsed_s',elapsed, ...
    'per_sample_s',elapsed/n_samples,'output_bytes',numel(Hsca)*16+numel(Htotal)*16);
end

function path=local_resolve_path(requested,test)
requested=lower(strtrim(requested));
if ~strcmp(requested,'auto'), path=requested; return; end
path='proper';
if isstruct(test) && isfield(test,'reject_proper_at_5pct') && test.reject_proper_at_5pct
    path='improper';
elseif isstruct(test) && isfield(test,'reject_proper_at_5pct_one_sided') && test.reject_proper_at_5pct_one_sided
    path='improper';
end
end

function rank=local_rank(choice,selected,c,augmented)
rank=selected;
choice=lower(strrep(strtrim(choice),'%',''));
if augmented, prefix='augmented_'; else, prefix=''; end
switch choice
    case 'selected', return
    case {'99','variance_99'}, rank=c.([prefix 'variance_99']);
    case {'99.9','999','variance_999'}, rank=c.([prefix 'variance_999']);
    case {'full','full_rank'}, rank=c.([prefix 'full']);
    otherwise, error('rank must be selected, full, 99.9, or 99.');
end
end

function value=local_condition_value(condition,field,default)
if isfield(condition,field), value=condition.(field); else, value=default; end
end
