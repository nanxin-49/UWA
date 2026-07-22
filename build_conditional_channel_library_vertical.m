function library=build_conditional_channel_library_vertical(models,options)
%BUILD_CONDITIONAL_CHANNEL_LIBRARY_VERTICAL Build a discrete wind-speed library.
%   MODELS is a struct array or cell array of validated one-condition models.
%   This first version intentionally rejects duplicate nodes and does not
%   construct any interpolation between wind speeds.

arguments
    models
    options (1,1) struct = struct()
end
if isstruct(models), models=num2cell(models); end
if ~iscell(models) || isempty(models), error('models must contain at least one condition model.'); end
geometry=local_option(options,'fixed_geometry',struct());
for jj=1:numel(models)
    models{jj}=upgrade_conditional_channel_phase_vertical(models{jj},geometry);
end
required={'frequency_axis','phase_reference_meta','H_direct_f','H_ref_coh_f','stats','condition'};
for jj=1:numel(models)
    for kk=1:numel(required)
        if ~isfield(models{jj},required{kk}), error('Model %d is missing %s.',jj,required{kk}); end
    end
end
f=models{1}.frequency_axis(:); phase_meta=models{1}.phase_reference_meta;
conditions=repmat(local_condition(models{1}),numel(models),1);
for jj=1:numel(models)
    m=models{jj};
    if numel(m.frequency_axis)~=numel(f) || any(abs(m.frequency_axis(:)-f)>1e-10)
        error('All condition nodes must use the same frequency axis.');
    end
    if ~strcmp(m.phase_reference_meta.target_reference,phase_meta.target_reference) || ...
            abs(m.phase_reference_meta.relative_delay_s-phase_meta.relative_delay_s)>1e-12
        error('All condition nodes must use the same phase reference and geometry delay.');
    end
    conditions(jj)=local_condition(m);
end
wind=[conditions.wind_speed_mps];
if any(~isfinite(wind)) || numel(unique(wind))~=numel(wind)
    error('Every node needs one unique finite condition.wind_speed_mps.');
end
[~,order]=sort(wind); conditions=conditions(order);
library=struct('schema_version','2.0.0-discrete', ...
    'code_revision',local_option(options,'code_revision','uncommitted'), ...
    'frequency_axis',f,'reference_delay',0, ...
    'nominal_reflection_delay_s',phase_meta.relative_delay_s, ...
    'phase_reference_meta',phase_meta, ...
    'fixed_geometry',local_option(options,'fixed_geometry',struct()), ...
    'fixed_pe_config',local_option(options,'fixed_pe_config',struct()), ...
    'conditions',conditions, ...
    'supported_wind_speeds_mps',[conditions.wind_speed_mps], ...
    'selection_policy','exact discrete node only; no silent nearest-node selection or interpolation', ...
    'limitations','No wind-speed interpolation, bubbles, Doppler, noise, or communication-chain integration.');
end

function c=local_condition(m)
s=m.stats; cond=m.condition;
c=struct('wind_speed_mps',local_field(cond,'wind_speed_mps',NaN), ...
    'wind_convention',local_field(cond,'wind_convention',''), ...
    'Hs_implied_m',local_field(cond,'Hs_implied_m',NaN), ...
    'raw_pm_meta',local_field(cond,'raw_pm_meta',cond), ...
    'H_direct_f',m.H_direct_f,'H_ref_coh_f',m.H_ref_coh_f, ...
    'mu_scatter_f',s.mu_scatter_f,'C_scatter_f',s.C_scatter_f, ...
    'P_scatter_f',s.P_scatter_f,'eigenvectors',s.eigenvectors, ...
    'eigenvalues',s.eigenvalues,'rank_selected',s.rank_selected, ...
    'properness_result',s.properness_test, ...
    'phase_reference_meta',m.phase_reference_meta, ...
    'validation',m.validation,'timing',m.timing,'model',m);
end

function value=local_option(s,name,default)
if isfield(s,name), value=s.(name); else, value=default; end
end

function value=local_field(s,name,default)
if isfield(s,name), value=s.(name); else, value=default; end
end
