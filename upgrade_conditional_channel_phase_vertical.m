function [model,migration] = upgrade_conditional_channel_phase_vertical(model,geometry_context)
%UPGRADE_CONDITIONAL_CHANNEL_PHASE_VERTICAL Upgrade legacy reduced models.
%   Project schema 1.x stored reflected receiver statistics in the raw PE
%   reduced-envelope reference. Schema 2.x stores the direct-referenced DSP
%   convention. Unknown external channel arrays are intentionally outside
%   this project-model migration API.

arguments
    model (1,1) struct
    geometry_context (1,1) struct = struct()
end

if isfield(model,'conditions')
    [model,migration] = local_upgrade_library(model,geometry_context);
    return
end
version = local_field(model,'schema_version','');
if startsWith(version,'2.')
    migration = struct('applied',false,'source_schema',version, ...
        'target_schema',version,'reason','already_current');
    return
end
if ~startsWith(version,'1.')
    error('upgrade_conditional_channel_phase_vertical:UnknownSchema', ...
        'Only project conditional schemas 1.x and 2.x are supported.');
end
required = {'frequency_axis','H_direct_f','H_ref_coh_f','stats'};
for ii=1:numel(required)
    if ~isfield(model,required{ii})
        error('upgrade_conditional_channel_phase_vertical:InvalidModel', ...
            'Legacy model.%s is required.',required{ii});
    end
end

[geometry,geometry_source] = local_geometry(model,geometry_context);
f = model.frequency_axis(:);
dummy = struct('direct_f',model.H_direct_f(:), ...
    'reflect_coh_f',model.H_ref_coh_f(:));
[deterministic,phase_meta] = apply_pe_channel_phase_reference_vertical( ...
    f,dummy,geometry,'direct_dsp');
d = phase_meta.reflect_dsp_factor_f(:);

model.H_direct_reduced_f = model.H_direct_f(:);
model.H_ref_coh_reduced_f = model.H_ref_coh_f(:);
model.H_direct_f = deterministic.direct_f;
model.H_ref_coh_f = deterministic.reflect_coh_f;
s = model.stats;
s.mu_scatter_reduced_f = s.mu_scatter_f;
s.C_scatter_reduced_f = s.C_scatter_f;
s.P_scatter_reduced_f = s.P_scatter_f;
s.mu_scatter_f = d.*s.mu_scatter_f;
s.C_scatter_f = d.*s.C_scatter_f.*conj(d.');
s.P_scatter_f = d.*s.P_scatter_f.*d.';
if isfield(s,'C_scatter_raw_f')
    s.C_scatter_raw_reduced_f = s.C_scatter_raw_f;
    s.C_scatter_raw_f = d.*s.C_scatter_raw_f.*conj(d.');
end
if isfield(s,'eigenvectors')
    s.eigenvectors_reduced = s.eigenvectors;
    s.eigenvectors = d.*s.eigenvectors;
end
F = numel(f);
T = [diag(real(d)),-diag(imag(d));diag(imag(d)),diag(real(d))];
if isfield(s,'augmented_mean')
    s.augmented_mean_reduced = s.augmented_mean;
    s.augmented_mean = T*s.augmented_mean;
end
if isfield(s,'augmented_covariance')
    s.augmented_covariance_reduced = s.augmented_covariance;
    s.augmented_covariance = T*s.augmented_covariance*T.';
    s.augmented_covariance = 0.5*(s.augmented_covariance+s.augmented_covariance.');
end
if isfield(s,'augmented_eigenvectors')
    s.augmented_eigenvectors_reduced = s.augmented_eigenvectors;
    s.augmented_eigenvectors = T*s.augmented_eigenvectors;
end
model.stats = s;
model.phase_reference_meta = phase_meta;
model.phase_reference_meta.source = 'legacy_project_model_migration';
model.phase_reference_meta.geometry_source = geometry_source;
model.nominal_reflection_delay_s = phase_meta.relative_delay_s;
model.reference_delay_s = 0;
model.schema_version = '2.0.0';
migration = struct('applied',true,'source_schema',version, ...
    'target_schema','2.0.0','geometry_source',geometry_source, ...
    'relative_delay_s',phase_meta.relative_delay_s, ...
    'formula','mu=D*mu; C=D*C*D^H; P=D*P*D^T', ...
    'augmented_transform_size',[2*F,2*F]);
model.phase_migration = migration;
end

function [library,migration] = local_upgrade_library(library,context)
version=local_field(library,'schema_version','');
if startsWith(version,'2.')
    migration=struct('applied',false,'source_schema',version, ...
        'target_schema',version,'reason','already_current');
    return
end
if ~startsWith(version,'1.')
    error('upgrade_conditional_channel_phase_vertical:UnknownLibrarySchema', ...
        'Only project conditional-library schemas 1.x and 2.x are supported.');
end
if isempty(fieldnames(context)) && isfield(library,'fixed_geometry')
    context=library.fixed_geometry;
end
node_migrations=cell(numel(library.conditions),1);
for ii=1:numel(library.conditions)
    if ~isfield(library.conditions(ii),'model')
        error('upgrade_conditional_channel_phase_vertical:MissingNodeModel', ...
            'Legacy library condition %d has no nested model.',ii);
    end
    [m,node_migrations{ii}]=upgrade_conditional_channel_phase_vertical( ...
        library.conditions(ii).model,context);
    library.conditions(ii).model=m;
    library.conditions(ii).H_direct_f=m.H_direct_f;
    library.conditions(ii).H_ref_coh_f=m.H_ref_coh_f;
    library.conditions(ii).mu_scatter_f=m.stats.mu_scatter_f;
    library.conditions(ii).C_scatter_f=m.stats.C_scatter_f;
    library.conditions(ii).P_scatter_f=m.stats.P_scatter_f;
    library.conditions(ii).eigenvectors=m.stats.eigenvectors;
    library.conditions(ii).phase_reference_meta=m.phase_reference_meta;
end
library.schema_version='2.0.0-discrete';
library.phase_reference_meta=library.conditions(1).model.phase_reference_meta;
library.reference_delay=0;
migration=struct('applied',true,'source_schema',version, ...
    'target_schema',library.schema_version,'node_migrations',{node_migrations});
library.phase_migration=migration;
end

function [geometry,source]=local_geometry(model,context)
candidate=context;
source='geometry_context';
if isempty(fieldnames(candidate)) && isfield(model,'fixed_geometry')
    candidate=model.fixed_geometry;
    source='model.fixed_geometry';
end
if isfield(candidate,'cfg'), candidate=candidate.cfg; source=[source '.cfg']; end
if all(isfield(candidate,{'z_tx','z_rx','c0'}))
    geometry=struct('z_tx',candidate.z_tx,'z_rx',candidate.z_rx, ...
        'z_surface',local_field(candidate,'z_surface',0),'c0',candidate.c0);
    return
end
wind=NaN;
if isfield(model,'condition'), wind=local_field(model.condition,'wind_speed_mps',NaN); end
if any(abs(wind-[5,8])<1e-12) && numel(model.frequency_axis)==64 && ...
        abs(model.frequency_axis(1)-4000)<1e-9 && abs(model.frequency_axis(end)-8000)<1e-9
    geometry=struct('z_tx',100,'z_rx',3,'z_surface',0,'c0',1500);
    source='registered_u5_u8_f64_validation_geometry';
    return
end
error('upgrade_conditional_channel_phase_vertical:MissingGeometry', ...
    ['Legacy project model needs reliable z_tx, z_rx, z_surface, and c0. ', ...
    'The deprecated signed reference_delay_s is not sufficient.']);
end

function value=local_field(s,name,default)
if isfield(s,name), value=s.(name); else, value=default; end
end
