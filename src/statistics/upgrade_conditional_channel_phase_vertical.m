function [model,migration] = upgrade_conditional_channel_phase_vertical(model,geometry_context)
%UPGRADE_CONDITIONAL_CHANNEL_PHASE_VERTICAL Upgrade legacy reduced models.
%   Migration requires explicit source_representation and geometry. A schema
%   number alone never identifies the Fourier convention. Native reduced uses
%   D*conj(z); historical carrier-only DSP uses D^2*conj(z).

arguments
    model (1,1) struct
    geometry_context (1,1) struct = struct()
end

if isfield(model,'conditions')
    [model,migration] = local_upgrade_library(model,geometry_context);
    return
end
version = local_field(model,'schema_version','');
if startsWith(version,'3.')
    assert_channel_dsp_metadata_vertical(model.phase_reference_meta);
    migration = struct('applied',false,'source_schema',version, ...
        'target_schema',version,'reason','already_current');
    return
end
if ~startsWith(version,'1.') && ~startsWith(version,'2.')
    error('upgrade_conditional_channel_phase_vertical:UnknownSchema', ...
        'Supported schemas are explicit-source 1.x/2.x or current 3.x.');
end
required = {'frequency_axis','H_direct_f','H_ref_coh_f','stats'};
for ii=1:numel(required)
    if ~isfield(model,required{ii})
        error('upgrade_conditional_channel_phase_vertical:InvalidModel', ...
            'Legacy model.%s is required.',required{ii});
    end
end

[geometry,geometry_source] = local_geometry(model,geometry_context);
if ~isfield(geometry_context,'source_representation')
    error('upgrade_conditional_channel_phase_vertical:MissingRepresentation', ...
        'Declare source_representation=native_pe_reduced or legacy_carrier_only_direct_dsp; schema alone is insufficient.');
end
source_representation = char(geometry_context.source_representation);
if ~ismember(source_representation,{'native_pe_reduced','legacy_carrier_only_direct_dsp'})
    error('upgrade_conditional_channel_phase_vertical:SourceRepresentation','Unsupported migration source.');
end
f = model.frequency_axis(:);
phase_d = exp(-1i*2*pi*f*((geometry.z_tx+geometry.z_rx-2*geometry.z_surface- ...
    (geometry.z_tx-geometry.z_rx))/geometry.c0));
old_direct = model.H_direct_f(:); old_coherent = model.H_ref_coh_f(:);
native_coherent = old_coherent;
native_inverse = ones(size(f));
if strcmp(source_representation,'legacy_carrier_only_direct_dsp')
    native_coherent = old_coherent./phase_d;
    native_inverse = 1./phase_d;
end
dummy = struct('direct_f',model.H_direct_f(:), ...
    'reflect_coh_f',native_coherent);
[deterministic,phase_meta] = apply_pe_channel_phase_reference_vertical( ...
    f,dummy,geometry,'direct_dsp');
d = phase_meta.reflect_dsp_factor_f(:);

model.H_direct_reduced_f = old_direct;
model.H_ref_coh_reduced_f = native_coherent;
model.H_direct_f = deterministic.direct_f;
model.H_ref_coh_f = deterministic.reflect_coh_f;
s = model.stats;
s.mu_scatter_reduced_f = native_inverse.*s.mu_scatter_f;
s.C_scatter_reduced_f = native_inverse.*s.C_scatter_f.*conj(native_inverse.');
s.P_scatter_reduced_f = native_inverse.*s.P_scatter_f.*native_inverse.';
a = d;
if strcmp(source_representation,'legacy_carrier_only_direct_dsp'), a=d.^2; end
s.mu_scatter_f = a.*conj(s.mu_scatter_f);
s.C_scatter_f = a.*conj(s.C_scatter_f).*conj(a.');
s.P_scatter_f = a.*conj(s.P_scatter_f).*a.';
if isfield(s,'C_scatter_raw_f')
    s.C_scatter_raw_reduced_f = native_inverse.*s.C_scatter_raw_f.*conj(native_inverse.');
    s.C_scatter_raw_f = a.*conj(s.C_scatter_raw_f).*conj(a.');
end
if isfield(s,'eigenvectors')
    s.eigenvectors_reduced = native_inverse.*s.eigenvectors;
    s.eigenvectors = a.*conj(s.eigenvectors);
end
F = numel(f);
% Real representation of a*conj(z): rotation(a)*diag(I,-I).
T = [diag(real(a)),diag(imag(a));diag(imag(a)),-diag(real(a))];
Rinv = [diag(real(native_inverse)),-diag(imag(native_inverse)); ...
    diag(imag(native_inverse)),diag(real(native_inverse))];
if isfield(s,'augmented_mean')
    s.augmented_mean_reduced = Rinv*s.augmented_mean;
    s.augmented_mean = T*s.augmented_mean;
end
if isfield(s,'augmented_covariance')
    s.augmented_covariance_reduced = Rinv*s.augmented_covariance*Rinv.';
    s.augmented_covariance = T*s.augmented_covariance*T.';
    s.augmented_covariance = 0.5*(s.augmented_covariance+s.augmented_covariance.');
end
if isfield(s,'augmented_eigenvectors')
    s.augmented_eigenvectors_reduced = Rinv*s.augmented_eigenvectors;
    s.augmented_eigenvectors = T*s.augmented_eigenvectors;
end
model.stats = s;
model.phase_reference_meta = phase_meta;
model.phase_reference_meta.source = 'legacy_project_model_migration';
model.phase_reference_meta.geometry_source = geometry_source;
model.nominal_reflection_delay_s = phase_meta.relative_delay_s;
model.reference_delay_s = 0;
model.schema_version = '3.0.0';
migration = struct('applied',true,'source_schema',version, ...
    'target_schema','3.0.0','geometry_source',geometry_source, ...
    'source_representation',source_representation, ...
    'relative_delay_s',phase_meta.relative_delay_s, ...
    'formula','z_new=A*conj(z_old); A=D (native), D^2 (carrier-only DSP)', ...
    'augmented_transform_size',[2*F,2*F]);
model.phase_migration = migration;
end

function [library,migration] = local_upgrade_library(library,context)
version=local_field(library,'schema_version','');
if startsWith(version,'3.')
    for ii=1:numel(library.conditions)
        assert_channel_dsp_metadata_vertical(library.conditions(ii).model.phase_reference_meta);
    end
    migration=struct('applied',false,'source_schema',version, ...
        'target_schema',version,'reason','already_current');
    return
end
if ~startsWith(version,'1.') && ~startsWith(version,'2.')
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
library.schema_version='3.0.0-discrete';
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
error('upgrade_conditional_channel_phase_vertical:MissingGeometry', ...
    ['Legacy project model needs reliable z_tx, z_rx, z_surface, and c0. ', ...
    'The deprecated signed reference_delay_s is not sufficient.']);
end

function value=local_field(s,name,default)
if isfield(s,name), value=s.(name); else, value=default; end
end
