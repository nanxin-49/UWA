function [result,meta] = run_cached_joint_kstat_pe_executor_vertical(cache,deltaG_xy_fm)
%RUN_CACHED_JOINT_KSTAT_PE_EXECUTOR_VERTICAL Propagate scatter batches only.

arguments
    cache (1,1) struct
    deltaG_xy_fm {mustBeNumeric}
end
sz=size(deltaG_xy_fm);
if numel(sz)<4, sz(4)=1; end
F=numel(cache.f_axis_hz); M=sz(4);
if sz(3)~=F, error('Frequency dimension mismatch.'); end
if sz(1)==cache.cfg.ny && sz(2)==cache.cfg.nx
    delta_pe=deltaG_xy_fm;
    mapping_applied=false;
elseif sz(1)==cache.pm_mapping.pm_grid(1) && sz(2)==cache.pm_mapping.pm_grid(2)
    delta_pe=deltaG_xy_fm(cache.pm_mapping.iy,cache.pm_mapping.ix,:,:);
    mapping_applied=true;
else
    error('Input spatial grid matches neither PM nor PE grid.');
end

timer=tic;
H_ref_sca_reduced_fm=complex(zeros(F,M));
for ii=1:F
    incident=cache.incident_surface_xy_f(:,:,ii);
    pages=double(reshape(delta_pe(:,:,ii,:),cache.cfg.ny,cache.cfg.nx,M)).*incident;
    psi_end=apply_forward_surface_to_receiver_vertical(cache,pages,ii);
    H_ref_sca_reduced_fm(ii,:)=reshape(psi_end(cache.iy_rx,cache.ix_rx,:),1,M);
end
elapsed_s=toc(timer);
H_direct_reduced_f=local_cache_field(cache,'H_direct_reduced_f','H_direct_f');
H_ref_coh_reduced_f=local_cache_field(cache,'H_ref_coh_reduced_f','H_ref_coh_f');
geometry=local_geometry(cache);
reduced=struct('direct_f',H_direct_reduced_f, ...
    'reflect_coh_f',H_ref_coh_reduced_f, ...
    'reflect_sca_fm',H_ref_sca_reduced_fm);
[dsp,phase_meta]=apply_pe_channel_phase_reference_vertical( ...
    cache.f_axis_hz(:),reduced,geometry,'direct_dsp');
[physical,~]=apply_pe_channel_phase_reference_vertical( ...
    cache.f_axis_hz(:),reduced,geometry,'absolute_physical');
H_total_reduced_fm=H_direct_reduced_f+H_ref_coh_reduced_f+H_ref_sca_reduced_fm;
component_error=dsp.total_fm-dsp.direct_f-dsp.reflect_coh_f-dsp.reflect_sca_fm;
result=struct('H_direct_f',dsp.direct_f, ...
    'H_ref_coh_f',dsp.reflect_coh_f, ...
    'H_ref_sca_fm',dsp.reflect_sca_fm, ...
    'H_total_fm',dsp.total_fm, ...
    'H_direct_reduced_f',H_direct_reduced_f, ...
    'H_ref_coh_reduced_f',H_ref_coh_reduced_f, ...
    'H_ref_sca_reduced_fm',H_ref_sca_reduced_fm, ...
    'H_total_reduced_fm',H_total_reduced_fm, ...
    'H_direct_physical_f',physical.direct_f, ...
    'H_ref_coh_physical_f',physical.reflect_coh_f, ...
    'H_ref_sca_physical_fm',physical.reflect_sca_fm, ...
    'H_total_physical_fm',physical.total_fm, ...
    'phase_reference_meta',phase_meta, ...
    'component_sum_max_abs_error',max(abs(component_error(:))));
meta=struct('elapsed_s',elapsed_s,'per_realization_s',elapsed_s/M, ...
    'batch_size',M,'mapping_applied',mapping_applied, ...
    'working_array_bytes_estimate', ...
        cache.cfg.nx*cache.cfg.ny*M*16*4+numel(delta_pe)*8, ...
    'memory_snapshot_bytes',local_memory_used_bytes());
end

function value=local_cache_field(cache,preferred,legacy)
if isfield(cache,preferred), value=cache.(preferred)(:); else, value=cache.(legacy)(:); end
end

function geometry=local_geometry(cache)
if isfield(cache,'phase_geometry')
    geometry=cache.phase_geometry;
else
    geometry=struct('z_tx',cache.cfg.z_tx,'z_rx',cache.cfg.z_rx, ...
        'z_surface',0,'c0',cache.cfg.c0);
end
end

function value=local_memory_used_bytes()
value=NaN; try, info=memory; value=info.MemUsedMATLAB; catch, end
end
