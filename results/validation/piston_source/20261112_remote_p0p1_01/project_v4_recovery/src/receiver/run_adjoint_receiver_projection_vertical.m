function [result, meta] = run_adjoint_receiver_projection_vertical( ...
    projection, deltaG_xy_fm)
%RUN_ADJOINT_RECEIVER_PROJECTION_VERTICAL Project K-stat residuals to one Rx.
%   DELTAG may be on the PE grid or the mapped PM grid. The input residual
%   already includes R0, so the stored a_pe weights contain no extra R0.

arguments
    projection (1,1) struct
    deltaG_xy_fm {mustBeNumeric}
end
required = {'f_axis_hz','a_pe_xy_f','H_direct_f','H_ref_coh_f','pm_mapping'};
for ii = 1:numel(required)
    if ~isfield(projection,required{ii})
        error('run_adjoint_receiver_projection_vertical:InvalidProjection', ...
            'projection.%s is required.',required{ii});
    end
end
if isa(deltaG_xy_fm,'gpuArray')
    error('run_adjoint_receiver_projection_vertical:CPURequired', ...
        'v1 does not support GPU input.');
end
sz = size(deltaG_xy_fm);
if numel(sz) < 4
    sz(4) = 1;
end
F = numel(projection.f_axis_hz);
M = sz(4);
pe_grid = projection.pm_mapping.pe_grid;
pm_grid = projection.pm_mapping.pm_grid;
if sz(3) ~= F
    error('run_adjoint_receiver_projection_vertical:FrequencyMismatch', ...
        'Input frequency dimension does not match projection.f_axis_hz.');
end
if sz(1) == pe_grid(1) && sz(2) == pe_grid(2)
    delta_pe = deltaG_xy_fm;
    mapping_applied = false;
elseif sz(1) == pm_grid(1) && sz(2) == pm_grid(2)
    delta_pe = deltaG_xy_fm(projection.pm_mapping.iy,projection.pm_mapping.ix,:,:);
    mapping_applied = true;
else
    error('run_adjoint_receiver_projection_vertical:GridMismatch', ...
        'Input grid matches neither the cached PE grid nor the mapped PM grid.');
end

timer = tic;
H_ref_sca_reduced_fm = complex(zeros(F,M));
for ii = 1:F
    weights = projection.a_pe_xy_f(:,:,ii);
    pages = double(reshape(delta_pe(:,:,ii,:),pe_grid(1),pe_grid(2),M));
    weighted = conj(weights).*pages;
    H_ref_sca_reduced_fm(ii,:) = reshape(sum(sum(weighted,1),2),1,M);
end
elapsed_s = toc(timer);
H_direct_reduced_f = local_projection_field(projection, ...
    'H_direct_reduced_f','H_direct_f');
H_ref_coh_reduced_f = local_projection_field(projection, ...
    'H_ref_coh_reduced_f','H_ref_coh_f');
geometry = local_geometry(projection);
reduced = struct('direct_f',H_direct_reduced_f, ...
    'reflect_coh_f',H_ref_coh_reduced_f, ...
    'reflect_sca_fm',H_ref_sca_reduced_fm);
[dsp,phase_meta] = apply_pe_channel_phase_reference_vertical( ...
    projection.f_axis_hz(:),reduced,geometry,'direct_dsp');
[physical,~] = apply_pe_channel_phase_reference_vertical( ...
    projection.f_axis_hz(:),reduced,geometry,'absolute_physical');
H_total_reduced_fm = H_direct_reduced_f+H_ref_coh_reduced_f+H_ref_sca_reduced_fm;
component_error = dsp.total_fm-dsp.direct_f-dsp.reflect_coh_f-dsp.reflect_sca_fm;
result = struct('H_direct_f',dsp.direct_f, ...
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
meta = struct('elapsed_s',elapsed_s, ...
    'per_realization_s',elapsed_s/M, ...
    'batch_size',M, ...
    'mapping_applied',mapping_applied, ...
    'input_precision',class(deltaG_xy_fm), ...
    'working_array_bytes_estimate',numel(delta_pe)*8+numel(H_ref_sca_reduced_fm)*16, ...
    'memory_snapshot_bytes',local_memory_used_bytes());
end

function value = local_projection_field(projection,preferred,legacy)
if isfield(projection,preferred), value=projection.(preferred)(:); else, value=projection.(legacy)(:); end
end

function geometry = local_geometry(projection)
if isfield(projection,'phase_geometry')
    geometry=projection.phase_geometry;
elseif isfield(projection,'phase_reference_meta')
    m=projection.phase_reference_meta;
    geometry=struct('z_tx',m.z_tx_m,'z_rx',m.z_rx_m, ...
        'z_surface',m.z_surface_m,'c0',m.c0_mps);
else
    error('run_adjoint_receiver_projection_vertical:MissingPhaseGeometry', ...
        'Projection needs phase_geometry for direct-DSP receiver assembly.');
end
end

function value = local_memory_used_bytes()
value = NaN;
try
    info = memory;
    value = info.MemUsedMATLAB;
catch
end
end
