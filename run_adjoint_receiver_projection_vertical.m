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
H_ref_sca_fm = complex(zeros(F,M));
for ii = 1:F
    weights = projection.a_pe_xy_f(:,:,ii);
    pages = double(reshape(delta_pe(:,:,ii,:),pe_grid(1),pe_grid(2),M));
    weighted = conj(weights).*pages;
    H_ref_sca_fm(ii,:) = reshape(sum(sum(weighted,1),2),1,M);
end
elapsed_s = toc(timer);
H_total_fm = projection.H_direct_f(:)+projection.H_ref_coh_f(:)+H_ref_sca_fm;
component_error = H_total_fm-projection.H_direct_f(:)- ...
    projection.H_ref_coh_f(:)-H_ref_sca_fm;
result = struct('H_direct_f',projection.H_direct_f(:), ...
    'H_ref_coh_f',projection.H_ref_coh_f(:), ...
    'H_ref_sca_fm',H_ref_sca_fm, ...
    'H_total_fm',H_total_fm, ...
    'component_sum_max_abs_error',max(abs(component_error(:))));
meta = struct('elapsed_s',elapsed_s, ...
    'per_realization_s',elapsed_s/M, ...
    'batch_size',M, ...
    'mapping_applied',mapping_applied, ...
    'input_precision',class(deltaG_xy_fm), ...
    'working_array_bytes_estimate',numel(delta_pe)*8+numel(H_ref_sca_fm)*16, ...
    'memory_snapshot_bytes',local_memory_used_bytes());
end

function value = local_memory_used_bytes()
value = NaN;
try
    info = memory;
    value = info.MemUsedMATLAB;
catch
end
end
