function [projection, meta] = build_adjoint_receiver_projection_vertical(cache, pm_spec)
%BUILD_ADJOINT_RECEIVER_PROJECTION_VERTICAL Build one exact kernel per frequency.
%   Validation-only v1 for the fixed cached uniform CPU-double path.

arguments
    cache (1,1) struct
    pm_spec (1,1) struct
end

required_cache = {'cfg','f_axis_hz','ix_rx','iy_rx','incident_surface_xy_f', ...
    'H_direct_f','H_ref_coh_f','pm_mapping'};
for ii = 1:numel(required_cache)
    if ~isfield(cache, required_cache{ii})
        error('build_adjoint_receiver_projection_vertical:InvalidCache', ...
            'cache.%s is required.', required_cache{ii});
    end
end
required_pm = {'nx','ny','dx_m','dy_m'};
for ii = 1:numel(required_pm)
    if ~isfield(pm_spec, required_pm{ii})
        error('build_adjoint_receiver_projection_vertical:InvalidPMSpec', ...
            'pm_spec.%s is required.', required_pm{ii});
    end
end
local_validate_scope(cache, pm_spec);

F = numel(cache.f_axis_hz);
q_surface_xy_f = complex(zeros(cache.cfg.ny, cache.cfg.nx, F));
a_pe_xy_f = complex(zeros(cache.cfg.ny, cache.cfg.nx, F));
per_frequency_s = zeros(F,1);
receiver_source = complex(zeros(cache.cfg.ny, cache.cfg.nx));
receiver_source(cache.iy_rx, cache.ix_rx) = 1;

total_timer = tic;
for ii = 1:F
    frequency_timer = tic;
    q_surface_xy_f(:,:,ii) = apply_adjoint_receiver_to_surface_vertical( ...
        cache, receiver_source, ii);
    a_pe_xy_f(:,:,ii) = conj(cache.incident_surface_xy_f(:,:,ii)) .* ...
        q_surface_xy_f(:,:,ii);
    per_frequency_s(ii) = toc(frequency_timer);
end
total_s = toc(total_timer);

projection = struct();
projection.kind = 'adjoint_pe_receiver_projection_uniform_v1';
projection.schema_version = '2.0.0';
projection.f_axis_hz = cache.f_axis_hz(:);
projection.q_surface_xy_f = q_surface_xy_f;
projection.a_pe_xy_f = a_pe_xy_f;
projection.H_direct_reduced_f = local_cache_field(cache, ...
    'H_direct_reduced_f','H_direct_f');
projection.H_ref_coh_reduced_f = local_cache_field(cache, ...
    'H_ref_coh_reduced_f','H_ref_coh_f');
projection.mu_scatter_f = complex(zeros(F,1));
projection.mu_total_reduced_f = projection.H_direct_reduced_f + ...
    projection.H_ref_coh_reduced_f;
projection.pm_mapping = cache.pm_mapping;
projection.receiver_grid_index = [cache.iy_rx, cache.ix_rx];
projection.c0_mps = cache.cfg.c0;
projection.reflect_coeff = cache.cfg.reflect_coeff;
projection.phase_geometry = struct('z_tx',cache.cfg.z_tx, ...
    'z_rx',cache.cfg.z_rx,'z_surface',0,'c0',cache.cfg.c0);
[deterministic_dsp,phase_meta] = apply_pe_channel_phase_reference_vertical( ...
    projection.f_axis_hz,struct('direct_f',projection.H_direct_reduced_f, ...
    'reflect_coh_f',projection.H_ref_coh_reduced_f), ...
    projection.phase_geometry,'direct_dsp');
projection.H_direct_dsp_f = deterministic_dsp.direct_f;
projection.H_ref_coh_dsp_f = deterministic_dsp.reflect_coh_f;
projection.H_direct_f = projection.H_direct_dsp_f;
projection.H_ref_coh_f = projection.H_ref_coh_dsp_f;
projection.mu_total_f = projection.H_direct_dsp_f+projection.H_ref_coh_dsp_f;
projection.phase_reference_meta = phase_meta;
projection.deltaG_definition = ...
    'deltaG_i=R0_i*exp(1i*alpha_i*eta)-R_coh_i; C_deltaG and P_deltaG include R0';
projection.weight_definition = 'a_pe_i=conj(psi_inc_i).*q_i; no extra R0 factor';
projection.mean_definition = ...
    'mu_total_f=H_direct_f+H_ref_coh_f+mu_scatter_f; mu_scatter_f=0';
projection.limitations = ['Uniform sound speed, CPU double, fixed transmitter and one ', ...
    'nearest-grid receiver, fixed PE/PM grids and frequency axis, no bubbles or Doppler.'];

array_bytes = numel(q_surface_xy_f)*16 + numel(a_pe_xy_f)*16 + ...
    numel(projection.H_direct_f)*16 + numel(projection.H_ref_coh_f)*16;
meta = struct('per_frequency_s',per_frequency_s, 'total_s',total_s, ...
    'mean_per_frequency_s',mean(per_frequency_s), ...
    'projection_array_bytes',array_bytes, ...
    'receiver_source_nonzero_count',nnz(receiver_source), ...
    'receiver_sampling','nearest_grid_point', ...
    'adjoint_kind','exact_discrete_conjugate_transpose', ...
    'memory_snapshot_bytes',local_memory_used_bytes());
projection.build_meta = meta;
end

function value = local_cache_field(cache,preferred,legacy)
if isfield(cache,preferred)
    value = cache.(preferred)(:);
else
    value = cache.(legacy)(:);
end
end

function local_validate_scope(cache, pm_spec)
if isfield(cache.cfg,'env_mode') && ~strcmpi(cache.cfg.env_mode,'uniform')
    error('build_adjoint_receiver_projection_vertical:UnsupportedEnvironment', ...
        'v1 requires env_mode=uniform.');
end
if isfield(cache.cfg,'enable_bubbles') && cache.cfg.enable_bubbles
    error('build_adjoint_receiver_projection_vertical:BubblesUnsupported', ...
        'v1 requires bubbles disabled.');
end
if ~isequal(cache.pm_mapping.pm_grid, [pm_spec.ny, pm_spec.nx]) || ...
        ~isequal(cache.pm_mapping.pe_grid, [cache.cfg.ny, cache.cfg.nx])
    error('build_adjoint_receiver_projection_vertical:MappingMismatch', ...
        'pm_spec and cached PM-to-PE mapping dimensions do not match.');
end
if abs(pm_spec.dx_m-cache.pm_mapping.dx_m) > 1e-12 || ...
        abs(pm_spec.dy_m-cache.cfg.yw/cache.cfg.ny) > 1e-12
    error('build_adjoint_receiver_projection_vertical:SpacingMismatch', ...
        'PM and PE spacing must match the cached central-crop mapping.');
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
