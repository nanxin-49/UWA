function psi_rx_xy_m = apply_forward_surface_to_receiver_vertical( ...
    cache, psi_surface_xy_m, frequency_index)
%APPLY_FORWARD_SURFACE_TO_RECEIVER_VERTICAL Cached uniform PE plane operator.
%   Validation-only CPU-double implementation. The output is the complete
%   receiver-depth plane; receiver sampling remains the caller's operation.

arguments
    cache (1,1) struct
    psi_surface_xy_m {mustBeNumeric}
    frequency_index (1,1) double {mustBeInteger,mustBePositive}
end

local_validate_cache(cache, frequency_index);
if isa(psi_surface_xy_m, 'gpuArray') || ~isa(psi_surface_xy_m, 'double')
    error('apply_forward_surface_to_receiver_vertical:CPUdoubleRequired', ...
        'psi_surface_xy_m must be a CPU double array.');
end
sz = size(psi_surface_xy_m);
if numel(sz) < 3
    sz(3) = 1;
end
if sz(1) ~= cache.cfg.ny || sz(2) ~= cache.cfg.nx || numel(sz) > 3
    error('apply_forward_surface_to_receiver_vertical:GridMismatch', ...
        'Input must have size ny-by-nx-by-M on the cached PE grid.');
end

fr = cache.surface_rx_fr(:,:,frequency_index);
screen = cache.surface_rx_screen(:,:,frequency_index);
n_step = cache.surface_rx_nstep(frequency_index);
psi_k = fft2(psi_surface_xy_m);
for step_index = 1:n_step
    psi_k = fr .* fft2(screen .* ifft2(fr .* psi_k));
end
psi_rx_xy_m = ifft2(psi_k);
end

function local_validate_cache(cache, frequency_index)
required = {'cfg','f_axis_hz','surface_rx_fr','surface_rx_screen','surface_rx_nstep'};
for ii = 1:numel(required)
    if ~isfield(cache, required{ii})
        error('apply_forward_surface_to_receiver_vertical:InvalidCache', ...
            'cache.%s is required.', required{ii});
    end
end
if frequency_index > numel(cache.f_axis_hz)
    error('apply_forward_surface_to_receiver_vertical:FrequencyIndex', ...
        'frequency_index exceeds the cached frequency count.');
end
if ~(isscalar(cache.surface_rx_nstep(frequency_index)) && ...
        isfinite(cache.surface_rx_nstep(frequency_index)) && ...
        cache.surface_rx_nstep(frequency_index) >= 1)
    error('apply_forward_surface_to_receiver_vertical:InvalidStepCount', ...
        'Cached surface-to-receiver step count must be positive.');
end
end
