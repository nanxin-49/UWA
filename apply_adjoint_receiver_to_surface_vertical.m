function psi_surface_xy_m = apply_adjoint_receiver_to_surface_vertical( ...
    cache, psi_rx_xy_m, frequency_index)
%APPLY_ADJOINT_RECEIVER_TO_SURFACE_VERTICAL Exact cached PE discrete adjoint.
%   This is the conjugate transpose of the implemented discrete forward
%   plane operator. It is not reciprocal back-propagation or inversion.

arguments
    cache (1,1) struct
    psi_rx_xy_m {mustBeNumeric}
    frequency_index (1,1) double {mustBeInteger,mustBePositive}
end

local_validate_cache(cache, frequency_index);
if isa(psi_rx_xy_m, 'gpuArray') || ~isa(psi_rx_xy_m, 'double')
    error('apply_adjoint_receiver_to_surface_vertical:CPUdoubleRequired', ...
        'psi_rx_xy_m must be a CPU double array.');
end
sz = size(psi_rx_xy_m);
if numel(sz) < 3
    sz(3) = 1;
end
if sz(1) ~= cache.cfg.ny || sz(2) ~= cache.cfg.nx || numel(sz) > 3
    error('apply_adjoint_receiver_to_surface_vertical:GridMismatch', ...
        'Input must have size ny-by-nx-by-M on the cached PE grid.');
end

fr_h = conj(cache.surface_rx_fr(:,:,frequency_index));
screen_h = conj(cache.surface_rx_screen(:,:,frequency_index));
n_step = cache.surface_rx_nstep(frequency_index);
q_k = fft2(psi_rx_xy_m);
for step_index = n_step:-1:1 %#ok<NASGU>
    q_k = fr_h .* fft2(screen_h .* ifft2(fr_h .* q_k));
end
psi_surface_xy_m = ifft2(q_k);
end

function local_validate_cache(cache, frequency_index)
required = {'cfg','f_axis_hz','surface_rx_fr','surface_rx_screen','surface_rx_nstep'};
for ii = 1:numel(required)
    if ~isfield(cache, required{ii})
        error('apply_adjoint_receiver_to_surface_vertical:InvalidCache', ...
            'cache.%s is required.', required{ii});
    end
end
if frequency_index > numel(cache.f_axis_hz)
    error('apply_adjoint_receiver_to_surface_vertical:FrequencyIndex', ...
        'frequency_index exceeds the cached frequency count.');
end
end
