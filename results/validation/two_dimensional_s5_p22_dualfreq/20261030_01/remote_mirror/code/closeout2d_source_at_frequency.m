function src = closeout2d_source_at_frequency(frozen_src,frequency_hz,c0_mps)
%CLOSEOUT2D_SOURCE_AT_FREQUENCY Hold the frozen Tx field; recompute only kz.
% Validation-only candidate for the fixed 2-D Gaussian-source S5 study.
arguments
    frozen_src (1,1) struct
    frequency_hz (1,1) double {mustBeFinite,mustBePositive}
    c0_mps (1,1) double {mustBeFinite,mustBePositive}
end
required={'kx','kz','coeff','n_original','x0_m','z_tx_m'};
assert(all(isfield(frozen_src,required)),'Frozen source schema mismatch');
assert(numel(frozen_src.kx)==frozen_src.n_original && ...
    numel(frozen_src.coeff)==frozen_src.n_original);
src=frozen_src;
k0=2*pi*frequency_hz/c0_mps;
src.kz=sqrt(complex(k0^2-src.kx(:).^2,0));
src.frequency_hz=frequency_hz;
src.c0_mps=c0_mps;
src.s5_source_semantics='frozen Tx spatial field and Fourier coefficients; kz recalculated';
end
