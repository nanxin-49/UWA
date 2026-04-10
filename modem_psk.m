function varargout = modem_psk(mode, varargin)
%MODEM_PSK Utility wrapper for PSK modulation, demodulation and error metrics.

if nargin < 1
    error('Usage: modem_psk(mode, ...)');
end

switch lower(mode)
    case 'modulate'
        [symbols, bits_used] = local_psk_modulate(varargin{:});
        varargout{1} = symbols;
        varargout{2} = bits_used;
    case 'demodulate'
        bits_hat = local_psk_demodulate_coherent(varargin{:});
        varargout{1} = bits_hat;
    case 'error_rate'
        [ber, ser] = local_psk_error_rate(varargin{:});
        varargout{1} = ber;
        varargout{2} = ser;
    otherwise
        error('Unsupported mode: %s', mode);
end
end

function [symbols, bits_used] = local_psk_modulate(bits_in, M)
%LOCAL_PSK_MODULATE Map input bits to M-PSK symbols with unit average power.

if nargin < 2
    M = 4;
end

if ~isscalar(M) || M < 2 || abs(round(log2(M)) - log2(M)) > 1e-12
    error('M must be a power of 2 and >= 2.');
end

k = round(log2(M));
bits_in = bits_in(:);
n_sym = floor(numel(bits_in) / k);
if n_sym < 1
    error('Not enough bits for one symbol.');
end

bits_used = bits_in(1:(n_sym*k));
bit_mat = reshape(bits_used, k, []).';
sym_idx = local_bits_to_index(bit_mat);

phase = 2*pi*double(sym_idx)/M;
symbols = exp(1i*phase);
end

function bits_hat = local_psk_demodulate_coherent(rx_symbols, M)
%LOCAL_PSK_DEMODULATE_COHERENT Coherent hard-decision demod for M-PSK symbols.

if nargin < 2
    M = 4;
end
if ~isscalar(M) || M < 2 || abs(round(log2(M)) - log2(M)) > 1e-12
    error('M must be a power of 2 and >= 2.');
end

k = round(log2(M));
rx_phase = mod(angle(rx_symbols), 2*pi);
idx_hat = mod(round((rx_phase/(2*pi))*M), M);
bit_mat_hat = local_index_to_bits(idx_hat, k);
bits_hat = reshape(bit_mat_hat.', [], 1);
end

function [ber, ser] = local_psk_error_rate(bits_tx, bits_rx, M)
%LOCAL_PSK_ERROR_RATE Compute BER and SER for M-PSK.

if nargin < 3
    M = 4;
end
k = round(log2(M));

n_bits = min(numel(bits_tx), numel(bits_rx));
bits_tx = bits_tx(1:n_bits);
bits_rx = bits_rx(1:n_bits);

ber = mean(bits_tx ~= bits_rx);

n_sym = floor(n_bits / k);
if n_sym < 1
    ser = NaN;
    return
end

sym_err = false(n_sym, 1);
for ii = 1:n_sym
    s = (ii-1)*k + 1;
    e = ii*k;
    sym_err(ii) = any(bits_tx(s:e) ~= bits_rx(s:e));
end
ser = mean(sym_err);
end

function idx = local_bits_to_index(bit_mat)
%LOCAL_BITS_TO_INDEX Convert binary row vectors (left-msb) to decimal index.

[n_row, n_col] = size(bit_mat);
weights = 2.^((n_col-1):-1:0);
idx = zeros(n_row, 1);
for ii = 1:n_row
    idx(ii) = sum(double(bit_mat(ii, :)) .* weights);
end
end

function bit_mat = local_index_to_bits(idx, k)
%LOCAL_INDEX_TO_BITS Convert decimal index to binary row vectors (left-msb).

idx = double(idx(:));
n = numel(idx);
bit_mat = zeros(n, k);
for bb = 1:k
    shift = k - bb;
    bit_mat(:, bb) = mod(floor(idx / (2^shift)), 2);
end
end
