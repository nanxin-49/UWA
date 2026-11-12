function [rx_noisy, noise, meta] = noise_inject_vertical(rx_clean, noise_cfg, signal_ref)
%NOISE_INJECT_VERTICAL Add configurable noise to baseband received symbols.
%
% Signature:
%   [rx_noisy, noise, meta] = noise_inject_vertical(rx_clean, noise_cfg, signal_ref)
%
% noise_cfg fields:
%   enable_noise    : logical, default false
%   model           : 'awgn' (default)
%   snr_db          : optional SNR in dB
%   ebn0_db         : optional Eb/N0 in dB
%   bits_per_symbol : used with ebn0_db, default 1
%   seed            : RNG seed, default 12345
%   custom_noise_fn : optional function handle for custom models

if nargin < 1
    error('rx_clean is required.');
end
if nargin < 2 || isempty(noise_cfg)
    noise_cfg = struct();
end
if nargin < 3 || isempty(signal_ref)
    signal_ref = rx_clean;
end

if ~isnumeric(rx_clean)
    error('rx_clean must be numeric.');
end
if ~isnumeric(signal_ref)
    error('signal_ref must be numeric.');
end

if ~isfield(noise_cfg, 'enable_noise') || isempty(noise_cfg.enable_noise)
    noise_cfg.enable_noise = false;
end
if ~isfield(noise_cfg, 'model') || isempty(noise_cfg.model)
    noise_cfg.model = 'awgn';
end
if ~isfield(noise_cfg, 'snr_db')
    noise_cfg.snr_db = [];
end
if ~isfield(noise_cfg, 'ebn0_db')
    noise_cfg.ebn0_db = [];
end
if ~isfield(noise_cfg, 'bits_per_symbol') || isempty(noise_cfg.bits_per_symbol)
    noise_cfg.bits_per_symbol = 1;
end
if ~isfield(noise_cfg, 'seed') || isempty(noise_cfg.seed)
    noise_cfg.seed = 12345;
end
if ~isfield(noise_cfg, 'custom_noise_fn')
    noise_cfg.custom_noise_fn = [];
end

if isstring(noise_cfg.model)
    noise_cfg.model = char(noise_cfg.model);
end
if ~ischar(noise_cfg.model)
    error('noise_cfg.model must be a string or char.');
end

meta = struct();
meta.enabled = logical(noise_cfg.enable_noise);
meta.model = lower(noise_cfg.model);
meta.seed = noise_cfg.seed;
meta.target_snr_db = NaN;
meta.effective_snr_db = Inf;
meta.noise_power = 0;
meta.signal_power_ref = mean(abs(signal_ref(:)).^2);

if ~meta.enabled
    noise = complex(zeros(size(rx_clean)));
    rx_noisy = rx_clean;
    return
end

if ~isempty(noise_cfg.custom_noise_fn)
    if ~isa(noise_cfg.custom_noise_fn, 'function_handle')
        error('custom_noise_fn must be a function handle.');
    end
    [rx_noisy, noise, custom_meta] = noise_cfg.custom_noise_fn(rx_clean, noise_cfg, signal_ref);
    if ~isnumeric(rx_noisy) || ~isequal(size(rx_noisy), size(rx_clean))
        error('custom_noise_fn must return rx_noisy with same size as rx_clean.');
    end
    if ~isnumeric(noise) || ~isequal(size(noise), size(rx_clean))
        error('custom_noise_fn must return noise with same size as rx_clean.');
    end
    if exist('custom_meta', 'var')
        meta.custom = custom_meta;
    end
    n_pow = mean(abs(noise(:)).^2);
    if n_pow > 0
        meta.effective_snr_db = 10*log10(max(mean(abs(rx_clean(:)).^2), eps) / n_pow);
    end
    return
end

switch meta.model
    case 'awgn'
        if ~isempty(noise_cfg.snr_db)
            snr_target_db = noise_cfg.snr_db;
        elseif ~isempty(noise_cfg.ebn0_db)
            snr_target_db = noise_cfg.ebn0_db + 10*log10(noise_cfg.bits_per_symbol);
        else
            error('noise_cfg.snr_db or noise_cfg.ebn0_db must be provided when noise is enabled.');
        end

        signal_power = mean(abs(signal_ref(:)).^2);
        if signal_power <= 0
            error('signal_ref power must be positive for AWGN injection.');
        end

        noise_power = signal_power / (10^(snr_target_db/10));
        rng(noise_cfg.seed, 'twister');
        noise = sqrt(noise_power/2) .* ...
                (randn(size(rx_clean)) + 1i*randn(size(rx_clean)));
        rx_noisy = rx_clean + noise;

        meta.target_snr_db = snr_target_db;
        meta.noise_power = noise_power;
        meta.effective_snr_db = 10*log10(max(mean(abs(rx_clean(:)).^2), eps) / ...
                                         max(mean(abs(noise(:)).^2), eps));
    otherwise
        error('Unsupported noise model: %s', meta.model);
end

end
