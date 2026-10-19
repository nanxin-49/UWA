function output = black_phase_doppler(channel)
%BLACK_PHASE_DOPPLER Fit continuous phase slope in non-overlapping 0.1 s bins.
fs_delay = channel.params.fs_delay;
M = size(channel.h_hat,2);
if isfield(channel,'phi_hat')
    field = 'phi_hat';
elseif isfield(channel,'theta_hat')
    field = 'theta_hat';
else
    output = struct('available',false,'field','none','time_s',[], ...
        'frequency_hz',zeros(0,M),'bulk_frequency_hz',NaN,'block_s',0.1);
    return
end
phase = channel.(field);
assert(size(phase,1)==M && all(isfinite(phase(:))), 'black:PhaseShape', 'Invalid phase trajectory.');
N = round(0.1*fs_delay);
assert(N >= 3,'black:PhaseRate','Too few phase samples per block.');
B = floor(size(phase,2)/N);
assert(B >= 1,'black:PhaseDuration','Phase record too short.');
t = (0:N-1).'/fs_delay;
t = t-mean(t);
frequency_hz = zeros(B,M);
for mm = 1:M
    blocks = reshape(phase(mm,1:N*B),N,B);
    % Centering the values improves conditioning without wrapping the phase.
    blocks = blocks - mean(blocks,1);
    frequency_hz(:,mm) = (t.'*blocks).' / sum(t.^2) / (2*pi);
end
bulk = 0;
if isfield(channel,'f_resamp')
    bulk = (1/channel.f_resamp-1)*channel.params.fc;
    frequency_hz = frequency_hz + bulk;
end
output = struct('available',true,'field',field, ...
    'time_s',((0:B-1)*N+(N-1)/2).'/fs_delay,'frequency_hz',frequency_hz, ...
    'bulk_frequency_hz',bulk,'block_s',N/fs_delay);
end
