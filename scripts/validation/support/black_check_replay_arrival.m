function evidence = black_check_replay_arrival(y,baseband,channel,cfg,start)
%BLACK_CHECK_REPLAY_ARRIVAL Block matched-filter arrival vs stored taps/drift.
% Eight disjoint probe blocks reduce sensitivity to phase change over the probe.
fs_delay = channel.params.fs_delay;
[p,q] = rat(fs_delay/cfg.fs_passband);
t = (0:size(y,1)-1).'/cfg.fs_passband;
rx = resample(y.*exp(-2j*pi*channel.params.fc*t),p,q);
tx = resample(baseband,p,q);
B = 8;
L = floor(numel(tx)/B);
profile_length = size(rx,1)-B*L+1;
assert(profile_length>0,'black:CorrelationLength','Invalid matched-filter support.');
profiles = zeros(profile_length,size(y,2));
for bb = 0:B-1
    piece = tx(bb*L+(1:L));
    nfft = 2^nextpow2(size(rx,1)+L-1);
    kernel_f = fft(conj(flipud(piece)),nfft);
    for mm = 1:size(y,2)
        correlation = ifft(fft(rx(:,mm),nfft).*kernel_f);
        power = abs(correlation(L:size(rx,1))).^2;
        profiles(:,mm) = profiles(:,mm)+power(bb*L+(1:profile_length))/B;
    end
end
time_first = (start-1)/fs_delay;
time_last = time_first+numel(tx)/fs_delay+2*cfg.guard_s;
snapshot_time = (0:size(channel.h_hat,3)-1)/channel.params.fs_time;
snapshots = snapshot_time>=time_first & snapshot_time<=time_last;
pdp = mean(abs(channel.h_hat(:,:,snapshots)).^2,3);
expected = zeros(1,size(y,2));
observed = zeros(size(expected));
drift = zeros(size(expected));
if isfield(channel,'phi_hat')
    ids = start:min(size(channel.phi_hat,2),start+ceil((time_last-time_first)*fs_delay));
    drift = mean(channel.phi_hat(:,ids),2).'/(2*pi*channel.params.fc)*fs_delay;
end
if isfield(channel,'f_resamp')
    assert(abs(channel.f_resamp-1)<1e-12,'black:ArrivalBulk', ...
        'Nonunit f_resamp requires an explicit arrival-time coordinate adjustment.');
end
for mm = 1:numel(expected)
    [~,k] = max(pdp(:,mm));
    expected(mm) = cfg.guard_s*fs_delay + local_peak(pdp(:,mm),k) - drift(mm);
    [~,k] = max(profiles(:,mm));
    observed(mm) = local_peak(profiles(:,mm),k);
end
error_samples = observed-expected;
evidence = struct('passed',all(abs(error_samples)<=2),'tolerance_samples',2, ...
    'expected_samples',expected,'observed_samples',observed,'error_samples',error_samples, ...
    'error_us',error_samples/fs_delay*1e6,'drift_samples',drift, ...
    'profile_delay_s',(0:profile_length-1).'/fs_delay,'matched_filter_power',profiles, ...
    'method','8-block incoherent matched filter; sub-sample parabolic peak; stored PDP and mean phi drift');
fprintf('Arrival errors (fs_delay samples): %s\n',mat2str(error_samples,5));
end

function sample = local_peak(p,k)
sample = k-1;
if k>1 && k<numel(p)
    curvature = p(k-1)-2*p(k)+p(k+1);
    if curvature~=0, sample=sample+0.5*(p(k-1)-p(k+1))/curvature; end
end
end
