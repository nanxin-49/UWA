function [metrics, pdp] = black_delay_metrics(h_hat, fs_delay)
%BLACK_DELAY_METRICS Thresholded PDP moments and central 90%% delay interval.
assert(isnumeric(h_hat) && all(isfinite(h_hat(:))), 'black:FiniteTaps', 'Invalid channel taps.');
assert(fs_delay > 0, 'black:DelayRate', 'Positive delay sampling rate required.');
pdp = mean(abs(h_hat).^2, 3);
tau = (0:size(h_hat,1)-1).' / fs_delay;
M = size(h_hat,2);
values = zeros(M,7);
for mm = 1:M
    p = pdp(:,mm);
    assert(sum(p) > 0, 'black:ZeroChannel', 'Zero-energy channel.');
    mask = p >= max(p)*10^(-25/10);
    retained = sum(p(mask))/sum(p);
    p(~mask) = 0;
    p = p/sum(p);
    mean_delay = sum(tau.*p);
    rms_delay = sqrt(sum((tau-mean_delay).^2.*p));
    cumulative = cumsum(p);
    low = tau(find(cumulative >= 0.05,1));
    high = tau(find(cumulative >= 0.95,1));
    [~, main_tap] = max(pdp(:,mm));
    values(mm,:) = [mean_delay,rms_delay,low,high,high-low,retained,main_tap];
end
metrics = array2table([(1:M).', values], 'VariableNames', ...
    {'receiver','mean_delay_s','rms_delay_s','energy90_low_s','energy90_high_s', ...
    'energy90_width_s','retained_energy_fraction','main_tap_index'});
end
