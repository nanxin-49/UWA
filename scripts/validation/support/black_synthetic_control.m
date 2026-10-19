function evidence = black_synthetic_control()
%BLACK_SYNTHETIC_CONTROL Independent analytic checks of units and phase sign.
test_channel.params = struct('fs_delay',10000,'fs_time',125,'fc',18000);
test_channel.h_hat = zeros(41,1,4);
test_channel.h_hat(1,1,:) = 1;
test_channel.h_hat(41,1,:) = 0.5;
[metrics,~] = black_delay_metrics(test_channel.h_hat,10000);
% Powers 1 and 1/4, path separation 4 ms: mean .8 ms and RMS 1.6 ms.
assert(abs(metrics.mean_delay_s-0.0008) < 1e-14, 'black:ControlMean','Delay mean control failed.');
assert(abs(metrics.rms_delay_s-0.0016) < 1e-14, 'black:ControlRms','Delay RMS control failed.');
assert(abs(metrics.energy90_width_s-0.004) < 1e-14, 'black:ControlWindow','Delay interval control failed.');
t = (0:9999)/10000;
test_channel.phi_hat = 2*pi*(-3.25)*t + 1.2;
phase = black_phase_doppler(test_channel);
assert(max(abs(phase.frequency_hz+3.25)) < 1e-10, 'black:ControlDoppler','Doppler sign/units failed.');
test_channel.f_resamp = 1/(1+2/18000);
phase_bulk = black_phase_doppler(test_channel);
assert(max(abs(phase_bulk.frequency_hz+1.25)) < 1e-9, 'black:ControlBulk','Bulk Doppler control failed.');
tone = exp(2j*pi*10*(0:999).'/125);
[power,frequency] = pwelch(tone,hann(500,'periodic'),250,500,125,'centered');
[~,index] = max(power);
assert(abs(frequency(index)-10)<1e-12,'black:ControlSpectrum','Centered Welch sign/units failed.');
evidence = struct('passed',true,'mean_delay_s',metrics.mean_delay_s, ...
    'rms_delay_s',metrics.rms_delay_s,'known_frequency_hz',-3.25, ...
    'frequency_error_hz',max(abs(phase.frequency_hz+3.25)),'welch_peak_hz',frequency(index));
end
