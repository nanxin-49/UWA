function receipt = validate_two_dimensional_s5_prerequisites(output_dir)
%VALIDATE_TWO_DIMENSIONAL_S5_PREREQUISITES Solver-free phase/CIR controls.
% This does not qualify any PE, BIE, or FullFK broadband field.
arguments
    output_dir (1,:) char
end
assert(~isfolder(output_dir) && ~isfile(output_dir),'Use a fresh result directory');
mkdir(output_dir);
f=(4000+(0:63)'*62.5); N=numel(f);
tau_d=97/1500; tau_r=103/1500; dt=tau_r-tau_d;
g=0.3*exp(1i*0.6);
H_direct=ones(N,1);
H_reflect=g*exp(-1i*2*pi*f*dt);
H_total=H_direct+H_reflect;
P_direct=conj(H_direct).*exp(1i*2*pi*f*tau_d);
P_reflect=conj(H_reflect).*exp(1i*2*pi*f*tau_d);
D=conj(P_direct).*exp(1i*2*pi*f*tau_d);
R=conj(P_reflect).*exp(1i*2*pi*f*tau_d);
mapping_error=norm([D-H_direct;R-H_reflect])/norm([H_direct;H_reflect]);
roundtrip_error=norm([conj(D).*exp(1i*2*pi*f*tau_d)-P_direct; ...
    conj(R).*exp(1i*2*pi*f*tau_d)-P_reflect])/norm([P_direct;P_reflect]);
algebra_error=norm(H_total-D-R)/norm(H_total);
cir=build_channel_cir_vertical(H_total,f,struct('input_reference','direct_dsp'));
h=cir.h_tau;
fft_error=norm(fft(h)-H_total)/norm(H_total);
expected=complex(zeros(N,1)); expected(1)=1; expected(17)=g;
cir_error=norm(h-expected)/norm(expected);
delay_error_s=abs(cir.delay_axis_s(17)-dt);
amplitude_error=max(abs([h(1)-1;h(17)-g]));
H_plane=H_direct-exp(-1i*2*pi*f*dt);
plane_cir=build_channel_cir_vertical(H_plane,f, ...
    struct('input_reference','direct_dsp'));
h_plane=plane_cir.h_tau;
expected_plane=complex(zeros(N,1)); expected_plane(1)=1;
expected_plane(17)=-1;
plane_cir_error=norm(h_plane-expected_plane)/norm(expected_plane);
assert(mapping_error<=1e-10 && roundtrip_error<=1e-10 && ...
    algebra_error<=1e-10 && fft_error<=1e-10 && cir_error<=1e-10);
assert(delay_error_s<=1e-12 && amplitude_error<=1e-10 && ...
    plane_cir_error<=1e-10);
metrics=table(mapping_error,roundtrip_error,algebra_error,fft_error, ...
    cir_error,plane_cir_error,delay_error_s,amplitude_error, ...
    'VariableNames',{'native_dsp_relative','native_roundtrip_relative', ...
    'component_sum_relative','fft_roundtrip_relative','two_path_cir_relative', ...
    'plane_cir_relative','two_path_delay_error_s','two_path_amplitude_error'});
receipt=struct('status','ALGEBRA_CONTROLS_PASS_ONLY','matlab_version',version, ...
    'frequency_min_hz',f(1),'frequency_max_hz',f(end), ...
    'frequency_step_hz',62.5,'independent_bins',N, ...
    'delay_window_s',1/62.5,'delay_bin_s',1/(N*62.5), ...
    'native_convention','exp(-i omega t)','dsp_formula', ...
    'conj(P_native)*exp(i 2 pi f tau_direct); shared origin', ...
    'acoustic_solver_calls',0,'metrics',table2struct(metrics));
save(fullfile(output_dir,'algebra_controls.mat'),'f','H_direct', ...
    'H_reflect','H_total','H_plane','P_direct','P_reflect', ...
    'h','h_plane','expected','expected_plane','metrics','receipt');
writetable(metrics,fullfile(output_dir,'algebra_metrics.csv'));
fid=fopen(fullfile(output_dir,'receipt.json'),'w');
assert(fid>0); cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'%s',jsonencode(receipt,PrettyPrint=true));
disp(receipt);
end
