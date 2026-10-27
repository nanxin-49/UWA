function receipt = validate_two_dimensional_s5_nine_point_cir_controls(output_dir)
%VALIDATE_TWO_DIMENSIONAL_S5_NINE_POINT_CIR_CONTROLS Analytic CIR controls.
arguments
    output_dir (1,:) char
end
assert(~isfolder(output_dir) && ~isfile(output_dir), ...
    'Use a fresh result directory.');

f = (4000+(0:63)'*62.5);
N = numel(f);
df = f(2)-f(1);
dt = 1/(N*df);
c0 = 1500;
ztx = 100;
xv = [-0.1953125,0,0.1953125];
zv = [2.875,3,3.125];
[xg,zg] = meshgrid(xv,zv);
xrx = xg(:);
zrx = zg(:);
ids = strings(9,1);
for iz = 1:3
    for ix = 1:3
        ids((iz-1)*3+ix) = sprintf('P%d%d',iz,ix);
    end
end
range_d_m = hypot(xrx,ztx-zrx);
range_r_m = hypot(xrx,ztx+zrx);
tau_d = range_d_m/c0;
tau_r = range_r_m/c0;
tau_ref = min(tau_d);
arrival_d = tau_d-tau_ref;
arrival_r = tau_r-tau_ref;
assert(max(arrival_r)<1/df);

alpha = [-1,0.3*exp(1i*0.6)];
case_names = ["flat_image_path_control","known_dual_path_control"];
Hd = complex(zeros(N,9,2));
Hr = Hd;
Ht = Hd;
hd = Hd;
hr = Hd;
ht = Hd;
expected_hd = Hd;
expected_hr = Hd;
expected_ht = Hd;
names = strings(18,1);
receiver_ids = strings(18,1);
numeric = zeros(18,13);

for ic = 1:2
    for ir = 1:9
        row = (ic-1)*9+ir;
        pd = exp(1i*2*pi*f*tau_d(ir));
        pr = alpha(ic)*exp(1i*2*pi*f*tau_r(ir));
        Hd(:,ir,ic) = conj(pd).*exp(1i*2*pi*f*tau_ref);
        Hr(:,ir,ic) = conj(pr).*exp(1i*2*pi*f*tau_ref);
        Ht(:,ir,ic) = Hd(:,ir,ic)+Hr(:,ir,ic);
        cd = build_channel_cir_vertical(Hd(:,ir,ic),f, ...
            struct('input_reference','direct_dsp'));
        cr = build_channel_cir_vertical(Hr(:,ir,ic),f, ...
            struct('input_reference','direct_dsp'));
        ct = build_channel_cir_vertical(Ht(:,ir,ic),f, ...
            struct('input_reference','direct_dsp'));
        ed = local_expected_cir(1,arrival_d(ir),f(1),df,N);
        er = local_expected_cir(conj(alpha(ic)),arrival_r(ir),f(1),df,N);
        et = ed+er;
        hd(:,ir,ic) = cd.h_tau;
        hr(:,ir,ic) = cr.h_tau;
        ht(:,ir,ic) = ct.h_tau;
        expected_hd(:,ir,ic) = ed;
        expected_hr(:,ir,ic) = er;
        expected_ht(:,ir,ic) = et;
        id = round(arrival_d(ir)/dt)+1;
        irf = round(arrival_r(ir)/dt)+1;
        [~,ipd] = max(abs(cd.h_tau));
        [~,ipr] = max(abs(cr.h_tau));
        numeric(row,:) = [ ...
            norm(conj(Hd(:,ir,ic)).*exp(1i*2*pi*f*tau_ref)-pd)/norm(pd), ...
            norm(conj(Hr(:,ir,ic)).*exp(1i*2*pi*f*tau_ref)-pr)/norm(pr), ...
            norm(Ht(:,ir,ic)-Hd(:,ir,ic)-Hr(:,ir,ic))/norm(Ht(:,ir,ic)), ...
            norm(fft(cd.h_tau)-Hd(:,ir,ic))/norm(Hd(:,ir,ic)), ...
            norm(fft(cr.h_tau)-Hr(:,ir,ic))/norm(Hr(:,ir,ic)), ...
            norm(fft(ct.h_tau)-Ht(:,ir,ic))/norm(Ht(:,ir,ic)), ...
            norm(cd.h_tau-ed)/norm(ed),norm(cr.h_tau-er)/norm(er), ...
            norm(ct.h_tau-et)/norm(et),abs((ipd-1)*dt-arrival_d(ir)), ...
            abs((ipr-1)*dt-arrival_r(ir)),abs(cd.h_tau(id)-ed(id)), ...
            abs(cr.h_tau(irf)-er(irf))];
        names(row) = case_names(ic);
        receiver_ids(row) = ids(ir);
    end
end

metrics = table(names,receiver_ids,numeric(:,1),numeric(:,2),numeric(:,3), ...
    numeric(:,4),numeric(:,5),numeric(:,6),numeric(:,7),numeric(:,8), ...
    numeric(:,9),numeric(:,10),numeric(:,11),numeric(:,12),numeric(:,13), ...
    'VariableNames',{'case_name','receiver_id','mapping_direct_relative', ...
    'mapping_reflect_relative','component_closure_relative','fft_direct_relative', ...
    'fft_reflect_relative','fft_total_relative','cir_direct_relative', ...
    'cir_reflect_relative','cir_total_relative','direct_peak_error_s', ...
    'reflect_peak_error_s','direct_amplitude_error','reflect_amplitude_error'});
max_complex_error = max(numeric(:,1:9),[],'all');
max_arrival_error_s = max(numeric(:,10:11),[],'all');
max_amplitude_error = max(numeric(:,12:13),[],'all');
gate_complex = max_complex_error<=1e-10;
gate_arrival = max_arrival_error_s<=dt/2;
gate_amplitude = max_amplitude_error<=1e-10;
all_gates_pass = gate_complex && gate_arrival && gate_amplitude;
if all_gates_pass
    status = 'NINE_POINT_ANALYTIC_CIR_CONTROLS_PASS';
else
    status = 'NINE_POINT_ANALYTIC_CIR_CONTROLS_FAIL';
end
cfg = struct('frequency_hz',f,'receiver_x_m',xrx,'receiver_z_m',zrx, ...
    'tx_depth_m',ztx,'sound_speed_mps',c0,'tau_ref_s',tau_ref, ...
    'frequency_count',N,'delta_f_hz',df,'delay_sample_spacing_s',dt, ...
    'control_model','normalized point paths with exact direct/image ranges', ...
    'direct_path_gain',1, ...
    'reflection_coefficients',alpha,'native_to_dsp', ...
    'conj(P_native)*exp(i*2*pi*f*tau_ref)', ...
    'cir_window','none','zero_padding_factor',1);

receipt = struct('status',status, ...
    'matlab_version',version,'acoustic_solver_calls',0, ...
    'frequency_min_hz',f(1),'frequency_max_hz',f(end), ...
    'frequency_step_hz',df,'frequency_count',N,'delay_window_s',1/df, ...
    'delay_resolution_s',1/(f(end)-f(1)),'delay_sample_spacing_s',dt, ...
    'tau_ref_s',tau_ref,'receiver_count',9,'analytic_case_count',2, ...
    'max_complex_relative_error',max_complex_error, ...
    'max_arrival_error_s',max_arrival_error_s, ...
    'max_amplitude_error',max_amplitude_error, ...
    'failed_check_count',3-(gate_complex+gate_arrival+gate_amplitude), ...
    'gates',struct('complex_error_pass',gate_complex, ...
    'arrival_tolerance_pass',gate_arrival,'amplitude_pass',gate_amplitude));
mkdir(output_dir);
save(fullfile(output_dir,'nine_point_cir_controls.mat'),'f','xrx','zrx', ...
    'ids','range_d_m','range_r_m','tau_d','tau_r','tau_ref','arrival_d','arrival_r','alpha', ...
    'Hd','Hr','Ht','hd','hr','ht','expected_hd','expected_hr', ...
    'expected_ht','metrics','receipt','cfg');
writetable(metrics,fullfile(output_dir,'nine_point_cir_metrics.csv'));
fid = fopen(fullfile(output_dir,'receipt.json'),'w');
assert(fid>0,'Cannot create receipt.json.');
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s',jsonencode(receipt,PrettyPrint=true));
disp(receipt);
if ~all_gates_pass
    error('Nine-point CIR controls failed; the saved result directory contains the metrics.');
end
end

function h = local_expected_cir(a,tau,f0,df,N)
n = (0:N-1).';
k = 0:N-1;
h = a*exp(-1i*2*pi*f0*tau).* ...
    (sum(exp(1i*2*pi*(n/N-tau*df).*k),2)/N);
end
