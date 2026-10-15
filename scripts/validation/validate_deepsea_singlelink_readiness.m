function output=validate_deepsea_singlelink_readiness()
%VALIDATE_DEEPSEA_SINGLELINK_READINESS Deterministic controls, never PM inference.
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
outdir=fullfile(root,'results','validation','deepsea_two_track_20261014');
file=fullfile(outdir,'singlelink_readiness.mat');
assert(~isfile(file),'Readiness:Exists','Preserve prior run.');
timer=tic;[before,sys]=memory;oldrng=rng;cleaner=onCleanup(@()rng(oldrng));
data_seed=2610141;noise_seed=2610142;
meta=struct('physics_tier',"analytic_control",'scope',"deterministic controls only", ...
    'G0',"PASS",'G1_target_PM',"OPEN",'WP2_formal',false, ...
    'source',"unit point Green exp(ikR)/(4pi R)",'time_convention',"exp(-i omega t)", ...
    'native_to_absolute_DSP',"conjugation",'DSP_reference',"exp(i2pi f tau_d)", ...
    'c_mps',1500,'tx_m',[0 0 1500],'rx_m',[0 0 5], ...
    'parent_surface',"flat analytic; no random coefficients",'support',"infinite flat image", ...
    'data_seed',data_seed,'noise_seed',noise_seed, ...
    'plan',"reports/deep_sea_channel_research_plan.md", ...
    'version_hash_manifest',"results/validation/deepsea_two_track_20261014/run_manifest.json");
deepsea_readiness_guard(meta,"readiness");
responses=struct([]);rows=struct([]);
for q=[0 .1]
    tx=[0 q*1505 1500];rx=[0 0 5];Ld=norm(tx-rx);Lr=norm([tx(1:2) -tx(3)]-rx);
    tau=Ld/1500;delay=(Lr-Ld)/1500;
    for N=[256 512]
        f=linspace(4000,8000,N+1).';df=4000/N;
        Hd=exp(1i*2*pi*f*tau)/(4*pi*Ld);
        Hr=-exp(1i*2*pi*f*(Lr/1500))/(4*pi*Lr);
        native=Hd+Hr;absolute=conj(native);
        ref=exp(1i*2*pi*f*tau);relative=absolute.*ref;
        direct=conj(Hd).*ref;reflect=conj(Hr).*ref;
        hd=ifft(ifftshift(direct(1:N)));hr=ifft(ifftshift(reflect(1:N)));h=hd+hr;
        recovered=fftshift(fft(h));roundtrip=norm(recovered-relative(1:N))/norm(relative(1:N));
        component=max([norm(fftshift(fft(hd))-direct(1:N))/norm(direct(1:N)), ...
            norm(fftshift(fft(hr))-reflect(1:N))/norm(reflect(1:N))]);
        parseval=abs(sum(abs(h).^2)-sum(abs(relative(1:N)).^2)/N)/sum(abs(h).^2);
        measured=-angle(sum(reflect(2:end).*conj(reflect(1:end-1))))/(2*pi*df);
        % Half-bin frequencies are computed independently from path geometry.
        mid=(f(1:end-1)+f(2:end))/2;
        exact=1/(4*pi*Ld)-exp(-1i*2*pi*mid*delay)/(4*pi*Lr);
        interp=(relative(1:end-1)+relative(2:end))/2;
        interp_error=norm(interp-exact)/norm(ones(N,1)/(4*pi*Ld));
        absolute_formula=conj(Hd)+conj(Hr);
        convention=norm(absolute-absolute_formula)/norm(Hd);
        assert(max([roundtrip,component,parseval,convention])<1e-10);
        assert(abs(measured-delay)<1e-10);
        if q==0,assert(abs(delay-1/150)<1e-14);end
        % Padding time samples changes frequency interpolation only.
        padded=fft([h;zeros(N,1)]);
        zeropad_match=norm(padded(1:2:end)-fft(h))/norm(fft(h));
        assert(zeropad_match<1e-10);
        r=struct('q',q,'N_intervals',N,'N_source_points',N+1,'df_hz',df, ...
            'window_s',1/df,'dt_s',1/4000,'delay_s',delay, ...
            'interference_period_hz',1/delay,'roundtrip',roundtrip, ...
            'component_roundtrip',component,'parseval',parseval, ...
            'delay_error_s',abs(measured-delay),'midpoint_linear_interpolation_error',interp_error, ...
            'zero_padding_match',zeropad_match,'minimum_gain',min(abs(relative)), ...
            'spectral_condition',max(abs(relative))/min(abs(relative)));
        if isempty(rows),rows=r;else,rows(end+1)=r;end %#ok<AGROW>
        item=struct('meta',meta,'q',q,'tx_m',tx,'f_axis_hz',f,'tau_direct_s',tau, ...
            'direct_native',Hd,'reflect_native',Hr,'total_native',native, ...
            'absolute_DSP',absolute,'relative_DSP',relative,'h_direct_bb',hd, ...
            'h_reflect_bb',hr,'h_total_bb',h,'midpoint_hz',mid,'midpoint_exact',exact, ...
            'transform',"first N bins; right endpoint excluded; ifft(ifftshift); rectangular band; dt=.25ms");
        item.meta.tx_m=tx;
        if isempty(responses),responses=item;else,responses(end+1)=item;end %#ok<AGROW>
    end
end
% Same frequency range: true refinement doubles the periodic window, not resolution.
assert(rows(2).window_s==2*rows(1).window_s && rows(2).dt_s==rows(1).dt_s);
assert(rows(2).midpoint_linear_interpolation_error<rows(1).midpoint_linear_interpolation_error);
rng(data_seed,'twister');bits=randi([0 1],500000,1);
s=((1-2*bits(1:2:end))+1i*(1-2*bits(2:2:end)))/sqrt(2);
rng(noise_seed,'twister');awgn=struct([]);
for eb=[0 4 8]
    gamma=10^(eb/10);noise=(randn(size(s))+1i*randn(size(s)))/sqrt(4*gamma);
    decoded=local_demod(s+noise);errors=sum(decoded~=bits);ber=errors/numel(bits);
    theory=.5*erfc(sqrt(gamma));stderr=sqrt(theory*(1-theory)/numel(bits));
    assert(abs(ber-theory)<=5*stderr+1/numel(bits));
    row=struct('EbN0_dB',eb,'EsN0_dB',eb+10*log10(2),'bits',numel(bits), ...
        'errors',errors,'BER',ber,'theory',theory,'standard_error',stderr, ...
        'expected_complex_noise_power',1/(2*gamma),'observed_complex_noise_power',mean(abs(noise).^2));
    if isempty(awgn),awgn=row;else,awgn(end+1)=row;end %#ok<AGROW>
end
assert(all(local_demod((.4+.7i)*s/(.4+.7i))==bits));
% Analytic flat two-path block controls with a full finite-DFT impulse CP.
flat_errors=zeros(1,numel(responses));cp_errors=flat_errors;
for j=1:numel(responses)
    H=responses(j).relative_DSP(1:end-1);N=numel(H);symbols=s(1:N/2);
    % Single-carrier QPSK impulses at 2 samples/symbol, Fs4k/Rs2k.
    transmitted=zeros(N,1);transmitted(1:2:end)=symbols;h=responses(j).h_total_bb;
    received=ifft(ifftshift(H).*fft(transmitted));
    equalized=ifft(fft(received)./ifftshift(H));equalized=equalized(1:2:end);
    flat_errors(j)=sum(local_demod(equalized)~=local_demod(symbols));
    cp=[transmitted(end-(N-2):end);transmitted];linear=conv(cp,h);
    retained=linear(N:2*N-1);cp_errors(j)=norm(retained-received)/norm(received);
    assert(flat_errors(j)==0 && cp_errors(j)<1e-10);
end
% Three frozen accepted points remain three scalar fixtures, never a dual CIR.
path=fullfile(root,'results','validation','deepsea3d_g1_frequency','frequency_summary.csv');
T=readtable(path,TextType='string');assert(height(T)==3 && all(T.status=="CONDITIONAL_G1_PASS"));
dualmeta=meta;dualmeta.physics_tier="compact_dual_conditional";
dualmeta.scope="alpha2.5 fixed y-offset; only frozen 4/6/8kHz scalar gains";
dualmeta.tx_m=[0 150.5 1500];dualmeta.certified_frequencies_hz=[4000 6000 8000];
dualmeta.parent_surface="fixed dual; amplitudes .05/.025m; K(.4,0)/(.2,.6); phase .2/.7";
dualmeta.support="C2 plateau8/support12";dualmeta.source_evidence=string(path);
dual_native=complex(T.reference_reflect_real,T.reference_reflect_imag);
dual_total=zeros(3,1);dual_dsp=dual_total;dual_errors=dual_total;
for j=1:3
    f=T.frequency_hz(j);dualmeta.frequency_hz=f;deepsea_readiness_guard(dualmeta,"readiness");
    Ld=norm(dualmeta.tx_m-dualmeta.rx_m);tau=Ld/1500;
    dual_total(j)=dual_native(j)+exp(1i*2*pi*f*tau)/(4*pi*Ld);
    dual_dsp(j)=conj(dual_total(j))*exp(1i*2*pi*f*tau);
    equalized=(dual_dsp(j)*s)/dual_dsp(j);
    dual_errors(j)=sum(local_demod(equalized)~=bits);assert(dual_errors(j)==0);
end
guard_rejections=0;
for purpose=["training","statistical_acceptance","PM_BER"]
    try
        deepsea_readiness_guard(meta,purpose);
    catch err
        assert(strcmp(err.identifier,'Readiness:Provenance'));guard_rejections=guard_rejections+1;
    end
end
pm=meta;pm.physics_tier="PM_unvalidated_diagnostic";
try
    deepsea_readiness_guard(pm,"readiness");
catch err
    assert(strcmp(err.identifier,'Readiness:Provenance'));guard_rejections=guard_rejections+1;
end
invaliddual=dualmeta;invaliddual.frequency_hz=5000;
try
    deepsea_readiness_guard(invaliddual,"readiness");
catch err
    assert(strcmp(err.identifier,'Readiness:FrequencyScope'));guard_rejections=guard_rejections+1;
end
assert(guard_rejections==5);
[after,~]=memory;
output=struct('status','DETERMINISTIC_READINESS_PASS','meta',meta,'CIR_checks',rows, ...
    'AWGN',awgn,'known_scalar_bit_errors',0,'flat_two_path_bit_errors',flat_errors, ...
    'flat_CP_discrepancy',cp_errors,'dual_scalar_bit_errors',dual_errors, ...
    'provenance_rejections',guard_rejections,'G1_target_PM','OPEN','WP2_formal',false, ...
    'wall_s',toc(timer),'process_before_MiB',before.MemUsedMATLAB/2^20, ...
    'process_after_MiB',after.MemUsedMATLAB/2^20,'available_before_MiB',sys.PhysicalMemory.Available/2^20, ...
    'peak_memory_MiB',NaN,'peak_kind','unknown; snapshots only', ...
    'receiver_scope','perfect CSI/timing; QPSK symbol AWGN and block CP/frequency equalization; no RRC/FSE/pilot', ...
    'optional_2D_3D','deferred; no new gate');
save(file,'output','responses','bits','dualmeta','dual_native','dual_total','dual_dsp','-v7');
writetable(struct2table(rows),fullfile(outdir,'CIR_checks.csv'));
writetable(struct2table(awgn),fullfile(outdir,'AWGN_checks.csv'));
fid=fopen(fullfile(outdir,'singlelink_readiness.json'),'w');assert(fid>=0);
closefile=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(output,PrettyPrint=true));
fprintf('B %s: AWGN errors %s; flat/dual zero errors; guard rejections%d; wall%.3fs\n', ...
    output.status,mat2str([awgn.errors]),guard_rejections,output.wall_s);
end
function bits=local_demod(symbols)
bits=zeros(2*numel(symbols),1);bits(1:2:end)=real(symbols)<0;bits(2:2:end)=imag(symbols)<0;
end
