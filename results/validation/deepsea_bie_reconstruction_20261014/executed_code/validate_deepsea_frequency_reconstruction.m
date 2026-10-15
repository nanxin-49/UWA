function output=validate_deepsea_frequency_reconstruction()
%VALIDATE_DEEPSEA_FREQUENCY_RECONSTRUCTION Independent off-grid two-path tests.
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
outdir=fullfile(root,'results','validation','deepsea_bie_reconstruction_20261014');
file=fullfile(outdir,'frequency_reconstruction.mat');assert(~isfile(file),'Reconstruction:Exists','No overwrite.');
timer=tic;[mem0,sys]=memory;assert(128<.7*sys.PhysicalMemory.Available/2^20);
base=load(fullfile(root,'results','validation','deepsea_two_track_20261014','singlelink_readiness.mat'),'responses','output');
meta=struct('physics_tier',"analytic_control",'known_path_components',true, ...
    'G1_target_PM',"OPEN",'WP2_formal',false,'c_mps',1500, ...
    'scope',"flat image and known geometric paths; no identification from total H");
rows=struct([]);saved=struct([]);sensitivity=struct([]);
for j=1:numel(base.responses)
    r=base.responses(j);f=r.f_axis_hz;N=numel(f)-1;df=4000/N;
    tx=r.tx_m;rx=r.meta.rx_m;Ld=norm(tx-rx);Lr=norm([tx(1:2),-tx(3)]-rx);
    tau_d=Ld/1500;tau=(Lr-Ld)/1500;scale=1/(4*pi*Ld);
    paths=[conj(r.direct_native),conj(r.reflect_native)].*exp(1i*2*pi*f*tau_d);
    testf=reshape(f(1:end-1)+df*[.25 .5 .75],[],1);testf=sort(testf);
    % Include exact interference minima and very close side points.
    minima=(ceil(4000*tau):floor(8000*tau)).'/tau;
    minima=reshape(minima+[-1e-6 0 1e-6],[],1);testf=unique([testf;minima]);
    exact=scale-exp(-1i*2*pi*testf*tau)/(4*pi*Lr);
    linear=interp1(f,r.relative_DSP,testf,'linear');
    reconstructed=deepsea_delay_separated_reconstruct(f,paths,[0 tau],testf,meta);
    e=norm(reconstructed-exact)/(scale*sqrt(numel(testf)));
    el=norm(linear-exact)/(scale*sqrt(numel(testf)));
    mag=max(abs(abs(reconstructed)-abs(exact)))/scale;
    phase=angle(reconstructed.*conj(exact));phase_weighted=sqrt(sum(abs(exact).^2.*phase.^2)/sum(abs(exact).^2));
    deep=abs(exact)<.05*scale;assert(any(deep));
    deep_error=max(abs(reconstructed(deep)-exact(deep))./abs(exact(deep)));
    phase_max=max(abs(phase));total_relative=max(abs(reconstructed-exact)./abs(exact));
    mid=(f(1:end-1)+f(2:end))/2;
    mid_exact=scale-exp(-1i*2*pi*mid*tau)/(4*pi*Lr);
    mid_linear=interp1(f,r.relative_DSP,mid,'linear');
    old=base.output.CIR_checks(j).midpoint_linear_interpolation_error;
    baseline_match=abs(norm(mid_linear-mid_exact)/(scale*sqrt(N))-old);
    assert(baseline_match<1e-12);
    % Refinement to all half bins produces twice the physical transform window.
    f_dense=linspace(4000,8000,2*N+1).';
    dense=deepsea_delay_separated_reconstruct(f,paths,[0 tau],f_dense,meta);
    exact_dense=scale-exp(-1i*2*pi*f_dense*tau)/(4*pi*Lr);
    h=ifft(ifftshift(dense(1:end-1)));h_exact=ifft(ifftshift(exact_dense(1:end-1)));
    h_linear=ifft(ifftshift(interp1(f,r.relative_DSP,f_dense(1:end-1),'linear')));
    cir_error=norm(h-h_exact)/norm(h_exact);cir_linear_error=norm(h_linear-h_exact)/norm(h_exact);
    dt=1/4000;t=(0:2*N-1).'*dt;
    % Rectangular band tails, saved rather than cropped or called extra paths.
    tail=(t>.020 & t<2/df-.005);
    tail_exact=sum(abs(h_exact(tail)).^2)/sum(abs(h_exact).^2);
    tail_error=abs(sum(abs(h(tail)).^2)/sum(abs(h).^2)-tail_exact);
    assert(max([e,mag,phase_weighted,deep_error,total_relative,cir_error])<1e-3);
    for delay_error=[0 1e-6 1e-4 1e-3]
        shifted=deepsea_delay_separated_reconstruct(f,paths,[0 tau+delay_error],testf,meta);
        z=struct('q',r.q,'source_points',N+1,'delay_error_s',delay_error, ...
            'complex_direct_scale',norm(shifted-exact)/(scale*sqrt(numel(testf))));
        if isempty(sensitivity),sensitivity=z;else,sensitivity(end+1)=z;end %#ok<AGROW>
    end
    row=struct('q',r.q,'source_points',N+1,'df_hz',df,'query_points',numel(testf), ...
        'linear_complex',el,'separated_complex',e,'magnitude_max_direct_scale',mag, ...
        'phase_weighted_rad',phase_weighted,'phase_max_rad',phase_max, ...
        'total_relative_max',total_relative,'deep_fade_relative_max',deep_error, ...
        'minimum_exact_gain',min(abs(exact)),'CIR_error',cir_error,'linear_CIR_error',cir_linear_error, ...
        'tail_energy_exact',tail_exact,'tail_energy_error',tail_error,'baseline_match',baseline_match, ...
        'window_s',2/df,'dt_s',dt,'passed',true);
    item=struct('q',r.q,'source_f_hz',f,'path_response',paths,'delay_s',[0 tau], ...
        'query_hz',testf,'exact',exact,'linear',linear,'reconstructed',reconstructed, ...
        'f_dense_hz',f_dense,'dense_response',dense,'h',h,'h_exact',h_exact,'h_linear',h_linear);
    if isempty(rows),rows=row;saved=item;else,rows(end+1)=row;saved(end+1)=item;end %#ok<AGROW>
end
rejections=0;
for tier=["compact_dual_conditional","PM_unvalidated_diagnostic"]
    rejected=meta;rejected.physics_tier=tier;
    try
        deepsea_delay_separated_reconstruct(f,paths,[0 tau],f,rejected);
    catch err
        assert(strcmp(err.identifier,'Reconstruction:Scope'));rejections=rejections+1;
    end
end
assert(rejections==2);[mem1,~]=memory;
output=struct('status','PASS_ANALYTIC_KNOWN_PATHS_ONLY','rows',rows,'sensitivity',sensitivity, ...
    'meta',meta,'scope_rejections',rejections,'generic_total_only_reconstruction',false, ...
    'dual_intermediate_certified',false,'noisy_flat_two_path_next_experiment_suitable',true, ...
    'G1_target_PM','OPEN','WP2_formal',false,'wall_s',toc(timer), ...
    'memory_before_MiB',mem0.MemUsedMATLAB/2^20,'memory_after_MiB',mem1.MemUsedMATLAB/2^20);
save(file,'output','saved','-v7');writetable(struct2table(rows),fullfile(outdir,'reconstruction.csv'));
writetable(struct2table(sensitivity),fullfile(outdir,'delay_sensitivity.csv'));
fid=fopen(fullfile(outdir,'frequency_reconstruction.json'),'w');assert(fid>=0);closer=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(output,PrettyPrint=true));
fprintf('Reconstruction PASS, max complex %.3g, max fade-relative %.3g, wall %.3fs\n',max([rows.separated_complex]),max([rows.deep_fade_relative_max]),output.wall_s);
end
