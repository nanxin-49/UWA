function receipt = diagnose_pe_screen_p22_aliasing(run_id,parent_run_id)
%DIAGNOSE_PE_SCREEN_P22_ALIASING Compute independent 257-grid odd bins.
arguments
    run_id (1,:) char
    parent_run_id (1,:) char
end
assert(~isempty(regexp(run_id,'^[A-Za-z0-9_]+$','once')));
assert(~isempty(regexp(parent_run_id,'^[A-Za-z0-9_]+$','once')));
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
setup_vertical_project;
parent_dir=fullfile(root,'results','validation','pe_screen_p22_comm',parent_run_id);
outdir=fullfile(root,'results','validation','pe_screen_p22_aliasing',run_id);
assert(~isfolder(outdir) && ~isfile(outdir),'Run ID already exists');
p=load(fullfile(parent_dir,'pe_screen_p22_broadband.mat'));
assert(numel(p.f129)==129 && all(p.f129==4000+(0:128)'*31.25));
source_path=fullfile(root,'results','comparisons', ...
    'two_dimensional_channel_closeout','20261026_s3a_seed260003_01', ...
    'input','frozen_source.mat');
v=load(source_path,'src'); frozen=v.src;
assert(frozen.n_original==984 && frozen.z_tx_m==100);
N=984; W=192.1875; c0=1500; ztx=100; zrx=3;
x=-W/2+(0:N-1)'*W/N;
source_tx=ifft(frozen.coeff(:));
assert(max(abs(source_tx-exp(-.5*(x/.3).^2)))<1e-12);
[eta_weak,~]=closeout_s2b_surface(x,.01);
f257=4000+(0:256)'*15.625;
assert(isequal(f257(1:2:end),p.f129));
H_direct=complex(zeros(257,2)); H_reflect=H_direct;
H_direct(1:2:end,:)=p.H_direct;
H_reflect(1:2:end,:)=p.H_reflect;
seconds=zeros(257,2); direct_error=zeros(257,2);
tau0=(ztx-zrx)/c0;
tick_all=tic;
for isurf=1:2
    eta=zeros(N,1);
    if isurf==2, eta=eta_weak; end
    for ii=2:2:256
        f=f257(ii);
        src=closeout2d_source_at_frequency(frozen,f,c0);
        cfg=struct('frequency_hz',f,'c0_mps',c0,'xw_m',W, ...
            'nx',N,'z_tx_m',ztx,'z_rx_m',zrx,'sigma_src_m',.3, ...
            'surface_reflect_coeff',-1,'step_m',.05,'x_rx_m',0, ...
            'reflection_model','model0_normal', ...
            'source_field_override',source_tx,'eta_z_m',eta);
        tick=tic; pe=run_pe_physical_surface_model0_validation(cfg);
        seconds(ii,isurf)=toc(tick);
        phase=exp(1i*2*pi*f*tau0);
        ud=pe.direct_receiver*phase;
        ua=fk_single_receiver_incident(src,0,zrx);
        direct_error(ii,isurf)=abs(ud-ua)/abs(ua);
        H_direct(ii,isurf)=conj(ud)*phase;
        H_reflect(ii,isurf)=conj(pe.reflected_native(pe.ix_rx))*phase;
    end
end
H_total=H_direct+H_reflect;
h257=ifft(H_total(1:256,:),[],1);
hd257=ifft(H_direct(1:256,:),[],1);
hr257=ifft(H_reflect(1:256,:),[],1);
alias65_129=zeros(1,2); alias129_257=zeros(1,2);
tail16_32=zeros(1,2); tail32_64=zeros(1,2);
direct_tail16_32=zeros(1,2); reflect_tail16_32=zeros(1,2);
fft_error=zeros(1,2); max_odd_direct_error=max(direct_error,[],'all');
for j=1:2
    a=p.h65(:,j)-p.h129(1:64,j)-p.h129(65:128,j);
    alias65_129(j)=norm(a)/norm(p.h129(:,j));
    b=p.h129(:,j)-h257(1:128,j)-h257(129:256,j);
    alias129_257(j)=norm(b)/norm(h257(:,j));
    den=sum(abs(h257(:,j)).^2);
    tail16_32(j)=sum(abs(h257(65:128,j)).^2)/den;
    tail32_64(j)=sum(abs(h257(129:end,j)).^2)/den;
    direct_tail16_32(j)=sum(abs(hd257(65:128,j)).^2)/sum(abs(hd257(:,j)).^2);
    reflect_tail16_32(j)=sum(abs(hr257(65:128,j)).^2)/sum(abs(hr257(:,j)).^2);
    fft_error(j)=norm(fft(h257(:,j))-H_total(1:256,j))/norm(H_total(1:256,j));
end
% This is the original unwindowed 65/129 gate, deliberately unchanged.
original_refinement=max(p.receipt.cir_refinement_relative);
gates=struct('new_direct_source',max_odd_direct_error<=1e-10, ...
    'nested_grid',isequal(H_total(1:2:end,:),p.H_total), ...
    'alias_identity_65_129',max(alias65_129)<=1e-12, ...
    'alias_identity_129_257',max(alias129_257)<=1e-12, ...
    'fft_roundtrip',max(fft_error)<=1e-10, ...
    'original_frequency_refinement',original_refinement<=1e-3);
failed=sum(~cell2mat(struct2cell(gates)));
status='ALIASING_DIAGNOSIS_PASS_CIR_GATE_FAILED';
if failed==0, status='PE_SCREEN_CIR_GATE_PASS'; end
receipt=struct('status',status,'matlab_version',version, ...
    'parent_run_id',parent_run_id,'source_path',source_path, ...
    'frequency_count',257,'new_pe_calls',256, ...
    'frequency_spacing_hz',15.625,'time_window_s',.064, ...
    'time_origin_s',tau0,'original_refinement_relative',original_refinement, ...
    'alias_identity_65_129',alias65_129,'alias_identity_129_257',alias129_257, ...
    'tail_energy_16_32ms',tail16_32,'tail_energy_32_64ms',tail32_64, ...
    'direct_tail_energy_16_32ms',direct_tail16_32, ...
    'reflect_tail_energy_16_32ms',reflect_tail16_32, ...
    'max_new_direct_source_relative',max_odd_direct_error, ...
    'fft_roundtrip_relative',fft_error,'pe_seconds',sum(seconds,'all'), ...
    'elapsed_seconds',toc(tick_all),'failed_check_count',failed,'gates',gates, ...
    'bie_calls',0,'fullfk_calls',0,'three_dimensional_calls',0);
mkdir(outdir);
save(fullfile(outdir,'pe_screen_p22_257.mat'),'receipt','f257', ...
    'H_direct','H_reflect','H_total','h257','hd257','hr257','seconds', ...
    'direct_error','-v7.3');
t=table(repmat(f257,2,1),repelem(["flat";"weak_sine"],257), ...
    real(H_direct(:)),imag(H_direct(:)),real(H_reflect(:)),imag(H_reflect(:)), ...
    real(H_total(:)),imag(H_total(:)),seconds(:), ...
    'VariableNames',{'frequency_hz','surface','H_direct_real','H_direct_imag', ...
    'H_reflect_real','H_reflect_imag','H_total_real','H_total_imag','pe_seconds'});
writetable(t,fullfile(outdir,'frequency_response_257.csv'));
delay_s=(0:255)'*.00025;
td=table(repmat(delay_s,2,1),repelem(["flat";"weak_sine"],256), ...
    real(h257(:)),imag(h257(:)),abs(h257(:)).^2, ...
    'VariableNames',{'delay_s','surface','cir_real','cir_imag','pdp'});
writetable(td,fullfile(outdir,'cir_pdp_257.csv'));
fid=fopen(fullfile(outdir,'receipt.json'),'w'); assert(fid>0);
cleaner=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s',jsonencode(receipt,PrettyPrint=true));
disp(receipt);
end
