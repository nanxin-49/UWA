function receipt = run_pe_screen_p22_qpsk(run_id)
%RUN_PE_SCREEN_P22_QPSK Known-CSI QPSK through saved P22 PE CIR.
arguments
    run_id (1,:) char
end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
setup_vertical_project;
outdir=fullfile(root,'results','communication','pe_screen_p22_qpsk',run_id);
config_path=fullfile(outdir,'config.json');
assert(isfile(config_path) && ~isfile(fullfile(outdir,'qpsk_results.mat')));
cfg=jsondecode(fileread(config_path));
assert(strcmp(cfg.status,'FROZEN_BEFORE_QPSK_RUN'));
assert(cfg.M==4 && cfg.symbol_rate_hz==4000 && ...
    cfg.channel_tap_count==256 && cfg.block_symbols==512 && ...
    cfg.cyclic_prefix_symbols==255 && cfg.frame_bits==1024 && ...
    cfg.frames_per_snr==256 && cfg.bit_count_per_snr==262144);
gate=jsondecode(fileread(fullfile(root,cfg.cir_gate)));
assert(gate.common_window_gate_pass && ...
    strcmp(gate.status,'PE_SCREEN_COMMON_WINDOW_CIR_INDEPENDENT_AUDIT_PASS'));
q=load(fullfile(root,cfg.channel_source),'h257','H_total','f257');
assert(isequal(size(q.h257),[256,2]) && numel(q.f257)==257);
assert(max(abs(fft(q.h257,[],1)-q.H_total(1:256,:)),[],'all')<1e-12);
guard=cfg.common_causal_guard_samples;
channels=circshift(q.h257,guard,1);
nblock=cfg.block_symbols; nframe=cfg.frames_per_snr;
assert(cfg.cyclic_prefix_symbols>=size(channels,1)-1);
rng(cfg.bits_seed,'twister');
bits=randi([0,1],cfg.bit_count_per_snr,1);
bits_natural=local_gray_to_natural(bits);
[tx,used]=modem_psk('modulate',bits_natural,cfg.M);
assert(isequal(used,bits_natural));
x=reshape(tx,nblock,nframe);
cp=cfg.cyclic_prefix_symbols;
tx_cp=[x(end-cp+1:end,:);x];
nsnr=numel(cfg.snr_db);
eq_symbols=complex(zeros(nblock,nframe,2,nsnr));
bits_hat=false(numel(bits),2,nsnr);
noiseless_error=zeros(1,2); noiseless_bit_errors=zeros(1,2);
convolution_error=zeros(1,2); circular_error=zeros(1,2);
ber=zeros(2,nsnr); per=ber; evm=ber;
bit_errors=zeros(2,nsnr); frame_errors=bit_errors;
ber_low=ber;ber_high=ber;per_low=ber;per_high=ber;
effective_snr=ber; noise_power=ber;
tick=tic;
for j=1:2
    h=channels(:,j);
    Hblock=fft(h,nblock);
    yfull=filter(h,1,tx_cp,[],1);
    yc=yfull(cp+1:cp+nblock,:);
    check=conv(tx_cp(:,1),h);
    convolution_error(j)=norm(yfull(:,1)-check(1:size(yfull,1)))/norm(check(1:size(yfull,1)));
    circular_error(j)=norm(fft(yc,[],1)-fft(x,[],1).*Hblock)/norm(fft(yc,[],1));
    x0=ifft(fft(yc,[],1)./Hblock,[],1);
    noiseless_error(j)=norm(x0-x)/norm(x);
    nohat=local_natural_to_gray(modem_psk('demodulate',x0(:),4));
    noiseless_bit_errors(j)=sum(nohat~=bits);
    for k=1:nsnr
        noise_cfg=struct('enable_noise',true,'model','awgn', ...
            'snr_db',cfg.snr_db(k),'seed',cfg.noise_seed_base+100*k);
        [yn,~,nm]=noise_inject_vertical(yc,noise_cfg,yc);
        W=conj(Hblock)./(abs(Hblock).^2+nm.noise_power);
        xe=ifft(fft(yn,[],1).*W,[],1);
        eq_symbols(:,:,j,k)=xe;
        bh=local_natural_to_gray(modem_psk('demodulate',xe(:),4));
        bits_hat(:,j,k)=logical(bh);
        errors=bh~=bits;
        bit_errors(j,k)=sum(errors);
        frame_errors(j,k)=sum(any(reshape(errors,cfg.frame_bits,nframe),1));
        ber(j,k)=bit_errors(j,k)/numel(bits);
        per(j,k)=frame_errors(j,k)/nframe;
        evm(j,k)=sqrt(mean(abs(xe(:)-tx).^2));
        [ber_low(j,k),ber_high(j,k)]=local_interval(bit_errors(j,k),numel(bits));
        [per_low(j,k),per_high(j,k)]=local_interval(frame_errors(j,k),nframe);
        effective_snr(j,k)=nm.effective_snr_db;
        noise_power(j,k)=nm.noise_power;
    end
end
gates=struct('cir_qualification',gate.common_window_gate_pass, ...
    'linear_convolution',max(convolution_error)<=1e-12, ...
    'cyclic_prefix',max(circular_error)<=1e-12, ...
    'noiseless_equalization',max(noiseless_error)<=1e-10, ...
    'noiseless_zero_bit_errors',all(noiseless_bit_errors==0), ...
    'sufficient_bits',numel(bits)>=100000);
failed=sum(~cell2mat(struct2cell(gates)));
status='PE_SCREEN_QPSK_SIMULATION_FAILED';
if failed==0, status='PE_SCREEN_COMM_SIMULATION_PASS'; end
receipt=struct('status',status,'matlab_version',version, ...
    'scope','P22 PE phase-screen simulation; not formal S5/S6 physics acceptance', ...
    'channel_cir_gate',gate.status,'symbol_rate_hz',cfg.symbol_rate_hz, ...
    'block_symbols',nblock,'cyclic_prefix_symbols',cp, ...
    'guard_samples',guard,'bit_count_per_snr',numel(bits), ...
    'snr_db',cfg.snr_db,'noiseless_equalization_relative',noiseless_error, ...
    'noiseless_bit_errors',noiseless_bit_errors, ...
    'linear_convolution_relative',convolution_error, ...
    'cyclic_prefix_relative',circular_error,'bit_errors',bit_errors, ...
    'frame_errors',frame_errors,'ber',ber,'per',per,'evm_rms',evm, ...
    'effective_snr_db',effective_snr,'noise_power',noise_power, ...
    'elapsed_seconds',toc(tick),'failed_check_count',failed,'gates',gates, ...
    'new_acoustic_solver_calls',0,'bie_calls',0,'fullfk_calls',0, ...
    'three_dimensional_calls',0);
save(fullfile(outdir,'qpsk_results.mat'),'cfg','receipt','bits','tx', ...
    'channels','eq_symbols','bits_hat','ber','per','evm', ...
    'bit_errors','frame_errors','ber_low','ber_high','per_low','per_high', ...
    'effective_snr','noise_power','-v7.3');
surface=repmat(["flat";"weak_sine"],nsnr,1);
snr_db=repelem(cfg.snr_db(:),2);
summary=table(surface,snr_db,bit_errors(:),ber(:),ber_low(:),ber_high(:), ...
    frame_errors(:),per(:),per_low(:),per_high(:),evm(:),effective_snr(:), ...
    'VariableNames',{'surface','snr_db','bit_errors','ber','ber_low95', ...
    'ber_high95','frame_errors','per','per_low95','per_high95', ...
    'evm_rms','effective_snr_db'});
writetable(summary,fullfile(outdir,'performance.csv'));
fid=fopen(fullfile(outdir,'receipt.json'),'w'); assert(fid>0);
cleaner=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s',jsonencode(receipt,PrettyPrint=true));
disp(receipt);
end

function natural=local_gray_to_natural(gray)
natural=gray;
natural(2:2:end)=xor(gray(1:2:end),gray(2:2:end));
end

function gray=local_natural_to_gray(natural)
gray=natural;
gray(2:2:end)=xor(natural(1:2:end),natural(2:2:end));
end

function [lo,hi]=local_interval(k,n)
if k==0
    lo=0; hi=1-0.05^(1/n); return
end
p=k/n; z=1.95996398454005;
center=(p+z^2/(2*n))/(1+z^2/n);
half=z*sqrt(p*(1-p)/n+z^2/(4*n^2))/(1+z^2/n);
lo=max(0,center-half); hi=min(1,center+half);
end
