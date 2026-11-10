function audit = audit_pe_screen_p22_qpsk(run_id)
%AUDIT_PE_SCREEN_P22_QPSK Reopen and replay saved known-CSI simulation.
arguments
    run_id (1,:) char
end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
outdir=fullfile(root,'results','communication','pe_screen_p22_qpsk',run_id);
s=load(fullfile(outdir,'qpsk_results.mat'));
t=readtable(fullfile(outdir,'performance.csv'));
cfg=jsondecode(fileread(fullfile(outdir,'config.json')));
gate=jsondecode(fileread(fullfile(root,cfg.cir_gate)));
assert(gate.common_window_gate_pass);
q=load(fullfile(root,cfg.channel_source),'h257');
assert(height(t)==2*numel(cfg.snr_db));
assert(isequal(s.channels,circshift(q.h257,cfg.common_causal_guard_samples,1)));
assert(isequal(s.cfg,cfg));
assert(numel(s.bits)==cfg.bit_count_per_snr);
x=reshape(s.tx,cfg.block_symbols,cfg.frames_per_snr);
cp=cfg.cyclic_prefix_symbols;
tx_cp=[x(end-cp+1:end,:);x];
max_eq_replay_error=0;
max_csv_difference=0;
for k=1:numel(cfg.snr_db)
    for j=1:2
        row=(k-1)*2+j;
        assert(t.snr_db(row)==cfg.snr_db(k));
        if j==1, assert(strcmp(string(t.surface(row)),"flat"));
        else, assert(strcmp(string(t.surface(row)),"weak_sine")); end
        h=s.channels(:,j);
        full=filter(h,1,tx_cp,[],1);
        yc=full(cp+1:cp+cfg.block_symbols,:);
        H=fft(h,cfg.block_symbols);
        nc=struct('enable_noise',true,'model','awgn', ...
            'snr_db',cfg.snr_db(k),'seed',cfg.noise_seed_base+100*k);
        [yn,~,nm]=noise_inject_vertical(yc,nc,yc);
        weights=conj(H)./(abs(H).^2+nm.noise_power);
        eq=ifft(fft(yn,[],1).*weights,[],1);
        max_eq_replay_error=max(max_eq_replay_error, ...
            max(abs(eq-s.eq_symbols(:,:,j,k)),[],'all'));
        b=local_natural_to_gray(modem_psk('demodulate',eq(:),4));
        assert(isequal(logical(b),s.bits_hat(:,j,k)));
        errors=b~=s.bits;
        bit_errors=sum(errors);
        frame_errors=sum(any(reshape(errors,cfg.frame_bits,[]),1));
        ber=bit_errors/numel(errors);
        per=frame_errors/cfg.frames_per_snr;
        evm=sqrt(mean(abs(eq(:)-s.tx).^2));
        assert(t.bit_errors(row)==bit_errors && t.frame_errors(row)==frame_errors);
        max_csv_difference=max(max_csv_difference,max(abs([ ...
            t.ber(row)-ber,t.per(row)-per,t.evm_rms(row)-evm, ...
            t.effective_snr_db(row)-nm.effective_snr_db])));
        assert(t.ber_low95(row)<=ber && t.ber_high95(row)>=ber);
        assert(t.per_low95(row)<=per && t.per_high95(row)>=per);
    end
end
assert(max_eq_replay_error<=1e-12 && max_csv_difference<=1e-12);
assert(strcmp(s.receipt.status,'PE_SCREEN_COMM_SIMULATION_PASS'));
assert(all(s.receipt.noiseless_bit_errors==0));
audit=struct('status','PE_SCREEN_QPSK_INDEPENDENT_AUDIT_PASS', ...
    'max_equalized_replay_abs_error',max_eq_replay_error, ...
    'max_csv_metric_difference',max_csv_difference, ...
    'bit_count_per_snr',cfg.bit_count_per_snr, ...
    'snr_point_count',numel(cfg.snr_db), ...
    'surface_count',2,'noiseless_bit_errors',s.receipt.noiseless_bit_errors, ...
    'cir_gate_status',gate.status);
figdir=fullfile(outdir,'figures');
assert(~isfolder(figdir)); mkdir(figdir);
f=figure('Visible','off','Color','w');
tiledlayout(3,1);
names={'Flat','Weak sine'};
values={s.ber,s.per,s.evm};
labels={'BER','PER','EVM RMS'};
for k=1:3
    nexttile;
    for j=1:2
        y=values{k}(j,:);
        if k<3, semilogy(cfg.snr_db,max(y,1e-7),'-o','LineWidth',1.2); 
        else, plot(cfg.snr_db,y,'-o','LineWidth',1.2); end
        hold on;
    end
    grid on; xlabel('Post-channel SNR (dB)'); ylabel(labels{k});
    legend(names,'Location','best');
end
savefig(f,fullfile(figdir,'ber_per_evm.fig'));
exportgraphics(f,fullfile(figdir,'ber_per_evm.png'),'Resolution',160);
close(f);
fid=fopen(fullfile(outdir,'audit.json'),'w'); assert(fid>0);
cleaner=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s',jsonencode(audit,PrettyPrint=true));
disp(audit);
end

function gray=local_natural_to_gray(natural)
gray=natural;
gray(2:2:end)=xor(natural(1:2:end),natural(2:2:end));
end
