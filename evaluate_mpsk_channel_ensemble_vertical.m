function result=evaluate_mpsk_channel_ensemble_vertical(channel_set,cfg)
%EVALUATE_MPSK_CHANNEL_ENSEMBLE_VERTICAL Common MPSK Monte Carlo evaluator.
%   CHANNEL_SET.H_fm is F-by-L with CHANNEL_SET.f_axis_hz. Alternatively,
%   CHANNEL_SET.h_tm is Ntap-by-L. Noise is generated only after channel
%   convolution. Perfect channel knowledge is used for peak sync and MMSE.

arguments
    channel_set (1,1) struct
    cfg (1,1) struct = struct()
end
c=local_defaults(cfg); k=round(log2(c.M));
if isfield(channel_set,'H_fm'), L=size(channel_set.H_fm,2);
elseif isfield(channel_set,'h_tm'), L=size(channel_set.h_tm,2);
else, error('channel_set requires H_fm or h_tm.'); end
rng(c.bits_seed,'twister'); bits_tx=randi([0,1],c.n_sym*k,1);
[tx,bits_tx]=modem_psk('modulate',bits_tx,c.M); N=numel(tx); E=numel(c.EbN0_dB_list);
bit_err_pre=zeros(L,E); bit_err_eq=zeros(L,E); sym_err_pre=zeros(L,E); sym_err_eq=zeros(L,E);
noiseless_bit_err=zeros(L,1); noiseless_sym_err=zeros(L,1); effective_snr=zeros(L,E);
tap_count=zeros(L,1); peak_index=zeros(L,1); prepeak=zeros(L,1); closure=zeros(L,1); energy_kept=zeros(L,1);
timer=tic;
for mm=1:L
    if isfield(channel_set,'H_fm')
        input=struct('H_f',channel_set.H_fm(:,mm),'f_axis_hz',channel_set.f_axis_hz);
        if isfield(channel_set,'phase_reference_meta')
            input.phase_reference_meta=channel_set.phase_reference_meta;
        end
    else
        input=struct('h_t',channel_set.h_tm(:,mm));
    end
    [h,t]=build_communication_taps_vertical(input,struct('symbol_rate_hz',c.symbol_rate_hz, ...
        'n_fft',c.tap_fft_size,'f_ref_hz',c.f_ref_hz,'tap_energy_ratio',c.tap_energy_ratio, ...
        'circular_peak_align',c.circular_peak_align));
    if c.circular_peak_align
        [rx_clean,h_eq]=local_apply_causal_head(tx,h);
    else
        [rx_clean,h_eq]=local_apply_peak_sync(tx,h);
    end
    [be,se]=local_error_counts(bits_tx,modem_psk('demodulate',local_mmse_equalize(rx_clean,h_eq,c.mmse_reg_eps),c.M),c.M);
    noiseless_bit_err(mm)=be; noiseless_sym_err(mm)=se;
    tap_count(mm)=t.tap_count; peak_index(mm)=t.peak_index; prepeak(mm)=t.discarded_pre_peak_energy_fraction_if_peak_sync;
    closure(mm)=t.frequency_ifft_closure_relative_error; energy_kept(mm)=t.tap_energy_kept;
    for ee=1:E
        nc=struct('enable_noise',true,'model',c.noise_model,'ebn0_db',c.EbN0_dB_list(ee), ...
            'bits_per_symbol',k,'seed',c.noise_seed_base+100000*mm+ee,'custom_noise_fn',[]);
        if strcmpi(c.ebn0_reference,'rx_clean'), ref=rx_clean; else, ref=tx; end
        [rx,~,nm]=noise_inject_vertical(rx_clean,nc,ref); effective_snr(mm,ee)=nm.effective_snr_db;
        [bit_err_pre(mm,ee),sym_err_pre(mm,ee)]=local_error_counts(bits_tx,modem_psk('demodulate',rx,c.M),c.M);
        eq=local_mmse_equalize(rx,h_eq,c.mmse_reg_eps);
        [bit_err_eq(mm,ee),sym_err_eq(mm,ee)]=local_error_counts(bits_tx,modem_psk('demodulate',eq,c.M),c.M);
    end
end
per_channel_ber=bit_err_eq/(N*k); per_channel_ser=sym_err_eq/N;
per_channel_ber_pre=bit_err_pre/(N*k); per_channel_ser_pre=sym_err_pre/N;
ber=sum(bit_err_eq,1)/(L*N*k); ser=sum(sym_err_eq,1)/(L*N);
ber_pre=sum(bit_err_pre,1)/(L*N*k); ser_pre=sum(sym_err_pre,1)/(L*N);
result=struct('name',local_field(channel_set,'name','external'), ...
    'wind_speed_mps',local_field(channel_set,'wind_speed_mps',NaN), ...
    'channel_count',L,'cfg',c,'BER',ber,'SER',ser,'BER_pre_equalizer',ber_pre, ...
    'SER_pre_equalizer',ser_pre,'per_channel_BER',per_channel_ber, ...
    'per_channel_SER',per_channel_ser,'per_channel_BER_pre_equalizer',per_channel_ber_pre, ...
    'per_channel_SER_pre_equalizer',per_channel_ser_pre, ...
    'BER_wilson95',local_wilson(sum(bit_err_eq,1),L*N*k), ...
    'SER_wilson95',local_wilson(sum(sym_err_eq,1),L*N), ...
    'BER_cluster95',local_cluster_ci(per_channel_ber,c.bootstrap_count,c.bootstrap_seed+1), ...
    'SER_cluster95',local_cluster_ci(per_channel_ser,c.bootstrap_count,c.bootstrap_seed+2), ...
    'effective_snr_db_mean',mean(effective_snr,1), ...
    'noiseless_BER',sum(noiseless_bit_err)/(L*N*k), ...
    'noiseless_SER',sum(noiseless_sym_err)/(L*N), ...
    'tap_meta',struct('tap_count',tap_count,'peak_index',peak_index, ...
        'prepeak_energy_fraction',prepeak, ...
        'actual_discard_fraction',prepeak.*double(~c.circular_peak_align), ...
        'ifft_closure_error',closure,'energy_kept',energy_kept), ...
    'elapsed_s',toc(timer),'noise_stored_in_channel',false, ...
    'equalizer','perfect-known-channel frequency-domain MMSE after peak_sync', ...
    'synchronization',local_sync_name(c.circular_peak_align));
end

function c=local_defaults(c)
d=struct('M',4,'n_sym',8000,'EbN0_dB_list',0:2:20,'symbol_rate_hz',1000, ...
    'tap_fft_size',2048,'f_ref_hz',6000,'tap_energy_ratio',0.999, ...
    'mmse_reg_eps',1e-6,'bits_seed',910001,'noise_seed_base',920000, ...
    'noise_model','awgn','ebn0_reference','rx_clean','bootstrap_count',2000,'bootstrap_seed',930001);
d.circular_peak_align=true;
n=fieldnames(d); for ii=1:numel(n), if ~isfield(c,n{ii}), c.(n{ii})=d.(n{ii}); end, end
if abs(round(log2(c.M))-log2(c.M))>1e-12, error('M must be a power of two.'); end
c.EbN0_dB_list=c.EbN0_dB_list(:).';
end
function name=local_sync_name(tf)
if tf, name='shortest circular 99.9% energy window, unwrapped to causal head';
else, name='legacy unshifted IFFT prefix, then peak_sync'; end
end

function [rx,h_eq]=local_apply_peak_sync(tx,h)
e=abs(h).^2; [~,p]=max(e); h_eq=h(p:end); full=conv(tx,h_eq,'full'); rx=full(1:numel(tx));
end
function [rx,h_eq]=local_apply_causal_head(tx,h)
h_eq=h(:); full=conv(tx,h_eq,'full'); rx=full(1:numel(tx));
end
function y=local_mmse_equalize(x,h,reg)
N=numel(x); L=min(numel(h),N); hp=[h(1:L);zeros(N-L,1)]; H=fft(hp); y=ifft(fft(x).*conj(H)./(abs(H).^2+reg));
end
function [be,se]=local_error_counts(a,b,M)
n=min(numel(a),numel(b)); a=a(1:n); b=b(1:n); be=sum(a~=b); k=round(log2(M)); ns=floor(n/k);
aa=reshape(a(1:ns*k),k,[]); bb=reshape(b(1:ns*k),k,[]); se=sum(any(aa~=bb,1));
end
function ci=local_wilson(err,n)
p=err/n; z=1.95996398454005; den=1+z^2/n; ctr=(p+z^2/(2*n))/den;
half=z*sqrt(p.*(1-p)/n+z^2/(4*n^2))/den; ci=[max(ctr-half,0);min(ctr+half,1)];
end
function ci=local_cluster_ci(x,B,seed)
rng(seed,'twister'); [L,E]=size(x); means=zeros(B,E);
for bb=1:B, idx=randi(L,L,1); means(bb,:)=mean(x(idx,:),1); end
means=sort(means,1); lo=max(1,round(.025*(B-1)+1)); hi=min(B,round(.975*(B-1)+1)); ci=[means(lo,:);means(hi,:)];
end
function v=local_field(s,n,d), if isfield(s,n), v=s.(n); else, v=d; end, end
