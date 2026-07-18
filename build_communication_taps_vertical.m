function [h_bb,meta]=build_communication_taps_vertical(channel_input,options)
%BUILD_COMMUNICATION_TAPS_VERTICAL Convert external H(f) or h(t) to symbol taps.
%   The H(f) path reproduces the current communication demo convention:
%   shift around f_ref, interpolate on the symbol-rate FFT grid, IFFT, then
%   keep the causal prefix containing the requested energy. Zero padding is
%   interpolation only and is not reported as physical delay resolution.

arguments
    channel_input (1,1) struct
    options (1,1) struct = struct()
end
o=local_defaults(options);
if isfield(channel_input,'H_f') && ~isempty(channel_input.H_f)
    if ~isfield(channel_input,'f_axis_hz'), error('H_f input requires f_axis_hz.'); end
    f=double(channel_input.f_axis_hz(:)); H=double(channel_input.H_f(:));
    if numel(f)~=numel(H) || numel(f)<2, error('f_axis_hz and H_f must have equal length >=2.'); end
    df=diff(f); if any(df<=0) || max(abs(df-mean(df)))>1e-9*max(abs(mean(df)),1)
        error('External H(f) frequency axis must be strictly increasing and equally spaced.');
    end
    [~,iref]=min(abs(f-o.f_ref_hz)); fc=f(iref); f_rel=f-fc;
    fbb=((0:o.n_fft-1).'-floor(o.n_fft/2))*(o.symbol_rate_hz/o.n_fft);
    if min(fbb)<min(f_rel)-1e-9 || max(fbb)>max(f_rel)+1e-9
        error('H(f) does not cover the requested symbol-rate baseband.');
    end
    Hbb=interp1(f_rel,H,fbb,'linear');
    h_full=ifft(ifftshift(Hbb));
    closure=max(abs(fftshift(fft(h_full))-Hbb))/max(max(abs(Hbb)),eps);
    source='frequency_response';
    physical_resolution_s=1/(f(end)-f(1));
    maximum_unambiguous_delay_s=1/mean(df);
elseif isfield(channel_input,'h_t') && ~isempty(channel_input.h_t)
    h_full=double(channel_input.h_t(:)); Hbb=fftshift(fft(h_full));
    fbb=((0:numel(h_full)-1).'-floor(numel(h_full)/2))* ...
        (o.symbol_rate_hz/numel(h_full));
    closure=0; fc=o.f_ref_hz; iref=NaN; source='time_response';
    physical_resolution_s=1/o.symbol_rate_hz;
    maximum_unambiguous_delay_s=numel(h_full)/o.symbol_rate_hz;
else
    error('channel_input must contain nonempty H_f or h_t.');
end
[~,raw_peak]=max(abs(h_full).^2); if isempty(raw_peak), raw_peak=1; end
raw_energy=abs(h_full).^2; total=sum(raw_energy);
if o.circular_peak_align && total>0
    [window_start,window_length]=local_shortest_circular_window(raw_energy,o.tap_energy_ratio);
    idx=mod((window_start-1)+(0:window_length-1),numel(h_full))+1;
    h_ordered=h_full(idx); circular_shift=1-window_start;
else
    window_start=1; window_length=numel(h_full); circular_shift=0; h_ordered=h_full;
end
energy=abs(h_ordered).^2;
if total<=0
    h_bb=complex(0); kept=0;
else
    if o.circular_peak_align, n=numel(h_ordered);
    else, n=find(cumsum(energy)>=o.tap_energy_ratio*total,1,'first'); n=max(n,1); end
    h_bb=h_ordered(1:n); kept=sum(abs(h_bb).^2)/total;
end
[~,peak]=max(abs(h_bb).^2); if isempty(peak), peak=1; end
pre=sum(abs(h_bb(1:max(peak-1,0))).^2)/max(sum(abs(h_bb).^2),eps);
meta=struct('source',source,'symbol_rate_hz',o.symbol_rate_hz,'n_fft',numel(h_full), ...
    'f_ref_hz',fc,'idx_f_ref',iref,'f_bb_axis_hz',fbb,'H_baseband',Hbb, ...
    'h_full',h_full,'h_ordered',h_ordered,'tap_count',numel(h_bb),'tap_energy_kept',kept, ...
    'raw_circular_peak_index',raw_peak,'circular_shift_samples',circular_shift, ...
    'circular_window_start',window_start,'circular_window_length',window_length, ...
    'circular_peak_align',o.circular_peak_align, ...
    'peak_index',peak,'discarded_pre_peak_energy_fraction_if_peak_sync',pre, ...
    'frequency_ifft_closure_relative_error',closure, ...
    'physical_delay_resolution_s',physical_resolution_s, ...
    'maximum_unambiguous_delay_s',maximum_unambiguous_delay_s, ...
    'symbol_tap_spacing_s',1/o.symbol_rate_hz, ...
    'zero_padding_note','FFT grid density interpolates symbol-spaced taps; it does not improve 1/B physical resolution.');
end

function [best_start,best_length]=local_shortest_circular_window(power,target)
q=real(power(:)); q=q/max(sum(q),eps); N=numel(q); q2=[q;q];
best_start=1; best_length=N; right=0; running=0;
for left=1:N
    while right<left+N-1 && running<target
        right=right+1; running=running+q2(right);
    end
    if running>=target && right-left+1<best_length
        best_start=left; best_length=right-left+1;
    end
    running=running-q2(left);
    if right<left, right=left; running=0; end
end
end

function o=local_defaults(o)
d=struct('symbol_rate_hz',1000,'n_fft',2048,'f_ref_hz',6000,'tap_energy_ratio',0.999, ...
    'circular_peak_align',false);
n=fieldnames(d); for ii=1:numel(n), if ~isfield(o,n{ii}), o.(n{ii})=d.(n{ii}); end, end
validateattributes(o.symbol_rate_hz,{'numeric'},{'scalar','positive','finite'});
validateattributes(o.n_fft,{'numeric'},{'scalar','integer','positive'});
validateattributes(o.tap_energy_ratio,{'numeric'},{'scalar','>',0,'<=',1});
end
