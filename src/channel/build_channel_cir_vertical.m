function cir = build_channel_cir_vertical(H_f,f_axis_hz,options)
%BUILD_CHANNEL_CIR_VERTICAL Build a CIR with an explicit phase convention.
%   DIRECT_DSP uses MATLAB IFFT. ABSOLUTE_PHYSICAL uses FFT/N under the
%   project p(t)=real(P*exp(-i*omega*t)) convention. A time-origin shift is
%   common to the complete supplied channel; path-wise alignment is not
%   performed here.

arguments
    H_f {mustBeNumeric}
    f_axis_hz (:,1) double {mustBeFinite,mustBePositive}
    options (1,1) struct = struct()
end
options=local_options(options);
F=numel(f_axis_hz);
if size(H_f,1)~=F, error('The first dimension of H_f must match f_axis_hz.'); end
if F<2, error('At least two frequency samples are required.'); end
df=diff(f_axis_hz); df0=mean(df);
if max(abs(df-df0))>max(1e-10*abs(df0),1e-9)
    error('Frequency axis must be equally spaced.');
end
B=f_axis_hz(end)-f_axis_hz(1);
Nfft=options.zero_padding_factor*F;
w=local_window(F,options.window_type);
switch options.input_reference
    case {'direct_dsp','external_dsp_assumed'}
        phase=exp(-1i*2*pi*f_axis_hz(:)*options.time_origin_shift_s);
        transform='ifft';
        timer=tic; h=ifft(H_f.*(phase.*w),Nfft,1); elapsed=toc(timer);
    case 'absolute_physical'
        phase=exp(1i*2*pi*f_axis_hz(:)*options.time_origin_shift_s);
        transform='fft_over_nfft';
        timer=tic; h=fft(H_f.*(phase.*w),Nfft,1)/Nfft; elapsed=toc(timer);
    case 'legacy_reduced'
        error('build_channel_cir_vertical:UnalignedReducedInput', ...
            ['legacy_reduced components do not share a path phase reference. ', ...
            'Convert direct and reflected components before constructing a total CIR.']);
end
Tmax=1/df0;
delay=(0:Nfft-1).'*Tmax/Nfft;
cir=struct('h_tau',h,'h_physical_tau',h,'delay_axis_s',delay, ...
    'frequency_axis_hz',f_axis_hz(:),'input_reference',options.input_reference, ...
    'time_origin_shift_s',options.time_origin_shift_s, ...
    'common_time_shift_only',true,'transform',transform, ...
    'window_type',options.window_type,'window_values',w, ...
    'zero_padding_factor',options.zero_padding_factor,'delta_f_hz',df0, ...
    'bandwidth_hz',B,'physical_delay_resolution_s',1/B, ...
    'delay_axis_spacing_s',Tmax/Nfft,'maximum_unambiguous_delay_s',Tmax, ...
    'zero_padding_note','Zero padding interpolates delay samples; it does not improve 1/B physical resolution.', ...
    'transform_elapsed_s',elapsed,'ifft_elapsed_s',elapsed);
end

function options=local_options(options)
defaults=struct('input_reference','direct_dsp','time_origin_shift_s',0, ...
    'window_type','none','zero_padding_factor',1);
names=fieldnames(defaults);
for ii=1:numel(names)
    if ~isfield(options,names{ii}), options.(names{ii})=defaults.(names{ii}); end
end
options.input_reference=local_choice(options.input_reference, ...
    {'direct_dsp','external_dsp_assumed','absolute_physical','legacy_reduced'},'input_reference');
options.window_type=local_choice(options.window_type, ...
    {'none','rectangular','rect','hann','hanning','hamming'},'window_type');
validateattributes(options.time_origin_shift_s,{'numeric'},{'scalar','finite'});
validateattributes(options.zero_padding_factor,{'numeric'},{'scalar','integer','positive'});
end

function value=local_choice(value,allowed,name)
if isstring(value), value=char(value); end
if ~ischar(value), error('%s must be a char vector or string scalar.',name); end
value=lower(strtrim(value));
if ~ismember(value,allowed), error('Unsupported %s: %s.',name,value); end
end

function w=local_window(F,kind)
n=(0:F-1).';
switch kind
    case {'none','rectangular','rect'}, w=ones(F,1);
    case {'hann','hanning'}, w=0.5-0.5*cos(2*pi*n/(F-1));
    case 'hamming', w=0.54-0.46*cos(2*pi*n/(F-1));
end
end
