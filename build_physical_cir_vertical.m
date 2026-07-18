function cir=build_physical_cir_vertical(H_f,f_axis_hz,reference_delay_s,window_type,zero_padding_factor)
%BUILD_PHYSICAL_CIR_VERTICAL Convert an equally spaced physical H(f) to CIR.
%   The project PE convention has negative phase slope for positive delay.
%   REFERENCE_DELAY_S is signed in exp(-i*2*pi*f*reference_delay_s).

arguments
    H_f {mustBeNumeric}
    f_axis_hz (:,1) double {mustBeFinite,mustBePositive}
    reference_delay_s (1,1) double {mustBeFinite} = 0
    window_type char = 'none'
    zero_padding_factor (1,1) double {mustBeInteger,mustBePositive} = 1
end
F=numel(f_axis_hz);
if size(H_f,1)~=F, error('The first dimension of H_f must match f_axis_hz.'); end
if F<2, error('At least two frequency samples are required.'); end
df=diff(f_axis_hz); df0=mean(df);
if max(abs(df-df0))>max(1e-10*abs(df0),1e-9), error('Frequency axis must be equally spaced.'); end
B=f_axis_hz(end)-f_axis_hz(1); Nfft=zero_padding_factor*F;
w=local_window(F,window_type);
phase=exp(-1i*2*pi*f_axis_hz(:)*reference_delay_s);
weighted=H_f.*(phase.*w);
timer=tic; h=ifft(weighted,Nfft,1); elapsed=toc(timer);
Tmax=1/df0; delay=(0:Nfft-1).'*Tmax/Nfft;
cir=struct('h_physical_tau',h,'delay_axis_s',delay, ...
    'frequency_axis_hz',f_axis_hz(:),'reference_delay_s',reference_delay_s, ...
    'pe_phase_convention','negative phase slope for positive physical delay', ...
    'reference_phase_formula','exp(-1i*2*pi*f*reference_delay_s)', ...
    'window_type',lower(window_type),'window_values',w, ...
    'zero_padding_factor',zero_padding_factor,'delta_f_hz',df0, ...
    'bandwidth_hz',B,'physical_delay_resolution_s',1/B, ...
    'delay_axis_spacing_s',Tmax/Nfft,'maximum_unambiguous_delay_s',Tmax, ...
    'zero_padding_note','Zero padding interpolates delay samples; it does not improve 1/B physical resolution.', ...
    'ifft_elapsed_s',elapsed);
end

function w=local_window(F,kind)
n=(0:F-1).'; kind=lower(strtrim(kind));
switch kind
    case {'none','rectangular','rect'}, w=ones(F,1);
    case {'hann','hanning'}, w=0.5-0.5*cos(2*pi*n/(F-1));
    case 'hamming', w=0.54-0.46*cos(2*pi*n/(F-1));
    otherwise, error('window_type must be none, hann, or hamming.');
end
end
