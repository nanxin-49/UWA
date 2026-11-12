function cir=build_physical_cir_vertical(H_f,f_axis_hz,reference_delay_s,window_type,zero_padding_factor)
%BUILD_PHYSICAL_CIR_VERTICAL Legacy wrapper for direct-DSP CIR synthesis.
%   REFERENCE_DELAY_S is retained as one common time-origin shift. It must
%   not be used to align direct and reflected PE paths separately.

arguments
    H_f {mustBeNumeric}
    f_axis_hz (:,1) double {mustBeFinite,mustBePositive}
    reference_delay_s (1,1) double {mustBeFinite} = 0
    window_type char = 'none'
    zero_padding_factor (1,1) double {mustBeInteger,mustBePositive} = 1
end
options=struct('input_reference','direct_dsp', ...
    'time_origin_shift_s',reference_delay_s,'window_type',window_type, ...
    'zero_padding_factor',zero_padding_factor);
cir=build_channel_cir_vertical(H_f,f_axis_hz,options);
cir.reference_delay_s=reference_delay_s;
cir.legacy_api=true;
cir.reference_phase_formula='exp(-1i*2*pi*f*reference_delay_s), common shift only';
end
