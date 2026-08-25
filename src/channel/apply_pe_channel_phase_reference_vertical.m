function [components,meta] = apply_pe_channel_phase_reference_vertical( ...
    f_axis_hz,reduced_components,geometry,target_reference)
%APPLY_PE_CHANNEL_PHASE_REFERENCE_VERTICAL Align reduced PE path phases.
%   The PE march stores a reduced envelope for each longitudinal path. This
%   function restores a common receiver reference after propagation; it does
%   not modify the PE operator, surface screen, or adjoint weights.

arguments
    f_axis_hz (:,1) double {mustBeFinite,mustBePositive}
    reduced_components (1,1) struct
    geometry (1,1) struct
    target_reference = 'direct_dsp'
end

required_geometry = {'z_tx','z_rx','c0'};
for ii = 1:numel(required_geometry)
    if ~isfield(geometry,required_geometry{ii})
        error('apply_pe_channel_phase_reference_vertical:MissingGeometry', ...
            'geometry.%s is required.',required_geometry{ii});
    end
end
if ~isfield(geometry,'z_surface'), geometry.z_surface = 0; end
validateattributes(geometry.z_tx,{'numeric'},{'scalar','finite'});
validateattributes(geometry.z_rx,{'numeric'},{'scalar','finite'});
validateattributes(geometry.z_surface,{'numeric'},{'scalar','finite'});
validateattributes(geometry.c0,{'numeric'},{'scalar','finite','positive'});
if geometry.z_rx < geometry.z_surface || geometry.z_rx >= geometry.z_tx
    error('apply_pe_channel_phase_reference_vertical:InvalidGeometry', ...
        'Require z_surface <= z_rx < z_tx.');
end
if ~isfield(reduced_components,'direct_f')
    error('apply_pe_channel_phase_reference_vertical:MissingDirect', ...
        'reduced_components.direct_f is required.');
end

target_reference = local_choice(target_reference);
F = numel(f_axis_hz);
direct_reduced = local_frequency_array(reduced_components.direct_f,F,'direct_f');
if size(direct_reduced,2) ~= 1
    error('apply_pe_channel_phase_reference_vertical:DirectShape', ...
        'direct_f must contain one deterministic column.');
end

tau_direct = (geometry.z_tx-geometry.z_rx)/geometry.c0;
tau_reflect = ((geometry.z_tx-geometry.z_surface)+ ...
    (geometry.z_rx-geometry.z_surface))/geometry.c0;
relative_delay = tau_reflect-tau_direct;
physical_direct = exp(1i*2*pi*f_axis_hz*tau_direct);
physical_reflect = exp(1i*2*pi*f_axis_hz*tau_reflect);
dsp_direct = ones(F,1);
dsp_reflect = exp(-1i*2*pi*f_axis_hz*relative_delay);

switch target_reference
    case 'legacy_reduced'
        direct_factor = ones(F,1);
        reflect_factor = ones(F,1);
    case 'direct_dsp'
        direct_factor = dsp_direct;
        reflect_factor = dsp_reflect;
    case 'absolute_physical'
        direct_factor = physical_direct;
        reflect_factor = physical_reflect;
end

components = struct();
components.direct_f = direct_factor.*direct_reduced;
optional = {'reflect_coh_f','reflect_sca_fm','reflect_fm'};
for ii = 1:numel(optional)
    name = optional{ii};
    if isfield(reduced_components,name) && ~isempty(reduced_components.(name))
        value = local_frequency_array(reduced_components.(name),F,name);
        components.(name) = reflect_factor.*value;
    end
end
if ~isfield(components,'reflect_fm')
    if isfield(components,'reflect_coh_f') && isfield(components,'reflect_sca_fm')
        components.reflect_fm = components.reflect_coh_f+components.reflect_sca_fm;
    elseif isfield(components,'reflect_coh_f')
        components.reflect_fm = components.reflect_coh_f;
    elseif isfield(components,'reflect_sca_fm')
        components.reflect_fm = components.reflect_sca_fm;
    else
        components.reflect_fm = complex(zeros(F,1));
    end
end
components.total_fm = components.direct_f+components.reflect_fm;

meta = struct( ...
    'schema_version','1.0.0', ...
    'target_reference',target_reference, ...
    'physical_time_convention','p(t)=real(P(f)*exp(-1i*2*pi*f*t))', ...
    'dsp_synthesis_convention','MATLAB ifft; positive relative delay has negative frequency-phase slope', ...
    'reduced_operator_convention','PE stores exp(1i*d*(kz-k0)) reduced envelopes', ...
    'z_tx_m',geometry.z_tx,'z_rx_m',geometry.z_rx, ...
    'z_surface_m',geometry.z_surface,'c0_mps',geometry.c0, ...
    'tau_direct_s',tau_direct,'tau_reflect_s',tau_reflect, ...
    'relative_delay_s',relative_delay, ...
    'direct_factor_f',direct_factor,'reflect_factor_f',reflect_factor, ...
    'direct_dsp_factor_f',dsp_direct,'reflect_dsp_factor_f',dsp_reflect, ...
    'direct_physical_factor_f',physical_direct, ...
    'reflect_physical_factor_f',physical_reflect, ...
    'alignment_policy','Deterministic post-propagation reference conversion; no per-realization peak alignment.');
end

function value = local_frequency_array(value,F,name)
if ~isnumeric(value) || size(value,1) ~= F
    error('apply_pe_channel_phase_reference_vertical:FrequencyShape', ...
        '%s must have frequency as its first dimension (%d rows).',name,F);
end
value = double(value);
end

function value = local_choice(value)
if isstring(value), value = char(value); end
if ~ischar(value)
    error('apply_pe_channel_phase_reference_vertical:ReferenceType', ...
        'target_reference must be a char vector or string scalar.');
end
value = lower(strtrim(value));
if ~ismember(value,{'legacy_reduced','direct_dsp','absolute_physical'})
    error('apply_pe_channel_phase_reference_vertical:ReferenceValue', ...
        'target_reference must be legacy_reduced, direct_dsp, or absolute_physical.');
end
end
