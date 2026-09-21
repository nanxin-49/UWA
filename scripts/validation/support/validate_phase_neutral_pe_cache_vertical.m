function audit=validate_phase_neutral_pe_cache_vertical(cache,pm_spec,joint_model,expected)
%VALIDATE_PHASE_NEUTRAL_PE_CACHE_VERTICAL Verify reusable pre-phase artifacts.
arguments
    cache (1,1) struct
    pm_spec (1,1) struct
    joint_model (1,1) struct
    expected (1,1) struct
end
f=cache.f_axis_hz(:); ef=expected.f_axis_hz(:);
audit=struct('phase_neutral',true, ...
    'frequency_axis_max_abs_error_hz',local_axis_error(f,ef), ...
    'pe_shape_matches',cache.cfg.nx==expected.nx&&cache.cfg.ny==expected.ny, ...
    'geometry_matches',cache.cfg.z_tx==expected.z_tx&&cache.cfg.z_rx==expected.z_rx&&cache.cfg.c0==expected.c0, ...
    'pm_shape_matches',pm_spec.nx==expected.pm_nx&&pm_spec.ny==expected.pm_ny, ...
    'mapping_is_central',local_mapping(cache.pm_mapping,cache.cfg,pm_spec), ...
    'operator_is_uniform_cpu_double',strcmp(cache.cfg.env_mode,'uniform')&&~cache.cfg.enable_bubbles&&isa(cache.surface_rx_fr,'double'), ...
    'joint_frequency_matches',local_joint_frequency(joint_model,f), ...
    'reason','PE operator, PM spectrum/mapping, and joint spatial factor precede receiver phase-reference rotation');
audit.pass=audit.frequency_axis_max_abs_error_hz==0&&audit.pe_shape_matches&& ...
    audit.geometry_matches&&audit.pm_shape_matches&&audit.mapping_is_central&& ...
    audit.operator_is_uniform_cpu_double&&audit.joint_frequency_matches;
if ~audit.pass
    error('validate_phase_neutral_pe_cache_vertical:Mismatch', ...
        'Reusable PE/PM/joint cache does not match the formal configuration.');
end
end

function e=local_axis_error(a,b)
if numel(a)~=numel(b), e=Inf; else, e=max(abs(a-b)); end
end

function tf=local_mapping(m,cfg,pm)
iy=floor((pm.ny-cfg.ny)/2)+(1:cfg.ny);
ix=floor((pm.nx-cfg.nx)/2)+(1:cfg.nx);
tf=isequal(m.iy(:),iy(:))&&isequal(m.ix(:),ix(:));
end

function tf=local_joint_frequency(m,f)
if isfield(m,'f_axis_hz'), mf=m.f_axis_hz(:);
elseif isfield(m,'frequency_axis_hz'), mf=m.frequency_axis_hz(:);
else, tf=false; return; end
tf=numel(mf)==numel(f)&&max(abs(mf-f))==0;
end
