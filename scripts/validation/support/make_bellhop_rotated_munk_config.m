function cfg = make_bellhop_rotated_munk_config(direction)
%MAKE_BELLHOP_ROTATED_MUNK_CONFIG Validation-only rotation of frozen FFP model.
arguments
    direction (1,:) char {mustBeMember(direction,{'up','down'})} = 'up'
end
base=make_munk_5000m_ffp_config();
cfg=struct('frequency_hz',base.frequency_hz,'source_z_m',base.source_depth_m, ...
    'water_depth_m',base.water_depth_m,'munk_parameters',base.munk_parameters, ...
    'bottom_cp_mps',base.bottom_cp_mps,'bottom_rho_gcm3',base.bottom_rho_gcm3, ...
    'bottom_atten_db_wavelength',base.bottom_atten_db_wavelength, ...
    'direction',direction,'transverse_origin_m',1000,'artificial_halfwidth_m',1000, ...
    'ssp_spacing_m',5,'step_m',1,'beam_count',1001,'angle_limits_deg',[-1 1], ...
    'receiver_x_m',0,'run_type','A','source_geometry','R', ...
    'uniform_control',false,'uniform_c_mps',1500);
if strcmp(direction,'up')
    cfg.receiver_range_m=[50 1950 3950 4940];
    cfg.range_limit_m=4950;
else
    cfg.receiver_range_m=[5 10 25 49];
    cfg.range_limit_m=50;
end
cfg.source_model_note='Ideal omnidirectional source; native Bellhop R spreading is diagnostic until rotated 3-D spreading is qualified.';
cfg.phase_note='Native Bellhop phasor preserved; conversion to exp(-i omega t) FFP must be independently verified.';
end