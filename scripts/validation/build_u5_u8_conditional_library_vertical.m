%% Build the exact-node U=5/U=8 conditional receiver-channel library
clear; clc;
root=fileparts(fileparts(fileparts(mfilename('fullpath')))); addpath(root);
validation_run_meta=pe_phase_release_run_meta_vertical(root,struct( ...
    'frequency_axis_hz',linspace(4000,8000,64).', ...
    'seed_definition',struct('library_u5',2500001,'library_u8',2500002)));
u5_file=fullfile(root,'results','validation','u5_conditional_channel_f64', ...
    'u5_conditional_channel_model_f64_full.mat');
u8_file=fullfile(root,'results','validation','u8_conditional_channel_f64', ...
    'u8_conditional_channel_model_f64_full.mat');
if ~exist(u5_file,'file'), error('Missing validated U=5 full model: %s',u5_file); end
if ~exist(u8_file,'file'), error('Missing validated U=8 full model: %s',u8_file); end
u5=load(u5_file,'model'); u8=load(u8_file,'model');
assert_pe_phase_release_artifact_vertical(u5.model,validation_run_meta.run_id,'U=5 full model');
assert_pe_phase_release_artifact_vertical(u8.model,validation_run_meta.run_id,'U=8 full model');
fixed_geometry=struct('x_tx_m',0,'y_tx_m',0,'z_tx_m',100, ...
    'x_rx_m',0,'y_rx_m',0,'z_rx_m',3,'surface_depth_m',0);
fixed_pe_config=struct('x_width_m',50,'y_width_m',50,'nx',128,'ny',128, ...
    'sound_speed_mode','uniform','bubbles','off','doppler','off');
options=struct('code_revision',validation_run_meta.code_revision,'fixed_geometry',fixed_geometry, ...
    'fixed_pe_config',fixed_pe_config);
library=build_conditional_channel_library_vertical({u5.model,u8.model},options);
library.validation_run_meta=validation_run_meta;
out_dir=fullfile(root,'results','validation','conditional_channel_library');
if ~exist(out_dir,'dir'), mkdir(out_dir); end
library_file=fullfile(out_dir,'conditional_channel_library_u5_u8_f64.mat');
schema_version=library.schema_version; phase_reference_meta=library.phase_reference_meta;
save(library_file,'library','schema_version','phase_reference_meta','validation_run_meta','-v7.3');
[sample5,meta5]=sample_conditional_channel_library_vertical(library,5,16,2500001,struct('path','auto'));
[sample8,meta8]=sample_conditional_channel_library_vertical(library,8,16,2500002,struct('path','auto'));
unknown_rejected=false; unknown_message='';
try
    sample_conditional_channel_library_vertical(library,6,1,2500003); %#ok<NASGU>
catch ME
    unknown_rejected=true; unknown_message=ME.message;
end
result=struct('library_file',library_file,'supported_nodes',[5,8], ...
    'u5_sample_meta',meta5,'u8_sample_meta',meta8, ...
    'u5_component_error',max(abs(sample5.H_total_f-(u5.model.H_direct_f+u5.model.H_ref_coh_f+sample5.H_ref_sca_f)),[],'all'), ...
    'u8_component_error',max(abs(sample8.H_total_f-(u8.model.H_direct_f+u8.model.H_ref_coh_f+sample8.H_ref_sca_f)),[],'all'), ...
    'unknown_node_rejected',unknown_rejected,'unknown_node_message',unknown_message, ...
    'interpolation_implemented',false);
result.schema_version='2.0.0'; result.phase_reference_meta=library.phase_reference_meta;
result.validation_run_meta=validation_run_meta;
schema_version=result.schema_version; phase_reference_meta=result.phase_reference_meta;
save(fullfile(out_dir,'conditional_channel_library_u5_u8_smoke.mat'),'result', ...
    'schema_version','phase_reference_meta','validation_run_meta','-v7.3');
fprintf('Two-node library saved. Unknown U rejected=%d; component errors %.3e / %.3e.\n', ...
    unknown_rejected,result.u5_component_error,result.u8_component_error);
