%BOOTSTRAP_PROJECT Configure paths and output folder for project scripts.
%
% Scripts under scripts/ are intended to run from any current directory.
% This bootstrap keeps core APIs on the MATLAB path and routes legacy
% relative outputs into the appropriate results/ subfolder.

bootstrap_file = mfilename('fullpath');
project_root = fileparts(fileparts(bootstrap_file));

addpath(project_root);
addpath(fullfile(project_root, 'scripts'));
addpath(fullfile(project_root, 'scripts', 'validation'));
addpath(fullfile(project_root, 'scripts', 'comparisons'));
addpath(fullfile(project_root, 'scripts', 'experiments'));
addpath(fullfile(project_root, 'scripts', 'reporting'));

stack_info = dbstack('-completenames');
if numel(stack_info) >= 2
    script_file = stack_info(2).file;
else
    script_file = '';
end
[~, script_name] = fileparts(script_file);

script_category = 'experiments';
if contains(script_file, [filesep 'validation' filesep])
    script_category = 'validation';
elseif contains(script_file, [filesep 'comparisons' filesep])
    script_category = 'comparisons';
elseif strcmp(script_name, 'comm_main_vertical_psk') || startsWith(script_name, 'comm_')
    script_category = 'communication';
elseif strcmp(script_name, 'main_vertical') || strcmp(script_name, 'explain_main_vertical') || ...
        strcmp(script_name, 'visualize_surface_reflection_wavefield_vertical')
    script_category = 'visualization';
elseif strcmp(script_name, 'generate_ssa1_group_meeting_report_vertical')
    script_category = 'reports';
end

script_output_dir = fullfile(project_root, 'results', script_category);
if ~exist(script_output_dir, 'dir')
    mkdir(script_output_dir);
end
cd(script_output_dir);
