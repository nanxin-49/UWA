run(fullfile(fileparts(mfilename('fullpath')), 'scripts', 'bootstrap_project.m'));
% Compatibility entrypoint for the channel-only vertical demo.

clear
format compact

run(fullfile(fileparts(mfilename('fullpath')), 'main_vertical.m'));
