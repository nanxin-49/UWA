function out_dir = project_result_dir(varargin)
%PROJECT_RESULT_DIR Return and create a results subdirectory.
project_root = fileparts(fileparts(mfilename('fullpath')));
out_dir = fullfile(project_root, 'results', varargin{:});
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end
end
