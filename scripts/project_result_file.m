function file_path = project_result_file(varargin)
%PROJECT_RESULT_FILE Return a path under the project results directory.
if nargin < 1
    error('project_result_file requires at least one path component.');
end
parts = varargin;
file_name = parts{end};
dir_parts = parts(1:end-1);
file_path = fullfile(project_result_dir(dir_parts{:}), file_name);
end
