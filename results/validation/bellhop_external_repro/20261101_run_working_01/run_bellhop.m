function success = run_bellhop(bellhop_exe, output_dir, base_name)
% run_bellhop - Execute BELLHOP and check results
%
% Inputs:
%   bellhop_exe - path to bellhop executable
%   output_dir - directory containing .env file
%   base_name - base name of env file (without .env extension)
%
% Returns:
%   success - true if BELLHOP completed successfully

    old_dir = pwd;
    cd(output_dir);

    cmd = sprintf('"%s" %s', bellhop_exe, base_name);
    fprintf('Running: %s\n', cmd);

    [status, output] = system(cmd);

    cd(old_dir);

    % Check status
    if status ~= 0
        fprintf('WARNING: BELLHOP returned status %d\n', status);
        if ~isempty(output)
            fprintf('Output: %s\n', output);
        end
    end

    % Check print file
    prt_file = fullfile(output_dir, [base_name '.prt']);
    success = false;

    if isfile(prt_file)
        prt_text = fileread(prt_file);

        if contains(prt_text, 'FATAL ERROR', 'IgnoreCase', true)
            fprintf('ERROR: BELLHOP reported fatal error\n');
            % Print relevant error lines
            lines = splitlines(prt_text);
            for i = 1:length(lines)
                if contains(lines{i}, 'ERROR', 'IgnoreCase', true) || ...
                   contains(lines{i}, 'FATAL', 'IgnoreCase', true)
                    fprintf('  %s\n', lines{i});
                end
            end
        else
            fprintf('BELLHOP completed.\n');
            success = true;
        end
    else
        fprintf('WARNING: Print file not found\n');
    end
end
