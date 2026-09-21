function bellhop_exe = resolve_bellhop_exe_vertical()
%RESOLVE_BELLHOP_EXE_VERTICAL Resolve the external Bellhop executable.
%   The repository never bundles Bellhop. Prefer BELLHOP_EXE; otherwise use
%   the standard Windows binary location below BELLHOP_TOOLBOX_ROOT.

bellhop_exe = strtrim(getenv('BELLHOP_EXE'));
if ~isempty(bellhop_exe)
    return;
end

toolbox_root = strtrim(getenv('BELLHOP_TOOLBOX_ROOT'));
if isempty(toolbox_root)
    return;
end

candidate = fullfile(toolbox_root, 'windows-bin-20201102', 'bellhop.exe');
if isfile(candidate)
    bellhop_exe = candidate;
end
end
