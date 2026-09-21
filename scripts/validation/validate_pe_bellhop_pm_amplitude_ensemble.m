function ensemble = validate_pe_bellhop_pm_amplitude_ensemble(overrides)
%VALIDATE_PE_BELLHOP_PM_AMPLITUDE_ENSEMBLE Stage 3B PM ensemble entrypoint.
%   Runs the common Stage 3 paired PE/Bellhop workflow with independent
%   Gaussian Fourier coefficient amplitudes drawn from the canonical PM
%   spectral-density band.  This is validation-only; PE/Bellhop core physics
%   and communication code are untouched.
if nargin < 1 || isempty(overrides), overrides = struct(); end
if ~isstruct(overrides) || ~isscalar(overrides)
    error('overrides must be a scalar struct.');
end
script_dir = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(script_dir));
addpath(root);
addpath(script_dir);
addpath(fullfile(script_dir,'support'));
setup_vertical_project();
defaults = struct('ensemble_mode','amplitude','beam_count',5001, ...
    'minimum_seed_count',8, ...
    'output_dir',fullfile(root,'results','validation','pe_bellhop_pm_amplitude_ensemble'), ...
    'report_path',fullfile(root,'reports','pe_bellhop_pm_amplitude_ensemble_report.md'));
names = fieldnames(overrides);
for ii = 1:numel(names)
    % Other common Stage 3 controls are accepted by the delegated entrypoint.
    defaults.(names{ii}) = overrides.(names{ii});
end
ensemble = validate_pe_bellhop_pm_ensemble(defaults);
end
