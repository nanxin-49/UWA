function meta = assert_pe_phase_release_artifact_vertical(value,expected_run_id,label)
%ASSERT_PE_PHASE_RELEASE_ARTIFACT_VERTICAL Reject stale or mixed-run inputs.
if nargin<3, label='artifact'; end
if isstruct(value) && isfield(value,'validation_run_meta')
    meta=value.validation_run_meta;
elseif isstruct(value) && isfield(value,'run_id') && isfield(value,'code_fingerprint')
    meta=value;
else
    error('assert_pe_phase_release_artifact_vertical:MissingMetadata', ...
        '%s has no validation_run_meta.',label);
end
if ~strcmp(meta.run_id,expected_run_id)
    error('assert_pe_phase_release_artifact_vertical:RunIdMismatch', ...
        '%s belongs to run_id %s, expected %s.',label,meta.run_id,expected_run_id);
end
if ~isfield(meta,'phase_reference') || ~strcmp(meta.phase_reference,'direct_dsp')
    error('assert_pe_phase_release_artifact_vertical:PhaseReferenceMismatch', ...
        '%s is not marked direct_dsp.',label);
end
end
