function assert_channel_dsp_metadata_vertical(meta)
%ASSERT_CHANNEL_DSP_METADATA_VERTICAL Reject ambiguous project DSP contracts.
% Missing metadata is an external pass-through policy of the calling API,
% never evidence that a historical project array has been migrated.
arguments
    meta (1,1) struct
end
if ~isfield(meta,'target_reference') || ...
        ~ismember(char(meta.target_reference),{'direct_dsp','external_dsp_assumed'})
    error('assert_channel_dsp_metadata_vertical:Reference','DSP-ready input required.');
end
if strcmp(char(meta.target_reference),'external_dsp_assumed'), return; end
if ~isfield(meta,'phase_contract_version') || ...
        ~strcmp(meta.phase_contract_version,'2.0.0') || ...
        ~isfield(meta,'output_representation') || ...
        ~strcmp(meta.output_representation,'direct_dsp')
    error('assert_channel_dsp_metadata_vertical:LegacyContract', ...
        'Unversioned/old project direct_dsp is not native-to-DSP v2. Migrate using explicit source representation and geometry.');
end
end
