function result = validate_bellhop_shd_receiver_range_pairing_vertical(options)
%VALIDATE_BELLHOP_SHD_RECEIVER_RANGE_PAIRING_VERTICAL Freeze 103/97 m pairing.
% This lightweight regression reads the saved curved-wall audit SHD files.
% It proves that the rotated reflected field uses 103 m, while the native
% total and direct fields both use 97 m.  The adjacent 102 m column is kept
% as a deliberate negative control for the historical 2.095-rad error.
arguments
    options.audit_dir (1,:) char = ''
    options.output_dir (1,:) char = ''
    options.fail_on_check (1,1) logical = true
end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
setup_vertical_project;
if isempty(options.audit_dir)
    options.audit_dir=fullfile(root,'results','validation','bellhop_curved_wall_beam_frame_audit');
end
if isempty(options.output_dir), options.output_dir=options.audit_dir; end
if ~exist(options.output_dir,'dir'), mkdir(options.output_dir); end

rotated_file=local_evidence_file(fullfile(options.audit_dir,'cases','rotated'), ...
    'A_0.25_N_161_S_0.1_B_5001.shd','*B5001_rotated.shd');
native_file=local_evidence_file(fullfile(options.audit_dir,'cases','native'), ...
    'A_0.25_N_161_S_0.1_B_5001.shd','*B5001_native.shd');
direct_file=local_evidence_file(fullfile(options.audit_dir,'cases','direct'), ...
    'A_0.25_N_161_S_0.1_B_5001.shd','*B5001_direct.shd');
rotated=read_bellhop_shd_unfolded_vertical(rotated_file);
native=read_bellhop_shd_unfolded_vertical(native_file);
direct=read_bellhop_shd_unfolded_vertical(direct_file);

[p_rotated,rot_meta]=select_bellhop_shd_pressure_at_range_vertical(rotated,103,0);
[p_rotated_wrong,wrong_meta]=select_bellhop_shd_pressure_at_range_vertical(rotated,102,0);
[p_native_total,native_meta]=select_bellhop_shd_pressure_at_range_vertical(native,97,0);
[p_native_direct,direct_meta]=select_bellhop_shd_pressure_at_range_vertical(direct,97,0);
p_native_reflected=p_native_total-p_native_direct;

correct_phase=angle(p_rotated*conj(p_native_reflected));
wrong_phase=angle(p_rotated_wrong*conj(p_native_reflected));
correct_tl=20*log10(abs(p_rotated)/abs(p_native_reflected));
correct_relative=abs(p_rotated-p_native_reflected)/abs(p_native_reflected);
checks=struct;
checks.explicit_ranges=rot_meta.actual_range_m==103 && wrong_meta.actual_range_m==102 && ...
    native_meta.actual_range_m==97 && direct_meta.actual_range_m==97;
checks.expected_columns=rot_meta.range_index==2 && wrong_meta.range_index==1 && ...
    native_meta.range_index==1 && direct_meta.range_index==1;
checks.correct_phase=abs(correct_phase)<=1e-3;
checks.negative_control=abs(wrong_phase)>=1.0 && abs(wrong_phase-correct_phase)>=1.0;
checks.all=all(structfun(@(x)logical(x),checks));

summary=table(103,97,rot_meta.range_index,native_meta.range_index,direct_meta.range_index, ...
    correct_phase,wrong_phase,correct_tl,correct_relative,checks.all, ...
    'VariableNames',{'rotated_range_m','native_range_m','rotated_range_index', ...
    'native_total_range_index','native_direct_range_index','correct_phase_difference_rad', ...
    'wrong_102m_phase_difference_rad','tl_difference_db','complex_relative_error','passed'});
csv_file=fullfile(options.output_dir,'receiver_range_pairing_regression.csv');
writetable(summary,csv_file);
result=struct('checks',checks,'summary',summary,'rotated_meta',rot_meta, ...
    'wrong_102m_meta',wrong_meta,'native_meta',native_meta,'direct_meta',direct_meta, ...
    'files',struct('rotated_shd',rotated_file,'native_shd',native_file, ...
    'direct_shd',direct_file,'summary_csv',csv_file));
if options.fail_on_check && ~checks.all
    error('BellhopRangeSelection:RegressionFailed', ...
        'SHD receiver-range pairing regression failed; inspect %s.',csv_file);
end
disp(summary);
end

function path=local_evidence_file(folder,preferred_name,legacy_pattern)
preferred=fullfile(folder,preferred_name);
if isfile(preferred), path=preferred; return; end
pattern=fullfile(folder,legacy_pattern);
files=dir(pattern);
if numel(files)~=1
    error('BellhopRangeSelection:EvidenceFile', ...
        'Expected exactly one saved SHD matching %s; found %d.',pattern,numel(files));
end
path=fullfile(files(1).folder,files(1).name);
end
