% Summarize current SSA surface-model applicability from compact result files.
% This script reads existing validation/calibration outputs only; it does not
% run PE/WAPE or the communication chain.

clear
format compact

result_file = getenv('SSA_APPLICABILITY_RESULT_FILE');
if isempty(result_file)
    result_file = 'summarize_ssa_applicability_vertical_result.mat';
end

rows = struct([]);

ssa2 = local_load_table_if_available( ...
    'compare_ssa1_ssa2_coherent_reflection_vertical_result.mat', 'summary_table');
padding = local_load_table_if_available( ...
    'compare_ssa_conv_padding_vertical_result.mat', 'pair_summary_table');
scale_cal = local_load_table_if_available( ...
    'calibrate_ssa_scatter_scale_vertical_result.mat', 'calibration_table');
freq_corr = local_load_table_if_available( ...
    'validate_ssa_frequency_correlation_vertical_result.mat', 'summary_table');

Hs_values = local_union_numeric_columns({ssa2, scale_cal}, 'Hs_target');
f_values = local_union_numeric_columns({ssa2, scale_cal}, 'f_hz');
if isempty(Hs_values)
    Hs_values = [0.05, 0.2, 0.5, 1.0];
end
if isempty(f_values)
    f_values = [4000, 6000, 8000, 10000];
end

for ih = 1:numel(Hs_values)
    for ifq = 1:numel(f_values)
        Hs = Hs_values(ih);
        f_hz = f_values(ifq);
        row = struct();
        row.Hs_target = Hs;
        row.f_hz = f_hz;
        row.recommended_kernel = "ssa1_geometry";
        row.recommended_scale = local_lookup_numeric(scale_cal, Hs, f_hz, 'recommended_scale', NaN);
        row.scale_residual = local_lookup_numeric(scale_cal, Hs, f_hz, 'best_objective', NaN);
        loss_delta = local_lookup_numeric(ssa2, Hs, f_hz, 'coherent_loss_delta_db', NaN);
        row.ssa2_loss_delta_db = loss_delta;
        row.ssa2_attention = abs(loss_delta) > 0.1;
        row.padding_difference_attention = local_padding_attention(padding, Hs, f_hz);
        row.wideband_frequency_correlation_risk = "needs_review_for_random_wideband_runs";
        row.notes = local_make_note(row);
        rows = local_append_struct(rows, row);
    end
end

applicability_table = struct2table(rows);
report_meta = struct();
report_meta.script = mfilename;
report_meta.created_at = char(datetime('now'));
report_meta.sources = struct( ...
    'ssa2_result', exist('compare_ssa1_ssa2_coherent_reflection_vertical_result.mat', 'file') == 2, ...
    'padding_result', exist('compare_ssa_conv_padding_vertical_result.mat', 'file') == 2, ...
    'scale_calibration_result', exist('calibrate_ssa_scatter_scale_vertical_result.mat', 'file') == 2, ...
    'frequency_correlation_result', exist('validate_ssa_frequency_correlation_vertical_result.mat', 'file') == 2);
if ~isempty(freq_corr)
    report_meta.frequency_correlation_summary = freq_corr;
else
    report_meta.frequency_correlation_summary = table();
end
report_meta.limitations = ['This table is a compact engineering guide from reduced-grid diagnostics. ', ...
    'It is not an experimental scattering calibration or a complete SSA/NLSSA validity map.'];

save(result_file, 'applicability_table', 'report_meta');
writetable(applicability_table, 'summarize_ssa_applicability_vertical_table.csv');
disp(applicability_table)
fprintf('Saved %s\n', result_file);

function T = local_load_table_if_available(file_name, var_name)
if exist(file_name, 'file') ~= 2
    T = table();
    return
end
info = whos('-file', file_name);
if ~any(strcmp({info.name}, var_name))
    T = table();
    return
end
S = load(file_name, var_name);
if isfield(S, var_name)
    T = S.(var_name);
else
    T = table();
end
end

function values = local_union_numeric_columns(tables, column_name)
values = [];
for ii = 1:numel(tables)
    T = tables{ii};
    if ~isempty(T) && any(strcmp(T.Properties.VariableNames, column_name))
        values = [values; T.(column_name)(:)]; %#ok<AGROW>
    end
end
values = unique(values(isfinite(values))).';
end

function value = local_lookup_numeric(T, Hs, f_hz, field_name, default_value)
value = default_value;
if isempty(T) || ~all(ismember({'Hs_target', 'f_hz', field_name}, T.Properties.VariableNames))
    return
end
idx = abs(T.Hs_target - Hs) < 1e-12 & abs(T.f_hz - f_hz) < 1e-9;
if any(idx)
    candidate = T.(field_name)(find(idx, 1, 'first'));
    if isnumeric(candidate) && isfinite(candidate)
        value = candidate;
    end
end
end

function attention = local_padding_attention(T, Hs, f_hz)
attention = false;
if isempty(T)
    return
end
names = T.Properties.VariableNames;
if any(strcmp(names, 'Hs_target'))
    Hs_col = T.Hs_target;
elseif any(strcmp(names, 'sea_hs_target'))
    Hs_col = T.sea_hs_target;
else
    return
end
idx = abs(Hs_col - Hs) < 1e-12;
if any(strcmp(names, 'f_hz'))
    idx = idx & abs(T.f_hz - f_hz) < 1e-9;
end
if ~any(idx)
    return
end
candidate_fields = {'periodic_zero_padded_rel_diff', 'relative_difference', ...
    'P_sca_raw_relative_difference', 'abs_h_reflect_relative_difference', ...
    'rel_diff_H_reflect_f', 'rel_diff_H_f', 'rel_diff_abs_h_reflect', ...
    'rel_diff_E_sca_raw_over_E_inc'};
for ii = 1:numel(candidate_fields)
    if any(strcmp(names, candidate_fields{ii}))
        v = T.(candidate_fields{ii})(idx);
        attention = any(isfinite(v) & abs(v) > 0.1);
        return
    end
end
end

function note = local_make_note(row)
parts = strings(0);
if isfinite(row.recommended_scale)
    parts(end+1) = sprintf('scale=%g from reduced Kirchhoff calibration', row.recommended_scale); %#ok<AGROW>
else
    parts(end+1) = "scale calibration not available";
end
if row.ssa2_attention
    parts(end+1) = "SSA2 coherent delta exceeds 0.1 dB";
else
    parts(end+1) = "SSA2 coherent delta small or unavailable";
end
if row.padding_difference_attention
    parts(end+1) = "periodic/zero_padded difference needs review";
end
note = strjoin(parts, '; ');
end

function rows = local_append_struct(rows, row)
if isempty(rows)
    rows = row;
else
    rows(end+1) = row; %#ok<AGROW>
end
end
