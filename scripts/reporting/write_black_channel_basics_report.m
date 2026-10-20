function write_black_channel_basics_report(analysis,results,filename)
%WRITE_BLACK_CHANNEL_BASICS_REPORT Render the plain-language BLACK briefing.
% Numeric values come from the saved analysis/status; prose lives in a template.
[folder,~,~] = fileparts(filename);
if ~isfolder(folder), mkdir(folder); end
reporting_dir = fileparts(mfilename('fullpath'));
template_file = fullfile(reporting_dir,'templates','black_channel_basics_report_zh.md');
text = fileread(template_file);
artifact_prefix = local_relative_prefix(folder,results.run_dir);
figure_names = {'black_01_time_delay.png','black_02_power_delay.png','black_03_doppler.png'};
if all(cellfun(@(name) isfile(fullfile(folder,name)),figure_names))
    figure_prefix = '';
else
    figure_prefix = artifact_prefix;
end
m = analysis.metrics;
values = struct;
values.FIGURE_PREFIX = figure_prefix;
values.ARTIFACT_PREFIX = artifact_prefix;
values.DURATION_S = sprintf('%.3f',analysis.duration_s);
values.RMS_RANGE_MS = sprintf('%.3f–%.3f',min(m.rms_delay_s)*1e3,max(m.rms_delay_s)*1e3);
values.RMS_RANGE_PRECISE_MS = sprintf('%.4f–%.4f',min(m.rms_delay_s)*1e3,max(m.rms_delay_s)*1e3);
values.ENERGY90_RANGE_MS = sprintf('%.3f–%.3f',min(m.energy90_width_s)*1e3,max(m.energy90_width_s)*1e3);
values.RETAINED_RANGE_PERCENT = sprintf('%.2f–%.2f',min(m.retained_energy_fraction)*100,max(m.retained_energy_fraction)*100);
values.PHASE_FIELD = analysis.phase.field;
if analysis.phase.available
    values.FD_RANGE_HZ = sprintf('%.2f 至 %+.2f',min(analysis.phase.frequency_hz(:)),max(analysis.phase.frequency_hz(:)));
else
    values.FD_RANGE_HZ = '无法恢复（缺少相位字段）';
end
values.DATA_DOI = analysis.manifest.doi;
values.DATA_VERSION = analysis.manifest.dataset_version;
values.DATA_LICENSE = analysis.manifest.license;
values.SOURCE_REVISION = analysis.manifest.source_revision;
values.FS_DELAY_KHZ = sprintf('%.0f',analysis.params.fs_delay/1e3);
values.DELAY_STEP_US = sprintf('%.0f',1e6/analysis.params.fs_delay);
values.FS_TIME_HZ = sprintf('%.0f',analysis.params.fs_time);
values.SNAPSHOT_STEP_MS = sprintf('%.0f',1e3/analysis.params.fs_time);
values.WELCH_BIN_HZ = sprintf('%.3f',1/analysis.cfg.welch_window_s);
values.NYQUIST_HZ = sprintf('%.1f',analysis.params.fs_time/2);
values.BULK_HZ = sprintf('%.6f',analysis.phase.bulk_frequency_hz);
values.STORED_SPAN_MS = sprintf('%.2f',(max(analysis.delay_s)-min(analysis.delay_s))*1e3);
values.POST_MAIN_MS = sprintf('%.1f',analysis.delay_s(end)*1e3-0.8);
values.STATUS = results.status;
values.REPEAT_ERROR = sprintf('%.3g',results.repeat_relative_error);
values.SNR_DB = sprintf('%.2f',results.measured_snr_db);
values.ARRIVAL_ERROR_SAMPLES = sprintf('%.6f',max(abs(results.arrival.error_samples)));
values.ARRIVAL_ERROR_US = sprintf('%.2f',max(abs(results.arrival.error_samples))/analysis.params.fs_delay*1e6);
values.MATLAB_VERSION = results.matlab_version;
values.RUN_DIR = strrep(results.run_dir,'\','/');
values.BUF_LEFT = sprintf('%.4f',analysis.cfg.buffer_left);
values.BUF_RIGHT = sprintf('%.4f',analysis.cfg.buffer_right);
metrics_table = sprintf(['| 阵元 | RMS 时延 (ms) | 90%% 区间 (ms) | 窗口内保留能量 (%%) | 公共频移中位数 (Hz) | 公共频移范围 (Hz) |\n' ...
    '|---|---:|---|---:|---:|---|\n']);
for mm=1:height(m)
    row = m(mm,:);
    metrics_table = [metrics_table,sprintf('| %d | %.4f | %.3f–%.3f | %.3f | %.4f | %.4f–%.4f |\n', ...
        row.receiver,row.rms_delay_s*1e3,row.energy90_low_s*1e3,row.energy90_high_s*1e3, ...
        row.retained_energy_fraction*100,row.common_fd_median_hz,row.common_fd_min_hz,row.common_fd_max_hz)]; %#ok<AGROW>
end
values.METRICS_TABLE = strtrim(metrics_table);
names = fieldnames(values);
for ii=1:numel(names)
    text = strrep(text,['{{' names{ii} '}}'],values.(names{ii}));
end
assert(isempty(regexp(text,'\{\{[A-Z_]+\}\}','once')), ...
    'black:ReportTemplate','Unresolved report template token.');
fid = fopen(filename,'w','n','UTF-8');
assert(fid>=0,'black:ReportOpen','Cannot write report.');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid,'%s',text);
end

function prefix = local_relative_prefix(from_dir,to_dir)
from_file = java.io.File(from_dir);
to_file = java.io.File(to_dir);
from_path = from_file.getCanonicalFile().toPath();
to_path = to_file.getCanonicalFile().toPath();
relative = char(from_path.relativize(to_path).toString());
if isempty(relative), prefix = '';
else, prefix = [strrep(relative,'\','/') '/']; end
end
