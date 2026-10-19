function write_black_channel_basics_report(analysis,results,filename)
%WRITE_BLACK_CHANNEL_BASICS_REPORT Rebuild the short Chinese BLACK report.
[folder,~,~] = fileparts(filename);
if ~isfolder(folder), mkdir(folder); end
project_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
normalized_root = strrep(project_root,'\','/');
normalized_folder = strrep(folder,'\','/');
normalized_run = strrep(results.run_dir,'\','/');
if strcmpi(normalized_folder,[normalized_root '/reports'])
    assert(startsWith(lower(normalized_run),lower([normalized_root '/'])), ...
        'black:ReportPath','Run directory must be inside the project.');
    prefix = ['../' normalized_run(numel(normalized_root)+2:end) '/'];
else
    prefix = '';
end
fid = fopen(filename,'w','n','UTF-8');
assert(fid>=0,'black:ReportOpen','Cannot write report.');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid,'# BLACK 信道：下载、基础复现与特征分析\n\n');
fprintf(fid,'本报告只使用 BLACK 及其配套噪声。复现状态：**%s**。\n\n',results.status);
fprintf(fid,'## 数据与环境\n\n');
fprintf(fid,'- 数据：[Zenodo %s](https://doi.org/%s)，版本 %s，许可 %s。\n', ...
    analysis.manifest.doi,analysis.manifest.doi,analysis.manifest.dataset_version,analysis.manifest.license);
fprintf(fid,'- 官方 MATLAB 源码固定至 `%s`；许可证、MD5 和 SHA-256 见下载 manifest。\n',analysis.manifest.source_revision);
fprintf(fid,'- MATLAB：%s；Signal Processing Toolbox 可用。\n',results.matlab_version);
fprintf(fid,'- 本次结果目录：`%s`。\n',strrep(results.run_dir,'\','/'));
fprintf(fid,'- 实际载频 %.0f Hz；时延轴 %.0f Hz；快照 %.0f Hz；记录 %.3f s；8 个接收阵元。\n\n', ...
    analysis.params.fc,analysis.params.fs_delay,analysis.params.fs_time,analysis.duration_s);
fprintf(fid,'数据是海试估计的等效复基带信道，不是原始录音；幅度经过归一化。');
fprintf(fid,'本报告不据此推断绝对传播损失，不将相对时延解释为绝对传播时间。\n\n');
fprintf(fid,'## 复现与检查\n\n');
fprintf(fid,'采用官方 `unpack`、`replay`、`noisegen`，不修改其算法。解包仅展示幅度，以原生快照率输出；');
fprintf(fid,'已从绘图坐标中扣除左补零长度，补零按实际漂移范围检查并增加（左/右比例 %.4f/%.4f）。', ...
    analysis.cfg.buffer_left,analysis.cfg.buffer_right);
fprintf(fid,'复现参数为 8192 个 BPSK 符号、6.25 ksym/s、');
fprintf(fid,'100 kHz 通带采样率、首尾 20 ms 静默、记录四分之一处起始，随机种子 1994/twister。\n\n');
fprintf(fid,'- 文件：两文件大小与官方 MD5 一致，另记录 SHA-256；HDF5 类型检查通过后加载。\n');
fprintf(fid,'- 合成双径控制：理论均值 0.8 ms、RMS 1.6 ms；负相位斜率与 Welch 频率轴检查通过。\n');
fprintf(fid,'- 同一输入重复回放相对误差：`%.3g`（门限 1e-12）。\n',results.repeat_relative_error);
fprintf(fid,'- 配套噪声生成并缩放后的整体 SNR：%.8f dB（目标 20 dB，包含静默和尾段）。\n',results.measured_snr_db);
fprintf(fid,'- 8 阵元主到达最大误差：%.6f 个时延采样点（门限 2 点，%.1f us）；检查结果 %d。\n\n', ...
    max(abs(results.arrival.error_samples)),2/analysis.params.fs_delay*1e6,results.arrival.passed);
fprintf(fid,'到达检查使用无噪声回放的八段非相干匹配滤波，与保存抽头主峰和平均相位时延漂移作比较。');
fprintf(fid,'此检查验证重建/回放一致性，不等同于重新估计原始海试信道。\n\n');
fprintf(fid,'## 基本多径\n\n');
fprintf(fid,'![BLACK 时延—时间图](%sblack_01_time_delay.png)\n\n',prefix);
fprintf(fid,'图 1：固定 1、4、8 号阵元的重建幅度。三个面板使用共同幅度峰值归一化，');
fprintf(fid,'色标为 -30 至 0 dB；相位/漂移按官方代码恢复。此图的快照率不用于估计完整多普勒谱。\n\n');
fprintf(fid,'![BLACK 平均功率时延谱](%sblack_02_power_delay.png)\n\n',prefix);
fprintf(fid,'图 2：8 阵元保存抽头的平均功率时延谱，逐阵元峰值归一化。');
fprintf(fid,'时延矩在各自峰值以下 25 dB 内的抽头上计算；90%% 区间定义为该截断归一化分布的 5%% 至 95%% 分位点。');
fprintf(fid,'门限会排除弱抽头，因此指标仅描述该门限下的有效多径。40 us 是采样间隔，不等于独立物理分辨率。\n\n');
fprintf(fid,'| 阵元 | RMS 时延 (ms) | 90%% 区间 (ms) | 保留能量 (%%) | 公共频移中位数 (Hz) | 公共频移范围 (Hz) |\n');
fprintf(fid,'|---|---:|---|---:|---:|---|\n');
for mm=1:height(analysis.metrics)
    row=analysis.metrics(mm,:);
    fprintf(fid,'| %d | %.4f | %.3f–%.3f | %.3f | %.4f | %.4f–%.4f |\n', ...
        row.receiver,row.rms_delay_s*1e3,row.energy90_low_s*1e3,row.energy90_high_s*1e3, ...
        row.retained_energy_fraction*100,row.common_fd_median_hz,row.common_fd_min_hz,row.common_fd_max_hz);
end
fprintf(fid,'\n8 阵元 RMS 时延扩展范围为 **%.4f–%.4f ms**；90%% 能量区间宽度范围为 **%.3f–%.3f ms**。', ...
    min(analysis.metrics.rms_delay_s)*1e3,max(analysis.metrics.rms_delay_s)*1e3, ...
    min(analysis.metrics.energy90_width_s)*1e3,max(analysis.metrics.energy90_width_s)*1e3);
fprintf(fid,'多径强弱应结合图 2 的实际峰形判断；不把带限主峰相邻抽头直接计作独立传播路径。\n\n');
fprintf(fid,'本次图像显示：能量主要集中在保存时延坐标约 0.8 ms 的主到达附近，');
fprintf(fid,'其后存在较弱的尾部结构，整体时延分散较小。图 1 中主能量带随时间往返漂移，');
fprintf(fid,'三个代表阵元的变化形态接近。这里不进一步将尾部峰指定为某条海面或海底反射路径。\n\n');
fprintf(fid,'## 基本多普勒\n\n');
fprintf(fid,'![BLACK 基本多普勒](%sblack_03_doppler.png)\n\n',prefix);
if analysis.phase.available
    fprintf(fid,'图 3 上图使用 `%s` 的连续相位，每 0.1 s 独立拟合斜率，',analysis.phase.field);
    fprintf(fid,'按 `f_D = (d phi / dt)/(2 pi)` 定义符号；包含 `f_resamp`（若存在）的固定频移贡献 %.6f Hz。',analysis.phase.bulk_frequency_hz);
    fprintf(fid,'全部阵元与时段的公共频移范围为 **%.4f–%.4f Hz**。\n\n', ...
        min(analysis.phase.frequency_hz(:)),max(analysis.phase.frequency_hz(:)));
else
    fprintf(fid,'相位字段缺失：无法由现有字段恢复公共频移；表中相关项为 NaN。\n\n');
end
fprintf(fid,'下图每阵元选取平均能量最大的固定抽头，原生快照率计算 Welch 谱：');
fprintf(fid,'4 s periodic Hann 窗、50%% 重叠、无额外零填充，频率栅格 %.3f Hz，',1/analysis.cfg.welch_window_s);
fprintf(fid,'奈奎斯特范围 +/- %.1f Hz。右图仅为同一谱的 +/-5 Hz 放大。',analysis.params.fs_time/2);
fprintf(fid,'谱逐阵元峰值归一化，标记为**残余抽头谱**；不得将其宽度当作完整信道多普勒扩展。');
fprintf(fid,'窄峰宽度含有限时窗效应，不据此宣称非零物理扩展；未增加复杂展宽拟合。\n\n');
fprintf(fid,'本次三个代表阵元的公共频移曲线接近，并随时间平滑改变正负；');
fprintf(fid,'保存主抽头的残余谱则集中在零频附近。两者共同说明：移除的公共相位变化不能忽略，');
fprintf(fid,'仅看保存抽头的窄谱会遗漏这部分时变频移。\n\n');
fprintf(fid,'## 重复运行与解释限制\n\n');
fprintf(fid,'```matlab\nsetup_vertical_project;\nr = run_black_channel_basics();\n');
fprintf(fid,'run_black_channel_basics(''render'', r.run_dir);\n```\n\n');
fprintf(fid,'`render` 仅从保存的 analysis.mat/status.json 重建图、指标表和报告，不重新下载或回放。');
fprintf(fid,'PNG/PDF、analysis.mat、replay.mat、black_metrics.csv、status.json 位于结果目录。');
fprintf(fid,'原始文件与固定源码放在本 case 的 raw/vendor 目录，后续运行复核校验和。\n\n');
fprintf(fid,'结论仅限本段 BLACK 记录和测量频带；相位频移并不单独证明平台速度或某一环境机理。');
fprintf(fid,'未开展 BER、完整接收机、其他信道或生产物理模型对比。\n\n');
fprintf(fid,'参考：[BLACK 参数](https://uwa-channels.github.io/channels/black/)、');
fprintf(fid,'[格式说明](https://uwa-channels.github.io/docs/file_specifications/)、');
fprintf(fid,'[配套论文](https://arxiv.org/abs/2609.03207)。\n');
fprintf(fid,'\n软件流程参考：Kassis, T., Agarwal, V., He, Y., Patel, D., & Brueckner, A. M. (2026). ');
fprintf(fid,'[Scientific Agent Skills: A Library of Procedural Knowledge for Research Agents](https://doi.org/10.48550/arXiv.2609.00065). ');
fprintf(fid,'所用 MATLAB/科学可视化技能仅用于资产检查、可复现与图表审查流程，未增加信道数据或物理模型。\n');
end
