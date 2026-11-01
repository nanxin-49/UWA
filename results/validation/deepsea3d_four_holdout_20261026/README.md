# 四频点固定 dual 盲测：事前冻结与执行索引

2026-10-26。最终状态：`COMPACT_DUAL_FREQUENCY_HOLDOUT_SCREENING_PASS`。四个新频点的 FE 原生真值均独立取得且过原数值资格门；下列预测仍是FE运行前冻结版本。[正式报告](../../../reports/deepsea3d_four_frequency_holdout_2026-10-26.md)先给参考资格与盲测误差，再给限定GO判断。

- [预测锁](prediction_lock.json)、[八个原生复预测](predictions_pre_FE.csv)：A 使用旧 4/6/8 kHz 的 `Q=R/B*`；B 使用已有 4/5/6/7/8 kHz 的分段复线性 `Q`。`η*=0.06251088653876558 m`，`B*` 为固定水平高度的解析像源反射。
- [冻结脚本](../../../scripts/validation/support/freeze_deepsea3d_four_holdout.py)只读取旧五频点与旧锁，不读取新频点FE。CSV SHA-256：`3c3812063728212adde4304f4ccd49fefcb4f225bf769f70772b053df717c2ac`。
- [原资格门](../../../reports/deepsea3d_g1_frequency_protocol.md)和[前次 5/7 kHz 证据](../../../reports/deepsea_dual_5k7k_independent_validation_2026-10-26.md)约束 FE 阶数、真实边界、离节点、网格/域、M32/M40 尾项、出射与完整 FullFK 五变体。两个预测器的反射相对复误差筛查门均为严格小于 `1e-3`。参考未合格者不进入盲测通过数。
- 本任务仅用 `/home/misc/cwx/deepsea3d_four_holdout_20261026/` 独立远端目录；输入以逐文件 SHA 校验。旧 5/7 kHz 完整远端进程墙钟 27/33 s。新四点各预留 85 s（外部 timeout 含 5 s 收尾），另留 40 s 后续 MATLAB 审计，初始保守计划 `380/400 s`，每点低于 `120 s`。启动前按 1520 MiB 新增工作集估计检查物理 `MemAvailable - 1520 MiB >= 30% MemTotal`；不用 swap 放行。
- 若任一参考或资源门失败，保留原日志/检查点并停止后续新物理计算，不放宽门。只有四点均合格且保存原生复场、SHA、成本和独立审计后，才可能作完整 GO/STOP 结论。

## 完成后的证据

- [原生四频点 MAT](four_holdout.mat)，`f4500/`、`f5500/`、`f6500/`、`f7500/` 中每点原生FE、部分和、真实边界、M32/M40、FullFK五配置、成本、复场 CSV 及独立频点检查点。
- [九点原生场与参考资格 CSV](nine_frequency_qualification.csv)、[四点A/B盲测指标 CSV](holdout_blind_metrics.csv)。
- [反射复响应与重建误差 PNG](reflection_complex_and_reconstruction.png)、[预测差与FE敏感性 PNG](prediction_vs_reference_sensitivity.png)；同名 SVG 是可编辑矢量图。图上九点仅为离散样本，未连接成连续频响。
- [独立审计](independent_audit.json)、[完整执行检查点](checkpoint.json)、[远端回传逐文件SHA](remote_output.sha256)。`output_bundle.tar` 的SHA-256为`f995a565f5e528df74f3bb55f5dfecc458b32ba210bb7b378e215842f3dca875`；106/106文件回传一致。
- `run_f4500_01/`、`run_f5500_01/`、`run_f6500_01/`、`run_f7500_01/` 和 `run_audit_01/` 保留完整stdout/stderr、真实退出码、墙钟、物理RAM快照与峰值RSS。五个新进程合计141/400s，无重试。
- 可复核入口：`scripts/validation/validate_deepsea3d_four_holdout_frequency.m`、`run_deepsea3d_four_holdout_frequency.m`、`support/freeze_deepsea3d_four_holdout.py`、`support/package_deepsea3d_four_holdout.py`、`support/audit_deepsea3d_four_holdout_saved.py`。只对本次冻结问题有效，不授权重算或扩频。

二维 `CORE_COMPLETE`、旧五点 `CONDITIONAL_G1_PASS`、Q1 BIE `STOP`、`G1_target_PM=OPEN`、`WP2_formal=false` 均保持原样。此包不授权连续频响、CIR、BER、PM 散射或生产 API 改动。
