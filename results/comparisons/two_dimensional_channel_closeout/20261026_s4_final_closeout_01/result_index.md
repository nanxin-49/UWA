# S4 二维核心收尾：正式结果与可追溯索引

状态：`CORE_COMPLETE`（S3A固定PM `3/3`）、`S3B_NOT_RUN_BUDGET_LIMITED`、`S5_NOT_READY`。权威解释见[唯一主报告](../../../../reports/two_dimensional_channel_closeout_report.md)顶部S4摘要。本包只有保存结果索引、只读核查和阶段检查点；新增声学场为0。以下相对路径均从仓库根目录解释，`cash/`不在证据集内。

| 核心主张 | 正式MAT/CSV与配置入口 | SHA/独立审计 | 图件 |
|---|---|---|---|
| S2a三核心方法同源、平面表示；Bellhop源未闭合 | `results/comparisons/two_dimensional_channel_closeout/20261023_s2a_final/data/source_qualification.csv`、`source_manifest.csv`、`plane_control.csv` | 同包`audit/independent_audit.json`；本轮关键源资格CSV SHA `14275e6e2a8cfdc877ef3610231297c5c303e5dc969c867c33b118c2a33a2a1a` | 同包`figures/F02_flat_control.png` |
| S2b平面/三正弦九点与411线；弱中FullFK九点未分辨 | `results/comparisons/two_dimensional_channel_closeout/20261023_s2b_1305/formal_reference_v2/fields_core.csv`、`fields_line.csv`、`array_metrics.csv`、`line_metrics.csv`、`resolutions.csv`及`final_certificate.mat` | 同包`provenance/final_output_hashes.csv`本轮247/247；`audit/independent_audit.json`、`audit/resolution_audit.json`；配置/代码SHA在`provenance/` | 同包`figures/F01–F07`，以真实文件名为准 |
| S2c-A空间支持与BIE物理响应 | `results/comparisons/two_dimensional_channel_closeout/20261023_s2c_a_1642/mechanism_baseline/`中的保存场和支持指标；`protocol/`冻结设置 | 同包`provenance/final_output_hashes.csv`本轮203/203；`audit/independent_audit.json` | 同包`mechanism_atlas.md`，M00–M03 |
| K_HIGH：BIE物理变化、两法模型残差及未分辨五点变化 | `results/comparisons/two_dimensional_channel_closeout/20261023_s2c_b_khigh_02/analysis/raw_fields.csv`、`support_metrics.csv`、`paired_changes.csv`及`protocol/` | 同包`provenance/P5_final_artifact_hashes_v6.csv`本轮132/132；`audit/final_independent_audit_v4.json` | 同包`mechanism_atlas.md`，正式M04a–d |
| PHI_PLUS90：真场齐全、参考受混合平台求解限制 | `results/comparisons/two_dimensional_channel_closeout/20261023_s2c_b_phi_plus90_01/analysis/support_metrics.csv`、`raw_fields.csv`及`protocol.md`；远端迁移见`remote_recovery/remote_return_manifest.json` | 同包`audit/final_independent_audit_phi_v2.json`、`audit/remote_saved_density_boundary_20261025.json`；本轮关键CSV、最终审计及成本文件SHA见`s4_saved_evidence_audit.json` | 同包`mechanism_atlas.md`及`figures_final_complete_v2/`的M05a–e |
| S3A seed260001/260002/260003：三方法九点与整线六行对照 | `results/comparisons/two_dimensional_channel_closeout/20261025_s3a_seed260001_01/`、`20261025_s3a_seed260002_01/`、`20261026_s3a_seed260003_01/`各自`analysis_v2/`或`analysis/`、`input/`、`checkpoint.json`；统一表在第三包`three_seed_closeout/three_seed_comparison.csv` | 三包各自`audit/SHA256SUMS*.txt`本轮正式数据与冻结快照24/24、23/23、130/130；各自`audit/independent_saved_data_audit_v3.json`，第三包`three_seed_closeout/independent_audit.json` | 第三包`three_seed_closeout/F08_S3A_three_seed_summary.png`；各seed F01–F03 |

SHA限定：seed260001的旧清单另有8条指向后续演进的活脚本路径，本轮未把这些现行字节与历史冻结SHA等同；24条冻结输入/结果全匹配。其余两seed按各自最终清单全匹配。上述现存正式包未移动或覆盖。S4的只读检查重新计算文件SHA、读取最终审计状态与支持/分母字段；原MAT/CSV数值重算、FIG重开及PNG逐像素重绘采用各阶段已经通过的独立审计，不在S4重跑。

成本与退出：S3A保守扣记2010/3000 s；共享累计5411.424104/7506.601 s，余额2095.176896 s。S3B新增五个按1100 s/seed需5500 s，缺3404.823104 s；未抽取后五个seed。S5宽带源、时间原点、频率采样与预算门未过。S4阶段门和未闭合项逐项记录于`checkpoint.json`。
