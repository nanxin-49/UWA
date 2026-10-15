# 2026-10-14 双轨检查点产物

正式依据：[研究计划](../../../reports/deep_sea_channel_research_plan.md)。[协议](../../../reports/deepsea_two_track_protocol_2026-10-14.md)、[A报告](../../../reports/deepsea3d_preconditioner_decision_report.md)、[B报告](../../../reports/deepsea_singlelink_readiness_report.md)。G1 OPEN/WP2=false。

- `frozen_before.json`：原src/scripts及全部既有results/validation文件SHA256。scripts/README导航是允许变更项，物理代码及原结果逐一检查不变。
- `run_manifest.json`：耗时求解前冻结的预算、方法、代码哈希和物理/数据seed；`final_audit.json`记录结束状态和最终哈希。计划执行后新增状态注记的哈希与执行前版本分开保留。
- `preconditioner_decision.mat/json`：小规模制造解与RAM预检；不是正式PM结果。
- `right_flat/dual/quarter/pm.mat/json`：正式一次干预，MAT保留完整W、原生压力、配置及原始残差；PM未通过，不能进入统计训练或PM BER。
- `A_comparison.csv`、`data_audit.json`：冻结算子重算、同海面/网格检查及停止状态。
- `singlelink_readiness.mat/json`：解析宽带三类复场/CIR、发射bit、seed、离散dual复场；`CIR_checks.csv`/`AWGN_checks.csv`保留正负检查。513点线性插值未达1e-3筛查。
- `BIE_cost_estimate.csv/json`：非对角算术微基准及稠密/分块容量估计，没有BIE求解。
- `execution_log.txt`、`regression_tests.csv`、`code_analysis.json`：真实执行、标准入口遮蔽失败及后续4/4回归、静态检查。
- `plan_merge.diff`、`plan_structure_audit.json`：原计划完整备份到合并版本差异和文档核对；最后执行状态注记也纳入最终差异。

共享MATLAB执行，先 `addpath('E:\MISC\CARPE3D_matlab\Explain'); setup_vertical_project`，不依赖cwd。
在未有本轮结果的独立副本中，顺序调用：

```matlab
validate_deepsea3d_preconditioner_decision;
validate_deepsea3d_preconditioner_checkpoint('flat');
validate_deepsea3d_preconditioner_checkpoint('dual');
validate_deepsea3d_preconditioner_checkpoint('quarter');
validate_deepsea3d_preconditioner_checkpoint('pm'); % controls/RAM/预算门先通过
validate_deepsea_singlelink_readiness;
audit_deepsea_two_track;
estimate_deepsea3d_bie_cost;
suite=matlab.unittest.TestSuite.fromFile(which('test_deepsea3d_mapped_direct'));
runner=matlab.unittest.TestRunner.withTextOutput;
tests=runner.run(suite);
```

当前case已存在MAT，求解入口会拒绝覆盖。未来复跑应使用独立副本或明确的新输出case，并重新冻结预算/manifest，不删除或改写本轮失败证据。AWGN可用已保存bit、code与noise seed完整重放；不能把当前人工控制移交为PM样本。
