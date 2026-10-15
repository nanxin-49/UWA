# BIE最小验证与解析频率重建产物

[正式计划](../../../reports/deep_sea_channel_research_plan.md)、[预声明协议](../../../reports/deepsea_bie_reconstruction_protocol_2026-10-14.md)、[BIE报告](../../../reports/deepsea3d_bie_minimum_report.md)、[重建报告](../../../reports/deepsea_frequency_reconstruction_report.md)。G1_target_PM=OPEN/WP2_formal=false。

## 状态与证据

- `bie_minimum.mat/json`：六案Gaussian制造逆问题、完整密度、rhs、残差历史、独立边界/场、原始近跳跃失败、实际图面Cauchy制造数据、RAM时间采样。六案线性通过，但全部reference筛查失败。
- `bie_manufactured.csv`、`bie_near_jump.csv`、`bie_cauchy.csv`：通过与失败均保留，不是PM参考。
- `bie_window_diagnosis.mat/json`、`bie_window.csv`：固定h窗口、真实C2图面PV和法向/面元差分，仍未过真实边界门；没有未知dual或PM求解。
- `frequency_reconstruction.mat/json`：同257/513源点，quarter/midpoint/three-quarter及深衰落query、原失败linear、已知几何分量重建、dense频响/CIR。
- `reconstruction.csv`、`delay_sensitivity.csv`：数值筛查与错误时延诊断，严格限analytic_control。
- `data_audit.json`：重读原始MAT的门分类、独立反射分量/PDE/单射线出射检查；PDE通过不掩盖边界失败。
- `unit_tests.csv`、`regression_tests.csv`、`code_analysis.json`、`execution_log.txt`：6项新unit、4项冻结回归、静态检查及执行原始记录。
- `frozen_before.json`、`run_manifest.json`、`executed_code/`、`final_audit.json`：全部既有源/结果/报告哈希、执行前代码/预算、执行源快照及最终审计。源快照不在setup路径中；最终代码只清理两处多余分析器抑制注释，计算内容不变。
- `state.json`：当前双轨状态；此前多项式失败、通用插值未通过仍保留。
- `pm_representation_geometry.json`：冻结PM全域下界来源/哈希，未来a=-4m的几何匹配提议；本轮a=-1m仅平面/dual控制，未运行PM。
- `completion_coverage.json`：逐项交付核对，区分完成的否决决策与仍未具备的PM参考资格。

## 复现

仅共享MATLAB。先`addpath('E:\MISC\CARPE3D_matlab\Explain'); setup_vertical_project`，不依赖当前目录。当前MAT已存在，求解入口会拒绝覆盖；需要复跑时使用独立仓库副本或新case输出路径，重新冻结manifest，不删除本轮负结果。审计在本case重建摘要，不改变原始MAT。

```matlab
validate_deepsea_frequency_reconstruction;
validate_deepsea_bie_minimum;
diagnose_deepsea_bie_windows;
audit_deepsea_bie_reconstruction;
runner=matlab.unittest.TestRunner.withTextOutput;
runner.run(matlab.unittest.TestSuite.fromFile(which('test_deepsea_bie_reconstruction')));
runner.run(matlab.unittest.TestSuite.fromFile(which('test_deepsea3d_mapped_direct')));
```

MATLAB R2025b Update6；所有新代码只用基础MATLAB，监测调用Windows自带.NET进程属性。无新toolbox依赖。BIE辅助平面与图面外区处理仍需独立验证后才可推广，图面制造恒等式不等于未知散射解或PM物理认证。
