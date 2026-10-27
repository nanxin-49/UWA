# 固定 dual 65 点网格：四频点独立验收结果索引

2026-10-27。状态 `COMPACT_DUAL_65GRID_HOLDOUT_SCREENING_PASS`，限固定三维单链路、4 个新增离散点。正式数值、局限和下一步见[中文报告](../../../reports/deepsea3d_65grid_four_holdout_2026-10-27.md)。不认证连续宽带、CIR、BER、原 PM。

- 预 FE [预测锁](prediction_lock.json)、[四点原生复预测](predictions_pre_FE.csv)、[独立锁审计](prediction_lock_audit.json)。锁 SHA `81f03e3973830c9a02d5441515c299c96c5424268399cb975465322eedcbe8d3`；预测 CSV SHA `badfe577751812388fc519c6872453f1de4dd7b2150799176bb4f0756d368324`。锁生成/审计均早于首个新 MATLAB 进程及真正 FE。
- [原生 MATLAB state](65grid_four.mat)及 `f4750/`、`f7250/`、`f4125/`、`f6562.5/` 保存各阶段 FE、M32、回放、原五配置 FK、真实边界、尾界、成本 CSV/JSON。旧九点只读，6 kHz 新解为零。
- [四点独立复算表](four_point_blind_metrics.csv)、[独立保存数据审计](independent_audit.json)、[MATLAB 工件断言](artifact_assertions.json)、[最终检查点](checkpoint.json)。四点 FE/FullFK 条件资格及 `e_R/e_Δ/a_D<1e−3` 各 4/4 通过。`independent_audit_v1.json` 与 `checkpoint_v1.json` 保存一次时间字段命名纠正前的快照，数值及门未变。
- 最终两图：[原生反射复平面 PNG](four_point_reflection_complex_v3.png)/[FIG](four_point_reflection_complex_v3.fig)，[预测复差与 FE 经验敏感性 PNG](prediction_vs_FE_sensitivity_v3.png)/[FIG](prediction_vs_FE_sensitivity_v3.fig)；[FIG 对象审计](figure_object_audit.json)重开并核对四点坐标、残差及三组离散误差标记。早期版本及失败进程留证。
- `run_f4750_01/` 是设置目录缺失、未开始 FE 的失败进程；`run_f4750_02/`、`run_f7250_01/`、`run_f4125_01/`、`run_f6562.5_01/` 为四个成功物理进程；`run_audit_01/`、`run_render_01/02/03/`、`run_figaudit_01/` 留完整 stdout/stderr、真实退出码、wall、前后 RAM/负载/磁盘与 RSS。95 s MATLAB 版本探测保守扣账但当时完整日志未持久保存。全部 MATLAB 进程累计 `341/720 s`，四物理点各 `<150 s`。
- [远端首批 133 文件 SHA 清单](output.sha256)、[最终图件 SHA](figure_output_v3.sha256)、[冻结输入 20 文件 SHA](remote_input.sha256)、两批回传 bundle 均已验证。大 MAT 与 bundle 是本地结果，不属于生产代码依赖。

复核代码在 `scripts/validation/validate_deepsea3d_65grid_four_frequency.m`、`run_deepsea3d_65grid_four_frequency.m`、`audit_deepsea3d_65grid_four_figures.m`，只读锁/回传审计及打包脚本在 `scripts/validation/support/`，最终绘图入口为 `scripts/reporting/render_deepsea3d_65grid_four_v3.m`。此已完成目录不可重跑或覆盖；新物理频点、129 点、CIR 或通信须另立授权。
