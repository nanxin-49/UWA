# 三维 compact-dual 5/7 kHz 独立频点结果索引

2026-10-26；状态 `COMPLETE_DISCRETE_TWO_POINT_SCREEN_PASS`。固定 5/7 kHz 物理参考与 FullFK 门均为 `CONDITIONAL_G1_PASS`，冻结重建两个中点均低于 `1e-3`。只支持五个离散频点和两个局部重建检查。`G1_target_PM=OPEN`，`WP2_formal=false`。

- [正式报告](../../../reports/deepsea_dual_5k7k_independent_validation_2026-10-26.md)；[总体检查点 v3](checkpoint_v3.json)；[独立审计 v5](independent_audit_v5.json)。早期审计/检查点保留为演进记录。
- 原生复数及全部 FE/FK 变体：[five_seven.mat](five_seven.mat)、[5 kHz 表](f5000/reference_runs.csv)、[7 kHz 表](f7000/reference_runs.csv)；旧 4/6/8 kHz MAT 位于相邻旧目录，没有改写。
- [五频点复数/误差 CSV](five_frequency_native_summary_v2.csv)；[正式 PNG](five_frequency_reflection_and_error_v2.png)和[可编辑 FIG](five_frequency_reflection_and_error_v2.fig)。v1 图及日志保留为版式演进证据。
- [跨平台诊断](compatibility_6000.mat)；[FIG 坐标审计](figure_object_audit_v3.json)；[远端输出 SHA 清单](output.sha256)覆盖 115 项回传文件；[输入 SHA/长度清单](input_manifest.json)覆盖 17 项冻结输入。
- `run_diag6000_01/`、`run_f5000_01/`、`run_f7000_01/`、`run_render_01/`、`run_render2_01/`、`run_figaudit_01/`、`run_figaudit2_01/` 各有完整 `stdout.log`、`stderr.log`、真实退出码、进程墙钟、前后物理 RAM/负载/磁盘记录及检查点。`f5000/`、`f7000/` 各有参考/FullFK/尾界/物理指标和成本表。
- [远端启动脚本](launch_ubuntu_matlab.sh)用于复核执行形式；不应在旧输出目录重复运行。新作业必须重新授权、用新独立目录并重新过资源门。

保存 MATLAB 校验入口：`scripts/validation/validate_deepsea3d_5k7k_frequency.m`、`run_deepsea3d_5k7k_frequency.m`、`check_deepsea3d_ubuntu_compatibility.m`、`audit_deepsea3d_5k7k_figure.m`；报告图入口为 `scripts/reporting/render_deepsea3d_5k7k_frequency.m`；本地外部独立复算为 `scripts/validation/support/audit_deepsea3d_5k7k_saved.py`。原三维 FE/FullFK 验证算子与生产 API 未改。
