# S2c-A 正弦空间机理与受控试验冻结

状态：`S2C_A_DIAGNOSTIC_PROTOCOL_COMPLETE`。新增声学计算 **0**。正式输入仍为S2b `formal_reference_v2`，S0–S2b原场、指标、状态及版本均未改动。

先读 [中文可预览图册](mechanism_atlas.md)：M00四联总览、M01三档几何、三档M02a/b物理响应/模型误差、三档M03全线/局部/真实九点，共11组FIG+PNG。

- [保存复场与绘图MAT](mechanism_baseline/diagnostic_data.mat)：417唯一接收位置×4案例×3方法，411线加6离深度点；没有新场或插值。
- [来源及版本](mechanism_baseline/source_index.csv)、[原生线CSV](mechanism_baseline/fields_line.csv)、[原九点CSV](mechanism_baseline/fields_core.csv)。
- [BIE物理响应](mechanism_baseline/BIE_physical_response.csv)：1668行；[模型逐点残差](mechanism_baseline/model_point_diagnostics.csv)：3336行；[分支持指标](mechanism_baseline/support_metrics.csv)：40行。
- [几何](mechanism_baseline/geometry.csv)、[真实接收x投影几何](mechanism_baseline/geometry_at_receivers.csv)、[保存边界照射](mechanism_baseline/saved_boundary_illumination.csv)、[描述性相关](mechanism_baseline/spatial_associations.csv)。几何与接收x并列不代表局部反射贡献。
- [冻结协议](protocol/design_freeze.json)、[候选解析几何MAT](protocol/candidate_geometry.mat)、[几何统计](protocol/candidate_geometry_stats.csv)、[配置哈希](protocol/configuration_hashes.csv)、[成本/RAM门](protocol/cost_gate.md)。候选均没有声学场。
- [独立审计](audit/independent_audit.json)、[图像检查](audit/visual_review.csv)、[代码检查](audit/code_analyzer.json)、[最终证书](completion_certificate.json)。

`protocol/`保存原离散源谱、411位置与q=1/411、原九点、13个冻结PE倍网格探针、解析候选、实际含C2的统计及已有入射强度加权几何。`provenance/`保存输入/数值参考/代码/文档/输出SHA及源码快照。`execution_attempts/`保留本轮未完成导出和审计导入/文本舍入修正过程，均为零新场；其图与pending记录不是正式交付，不从中挑值。首次失败目录已收纳进本唯一run_id。

只读重绘（目标必须不存在；使用已有MATLAB会话，setup不改变cwd）：

```matlab
setup_vertical_project;
r = 'E:/MISC/CARPE3D_matlab/Explain/results/comparisons/two_dimensional_channel_closeout/20261023_s2c_a_1642';
render_two_dimensional_closeout_s2c_a(fullfile(r,'mechanism_baseline/diagnostic_data.mat'),fullfile(r,'unused_redraw'));
```

`audit_two_dimensional_closeout_s2c_a(r)`从CSV独立复算并重绘；它拒绝覆盖已有audit，因此再审计须先复制本包到新的审核目录并保留原路径来源索引，不能在本正式包上直接重复覆盖。`run_two_dimensional_closeout_s2c_a(root,new_run_id)`可重新生成保存场诊断和解析几何基础协议；当前完整事前执行门以本包已哈希的`design_freeze.json`、13探针CSV与`cost_gate.md`为准，它们不是声学计算入口。

强档PE全线/中心11点/九点复L2分别为6.2911%/0.02098%/0.008051%。中心高度零但斜率非零，BIE物理相位也发生变化；空间关联仅是候选机制。弱/中档FullFK九点仍未在经验5倍规则下分辨。保持`G1_target_PM=OPEN`、`WP2_formal=false`。

共享成本仍扣记1408.660105/7200 s、余5791.339895 s。建议S2c-B子额度≤1800 s且预留S3A规划3000 s；两个候选完整参考可能超额，建议先一个K_HIGH代表案例。S2c-B/S3均须用户单独授权。本阶段结束，不启动下一物理计算，不提交推送。
