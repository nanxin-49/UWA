# PE–Bellhop验证脚本归档清单

归档日期：2026-09-01  
归档位置：`cash/pe_bellhop_validation_scripts_2026-09-01/`

## 归档原因

本批文件属于阶段性的PE–Bellhop交叉验证、正常坐标近垂直/epsilon诊断和历史绘图入口。当前验证工作准备在新的窗口中重新组织，旧入口包含已经结束、退役或仅用于诊断的路径；继续留在活动脚本目录会与PE生产主线及后续验证方案混淆。

此次归档是可恢复的文件隔离，不代表否定其中的历史结论，也不修改PE marching、FFT约定、Gaussian源定义、海面反射实现或通信主线。`cash/`为Git忽略的本地隔离区，后续正常检索不读取其内容。

## 文件清单

### 历史 PE–Bellhop入口

- `scripts/validation/validate_pe_bellhop_flat_surface_vertical.m`
- `scripts/validation/validate_pe_bellhop_flat_surface_matrix_vertical.m`
- `scripts/validation/validate_pe_bellhop_flat_surface_current_vertical.m`
- `scripts/validation/rerun_pe_bellhop_point_source_postbudget_vertical.m`
- `scripts/validation/support/pe_bellhop_current_run_meta_vertical.m`

### 已退役的正常坐标/epsilon路径

- `scripts/validation/validate_bellhop_rough_surface_near_vertical_limit.m`
- `scripts/validation/validate_bellhop_100m_boundary_source_vertical.m`
- `scripts/validation/validate_bellhop_100m_arrival_field_audit_vertical.m`

### 历史绘图入口

- `scripts/reporting/generate_pe_bellhop_comparison_atlas_vertical.m`
- `scripts/reporting/generate_bellhop_flat_surface_visuals_vertical.m`
- `scripts/reporting/plot_bellhop_rough_surface_near_vertical_rays.m`
- `scripts/reporting/plot_bellhop_rough_surface_full_paths.m`

共归档12个MATLAB脚本（含1个专用helper）。本次不移动报告、结果数据、通用角谱/Weyl辅助函数、当前离底1 m Bellhop代表绘图入口或PE生产代码。
