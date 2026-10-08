# FK 优化：汇报可视化

日期：2026-10-08。仅整理已有验证证据，不推进主线研究。

面向组会/PPT，用 11 张白底宽屏图讲清“物理图像 → 算法精度与速度 → 物理适用域 → 跨频一致性”。每图独立可用；术语统一为 Incident、Reflected、FullFK、FastFK、BIE。PNG 可直接插入 PPT，PDF 为矢量备份。

## 数据来源

- [Prototype](spectral_fast_kirchhoff_prototype_report.md)：`results/validation/spectral_fast_kirchhoff_prototype/` 的 `all_metrics.csv`、`scaling.csv`、`weak_fields.mat`、`strong_height_fields.mat`。
- [PM pilot](spectral_fast_kirchhoff_pm_pilot_report.md)：`results/validation/spectral_fast_kirchhoff_pm_pilot/` 的 `field_metrics.csv`、`seed_260001_fields.mat`。
- [物理边界](full_kirchhoff_pm_validity_boundary_report.md)：`results/validation/full_kirchhoff_pm_validity_boundary/` 的 `validity_metrics.csv`、`case_data/cutoff_n12_receiver_fields.csv`。
- [跨频一致性](spectral_fast_kirchhoff_multifrequency_consistency_report.md)：`results/validation/spectral_fast_kirchhoff_multifrequency_consistency/` 的 `dense_metrics.csv`、`bie_anchors.csv`、`phase_group_delay.csv`、主 `reflected_cir.csv`；相位语义以 [vertical_comm_guide](../vertical_comm_guide.md) 为准。

没有重测速度，没有运行 PE、FK、BIE 或随机 realization。弱海面 Incident 只按原 prototype 的离散角谱公式从保存系数重建；PM Incident 直接读取保存的 exact input。重建剖面保存在输出根目录的两个 `*_incident_source.csv`。`visualization_audit.json` 记录关键数值，`source_manifest.csv` 记录输入 SHA256。

## 图片与讲述用途

全部文件位于 `results/visualization/fk_optimization_ppt/figures/`，每个下列 PNG 均有同名 PDF。

|图|文件|一句话用途|
|---|---|---|
|1|[ppt_overview_roadmap.png](../results/visualization/fk_optimization_ppt/figures/ppt_overview_roadmap.png)|开场先说明工作链和四阶段结论。|
|2|[ppt_geometry_setup.png](../results/visualization/fk_optimization_ppt/figures/ppt_geometry_setup.png)|解释上行入射、粗糙海面反射和 Rx 平面。|
|3|[ppt_incident_field.png](../results/visualization/fk_optimization_ppt/figures/ppt_incident_field.png)|直观看到到达海面的入射幅度、相位及其压缩一致性。|
|4|[ppt_reflected_field_comparison.png](../results/visualization/fk_optimization_ppt/figures/ppt_reflected_field_comparison.png)|用大幅接收场剖面展示 FullFK/FastFK 重合及独立 BIE 残差。|
|5|[ppt_incident_vs_reflected.png](../results/visualization/fk_optimization_ppt/figures/ppt_incident_vs_reflected.png)|同一 canonical 海面、同一幅度基准下展示反射前后的场变化。|
|6|[ppt_prototype_accuracy_speed.png](../results/visualization/fk_optimization_ppt/figures/ppt_prototype_accuracy_speed.png)|说明确定性案例的误差、rank 收敛和 141 倍热启动收益。|
|7|[ppt_pm_pilot_algorithm_vs_physics.png](../results/visualization/fk_optimization_ppt/figures/ppt_pm_pilot_algorithm_vs_physics.png)|区分“压缩成功”与“高截止物理门失败”。|
|8|[ppt_pm_validity_boundary.png](../results/visualization/fk_optimization_ppt/figures/ppt_pm_validity_boundary.png)|说明单 realization 的条件边界，且高度/Kmax 单独不足以描述适用域。|
|9|[ppt_multifrequency_consistency.png](../results/visualization/fk_optimization_ppt/figures/ppt_multifrequency_consistency.png)|展示频带误差、三个 BIE 锚点、相位/群时延和 diagnostic CIR。|
|10|[ppt_rank_evolution.png](../results/visualization/fk_optimization_ppt/figures/ppt_rank_evolution.png)|说明统一 rank 从 3/6、3/8 到 3/12 的有限阶梯选择。|
|11|[ppt_summary_dashboard.png](../results/visualization/fk_optimization_ppt/figures/ppt_summary_dashboard.png)|最后汇总成果、速度和仍需保留的条件。|

推荐完整顺序：**1 → 2 → 3 → 4 → 5 → 6 → 7 → 8 → 9 → 10 → 11**。若限制为十页，可将图 5 作为备选。

若只保留六张：**1、3、4、6、8、9**。依次讲全局、入射、反射、精度/速度、物理边界、宽带一致性；图 8 已包含物理误差，讲述时口头强调它不同于 FastFK–FullFK。

## 阅读边界与归一化

- 本套几何和场图属于二维 Tx100/Rx3 研究夹具；不代表 Tx99 的 Bellhop 安装验证，底部 Tx 符号也不表示加入了海底边界模型。传播箭头只是示意。
- 图 3 显示 exact FullFK/BIE 共享输入与 FastFK rank3 输入；**不是新增的 PE 独立一致性验收**。不把压缩后的入射场冒称 exact 输入。
- 图 3 各案例幅度以自身 exact Incident 峰值归一；图 4 每列共用 FullFK 核心峰值；图 5 前后共用 Incident 峰值；图 9 CIR 共用 FullFK CIR 峰值。没有幅度拟合、相位对齐或峰位移动。
- 场图显示约 ±40 m 的核心孔径，省略的外部点仍在原始文件中。相位差沿用保存 M99 / Full-only −40 dB 掩膜；原相位剖面使用 unwrap，只为展示，不进行定相位。
- 所有 E 为相对加权复 L2；算法误差以 FullFK 为参照，物理误差以 BIE 为参照。rank 图的 outgoing 阶梯固定 incoming=4，不能与最终选定 3/6 混为同一实验。
- 图 8 端点误差棒为三个数值变化量之和，约 6.8e-6，**不是统计置信区间**。其余多数点是 baseline-only；连接线仅指示已采样 continuation 顺序，不给出插值定理、拟合公式或通用阈值。
- 4–8 kHz 一致性仅针对固定 canonical seed260001/C2 海面、65 点和统一 3/12。三个 FullFK–BIE 锚点的物理结论限于冻结 M99 核心。
- CIR 使用已验证 native-to-DSP 全包络符号桥接和反射分量；65 点矩形 IFFT 的 16 ms 周期窗**未认证所有长尾无混叠**。图中约 3.69 ms 的采样峰不等同于约 3.595 ms 的相位群时延。
- 141x 和 89–112x 是已有 kernel warm wall-time 比值，不能等同于端到端生产加速；跨频最低 warm/cold 为 41.94x/9.83x。所有时间依赖当时硬件/会话。

## 复现与检查

`setup_vertical_project; generate_fk_optimization_ppt_visuals();`

入口：[scripts/reporting/generate_fk_optimization_ppt_visuals.m](../scripts/reporting/generate_fk_optimization_ppt_visuals.m)。自动解析根路径，可从任意 MATLAB cwd 调用；默认重绘本报告输出，不覆盖验证数据。白底、16:9、大字号，PNG 和矢量 PDF 成对导出。Code Analyzer、导出文件尺寸、关键指标复核和人工版面检查见输出 `qa_summary.md`。

仅新增报告/可视化入口及导航，生产模型、接口和方法公式均未改变；因此不需要标量信道/宽带通信回归，也不修改 `vertical_comm_guide.md` 的物理方法正文。`PROJECT_CONTEXT.md` 只追加可视化导航与主线暂停说明。

工作流使用 scientific-visualization 与 matlab 技能；归因：Kassis T., Agarwal V., He Y., Patel D., Brueckner A. M. (2026), [Scientific Agent Skills](https://doi.org/10.48550/arXiv.2609.00065)（当前 arXiv 记录已核对）。此引用为工作流来源，不是声学结论依据。
