# R2两Seed保存数据独立审查、图件稳定化与最终封存

2026-11-17。**seed260002和260003均已在新的独立工作区完成全套保存数值审计、严格FIG/PNG审计及分别封存；R1+R2三Seed配对汇总通过。** R2本次限定范围完成，不自动进入R3或宽带通信。

本轮只使用现有本地MATLAB MCP R2025b Update6，临时单线程、结束恢复原路径/线程；BIE、PS/FullFK及远端/正式执行入口由六个哨兵遮蔽。新增声学求解0、SSH连接0、远程MATLAB启动0、Git提交/推送0。物理公式、原数值门、PM系数、411坐标和原复场不改。

交付根目录：`results/comparisons/piston_2d_pm/r2_saved_closeout_20261117_01/`。原三份回收的287/229/239项SHA及输入/阶段/原最终链在前后两次检查中都有效。seed260002原FAILED_STOPPED、退出1、缺原final seal及144170像素失败证据保留；新本地封存是新增证据链，不追改原运行。

## 正式指标与资格

| Seed | 模型 | 复场相对L2 (%) | 5倍经验阈值 (%) | 相位项占平方残差 (%) | 复残差峰x (m) | 资格 |
|---|---|---:|---:|---:|---:|---|
| 260001 | PS | 15.92349807 | 0.000189577 | 99.769088 | 21.0937500 | 既有R1正式审查＋本轮SHA复核 |
| 260001 | FullFK | 0.07472082 | 0.000156441 | 99.943711 | -18.5546875 | 既有R1正式审查＋本轮SHA复核 |
| 260002 | PS | 15.70519754 | 0.000140009 | 99.918940 | 33.7890625 | 本轮完整独立审计＋新封存 |
| 260002 | FullFK | 0.04906342 | 0.000128004 | 99.959181 | -26.5625000 | 本轮完整独立审计＋新封存 |
| 260003 | PS | 18.02830013 | 0.000179196 | 99.875852 | -39.6484375 | 本轮完整独立审计＋新封存 |
| 260003 | FullFK | 0.05568423 | 0.000148500 | 99.971724 | 7.0312500 | 本轮完整独立审计＋新封存 |

六项模型差均超过原5倍经验分辨阈值。三个实现上PS误差为15.7052%–18.0283%，完整FullFK为0.0490634%–0.0747208%；PS/FullFK误差比分别213.107、320.100、323.760。主要误差仍是相位项。完整FullFK小而可分辨，不能写成严格一致；这只是三个固定二维PM实现，不是总体置信区间或实测3-D换能器结论。

- [正式六行汇总表](../results/comparisons/piston_2d_pm/r2_saved_closeout_20261117_01/three_seed_summary/qualified_metrics.csv)：包含formal_pass、qualification、historical_remote_status、幅相/峰位/逐点有效数及相位能量。
- [资格登记](../results/comparisons/piston_2d_pm/r2_saved_closeout_20261117_01/three_seed_summary/qualification_registry.csv)：逐Seed绑定正式审计/封存路径与SHA。R1复用2026-11-16正式审查并复核来源SHA，本轮没有重做R1声学审计。
- [原保存指标快照](../results/comparisons/piston_2d_pm/r2_saved_closeout_20261117_01/three_seed_summary/original_stored_metrics.csv)：六行都标ORIGINAL_SAVED_SNAPSHOT_NOT_A_NEW_AUDIT，仅用于来源对照，不以原回执代替本轮独立审计。seed260003已从上轮“仅原回执”升级为本轮完整审计；当前正式主表没有未审计的R2行。

## 独立数值审计覆盖

两Seed都重新审计已有C0–C4保存密度，不求解新密度。每层核对应源/rhs、GL8几何/法向/弧长/窗权、密度接收重建、MAT/CSV、分量闭合、完整平台中点支持及残差，另做19点GL112独立边界重建。中点数保持340/408/408/476/476，未减点。综合审计另核独立源/梯度、PS直接Green、完整FullFK双项、411点复场与幅相。

| 综合独立指标 | seed260002 | seed260003 | 原门 |
|---|---:|---:|---:|
| BIE密度接收相对重建差 | 1.8084275e-14 | 1.6787021e-14 | <1e-11 |
| MAT/CSV最大绝对差 (Pa) | 6.9095486e-14 | 7.0340142e-14 | <1e-11 Pa |
| 独立边界相对残差 | 6.8124865e-08 | 1.1337542e-07 | <1e-4 |
| PS直接Green相对重建差 | 8.0287576e-09 | 2.0160562e-08 | <1e-4 |
| 完整FullFK双项相对重建差 | 1.1630512e-14 | 1.0669201e-14 | <1e-11 |

每Seed五层审计、综合审计、四敏感性/参考选择/不确定性重建及411点相位掩码检查均0失败，三个Seed均411/411点满足固定相位掩码。

| BIE敏感性（全线/原四段最坏复L2） | seed260002 | seed260003 |
|---|---:|---:|
| 空间10→12 | 1.270856987e-07 | 1.384817036e-07 |
| 空间12→14 | 1.423734211e-07 | 1.603244197e-07 |
| 自求积64→96 | 9.348240446e-08 | 1.024770654e-07 |
| 窗口60→70 | 2.271371792e-14 | 1.990009276e-14 |

四项均保留原<1e-4门。空间差分别上升12.02946%和15.77300%，原获批v2继续将其作为诊断，未宣称原v1严格单调通过；R1原+8.18133%失败历史也保留。参考始终为C3，C4仅验证窗宽。

误差公式不变：`E=||u_model-u_C3||2/||u_C3||2`；`bunc=|C0-C1|+|C1-C3|+|C2-C1|+|C1-C4|`；分辨条件仍为`E>5||bunc+model_unc||2/||u_C3||2`。PS/FullFK原数值敏感性贡献均保留。逐点残差统一除以本Seed整线BIE RMS，不用局部接近零声压作分母。相位掩码仍是`|BIE|>=0.01*RMS`，幅相平方残差分解仍为`|u-b|²=(|u|-|b|)²+4|u||b|sin²(Δphi/2)`。经验分辨率不是严格误差界或统计置信区间。

## 图件稳定化与严格像素门

采用新增的保存数据绘图器：独立坐标轴，固定像素位置/尺寸、坐标限、刻度/标签、指数显示和图例位置；图例在数据区域外，不覆盖误差峰。FIG重开时显式恢复1800×1400原画布，再进行drawnow和导出。本机实测openfig会把保存画布从1800×1400缩成1800×1041；恢复原尺寸后像素差归零。新流程同时避免tiledlayout自动布局，未改变误差数组或PNG比较标准。原远端F02内部首次变化仍只保留历史诊断，不将本机画布证据冒充远端过程的完整复现。

两Seed各核原FIG12条曲线及新FIG12条曲线；汇总新FIG另核12条，共60条语义曲线检查。每条都用保存复场和保存RMS/掩码逐点核对X/Y。4张单Seed误差图加1张三Seed误差图，每张独立重开两次并重导出，**10次PNG比较的尺寸相同、不同通道值数均为0**，且显式布局前后完全一致。没有用复制PNG、放宽像素门、平滑/插值/修正误差来通过。交付5组PNG/FIG，审计重导图和失败候选另行保留。

![三个Seed的统一口径误差图](../results/comparisons/piston_2d_pm/r2_saved_closeout_20261117_01/three_seed_summary/figures/F05_three_seed_errors.png)

- seed260002：[PS/FullFK幅相误差](../results/comparisons/piston_2d_pm/r2_saved_closeout_20261117_01/seed260002/work/figures_final/F02_pointwise_errors.png)；[FullFK独立纵轴](../results/comparisons/piston_2d_pm/r2_saved_closeout_20261117_01/seed260002/work/figures_final/F04_FullFK_pointwise_errors.png)。
- seed260003：[PS/FullFK幅相误差](../results/comparisons/piston_2d_pm/r2_saved_closeout_20261117_01/seed260003/work/figures_final/F02_pointwise_errors.png)；[FullFK独立纵轴](../results/comparisons/piston_2d_pm/r2_saved_closeout_20261117_01/seed260003/work/figures_final/F04_FullFK_pointwise_errors.png)。

5张最终PNG均已实际目视检查，坐标、单位、图例、峰值呈现与保存数据相符；对应FIG可编辑。未生成未请求的PDF。

## 跨Seed空间分布

PS的Q4能量占比分别41.9369%、68.8164%、8.5851%，seed260003的Q1占71.2867%；PS线内复残差峰分别在x=21.09375、33.7890625、−39.6484375 m。FullFK线内峰在−18.5546875、−26.5625、7.03125 m，其能量最集中区依次为Q2、Q1、Q3（约33.64%、36.74%、60.43%）。因此模型排序保持，但误差热点与能量分布随Seed明显变化，不能固定热点位置或用411点当411个独立随机重复。

[逐点2466行](../results/comparisons/piston_2d_pm/r2_saved_closeout_20261117_01/three_seed_summary/point_errors.csv)、[互斥分段能量24行](../results/comparisons/piston_2d_pm/r2_saved_closeout_20261117_01/three_seed_summary/disjoint_energy_regions.csv)、[原共享端点区间范数24行](../results/comparisons/piston_2d_pm/r2_saved_closeout_20261117_01/three_seed_summary/shared_endpoint_region_norms.csv)、[三Seed排序表](../results/comparisons/piston_2d_pm/r2_saved_closeout_20261117_01/three_seed_summary/model_ordering.csv)均已保存。能量分区不重叠，范数/收敛区间保留原共享端点；各Seed使用自己的BIE尺度，不相干平均不同海面的复场。

## 封存、失败保留与执行记录

| Seed | 本轮正式状态 | 封存绑定文件数 | committed.json SHA-256 |
|---|---|---:|---|
| 260002 | FORMAL_LOCAL_SAVED_CLOSEOUT_PASS | 125 | `599940b3388939b582da7a6b803a5155d07076f23458b2211c1d722cd5b2e9cc` |
| 260003 | FORMAL_LOCAL_SAVED_CLOSEOUT_PASS | 84 | `74e31f094d0ac1d8a4f082f37bd2e92d489f94a1124b53d7202f6f1be8919009` |

每Seed的`closeout_receipt.json`、`numeric_receipt.json`、`figures_receipt.json`、`code_versions.json`与`committed.json`独立保存于`results/comparisons/piston_2d_pm/r2_saved_closeout_20261117_01/seed<seed>/`。签名绑定原manifest、所有输入/场/审计/图件/代码及历史失败证据。任务总封存在根`committed.json`，摘要表的资格登记引用各Seed封存SHA。原seed260002仍是FAILED_STOPPED/退出1，新本地完成不改变原身份。

成功数值审计耗时：seed260002 59.700834 s，seed260003 60.401083 s；最终图件审计分别28.510568/29.525169 s；三Seed汇总22.617891 s。这些是各阶段函数计时，不是全部开发/失败尝试的总墙钟。新入口Code Analyzer均0消息。本地审计在阶段边界保留30%启动可用RAM门，不声称有远端作业式1秒采样序列。

原远端像素失败、此前诊断及本轮候选失败全部保留。本轮新增审计最初直接重算RMS产生7.10543e-15的尺度舍入差，导致归一化曲线最多1.66533e-16差异；改为先按原1e-12尺度门独立核验，再使用原保存尺度做逐点精确一致性，未修改任何数组或容差。固定坐标轴的首候选仍因FIG重开缩小画布造成904148通道值差，查清后加入原画布恢复；内部图例候选也保留，最终图例移到轴外。

汇总包装器的空日志SHA及局部函数作用域错误在生成正式汇总前被拦截，已修正为复用原获批v2空文件SHA原语和独立局部函数；两Seed封存未受影响。候选代码/失败回执保存在`failed_attempts/`、`preliminary_layouts/`和`summary_debug/`。MCP捕获输出时部分diary可为空，原样封存且不当作PASS依据；完整数值/图件JSON、MCP工具回执、汇总stdout和SHA检查提供实际完成证据。

## 变更、限制与停止

新增`closeout_r2_saved_seed`、`audit_r2_final_figures`、`render_r2_stable_saved`、`summarize_r2_qualified_closeout`及准备/封存工具；复用冻结R1独立数值审计，不改生产实现、公共API、原运行包或数值门。对应脚本索引、任务表、PROJECT_CONTEXT及进度源同步更新；无生产/通信行为变化，因此不运行通信回归、无需改vertical_comm_guide。

本次只关闭固定0.10 m理想二维条带、4 kHz、三个既定PM实现和411真实接收点的有限复核与保存交付。不能外推其他频率、海况、真实3-D有效孔径、CIR或通信性能；R3模型决策仍待单独明确范围。本轮完成后停止，不SSH、不继续求解、不自动进入R3/宽带通信。
