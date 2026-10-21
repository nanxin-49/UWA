# 垂直水声 PE/WAPE 信道与通信项目

本仓库用于研究海底发射机到近海面接收机的垂直水声传播。公共模型以 PE/WAPE 计算直达与海面反射频响，并提供粗糙海面、joint-kstat、接收端精确离散伴随、二阶统计和 MPSK 通信验证工具。

## 快速开始

在 MATLAB 中先初始化一次路径：

```matlab
cd('<repository-root>')
setup_vertical_project;
```

最短调用：

```matlab
paramsV = struct('f0',6000,'enable_wideband',false, ...
    'enable_surface_reflection',true,'show_figures',false);
output = vertical_channel_model(paramsV);
```

也可以运行兼容入口：

```matlab
main_vertical              % 信道演示
explain_main_vertical      % main_vertical 的兼容别名
comm_main_vertical_psk     % MPSK 通信演示
```

这些名称保持不变；实际演示脚本位于 `examples/`，核心实现位于 `src/`。

BLACK 实测信道的独立基础分析入口：`setup_vertical_project; r = run_black_channel_basics();`。
仅下载 BLACK 及其配套噪声，调用固定版本官方代码解包/短段回放，输出多径与基本多普勒图表；
`run_black_channel_basics('render', r.run_dir)` 从保存结果重绘。
见 [BLACK 场景与特征通俗汇报报告](reports/black_channel_basics_report.md)和[脚本索引](scripts/README.md)。
数据、源码快照和每次运行结果位于 `results/validation/black_channel_basics/`。
可直接复制进笔记库的 [BLACK Obsidian 阅读包](results/visualization/black_channel_basics/20261019_obsidian/BLACK深海垂直信道汇报.md)包含报告与三张原图；图片使用相对 wiki 引用，不含数值数据副本。

## 公共入口与目录

2026-10-21 [PE/BIE 高度相位符号审计](reports/pe_bie_height_phase_sign_audit.md)
解释了 TR01 幅值近似镜像：旧固定比较的高度符号与单边共轭未一致衔接；旧图册指标需连同此更正阅读。
[修正后的三方接收线与二维场图册](reports/pe_bellhop_bie_corrected_visualization_report.md)
使用同一保存入射场重算负号PE相位屏并与BIE/Bellhop原生复场比较，旧图册保留作历史参照。

两套图册的[通俗展示汇总](reports/pe_bellhop_bie_visualization_summary.md)与
[配套数据说明](reports/pe_bellhop_bie_visualization_data_guide.md)：
按“声场外形—幅相差异—定量结论”解读全部11张图，区分原始二维场与接收线粗糙/平面比。

Bellhop–BIE 与 PE–Bellhop–BIE 固定二维对比图册：
[图册、数据和复现命令](reports/pe_bellhop_bie_comparison_visualization_report.md)。
入口 `prepare_pe_bellhop_bie_visualization_data` 复用PE/BIE并补算Bellhop，
`generate_pe_bellhop_bie_comparison_visuals(run_dir)` 从完整复数数据独立重绘7张PNG/PDF。
另见[原始x–z声场扩展](reports/pe_bellhop_bie_xz_visualization_report.md)：
强海面411×45真实网格的幅值、相位、幅相合成与原始剖面；保留原生源标度，不用比值替代原场。

FK 优化汇报图集（仅重绘已有验证数据）：
[图集与展示顺序](reports/fk_optimization_ppt_visualization_report.md)。执行
`setup_vertical_project; generate_fk_optimization_ppt_visuals();` 输出 11 张 PNG/PDF，
包含 Incident/Reflected 场图，并分别展示算法保真和物理适用域。

| 路径 | 用途 |
|---|---|
| `vertical_channel_model.m` | 稳定公共 API 包装；初始化路径并调用信道实现 |
| `main_vertical.m` / `explain_main_vertical.m` | 信道演示兼容入口 |
| `comm_main_vertical_psk.m` | 通信演示兼容入口 |
| `setup_vertical_project.m` | 幂等加入 `src/`、示例和必要脚本目录，不改变工作目录 |
| `src/` | 可复用生产实现，按 channel/propagation/surface/receiver/statistics/communication/bubble 划分 |
| `examples/` | 演示脚本实际内容 |
| `scripts/` | 验证、比较、实验和报告入口 |
| `reports/` | 一次性验证报告和当前状态报告 |
| `results/` | MAT、CSV、图片、视频及检查点 |
| `docs/history/` | 非权威历史规格与旧说明，仅供追溯 |
| `cash/` | 已退役数据的本地隔离区；不参与正常读取/检索/运行，Git 忽略全部内容 |

模块和依赖方向见 `src/README.md`，脚本命令见 `scripts/README.md`。

## 当前能力与边界

| 项目 | 当前状态 |
|---|---|
| uniform CPU-double PE/WAPE，直达与海面反射 | 可用；公共默认相位参考为 `direct_dsp` |
| Kirchhoff spatial/k-domain 与 joint-kstat | 可用；公共默认 surface model 未改变 |
| cached forward、精确离散伴随、PM 周期 FFT `C/P` | 固定环境、单最近网格点接收机原型已验证 |
| U=5/U=8 条件统计生成器 | 已有 F=64 验证证据；适用于已登记固定条件 |
| bubble、SSA 诊断分支 | 已实现，但适用范围和验证强度不同，见技术指南 |
| layered、Doppler、GPU、多接收机、插值接收的伴随统计链 | 尚未纳入当前原型 |
| 验证专用 Full/Fast FK 跨频一致性 | 固定 PM、4–8 kHz/F65、统一3/12 rank PASS；包含原生物理到IFFT桥接，未生产集成 |

当前公共窗口参数仍保持兼容默认值，并未自动切换到近期收敛验证的推荐窗口。正式研究运行前请阅读技术指南中的窗口建议和证据边界。

## 文档导航

| 文档 | 主要读者 | 内容 |
|---|---|---|
| `README.md` | 初次进入仓库的人 | 快速入口、目录和导航 |
| `vertical_comm_guide.md` | 研究与开发人员 | 当前物理、接口、结果语义和限制的权威说明 |
| `PROJECT_CONTEXT.md` | 后续模型与维护任务 | 顶部当前索引和完整追加式项目日志 |
| `scripts/README.md` | 运行验证的人 | 脚本入口、环境变量、成本和输出位置 |
| `reports/current_project_status_20260820.md` | 项目审阅者 | 已完成、部分完成、开放问题和推荐配置差异 |
| `reports/spectral_fast_kirchhoff_multifrequency_consistency_report.md` | FK扩频研究人员 | 固定海面、统一rank、BIE锚点、绝对相位与诊断CIR |
| [Deep-sea research handoff](reports/deep_sea_channel_research_background_for_new_ai.md) | 新的研究评估者 | 英文自包含需求、历史证据、单点计时转向和待论证计划 |
| [External-reply reassessment](reports/deep_sea_channel_external_reply_reassessment.md) | 研究路线讨论者 | 英文外部建议核查、公开深海部署案例与建议测试范围；未改变生产模型 |
| [Official deep-sea research plan](reports/deep_sea_channel_research_plan.md) | 后续研究执行者 | 唯一正式研究计划；2026-10-15决策B：固定dual仅4/6/8kHz离散交付，停止参考连续修补，G1 OPEN/WP2=false；备份和v3仅版本记录 |
| [二维点阵对比收尾计划](two_dimensional_channel_closeout_plan.md) | Codex与二维研究执行者 | 2026-10-20专项计划；S0–S1保存数据首轮审查已完成并登记缺口，S2–S6未执行；不替代深海总计划 |
| [二维计划审查与首轮Goal提示词](reports/two_dimensional_channel_closeout_plan_review.md) | 首轮执行者 | 只启动S0–S1的完成契约；明确保存数据后处理、缺失点/平面分母、部分集统计及环境阻塞边界 |
| [二维 S0–S1 保存数据审查报告](reports/two_dimensional_channel_closeout_report.md) | 研究审阅者与下一阶段执行者 | 三档旧接收线重建、单点/部分点阵与原始场审查、图组与独立重绘证据、S2 最小补算清单 |
| [固定dual物理决策](reports/deepsea_dual_broadband_decision_2026-10-15.md) / [宽带成本](reports/deepsea_dual_broadband_cost_2026-10-15.md) | 研究评估者 | 反射幅度约降6.43%，总场干涉变化清楚；稀疏候选未获新独立频点认证，5/7kHz被资源门拒绝，未进行宽带扫描 |
| `reports/` | 验证审阅者 | 单次实验配置、结果和结论 |

历史文件原文保存在 `docs/history/`，不得作为当前 API、默认值或验证状态的唯一依据。

## 结果存放规则

- 验证结果：`results/validation/<case>/`
- 可视化结果：`results/visualization/<case>/`
- 一次性报告：`reports/`
- 大型正式结果不放回根目录，也不由目录整理任务自动删除或重算。
- `cash/` 仅在用户明确要求归档操作时访问；禁止日常自动读取或恢复。用途及约束见 `AGENTS.md`（`AGENT.md` 为说明入口）。

公共能力基线日期为 2026-08-20；最新研究验证索引见 `PROJECT_CONTEXT.md`（FK跨频验证更新于2026-10-08）。目录迁移仅改变文件位置和初始化方式，没有改变 PE marching、相位公式、海面统计公式或公共默认配置。

Bellhop 是外部验证依赖，不随仓库分发。运行 Bellhop 验证时可设置
`BELLHOP_EXE` 指向可执行文件，或设置 `BELLHOP_TOOLBOX_ROOT` 指向外部
Acoustics Toolbox 根目录；脚本会检查其 `windows-bin-20201102/bellhop.exe`。
缺少外部依赖时，验证入口应停在依赖检查，不得伪造运行成功。

## 深海图12(b)原生Bellhop复现尝试（2026-10-10）

独立入口 `setup_vertical_project; validate_bellhop_fig12b('full');` 建立5000 m
Munk、Tx4950 m、10 kHz的34个有符号距离声场/声线案例，结果在
`results/validation/bellhop_fig12b/`。原生点声源把r=0整列置零；全角近轴
幅值与声束数不稳定，本轮仅完成部分定性复现，不能作为已验证通信信道。
见[中文报告及补充假设](reports/bellhop_fig12b_reproduction_report.md)。本案例独立
于下述100 m安装验证，没有改变生产参数或公共API。

随后按用户限定完成三种追加方案（高斯点源、Cerveny点源、高斯线源），
每种含声束/步长及均匀介质控制，共15次运行，仍未可信复现，现已停止搜索。
入口`validate_bellhop_fig12b_three_attempts('full')`；见
[三次尝试与停止报告](reports/bellhop_fig12b_three_attempts_stop_report.md)。

## PE--Bellhop 展开坐标验证

正常坐标粗糙海面实验现按用户决定采用 **100 m 水深、物理源离底1 m（Tx99 m）**，
Rx仍为3 m。当前独立绘图入口为 `generate_bellhop_100m_source1m_visuals_vertical`，
先做4 kHz/89.5°代表案例，显示全部有限范围内的接收贡献和完整发射扇区。
见 [执行计划](reports/bellhop_100m_rebaseline_plan.md)。旧500 m结果已隔离；
[边界源实测报告](reports/bellhop_100m_boundary_source_report.md)保留为历史诊断，
不再继续微小离底极限扫描。没有修改PE生产默认收发位置。
两类图与范围说明见 [离底1 m展示报告](reports/bellhop_100m_source1m_visuals_report.md)。
若只查看与当前PE拓扑相同的路径，运行
`generate_bellhop_pe_equivalent_ray_paths_vertical`；该入口仅绘制一条直达和
一条“一次海面、零海底”Bellhop轨迹，不把含底或高阶多次反射纳入图片。
以下镜像展开证据仍在其原限定环境内保留。

当前新增的验证入口
`scripts/validation/validate_pe_bellhop_unfolded_flat_gaussian_vertical.m`
采用镜像展开：直达路径对应 Bellhop 97 m，平面压力释放反射路径对应
103 m 的镜像接收端并乘以 `-1`。该入口使用现有生产 Gaussian 源的角谱
指向性生成 `.sbp`，不修改 PE 核心或生产默认窗口；正式输出位于
`results/validation/pe_bellhop_unfolded_flat_gaussian/`。
