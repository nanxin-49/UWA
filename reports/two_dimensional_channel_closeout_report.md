# 二维收尾 S0–S1：保存数据审查报告

日期：2026-10-20（Asia/Shanghai）。阶段状态：**`S0_S1_COMPLETE_WITH_DECLARED_GAPS`**。本报告只完成 S0–S1；未执行 S2–S6。`G1_target_PM=OPEN`、`WP2_formal=false`。

> **2026-10-21解释更正（原结果保留）：**[高度相位符号审计](pe_bie_height_phase_sign_audit.md)已定位固定比较调用链的几何高度/PE屏相位错配及旧单边共轭问题。本报告的公式、复数、图像复现通过及S0–S1完成状态保持有效；下文8.0908%、相位贡献和点阵局部值均是旧分支算术结果，不能继续无条件解释为同一物理边界和表示下的纯PE模型误差。原生PE/BIE的差别不能只归因于后处理表示。旧PE粗糙场也不能直接复用于修正基线。当前已有[修正图册v2](pe_bellhop_bie_corrected_visualization_report.md)，下一步按[主计划S2a补齐约定控制与基线资格核验，再S2b点阵补齐](../two_dimensional_channel_closeout_plan.md)执行；旧缺口表作为历史记录保留。

本轮固定 4 kHz、c=1500 m/s、Tx=(0,100) m、Gaussian σ=0.3 m；接收点为 x=[−0.1953125,0,0.1953125] m、z=[2.875,3,3.125] m 的 P11–P33。PE 为 Model-0，BIE 用保存的数值参考，Bellhop 为 X 线源、internal-wall 分支。所有处理只读已有复数场；新增 PE/角谱/BIE/FullFK/Bellhop 物理场计算为 **0**，没有从保存密度评价新点。没有改变生产求解器或公共接口，也没有提交、推送或读取 `cash/`。

## 结论

**旧接收线公式未发现错误。** 对三档正弦、三组方法对照，使用原 411 点 M99、原入射能量权重和原 G 表示，从保存完整复数重新计算的九组复 L2、相位 RMS 与幅度比 RMS，均在 `atol=1e-12, rtol=1e-10` 内重建。独立从原 CSV 重算的最大绝对差为 `9.06e-16`。强正弦 PE–BIE 的旧值 `E=0.0809081077437715`、相位 RMS `0.0803363527071867 rad`、幅度比 RMS `0.0869963865863484 dB` 均重现；没有用报告的舍入表值作高精度断言。

**旧评价对象是整条 z=3 m 接收线上的粗糙/平面反射比 G，不能直接代表单点或完整 3×3 点阵，更不是原始绝对压力、总信道或 BER。** 当前中心 P22 的强正弦 G 逐点相对复差为 `0.0779531`，相位差 `−0.0779746 rad`，幅度比 `−0.0003864 dB`。同深度三点 P21/P22/P23 的 `E_subset=0.0780765`，最大逐点相对差 `0.0781852` 位于 P21；与整线 `0.0809081` 量级接近，但这只是现有三点的观察，不能填补另外六点。整线强正弦逐点相对差的等点数均值/中位数/最大值分别为 `0.0979288/0.0777888/0.373721`，最大处 x=40.0390625 m；它们与入射能量加权复 L2 是不同定义。

强正弦旧线复差平方的幅度/相位代数贡献为 `1.5324%/98.4676%`；三点局部 G 为 `0.0010%/99.9990%`。三个低 K 正弦在线上都以相位项为主。较大的复误差不等于同等幅值误差；图 F03 的大尺度、幅度主导色图可显得相似，却不展示这种相位残差。贡献比例是恒等式的平方项份额，不是声能损耗比例。

## 来源、可比性与数值结果

输入哈希、存在性、源/表示/单位语义在本次[来源表](../results/comparisons/two_dimensional_channel_closeout/20261020_2330_s0s1/source_manifest.csv)；[协议](../results/comparisons/two_dimensional_channel_closeout/20261020_2330_s0s1/protocol.json)、[9点坐标](../results/comparisons/two_dimensional_channel_closeout/20261020_2330_s0s1/receiver_grid.csv)、[案例表](../results/comparisons/two_dimensional_channel_closeout/20261020_2330_s0s1/case_manifest.csv)、[可用性长表](../results/comparisons/two_dimensional_channel_closeout/20261020_2330_s0s1/data_availability.csv)与[成本表](../results/comparisons/two_dimensional_channel_closeout/20261020_2330_s0s1/cost_plan.json)构成 S0 冻结。坐标与保存轴逐点精确匹配，误差小于 1e-9 m；没有插值或换点。E01 的 G 定义为 rough/flat，BIE/Bellhop 比值沿用已经保存的共轭表示；E02 的 PE 恢复解析载波，BIE/Bellhop 原生场保持自身表示。Bellhop 原生源幅度未与 Gaussian 峰值统一，因此不作 Bellhop 原始场绝对误差排名。

| PE–BIE G 案例 | 411点原权重复 L2 | 原权重相位 RMS/rad | 幅度比 RMS/dB | P21/P22/P23 的 E_subset | 三点均值/中位数/最大逐点相对差 |
|---|---:|---:|---:|---:|---:|
| 弱 A=.01 m | 0.00315757147 | 0.00311763566 | 0.00434837965 | 0.000195790889 | 0.000195788 / 0.000195580 / 0.000197248 (P23) |
| 中 A=.05 m | 0.01608996037 | 0.01589449468 | 0.02174179797 | 0.004884690803 | 0.004884669 / 0.004892021 / 0.004897646 (P23) |
| 强 A=.20 m | 0.08090810774 | 0.08033635271 | 0.08699638659 | 0.07807649395 | 0.078076426 / 0.078090945 / 0.078185189 (P21) |

另外两组 BH–BIE 与 PE–BH 的全部三档重建、逐点和幅相分解保留在[整线汇总](../results/comparisons/two_dimensional_channel_closeout/20261020_2330_s0s1/line_metrics.csv)、[整线逐点](../results/comparisons/two_dimensional_channel_closeout/20261020_2330_s0s1/line_point_metrics.csv)、[点阵逐点](../results/comparisons/two_dimensional_channel_closeout/20261020_2330_s0s1/point_metrics.csv)和[点阵汇总](../results/comparisons/two_dimensional_channel_closeout/20261020_2330_s0s1/point_summary.csv)。三方点阵指标使用共同的 P21/P22/P23 支持，不跨不同点集排名。每一行含实际复数、幅度、绝对/相对差、有符号 dB、wrapped 相位、代数分解、点号和有效性；汇总按指标列出有效数及最大点号。`E_array` 仅留给九个可比复数都齐全的情况；当前 G 只有 `E_subset`。

强正弦 E02 真实二维网格的 9 个点，PE/BIE/Bellhop 原生反射复数都已提取到[原生场样本表](../results/comparisons/two_dimensional_channel_closeout/20261020_2330_s0s1/native_field_samples.csv)。PE/BIE 原生**幅度**具有同一 Gaussian 源标度；九点幅度相对差均值/中位数/最大值为 `0.000215237/0.000254523/0.000375290`，最大在 P31；幅度 L2 为 `0.000249258`。这不是九点可比复相位或九点 G 误差。z=3 m 的 P21/P22/P23 另有保存平面参考，按既有 E01 共轭/载波约定建立比较表示的反射复场，`E_subset=0.0780764884`，与 G 的三点结果闭合。对 E02 原生 PE 与原生 BIE 粗糙场直接作复相位差会混入表示约定；未将该值列为模型误差。其他六点虽有粗糙原始场，缺对应平面分母，G 仍标缺失。

零参考时相对差、幅度比和整体归一误差无定义；零候选时相位无定义且 dB 为 `−Inf`；双方为零时仅绝对差定义。近零但非零值保留原始比值。本轮没有同源、同量、同点的参考不确定性界，`reference_uncertainty=unknown`，没有借用别处门限。9 点是同一海面的诊断位置，不是九个独立海面或置信区间。

## 公式、图像与独立重绘

逐点计算 `d=|A−B|`、`e=d/|B|`、`ΔA=20log10(|A|/|B|)`、`Δφ=arg(A conj(B))`；整体复差为 `sqrt(sum(w|A−B|²)/sum(w|B|²))`。检查的精确恒等式为 `|A−B|²=(|A|−|B|)²+4|A||B|sin²(Δφ/2)`。同场、纯幅差、纯相差、零候选、零参考、双方零、近零参考、缺失样本与分解闭合九个控制全部通过。整线最大逐点分解余量不超过 `5.56e-17`。公式测试只用人工复数数组，未计算物理场。

[F01 几何](../results/comparisons/two_dimensional_channel_closeout/20261020_2330_s0s1/figures/F01.png)、[F03 原始二维场](../results/comparisons/two_dimensional_channel_closeout/20261020_2330_s0s1/figures/F03.png)、[F04 幅相接收线](../results/comparisons/two_dimensional_channel_closeout/20261020_2330_s0s1/figures/F04.png)、[F05 残差线](../results/comparisons/two_dimensional_channel_closeout/20261020_2330_s0s1/figures/F05.png)、[F05 点阵缺失图](../results/comparisons/two_dimensional_channel_closeout/20261020_2330_s0s1/figures/F05_grid.png)、[F06 五点对照](../results/comparisons/two_dimensional_channel_closeout/20261020_2330_s0s1/figures/F06.png)、[F07 汇总](../results/comparisons/two_dimensional_channel_closeout/20261020_2330_s0s1/figures/F07.png)均实际生成并逐张检查；每张另有 PDF。完整图注与检查点在[图注](../results/comparisons/two_dimensional_channel_closeout/20261020_2330_s0s1/figures/captions.md)。F02 涉及尚缺的匹配入射/平面多深度控制，本轮未补算，故未冒充生成；F08–F10 属后续阶段。

两次只读渲染均以[保存的 plot_data.mat](../results/comparisons/two_dimensional_channel_closeout/20261020_2330_s0s1/plot_data/plot_data.mat)为唯一绘图输入，写入独立目录；七张 PNG 逐像素相同。[独立审计](../results/comparisons/two_dimensional_channel_closeout/20261020_2330_s0s1/audit/independent_audit.json)还从原 E01/E02 CSV 复核九组整线数值、3699 条整线逐点、99 条点阵记录、27 条原生场记录。CSV 与 MAT 整线指标最大差 `9.06e-16`，原始来源 CSV 与点阵 G/原生复数的最大差为 `0`。[公式控制](../results/comparisons/two_dimensional_channel_closeout/20261020_2330_s0s1/audit/formula_controls.csv) 9/9 通过。五个新增 MATLAB 文件的 Code Analyzer 均为零消息；本轮没有生产物理或公共接口修改，因此未运行生产标量/宽带回归，`vertical_comm_guide.md`亦无需改写。

## 缺口及 S2 最小补算建议（未执行）

**本节为2026-10-20旧分支缺口清单，后续执行以主计划2026-10-21修订为准。**旧强正弦PE场不能直接复用；当前修正图册v2已有新版强正弦二维场和三条接收线，应先核验复用。若其通过S2a资格核验，S2b名义缺口为18个平面分母、24个BIE/Bellhop弱中粗糙值及12个修正PE弱中离深度值；按其他合格现值继续去重。这些是点值覆盖数而非作业数。禁止沿用下方旧建议中的PE正号反射起始场及针对旧PE匹配的BIE/Bellhop单边比值共轭。

[逐条缺失 G 点](../results/comparisons/two_dimensional_channel_closeout/20261020_2330_s0s1/missing_g_samples.csv)记录 3 案例 ×3 方法 ×6 离开 z=3 m 的点，共 54 个缺失 G 值；[平面分母缺失表](../results/comparisons/two_dimensional_channel_closeout/20261020_2330_s0s1/missing_flat_denominators.csv)记录每方法六个离深度点，共 18 个独特复分母。强正弦三方法的 18 个离深度粗糙复场已在 E02 保存；弱/中正弦离深度粗糙场尚缺 36 个复数值。这些计数是**所需点值下界**，不是作业/求解次数。S2 授权后最小值清单见[s2_minimum_recompute.csv](../results/comparisons/two_dimensional_channel_closeout/20261020_2330_s0s1/s2_minimum_recompute.csv)：先补共享平面 18 值、弱中粗糙 36 值，复用强粗糙现值；按 PE 载波、BIE/Bellhop 比值共轭、X 源和相位/幅度单位逐点核对，再以匹配点的空间、求积、波束等敏感性建立数值参考不确定性。BIE 是否能仅评价保存密度、Bellhop 是否需新作业，应在 S2 前做可追溯预算；本轮均未触发。完整九点 G 与原始复场相位评价、Bellhop 原始场定标在此之前保持未完成。

E04 是 σ=1 m 入射/反射诊断，不直接替代本轮 σ=.3 m 参考；E05/E06 只给固定 PM 及条件性 FullFK 信息，未新增 PM；E07 的单接收点计划不外推到点阵；E08 为 Bellhop R 源，不能合并到 X 样本；E09 的 CIR/通信控制不构成本轮宽带信道证据。因此没有 PM 认证、完整九点物理通过、H(f)/CIR、EVM/BER 或三维结论。

## 复算入口与变更

只在已有 MATLAB MCP 会话中执行，入口可从任意 MATLAB 当前目录运行，使用绝对项目根且拒绝覆盖已有目录：

```matlab
addpath('E:\MISC\CARPE3D_matlab\Explain');
setup_vertical_project;
test_two_dimensional_closeout_metrics();
r = run_two_dimensional_closeout_s0_s1('your_new_unique_run_id');
audit_two_dimensional_closeout_s0_s1(r.run_dir);
render_two_dimensional_closeout_s0_s1(fullfile(r.run_dir,'plot_data','plot_data.mat'), ...
    fullfile(r.run_dir,'audit','another_new_redraw'));
```

本次实际入口、MATLAB 版本/cwd、执行时间、先前失败及修复见[执行记录](../results/comparisons/two_dimensional_channel_closeout/20261020_2330_s0s1/logs/execution.md)。新增代码限于 `scripts/comparisons/run_two_dimensional_closeout_s0_s1.m`、`scripts/reporting/render_two_dimensional_closeout_s0_s1.m`、`scripts/validation/support/two_dimensional_closeout_metrics.m`、`scripts/validation/support/audit_two_dimensional_closeout_s0_s1.m` 与 `scripts/validation/test_two_dimensional_closeout_metrics.m`。除报告和导航外没有修改物理、通信方法或公共输出。S0–S1 成果完成于已保存且可比的数据范围；缺口明确保留给 S2，未自动进入下一阶段。
