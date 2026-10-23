# 二维收尾主报告：S0–S1、S2a 与 S2b 点阵及参考限制

日期：2026-10-20；S2a/S2b 更新：2026-10-23（Asia/Shanghai）。当前 S2b 状态：**`S2B_COMPLETE_REFERENCE_QUALIFIED`**，限固定离散源、海面、411点接收线/原九点及已声明的经验数值分辨率，不等于模型物理PASS。P0、108/108真实点值、当前版本BIE再资格、FullFK完整外场/求积控制、图表与独立CSV/重绘审计均已完成。最新证据在同一run_id下的[formal_reference_v2](../results/comparisons/two_dimensional_channel_closeout/20261023_s2b_1305/formal_reference_v2/README.md)，详见文末“参考闭合补证”附录；前一“参考受限”阶段作为历史快照保留。历史 `S0_S1_COMPLETE_WITH_DECLARED_GAPS`、`S2A_CONVENTION_SOURCE_BASELINE_COMPLETE` 保留；S3–S6未执行。`G1_target_PM=OPEN`、`WP2_formal=false`。

**S2a权重解释补正：**以下S2a表及其原图保留历史入射能量加权口径；新空间主指标为同一M99的411点等权，另存本次独立结果包，不能将旧4.39622%改标签。强正弦PE–BIE原生反射场复L2的新空间值为6.29112453%，旧入射加权值为4.39621797%；都不是旧混合约定G的8.0908%。本次没有改写S2a历史MAT/CSV/图。

## S2a 实绩（2026-10-23）

本轮只执行[主计划 S2a](../two_dimensional_channel_closeout_plan.md)的实际入口、声源资格和平面基线审查。交付目录为 [20261023_s2a_final](../results/comparisons/two_dimensional_channel_closeout/20261023_s2a_final/)。主量是同源、原生 `exp(−iωt)` 表示下的反射复场 `u_flat/u_rough`，单位为**单位源响应**，不是 Pa 或绝对 TL；`G=u_rough/u_flat` 和 `δu=u_rough−u_flat` 另存为辅助诊断。固定 4 kHz、c=1500 m/s、Tx=(0,100) m、Rx z=3 m、Gaussian σ=.3 m，接收线使用原 M99 的 411 点和入射能量权重。没有拟合复系数，也没有以中心幅比定源强。

### 实际源、约定与资格

物理源是在 z=100 m 发射平面的二维 Gaussian `exp[−x²/(2σ²)]`，峰值为 1；`N=984`、周期宽度 `W=192.1875 m`、`x0=−96.09375 m`、`Δkx=2π/W`，离散复谱 `coeff=fft(source)` 并以 `1/N` 反变换。`kz=sqrt(k²−kx²)` 保留倏逝分量，入射表达为 `Σ coeff·exp[i kx(x−x0)]·exp[−i kz(z−100)]/N`；导数分别乘 `i kx` 和 `−i kz`。BIE 实际 `incident_fn` 使用此离散源；FullFK 用相同谱计算真实 `z=eta_z` 边界上的入射场与 `∂n u_inc`，入水法向 `n=(−eta_z′,1)/sqrt(1+eta_z′²)`。PE 保存初始数组与解析 Gaussian 逐点完全相同，保存海面入射包络恢复 `exp(+ik·100)` 载波后与角谱相对差 `4.56×10⁻¹³`。恢复 `exp[+ik(100+z)]` 的 PE 平面反射场与 BIE 保存原生平面场相对差 `4.70×10⁻¹³`，BIE 平面与独立镜像源角谱相对差 `3.53×10⁻¹⁶`。

[源资格表](../results/comparisons/two_dimensional_channel_closeout/20261023_s2a_final/data/source_qualification.csv)将 PE–BIE、FullFK–BIE、PE–FullFK 列为 `SOURCE_MATCHED`，此标签只证实上述固定源和表示的输入合同，**不表示反射模型精度通过**。PE 的物理高度 `eta_z` 在[验证适配入口](../scripts/validation/support/run_pe_physical_surface_model0_validation.m)中仅一次转换为旧验证 helper 的向上高程 `eta_up=−eta_z`，所以实际 Model-0 屏为 `−exp(−2ik eta_z)`；传播核沿用保存比较的 `exp[iz(kz−k)]`，统一恢复原生载波后才作比值。旧 `+2ik eta_z` 和单边共轭分支仍标历史诊断。Bellhop 的 `X` 角度指向性、原生 SHD 共轭转换可追踪，但源强、整体相位及二维源单位到 Gaussian 谱的物理映射仍缺，列为 `SOURCE_UNRESOLVED`；未以固定平面最小二乘 α、中心幅比或 0.27% 旧拟合残差升级资格。原生、物理映射和可选固定拟合严格分层，本轮**未应用拟合**。

[继承 E01–E09 哈希复核](../results/comparisons/two_dimensional_channel_closeout/20261023_s2a_final/data/e01_e09_source_recheck.csv)为 20/20 一致；[本次来源与代码哈希](../results/comparisons/two_dimensional_channel_closeout/20261023_s2a_final/data/source_manifest.csv)包括 S0–S1 审计、E10 符号审计、E11 修正 v2、E12 FullFK 及本次入口。E11 原包的保存数据/像素审计仍为通过。E12 保存弱/强海面 FullFK 入射值、法向导数与本次同源重算相对差至多 `3.43×10⁻¹⁶`，法向逐点相同。其 4097 点完整双项积分保留真实边界距离、`arc_length` 与 `local_radius_m=Inf`；`receiver_field` 为反射分量，不混入直达。平面和弱档用当前[完整积分入口](../scripts/validation/support/solve_full_kirchhoff_2d.m)重算，与 E12 保存线相对差分别为 `1.31×10⁻¹⁵` 和 `1.31×10⁻¹⁵`，故强档保留 E12 保存线；中档为本轮唯一新增的 FullFK 正弦线。平面、弱档、中档的两个积分项单独保存，和总场严格闭合；强档历史包未保存分项，分项标缺失，不伪称已复核。

[版本复用/失效表](../results/comparisons/two_dimensional_channel_closeout/20261023_s2a_final/data/version_invalidations.csv)另列旧正号PE与单边共轭仅作诊断、修正v2 PE线和E12 FullFK核验复用、Bellhop待源映射。

### 平面控制与接收线数值

正入射常数源在真实[适配—屏—传播—比值入口](../scripts/comparisons/run_two_dimensional_closeout_s2a.m)上运行 `eta_z=0,+0.01,−0.01 m`。独立参考由路径 `(100−h)+(3−h)=103−2h m` 推出 `G_h=exp(−2ikh)`，不是再写一次屏公式；三档最大复差 `1.75×10⁻¹³`，逐项结果和耗时在[平面控制表](../results/comparisons/two_dimensional_channel_closeout/20261023_s2a_final/data/plane_control.csv)。同一入口对保存 Gaussian 运行平面和三档正弦，平面场及修正 v2 的三个 G 线都闭合。Gaussian 整面平移并未被强制满足正入射常数比。

| 411 点原权重、原生反射场 | PE–BIE 复 L2 | FullFK–BIE 复 L2 | PE–BIE 相位 RMS/rad |
|---|---:|---:|---:|
| 平面 | 3.56×10⁻¹³ | 0.000185546 | 1.94×10⁻¹³ |
| 弱正弦 A=.01 m | 0.00217475 | 0.000185488 | 0.00311514 |
| 中正弦 A=.05 m | 0.0108807 | 0.000185735 | 0.0155852 |
| 强正弦 A=.20 m | 0.0439622 | 0.000189670 | 0.0629360 |

全部 PE–FullFK、幅度、dB、相位、两平方项及有效数见[12 行指标](../results/comparisons/two_dimensional_channel_closeout/20261023_s2a_final/data/line_metrics.csv)；[4,932 行原始复场](../results/comparisons/two_dimensional_channel_closeout/20261023_s2a_final/data/line_fields.csv)、[3,699 行 G/δu](../results/comparisons/two_dimensional_channel_closeout/20261023_s2a_final/data/G_delta_u.csv)、[9 行平面—增量—交叉项](../results/comparisons/two_dimensional_channel_closeout/20261023_s2a_final/data/increment_decomposition.csv)保留全精度。强正弦 PE–BIE 的 `0.0439622` 是原生 `u_rough` 加权相对复差，不能与旧分支 G 的 `0.0809081` 或修正 G 的 `0.0629113` 串成“误差改善”序列。恒等式 `e_rough=e_flat+e_delta` 的复内积交叉项独立保存，不能把两平方项当作唯一因果份额。

平面 PE–BIE 几乎闭合证明源/载波路径一致；FullFK–BIE 平面差 `1.86×10⁻⁴` 是本有限边界/求积/模型设置下的数值差，不能称严格误差上界。E12 的 2049→4097 点弱/强线相对变化约 `3.72×10⁻⁵`，中档本轮仅 4097 点；BIE 的匹配点参考不确定性和 FullFK 有限窗口敏感性仍未充分定界。因此上述差异是**已声明支持上的观察**，不是模型物理精度 PASS；`reference_uncertainty=unknown`。

### 图像、独立审计与 S2b 缺口

[F02 平面幅相与残差](../results/comparisons/two_dimensional_channel_closeout/20261023_s2a_final/figures/F02_flat_control.png)、[F04 三档原生反射线](../results/comparisons/two_dimensional_channel_closeout/20261023_s2a_final/figures/F04_rough_native_lines.png)、[F05 原生复差与辅助 G/δu](../results/comparisons/two_dimensional_channel_closeout/20261023_s2a_final/figures/F05_native_residuals.png)均实际渲染、逐张检查；同名 `.fig` 可编辑。图中 PE/BIE/FullFK 使用同一单位源标度，只有 z=3 m 的 411 个真实样本，未把线插成二维场。图注及限制见[图注](../results/comparisons/two_dimensional_channel_closeout/20261023_s2a_final/figures/captions.md)。[独立审计](../results/comparisons/two_dimensional_channel_closeout/20261023_s2a_final/audit/independent_audit.json)核对 19 个本次哈希、4,932 条复场、12 条指标、3,699 条 G/δu、9 条增量分解；中档 FullFK 的三个接收点从保存求积数组直接重新积分，最大差 `2.83×10⁻¹⁷`。三张 PNG 从 [plot_data.mat](../results/comparisons/two_dimensional_channel_closeout/20261023_s2a_final/data/plot_data.mat) 独立重绘逐像素一致，三张 `.fig` 均重开成功。六个改动/新增 MATLAB 文件的 Code Analyzer 为零消息；修改的旧 PE helper 默认 Gaussian 与显式相同 Gaussian 在 128 点回归中反射场和入射场最大差均为 0。

[S2b 108 格覆盖/缺失长表](../results/comparisons/two_dimensional_channel_closeout/20261023_s2a_final/data/s2b_coverage_gap.csv)按方法×平面/三档×九点×源版本保存：核心已存在 48 格，缺 **60 个复场点值**，其中 PE 18、BIE 18、FullFK 24；这是点值覆盖，不是 60 次求解。[分案成本表](../results/comparisons/two_dimensional_channel_closeout/20261023_s2a_final/data/s2b_cost_tasks.csv)列明每种方法/案例。强正弦 PE/BIE 已有九点真实原始场，其他离深度平面分母与弱中粗糙场依旧缺失；FullFK 本轮只有 z=3 m 线。Bellhop 源映射仍未闭合，作为条件扩展另有名义 18 个离深度原生点值缺口，不计入核心 108 格，也不阻塞合格子集。S2b 前先定 BIE 新密度评价/求解必要性及参考预算，按保存现值去重；不得由本轮直接开算。

本次新增声学计算的计时与内存门见[成本记录](../results/comparisons/two_dimensional_channel_closeout/20261023_s2a_final/data/cost.json)：代表 FullFK 作业约 `0.40 s`，弱档复核约 `0.37 s`，中档新线约 `0.36 s`，均低于单作业 900 s；PE 与平面控制合计亦远低于 7200 s。作业前 Windows 物理可用量约 `1.719 GB`、保守增量估计约 `0.409 GB`，保留 30% 的门通过；估计不是测量峰值。新增 BIE 求解/密度评价、Bellhop、完整 x–z 场、PM、宽带与 FastFK 调用均为零。生产 API 未改；`vertical_comm_guide.md` 因此不需改写。

复算命令（仅在已有 MATLAB MCP 会话中执行；运行名必须未使用，并在运行前从 Windows `\Memory\Available Bytes` 获取当时可用物理字节数）：

```matlab
setup_vertical_project;
r = run_two_dimensional_closeout_s2a('new_unique_id', available_physical_bytes);
a = audit_two_dimensional_closeout_s2a(r.run_dir);
```

`render_two_dimensional_closeout_s2a(fullfile(r.run_dir,'data','plot_data.mat'),new_figure_directory)` 是纯保存数据重绘入口。MATLAB 当前目录始终保持项目根目录；两次诊断运行 `20261023_s2a_a`、`20261023_s2a_b` 及其修复痕迹保留，最终证据仅以 `20261023_s2a_final` 为准。没有读取 `cash/`，没有提交或推送。

## S0–S1 历史审查（2026-10-20，原结论和数据保持）

> **2026-10-23计划修订时的历史说明（当时未执行新实验）：**[唯一主计划](../two_dimensional_channel_closeout_plan.md)改为PE/BIE/完整Kirchhoff核心三方法，Bellhop保留条件扩展；主量为同源反射复场，平面控制并列，G/δu辅助。当时把S2a声源资格与平面控制列为下一步，并指出旧54个缺值清单不涵盖新增FullFK需求。S2a现已按本报告上节执行；本报告下文仍是历史PE/BIE/Bellhop及G的S0–S1实绩，不改标签为新版三方法结果，也不重置其完成状态。

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

## 附录 S2b：实绩与参考限制（2026-10-23）

> 本节是同一Goal早期的`S2B_DATA_COMPLETE_REFERENCE_LIMITED`快照，包含当时真实阻塞、失败和成本；其“Goal未完成/无新BIE求解/下一单项”等时间状态现由文末参考闭合补证覆盖。保留历史数据与判断，不追溯改写为当时已通过。最终数据/图件为`formal_reference_v2`，早期有限积分版仍保留。

本次按用户授权创建并执行单个Goal，未设置token预算。唯一run_id为[20261023_s2b_1305](../results/comparisons/two_dimensional_channel_closeout/20261023_s2b_1305/)；run_id是目录标签，不表示精确起始时刻，实际作业时间见回执。全部声学运行经已有MATLAB R2025b Update 6 MCP会话，cwd保持本仓库，无shell MATLAB、无提交/推送、无cash访问。仅新增比较/审计/绘图工作流，未改生产求解器或公共API；不运行范围以外的通信回归，vertical_comm_guide.md的实现公式无变更。

### P0：相同原生场、双权重及零新场审计

只读S2a正式包的s2a_data.mat和line_fields.csv，使用相同4932条复数（4案例×3方法×411点），未调用PE传播、BIE密度求值、角谱场求值或FullFK积分。空间主权重q_i=1/411，support_id=M99_z3_411；历史w_i为保存d.weights的M99归一化值。24行表中分别保留spatial_uniform_411及legacy_incident_energy_M99各12行，附分母、n_total/n_valid、零幅值、近零、非有限dB、复L2、幅度L2、dB/相位RMS及幅相平方项。near_null阈值只是1e-12×本参考最大幅值的数值尺度标记，不能冒充参考误差界；reference_uncertainty仍有限定。

| 案例 | PE–BIE空间主复L2 | PE–BIE历史w复L2 | FullFK–BIE空间主复L2 | FullFK–BIE历史w复L2 |
|---|---:|---:|---:|---:|
| flat | 3.85502826e-13 | 3.56215625e-13 | 0.000253643479 | 0.000185545534 |
| weak .01 m | 0.003115550552 | 0.002174746032 | 0.000253489884 | 0.000185487806 |
| medium .05 m | 0.015586811401 | 0.010880734386 | 0.000253704239 | 0.000185734864 |
| strong .20 m | 0.062911245315 | 0.043962179741 | 0.000257000813 | 0.000189670369 |

权重改变了离散空间统计，不改变任何物理场或点位。新u空间复差数值与其他口径的G误差接近也不能互换名称；仍须按各自实际场量、分母和权重解释。全部三方法对见[data/line_metrics_both_weights.csv](../results/comparisons/two_dimensional_channel_closeout/20261023_s2b_1305/data/line_metrics_both_weights.csv)，权重向量见line_weights.csv；原始来源CSV副本、MAT和P0审计在同目录。CSV/MAT复场最大差5.02e-16，历史指标最大复现差4.34e-19，独立笛卡尔平方/余弦恒等式审计最大差1.12e-16。同场、纯幅差、纯相差、相消、零候选/零参考/双方零、近零及缺失共9组人工控制通过，新增物理场计算为0。

### P1：本地来源、四案例参考、资源和冻结工作清单

先核验本地正式包及S2a冻结的19项输入/代码哈希，再固定弱/中/强参考文件哈希，合计22项。弱/强密度在已有强场包source_data内；中档在已登记results/archive/cl-verification-worktree-unmanaged-20260930内，未读取cash，也未把GitHub报告当作大数据已同步的证据。各参考的4 kHz、1500 m/s、N984/W192.1875、σ=.3 m、Tx100 m、A=.01/.05/.20、C2 42/50 m、物理高度及入水法向逐项核验；保存z=3 m参考线与S2a复场一致。其密度可评价全部固定九点：z>=2.875 m低于粗糙海面，满足现有评价器限制。

- 平面：同一冻结离散复谱的镜像角谱解，无需平面BIE密度或边界重解；名义缺6个离深度点用解析平面参考评价。
- 弱档：保存空间8/10/12 ppw、窗口60/70/80 m和自求积64/96层级；6个缺点可直接用12 ppw密度求值。
- 中档：保存空间8/10/12 ppw、窗口60/70/80 m；6个缺点可用12 ppw密度求值。缺独立自求积层级。
- 强档：保存空间8/10/12/14 ppw及窗口60/70/80 m；九点现值保留，保存密度各层级在九点和接收线最差/最小参考幅值探针上重评，检查原值。缺独立自求积层级。

保存求解对象的schema/formulation/时间约定/法向/数值设置及当前冻结求解代码SHA可追踪，但原执行时源码SHA未保存，不能以当前代码哈希倒填历史执行版本。完整密度数字快照、节点、法向、窗口和配置见data/p1.mat；清单见bie_density_inventory.csv。保存对象中的字符编码保留原值；数值字段按HDF5对象引用读取，不执行原MAT内可能存在的函数句柄。

P1当时Windows物理可用RAM=1,351,811,072 bytes；保存密度评价增量初估64 MiB，实际各批统一保守上调至96 MiB。边界重解保守增量估计1,418,723,328 bytes，超过当时可用量70%，故该门失败，未启动新BIE边界解。这里是估计峰值，不是实测峰值；重解耗时未测，未假称可在900 s内完成。所有后续低内存批次前均重新采Windows Memory Available Bytes，逐批回执记录测量值、估计、时间和累计预算；不使用MATLAB页文件/虚拟容量替代物理RAM。后来的内存短暂回升不自动授权新的高内存求解作业。

S2a三个保存运行a/b/final实计合计9.2915103 s，均扣除，包括日志曾描述为停止的诊断a（其本地cost.json实际存在且记录物理成本，按文件扣记）。另为原未计时的预输出尝试/审计积分保守预留120 s，不将其冒充实测耗时。冻结作业清单见data/frozen_jobs.md：先强档BIE代表评价，再flat/bie1/bie2/pe1–4/fk1–4串行批次；每批900 s、累计7200 s、至少30%当时可用物理RAM保留。代表作业成功0.1216 s后才继续。完整冻结P1门见p1_gate.json。

### P2：108/108真实点值、复用及成本

名义60个缺值全部补齐：PE18，BIE18（其中12个粗糙密度场评价及6个平面镜像场评价），FullFK24；并非60次边界求解。保留48个原有值，未用复查值覆盖它们；重查最大绝对差7.83e-16。三方法×四案例各9点均有真实值，坐标严格为原P11–P33，无插值、共轭拟合、重选点或换源。全部108格源/单位/版本/原生与映射值/逐格来源见[data/coverage_after.csv](../results/comparisons/two_dimensional_channel_closeout/20261023_s2b_1305/data/coverage_after.csv)及fields.csv，MAT见s2b_data.mat；所有数据均为反射单位源响应，不是Pa、绝对TL或总信道。

PE从保存入射场重新构造负号物理高度屏，不复用旧正号反射起始场；FullFK采用完整两项、真实C2边界/法向、弧长求积、local_radius_m=Inf，不调用FastFK或局部截断。主FullFK保持4097/L60与S2a版本一致；更细求积或大窗口结果单列数值敏感性，未挑选较小误差替换主场。

| 九点等权、未经rough/flat归一化 | PE–BIE E_array | FullFK–BIE E_array | PE–FullFK E_array |
|---|---:|---:|---:|
| flat | 1.52177981e-13 | 0.000209724318 | 0.000209724931 |
| weak | 0.000001685861 | 0.000209666308 | 0.000209640275 |
| medium | 0.000012799862 | 0.000209654998 | 0.000207019495 |
| strong | 0.000080510858 | 0.000209688033 | 0.000205758686 |

九点仅覆盖x约±.1953 m、z=2.875–3.125 m的微小邻域，不能替代整条411点线。强档PE–BIE九点最大逐点相对差0.000145868625位于P33；九点E_array约0.00805%，与整线空间6.29112%的差距反映支持不同，不是方法精度突然提升。全部绝对/相对均值、中位数、最差点、幅相平方项、dB/相位和每项分母见array_metrics.csv、point_metrics.csv。九点和P0线均无缺失/零参考/尺度近零点；中强参考误差仍未知的部分不因此消失。G、δu及平面—增量—交叉项另存diagnostics.csv、increment.csv，仅作辅助。

12个成功新增声学批次总实测13.9342315 s；加三次S2a实测后为23.2257418 s。首次代表作业在密度评价后因cell receiver_ids触发struct展开，结果封装失败；修正一处类型后仅重试一次，声学输入未变。失败未取得核内计时，独立保守扣记30 s并保存bie3_failed记录。总预算扣记=23.2257418+120+30=173.2257418 s，剩余7026.7742582 s；不把150 s保守扣记称为实测。成功单批最大2.4232 s，远低于900 s。第二次非物理汇总因char表列不同长度报错，改string后重跑纯保存数据汇总；旧失败与修正均记于logs/execution.md。绘图/审计/I/O不计声学核，图形时间未单独仪器化，不报告伪造耗时。

### P3：匹配点数值检查及未闭合问题

BIE以相邻最终空间层级变化、窗口60→70与70→80的最大变化、可用自求积变化逐点记录绝对复场差，并保存每个层级真实复数到sensitivity_fields.csv（总2250行包含PE/FK层级）。已知分量之和的最大值约5.45e-8单位源响应，只是经验敏感性包络，不是严格误差界；中、强档自求积列为NaN，未借用弱档数值填补。平面参考是本离散源的解析镜像解，粗糙参考仍需执行版本和独立自求积控制补证，不能用线性残差小替代接收场误差认证。

PE步长.05→.025 m的411点相对变化约4.16e-13–4.24e-13；均匀介质平方根算子各步可交换，这主要是数值舍入一致性。横向984→1968使用同一冻结N984复谱零填充、相同W，并验证原节点源值一致；不能换成另一个离散Gaussian源来冒称同源收敛。匹配九点及最差/最小参考探针的变化从平面1.71e-13到强档3.33e-6。该检查不能推导未测全场的统一误差界。

FullFK L60的2049→4097相对变化约3.62e-5–3.86e-5，4097→8193降至8.73e-6–9.31e-6；但固定细网距120/8192 m、L60→L90（8193→12289点，物理C2支撑仍42/50）变化约2.52e-4–2.69e-4。原P1文字的L80不适合整数个相同细网距，故在第一作业前的frozen_jobs.md明确冻结L90。窗口变化与主FullFK–BIE九点差约2.10e-4、411点空间差约2.54e-4相当，甚至更大；**约0.019%的历史差异绝非已认证模型精度**，有限窗口影响尚未消除。这是实测的不确定性结果，没有通过缩点、调源或替换主场取得PASS。全部设置与绝对/相对变化见numerical.csv，各作业保留完整积分项、权重和几何。

因此本阶段为 **S2B_DATA_COMPLETE_REFERENCE_LIMITED**，不是S2B_COMPLETE_REFERENCE_QUALIFIED，也不是物理PASS。数据/图表工作完成，但本Goal未标完成；本轮按用户要求结束，不继续S3–S6、Bellhop源映射或通信，不自行设自动后续。下一次可先执行一个低内存单项：只针对平面、原九点及已冻结最差/最小参考探针，以同一离散源和解析镜像参考进一步审查FullFK有限窗口尾项；须先冻结窗口/源周期解释和成本门，不能盲目增窗越过源周期边界。BIE中强自求积控制和原执行源码证据作为独立待办，在资源/版本门通过前不强启重解。

### 图件、独立审计和复算入口

最终图件为[figures_v2](../results/comparisons/two_dimensional_channel_closeout/20261023_s2b_1305/figures_v2/)内F01–F07各FIG+PNG，无PDF。F01坐标/深度和九点几何；F02平面幅相与绝对/相位残差；F03仅使用已有强档PE/BIE真实411×45幅相场，FullFK连续二维场明确缺失；F04三档未归一化原生反射线；F05九点逐格残差（独立标注色标，不跨色标直接比较颜色）；F06五个代表点；F07空间线/历史线/九点等权指标并列。全部7张已实际检查；v2修复Tx标记贴边和F05标题拥挤，初版图与初次审计仍保留。包内captions.md注明相位wrap、幅值单位及支持；没有用九点或线插值出细声场。

[最终独立审计](../results/comparisons/two_dimensional_channel_closeout/20261023_s2b_1305/audit_v2/independent_audit.json)从CSV独立复算108复场、108逐点、12点阵、24双权重线指标、G/δu和交叉恒等式；复数用笛卡尔差及atan2计算，不调用主指标函数。MAT/CSV复场最大差5.01e-16，独立点阵指标最大差1.51e-14，线指标1.32e-15，G/δu闭合9.05e-16，交叉项7.90e-19；2250行敏感性复数一致，BIE敏感性从CSV重算最大差8.53e-16。22来源哈希通过；7张PNG只读plot_data.mat重绘后像素最大差0，7张FIG均重开。审计和重绘新增物理场计算为0。12个新增MATLAB文件最终Code Analyzer均零消息；初次仅有4条样式/性能提示，修正后复核记录见logs/code_analyzer_final.txt。最终cwd核对以分隔符归一化比较，避免把正反斜杠差异误判为目录变化。

仅保存数据重绘/复算（不得指定已有目录）：

```matlab
setup_vertical_project;
out = fullfile('E:/MISC/CARPE3D_matlab/Explain/results/comparisons/two_dimensional_channel_closeout','20261023_s2b_1305');
render_two_dimensional_closeout_s2b(fullfile(out,'data/plot_data.mat'),fullfile(out,'my_new_figures'));
a = audit_two_dimensional_closeout_s2b(out,'my_new_figures','my_new_audit');
```

新的物理补算必须重新执行P1门；不能把上述纯重绘命令当作新求解授权。工作流源文件与哈希、全部输入和S2a保存产物的完整性清单置于provenance/；data/qualification.json补充FullFK窗口限制，data/summary.json保留首次汇总快照。PROJECT_CONTEXT、README和scripts/README同步当前状态，主计划只追加执行状态说明，保留本轮开始时已有的计划修改。G1_target_PM=OPEN、WP2_formal=false及S2a历史状态/8.0908%、6.2911%、4.39622%原口径保持。

## 附录 S2b：参考闭合补证与最终退出（2026-10-23，同一Goal续行）

**最终阶段为`S2B_COMPLETE_REFERENCE_QUALIFIED`。**这是固定4 kHz、c=1500、σ=.3 m、Tx=(0,100)m、原C2平面/弱中强正弦、M99的411真实样本和P11–P33支持上的数值参考/交付完成，不等于PE或FullFK通过任意通用物理精度门。前节`REFERENCE_LIMITED`是同一Goal早期的真实快照，不删除也不改成当时已通过。全部新增补证仍在唯一run_id下：[正式数值包](../results/comparisons/two_dimensional_channel_closeout/20261023_s2b_1305/formal_reference_v2/README.md)。

### 1. 先重开资源门，未绕过原阻塞

续行先核对真实本地状态与现有MATLAB会话。原先1.35 GB物理可用时拒绝高内存BIE的决定仍有效；后续物理可用量改变后，逐作业重新采Windows Memory Available Bytes并独立通过门，没有拿页文件/虚拟内存替代。完整事前决策在reference_followup_01/bie_gate_extension.md；新回执保留actual available、conservative peak estimate、先验耗时估计、实耗、节点/求积阶、当前源数组及原BIE求解器SHA。原`solve_helmholtz_bie_halfplane_vertical.m`和生产求解器未修改。

| 新BIE作业 | 节点 / ppw / 自求积 | 实测秒 | 用途 |
|---|---|---:|---|
| medium pilot2 | 640 / 2 / 96 | 9.6831 | 只估成本，不作为可信参考或新增正式点 |
| medium pilot4 | 1280 / 4 / 96 | 30.2089 | 只估成本；另测非声学LU时间，不占声学核 |
| medium96 | 3840 / 12 / 96 | 242.9794 | 当前代码复现中档整密度 |
| medium64 | 3840 / 12 / 64 | 168.6432 | 中档同网格自求积对照 |
| strong96 | 4480 / 14 / 96 | 289.7415 | 强档原生14ppw整密度复现 |
| strong64 | 4480 / 14 / 64 | 196.0132 | 强档原生网格自求积对照 |
| weak64 | 3840 / 12 / 64 | 276.6899 | 当前代码复现已有弱档64阶整密度 |

两级成本代表先行，3840/qa96保守估计702.32 s；qa128估计超过1600 s，在开算前拒绝。一次明确的成本修订采用已有弱档相同的64/96对照，没有改物理源、海面、点位或误差判据。3840实测242.98 s后才将强档原生4480/qa96估计更新为669.60 s；强档新物理RAM2,032,869,376 bytes通过1,418,723,328 bytes保守增量/70%门。每个随后作业重新读RAM；没有复用前一作业的内存容量来宣称新门通过。全部新边界解最大实耗289.742 s，低于900 s。参考评价/积分闭合批次预检物理RAM837,144,576 bytes，192 MiB保守增量门通过，实耗21.4751801 s。内存增量是估计，不宣称测得峰值。

### 2. 未倒填历史哈希：用当前已知代码作数值再资格

当前未改BIE求解器SHA为`76f1a9d52c898d62507ad7d971813c79b1f09cd9ced4f6dbce556ed11c349fa9`。弱64、中96、强原生14ppw96的**整个边界密度**与对应历史保存层逐点完全一致，三案相对差和最大差均为0，而不是只比中心或只比一条接收线。因此原物理源、法向/高度、h=-2半平面Green、窗和相位数值合同在当前已知代码下重新得到实际支持。历史执行时SHA仍写未知，没有把当前SHA追填成旧执行记录。

在411线+6个离深度点的417个唯一位置重新评价所需64/96密度，并和保存空间/窗口层级严格同点配对。三档自求积变化的支持整体相对量约`3.08e-8/3.09e-8/2.84e-8`；中、强九点/原探针同网格对照约`3.09e-8/2.74e-8`。空间最后一层、窗口相邻变化及自求积变化分别保存，经验分辨率估计为三种绝对变化之和；大约九点`1.51e-7–1.53e-7`、线`2.55e-7–2.63e-7`相对量。它不是严格数学误差界或统计置信区间；线性残差/边界残差也未被直接当作接收场前向误差。

逐点原始层级复数、幅值和误差组成见formal_reference_v2/sensitivity_fields.csv、reference.csv、density_checks.csv及reference_data.mat。当前整密度复现解决“仅凭当下代码标签升级旧数据”的缺口；完整历史执行元数据缺失仍作为溯源限制保留，但不再将已有直接重解再资格的数值一概称为版本未知参考。

### 3. FullFK外场补全：同一完整双项积分的等价计算

原L60→L90变化约2.6e-4说明有限积分尾项尚重要，不是可直接忽略的模型差。对每个固定角谱模态，平面源法向导数是`−i kz u_inc`；出射Green的Weyl积分给出完整Kirchhoff两项各一半镜像反射。因此无限平面FullFK恰为同源解析镜像`u_flat,analytic`。物理C2海面在`|x|>=50 m`处严格平坦，粗糙和平面的真实法向、弧长、入射值/导数及两项Green被积项完全相同，故

`u_FullFK,complete = u_flat,analytic + FK_rough([-L,L]) − FK_flat([-L,L]), L>50 m`。

这不是拟合α、rough/flat归一化、逐点幅相修正、FastFK或扔掉物理外场。被积项之差在外部**恒等为0**，完整平面外部贡献解析计入；它只是原完整Kirchhoff近似的数值等价改写，不是新散射模型。原有限积分值、分项和全部旧图/指标保持不变，候选数组、推导及L60/L90对照在reference_followup_01内。

随后用独立实现直接积分粗糙—平面两项被积项之差，加回解析平面。使用真实向下几何、入水法向、解析弧长Jacobian、Gauss8×400/800面板，`local_radius_m=Inf`；没有receiver局部截断。原来两次有限积分相减+解析尾项与直接差积分对照，417支持最大差依次约`1.50e-12/8.28e-12/8.22e-11`；直接Gauss400→800相对变化不超过`1.5e-12`，两项和闭合不超过`1.56e-17`。旧保存数组解析补全后L60/L90九点最大差`3.93e-17`。这排除了将窗口误差当成FullFK–BIE物理差异的旧解释。

**正式FullFK数值版本**在formal_reference_v2中明确标为`FullFK_complete_analytic_flat_exterior_GL8_v1`，与旧有限L60版分开。初期48现值+60补值的基线仍完整保留；正式PE/BIE72格原值逐点不改，FullFK36格为有来源的独立完整积分数值完成。正式411线FullFK也保存为完整原生场；不只给G或δu。P0的S2a原始FullFK双权重补正表仍是上级data/line_metrics_both_weights.csv，绝不拿新FullFK值冒充P0复算结果。

### 4. 当前主反射复场结果、分辨率与局限

| 原九点等权 E_array | PE–BIE | 完整FullFK–BIE | PE–完整FullFK |
|---|---:|---:|---:|
| flat | 1.52177981e-13 | 4.68282779e-17 | 1.52180753e-13 |
| weak | 1.68586098e-6 | 1.29144435e-7 | 1.67618311e-6 |
| medium | 1.27998616e-5 | 2.66685007e-7 | 1.27734160e-5 |
| strong | 8.05108583e-5 | 1.17048429e-6 | 7.97260855e-5 |

| M99的411点空间等权 E_line | PE–BIE | 完整FullFK–BIE | PE–完整FullFK |
|---|---:|---:|---:|
| flat | 3.85502826e-13 | 3.28614041e-16 | 3.85504548e-13 |
| weak | 0.003115550552 | 2.16654107e-6 | 0.003113832480 |
| medium | 0.015586811401 | 1.07327689e-5 | 0.015578221900 |
| strong | 0.062911245315 | 4.29152410e-5 | 0.062876932200 |

所有量都是同源、未作rough/flat归一化的原生反射单位源响应；原九点坐标、频率、σ、seed/确定性海面、相位和表示不变，无单边共轭或拟合。fields_core.csv逐格含原生/映射实虚部、源/单位/相位/数值版本；fields_line.csv含4932行真实复场，array_metrics.csv的12行显式E_array及最差点，line_metrics.csv空间/入射权重各12行，point_metrics_detailed.csv给出108行候选/参考幅值、分母、dB/相位、平方项、有效/近零标志。当前支持无缺失、零参考或尺度/经验近零标志，不能把“无近零”解释成参考误差为0。

数值分辨率按事前规则明确列在resolutions.csv：BIE空间/窗口/64–96变化分别计入；含PE的方法对另计半步长/倍网格实测变化，后者只在预冻结探针集测量，最大值向支持保守携带，不宣称已做整线倍网格的严格上界；含FullFK的方法对另计Gauss加密，平面另有固定舍入尺度。普通模型差异与5倍该经验估计比较，**不是统计显著性或通用物理容差**。平面差异均未分辨；弱、中档FullFK–BIE九点未达到5倍判据，不能据此声称1e-7级物理精度。强档九点与三档411线的FullFK–BIE差异达到本固定规则；PE三档在两种支持均可分辨。结论只限实际支持，不延伸到PM/3D/宽带/通信。

### 5. 真实成本、独立审计和完成性核验

七次新BIE边界解总实计约1213.9591831 s；加初期全部成功声学批次/三次S2a23.2257418 s、参考评价/完整积分闭合21.4751801 s，**总实计1258.660105 s**。另保留旧未计时/失败保守扣记150 s，**累计扣记1408.660105/7200 s，余5791.339895 s**。最大单作业289.742 s。没有把图形/I-O的未单独仪器时间虚构为已测，也没有把旧成本或失败从预算中删除；所有新边界求解和密度新点评价均计入。

formal_reference_v2/audit/independent_audit.json从CSV独立复核108核心场、4932线场、24线指标、12阵列指标、15429敏感性复数及G/δu/交叉恒等式；不调用主指标函数重算指标。MAT/CSV场最大差5.02e-16、独立指标1.23e-14、BIE分辨率组成1.51e-15、FullFK数值对照3.74e-16、G/δu9.05e-16、交叉项2.16e-19。resolution_audit.json另从原PE敏感性CSV和本次参考CSV独立复算24行总分辨率及全部5倍分类，最大差9.97e-16；逐点幅相平方分解最大余量2.87e-22。审计新增物理场为0。

正式F01–F07各FIG/PNG在formal_reference_v2/figures；全部7张逐张实际查看，坐标/单位/相位wrap/色标/截断与局部轴放大已检查；FIG全部重开，保存MAT重绘PNG像素最大差0。F03仅已有强档PE/BIE411×45真场，FullFK连续x-z未计算而明确缺失，未由九点/线插图。18个新增MATLAB文件最终Code Analyzer零消息；原S2a的30个保存产物完整性复核及最终源码/输入/输出哈希见provenance，git diff --check通过。仅新增比较/参考/审计/绘图代码，公共API/生产求解器未改，未作范围外的通信回归或改写vertical_comm_guide.md实现。

完成条件逐条对应证据见formal_reference_v2/completion_requirements.csv和final_certificate.mat/summary.json。原MAT内的pending_audit字段是审计前不可变快照，最终证书与独立审计覆盖其当前状态，不回写物理数组。下次只需纯保存数据重绘时：

```matlab
setup_vertical_project;
base = fullfile('E:/MISC/CARPE3D_matlab/Explain/results/comparisons/two_dimensional_channel_closeout','20261023_s2b_1305');
render_two_dimensional_closeout_s2b(fullfile(base,'formal_reference_v2/plot_data.mat'),fullfile(base,'unused_redraw_directory'));
```

当前保留的问题是：历史执行时SHA无法追补；经验误差估计不是严格界；弱中档九点FullFK差异尚未分辨；PE倍网格仅有冻结探针支持；FullFK连续x-z图未计算；源周期/离散Gaussian及紧支撑平坦外场假设不能默认跨场景。它们均已限定本次结论，没有制造PASS。若以后需要更小差异的物理判断，应另定更严格参考及资源预算；本Goal不继续。S2a状态/旧8.0908%、6.2911%、4.39622%保持原口径，`G1_target_PM=OPEN`、`WP2_formal=false`，不自动执行S3–S6/Bellhop源映射/通信，不提交推送。
