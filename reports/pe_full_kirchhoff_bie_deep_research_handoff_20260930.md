# PE–Full Kirchhoff–BIE 粗糙海面反射：深度研究交接报告

日期：2026-09-30  
性质：现有代码、报告及保存结果的证据综合；本文没有运行求解器，也没有引入新的物理模型或验证结果。  
研究范围：垂直水声信道中，局部 PE 海面反射近似与二维 Full Kirchhoff、Helmholtz BIE 的差异，以及低成本非局部反射算子的可行性。

## 1. 给研究者的结论摘要

在已经验证的**二维、单横向维、4 kHz、光滑确定性正弦海面**范围内，使用统一的场表示后，Full Kirchhoff（下文简称 FK）与独立 Helmholtz BIE 的归一化复误差约为 `2.54e-4` 至 `3.32e-4`；同一批参数扫描中的局部 PE 相位屏对 BIE 误差为 `1.58e-3` 至 `8.09e-2`。所以，此范围内最值得研究的是 PE 局部反射近似丢失了哪些传播路径和相干空间耦合，而不是预设 FK 已失效。此判断不能推广为任意粗糙度、频率、三维表面或随机 PM 海况的定理。

入射场在 `z=0` 与真实粗糙表面上，PE 与 FK 的固定掩膜复误差约 `2.6e-13`。较大的反射场差异发生在共同入射场之后。`k_z` 感知的局部 Model-1 能在宽角源、低表面波数条件下减少误差，但强高度或高表面波数的残差仍明显。接收深度 `z=1 m` 的历史结果存在近边界 BIE 计算和 FK 积分收敛问题，必须保留 **UNVERIFIED** 标签。

Reduced Kirchhoff 的已测版本是对完整 FK 积分按 `|x_r-x_s|<L` 截断，保留真实曲面、法向、弧长测度和 Green 核。它能在部分有限孔径下优于 PE，但 `L<=64 m` 没有在三个代表案例同时达到 `E_c<1e-3`。后续 11 个扫描案例都首次在 `L=128 m` 通过这一门槛，而该孔径已覆盖相关曲面支撑，因此还没有证据支持一个随海面几何变化且真正局域的最小孔径规律。

## 2. 项目主线与本研究的位置

项目公开入口为 `vertical_channel_model(paramsV)`；主要实现由 `src/channel/vertical_channel_model_impl.m`、`src/propagation/vertical_wape_propagator.m` 与 `src/surface/pm_surface_boundary_model.m` 组成。主线把信道分为直达与海面反射分量，并输出 `H_f`、`h_direct`、`h_reflect`、`h_total` 等字段；`H_f=H_direct_f+H_reflect_f`，参考频点满足 `h_total=H_f(idx_f_ref)`。通信演示 `comm_main_vertical_psk.m` 消费频域信道。坐标采用海面 `z=0`、正 `z` 向下，Tx 比 Rx 深。

本报告的 PE–FK–BIE 工作是**验证支线**。二维单横向 PE 辅助程序 `scripts/validation/support/run_pe_1d_surface_reflection_validation.m` 中，Model-0 使用局部 `R_0 exp(+i 2 k_0 eta(x))` 相位屏，Model-1 对入射角谱各分量使用 `R_0 exp(+i 2 k_z(k_x) eta(x))`。FK 与 BIE 分别由 `scripts/validation/support/solve_full_kirchhoff_2d.m` 和 `solve_helmholtz_bie_halfplane_vertical.m` 支撑。这些结果没有把 Model-1、FK 或 Reduced Kirchhoff 接入生产 WAPE 海面边界，也不直接验证三维生产信道的通信误码率。

项目现行背景参见本 checkout 的 `PROJECT_CONTEXT.md`、`vertical_comm_guide.md` 与 `reports/pe_surface_operator_improvement_final_report.md`。本研究所综合的较新验证文件保存在另一 checkout，路径及状态见下节。两个 checkout 的报告集不同，不能把后一 checkout 的未合并文件说成当前 checkout 已有的生产项目报告。

## 3. 证据位置与版本边界

以下相对路径均以 `C:\Users\ASUS\.codex\worktrees\baf5\Explain` 为根。读取时该 worktree 位于 `cl-verfication` 分支、HEAD `ee21d3dee73d3ea0ef385f5a0247eefddf61704c`，工作树另有未提交改动；其中 receiver-depth audit 和 incident/reflection covariance audit 报告为未跟踪文件。本文记录的是 **2026-09-30 当时磁盘文件的内容**，不是该 HEAD 单独可复现的快照。五份报告所指向的 CSV 均已检查存在。

| 主题 | 报告 | 直接保存数据 |
|---|---|---|
| FK–BIE 参数扫描 | `reports/pe_bie_full_kirchhoff_parameter_sweep_report.md` | `results/validation/pe_bie_full_kirchhoff_parameter_sweep/pe_bie_full_kirchhoff_parameter_sweep.csv` |
| Reduced Kirchhoff Stage-1 | `reports/pe_bie_full_kirchhoff_reduced_operator_stage1_report.md` | `results/validation/pe_bie_full_kirchhoff_reduced_stage1/pe_bie_full_kirchhoff_reduced_stage1_summary.csv` |
| Gaussian `sigma` 与 `k_z`-aware | `reports/pe_bie_full_kirchhoff_gaussian_sigma_kzaware_comparison_report.md` | `results/validation/pe_bie_full_kirchhoff_gaussian_sigma_kzaware/gaussian_sigma_kzaware_master.csv` |
| 接收深度数值审查 | `reports/pe_bie_receiver_depth_numerical_audit_report.md` | `results/validation/pe_bie_receiver_depth_numerical_audit/receiver_depth_numerical_audit.csv` |
| 入射/反射场一致性审查 | `reports/pe_bie_incident_reflection_covariance_report.md` | `results/validation/pe_bie_incident_reflection_covariance/incident_reflection_covariance.csv` |

必要的先验与后续证据还包括：`reports/pe_bie_full_kirchhoff_convention_audit_report.md`（符号与共轭约定）、`reports/pe_bie_full_kirchhoff_validation_report.md`（平面 Gate-0 与三案例）、`reports/pe_bie_error_decomposition_report.md`（幅度/相位拆分）、`reports/pe_bie_full_kirchhoff_gaussian_sigma_sensitivity_report.md`、`reports/pe_bie_full_kirchhoff_receiver_depth_sensitivity_report.md`、`reports/pe_bie_full_kirchhoff_finite_aperture_report.md` 和 `reports/pe_bie_full_kirchhoff_adaptive_aperture_stage1_report.md`。全部路径均相对于上述 `baf5` checkout；不要读取仓库根目录的 `cash/` 作为当前证据。

## 4. 几何、场表示、指标与不能混用的案例标签

典型固定条件为 `f=4000 Hz`、`c=1500 m/s`、Tx 深度 `100 m`、压力释放 Dirichlet 海面，正弦形状采用指定 C2 截断。参数扫描固定 Rx 深度 `3 m`、Gaussian 源宽 `sigma=0.3 m`；`sigma` 研究把源宽改为 `0.3,0.5,1,2,5 m`；深度和入射/反射一致性审查固定 `sigma=1 m`。只在相同 `A,K,sigma,z_rx`、掩膜和场约定下比较误差。

权威几何是 `z_s=+eta(x)`，水侧单位法向为 `(-eta'(x),1)/sqrt(1+eta'(x)^2)`。FK/BIE 原生求解采用 `exp(-i omega t)`、出射二维 Green 函数 `G=(i/4)H_0^(1)(kR)`。比较使用粗/平反射场比值，并将 FK/BIE **原生比值各共轭一次**映射到冻结的 PE 比较表示。旧的“大相位异常”来自原生 FK 比值与已共轭的 BIE 比值直接比较；不能再次共轭已保存的 `G_FK` 或 `G_BIE`，也不能逐案例选符号、拟合复常数或移动相位零点。

主要复误差 `E_c(a,b)` 的实现为 `sqrt(sum(w*|a-b|^2)/sum(w*|b|^2))`，其中 `b` 是表中参照场，权重在有效掩膜上重新归一。相位 RMS 使用 `angle(a*conj(b))` 的加权平方平均；“magnitude L2”使用 `|a|-|b|`，与复误差共用参照场能量。典型掩膜为 Stage-0 的 99% 入射能量足迹、有限值、比较双方均高于这对场共同峰值的 `-40 dB`。`sigma` 扫描采用各源宽自己的 M99 足迹；接收深度审查固定 `sigma=1` 的 M99 足迹和权重。**成对掩膜与四场共同掩膜不是同一指标**；不能直接把不同报告中的末位数字差别解释为物理变化。

尤其注意 `weak_low_K` 不是跨文件唯一案例 ID。FK 参数扫描和 Reduced Stage-1 的弱低 K 案例对应 `A=0.01 m,K=0.10 rad/m`，其 PE/BIE 误差约 `0.00316`。`sigma+k_z`、receiver-depth 和 incident/reflection 审查中的同名案例对应 `A=0.05 m,K=0.10 rad/m`，在 `sigma=0.3 m` 时 PE/BIE 误差 `0.01609`。深度审查改为 `sigma=1 m`，`z=3 m` 的 Model-1/BIE 误差为 `0.0102641`。研究数据索引必须以完整数值四元组 `(A,K,sigma,z_rx)` 为主键。

## 5. 综合证据链

### 5.1 FK 与 BIE 的受控闭合范围

参数扫描以 `K=0.10 rad/m` 固定高度 `A=0.005–0.20 m`，以及 `A=0.02 m` 固定波数 `K=0.05–0.70 rad/m`。FK/BIE 的 `E_c` 为 `0.000253613–0.00033158`，FK 相位 RMS 为 `0.000178554–0.000279571 rad`；2049 到 4097 点的 FK 积分改变量最大复 L2 为 `3.72e-5`。同一比较中 PE/BIE 的 `E_c` 为 `0.00157783–0.0809081`，且在 `A=0.02,K=0.70` 达 `0.0484307`。平面 Gate-0 的 FK/平面参考复 L2 为 `0.000185546`。所以已测扫描没有发现 FK 失效边界；FK/BIE 的小残差不能被称为零，也不能据此确定外推极限。

### 5.2 有限角度修正与剩余非局部误差

在 `A=0.05,K=0.10,z=3 m`，`sigma=0.3 m` 的 `theta99=21.2801°`；Model-0/BIE `E_c=0.01609`，Model-1/BIE `0.00398989`，改善约 `4.03` 倍。`sigma=0.5 m` 改善约 `2.80` 倍；`sigma>=1 m` 时此案例 Model-1 与 Model-0 的复误差接近。强高度低 K (`A=0.20,K=0.10`) 在 `sigma=0.3` 时从 `0.0809081` 降到 `0.0507808`，但 `sigma=1` 时两者仍约 `0.0713`。弱高度高 K (`A=0.02,K=0.47`) 在 `sigma=0.3` 仅从 `0.0229061` 降到 `0.0215668`，较窄源时仍约 `0.025–0.029`。对应 FK/BIE 在这些源宽和案例中保持约 `5e-6–2.7e-4`。结论是有限角度解释部分低 K 偏差；源变窄不保证 PE 接近 BIE，高 K 幅度残差尤其如此。

### 5.3 误差从哪里开始

入射场审查在三种表面上给出 PE/FK 参考平面误差 `2.62103e-13`，真实表面误差约 `2.63e-13–2.66e-13`。BIE 使用相同角谱构造作为源项，但没有另存独立入射边界迹，所以不能称“三种独立入射解已数值一致”。在 `sigma=1,z=3 m` 时，三案例的 PE/FK 反射场误差分别约 `0.0102641`、`0.0713035`、`0.0252701`，而 FK/BIE 分别约 `8.20e-6`、`3.28e-5`、`9.68e-5`。这把差异定位到共同入射构造之后的反射/空间频谱重建链条；不能仅凭此进一步把全部误差唯一归因于某个单独的算子项。

### 5.4 Reduced Kirchhoff 的实际边界

Stage-1 消融显示，保留真实表面高度但把曲面法向改成平法向，三个代表案例的 FK 相对误差仅 `2.12e-6`、`1.32e-4`、`9.81e-5`；仅把 `ds` 改为 `dx` 则为 `3.13e-7`、`1.25e-4`、`2.70e-5`。把物理海面整体换到 `z=0` 的误差却约为 `0.230`、`1.217`、`0.464`。这说明真实高度/路径相位至关重要，但 **PE 已有局部 `2k eta` 高度相位**；“平面消融很差”本身不能证明 PE 完全遗漏高度项。问题是局部代理未能重现完整的曲面至接收机路径相位及相干积分。

以“把 PE/FK 误差至少减半”为宽松目标，接收机中心截断所需的最小已测 `L` 对弱低 K (`A=0.01`)、强高 A (`A=0.20`)、弱高 K (`K=0.47`) 分别为 `64,8,32 m`。若目标改成所有三案例 Reduced-Kirchhoff/BIE `E_c<1e-3`，`L<=64 m` 全部失败，完整 FK 才同时通过。后续 11 个唯一参数点的 `L=4,8,16,32,64,128,Inf m` 数据库在 `L=128 m` 才全部首过 `1e-3`，且该长度实质覆盖固定足迹下完整相关支撑。因此“宽松改善所需支撑因案例而异”与“严格误差门槛未得到有效短孔径”并不矛盾。低秩或局部平稳核近似是研究候选，尚无已验证的快速生成器、通用小卷积核或生产实现。

### 5.5 接收深度：可信趋势与未验证边界

深度审查保持 `sigma=1` 的固定 117 个 M99 样本、权重和离散 `z_rx=1,2,3,5,10,20,30 m`，不插值。`z=3→30 m` 的 Model-1/BIE `E_c`：低 K 中高度 (`A=0.05`) 从 `0.0102641` 到 `0.0676774`；强高度 (`A=0.20`) 从 `0.0713033` 到 `0.545856`；高 K (`A=0.02,K=0.47`) 从 `0.0252671` 到 `0.187208`。相应 FK/BIE 在 `z>=2 m` 保持约 `1e-5–1e-4`。30 m 处强高度案例是相位贡献较大，低 K 中高度案例则幅度贡献较大，不能使用单一幅度修正解释。

此处证实的是固定归一化比较下的可复现随深度变化；还须结合平面解析镜像场控制、PE 包络与解析/FK 载波约定，才能把趋势提升为物理“接收深度定律”。`z=1 m` 的 FK 标称/细化差达到 `0.172024`，保存的 BIE 评估对亚米距离用过 `max(r,1)`；其 PE/BIE 与 FK/BIE 巨大值 **UNVERIFIED**，不能纳入模型拟合或门槛判断。

## 6. 研究假设、反证要求与优先顺序

1. **优先检验完整路径相位及非局部核。** 用同一入射角谱、同一真实海面、同一平面参考和同一接收网格，把 FK 核分解为曲面位置、传播 Green 核、入射/出射法向导数和测度。分解必须能从原核重构完整 FK；任何近似都要在留出的 `A,K,sigma,z` 案例上盲测。
2. **分开有限角与表面空间耦合。** 以 `theta95/theta99` 而非仅 `sigma` 表示源角宽；比较 Model-0、Model-1、FK、BIE 的复数、幅度、相位和角谱误差。宽角修正已经解释的部分不可再次计为非局部收益。
3. **避免把能量极弱频谱箱的相位当主证据。** 入射/反射审查的 `|k_x|/k` 分带相位统计未加权，低振幅箱可能相位不稳。需另给能量加权频谱误差、频带贡献、共同掩膜和对掩膜阈值的敏感性。
4. **先修复近边界数值控制。** 若研究 `z=1 m`，应去掉 BIE 亚米距离裁剪，处理真实近奇异积分，并证明 FK 网格收敛，再重新审查平面与粗糙控制。当前 `z>=2 m` 可用于已核对的离散深度问题。
5. **给降阶方法独立复杂度目标。** 对低秩、局部平稳或其他近似同时报告 `E_c`、相位/幅度、算子装配和每接收点计算成本、内存及参数外推表现。`L=128 m` 在现有足迹下近乎完整积分，不应被包装成已实现的局部快速法。

现有证据**不支持**宣称：FK 对所有海况等于 BIE；已找出 FK 失效边界；Model-1 已取代生产海面算子；receiver-depth 增长全由非局部散射造成；`z=1 m` 结果可信；现有 Reduced Kirchhoff 已达到严格目标且有通用自适应孔径；或这些二维单频结果直接给出宽带通信链 BER 预测。

## 7. 可直接交给 GPT 深度研究的任务书

> 请从上列源文件和保存 CSV/MAT 出发，研究二维垂直水声链路粗糙海面反射中，局部 PE 相位屏、`k_z` 感知相位屏、完整 FK 积分与 Helmholtz BIE 的关系。先核对每个数据点的 `(A,K,sigma,z_rx)`、曲面、源、平面参考、原生/共轭表示、掩膜与权重；不要使用 `cash/`。以已测范围的 FK/BIE 数值闭合作为受限基准，提出能减少 PE 残差且保持物理可解释性的非局部算子简化。逐项列出假设、推导、可能的反例、与现存代码的接口关系、离散误差和计算成本。明确区分“从已有数据推断”“文献理论支持”“需要新计算验证”；不要把后两者写成项目已完成结论。

建议交付顺序：

1. **证据复核表**：逐数据点列出文件/字段、四元组、比较场与掩膜、基线和不确定性；解释 `weak_low_K` 标签复用与 `z=1` 排除。
2. **理论推导**：给出 FK 边界积分到局部相位屏所需的具体近似步骤，定位可能产生高 K 幅度偏差、强高度相位偏差与随深度变化的项；用代码中的法向和 Green 约定核对符号。
3. **候选降阶算子**：至少比较“硬截断孔径”“低秩/分离核”“局部平稳或角谱耦合核”，说明是否保留真实表面路径相位、法向导数和非局部相干耦合。
4. **最小验证矩阵**：含平面 Gate-0、`A=0.01/0.05/0.20` 低 K、`K=0.47/0.70`、`sigma=0.3/1/5`、`z=2/3/10/30 m` 的适当留出组合；只在修正近边界求解后加入 `z=1 m`。记录 `E_c`、幅度、相位、场图、频谱加权误差、数值收敛与成本。
5. **结论等级**：每条结论标为 `VERIFIED_IN_SAVED_CASES`、`INFERENCE` 或 `UNVERIFIED`，并列明不能外推的参数和物理条件。

## 8. 文档维护说明

这是跨 checkout 的综合研究交接文件，不改变 MATLAB 接口、生产模型或保存结果。当前 checkout 的 `PROJECT_CONTEXT.md`、`vertical_comm_guide.md` 和 README 不需要因本文而把验证支线改写成生产能力。若未来将某个降阶反射算子接入 `src/`，应再按项目规则更新主线文档、公共接口说明和回归验证报告。
