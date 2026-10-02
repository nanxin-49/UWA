# Full Kirchhoff spectral separability / low-rank audit

日期：2026-09-30。状态：**VERIFIED_IN_SAVED_CASES；建议继续开发 Spectral Fast FK prototype**。本阶段已完成并停止；没有开发生产 Fast FK。

## 结论与版本边界

真实入射能量加权、保留全部已保存传播频谱时，15 个受控案例达到矩阵相对 Frobenius 误差 `1e-2 / 1e-3 / 1e-4` 分别只需 rank **1–2 / 1–2 / 1–3**。强高度宽角案例为 **2 / 2 / 3**。三个完整 FK 运算检查中，加权 rank 2 的接收场粗/平比值误差全部小于 `1e-3`。固定物理区间加密没有引起 rank 增长。

本地基线是干净的 main、HEAD `02a29d92b84eab9f07fb115969022dff4960788e`，已通过 `e78c24d` 整合研究支线。GitHub 查询显示 [PR #3](https://github.com/nanxin-49/UWA/pull/3) 仍为 Draft/open、merged=false，head SHA 与本地 HEAD 相同；因此“本地已整合”和“GitHub PR 已合并”必须区分。没有同步、提交、推送或改变 PR。

结论仅涉及二维单横向、4 kHz、c=1500 m/s、Tx=100 m、Rx=3 m、C2 截断确定性正弦面及保存 Gaussian 入射谱。没有读取 cash/，没有重跑 PE、BIE、PM Monte Carlo、宽带数据库或 z=1 m 结果。交接报告中的旧 worktree 路径现在是历史 provenance；本轮使用主线结果目录。

## 数学对象、范数与权重

令 `kz_j=sqrt(k^2-kx_j^2)`、`delta_j=k-kz_j`，审查 `C_mj=exp(i delta_j eta_m)`，其中 `eta_m` 均匀覆盖 `[-A,A]`，默认 257 点。完整高度相位为 `exp(-i kz_j eta_m)=exp(-i k eta_m) C_mj`。强高度载波完整保留，绝不进行 small-height Taylor expansion。左乘单位模对角矩阵不改变奇异值，因此完整高度矩阵和残余矩阵的 SVD rank 理论上相同；分解的价值是揭示控制耦合的是 delta 的跨度，而不是去掉某个物理大相位。

未加权矩阵使用 `B=C/sqrt(n_kx)`；加权使用 `B=C diag(sqrt(w_j))`，`w_j=|co_j|^2/sum_j |co_j|^2`。使用**能量的平方根**乘列，不是直接以能量乘列。误差定义为 `||B-B_r||_F/||B||_F=sqrt(sum_{j>r} s_j^2/sum_j s_j^2)`，报告满足门槛的最小正整数 rank。eta 方向为均匀测度；它不是表面实际高度分布的概率加权。该范数是平均矩阵误差，不是最坏入射算子范数或每个接收点的误差上界。尾能量从小奇异值直接累加，避免 total-minus-cumsum 的抵消。

三个频域：`effective99` 为真实谱的最小对称 99% 能量区间；`saved_propagating` 保留所有已保存传播箱；`full_propagating_control` 为 1025 点均匀 `[-k,k]`、仅未加权数学对照。95%/99% 分位数按 |kx| 排序，边界相同的正负箱全部保留，实际能量可略大于 99%。没有对入射场加 M99/Hann 窗来生成权重：优先使用保存源的原始复角谱系数，避免窗口展宽。已有 covariance CSV 仅作为交叉背景，其 Hann 谱不是这里的权重。

`saved_propagating` 的源网格 N=984、宽度 192.1875 m；所有这些箱都为传播波，范围约 `[-0.960,0.958] k`。没有丢弃已保存尾部，但未采样的近掠射区间不等于已实测工程谱。数学 full-propagating 对照补充这一差别。截掉 1% 能量一般可产生约 0.1 的谱 L2 误差，所以 effective99 rank 不能单独作为严格场误差保证。

## 数据复用与实际角宽

- `results/validation/pe_bie_full_kirchhoff_parameter_sweep/pe_bie_full_kirchhoff_parameter_sweep.mat`
- `results/validation/pe_bie_full_kirchhoff_gaussian_sigma_sensitivity/gaussian_sigma_sensitivity.mat`
- `results/validation/pe_bellhop_controlled_comparison/stage0/stage0_validation.mat`

复用完整 sweep 的 11 个唯一 A–K 点；相同 `(A=.02,K=.1)` 重复行去重。额外仅对 strong-height A=.2、K=.1 复用 sigma=.5/1/2/5 的数值源谱，不重跑原 sweep。保存 MAT 中过期函数句柄只加载而不调用；代码仅定向抑制已确认的 `UnresolvedFunctionHandle` 路径警告，读取数值 co/kx/kz。sigma=.3 的保存 co 与原 Gaussian 离散 FFT 闭合在 `1e-12` 内。

|sigma (m)|theta95 (deg)|theta99 (deg)|kx95/k|kx99/k|99% 区间实际能量|
|---:|---:|---:|---:|---:|---:|
|0.3|15.96928|21.28005|0.275122|0.3629268|0.990314569|
|0.5|9.546829|12.50855|0.1658537|0.2165854|0.990052586|
|1|4.700731|6.273138|0.08195122|0.1092683|0.991011819|
|2|2.348387|3.131865|0.04097561|0.05463415|0.991619115|
|5|0.8944095|1.229858|0.01560976|0.02146341|0.992288206|

## 参数与 rank 结果

各 rank 元组依次对应 `(1e-2,1e-3,1e-4)`。U99/W99 为有效区间未加权/加权，Uall/Wall 为全部保存频谱未加权/加权，Ufull 为完整传播区间数学对照。每行完整 theta95/99、样本数、保留能量、耗时和小角矩阵误差均见 `rank_summary.csv`。

|A (m)|K (rad/m)|sigma (m)|kA|AK|AK²/k|K/k|theta99 (deg)|chi|U99|W99|Uall|Wall|Ufull|
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---|---|---|
|0.005|0.1|0.3|0.0837758|0.0005|2.98416e-06|0.00596831|21.2801|0.00571203|1/1/2|1/1/2|1/2/2|1/1/2|2/2/2|
|0.01|0.1|0.3|0.167552|0.001|5.96831e-06|0.00596831|21.2801|0.0114241|1/2/2|1/2/2|2/2/3|1/2/2|2/2/3|
|0.02|0.1|0.3|0.335103|0.002|1.19366e-05|0.00596831|21.2801|0.0228481|1/2/2|1/2/2|2/2/3|1/2/2|2/2/3|
|0.05|0.1|0.3|0.837758|0.005|2.98416e-05|0.00596831|21.2801|0.0571203|1/2/2|1/2/2|2/3/3|1/2/2|2/3/4|
|0.1|0.1|0.3|1.67552|0.01|5.96831e-05|0.00596831|21.2801|0.114241|2/2/3|2/2/2|3/3/4|2/2/3|3/4/4|
|0.2|0.1|0.3|3.35103|0.02|0.000119366|0.00596831|21.2801|0.228481|2/2/3|2/2/3|3/4/5|2/2/3|4/5/6|
|0.02|0.05|0.3|0.335103|0.001|2.98416e-06|0.00298416|21.2801|0.0228481|1/2/2|1/2/2|2/2/3|1/2/2|2/2/3|
|0.02|0.2|0.3|0.335103|0.004|4.77465e-05|0.0119366|21.2801|0.0228481|1/2/2|1/2/2|2/2/3|1/2/2|2/2/3|
|0.02|0.3|0.3|0.335103|0.006|0.00010743|0.0179049|21.2801|0.0228481|1/2/2|1/2/2|2/2/3|1/2/2|2/2/3|
|0.02|0.47|0.3|0.335103|0.0094|0.00026368|0.0280511|21.2801|0.0228481|1/2/2|1/2/2|2/2/3|1/2/2|2/2/3|
|0.02|0.7|0.3|0.335103|0.014|0.000584894|0.0417782|21.2801|0.0228481|1/2/2|1/2/2|2/2/3|1/2/2|2/2/3|
|0.2|0.1|0.5|3.35103|0.02|0.000119366|0.00596831|12.5086|0.0795412|2/2/2|1/2/2|3/4/5|1/2/3|4/5/6|
|0.2|0.1|1|3.35103|0.02|0.000119366|0.00596831|6.27314|0.020065|1/2/2|1/2/2|3/4/5|1/2/2|4/5/6|
|0.2|0.1|2|3.35103|0.02|0.000119366|0.00596831|3.13187|0.00500497|1/1/2|1/1/2|3/4/5|1/1/2|4/5/6|
|0.2|0.1|5|3.35103|0.02|0.000119366|0.00596831|1.22986|0.000771963|1/1/2|1/1/1|3/4/5|1/1/1|4/5/6|

### 参数控制量解释

`chi=A[k-kz(theta99)]=kA(1-cos(theta99))`；小角估计为 `kA theta99^2/2`（theta 用弧度），同时保留 `delta≈kx²/(2k)` 的真实矩阵误差。固定 A=.2 和 k 时，sigma=.3→5 使 chi=.22848→.00077196、全谱加权 rank 从 2/2/3 降至 1/1/1，而 kA 始终为 3.35103。因此 **chi 连同真实能量分布比单独 kA 更能解释 rank**。相同 chi 也不保证相同加权 rank，因为尾部谱能量及范数权重仍有影响。小角式只作诊断，SVD 与算子检查均使用精确 kz。

在均匀 eta 区间矩阵的定义中，K 不出现。因此固定 A、sigma 的 K=.05–.70 扫描 rank 完全相同，是数学构造的必然结果，不能解读为 K 对完整 FK 不重要。AK、AK²/k、K/k 会通过实际 eta(x)、法向、散射出射带宽和 Green 核影响完整运算。矩阵层面没有强行拟合经验公式；K=.47 的真实曲面检查提供有限的算子支持，K=.70 尚未做完整 reconstruction 检查。

小角诊断，固定 A=.2、K=.1；末列为 effective99 域的能量加权相对 Frobenius 误差 `||C_small-C||/||C||`：

|sigma|chi exact|chi small-angle|small-angle matrix error|
|---:|---:|---:|---:|
|0.3|0.22848128|0.23112594|0.00059175264|
|0.5|0.079541161|0.079857837|7.3699629e-05|
|1|0.020064996|0.020085052|4.6733741e-06|
|2|0.0050049687|0.0050062151|2.9630702e-07|
|5|0.0007719626|0.00077199224|7.7061501e-09|

宽角端的小角近似自身误差不保证达到 1e-4，不能以低秩成功替代精确 kz；实际原型应保留精确谱相位。

## 网格稳定性

固定 A=.2、sigma=.3 或 1，保持物理 eta/kx 区间不变，把 `(n_eta,n_kx)` 从 `(129,257)` 加密到 `(257,513)`、`(513,1025)`。权重由保存能量谱线性插值，并在各网格归一；没有生成新源场或重新求解。这是离散积分/采样密度稳定性，而不是新的物理 Nyquist 带宽或真实传播网格收敛。

|sigma|domain|weighted|rank tuples at the three densities|
|---:|---|---|---|
|0.3|effective99|0|2/2/3 2/2/3 2/2/3 |
|0.3|effective99|1|2/2/3 2/2/3 2/2/3 |
|0.3|saved_propagating|0|3/4/5 3/4/5 3/4/5 |
|0.3|saved_propagating|1|2/2/3 2/2/3 2/2/3 |
|1|effective99|0|1/2/2 1/2/2 1/2/2 |
|1|effective99|1|1/2/2 1/2/2 1/2/2 |
|1|saved_propagating|0|3/4/5 3/4/5 3/4/5 |
|1|saved_propagating|1|1/2/2 1/2/2 1/2/2 |

所有测试 rank 均稳定，未出现随采样数显著增长。

## 完整 FK 运算 sanity check

使用参数扫描中 weak `(A=.01,K=.1)`、strong-height `(.2,.1)`、high-K `(.02,.47)`，sigma=.3、Rx=3 m。保持原保存曲面 4097 点、完整 [-60,60] m 支撑、解析 C2 法向、原弧长梯形测度、全部 984 个传播谱箱与固定 M99 接收 411 点/原能量权重。没有额外 -40 dB 成对筛选，没有逐例拟合、归一相位或载波校准。

以 uniform-eta SVD 的右奇异子空间 V_r 投影实际曲面耦合：`C_r(eta(x))=[C(eta(x)) diag(sqrt(w)) V_r] V_r^H diag(1/sqrt(w))`。这是同一右子空间到保存高度样本的直接延拓，不是插值左奇异向量。显式计算 `u=sum_j co_j exp(i kx_j(x-x0)) exp(i kz_j zTx) exp(-ik eta) C_r/n`；ux、uz 使用同一个 C_r，并分别完整保留 `i kx_j`、`-i kz_j` 乘子；物理坐标导数不对曲面限制后的 u(x,eta(x)) 做切向求导。法向导数为 `n_x ux+n_z uz`。

先调用原 `solve_full_kirchhoff_2d` 得到未截断基线；再只在本脚本装配相同的全支撑 Green/导数线性收缩，并确认与原 solver 相对差 `<1e-11`，避免每个 rank 重算 Bessel 核。以近似后的 u/dn 作用到该完整收缩，出射几何/Green 完全不近似。场比值用固定 pressure-release 平面参考，`G=conj(u_ref/flat)` 一次共轭，未再次共轭已保存 G。

接收原生场和粗/平比值均使用 `sqrt(sum(w |a-b|²)/sum(w |b|²))`；没有相位重定。矩阵误差、入射场误差、法向导数误差、receiver_native_l2、receiver_ratio_l2 都保存在 operator CSV。

|case|saved incident closure|saved FK ratio closure|weighted r=1 FK error|r=2 matrix error|r=2 FK error|r=3 FK error|r=4 FK error|
|---|---:|---:|---:|---:|---:|---:|---:|
|height_A0.01_K0.1|9.504e-14|6.03e-14|0.000965946|2.32036e-06|1.66368e-06|3.51296e-09|8.74459e-12|
|height_A0.2_K0.1|9.531e-14|6.294e-14|0.0193193|0.000922052|0.000661285|2.76674e-05|1.3473e-06|
|wavenumber_A0.02_K0.47|9.554e-14|6.278e-14|0.00291836|9.28097e-06|7.17671e-06|2.64516e-08|1.4696e-10|

场误差随 rank 下降，与矩阵截断误差趋势一致；加权 rank 2 所有案例通过 1e-3，rank 3 所有案例通过 1e-4。强高度 r=2 的 grid Frobenius 误差约 9.22e-4，而实际高度样本 Frobenius 误差约 1.076e-3、FK ratio 误差约 6.61e-4，说明趋势一致不意味着数值相等或严格场误差上界。未加权强高度 r=4 约 3.80e-4、r=5 约 1.96e-5；完整 CSV 包括两种加权方式的 r=1–8。机器精度附近不再解释微小非单调变化。

## CPU、内存与复现

MATLAB `25.2.0.3312555 (R2025b) Update 6`，共享 MCP 会话，标准 MATLAB FFT/SVD/Bessel/graphics，无新增工具箱。完整最终一次运行（加载、矩阵 audit、网格 audit、三个算子检查、PNG/PDF 导出；不含最后保存/report 写入）wall **33.906 s**，进程 CPU **122.234 s**。CPU 包括可能的多线程计算，因此可大于 wall；它不是实时性能或生产速度比较。未运行重复基准、未测峰值 RSS、未证明加速比。

主 audit 残余矩阵 complex-double payload 为 0.090–4.020 MiB；最大 U/S/V payload 合计 6.035 MiB。每个 sanity 的实际曲面耦合矩阵 61.515 MiB，两块 FK Green/导数收缩矩阵合计 51.388 MiB。还有 phase、C_r、Bessel 临时数组及 SVD 工作空间，所以这些是**数组存储核算，不是进程峰值**。原型仍有 dense 矩阵和 O(N_surface N_kx)、O(N_receiver N_surface) 工作；低秩重构的误差检查不是 FFT 性能实现。

运行：先将项目根加入 MATLAB path 并执行 `setup_vertical_project`，再调用 `audit=validate_full_kirchhoff_spectral_rank_audit()`。入口通过 mfilename 定位根目录，不改变 cwd，输出路径明确；重新运行会更新本 audit 目录和本报告。数值数据无 RNG。输入和代码 SHA-256 见目录中的 `provenance.json`（由本轮仓库工具记录，不属于求解器）。

仅从本 audit 保存数据重绘/更新报告：`validate_full_kirchhoff_spectral_rank_audit(true)`；此模式不重算 SVD/FK，报告中的运行时间仍对应保存的完整 numerical run，排除后续纯渲染时间。

## 图与数据索引

![incident_energy_ranges](../results/validation/full_kirchhoff_spectral_rank_audit/figures/incident_energy_ranges.png)

原保存源的归一 bin 能量与累积能量；95/99% 角宽不由全传播区间代替。

![singular_value_decay](../results/validation/full_kirchhoff_spectral_rank_audit/figures/singular_value_decay.png)

三个案例有效/全保存谱奇异值衰减，以及强高度完整传播数学对照；机器精度以下不作物理解读。

![rank_dimensionless_parameters](../results/validation/full_kirchhoff_spectral_rank_audit/figures/rank_dimensionless_parameters.png)

15 个案例全保存谱加权 1e-3 rank 对各无量纲量；点重叠表示离散 rank 相同，颜色为 sigma。

![operator_reconstruction](../results/validation/full_kirchhoff_spectral_rank_audit/figures/operator_reconstruction.png)

真实完整 FK 运算误差随 rank 下降；矩阵与场使用各自明确的相对范数。

- `spectral_rank_audit.mat`: rank/grid/spectrum/奇异值/算子结果和运行统计，无旧 worktree 函数句柄。
- `rank_summary.csv`: 75 行，两种真实域×两种权重、一个数学对照，全部参数/误差/成本。
- `grid_convergence.csv`: 24 行，两个域/两种源宽/三种密度/两种权重。
- `incident_spectrum_01..15.csv`: 原始 kx bin、归一能量权重和传播标志。
- `operator_02/06/11.csv` 与对应 MAT: 三个真实曲面运算检查，包含 r=1–8 的全部原生重构场。
- `figures/`: 四组 PNG 和 vector PDF。

## 下一阶段判断与停止点

**继续开发 Spectral Fast FK prototype。** 依据是本范围内 1e-3 有效加权 rank=1–2、固定物理区间加密 rank 稳定，以及三类完整积分运算检查通过。这个判断是下一阶段研究建议，并不授权本阶段继续编码。

**INFERENCE**：入射 residual 分离可将 sum_j 展成少量可 FFT 的谱乘子与空间系数，具备 O(r N log N) 候选结构。**UNVERIFIED**：整个非局部 FK 的出射 Green 表示、散射角带宽、倏逝尾部、FFT padding/periodization/有限支撑处理，以及装配成本摊销仍需下一阶段证明。入射 theta99 不自动界定出射 theta99；K 会改变出射带宽。未经这些检查不能把当前结果称为完整 Fast FK 或全核低秩证明。

真实 PM 海面、比 sigma=.3 更宽的工程角谱、更大 kA、更高频率、三维和宽带均未验证。均匀 eta 区间在给定 A 下能覆盖高度，但不能自动验证 PM 的导数、散射和出射谱。若后续这些范围所需 rank 很高、随物理带宽/网格增长或出射重构不稳定，应明确停止全局分离路线并比较 block-Fresnel、locally stationary approximation、SSA/SSA2；本次没有这些失败证据，也没有对备选方法作优劣认证。

## 改动与文档纪律

新增一个 validation-only 入口、独立结果目录及本报告，更新 scripts/README、PROJECT_CONTEXT 与 vertical_comm_guide 的验证边界说明。未改变 production paramsV/output/results 接口、PE/WAPE、BIE、solve_full_kirchhoff_2d 或 local phase screen。没有新增物理公式实现，只审查现有谱耦合的数学等价分解。当前不涉及生产信道回归，验证为保存源/场闭合、原 FK contraction 等价性、rank 密度稳定性和真实重构误差；未运行 scalar-channel/wideband-communication 回归，也不声称它们本轮通过。
