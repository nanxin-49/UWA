# Research Plan for a Deep-to-Near-Surface Acoustic Communication Channel

> **固定dual 65点网格四新频点验收 — 2026-10-27。** 先锁五锚点固定`η*`复线性预测，再在Ubuntu R2023b按4750/7250Hz常规、4125/6562.5Hz相消压力顺序各取得条件性FE/原五配置FullFK参考；四点的`e_R/e_Δ/a_D`均严格低于事前`1e−3`门，最大分别`5.932e−5/3.013e−5/5.515e−5`。4125Hz相消指标`17.05`、`e_H=2.654e−4`，以直达标度并列报告；4750/4125/6562.5Hz预测差大于五倍FE经验敏感性。状态仅`COMPACT_DUAL_65GRID_HOLDOUT_SCREENING_PASS`，表示固定案例四个新离散频点对65点**预测输出**网格的有限筛查支持，不认证连续宽带、无混叠CIR、BER或原PM。[正式报告](deepsea3d_65grid_four_holdout_2026-10-27.md)、[原MAT/CSV/锁/SHA/审计/检查点/图](../results/validation/deepsea3d_65grid_four_20261027/README.md)。总保守MATLAB进程账`341/720 s`含两次非物理失败和图件审计；不结转余量，不自动启动时域窗/极值、附加频点或129点。Q1 BIE`STOP`、二维S4`CORE_COMPLETE`、`G1_target_PM=OPEN`、`WP2_formal=false`不变；下条2026-10-27网格设计为本验收的事前方案。

> **固定dual宽带网格设计 — 2026-10-27。** 已在不新增声场求解的条件下复核九个条件性 FE 频点、固定几何/native→DSP 时延与数值成本，优先设计 4–8 kHz **65 点、62.5 Hz** 的预测输出网格（相对周期窗 16 ms）；129 点仅为后续候选。下一轮建议先锁 4750/7250 Hz 常规保留及 4125/6562.5 Hz 相消压力点的复预测和 SHA，再申请不超过 720 s 的四点独立物理验证；本文不授权执行。四旧半千赫盲测已参与本次设计，不能再作为新盲测。[唯一新方案](deepsea_compact_dual_broadband_grid_design.md)。连续宽带/CIR/BER 未认证，Q1 BIE `STOP`、二维 S4 `CORE_COMPLETE`、`G1_target_PM=OPEN`、`WP2_formal=false` 不变；以下 2026-10-26 GO 只表示当时允许设计。

> **最新固定dual跨频盲测 — 2026-10-26。** A(4/6/8 kHz)与B(4/5/6/7/8 kHz)的固定几何因子复线性预测，在新4500/5500/6500/7500 Hz四个预锁定独立FE/FullFK条件性参考上均通过原`1e−3`反射复筛查；A/B最大误差分别`5.193e−5/3.050e−5`。5500/6500 Hz预测差高于五倍FE经验敏感性，4500/7500 Hz无法判定更细精度。7500 Hz粗Gauss门`9.892e−4<1e−3`但裕量小。标记`COMPACT_DUAL_FREQUENCY_HOLDOUT_SCREENING_PASS`，仅GO设计固定dual宽带网格，不执行更多频点或认证连续宽带/CIR。[正式报告](deepsea3d_four_frequency_holdout_2026-10-26.md)、[原结果与审计](../results/validation/deepsea3d_four_holdout_20261026/README.md)。Q1 BIE `STOP`、二维S4 `CORE_COMPLETE`、`G1_target_PM=OPEN`、`WP2_formal=false`保持；下文五频点结论是此前阶段快照。

> **最新三维目标 PM 参考方法决策 — 2026-10-26：STOP 当前 Q1 BIE 路线。** 按单次限额检查点在旧 6 kHz 平面制造问题上测试连续节点 Q1 密度及匹配局部 Duffy/Gauss 求积；最佳独立离节点边界/离面复场误差 `0.0250587454/0.0409758501`，均远超原 `1e−4` 门。基础求积 8→12 阶复场变化约 `5e−9`，新鲜线性残差约 `1e−15`；高频网格分辨与未知散射密度的无限平面外区仍未闭合。三个正式作业加版本探测保守 `53/600 s`，49 项回传 SHA 和保存复数独立审计通过。**不进入 compact-dual 或原 PM BIE；不转试另一求解器。**[本次决策与原数值](deepsea3d_bie_q1_decision_2026-10-26.md)、[证据索引](../results/validation/deepsea_bie_q1_20261026/README.md)。先前下文“高阶 BIE 优先”的建议现为历史提案，不再构成自动执行授权。原五频点 compact-dual 条件性结论、二维 S4 `CORE_COMPLETE` 不变；`G1_target_PM=OPEN`、`WP2_formal=false`。

> **Latest authoritative compact-dual decision — 2026-10-26.** 在冻结点源、C2 dual 海面与旧 4/6/8 kHz 证据不变的条件下，Ubuntu R2023b 完成事前锁定的 5/7 kHz 独立 FE 参考和随后完整 FullFK：两个新频点均 `CONDITIONAL_G1_PASS`；固定 `Q=R/B*` 复线性候选对独立 FE 的反射复误差分别 `6.766e-5/2.954e-5`，过原 `1e-3` 筛查。但两点绝对预测差与保存 FE 变体的经验敏感性同量级，不能据此声称更细精度、连续 4–8 kHz 响应、CIR 或带内衰落位置得到认证。[本轮三维报告](deepsea_dual_5k7k_independent_validation_2026-10-26.md)和[结果/审计](../results/validation/deepsea3d_5k7k_20261026/README.md)为证据。下一独立研究优先增强可信频率证据，先提高参考分辨率把握；本 Goal 停止，不自动扩频。目标原振幅 PM 参考仍未认证，`G1_target_PM=OPEN`、`WP2_formal=false`；不重开历史 FE/TFE/BIE 修补、PM 集合或通信任务。

> **Historical decision — 2026-10-15 / 当时仅 4/6/8 kHz 离散交付。** FE/TFE/BIE连续修补任务链停止；以下旧条目中的“next”“active”“高阶BIE优先”均为历史建议，不再授权执行。固定Tx(0,150.5,1500)/Rx(0,0,5)、dual α2.5/C2 support12的4/6/8kHz冻结参考显示反射幅度约降6.43%，总场相对平面分别−4.99/+2.38/+12.82dB；效应明显大于已观察数值敏感性。以固定镜面点海面高度分离相位的稀疏候选有探索性支持，但当时没有新的独立频点认证。已锁定5/7kHz测试在当时执行前被 Windows 内存门拒绝：物理可用639.199MiB，30%余量后447.439MiB，小于保守增量估计1520MiB；当时新增物理点0、核0s。时间外推没有证明宽带不可负担；当前完整宽带点数及资格未知，不生成认证CIR/带内衰落位置。G1_target_PM=OPEN、WP2_formal=false。见[当时中文物理决策](deepsea_dual_broadband_decision_2026-10-15.md)、[完整成本评估](deepsea_dual_broadband_cost_2026-10-15.md)、[协议](deepsea_dual_broadband_protocol_2026-10-15.md)。下文原公式、WP0–WP6、门槛及历史结果完整保留。

> **Canonical plan / 唯一正式计划：**本文件 `reports/deep_sea_channel_research_plan.md` 是后续研究唯一正式依据。v3 内容于 2026-10-14 核对后合并；[原计划完整备份](deep_sea_channel_research_plan_backup_2026-10-14.md)和[原始 v3 版本记录](deep_sea_channel_research_plan_v3_2026-10-14.md)均只读留存，不作为第二份正式计划。配套[双轨执行记录](codex_goal_two_track_2026-10-14.md)受本计划及用户指令约束。见[合并核对报告](deep_sea_channel_plan_merge_audit_2026-10-14.md)。本轮始终保持 `G1 OPEN`、`WP2=false`；冻结求解器和历史结果不改。

Original plan: 2026-10-09 (Asia/Shanghai). Execution-policy revision: 2026-10-14 (Asia/Shanghai). Evidence-ladder and Goal-mode execution revision (v3): 2026-10-14.

Status: **two-track, evidence-gated research plan**, not an approved engineering operating envelope. G0/WP0 analytic controls are complete; finite compact dual-surface cases have **conditional, case-specific G1 passes**. The original-amplitude finite-band PM target remains **INCONCLUSIVE / G1 OPEN**, and **formal WP2=false**. Production physics and frozen baselines remain unchanged. Unperformed experiments and acceptance levels remain proposals.


> **Executed checkpoint — 2026-10-14:** The single fixed right-Neumann intervention P=2I−B passed controls but worsened the original-PM linear/true PDE/DtN residuals at the unchanged 96-inner-step cap. This intervention is closed; do not read Section 10 as authorization to repeat it. Next Stream A work requires a separately budgeted matched 3D BIE/BEM window/near-singular minimum feasibility checkpoint. Stream B deterministic controls passed, but simple broadband linear interpolation and a complete practical receiver remain unaccepted. [A report](deepsea3d_preconditioner_decision_report.md), [B report](deepsea_singlelink_readiness_report.md). `G1_target_PM=OPEN`, `WP2_formal=false` remain unchanged.

> **Active follow-up — 2026-10-14:** The user authorized a bounded BIE/BEM minimum feasibility validation and delay-separated reconstruction of the analytic two-path channel. The tested polynomial preconditioner is closed as FAILED; deterministic communication controls passed, while generic linear interpolation remains unaccepted. Follow the [new protocol](deepsea_bie_reconstruction_protocol_2026-10-14.md); do not repeat the old intervention or run large original-PM BIE. `G1_target_PM=OPEN`, `WP2_formal=false`.

> **Follow-up outcome — 2026-10-14 (latest decision):** The minimum BIE experiment is closed with reference rejection: six manufactured inversions pass linear checks but fail independent boundary/field gates; at 6kHz the finest P0 trace error is 3.385%, and actual-graph manufactured PV error is 3.616e-4. Do not invest in a large original-PM solve yet. Known-path analytic two-path delay separation passes independent off-grid reconstruction (complex error ≤6.22e-12); generic total-only interpolation and dual intermediate frequencies remain unaccepted. Next single research checkpoint should test higher-order density/near quadrature on the same 6kHz flat manufactured boundary before rough scattering. [BIE report](deepsea3d_bie_minimum_report.md), [reconstruction report](deepsea_frequency_reconstruction_report.md). Earlier next-task text is historical and does not authorize rerunning either closed intervention. `G1_target_PM=OPEN`, `WP2_formal=false`.

## 执行摘要（中文；v3，供项目决策和 Codex Goal 模式使用）

**最终目标不变：**针对深海 Tx（初始 1500 m）到近水面单 Rx（初始 5 m），在 4–8 kHz 频带内建立物理可信、计算成本可核算的**宽带复信道**，最后以单链路通信性能验证。静态、均匀声速、压力释放海面是第一阶段理想化问题；动态海面、多普勒和实际海洋传播在后续阶段处理。

**为什么不重做二维：**二维 FullFK–BIE 的受控物理对照已经开展，既有通过案例也有高波数失败边界；它们应保留为回归证据。二维通常对应线源与一维粗糙剖面，三维对应点源与二维粗糙表面，Green 函数和出平面散射机制不同，二维通过**不能自动认证三维**。可做一次低优先级、条件触发的**匹配线源 + y 不变表面**二维—三维极限检查，核对归一化、法向量和符号；这只是公式/实现一致性诊断，**不是新的强制 G1 门**，也不能取代三维独立参考。

**当前证据分层（2026-10-26）：**`G0=PASS`；已验证紧支撑 dual 海面的**特定**振幅、方位、几何和 4/5/6/7/8 kHz 五个离散单点为 `G1_compact_dual=CONDITIONAL_PASS`，其中 5/7 kHz 是本轮独立新增，两个锁定重建中点仅过局部筛查；原振幅有限带 PM（固定 seed261013、C2 support12）的独立 FE/TFE/直接参考仍不足以认证，`G1_target_PM=OPEN`、`WP2_formal=false`。**不得把同一近似方法的网格收敛当作独立物理认证。**

**双轨推进：**

| 轨道 | 当前允许做什么 | 绝不能推导出什么 |
|---|---|---|
| **A：三维目标 PM 物理认证** | 复用已建矩阵自由 `t=1` 展平 Helmholtz 直接求解器；只做**一次预声明且有资源上限**的预条件器改进，随后检查真实 PDE/DtN/Dirichlet 残差；失败则转为**匹配三维 BIE/BEM 最小原型的技术和成本评估** | 未完成独立参考、网格/纵向/计算域及支撑解释之前，不得宣称 FullFK–PM 物理准确，不能放行正式 WP2 |
| **B：通信基础设施准备** | 用解析平面与已通过的 dual 案例验证 native/DSP 复相位、频率网格、时延/CIR、AWGN、已知信道与解析两路径通信，以及数据证据标签 | 不得据此拟合 PM 统计模型、宣称 PM BER/PER 或训练 PM 路径/高斯发生器 |

**执行顺序：**A 的一次预条件器检查点与 B 的确定性通信准备可独立验收。A 线性求解若改善，仍需真实物理残差和独立离散收敛；A 若超预算/失败，**停止反复调参**并给出 BIE 方案与成本边界。两个轨道都不能自动改变 `G1_target_PM` 和 `WP2_formal`。

**海面解释不可省略：**128 m 父海面、PM 谱截止、固定随机系数、C2 support12、声学孔径、求解器盒长及求积网格是不同层次。support12→support24 **改变物理海面**，不等于数值网格加密。有限支撑样本的条件性通过，也不等于无限随机 PM 海况通过。

**下一次 Codex 长任务**以本计划第 10 节的 A/B 决策为准，使用配套 `codex_goal_two_track_2026-10-14.md`。它是可连续执行的研究任务，不是授权无限加大 Krylov、海面样本或计算域。

---

## Decision update — 2026-10-14: parallel progress without unvalidated PM inference

**Research destination is unchanged:** a credible, efficient complex broadband model of a *single* deep-to-near-surface acoustic communication link, with communication-relevant validation and total-cost accounting. The two work streams may proceed in parallel, but their evidence and authorization must not be mixed.

| Stream | Authorized near-term work | Gate / prohibited inference |
|---|---|---|
| **A — target-PM physical certification (WP1/G1)** | **One budget-capped preconditioner improvement** for the already-built matrix-free, non-Taylor, `t=1` mapped Helmholtz direct solver; if it fails, evaluate a small *matched* 3D BIE/BEM feasibility alternative | PM physical validity remains OPEN. No PM FullFK accuracy claim, new PM ensemble, or WP2 release without independent-reference and surface-support evidence |
| **B — non-inferential single-link readiness** | Validate native→DSP and phase/delay/CIR metadata, deterministic flat and *certified compact dual* single-link controls, frequency/discretization guards, reproducible interfaces, and AWGN/known-channel/two-path communication sanity checks. A dimensionally matched 2D→3D limit check is an **optional, low-priority diagnostic**, not a prerequisite. | This is **WP2-preparation, not WP2**. Do not fit PM Gaussian/path statistics, publish PM performance claims, expand PM seeds, or optimize FastFK on unvalidated PM responses |

**Recorded evidence (not a new acceptance):** the `t=1` direct-solve pilot already exists. Plane/dual/quarter-amplitude controls passed; original-PM restart 8/12 runs at their recorded 96-inner-step cap were unresolved, and the best saved true PDE/DtN and true-boundary residuals did not meet the original checks. Do **not** re-run that pilot as a new prototype task. One *bounded* preconditioner intervention is a decision checkpoint, not a commitment to an indefinite solver program.

**Stream A decision rule:** predeclare solver/test configuration and runtime, iteration, memory and residual ceilings from recorded costs. If a preconditioner improves linear solves, require true physical PDE, transparent-boundary (DtN), real-surface Dirichlet, independent horizontal/vertical resolution, and computational-domain convergence **before** any matched FullFK–PM comparison. GMRES residual or receiver pressure stability alone is insufficient. If the capped intervention fails, stop iterative/grid escalation and issue a costed decision for a *matched* 3D BIE/BEM minimum pilot rather than silently cycling through Taylor order, padding, or new preconditioners.

**Surface-definition constraint:** the frozen candidate is a **finite-band PM parent realization**, not a certified unbounded PM ocean. Preserve the parent spectrum/coefficient law and seed, and keep spectral cutoff, C2 support embedding, acoustic integration aperture, solver-box extent, and quadrature resolution as **distinct factors**. Enlarging support12 to support24 alters the modeled physical sea; it is not ordinary numerical convergence. A reliable reference plus a justified support/noncompact treatment is a separate necessary part of PM G1 evidence.

**Authorization ledger:** `G0=PASS`; `G1_compact_dual=CONDITIONAL_PASS` only for tested cases; `G1_target_PM=OPEN`; `WP2_formal=false`; `WP3_PM_training=false`; `WP4_baseline_controls=allowed` with analytic/accepted-fixture labels only. Allow a future **explicitly narrowed** physical pilot only after its domain and reference are independently accepted; an infrastructure test or one diagnostic PM calculation never counts as WP2. Keep physical-model error, numerical/reference uncertainty, algorithm fidelity and statistical/communication adequacy in separate records.

### Archived execution evidence — 2026-10-12 to 2026-10-13

Stage-specific suggestions below are retained as an audit trail; when they conflict with the 2026-10-14 decision and the newest completed results, **this decision record supersedes the old proposed next tasks**. Older `INCONCLUSIVE` notes remain accurate for their then-unfinished stage, not for subsequently completed discrete cases.

2026-10-13 matrix-free t=1 mapped Helmholtz minimum pilot completed:
[report](deepsea3d_mapped_direct_report.md), [protocol](deepsea3d_mapped_direct_protocol.md).
Exact source-on-graph Dirichlet, unchanged dealiased operators and outgoing DtN, L0-preconditioned GMRES;
no Taylor expansion. Plane/dual/quarter controls pass, dual increment error1.0935e-5 to frozen reference.
Original PM restart8/12 each capped96 inner steps remain linearly unresolved; best true PDE/DtN
1.329e-3/1.273e-3 and true boundary1.04033e-4 fail original gates. Pressure stability is insufficient.
This improves Taylor32 residuals substantially but does not certify physical error. Stop iteration/grid expansion;
next bounded work should address preconditioning before formal convergence checks, with matched3D BIE fallback.
PM INCONCLUSIVE, G1 OPEN/WP2=false, no support24/FullFK–PM. Kernel wall464.135s, frozen312 files unchanged.

2026-10-13 bounded TFE dealiasing diagnosis completed:
[report](deepsea3d_tfe_dealiasing_report.md), [protocol](deepsea3d_tfe_dealiasing_protocol.md).
Cubic-exact2n padding and matched symmetric Fourier space, independent3n/4n true geometry;
dual/quarter local controls pass. Original PM N32 last coefficient norm improves2.363x,
but true PDE/DtN worsen .410575→.475266/.287852→.326480.99.60% of residual energy is retained.
Dealiasing does not rescue the current Taylor sum; do not stack further orders. Next preferred
bounded reference is a direct matrix-free t=1 mapped Helmholtz solve, with matched3D BIE fallback.
PM INCONCLUSIVE, G1 OPEN, WP2=false; no support24/FullFK–PM. Frozen284 files unchanged.
Kernel wall486.343s includes initial invalid comparator and saved-field independent probe.

2026-10-13 independent 3D TFE minimum prototype completed:
[report](deepsea3d_tfe_prototype_report.md), [protocol](deepsea3d_tfe_prototype_protocol.md).
Mapped Helmholtz with exact outgoing Fourier DtN, independently checked against manufactured modes
and the certified dual case: all dual numerical gates pass, fine increment differs6.308e-7 from frozen FE.
The same frozen PM seed261013/support12 fails at N40; h.125 still has physical PDE/DtN residuals .120/.173.
Quarter-amplitude local diagnostics pass but are not the target or a formally converged reference.
No support24 or FullFK–PM comparison; reference INCONCLUSIVE, full G1 OPEN, WP2=false.
Next bounded reference work may address dealiasing/TFE truncation or a matched3D BIE prototype.
Sixteen rough solves plus one plane control, saved kernel wall278.451s; frozen263 files unchanged.

2026-10-13 fixed-realization 3D PM feasibility completed:
[report](deepsea3d_pm_feasibility_report.md), [protocol](deepsea3d_pm_feasibility_protocol.md).
One isotropic finite-band Hm0=.5m parent (128m, kmin2pi/128, kmax1rad/m), seed261013,
6kHz, Tx1500/Rx5, q=.10; fixed alpha determines finite-band Tp3.638671s. No realized-RMS rescaling.
Nine bounded ordinary-FE probes fail order/source/true-boundary gates at N20; receiver mesh/domain
stability does not make the reference credible. Exact finite-mode exponential trace confirms failure.
Changing C2 support changes the physical sea; cutoff/support pressure differences remain unqualified.
No FullFK comparison or new conditional G1. Single-case INCONCLUSIVE, full G1 OPEN, WP2=false.
Next need matched TFE/3D BIE and justified noncompact support, then original physical/numerical gates;
no automatic WP2 release from the dual/frequency results below. Frozen solvers/production/data unchanged.

2026-10-13 fixed-y 4/8kHz conditional physics checks completed:
[report](deepsea3d_g1_frequency_report.md), [protocol](deepsea3d_g1_frequency_protocol.md).
Geometry/alpha2.5/dual/C2 frozen, 6kHz reused. Independent references require N16 at4kHz,
N20 at8kHz; mesh/domain/true-boundary/M32-to40 and cancellation controls pass.
Initial8kHz Gauss h.05 failed1e-3; one documented h.025 refinement passed without enlarging FE.
All three discrete points CONDITIONAL_G1_PASS. Reflected/increment complex errors decrease
at the primary settings, but order-sensitivity ranges overlap6-to8kHz: full complex-error
monotonicity is unresolved. Phase-error decrease and increment magnitude-error increase6-to8
are separable in observed variants. This does not certify a continuous band or generic trend.
Full G1 target OPEN, WP2=false; frozen solvers, production and historical data unchanged.

2026-10-13 fixed-sea q=.10 azimuth comparison completed:
[report](deepsea3d_g1_azimuth_report.md), [protocol](deepsea3d_g1_azimuth_protocol.md).
Tx changed from(150.5,0,1500) to(0,150.5,1500), Rx unchanged; no sea rotation/translation.
Flat backgrounds match; y-direction independent FE and FK numerical guards passed.
CONDITIONAL_G1_PASS for y: reflected/increment errors2.843e-4/1.375e-4.
Independent reference directional complex changes5.111%/2.467% are distinguished from
model-error reductions1.920%/1.762% and observed numerical/reference variations.
Empirical model-error ranges remain disjoint after saved lower-order pressure sensitivity.
No universal azimuth trend follows. Full target G1 OPEN, WP2=false; frozen data/solvers unchanged.

2026-10-13 q=.10 geometry-only continuation completed:
[report](deepsea3d_g1_q10_report.md), [protocol](deepsea3d_g1_q10_protocol.md).
Tx=(150.5,0,1500), Rx=(0,0,5), fixed alpha2.5/6kHz/dual/C2 and surface center(0,0).
WP0 plane controls were reused without reusing old moved-surface rough cases.
New N12/N16 precheck required N20; independent mesh/domain/true-boundary gates passed.
FullFK numerical guards passed, reflected/increment native errors2.899e-4/1.399e-4:
CONDITIONAL_G1_PASS only at this discrete geometry. Errors decreased2.03%/2.94% versus
frozen q0, with disjoint empirical numerical-variant ranges; total relative error slightly increased.
The complete G1 target domain stays OPEN, WP2=false; no amplitude scan/PM/FastFK/SSP/dynamics.

2026-10-13 alpha2.5 domain and physical check completed:
[report](deepsea3d_alpha25_domain_report.md), [protocol](deepsea3d_alpha25_domain_protocol.md).
Frozen L24-to32/h=.0625/N20 domain change3.242e-7 and true half-cell residual2.977e-7
pass the original1e-4 gates; L48 was not needed. Only after this, matched FullFK
mesh/Gauss/aperture controls passed; reflected/increment native complex errors
2.959e-4/1.442e-4 give CONDITIONAL_G1_PASS at alpha2.5 for this discrete case.
This supersedes the alpha2.5 physical INCONCLUSIVE status in the earlier stage records below.
The complete G1 target domain remains OPEN, WP2=false; no larger-amplitude/PM/FastFK/WP2 work.
Frozen solvers, production and historical acceptance data remain unchanged.

2026-10-13 alpha2.5 fine-reference check completed: [report](deepsea3d_alpha25_fine_report.md).
L24/h=.0625 N16/N20 passed order/grid/true-boundary controls; no observed FE instability.
Next a bounded domain check, then matched FullFK comparison if the reference remains stable.
No domain extension or FK solve was run in this stage; alpha2.5 physical G1 is INCONCLUSIVE,
the complete domain remains OPEN and WP2 is not authorized.

2026-10-13 alpha2.5 bounded precheck completed: [report](deepsea3d_alpha25_preflight_report.md).
N12 is insufficient; N16-to20 change and true half-cell boundary residual pass the existing
coarse precheck. Ordinary FE has no observed catastrophic instability on this configuration.
Next only a controlled finer-grid reference check, then domain/FK work if stable.
No formal alpha2.5 G1 pass or WP2 authorization follows; no full validation was started.

2026-10-13 amplitude-only continuation: [report](deepsea3d_g1_amplitude_report.md).
Frozen alpha1 reused; new .5/1.5/2 each passed independent order/mesh/domain/true-boundary
controls and conditional pointwise G1 at unchanged geometry/spectrum/phase/support.
No credible physical failure or FE instability was found. This certifies four sampled points,
not a continuous interval or the complete G1 domain; WP2 remains blocked by coverage gaps.
If pursuing larger amplitude, first probe alpha2.5 order stability at bounded cost.

Latest continuation: [finite-roughness single-case report](deepsea3d_g1_finite_report.md)
and [protocol](deepsea3d_g1_finite_protocol.md). An independent order16 Dirichlet FE reference
with grid/domain/order and true-boundary residual controls conditionally passes one
Tx1500/Rx5, q0, 6 kHz compact dual-direction surface, amplitudes (.02,.01) m.
FullFK reflected/increment errors are 1.36e-4/1.06e-4; finite-amplitude effects are substantial.
This closes the first-order-only gap for this case, not the full target domain G1.
WP2 remains unauthorized. Next use one controlled physical-domain extension; do not require
a dense BIE or repeat this completed single case as a blanket prerequisite.

[Chinese stage report](deepsea3d_wp0_wp1_report.md) and [frozen protocol](deepsea3d_wp0_wp1_protocol.md) supersede the initial execution order where specified here.

- The 2D PM physical-domain audit is complete (15 controlled inputs, four guarded cases). Do not repeat it as a blanket prerequisite. Its conditional boundary is not a universal cutoff or 3D acceptance criterion.
- Newer 2D single-Rx FastIn/composite prototypes are conditionally frozen only at original Rx(0,3), with cancellation-position failures. They add algorithm fidelity, not BIE physics; no acceleration/rank work is required for this 3D stage.
- WP0/G0 is PASS for Tx1500/Rx5, q=0/.10, unit Green point source, image fields, shifted planes, rotations and native/DSP algebra.
- Minimum WP1 now has the full two-term single-Rx 3D Kirchhoff integral with analytic gradients and exact infinite-flat background subtraction for compact C2 perturbations. Numerical gates pass at 6 kHz, then 4/8 kHz. Millimetric weak/dual-direction controls do not represent the future Hm0=.5/1 m PM pilot.
- One independent first-order Dirichlet shape-reference formulation is implemented: 12 amplitude/geometry comparisons, Gauss grid/support guards and Fourier cross-check. This is a weak-limit check, not finite-roughness BIE certification. Full G1 remains OPEN; these results do not authorize WP2/PM ensembles.
- The earlier next-step finite-amplitude reference is now completed conditionally for the case above. For harder cases test FE stability anew, or use TFE/matched half-space BIE. An open BEM patch or a 1% FK–first-order discrepancy still does not establish full G1.
- Section 9 numerical 1e-3 and physical screening 1% proposals remain distinct from the historical 2D physical 1e-3 gate. No tolerance was relaxed to pass this campaign.

## 1. Recommended direction and intended contribution

**Develop a receiver-focused, physically informed model of the complex broadband response of one deep-to-near-surface link. Maintain independent physical certification for the claimed sea-state domain while progressing deterministic single-link and communication infrastructure on already accepted controls; compare statistical and direct-sampling routes only when trustworthy target-domain data are available.**

The existing physical candidate is an analytic homogeneous-medium incident field plus a direct single-receiver 3D rough-surface Kirchhoff integral. Analytic controls and several compact dual-surface comparisons have passed, but this does **not** validate the finite-band, original-amplitude PM target. It remains a conditional physical evaluator, not an unqualified truth source or a production-PE replacement.

After an appropriate physical evidence gate, compare direct physical sampling (P0), a conditional joint-frequency Gaussian process (S1), and a random-path-plus-residual model (S2). S2 is a mechanism-oriented **hypothesis, not a preselected winner**; neither it nor S1 is authorized for PM fitting today. Direct evaluation is a valid final choice when its end-to-end cost and accuracy suffice.

The scientific contribution should be a **validated domain of applicability, a justified stochastic representation, and its communication-relevant accuracy/cost tradeoff**. A useful negative finding, such as a physically meaningful regime where few-path or Gaussian models fail, is part of that contribution. Publication novelty must be established against primary literature before claiming a new model; hybrid physical/statistical channel modeling already has precedents. [E4]

The main research questions are:

1. Which rough-surface mechanisms control the single-link broadband response in the proposed deep, near-vertical geometry?
2. Can their effects be represented by a small set of jointly distributed path variables, or is a richer joint-frequency process necessary?
3. Over which geometry, surface spectrum, and receiver-depth ranges does the selected physical/statistical approximation meet the chosen communication tolerance?
4. When do offline construction and condition updates pay for themselves relative to direct physical sampling?
5. Which state variables must be retained to extend the static model to consistent time evolution and Doppler?

## 2. Requirements, assumptions, and scope boundaries

| Category | Current specification | Status |
|---|---|---|
| Application | Deep transmitter to a cable-deployed receiver near a vessel or surface platform; one Tx/Rx link | User requirement |
| Depth | Tx deeper than 1000 m; Rx 3–10 m | User requirement; maximum Tx depth and bottom clearance unknown |
| Geometry | Near vertical, with horizontal displacement permitted | User requirement; operational displacement distribution unknown |
| Spatial physics | 3D propagation and a rough surface with two horizontal coordinates | User requirement; reductions need demonstrated error bounds |
| First stage | Static surface, uniform sound speed | User requirement |
| Frequency | Start with 4–8 kHz | User preference; changing band requires a scientific reason |
| Later capability | Evolving sea, time correlation, and Doppler | User requirement |
| Sea model | PM is a reasonable initial candidate | User preference; directionality, cutoff, and later alternatives open |
| Source | Transducer model and directivity unknown | Explicit uncertainty |
| Deployment | Cable deployment confirmed; hull clearance unspecified | Explicit uncertainty |
| Computing | Limited offline work allowed; count construction, validation, updates, and online sampling | User requirement; numerical budgets unknown |
| Acceptance | Communication performance plus necessary physical/statistical checks | User requirement; proposed tolerances appear in Section 9 |

The initial physical study is a uniform, lossless water half-space with a pressure-release surface, fixed endpoints, and an explicitly normalized point source. It isolates free-surface effects. It is not a calibrated shipboard link budget. Constant sound speed does not imply zero absorption in a real ocean; absorption is omitted only in this first controlled model.

Before making engineering range, transmit-power, or outage claims, add or bound frequency-dependent absorption, source/receiver responses, ambient and platform noise, bottom-related paths, and sound-speed structure as appropriate. Tx depth is not water depth. A source near the bottom may require a bottom-reflected contribution, so bottom clearance is a high-priority future input.

Hull scattering, wakes, bubbles, and motion are not declared negligible. Initially excluding them defines the idealized study domain; later deployment sensitivity determines whether the idealized result transfers. Unknown engineering information should not prevent analytic controls, but it limits the final claim.

## 3. What existing work changes in this plan

| Existing result | Appropriate use now | Consequence |
|---|---|---|
| PE/WAPE, receiver adjoint, cross-frequency statistics, and MPSK tools | Reuse conventions, interfaces, and independently verified numerical components | Do not assume that old local-screen statistics are statistics of a new rough-boundary solver |
| Completed FullFK–BIE deterministic/PM-domain comparisons | Reuse 2D regressions, conditional passes and guarded failures | Do not repeat the completed audit; no automatic physical transfer to 3D |
| Fixed-PM 4–8 kHz FastFK study, ranks 3/12 | Algorithm-fidelity and phase-convention regression | No reason to retain those ranks, sample counts, or surface cutoff in the new model |
| Single-receiver timing reversal | Match computation to the requested output | Suspend dense-field FastFK optimization as the main research activity |
| 2026-10-12 single-Rx conditional freeze | Reuse algorithm contracts and cancellation diagnostics | Not additional BIE physics or a 3D speed forecast |
| Finite-band CIR agreement | Reuse native-to-DSP conversion tests | Re-establish delay support and frequency-grid convergence |
| 3D WP0 and compact dual-surface finite-amplitude controls, 4/6/8 kHz and two q=.10 azimuths | Reuse frozen image/phase checks and matched high-order Dirichlet physical-reference cases | Conditional G1 **only for named compact surfaces, geometries, and frequencies**, not PM or a continuous band |
| Fixed Hm0=.5 m PM feasibility / FE–TFE–dealiasing / mapped t=1 direct pilot | Preserve failures, residual budgets, and exact surface/solver definitions | Original PM remains INCONCLUSIVE: Taylor/TFE residuals failed; direct-pilot linear and true-residual checks also failed; one bounded preconditioner decision, then matched BIE feasibility if needed |
| External AI suggestions | Retain useful questions and candidate methods | Recompute geometry, roughness moments, and budgets; do not inherit unsupported thresholds |

The measured fixed 4 kHz, 2D Tx100/Rx3 case gave 0.98571 ms for FullFK single-point evaluation versus 36.0479 ms for FastFK whole-line evaluation followed by extraction. Including preparation gave 111.55524 versus 214.69579 ms. FullFK was therefore approximately 36.57 and 1.925 times faster, respectively; pointwise complex difference was about 0.0177%. These are workload-specific measurements, not projected 3D timings. [L2]

In that implementation, most FullFK cost was incident-field reconstruction on the actual surface. An analytic 3D point-source field changes this preparation problem substantially, but it also changes the source model from the old 2D Gaussian fixture. Any resulting cost reduction must be reported as a new formulation/workload comparison, not an equivalent acceleration of the old fixture.

Continue to preserve the production WAPE model and all historical references. Research alternatives should be introduced through validation workflows and only integrated into public interfaces after their domain and semantics are established.

### 3.1 Evidence ladder: why completed 2D validation precedes, but cannot replace, 3D G1

The historical 2D FullFK–BIE campaign is a **completed, bounded reference study**. It checked the complete two-term Kirchhoff representation against an independent Dirichlet BIE over declared deterministic and finite-band PM cases, including a documented high-wavenumber limit where the approximation missed its 2D acceptance gate. That evidence justifies keeping the 2D formula, code and known failures as frozen regression assets. It does **not** certify the 3D point-source problem, the two-coordinate surface geometry, azimuthal scattering, or the original-amplitude target PM population. Conversely, lack of a 3D PM reference does **not** invalidate the completed 2D conclusions.

The evidence ladder is therefore:

1. **2D physical mechanism:** matched 2D FullFK–BIE successes **and** observed failures; independently evaluated, no blanket correctness claim.
2. **3D analytic/controlled reference:** flat Green/image checks, compact dual-surface FE references and numerical guards; conditional, named-case G1 passes.
3. **3D target PM independent reference:** **OPEN** until a non-Kirchhoff solver satisfies true PDE, pressure-release boundary, outgoing condition, resolution/box convergence, and a defensible finite-support/cutoff treatment. Only then run a matched FullFK–PM physical error comparison.
4. **Single-link statistical/communication inference:** formal PM mechanism tests, P0/S1/S2 comparisons and communication claims require their own physical-domain and evidence gates; independent deterministic DSP readiness may proceed before that.

**Optional cross-dimensional regression (not a new WP1 stop gate).** If it can be implemented cheaply from existing validation code, verify a truly *matched* 2D limit: a y-translation-invariant surface and matching **line-source** excitation, phase convention, normals, boundary model, radiation condition, and response normalization. The y-integrated outgoing 3D Green kernel must reproduce the 2D Hankel kernel (with justified oscillatory/infinite-aperture treatment). A 3D **point source** placed above an eta(x)-only sea is **not** automatically the old 2D Gaussian-source experiment; comparing those directly would create a false failure. First audit the analytic reduction and matched fixtures; implement only if the comparison can be made without a new large numerical solver. Report agreement or a documented deferral as a formula/implementation diagnostic, never as a target-PM physics pass. Do not let this optional check displace the capped preconditioner decision or deterministic single-link readiness.

## 4. Proposed physical and statistical architecture

### 4.1 Physical sample generation

For a source at position x_t and a surface point s, start from the native outgoing 3D Green function

\[
G_k(s,x_t)=\frac{e^{ikR_t}}{4\pi R_t},\qquad R_t=\lVert s-x_t\rVert,\quad k=2\pi f/c,
\]

with time convention exp(-i omega t). Its analytic spatial derivative is

\[
\nabla_sG_k=G_k\left(ik-\frac{1}{R_t}\right)\frac{s-x_t}{R_t}.
\]

A declared source-strength spectrum multiplies this field. A directional source later requires the corresponding angular weighting and derivatives, not an undocumented beam factor. The same source normalization must be used for the direct and reflected responses.

Compute the receiver-directed surface integral using the actual two-coordinate surface, consistent normals, surface measure, and spherical source/receiver distances. Derive the chosen 3D boundary approximation before implementation. Do not copy the old Hankel kernel or assume that changing a 2D kernel to a 3D kernel certifies the entire formulation.

At fixed quadrature, direct single-point incident construction and integration can be organized with work proportional to the number of surface samples per frequency, rather than receiver-count times surface-count. Oscillatory quadrature, accuracy-driven surface size, and reference-solver complexity remain major unknowns. The proposed first implementation should prioritize transparent convergence and streamed storage over compression.

Use a hierarchy of references:

| Level | Reference | What it can establish |
|---|---|---|
| A | Exact flat pressure-release image solution; uniformly shifted flat surface; coordinate rotations | Normalization, geometry, signs, conventions, and simple limits |
| B | Refinement of the same 3D approximate integral | Its numerical convergence only |
| C | Independently formulated 3D boundary solution on manageable matched cases, plus weak-roughness limiting checks | Physical approximation error under those matched conditions |
| D | Eventually, measured channel responses with known geometry/environment | Transfer to actual deployment |

A converged 3D Kirchhoff implementation is not Level C evidence for itself. Agreement of two Kirchhoff implementations is algorithm fidelity. Existing 2D BIE results remain useful regressions but cannot replace a three-dimensional physical check. The accepted compact dual cases provide **case-specific** Level C evidence from independently checked high-order Dirichlet references; they do not transfer to the original-amplitude PM sea. If a matched 3D PM reference is infeasible, report the remaining physical uncertainty and limit the claim; do not label the reference gap as a pass. Deterministic DSP/communication infrastructure may still be tested with separately labeled Level A/C controls.

### 4.2 Statistical output and candidate generators

All stochastic and communication quantities should use one documented DSP frequency convention. Preserve the native physical response separately. With a declared direct travel-time reference tau_d, use a response such as

\[
\widetilde H(f;\xi)=H_{\rm DSP}(f;\xi)e^{i2\pi f\tau_d},
\]

where xi contains the random sea state realization and, in later stages, motion. The direct reference must be prescribed by geometry; no realization-specific fitted alignment is permitted in physical-solver comparisons.

**Only after the relevant physical-domain gate and dataset are established**, compare the following candidate routes. Listing them is a research design, **not approval to fit statistical models to uncertified PM data**, and no route is guaranteed to win:

| Route | Representation | Reason to retain it |
|---|---|---|
| P0: direct physical sampling | Generate a surface and evaluate the single-link response | Accuracy/cost baseline and a possible final solution |
| S1: conditional joint-frequency Gaussian | Estimate mean, covariance C, and, where needed, pseudo-covariance P; sample a regularized joint model | Simple baseline; tests whether second-order structure is enough |
| S2: random paths plus residual | Joint random path delays, amplitudes, and phases, with a conditional joint-frequency residual | Preferred mechanism-based hypothesis if individual arrivals are stable and identifiable |

An illustrative S2 form is

\[
\widetilde H(f;\xi)=A_d(f)+\sum_{p=1}^{P}A_p(f;\xi)e^{-i2\pi f\Delta\tau_p(\xi)}+R(f;\xi).
\]

The amplitudes must have a constrained, smooth frequency representation; otherwise arbitrary complex amplitudes can absorb any delay and make the decomposition unidentifiable. Start with one reflected cluster, then add complexity only if held-out evidence justifies it. Do not require that every rough realization has one specular point, and do not interpret sub-resolution fitted paths as uniquely observed physical arrivals.

Estimate path variables jointly and retain dependencies between delays, amplitudes, and residual statistics. If a phase/delay alignment helps fit a conditional residual, save the transformation and sample its distribution so the final generator restores the original phase and timing. Never use such alignment to hide physical-reference disagreement.

A small ensemble reflected mean does not imply a weak reflected realization. Evaluate reflected-only coherence, reflected power, and direct-to-reflected power separately. A nonzero mean does not require nonzero centered P. Joint proper/improper modeling should follow the actual covariance definitions. [E2]

Add a non-Gaussian mixture, amplitude-dependent residual, or other richer family only if both simpler candidates fail important held-out statistics or communication criteria. Deep learning and extensive surrogate training are deferred until a measured need and a sufficient independent dataset exist.

## 5. Scenario specification and sea-surface design

### 5.1 Parameter priorities

The following are research settings, not user-approved installation limits or published offset quantiles:

| Role | Geometry and sea setting | Purpose |
|---|---|---|
| First analytic control | Tx1500 m, Rx5 m, c1500 m/s, flat surface, q=0 and 0.10 | Establish deep geometry and absolute conventions |
| First numerical roughness controls | Same geometry; compact weak sinusoid and dual-direction perturbations | Executed WP1; amplitudes .001 and (.002,.001) m, not PM |
| Target PM mechanism pilot (not yet authorized) | Hm0=0.5 and 1 m crossed with q=0 and 0.10 | **Formal WP2=false** pending matched PM physical-reference, spectrum/support and numerical-domain acceptance; not an instruction to generate 32 cases now |
| Nearby validation | Tx3000 m; Rx3 and 10 m; q up to 0.25; Hm0 up to 2 m | Test transfer one controlled extension at a time, then selected interactions |
| Later stress cases | Tx6000 m, q=0.50, Hm0=3 m | Extend only after the core domain and cost are understood |

Here q=rho/(z_t+z_r), rho is horizontal displacement, and the flat-reflection normal angle is arctan(q). q=0.10 corresponds to 5.71 degrees and rho=150.5 m at Tx1500/Rx5. q=0.25 corresponds to 14.04 degrees. The optional q=0.50 case is a wider-angle stress test, not part of an asserted near-vertical operating specification.

Published deep-vehicle deployments and descriptive sea-state sources motivating these choices are documented in the companion reassessment. They support realistic application context, not the exact proposed grid. [L4]

### 5.2 Define the spectrum before defining the mesh

Start with a fully developed PM candidate and define its spectral measure and units explicitly. PM is a model with assumptions, including wind conventions; a label alone does not specify a realization. [E1]

The surface specification must include Hm0=4 sqrt(m0), peak period, directional distribution, long-wave and short-wave limits, coefficient law, Fourier normalization, and seed. For a fixed PM formulation, height and period cannot be chosen independently without changing the model family. A useful implementation path is to select target Hm0, solve for the corresponding parameter of that declared PM family, and report its implied period; any finite-band adjustment must be explicit.

Use a normalized isotropic directional spectrum as a reproducibility control. Then compare a declared directional spectrum with its mean wave direction parallel and perpendicular to horizontal offset. The spreading function and concentration are proposed sensitivity settings until supported by a target environment; none should be called a measured sea condition.

Separate three lengths:

1. The parent surface domain, large enough to represent long waves and local slope/curvature.
2. The acoustic integration support, determined by contribution and convergence checks.
3. The quadrature spacing, determined by both acoustic oscillation and resolved surface geometry.

Long-wave modes larger than the integration patch can affect its local height and tilt. Do not remove them merely by generating a small periodic patch. Conversely, a large illumination footprint is not automatically the required integration aperture. Use a common parent realization when changing the integration window or numerical mesh.

The frozen target-PM feasibility specimen uses seed261013, an isotropic finite-band 128 m parent surface with kmin=2pi/128 and kmax=1 rad/m, Hm0=0.5 m (finite-band model specification), 6 kHz, Tx1500/Rx5 and q=.10. Its C2 support12 embedding is a *separate modeled physical boundary* and has not passed an independent finite-amplitude PM reference. A support24 case is **a changed surface**, not a finer solve of support12. Track parent spectral cutoff, realization, C2 envelope/support, acoustic integration support, FE/TFE/solver box and horizontal/vertical/quadrature grids separately; do not classify physical support or cutoff changes as pure discretization convergence.

The PM gravity-wave tail does not provide finite slope/curvature moments without a short-wave limit. Choose and document the physical high-wavenumber treatment before statistical production. Run nested spectral-band checks with common low-frequency coefficients. If omitted bands materially change receiver observables, extend the spectrum or narrow the model's claimed domain. Numerical mesh refinement must not silently add physical roughness.

Record actual height/slope moments and receiver clearance for every realization. If a surface crosses a fixed receiver, classify it as outside the wet-receiver idealization; report its occurrence instead of silently resampling and changing the ensemble.

## 6. Work packages, deliverables, and decision gates

### WP0 — Lock definitions and establish analytic controls

**Question:** Is the proposed single-link response represented consistently before roughness is introduced?

Specify coordinates, source units, native/DSP conversion, direct time reference, source/receiver response assumptions, and the surface spectrum. Obtain the flat 3D solution

\[
H_{\rm native}(f)=\frac{e^{ikL_d}}{4\pi L_d}-\frac{e^{ikL_r}}{4\pi L_r},\quad
L_d=\sqrt{\rho^2+(z_t-z_r)^2},\quad L_r=\sqrt{\rho^2+(z_t+z_r)^2}.
\]

Check direct, reflected, and total fields separately. For vertical geometry, reflected excess delay is 2 z_r/c; Rx5 m at c1500 m/s gives 6.6667 ms and a 150 Hz two-path interference period. This persists when Tx depth changes in the vertical, uniform, flat control.

Include a uniformly shifted flat reflector to show how random phase/delay can suppress the ensemble mean without creating diffuse multipath. Add horizontal coordinate rotation tests so that isotropic-source, isotropic-surface results do not depend on axis naming.

**Deliverable D0:** a versioned problem specification, independent analytic test cases, and a phase/CIR metadata contract.

**Gate G0:** analytic algebra and numerical representations agree within the proposed numerical tolerances; no convention fitting is used. Failure stops rough-surface interpretation. **Current status: PASS on the documented flat, shifted-plane, rotation and native/DSP controls**; repeat only when conventions or source geometry change materially.

### WP1 — Build a small receiver-focused 3D physical reference campaign

**Question:** Can the selected rough-surface approximation produce a defensible response at an acceptable cost?

Reuse the completed receiver-focused 3D candidate and its controls rather than restarting WP1. Preserve the completed 2D BIE physical audits as regressions, and treat the matched line-source 2D→3D limit of Section 3.1 as optional, **not** as a prerequisite to continue the 3D PM reference. Compact perturbations use exact infinite-flat background plus rough/flat full-density differences, with raw finite-aperture flat errors separately recorded. The validated compact dual family includes controlled finite amplitudes, q=0 and q=.10 x/y orientations, and 4/6/8 kHz **discrete** comparisons to independent high-order references. These passes do not establish population, continuous-band or PM reliability.

**Current gap:** the fixed original-amplitude PM realization has *no credible independent finite-amplitude reference*. Ordinary FE and Taylor/TFE checks failed; Fourier dealiasing did not restore the TFE physical residual. The mapped non-Taylor `t=1` Helmholtz pilot passed control cases but its capped PM GMRES/true PDE, DtN and real-boundary checks remained unresolved. No matched FullFK–PM comparison was accepted. Do not reproduce completed dual or failed PM experiments as new progress.

One **budgeted preconditioner improvement** on the existing mapped direct solver is now the first Stream A decision. Freeze the specimen and reference operator; diagnose operator/preconditioner/GMRES consistency and record true residuals independently of Krylov stopping. A linear solve passing is necessary but not sufficient: then check real-surface boundary handling, horizontal and vertical mesh refinement, transparent boundary, box/support sensitivity and the original physical-reference conditions. If it cannot pass within the predeclared cap, stop this solver route and evaluate a minimum *matched* 3D BIE/BEM pilot with feasibility and memory budgeting before implementation. Never use a free-space closed-body BEM result as an unexamined reference for an infinite rough half-space.

The original reference-campaign size (up to 12 cases) is **not presently authorized for target PM**. A later bounded multi-case campaign requires an accepted reference and separately justified spectral cutoff, finite embedding/noncompact support and integration aperture. Numerical mesh refinement must hold the *same* physical sea fixed.

**Deliverable D1:** convergence report, measured cost breakdown, and a case-by-case applicability table including failures and reference limitations.

**Gate G1 — explicitly scoped, not a single project-wide binary claim:**

- **Already accepted:** `G1_compact_dual=CONDITIONAL_PASS` only for the archived, independently referenced compact examples. These may support deterministic control and method studies, not PM population claims.
- **Still open:** `G1_target_PM=OPEN` and `WP2_formal=false`. A valid direct mapped Helmholtz or matched BIE result first needs independent true PDE/Dirichlet/DtN checks, mesh/vertical/box convergence and an explicit support/cutoff interpretation; **then** compare FullFK versus reference on exactly matched PM surfaces. A single conditional success would not automatically certify the original four-condition PM pilot.
- **Numerical/reference failure:** diagnose the reference within a fixed budget; keep `INCONCLUSIVE` if it cannot be made reliable. Do not label failed numerical checks as FullFK physical failure.
- **Converged physical disagreement:** use a justified alternative (e.g. SSA/two-scale) or explicitly restrict the accepted sea-state domain; increasing FastFK rank cannot repair physical approximation error.
- **Adequate physical agreement:** authorize only a *declared and independently supported* physical subdomain. Expand to the formal WP2 pilot only when coverage of its planned frequency, geometry and surface family has been justified and recorded in a new approval decision.

### Cross-track readiness — allowed while target-PM G1 is open (not WP2)

Allowed independent deliverables: native/DSP conversion and shared frequency/delay/CIR metadata; deterministic flat-image and independently accepted dual-surface regression signals; numerically checked phase/CIR windows and frequency grids; validation infrastructure, dataset schema/provenance/flags; AWGN, known-channel and flat two-path QPSK receiver controls; and resource accounting for those **named workloads**. An analytically matched 2D-line-source to y-invariant 3D regression may be added **only if low-cost and non-blocking**, with `cross_dim_formula_diagnostic` provenance; a point-source-versus-old-Gaussian comparison is invalid. This check is not a new G1 acceptance gate. Clearly label whether a result is `analytic_control`, `compact_dual_conditional`, or `PM_unvalidated_diagnostic` and never mix these tiers in physical/statistical acceptance.

Do **not** collect an enlarged original-amplitude PM ensemble, fit PM S1/S2 distributions, calibrate residuals from PM numerical artifacts, tune FastFK on unvalidated PM, infer PM outages/BER, or treat a deterministic signal-processing test as PM communication validation. A diagnostic PM output may be retained with a non-acceptance label but must not enter training, model-selection, held-out testing or reported physical statistics. Stream B deliverables have their own readiness checklist and do not change G1, G2, G3, G4 or `WP2_formal`.

### WP2 — Run a small paired mechanism experiment (formal; currently blocked)

**Question:** Does roughness primarily randomize a small number of path variables, or does it require a richer frequency process?

**Do not execute this package while `WP2_formal=false`.** The four conditions and counts below remain a contingent experiment design, not permission to generate target-PM training/validation data. When a later G1 decision explicitly releases a defined sea-state domain, use the four pilot conditions in Section 5 only if *each* is supported, or formally narrow the pilot and its claims. The originally proposed crossed conditions are Hm0=0.5/1 m and q=0/0.10 at Tx1500/Rx5. Start with eight independent base coefficient seeds only after release; pair compatible conditions with common standardized coefficients, preserving actual spectra so that coupling is not confused with identical physical seas.

This gives 32 condition-realization responses. At the initial 4/6/8 kHz screening stage it represents 96 scalar-frequency evaluations per physical candidate, excluding controls, refinements, and independent-reference checks. These counts are budgeting definitions, not a statistically powered final sample. Use the same realizations across frequency; never redraw a surface independently at each frequency.

After the sparse checks, compute a genuinely converged broadband response for selected realizations and then the pilot set if cost permits. A first trial grid is 257 frequencies over 4–8 kHz, with spacing 15.625 Hz and periodic transform window 64 ms. Compare with 513 frequencies at fixed band for a 128 ms window and evaluate intermediate frequencies to test interpolation. Neither grid is accepted a priori. Fixed-band refinement increases the unaliased window, while nominal delay resolution remains about 1/B=0.25 ms; zero-padding only interpolates.

Use geometric path-delay estimates and observed tail energy to enlarge the window when needed. A small energy tail can still affect deep fades, so also test its communication effect. Save transform endpoints, bin ordering, taper, amplitude normalization, reference time, and carrier definition. Do not infer an absent reflected component from its falling outside a periodic CIR window.

Extract reflected power, direct/reflected ratio, reflected-only coherent fraction, phase and delay variation, frequency covariance and pseudo-covariance, band-limited delay structure, and fades. Interpret peak count with the finite bandwidth in mind.

**Deliverable D2 (deferred):** paired effect plots and a comparison of P0/S1/S2 feasibility with uncertainty and counterexamples, limited to the newly authorized physical domain. Readiness controls from Stream B are not D2 physical statistics.

**Gate G2 (not evaluated):** proceed only when the chosen observables and frequency representation are numerically stable *and the physical inputs are accepted for their claimed domain*. If path extraction is unstable or residual structure remains strong, favor direct sampling, S1, or a richer joint model instead of forcing a few-path interpretation.

### WP3 — Fit and validate the conditional static generator

**Question:** Which representation generalizes to unseen surfaces and nearby conditions at the lowest useful cost?

**WP3 is not authorized while formal WP2 and G2 remain blocked.** Do not fit statistical PM generators to unvalidated physical outputs. Once an authorized physical pilot passes G2, use it for design and tuning only. A provisional first development batch is 32 new independent base seeds, followed by 32 disjoint validation seeds across the selected core conditions. These are batch sizes for cost and precision assessment, not a guarantee of adequate covariance or tail estimation. Reserve a separate final test set after selecting the model and tolerances.

At hundreds of frequencies, small ensembles cannot identify an unrestricted full covariance. Use shrinkage, a reduced basis, or a structured kernel selected on development data. Verify positive semidefiniteness, pseudo-covariance compatibility, and sensitivity to basis dimension. Do not let low sample rank masquerade as physical low dimensionality.

For S2, regularize path count and spectral amplitude smoothness and preserve joint path/residual dependencies. For S1, assess proper and improper variants; retain P only where its estimated structure is stable and useful. Compare amplitude distributions, joint-frequency functionals, delay statistics, fades, and held-out predictive scores. Paired complex errors apply to methods supplied the same realization; an independently sampled statistical generator should be compared at the distribution/functional level, not by arbitrary sample pairing.

Test interpolation in operating conditions separately from independent-seed generalization. Hold out at least one interior condition, such as an intermediate offset, and one boundary condition. Do not claim extrapolation to Tx6000 m, strong seas, or an unseen directional family from interpolation evidence.

**Deliverable D3:** frozen generator, model card, training/validation manifest, uncertainty report, and model-selection record.

**Gate G3:** the simplest model satisfying the physical-domain restrictions and held-out statistical/communication criteria wins. If a richer model does not show a stable advantage, retain the simpler one. If all statistical candidates fail, keep direct physical sampling and diagnose the missing structure before expanding model complexity.

### WP4 — Validate communication performance and total cost

**Question:** Does replacing physical channel sampling change relevant communication conclusions within the declared tolerance?

Use single-carrier QPSK as the primary static benchmark because the repository already has MPSK infrastructure. Propose fc=6 kHz, symbol rate 2000 symbols/s, root-raised-cosine roll-off 0.25, and uncoded raw rate 4000 bit/s. The nominal occupied band is about 4.75–7.25 kHz, within the preferred band. These are research choices, not specifications for the unknown transducer. Add a wider-band waveform sensitivity test only after the first benchmark is trustworthy.

Establish two receiver modes: a perfect-CSI/timing control and a practical pilot-based channel estimator with a documented equalizer and synchronization procedure. Derive sample-rate, pulse-filter, noise-variance, Es/N0 and Eb/N0 mappings explicitly. **The AWGN, analytic flat two-path, and noiseless known-channel tests may proceed now as Stream B controls**, without any target-PM simulated dataset or PM performance interpretation. The old processing-related BER floor must be resolved or the benchmark rebuilt independently; it must not serve as the channel acceptance tolerance.

For a concrete starting configuration, propose 1024 payload QPSK symbols per packet, a 128-symbol known pilot, two complex samples per symbol, and a 10-symbol pulse-filter span. Start with a fractionally spaced linear MMSE equalizer; determine and freeze its length from the accepted channel-memory and pulse-support study, then check a longer filter. Use a regularized pilot-based impulse-response estimate in the practical receiver. Set guard duration from the converged channel support rather than a historical CIR window. If linear equalization is an inadequate benchmark under relevant spectral nulls, evaluate a documented DFE extension and keep its settings identical across generator comparisons. Report net rate including pilot and guard overhead. These packet choices belong to the explicitly static first study; the later dynamic study must test variation within a packet rather than assume it remains frozen.

**Deferred WP4 acceptance experiments (not Stream B controls):**

1. **Generator fidelity:** only with an accepted physical ensemble, and using the same receiver algorithm and noise policy, compare performance distributions under physical and generated channel ensembles.
2. **Receiver-model mismatch:** only on held-out *accepted* physical channels, transmitted data, and pilot/noise realizations, evaluate receiver designs or priors; keep a perfect-CSI control to separate estimation losses.

No PM BER/PER, fading, outages, or generator-equivalence claim is authorized by an AWGN, image, or compact-dual sanity check. Stream B may verify software interfaces, phase reference and waveform algebra but must label such results **pipeline readiness**, not physical or communication acceptance for the PM target.

For fade and outage tests, retain channel gain under a fixed transmitted-energy/noise convention. A per-realization normalization experiment may isolate frequency-selective distortion, but must be labeled secondary because it removes part of the fading effect.

Measure surface generation, source construction, scattering, broadband assembly, fitting, validation, condition updates, and sample generation separately and end to end. Include first-call startup and warmed execution as separate categories. No claimed single-link speedup may use an unnecessary dense-field baseline.

**Deliverable D4:** BER/SNR and packet-error comparisons with uncertainty, fixed-energy fade/outage comparisons, and complete cost/break-even curves.

**Gate G4:** acceptance requires the proposed communication tolerances with adequate Monte Carlo precision, plus no unresolved physical/statistical failure in the claimed domain. Speed alone cannot compensate for incorrect tails or communication behavior.

### WP5 — Establish the domain and deployment sensitivity

Expand selected factors after the core model is stable: Tx3000 m, Rx3/10 m, q up to 0.25, Hm0 up to 2 m, and wave direction relative to the offset. Prioritize interactions indicated by WP2 rather than an exhaustive Cartesian product. Use a small space-filling or deliberately selected boundary design once the range is frozen.

Then assess source directivity, receiver pattern, absorption, plausible bottom clearance, and representative sound-speed profiles. These additions may change the preferred propagation tool; PE, rays, or a hybrid can be appropriate here after source/phase consistency is checked. The first-stage analytic homogeneous propagator is a deliberate controlled starting point, not a universal replacement for PE.

If deployment geometry or measured channel data become available, compare trends and distributions under matched conditions. Without them, the strongest defensible deliverable is a simulation-validated idealized model with an explicit domain, not a sea-trial-validated shipboard model.

**Deliverable D5:** applicability map, sensitivity ranking, update-cost map, and clear boundaries between validated, interpolated, stress-tested, and untested conditions.

### WP6 — Extend to temporal evolution and Doppler

Represent an initial linear evolving sea using consistent modal coefficients, for example

\[
\eta(x,y,t)=\sum_m a_m\cos(\boldsymbol K_m\cdot(x,y)-\Omega_m t+\phi_m),
\qquad \Omega_m^2=gK_m
\]

in the deep-water gravity-wave regime. Finite-depth dispersion uses g K tanh(K h); the relevant water depth and spectral regime must be checked. Capillary or other short-wave extensions require their corresponding dispersion model. [E3]

Retain coefficients across time and evolve phases; add source/receiver trajectories separately. Compare deterministic single-wave and constant-velocity controls before a random time series. For slow motion, a path has a Doppler contribution approximately -(f/c) dL/dt under the chosen DSP convention. This is a diagnostic relation, not a proof of a frozen-boundary approximation.

At kilometer depths, acoustic flight time can be appreciable relative to surface evolution. Define whether t denotes transmit, boundary-interaction, or receive time. Test the error of frozen-surface snapshots against a retarded-interaction or suitable moving-boundary reference on simple controls before assuming the whole path sees one instantaneous sea surface. Refining the snapshot interval alone does not establish that approximation.

Validate temporal correlation, path-dependent Doppler, amplitude–Doppler dependence, and communication performance under temporal sampling refinement. Include common Doppler/time scaling from platform motion and residual path variation where relevant. A short coherence time does not justify independent redraws; a small delay spread does not establish WSSUS.

Begin with independent continuous records spanning many dominant wave periods; extend record length based on effective independent information. Split validation by record and seed, not by neighboring snapshots. Fit a temporal stochastic state model only after the evolving physical reference is adequate.

**Deliverable D6:** a continuous-time-compatible generator, temporal validity domain, and a justified update rule.

**Gate G6:** static fidelity must survive temporal extension; correlation and Doppler must agree within predefined uncertainty. If the frozen-boundary or temporal state model fails, retain the static model as the completed first-stage result while revising the dynamic formulation.

## 7. Experimental units, replication, and uncertainty

The independent static unit is the parent surface realization at a specified environmental condition. Frequency bins, receiver samples, repeated timing calls, and repeated packets on one surface are not additional independent sea realizations. For time-varying work, independent records/seeds are the primary units; snapshots within a record are correlated.

Tag every artifact with a **physics-evidence tier**, solver/reference identity and version, surface-definition and parent-coefficient hash, cutoff, C2 embedding/support, numerical aperture, source/receiver/frequency convention, and gates. `analytic_control`, `compact_dual_conditional`, `cross_dim_formula_diagnostic`, and `PM_unvalidated_diagnostic` must remain distinguishable; only data carrying accepted physical-domain scope may be used for later statistical fitting or PM communication conclusions.

Use common random inputs for paired algorithm comparisons and compatible factor changes. Randomize method order within a timing block and block by machine/session. Keep numerical refinement on exactly the same physical surface. Separate seeds for surface coefficients, communication data/noise, and run ordering; record all streams in a manifest. No schedule or new experiment is generated by this planning document.

Pilot batches are exploratory. Estimate uncertainty at the parent-realization level using an appropriate clustered analysis; paired differences are preferable where pairing is real. Repeated packets improve conditional BER estimation, not environmental coverage. Regularized covariance comparisons must acknowledge finite-sample bias.

After the pilot, freeze primary metrics, thresholds, model choices, and a final test size using pilot variability and the available budget. Do not repeatedly inspect one final test set and stop at a favorable result. If adaptive sampling is necessary, predefine the stopping rule or use suitable sequential intervals. Rare-outage claims require enough independent channel events; 32 or 64 surfaces cannot certify extreme tail probabilities. Report unresolved precision as inconclusive, not as a pass.

## 8. Cost model and resource decisions

For one operating condition and M requested channel realizations, compare

\[
T_{\rm physical}(M)=T_{\rm setup,P}+M T_{\rm sample,P},
\]

\[
T_{\rm statistical}(M)=T_{\rm data}+T_{\rm fit}+T_{\rm validation}+T_{\rm update}+M T_{\rm sample,S}.
\]

If sample costs are stable and T_sample,P exceeds T_sample,S, the break-even count is

\[
M_* = \max\left(0,\frac{T_{\rm data}+T_{\rm fit}+T_{\rm validation}+T_{\rm update}-T_{\rm setup,P}}
{T_{\rm sample,P}-T_{\rm sample,S}}\right).
\]

Use matched scopes and avoid double-counting shared validation costs. Report both research-development cost and recurring deployment cost, and how the result changes when conditions require rebuilding. If the denominator is nonpositive, the statistical route has no online speed advantage under those measurements.

Report total time, amortized time, memory peak, batch size, hardware, precision, startup, and uncertainty. Benchmark workloads of 1, 100, and 10,000 generated channels where feasible; these are comparison workloads, not user throughput requirements. A target such as tenfold warmed sampling speedup may be useful as a research aspiration, but total break-even and adequate fidelity determine usefulness.

The original 2D timing cannot budget 3D BIE. Before any expensive batch, run one candidate and one independent-reference case, forecast memory and total cost, and choose a bounded checkpointed campaign. Stop scaling when the reference is unaffordable, accuracy is inadequate, or the intended generation volume falls below break-even. Then reduce the experimental scope or change the formulation explicitly rather than silently relaxing convergence.

**Immediate Stream A budget decision:** use the saved `t=1` mapped-solver pilot as baseline, state one preconditioner change, define maximum Krylov work/peak-memory/wall-time and independent physical-residual limits *before running*, and compare like-for-like true residual histories. If the single bounded intervention fails, close that intervention and produce a small matched BIE feasibility estimate instead of authorizing additional unchecked solver trials. Stream B cost accounting applies only to deterministic validated controls and infrastructure; it cannot be extrapolated to target-PM generator throughput.

## 9. Proposed acceptance criteria

These are initial design tolerances proposed by this plan, not validated universal thresholds or previously confirmed user requirements. Freeze or revise them after an *authorized* pilot and communication sensitivity study, before final testing. Report all failures as well as passes. **Do not reinterpret a single-frequency, compact-dual conditional pass as the target-PM band or population meeting the provisional 1% criterion.** WP0 algebraic, numerical, independent-physics, PM-support and communications gates are distinct; no error tolerance is relaxed by the dual-track policy.

| Layer | Proposed initial criterion | Interpretation and safeguard |
|---|---|---|
| Analytic identities / convention conversion | Approximately 1e-10 normalized discrepancy for well-conditioned double-precision algebraic controls | Does not impose this tolerance on a discretized boundary integral |
| Numerical integration / broadband representation | Initially target 1e-3 change in normalized band response under refinement | Separate direct and reflected components; also track phase/delay and communication changes |
| Independent physical comparison | Use 1% band-relative reflected-field error as a provisional screening level | Not sufficient alone; inspect spectral structure, fades, and communication sensitivity; reference uncertainty should be substantially smaller |
| Statistics | Initial screening: mean reflected power within 0.5 dB and stable covariance/functionals within roughly 10% | Use uncertainty-aware metrics and an absolute floor when a denominator is near zero; tail failures can reject a model despite passing averages |
| Primary communication | Absolute SNR shift no greater than 0.5 dB at uncoded QPSK BER=1e-3 | Proposed criterion; report intervals and require uncertainty appreciably below 0.5 dB; do not extrapolate beyond the measured crossing |
| Packet-error sensitivity | Absolute PER difference no greater than 0.02 near a declared reference PER of 0.10, using a frozen packet size | Secondary proposal; specify overhead/conditioning and cluster uncertainty by channel |
| Cost | Document break-even and useful throughput under measured scopes | No predetermined single-link speedup is treated as achieved |
| Time variation | Tolerances set from Doppler/correlation sensitivity before final time-series testing | No universal Hz, coherence-time, or update-interval value is assumed |

Normalize same-realization complex errors carefully. Report component-relative metrics and an absolute scale based on the declared direct free-field band norm. Never divide only by the total field at deep interference nulls. Conversely, a small total-field error dominated by the direct path must not hide a wrong reflected component.

Phase and delay estimates near zeros require amplitude-aware interpretation; do not discard those frequencies from all tests, because the nulls may matter to communication. Use integrated metrics, raw response plots, and communication checks together. Demonstrate that numerical, reference, and Monte Carlo uncertainty do not consume the claimed tolerance.

## 10. Conditional schedule and the next executable task

The estimates below are effort bands for planning, not measured runtimes or a delivery promise. Reference development is the main uncertainty. Gates, not calendar dates, govern progression.

| Work package | Indicative focused effort | Prerequisite |
|---|---|---|
| WP0 | Completed (G0 PASS) | Retain controls; no new work unless contracts change |
| WP1/G1 — target-PM reference | **Active, one bounded decision checkpoint**; duration/cost to be measured, not promised | Existing failed mapped-Helmholtz direct pilot; frozen PM/spec/support; predeclared solver-budget and stop rule |
| Stream B — deterministic readiness | **Can start now**; separate work log/budget | Analytic flat + accepted compact-dual controls; no PM statistical inputs |
| Optional 2D→3D matched-limit audit | Low-priority, bounded regression only; may be deferred without changing any gate | Analytically matched infinite line source, y-invariant sea, same radiation/normalization; no added PM/3D physics claim |
| WP2 — formal PM mechanism pilot | Deferred; original 1–2-week estimate is not authorization | **New scope-specific G1 approval**, frozen spectrum/support and accepted independent reference; `WP2_formal=false` now |
| WP3 — statistical generator | Deferred; original development estimates remain provisional | G2 and accepted physical-data provenance |
| WP4 — receiver controls | AWGN/flat/two-path baseline work allowed; PM generator/communication validation deferred | Stream B contract for controls; WP4 acceptance additionally requires accepted physical/statistical inputs |
| WP5 | Deferred | Accepted WP3/WP4 and defined applicability domain |
| WP6 | Deferred separate work | Accepted static model and dynamics specification |

**Next executable work is two bounded, non-overlapping tasks rather than redoing WP0/WP1:**

1. **A / reference decision (first priority for target-PM certification):** review the recorded mapped `t=1` PM failure; make and test *one* predeclared preconditioner change with caps on solver steps, memory, and wall time. Check true PDE/DtN/Dirichlet residuals independently of the linear solver. If it still fails, stop escalation and prepare a costed, mathematically matched 3D BIE/BEM minimum-pilot decision. If it passes, proceed in gated order to real boundary, horizontal/vertical mesh, domain and surface-support evidence **before** matched FullFK–PM physics. Do not run PM ensembles, support24 extrapolation, or WP2 during this task.
2. **B / parallel readiness (does not wait for target-PM G1):** freeze native/DSP, phase/delay/CIR and frequency contracts; build reproducible outputs and physics-tier labels using image and accepted dual controls; pass AWGN, noiseless known-channel and flat two-path modem sanity tests; record end-to-end costs for these exact deterministic workloads. **Only if it is genuinely cheap after the mandatory controls**, optionally add the matched 2D line-source/y-invariant-3D formula regression from Section 3.1; otherwise document deferral. No PM parameter fitting, PM BER/outage result, FastFK optimization, or formal WP2.
3. **Decision review:** report each stream separately. Stream A yields `reference_ready` / `INCONCLUSIVE`, residual and support evidence and a next-method stop/go choice. Stream B yields only deterministic **infrastructure readiness**. Neither authorizes formal WP2 unless a later, explicit, target-domain G1 approval is documented.

The next Codex Goal-mode execution shall use the accompanying `codex_goal_two_track_2026-10-14.md`: continue through preflight, the **single** capped preconditioner decision, independent residual/audit work, B-stream deterministic tests, and the final separated reports without requiring step-by-step approvals. It must stop expensive Stream A exploration at its declared cap and must not interpret autonomy as permission to expand the PM physics domain.

The original 32-response and broadband statistical experiments remain **future gated designs**, not the current task. The ultimate objective is still a defensible *single-link communication channel*; once physical data exist, select P0, S1 or S2 by accuracy, cost and communication evidence, not by prior preference for FastFK or any one statistical family.

## 11. Deliverables and repository organization

The static research is complete only when D0–D5 provide a reproducible physical domain, a selected generator or direct-sampling solution, held-out statistical and communication evidence, and full cost accounting. D6 is a distinct later deliverable. **Stream B readiness is an interim engineering/research artifact, not D2–D5 acceptance**; a frozen compact-dual numerical/communication control must never be presented as validation of the target finite-band PM population.

For future implementation, follow existing repository responsibility boundaries:

- Independent reference solvers and validation-only helpers: `scripts/validation/support/`.
- Invariants, analytic controls, and acceptance workflows: `scripts/validation/`.
- Controlled route comparisons: `scripts/comparisons/`.
- Surface/condition ensembles: `scripts/experiments/`.
- Saved-data figures and summaries: `scripts/reporting/`.
- Validated reusable components, when promoted: appropriate `src/surface/`, `src/receiver/`, `src/statistics/`, or `src/channel/` modules.
- Case-specific data: `results/validation/`, `results/comparisons/`, or `results/experiments/`; conclusions in `reports/`.

Persist raw native complex fields, DSP responses, frequency axes, phase metadata, realized surface coefficients, source definition, reference settings, numerical guards, seed manifests, timing scopes, fitted parameters (only once fitting is authorized), and train/test membership (only for authorized datasets). Preserve one explicit `cross_dim_formula_diagnostic` label whenever an optional dimensional-limit test is run; its output is not a 3D PM reference. Track analytic/convention checks, algorithm fidelity, physical validity, reference credibility, surface-support status, and statistical/communication adequacy **separately**. Emit `G1_target_PM`, `WP2_formal`, physics-tier provenance and exact support/cutoff metadata on every relevant research artifact.

Preserve the public `vertical_channel_model` interface and `H_f = H_direct_f + H_reflect_f`. New production integration is a separate stage requiring the repository's scalar-frequency and wideband regressions. Use the shared MATLAB MCP session for execution. Do not recover obsolete evidence from `cash/`.

## 12. Unresolved inputs and how the plan handles them

| Missing input | Current treatment | Needed before |
|---|---|---|
| Source model, calibration, directivity | Ideal normalized point-source controls; explicit sensitivity later | Device-specific received power and range claims |
| Hull clearance and platform geometry | Idealized hull-free domain, with no transfer claim | Shipboard validation |
| Maximum depth and offsets | Proposed nested research settings, not engineering limits | Final applicability envelope |
| Sea-state distribution and short-wave model | Explicit finite-band PM candidate, cutoff/direction sensitivity; no unbounded-ocean certificate | Population weighting or ocean-wide claims |
| Target-PM independent reference and linear convergence | Ordinary FE/Taylor-TFE and bounded mapped `t=1` pilot did not pass original full-amplitude PM checks; one capped preconditioner step, then matched BIE feasibility if necessary | FullFK–PM physical acceptance, formal WP2 and PM-based statistical inference |
| Finite support versus numerical convergence | Parent PM seed/cutoff and C2 support12 are frozen; support24 changes the modeled physical sea | Noncompact PM claims, support-independent physical accuracy and WP2 domain approval |
| Water depth and source bottom clearance | No bottom in first half-space study | Claims about a bottom-mounted installation |
| Computing and update budgets | Measured cost curves and break-even workloads | Optimization target and deployment choice |
| Measured channel data | Independent numerical reference plus explicit measurement gap | Field-validated model claim |
| Required modem and tolerance | Proposed QPSK benchmark and tolerances above | Final engineering acceptance |

This uncertainty register does not require another clarification round before analytic work. It prevents provisional choices from becoming hidden assumptions in later conclusions.

## 13. Evidence, references, and planning assistance

### Local evidence

- **L1.** [Research background and handoff](deep_sea_channel_research_background_for_new_ai.md): confirmed requirements and full historical evidence register.
- **L2.** [Single-receiver timing report](fast_fk_single_receiver_timing_report.md): fixed 4 kHz timings, source fixture, repetition protocol, and exclusions.
- **L3.** [Multifrequency consistency report](spectral_fast_kirchhoff_multifrequency_consistency_report.md): ranks 3/12, phase conventions, BIE anchors, and diagnostic CIR limitations.
- **L4.** [External-reply reassessment](deep_sea_channel_external_reply_reassessment.md): corrected geometry, PM arithmetic, statistical interpretation, published deployment examples, and scenario provenance.
- **L5.** [Technical guide](../vertical_comm_guide.md), [2D FullFK reference](../scripts/validation/support/solve_full_kirchhoff_2d.m), and [WAPE implementation](../src/propagation/vertical_wape_propagator.m): existing model definitions, numerical conventions, and boundaries.
- **L6.** [FullFK PM validity boundary](full_kirchhoff_pm_validity_boundary_report.md): conditional physical approximation limits; not a universal cutoff or 3D certification.
- **L7.** [3D compact dual/frequency certification](deepsea3d_g1_frequency_report.md), [azimuth](deepsea3d_g1_azimuth_report.md), and [finite-amplitude domain](deepsea3d_alpha25_domain_report.md): conditional per-case G1 evidence, not a PM certificate.
- **L8.** [3D PM feasibility](deepsea3d_pm_feasibility_report.md), [Taylor/TFE dealiasing](deepsea3d_tfe_dealiasing_report.md), and [mapped `t=1` direct-solve pilot](deepsea3d_mapped_direct_report.md): target-PM reference still fails numerical acceptance; current decision starts from this recorded endpoint, not from a blank implementation.

### Primary and institutional sources

- **E1.** Pierson and Moskowitz (1964), [A proposed spectral form for fully developed wind seas based on the similarity theory of S. A. Kitaigorodskii](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/JZ069i024p05181). Publisher abstract/record checked; supports the fully developed-sea and wind-convention context, not a new numeric cutoff in this plan.
- **E2.** Neeser and Massey (1993), [Proper Complex Random Processes with Applications to Information Theory](https://www.isiweb.ee.ethz.ch/archive/massey_pub/pdf/BI437.pdf). Original author-hosted paper; used for centered covariance/pseudo-covariance semantics, not a claim that the target channel is Gaussian.
- **E3.** MIT Marine Hydrodynamics, [Solution of the Dispersion Relationship](https://web.mit.edu/fluids-modules/www/potential_flows/LecturesHTML/lec19bu/node5.html). Institutional teaching source for the linear gravity-wave dispersion relation and its deep-water limit; not a validation of the proposed moving acoustic boundary model.
- **E4.** Qarabaqi and Stojanovic (2013), [Statistical Characterization and Computationally Efficient Modeling of a Class of Underwater Acoustic Communication Channels](https://doi.org/10.1109/JOE.2013.2278787), IEEE Journal of Oceanic Engineering 38(4), 701–717. The accessible abstract describes a physical/statistical channel formulation. The author-hosted full PDF could not be fetched during this turn, so no detailed theorem, validity region, or numeric performance claim is imported. Read the full paper before claiming novelty or adopting its assumptions for this geometry.

Sources were checked on 2026-10-09. Equations for free-space propagation, image geometry, proposed cost accounting, and planning values are explicitly given above; none is a newly measured project result.

### Planning assistance

The local experimental-design skill informed independent-unit definitions, paired seeds, timing blocks, and separation of pilot/development/final validation. Its provenance is Kassis, T., Agarwal, V., He, Y., Patel, D., and Brueckner, A. M. (2026), [Scientific Agent Skills: A Library of Procedural Knowledge for Research Agents](https://doi.org/10.48550/arXiv.2609.00065). The arXiv author/year record was checked; the current record is revised September 2026. This is a workflow acknowledgment, not evidence of acoustic validity. No acoustic-formula-audit skill workflow was used for this plan.
