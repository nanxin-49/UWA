# UWA 垂直水声信道 MATLAB 仓库技术审查报告

审查日期：2026-09-23
审查方式：只读源码审查 + 有效文档核对（未运行任何代码，未读取 `cash/`，未分析 results 数据文件）。
依据优先级：`src/` 源码 > `vertical_comm_guide.md`（权威技术指南，2026-09-03）> `reports/` 单次验证报告 > `docs/history/`（不作为依据）。

---

## 1. 项目整体目标

研究**海底发射机（Tx，z=100 m）到近海面接收机（Rx，z=3 m）的垂直水声信道**，标准几何为均匀介质 c0=1500 m/s、名义直达/海面反射相对时延 4 ms（`vertical_comm_guide.md` §1）。项目用 PE/WAPE 分步角谱法逐频计算直达与海面反射频响 H(f)=H_direct+H_reflect，在其上构建四层研究能力：

1. **物理仿真链**：PE/WAPE 传播 + 4 种海面边界模型（显式 Kirchhoff、k 域 Kirchhoff、joint-kstat 统计屏、SSA 统计核）+ 气泡有效介质（默认关闭）；
2. **统计加速链**：cached forward + 精确离散伴随 + 解析跨频协方差/伪协方差 C/P 收缩，免除逐 realization 的 PE；
3. **条件统计信道生成器**：U=5/U=8 两个固定风速节点、F=64（4–8 kHz）的复高斯（proper/improper）信道采样；
4. **通信验证链**：QPSK/MPSK + AWGN + 已知信道 MMSE 的 Monte Carlo BER/SER 评估，用于验证统计生成器与物理链的一致性。

外部基准使用 Bellhop 2020（validation-only internal-wall 二进制，`scripts/validation/support/bellhop_internal_*/`）和独立 2-D Helmholtz BIE 求解器（`scripts/validation/support/solve_helmholtz_bie_halfplane_vertical.m`）。

---

## 2. 仓库结构

| 目录 | 功能 | 规模 |
|---|---|---|
| 根目录 | 4 个兼容入口包装 + `setup_vertical_project.m` + 3 份权威文档 | ~10 文件 |
| `src/` | 可复用生产实现，按 channel/propagation/surface/receiver/statistics/communication/bubble 划分 | 40 个 .m，约 8600 行 |
| `examples/` | 演示脚本实际内容（入口在根目录） | 2 个 |
| `scripts/validation/` | 验证入口（PE 自检、伴随、Bellhop、BIE、窗口收敛、条件模型、通信） | ~90 个 .m |
| `scripts/validation/support/` | Bellhop/Weyl/Li2009/properness 等验证专用辅助 + 5 套 internal-wall Bellhop 源码 | ~30 个 |
| `scripts/comparisons/` | 海面模型、气泡模型、SSA/Kirchhoff 对比 | 8 个 |
| `scripts/experiments/` | Monte Carlo、标定、参数扫描 | 9 个 |
| `scripts/reporting/` | 图册、报告绘图 | 16 个 |
| `reports/` | 一次性验证报告（Markdown） | ~90 份 |
| `results/` | MAT/CSV/PNG/MP4 输出（本次未分析） | — |
| `docs/history/`, `docs/old/` | 非权威历史规格与旧代码 | — |
| `cash/` | 退役数据隔离区，Git 忽略，正常检索不读取（本次遵守） | — |

关键文档：`README.md`（导航）、`vertical_comm_guide.md`（物理/接口/验证边界权威说明）、`PROJECT_CONTEXT.md`（追加式项目日志，顶部为当前状态索引）、`scripts/README.md`（验证运行索引）、`reports/current_project_status_20260820.md`（状态快照）。

---

## 3. 核心文件索引

### P0：核心算法

| 文件 | 类别 | 作用 | 重要性 |
|---|---|---|---|
| `src/propagation/vertical_wape_propagator.m` (1247 行) | propagation | PE/WAPE 分步角谱推进；直达与 Tx→Surface→Rx 两段链；1/R 自检；全输出组装 | P0 |
| `src/surface/pm_surface_boundary_model.m` (2110 行) | surface | PM 谱生成与海面边界算子总入口：kirchhoff_spatial / kirchhoff_kdomain / kirchhoff_kstat / ssa_stat_kernel 四分支 | P0 |
| `src/channel/vertical_channel_model_impl.m` (792 行) | channel | 公共 API 实现：默认参数、~100 项输入校验、调用 propagator、输出语义组装 | P0 |
| `src/channel/apply_pe_channel_phase_reference_vertical.m` | channel | 约化 PE 包络 → direct_dsp / absolute_physical 相位转换（Δτ0 相对时延因子） | P0 |
| `src/receiver/build_adjoint_receiver_projection_vertical.m` | receiver | 精确离散伴随接收核构建（q_i、a_pe_i、PM 映射） | P0 |
| `src/receiver/apply_adjoint_receiver_to_surface_vertical.m` / `apply_forward_surface_to_receiver_vertical.m` | receiver | 伴随算子（conj 转置反向 march）与共享 cached 前向算子 | P0 |
| `src/receiver/contract_kstat_receiver_stats_vertical.m` | receiver | 接收端 C_H/P_H 二阶统计收缩（dense/FFT 两法，direct-DSP 转换） | P0 |
| `src/statistics/build_kirchhoff_kstat_joint_model_streaming_vertical.m` | statistics | 跨频联合 C/P 的级数展开流式构建（631 s→6.5 s、6.6 GiB→1.7 GiB 的优化版） | P0 |
| `src/statistics/build_kirchhoff_kstat_joint_model_vertical.m` | statistics | 联合因子模型原版（F×F×N_K 全张量，作基准） | P0 |
| `src/communication/build_communication_taps_vertical.m` | communication | H(f)/h(t) → 符号率 taps（f_ref 平移、插值、IFFT、最短循环能量窗） | P0 |
| `src/communication/evaluate_mpsk_channel_ensemble_vertical.m` | communication | 统一 MPSK Monte Carlo 评估（同步、MMSE、Wilson/聚类 bootstrap CI） | P0 |
| `src/communication/modem_psk.m` | communication | M-PSK 调制/相干解调/BER/SER | P0 |

### P1：重要模块

| 文件 | 类别 | 作用 | 重要性 |
|---|---|---|---|
| `src/propagation/gaussian_source_initial_field_vertical.m` | propagation | 生产 Gaussian 初场 exp(-ρ²/2σ²)，σ∈[0.2,0.5] m | P1 |
| `src/surface/raw_pm_spectrum_grid_vertical.m` | surface | 可审计 raw-PM 谱网格（含无限域方差、覆盖比诊断） | P1 |
| `src/surface/sample_raw_pm_surface_vertical.m` | surface | 显式 PM 海面采样 | P1 |
| `src/bubble/bubble_environment_vertical.m` / `bubble_effective_medium.m` / `bubble_hall_spectrum.m` | bubble | 气泡有效介质（c_eff、α）注入 PE 相位屏；Hall1D 谱 | P1 |
| `src/channel/build_channel_cir_vertical.m` | channel | H(f)→CIR（direct_dsp 用 IFFT，physical 用 FFT/N，显式相位约定） | P1 |
| `src/statistics/estimate_conditional_channel_stats_vertical.m` / `sample_conditional_channel_vertical.m` / `build_conditional_channel_library_vertical.m` / `sample_conditional_channel_library_vertical.m` | statistics | 条件统计模型估计/采样/两节点库 | P1 |
| `src/statistics/sample_kirchhoff_kstat_joint_frequency_vertical.m` / `sample_kirchhoff_kstat_factor_model_vertical.m` | statistics | 联合因子模型采样 | P1 |
| `src/communication/noise_inject_vertical.m` | communication | 信道卷积后 AWGN 注入（唯一噪声入口） | P1 |
| `src/receiver/build_cached_joint_kstat_pe_executor_vertical.m` / `run_cached_joint_kstat_pe_executor_vertical.m` / `generate_cached_kstat_ensemble_vertical.m` / `generate_cached_kdomain_ensemble_vertical.m` | receiver | cached 前向批执行器与 ensemble 生成 | P1 |
| `src/statistics/build_surface_empirical_channel_model_vertical.m` / `sample_surface_empirical_channel_vertical.m` | statistics | 经验信道模型（realization→统计生成器） | P1 |

### P2：验证/实验脚本（按目录汇总）

| 目录 | 代表入口 | 重要性 |
|---|---|---|
| `scripts/validation/` (~90) | `validate_pe_as_freefield_vertical.m`、`validate_adjoint_pe_receiver_projection_vertical.m`、`validate_pe_bellhop_pm_*`（Stage 0–4）、`validate_pe_bellhop_controlled_comparison.m`、`validate_pe_bellhop_helmholtz_bie_reference.m`、`validate_u5/u8_conditional_channel_generator_vertical.m`、`validate_two_node_communication_vertical.m`、`validate_reflected_chain_window_convergence_vertical.m`、`validate_li2009_explicit_surface_vertical.m`、`validate_2k_*.m` | P2 |
| `scripts/comparisons/` (8) | `compare_ssa1_ssa2_*`、`compare_kirchhoff_kdomain_kstat_*`、`comm_compare_bubble_models_*` | P2 |
| `scripts/experiments/` (9) | `monte_carlo_comm_surface_psk_vertical.m`、`calibrate_hall1d_bubble_vertical.m`、`sweep_*` | P2 |
| `scripts/reporting/` (16) | `generate_pe_propagation_atlas_vertical.m`、`plot_u5_f64_conditional_validation_vertical.m` | P2 |

### P3：辅助文件

| 文件 | 作用 | 重要性 |
|---|---|---|
| `vertical_channel_model.m` / `main_vertical.m` / `explain_main_vertical.m` / `comm_main_vertical_psk.m` | 根目录兼容入口（薄包装，保持历史调用名） | P3 |
| `setup_vertical_project.m` | 幂等路径初始化（不动当前目录） | P3 |
| `scripts/bootstrap_project.m`、`scripts/project_result_dir.m`、`scripts/project_result_file.m` | 旧脚本引导与结果路径工具 | P3 |
| `examples/main_vertical_demo.m`、`examples/comm_main_vertical_psk_demo.m` | 演示内容 | P3 |

---

## 4. 算法流程（核心调用链）

三个入口都是薄包装：`vertical_channel_model.m:6-7` 先 `setup_vertical_project()` 再转调 `src/channel/vertical_channel_model_impl.m`；`main_vertical.m:3` 与 `comm_main_vertical_psk.m:3` 分别 `run` examples 下的演示脚本。

```
参数输入 paramsV (频率/几何/网格/海面/气泡 ~100 字段)
        │
        ▼
vertical_channel_model.m → setup_vertical_project() (路径)
        │
        ▼
src/channel/vertical_channel_model_impl.m
  local_prepare_config: 默认值合并 + 逐字段校验 + dx/dz/numstep 派生
        │
        ▼
src/propagation/vertical_wape_propagator.m   ←── 核心推进
  ├─ src/propagation/gaussian_source_initial_field_vertical.m  (初场)
  ├─ local_absorption_profile (边缘 sponge α(x,y))
  ├─ 频率轴解析 local_resolve_frequency_axis (宽带自适应 Δf≈0.9/Δτ)
  │
  │  逐频 f ∈ f_axis:
  │  ① 直达: z_tx→z_rx marching（Strang 分步角谱）
  │       psi_k = fr .* fft2( screen .* ifft2( fr .* psi_k ) )   [L189]
  │       fr    = exp(-i·ds/2 · κ²/(k0+√(k0²-κ²)))               [L148-149]
  │       screen= exp(-i·k0·ds·(U(x,y) - i·α_total/k0))          [L1119]
  │  ② 反射: local_march_field(Tx→z=0)                            [L273-275]
  │       → src/surface/pm_surface_boundary_model.m                [L276-277]
  │          ├ kirchhoff_spatial : ψ_ref = R0·e^{i·k0(cos_i+cos_r)η} ⊙ ψ_inc
  │          ├ kirchhoff_kdomain : 同一屏的 FFT 隐式卷积
  │          ├ kirchhoff_kstat   : 统计屏 R_coh + 散射功率谱 P_sca
  │          └ ssa_stat_kernel   : SSA1/SSA2 统计核 + 能量限额
  │       → local_march_field(z=0→Rx) → h_reflect                   [L279-282]
  │  ③ 1/R 中心振幅自检（斜率∈[-1.05,-0.95]，RMS≤10%）             [L755-781]
  │
  ▼
src/channel/apply_pe_channel_phase_reference_vertical.m
  H_reflect = e^{-i2πfΔτ0}·H_ref^red   (direct_dsp 默认)         [L51]
  H_f = H_direct + H_reflect                                       [L86]
        │
        ▼
  output（H_direct_f/H_reflect_f/H_f + reduced/physical/元数据）
        │
        ├──────────────────────────────┐
        ▼                              ▼
  通信链 src/communication/         研究加速链（固定环境，验证专用）
  build_communication_taps_vertical   cached forward + 精确伴随 src/receiver/
   (f_ref 平移→插值→IFFT→             build_adjoint_receiver_projection_vertical
    最短循环 99.9% 能量窗→因果 taps)     → q_i = A_i^H r, a_i = conj(ψ_inc)⊙q_i
        → conv(tx,h)                contract_kstat_receiver_stats_vertical
        → noise_inject_vertical       → C_H(i,j)=ã_i^H C_δG,ij ã_j (解析二阶统计)
        → 频域 MMSE (已知信道)      src/statistics/ 条件模型链
        → modem_psk 解调               build_kirchhoff_kstat_joint_model_streaming
        → evaluate_mpsk_channel_       → estimate/sample_conditional_channel
          ensemble_vertical             → U=5/U=8 两节点库 → 大规模通信 Monte Carlo
```

方向约定：z 向下为正、海面 z=0、向上传播即 z 减小；PE 推进的是相对参考载波的复包络（`vertical_comm_guide.md` §1-2）。依赖方向由 `src/README.md` 固定：channel→propagation→{surface,bubble}，receiver↔statistics→communication。

---

## 5. PE/WAPE 传播分析

**解决问题**：均匀（或分层预留）介质中垂直向上 97 m 路径的单频/宽带复包络传播，输出接收点频响与诊断场。

**数学模型**（`vertical_wape_propagator.m`）：
- 一步推进为对称 Strang 分裂：半步横向谱因子 `fr`、中间介质相位屏 `screen`、再半步 `fr`（`vertical_wape_propagator.m:189`）。`fr = exp(-i·(ds/2)·κ²/(k0+√(k0²-κ²)))`（L148-149）正是精确单程平方根算子 `exp(i·ds·(kz-k0))` 的恒等变形（kz-k0 = -κ²/(k0+kz)），因此横向绕射部分**不是近似 Padé，而是精确角谱因子**；近似的来源只在介质屏与分裂误差。
- 介质屏（L1119）：`screen = exp(-i·k0·ds·(U - i·α/k0))`，U=(c_eff-c0)/c_eff 为折射指数扰动（含气泡有效声速），α 为 sponge+气泡吸收。
- 直达与反射各自存储**约化包络**（消去确定性载波延迟），载波相位只在传播完成后由 `apply_pe_channel_phase_reference_vertical.m` 统一转换；名义 Δτ0=(z_tx+z_rx-2z_s)/c0 = 4 ms。
- 横向网格为周期 FFT 域，边缘用立方 sponge 抑制绕回；**sponge 不是无限孔径替代品**（指南 §3 明确警告）。
- 内建 1/R 几何自检（L755-781）作为均匀介质直达路径的运行时回归门。

**关键函数**：`vertical_wape_propagator.m` 主循环 L117-344；`local_march_field` L783-857（反射两段共用）；`local_phase_screen` L1086-1120；`local_resolve_frequency_axis` L1161-1191。

**数值验证状态**：与独立一步精确角谱对比达双精度（`reports/pe_1d_validation_bridge_report.md`：复误差 ~1.1e-12）；1-D split-step、exact-AS、生产 PE `ky=0` 投影三方闭合（`scripts/validation/validate_pe_1d_validation_bridge.m`）；相位符号独立审计通过（`reports/pe_phase_reference_release_candidate_report.md`）。

**已知限制**：公共默认窗口 50 m+sponge ≠ 验证收敛配置 192.1875 m/no-sponge（详见 §10）；`env_mode='layered'` 仅留接口与测试剖面（`local_sound_speed` L725-753），无验证链；GPU 分支存在但未纳入验证。

---

## 6. 海面反射模型分析

**解决问题**：z=0 压力释放（R0=-1）粗糙海面对上行场的反射，统一入口 `src/surface/pm_surface_boundary_model.m`。

**PM 海面统计**（`local_pm_spectrum` L514-543）：E1D(K)=α/(2K³)·exp(-βg²/(U⁴K²))，α=8.10e-3、β=0.74；Φ2D=E1D/(2πK)；无限域方差 σ²=αU⁴/(4βg²)。显式海面由谱幅度×复高斯白噪声 IFFT 合成，默认按目标 Hs 重标（`target_hs` 模式）。

**四个边界分支**：

1. **kirchhoff_spatial**（公共默认，L290-339）：显式 η(x,y) 相位屏。相位 `delta_phi = k0·(cos_i+cos_r)·η`，法向入射退化为 2k0η；反射场 `ψ_ref = R0·e^{iΔφ}⊙ψ_inc`。历史上曾误用 2k0(cos_i+cos_r)η（法向 4k0η），已修正并在 `meta.kirchhoff_phase_screen_legacy_note` 留档（L489-490）。
2. **kirchhoff_kdomain**（L1955-2011）：同一相位屏在 k 域按 FFT 循环卷积施加，与 spatial 分支可用 `boundary_check_equivalence` 互检到 1e-10。
3. **kirchhoff_kstat**（L771-956）：不生成显式海面，直接建随机残差统计。核心随机变量 δG_i=R0·e^{iα_iη}-R_coh,i，α=2k0；相干反射 R_coh=R0·exp(-0.5·(2k0)²σ²)（压力释放即 -exp(-2k0²σ²)，L902）；散射功率谱 P_sca=|R0|²·(S_δG ∘ 卷积 |Ψ_inc|²)，其中 C_δG(ρ)=exp(α²(C_η(ρ)-σ²))-exp(-α²σ²)（L801）。能量闭合 |G_mean|²+ΣS_δG=1 被显式审计（L813-816）。
4. **ssa_stat_kernel**（L578-769）：一阶 Dirichlet SSA 统计核 G_SSA1=4γ(K)γ(K')，含 SSA2 (Broschat coherent) 相干项；散射功率受入射能量限额约束（L624-657）。**定位为诊断分支**，不含 NLSSA/多次散射/阻抗边界（meta.limitations L953-955）。

**跨频联合统计**（`src/statistics/build_kirchhoff_kstat_joint_model_vertical.m:32-47`）：
C_ij=R0_i·conj(R0_j)·[exp(-½(α_i²+α_j²)σ²+α_iα_jC_η) - exp(-½(α_i²+α_j²)σ²)]，P_ij 同型但指数取 -α_iα_jC_η——同一海面在多频上的联合协方差/伪协方差。流式版用级数展开 C_G(ρ)=Σ b_n b_n^H t(ρ)^n 把构建时间压到 1%（`reports/u8_joint_optimization_two_node_library_report.md`）。

**外部基准裁定（最重要的物理结论）**：
- PE–Bellhop 受控比较（`reports/pe_bellhop_pm_complete_execution_summary.md`）：数值/几何链 Stage 0A–0E 全 PASS，但 PE Kirchhoff 屏与 Bellhop local-specular 的 reflected-only 差异稳定为 ~0.31 dB / -2.25 rad，判定为**反射模型差异而非数值失败**。
- **Helmholtz BIE 第三参考**（`reports/pe_bellhop_helmholtz_bie_final_report.md`，R3–R5）：A=0.01→0.20 m 时 PE–BIE 误差 E_G 从 0.0032 增至 0.081，而 Bellhop–BIE 仅 8e-6→1e-4——**Bellhop 显著更接近独立 Helmholtz 参考**。
- 海面算子改进工作流（`reports/pe_surface_operator_improvement_final_report.md`）：G2 kz-aware 局部修正有改进，G3 角度/斜率修正失败，最终结论 **NONLOCAL_EFFECT_REQUIRED**——剩余差异需要非局域海面算子，canonical fixed-PM 工况（RMS 高度 ~0.19 m）被判为 Region III（PE Kirchhoff 适用域外，`reports/pe_bellhop_controlled_comparison_final_report.md`）。

---

## 7. 统计信道模型分析

**解决问题**：免 PE 的接收端二阶统计与大规模信道采样。

**三层结构**：

1. **精确离散伴随**（`src/receiver/`）：cached forward 把固定频率/网格/环境的 Surface→Rx 步进因子缓存，对许多海面 realization 复用；伴随算子是严格的共轭转置反向 march（指南 §5 给出 `q_k = conj(fr).*fft2(conj(screen).*ifft2(conj(fr).*q_k))`）。接收灵敏度核 q_i=A_i^H r，权重 a_i=conj(ψ_inc,i)⊙q_i（`build_adjoint_receiver_projection_vertical.m:37-41`）。验证到双精度（内积检验 7.8e-15，`reports/adjoint_pe_receiver_projection_feasibility_report.md`）。
2. **解析 C/P 收缩**（`contract_kstat_receiver_stats_vertical.m`）：C_H(i,j)=ã_i^H C_δG,ij ã_j、P_H(i,j)=ã_i^H P_δG,ij ã_j*；PM 大网格→PE 小网格用中央裁剪 E 的伴随零嵌（L20-36 有位级嵌入校验）；direct-DSP 统计经 D=diag(e^{-i2πfΔτ0}) 转换（L107-109）。F=64/512 realization 的样本统计落在 split-sample floor 内（0.851<1.25，`reports/pe_phase_reference_release_candidate_report.md`）。
3. **条件统计生成器**（`src/statistics/`）：从验证过的接收端样本估计 μ、C、P，proper 路径走复 EVD（U·√d·z），improper 路径走实增广 EVD，按 properness null 检验自动选择（`sample_conditional_channel_vertical.m:17-34`）；两节点库明确**拒绝未登记风速、不插值**（`build_conditional_channel_library_vertical.m`）。U=5 F=64 验证：joint-kstat 散射项不能拒绝 proper 复高斯假设；相对 held-out kdomain 的 PDP/LFM 相关 0.952/0.984（`reports/u5_conditional_channel_generator_f64_report.md`）。

**角色边界**（指南 §6、§10）：条件生成器承担大规模通信 Monte Carlo，**不替代 PE 物理回归**；schema 2.x 显式记录 target_reference/几何/频率轴，旧 schema 1.x 只有在几何可靠时才允许迁移。

---

## 8. 通信链分析

**解决问题**：给定 H(f) 或 h(t)，评估 M-PSK 在该信道下的 BER/SER 统计。

**链路**（`evaluate_mpsk_channel_ensemble_vertical.m` + `build_communication_taps_vertical.m`）：
1. H(f)（要求 direct_dsp，显式拒绝 legacy_reduced/absolute_physical 输入，`build_communication_taps_vertical.m:15-21`）平移到 f_ref，线性插值到符号率 FFT 网格，IFFT 得 h_full；
2. **最短循环 99.9% 能量窗**解绕为因果 taps（`local_shortest_circular_window`，滑窗算法）——取代旧 peak-sync（后者在 U=8 会丢弃 61.8% 峰前能量，`reports/two_node_statistical_channel_communication_validation_report.md` §4）；
3. 线性卷积取因果头 → **噪声只在卷积后注入**（`noise_inject_vertical.m`，AWGN，EbN0 以 rx_clean 为参考）→ 已知信道频域 MMSE（正则 1e-6）→ 相干硬判决；
4. 统计：Wilson 95% CI + 按信道 realization 聚类的 bootstrap CI。

**F 与时延窗语义**（指南 §7）：F=64 覆盖 4–8 kHz（df=63.5 Hz，无模糊窗 15.75 ms）；F=65/62.5 Hz 用于 4 ms 相位验收；F=9 smoke 的 2 ms 窗会混叠，禁用于相位结论。通信 taps 间隔 1 ms 是符号率投影，不是物理时延分辨率（代码与报告均显式声明）。

**核心验证结果**（`reports/two_node_statistical_channel_communication_validation_report.md`）：四种信道来源（kdomain+cached PE、joint-kstat+cached PE、统计 full、统计 99.9%）在 U=5/U=8、全部 Eb/N0 点的聚类 95% CI 均重叠；U=8/U=5 平均 BER 比 1.35–2.17，物理趋势一致。**但无噪声 BER 底为 0.2%–0.5%**（§7），源于循环 FFT 逆滤波近似线性卷积、无保护间隔、0.999 截断与固定正则——高 SNR 绝对性能不可作为物理极限。

---

## 9. 已实现和已验证内容（证据表）

状态分级沿用仓库自身语义：PASS（满足预设门限）/ CONVERGED / PASS_WITH_LIMITS / PRELIMINARY / OPEN。

| 结论 | 验证脚本 | 证据 | 状态 |
|---|---|---|---|
| PE 分步角谱实现 = 独立一步精确角谱；1-D bridge 三方闭合 | `scripts/validation/validate_pe_as_freefield_vertical.m`、`validate_pe_1d_validation_bridge.m` | 复误差 ~1.1e-12；平面压力释放符号误差 2.6e-13（`reports/pe_1d_validation_bridge_report.md`） | **已验证** |
| direct_dsp 载波相位统一直达/反射；F=65 名义 4 ms 时延 | release candidate 工作流（`scripts/README.md` L115-160） | run `phase_rc_20260722_174945` 全 13 项 PASS，4 ms 时延误差 3.47e-15（`reports/pe_phase_reference_release_candidate_report.md`） | **已验证** |
| 精确离散伴随 = cached forward 的共轭转置；投影；dense=FFT C/P | `validate_adjoint_pe_receiver_projection_vertical.m` | 内积 7.8e-15；PM 零嵌位级一致；F=9/4096 与 F=64/512 样本统计在 split floor 内（`reports/adjoint_pe_receiver_projection_feasibility_report.md`） | **已验证**（限 uniform/CPU-double/单最近网格点 Rx/无气泡 Doppler） |
| 固定海面完整反射链 192.1875 m/no-sponge 收敛 | `validate_reflected_chain_window_convergence_vertical.m` | 50 m 窗反射 TL 偏 1.15 dB；192 m 边缘 -52.3 dB PASS（`reports/pe_reflected_chain_window_validation_report.md`） | **已验证**（4 kHz + 3–5 kHz 固定面） |
| 随机海面 4 kHz 下 192.1875 m 窗口鲁棒性 | `validate_random_surface_window_robustness_vertical.m` | 3 海况×5 seed = 15/15 严格通过；160 m 边缘 0/15（`reports/pe_random_surface_window_robustness_report.md`） | **已验证**（仅 4 kHz；宽带成对未完成） |
| U=5/U=8 F=64 条件统计生成器 | `validate_u5/u8_conditional_channel_generator_vertical.m` | properness、秩选择、10k 样本；PDP/LFM 相关 ~0.95（`reports/u5_conditional_channel_generator_f64_report.md`、`u8_joint_optimization_two_node_library_report.md`） | **已验证**（仅两个固定风速节点） |
| 统计生成器复现物理链 BER/SER 趋势 | `validate_two_node_communication_vertical.m` | 四来源全部 Eb/N0 点聚类 CI 重叠 1.0；U=8 更差的物理趋势一致（`reports/two_node_statistical_channel_communication_validation_report.md`） | **已验证**（受 0.2–0.5% 无噪声 BER 底限制） |
| PE–Bellhop 数值/几何链（canonical mapper、flat source、常数高度符号、数值预算） | `validate_pe_bellhop_pm_stage0_*` 等 | wall residual ~8e-15 m；PE/Bellhop 扫描均在 0.1 dB/0.02 rad 预算内（`reports/pe_bellhop_pm_complete_execution_summary.md` Stage 0A–0E） | **已验证** |
| PE Kirchhoff 屏 vs Bellhop local-specular 一致性 | `validate_pe_bellhop_pm_stage1_tier1.m` 及后续 | 稳定差异 0.3104 dB / -2.2482 rad（X 源重跑 0.31037 dB），判为模型差异 | **部分验证**（差异已定征，一致性不成立） |
| PE Kirchhoff 屏 vs 独立 Helmholtz BIE | `validate_pe_bellhop_helmholtz_bie_reference.m` | A=0.01→0.20 m：PE–BIE E_G 0.0032→0.081，Bellhop–BIE 8e-6→1e-4；**Bellhop 更近**（`reports/pe_bellhop_helmholtz_bie_final_report.md`） | **部分验证**（结论为负面：PE 屏粗糙度适用性受限） |
| 受控比较：差异机制 = 空间相关失真；canonical PM ∈ Region III | `validate_pe_bellhop_controlled_comparison.m` Stage 0–7 | 最终状态 PHASE_MECHANISM_IDENTIFIED（`reports/pe_bellhop_controlled_comparison_final_report.md`） | **部分验证** |
| PE 海面算子局部改进可行性（G0–G4） | `validate_pe_surface_normal_approximation_audit.m` 等 | G2 kz-aware 有改进；G3 失败；NONLOCAL_EFFECT_REQUIRED（`reports/pe_surface_operator_improvement_final_report.md`） | **部分验证** |
| PE–Bellhop PM 统计差异稳定性（Stage 3/4） | `validate_pe_bellhop_pm_amplitude_ensemble.m`、`validate_pe_bellhop_pm_model_discrepancy_statistics.m` | M=50：mean ΔTL -0.0075 dB；bootstrap 半宽 0.258>0.25 dB 门限 → PRELIMINARY（`reports/pe_bellhop_pm_model_discrepancy_statistics_report.md`） | **部分验证** |
| Li2009 粗糙海面旅行时间统计复现 | `validate_li2009_explicit_surface_vertical.m` | 首阈值时延标准差 +126% 与论文趋势一致；峰值展宽未复现；**横向网格收敛未通过**（`reports/li2009_explicit_surface_validation_report.md`） | **部分验证** |
| SSA1 / 2k 相位与海谱诊断 | `validate_2k_phase_approx.m`、`validate_2k_ocean_spectra.m` | 只验证相位系数近似与最低阶趋势 | **部分验证**（仓库自定位为诊断，非物理验收） |
| 气泡 Hall1D | `calibrate_hall1d_bubble_vertical.m` 等 | 部分复现；横向网格收敛与独立交叉验证未闭环（`vertical_comm_guide.md` §4） | **部分验证** |
| 随机海面 3–5 kHz 宽带成对收敛 / ensemble 最差群时延 | `validate_random_surface_window_robustness_wideband_vertical.m` | 运行被停止，部分 checkpoint 不作资格证据 | **尚未验证** |
| layered / Doppler / GPU / 多接收机 / 插值接收 / 气泡伴随统计 | 无验证入口 | 代码接口存在（如 `local_sound_speed` layered 分支、`rx_position_fn`/`doppler_fn` 参数），验证链明确不覆盖（指南 §9） | **尚未验证** |
| 风速插值条件模型 | `build_conditional_channel_library_vertical.m` 显式拒绝 | 设计上禁止 | **尚未验证**（by design） |
| Bellhop 绝对幅度 / 点源无限孔径 | `validate_bellhop_freefield_normalization_vertical.m` 等 | ~0.26 dB backward-range 偏移保留为 diagnostic 未拟合（多份报告） | **开放** |

---

## 10. 存在的问题

按风险从高到低排列：

### R1（根本性）：海面反射物理模型与整条统计链的基准耦合
整条 joint-kstat → 伴随统计 → 条件生成器 → 通信验证链都建立在 Kirchhoff 统计相位屏（δG=R0·e^{2ik0η}-R_coh）之上，而独立 Helmholtz BIE 仲裁表明该屏在 A≥0.05 m 时系统性偏离精确解（E_G 0.016→0.081，Bellhop 同条件 4e-5→1e-4），canonical PM 工况已入 Region III（`reports/pe_bellhop_controlled_comparison_final_report.md`）。数值实现本身可信（伴随、统计收缩均达双精度），**但物理基准的绝对准确性在高海况下不成立**——以此链产生的绝对 BER/时延统计会继承该偏差。仓库对此有清醒的自我声明（指南 §9 区分 PASS 与“仅诊断”），但尚未把适用域边界写进公共入口的输出语义。

### R2（工程）：公共默认配置 ≠ 验证收敛配置
公共默认 `xw=yw=50 m + sponge(0.12, 0.15)`（`vertical_channel_model_impl.m:122-140`）与严格验证推荐（直达 160 m、完整反射链 192.1875 m、均 no-sponge）不一致：50 m 窗的反射 TL 偏约 1.15 dB、边缘能量 11.8%（`reports/pe_reflected_chain_window_validation_report.md`）。文档反复警告，但代码未切换、输出中也无强制标注——新使用者极易用默认参数产生带系统偏差的结果。注意 192 m 窗在公共入口需 `validation_allow_extended_window=true`（窗口上限默认 100 m，`vertical_channel_model_impl.m:270-274`），即**推荐配置在公共路径上当前不可直达**。

### R3（通信）：无噪声 BER 底 0.2%–0.5%
循环 FFT 逆滤波近似线性卷积 + 帧无保护间隔 + 0.999 能量截断 + 固定 MMSE 正则，导致高 Eb/N0 的 BER 曲线存在算法底（`reports/two_node_statistical_channel_communication_validation_report.md` §7），两节点报告自己已把“改通信算法”排在“扩风速”之前。

### R4（证据规模）：统计证据仍是小样本
Stage 3/4 PE–Bellhop ensemble 仅 8–50 seed 且 bootstrap 半宽 0.258 dB 超 0.25 dB 门限（PRELIMINARY）；随机海面宽带成对收敛未完成，ensemble 最差群时延不可报告；条件模型只有 U=5/U=8 两点。

### R5（环境单一）：uniform-only
layered（含 erfc 测试剖面）、Doppler、GPU、多接收机、插值接收、气泡在伴随/统计链中全部不覆盖；`rx_position_fn`/`doppler_fn` 是空接口（`vertical_wape_propagator.m:22-61` 只在 t=0 求值一次）。

### R6（小问题 / 卫生）
- `examples/main_vertical_demo.m:25` 用 `z_rx=4`，与权威标准几何 z_rx=3（comm demo 与指南）不一致；
- 根目录残留 `psk_comm_result.mat`（2.1 MB 演示输出）、`A_0.25_N_161_S_0.prt`、`.prt`、`tmp/`，与 README“大型结果不放回根目录”的自定规则不完全一致；
- `scripts/README.md` 中大量入口已归档到 `cash/`，链接指向归档路径，活动/退役边界需要读者自行分辨（README 有声明，但认知负担高）。

---

## 11. 下一步研究建议

按“性价比 × 解除的风险”排序：

1. **关闭 R2（低成本、确定性收益）**：完成公共默认窗口向 192.1875 m/no-sponge 的迁移（或至少在 `output.config` 之外增加强制告警字段/日志），并放宽 `validation_allow_extended_window` 对推荐配置的可达性。这一步不改任何物理公式，只消除已量化的 1.15 dB 系统偏差。
2. **完成随机海面宽带成对收敛**（R4）：跑完 3–5 kHz 候选–256 m 参考成对集合，取得 ensemble 最差群时延上界——这是宽带通信结论（群时延、时延扩展）资格化的前置条件，且基础设施（checkpoint、resumable 框架）已就绪。
3. **通信基线修复**（R3）：引入循环前缀/overlap-save 线性卷积基线、帧保护间隔、自适应正则，把无噪声 BER 压到 0；随后高 SNR 结论才有意义。两节点报告 §8.6 已给出同样排序。
4. **非局域海面算子研究**（R1 的正面进攻，已锁定方向）：G4 结论 NONLOCAL_EFFECT_REQUIRED + BIE 参考在手，可 (a) 研究 kz-aware 谱域非局域屏（G2 已证局部 kz 修正方向正确），或 (b) 以 Bellhop/BIE 建立 Kirchhoff 屏的 (Hs, k0σ, 坡度) 适用域边界图（Stage 4 validity map 已有雏形），把统计链的物理声明显式限制在 Region I/II。
5. **以上完成后再扩展**：U=10 等新风速节点、宽带条件模型、layered 环境验证、实测数据对比。在 R1/R3 未解前扩展节点只会复制带偏差的基准。

---

## 附：四个重点问题的直接回答

**Q1 当前 UWA 项目真正实现了什么？**
一条数值实现质量很高的垂直信道仿真与统计研究链：均匀介质 PE/WAPE（精确角谱分裂）两段传播 + 4 种海面边界模型 + 气泡接口 + 精确伴随接收统计 + U=5/U=8 条件信道生成器 + MPSK 通信 Monte Carlo。代码防御性极强（约 100 项输入校验、全链路元数据、seed 可复现、缓存指纹防混用），内部数值一致性验证到双精度级别。物理上，直达路径与平面/弱粗糙海面有 Bellhop/BIE 外部基准支持；粗糙海面反射的 Kirchhoff 相位屏在高海况的绝对物理准确性被独立 BIE 参考否定，项目已自我识别并定位为非局域算子问题（§6、§9）。

**Q2 哪些结果有可靠验证？**
可靠（双精度级内部一致性 + 明确门限）：PE=角谱实现、direct_dsp 相位语义与 4 ms 名义时延、精确伴随/投影/C-P 收缩、192.1875 m 窗口收敛（固定面 CONVERGED、4 kHz 随机面 15/15）、U=5/U=8 条件生成器、两节点通信趋势复现（CI 级）、PE–Bellhop Stage 0 数值链。受限：PE Kirchhoff 粗糙面绝对物理（Region III）、宽带随机海面群时延（未完成）、高 SNR BER（算法底）、一切非 uniform 环境与风速插值（§9 证据表）。

**Q3 最大的技术风险是什么？**
海面反射物理基准与统计链的耦合（R1）：统计信道→通信结论的“物理正确性”上限由 Kirchhoff 统计屏决定，而该屏在目标 PM 海况被 BIE 裁定为系统性偏离（Bellhop 更近 2–3 个数量级）。若不先把适用域声明进输出语义或改进算子，基于该链的绝对通信性能结论会整体继承这一偏差。次级风险是默认窗口未迁移（1.15 dB）与通信 BER 底（0.2–0.5%）。

**Q4 下一步最值得推进的方向是什么？**
依次：(1) 默认窗口迁移（消除已知系统偏差，零物理风险）；(2) 随机海面宽带成对收敛（补齐宽带资格证据）；(3) 通信基线修复（CP/线性卷积均衡，消除误码底）；(4) 非局域海面算子或 Kirchhoff 屏适用域边界图（正面解决 R1）。风速节点扩展应排在这些之后。

---

*本报告基于 2026-09-23 的仓库快照（main @ 361bcec）。所有行号引用以该快照为准。*
