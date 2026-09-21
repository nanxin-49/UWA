# PE--Bellhop controlled rough-surface comparison Goal audit

审核日期：2026-09-10  
审核状态：**GOAL_MODE_READY_WITH_MINOR_CHANGES**

本轮仅完成源码、权威报告和可执行条件审查；没有运行新的 PE/Bellhop
case，没有修改 production PE、Bellhop `Reflect2D`、`InfluenceGeoHatCart`
或 SHD selector，也没有使用子 Agent。

## 1. 审核结论

新的 controlled comparison 可以执行，但必须先冻结以下小改动：

1. Bellhop comparator 固定为 2-D Cartesian line source `X`。历史 `R`
   的约 `0.26065 dB` range normalization 已被定位，不能重新引入。
2. 物理 PE receiver 坐标 `x` 在 post-reflection proper rotation chart 中
   映射为 Bellhop receiver depth `z'_BH=-x`。写入 ENV 前按 Bellhop 要求将
   depths 排序，读回后按保存的逆置换恢复 PE `x` 顺序；range 固定写
   `[102,103] m` 并显式选择 `103 m`。
3. `eta(x)=A cos(Kx)` 不使用只接受 `r=R0-A sin(Kz)` 的旧 sinusoidal
   validation runner。直接复用已通过参数化几何修复的 generic PM
   internal-wall binary，输入点列 `Gamma(s)=[R0-A cos(Ks),s]`。这不改变
   任何 Bellhop 核心反射或 influence 公式。
4. 5001-beam weak-sinusoid 场已有约 `3.16 dB` 的未收敛先例。新的第一条
   weak sinusoid 必须做 `10001 -> 20001` beams；失败时只允许一次
   `20001 -> 40001` escalation，仍失败则停止。
5. 相关性必须拆成复 coherence。`|c|` 对 global phase 不敏感，
   `arg(c)` 才是 global phase offset，`Re(c)` 是 phase-sensitive raw
   correlation；不能把 `|c|` 的 phase-aligned 前后变化当作归因证据。
6. 当前 Stage 2/3/4 数值报告保留历史 point-source `R` provenance。
   它们可用于数值量级和 M=50 机制解释，不能作为新的 `X` receiver-line
   结果，也不得重写标签。

完成上述冻结后，不存在需要重写核心物理的阻塞，故不是
`GOAL_MODE_BLOCKED`。

## 2. 当前真实状态

| 项目 | 当前证据 | 本 Goal 中的处理 |
|---|---|---|
| 1T PE propagation | `pe_1d_validation_bridge_report.md`：与 exact AS/production `ky=0` 达双精度量级 | 直接复用，不重验 marching |
| incident source mapping | `pe_bellhop_incident_field_comparison_report.md`：X-source，13 gates PASS；M99 L2 `0.00610485` | 作为 source/free-field floor；不重拟合 `.sbp` |
| Bellhop internal wall | `bellhop_internal_wall_implementation_report.md`：参数化 `t/n/kappa`、single `Reflect2D`、proper rotation、receiver selector 已通过 | 直接复用 generic parametric wall binary |
| R/X source audit | `bellhop_source_geometry_rx_audit_report.md`：X 将 native/internal `-0.26065 dB` 降至 `4.24e-6 dB` | 所有新 case 固定 X |
| fixed-PM numerical budget | `pe_bellhop_pm_numerical_error_budget.md`：PE/BH 扫描 PASS | 保留为先验预算；因其为 pre-X/轴上证据，不替代新 spatial convergence |
| fixed-PM Tier-1 | X seed-260001 为 `0.310371 dB / -2.248201 rad` | 仅作新 Stage 6 轴上 sanity check |
| dimensionality | 1T->2T 为 `-0.15886 dB / +0.01251 rad` | 报告但不从新结果中数值相减 |
| 4/6/8 kHz、M=50 | 现存正式值为历史 R provenance | 不重跑；只在最终解释阶段引用并标注限制 |

现有证据没有回答新的核心问题：`A -> 0` 时，PE Kirchhoff phase screen 与
Bellhop local-specular reflected **receiver-line complex field** 是否回到共同
flat/numerical floor。因此旧 fixed-PM、frequency 和 ensemble 计算不能代替
本 Goal 的 Stage 0--1。

## 3. 已冻结的可比性定义

### 3.1 坐标与 receiver

- PE validation geometry：source-to-nominal-surface distance `100 m`，
  surface-to-receiver distance `3 m`，flat image distance `103 m`。
- Bellhop incident chart：source `(r,z)=(0,0)`，wall
  `Gamma(s)=[100-eta(s),s]`。
- wall hit 后仅执行一次原生 `Reflect2D`，再作
  `T(r,z)=(200-r,-z)`；`p/q/tau/Amp/Phase` 不变。
- physical receiver `(x,3 m)` 对应 mapped receiver `(r',z')=(103,-x)`。
- Bellhop receiver depths 必须排序写入，结果必须通过明确的坐标表恢复到
  PE x-grid；禁止通过误差最小化选择正负号、range 或 column。

### 3.2 source 与相位

- 4 kHz，`c=1500 m/s`，Gaussian `sigma=0.3 m`。
- `.sbp` 固定为
  `D(theta)=cos(theta) exp[-(k sigma sin(theta))^2/2]`，2401 samples，
  `[-30,30] deg`，`-120 dB` floor。
- Bellhop run type 固定 coherent `C`，source geometry 固定 `X`。
- Bellhop raw/conjugate 使用既有 normalization audit 的固定 spatial sign；
  不允许逐 case 选择、不应用历史 point-source global constant。
- PE/Bellhop flat-normalized ratio 自动消除共同 carrier；不得额外加入或删去
  carrier phase。

### 3.3 flat 与 rough 定义

对于同一 receiver-line sample：

```text
G_PE(x) = Href_PE,rough(x) / Href_PE,flat(x)
G_BH(x) = Href_BH,rough(x) / Href_BH,flat(x)
```

PE flat/rough 必须共享全部非 surface 配置；Bellhop flat/rough 必须使用同一
generic internal-wall executable、source、beam fan、beam count、step、receiver
array 和 SHD selector。flat denominator 不能来自不同 source type 或旧 cache。

95%/99% footprints 只由 Stage-0 flat exact-AS field at `103 m` 构造，之后冻结；
rough result 不得参与 mask 选择。pointwise phase 只在 M95 且双方 flat
denominator 均高于 axis `-40 dB` 的样本上统计。

## 4. 删除的重复或无效工作

后续不再执行：

- PE marching、1T bridge、internal-wall `Reflect2D`/curvature/rotation 的重复审计；
- 89/89.5/89.9-degree native ATI 对照；本 Goal 直接使用 exact-90 internal wall；
- 6/8 kHz 扩展、更多 PM seeds 或新的 Monte Carlo；
- 旧 `R` cases 的重新解释、经验 amplitude correction 或 `.sbp` refit；
- total field、communications、arrivals/eigenray ordering；
- 在 weak-limit 失败后继续增大 A/K。

## 5. 数值 floor 与预注册门槛

已有 X incident audit 的量级为：M99 complex L2 `0.00610485`、M99 phase RMS
`0.00610483 rad`、M95 phase P95 `0.0142931 rad`、M95 TL P95
`0.00190364 dB`。这些是先验量级，不直接替代 reflected Stage-0 measurement。

Stage 0 的 hard gates 冻结为：

| check | limit |
|---|---:|
| PE flat reflected vs exact AS, normalized complex error | `1e-10` |
| PE outer-5% reflected energy | `1e-5` |
| receiver coordinate residual | `1e-6 m` |
| BH generic-flat vs official/free-field `-P_X(103 m)`, M99 L2 | `2e-3` |
| 5001->10001 beam M99 L2 | `2e-3` |
| beam M95 phase/TL RMS | `5e-3 rad / 0.02 dB` |
| main PE--BH flat M99 normalized L2 | `0.02` |
| main PE--BH M99 phase RMS | `0.02 rad` |
| main M95 phase/TL absolute P95 | `0.05 rad / 0.10 dB` |
| complex coherence magnitude `|c|` | `>=0.9995` |
| global phase `|arg(c)|` | `<=0.02 rad` |
| half-dx shared-node M99 L2 | `0.005` |
| NaN/Inf | `0` |

Stage 0 通过后记录而不回调阈值：

```text
F_E   = max(flat PE-BH L2, BH-chain L2, beam L2, half-dx L2)
F_phi = max(flat PE-BH phase RMS, beam phase RMS, half-dx phase RMS)
F_TL  = max(flat PE-BH TL RMS, beam TL RMS, half-dx TL RMS)
```

weak-limit acceptance bands 在查看 sinusoid 结果前固定为：

```text
T_E   = F_E   + max(0.005, 0.5 F_E)
T_phi = F_phi + max(0.005 rad, 0.5 F_phi)
T_TL  = F_TL  + max(0.01 dB, 0.5 F_TL)
```

并始终要求 `|c|>=0.9995`。这些量只用于判定是否回到 comparison floor，
不能用来校准任何场。

## 6. 最小 preparation

当前 worktree 不含 git-ignored validation binaries/results，但原项目目录中的
官方 Bellhop 2020、canonical coefficients 和 M=50 artifacts 可用。Goal 执行前：

1. 从当前源码构建 flat 和 generic PM validation binaries到本 Goal 的 result
   tree；不要复制旧 case cache。
2. 仅在 Stage 6 到达后复制 authoritative canonical coefficient/master CSV，
   并核对 coefficient SHA-256
   `1f7eda465e4ae85b8ac038310edf053f2d062563a3b943bfd83bd59d07015f67`。
3. Stage 7 只读既有 M=50 CSV/MAT，记录来源和 SHA-256；不重跑 seeds。
4. 所有新 case root 必须包含 source `X`、A、K、profile N、beam count 和 step，
   cache fingerprint 至少覆盖这些字段及 executable/profile hashes。

这些属于可重复的 validation preparation，不改变研究模型，故最终判定为：

```text
GOAL_MODE_READY_WITH_MINOR_CHANGES
```

后续唯一执行规范为
`reports/pe_bellhop_controlled_comparison_GOAL.md`；旧设计报告继续作为历史证据，
不再决定阶段顺序或验收门槛。
