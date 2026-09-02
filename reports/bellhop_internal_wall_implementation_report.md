# Bellhop 2020 参数化 internal-wall 实现与验证总报告

## 1. 文档状态与结论

本报告是当前 Bellhop 2020 validation-only internal-wall 路线的唯一权威总报告，整合并取代 feasibility、flat、tilted、sinusoidal、beam/frame audit、旧 PM 尝试、PM 重设计审查、参数化源码审计和几何修复回归等阶段报告。

最终状态：**PASS_WITH_LIMITS**。

- flat、tilted、smooth sinusoidal 和严格 vertical-tangent curved wall 已通过。
- curved internal wall 已从 rotated chart 的 `z(r)`、`dz/dr`、`Dss`、`Delta z/Delta r`、`1/Delta r` 和 range 单调性假设中解耦。
- 物理反射只调用一次 Bellhop 2020 原生 2D `Reflect2D`；`InfluenceGeoHatCart`、Bellhop 原生 `p/q` 公式、PE 和通信代码均未修改。
- 严格 90 度下的旧 native ATI 与 rotated returned branch 不是同一个物理问题，因此旧 PM hard comparison 已废止。fixed-band-limited PM 的有限角度 Bellhop-only POC 尚未实现，不能据此宣称 PM rough wall 已通过。

本报告只覆盖独立 validation binary。官方 Bellhop 2020 executable 及其源码基线未被替换。

### 1.1 被整合的阶段报告

下表记录本报告吸收的全部历史报告及其最终处置。文件已归档，不应再被单独引用为当前结论。

| 归档报告 | 被吸收的有效结论 | 当前处置 |
|---|---|---|
| `bellhop_internal_rough_wall_feasibility_review.md` | 原生 `Reflect2D` + 反射后 proper rotation 的 validation-only 路线在限定范围内可行。 | 由已完成 POC 取代。 |
| `bellhop_internal_flat_wall_poc_report.md` | flat wall 精确恢复 `-P_BH(103 m)`。 | 保留为本报告的 flat 回归结论。 |
| `bellhop_internal_tilted_wall_poc_report.md` | tilted geometry/state 通过；backward-range amplitude 只作 diagnostic。 | 保留为本报告的 tilted 回归结论。 |
| `bellhop_internal_sinusoidal_wall_poc_report.md` | nonzero-curvature 反射链通过；早期 coherent phase 解释需后续审计修正。 | 被 frame audit 和几何修复结果覆盖。 |
| `bellhop_curved_wall_beam_frame_audit_report.md` | `2.095 rad` 来自 SHD receiver-column 误配，不是 frame covariance 缺陷。 | 结论并入本报告；一次性 logging 源码退役。 |
| `bellhop_internal_pm_wall_validation_report.md` | strict-90 native/rotated hard comparison 不是等价物理问题。 | 该 PM 方案与入口退役。 |
| `bellhop_pm_wall_validation_redesign_review.md` | 下一步应使用同一 fixed-band-limited realization 和有限角度全系统旋转。 | 作为下一阶段最小 POC 设计并入。 |
| `bellhop_internal_wall_parameterization_source_audit_report.md` | pre-fix curved adapter 仍含 range-based curvature，结论为需要修复。 | 问题已修复，审计被回归报告取代。 |
| `bellhop_internal_wall_geometry_fix_regression_report.md` | 参数曲率、单位 frame、finite support 和 vertical-tangent 回归通过。 | 构成本报告当前 `PASS_WITH_LIMITS` 的主要证据。 |

## 2. 当前实现

当前链路为：

```text
sampled parametric wall Gamma(s)=[r(s),z(s)]          [custom]
  -> finite-support segment intersection             [custom]
  -> ReduceStep2D stops at the true wall hit          [native hook + custom limiter]
  -> unit segment/node/hit tangent                    [custom]
  -> TOP-semantic unit normal n=(t_z,-t_r)            [custom]
  -> signed turning-angle / arc-length curvature      [custom]
  -> native Reflect2D once                            [Bellhop 2020 native]
  -> native RN/RM and p/q reflection update           [Bellhop 2020 native]
  -> one pressure-release pi phase change             [Bellhop 2020 native]
  -> proper half-turn: r'=2R0-r, z'=-z, t'=-t         [custom chart map]
  -> transformed post-wall branch isolation           [custom]
  -> InfluenceGeoHatCart                              [Bellhop 2020 native, unchanged]
```

### 2.1 参数化几何

墙由有序、相邻点不重合的参数点列定义。曲率不再先构造 `dz/dr`，而是由相邻单位切向的有符号转角除以局部弧长尺度得到：

```text
t_i       = Delta Gamma_i / |Delta Gamma_i|
Delta phi = atan2(cross(t_i,t_{i+1}), dot(t_i,t_{i+1}))
kappa     = Delta phi / Delta s
```

命中点切向插值后重新单位化；法向由单位切向按 TOP 标架一次性构造，不再独立插值。输入允许 `dr/ds>0`、`dr/ds=0` 和 `dr/ds<0`，只要求参数顺序明确、相邻段非退化且曲线局部正则。命中真实 profile support 以外的射线明确拒绝，不采用 native ATI 风格的虚假无限延拓。

### 2.2 反射后 proper rotation

固定变换

```text
r' = 2 R0 - r
z' = -z
```

是 determinant `+1` 的 180 度刚性旋转，只改变反射后分支的坐标、切向和射线法向表示。它不重新计算墙几何或 curvature，不再次调用 `Reflect2D`，不再次施加 reflection coefficient / pressure-release phase，也不修改已保存的 `p/q`、`tau`、`Amp`、`Phase` 或 KMAH 状态。

## 3. 源码布局

保留的 current validation 入口与 overlay：

| 目的 | 入口/源码 |
|---|---|
| flat wall 回归 | `scripts/validation/validate_bellhop_internal_flat_wall_poc.m`; `scripts/validation/support/bellhop_internal_flat_wall_poc/` |
| tilted wall 回归 | `scripts/validation/validate_bellhop_internal_tilted_wall_poc.m`; `scripts/validation/support/bellhop_internal_tilted_wall_poc/` |
| sinusoidal curved-wall 回归 | `scripts/validation/validate_bellhop_internal_sinusoidal_wall_poc.m`; `scripts/validation/support/bellhop_internal_sinusoidal_wall_poc/` |
| vertical-tangent curved-wall 回归 | `scripts/validation/validate_bellhop_internal_vertical_tangent_poc.m` |
| 通用 parametric/未来 PM overlay | `scripts/validation/support/run_bellhop_internal_pm_wall_poc_vertical.m`; `scripts/validation/support/bellhop_internal_pm_wall_poc/` |
| SHD receiver-range 配对回归 | `scripts/validation/validate_bellhop_shd_receiver_range_pairing_vertical.m` |

参数化曲率修复位于三个 curved validation overlay 的 `Step.f90`：

- `scripts/validation/support/bellhop_internal_sinusoidal_wall_poc/Step.f90`
- `scripts/validation/support/bellhop_internal_pm_wall_poc/Step.f90`
- historical beam/frame overlay 曾同步应用相同修复；该一次性插桩现已归档，不再是 current source。

## 4. 验证结果

### 4.1 回归总览

| 验证 | 结果 | 关键数值与解释 |
|---|---:|---|
| flat `r=100 m` | PASS | 4 kHz，step `0.2/0.1/0.05 m`、beam `2001/5001/10001`；全部恢复 `H_wall=-P_BH(103 m)`；最大墙残差 `3.54e-12 m`，最大 travel-time error `9.71e-16 s`，rotation 前后 `p/q` 误差为 0。 |
| tilted `r=100+0.005z` | PASS | 交点约 `1e-11 m`、specular direction error `6.7e-16`；一次 pi phase、`Amp` 与 `p/q` 无人工修正；native backward-range 场约 `0.286 dB` 幅度差仅作 diagnostic。 |
| smooth sinusoidal | PASS | 36-case step/beam/profile scan；最大交点残差 `3.54e-12 m`，tangent/normal error `3.64e-6`，direction error `8.89e-16`，正确 receiver pairing 后 phase difference 不超过 `6.71e-5 rad`。 |
| vertical tangent, circular wall | PASS | radius `100 m`，samples `65/129/257`；命中点严格 `dr/ds=0`，`u dot n=1`；无 NaN/Inf，signed kappa 从 `-0.010000456` 收敛至 `-0.010000028 1/m`，RN error 从 `6.08e-10` 收敛至 `3.80e-11`。 |

### 4.2 vertical-tangent 关键结论

圆弧命中处墙切向严格竖直，但曲线正则且曲率有限。三档采样均得到单位正交 frame、有限 curvature、正确镜面方向、稳定 `p/q` kick、一次 pressure-release phase change 和正向 transformed range。由此证明 `dr/ds=0` 本身不是当前 internal-wall adapter 的数值或理论奇异点。

真正仍需警惕的是物理 grazing：`u dot n -> 0`。它会使 Bellhop 原生 reflection/beam geometry 敏感，与墙在坐标图中是否近竖直是两个不同问题。

### 4.3 beam/frame 与 receiver-column 审计结论

曾报告的约 `-2.095 rad` native/rotated phase offset 不是 reflection-frame covariance 错误，而是把 rotated SHD 的 `102 m` 首列误当成 `103 m` 目标 receiver。固定显式 receiver-range selector，并对 native total/direct 使用同一列后，反射场 phase difference 约为 `1.02e-5 rad`。`Reflect2D` 出口的 `RN/RM/p/q/frame` 与 proper rotation 保存的 beam state 没有发现需修正之处。

native backward-range Cartesian influence 仍有约 `-0.2606` 至 `-0.286 dB` 的稳定幅度差。它保持为 point-source/range-chart diagnostic，不作拟合、renormalization 或 hard gate。

## 5. PM 阶段的当前状态

旧 strict-90-degree 方案把 native `z=eta(r)` 与 rotated internal-wall returned branch 当作等价场景，但两者在 native range chart 中不是同一个物理问题；同时旧 density sweep 会随 `N` 改变 PM realization 的高波数内容。该方案及其 runner 已退役，不能作为当前失败或通过证据。

下一步最小 POC 应当：

1. 只生成一次 fixed-seed、fixed-`Kmax`、高分辨率 band-limited `eta_ref`；所有 density cases 均从它采样/插值。
2. 先用全系统 proper rotation `89 deg`，统一变换 Tx、Rx、source direction 和同一条 rough profile，使 native ATI 与 internal wall 表示完全相同的物理几何；通过后再测试 `89.5 deg`，`89.9 deg` 只作极限诊断。
3. 推荐起始参数：seed `260001`、`U=6 m/s`、profile length `160 m`、reference samples `4097`、candidate `Kmax=0.5 rad/m`（须先通过 slope/curvature/grazing preflight）、density `513/1025/2049`、ray step `0.1 m`、beam count `2001`、4 kHz、显式 line source `X`。
4. 先验收 intersection、`t/n/kappa`、direction、RN/RM、`p/q`、delay 和 reflection phase，再查看 receiver complex field。不得进入 PE 或 Monte Carlo，直到 Bellhop-only covariance 通过。

需要特别区分：internal adapter 使用参数曲线的连续几何曲率；native C-ATI 的有限采样实现使用 Bellhop 自身的 `Dss` 离散。有限 `N` 下不应要求两者 bitwise 相同，而应检验它们对同一连续曲线的收敛。若未来需要严格的“同一离散 ATI frame”审计，应作为单独 diagnostic mode，不应重新把通用 internal wall 绑定到 `z(r)`。

## 6. 已知限制

- 仅 Bellhop 2020、2D、uniform `c=1500 m/s`、pressure-release、一次 wall reflection、当前 Gaussian `.sbp`。
- internal wall 只有有限 profile support；profile 外命中被拒绝。
- sidecar 当前仍按明确的采样参数顺序输入；允许 range 折返，但不允许退化零长度 segment。
- 未验证 layered SSP、3D、动态海面、多次粗糙面反射、跨 seam 的 arrivals/eigenrays 或大规模随机统计。
- `InfluenceGeoHatCart` 核心公式没有修改；native backward-range absolute amplitude 仍是 diagnostic limitation。
- fixed-realization PM、PE rough-wall cross-validation 和 Monte Carlo 均未完成。

这些限制不会否定已完成的参数化几何和 vertical-tangent 结论，但会阻塞将当前结果直接外推为完整 PM/PE 粗糙面模型。

## 7. 归档说明

本次收尾将以下阶段性材料移入 `cash/bellhop_internal_wall_superseded_20260902/`，不再作为 active context：

- 九份被本报告取代的 feasibility/POC/audit/PM/redesign/geometry-fix 阶段报告；
- 已完成使命的 curved-wall beam/frame logging 入口、runner 和源码 overlay；
- 物理比较关系无效的 strict-90-degree PM validator、native PM runner 和旧 profile generator。

归档不删除其历史内容；current 数值证据仍保留在 `results/validation/`。flat、tilted、sinusoidal、vertical-tangent、SHD selector 以及通用 parametric-wall build/runner 均保留为有效回归或下一阶段基础。

## 8. 最终判断

**PASS_WITH_LIMITS**：当前 validation-only curved internal wall 已真正按参数曲线处理，`dr/ds=0` 回归和 flat/tilted/sinusoidal 回归通过；原生 `Reflect2D`、原生 beam-state 更新和未修改的 `InfluenceGeoHatCart` 保持在链路中。下一步可以开展“同一 fixed-band-limited PM realization + 有限角度全系统旋转”的 Bellhop-only POC，但尚不能进入 PE rough-wall 对比或宣称 PM validation 完成。
