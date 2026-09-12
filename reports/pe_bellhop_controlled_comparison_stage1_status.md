# PE--Bellhop controlled comparison — Stage 1 status

日期：2026-09-10  
权威执行入口：`scripts/validation/validate_pe_bellhop_controlled_comparison.m`  
权威 Goal：`reports/pe_bellhop_controlled_comparison_GOAL_revised_stage1x.md`

## 已执行配置

4 kHz、均匀 `c=1500 m/s`、X 线源、coherent C、固定弱正弦
`eta(s)=A cos(0.1 s)`。Stage 0 已 PASS。根据用户最新运行约束，后续
Bellhop case 统一使用 10,001 beams；不再启动 20,001-beam。

## Stage 1 结果

| A (m) | profile L2 (M99) | beam L2 (M99) | model L2 (M99) | model phase RMS (rad) | model TL RMS (dB) | rho | 状态 |
|---:|---:|---:|---:|---:|---:|---:|:---|
| 0.0100 | 1.2025e-6 | 1.5437e-6 | 0.47247 | 0.47920 | 0.004570 | 0.91923 | FAIL |
| 0.0050 | 6.0272e-7 | 1.4363e-6 | 0.23875 | 0.23960 | 0.002285 | 0.97928 | FAIL |
| 0.0025 | 3.0923e-7 | 0 (same 10,001-beam endpoint) | 0.11969 | 0.11980 | 0.001142 | 0.99479 | FAIL |

三组 case 的 profile/geometry guards 均通过；wall residual 约 `7.1e-15 m`，
pressure-release phase jump 误差约 `7.1e-15 rad`，最小 transformed `dr`
约 `0.0433 m`，center travel-time 误差约 `4.7e-15 s`，且无 NaN/Inf。

Stage-0 floor 给出的弱限值为 `T_E=0.0108994`、`T_phi=0.0108992 rad`、
`T_TL=0.0107979 dB`。即使 `A=0.0025 m`，`E_G`、phase RMS 和 `rho_shape`
仍未达到 floor gate；TL RMS 已通过。随着 A 减半，模型 L2/phase 近似减半，
数值收敛项保持在 `1e-6` 量级，说明当前阻塞是模型可比性而非 solver
收敛或 internal-wall 几何故障。

## 决策

按修订版 Goal：Stage 2--7 仍锁定，先执行只读 Stage 1X
phase/convention audit，再决定是否允许后续 Helmholtz/BEM branch。Stage 1X
报告为 `reports/pe_bellhop_controlled_comparison_stage1x_phase_convention_audit.md`；
其固定共轭诊断闭合三组 A，但不自动改写 Stage 1 的原始 PASS/FAIL。
A=0.01 与 A=0.005 的历史 20,001-beam endpoint 结果保留为诊断；A=0.0025
遵循最新约束仅使用 10,001 beams。

详细逐 case 文件：

- `results/validation/pe_bellhop_controlled_comparison/stage1/A0p01/`
- `results/validation/pe_bellhop_controlled_comparison/stage1/A0p005/`
- `results/validation/pe_bellhop_controlled_comparison/stage1/stage1_report.md`

本轮未修改 PE、Bellhop 核心、`Reflect2D`、`InfluenceGeoHatCart`、SHD
selector 或通信代码。

## Convention-only 回归（Stage 1X 后）

按修订 Goal 对 Stage 1X 的处理，主驱动只在比较层增加
`G_BH_comparison=conj(G_BH)`；Bellhop/PE 原始场和旧结果均保留不变。
Stage 0 已重新执行并继续 PASS。A=0.01 的 10,001-beam 主程序回归 PASS；
A=0.005 与 A=0.0025 使用已保存的同一 raw solver MAT 做 comparison-only
重建，未重新启动 Bellhop，且各自的单项 checks 重新计算后均 PASS。

| A (m) | corrected E_G | corrected phase RMS (rad) | corrected TL RMS (dB) | rho | 状态 |
|---:|---:|---:|---:|---:|:---|
| 0.0100 | 0.00409346 | 0.00405930 | 0.00456998 | 0.99999242 | PASS |
| 0.0050 | 0.00204123 | 0.00202407 | 0.00228491 | 0.99999810 | PASS |
| 0.0025 | 0.00101951 | 0.00101096 | 0.00114248 | 0.99999952 | PASS |

fixed 输出位于
`results/validation/pe_bellhop_controlled_comparison/stage1_convention_fixed/`；
原始 Stage 1 FAIL 结果位于 `stage1/`，作为 superseded-by-convention-fix
诊断保留。由于后两组是 comparison-only 重建，本节不宣称新增 Bellhop
solver 收敛证据；Stage 2--7 仍需按修订 Goal 另行解锁。

## Stage 1Y theoretical convention closure

Stage 1Y 已完成，且没有重新运行 PE/Bellhop。报告
`reports/pe_bellhop_controlled_comparison_stage1y_theoretical_convention_closure.md`
从 PE helper、Bellhop 2020 `influence.f90`/`ReflectMod.f90`、SHD reader、X
source normalization 以及既有 free-field 和 constant-height 数据独立闭合了
复数传播约定。正式比较层固定为：

```text
B_abs = conj(B_raw)
G_BH_comparison = conj(G_BH_abs) = G_BH_raw
```

这不是按误差择优，也不是幅相校正；PE screen、Bellhop `Reflect2D`、
`InfluenceGeoHatCart`、SHD selector 和通信代码均未修改。`eta=+0.05,0,-0.05 m`
的既有 constant-height audit 给出 Bellhop ratio phase residual 约
`2.56e-5 rad`，与理论 `exp(+i*2*k*eta)` 符号一致；incident PE--AS L2 为
`3.73e-13`，Bellhop--AS 10001-beam L2 为 `6.11e-3`，均为既有结果。

最终状态：

```text
Stage 1 raw comparison             = FAIL (永久保留)
Stage 1 convention-fixed regression= PASS
Stage 1Y                           = CONVENTION_THEORETICALLY_CLOSED
Stage 2                           = UNLOCKED at Stage1Y completion (see follow-up below)
```

## Stage 2 follow-up (2026-09-11)

Stage 2 fixed-`K` height sweep has since completed with X/C, `N=4097`, and the
user-mandated 10,001 beams. The authoritative report is
`reports/pe_bellhop_controlled_comparison_stage2_height_sweep_report.md`;
all geometry/finite/metric guards PASS and Stage 3 is unlocked. The measured
M99 `E_G` values for `A=[0.01,0.02,0.05,0.10,0.20] m` are
`[0.004093,0.008234,0.021036,0.044179,0.100485]`; corresponding phase RMS
values are `[0.004059,0.008166,0.020869,0.043877,0.100142] rad`.
Model discrepancy increases with height, but this is not a Stage-2 solver or
mapping failure under the Goal rules. Stage 3 has not been executed.

## Stage 3 follow-up (2026-09-12)

Stage 3 subsequently froze `A=0.02 m` as the largest Stage-2 Region-I height
and completed `K=[0.10,0.20,0.35,0.47] rad/m` at 10,001 beams. All hard
solver/mapping/geometry/finite guards pass; the sampled classifications are
`I/II/II/II`. See
`reports/pe_bellhop_controlled_comparison_stage3_slope_curvature_report.md`.
Stage 4 has not been executed.

## Stage 4 follow-up (2026-09-12)

The read-only Stage-4 validity map subsequently completed with zero solver
calls. It reports sampled Region-I transition brackets `A=(0.02,0.05] m` at
`K=0.10 rad/m` and `K=(0.10,0.20] rad/m` at `A=0.02 m`; no boundary is
interpolated or extrapolated. See
`reports/pe_bellhop_controlled_comparison_stage4_validity_map_report.md`.
