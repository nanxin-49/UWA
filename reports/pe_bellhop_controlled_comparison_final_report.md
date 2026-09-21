# PE--Bellhop controlled comparison — final report

日期：2026-09-13  
最终状态：**PHASE_MECHANISM_IDENTIFIED**

## Executive conclusion

Stage 0–7 已按冻结顺序完成。PE 与 Bellhop 在共同 weak-roughness limit 中具有
一致的主导 reflected complex-field behavior；随着粗糙高度和 slope/curvature
增加，差异表现为空间相关失真，而不是未处理的统一 complex convention、receiver
selection 或 Bellhop internal-wall 数值失败。canonical fixed-PM case 已进入
Region III，因此本研究不能据此宣称 PE Kirchhoff phase-screen 在目标 PM 工况下
具有绝对物理准确性。

## 1. Flat spatial comparison floor

Stage-0 flat receiver-line closure 为：

| metric | value |
|---|---:|
| normalized complex L2 `F_E` | 0.00589939216 |
| phase RMS `F_phi` | 0.00589924793 rad |
| TL RMS `F_TL` | 0.000797852579 dB |
| rho_shape | 0.999996312 |
| global phase | -0.00197233300 rad |
| aligned L2 | 0.00271602777 |

对应冻结门限为 `T_E=0.0108993922`、`T_phi=0.0108992479 rad`、
`T_TL=0.0107978526 dB`。internal-wall 与独立 free-field chain 的误差为零到
双精度量级；10,001-beam endpoint 的 beam-convergence L2 为 `5.85e-6`。

## 2–3. Fixed conjugation and theoretical closure

Bellhop SHD 输出与项目正向空间相量的绝对场定义不同。冻结定义为：

```text
B_raw           = SHD direct complex pressure
B_abs           = conj(B_raw)
G_BH_abs        = B_abs,rough / B_abs,flat
G_BH_comparison = conj(G_BH_abs) = G_BH_raw
```

这不是逐 case 择优或相位拟合。Stage 1Y 从 PE 的 `exp(+i(kz-k0)L)`、Kirchhoff
screen `-exp(+i2k eta)`、Bellhop coherent influence
`exp(-i(omega*tau-phase))`、SHD real/imag reader、X-source 实数归一化、单次
vacuum `+pi` 和 proper half-turn 分别闭合了符号，结论为
`CONVENTION_THEORETICALLY_CLOSED`。

## 4–5. Weak-limit convergence

固定 `K=0.10 rad/m` 时，convention-fixed 结果随 A 减小为：

| A (m) | E_G | phase RMS (rad) | TL RMS (dB) | rho_shape |
|---:|---:|---:|---:|---:|
| 0.0100 | 0.00409346 | 0.00405930 | 0.00456998 | 0.99999242 |
| 0.0050 | 0.00204123 | 0.00202407 | 0.00228491 | 0.99999810 |
| 0.0025 | 0.00101951 | 0.00101096 | 0.00114248 | 0.99999952 |

因此 `A->0` 时误差近似线性下降并趋向共同 flat/source floor；weak sinusoid
`A=0.01 m` 满足全部冻结门限。

## 6. Height, slope and curvature dependence

Stage 2 固定 K 增大 A 时，`E_G` 从 `0.00409` 单调增至 `0.10049`，说明差异
首先随 height/phase excursion 清晰增长。Stage 3 固定 `A=0.02 m` 增大 K 时，
`E_G=0.00823/0.00837/0.01414/0.02282`，而 TL RMS 对 K 更敏感并增至
`0.16344 dB`。由于改变 K 会同时改变 slope 与 curvature，本数据只能识别
“slope/curvature 耦合效应”，不能把二者因果贡献再单独拆分。

## 7. Phase mechanism

三个 Region-II 代表点的 global-phase alignment 仅降低
`7.85%/20.24%/7.94%` 的 `E_G`，均按冻结规则归类为 **spatial distortion**。
Bellhop 相对解析 stationary-path phase 的 RMS 仅
`1.7e-6/6.2e-6/1.3e-5 rad`，PE 对应残差为
`0.02087/0.10014/0.01290 rad`。因此剩余相位差不是统一 offset，主要来源是
Kirchhoff phase-screen 与 local-specular stationary/`Reflect2D` 的空间模型差异。

## 8. Sampled separation intervals

- 固定 `K=0.10 rad/m`：Region-I 到非 I 的采样括区为
  `A=(0.02,0.05] m`。
- 固定 `A=0.02 m`：括区为 `K=(0.10,0.20] rad/m`。

这些是 sampled brackets，不是插值或外推得到的精确 validity boundary。

## 9. Fixed-PM validity region

seed-260001、4 kHz、X/C、10,001-beam fixed-PM 的完整 M99 指标为：

```text
E_G       = 0.958268507
TL_RMS    = 1.27399043 dB
phase_RMS = 1.10771425 rad
phi0      = -0.798545974 rad
rho_shape = 0.771212836
E_aligned = 0.678503463
```

其分类为 **Region III**。同时 wall residual、incidence/grazing、curvature、
pressure-release phase、`p/q`、tau、positive post-range 和 field-finite guards
全部通过。轴上 sanity 为 `0.310370739 dB / -2.24820073 rad`，但轴上点不能代替
完整 receiver-line 判断。

## 10. Historical M=50 interpretation

M=50 的 provenance 明确为 **historical Bellhop point-source R**；没有重命名为
X，也没有经验 R→X 修正。原统计为 PE/BH mean powers
`1.004510/1.007007`、mean delta TL `-0.00747858 dB`、circular mean phase
`-0.375345 rad`。mean-power agreement 可以由 realization-level 正负幅度差平均
抵消解释，并不意味着复场闭合；非零圆均相位与 Stage-5 spatial distortion 以及
Stage-6 Region-III 结果一致，而不是剩余 convention 错误。

## 11. Is a BIE/BEM reference required?

对当前 Goal 的 comparability、weak-limit closure 和 phase-mechanism 定位，条件式
BIE/BEM branch **不需要触发**：Stage1X 的 fixed convention 已由 Stage1Y 理论证明。
若下一步要裁决 Region-III PM 工况中“PE 或 Bellhop 哪个更接近精确 Helmholtz 解”，
则需要独立 BIE/BEM reference；应先做 flat、weak sinusoid、one stronger sinusoid，
不能直接进入 PM ensemble。

## 12. Validity of the target PM PE model

当前证据支持：PE marching/FFT/source 链数值自洽，且在 weak roughness limit 与
Bellhop 主导复场行为一致。当前证据**不支持**宣称目标 fixed-PM 下的 PE rough
surface model 已被独立验证为物理准确，因为该点处于 Region III，且没有 exact
Helmholtz reference。它应继续被表述为一个稳定但模型适用性尚待独立裁决的
Kirchhoff phase-screen 近似。

## Final status

```text
PHASE_MECHANISM_IDENTIFIED
```

`COMMON_LIMIT_CONFIRMED` 的全部必要证据同时满足；最终采用更具信息量的
`PHASE_MECHANISM_IDENTIFIED`，因为 Stage5 已将剩余差异定位为空间反射模型失真。
未修改 PE/Bellhop 核心物理，未运行 20,001/40,001-beam case，未扩展 PM
Monte Carlo。
