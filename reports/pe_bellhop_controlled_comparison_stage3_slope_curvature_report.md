# PE--Bellhop controlled comparison — Stage 3 slope/curvature sweep

日期：2026-09-12  
状态：**PASS**

## 冻结配置

Stage 2 的 Region-I 硬门同时由 `E_G`、`E_aligned`、TL RMS、phase RMS 和
`rho_shape` 决定。满足全部门限的最大高度为 `A=0.02 m`，因此 Stage-3
manifest 将其一次性冻结；没有根据 K-sweep 结果改选 A。其余固定项为 4 kHz、
uniform `c=1500 m/s`、X source、coherent C、10,001 beams、profile `N=4097`、
相同 flat denominator、AS footprint 和 Stage-1Y comparison convention。

## 结果

| K (rad/m) | max slope | max/RMS curvature (1/m) | min radius (m) | E_G | phase RMS (rad) | TL RMS (dB) | phi0 (rad) | rho_raw | rho_shape | E_aligned | hit max |kappa| (1/m) | hit-kappa max error (1/m) | wall residual (m) | min post dr (m) | region | guards |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---:|:---:|
| 0.10 | 0.0020 | 0.000200 / 0.000140143 | 5000 | 0.00823397 | 0.00816578 | 0.00914010 | -0.00270573 | 0.999966 | 0.999970 | 0.00777803 | 0.000199999 | 3.87e-7 | 7.10e-15 | 0.0433501 | I | PASS |
| 0.20 | 0.0040 | 0.000800 / 0.000570535 | 1250 | 0.00836932 | 0.00757942 | 0.0308524 | -0.00218459 | 0.999965 | 0.999967 | 0.00807911 | 0.000799989 | 3.11e-6 | 7.11e-15 | 0.0434707 | II | PASS |
| 0.35 | 0.0070 | 0.002450 / 0.00172429 | 408.163 | 0.0141398 | 0.00942753 | 0.0915312 | -0.00472523 | 0.999900 | 0.999911 | 0.0133263 | 0.00244990 | 1.66e-5 | 7.12e-15 | 0.0429552 | II | PASS |
| 0.47 | 0.0094 | 0.004418 / 0.00311980 | 226.347 | 0.0228161 | 0.0129059 | 0.163441 | -0.00891677 | 0.999740 | 0.999779 | 0.0210043 | 0.00441769 | 4.04e-5 | 7.12e-15 | 0.0428687 | II | PASS |

所有点的 pressure-release phase、`p/q` rotation、wall intersection、travel
time、receiver mapping 和 finite-state guards 均通过。每条射线的 hit-point
curvature、`Tg/Th`、`RN/RM`、reflection 前后 `p/q`、tau 和 wall residual 已保存在
各 case MAT 中，并保留原始 `.iwdiag`。随 K 增大出现的模型差异按 Goal 不触发
Stage-3 停止，也没有通过拟合、相位扣除或物理公式修改予以消除。

## 产物

- `results/validation/pe_bellhop_controlled_comparison/stage3_slope_curvature_sweep/stage3_manifest.mat`
- `results/validation/pe_bellhop_controlled_comparison/stage3_slope_curvature_sweep/stage3_validation.mat`
- `results/validation/pe_bellhop_controlled_comparison/stage3_slope_curvature_sweep/stage3_metrics.csv`
- `results/validation/pe_bellhop_controlled_comparison/stage3_slope_curvature_sweep/stage3_report.md`
- 各 `slope_A0p02_K*_N4097_B10001/` 下的 `stage3_case.mat` 与 case report

## 判定

```text
Stage 3 = PASS
Stage 4 = COMPLETED / PASS (subsequent read-only validity map)
```

本阶段随后已进入 Goal 的 Stage 4 read-only validity map；见
`reports/pe_bellhop_controlled_comparison_stage4_validity_map_report.md`。
仍未启动 fixed-PM 或任何 Monte Carlo。
