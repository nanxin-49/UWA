# PE--Bellhop controlled comparison — Stage 5 phase attribution

日期：2026-09-12  
状态：**PASS**  
solver calls：`0`

## 方法

从 Stage-4 Region-II 样本中固定选取：首个离开 Region I 的高度点、最大高度点、
最大曲率点。分类严格使用 Goal 的既定规则：只有同时满足
`rho_shape>=0.995` 和 `E_aligned<=0.5 E_G` 才归为 global/coherent phase
dominated；否则按规则归为 spatial distortion。

额外诊断从已保存的场构造，不重新运行 PE/Bellhop：对每个 receiver 求解析正弦面
上的 stationary reflection point，核查镜面方向，并计算 phase gradient、
`2k eta` 和 angle-aware `k(cos_i+cos_o)eta`。这些诊断没有进入 production PE
screen，也没有用于校正场。

## 结果

| A (m) | K (rad/m) | E_G | E_aligned | alignment reduction | phi0 (rad) | rho_shape | phase RMS (rad) | TL RMS (dB) | phase-gradient RMS/max (rad/m) | PE/BH vs stationary RMS (rad) | `2k eta` vs angle RMS (rad) | classification | guards |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---|:---:|
| 0.05 | 0.10 | 0.0210357 | 0.0193849 | 7.85% | -0.00818967 | 0.999812 | 0.0208691 | 0.0228496 | 0.002153/0.005593 | 0.0208684 / 1.73e-6 | 0.0202833 | spatial distortion | PASS |
| 0.20 | 0.10 | 0.100485 | 0.0801428 | 20.24% | -0.0609644 | 0.996796 | 0.100142 | 0.0914359 | 0.008103/0.016091 | 0.100140 / 6.18e-6 | 0.0825710 | spatial distortion | PASS |
| 0.02 | 0.47 | 0.0228161 | 0.0210043 | 7.94% | -0.00891677 | 0.999779 | 0.0129059 | 0.163441 | 0.006549/0.026122 | 0.0129018 / 1.29e-5 | 0.00709721 | spatial distortion | PASS |

M99 内 stationary residual 最大 `1.40e-7`，解析镜面方向误差最大 `1.42e-7`，
均通过明确的 `1e-6` validation-only 几何容差。全局相位对齐没有消除一半以上的
复场误差，因此稳定的小 `phi0` 不是主因。Bellhop 与解析 stationary phase
高度一致，而 PE 残差与 PE–Bellhop phase RMS 同量级；这把剩余差异进一步定位为
空间相关的反射模型差异，符合 Kirchhoff phase-screen 与 local-specular
stationary/`Reflect2D` 描述不同，而不是新的 convention 或 receiver bookkeeping
问题。

## 产物

- `results/validation/pe_bellhop_controlled_comparison/stage5_phase_attribution/stage5_validation.mat`
- `results/validation/pe_bellhop_controlled_comparison/stage5_phase_attribution/stage5_attribution.csv`
- `results/validation/pe_bellhop_controlled_comparison/stage5_phase_attribution/stage5_report.md`

```text
Stage 5 = PASS
Stage 6 = COMPLETED / PASS (subsequent fixed-PM case)
```

Stage 6 后续结果见
`reports/pe_bellhop_controlled_comparison_stage6_fixed_pm_report.md`。
