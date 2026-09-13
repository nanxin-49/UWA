# PE--Bellhop controlled comparison — Stage 4 sampled validity map

日期：2026-09-12  
状态：**PASS**  
solver calls：`0`

Stage 4 只读取通过 numerical/mapping guards 的 Stage 1–3 结果，去重共享的
`(A,K)=(0.02,0.10)` case，并应用 Goal 中冻结的 Region I/II/III 判据。没有
重新运行 PE 或 Bellhop，也没有拟合、插值或修改任何场。

| surface | A (m) | K (rad/m) | E_G | E_aligned | TL RMS (dB) | phase RMS (rad) | phi0 (rad) | rho_raw | rho_shape | region | guards |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---:|:---:|
| sinusoid | 0.01 | 0.10 | 0.00409346 | 0.00389577 | 0.00456998 | 0.00405930 | -0.00125773 | 0.999992 | 0.999992 | I | PASS |
| sinusoid | 0.02 | 0.10 | 0.00823397 | 0.00777803 | 0.00914010 | 0.00816578 | -0.00270573 | 0.999966 | 0.999970 | I | PASS |
| sinusoid | 0.05 | 0.10 | 0.0210357 | 0.0193849 | 0.0228496 | 0.0208691 | -0.00818967 | 0.999779 | 0.999812 | II | PASS |
| sinusoid | 0.10 | 0.10 | 0.0441790 | 0.0388528 | 0.0457022 | 0.0438770 | -0.0211066 | 0.999024 | 0.999246 | II | PASS |
| sinusoid | 0.20 | 0.10 | 0.100485 | 0.0801428 | 0.0914359 | 0.100142 | -0.0609644 | 0.994944 | 0.996796 | II | PASS |
| sinusoid | 0.02 | 0.20 | 0.00836932 | 0.00807911 | 0.0308524 | 0.00757942 | -0.00218459 | 0.999965 | 0.999967 | II | PASS |
| sinusoid | 0.02 | 0.35 | 0.0141398 | 0.0133263 | 0.0915312 | 0.00942753 | -0.00472523 | 0.999900 | 0.999911 | II | PASS |
| sinusoid | 0.02 | 0.47 | 0.0228161 | 0.0210043 | 0.163441 | 0.0129059 | -0.00891677 | 0.999740 | 0.999779 | II | PASS |
| fixed PM | — | — | 0.958269 | 0.678503 | 1.27399 | 1.10771 | -0.798546 | 0.538113 | 0.771213 | III | PASS |

## Sampled validity brackets

- 固定 `K=0.10 rad/m`：Region I 的已采样高度为 `A=0.01,0.02 m`；非 I
  已采样高度为 `A=0.05,0.10,0.20 m`，因此转换括区为 `A=(0.02,0.05] m`。
- 固定 `A=0.02 m`：Region I 的已采样波数为 `K=0.10 rad/m`；非 I 已采样
  波数为 `K=0.20,0.35,0.47 rad/m`，因此转换括区为 `K=(0.10,0.20] rad/m`。
- 所有非 I 的正弦样本均为 Region II；正弦 sweep 本身没有 Region III 样本。

Stage 6 完成后，canonical fixed-PM case 已按同一门限加入 map 并落入 Region III；
这不改变上述正弦采样括区。

上述区间只描述实际采样点，不代表精确 validity boundary。

## 产物

- `results/validation/pe_bellhop_controlled_comparison/stage4_validity_map/stage4_validation.mat`
- `results/validation/pe_bellhop_controlled_comparison/stage4_validity_map/stage4_validity_map.csv`
- `results/validation/pe_bellhop_controlled_comparison/stage4_validity_map/stage4_report.md`

```text
Stage 4 = PASS
Stage 5 = COMPLETED / PASS (subsequent attribution audit)
```

Stage 5 后续已完成，见
`reports/pe_bellhop_controlled_comparison_stage5_phase_attribution_report.md`。
