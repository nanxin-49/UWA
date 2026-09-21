# PE--Bellhop controlled comparison — Stage 2 height-phase sweep

日期：2026-09-11  
状态：**PASS**

## 范围

Stage 1Y 已达到 `CONVENTION_THEORETICALLY_CLOSED`，因此本阶段按 Goal 解锁，
固定 `K=0.10 rad/m`、4 kHz、uniform `c=1500 m/s`、X source、coherent C、
10,001 beams、`N=4097`，按 `A=[0.01,0.02,0.05,0.10,0.20] m` 顺序执行。
主入口为 `scripts/validation/validate_pe_bellhop_controlled_comparison.m`
的 `stage2` 分支；其结果目录为
`results/validation/pe_bellhop_controlled_comparison/stage2_height_sweep/`。

本轮没有修改 PE/Bellhop 核心、Reflect2D、InfluenceGeoHatCart、SHD selector
或 source normalization。

## 已落盘结果

`A=0.01 m` 复用了已通过的 Stage-1 convention-fixed MAT；A=0.02/0.05/0.10/0.20
各运行一次固定配置的 Bellhop/PE case。所有 geometry、finite 和 metric guards
均通过：

| A (m) | `2kA` | max slope | max/RMS curvature (1/m) | E_G (M99) | phase RMS (rad) | TL RMS (dB) | `phi0` (rad) | `rho_raw` | `rho_shape` | E_aligned | geometry | finite |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---:|:---:|
| 0.01 | 0.335103216 | 0.001 | 0.0001 / 0.0000700716 | 0.00409345672 | 0.00405930291 | 0.00456997571 | -0.00125773333 | 0.999991624 | 0.999992415 | 0.00389576987 | PASS | PASS |
| 0.02 | 0.670206433 | 0.002 | 0.0002 / 0.000140143 | 0.00823396805 | 0.00816578444 | 0.00914010217 | -0.00270573185 | 0.999966108 | 0.999969768 | 0.00777803380 | PASS | PASS |
| 0.05 | 1.67551608 | 0.005 | 0.0005 / 0.000350355 | 0.0210357376 | 0.0208690864 | 0.0228496268 | -0.00818966942 | 0.999778741 | 0.999812270 | 0.0193849209 | PASS | PASS |
| 0.10 | 3.35103216 | 0.010 | 0.0010 / 0.000700690 | 0.0441790479 | 0.0438770230 | 0.0457022488 | -0.0211065624 | 0.999023639 | 0.999246207 | 0.0388527853 | PASS | PASS |
| 0.20 | 6.70206433 | 0.020 | 0.0020 / 0.00140122 | 0.100485227 | 0.100141661 | 0.0914359332 | -0.0609643815 | 0.994943993 | 0.996795793 | 0.0801428322 | PASS | PASS |

几何/数值 guard 明细（全部使用 10,001 beams）：

| A (m) | wall residual (m) | phase jump error (rad) | p/q rotation error | min post dr (m) | center tau error (s) |
|---:|---:|---:|---:|---:|---:|
| 0.01 | `7.10e-15` | `7.11e-15` | `0 / 0` | `0.043326` | `4.68e-15` |
| 0.02 | `7.10e-15` | `7.11e-15` | `0 / 0` | `0.043350` | `4.70e-15` |
| 0.05 | `3.54e-12` | `7.11e-15` | `0 / 0` | `0.043423` | `4.16e-17` |
| 0.10 | `3.54e-12` | `7.11e-15` | `0 / 0` | `0.043545` | `2.78e-17` |
| 0.20 | `3.52e-12` | `7.11e-15` | `0 / 0` | `0.043789` | `1.39e-17` |

逐 case 文件均位于
`results/validation/pe_bellhop_controlled_comparison/stage2_height_sweep/height_A*_K0p1_N4097_B10001/`；
汇总文件为 `stage2_validation.mat`、`stage2_metrics.csv` 和 `stage2_report.md`。

## 执行与 resume 说明

首次长任务曾超过 MCP 单次等待窗口；主驱动的 case-granular resume 机制避免了
重复 A=0.01/0.02/0.05 的已完成计算，并允许继续完成 A=0.10/0.20。若某个
`height_A*_K*_N4097_B10001/stage2_case.mat` 已完整存在且 A/K 匹配，后续调用
只读取该 case，不重复启动求解器；这保证恢复执行不会覆盖或重复已完成结果。

## 当前判定

```text
Stage 2 = PASS
Stage 3 = COMPLETED / PASS (见独立 Stage-3 报告)
```

模型差异随 A 增大而增大，但按 Goal 不作为 Stage2 停止条件；本报告不把它解释为
solver 或 internal-wall 失败。Stage3 后续已按 Goal 单独完成，见
`reports/pe_bellhop_controlled_comparison_stage3_slope_curvature_report.md`。
