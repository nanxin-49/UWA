# PE--Bellhop controlled comparison — Stage 6 fixed PM

日期：2026-09-13  
状态：**PASS**  
validity region：**III**

## Frozen case

- seed `260001`，4 kHz，uniform `c=1500 m/s`，X source，coherent C；
- Bellhop 10,001 beams，profile `N=4097`；
- canonical coefficient SHA-256：
  `1f7eda465e4ae85b8ac038310edf053f2d062563a3b943bfd83bd59d07015f67`；
- 复用 Stage-0 flat denominator、receiver line、M95/M99 masks 和已关闭的
  comparison convention；没有重生成、重归一化、平滑或拟合表面。

| E_G | TL RMS (dB) | phase RMS (rad) | phi0 (rad) | rho_raw | rho_shape | E_aligned |
|---:|---:|---:|---:|---:|---:|---:|
| 0.958268507 | 1.27399043 | 1.10771425 | -0.798545974 | 0.538113006 | 0.771212836 | 0.678503463 |

PM 统计：RMS/max height `0.186410/0.499456 m`，RMS/max slope
`0.053977/0.116674`，RMS/max curvature `0.0209072/0.0445944 1/m`，最小
曲率半径 `22.4243 m`，`2k sigma_eta=6.24667`。

Bellhop wall residual `7.27e-15 m`，minimum incidence `0.831901`，grazing
fraction `0`，hit maximum curvature `0.0425258 1/m`，minimum post-wall range
increment `0.0396983 m`，pressure-release phase error `7.11e-15 rad`，`p/q`
rotation error `0/0`。全部 numerical/mapping/geometry/finite guards PASS。

轴上 `Delta TL=0.310370739 dB`、`Delta phase=-2.24820073 rad`，精确复现旧 X
audit 的 sanity 值；它不作为整条 receiver-line 的替代。完整 M99 场落入 Region
III，说明该 PM realization 已超出当前光滑正弦采样得到的 Region-I/II 范围；这不是
Bellhop internal-wall 数值失败。

产物位于
`results/validation/pe_bellhop_controlled_comparison/stage6_fixed_pm/`。

