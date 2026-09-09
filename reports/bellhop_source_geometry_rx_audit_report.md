# Bellhop point-source R versus line-source X amplitude audit

状态：**PASS_WITH_LIMITS**。本审计仅改变 Bellhop `RunType(4)`：`R`（point source）或 `X`（line source）；Gaussian `.sbp`、step=0.05 m、beam 数、receiver、wall geometry、Reflect2D、p/q 与 InfluenceGeoHatCart 均未改变。

Bellhop 2020 源码静态核查与数值结果一致：`InfluenceGeoHatCart` 对 `R` 使用
`sqrt(abs(cos(alpha)))` 的 source factor，而 `ScalePressure` 对 `R` 额外使用
`1/sqrt(abs(r))`，`X` 则采用与 range 无关的 line-source scale。该源码仅被读取，
没有修改。

## Flat internal-wall regression

| source | beams | TL error (dB) | phase error (rad) | complex relative error |
|---|---:|---:|---:|---:|
| R | 5001 | 0 | 0 | 0 |
| R | 10001 | 0 | 0 | 0 |
| X | 5001 | 0 | 6.88168084e-17 | 6.88168084e-17 |
| X | 10001 | 0 | 6.88168084e-17 | 6.88168084e-17 |

Reference is `H_wall=-P_BH(103 m)` under the same source convention.

## Weak sinusoidal native ATI versus internal wall

A=0.25 m, K=-0.01 1/m, N=161; native receiver r=97 m and rotated-wall receiver r=103 m.

| source | beams | TL(native-wall) (dB) | phase (rad) | complex relative error |
|---|---:|---:|---:|---:|
| R | 5001 | 2.89572521 | 4.535663e-06 | 0.395681303 |
| R | 10001 | -0.260650435 | 2.7392258e-08 | 0.0295627058 |
| X | 5001 | 3.15637973 | 4.47710572e-06 | 0.438199013 |
| X | 10001 | 4.24400274e-06 | -3.91849339e-08 | 4.90177729e-07 |

The point-source range factor predicts `10 log10(97/103)=-0.260654904 dB`; measured R is -0.260650435 dB. X leaves 4.24400274e-06 dB. The absolute-error reduction is 0.260646191 dB.

## 5001 -> 10001 beam convergence

| case | source | TL change (dB) | phase change (rad) | complex change |
|---|---|---:|---:|---:|
| flat | R | 0 | 0 | 0 |
| flat | X | 0 | 0 | 0 |
| sinusoidal | R | 3.15637565 | 4.50827074e-06 | 0.30468561 |
| sinusoidal | X | 3.15637549 | 4.51629066e-06 | 0.304685597 |

The 5001-beam sinusoidal field is not converged; both R and X share the same approximately 3.16 dB 5001-to-10001 change. Conclusions therefore use the requested 10001-beam endpoint, and no parameter was tuned.

## Fixed-PM Tier-1, 4 kHz, seed=260001

| source | G_PE | G_BH | Delta TL (dB) | Delta phase (rad) | complex error |
|---|---:|---:|---:|---:|---:|
| R (authoritative prior) | -0.086869244+0.98806868i | -0.69033909-0.66286395i | 0.310433523 | -2.24820076 | 1.83664201 |
| X (this audit) | -0.086869244+0.98806868i | -0.6903441-0.66286873i | 0.310370739 | -2.24820073 | 1.83663521 |

## Conclusions

- The ~0.26 dB native/internal covariance bias is primarily attributable to the point-source R range normalization.
- Tier-1 should formally use X as the dimensionally matched line-source Bellhop comparator.
- After X, the original 0.310434 dB PE--Bellhop Tier-1 amplitude delta becomes 0.310371 dB.
- The reflection-model discrepancy conclusion does not require revision. X removes a Bellhop coordinate/source normalization artifact but does not remove the rough/flat cross-model residual.

No empirical scaling, `.sbp` refit, phase correction, or physics change was used.
