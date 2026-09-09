# Bellhop internal-wall visualization

本报告由已保存的 Bellhop validation-only SHD/diagnostic 数据生成；入射场 SHD 来自官方 matched-halfspace `C*X` run，绘图脚本本身不重新运行 PE 或 Bellhop。

- walls: `r=100+0.05z`, `r=100-2sin(0.04z)`; z support: [-32, 32] m; physical Rx: `(97,0)` m.
- inverse map: `r=200-r'`, `z=-z'`; dense SHD contains only the post-wall reflected branch and is plotted as `TL=-20 log10(abs(p))`.
- both TL figures use smooth interpolated Bellhop-style backgrounds; samples below 80 dB relative to the field maximum (or with no finite beam support) are masked as NaN and shown with the axes background.
- incident-only field uses the same dense receiver depth grid and extends past the wall support for context; the two wall curves are overlays, not reflecting boundaries in that run.
- selected representative rays: 8 per wall; gray dashed paths are native backward branches and colored paths are inverse-mapped proper-rotation branches.

## Ray diagnostics

| wall | alpha (deg) | hit r (m) | hit z (m) | residual (m) | inc u_r | inc u_z | ref u_r |
|---|---:|---:|---:|---:|---:|---:|---:|
| tilted | -14.004 | 98.768348 | -24.6330391 | -6.217e-15 | 0.970279 | -0.24199 | -0.989578 |
| tilted | -9.996 | 99.1264236 | -17.4715276 | -1.110e-16 | 0.98482 | -0.173579 | -0.997223 |
| tilted | -6 | 99.4772261 | -10.4554778 | -6.994e-15 | 0.994522 | -0.104528 | -0.999988 |
| tilted | -2.004 | 99.8253522 | -3.49295575 | -1.887e-15 | 0.999388 | -0.0349693 | -0.997892 |
| tilted | 0 | 100 | 0 | -3.539e-12 | 1 | 0 | -0.995012 |
| tilted | 2.004 | 100.17526 | 3.50519926 | 4.996e-15 | 0.999388 | 0.0349693 | -0.990916 |
| tilted | 6 | 100.528297 | 10.5659498 | 5.551e-16 | 0.994522 | 0.104528 | -0.979135 |
| tilted | 9.996 | 100.88911 | 17.78221 | -2.554e-15 | 0.98482 | 0.173579 | -0.962593 |
| sinusoidal | -13.998 | 101.697954 | -25.3523772 | -5.637e-15 | 0.970304 | -0.241888 | -0.946432 |
| sinusoidal | -10.002 | 101.310742 | -17.8674637 | -7.922e-16 | 0.984802 | -0.173683 | -0.956727 |
| sinusoidal | -6 | 100.822565 | -10.5968786 | 1.514e-16 | 0.994522 | -0.104528 | -0.96884 |
| sinusoidal | -1.998 | 100.278946 | -3.4983133 | 6.354e-15 | 0.999392 | -0.0348646 | -0.98144 |
| sinusoidal | 0 | 100 | 0 | -3.527e-12 | 1 | 0 | -0.987282 |
| sinusoidal | 1.998 | 99.7225924 | 3.47890445 | 3.636e-15 | 0.999392 | 0.0348646 | -0.992416 |
| sinusoidal | 6 | 99.1899468 | 10.4252835 | -2.817e-15 | 0.994522 | 0.104528 | -0.999147 |
| sinusoidal | 10.002 | 98.717118 | 17.4100444 | 4.883e-15 | 0.984802 | 0.173683 | -0.998649 |

![ray trajectories](../results/visualization/bellhop_internal_wall_visuals/internal_wall_ray_trajectories.png)

![reflected-only 2-D TL](../results/visualization/bellhop_internal_wall_visuals/internal_wall_tl_2d.png)

![incident-only 2-D TL](../results/visualization/bellhop_internal_wall_visuals/internal_wall_incident_tl_2d.png)
