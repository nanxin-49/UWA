# Bellhop internal-wall visualization

本报告由已保存的 Bellhop validation-only SHD/diagnostic 数据生成；入射场 SHD 来自官方 matched-halfspace coherent run，绘图脚本本身不重新运行 PE 或 Bellhop。

## 权威状态

本文件是当前唯一权威的 internal-wall 可视化报告，整合并取代早期小振幅墙、独立大振幅墙面预览以及 two-period full-fan 阶段报告。当前图件统一采用 `A=2.4 m`、`K=0.20944 1/m`、完整 `[-30,30] deg` 发射扇区和已标出的 Rx eigenray；旧报告与重复图件仅保存在可恢复归档中。

- walls: `r=100+0.05z`, `r=100-2.4sin(0.20944z)`; sampled wall support: [-75, 75] m; trajectory z window: [-65, 65] m; physical Rx: `(97,0)` m.
- inverse map: `r=200-r'`, `z=-z'`; dense SHD contains only the post-wall reflected branch. Both TL maps use the same axial source reference `p_ref=|p_inc(r=1 m,z=0)|=0.671841281` and are plotted as `TL=-20 log10(abs(p)/p_ref)`.
- This common reference only harmonizes the visualization; it does not renormalize the reflected branch or alter Bellhop physics. The reflected-only field remains a separate post-wall branch and must not be interpreted as a direct-plus-reflected total field.
- both TL figures use smooth interpolated Bellhop-style backgrounds; samples below 80 dB relative to the field maximum (or with no finite beam support) are masked as NaN and shown with the axes background.
- incident-only field uses the same dense receiver depth grid and extends past the wall support for context; the two wall curves are overlays, not reflecting boundaries in that run.
- TL axes include the complete sampled wall support. Regions outside the computed SHD receiver grid remain background/NaN; no field values are extrapolated to fill them.
- selected representative rays: tilted 9, sinusoidal 10; the latter includes the Rx eigenray near 1.23856 deg. Gray dashed paths are native backward branches and colored paths are inverse-mapped proper-rotation branches.

## 当前图件介绍与分析

两张二维 TL 图的背景都是 Bellhop 密集接收网格上的相干复声压；TL 仅表示压力幅度，不是单条声线能量或总场功率。入射图只包含 Tx 到墙前的场，反射图只包含一次 wall reflection 后、经 proper rotation 并逆映射回物理坐标的 reflected-only 分支。

### 统一参考与绝对幅度

两张图现在共同使用 `p_ref=|p_inc(r=1 m,z=0)|=0.671841281`，色标统一为 `[-20, 100] dB`。这一统一只改变显示参考，不对任何分支做幅度拟合、重标定或物理修正。

| field | max `|p|` (Bellhop units) | common-reference TL range (dB) | median TL (dB) | valid fraction |
|---|---:|---:|---:|---:|
| incident-only | 20 | [-29.475, 50.524] | 38.094 | 0.334 |
| tilted reflected-only | 0.00999501 | [36.550, 47.909] | 40.061 | 0.987 |
| sinusoidal reflected-only | 0.068677 | [19.809, 89.716] | 46.298 | 0.966 |

统一参考后，入射、tilted reflected-only 和 sinusoidal reflected-only 的最大原始幅度分别为 `20`、`0.00999501` 和 `0.068677`；两个反射图最大值相对入射图最大值分别为 `-66.025 dB` 和 `-49.284 dB`。这些是不同空间区域内的场最大值，只用于说明共同色标下的显示动态范围，不能直接解释为墙面的能量反射率。

### 图形形状、灰色区域与平滑性

入射图使用有限 Gaussian 发射扇区，灰色三角区域对应没有有效 beam support 的位置（当前入射网格有效比例约为 `33.41%`）。反射图的物理 range 来自 `r=200-r'`，因此正 range mapped branch 映射回物理坐标后，反射场自然显示在墙的左侧。tilted case 的有效比例约为 `98.66%`；sinusoidal case（`r=100-2.4sin(0.20944z)`）约为 `96.57%`。

背景的平滑性来自相干 Bellhop beam influence 在密集二维接收网格上的累积以及绘图时的面内插值；当前没有海底、多次反射或直达/反射总场干涉，因此不应期待标准多途 Bellhop 图中常见的强烈干涉条纹。图上的代表性声线（tilted 9 条、sinusoidal 10 条）只是几何叠加，不是 TL 背景中的全部 beam 或能量脊线。掩膜边界的少量台阶来自接收网格和 NaN support 边界，而不是传播算子台阶。

### 声线几何解释

声线不需要垂直撞击墙面；正确条件是入射、反射方向关于局部墙面法线镜像对称。保存的全射线几何审计给出 tilted/sinusoidal 最大镜面方向误差分别为 `6.697e-16` 和 `7.022e-16`，对应 `max|t dot n|` 分别为 `0.000e+00` 和 `0.000e+00`。因此轨迹图中看似斜入射的线是正确的 specular reflection，而不是法线或比例错误。

## Ray diagnostics

| wall | alpha (deg) | hit r (m) | hit z (m) | residual (m) | inc u_r | inc u_z | ref u_r |
|---|---:|---:|---:|---:|---:|---:|---:|
| tilted | -27.996 | 97.4107261 | -51.7854789 | 7.105e-15 | 0.88298 | -0.46941 | -0.9254 |
| tilted | -20.004 | 98.2122942 | -35.7541167 | -3.775e-15 | 0.939669 | -0.342086 | -0.969105 |
| tilted | -12 | 98.9483935 | -21.0321303 | 6.661e-16 | 0.978148 | -0.207912 | -0.994008 |
| tilted | -3.996 | 99.6519325 | -6.96135094 | -4.607e-15 | 0.997569 | -0.0696868 | -0.999545 |
| tilted | 0 | 100 | 0 | -3.539e-12 | 1 | 0 | -0.995012 |
| tilted | 3.996 | 100.350508 | 7.01015106 | 1.388e-15 | 0.997569 | 0.0696868 | -0.985642 |
| tilted | 12 | 101.074199 | 21.4839843 | 4.441e-16 | 0.978148 | 0.207912 | -0.95253 |
| tilted | 20.004 | 101.853994 | 37.0798749 | 3.109e-15 | 0.939669 | 0.342086 | -0.900859 |
| tilted | 27.996 | 102.730684 | 54.6136742 | -2.220e-15 | 0.88298 | 0.46941 | -0.831753 |
| sinusoidal | -27.996 | 97.6191681 | -51.8962909 | -5.825e-15 | 0.88298 | -0.46941 | -0.935183 |
| sinusoidal | -20.004 | 102.397367 | -37.2776894 | -2.529e-15 | 0.939669 | -0.342086 | -0.922644 |
| sinusoidal | -12 | 97.7543728 | -20.7783334 | 3.707e-15 | 0.978148 | -0.207912 | -0.989993 |
| sinusoidal | -3.996 | 102.393653 | -7.15287838 | -5.444e-15 | 0.997569 | -0.0696868 | -0.989831 |
| sinusoidal | 0 | 100 | 0 | -3.162e-12 | 1 | 0 | -0.596612 |
| sinusoidal | 1.236 | 98.9621736 | 2.1351699 | 1.347e-15 | 0.999767 | 0.0215706 | -0.675265 |
| sinusoidal | 3.996 | 97.6243504 | 6.81971083 | -3.369e-15 | 0.997569 | 0.0696868 | -0.997354 |
| sinusoidal | 12 | 102.371187 | 21.7596675 | -3.334e-15 | 0.978148 | 0.207912 | -0.934346 |
| sinusoidal | 20.004 | 97.7876239 | 35.5995158 | 4.253e-15 | 0.939669 | 0.342086 | -0.999363 |
| sinusoidal | 27.996 | 102.223232 | 54.3439031 | 6.588e-15 | 0.88298 | 0.46941 | -0.993463 |

![ray trajectories](../results/visualization/bellhop_internal_wall_visuals/internal_wall_ray_trajectories.png)

![reflected-only 2-D TL](../results/visualization/bellhop_internal_wall_visuals/internal_wall_tl_2d.png)

![incident-only 2-D TL](../results/visualization/bellhop_internal_wall_visuals/internal_wall_incident_tl_2d.png)
