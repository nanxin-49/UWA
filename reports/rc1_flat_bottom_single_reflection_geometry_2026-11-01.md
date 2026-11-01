# RC1 平坦海底一次镜面反射：D/S/B 解析几何补充

日期：2026-11-01。此补充独立于既有 [D/S 解析图与孔径设计报告](rc1_target_deepsea_geometry_and_pe_aperture_design_2026-11-01.md)；原报告、原图和原始结果未改。仅运行 MATLAB 镜像法解析脚本 [`render_target_deepsea_dsb_geometry.m`](../scripts/reporting/render_target_deepsea_dsb_geometry.m)，没有运行 PE/WAPE、BIE、FullFK 或声压计算。

## 定义与结果

保持水深 `H=1500 m`、发射深度 `z_t=1490 m`、接收深度 `z_r∈{3,6,10} m`、水平偏移 `x=0:1:200 m`、等声速 `c=1500 m/s`。深度向下为正。D 为直达，S 为平坦压力释放海面一次镜面路径；新增 B 为**平坦海底一次镜面路径的几何线**，将声源镜像到 `z'_t=2H−z_t=1510 m`。路径长度和到达时间为

\[
L_D=\sqrt{x^2+(z_t-z_r)^2},\quad
L_S=\sqrt{x^2+(z_t+z_r)^2},\quad
L_B=\sqrt{x^2+(2H-z_t-z_r)^2},\qquad
\tau_j=L_j/c.
\]

图与 CSV 使用同一接收点的 D 到达为零点，列出 `τ_S−τ_D`、`τ_B−τ_D`，并保存 `τ_B−τ_S`。未对不同接收点重新作额外移峰，也未以轴向距离代替斜路径长度。

| 接收深度 | 偏移 | S 相对 D | B 相对 D | B 相对 S |
|---:|---:|---:|---:|---:|
| 3 m | 0 m | 4.000 ms | 13.333 ms | 9.333 ms |
| 3 m | 200 m | 3.964 ms | 13.216 ms | 9.251 ms |
| 6 m | 0 m | 8.000 ms | 13.333 ms | 5.333 ms |
| 6 m | 200 m | 7.929 ms | 13.215 ms | 5.287 ms |
| 10 m | 0 m | 13.333 ms | 13.333 ms | **0 ms** |
| 10 m | 200 m | 13.215 ms | 13.215 ms | **0 ms** |

当 `z_r=H−z_t=10 m` 时，S 的竖向镜像距离 `z_t+z_r` 与 B 的 `2H−z_t−z_r` 都是 1500 m。因此 `L_S=L_B`、`τ_S=τ_B` 对本模型的**所有水平偏移**成立，而不只是某一个采样点。0–200 m 的 201 个接收偏移的保存值最大差为 0 s。深度 3 和 6 m 时 B 晚于 S。

## 图表、复核和解释边界

- [D/S/B 相对时延图 PNG](../results/visualization/rc1_target_deepsea_dsb_geometry/20261101_01/F01_DSB_relative_delays.png)；[可编辑 FIG](../results/visualization/rc1_target_deepsea_dsb_geometry/20261101_01/F01_DSB_relative_delays.fig)。三个接收深度采用共同 0–14.5 ms 纵轴；10 m 面板用叠合虚线和空心点显式显示 S/B 重合。标题标明等声速、平坦边界和解析几何，并声明非 PE 或实测声压。
- [完整 CSV（603 行）](../results/visualization/rc1_target_deepsea_dsb_geometry/20261101_01/dsb_geometry_samples.csv)、[关键点 CSV（12 行）](../results/visualization/rc1_target_deepsea_dsb_geometry/20261101_01/dsb_geometry_key_points.csv)、[配置](../results/visualization/rc1_target_deepsea_dsb_geometry/20261101_01/config.json)。没有平滑、拟合、插值或不确定度估计。
- MATLAB Code Analyzer 对新增脚本报 0 项；解析运行回执为 `RC1_DSB_ANALYTIC_PASS`，603 行、1 对图件、10 m S/B 最大差 0 s。独立读取 CSV 与旧 D/S CSV 逐行对齐后，D/S 时延最大差均为 0 s；另按 B 镜像公式复算 CSV 最大绝对差约 `4.89×10⁻¹⁵ s`（十进制写出舍入）。FIG 已重开，含 3 个坐标轴和 10 条线对象；PNG 已目视检查。

**适用边界：**B 仅表示假设平坦海底存在单次镜面几何路径时的候选到达时刻。未指定或估计海底声阻抗、反射系数、粗糙度、损耗和相位；也没有计算任何路径的实际复声压。S/B 到达时间重合不代表声压幅相相同，不支持关于海底反射强度、可分辨多径数量或实际信道簇结构的推断。图不构成目标几何的 PE 资格或宽带 CIR 验证。
