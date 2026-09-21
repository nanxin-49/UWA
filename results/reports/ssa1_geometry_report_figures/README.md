# ssa1_geometry 汇报图说明

本文件夹保存本周 `ssa1_geometry` 相关汇报图。核心结论是：当前海面分支可表述为
**PM-spectrum-driven first-order pressure-release Dirichlet SSA statistical reflection/scattering model**。
它用 PM 高度谱直接生成统计散射功率，不生成具体海面高度 realization；`kirchhoff_spatial`
保留为具体海面 realization 相位屏对照。

## 核心公式

PM 高度谱归一化：

\[
\sigma_\eta=\frac{H_s}{4},
\qquad
\sum W_\eta(K_x,K_y)\Delta K_x\Delta K_y=\sigma_\eta^2.
\]

相干镜面反射：

\[
R_{\rm coh}
=
R_0
\exp\left[-\frac12(\gamma_i+\gamma_s)^2\sigma_\eta^2\right].
\]

法向镜面方向：

\[
\gamma_i=\gamma_s=k_0,\qquad
R_{\rm coh}=-\exp[-2k_0^2\sigma_\eta^2].
\]

一阶 pressure-release / Dirichlet SSA 几何因子：

\[
\gamma(K,f)=\sqrt{\max(k_0^2-|K|^2,0)},
\]

\[
G_{\rm SSA1}(K_s,K_i;f)=4\gamma(K_s,f)\gamma(K_i,f).
\]

非相干统计散射功率：

\[
P_{\rm sca}(K_s)
=
4C_{\rm norm}\gamma(K_s)
\left[
W_\eta*
\left(\gamma(K_i)|\Psi_{\rm inc}(K_i)|^2\right)
\right](K_s)
\Delta K_x\Delta K_y.
\]

能量约束：

\[
E_{\rm coh}+E_{\rm sca}\le E_{\rm inc}.
\]

其中 \(C_{\rm norm}\) 对应 `surface_ssa_scatter_scale`，当前是工程归一化参数，不是实验标定的绝对散射截面常数。

## 图像说明

### 1. `validate_ssa_abs_R_coh_vs_Hs.png`

展示不同频率下 \(|R_{\rm coh}|\) 随 \(H_s\) 的变化。

汇报重点：
- \(H_s=0\) 时 \(|R_{\rm coh}|=1\)，退化为平整自由海面。
- \(H_s\) 增大时相干镜面反射下降。
- 高频下降更明显，因为指数项含 \(k_0^2\)。

### 2. `validate_ssa_abs_R_coh_vs_frequency.png`

展示不同 \(H_s\) 下 \(|R_{\rm coh}|\) 随频率变化。

汇报重点：
- 固定非零 \(H_s\) 时，频率越高 coherent loss 越强。
- 该趋势直接对应
  \[
  R_{\rm coh}=-\exp[-2k_0^2\sigma_\eta^2].
  \]

### 3. `validate_ssa_coherent_energy_fraction_vs_Hs.png`

展示

\[
E_{\rm coh}/E_{\rm inc}
\]

随 \(H_s\) 的变化。

汇报重点：
- 粗糙度增强后，相干镜面能量占比下降。
- 这说明粗糙海面会削弱镜面相干反射。

### 4. `validate_ssa_scatter_budget_fraction_vs_Hs.png`

展示

\[
E_{\rm sca}^{\max}/E_{\rm inc}
\]

随 \(H_s\) 的变化。

汇报重点：
- \(H_s=0\) 时散射预算为 0。
- \(H_s\) 增大后，允许的统计散射预算上升。
- 它和 coherent energy 下降形成互补关系。

### 5. `validate_ssa_W_eta_variance_check.png`

展示目标方差 \(\sigma_\eta^2\) 与离散谱积分

\[
\sum W_\eta\Delta K_x\Delta K_y
\]

的对比。

汇报重点：
- 点落在 \(y=x\) 附近，说明 PM 谱按 `Hs_target` 缩放正确。
- 当前验证中最大相对误差约为 \(8.67\times10^{-15}\)。

### 6. `compare_ssa1_kirchhoff_abs_h_reflect_vs_Hs.png`

比较 `ssa1_geometry` 与 `kirchhoff_spatial` 的参考频点反射通道幅度统计。

汇报重点：
- 该图用于看反射通道幅度随海况变化的统计趋势。
- 不要求 SSA1 与 Kirchhoff realization 逐 seed 相等。

### 7. `compare_ssa1_kirchhoff_abs_R_coh_vs_Hs.png`

展示 SSA1 metadata 中的 \(|R_{\rm coh}|\) 随 \(H_s\) 下降趋势。

汇报重点：
- `kirchhoff_spatial` 没有解析 \(R_{\rm coh}\) metadata，因此相关曲线以 SSA1 为主。
- 该图用于支撑相干反射公式回归。

### 8. `compare_ssa1_kirchhoff_coherent_fraction_vs_Hs.png`

展示 SSA1 的 \(E_{\rm coh}/E_{\rm inc}\) 随 \(H_s\) 下降。

汇报重点：
- 这是相干镜面能量损失的能量审计版本。
- 与 \(|R_{\rm coh}|\) 图给出一致趋势。

### 9. `compare_ssa1_kirchhoff_scatter_budget_vs_Hs.png`

展示 SSA1 的 \(E_{\rm sca}^{\max}/E_{\rm inc}\) 随 \(H_s\) 增强。

汇报重点：
- 粗糙度越强，相干项下降后留给统计散射的能量预算越大。
- 能量限制保证不会出现非物理总反射放大。

### 10. `compare_ssa1_kirchhoff_rms_delta_k_vs_Hs.png`

展示 `ssa1_geometry` 与 `kirchhoff_spatial` 的反射角谱 RMS 横向波数展宽随 \(H_s\) 的变化。

汇报重点：
- Kirchhoff 是具体海面 realization 相位屏模型，反射谱展宽来自随机相位屏。
- SSA1 是统计散射模型，展宽来自 PM 谱驱动的散射功率。
- 当前 16-seed 对照中，Kirchhoff RMS delta-k 从 \(H_s=0\) 到 \(H_s=1.0\) 均明显增大：
  - 4000 Hz: 增量约 \(3.4532\ {\rm rad/m}\)
  - 6000 Hz: 增量约 \(3.2032\ {\rm rad/m}\)
  - 8000 Hz: 增量约 \(3.1521\ {\rm rad/m}\)
  - 10000 Hz: 增量约 \(3.1419\ {\rm rad/m}\)

### 11. `compare_ssa1_kirchhoff_phase_variance_vs_Hs.png`

展示反射通道相位 circular variance 随 \(H_s\) 的变化。

汇报重点：
- 该图用于观察粗糙海面引入的反射相位随机性。
- 有限 seed 下不要求每个相邻海况严格单调，只看统计波动和整体趋势是否与粗糙度增强相容。

## 数据补充

### SSA1 物理趋势验证

数据文件：`../validate_ssa_physical_trends_vertical_result.mat`

覆盖：
- 模型：`pm_convolution`、`ssa1_geometry`
- \(H_s=[0,0.05,0.2,0.5,1.0]\)
- \(f=[4000,6000,8000,10000]\) Hz
- 共 40 个 reduced scalar case
- 74 个检查项，失败数为 0

关键数值：
- coherent reflection 公式回归最大误差：0
- \(W_\eta\) 离散归一化最大相对误差：约 \(8.67\times10^{-15}\)
- 最大能量守恒误差：0
- \(|R_{\rm coh}|\) 范围：约 \(5.60\times10^{-96}\) 到 1

### SSA1 与 Kirchhoff 16-seed 统计对照

数据文件：`../compare_ssa1_kirchhoff_surface_statistics_vertical_result.mat`

完整 reduced-grid 对照覆盖：
- 模型：`kirchhoff_spatial`、`ssa1_geometry`
- \(H_s=[0,0.05,0.2,0.5,1.0]\)
- \(f=[4000,6000,8000,10000]\) Hz
- seed：`12345+(0:15)`，共 16 个
- 网格：`nx=ny=64`
- 共 \(2\times5\times4\times16=640\) 个 reduced scalar case
- `run_table` 640 行，`summary_table` 40 行，`trend_table` 156 行

关键数值：
- hard checks 全部通过：`hard_fail_count=0`
- compatibility checks 全部观察到：`compat_fail_count=0`
- `compatibility_observed_fraction=1.000`
- 最大通道不变量误差：
  \[
  \max |H_f-H_{\rm direct}-H_{\rm reflect}|
  =
  1.55\times10^{-17}.
  \]
- \(H_s=0\) 时 SSA1 与 Kirchhoff 的 `h_reflect/h_total` 差异为舍入误差级，最大约 \(2.22\times10^{-16}\)。
- SSA1 最大能量守恒误差约 \(9.93\times10^{-16}\)。

## 汇报时建议强调

1. `ssa1_geometry` 已经加入一阶 pressure-release / Dirichlet SSA 几何因子：
   \[
   G_{\rm SSA1}=4\gamma_s\gamma_i.
   \]
2. 相干反射项已修正为文献一致的 normal 退化形式：
   \[
   R_{\rm coh}=-\exp[-2k_0^2\sigma_\eta^2].
   \]
3. \(H_s=0\) 时模型退化正确。
4. PM 谱离散归一化正确。
5. FFT 周期卷积与 dense 显式求和一致。
6. 与 Kirchhoff realization 对照时，16-seed 统计趋势相容，但不声称逐 realization 等价。

## 当前限制

- 当前不是完整 SSA/NLSSA。
- 没有 second-order SSA coherent correction。
- 没有 T-matrix 求解器。
- 没有阻抗边界或 Neumann 边界。
- 没有多重散射。
- 没有倏逝散射注入。
- 没有实验标定的 bistatic scattering cross section。
- `surface_ssa_scatter_scale` / \(C_{\rm norm}\) 仍是工程归一化参数。
