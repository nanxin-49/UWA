# PE–Full Kirchhoff–BIE 可视化综合分析报告

状态：**validation-only 分析汇总；未修改 PE、BIE 或 Full Kirchhoff 生产物理**。

本报告统一解释以下现有结果：

- deterministic sinusoidal Full Kirchhoff parameter sweep；
- PE–Full Kirchhoff reduced-operator analysis；
- finite/adaptive aperture scans；
- Full-Kirchhoff operator visualization atlas。

所有正式指标均使用固定 Stage-0 M99 incident-energy footprint、共同 `-40 dB` 有效幅度门槛和保存的能量权重。Full Kirchhoff 与 BIE 使用冻结的 native-to-comparison conjugation convention。

## 1. Receiver-line 场分布

代表图：

- `results/validation/pe_bie_full_kirchhoff/figures/weak_low_K_reflected.png`
- `results/validation/pe_bie_full_kirchhoff/figures/strong_height_low_K_reflected.png`
- `results/validation/pe_bie_full_kirchhoff/figures/weak_high_K_reflected.png`

### weak/low-K

PE Model-0、kz-aware Model-1、Full Kirchhoff 和 BIE 在有效中心 receiver 区域基本重合。该图说明在小高度、长空间波长条件下，local phase-screen 已能保留主要场形；剩余差异主要是小的有限角度/归一化误差。

### strong-height/low-K

随着 `A` 增大，Model-0 和 Model-1 的幅度峰谷、局部相位和场形开始偏离 BIE。Full Kirchhoff 仍与 BIE 重合。该现象说明 PE 中的局部 `2*k*eta` 相位代理不能完全复现真实 surface-to-receiver path phase 和表面上的非局部相干积分。

### weak/high-K

PE 场出现明显的准周期幅度起伏，而 Full Kirchhoff/BIE 仍保持一致。这是非局域 coherent redistribution / spectral coupling 缺失的直接图像证据；它不符合简单全局幅值因子或单一 phase offset 的特征。

图两端低幅度区域的 Full-Kirchhoff 快速振荡不作为主要物理结论。正式误差只在 M99 与 `-40 dB` mask 内评价。

## 2. PE 与 Full Kirchhoff 的误差分布

代表图：

- `results/validation/pe_bie_full_kirchhoff_reduced_stage1/figures/weak_low_K_PS_FK_error.png`
- `results/validation/pe_bie_full_kirchhoff_reduced_stage1/figures/strong_height_low_K_PS_FK_error.png`
- `results/validation/pe_bie_full_kirchhoff_reduced_stage1/figures/weak_high_K_PS_FK_error.png`

误差不是白噪声，而是随 receiver 横向坐标变化的平滑、相关结构。典型趋势是：

- weak/low-K：误差较小，主要由 phase difference 构成；
- strong-height/low-K：phase residual 显著增大，magnitude residual 仍较小；
- weak/high-K：出现明显的周期性 magnitude difference，同时保留较小但非零的 phase difference。

因此 PE–Full Kirchhoff 差异是结构性模型差异，而不是随机数值噪声。

## 3. 入射场二维图

图件：

- `results/validation/pe_bie_full_kirchhoff_parameter_sweep/figures/weak_lowK_incident_xz.png`
- `results/validation/pe_bie_full_kirchhoff_parameter_sweep/figures/weak_lowK_incident_surface.png`

入射场从 `z=100 m` 的 Gaussian source 向上扩展，在 surface 附近仍具有有限横向 footprint。不同 surface point `x_s` 的 incident field 并非只对应同一横向 receiver；Full Kirchhoff 会把各 surface point 经 Green propagation 后的贡献重新相干叠加。

这解释了为什么 local phase screen 在强高度或高 K 下失效：它只保留局部边界相位，而没有显式保留完整的 source-point-to-receiver coupling。

## 4. Full Kirchhoff integrand kernel 图

图件位于：

- `results/visualization/pe_bie_full_kirchhoff_operator_atlas/figures/weak_low_K_integrand_kernel.png`
- `results/visualization/pe_bie_full_kirchhoff_operator_atlas/figures/strong_height_low_K_integrand_kernel.png`
- `results/visualization/pe_bie_full_kirchhoff_operator_atlas/figures/weak_high_K_integrand_kernel.png`

这里显示的不是单独 `|G|`，而是完整的离散积分密度：

`I(x_r,x_s) = -[G d_n(p_inc) + (d_n_s G) p_inc] w_s`。

图中的 magnitude 结构显示贡献集中在 receiver 与 surface source point 的几何耦合带附近，但并不严格局限于窄对角线。phase 图呈现高密度振荡条纹，说明远离局部几何邻域的 surface contributions 仍可能通过相干叠加改变最终复场。

因此，Full Kirchhoff 不能安全地被解释为一个非常短的、固定平移不变 convolution kernel。

## 5. 单 receiver 的海面逐点贡献

图件：

- `results/visualization/pe_bie_full_kirchhoff_operator_atlas/figures/weak_low_K_single_receiver_contribution.png`
- `results/visualization/pe_bie_full_kirchhoff_operator_atlas/figures/strong_height_low_K_single_receiver_contribution.png`
- `results/visualization/pe_bie_full_kirchhoff_operator_atlas/figures/weak_high_K_single_receiver_contribution.png`

代表 receiver 为 `x_r≈19.92 m`。图中同时给出：

- `|I(x_r,x_s)|`；
- contribution phase；
- real/imaginary contribution；
- 从 surface support 累积的 coherent sum。

三个 case 的 contribution cancellation ratio 均约为 `22.8`。也就是说，最终 receiver field 远小于所有 surface contributions 的绝对值总和，结果依赖精细的相位抵消，而不是少数局部点的简单叠加。

这也是 aperture 截断需要较大 support 的主要直观原因：被删去的远端项单独看可能不大，但其相位关系会改变总和。

## 6. Geometry、source、receiver 和 aperture

图件：

- `results/visualization/pe_bie_full_kirchhoff_operator_atlas/figures/weak_low_K_geometry_support_aperture.png`
- `results/visualization/pe_bie_full_kirchhoff_operator_atlas/figures/strong_height_low_K_geometry_support_aperture.png`
- `results/visualization/pe_bie_full_kirchhoff_operator_atlas/figures/weak_high_K_geometry_support_aperture.png`

这些图明确标出：

- rough surface support；
- Tx at `z=100 m`；
- receiver line at `z=3 m`；
- 代表 receiver；
- `L=4,16,64,128 m` 的 receiver-centered aperture。

由于 surface height 是厘米到分米量级，而 Tx/Rx 深度是百米量级，几何图中 roughness 本身会被压缩在 `z≈0` 附近。这是尺度差异，不代表 surface geometry 没有进入积分；kernel 图和 contribution 图才是观察其传播影响的主要工具。

## 7. Height sweep 与 K sweep

图件：

- `results/validation/pe_bie_full_kirchhoff_parameter_sweep/figures/height_sweep_complex_error.png`
- `results/validation/pe_bie_full_kirchhoff_parameter_sweep/figures/wavenumber_sweep_complex_error.png`
- 多指标版本：`results/visualization/pe_bie_full_kirchhoff_operator_atlas/figures/height_sweep_all_metrics.png`
- 多指标版本：`results/visualization/pe_bie_full_kirchhoff_operator_atlas/figures/wavenumber_sweep_all_metrics.png`

### Height sweep，固定 `K=0.10 rad/m`

PE–BIE complex L2 从：

- `1.58e-3` at `A=0.005 m`
- 增至 `8.09e-2` at `A=0.20 m`

Full Kirchhoff–BIE 保持约 `2.54e-4` 到 `2.57e-4`。因此 height 增大主要放大 PE 的 phase/path residual，而不是 Full Kirchhoff 的误差。

### K sweep，固定 `A=0.02 m`

PE–BIE complex L2 从约 `9.41e-3` 到 `4.84e-2`，在 `K=0.70 rad/m` 达到最大。Full Kirchhoff–BIE 只从约 `2.54e-4` 增至 `3.32e-4`，仍接近 BIE/integration floor。

高 K 的 PE 误差同时包含明显 magnitude redistribution，因此不能只归结为一个局部相位漏项。

## 8. Reduced Kirchhoff aperture 收敛

图件：

- `results/visualization/pe_bie_full_kirchhoff_operator_atlas/figures/all_cases_aperture_convergence.png`
- `results/validation/pe_bie_full_kirchhoff_adaptive_aperture_stage1/figures/representative_aperture_convergence.png`

典型趋势：

- `L=4 m`：complex L2 约 `6.4e-2`；
- `L=16 m`：约 `2.0e-2`；
- `L=64 m`：约 `1.5e-3`；
- `L=128 m`：约 `2.5e-4` 到 `3.3e-4`，进入 `1e-3` gate；
- `L=Inf`：与 Full Kirchhoff 完全一致。

有限 aperture 明显优于 local phase screen 的某些误差水平，但没有一个小的固定 aperture 在严格 `Ec<1e-3` 门槛下覆盖所有 case。自适应 Stage-1 的 11 个 case 全部首次在 `L=128 m` 通过，说明当前硬截断定义在该精度目标下实际上是近似全局的。

因此目前没有足够证据建立有意义的 `L_min=f(A,K)` 经验规律。

## 9. 综合物理判断

### Full Kirchhoff 与 BIE

在当前确定性正弦海面、4 kHz、one-transverse、均匀介质范围内，Full Kirchhoff 已经数值闭合到 Helmholtz BIE：

- 最大 Full Kirchhoff–BIE complex L2 约 `3.32e-4`；
- 2049→4097 surface integration 的最大 convergence complex L2 约 `3.72e-5`；
- operator atlas 的 integrand 回代相对误差为 `0`。

因此剩余 PE–BIE discrepancy 不能主要归因于 Full Kirchhoff 的实现误差。

### PE 误差来源

现有消融和图件共同支持以下排序：

1. surface height / complete path phase 是必要因素；
2. local phase-screen 对完整 surface-to-receiver path phase 的近似不足；
3. high-K 时还存在明显的 nonlocal coherent redistribution / spectral coupling；
4. curved normal 和 `ds` Jacobian 在当前测试坡度下不是主要误差来源。

kz-aware Model-1 在 weak/low-K 可显著减少 phase error，但在 strong-height 和 high-K 情况下仍不能闭合到 BIE，说明严格法向修正不是全部答案。

### 当前不能得出的结论

这些图还不能证明：

- 某个固定小 `L` 对所有海况都适用；
- 存在简单的平移不变卷积 kernel；
- 一个全局 amplitude correction 可以代替非局部积分；
- 当前确定性正弦结果可以直接外推到 PM random surface 或 3D。

## 10. 下一步建议

当前最合理的研究顺序是：

1. 以 operator atlas 为基础，测试低秩或 locally-stationary kernel compression；
2. 保留完整 geometry-dependent Green/normal structure，避免直接把 `u_FK-u_PE` 当 correction；
3. 先在 weak/low-K 和 weak/high-K 两个端点验证压缩误差，再考虑 strong-height；
4. 只有在压缩后的 reduced operator 通过同一 BIE mask/gate 后，才考虑接入更大范围验证。

## 最终结论

当前图件形成了一个闭合证据链：

`Full Kirchhoff ≈ BIE`，而 `PE local phase screen` 在 height 增大或 K 增大时逐步偏离；偏离结构主要来自完整非局部表面积分及其相干抵消，而不是单纯 amplitude normalization 或 Full Kirchhoff 数值不稳定。

有限 aperture 能改善 PE 误差，但严格精度目标要求的支持范围接近全局。因此下一阶段应优先研究 **low-rank / locally-stationary reduced Kirchhoff operator**，而不是继续拟合固定 aperture 规律。

## 11. Gaussian 源宽度与 angular-spectrum 结果

此前未纳入本汇总的图件：

- [angular spectrum versus sigma](../results/validation/pe_bie_full_kirchhoff_gaussian_sigma_sensitivity/figures/angular_spectrum_sigma.png)
- [source beam magnitude versus sigma](../results/validation/pe_bie_full_kirchhoff_gaussian_sigma_sensitivity/figures/source_beam_magnitude_sigma.png)
- [error versus sigma](../results/validation/pe_bie_full_kirchhoff_gaussian_sigma_sensitivity/figures/error_vs_sigma.png)
- [theta99 versus PE error](../results/validation/pe_bie_full_kirchhoff_gaussian_sigma_sensitivity/figures/theta99_vs_pe_error.png)
- [weak/low-K receiver fields](../results/validation/pe_bie_full_kirchhoff_gaussian_sigma_sensitivity/figures/weak_low_K_receiver_fields_sigma.png)
- [strong-height/low-K receiver fields](../results/validation/pe_bie_full_kirchhoff_gaussian_sigma_sensitivity/figures/strong_height_low_K_receiver_fields_sigma.png)
- [weak/high-K receiver fields](../results/validation/pe_bie_full_kirchhoff_gaussian_sigma_sensitivity/figures/weak_high_K_receiver_fields_sigma.png)

减小 sigma 会收窄角谱：`theta99` 从 `21.28 deg`（sigma=`0.3 m`）降到 `1.23 deg`（sigma=`5 m`）。在 weak/low-K case 中，PE phase RMS 从约 `0.0158 rad` 降到 `0.0040 rad`，说明有限入射角确实是 Model-0 误差的一部分。可是 PE complex error 没有单调降到 Full Kirchhoff/BIE floor；strong-height/low-K 仍保持约 `0.076` 的 complex error，weak/high-K 在窄角端的 magnitude error 反而上升到约 `0.0266`。

因此这些图只支持有限角度贡献的判断，不能把 strong-height 或 high-K residual 归因于角谱宽度。Full Kirchhoff–BIE 在同一 sigma sweep 中仍保持约 `1e-5` 到 `2.7e-4`。

## 12. Sign / convention audit 可视化

- [weak/low-K before-after](../results/validation/pe_bie_full_kirchhoff_convention_audit/figures/weak_low_K_before_after.png)
- [strong-height/low-K before-after](../results/validation/pe_bie_full_kirchhoff_convention_audit/figures/strong_height_low_K_before_after.png)
- [weak/high-K before-after](../results/validation/pe_bie_full_kirchhoff_convention_audit/figures/weak_high_K_before_after.png)
- [weak/low-K delta-phi](../results/validation/pe_bie_full_kirchhoff_convention_audit/figures/weak_low_K_delta_phi_vs_4keta.png)
- [strong-height delta-phi](../results/validation/pe_bie_full_kirchhoff_convention_audit/figures/strong_height_low_K_delta_phi_vs_4keta.png)
- [weak/high-K delta-phi](../results/validation/pe_bie_full_kirchhoff_convention_audit/figures/weak_high_K_delta_phi_vs_4keta.png)

审计图显示旧的 native rough/flat ratio 与已经 conjugated 的 comparison field 混合比较时，会出现接近 `-4*k*eta` 的系统相位结构。统一 `z_s=+eta`、water normal、`exp(-i omega t)` outgoing Green function，并只应用一次 conjugation 后，phase residual 回到约 `1e-4 rad` 到 `2e-4 rad` 的 integration floor。

因此 surface-z、normal、Green derivative 和 field-definition mismatch 已被排除为当前 PE–BIE residual 的主来源。`z_s=-eta` 的 Case B 是镜像物理曲面，不是同一 authoritative case。

## 13. 物理项消融和 PE–Full Kirchhoff 趋势

- [weak/low-K ablations](../results/validation/pe_bie_full_kirchhoff_reduced_stage1/figures/weak_low_K_ablations.png)
- [strong-height/low-K ablations](../results/validation/pe_bie_full_kirchhoff_reduced_stage1/figures/strong_height_low_K_ablations.png)
- [weak/high-K ablations](../results/validation/pe_bie_full_kirchhoff_reduced_stage1/figures/weak_high_K_ablations.png)
- [PE–Full Kirchhoff parameter trends](../results/validation/pe_bie_full_kirchhoff_reduced_stage1/figures/PS_FK_parameter_trends.png)

将实际 surface height/path phase 替换为 flat surface `z=0` 会造成约 `0.23`、`1.22` 和 `0.46` 的 complex error；把 curved normal 换成 flat normal，或把 `ds` 换成 `dx`，误差通常只有 `1e-6` 到 `1e-4`。这说明 surface height/path phase 是必要物理项，而 normal 和 Jacobian 在当前坡度范围内不是主要缺失项。

PE–Full Kirchhoff 趋势图进一步显示：weak/low-K 的误差主要是 phase；strong-height/low-K 的 phase residual 随 height 增长；weak/high-K 同时出现相干幅度重分配。kz-aware Model-1 能明显削弱弱低 K 的有限角度 phase 误差，但不能处理完整的 nonlocal surface coupling。

## 14. Aperture reduction 的补充解释

- [finite-aperture reduction factor](../results/validation/pe_bie_full_kirchhoff_finite_aperture/figures/reduced_aperture_reduction_factor.png)
- [finite-aperture convergence](../results/validation/pe_bie_full_kirchhoff_finite_aperture/figures/reduced_aperture_convergence.png)
- [strong-height reduced versus Full error](../results/validation/pe_bie_full_kirchhoff_finite_aperture/figures/strong_height_low_K_reduced_full_error.png)
- [weak/high-K reduced versus Full error](../results/validation/pe_bie_full_kirchhoff_finite_aperture/figures/weak_high_K_reduced_full_error.png)
- [sparse A–K L_min map](../results/validation/pe_bie_full_kirchhoff_adaptive_aperture_stage1/figures/Lmin_sparse_AK_map.png)
- [L_min relationships](../results/validation/pe_bie_full_kirchhoff_adaptive_aperture_stage1/figures/Lmin_relationships.png)

当 `E_PE,BIE/E_RK,BIE>1` 时，有限 aperture 已经比 PE phase screen 更接近 BIE；但严格 `Ec<1e-3` gate 下，11 个 Stage-1 case 均首次在 `L=128 m` 通过。A–K map 和关系图呈现的是统一的 right-censoring，而不是已确定的 `L_min(A,K)` 规律。

这组图支持“有限支持有用”，但不支持“固定小 aperture 已足够”或“简单 aperture law 已确定”。

## 15. 当前可视化证据链

1. angular-spectrum 图说明窄角只消除部分低 K phase residual；
2. convention audit 图排除了 Full Kirchhoff 的比较符号错误；
3. ablation 图表明 surface height/path phase 是必要项，normal/Jacobian 是次要项；
4. kernel/contribution 图显示 receiver field 依赖广泛的相干抵消；
5. aperture 图说明达到严格误差门槛需要接近全局 support；
6. height/K sweep 图说明 PE residual 随海面高度和空间 wavenumber 增大，而 Full Kirchhoff 仍接近 BIE。

综合判断仍是：剩余差异主要来自 PE local phase-screen 对完整非局部 surface integral 的降阶误差。下一步应优先研究 geometry-dependent 的 low-rank 或 locally-stationary reduced operator。
