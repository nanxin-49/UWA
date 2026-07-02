# SSA1 组会汇报材料（2026-06-24）

本目录由 `generate_ssa1_group_meeting_report_vertical.m` 生成。核心模型未修改，新增内容仅用于汇报统计、可视化和适用性审计。

## 建议汇报顺序

1. PM 谱归一化：`03_PM_variance_normalization.png`。
2. 相干反射随海况和频率变化：`01`、`02`。
3. 相干能量与散射预算：`04`。
4. 复数集合平均得到的相干通道代理：`05`。
5. SSA1 与 Kirchhoff 的角谱展宽和相位统计：`06`、`07`。
6. 16-seed 稳定性：`08`。
7. 代表性二维角谱与径向谱：`09`、`10`。
8. 趋势相容性和适用性边界：`11`、`12`。

## 图表说明

- `01_Rcoh_vs_Hs.png`：验证粗糙度增大时解析相干反射下降。
- `02_Rcoh_vs_frequency.png`：验证固定海况下高频 coherent loss 更强。
- `03_PM_variance_normalization.png`：验证离散 PM 谱积分等于目标海面高度方差。
- `04_coherent_and_scatter_budget.png`：展示相干能量下降与可用散射预算上升的互补关系。
- `05_coherent_channel_proxy_compare.png`：使用复数集合平均比较 SSA1、Kirchhoff 和解析相干项。强粗糙度下 SSA1 的非零平台主要是 16-seed 随机散射残差，不应解释为解析相干项。
- `06_rms_k_broadening_increment.png`：扣除平整海面基线后的角谱展宽。Kirchhoff 展宽明显强于当前 SSA1。
- `07_phase_circular_variance_ci.png`：反射通道 circular variance 及 bootstrap 95% 区间，区间较宽说明 16 seeds 对相位高阶统计仍有限。
- `08_seed_stability.png`：展示 N=2、4、8、16 时 bootstrap 统计。RMS 展宽较稳定，相干代理和相位方差收敛更慢。
- `09_representative_angular_spectrum_maps.png`：6000 Hz 的集合平均二维反射角谱。Kirchhoff 粗糙海面接近全窗口扩散，SSA1 保留更集中的低横向波数分布。
- `10_radial_angular_spectrum_profiles.png`：二维角谱的径向积分版本，便于定量比较谱能量向高横向波数迁移。
- `11_trend_compatibility.png`：SSA1 散射预算与 Kirchhoff 展宽的描述性秩相关；只表示趋势同向，不表示数值等价。
- `12_ssa1_applicability_parameters.png`：展示 \(k_0\sigma_\eta\)、coherent loss 和有限网格 PM RMS slope。未使用无文献依据的硬阈值。

## 核心公式

$$\sum W_\eta\Delta K_x\Delta K_y=\sigma_\eta^2,\quad \sigma_\eta=H_s/4.$$

$$R_{\rm coh}=R_0\exp(-2k_0^2\sigma_\eta^2).$$

$$G_{\rm SSA1}=4\gamma_s\gamma_i.$$

$$E_{\rm coh}+E_{\rm sca}\le E_{\rm inc}.$$

Kirchhoff 的共同相干指标使用

$$L_{\rm coh}=|\mathbb{E}[h_{\rm reflect}]|/|h_{\rm flat}|,$$

而不是 $\mathbb{E}[|h_{\rm reflect}|]$。

## 验证记分表

| 指标 | 数值 | 验收条件 |
| --- | ---: | --- |
| coherent_formula_max_abs_error | 0 | <= 1e-12 |
| PM_variance_max_relative_error | 8.67362e-15 | <= 1e-10 |
| dense_vs_periodic_FFT_relative_error | 2.52468e-15 | <= 1e-10 |
| SSA1_max_energy_conservation_error | 9.93216e-16 | <= 1e-12 |
| Hs0_SSA1_Kirchhoff_max_channel_difference | 2.22051e-16 | <= 1e-10 |
| max_channel_invariant_error | 1.55158e-17 | <= 1e-10 |

## 趋势相容性

| 频率 / Hz | Spearman rho | Kirchhoff 最大海况展宽 / rad m^-1 |
| ---: | ---: | ---: |
| 4000 | 1.0000 | 3.4532 |
| 6000 | 0.9747 | 3.2032 |
| 8000 | 0.9747 | 3.1521 |
| 10000 | 0.9747 | 3.1419 |

这些相关系数只描述有限测试点上的趋势相容性，不证明两个模型等价。

## 结论边界

- 已验证公式退化、PM 归一化、FFT 实现、能量约束和 16-seed 统计趋势。
- Kirchhoff 的角谱展宽通常强于当前 SSA1，强海况反射幅度也存在明显差异。
- SSA1 在中强海况下可能因能量限制出现散射谱和反射统计饱和，不能将平台误读为真实海面散射已经不再增强。
- `surface_ssa_scatter_scale` 仍是工程归一化参数，不是实验标定截面常数。
- $H_s=1$ m 是压力测试点；离散 RMS slope 仅作诊断，尚未证明其满足一阶 SSA 适用条件。
- 尚未验证斜入射 benchmark、绝对散射截面和物理跨频率相关模型。
