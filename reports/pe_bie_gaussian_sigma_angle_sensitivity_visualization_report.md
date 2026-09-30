# Gaussian 入射角谱与 PE/FK/BIE 误差可视化报告

状态：**只做已有数据核查、指标重算与绘图；没有执行 PE、Full Kirchhoff 或 BIE 求解，也没有修改任何物理模型。**

## 1. 实验设置和数据来源

读取 `results/validation/pe_bie_full_kirchhoff_gaussian_sigma_sensitivity/gaussian_sigma_sensitivity.mat` 与 `results/validation/pe_bie_full_kirchhoff_gaussian_sigma_kzaware/gaussian_sigma_kzaware_comparison.mat`。预期 3 个海面场景 × 5 个 sigma 共 15 组；MAT 中已找到全部 15 组四种场（Model-0、Model-1、Full Kirchhoff、BIE）的 receiver 复场。

环境为 2D/one-transverse，4 kHz，均匀 `c=1500 m/s`，Tx=100 m、Rx=3 m、压力释放海面，sigma=`0.3, 0.5, 1, 2, 5 m`。场景：weak_low_K (`A=0.05 m, K=0.10 rad/m`)、strong_height_low_K (`A=0.20 m, K=0.10 rad/m`)、weak_high_K (`A=0.02 m, K=0.47 rad/m`)。

这些 MAT 中的接收复场是无量纲 rough/flat reflected-response ratio，并按既有比较约定处理 FK/BIE 的共轭；不是未归一化的绝对声压。源数据与原仿真结果均保持未改。

### 评价 mask 与指标定义

已有 sensitivity / kz-aware 报告使用 sigma-specific M99 与逐模型对 BIE 的 pairwise -40 dB 门槛；不同模型 Ec 可能来自不同的接收点集合，不应横向视为完全同采样区域。本报告另外从已存复场重算统一指标：每个 sigma 使用其已存 M99 mask，取所有四个复场均有限且各自幅度不低于四场共同峰值的 -40 dB 的交集，并在交集上重新归一化原 M99 energy weights。mask 随 sigma 仍可能变化，故跨 sigma 是各自有效 footprint 上的一致口径，不是固定空间区域的严格比较。

统一复场误差 `E_c=sqrt(sum(w|u-u_BIE|^2)/sum(w|u_BIE|^2))`；幅度误差 `E_a=sqrt(sum(w(|u|-|u_BIE|)^2)/sum(w|u_BIE|^2))`；phase RMS 为 `sqrt(sum(w angle(u conj(u_BIE))^2))`。完整的新旧 15 组指标见 CSV。

图件均导出为矢量 PDF 和 300 dpi PNG，存于 `results/visualization/pe_bie_gaussian_sigma_angle_sensitivity/figures`。

## 2. 入射角谱随 sigma 的变化

![01_incident_angle_quantiles](../results/visualization/pe_bie_gaussian_sigma_angle_sensitivity/figures/01_incident_angle_quantiles.png)

**图1：同一 Gaussian 源的传播角谱 theta90、theta95、theta99 累积能量角度分位数随 sigma 增大而下降；三个海面场景使用相同源谱，因此角度分位数相同。**

[PDF 矢量版](../results/visualization/pe_bie_gaussian_sigma_angle_sensitivity/figures/01_incident_angle_quantiles.pdf)

数值从宽角端至窄角端，theta99 为 `21.28` 到 `1.23 deg`；theta90 为 `13.43` 到 `0.7826 deg`。角度统计只计传播谱分量，曲线没有插值或平滑。

## 3. Model-0、Model-1、Full Kirchhoff 总体精度

![02_complex_error_by_sigma](../results/visualization/pe_bie_gaussian_sigma_angle_sensitivity/figures/02_complex_error_by_sigma.png)

**图2（核心）：三联图比较三个场景的 Model-0、kz-aware Model-1 和 Full Kirchhoff 相对 BIE 的统一 mask 复场误差；共用 log y 轴和一致的模型色/线型。**

[PDF 矢量版](../results/visualization/pe_bie_gaussian_sigma_angle_sensitivity/figures/02_complex_error_by_sigma.pdf)

在 weak_low_K、sigma=0.3 m，统一 mask 下 Model-1 相对 Model-0 的 Ec 改善为 `4.03` 倍；三场景、sigma=5 m 时 Ec(M0)/Ec(M1) 分别为 `0.997 / 0.999 / 1`。Full Kirchhoff 的同口径误差范围为 `5.37e-06` 至 `0.000271`。

## 4. kz-aware 改善倍数

![03_kz_aware_improvement_vs_theta99](../results/visualization/pe_bie_gaussian_sigma_angle_sensitivity/figures/03_kz_aware_improvement_vs_theta99.png)

**图3：改善倍数定义为同一 common mask 上 Ec(Model-0,BIE)/Ec(Model-1,BIE)；水平线 1 表示没有改善，低于 1 表示 Model-1 略差。**

[PDF 矢量版](../results/visualization/pe_bie_gaussian_sigma_angle_sensitivity/figures/03_kz_aware_improvement_vs_theta99.pdf)

在宽角 sigma=0.3 m，三个 case 的改善倍数依次为 `4.03, 1.59, 1.06`；在 sigma=5 m 分别为 `0.997, 0.999, 1`。宽角下收益具有场景依赖性；缩窄角谱后多数比值趋近 1，不应据此推断误差全都消失。

## 5. 相位与幅度误差分解

![04_phase_error_by_sigma](../results/visualization/pe_bie_gaussian_sigma_angle_sensitivity/figures/04_phase_error_by_sigma.png)

**图4：三联图展示两种 PE 与 Full Kirchhoff 的加权相位 RMS；对数纵轴突出相位误差数量级。**

[PDF 矢量版](../results/visualization/pe_bie_gaussian_sigma_angle_sensitivity/figures/04_phase_error_by_sigma.pdf)

![05_magnitude_error_by_sigma](../results/visualization/pe_bie_gaussian_sigma_angle_sensitivity/figures/05_magnitude_error_by_sigma.png)

**图5：三联图展示幅度相对 L2；重点可检查 weak_high_K 随 sigma 增大后的幅度误差，而不把相位改善误认为幅度也改善。**

[PDF 矢量版](../results/visualization/pe_bie_gaussian_sigma_angle_sensitivity/figures/05_magnitude_error_by_sigma.pdf)

weak_low_K 中随角谱变窄，Model-0/Model-1 phase RMS 整体显著下降，支持有限角度项的存在。strong_height_low_K 仍保留较大残差。weak_high_K 在窄角端仍有可见幅度误差；因此其剩余差异不能仅用入射角宽度解释。该组数据本身不能唯一判定残差是局部边界近似、有限窗/离散误差或其他模型因素。

## 6. 代表性接收场分析

![07_receiver_weak_low_K_sigma_0.3](../results/visualization/pe_bie_gaussian_sigma_angle_sensitivity/figures/07_receiver_weak_low_K_sigma_0.3.png)

**代表性接收线图：幅度和 wrapped phase 对比四种方法，下排给出 Model-0/BIE、Model-1/BIE 的逐点幅度(dB)与 wrapped phase 差；只画 common M99/-40 dB 有效点，坐标 gaps 表示被 mask 的样本。(weak_low_K, sigma=0.3 m)。**

[PDF 矢量版](../results/visualization/pe_bie_gaussian_sigma_angle_sensitivity/figures/07_receiver_weak_low_K_sigma_0.3.pdf)

![07_receiver_strong_height_low_K_sigma_1](../results/visualization/pe_bie_gaussian_sigma_angle_sensitivity/figures/07_receiver_strong_height_low_K_sigma_1.png)

**代表性接收线图：幅度和 wrapped phase 对比四种方法，下排给出 Model-0/BIE、Model-1/BIE 的逐点幅度(dB)与 wrapped phase 差；只画 common M99/-40 dB 有效点，坐标 gaps 表示被 mask 的样本。(strong_height_low_K, sigma=1.0 m)。**

[PDF 矢量版](../results/visualization/pe_bie_gaussian_sigma_angle_sensitivity/figures/07_receiver_strong_height_low_K_sigma_1.pdf)

![07_receiver_weak_high_K_sigma_5](../results/visualization/pe_bie_gaussian_sigma_angle_sensitivity/figures/07_receiver_weak_high_K_sigma_5.png)

**代表性接收线图：幅度和 wrapped phase 对比四种方法，下排给出 Model-0/BIE、Model-1/BIE 的逐点幅度(dB)与 wrapped phase 差；只画 common M99/-40 dB 有效点，坐标 gaps 表示被 mask 的样本。(weak_high_K, sigma=5.0 m)。**

[PDF 矢量版](../results/visualization/pe_bie_gaussian_sigma_angle_sensitivity/figures/07_receiver_weak_high_K_sigma_5.pdf)

相位差采用 `angle(u_model conj(u_BIE))`，限制在 `[-pi,pi]`；没有对曲线做 unwrap，以免把真实的 wrapped residual 改写成累计相位。每个图标题列出 mask 样本数和接收横向范围。

## 7. 全部 15 组总结热力表

![06_all_case_metric_heatmap](../results/visualization/pe_bie_gaussian_sigma_angle_sensitivity/figures/06_all_case_metric_heatmap.png)

**图6：左侧颜色表示 log10(Ec)，单元格标注原始 Ec；右侧以独立色阶显示 improvement，单元格标注 M0/M1 倍数。每行是一个场景和 sigma。**

[PDF 矢量版](../results/visualization/pe_bie_gaussian_sigma_angle_sensitivity/figures/06_all_case_metric_heatmap.pdf)

完整数据可读 `results/visualization/pe_bie_gaussian_sigma_angle_sensitivity/gaussian_sigma_common_mask_summary.csv`；其中同时保留 common-mask 新指标与旧 pairwise-mask 指标，并有 `n_common` 及 mask x 边界供审计。

## 8. 结论与适用范围

1. 宽角谱端 weak_low_K 的 kz-aware 改善约 `4.03` 倍；strong-height 约 `1.59` 倍，weak-high-K 约 `1.06` 倍，体现明显场景依赖。
2. sigma 增大时传播角谱变窄，Model-0 中有限入射角相关的相位误差下降；两种 PE 的总误差逐渐接近，表示 kz 修正的增益收缩到很小，并不代表它们都与 BIE 完全一致。
3. sigma=5 m 时，common-mask Ec(M0/M1/FK) 的逐场景数值如下：

| case | M0 | M1 | Full Kirchhoff |
|---|---:|---:|---:|
| weak_low_K | 0.01133 | 0.01136 | 7.244e-06 |
| strong_height_low_K | 0.0764 | 0.07648 | 2.897e-05 |
| weak_high_K | 0.02885 | 0.02886 | 9.076e-05 |

4. weak_low_K 的宽角误差主要有可被 kz-aware 削弱的相位部分；strong-height 仍有相位残差；weak_high_K 窄角端的幅度残差仍不可忽略。
5. 现有数据足以支持开展 receiver-depth sensitivity test：已有 15 个源宽度/场景的完整复场和明确 mask，可将接收深度作为下一项受控维度；但需要锁定统一指标定义并逐深度记录有效 mask，避免跨深度误读为同一采样区域。该建议是可行性判断，不表示本轮已执行 depth sweep。

### 原有指标与本轮统一指标的区分

原始报告的 pairwise mask 指标来自原始 MAT 中 `metrics_Model0_BIE`、`metrics_Model1_BIE`、`metrics_FK_BIE`；本报告的曲线、表和改善倍数全部由保存的四种原始复场在同一个四场交集 mask 上重算。CSV 中的 `old_*_pairwise` 列仅供数值核对，不能与 common-mask 列混为同一评价口径。原始报告的 mask 定义为：Sigma-specific 99-percent source-energy footprint (strict-normal footprint when available; source-path footprint for sigma=5) plus pairwise -40 dB threshold; weights renormalized on that footprint

## 文件

- 绘图脚本：`scripts/reporting/plot_gaussian_sigma_angle_sensitivity_atlas.m`
- 完整 15 组 CSV：`results/visualization/pe_bie_gaussian_sigma_angle_sensitivity/gaussian_sigma_common_mask_summary.csv`
- 可审计 MAT：`results/visualization/pe_bie_gaussian_sigma_angle_sensitivity/gaussian_sigma_angle_sensitivity_visualization.mat`
- 数据源：`results/validation/pe_bie_full_kirchhoff_gaussian_sigma_sensitivity/gaussian_sigma_sensitivity.mat`、`results/validation/pe_bie_full_kirchhoff_gaussian_sigma_kzaware/gaussian_sigma_kzaware_comparison.mat`
