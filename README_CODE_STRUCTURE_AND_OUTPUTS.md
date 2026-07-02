# 核心代码结构与输出结果说明

本文档面向需要快速读懂当前项目的人：先说明代码文件怎么分层，再说明主要脚本做什么、输出文件里有什么、这些结果应该如何解释。

当前仓库是一个垂直水声信道与 QPSK 通信仿真项目。主线是用 PE/WAPE 风格传播模型得到频域信道 `H_f`，再在此基础上做粗糙海面 Kirchhoff 相位屏诊断、Monte Carlo 统计、通信性能统计，以及基于 Monte Carlo 结果的经验随机信道生成器原型。

## 1. 推荐阅读顺序

1. 先读 `vertical_channel_model.m`：这是最稳定的信道 API 入口。
2. 再读 `vertical_wape_propagator.m`：这是 direct path、surface reflection 和频率循环所在位置。
3. 再读 `pm_surface_boundary_model.m`：这是粗糙海面 Kirchhoff 边界屏、k-domain 接口和谱诊断所在位置。
4. 需要通信链路时读 `comm_main_vertical_psk.m`、`monte_carlo_comm_surface_psk_vertical.m` 和 `monte_carlo_comm_surface_psk_representative_vertical.m`。
5. 需要统计信道生成器时读 `build_surface_empirical_channel_model_vertical.m`、`sample_surface_empirical_channel_vertical.m` 和 `validate_surface_empirical_channel_generator_vertical.m`。

## 2. 核心运行链路

典型信道调用是：

```matlab
paramsV = struct();
paramsV.enable_wideband = true;
paramsV.f_band_hz = [4000 8000];
paramsV.Nf_min = 16;
paramsV.Nf_max = 16;
paramsV.f_ref_hz = 6000;
paramsV.enable_surface_reflection = true;
paramsV.surface_boundary_model = 'kirchhoff_kdomain';
paramsV.show_figures = false;

channel = vertical_channel_model(paramsV);
```

数据流可以概括为：

```text
paramsV
  -> vertical_channel_model.local_prepare_config
  -> vertical_wape_propagator frequency loop
  -> direct path H_direct_f
  -> optional surface reflection H_reflect_f
  -> H_f = H_direct_f + H_reflect_f
  -> h_total = H_f(idx_f_ref)
```

通信链路再把 `H_f` 转成 baseband 频响和 tap，再做发射、信道卷积、加噪、同步/均衡、解调和 BER/SER 统计。当前 D2 后通信接收窗口采用明确的峰值同步策略，Eb/N0 噪声参考采用接收端 clean signal 功率口径。

## 3. 核心代码文件

### 3.1 信道传播核心

`vertical_channel_model.m`

- 公共信道 API。
- 负责把用户输入 `paramsV` 变成运行配置 `cfg`。
- 负责默认值、参数校验和输出结构整理。
- 必须保持对外输出字段兼容。

`vertical_wape_propagator.m`

- 传播核心和频率循环。
- 生成 `H_direct_f`、`H_reflect_f`、`H_f`。
- 保持核心不变量：`H_f = H_direct_f + H_reflect_f`。
- 反射开启时走 `tx -> surface -> boundary operator -> rx`。
- 宽带运行中，粗糙海面诊断主要保存在参考频点 `idx_f_ref` 的 metadata 中，避免每个频点都保存大型诊断。

`pm_surface_boundary_model.m`

- 粗糙海面 Kirchhoff 相位屏模块。
- 支持旧默认空间域模型 `surface_boundary_model='kirchhoff_spatial'`。
- 支持 k-domain 接口 `surface_boundary_model='kirchhoff_kdomain'`。
- k-domain 第一版不构造四维稠密矩阵，而是用 FFT 乘积-卷积等价实现：

```matlab
G_xy = surface_reflect_coeff .* exp(1i * delta_phi);
Psi_inc_k = fft2(psi_inc_xy);
Psi_ref_k = fft2(G_xy .* ifft2(Psi_inc_k));
psi_ref_xy = ifft2(Psi_ref_k);
```

- 该写法表达的是离散周期横向网格上的隐式边界算子：

```text
Psi_ref(K) = B_xi[Psi_inc](K)
B_xi(K,K') proportional to G_hat(K-K')
```

- 当前实现仍是同一个 Kirchhoff 相位屏的接口重写，不是严格 T 矩阵、SSA/NLSSA，也不是新的粗糙面散射理论。

### 3.2 通信链路核心

`comm_main_vertical_psk.m`

- 单次端到端通信 demo。
- 消费 `vertical_channel_model(paramsV)` 返回的 `H_f`。
- 构造 baseband 频响和 tap，完成 QPSK 调制、信道、噪声、均衡、解调。

`modem_psk.m`

- PSK 调制/解调 helper。
- 当前通信统计默认使用 QPSK，即 `M=4`。

`noise_inject_vertical.m`

- AWGN 注入 helper。
- D2 后通信统计明确区分 Eb/N0 噪声参考功率口径，默认使用接收端 clean signal 功率。

`validate_comm_link_minimal_vertical.m`

- D1 最小闭环验证脚本。
- 用单位信道、单抽头复增益、已知短多径和 PE 生成 tap 检查调制、加噪、均衡和解调是否工作。

### 3.3 粗糙海面诊断与 Monte Carlo

`sweep_surface_boundary_redistribution_vertical.m`

- C2.5 频率、海况趋势诊断脚本。
- 重点输出入射谱到反射谱的实际谱再分布指标。

`monte_carlo_surface_channel_vertical.m`

- C3 固定海况、多 `sea_seed` 的 Monte Carlo 信道统计。
- 输出 `H_f`、参考频点信道、反射/直达比、C2/C2.5 诊断和 tap 摘要的经验统计。

`sweep_monte_carlo_surface_channel_vertical.m`

- C3.5 多海况 Monte Carlo 扫描。
- 默认扫描 `sea_hs_target = [0.05, 0.5, 1.0]`、`sea_wind_speed = [3, 5, 8, 12]`。
- 当前阶段已使用 `mc_count=16` 生成 12 个海况条件、192 次传播的统计结果。

`plot_c35_core_heatmaps_vertical.m`

- 从 C3.5 结果文件生成 5 类核心图：
  - `|H(f_ref)|` 均值热力图。
  - `|H(f_ref)|` 标准差热力图。
  - `reflect_rms_delta_k` 热力图。
  - `reflect_high_k_fraction` 热力图。
  - `tap_rms_delay` 热力图。

`monte_carlo_comm_surface_psk_vertical.m`

- 通信性能 Monte Carlo 统计脚本。
- 默认比较 weak/strong 海况下的 `direct_only` 与 `direct_plus_reflect`。
- 输出 BER/SER 曲线的均值、方差和分位数。

`monte_carlo_comm_surface_psk_representative_vertical.m`

- 代表性海况通信 Monte Carlo 脚本。
- 当前用于 weak、mid、strong 三个代表海况的通信性能统计。

### 3.4 C4 经验随机信道生成器

`build_surface_empirical_channel_model_vertical.m`

- 读取 C3/C3.5 Monte Carlo `.mat` 结果。
- 构建经验统计模型表。
- 不调用 `vertical_channel_model`，不运行 PE/WAPE。

`sample_surface_empirical_channel_vertical.m`

- 从经验模型中抽样。
- 支持最近邻海况匹配。
- 当前支持三类模式：
  - `summary`：低维摘要样本。
  - `wideband_hf`：从已有 Monte Carlo 样本中 bootstrap 宽带 `H_f`。
  - `tap_level`：从已有样本抽取 tap-level 摘要。

`validate_surface_empirical_channel_generator_vertical.m`

- 验证 C4 生成器。
- 读取 C3.5 结果作为统计来源，生成快速样本并比较均值、标准差和分位数。
- 不运行完整 PE 传播。

## 4. 主要输出字段解释

`vertical_channel_model(paramsV)` 的核心输出字段如下。

`output.H_f`

- 频域总信道。
- 每个频率点都满足：

```text
H_f = H_direct_f + H_reflect_f
```

`output.H_direct_f`

- 直达传播分量。
- 在 uniform/no-bubble 的设置下，跨不同 `sea_seed` 应保持 roundoff 级一致。

`output.H_reflect_f`

- 海面反射分量。
- 当 `enable_surface_reflection=false` 时应为 0。

`output.f_axis`

- 宽带频率轴，单位 Hz。
- 标量单频时通常只有一个频点。

`output.idx_f_ref`

- 参考频点在 `f_axis` 中的索引。
- `h_total`、`h_direct`、`h_reflect` 对应这个频点。

`output.h_total`

- 参考频点复数总信道：

```text
h_total = H_f(idx_f_ref)
```

`output.h_direct`

- 参考频点直达分量：

```text
h_direct = H_direct_f(idx_f_ref)
```

`output.h_reflect`

- 参考频点反射分量：

```text
h_reflect = H_reflect_f(idx_f_ref)
```

`output.roughness_meta`

- 粗糙海面与边界算子 metadata。
- 反射关闭时 `roughness_meta.enabled=false`，但新增诊断字段仍保留并标记为未执行。

## 5. `roughness_meta` 诊断字段

### 5.1 边界模型 metadata

`boundary_model`

- 当前边界模型。
- 默认是 `kirchhoff_spatial`。
- 可选 `kirchhoff_kdomain`。

`boundary_operator_form`

- 说明当前是空间域乘法，还是 k-domain 隐式卷积算子接口。

`boundary_dense_matrix_used`

- 当前应为 `false`。
- 第一版 k-domain 不构造四维稠密矩阵。

`boundary_fft_convention`

- MATLAB FFT 约定：
  - `fft2` 前向不归一化。
  - `ifft2` 包含 `1/(nx*ny)`。
  - 横向波数按周期网格解释。

`boundary_equivalence_error`

- 当启用等价性检查时，记录 spatial 与 kdomain 两条路径的 roundoff 级误差。

### 5.2 C2 屏函数耦合诊断

`boundary_coupling_diagnostics`

- 只分析边界屏函数：

```text
G_xy = R_s exp(i delta_phi)
G_hat_k = fft2(G_xy)
P_k = |G_hat_k|^2
```

- 典型指标：
  - `total_energy`
  - `zero_energy`
  - `nonzero_power_fraction`
  - `rms_delta_k_rad_per_m`
  - `energy_radius_90_rad_per_m`
  - `nonzero_energy_radius_90_rad_per_m`

这些指标说明屏函数本身在横向波数域的扩展程度。`nonzero_power_fraction` 越高、`rms_delta_k` 越大，说明 Kirchhoff 屏函数越不均匀，隐式 `G_hat(K-K')` 卷积核中潜在非对角耦合越强。但它们不是严格散射矩阵元，也不能解释为真实散射截面。

### 5.3 C2.5 入射加权谱再分布诊断

`boundary_redistribution_diagnostics`

- 分析当前入射谱经过粗糙屏后实际变成怎样的反射谱：

```text
Psi_inc_k = fft2(psi_inc_xy)
Psi_ref_k = fft2(psi_ref_xy)
Psi_flat_ref_k = fft2(R_s .* psi_inc_xy)
```

- 典型指标：
  - `incident_rms_delta_k_rad_per_m`
  - `reflect_rms_delta_k_rad_per_m`
  - `rms_delta_k_increase_rad_per_m`
  - `incident_centroid_kx_rad_per_m`
  - `reflect_centroid_kx_rad_per_m`
  - `centroid_shift_mag_rad_per_m`
  - `reflect_high_k_fraction`
  - `rough_vs_flat_rms_delta_k_increase_rad_per_m`

这类指标比 C2 更接近“当前入射场实际被怎样重分布”，因为它同时依赖 `G_xy` 和 `psi_inc_xy`。它仍然只是 Kirchhoff 相位屏下的谱诊断，不是 T 矩阵、SSA/NLSSA 或真实粗糙面散射截面。

## 6. C3.5 结果如何看

主要结果文件：

```text
sweep_monte_carlo_surface_channel_vertical_result.mat
```

核心结构：

`base_params`

- 本次 Monte Carlo 扫描使用的 reduced-grid 基础参数。

`sea_hs_values`

- 扫描的有效波高列表。

`sea_wind_values`

- 扫描的风速列表。

`seed_list`

- 每个海况使用的 `sea_seed` 列表。

`run_summary_table`

- 每行对应一个 `(Hs, wind, seed)` realization。
- 适合排查单次运行的 invariant error、诊断量、tap 指标。

`condition_summary_table`

- 每行对应一个 `(Hs, wind)` 条件。
- 已对该条件下多个 seed 做均值、标准差、分位数等聚合。

`condition_results`

- 每个海况条件的 compact Monte Carlo 结果。
- 包含小型 `H_f` 样本矩阵和统计量，不保存完整空间场。

`sweep_stats`

- 按 `[numel(Hs) x numel(wind)]` 整理的趋势矩阵。
- 适合画热力图。

当前 `mc_count=16` reduced C3.5 结果：

- 海况条件数：12。
- 总传播次数：192。
- `max(abs(H_f-(H_direct_f+H_reflect_f))) = 1.55e-17`。
- direct path 跨 seed 最大漂移约 `1.39e-17`。

5 类核心图文件：

```text
c35_core_abs_H_ref_mean_heatmap.png
c35_core_abs_H_ref_std_heatmap.png
c35_core_reflect_rms_delta_k_mean_heatmap.png
c35_core_reflect_high_k_fraction_mean_heatmap.png
c35_core_tap_rms_delay_mean_heatmap.png
```

读图方式：

- `|H(f_ref)|` 均值热力图：看参考频点平均信道强度。
- `|H(f_ref)|` 标准差热力图：看海面 realization 引入的信道起伏强度。
- `reflect_rms_delta_k` 热力图：看反射谱横向波数扩展程度。
- `reflect_high_k_fraction` 热力图：看反射谱中高横向波数能量比例。
- `tap_rms_delay` 热力图：看等效 baseband tap 的时延扩展。

这些图描述的是 reduced-grid、有限 seed、当前 Kirchhoff 相位屏模型下的经验趋势。

## 7. 通信 Monte Carlo 结果如何看

主要脚本：

```text
monte_carlo_comm_surface_psk_vertical.m
monte_carlo_comm_surface_psk_representative_vertical.m
```

主要结果文件：

```text
monte_carlo_comm_surface_psk_representative_result.mat
```

核心结构：

`sea_conditions`

- 通信统计使用的海况列表。
- 当前代表性运行使用 weak、mid、strong 三个海况。

`scenarios`

- 传播场景。
- `direct_only`：关闭海面反射。
- `direct_plus_reflect`：打开海面反射并使用 `kirchhoff_kdomain`。

`EbN0_dB_list`

- BER/SER 曲线横轴。

`run_summary_table`

- 每行对应一个 `(sea_condition, scenario, seed)` 的信道运行。

`curve_summary_table`

- 每行对应一个 `(sea_condition, scenario, EbN0)` 的 BER/SER 聚合结果。

`curve_stats`

- BER、SER、effective SNR 的均值、方差、标准差和分位数。

当前代表性通信 Monte Carlo 结果：

- 海况：weak `(Hs=0.05, wind=5)`、mid `(Hs=0.5, wind=8)`、strong `(Hs=1.0, wind=12)`。
- 场景：`direct_only` 与 `direct_plus_reflect`。
- 总信道运行：96。
- 最大 `H_f` 不变量误差：`1.39e-17`。
- direct-only 下最大反射分量：0。

结果解读重点：

- direct-only 是通信链路基线。
- direct_plus_reflect 显示海面反射和多径 tap 对 BER/SER 的影响。
- 中高 Eb/N0 下 BER/SER 不一定严格单调，因为默认符号数和 seed 数有限，零误码或个别非单调点属于有限样本现象。
- 当前通信统计是 reduced 实验，不是闭式 BER 公式。

## 8. C4 经验随机信道生成器如何看

C4 的目标是从已有 Monte Carlo 样本中快速生成经验样本，而不是替代 PE/WAPE 传播。

主要输入文件：

```text
sweep_monte_carlo_surface_channel_vertical_result.mat
```

主要输出文件：

```text
surface_empirical_channel_generator_validation_result.mat
```

支持的输入参数：

- `sea_hs_target`
- `sea_wind_speed`
- `f_ref_hz`
- `n_samples`
- `sample_mode`
- `rng_seed`

当前 v1 采样策略：

- 最近邻匹配海况。
- 从匹配海况的 Monte Carlo seed 样本中有放回 bootstrap。
- 不做海况插值。

可生成的样本类型：

`summary`

- `|H(f_ref)|`
- `angle(H(f_ref))`
- `h_reflect/h_direct`
- `reflect_rms_delta_k`
- `reflect_high_k_fraction`
- `tap_rms_delay_symbols`

`wideband_hf`

- 从 C3.5 中已有的 `H_f` 样本矩阵 bootstrap。
- 可用于快速得到低成本的宽带频响样本。
- 仍然只是在已有 Monte Carlo 样本上的经验重采样。

`tap_level`

- 抽取 tap 摘要样本。
- 当前保存 tap 指标，不保存完整空间场。

C4 限制：

- 只覆盖已有 Monte Carlo 扫描范围。
- 有限 seed 下的 bootstrap 会复用已有样本，不会创造新的物理 realization。
- 不是严格闭式统计信道模型。
- 不是 T 矩阵、SSA/NLSSA 或新的散射理论。
- 不应把 C4 输出当作 PE/WAPE 传播结果的替代品。

## 9. 图像输出说明

当前目录中的 PNG 图大致分为六类：C3.5 多海况信道统计图、C3 固定海况 Monte Carlo 图、C2.5 谱再分布趋势图、通信 BER/SER 图、C4 经验生成器验证图，以及历史 bubble/校准图。读图时要先看文件名前缀，确认图来自哪个脚本和哪一类实验。

### 9.1 C3.5 多海况信道热力图

核心文件：

```text
c35_core_abs_H_ref_mean_heatmap.png
c35_core_abs_H_ref_std_heatmap.png
c35_core_reflect_rms_delta_k_mean_heatmap.png
c35_core_reflect_high_k_fraction_mean_heatmap.png
c35_core_tap_rms_delay_mean_heatmap.png
```

这些图由 `plot_c35_core_heatmaps_vertical.m` 从 `sweep_monte_carlo_surface_channel_vertical_result.mat` 读取 compact statistics 生成，不会重新运行 PE/WAPE。

坐标含义：

- 横轴：`sea_wind_speed`，单位 m/s。
- 纵轴：`sea_hs_target`，单位 m。
- 每个格点：一个固定海况 `(Hs, wind)` 下，对多个 `sea_seed` 的 Monte Carlo 聚合值。
- 颜色和格点文字：该统计量的数值。

各图含义：

`c35_core_abs_H_ref_mean_heatmap.png`

- 显示参考频点 `|H(f_ref)|` 的 Monte Carlo 均值。
- 颜色越大，说明该海况下参考频点平均信道幅度越强。
- 它反映的是 direct 与 reflect 合成后的总信道幅度，不单独表示反射强度。

`c35_core_abs_H_ref_std_heatmap.png`

- 显示参考频点 `|H(f_ref)|` 跨 `sea_seed` 的标准差。
- 颜色越大，说明不同随机海面 realization 造成的信道幅度起伏越强。
- 该图适合观察信道随机性，而不是平均链路预算。

`c35_core_reflect_rms_delta_k_mean_heatmap.png`

- 显示 C2.5 反射谱 `|Psi_ref(K)|^2` 的 rms 横向波数均值，单位 rad/m。
- 颜色越大，说明反射后的角谱更宽，当前入射场被 Kirchhoff 粗糙屏重分布到更大的横向波数范围。
- 这仍是相位屏谱诊断，不是严格散射截面。

`c35_core_reflect_high_k_fraction_mean_heatmap.png`

- 显示反射谱中高横向波数能量占比的均值。
- 默认高 k 阈值来自入射谱 90% 能量半径。
- 颜色越大，说明反射后有更多能量进入入射主谱宽之外的横向波数区域，潜在非对角谱再分布更明显。

`c35_core_tap_rms_delay_mean_heatmap.png`

- 显示 baseband tap RMS delay 的 Monte Carlo 均值，单位 symbol。
- 颜色越大，说明等效时域信道能量分布更分散，多径/频率选择性更明显。
- 该指标受 tap 构造、频带宽度、同步策略和 reduced-grid 频率采样共同影响，不能单独解释为真实物理时延扩展。

完整 C3.5 sweep 脚本还会生成同类文件：

```text
sweep_monte_carlo_surface_channel_vertical_abs_H_ref_mean_heatmap.png
sweep_monte_carlo_surface_channel_vertical_abs_H_ref_std_heatmap.png
sweep_monte_carlo_surface_channel_vertical_reflect_rms_delta_k_mean_heatmap.png
sweep_monte_carlo_surface_channel_vertical_reflect_high_k_fraction_mean_heatmap.png
sweep_monte_carlo_surface_channel_vertical_rms_delta_k_increase_mean_heatmap.png
sweep_monte_carlo_surface_channel_vertical_tap_rms_delay_mean_heatmap.png
sweep_monte_carlo_surface_channel_vertical_reflect_high_k_fraction_vs_wind.png
sweep_monte_carlo_surface_channel_vertical_reflect_rms_delta_k_vs_wind.png
```

其中 `rms_delta_k_increase_mean_heatmap` 表示反射谱 rms 横向波数相对入射谱的增量；`*_vs_wind.png` 是在不同 `Hs` 曲线下观察指标随风速变化的线图。带 `smoke_` 的文件是低成本 smoke test 结果，只用于确认脚本能跑通，不建议用于阶段性物理结论。

### 9.2 C3 固定海况 Monte Carlo 图

核心文件：

```text
monte_carlo_surface_channel_vertical_H_f_mean_magnitude.png
monte_carlo_surface_channel_vertical_abs_h_ref_hist.png
monte_carlo_surface_channel_vertical_phase_h_ref_hist.png
monte_carlo_surface_channel_vertical_reflect_rms_delta_k_hist.png
monte_carlo_surface_channel_vertical_reflect_high_k_fraction_hist.png
monte_carlo_surface_channel_vertical_tap_rms_delay_hist.png
```

`H_f_mean_magnitude.png`

- 横轴是频率，纵轴是 `|H_f|`。
- 通常显示 Monte Carlo 均值曲线和一倍标准差包络。
- 用于观察整个频带内平均信道幅度和 realization 起伏。

`abs_h_ref_hist.png`

- 参考频点 `|H(f_ref)|` 的直方图。
- 用于观察固定海况、多 seed 下的幅度分布。

`phase_h_ref_hist.png`

- 参考频点 `angle(H(f_ref))` 的主值相位直方图。
- 相位是圆周变量，均值和方差应使用 circular statistics 解读，不能直接当普通线性变量。

`reflect_rms_delta_k_hist.png`

- C2.5 反射谱 rms 横向波数的样本分布。
- 分布越偏向大值，说明该海况下谱展宽更强。

`reflect_high_k_fraction_hist.png`

- 反射高 k 能量占比的样本分布。
- 用于观察不同随机海面导致的非主谱宽能量比例变化。

`tap_rms_delay_hist.png`

- baseband tap RMS delay 的样本分布。
- 用于观察通信等效信道的时延扩展随机性。

### 9.3 C2.5 谱再分布趋势图

核心文件：

```text
sweep_surface_boundary_redistribution_vertical_Hs_rms_broadening.png
sweep_surface_boundary_redistribution_vertical_Hs_high_k_fraction.png
sweep_surface_boundary_redistribution_vertical_wind_rms_broadening.png
sweep_surface_boundary_redistribution_vertical_frequency_rms_broadening.png
```

这些图来自 `sweep_surface_boundary_redistribution_vertical.m`，用于观察单次或 reduced sweep 下的入射加权谱再分布趋势。

`Hs_rms_broadening.png`

- 横轴通常是 `Hs`。
- 纵轴是反射谱 rms 横向波数相对入射谱或平整海面对照的展宽指标。
- 用于观察粗糙度增强时谱宽是否增加。

`Hs_high_k_fraction.png`

- 横轴通常是 `Hs`。
- 纵轴是反射高 k 能量占比。
- 用于观察更强海况是否把更多能量重分布到入射主谱之外。

`wind_rms_broadening.png`

- 横轴通常是风速。
- 纵轴是谱展宽指标。
- 用于观察 PM 海面参数变化下的趋势。

`frequency_rms_broadening.png`

- 横轴通常是频率。
- 纵轴是谱展宽指标。
- 用于观察同一海况下不同声频率的 Kirchhoff 相位屏调制强弱变化。

这些图只说明当前相位屏模型中的谱再分布趋势，不代表严格散射截面随 `Hs`、风速或频率的闭式规律。

### 9.4 通信 BER/SER 图

核心文件：

```text
monte_carlo_comm_surface_psk_vertical_BER_weak_strong_direct_plus_reflect.png
monte_carlo_comm_surface_psk_vertical_SER_weak_strong_direct_plus_reflect.png
monte_carlo_comm_surface_psk_vertical_BER_weak_direct_vs_reflect.png
monte_carlo_comm_surface_psk_vertical_SER_weak_direct_vs_reflect.png
monte_carlo_comm_surface_psk_vertical_BER_strong_direct_vs_reflect.png
monte_carlo_comm_surface_psk_vertical_SER_strong_direct_vs_reflect.png
```

代表性三海况脚本对应：

```text
monte_carlo_comm_surface_psk_representative_BER_representative_direct_plus_reflect.png
monte_carlo_comm_surface_psk_representative_SER_representative_direct_plus_reflect.png
monte_carlo_comm_surface_psk_representative_BER_weak_direct_vs_reflect.png
monte_carlo_comm_surface_psk_representative_BER_mid_direct_vs_reflect.png
monte_carlo_comm_surface_psk_representative_BER_strong_direct_vs_reflect.png
monte_carlo_comm_surface_psk_representative_SER_weak_direct_vs_reflect.png
monte_carlo_comm_surface_psk_representative_SER_mid_direct_vs_reflect.png
monte_carlo_comm_surface_psk_representative_SER_strong_direct_vs_reflect.png
```

读图方式：

- 横轴：`Eb/N0`，单位 dB。
- 纵轴：BER 或 SER，通常用对数坐标更容易观察。
- 曲线：不同海况或不同传播场景。
- 阴影/误差包络若存在：跨 `sea_seed` 的分位数或统计离散程度。

`*_direct_plus_reflect.png`

- 比较不同海况在包含海面反射时的通信性能。
- 曲线越低，误码性能越好。
- weak/mid/strong 的差异反映当前 reduced 信道和同步/均衡策略下，粗糙海面反射对通信链路的经验影响。

`*_direct_vs_reflect.png`

- 在同一海况下比较 `direct_only` 与 `direct_plus_reflect`。
- 若 `direct_plus_reflect` 曲线高于 `direct_only`，说明反射分量、多径和频率选择性在当前接收策略下增加了误码。
- 若两者接近，说明该海况和 reduced 参数下反射影响较弱，或被同步/均衡部分补偿。

带 `smoke_` 的通信图是 smoke test 结果，通常 seed 更少，只用于检查脚本和数据结构，不应用作最终性能结论。

### 9.5 D1 通信链路最小闭环验证图

核心文件：

```text
validate_comm_link_minimal_vertical_simple_channel_BER.png
validate_comm_link_minimal_vertical_simple_channel_SER.png
validate_comm_link_minimal_vertical_pe_h_bb_current_taps.png
validate_comm_link_minimal_vertical_pe_h_bb_aligned_taps.png
```

`simple_channel_BER.png` 和 `simple_channel_SER.png`

- 用单位信道、单抽头复增益和已知短多径测试基础通信链路。
- 目标是确认 QPSK 调制、AWGN、均衡和解调在简单信道下能给出随 Eb/N0 改善的 BER/SER 曲线。

`pe_h_bb_current_taps.png`

- 显示 PE 生成的 baseband tap 在当前窗口下的幅度。
- 用于诊断主峰是否偏离接收窗口起点，是否可能导致采样/均衡问题。

`pe_h_bb_aligned_taps.png`

- 显示诊断性峰值对齐后的 tap 幅度。
- D1 用它证明峰值同步能改善通信闭环；D2 随后把接收窗口固定为明确的峰值同步策略。

### 9.6 C4 经验生成器验证图

核心文件：

```text
surface_empirical_channel_generator_abs_h_ref_hist.png
surface_empirical_channel_generator_phase_h_ref_rad_hist.png
surface_empirical_channel_generator_reflect_direct_abs_hist.png
surface_empirical_channel_generator_redistribution_reflect_rms_delta_k_rad_per_m_hist.png
surface_empirical_channel_generator_redistribution_reflect_high_k_fraction_hist.png
surface_empirical_channel_generator_tap_rms_delay_symbols_hist.png
```

这些图来自 `validate_surface_empirical_channel_generator_vertical.m`，用于比较原始 Monte Carlo 样本和 C4 bootstrap 快速样本。

`abs_h_ref_hist.png`

- 比较 `|H(f_ref)|` 的真实 Monte Carlo 分布与快速样本分布。
- 两者接近说明 bootstrap 能复现参考频点幅度统计。

`phase_h_ref_rad_hist.png`

- 比较 `angle(H(f_ref))` 分布。
- 相位主值存在 `-pi/pi` 跳变，读图时应按圆周变量理解。

`reflect_direct_abs_hist.png`

- 比较 `|h_reflect/h_direct|` 分布。
- 用于观察快速样本是否复现反射相对直达强度的统计。

`redistribution_reflect_rms_delta_k_rad_per_m_hist.png`

- 比较反射谱 rms 横向波数诊断量分布。
- 用于观察 C4 是否复现谱展宽统计。

`redistribution_reflect_high_k_fraction_hist.png`

- 比较反射高 k 能量占比分布。
- 用于观察 C4 是否复现非主谱宽能量比例。

`tap_rms_delay_symbols_hist.png`

- 比较 tap RMS delay 分布。
- 用于观察 C4 是否复现通信等效时延扩展摘要。

C4 图只验证“经验统计复现”能力。它不证明 C4 生成了新的物理海面，也不说明 C4 能外推到未扫描海况。

### 9.7 历史 bubble 与校准图

当前目录还包含若干 bubble model、Hall 1D 校准和敏感性分析图，例如：

```text
compare_bubble_models_vertical_Figure41_H_magnitude.png
compare_bubble_models_vertical_Figure42_H_phase.png
compare_bubble_models_vertical_Figure43_delta_TL.png
compare_bubble_models_vertical_Figure44_components.png
compare_bubble_models_vertical_Figure45_bubble_meta.png
comm_compare_bubble_models_vertical_BER_vs_EbN0.png
comm_compare_bubble_models_vertical_SER_vs_EbN0.png
calibrate_hall1d_bubble_vertical_*.png
sweep_hall1d_sensitivity_vertical_*.png
```

这些图属于历史 bubble / Hall 1D 方向：

- `H_magnitude` / `H_phase`：比较不同 bubble 模型下的信道幅度和相位。
- `delta_TL`：比较 bubble 引起的传输损失变化。
- `components`：拆分 direct、reflect 或不同模型贡献。
- `bubble_meta`：展示 bubble 模型内部参数或诊断摘要。
- `BER_vs_EbN0` / `SER_vs_EbN0`：比较 bubble 模型对通信性能的影响。
- `calibrate_*`：校准 Hall 1D bubble 模型参数与目标量的匹配情况。
- `sweep_hall1d_sensitivity_*`：敏感性分析，观察某个 bubble 参数变化对最大吸收、损失或通信指标的影响。

这些历史图和当前 C1-C4 粗糙海面 Monte Carlo 主线相关，但不是同一组统计实验。引用阶段性结果时应优先使用 C3.5、代表性通信 Monte Carlo 和 C4 图。

### 9.8 基础 vertical demo 图

基础 demo 还会生成：

```text
vertical_upward_4k_uniform_Figure11_xy.png
vertical_upward_4k_uniform_Figure12_xz.png
vertical_upward_4k_uniform_Figure13_yz.png
vertical_upward_4k_uniform_Figure14_1overR.png
vertical_upward_4k_uniform_Figure15_surface_kirchhoff.png
```

这些图用于检查单次传播场：

- `xy`：某个深度或输出平面上的横向场分布。
- `xz` / `yz`：垂直切片场分布。
- `1overR`：uniform medium 中的 `1/R` 行为检查。
- `surface_kirchhoff`：粗糙海面 Kirchhoff 相位屏或相关反射诊断图。

它们更适合做传播模型 sanity check，不是 Monte Carlo 统计结果。

## 10. 常见运行入口

单次信道：

```matlab
channel = vertical_channel_model(paramsV);
```

单次通信 demo：

```matlab
comm_main_vertical_psk
```

C3.5 多海况信道 Monte Carlo：

```matlab
sweep_monte_carlo_surface_channel_vertical
```

生成 C3.5 核心热力图：

```matlab
plot_c35_core_heatmaps_vertical
```

代表海况通信 Monte Carlo：

```matlab
monte_carlo_comm_surface_psk_representative_vertical
```

C4 经验生成器验证：

```matlab
validate_surface_empirical_channel_generator_vertical
```

## 11. 输出文件与版本控制

本仓库 `.gitignore` 当前采用源码优先策略：

- MATLAB 源码 `*.m` 会纳入版本控制。
- Markdown 文档 `*.md` 会纳入版本控制。
- 大型或可再生成结果如 `.mat`、`.png` 默认忽略。

因此，结果文件通常留在本地工作目录中，用文档记录关键参数和代表性数值，不直接提交大型仿真产物。

## 12. 当前阶段性结论

当前代码已经形成 realization-based 仿真平台：

- 可以生成 PE/WAPE 风格的 direct 与 reflected 信道。
- 可以在 Kirchhoff 相位屏上使用 spatial 或 kdomain 边界接口。
- 可以输出屏函数耦合诊断和入射加权谱再分布诊断。
- 可以做固定海况和多海况 Monte Carlo 信道统计。
- 可以做 reduced QPSK 通信性能 Monte Carlo。
- 可以从已有 Monte Carlo 统计结果构建轻量经验随机信道生成器原型。

下一阶段的重点已经从“单个 realization 是否可算、诊断是否可解释”转向“如何基于经验统计构建更有用、更可控的随机信道生成器”。

## 2026-07-02 Addendum: `kirchhoff_kstat`

`surface_boundary_model='kirchhoff_kstat'` adds a Kirchhoff statistical phase-screen branch under `pm_surface_boundary_model.m`.

Core behavior:
- It does not synthesize an explicit sea surface `eta(x,y)`.
- It uses the PM height spectrum to form `C_eta`, `C_deltaG`, and `S_deltaG`.
- It returns a coherent reflected spectrum `Psi_coh_k=R_coh*Psi_inc_k`.
- It forms incoherent scatter power with `P_sca=|R0|^2*circconv(S_deltaG,abs(Psi_inc_k).^2)*dkx*dky/(2*pi)^2`.
- If `surface_kstat_random_scatter=true`, it adds `sqrt(P_sca).*CN(0,1)` as a random reflected spectrum realization.

New public configuration fields:
- `surface_kstat_random_scatter`, default `true`.
- `surface_kstat_seed_offset`, default `200000`.
- `surface_kstat_conv_padding`, default `periodic`, allowed `periodic` or `zero_padded`.
- `surface_kstat_trusted_angle_deg`, default `NaN`.

Metadata:
- `output.roughness_meta.kirchhoff_kstat_meta` stores `sigma_eta2_m2`, `alpha_rad_per_m`, `G_mean`, `R_coh`, `S_deltaG`, `P_sca`, seed fields, FFT normalization notes, phase-screen energy closure, and propagation-window energies.
- Disabled paths still include the `kirchhoff_kstat_meta` struct with `enabled=false`.

Model boundary:
- This branch is a Kirchhoff / Gaussian statistical phase-screen generator.
- It conserves full phase-screen energy through `abs(<G>)^2 + int S_deltaG/(2*pi)^2`.
- The propagation window `K_h<=k0` is diagnostic only and is not renormalized.
- SSA1 remains a weak-roughness/small-angle reference check, not the kstat generation formula.

## 13. 使用结果时的注意事项

1. reduced-grid 结果适合趋势分析和代码验证，不等同于最终高分辨率物理结论。
2. Kirchhoff k-domain 接口是相位屏的波数域接口重写，不是新的粗糙面散射理论。

## 2026-07-02 Addendum: Raw-PM Wind Comparison

- New config: `surface_roughness_scale_mode='target_hs'|'raw_pm'`.
- Default remains `target_hs`, so existing runs still scale PM roughness to `sea_hs_target`.
- `raw_pm` leaves the PM spectrum/realization amplitude determined by `sea_wind_speed`; `sea_hs_target` is retained only as metadata.
- `pm_surface_boundary_model.m` now records `roughness_scale_mode`, `sigma_eta_raw_m`, `Hs_raw_m`, `Hs_target_m`, `scale_factor`, `pm_variance_raw_discrete_m2`, and `pm_variance_raw_continuous_m2`.
- `scripts/comparisons/compare_kirchhoff_kdomain_kstat_wind_vertical.m` compares `kirchhoff_kdomain` and `kirchhoff_kstat` for wind speeds `[3,5,8,10,12,15]` with 32 seeds by default.
- Outputs are written only under `results/comparisons/`: MAT result, CSV summary, and five PNG figures for raw `Hs`, coherent reflection, reflected tap magnitude, incoherent energy, and propagating-window energy.
- Interpretation: explicit `kirchhoff_kdomain` uses the project's existing discrete PM realization convention, while `kirchhoff_kstat` uses the continuous `(2*pi)^-2` phase-screen convention. The comparison table exposes both raw PM variance audits instead of silently forcing them to match.
3. C2/C2.5 诊断指标描述的是当前离散模型中的谱扩展和谱再分布，不是严格散射截面。
4. C3/C3.5 Monte Carlo 是有限 seed 的经验统计，不是闭式随机信道模型。

5. C4 经验生成器依赖已有 Monte Carlo 数据，不能外推到未扫描海况。
6. 通信 BER/SER 是有限符号数和有限 AWGN realization 下的经验结果，高 Eb/N0 下出现 0 BER 是正常有限样本现象。
