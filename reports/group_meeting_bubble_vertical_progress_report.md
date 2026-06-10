# 风生气泡层垂直 PE 水声信道建模与通信验证研究进度汇报
本文中的数值均来自当前目录下已有结果文件，主要包括 `compare_bubble_models_vertical_result.mat`、`sweep_hall1d_sensitivity_vertical_result.mat`、`calibrate_hall1d_bubble_vertical_result.mat`、`diagnose_comm_chain_vertical_result.mat` 和 `comm_compare_bubble_models_vertical_result.mat`。当前结果均为 reduced-grid 验证结果，不是最终高分辨率生产仿真结果。

## 1. 研究进展总述

目前项目已经在原有垂直 PE 水声信道框架上完成了风生气泡层建模和初步通信验证。当前已完成的功能包括：

- 垂直 PE 信道主链路：保留原有高斯声源、垂直 WAPE 传播、直达路径、PM 粗糙海面、Kirchhoff 反射、`H_direct_f + H_reflect_f = H_f` 的信道输出结构。
- 近海面气泡层接入：气泡层作为频率相关、深度相关的等效介质进入 PE 相位屏，最终表现为有效声速 `c_eff` 和物理衰减 `alpha_bub`。
- Level0 经验气泡层：实现了指数深度衰减的经验衰减项和可选声速扰动，用于验证符号约定、衰减方向、直达/反射路径一致性。
- Hall1D 平均气泡层：实现了 Hall-Novarini-type 的一维平均气泡数密度谱、空化率 `beta` 计算和复等效介质转换。
- 风速解耦与 Hall1D 校准：将控制 PM 粗糙海面的 `sea_wind_speed` 与控制 Hall1D 气泡数密度的 `bubble_wind_speed` 解耦，避免风速扫描时同时改变海面粗糙度和气泡浓度。
- 信道频响对比：已有 no bubble、Level0、Hall1D default、Hall1D strong 等场景的频率响应、相位、相对传输损失变化和直达/反射分量对比。
- 通信验证：已有独立通信比较脚本，当前正式使用 `scalar_h_total + pilot-LS` 模式比较 BER/SER；wideband `h_bb` 路径已诊断为 tap/delay mismatch，目前仅作为 diagnostic-only，不作为正式 BER/SER 结论。

当前尚未完成的内容包括：非均匀 plume 气泡云、完整 wideband multipath 接收机、OFDM 或训练序列宽带 LS 信道估计、高分辨率仿真以及论文级最终图表。

## 2. 研究目标与当前模型框架

当前研究目标是在垂直 PE 水声信道中引入风生气泡层，分析其对上行垂直传播信道频率响应、直达/海面反射分量以及通信 BER/SER 的影响。坐标约定为海面 `z=0`，正 `z` 向下，发射端位于较深位置，接收端位于近海面。

原始信道由两部分组成：

- 直达路径：声场从 `z_tx` 向上传播到 `z_rx`，得到 `H_direct_f`。
- 海面反射路径：声场先传播到海面，经过 PM 粗糙海面和 Kirchhoff 反射模块，再传播到接收深度，得到 `H_reflect_f`。

总信道定义为：

```text
H_f(f) = H_direct_f(f) + H_reflect_f(f)
h_total = H_f(idx_f_ref)
```

气泡层进入传播方程的方式是改变局部等效介质。原背景声速为 `c_bg(z)`，加入气泡后变为 `c_eff(x,y,z,f)`；原数值边界吸收为 `alpha_sponge(x,y)`，加入气泡后总衰减为：

```text
alpha_total(x,y,z,f) = alpha_sponge(x,y) + alpha_bub(x,y,z,f)
```

这里 `alpha_bub` 是物理气泡衰减，单位为 Np/m。它与数值 sponge 吸收相加后进入 PE 相位屏。正的 `alpha_bub` 必须导致幅度衰减，这是 Hall1D 复介质符号约定验证的关键。

## 3. Level0 经验气泡层结果

Level0 模型的主要作用不是作为最终物理模型，而是验证气泡层接入 PE 相位屏后的基本方向是否正确，包括：正衰减是否导致幅度下降、声速扰动是否主要引起相位变化、直达路径和海面反射路径是否都受到气泡层影响。

Level0 使用的基本形式为：

```text
alpha_bub(z,f) = alpha0 * exp(-z/Lb) * (f/f_ref)^p_alpha
delta_c_bub(z,f) = delta_c0 * exp(-z/Lc) * (f/f_ref)^p_c
c_eff(z,f) = c_bg(z) + delta_c_bub(z,f)
```

在当前 `compare_bubble_models_vertical_result.mat` 中，Level0 场景使用 `bubble_alpha0_np_per_m=0.02`、`bubble_layer_decay_m=20`、`bubble_delta_c0_mps=0`。因此该场景主要体现近海面指数衰减层的影响。

| 场景 | abs(h_total) | 相位 rad | max Delta TL dB | max alpha_bub Np/m | invariant error |
|---|---:|---:|---:|---:|---:|
| no_bubble | 0.056546 | -1.3376 | 0 | 0 | 0 |
| level0_empirical | 0.038889 | -1.3523 | 3.3295 | 0.017214 | 0 |

Level0 的 `abs(h_total)` 相对 no bubble 约为 `0.6878`，说明加入正衰减后接收处信道幅度下降，方向与预期一致。该结果来自 reduced-grid 频响对比，不是最终生产尺度结果。

**图 1：频率响应幅度对比**

这张图比较 no bubble、Level0、Hall1D default 和 Hall1D strong 在 4-8 kHz 频带内的 `|H(f)|`。Level0 曲线相对 no bubble 整体降低，说明经验气泡衰减层对信道幅度有可见影响。

![频率响应幅度对比](../compare_bubble_models_vertical_Figure41_H_magnitude.png)

主要观察：Level0 的幅度衰减清晰可见；Hall1D default 与 no bubble 接近；Hall1D strong 显示更明显变化。

**图 2：频率响应相位对比**

这张图比较各场景的展开相位。Level0 当前主要设置为纯衰减，因此相位变化较小；Hall1D strong 因为等效声速变化和频散影响，相位变化更明显。

![频率响应相位对比](../compare_bubble_models_vertical_Figure42_H_phase.png)

主要观察：相位变化不仅由衰减决定，也受到 `c_eff` 频率相关变化和直达/反射相干叠加影响。

**图 3：相对 no bubble 的 Delta TL**

这里的 `Delta TL` 定义为：

```text
Delta TL(f) = -20 log10( |H_case(f)| / max(|H_no_bubble(f)|, eps) )
```

正值表示相对 no bubble 出现额外传输损失。

![Delta TL 对比](../compare_bubble_models_vertical_Figure43_delta_TL.png)

主要观察：Level0 出现正的额外传输损失；某些频点上 Delta TL 不完全单调，原因是总信道由直达和反射分量相干叠加。

## 4. Hall1D 平均气泡层结果

Hall1D 是 Hall-Novarini-type 气泡谱的一维平均层简化实现。模型先根据半径 `a`、深度 `z` 和气泡风速 `U10` 计算气泡数密度，再计算空化率 `beta`，最后通过复等效介质公式得到 `c_eff` 和 `alpha_bub`。

Hall-type 数密度形式为：

```text
N(a,z,U10) = p0 * D(z,U10) * G(a,z) * (U10/13)^3 * bubble_strength_scale
```

其中 `D(z,U10)` 描述近海面随深度衰减，`G(a,z)` 描述半径谱。随后计算：

```text
beta(z) = integral (4*pi/3) * a^3 * N(a,z,U10) da
```

复等效介质转换后得到 `c_eff` 和 `alpha_bub`。当前实现已经修正符号约定，使正常 Hall1D 测试中 `alpha_bub` 为有限非负值，不依赖负衰减裁剪作为正常路径。

| 场景 | abs(h_total) | 相位 rad | max Delta TL dB | max alpha_bub Np/m | max beta | invariant error |
|---|---:|---:|---:|---:|---:|---:|
| hall1d_default | 0.056361 | -1.3329 | 0.10545 | 6.6675e-08 | 4.1242e-11 | 0 |
| hall1d_strong | 0.039869 | 1.1151 | 4.1475 | 0.054023 | 4.1242e-05 | 0 |

当前结果表明，在 `U10=5 m/s`、`z_rx=3 m`、4-8 kHz 条件下，Hall1D default 气泡浓度较低，`max alpha_bub` 仅为 `6.6675e-08 Np/m`，因此信道影响很弱。将 `bubble_strength_scale` 放大到 `1e6` 后，气泡衰减和相位影响显著增强。

**图 4：气泡元数据汇总**

这张图汇总不同气泡模型下的 `alpha_bub`、`c_eff` 或 `beta` 等诊断信息，用于判断模型是否确实产生物理量变化。

![气泡元数据汇总](../compare_bubble_models_vertical_Figure45_bubble_meta.png)

主要观察：Hall1D default 的气泡衰减量级很小，而 strong 场景中 `alpha_bub` 和 `beta` 明显增大。这支持“默认 Hall1D 在当前场景下较弱”的结论。

## 5. 参数敏感性与风速解耦校准

### 5.1 Hall1D 敏感性扫描

`sweep_hall1d_sensitivity_vertical_result.mat` 对 `bubble_strength_scale`、`sea_wind_speed`、`bubble_delta_const`、半径网格上限和接收深度进行了 reduced-grid 扫描。需要注意：该脚本中的风速扫描仍然是耦合扫描，即同一个风速同时影响 PM 粗糙海面和 Hall1D 气泡密度，因此该风速结果不能简单解释为纯气泡浓度变化。

**图 5：bubble_strength_scale 对 max Delta TL 的影响**

这张图比较不同 `bubble_strength_scale` 下最大额外传输损失。

![strength scale 敏感性](../sweep_hall1d_sensitivity_vertical_strength_scale.png)

主要观察：增大 `bubble_strength_scale` 会显著增强 `alpha_bub` 和相位影响，但由于直达/反射干涉，`max Delta TL` 不严格单调。

| bubble_strength_scale | max Delta TL dB | mean Delta TL dB | max phase diff rad | max alpha_bub Np/m | max beta |
|---:|---:|---:|---:|---:|---:|
| 1 | 0.10545 | 0.0027207 | 0.01146 | 6.6675e-08 | 4.1242e-11 |
| 1e2 | 4.7301 | -0.052242 | 0.89569 | 6.6673e-06 | 4.1242e-09 |
| 1e4 | 4.0803 | 0.42012 | 0.52066 | 0.00066501 | 4.1242e-07 |
| 1e6 | 4.1475 | 0.57646 | 3.3634 | 0.054023 | 4.1242e-05 |

**图 6：耦合风速扫描**

这张图比较不同 `sea_wind_speed` 下的最大额外传输损失。

![wind speed 敏感性](../sweep_hall1d_sensitivity_vertical_wind_speed.png)

主要观察：风速升高通常会增加气泡浓度，但该图同时包含 PM 海面粗糙度变化，因此趋势不是纯气泡效应。这个问题直接推动了后续 `bubble_wind_speed` 与 `sea_wind_speed` 的解耦。

**图 7：半径网格上限诊断**

这张图比较气泡半径网格上限从 1 mm 扩展到 2 mm 和 3 mm 后的结果变化。

![radius grid max 敏感性](../sweep_hall1d_sensitivity_vertical_radius_grid_max.png)

主要观察：在当前 4-8 kHz、近海面 reduced-grid 场景下，半径上限从 1 mm 扩展到 2-3 mm 对 `max Delta TL` 和 `max alpha_bub` 的影响很小。当前 1 mm 默认上限对该频段的 reduced-grid 诊断基本足够，但更低频或不同风速条件仍需重新检查。

**图 8：接收深度诊断**

这张图比较 `z_rx=1, 3, 5 m` 时的气泡影响。

![receiver depth 敏感性](../sweep_hall1d_sensitivity_vertical_receiver_depth.png)

主要观察：`z_rx=1 m` 时气泡影响明显增强，`max Delta TL=10.504 dB`；`z_rx=3 m` 和 `z_rx=5 m` 的影响较弱。这符合近海面气泡层随深度快速衰减的物理直觉。

### 5.2 风速解耦与校准

为避免风速扫描同时改变 PM 粗糙海面和气泡浓度，当前已加入 `bubble_wind_speed`。当 `bubble_wind_speed` 为空时，Hall1D 沿用 `sea_wind_speed`，保持向后兼容；当给定 `bubble_wind_speed` 时，PM 海面仍使用 `sea_wind_speed`，Hall1D 气泡数密度使用 `bubble_wind_speed`。

`calibrate_hall1d_bubble_vertical_result.mat` 中固定 `sea_wind_speed=5 m/s`，分别扫描 `bubble_wind_speed` 和 `bubble_strength_scale`。

**图 9：固定 PM 风速下 bubble_wind_speed 对 Delta TL 的影响**

![bubble wind calibration](../calibrate_hall1d_bubble_vertical_max_delta_TL_vs_bubble_wind_speed.png)

主要观察：在 PM 海面风速固定为 5 m/s 时，增大 `bubble_wind_speed` 会增强气泡影响。该图比耦合风速扫描更适合解释“气泡风速 forcing”本身。

| bubble_wind_speed m/s | max Delta TL dB | mean Delta TL dB | max phase diff rad | max alpha_bub Np/m | max beta |
|---:|---:|---:|---:|---:|---:|
| 3 | 0.022681 | 0.0005896 | 0.0024441 | 1.4402e-08 | 8.9083e-12 |
| 5 | 0.10545 | 0.0027207 | 0.01146 | 6.6675e-08 | 4.1242e-11 |
| 8 | 0.49784 | 0.012371 | 0.055913 | 7.0097e-07 | 4.3359e-10 |
| 12 | 3.3519 | 0.054623 | 0.39331 | 6.3338e-05 | 3.9187e-08 |
| 15 | 5.5032 | -0.0070858 | 0.65144 | 0.00030202 | 1.8704e-07 |

**图 10：bubble_strength_scale 校准**

![strength calibration](../calibrate_hall1d_bubble_vertical_max_delta_TL_vs_strength_scale.png)

主要观察：在 `sea_wind_speed=5 m/s`、`bubble_wind_speed=8 m/s` 下，`bubble_strength_scale=1e2` 对应 `max Delta TL=3.0097 dB`，可作为当前代表性校准场景之一。

| bubble_strength_scale | max Delta TL dB | mean Delta TL dB | max phase diff rad | max alpha_bub Np/m | max beta |
|---:|---:|---:|---:|---:|---:|
| 1 | 0.49784 | 0.012371 | 0.055913 | 7.0097e-07 | 4.3359e-10 |
| 1e2 | 3.0097 | 0.016118 | 0.63514 | 7.0078e-05 | 4.3359e-08 |
| 1e3 | 4.0955 | 0.42156 | 0.51778 | 0.00069905 | 4.3359e-07 |
| 1e4 | 4.0875 | 0.43852 | 0.77591 | 0.0068245 | 4.3359e-06 |
| 1e5 | 4.1536 | 0.59221 | 3.9231 | 0.056301 | 4.3359e-05 |
| 1e6 | 4.4827 | 1.3745 | 17.108 | 0.27492 | 0.00043359 |

**图 11：目标 Delta TL 匹配**

![target matching](../calibrate_hall1d_bubble_vertical_target_matching.png)

| 目标 max Delta TL dB | 匹配 scale | bubble_wind_speed m/s | 实际 max Delta TL dB | 绝对误差 dB |
|---:|---:|---:|---:|---:|
| 0.5 | 1 | 8 | 0.49784 | 0.0021573 |
| 1.0 | 1 | 8 | 0.49784 | 0.50216 |
| 3.0 | 1e2 | 8 | 3.0097 | 0.0097334 |

该表说明当前参数网格中，`bubble_strength_scale=1e2`、`bubble_wind_speed=8 m/s` 可近似匹配 3 dB 最大额外传输损失。但该校准仍是数值诊断校准，后续需要和文献或实验测量进行物理标定。

## 6. 气泡模型对信道频响的影响

当前频响对比脚本比较了四类场景：no bubble、Level0 empirical、Hall1D default 和 Hall1D strong。其目的不是通信 BER/SER，而是直接观察气泡模型对 `H_f`、直达分量和反射分量的影响。

**图 12：总信道幅度响应**

![总信道幅度响应](../compare_bubble_models_vertical_Figure41_H_magnitude.png)

图中可以看到 Level0 和 Hall1D strong 相比 no bubble 有更明显的幅度变化，Hall1D default 基本接近 no bubble。

**图 13：总信道相位响应**

![总信道相位响应](../compare_bubble_models_vertical_Figure42_H_phase.png)

图中相位变化说明气泡层不仅可能带来幅度衰减，也可能通过等效声速和频散影响相位。Hall1D strong 的相位变化较明显。

**图 14：相对 no bubble 的 Delta TL**

![相对传输损失变化](../compare_bubble_models_vertical_Figure43_delta_TL.png)

Delta TL 不是所有频点都严格单调，主要原因是 `H_f` 是直达波和海面反射波的相干叠加。某一模型增强物理衰减后，局部频点的相长/相消干涉关系也会变化。

**图 15：直达与反射分量对比**

![直达与反射分量](../compare_bubble_models_vertical_Figure44_components.png)

这张图用于观察气泡层是否同时影响直达路径和反射路径。当前传播实现中，直达主循环和反射路径传播段都使用 bubble-aware 相位屏，因此反射分量同样会受到近海面气泡层影响。

## 7. 通信验证当前结果

当前通信比较脚本的正式结果采用 `scalar_h_total + pilot_ls`。该模式把每个场景在参考频点的总信道 `h_total` 作为标量窄带信道，发送端插入固定 pilot，接收端用 LS 估计标量信道，再对数据符号均衡并统计 BER/SER。

基本流程为：

```text
tx_frame = [pilot; data]
rx_clean = tx_frame * h_true
h_hat = sum(conj(x_pilot) .* y_pilot) / sum(abs(x_pilot).^2)
rx_eq_data = y_data / h_hat
```

BER/SER 只在 data symbols 上统计，不计入 pilot。该模式比 perfect CSI 更接近实际接收机，但仍然是标量窄带等效信道，不是 full wideband multipath receiver。

wideband `h_bb` 路径已经通过 `diagnose_comm_chain_vertical_result.mat` 诊断：ideal flat channel 和 scalar PE channel 的 BER 随 Eb/N0 下降，但原 wideband `h_bb` 路径在当前实现下 BER 接近随机判决。诊断结果显示 `tap_count_original=2000`、`conv_same_delay_est=1000`，存在 tap/delay mismatch。因此 wideband `h_bb` 当前只作为 diagnostic-only，不作为正式 BER/SER 结论。

**图 16：BER vs Eb/N0**

![BER 曲线](../comm_compare_bubble_models_vertical_BER_vs_EbN0.png)

主要观察：在 scalar `h_total` + pilot-LS 模式下，no bubble 的 BER 随 Eb/N0 增大而下降。Level0 和校准 Hall1D 因为信道幅度降低，在相同发射端 Eb/N0 定义下 BER/SER 更差。

**图 17：SER vs Eb/N0**

![SER 曲线](../comm_compare_bubble_models_vertical_SER_vs_EbN0.png)

主要观察：SER 趋势与 BER 一致。当前结果反映的是相同发射功率下的通道衰减影响，而不是等接收 SNR 条件下的调制抗噪性能。

**图 18：pilot-LS 信道估计相对误差**

![信道估计相对误差](../comm_compare_bubble_models_vertical_channel_est_error_vs_EbN0.png)

主要观察：随着 Eb/N0 增大，`h_hat` 相对误差下降，说明 pilot-LS 标量估计行为符合预期。弱信道场景下相对误差更大，因为接收 pilot 幅度更低。

| 场景 | abs(h_true) | BER 0 dB | BER 10 dB | BER 20 dB | SER 0 dB | SER 10 dB | SER 20 dB | h_hat rel err 0 dB | h_hat rel err 10 dB | h_hat rel err 20 dB | effective SNR 0/10/20 dB | invariant error |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|
| no_bubble | 0.056546 | 0.4785 | 0.4455 | 0.2815 | 0.7180 | 0.6475 | 0.3900 | 1.0823 | 0.50676 | 0.17804 | -21.906 / -11.912 / -2.0451 | 8.0553e-18 |
| level0_empirical | 0.038889 | 0.48225 | 0.4605 | 0.36575 | 0.7255 | 0.6810 | 0.5175 | 1.5737 | 0.73685 | 0.25887 | -25.157 / -15.163 / -5.2966 | 7.7001e-18 |
| hall1d_default | 0.056361 | 0.4785 | 0.4450 | 0.2820 | 0.7185 | 0.6470 | 0.3910 | 1.0858 | 0.50842 | 0.17862 | -21.934 / -11.940 / -2.0736 | 7.3598e-18 |
| hall1d_calibrated_3dB | 0.039987 | 0.48275 | 0.4605 | 0.35975 | 0.7265 | 0.6795 | 0.5085 | 1.5305 | 0.71662 | 0.25176 | -24.915 / -14.921 / -5.0548 | 7.4534e-18 |
| hall1d_high_bubble_wind | 0.050570 | 0.47675 | 0.44675 | 0.30975 | 0.7170 | 0.6560 | 0.4310 | 1.2102 | 0.56664 | 0.19908 | -22.876 / -12.882 / -3.0152 | 8.1827e-18 |

需要强调：上述通信结果是 reduced-grid、标量信道、发射端 Eb/N0 参考下的阶段性验证。若后续要比较“相同接收 SNR 下不同气泡模型的调制性能”，需要增加 receiver-SNR reference mode。

## 8. 当前主要结论

- 气泡层框架已经接入垂直 PE 主传播链路，并保持 no bubble 和 off 模式的基线输出不变。
- Level0 经验层验证了正衰减进入 PE 相位屏后会导致信道幅度下降，衰减方向和符号约定正确。
- Hall1D default 在当前 `U10=5 m/s`、`z_rx=3 m`、4-8 kHz reduced-grid 场景下影响很弱，主要原因是近海面气泡浓度低且随深度快速衰减。
- 将 `sea_wind_speed` 和 `bubble_wind_speed` 解耦后，可以更清楚地观察气泡 forcing 对信道的影响，而不混入 PM 海面粗糙度变化。
- 校准 Hall1D 可以产生可见信道变化；当前代表性校准场景 `sea_wind_speed=5 m/s`、`bubble_wind_speed=8 m/s`、`bubble_strength_scale=1e2` 对应 `max Delta TL=3.0097 dB`。
- scalar `h_total + pilot-LS` 通信比较能够反映气泡导致的信道幅度降低和 BER/SER 恶化趋势。
- full wideband multipath receiver 仍是后续工作；当前 wideband `h_bb` BER/SER 只作为诊断结果，不作为正式通信结论。

## 9. 当前不足与下一步工作

当前不足：

- 非均匀 plume 气泡云尚未实现。
- Hall1D 参数仍需与文献或实测气泡衰减数据进行物理校准。
- 当前正式通信验证是 scalar `h_total + pilot-LS`，不是 full wideband multipath LS/OFDM 接收机。
- wideband `h_bb` 路径存在 tap/delay mismatch，当前仅用于诊断。
- 当前图表和数值主要是 reduced-grid 结果，需要后续高分辨率生产仿真验证。
- 当前通信比较基于相同发射端 Eb/N0；如果要分离“传播损耗”和“调制抗噪性”，需要增加 receiver-SNR reference mode。

建议下一步：

- 优先基于文献或实验数据对 Hall1D 的 `bubble_strength_scale`、`bubble_wind_speed` 和阻尼参数进行校准。
- 增加 receiver-SNR reference mode，区分等发射功率和等接收 SNR 两类通信比较。
- 在 scalar pilot-LS 结果稳定后，再设计 full wideband receiver，例如 OFDM 或训练序列 LS。
- 在平均气泡层结论清楚后，再引入非均匀 plume 气泡云。
- 选择代表性场景进行高分辨率生产仿真，并生成论文级频响、Delta TL 和通信性能图。

## 10. 本报告使用的结果文件与图像

使用的结果文件：

- `compare_bubble_models_vertical_result.mat`
- `sweep_hall1d_sensitivity_vertical_result.mat`
- `calibrate_hall1d_bubble_vertical_result.mat`
- `diagnose_comm_chain_vertical_result.mat`
- `comm_compare_bubble_models_vertical_result.mat`

引用的图像：

- `../compare_bubble_models_vertical_Figure41_H_magnitude.png`
- `../compare_bubble_models_vertical_Figure42_H_phase.png`
- `../compare_bubble_models_vertical_Figure43_delta_TL.png`
- `../compare_bubble_models_vertical_Figure44_components.png`
- `../compare_bubble_models_vertical_Figure45_bubble_meta.png`
- `../sweep_hall1d_sensitivity_vertical_strength_scale.png`
- `../sweep_hall1d_sensitivity_vertical_wind_speed.png`
- `../sweep_hall1d_sensitivity_vertical_radius_grid_max.png`
- `../sweep_hall1d_sensitivity_vertical_receiver_depth.png`
- `../calibrate_hall1d_bubble_vertical_max_delta_TL_vs_bubble_wind_speed.png`
- `../calibrate_hall1d_bubble_vertical_max_delta_TL_vs_strength_scale.png`
- `../calibrate_hall1d_bubble_vertical_target_matching.png`
- `../comm_compare_bubble_models_vertical_BER_vs_EbN0.png`
- `../comm_compare_bubble_models_vertical_SER_vs_EbN0.png`
- `../comm_compare_bubble_models_vertical_channel_est_error_vs_EbN0.png`

待补充验证：

- 高分辨率生产网格下的最终频响和通信结果。
- full wideband multipath receiver 的 BER/SER。
- plume 非均匀气泡云模型结果。
- 与文献或实验数据对齐后的 Hall1D 物理标定结果。
