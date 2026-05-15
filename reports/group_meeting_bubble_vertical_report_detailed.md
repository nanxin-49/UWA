# 垂直 PE 水声信道中风生气泡层建模与通信验证详细阶段报告

## 1. 当前工作目标

本项目的基础程序是一个垂直方向水声信道仿真框架，主传播模型为垂直 WAPE/PE 型声场推进。当前阶段的工作目标是在已有垂直 PE 信道中加入风生近海面气泡层，使气泡通过传播介质参数进入声场计算，而不是作为接收端噪声或最终信道的经验标量损耗。

当前气泡层建模目标可以概括为：

```text
风速、深度、频率、气泡半径谱
-> 气泡数密度 n(a,z,U10)
-> 等效介质参数 c_eff(z,f), alpha_bub(z,f)
-> PE/WAPE 相位屏
-> H_direct_f, H_reflect_f, H_f
-> 通信链路 BER/SER 对比
```

本报告只总结当前已经实现并验证的内容。非均匀 plume 气泡云、full wideband multipath 接收机、OFDM 或宽带 LS 信道估计均未实现。

## 2. 代码文件结构与核心作用

下表按当前项目角色对主要 MATLAB 文件和文档进行分类。这里的“核心文件”指会直接影响主传播物理模型、公开 API 或当前正式通信比较结果的文件；诊断脚本虽然重要，但不属于主传播模型本体。

| 文件名 | 类型 | 是否核心文件 | 主要作用 | 主要输入 | 主要输出 | 与其他文件的调用关系 |
|---|---|---:|---|---|---|---|
| `CARPE3D_vertical.m` | 核心 API / 配置边界 | 是 | 公开信道 API；合并 `paramsV` 默认值；执行参数校验；调用传播核心；打包输出结构 | `paramsV` 用户参数结构体 | `output`，含 `H_f`、`h_total`、`bubble_meta`、`config` 等 | 被 `main_vertical.m`、`compare_bubble_models_vertical.m`、`comm_compare_bubble_models_vertical.m` 等调用；内部调用 `propWAPE_vertical(cfg)` |
| `propWAPE_vertical.m` | 核心传播 | 是 | 真正的垂直 WAPE/PE 传播核心；构造网格、频率轴、直达路径、反射路径和最终频域信道 | 已校验的 `cfg` | `h_direct`、`h_reflect`、`h_total`、`H_direct_f`、`H_reflect_f`、`H_f`、`f_axis`、`bubble_meta` | 被 `CARPE3D_vertical.m` 调用；内部调用 `bubble_environment_vertical.m` 和 `pm_surface_kirchhoff_module.m` |
| `pm_surface_kirchhoff_module.m` | 核心传播 / 粗糙海面 | 是 | 生成 PM 粗糙海面并施加 Kirchhoff 相位反射 | 入射海面场、横向波数、海面参数、风速等 | 反射后的场、海面起伏、相位扰动、粗糙面元数据 | 被 `propWAPE_vertical.m` 的反射路径调用 |
| `bubble_environment_vertical.m` | 气泡模型统一接口 | 是 | PE 相位屏调用的统一气泡环境函数；根据模型返回 `c_eff_xy` 与 `alpha_bub_xy` | `x,y,z_curr,f_hz,c_bg,cfg` | `c_eff_xy`、`alpha_bub_xy`、单步 `meta` | 被 `propWAPE_vertical.m/local_phase_screen` 调用；Hall1D 时调用 `bubble_hall_spectrum.m` 和 `bubble_effective_medium.m` |
| `bubble_hall_spectrum.m` | 气泡物理模型 | 是 | 计算 Hall 型一维平均气泡半径谱、空化率 `beta` 和谱元数据 | 半径网格、深度、气泡风速、`cfg` | `n_a`、`beta`、`spec_meta` | 被 `bubble_environment_vertical.m` 在 `bubble_model='hall1d'` 时调用 |
| `bubble_effective_medium.m` | 气泡等效介质 | 是 | 将气泡谱转换为等效声速 `c_eff` 和气泡衰减 `alpha_bub` | 半径网格、`n_a`、深度、频率、背景声速、`cfg` | `c_eff`、`alpha_bub`、`em_meta` | 被 `bubble_environment_vertical.m` 在 `bubble_model='hall1d'` 时调用 |
| `comm_compare_bubble_models_vertical.m` | 当前通信比较脚本 | 是 | 当前主气泡通信比较脚本；默认使用 `scalar_h_total + pilot_ls`；保留 `wideband_diagnostic` 路径 | 场景参数、信道输出、QPSK 数据和 pilot | 各场景 BER/SER、`h_hat` 误差、有效 SNR、图和 `.mat` 结果 | 调用 `CARPE3D_vertical.m`、`modem_psk.m`、`noise_inject_vertical.m` |
| `comm_main_vertical_psk.m` | 原始通信 demo | 否，参考入口 | 原始端到端 MPSK demo；曾作为 `H_f` 消费者和通信链路参考 | 内部脚本参数 | BER/SER、通信结果 | 调用 `CARPE3D_vertical.m`、`modem_psk.m`、`noise_inject_vertical.m`；当前不作为气泡通信主比较脚本 |
| `modem_psk.m` | 通信基础模块 | 是，通信基础 | MPSK/QPSK 调制、解调和误码率统计 | bit 序列或接收符号、调制阶数 `M` | PSK 符号、判决 bit、BER/SER | 被 `comm_compare_bubble_models_vertical.m`、`comm_main_vertical_psk.m` 和诊断脚本调用 |
| `noise_inject_vertical.m` | 通信噪声模块 | 是，通信基础 | AWGN 或自定义噪声注入；不包含气泡效应 | `rx_clean`、噪声配置、`signal_ref` | `rx_noisy`、噪声、噪声元数据 | 被通信脚本调用；气泡不能在这里实现 |
| `compare_bubble_models_vertical.m` | 诊断脚本 | 否 | reduced-grid 信道频响对比：no bubble、Level 0、Hall1D default、Hall1D strong | 内置 reduced-grid 参数 | `compare_bubble_models_vertical_result.mat` 和 PNG | 调用 `CARPE3D_vertical.m` |
| `sweep_hall1d_sensitivity_vertical.m` | 诊断脚本 | 否 | Hall1D 参数敏感性扫描，包括强度、风速、阻尼、半径上限、接收深度 | 内置 reduced-grid 扫描参数 | `sweep_hall1d_sensitivity_vertical_result.mat` 和 PNG | 调用 `CARPE3D_vertical.m` |
| `calibrate_hall1d_bubble_vertical.m` | 诊断 / 校准脚本 | 否 | 解耦 `bubble_wind_speed` 后进行 Hall1D 校准和目标 ΔTL 匹配 | 内置 reduced-grid 校准参数 | `calibrate_hall1d_bubble_vertical_result.mat`、校准表和 PNG | 调用 `CARPE3D_vertical.m` |
| `diagnose_comm_chain_vertical.m` | 诊断脚本 | 否 | 诊断 modem/noise、标量 PE 信道和 wideband `h_bb` 路径问题 | 内置 reduced-grid 通信参数 | `diagnose_comm_chain_vertical_result.mat` | 调用 `CARPE3D_vertical.m`、`modem_psk.m`、`noise_inject_vertical.m` |
| `main_vertical.m` | 原始信道演示脚本 | 否，参考入口 | 单次垂直信道仿真入口，使用较大默认网格并保存 `vertical_upward_4k_uniform.mat` 与 Figure11-15 | 内置 `paramsV` | `simulata_vertical`、`paramsV`、信道图 | 调用 `CARPE3D_vertical.m` |
| `AGENTS.md` | 工程文档 | 否 | 项目协作约束、数值一致性要求和文件修改规则 | 无 | 文档约束 | 指导后续开发 |
| `BUBBLE_EXTENSION_SPEC.md` | 工程 / 模型规格文档 | 否 | 气泡扩展公式、接口和分阶段实现计划 | 无 | 文档规格 | 指导气泡模型实现 |
| `PROJECT_CONTEXT.md` | 工程文档 | 否 | 记录项目代码结构、执行路径、数据流和扩展接口 | 无 | 文档说明 | 作为代码维护和交接参考 |
| `vertical_comm_guide.md` | 工程文档 / 历史说明 | 否 | 原始通信模块说明和公式补充；当前存在编码显示问题，需谨慎参考 | 无 | 文档说明 | 可作为历史参考，不作为当前实现的唯一依据 |

需要特别说明：

- `CARPE3D_vertical.m` 是公开 API 和配置校验边界，不直接承担 PE 数值推进。
- `propWAPE_vertical.m` 是真正的 PE/WAPE 传播核心。
- `bubble_environment_vertical.m` 是传播核心调用的气泡统一接口，不是独立后处理。
- `bubble_hall_spectrum.m` 和 `bubble_effective_medium.m` 是 Hall1D 物理模型的两个 helper。
- `comm_compare_bubble_models_vertical.m` 是当前气泡通信比较的主脚本。
- `comm_main_vertical_psk.m` 是原始通信 demo，目前主要作为参考，不再作为气泡通信比较主入口。
- `main_vertical.m` 是当前仓库中实际存在的原始信道演示入口；若只想快速看气泡模型对比，优先运行 reduced-grid 诊断脚本而不是该大网格入口。

## 3. 主程序数据流

当前主信道数据流从用户参数 `paramsV` 开始，经 `CARPE3D_vertical` 校验为运行时配置 `cfg`，再进入 `propWAPE_vertical` 执行频率循环和 PE 传播。

文字流程如下：

1. 用户或脚本构造 `paramsV`。
2. `CARPE3D_vertical(paramsV)` 调用 `local_prepare_config(paramsV)`。
3. `local_prepare_config` 合并默认值、校验网格、声源、接收器、气泡参数、反射参数和通信相关输出控制。
4. 得到 `cfg` 后，`CARPE3D_vertical` 调用 `propWAPE_vertical(cfg)`。
5. `propWAPE_vertical` 构造横向网格、波数网格、频率轴 `f_axis`。
6. 对每个频率执行直达路径 WAPE 推进，得到 `H_direct_f(ifq)`。
7. 若启用海面反射，先从发射深度推进到海面，调用 `pm_surface_kirchhoff_module` 施加粗糙海面反射，再从海面推进到接收深度，得到 `H_reflect_f(ifq)`。
8. 对每个频率计算 `H_f(ifq) = H_direct_f(ifq) + H_reflect_f(ifq)`。
9. 在参考频点 `idx_f_ref` 处提取 `h_direct`、`h_reflect` 和 `h_total`。
10. `CARPE3D_vertical` 将频域和标量结果打包到 `output`。
11. 通信脚本读取 `channel.H_f` 或 `channel.h_total` 计算 BER/SER。

ASCII 流程图：

```text
paramsV
  |
  v
CARPE3D_vertical(paramsV)
  |
  +--> local_prepare_config(paramsV)
  |       |
  |       v
  |      cfg
  |
  v
propWAPE_vertical(cfg)
  |
  +--> build x-y grid, kx/ky/kappa2, sponge
  |
  +--> build f_axis and idx_f_ref
  |
  +--> for each f in f_axis
  |       |
  |       +--> direct WAPE march: z_tx -> z_rx
  |       |       |
  |       |       +--> local_phase_screen
  |       |               |
  |       |               +--> bubble_environment_vertical
  |       |
  |       +--> optional reflection path
  |               |
  |               +--> local_march_field: z_tx -> 0
  |               +--> pm_surface_kirchhoff_module
  |               +--> local_march_field: 0 -> z_rx
  |
  +--> H_direct_f, H_reflect_f
  |
  +--> H_f = H_direct_f + H_reflect_f
  |
  +--> h_direct = H_direct_f(idx_f_ref)
  +--> h_reflect = H_reflect_f(idx_f_ref)
  +--> h_total = H_f(idx_f_ref)
  |
  v
output/channel struct
  |
  +--> diagnostic scripts use H_f and bubble_meta
  +--> communication scripts use h_total or diagnostic H_f -> h_bb
```

三类信道量的区别：

| 名称 | 含义 | 当前计算方式 |
|---|---|---|
| direct path | 发射机到接收机的直达上行传播 | 从 `z_tx` 直接推进到 `z_rx` |
| reflected path | 发射机到海面，再经粗糙海面反射到接收机 | `z_tx -> 0`，Kirchhoff 反射，`0 -> z_rx` |
| total channel | 接收端总复信道 | `H_f = H_direct_f + H_reflect_f` |

## 4. 如何运行代码查看结果

本节按“想看什么结果”列出推荐运行入口。除非需要原始大网格演示，建议优先运行 reduced-grid 诊断脚本；这些脚本运行成本更低，且结果文件已经被 `.gitignore` 忽略。

| 目的 | 推荐运行文件 | MATLAB 命令 | 主要输出文件 | 是否推荐作为当前正式结果 |
|---|---|---|---|---|
| 原始单次垂直信道演示 | `main_vertical.m` | `main_vertical` | `vertical_upward_4k_uniform.mat`、`vertical_upward_4k_uniform_Figure11_xy.png` 到 `Figure15_surface_kirchhoff.png` | 参考入口；默认 `1024 x 1024`，运行成本较高 |
| 气泡模型频响对比 | `compare_bubble_models_vertical.m` | `compare_bubble_models_vertical` | `compare_bubble_models_vertical_result.mat`、`compare_bubble_models_vertical_Figure41_H_magnitude.png` 等 | 是，reduced-grid 信道对比结果 |
| Hall1D 参数敏感性 | `sweep_hall1d_sensitivity_vertical.m` | `sweep_hall1d_sensitivity_vertical` | `sweep_hall1d_sensitivity_vertical_result.mat`、strength/wind/damping/radius/depth PNG | 是，reduced-grid 敏感性结果 |
| Hall1D 风速解耦与校准 | `calibrate_hall1d_bubble_vertical.m` | `calibrate_hall1d_bubble_vertical` | `calibrate_hall1d_bubble_vertical_result.mat`、target matching 和校准 PNG | 是，reduced-grid 校准结果 |
| 通信链路问题诊断 | `diagnose_comm_chain_vertical.m` | `diagnose_comm_chain_vertical` | `diagnose_comm_chain_vertical_result.mat` | 诊断用；用于说明 wideband `h_bb` 当前不可靠 |
| 当前气泡通信对比 | `comm_compare_bubble_models_vertical.m` | `comm_compare_bubble_models_vertical` | `comm_compare_bubble_models_vertical_result.mat`、BER/SER/channel estimation PNG | 是，当前正式通信比较为 `scalar_h_total + pilot_ls` |
| 原始 MPSK 通信 demo | `comm_main_vertical_psk.m` | `comm_main_vertical_psk` | `psk_comm_result.mat` 及通信图 | 参考入口；不是当前气泡通信主比较脚本 |

如果在 Windows PowerShell 中从外部调用 MATLAB，可使用以下形式：

```powershell
& 'D:\Matlab2025b\bin\matlab.exe' -batch "cd('E:/MISC/CARPE3D_matlab/Explain'); comm_compare_bubble_models_vertical;"
```

当前不同运行入口的定位如下：

- 想看传播模型是否工作：运行 `compare_bubble_models_vertical`。
- 想看 Hall1D 哪些参数最敏感：运行 `sweep_hall1d_sensitivity_vertical`。
- 想看固定 PM 海面后气泡风速和强度如何校准：运行 `calibrate_hall1d_bubble_vertical`。
- 想解释为什么旧 wideband BER/SER 接近随机：查看或运行 `diagnose_comm_chain_vertical`。
- 想比较气泡对 QPSK BER/SER 的影响：运行 `comm_compare_bubble_models_vertical`，但需明确它当前使用的是标量 `h_total` pilot-LS，不是 full wideband receiver。

## 5. PE 相位屏与气泡接入公式

参考波数为：

```math
k_0 = \frac{2\pi f}{c_0}.
```

原始无气泡相位屏为：

```math
S_0(x,y,z,f)
=
\exp\left[
-i k_0 \Delta s
\left(
\frac{c_{\rm bg}(z)-c_0}{c_{\rm bg}(z)}
- i\frac{\alpha_{\rm sponge}(x,y)}{k_0}
\right)
\right].
```

加入气泡后的相位屏为：

```math
S_{\rm bub}(x,y,z,f)
=
\exp\left[
-i k_0 \Delta s
\left(
\frac{c_{\rm eff}(x,y,z,f)-c_0}{c_{\rm eff}(x,y,z,f)}
- i\frac{\alpha_{\rm sponge}(x,y)+\alpha_{\rm bub}(x,y,z,f)}{k_0}
\right)
\right].
```

其中：

| 符号 | 含义 |
|---|---|
| `alpha_sponge(x,y)` | 数值边界吸收，用于抑制横向边界反射 |
| `alpha_bub(x,y,z,f)` | 物理气泡衰减，单位 Np/m |
| `c_bg(z)` | 背景声速 |
| `c_eff(x,y,z,f)` | 含气泡后的等效声速 |

两类衰减在相位屏中相加：

```math
\alpha_{\rm total}(x,y,z,f)
=
\alpha_{\rm sponge}(x,y)
+
\alpha_{\rm bub}(x,y,z,f).
```

正的 `alpha_bub` 必须导致幅度衰减。由相位屏可见，正的总衰减对应近似 `exp(-alpha_total*ds)` 的幅度因子。因此 Hall1D 有效介质公式中的虚部符号非常重要；若符号反，会出现负衰减，需要裁剪。当前实现中已经通过符号验证选择 `-1i*d` 的阻尼符号，使正常路径下 `alpha_bub = omega*imag(q)` 为非负。

代码实现上，`propWAPE_vertical.m` 中的 `local_phase_screen` 负责构造相位屏。该 helper 调用：

```matlab
[c_eff_xy, alpha_bub_xy, bubble_step_meta] = ...
    bubble_environment_vertical(x, y, z_curr, f_hz, c_bg, cfg);
```

然后执行：

```matlab
alpha_total_xy = alpha_xy + alpha_bub_xy;
U_real_xy = (c_eff_xy - cfg.c0) ./ c_eff_xy;
screen = exp(-1i * k0 * ds * (U_real_xy - 1i * alpha_total_xy / k0));
```

`local_phase_screen` 同时被两个地方调用：

- 主 direct propagation loop；
- `local_march_field`。

`local_march_field` 很关键，因为 reflected path 的两段传播都依赖它：

- `z_tx -> 0`；
- `0 -> z_rx`。

因此当前气泡并不是只接入直达路径，而是同时影响直达路径和海面反射路径。

## 6. Level 0 经验气泡层

Level 0 是经验模型，主要用于验证气泡相位屏链路是否正确，而不是最终物理模型。它可以独立测试衰减效应和声速扰动效应。

经验衰减公式：

```math
\alpha_{\rm bub}(z,f)
=
\alpha_0
\exp\left(-\frac{z}{L_b}\right)
\left(\frac{f}{f_{\rm ref}}\right)^{p_\alpha}.
```

经验声速扰动公式：

```math
\Delta c_{\rm bub}(z,f)
=
\Delta c_0
\exp\left(-\frac{z}{L_c}\right)
\left(\frac{f}{f_{\rm ref}}\right)^{p_c}.
```

等效声速：

```math
c_{\rm eff}(z,f) = c_{\rm bg}(z) + \Delta c_{\rm bub}(z,f).
```

对应配置字段：

| 字段 | 含义 |
|---|---|
| `bubble_alpha0_np_per_m` | 近海面参考气泡衰减强度，单位 Np/m |
| `bubble_layer_decay_m` | 衰减层 e-folding 深度 |
| `bubble_f_ref_hz` | 频率标度参考值 |
| `bubble_alpha_freq_exp` | 衰减频率幂指数 |
| `bubble_delta_c0_mps` | 近海面声速扰动强度 |
| `bubble_sound_speed_decay_m` | 声速扰动 e-folding 深度 |
| `bubble_sound_speed_freq_exp` | 声速扰动频率幂指数 |
| `bubble_apply_attenuation` | 是否启用经验衰减 |
| `bubble_apply_sound_speed` | 是否启用经验声速扰动 |

Level 0 的作用：

- 验证正衰减是否导致 `|H|` 下降；
- 验证声速扰动是否主要体现为相位变化；
- 验证 direct path 与 reflected path 是否都接入相同气泡相位屏；
- 验证 `H_f = H_direct_f + H_reflect_f` 不变量未被破坏。

当前可用 reduced-grid 结果来自 `compare_bubble_models_vertical_result.mat`：


| 场景 | <code>&#124;h_total&#124;</code> | `phase(h_total)` rad | max ΔTL dB | max `alpha_bub` Np/m | invariant error |
|---|---:|---:|---:|---:|---:|
| no_bubble | 0.056546 | -1.3376 | 0 | 0 | 0 |
| level0_empirical | 0.038889 | -1.3523 | 3.3295 | 0.017214 | 0 |

由当前结果可计算：

```text
|h_total_level0| / |h_total_no_bubble| = 0.038889 / 0.056546 = 0.6878
```

缺失或待补充：

| 指标 | 状态 |
|---|---|
| direct-only Level 0 `|H_bub|/|H_0|` | 当前结果文件未找到 |
| reflected-only Level 0 分量比值 | 当前结果文件未找到 |
| sound-speed-only 相位变化 | 当前结果文件未找到 |

## 7. Hall1D 平均气泡层公式与实现

Hall1D 是当前实现的物理型一维平均气泡层。该模型仍然是横向均匀的，即当前未实现 plume 或水平非均匀气泡云。

Hall 型数密度公式：

```math
N(a,z,u_{10})
=
p_0
D(z,u_{10})
G(a,z)
\left(\frac{u_{10}}{13}\right)^3
\cdot
\texttt{bubble\_strength\_scale}.
```

其中：

```math
p_0 = 1.6\times 10^{10}\ {\rm m^{-4}}.
```

深度衰减：

```math
D(z,u_{10}) = \exp\left[-\frac{z}{L(u_{10})}\right].
```

```math
L(u_{10}) =
\begin{cases}
0.4, & u_{10} \le 7.5, \\
0.4 + 0.115(u_{10}-7.5), & u_{10} > 7.5.
\end{cases}
```

半径谱函数：

```math
G(a,z) =
\begin{cases}
0, & a < 10\ \mu{\rm m}, \\
\left(\frac{a_{\rm ref}}{a}\right)^4,
& 10\ \mu{\rm m} \le a \le a_{\rm ref}, \\
\left(\frac{a_{\rm ref}}{a}\right)^{\chi},
& a_{\rm ref} < a \le 1000\ \mu{\rm m}, \\
0, & a > 1000\ \mu{\rm m}.
\end{cases}
```

```math
a_{\rm ref}(z)
=
54.4 + 1.984\times 10^{-6} z
\quad
{\rm micrometers}.
```

```math
\chi(z) = 4.37 + \left(\frac{z}{2.55}\right)^2.
```

空化率：

```math
\beta(z)
=
\int
\frac{4\pi}{3}
a^3
N(a,z,u_{10})
\,da.
```

如果 `beta > bubble_beta_max`，代码会按比例缩放 `n_a`，使 `beta` 不超过上限，并在元数据中记录 warning。

共振半径：

```math
a_{\rm res}(f,z)
=
\frac{1}{2\pi f}
\sqrt{
\frac{3\gamma P_0(z)}{\rho_w}
}.
```

静水压力：

```math
P_0(z) = P_{\rm atm} + \rho_w g z.
```

当前使用的复声速形式：

```math
\frac{1}{\tilde c^2}
=
\frac{1}{c_{\rm bg}^2}
+
\frac{1}{\pi f^2}
\int
\frac{a n_a}
{\left(\frac{a_{\rm res}}{a}\right)^2 - 1 - i d}
\,da.
```

复慢度：

```math
q = \frac{1}{\tilde c}.
```

输出介质参数：

```math
c_{\rm eff} = \frac{1}{\Re(q)}.
```

```math
\alpha_{\rm bub} = \omega \Im(q),
\quad
\omega = 2\pi f.
```

实现说明：

- `bubble_hall_spectrum.m` 计算 `n_a`、`beta` 和谱元数据。
- `bubble_effective_medium.m` 计算 `a_res`、复声速、`c_eff` 和 `alpha_bub`。
- `bubble_environment_vertical.m` 将一维标量 `c_eff` 和 `alpha_bub` 扩展为 `[numel(y), numel(x)]` 的矩阵。
- 阻尼分母使用 `-1i*d`。这是经过符号约定检查后的当前实现，使正的 `alpha_bub` 与 PE 相位屏中的幅度衰减一致。
- `negative_attenuation_clipped` 只保留为安全 fallback，不应是标准 reduced Hall1D 测试的正常路径。

当前 `compare_bubble_models_vertical_result.mat` 中 Hall1D reduced-grid 结果：

| 场景 | `h_total` | `phase(h_total)` rad | max ΔTL dB | max `alpha_bub` Np/m | max beta |
|---|---:|---:|---:|---:|---:|
| hall1d_default | 0.056361 | -1.3329 | 0.10545 | 6.6675e-08 | 4.1242e-11 |
| hall1d_strong, `bubble_strength_scale=1e6` | 0.039869 | 1.1151 | 4.1475 | 0.054023 | 4.1242e-05 |

结论：默认 Hall1D 在当前 `U10=5 m/s`、`z_rx=3 m`、4 到 8 kHz reduced-grid 场景中较弱；增大 `bubble_strength_scale` 后可产生明显幅度和相位变化。

## 8. 风速解耦逻辑

当前模型中有两个风速概念：

| 参数 | 控制对象 |
|---|---|
| `sea_wind_speed` | PM 粗糙海面谱和 Kirchhoff 反射 |
| `bubble_wind_speed` | Hall1D 气泡数密度谱 |

向后兼容规则：

- 如果 `bubble_wind_speed=[]` 或未提供，Hall1D 使用 `sea_wind_speed`；
- 如果提供 `bubble_wind_speed`，Hall1D 使用 `bubble_wind_speed`；
- PM 海面粗糙度始终使用 `sea_wind_speed`。

该解耦很重要。早期 wind sweep 中，`sea_wind_speed` 同时改变海面粗糙度和气泡浓度，导致频响变化中混有粗糙面散射和气泡介质两种因素。引入 `bubble_wind_speed` 后，可以固定 PM 海面条件，只扫描气泡浓度。

`calibrate_hall1d_bubble_vertical_result.mat` 中固定：

```text
sea_wind_speed = 5 m/s
```
```math
\Delta TL_{\rm bub}(f)
=
-20\log_{10}
\left(
\frac{|H_{\rm bub}(f)|}{|H_0(f)|}
\right)
```
然后单独扫描 `bubble_wind_speed`：

| `bubble_wind_speed` m/s | max ΔTL dB | mean ΔTL dB | max phase diff rad | max `alpha_bub` | max beta |
|---:|---:|---:|---:|---:|---:|
| 3 | 0.022681 | 0.0005896 | 0.0024441 | 1.4402e-08 | 8.9083e-12 |
| 5 | 0.10545 | 0.0027207 | 0.01146 | 6.6675e-08 | 4.1242e-11 |
| 8 | 0.49784 | 0.012371 | 0.055913 | 7.0097e-07 | 4.3359e-10 |
| 12 | 3.3519 | 0.054623 | 0.39331 | 6.3338e-05 | 3.9187e-08 |
| 15 | 5.5032 | -0.0070858 | 0.65144 | 0.00030202 | 1.8704e-07 |

在 `bubble_wind_speed=8 m/s` 下扫描强度：

| `bubble_strength_scale` | max ΔTL dB | mean ΔTL dB | max phase diff rad | max `alpha_bub` | max beta |
|---:|---:|---:|---:|---:|---:|
| 1 | 0.49784 | 0.012371 | 0.055913 | 7.0097e-07 | 4.3359e-10 |
| 100 | 3.0097 | 0.016118 | 0.63514 | 7.0078e-05 | 4.3359e-08 |
| 1000 | 4.0955 | 0.42156 | 0.51778 | 0.00069905 | 4.3359e-07 |
| 10000 | 4.0875 | 0.43852 | 0.77591 | 0.0068245 | 4.3359e-06 |
| 1e5 | 4.1536 | 0.59221 | 3.9231 | 0.056301 | 4.3359e-05 |
| 1e6 | 4.4827 | 1.3745 | 17.108 | 0.27492 | 0.00043359 |

目标匹配表：

| 目标 max ΔTL dB | 匹配 `bubble_strength_scale` | `bubble_wind_speed` m/s | 匹配 max ΔTL dB | 绝对误差 dB |
|---:|---:|---:|---:|---:|
| 0.5 | 1 | 8 | 0.49784 | 0.0021573 |
| 1.0 | 1 | 8 | 0.49784 | 0.50216 |
| 3.0 | 100 | 8 | 3.0097 | 0.0097334 |

因此当前代表性校准场景为：

```text
sea_wind_speed = 5 m/s
bubble_wind_speed = 8 m/s
bubble_strength_scale = 1e2
max_delta_TL_dB = 3.0097 dB
```

## 9. 核心验证结果

本节只使用当前已有结果文件，未重新运行大规模仿真。所有数值均为 reduced-grid 或诊断结果。

### 8.1 disabled/off 回归

生成脚本：

```text
calibrate_hall1d_bubble_vertical.m
```

结果文件：

```text
calibrate_hall1d_bubble_vertical_result.mat
```

验证表：

| 检查项 | rel `H_f` | rel `H_direct_f` | rel `H_reflect_f` | invariant error | pass |
|---|---:|---:|---:|---:|---|
| disabled override unchanged | 0 | 0 | 0 | 8.0553e-18 | true |
| off mode unchanged | 0 | 0 | 0 | 8.0553e-18 | true |
| level0 ignores `bubble_wind_speed` | 0 | 0 | 0 | 7.7001e-18 | true |
| hall1d legacy equals explicit sea wind | 0 | 0 | 0 | 7.3598e-18 | true |

结论：新增 `bubble_wind_speed` 没有破坏 disabled/off、Level0 和 Hall1D 旧行为。

### 8.2 Level0 验证

生成脚本：

```text
compare_bubble_models_vertical.m
```

结果文件：

```text
compare_bubble_models_vertical_result.mat
```

| 场景 | `h_total` | max ΔTL dB | max `alpha_bub` | invariant error |
|---|---:|---:|---:|---:|
| no_bubble | 0.056546 | 0 | 0 | 0 |
| level0_empirical | 0.038889 | 3.3295 | 0.017214 | 0 |

Level0 direct-only、reflected-only 和 sound-speed-only 的细分数值：当前结果文件未找到。

### 8.3 Hall1D 验证

生成脚本：

```text
compare_bubble_models_vertical.m
```

结果文件：

```text
compare_bubble_models_vertical_result.mat
```

| 场景 | `|h_total|` | max ΔTL dB | max `alpha_bub` | max beta | invariant error |
|---|---:|---:|---:|---:|---:|
| hall1d_default | 0.056361 | 0.10545 | 6.6675e-08 | 4.1242e-11 | 0 |
| hall1d_strong | 0.039869 | 4.1475 | 0.054023 | 4.1242e-05 | 0 |

默认 Hall1D 在当前场景下很弱；强度放大后可产生可见变化。

### 8.4 敏感性扫描

生成脚本：

```text
sweep_hall1d_sensitivity_vertical.m
```

结果文件：

```text
sweep_hall1d_sensitivity_vertical_result.mat
```

核心结论：

- `bubble_strength_scale` 能显著改变 `alpha_bub` 和相位响应；
- `sea_wind_speed` 扫描在未解耦时同时改变 PM 粗糙海面和气泡浓度，结果存在混合因素；
- `bubble_delta_const` 在当前扫描范围内对 max ΔTL 的影响较小；
- 半径上限从 `1 mm` 扩展到 `2 mm` 或 `3 mm` 对当前 4 到 8 kHz reduced case 影响很小；
- `z_rx=1 m` 的气泡影响明显强于 `z_rx=3 m` 和 `z_rx=5 m`。

图：

![Hall1D strength sweep](../sweep_hall1d_sensitivity_vertical_strength_scale.png)

![Hall1D wind speed sweep](../sweep_hall1d_sensitivity_vertical_wind_speed.png)

![Hall1D damping sweep](../sweep_hall1d_sensitivity_vertical_damping_const.png)

![Hall1D radius grid max sweep](../sweep_hall1d_sensitivity_vertical_radius_grid_max.png)

![Hall1D receiver depth sweep](../sweep_hall1d_sensitivity_vertical_receiver_depth.png)

![Hall1D max alpha summary](../sweep_hall1d_sensitivity_vertical_max_alpha_summary.png)

### 8.5 解耦校准结果

生成脚本：

```text
calibrate_hall1d_bubble_vertical.m
```

结果文件：

```text
calibrate_hall1d_bubble_vertical_result.mat
```

核心匹配结果：

| 目标 max ΔTL dB | 匹配强度 | 气泡风速 m/s | 匹配 max ΔTL dB |
|---:|---:|---:|---:|
| 0.5 | 1 | 8 | 0.49784 |
| 1.0 | 1 | 8 | 0.49784 |
| 3.0 | 100 | 8 | 3.0097 |

图：

![Calibration delta TL vs bubble wind](../calibrate_hall1d_bubble_vertical_max_delta_TL_vs_bubble_wind_speed.png)

![Calibration alpha vs bubble wind](../calibrate_hall1d_bubble_vertical_max_alpha_vs_bubble_wind_speed.png)

![Calibration delta TL vs strength](../calibrate_hall1d_bubble_vertical_max_delta_TL_vs_strength_scale.png)

![Calibration alpha vs strength](../calibrate_hall1d_bubble_vertical_max_alpha_vs_strength_scale.png)

![Calibration target matching](../calibrate_hall1d_bubble_vertical_target_matching.png)

### 8.6 信道频响对比

生成脚本：

```text
compare_bubble_models_vertical.m
```

结果文件：

```text
compare_bubble_models_vertical_result.mat
```

| 场景 | `|h_total|` | `phase(h_total)` rad | max ΔTL dB | max `alpha_bub` | max beta |
|---|---:|---:|---:|---:|---:|
| no_bubble | 0.056546 | -1.3376 | 0 | 0 | NaN |
| level0_empirical | 0.038889 | -1.3523 | 3.3295 | 0.017214 | NaN |
| hall1d_default | 0.056361 | -1.3329 | 0.10545 | 6.6675e-08 | 4.1242e-11 |
| hall1d_strong | 0.039869 | 1.1151 | 4.1475 | 0.054023 | 4.1242e-05 |

图：

![Channel magnitude](../compare_bubble_models_vertical_Figure41_H_magnitude.png)

![Channel phase](../compare_bubble_models_vertical_Figure42_H_phase.png)

![Bubble excess TL](../compare_bubble_models_vertical_Figure43_delta_TL.png)

![Direct and reflected components](../compare_bubble_models_vertical_Figure44_components.png)

![Bubble metadata summary](../compare_bubble_models_vertical_Figure45_bubble_meta.png)

### 8.7 通信 pilot-LS 结果

生成脚本：

```text
comm_compare_bubble_models_vertical.m
```

结果文件：

```text
comm_compare_bubble_models_vertical_result.mat
```

当前正式通信比较模式：

```text
comm_mode = scalar_h_total
csi_mode = pilot_ls
pilot_len = 64
```

| 场景 | `h_true` | BER 0/10/20 dB | SER 0/10/20 dB | `h_hat` rel err 0/10/20 dB | effective SNR 0/10/20 dB | invariant error |
|---|---:|---|---|---|---|---:|
| no_bubble | 0.056546 | 0.4785 / 0.4455 / 0.2815 | 0.7180 / 0.6475 / 0.3900 | 1.0823 / 0.50676 / 0.17804 | -21.906 / -11.912 / -2.0451 | 8.0553e-18 |
| level0_empirical | 0.038889 | 0.48225 / 0.4605 / 0.36575 | 0.7255 / 0.6810 / 0.5175 | 1.5737 / 0.73685 / 0.25887 | -25.157 / -15.163 / -5.2966 | 7.7001e-18 |
| hall1d_default | 0.056361 | 0.4785 / 0.4450 / 0.2820 | 0.7185 / 0.6470 / 0.3910 | 1.0858 / 0.50842 / 0.17862 | -21.934 / -11.940 / -2.0736 | 7.3598e-18 |
| hall1d_calibrated_3dB | 0.039987 | 0.48275 / 0.4605 / 0.35975 | 0.7265 / 0.6795 / 0.5085 | 1.5305 / 0.71662 / 0.25176 | -24.915 / -14.921 / -5.0548 | 7.4534e-18 |
| hall1d_high_bubble_wind | 0.050570 | 0.47675 / 0.44675 / 0.30975 | 0.7170 / 0.6560 / 0.4310 | 1.2102 / 0.56664 / 0.19908 | -22.876 / -12.882 / -3.0152 | 8.1827e-18 |

图：

![BER vs EbN0](../comm_compare_bubble_models_vertical_BER_vs_EbN0.png)

![SER vs EbN0](../comm_compare_bubble_models_vertical_SER_vs_EbN0.png)

![Channel estimation error](../comm_compare_bubble_models_vertical_channel_est_error_vs_EbN0.png)

![Communication delta TL](../comm_compare_bubble_models_vertical_delta_TL.png)

![Communication channel magnitude](../comm_compare_bubble_models_vertical_H_magnitude.png)

![Communication channel phase](../comm_compare_bubble_models_vertical_H_phase.png)

## 10. 通信模块原理与当前实现

当前通信验证重点是避免使用已知存在时延和 tap 长度问题的 wideband `h_bb` BER/SER 结果作为正式结论。因此当前正式通信比较使用标量等效信道 `h_total`。

### 9.1 perfect CSI scalar benchmark

perfect CSI 是理想基准：

```matlab
rx_clean = tx_symbols * h_true;
rx_eq = rx_noisy / h_true;
```

该模式用于确认在已知真实信道的情况下，气泡造成的标量信道衰减如何影响误码率。它是理想基准，不代表实际接收机需要估计信道。

### 9.2 pilot-LS scalar mode

当前默认模式是 `scalar_h_total + pilot_ls`。帧结构：

```text
tx_frame = [pilot_symbols; data_symbols]
```

其中：

```text
pilot_len = 64
pilot_symbols = 1 + 0i
```

标量信道：

```matlab
rx_clean_frame = tx_frame * h_true;
```

LS 信道估计：

```math
\hat h =
\frac{
\sum \operatorname{conj}(x_{\rm pilot}) y_{\rm pilot}
}{
\sum |x_{\rm pilot}|^2
}.
```

数据均衡：

```matlab
rx_eq_data = y_data / h_hat;
```

BER/SER 只对 data bits 计算，pilot symbols 不计入误码统计。该模式比 perfect CSI 更接近实际接收机，但仍是标量信道估计，不是宽带多径估计。

### 9.3 wideband diagnostic

保留的 `wideband_diagnostic` 路径为：

```text
H_f -> H_baseband -> h_bb -> conv(tx_symbols,h_bb,'same') -> MMSE equalization
```

`diagnose_comm_chain_vertical_result.mat` 显示该路径当前不适合作为正式 BER/SER 结果：

| 指标 | 数值 |
|---|---:|
| `H_nonzero_bins` | 2000 |
| `h_full_main_tap_index` | 1 |
| `h_full_energy_first_10` | 0.99199 |
| `h_full_energy_first_100` | 0.99282 |
| `tap_count_original` | 2000 |
| `conv_same_delay_est` | 1000 |
| original wideband BER@0/10/20 dB | 0.49575 / 0.511 / 0.505 |

诊断结论：

- ideal flat channel 可以正常随 Eb/N0 降低 BER；
- scalar PE `h_total` 信道 BER 随 Eb/N0 改善；
- wideband `h_bb` 路径因长 tap 和 `conv(...,'same')` 时延模型不匹配，BER 接近随机判决水平；
- 因此当前正式通信比较不使用 wideband `h_bb` 作为 BER/SER 结论。

未来 full wideband receiver 需要重新设计为 OFDM 或训练序列 LS/MMSE 接收机。

## 11. 代码实现逻辑与关键函数调用关系

### 10.1 `CARPE3D_vertical.m`

伪代码：

```text
function output = CARPE3D_vertical(paramsV)
    cfg = local_prepare_config(paramsV)
    [h_direct, h_reflect, h_total,
     f_axis, H_direct_f, H_reflect_f, H_f,
     idx_f_ref, bubble_meta] = propWAPE_vertical(cfg)

    output.h_direct = h_direct
    output.h_reflect = h_reflect
    output.h_total = h_total
    output.H_direct_f = H_direct_f
    output.H_reflect_f = H_reflect_f
    output.H_f = H_f
    output.f_axis = f_axis
    output.idx_f_ref = idx_f_ref
    output.bubble_meta = bubble_meta
    output.config = cfg
end
```

关键点：

- 所有新增物理参数通过 `paramsV` 进入；
- `local_prepare_config` 是默认值和参数校验边界；
- 旧输出字段保持存在；
- 新元数据通过 `output.bubble_meta` 添加。

### 10.2 `propWAPE_vertical.m`

伪代码：

```text
build x,y grid
build kx, ky, kappa2
build sponge alpha_xy
build f_axis and idx_f_ref

for each f_hz in f_axis
    initialize Gaussian source

    direct path:
        for z_tx -> z_rx
            screen = local_phase_screen(...)
            PE march one step
        H_direct_f(ifq) = field at receiver

    if enable_surface_reflection
        psi_surface = local_march_field(z_tx -> 0)
        psi_ref = pm_surface_kirchhoff_module(psi_surface)
        psi_rx = local_march_field(0 -> z_rx)
        H_reflect_f(ifq) = field at receiver
    else
        H_reflect_f(ifq) = 0

    H_f(ifq) = H_direct_f(ifq) + H_reflect_f(ifq)
end

h_direct = H_direct_f(idx_f_ref)
h_reflect = H_reflect_f(idx_f_ref)
h_total = H_f(idx_f_ref)
```

关键点：

- `local_phase_screen` 同时服务 direct loop 和 `local_march_field`；
- `local_march_field` 被 reflection path 的两段传播使用；
- 因此气泡接入对 direct/reflected 一致。

### 10.3 `bubble_environment_vertical.m`

伪代码：

```text
if enable_bubbles=false or bubble_model='off'
    c_eff_xy = c_bg * ones(ny,nx)
    alpha_bub_xy = zeros(ny,nx)
elseif bubble_model='level0_empirical'
    compute alpha_bub(z,f)
    compute delta_c_bub(z,f)
    c_eff = c_bg + delta_c_bub
    expand to xy matrices
elseif bubble_model='hall1d'
    U10 = bubble_wind_speed if provided, otherwise sea_wind_speed
    [n_a,beta] = bubble_hall_spectrum(...)
    [c_eff,alpha_bub] = bubble_effective_medium(...)
    expand to xy matrices
else
    error or future model
end
```

### 10.4 `comm_compare_bubble_models_vertical.m`

伪代码：

```text
define scenarios:
    no_bubble
    level0_empirical
    hall1d_default
    hall1d_calibrated_3dB
    hall1d_high_bubble_wind

generate common data bits
generate deterministic pilots

for each scenario
    channel = CARPE3D_vertical(paramsV)
    h_true = channel.h_total

    if comm_mode='scalar_h_total' and csi_mode='pilot_ls'
        tx_frame = [pilot; data]
        rx_clean = tx_frame * h_true
        add AWGN with deterministic seed
        estimate h_hat from pilots
        equalize data with h_hat
        demodulate data only
        compute BER/SER
    elseif csi_mode='perfect'
        equalize by h_true
    elseif comm_mode='wideband_diagnostic'
        build H_baseband and h_bb
        diagnostic only
    end
end
```

## 12. 核心文件与非核心文件区分

核心物理文件：

| 文件 | 说明 |
|---|---|
| `propWAPE_vertical.m` | PE/WAPE 传播核心 |
| `bubble_environment_vertical.m` | 气泡介质统一接口 |
| `bubble_hall_spectrum.m` | Hall1D 气泡半径谱 |
| `bubble_effective_medium.m` | 气泡等效介质 |
| `pm_surface_kirchhoff_module.m` | PM 粗糙海面和 Kirchhoff 反射 |

核心 API / 配置：

| 文件 | 说明 |
|---|---|
| `CARPE3D_vertical.m` | 公开信道 API 和参数校验边界 |

当前通信比较核心：

| 文件 | 说明 |
|---|---|
| `comm_compare_bubble_models_vertical.m` | 当前气泡通信比较主脚本 |
| `modem_psk.m` | PSK 调制解调和误码率 |
| `noise_inject_vertical.m` | AWGN 注入 |

原始或参考 demo：

| 文件 | 说明 |
|---|---|
| `comm_main_vertical_psk.m` | 原始通信 demo，目前不是主气泡比较脚本 |
| `main_vertical.m` | 原始信道说明 demo，默认大网格运行成本较高 |

诊断脚本：

| 文件 | 说明 |
|---|---|
| `compare_bubble_models_vertical.m` | 信道频响对比 |
| `sweep_hall1d_sensitivity_vertical.m` | Hall1D 参数敏感性 |
| `calibrate_hall1d_bubble_vertical.m` | 风速解耦和强度校准 |
| `diagnose_comm_chain_vertical.m` | 通信链路问题诊断 |

诊断脚本对验证很重要，但不属于主传播模型的一部分。后续修改物理模型时，应优先保证核心物理文件和核心 API 的行为稳定，再更新诊断脚本。

## 13. 当前不足与下一步

当前不足：

- plume 非均匀气泡云尚未实现；
- Hall1D 参数仍需要文献或实测数据进行物理校准；
- 当前通信正式比较是 `scalar_h_total + pilot_ls`，不是 full wideband multipath LS/OFDM 接收机；
- wideband `h_bb` 路径当前只作为诊断，因为 tap 长度和时延模型导致 BER 接近随机；
- 目前大部分结果是 `nx=128`、`ny=128` reduced-grid 验证结果，不是最终高分辨率生产仿真；
- 如果要比较相同接收端 SNR 条件，需要增加 receiver-SNR reference mode；
- direct-only Level0、reflected-only Level0 和 sound-speed-only 验证的分项数值在当前结果文件中未找到。

建议下一步：

| 优先级 | 工作 |
|---:|---|
| 1 | 用文献或实测气泡衰减数据校准 Hall1D 参数 |
| 2 | 增加接收端 SNR 参考模式，区分传播损耗和接收噪声条件 |
| 3 | 设计 full wideband 通信接收机，优先考虑 OFDM 或训练序列 LS |
| 4 | 在平均层通信结果稳定后实现 plume 非均匀气泡云 |
| 5 | 生成高分辨率生产仿真和论文级图件 |

## 14. 本报告使用的结果文件

| 文件 | 用途 |
|---|---|
| `compare_bubble_models_vertical_result.mat` | 信道频响、Level0、Hall1D default/strong 对比 |
| `sweep_hall1d_sensitivity_vertical_result.mat` | Hall1D 敏感性扫描 |
| `calibrate_hall1d_bubble_vertical_result.mat` | 风速解耦、强度校准和回归验证 |
| `diagnose_comm_chain_vertical_result.mat` | 通信链路诊断和 wideband `h_bb` 问题定位 |
| `comm_compare_bubble_models_vertical_result.mat` | 当前 `scalar_h_total + pilot_ls` 通信结果 |

## 15. 本报告引用的图文件

| 图文件 | 状态 |
|---|---|
| `../compare_bubble_models_vertical_Figure41_H_magnitude.png` | 已引用 |
| `../compare_bubble_models_vertical_Figure42_H_phase.png` | 已引用 |
| `../compare_bubble_models_vertical_Figure43_delta_TL.png` | 已引用 |
| `../compare_bubble_models_vertical_Figure44_components.png` | 已引用 |
| `../compare_bubble_models_vertical_Figure45_bubble_meta.png` | 已引用 |
| `../sweep_hall1d_sensitivity_vertical_strength_scale.png` | 已引用 |
| `../sweep_hall1d_sensitivity_vertical_wind_speed.png` | 已引用 |
| `../sweep_hall1d_sensitivity_vertical_damping_const.png` | 已引用 |
| `../sweep_hall1d_sensitivity_vertical_radius_grid_max.png` | 已引用 |
| `../sweep_hall1d_sensitivity_vertical_receiver_depth.png` | 已引用 |
| `../sweep_hall1d_sensitivity_vertical_max_alpha_summary.png` | 已引用 |
| `../calibrate_hall1d_bubble_vertical_max_delta_TL_vs_bubble_wind_speed.png` | 已引用 |
| `../calibrate_hall1d_bubble_vertical_max_alpha_vs_bubble_wind_speed.png` | 已引用 |
| `../calibrate_hall1d_bubble_vertical_max_delta_TL_vs_strength_scale.png` | 已引用 |
| `../calibrate_hall1d_bubble_vertical_max_alpha_vs_strength_scale.png` | 已引用 |
| `../calibrate_hall1d_bubble_vertical_target_matching.png` | 已引用 |
| `../comm_compare_bubble_models_vertical_BER_vs_EbN0.png` | 已引用 |
| `../comm_compare_bubble_models_vertical_SER_vs_EbN0.png` | 已引用 |
| `../comm_compare_bubble_models_vertical_channel_est_error_vs_EbN0.png` | 已引用 |
| `../comm_compare_bubble_models_vertical_delta_TL.png` | 已引用 |
| `../comm_compare_bubble_models_vertical_H_magnitude.png` | 已引用 |
| `../comm_compare_bubble_models_vertical_H_phase.png` | 已引用 |

## 16. 缺失或待补充验证项

| 项目 | 状态 |
|---|---|
| direct-only Level0 `|H_bub|/|H_0|` | 当前结果文件未找到 |
| reflected-only Level0 分量比值 | 当前结果文件未找到 |
| sound-speed-only Level0 相位差 | 当前结果文件未找到 |
| full wideband multipath LS/OFDM 通信结果 | 尚未实现 |
| plume 非均匀气泡云结果 | 尚未实现 |
| 高分辨率生产仿真结果 | 待补充验证 |
