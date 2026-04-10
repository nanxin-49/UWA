# 海底仪器到浮标水听器 MPSK 通信模块说明

## 1. 场景与目标
本模块用于仿真以下通信链路：

- 发射端：海底仪器（默认 `z_tx = 100 m`）
- 接收端：浮标下悬挂水听器（当前固定 `z_rx = 3 m`）
- 调制方式：MPSK（默认 QPSK, `M = 4`）
- 信道：上行 PE 直达路径 + 可选粗糙海面反射路径

核心开关：

- `enable_surface_reflection = false`：仅直达，`h_total = h_direct`
- `enable_surface_reflection = true`：直达 + 反射，`h_total = h_direct + h_reflect`

---

## 2. 物理模型与默认参数

### 2.1 坐标与传播约定
- `z` 轴向下为正。
- 传播方向为上行：从 `z_tx` 推进到 `z_rx`。
- 发射与接收横向坐标位于 `(x,y)` 平面网格内。

### 2.2 默认物理参数（主链路）
| 参数 | 默认值 | 含义 |
|---|---:|---|
| `f0` | `4000 Hz` | 载频 |
| `c0` | `1500 m/s` | 参考声速 |
| `z_max` | `100 m` | 水深上限 |
| `z_tx` | `100 m` | 发射深度 |
| `z_rx` | `3 m` | 接收深度（当前阶段固定） |
| `xw, yw` | `50 m, 50 m` | 横向计算域宽度 |
| `nx, ny` | `1024, 1024` | 横向网格数 |
| `stepz_lamb` | `0.5` | 垂向步长（波长单位） |
| `sigma_src_m` | `0.3 m` | 初始高斯源宽度 |
| `taper_ratio` | `0.12` | 横向海绵层比例 |

### 2.3 粗糙海面（PM 谱）默认参数
| 参数 | 默认值 | 含义 |
|---|---:|---|
| `U` | `5 m/s` | 风速 |
| `Hs_target` | `0.5 m` | 目标有效波高（标定后） |
| `seed` | `12345` | 随机种子 |
| `g` | `9.81` | 重力加速度 |
| `alpha_PM` | `8.10e-3` | PM 常数 |
| `beta_PM` | `0.74` | PM 常数 |

PM 波数谱采用：

\[
W(K)=\frac{\alpha_{PM}}{2K^3}\exp\left[-\beta_{PM}\frac{g^2}{U^4K^2}\right], \quad K=\sqrt{K_x^2+K_y^2}
\]

`K=0` 处能量置 0，避免奇点。

### 2.4 海面反射相位畸变（Kirchhoff 近似）
\[
\Delta\phi = 2\pi \frac{2\xi(x,y)}{\lambda_0}
\]
\[
\psi_{ref} = \psi_{inc}\cdot e^{j\Delta\phi}
\]

---

## 3. 文件与函数说明

## 3.1 `explain_main_vertical.m`
用途：单次信道求解入口脚本（会保存图与 `.mat`）。

关键点：
- 设置 `paramsV`（含几何、反射、海况等）
- 调用 `CARPE3D_vertical(paramsV)`
- 导出 Figure11~Figure15

## 3.2 `CARPE3D_vertical.m`
用途：参数检查、求解器调度、结果打包、绘图。

输入：
- `paramsV`（结构体）

输出（核心字段）：
- `h_direct`, `h_reflect`, `h_total`
- `rx_state_used`（实际使用的接收点状态）
- `fd_hz_used`（当前多普勒占位输出，默认 0）
- `surface_elevation`, `psi_ref`, `delta_phi`
- `direct_to_reflect_db`, `phase_diff_rad`

参数约束（关键）：
- `z_tx > z_rx >= 0`
- 接收点在网格内
- `nx, ny <= 2048`
- `xw, yw <= 100 m`

## 3.3 `propWAPE_vertical.m`
用途：上行 PE 核心推进。

主要流程：
- 在 `(x,y)` 截面构造初场
- split-step WAPE 推进至接收深度
- 提取 `h_direct`
- 可选反射分支：
  - 先推进到海面得到 `psi_surface_inc`
  - PM + Kirchhoff 得到 `psi_ref`
  - 从海面再推进回接收深度得到 `h_reflect`
  - 合成 `h_total`

未来接口（已预留）：
- `rx_position_fn(t_s, state)`：动态接收点
- `doppler_fn(t_s, tx_state, rx_state, env_state)`：多普勒频移

## 3.4 `pm_surface_kirchhoff_module.m`
用途：粗糙海面生成 + 反射相位畸变 + sanity plot。

输入：
- `psi_inc, KX, KY, x, y, xw, yw, lambda0, pm_cfg`

输出：
- `surface_elevation, delta_phi, psi_ref, meta`

`pm_cfg`：
- `U`, `Hs_target`, `seed`, `show_figure`

## 3.5 `noise_inject_vertical.m`
用途：统一噪声注入接口（可开关、可替换模型）。

签名：
- `[rx_noisy, noise, meta] = noise_inject_vertical(rx_clean, noise_cfg, signal_ref)`

`noise_cfg` 字段：
- `enable_noise`：是否加噪
- `model`：默认 `'awgn'`
- `snr_db` 或 `ebn0_db`
- `bits_per_symbol`
- `seed`
- `custom_noise_fn`（自定义噪声函数句柄）

行为：
- 关闭噪声：直接透传
- 开启噪声：按配置注入，返回 `effective_snr_db`

## 3.6 `comm_main_vertical_psk.m`
用途：端到端 MPSK 通信主脚本。

流程：
- 构造 `paramsV`（固定 `z_rx = 3`）
- 运行两工况：
  - `direct_only`
  - `direct_plus_reflect`
- 调制 -> 信道 -> 噪声注入 -> 相干均衡 -> 解调
- 输出 BER/SER 表格与对比图（Figure31/32）
- 保存 `psk_comm_result.mat`

## 3.7 `modem_psk.m`
用途：MPSK 调制/解调与误码统计工具函数。

支持模式：
- `modulate`
- `demodulate`
- `error_rate`

---

## 4. 常用运行方式

### 4.1 单次信道仿真
```matlab
explain_main_vertical
```

### 4.2 通信闭环仿真
```matlab
comm_main_vertical_psk
```

---

## 5. 结果字段速查

`CARPE3D_vertical` 输出关键字段：

- `h_direct`：直达复通道
- `h_reflect`：反射复通道（关闭反射时为 0）
- `h_total`：总复通道
- `rx_amplitude`, `rx_phase_rad`
- `direct_to_reflect_db`
- `phase_diff_rad`
- `rx_state_used`
- `fd_hz_used`

`comm_main_vertical_psk` 保存的 `results(ss)`：

- `name`（工况名）
- `BER`, `SER`
- `effective_snr_db`
- `h_direct`, `h_reflect`, `h_total`
- `fd_hz_used`, `rx_state_used`

---

## 6. 当前假设与后续扩展

当前假设：
- 窄带等效复通道（单载频）
- 接收点固定 3m
- 多普勒占位输出，默认 `fd_hz_used = 0`
- 噪声默认 AWGN

后续扩展建议：
- 通过 `rx_position_fn` 引入浮标随浪运动轨迹
- 通过 `doppler_fn` 注入时变频偏
- `custom_noise_fn` 接入海洋环境噪声模型
- 扩展到更高阶 MPSK 与帧同步/导频估计

