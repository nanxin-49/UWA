# 完整 PE -> 海面 -> PE 链路窗口收敛验证报告

状态：**CONVERGED**。本次只增加验证接口和诊断，不修改 PE marching、FFT 约定、Gaussian 源、Kirchhoff 反射公式或生产默认窗口/sponge。

## 1. 环境与方法

- 均匀介质，`c=1500 m/s`；Tx=`(0,0,100 m)`，Rx=`(0,0,3 m)`。
- 生产 Gaussian 初场，`sigma_src_m=0.3 m`；4 kHz 单频先验收，再运行 3--5 kHz、33 点宽带。
- Kirchhoff spatial 海面，`U=5 m/s`、目标 `Hs=0.5 m`、反射系数 `-1`；关闭气泡、Doppler、通信处理和 Monte Carlo。
- `dx=dy=50/256=0.1953125 m`。在名义 384 m、`1966 x 1966` 最大网格上用 seed 12345 生成一次 PM 海面，并且只归一化一次；各窗口从同一主海面中心逐点裁剪，不再重标 `Hs`。
- 单频窗口全部关闭 sponge。比较 `H_reflect_physical_f`，以及海面入射场、海面反射场、接收深度反射场和上/下行能量轨迹。

## 2. 4 kHz 单频结果

| 名义窗口 | 实际窗口 | 裁剪区实际 Hs | 反射 TL | 反射相位 | 最大外围 5% 能量占比 | 最差边界幅度 | 边缘门槛 |
|---:|---:|---:|---:|---:|---:|---:|:---:|
| 50 m | 50.0000 m | 0.54706 m | 39.119 dB | 2.2411 rad | 1.1839e-1 | -1.23 dB | FAIL |
| 128 m | 128.1250 m | 0.52340 m | 37.971 dB | 2.2439 rad | 1.0681e-3 | -28.86 dB | FAIL |
| 160 m | 160.1563 m | 0.52193 m | 37.960 dB | 2.2444 rad | 1.1769e-4 | -41.14 dB | FAIL |
| 192 m | 192.1875 m | 0.51434 m | 37.961 dB | 2.2438 rad | 1.7699e-5 | -52.29 dB | PASS |

| 窗口对 | 反射 TL 差 | 相位差 | 入射中心复场 L2 | 接收中心复场 L2 | 综合结果 |
|---|---:|---:|---:|---:|:---:|
| 50 -> 128 m | 1.1479 dB | -0.00281 rad | 13.263% | 14.982% | FAIL |
| 128 -> 160 m | 0.01171 dB | -0.000573 rad | 0.0263% | 0.0918% | 仅边缘门槛失败 |
| 160 -> 192 m | -0.00154 dB | 0.000612 rad | 0.00654% | 0.0455% | PASS |

160 m 相对 192 m 的接收复场、TL 和相位已经高度一致，但外围 5% 能量占比 `1.1769e-4` 略高于预先固定的 `1e-4` 门槛。因此，**160 m 是接收响应近似已收敛的低成本候选，但不是本轮严格合格窗口**；严格最小窗口为实际 `192.1875 m`。

阶段定位显示污染首先出现在上行 100 m PE：50/128/160/192 m 的上行最大外围 5% 能量占比分别约为 `1.1839e-1 / 7.8586e-4 / 7.4853e-5 / 8.2650e-6`。海面相位屏保持幅度不变量，随后 3 m 下行传播仅产生较小的额外边缘扩展。Kirchhoff 相位屏幅度不变量相对误差小于 `1e-16`。

## 3. 3--5 kHz、33 点结果

宽带以 `192.1875 m + no sponge` 为参考：

| 配置 | 最大反射 TL 差 | 最大反射相位差 | 最大内部群时延差 | 最大总信道 TL 差 | 最大总信道相位差 | 结果 |
|---|---:|---:|---:|---:|---:|:---:|
| 192.1875 m / no sponge | 0 | 0 | 0 | 0 | 0 | PASS |
| 160.1563 m / no sponge | 0.01921 dB | 0.002382 rad | 4.05 us | 0.03079 dB | 0.004021 rad | PASS（宽带响应） |
| 50 m / ratio 0.12 / alpha 0.15 | 1.3203 dB | 0.17782 rad | 70.94 us | 8.4134 dB | 0.64120 rad | FAIL |

所有频点均满足 `H_f=H_direct_f+H_reflect_f`，最大闭合残差不超过 `3.9e-18`。正式宽带计算采用 GPU double；4 kHz 与保存的 CPU double 参考做过逐复数值审计，相对差为 0。

## 4. 直接结论

1. **反射路径严格所需窗口为实际 192.1875 m。** 这是同时满足接收响应、中心二维复场和传播全过程边缘门槛的最小已测窗口。
2. **160 m 对最终接收 `H_reflect(f)` 已足够精确，但对预注册边缘门槛不够。** 它在宽带的最大反射 TL/相位/群时延差只有 `0.01921 dB / 0.002382 rad / 4.05 us`，但 4 kHz 外围 5% 能量超限约 17.7%。不能把它标为严格通过。
3. **推荐完整链路使用 192.1875 m + sponge off。** 该配置无 sponge 已通过边缘门槛，因此按固定决策规则不启动弱 sponge 扫描。当前 50 m/default sponge 明显改变反射和总信道，不适合作为已经验证的生产配置。
4. 本轮没有自动修改生产默认值；正式迁移仍需单独评审计算成本和数据平台容量。

## 5. 产物与限制

正式 MAT/CSV、逐阶段能量表、检查点及六张图位于 `results/validation/pe_reflected_chain_window/`。核心文件为 `reflected_chain_window_validation.mat`、`reflected_chain_wideband_validation.mat`、`single_frequency_cases.csv`、`reference_window_metrics.csv`、`stage_energy_traces.csv` 和 `wideband_summary.csv`。

结论仅适用于本次固定 Gaussian 源、固定 seed 海面、单次 realization、近轴接收点和 3--5 kHz 范围；它不是随机海面 ensemble 或 Monte Carlo 收敛结论。弱 sponge 分支未执行，因为 no-sponge 候选已满足全部严格门槛。

## Wideband qualification

| case | width (m) | max reflect TL (dB) | max reflect phase (rad) | max internal GD (us) | max total TL (dB) | max total phase (rad) | pass |
|---|---:|---:|---:|---:|---:|---:|:---:|
| reference_no_sponge | 192.187500 | 0 | 0 | 0 | 0 | 0 | PASS |
| nominal160_no_sponge | 160.156250 | 0.0192053 | 0.00238167 | 4.05067 | 0.0307945 | 0.00402087 | PASS |
| production_default | 50.000000 | 1.32035 | 0.177821 | 70.9435 | 8.41341 | 0.641202 | FAIL |

The physical-total and communication `H_f` comparison columns are both saved in `wideband_reflected_total_comparison.csv`; their window errors agree under the common direct-DSP phase reference. Recommended configuration: W=192.188 m, sponge=off.
