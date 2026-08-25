# 展开坐标 PE--Bellhop 平面海面 Gaussian 验证

状态：`passed=false`。主方法将原竖直传播展开为 Bellhop 水平自由场，并在展开接收端施加压力释放系数 `-1`。

## 环境

- `c=1500 m/s`，Tx 深度 `100 m`，Rx 深度 `3 m`，Gaussian `sigma=0.3 m`。
- 直达展开距离 `97 m`，反射展开距离 `103 m`；PE 窗口 `32 m`、`N=128`、sponge off。
- 频率轴 `3` 点：`4--8 kHz`；Bellhop 使用 `10001` 条波束，5001 条用于收敛审计。

## 结果摘要

| check | value | limit | pass |
|---|---:|---:|:---:|
| bellhop_beam_tl | 8.2000248e-07 | 0.1 | 1 |
| bellhop_beam_phase | 4.8354761e-06 | 0.02 | 1 |
| spatial_tl | 11.2276 | 0.25 | 0 |
| spatial_phase | 2.9580809 | 0.05 | 0 |
| spatial_complex_l2 | 2.4364961 | 0.02 | 0 |
| wideband_tl | 1.8027395 | 0.25 | 0 |
| wideband_phase | 0.94209107 | 0.05 | 0 |
| group_delay | 3.9288218e-05 | 2e-05 | 0 |
| pdp_delay | 6.25e-05 | 0.000125 | 1 |
| unfolded_pe_identity | 2.5103905e-13 | 1e-10 | 1 |
| native_auxiliary | 0 | 0 | 1 |
| H_closure | 3.469447e-18 | 1e-12 | 1 |

至少一个门槛未通过；失败项只作为诊断，不修改阈值，也不反向改变 PE 主线。

本结论不覆盖粗糙海面、深度相关声速、海底、气泡、Doppler 或实际换能器。更换换能器后必须用测量/建模的幅相指向性重新生成 `.sbp`，并重新确认 PE 横向窗口和边缘能量。

![summary](../results/validation/pe_bellhop_unfolded_flat_gaussian/pe_bellhop_unfolded_flat_gaussian.png)
