# PE--Bellhop 无反射自由场验证报告

状态：`passed=false`。本验证不含海面或海底反射、粗糙面、气泡、Doppler、随机信道或通信处理。

## 环境与数值设置

- 均匀声速：`1500 m/s`；频率：`[3000 4000 5000] Hz`。
- 虚拟源到初始面 `s0=[5 10 20] m`；PE 推进 `L=[20 40 70 100] m`；横向偏移 `[0 0.5 1 2] m`。
- PE 正式物理比较：`256 x 256`，横向宽度 `64.0 m`，`stepz_lamb=0.5`，`alpha_max=0.15 Np/m`。
- PE--AS hard gate 的 sponge 完全关闭；Bellhop 使用匹配上下半空间 free-space 构造。

## 四层结论

1. Bellhop normalization audit：通过；`|p|R` 最大均值偏差 `3.68283e-05`，统一 Green 转换为 `1/(4*pi)`，不做逐点拟合。
2. PE--独立一步 AS：通过；最大全场相对 L2 误差 `1.4287e-13`。
3. PE--解析球面波：未通过；最大 TL 误差 `18.2012 dB`。
4. PE--Bellhop：未通过；最大 TL 差 `18.2012 dB`，相位 RMS `1.38356 rad`；Bellhop--解析最大 TL 差 `1.99322e-06 dB`。

## 双相位参考

`H_plane = Psi(L)*exp(i*k*L)` 用于 PE--AS；`H_source = exp(i*k*s0)*H_plane` 用于解析球面波和 Bellhop。公共 PE 输出语义未改变。

## 解释与判定

PE--AS 已达到浮点误差量级，因此当前结果不支持“生产 PE 自由场平方根传播公式写错”的假设。失败集中在初始球面波被截断到有限横向平面并由周期 FFT 延拓；窗口从 32 m 增至 64 m 虽改善误差，但尚未收敛。下一步应独立验证无限孔径/Weyl 谱初始化或扩大无 sponge 孔径，而不是修改 PE marching 主线。

代表性 `R=80 m` 群时延：PE 误差 `3765.72 us`，Bellhop 误差 `0.00126667 us`。

## 自动检查

| check | value | limit | passed |
|---|---:|---:|:---:|
| pe_as_gate | 0 | 0 | 1 |
| bellhop_normalization_gate | 0 | 0 | 1 |
| pe_bellhop_tl | 18.201226 | 0.5 | 0 |
| pe_bellhop_phase | 1.3835644 | 0.1 | 0 |
| pe_group_delay | 0.00376572 | 5e-06 | 0 |
| bellhop_group_delay | 1.2666667e-09 | 5e-06 | 1 |
| bellhop_analytic_tl | 1.9932179e-06 | 0.5 | 1 |
| bellhop_analytic_phase | 4.3381516e-05 | 0.1 | 1 |

![自由场比较](../results/validation/pe_bellhop_freefield/formal/pe_bellhop_freefield_comparison.png)
