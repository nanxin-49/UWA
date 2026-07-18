# U=5/U=8 两节点统计信道通信性能验证报告

## 1. 结论摘要

U=5 和 U=8 条件统计信道库能够在当前 QPSK、接收端 AWGN、同步和已知信道 MMSE 设置下复现 `kirchhoff_kdomain + cached PE` 的 BER/SER 统计趋势。所有四种信道来源在全部 Eb/N0 点的按信道 realization 聚类 95% 置信区间均有重叠。

两节点库可以用于后续固定 U=5/U=8 的通信算法对比研究。当前更主要的限制不是统计信道，而是通信基线：现有线性卷积截取配合循环频域 MMSE 在无噪声时仍产生约 `0.2%--0.5%` BER，因此高 Eb/N0 曲线存在算法误码底。下一阶段应优先改进同步、块边界和均衡器，再扩展更多风速节点。

## 2. 实现范围

新增接口：

- `build_communication_taps_vertical`：将外部等间隔 H(f) 或 h(t) 转换为符号间隔 taps；
- `evaluate_mpsk_channel_ensemble_vertical`：统一调制、卷积、接收端噪声、同步、MMSE、解调和置信区间；
- `sample_conditional_channel_rank_pair_vertical`：使用共同潜在高斯系数构造 full/低秩配对信道；
- `validate_two_node_communication_vertical`：完成两节点四来源通信 Monte Carlo；
- `comm_main_vertical_psk` 增加默认关闭的 `COMM_EXTERNAL_CHANNEL_FILE`，支持 MAT 文件中的 H(f) 或 h(t)。环境变量为空时原 PE 路径不变。

噪声只由 `noise_inject_vertical` 在信道卷积之后加入。所有保存的 H(f)、物理 CIR 和通信 taps 均不含噪声。

## 3. 验证配置

| 项目 | 设置 |
|---|---|
| 风速节点 | 5、8 m/s |
| 信道来源 | kdomain+cached PE、joint-kstat+cached PE、统计 full、统计 99.9% |
| 新通信信道数 | 每风速、每来源 32 条 |
| 配对秩测试 | 每风速 128 对，共享潜变量 |
| 频率轴 | 4--8 kHz，F=64 |
| 调制 | QPSK |
| 每条信道符号数 | 8000 |
| 每条 BER 曲线总比特 | 512,000 |
| Eb/N0 | 0:2:20 dB |
| 符号率 | 1 ksym/s |
| 噪声参考 | `rx_clean`，与现有通信入口一致 |
| 同步 | 最短循环 99.9% 能量窗解绕为因果 taps |
| 均衡 | 已知有效 taps 的频域 MMSE |
| 置信区间 | 2000 次按信道聚类 bootstrap；同时保存 Wilson CI |

通信 bits seed 和同一 channel-index/Eb 点的噪声 seed 在所有来源之间相同；各信道来源使用与训练、接收端验证均不重叠的新 seed。

## 4. H(f)、物理 CIR 与通信 taps

F=64、4 kHz 带宽给出物理分辨率 `1/B=0.25 ms`，无模糊时延 `1/df=15.75 ms`。通信 taps 的间隔为 `1 ms`，它是 1 ksym/s 接收机的符号间隔投影，不是新的物理分辨率。

直接从未移位 IFFT 索引 1 累积能量会因循环时延和带限旁瓣对所有样本保留 2048 taps；旧 peak-sync 在 U=8 smoke 中还可能丢弃 61.8% 峰前能量。验证路径因此使用最短循环 99.9% 能量区间并解绕到因果时间轴。正式样本的中位 taps 数约为：

| U | kdomain | joint | stats full | stats 99.9% |
|---:|---:|---:|---:|---:|
| 5 | 28.5 | 28.5 | 67.0 | 28.5 |
| 8 | 65.5 | 77.0 | 132.5 | 98.5 |

实际同步丢弃能量为 0，保留能量不低于 0.999，H↔IFFT 数值闭合误差约 `1e-15`。零填充仅用于表示和插值，没有被解释为时延分辨率提升。

## 5. BER/SER 结果

### 5.1 代表性 BER

| U | 来源 | 0 dB | 8 dB | 20 dB |
|---:|---|---:|---:|---:|
| 5 | kdomain PE | 0.13675 | 0.00849 | 0.00222 |
| 5 | joint cached PE | 0.14817 | 0.01658 | 0.00378 |
| 5 | stats full | 0.13454 | 0.00574 | 0.00333 |
| 5 | stats 99.9% | 0.14256 | 0.01145 | 0.00269 |
| 8 | kdomain PE | 0.18957 | 0.03793 | 0.00479 |
| 8 | joint cached PE | 0.17846 | 0.02820 | 0.00343 |
| 8 | stats full | 0.18764 | 0.01908 | 0.00530 |
| 8 | stats 99.9% | 0.20989 | 0.05181 | 0.00733 |

点估计无需逐点相同；32 个信道 realization 的聚类 CI 在所有 Eb/N0 点均相交。完整 BER、SER、均衡前结果及上下界保存在 CSV/MAT 中。

### 5.2 曲线距离

| U | 对比 | BER log10 RMSE | SER log10 RMSE | BER/SER CI overlap |
|---:|---|---:|---:|---:|
| 5 | joint vs kdomain | 0.242 | 0.241 | 1.0 / 1.0 |
| 5 | stats full vs kdomain | 0.127 | 0.122 | 1.0 / 1.0 |
| 5 | stats full vs joint | 0.245 | 0.240 | 1.0 / 1.0 |
| 8 | joint vs kdomain | 0.086 | 0.088 | 1.0 / 1.0 |
| 8 | stats full vs kdomain | 0.194 | 0.189 | 1.0 / 1.0 |
| 8 | stats full vs joint | 0.159 | 0.156 | 1.0 / 1.0 |

因此没有发现 joint-kstat 与条件统计生成器之间超出当前有限信道样本不确定度的系统性能偏差。

### 5.3 U=5 与 U=8 趋势

按全部 Eb/N0 点平均，U=8/U=5 BER 比为：

- kdomain PE：`1.93`；
- joint cached PE：`1.35`；
- stats full：`1.69`；
- 独立 stats 99.9% 样本：`2.17`。

所有来源均给出 U=8 较差的趋势。趋势来自更强且更高秩的散射频率选择性；因为噪声按 `rx_clean` 标定，结果不混入简单接收功率差异。

## 6. Full rank 与低秩

主四来源的独立 32 信道比较中，full/99.9% 的 CI 全部重叠。为去除两批随机信道的 Monte Carlo 波动，又生成了 128 对共享潜在系数的 full/99.9% 信道：

- U=5：rank 34 对 rank 7；0--18 dB 的 BER 差 CI 包含 0，20 dB 低秩增加 `1.08e-4`，95% CI `[2.15e-5,2.04e-4]`；
- U=8：rank 42 对 rank 14；全部 Eb/N0 的 BER 差 CI 均包含 0。

99.9% 低秩在当前通信基线下没有工程上明显退化，可用于压缩实验；full rank 仍应保留为默认物理基准。99% 模型未在本阶段升格为通信默认压缩方案。

## 7. 均衡前后与误码底

未均衡 BER 约为 0.50，说明频率选择性和随机复相位使直接硬判决不可用。MMSE 后 0 dB BER 降至约 0.13--0.21，高 Eb/N0 降至约 0.002--0.007。

无噪声 BER 仍为：

- U=5：0.00224（kdomain）、0.00372（joint）、0.00328（stats full）；
- U=8：0.00462、0.00330、0.00536。

它与 20 dB 误码底接近，表明误码底属于当前通信处理，而非噪声或统计信道生成。主要原因是用循环 FFT 逆滤波近似有限长度线性卷积、帧首尾无保护间隔、0.999 能量截断以及固定 MMSE 正则化。绝对高 SNR 性能研究前应加入 overlap-save/循环前缀或线性均衡基线，并分别处理训练、同步和边界瞬态。

## 8. 六个核心问题的回答

1. **统计生成器能否复现 kdomain BER/SER？** 能。在两个风速、全部 Eb/N0 点的聚类 95% CI 均重叠，曲线 log-RMSE 为 0.127/0.122（U=5）和 0.194/0.189（U=8）。
2. **joint 与统计生成器是否有明显偏差？** 未发现超出当前有限信道样本置信范围的明显偏差；CI overlap 为 1.0。
3. **风速趋势是否合理？** 是。四种来源都预测 U=8 BER 高于 U=5，stats full 的平均比值为 1.69。
4. **低秩影响？** 99.9% 低秩总体无明显退化；U=5 仅 20 dB 出现很小但可检出的差值。full rank 仍为默认。
5. **两节点库能否正式用于通信研究？** 可以用于 U=5/U=8 固定节点的相对算法比较、Monte Carlo 和回归。不能把当前 MMSE 的高 SNR 误码底当作信道物理极限。
6. **下一步扩风速还是改通信算法？** 优先改进通信算法和帧/均衡结构。两个节点已经足以暴露同步、边界和均衡器问题；建立可信的无噪声零误码基线后，再扩展 U=10。

## 9. 产物与回归

- 完整 MAT：`results/validation/two_node_communication/two_node_communication_validation_full.mat`；
- 汇总 CSV：`two_node_comm_summary_full.csv`；
- U=5/U=8 BER/SER、均衡前后 PNG；
- fresh channel checkpoint：`two_node_comm_ensembles_full.mat`。

`comm_main_vertical_psk` 的外部 H(f) 和 h(t) smoke 均通过；清空环境变量后的原 direct-only/direct-plus-reflect reduced-Nf 回归也通过。公共传播默认值、噪声注入位置、两节点模型和条件库均未改变。
