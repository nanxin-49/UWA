# 独立伴随 PE 接收投影可行性报告

日期：2026-07-16

## 结论

当前 cached uniform PE 在 CPU double、固定网格、固定发射机和单个最近网格点接收机条件下具有可实现到双精度舍入误差的精确离散伴随。伴随接收投影可以在固定环境中完全替代每条 realization 的 surface-to-receiver PE；PM 中央裁剪后的解析 FFT 协方差和伪协方差也与稠密基准及 realization 路径一致。

建议把本原型作为独立研究与验证模块接入主分支，但暂不接入公共默认传播路径、默认 surface model 或通信入口。v1 不支持 layered medium、bubble、Doppler、GPU、多接收机或插值接收。

## 实现范围与接口

新增验证专用接口：

- `apply_forward_surface_to_receiver_vertical`：共享 cached surface-to-receiver 前向算子，支持 `ny x nx x M` CPU double pages。
- `apply_adjoint_receiver_to_surface_vertical`：上述离散算子的共轭转置，反向遍历深度步。
- `build_adjoint_receiver_projection_vertical`：构造每频接收核 `q_surface_xy_f`、PE 权重 `a_pe_xy_f` 及 PM 映射元数据。
- `run_adjoint_receiver_projection_vertical`：对同一批 joint-kstat realization 直接计算接收端散射响应。
- `contract_kstat_receiver_stats_vertical`：使用 `dense` 或 `fft` 方法收缩接收端 `C_H/P_H`。

`run_cached_joint_kstat_pe_executor_vertical` 仅把原内联前向循环改为调用共享 primitive，输入、输出字段及数值结果不变。公共 `vertical_channel_model`、`vertical_wape_propagator`、默认 surface model 和通信入口均未修改。

## 数学口径

随机变量沿用现有 joint model：

\[
\delta G_i=R_{0,i}\exp(i\alpha_i\eta)-R_{{\rm coh},i}.
\]

当前 `C_deltaG/P_deltaG` 已包含 `R0`，所以新权重不再次乘 `R0`。令

\[
a_{{\rm PE},i}=\operatorname{conj}(\psi_{{\rm inc},i})\odot q_i,
\]

PM 权重由现有中央裁剪算子伴随零嵌入：

\[
\widetilde a_i=E^Ha_{{\rm PE},i}.
\]

接收端统计为

\[
C_H(i,j)=\widetilde a_i^HC_{\delta G,ij}\widetilde a_j,
\qquad
P_H(i,j)=\widetilde a_i^HP_{\delta G,ij}\widetilde a_j^*.
\]

总均值单独构造为

\[
\mu_{H,{\rm total}}=H_{\rm direct}+H_{{\rm ref,coh}}+\mu_{H,{\rm sca}},
\qquad \mu_{H,{\rm sca}}=0.
\]

FFT 收缩采用未 shift 的 MATLAB FFT 排列，lag 零点在 `(1,1)`，没有经验性 `nx*ny` 补偿。

## 验证结果

### 精确伴随、投影和稠密基准

| 检查 | 最大相对误差 | 门限 | 结果 |
|---|---:|---:|---|
| 4/6/8 kHz 离散伴随内积，5 组随机输入/频率 | `7.77e-15` | `1e-10` | 通过 |
| 随机场、平整场、explicit Kirchhoff、joint-kstat 投影 | `5.10e-15` | `1e-10` | 通过 |
| 原 cached 循环与共享 forward primitive | `0` | `1e-12` | 通过 |
| PM 16² / PE 8²，dense 与 FFT 的 `C_H` | `6.18e-16` | `1e-10` | 通过 |
| PM 16² / PE 8²，dense 与 FFT 的 `P_H` | `9.35e-16` | `1e-10` | 通过 |

零嵌入中央区域与 PE 权重位级一致，中央区域外严格为零。稠密实现逐频率对流式构造，仅作为不超过 32² PM 点的验证基准。

### 宽带 realization 与统计

| 指标 | F=9，4096 条 | F=64，512 条 |
|---|---:|---:|
| forward/projection 最大相对误差 | `2.52e-13` | `1.07e-13` |
| covariance relative error | `0.02418` | `0.08866` |
| covariance error / split floor | `0.593` | `0.851` |
| pseudo-covariance error / `||C||` | `0.02569` | `0.06036` |
| pseudo error / split floor | `0.363` | `0.577` |
| 每频功率相对误差 | `0.01403` | `0.05864` |
| reflected-only PDP correlation | `0.999876` | `0.999784` |
| reflected-only LFM correlation | `0.999953` | `0.998688` |
| reflected-only matched-filter correlation | `0.999971` | `0.999774` |
| 最大散射均值标准误差比 | `1.388` | `1.782` |

两种尺度的解析—样本误差均低于对应 split-sample floor 的 1.25 倍。理论 `P` 接近零，因此传统 `||Phat-P||/||P||` 数值会被极小分母放大；验收使用 `C` 归一化误差和相同 realization 数目的 split floor，同时在 MAT 结果中保留传统指标。

## 性能与内存

F=9、PE 64²、PM 128²：

- 核构造：`0.0283 s`。
- 4096 条 cached forward 总计：`81.96 s`，约 `0.0200 s/条`。
- 4096 条投影总计：`1.661 s`，约 `0.000405 s/条`，批处理加速约 `49.4x`。
- 解析 `C/P`：`0.0346 s`。

F=64、PE 128²、PM 256²：

- 核构造：`0.497 s`。
- 512 条 cached forward 总计：`379.15 s`，约 `0.741 s/条`。
- 512 条投影总计：`10.76 s`，约 `0.0210 s/条`，批处理加速约 `35.2x`。
- 解析 `C/P`：`7.10 s`。
- `q+a_PE`：约 `32 MiB`；PM 嵌入/FFT 权重：约 `128 MiB`。
- 若错误保存完整空间 `F² C/P`，估算约 `563 TB`；可行实现不保存这些块，只保留权重、一个流式 lag kernel 和最终 `F x F` 结果。

F=64 的主要解析瓶颈是 2080 个上三角频率对对应的 PM lag FFT 与收缩，即 `F²` 频率对，而不是约 0.5 秒的一次性核构造。投影的在线成本远小于 cached forward PE。

## 公共回归

- `validate_cached_pe_public_consistency_vertical`：通过，cached/public total 相对误差 `7.04e-16`。
- `validate_public_channel_modes_vertical`：通过；`direct_only` 反射严格为零，`direct_plus_reflect` 分量闭合误差 `1.73e-18`。
- `validate_comm_link_minimal_vertical`：通过；reduced-grid QPSK peak-sync 无噪声 BER/SER 为零。
- 公共输出 `H_f=H_direct_f+H_reflect_f` 以及 `h_direct/h_reflect/h_total` 与 `idx_f_ref` 语义保持不变。

运行中出现的 BELLHOP 路径警告来自已有环境配置，不影响本次 uniform PE 验证。

## 最终判断

1. 当前 cached uniform PE 具有可实现到数值精度的精确离散伴随。
2. 对固定环境、固定单接收机和当前最近网格点采样，投影核可以完全替代每条 realization 的 surface-to-receiver PE。
3. PM 中央裁剪经 `E^H` 零嵌入后，解析 `C_H/P_H` 与稠密定义及 realization 路径一致。
4. F=64、PE 128²、PM 256² 的主要解析瓶颈是 `F²` 频率对的 PM FFT/收缩；核构造不是主要瓶颈。
5. 建议把独立原型、测试和报告接入主分支；在扩展范围得到验证前，不建议改变公共默认传播器或通信入口。
6. 推荐角色分工：伴随投影负责固定环境下快速且物理精确的接收端 realization；解析 FFT `C/P` 负责无需接收端 Monte Carlo 的二阶统计；cached forward PE 继续作为高可信回归 oracle；当前条件统计生成器在解析/样本验证后承担大规模通信 Monte Carlo，但不替代物理回归。

## 产物与剩余限制

验证入口为 `scripts/validation/validate_adjoint_pe_receiver_projection_vertical.m`。完整 MAT、文本摘要及误差、统计、特征值、时间和内存图位于 `results/validation/adjoint_pe_receiver_projection/`。

本实现是现有离散 marching operator 的 exact conjugate transpose，不是互易反传、逆传播或损耗补偿。每个新传播环境、接收位置或网格都必须重新构造核。layered、bubble、Doppler、GPU、多接收机和插值接收仍需分别定义并验证对应离散伴随。
