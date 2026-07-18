# Cached joint-kstat + PE 接收端验证报告

日期：2026-07-11

## 1. 结论

本阶段已经实现独立的 cached executor，并在固定 `U=5 m/s`、`raw_pm`、4--8 kHz、`F=32` 条件下完成指定规模的接收端验证。公共 `vertical_channel_model`、`vertical_wape_propagator`、默认配置和通信链路均未修改。

joint-frequency kstat 在接收端反射随机项上的表现显著优于 independent kstat：相对同一海面 `kirchhoff_kdomain` 测试集，协方差误差由 `0.9750` 降到 `0.2704`，reflected-only PDP correlation 为 `0.9714`，LFM matched-filter correlation 为 `0.9917`。`kdomain` 自身 128 样本训练集与 64 样本测试集之间的协方差差异为 `0.2610`，因此 joint 的 `0.2704` 已接近当前有限样本比较的误差底限，不能把它全部归因于模型偏差。

结论分级如下：

- cached PE executor：可直接用于下一轮单条件研究；double 同输入与公共路径总响应最大误差 `7.04e-16`。
- joint-frequency kstat + cached PE：已满足 covariance 优于 independent、PDP、LFM、分量求和、PM 映射和实际加速要求。
- 接收端条件统计生成器：**尚不建议立即宣告全部前置条件通过**。当前 64 样本测试集的样本 `||P||_F/||C||_F` 约为 `0.26`，很可能包含高维有限样本偏差，但 proper-complex-Gaussian 零假设 Monte Carlo 校准尚未执行完成；在该项完成前不能决定是否删除伪协方差。

## 2. 实现范围

新增验证接口：

- `build_kirchhoff_kstat_joint_model_vertical.m`：按跨频协方差和伪协方差构造 `K/-K` 增广谱因子；缓存截断单精度因子。
- `sample_kirchhoff_kstat_factor_model_vertical.m`：从缓存因子批量生成 joint 或 independent 随机反射相位屏。
- `build_cached_joint_kstat_pe_executor_vertical.m`：缓存发射端到海面的入射场、直达响应、相干反射接收响应、surface-to-receiver 核、介质屏、吸收层和 PM/PE 中央裁剪映射。
- `run_cached_joint_kstat_pe_executor_vertical.m`：每个 realization 只传播随机反射项，并输出 `H_direct_f`、`H_ref_coh_f`、`H_ref_sca_fm` 和 `H_total_fm`。
- `validate_cached_joint_kstat_pe_receiver_vertical.m`：执行三方训练/测试统计、CIR、PDP、LFM、分布和计时验证。
- `validate_cached_pe_public_consistency_vertical.m`：使用 double 边界输入独立验证 cached/public 数值一致性。
- `analyze_receiver_properness_null_vertical.m`：以训练协方差生成 proper 复高斯零假设，校准有限样本 `P/C`；该脚本已实现，本轮因额外 MATLAB 启动额度限制未执行。

cached executor 当前有意限制为：均匀声速、CPU double PE、无气泡、无多普勒。它是独立验证路径，不替代公共传播器。

## 3. 固定配置与随机样本

| 项目 | 设置 |
|---|---:|
| 风速 | 5 m/s |
| 海况 | `raw_pm` |
| 频带 | 4--8 kHz |
| 频点 | 32，等间隔 |
| PM 网格 | 100 m / 256 x 256 |
| PE 网格 | 50 m / 128 x 128 |
| 映射 | 相同 `dx=dy=0.390625 m`，中心裁剪，无插值、无去均值 |
| 发射/接收深度 | 100 m / 3 m |
| 训练/测试样本 | 128 / 64 |
| 主 batch size | 4 |
| 介质、气泡、多普勒 | uniform、off、off |

三种模式的训练/测试 seed 域完全分离。显式海面分别使用 `100001...100128` 和 `200001...200064`；joint 与 independent 使用互不重叠的 batch seed 段 `300xxx`、`400xxx`、`500xxx`、`600xxx`。

## 4. 分量与数值一致性

执行器采用

\[
H_{\rm total}=H_{\rm dir}+H_{\rm ref,coh}+H_{\rm ref,sca}.
\]

全部训练/测试样本的分量求和最大绝对误差为 `2.78e-17`，优于 `1e-10` 门槛。用户任务描述中的横线按分量列举符解释；项目的稳定物理接口仍保持加法不变量。

公共路径一致性使用相同 50 m/128² PM/PE 网格、相同显式海面和 double `deltaG`：

| 分量 | cached/public 最大绝对误差 |
|---|---:|
| direct | 舍入误差量级 |
| reflect | `7.04e-16` 量级 |
| total | `7.04e-16` |

主验证为了控制 100 m/256² joint batch 内存，将随机边界场保存为 single；若把同一显式边界也先量化成 single，公共/cached 总响应差为 `1.59e-9`。这不是 PE executor 不一致，而是输入精度误差。因此生产验证应明确记录边界场精度，严格一致性回归必须使用 double。

## 5. 接收端频域统计

全部指标只比较 `H_ref_sca(f)`，不使用强直达项。

| 指标（64 样本独立测试集） | joint | independent |
|---|---:|---:|
| `epsilon_C` | **0.2704** | 0.9750 |
| correlation-matrix relative error | **0.2699** | 0.9790 |
| adjacent-correlation RMSE | **0.0110** | 0.9972 |
| `epsilon_mu` | 1.536 | 1.492 |

散射项理论均值接近零，导致 `epsilon_mu` 的分母很小，两个大于 1 的数值不适合作为主要通过判据。应同时报告绝对均值和每频均值的标准误；后续统计生成器应把 `H_dir+H_ref,coh` 作为确定性条件均值，把随机散射样本均值作为零均值诊断。

训练/测试协方差相对差异分别为：kdomain `0.2610`、joint `0.1928`、independent `0.6650`。相关矩阵热图、邻频行为和特征值谱均显示 joint 保留了主要跨频结构，而 independent 基本破坏该结构。

## 6. CIR、PDP 与 LFM

PE 采用负群时延相位约定。验证脚本把任务公式中的有符号参考时延设为

\[
\tau_{\rm ref,signed}=-(z_{\rm tx}+z_{\rm rx})/c_0,
\]

从而 `exp(-i2*pi*f*tau_ref_signed)` 等价于移除 PE 的负相位斜率。频率设置给出真实时延分辨率 `1/B=0.25 ms`，离散 IFFT bin 为约 `0.242 ms`，最大无模糊时延 `1/df=7.75 ms`。零填充未被当作提高真实分辨率的手段。

| 指标 | kdomain | joint | independent |
|---|---:|---:|---:|
| PDP correlation（相对 kdomain） | 1 | **0.9714** | -0.0445 |
| LFM envelope correlation | 1 | **0.9917** | 0.5721 |
| PDP 主峰相对时延 | 0.969 ms | 1.211 ms | 0.484 ms |
| realization 平均 RMS delay | 2.521 ms | 2.817 ms | 2.378 ms |
| realization 平均 99% 尾长 | 7.334 ms | 7.497 ms | 7.444 ms |

joint 主峰相对 kdomain 偏移一个 IFFT bin。PDP/LFM 形状通过门槛，但 99% 能量窗逼近 `Tmax`，说明 F=32 仍是验证配置，而不是最终通信 CIR 配置。进入通信统计生成前应扩大无模糊时延窗（减小 `df`）并保持 4 kHz 总带宽；仅零填充不能解决混叠。

## 7. 伪协方差与分布

64 样本测试集观测到：

- kdomain `||P||_F/||C||_F = 0.2696`；
- joint `||P||_F/||C||_F = 0.2558`。

这两个值相近，而边界解析模型在同一 U=5 条件下的理论比值约为 `2.09e-10`。对于 F=32、L=64，高维样本伪协方差即使在真正 proper 的零假设下也不会为零，因此不能用 `0.26` 直接证明接收端 improper，也不能直接忽略 `P`。本阶段保留增广协方差生成能力，最终决定依赖已新增的 properness-null Monte Carlo、更多独立测试样本或解析传播后的 `A P A^T` 检查。

实部、虚部、幅值、相位、偏度、峰度和中心频点 QQ 图已保存在 MAT 与 PNG 输出中。joint-kstat 本身是满足二阶统计的复高斯场；显式同一海面经过非线性指数相位屏后可保留高阶非高斯特征，因此二阶统计生成器仍需用 QQ/峰度/通信性能检验其损失。

## 8. raw PM 映射

192 个独立显式 PM 海面给出的中央裁剪能量比为：

- 均值 `0.99237`，即平均有符号偏差 `-0.763%`；
- 95% 置信区间 `[0.97671, 1.00803]`；
- 置信区间半宽 `1.566%`。

因此 U=5 下 100 m/256² PM 与 50 m/128² PE 可按同 `dx` 中央裁剪解耦，且本轮置信区间满足预设 2% 稳定性判据。该结论不能直接外推到高风速；U=12/15 仍应使用更大 PM 孔径或低/高 K 分解，不应扩大完整 PE 网格来追逐低波数峰。

## 9. 实际性能

测试为真实 wall time，不是 FFT 次数投影。

| 项目 | 时间/内存 |
|---|---:|
| cached PE 固定项构建 | 10.92 s |
| joint K/-K 因子构建 | 32.85 s |
| cached PE 数组 | 20.0 MiB |
| joint 因子与 independent 谱 | 360.0 MiB |
| MATLAB 峰值内存快照 | 3.018 GiB |
| 公共 kdomain 单 realization，F=32 | 43.71 s |
| cached surface-to-rx 单 realization（边界已给定） | 0.233 s |
| 同口径传播加速 | 187.8x |

joint 随机场生成加 cached PE 的累计训练实测：

| L | 总时间 | 每 realization |
|---:|---:|---:|
| 16 | 9.527 s | 0.595 s |
| 64 | 37.222 s | 0.582 s |
| 128 | 74.141 s | 0.579 s |

固定 L=16 的 batch 实测：batch 1/2/4/8 总时间分别为 `11.406/10.296/9.130/7.608 s`，即每样本约 `0.713/0.643/0.571/0.476 s`。当前机器上 batch 8 最快，但内存随 batch 增长；默认 batch 4 是速度与内存的保守折中。

缓存可复用的量已经实际落地：tx-to-surface 入射场、direct、coherent reflection、surface-to-rx diffraction factor、步进数、吸收屏、介质项、PM/PE 映射。每条样本只执行随机边界生成和 surface-to-receiver 推进。

## 10. 通过条件审查

| 条件 | 状态 |
|---|---|
| joint covariance 明显优于 independent | 通过 |
| reflected-only PDP correlation >= 0.9 | 通过，0.9714 |
| LFM correlation >= 0.9 | 通过，0.9917 |
| component sum <= 1e-10 | 通过，2.78e-17 |
| PM 映射能量误差有稳定 CI | 通过，CI 半宽 1.566% |
| cached/public 同输入一致 | 通过，double 误差 7.04e-16 |
| cached 实际加速明确 | 通过，传播 187.8x；含 joint 生成约 75x |
| P 是否可忽略有接收端依据 | **未完成**；已有观测，缺有限样本零假设校准 |

因此当前状态是“需要一个小范围统计校准后进入 `mu,C,P` 生成器”，不是需要重构 PE 或跨频模型。主要剩余阻碍是伪协方差的有限样本判定，其次是 F=32 的时延窗接近散射尾部长度。

## 11. 下一阶段最小实现计划

1. 先运行 `analyze_receiver_properness_null_vertical.m`，以 L=64 和训练协方差做至少 2000 次 proper-null Monte Carlo；验收为报告 observed P/C、95% null interval 和 upper-tail p-value。
2. 将独立 kdomain 测试集扩到 128 或 256，但不重建多风速库；检查 `epsilon_C` 是否随样本数下降，以及 joint-kdomain 误差是否仍接近 kdomain split-sample floor。
3. 保持带宽 4 kHz，增加 F 到 64，令 `Tmax` 约翻倍；验收为 99% 散射尾部不再贴近无模糊窗边界，且 joint PDP/LFM correlation 仍大于 0.9。
4. 建立单一 U=5 接收端生成器：保存 `H_dir`、`H_ref,coh`、散射 `mu_sca,C_sca,P_sca`、特征值谱、频率轴、参考时延和训练 seed；使用未参与估计的 kdomain 测试集验证。
5. 若 properness 不能拒绝且 upper confidence bound 足够小，可提供 proper 快速路径；否则使用增广协方差生成 `[H;conj(H)]`，不得静默丢弃 `P`。

在以上 1--3 完成前，不建立多风速统计库，不修改公共默认值，不接入通信主链。

