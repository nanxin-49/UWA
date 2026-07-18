# U=5 m/s 条件统计信道生成器 F=64 验证报告

日期：2026-07-12

## 1. 最终结论

本阶段已完成固定 U=5 m/s、raw PM、4--8 kHz 条件下的接收端统计信道生成器。新增功能均为独立模块；公共 `vertical_channel_model`、`vertical_wape_propagator`、默认海面模型和通信链路未修改。

结论如下：

1. joint-kstat 接收端散射项按任务规定的 central 95% null interval 规则不能拒绝 proper complex Gaussian 假设；允许使用 proper 快速路径，同时保留 augmented-real improper 路径。
2. F=64 将 `Tmax` 扩大到 15.75 ms。显式 kdomain 的最短圆周连续 99% 能量窗为 4.657 ms，即 `0.296*Tmax`，解决了 F=32 时的时延窗疑虑。
3. joint-kstat 相对 held-out kdomain 的 reflected-only PDP/LFM correlation 为 0.9525/0.9839；independent 仅为 0.0048/0.4104。
4. 统计生成器 full-rank 模型相对 held-out kdomain 的 covariance 误差为 0.3922，PDP/LFM correlation 为 0.9548/0.9858。显式 kdomain 训练/测试 split 自身的 covariance 差异为 0.2800，因此剩余误差包含有限样本误差和 joint-to-kdomain 模型差异。
5. held-out 综合选择为 full rank（数值秩34）。99.9%/99% 模型仅需7/5阶，PDP/LFM几乎不下降，但 covariance 误差略高；当前正式推荐 full rank，低秩作为可选压缩路径。
6. 10,000 条 H+CIR 在无 PE 条件下冷启动生成0.178 s，MAT落盘1.007 s，总计1.185 s；文件28.23 MiB。

该单条件模型已经达到进入“多风速条件库原型”的技术条件，但多风速之前应先优化 F=64 空间联合因子构建的时间和峰值内存，并在每个风速节点重复 properness 检验。

## 2. 新增独立接口

- `estimate_conditional_channel_stats_vertical.m`：估计 `mu_scatter_f`、`C_scatter_f`、`P_scatter_f`、复协方差EVD、实增广协方差EVD、负特征值能量、累计能量和候选秩。
- `sample_conditional_channel_vertical.m`：支持 `proper`、`improper` 和 `auto` 路径；可选择 full、99.9%、99% 或已选秩；按固定分量重构总信道。
- `build_physical_cir_vertical.m`：检查等间隔频率轴、应用有符号参考时延和窗、执行IFFT，记录真实分辨率、插值间隔和无模糊时延。
- `properness_null_test_vertical.m`：从训练协方差生成 proper 零假设样本，输出均值、中位数、90/95/99区间、单侧上界、p-value 和决策。
- `generate_cached_kdomain_ensemble_vertical.m`、`generate_cached_kstat_ensemble_vertical.m`：可复用的显式海面及 joint/independent cached PE ensemble 接口。
- `validate_u5_conditional_channel_generator_vertical.m`：支持 `smoke` 和 `full` 两级验证及阶段缓存。
- `generate_u5_f64_sample_bundle_vertical.m`：保存10,000条带完整 H/CIR/metadata 的样本包。

模型保存字段包括：条件元数据、频率轴、参考时延、`H_direct_f`、`H_ref_coh_f`、`mu/C/P`、复与实增广特征分解、候选秩、properness、训练seed、验证和计时。噪声未写入 H 或 h。

## 3. 固定配置

| 项目 | 值 |
|---|---:|
| 风速 | 5 m/s |
| 海况/表面 | `raw_pm` / joint-frequency Kirchhoff-kstat |
| 频率 | 4--8 kHz，F=64 |
| `df` | 63.4921 Hz |
| `1/B` | 0.25 ms |
| `Tmax=1/df` | 15.75 ms |
| PM网格 | 100 m / 256² |
| PE网格 | 50 m / 128² |
| 映射 | 同dx中央裁剪 |
| 训练/测试 | 128 / 64 |
| 统计生成 | 10,000 |
| 气泡/多普勒 | off / off |

192个显式海面中央裁剪能量比均值为0.9913，95% CI `[0.9738,1.0088]`。raw-PM独立审计中U=5、100m/256²离散/无限域方差比为0.99833。

## 4. Properness-null Monte Carlo

每种模式使用训练协方差和 `Ltest=64`，执行2000次 proper complex Gaussian Monte Carlo。

| 模式 | observed P/C | null mean | null median | central 95% interval | upper-tail p | 规定规则结论 |
|---|---:|---:|---:|---:|---:|---|
| kdomain | 0.3327 | 0.2250 | 0.2201 | [0.1304,0.3473] | 0.0405 | 不拒绝 proper；单侧5%边缘警告 |
| joint | 0.1911 | 0.2285 | 0.2233 | [0.1364,0.3476] | 0.7326 | 不拒绝 proper |
| independent | 0.6991 | 0.6250 | 0.6252 | [0.6000,0.6491] | 0.00050 | 超过99%上界，拒绝 proper |

joint 的90%、95%、99% central intervals分别为 `[0.1478,0.3253]`、`[0.1364,0.3476]`、`[0.1183,0.3927]`。kdomain虽按任务规定的central interval规则不拒绝，但单侧p=0.0405，说明显式非线性相位屏可能存在轻微impropriety，也可能仍是样本波动。生产模型默认 proper，但不能删除 improper 路径。

## 5. F=64 接收端三方验证

全部指标只使用 `H_ref_sca(f)`。

| 指标 | joint | independent |
|---|---:|---:|
| covariance relative error | 0.4159 | 0.9956 |
| correlation-matrix relative error | 0.4119 | 1.0000 |
| adjacent-correlation RMSE | 0.0232 | 1.0236 |
| PDP correlation | 0.9525 | 0.0048 |
| LFM correlation | 0.9839 | 0.4104 |

joint在接收端显著优于independent。covariance误差较F=32增大，主要因为维数从32增至64而测试样本仍为64；kdomain自身训练/测试差异已达0.2800。下一阶段若需要更精确的协方差比较，应优先增加测试样本，而不是改变跨频物理模型。

## 6. 物理 CIR 与时延窗

PE相位约定为正传播时延对应负频率相位斜率。本条件使用有符号参考时延

`tau_ref_signed=-(z_tx+z_rx)/c0`，

使任务公式 `exp(-i2*pi*f*tau_ref)` 移除公共传播斜率。

矩形频率截断产生跨越IFFT周期边界的带限旁瓣。旧的“把峰移到第一个bin后线性累积99%”会把峰前旁瓣错误记为接近`Tmax`的长尾。本阶段改用“包含99%总能量的最短圆周连续区间”，并以主峰为中心计算有符号mean/RMS。该修正不使用零填充，也不改变 H(f)。

显式kdomain结果：

- 主峰相对位置：9.844 ms；
- 主峰中心有符号mean delay：-0.005 ms；
- RMS delay spread：0.526 ms；
- 最短圆周99%能量窗：4.657 ms；
- `tail99/Tmax=0.2957`。

因此F=64满足建议的 `<0.8*Tmax`。零填充仍只被允许作为显示插值。

## 7. 条件统计模型与低秩选择

模型由128条joint-kstat+cached PE训练样本估计，协方差收缩参数为0。Hermitian化后进行EVD，小负特征值裁零并记录负能量。当前复协方差数值full rank为34；99.9%与99%累计能量秩为7和5。

| 模型 | 秩 | C误差 vs kdomain | C误差 vs joint test | PDP corr | LFM corr | 幅值KS |
|---|---:|---:|---:|---:|---:|---:|
| full | 34 | **0.3922** | 0.2670 | **0.9548** | **0.9858** | 0.1251 |
| 99.9% | 7 | 0.4000 | 0.2669 | 0.9501 | 0.9846 | **0.1229** |
| 99% | 5 | 0.3992 | **0.2636** | 0.9533 | 0.9848 | 0.1235 |

三种候选均通过PDP/LFM门槛。按脚本预先定义的“通过时域门槛后最小化held-out kdomain covariance误差”规则，选择full rank。99%模型仅5阶且时域性能相当，可在多风速库文件大小或抽样速度成为瓶颈时作为压缩候选，但不替代当前full基准。

中心频点full模型KS距离为：实部0.115、虚部0.117、幅值0.125。生成器K factor为0.0037，held-out kdomain为0.0205；散射均值接近零，K估计对64条测试样本较敏感。均值绝对误差为0.00735，频点mean/standard-error中位数1.48、最大2.44，未显示系统性大偏移，但后续需增加测试样本确认。

## 8. 性能、内存与break-even

| 项目 | 实测 |
|---|---:|
| F=64 fixed PE cache构建 | 32.80 s |
| joint K/-K因子构建 | 631.33 s |
| 128条joint训练ensemble | 176.20 s |
| 接收端统计估计 | 0.023 s |
| 总建模时间 | 840.35 s |
| 公共kdomain+PE单条F=64 | 578.00 s |
| cached joint-kstat+PE单条 | 1.3766 s |
| 峰值MATLAB内存快照 | 6.626 GiB |
| F=64空间joint/PE缓存文件 | 815.1 MB |
| 最终接收端模型文件 | 694 KB |

热运行批量H生成时间：100/1000/10000条分别为0.00571/0.00587/0.01345 s；10,000条仅H的摊销约1.35微秒/条。独立冷启动样本包测试中，10,000条H+CIR为0.1776 s，即17.8微秒/条；28.23 MiB MAT落盘另需1.007 s，总计1.185 s。全部统计生成过程PE调用数为0。

按公共全路径比较：

`M_break-even=840.35/(578.00-Tsample)=1.45`，即约2条信道后回本。

按已有cached joint+PE比较：

`M_break-even≈840.35/(1.3766-Tsample)=610`。

第一种衡量替代公共逐条PE的收益；第二种衡量相对已优化训练路径的收益。空间联合因子构建和6.6GiB峰值内存是进入多风速库前的主要工程风险；最终接收端模型本身很小。

## 9. 回归结果

- cached/public double同输入：总响应最大误差 `7.04e-16`，F=32公共/cached时间47.46/0.247 s。
- reduced direct-only：`H_f-H_direct_f=0`。
- reduced direct-plus-reflect：分量和误差0。
- independent kstat：固定seed重复误差0；换seed反射变化0.0862；公共字段完整。
- 原joint边界回归，32²/16 realizations：joint/independent covariance误差0.0288/0.9523，PDP 0.9981/0.3991，LFM 0.9973/0.4076。
- raw-PM网格审计：U=5、100m/256² capture 0.99833；32-seed中央裁剪平均能量相对误差2.04%。
- target_hs/kirchhoff_spatial由reduced flat reflected回归覆盖；公共传播和通信噪声位置未改。

## 10. 验收与下一阶段

本阶段十项通过标准全部具备实现或实测证据：properness有结论；F=64尾部不贴窗；joint与生成器PDP/LFM通过；proper/improper均运行；full/99.9/99有held-out选秩；10,000条H/CIR无PE；计时、内存和break-even完整。

建议进入多风速条件库的最小下一步：

1. 只增加一个风速节点U=8 m/s，保持F=64、同一PE网格和接收端模型schema，验证流程可复用性。
2. 为F=64 joint空间因子实现K域分块落盘或频率低秩构建，目标峰值内存低于4 GiB、构建时间低于当前631 s。
3. 每个风速节点独立执行properness-null；若显式kdomain超过central 99%上界，则切换到augmented-real抽样。
4. 将kdomain测试样本增至128，降低64维协方差比较的采样噪声；训练样本暂保持128，必要时再增至256。
5. 多风速第一版仍只使用离散节点，不做协方差插值，不接入正式PSK/OFDM链路。

