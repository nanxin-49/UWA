# 跨频 Kirchhoff-kstat 与 raw-PM 网格前置验证报告

日期：2026-07-11

## 结论摘要

本阶段完成了两个独立前置能力：

1. 可重复的 raw-PM 孔径、网格和映射审计；
2. 基于解析跨频协方差和伪协方差的 joint-frequency Kirchhoff-kstat 边界相位屏原型。

没有修改 `vertical_channel_model`、`vertical_wape_propagator`、`pm_surface_boundary_model`、公共默认配置或通信链路。当前 `independent` kstat 分支仍保持原状，新 joint 模型只通过独立测试接口运行。

对任务所列七个问题的直接回答：

1. **边界层可以逼近。** U=5 m/s、4--8 kHz、F=32、64 realization、64² 边界网格下，joint kstat 对同一显式高斯海面的跨频协方差相对误差为 0.67%，邻频相关 RMSE 为 `2.15e-4`；independent 对应为 95.26% 和 0.9899。
2. **U=5 节点不需要显著的非零伪协方差，但生成公式必须保留。** 解析 `||P||F/||C||F=2.09e-10`；显式和 joint 有限样本值分别为 0.0027 和 0.0034。此结论不能外推到其他 U、频带或较弱相位起伏。
3. **raw PM 网格与 PE 网格可以解耦。** 同采样间隔的大 PM 网格中央裁剪到 PE 窗口无需插值；U=5 的 100 m/256² PM 网格映射到 50 m/128² PE 网格时，32 seed 平均能量误差为 2.04%。但 joint 模型在大 PM 网格上的 `F² Nxy` 存储仍需低秩/分块优化。
4. **U=5 推荐 PM=100 m/256²，PE=50 m/128²（reduced 原型）。** 若 PE 使用 50 m/256²，则保持相同 `dx` 的 PM 网格应为 100 m/512²。
5. **高风速不能继续用 50 m 孔径。** 固定 `dx=50/128 m` 的建议起点为：U=8 用 100 m/256²，U=10 用 200 m/512²，U=12 用 300 m/768²，U=15 用 400 m/1024²。U=12/15 更适合进一步采用低/高波数分解，而不是扩大完整 PE 网格。
6. **可缓存项明确。** 发射端到海面的入射场、直达响应、相干反射响应、传播核和固定介质/吸收项可缓存；每 realization 的随机散射场及 surface-to-receiver 推进不可直接复用。
7. **具备进入接收端统计生成阶段的数学前提，但尚未具备生产条件。** 下一步可以做单 U 的 cached-PE 接收端 `mu,C,P` 原型；在多风速建库前还需实现缓存执行路径、joint 模型的大网格低秩/分块版本，并完成 PE 后 held-out 验证。

## 新增实现

### raw PM 接口

- `raw_pm_spectrum_grid_vertical.m`
  - 复现项目 PM 波数谱；
  - 输出 `K_min`, `K_peak`, `K_nyquist`, 离散/无限域方差、隐含 `Hs`、离散方差比和理想径向支持比；
  - 明确指出旧项目参数 `U` 尚未声明测风高度/平均约定。
- `sample_raw_pm_surface_vertical.m`
  - 复现当前显式海面 `sqrt(2)*real(ifft2(...))` 约定；
  - 输出单 realization 能量误差。
- `scripts/validation/validate_raw_pm_grid_coverage_vertical.m`
  - 扫描 U=`[3,5,8,10,12,15]`；
  - 扫描 50--400 m 孔径及固定孔径网格分辨率；
  - 审计 100 m PM 网格到 50 m PE 网格的中央裁剪映射。

### joint-frequency kstat 接口

- `sample_kirchhoff_kstat_joint_frequency_vertical.m`
  - `mode='joint'|'independent'`；
  - joint 模式使用解析 `C_Gij(rho)` 和 `P_Gij(rho)`；
  - 对每个 `K,-K` 谱对构造增广协方差并 EVD 抽样；
  - 自共轭 FFT bin 使用 2F 维实增广协方差；
  - 记录负特征值裁剪、FFT 数、内存和时间；
  - independent 模式只保留每频对角协方差，作为原工程模型对照。
- `scripts/validation/validate_kstat_joint_frequency_boundary_vertical.m`
  - 同一高斯 PM 海面的显式 Kirchhoff 相位屏；
  - independent kstat；
  - covariance+pseudocovariance joint kstat；
  - 输出相干项、非相干谱、相关矩阵、伪协方差、PDP 和 LFM 匹配滤波代理。

### PE 缓存审计

- `scripts/validation/audit_pe_caching_vertical.m`
  - 逐频计算步数与 FFT 次数；
  - 给出 L=`[1,16,64,128]` 的当前/缓存后成本；
  - 测量 128² 公共路径基线；
  - 给出 batch=8 的缓存内存预算。

## 1. raw PM 网格覆盖

### 1.1 公式

代码采用

\[
E(K)=\frac{\alpha}{2K^3}
\exp\left[-\frac{\beta g^2}{U^4K^2}\right],
\qquad
\Phi_{2D}(K_x,K_y)=\frac{E(K)}{2\pi K}.
\]

径向谱峰值和无限域方差为

\[
K_{\rm peak}=\sqrt{\frac{2\beta}{3}}\frac{g}{U^2},
\qquad
\sigma_{\eta,\infty}^2=\frac{\alpha U^4}{4\beta g^2}.
\]

审计同时输出两个容易混淆的量：

- `capture_ratio_discrete_to_infinite`：二维矩形网格 Riemann 和除以解析总方差，混合了支持截断和离散求积误差；
- `capture_ratio_radial_support_idealized`：仅按 `K_min--K_max` 连续径向积分得到的理想支持覆盖，用于识别“离散求和偶然补偿漏谱”。

因此不能只用一个接近 1 的离散方差比判断网格已经解析谱峰。

### 1.2 固定 dx 的孔径扫描

固定 `dx=50/128=0.390625 m`，代表只增加孔径和低波数分辨率而保持相同高波数 Nyquist：

| U (m/s) | 推荐审计节点 | `Kmin` | `Kpeak` | 离散/无限方差比 | 理想径向支持比 | `Hs_implied` (m) |
| ---: | --- | ---: | ---: | ---: | ---: | ---: |
| 3 | 50 m / 128² | 0.1257 | 0.7656 | 0.9889 | 0.9932 | 0.1909 |
| 5 | 100 m / 256² | 0.06283 | 0.2756 | 0.9983 | 0.9991 | 0.5328 |
| 8 | 100 m / 256² | 0.06283 | 0.1077 | 1.0015 | 0.9876 | 1.3661 |
| 10 | 200 m / 512² | 0.03142 | 0.06890 | 0.9919 | 0.9992 | 2.1243 |
| 12 | 300 m / 768² | 0.02094 | 0.04785 | 0.9933 | 0.9996 | 3.0613 |
| 15 | 400 m / 1024² | 0.01571 | 0.03062 | 0.9913 | 0.9967 | 4.7782 |

U=8 的 100 m 节点仍只有约 1.7 个 bin 到峰值，虽然总能量口径已经接近；正式库建议再检查 150/200 m。U=10 的 100 m 节点离散比为 1.062，说明粗离散发生 6.2% 过积分，不能接受。

### 1.3 PM 到 PE 的映射

本次实现的最小映射是：

1. PM 统计网格：100 m/256²；
2. PE 网格：50 m/128²；
3. 两者 `dx` 完全相同；
4. 截取大网格中央 128² 样本，不进行插值，不对局部窗口去均值。

32 seed 结果：

- PM 网格目标方差：约 `0.0177 m²`；
- 映射窗口平均能量：约 `0.0181 m²`；
- 平均能量相对误差：2.04%；
- 若对 PE 窗口去均值，额外损失约 0.29% 能量。

2.04% 与有限 32 seed 抽样误差同量级，说明该映射可作为 U=5 原型，但生产门槛应增加 seed 数并给出置信区间。绝对海面偏置不应随意删除，因为大尺度低 K 模式在较小 PE 窗口内本来就可能表现为缓变偏置/斜率。

### 1.4 三种网格方案判断

#### 直接扩大 PE 网格

物理上直接，代码改动最小，但不推荐作为主方案。二维 FFT 成本约为 `N² log N`：

- 128² 到 256²：单 FFT 约增加 4--5 倍；
- 128² 到 768²：约增加 40--50 倍；
- 128² 到 1024²：约增加 70--80 倍，复数组内存增加 64 倍。

U=12/15 若扩大完整 PE 网格，会把低 K 海面统计问题变成传播计算问题，代价过高。

#### 大孔径 PM 网格映射到 PE 网格

推荐 U=5--10 的第一选择。优点是 PE 网格和传播成本不变；同 `dx` 中央裁剪无需插值。限制是当前 full joint kstat 为每个 K 保存 `F x F` 谱矩阵，其内存按 `F²Nxy` 增长。F=32、256² 时，仅 covariance+pseudocovariance 单精度临时谱约 1.07 GB，尚未包含输出和 EVD 工作区。

#### 低/高波数分解

推荐 U=12/15 的研究路线。将高斯海面分为独立频带

\[
\eta=\eta_{low}+\eta_{high},\qquad
C_\eta=C_{low}+C_{high}.
\]

低 K 在大孔径粗网格生成并投影到 PE 坐标，高 K 在 PE 网格生成。必须先在高度域合并，或在 kstat 公式中先合并 `C_eta` 再计算指数；不能把两个相位屏的 `S_deltaG` 简单相加，因为 `exp(i alpha eta)` 是非线性的。

## 2. joint-frequency kstat 推导与实现

### 2.1 协方差

令

\[
G_i(\mathbf r)=R_{0,i}\exp(i\alpha_i\eta(\mathbf r)),
\qquad
\mu_i=R_{0,i}\exp(-\alpha_i^2\sigma_\eta^2/2),
\]

则零均值相位屏 `deltaG_i=G_i-mu_i` 的跨频协方差为

\[
C_{G,ij}(\rho)=R_{0,i}R_{0,j}^*
e^{-(\alpha_i^2+\alpha_j^2)\sigma_\eta^2/2}
\left[e^{\alpha_i\alpha_jC_\eta(\rho)}-1\right].
\]

### 2.2 伪协方差

对应的伪协方差为

\[
P_{G,ij}(\rho)=R_{0,i}R_{0,j}
e^{-(\alpha_i^2+\alpha_j^2)\sigma_\eta^2/2}
\left[e^{-\alpha_i\alpha_jC_\eta(\rho)}-1\right].
\]

负号来自 `E[exp(i alpha_i eta_1+i alpha_j eta_2)]`。省略该项等价于预先假设 proper/circular，而不是从海面模型推导。

### 2.3 K 域联合抽样

平稳场的伪协方差使 `K` 与 `-K` 耦合。对非自共轭谱对构造

\[
\mathbf y=
\begin{bmatrix}
\mathbf X(\mathbf K)\\
\mathbf X^*(-\mathbf K)
\end{bmatrix},
\quad
E[\mathbf y\mathbf y^H]=
\begin{bmatrix}
\mathbf S_C(\mathbf K)&\mathbf S_P(\mathbf K)\\
\mathbf S_P^H(\mathbf K)&\mathbf S_C^*(-\mathbf K)
\end{bmatrix}.
\]

每个增广矩阵 Hermitian 化、EVD、裁剪数值负特征值后抽样。自共轭 bin 使用 `[Re X; Im X]` 的 2F 维实协方差。64²/F=32 验证中，被裁负特征值绝对和/正特征值和为 `1.86e-8`，属于周期 FFT 和单精度临时谱产生的小数值残差。

## 3. 边界层验证结果

### 3.1 设置

- `U=5 m/s`；
- `raw_pm`；
- `f=linspace(4000,8000,32)` Hz；
- 边界网格 64²、孔径 50 m；
- 64 independent realizations；
- 显式参考为同一 `eta(x,y)` 用于全部频率；
- 边界接收代理为反射屏的 aperture mean；
- LFM/PDP 是边界层代理，不包含 PE 传播。

该设置用于证明联合二阶统计构造，不是高风速/大网格生产验证。

### 3.2 数值结果

| 指标 | independent kstat | joint kstat |
| --- | ---: | ---: |
| 跨频 covariance 相对误差 | 0.9526 | **0.0067** |
| 邻频相关 RMSE | 0.9899 | **2.15e-4** |
| PDP correlation | 0.4890 | **0.9837** |
| LFM matched-filter correlation | 0.5163 | **0.9842** |
| 非相干谱相对误差 | 0.1755 | 0.1781 |
| coherent max absolute error | 0.0055 | 0.0033 |

两种 kstat 的单频非相干谱误差相近是预期结果：joint 改变跨频联合结构，不改变单频对角功率目标。17.8% 主要是 64 realization 的二维谱抽样误差，需通过 M 收敛单独降低。

joint 显著改善 PDP 和匹配滤波，但尚未达到“接收端真实信道验证”，因为这里没有 incident-field weighting、surface-to-receiver PE 或接收深度效应。

### 3.3 伪协方差判断

- 解析 zero-lag `||P||F/||C||F = 2.09e-10`；
- 显式样本值 0.0027；
- joint 样本值 0.0034；
- `P` 的相对误差大于 1，是因为分母对应的真实 `P` 几乎为零，不应作为失败指标。

U=5 raw-PM 在 4--8 kHz 的相位方差较大，反射相位充分缠绕，二阶统计近似 proper。实现仍保留 `P`，因为低 U、低频、较小粗糙度或非压力释放复反射系数下它可能不可忽略。

### 3.4 运行成本

64²/F=32/M=64：

- 显式同一海面 ensemble：0.751 s；
- joint 谱构建：0.160 s；
- joint 抽样和 IFFT：2.063 s；
- joint 输出数组：134.2 MB；
- joint covariance+pseudocovariance 临时谱：67.1 MB。

该时间不含 PE。

## 4. PE 重复计算与缓存

### 4.1 可缓存性

| 量 | 可缓存 | 条件 |
| --- | --- | --- |
| 发射端到海面入射场 | 是 | 固定环境、源、频率、网格 |
| 直达信道 | 是 | 可在一次 tx-to-surface 推进经过 z_rx 时提取 |
| 相干反射场/响应 | 是 | 固定 U、环境和频率 |
| PE 传播核 `fr` | 是 | 固定频率、网格、dz、介质 |
| 网格、吸收层、固定介质屏 | 是 | 固定环境和频率 |
| 随机散射边界场 | 否 | realization 改变 |
| 随机 surface-to-receiver 推进 | 否 | 输入场改变；可按 realization batch 并行 |

气泡或时变介质若随 realization 改变，对应传播核/屏不再可缓存；若气泡参数只是环境节点的一部分，则可在该节点内部缓存。

### 4.2 步进和 FFT 计数

固定几何 z_tx=100 m、z_rx=3 m，F=32：

| L | 当前 PE steps | 缓存后 PE steps | step 降幅 | 当前 FFT 数 | 缓存后 FFT 数 | FFT 降幅 |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 51,247 | 27,182 | 1.89x | 102,910 | 55,644 | 1.85x |
| 16 | 819,950 | 38,927 | 21.1x | 1,646,600 | 80,574 | 20.4x |
| 64 | 3,279,800 | 76,511 | **42.9x** | 6,586,200 | 160,350 | **41.1x** |
| 128 | 6,559,600 | 126,620 | 51.8x | 13,172,000 | 266,720 | 49.4x |

缓存后固定成本包括：一次 tx-to-surface（途中提取 direct）、一次 coherent surface-to-rx，以及 joint `C/P` 谱构建；每 realization 只保留随机 surface-to-rx。

### 4.3 时间

当前未改动公共路径在 128²、4 个频点上的实测平均为 5.9005 s/realization，即约 1.475 s/频点。在频点近似线性缩放下：

- F=32 单 realization 当前路径约 47.2 s；
- L=64 约 3021 s，约 50.4 min。

按本机 128² FFT primitive 计时，L=64：

- 当前 FFT 工作投影 1456.8 s；
- 缓存后 FFT 工作投影 35.5 s。

实际公共路径还包含介质屏、数组分配、边界统计和 MATLAB 调度，当前实测约为纯 FFT 投影的两倍。由此缓存后端到端时间的保守预测约 35--75 s/64 条，加上 joint 模型构建；这不是实测优化时间。要给出正式加速比，下一阶段必须实现独立 cached executor 并与当前 public path 同机计时。

### 4.4 内存

128²/F=32、complex double：

- 入射场缓存：8.39 MB；
- 传播核：8.39 MB；
- 固定介质屏：8.39 MB；
- L=64 全联合边界场：536.9 MB；
- realization batch=8：67.1 MB；
- 推荐缓存和 batch 工作集合计约 92.3 MB，不含 MATLAB header/临时 FFT workspace。

因此应流式分批传播，不能把所有 L 条二维宽带边界场同时保留。

## 5. 输出文件

主要结果：

- `results/validation/validate_raw_pm_grid_coverage_vertical_result.mat`
- `results/validation/validate_raw_pm_grid_coverage_vertical_table.csv`
- `results/validation/validate_kstat_joint_frequency_boundary_vertical_result.mat`
- `results/validation/audit_pe_caching_vertical_result.mat`
- `results/validation/audit_pe_caching_vertical_cost_table.csv`

图：

- `kstat_joint_frequency_correlation_compare.png`
- `kstat_joint_frequency_adjacent_correlation.png`
- `kstat_joint_frequency_pseudo_covariance.png`
- `kstat_joint_frequency_pdp_compare.png`
- `kstat_joint_frequency_matched_filter_compare.png`
- `kstat_joint_frequency_coherent_reflection.png`

## 6. 验证与回归

已运行：

1. raw-PM 全孔径/网格扫描及 32-seed 映射审计；
2. joint 最小 4²/2-frequency 和 16²/4-frequency 数值测试；
3. U=5、F=32、64²、64-realization 边界统计验证；
4. 128²/4-frequency 当前 kstat 公共路径两次计时；
5. reduced 通信回归：256²、F=4、128 symbols，`direct_only` 和 `direct_plus_reflect` 均完成，现有通信链未改变。

公共路径回归仍满足：

- direct-only 的 `H_reflect_f=0`；
- direct-plus-reflect 正常运行；
- 噪声仍由通信接收端注入；
- 公共默认仍为 `target_hs + kirchhoff_spatial`。

## 7. 限制与进入下一阶段的条件

已满足：

- joint 公式和 K-domain 增广抽样正确工作；
- U=5 边界跨频 covariance/PDP/LFM 明显优于 independent；
- raw PM/PE 网格解耦有可执行映射；
- PE 重复项和理论加速空间已量化。

尚未满足：

- joint 边界场尚未接入 cached surface-to-receiver PE；
- 尚未验证 PE 后接收端 `H_ref,sca(f)` 的 covariance 和 pseudo-covariance；
- 大 PM 网格 full joint 的 `F²Nxy` 内存尚未降阶；
- 低/高 K 分解尚未实现；
- 优化后运行时间尚未实测；
- 高风速建议只通过谱积分审计，尚未通过 joint+PE。

因此下一阶段应是单 U=5 的 **cached joint-kstat + PE 接收端验证**，而不是直接建立多风速统计库。建议先实现：

1. 独立 cached executor，不改公共传播器；
2. PM=100 m/256² 到 PE=50 m/128² 的分批 joint 边界映射；
3. PE 后保存 `H_dir`, `H_ref,coh`, `H_ref,sca`, `H_total`；
4. 用未参与建模的同一海面 kdomain PE ensemble 比较接收端 `C/P`、PDP 和 LFM；
5. 实测当前/缓存路径的 wall time 和峰值内存。
