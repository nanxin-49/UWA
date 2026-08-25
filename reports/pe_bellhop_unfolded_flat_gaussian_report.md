# 展开坐标 PE--Bellhop 平面海面 Gaussian 交叉验证报告

更新时间：2026-08-25  
正式数据完成时间：2026-08-21 20:54  
验证入口：`scripts/validation/validate_pe_bellhop_unfolded_flat_gaussian_vertical.m`  
正式状态：`passed=false`（10 项通过，2 项未通过）

## 1. 执行摘要

本轮已经完成正式的 `984×984`、`4--8 kHz / 65 点` PE--Bellhop
平面海面对比。网络/会话中断没有终止后台 MATLAB；全部 PE、Bellhop、
空间剖面、波束收敛、PDP 和辅助到达计算均已结束，并保存了完整产物。

主要结论如下。

1. 生产 PE marching 与独立一步精确角谱传播在展开反射距离 `103 m` 上的
   最大相对误差为 `4.74e-13`。该结果继续支持“不修改 PE 平方根 marching
   算子”的决定。
2. 轴上接收机的宽带反射/直达比
   `Q(f)=H_reflect(f)/H_direct(f)` 与 Bellhop 高度一致：最大幅度差
   `4.94e-4 dB`，最大相位差 `1.83e-3 rad`，内部频点最大群时延差
   `62.4 ns`，PDP 相对峰时延完全一致。
3. 相对总信道 `1+Q(f)` 的最大 TL 差为 `0.00808 dB`，最大相位差为
   `0.0296 rad`，仍在原定幅相门槛内。
4. Bellhop `5001/10001` 波束结果已经充分收敛；最大差只有
   `1.38e-6 dB` 和 `9.66e-6 rad`。当前剩余偏差不是波束数量不足造成的。
5. 横向剖面的幅度一致，但外侧相位未通过严格门槛。最差点位于
   `8 kHz、横向偏移 19.53125 m`：TL 差 `0.0158 dB`，但相位差
   `0.0612 rad > 0.05 rad`，复场相对误差 `0.0613 > 0.02`。

因此不能把本轮写成“全部指标通过”。更准确的表述是：

> PE 的传播实现、平面压力释放反射的轴向相对响应、路径时延和 4--8 kHz
> 宽带信道均获得强一致性证据；Bellhop 与 PE 的远轴横向复场仍存在约
> 0.05--0.06 rad 的系统相位差，需要单独定位源指向性映射、2D/3D 场定义
> 与射线束相位近似，不应据此修改 PE marching。

## 2. 验证目标与范围

本轮只验证最简单、可复现的均匀介质和平面压力释放海面。目标是判断：

- PE marching 是否正确实现了当前精确平方根谱传播；
- 平面海面反射系数 `-1` 是否给出正确的直达/反射相对幅相和时延；
- 在同一 Gaussian 发射声场下，PE 与 Bellhop 是否给出一致的相对复声场；
- 剩余差异来自 PE、Bellhop 数值设置，还是两者的源/维度映射。

本轮不覆盖粗糙海面、随机海面、Kirchhoff 粗糙相位屏、SSA、气泡、
Doppler、海底作用、深度相关 SSP 或通信调制过程。验证层没有修改 PE
marching、FFT convention、Gaussian 源定义或生产默认窗口/sponge 配置。

## 3. 坐标展开与物理等价关系

原物理坐标中海面为 `z=0`，Tx 位于 `z_tx=100 m`，Rx 位于
`z_rx=3 m`。为了避免在 Bellhop 的“距离终点”人为放置边界，本轮采用
平面镜像展开：

\[
L_d=z_{tx}-z_{rx}=97\ \mathrm{m},
\qquad
L_r=z_{tx}+z_{rx}=103\ \mathrm{m}.
\]

Bellhop 的水平距离轴代表原模型的竖直传播方向，Bellhop 深度轴代表 PE
的一条横向截面。Bellhop 中只计算均匀自由场传播：

\[
H_d=P_{BH}(97\ \mathrm{m}),
\qquad
H_r=-P_{BH}(103\ \mathrm{m}).
\]

负号是压力释放平面海面的反射系数。Bellhop 上、下边界放在人工深度
`±1000 m`，并使用与水体匹配的半空间，避免额外边界反射。该展开只对
均匀介质和平面镜面边界成立，不能直接推广到粗糙海面或深度相关 SSP。

## 4. 正式计算配置

| 类别 | 设置 |
|---|---|
| 介质 | 均匀水体，`c=1500 m/s` |
| 物理位置 | Tx `(0,0,100 m)`；Rx `(0,0,3 m)` |
| 展开距离 | 直达 `97 m`；反射镜像 `103 m` |
| 发射场 | 生产 Gaussian，`sigma_src_m=0.3 m` |
| PE 网格 | `984×984`，宽度 `192.1875 m`，`dx=dy=0.1953125 m` |
| PE 纵向步长 | `stepz_lamb=0.5`；实际步长随频率变化 |
| 边界 | `sponge_ratio=0`，`alpha_max=0` |
| 频率 | `4000--8000 Hz`，65 点，间隔 `62.5 Hz` |
| 海面 | 平面压力释放，反射系数 `-1` |
| 关闭功能 | 粗糙面、随机散射、气泡、Doppler、通信处理 |
| Bellhop 正式波束 | `10001` |
| Bellhop 波束审计 | `5001/10001`，频率 `4/6/8 kHz` |
| Bellhop 步长 | `0.05 m` |
| Bellhop 角扇区 | `[-30,30] deg` |
| `.sbp` 采样 | 2401 个角度点，低于 `-120 dB` 截断 |
| 横向剖面 | `4/6/8 kHz`；偏移 `0, 1.953, 4.883, 9.766, 14.648, 19.531 m` |

Bellhop 的 Gaussian 指向性写为：

\[
D(\theta,f)=\cos\theta\,
\exp\left[-\frac{(k\sigma\sin\theta)^2}{2}\right],
\]

并在轴向归一化为 `0 dB`。验证不做距离或频率相关的逐点拟合。

## 5. 实现的验证链

### 5.1 Bellhop 自由场归一化

保存的独立自由场审计确认当前 Bellhop 输出满足 `1/R` 型归一化：

- 归一化均值误差：`1.27e-8`；
- 距离不变性误差：`8.32e-8`；
- `1/R` 斜率误差：`2.07e-8`；
- 空间相位约定误差：`1.97e-8 rad`；
- 到达时间误差：`3.43 ns`；
- Bellhop 到 Green 函数 `1/(4πR)` 的固定换算为
  `20log10(4π)=21.9842 dB`。

本次正式任务记录为 `audit_execution=validated_saved_fallback`。原因是当前
安装的旧 Bellhop 可执行文件在重新运行 ±180° 点源自由场审计时发生外部
heap corruption；此前已经通过的审计 MAT 被复用。正式 Gaussian 对比采用
相对量，未使用逐频或逐距离校准。该回退不影响相对 `Q(f)`，但新的
Bellhop 版本就绪后仍应重新执行 4/6/8 kHz 归一化审计。

### 5.2 PE 实现级检查

生产 PE 多步 marching 与验证层独立实现的一步精确离散角谱传播进行比较：

\[
P_{\Delta z}^{N}\Psi_0
\quad\leftrightarrow\quad
P_{103\,m}\Psi_0.
\]

反射端同时乘以 `-1`。65 个频点的最大相对误差为 `4.74e-13`，远低于
`1e-10` 门槛。这是实现级 hard check，检查 FFT 顺序、传播方向、步长累计
和载波处理；它不是两个独立物理模型之间的证据。

### 5.3 独立连续 Gaussian/Weyl 参考

轴向 PE 与独立连续 Weyl/Hankel Gaussian 角谱参考的反射/直达比最大差为：

- TL：`2.50e-4 dB`；
- phase：`6.31e-5 rad`。

这说明轴向 PE 结果不仅与自身离散 propagator 一致，也与独立连续角谱
自由场参考一致。

### 5.4 Bellhop 波束收敛

| 频率 | 最大 TL 差 | phase RMS |
|---:|---:|---:|
| 4 kHz | `8.20e-7 dB` | `4.84e-6 rad` |
| 6 kHz | `1.38e-6 dB` | `7.24e-6 rad` |
| 8 kHz | `1.19e-6 dB` | `9.66e-6 rad` |

全部远小于 `0.1 dB / 0.02 rad` 门槛。继续单纯增加波束数不太可能解释
当前横向 `0.061 rad` 的相位差。

### 5.5 轴向宽带相对信道

硬验收量定义为：

\[
Q(f)=\frac{H_{reflect}(f)}{H_{direct}(f)},
\qquad
H_{rel}(f)=1+Q(f).
\]

这种定义消除了固定源强归一化，但保留了反射/直达之间的真实幅度、相位
和路径时延关系。

| 指标 | 正式结果 | 门槛 | 状态 |
|---|---:|---:|:---:|
| `Q(f)` 最大 TL 差 | `4.94e-4 dB`（8 kHz） | `0.25 dB` | PASS |
| `Q(f)` 最大 phase 差 | `1.83e-3 rad`（8 kHz） | `0.05 rad` | PASS |
| 内部频点最大群时延差 | `62.4 ns`（4187.5 Hz） | `20 us` | PASS |
| PDP 峰时延差 | `0` | `0.125 ms` | PASS |
| `1+Q` 最大 TL 差 | `0.00808 dB`（7937.5 Hz） | `0.25 dB` | PASS |
| `1+Q` 最大 phase 差 | `0.0296 rad`（8 kHz） | `0.05 rad` | PASS |

PE 的 `|Q|` 位于 `0.941744--0.941801`，Bellhop 位于
`0.941748--0.941748`。两者的相对 PDP 峰均为当前相位约定下的 `-4 ms`；
其绝对值正好等于：

\[
\frac{L_r-L_d}{c}=\frac{103-97}{1500}=4\ \mathrm{ms}.
\]

### 5.6 横向复场剖面

横向剖面先分别以轴上值归一化，再比较 PE 与 Bellhop，因此表中差异不受
固定源强影响。最外侧 `19.53125 m` 的结果为：

| 频率 | direct TL 差 | reflect TL 差 | direct phase 差 | reflect phase 差 | 最大复场误差 |
|---:|---:|---:|---:|---:|---:|
| 4 kHz | `0.00142 dB` | `0.000837 dB` | `-0.00900 rad` | `-0.00751 rad` | `0.00900` |
| 6 kHz | `0.00549 dB` | `0.00440 dB` | `-0.0270 rad` | `-0.0228 rad` | `0.0270` |
| 8 kHz | `0.01582 dB` | `0.01274 dB` | `-0.06120 rad` | `-0.05167 rad` | `0.06128` |

误差随频率和横向偏移平滑增加：幅度仍非常一致，主要偏差表现为远轴相位
曲率差。轴上归一化点的误差为零，当前生产接收机又位于轴上，因此该失败
不会推翻轴向 `H(f)` 的通过结论；但它说明当前 Bellhop `.sbp` 映射尚未证明
能够重现 PE 的完整横向复声场。

### 5.7 原坐标辅助到达

保留的原坐标 `Tx=80 m、Rx=10 m、水平偏移 6 m、4 kHz` Bellhop arrival
案例识别出两条路径：

| 路径 | 时延 | 幅度 | 海面反射 | 海底反射 |
|---|---:|---:|---:|---:|
| 直达 | `46.838 ms` | `0.014234` | 0 | 0 |
| 一次海面反射 | `60.133 ms` | `0.011087` | 1 | 0 |

该案例只验证几何、反射相位和到达结构，不作为主幅度验收。

### 5.8 不变量与回归

- 65 个频点均满足
  `H_f=H_direct_f+H_reflect_f`，最大闭合误差 `3.88e-18`；
- 接收二维场中心与 `H_reflect_reduced_f` 的诊断一致性在此前缩减回归中通过；
- `5001/10001` Bellhop 波束检查通过；
- MATLAB `checkcode` 对新增验证入口无问题；
- 缩减通信回归使用 `256×256`、8/16 个频点、100 个符号、
  `Eb/N0={0,10} dB`，direct-only 与 direct-plus-reflect 均完成；
- 没有修改 `vertical_wape_propagator.m`、Gaussian 源定义或通信主线。

## 6. 自动检查汇总

| 检查 | 数值 | 门槛 | 结果 |
|---|---:|---:|:---:|
| Bellhop beam TL | `1.3815e-6 dB` | `0.1 dB` | PASS |
| Bellhop beam phase | `9.6636e-6 rad` | `0.02 rad` | PASS |
| 横向 TL | `0.015822 dB` | `0.25 dB` | PASS |
| 横向 phase | `0.061204 rad` | `0.05 rad` | **FAIL** |
| 横向复场 L2 | `0.061277` | `0.02` | **FAIL** |
| 宽带 `Q` TL | `4.9373e-4 dB` | `0.25 dB` | PASS |
| 宽带 `Q` phase | `0.0018271 rad` | `0.05 rad` | PASS |
| 群时延 | `6.2412e-8 s` | `2e-5 s` | PASS |
| PDP 峰时延 | `0 s` | `1.25e-4 s` | PASS |
| PE 展开恒等式 | `4.7395e-13` | `1e-10` | PASS |
| 原坐标辅助案例 | `0` | `0` | PASS |
| `H` 分量闭合 | `3.8790e-18` | `1e-12` | PASS |

程序按预定规则保存全部产物后以非零状态退出，没有放宽阈值。

## 7. 可以确认与不能确认的结论

### 已有充分证据支持

- 当前 PE 精确平方根谱 marching 的实现没有发现错误；
- 平面压力释放反射的 `-1` 符号、轴向相对幅度和额外 `4 ms` 路径时延正确；
- 对当前 Gaussian 源和轴上接收机，4--8 kHz 的 PE 与 Bellhop 相对信道
  `Q(f)` 高度一致；
- `192.1875 m + no sponge` 在本次平面展开验证中没有表现出此前小窗口的
  明显周期污染；
- 当前差异不是 Bellhop 波束数量不足，也没有证据要求修改 PE marching。

### 尚不能确认

- PE 与 Bellhop 在整个横向复场上达到 `2%` 以内的一致性；
- Bellhop `.sbp` 已严格等价于 3D PE 的径向 Gaussian 角谱；
- Bellhop 绝对压力能够在当前旧可执行文件上完成新的 4/6/8 kHz
  点源归一化审计；
- 结论能够推广到粗糙海面、非均匀 SSP、海底、气泡或实际换能器；
- 更换实际换能器后仍可直接沿用当前 `192 m` 窗口。

## 8. 剩余横向差异的判断

现有证据把问题范围收窄到了验证映射层，而不是 PE 核心：

1. PE--独立离散 AS 已达到 `4.74e-13`；
2. PE--连续 Weyl 轴向差异只有 `2.50e-4 dB / 6.31e-5 rad`；
3. Bellhop 5001/10001 波束差异远小于横向误差；
4. 横向误差主要是相位，并随频率、偏移平滑增加；
5. 直达和反射剖面表现出相同趋势。

因此优先怀疑以下因素，但目前不把任何一项写成已证实原因：

- Bellhop 二维射线束场与 PE 三维横向角谱场的源维度差异；
- 从 PE Gaussian 孔径到 Bellhop `.sbp` 的幅相映射并非完整等价；
- Bellhop Gaussian/ray beam 的远轴相位曲率近似；
- Bellhop 步长、beam type 或源图样插值对外侧相位的影响；
- 横向剖面验收半径是否应依据实际换能器有效能量区预先定义。

## 9. 建议的下一步执行顺序

### 优先级 A：用独立 AS 横向剖面定位归属

新增一个只读取现有正式配置的验证脚本，在 `4/6/8 kHz` 对 `97 m` 和
`103 m` 直接计算独立一步二维精确角谱横向剖面，并分别比较：

1. `PE ↔ 独立 AS`；
2. `独立 AS ↔ Bellhop`；
3. `PE ↔ Bellhop`。

这一测试只需 6 次一步 FFT，不需要重跑 65 点完整 PE。判定规则应预先固定：

- 如果 PE--AS 全横向复场接近浮点误差，而 AS--Bellhop 保留约
  `0.061 rad`，则差异归于 Bellhop/源映射层，继续禁止修改 PE marching；
- 如果 PE--AS 也在相同位置偏离，先检查剖面提取、网格坐标和验证实现，
  仍不直接修改生产 propagator。

### 优先级 B：做 Bellhop-only 横向相位审计

只针对最差的 `8 kHz、19.53125 m` 附近运行低成本 Bellhop 扫描：

- `step_m = 0.1, 0.05, 0.025`；
- 角扇区 `±20, ±30, ±45 deg`；
- `.sbp` 角度采样加密；
- 对可选 beam type 做对照；
- 保留 5001/10001 波束结果作为基线。

重点看相位曲率是否收敛，不再重复已经充分通过的轴向 65 点 PE。

### 优先级 C：审计 2D/3D Gaussian 源映射

分别建立：

- PE 的二维横向 Fourier/Weyl 表达；
- Bellhop 二维射线平面所需的线源或点源角度权重；
- 当前 `cos(theta)` Jacobian 的来源和适用维度。

比较解析远场和有限距离 Fresnel 相位，明确 `.sbp` 应表示幅度指向性、功率
指向性还是包含额外维度/Jacobian 的等效权重。这一步最可能解释“幅度几乎
一致、远轴相位缓慢偏离”的现象。

### 优先级 D：恢复新的 Bellhop 归一化审计

更换或单独部署稳定版本的 Bellhop 后，重新运行 `4/6/8 kHz` 点源自由场
归一化审计，验证 `1/R`、空间相位方向和固定 `4π` 换算。不得使用频率或
距离相关拟合。

### 优先级 E：实际换能器阶段

当 Gaussian 基准闭环后，再用实测或建模换能器的频率相关幅相指向性替换
`.sbp`/PE 初场，并重新验证：

- `4--8 kHz` 全带宽；
- 横向窗口与边缘能量；
- 是否仍适合 `192 m + no sponge`；
- 轴向 `Q(f)`、总信道、群时延和 PDP。

当前 `192 m` 结论只属于 `sigma=0.3 m` Gaussian 源，不能直接继承给新的
实际换能器。

## 10. 当前工程建议

1. 保持 PE marching、FFT convention、Gaussian 源和默认生产参数不变。
2. 对当前轴上、平面海面、Gaussian 源的 4--8 kHz 相对信道，可把本轮结果
   视为强交叉验证证据。
3. 正式验证状态继续保留 `passed=false`，直到横向 AS 三方对比解释
   `0.061 rad / 6.13%` 的偏差；不要事后放宽原门槛。
4. 下一次计算优先执行“独立 AS 横向剖面”，它成本最低、信息增益最大。
5. 在完成实际换能器源建模前，不据此自动修改生产窗口或重新引入 sponge。

## 11. 输出文件

- 主结果：`results/validation/pe_bellhop_unfolded_flat_gaussian/pe_bellhop_unfolded_flat_gaussian_validation.mat`
- 自动检查：`results/validation/pe_bellhop_unfolded_flat_gaussian/checks.csv`
- 宽带比较：`results/validation/pe_bellhop_unfolded_flat_gaussian/relative_wideband.csv`
- 横向剖面：`results/validation/pe_bellhop_unfolded_flat_gaussian/spatial_profiles.csv`
- 波束审计：`results/validation/pe_bellhop_unfolded_flat_gaussian/beam_convergence.csv`
- PE--AS 恒等式：`results/validation/pe_bellhop_unfolded_flat_gaussian/unfolded_identity.csv`
- Bellhop 正式场：65 组 `.env/.sbp/.shd/.prt`
- Bellhop 剖面：3 组 `.env/.sbp/.shd/.prt`
- Bellhop 波束审计：6 组 `.env/.sbp/.shd/.prt`
- 正式执行日志：`formal_unfolded_run.log`、`formal_unfolded_run.err.log`

![正式验证汇总图](../results/validation/pe_bellhop_unfolded_flat_gaussian/pe_bellhop_unfolded_flat_gaussian.png)
