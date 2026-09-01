# PE–Bellhop 平面海面完整对比思路与结果汇总

更新日期：2026-08-26  
文档性质：当前验证链的总览与结论索引  
适用范围：均匀声速、平面压力释放海面、生产 Gaussian 源、直达与一次海面反射

## 1. 文档目的

仓库中已有多份 PE–Bellhop 报告，但它们分别记录不同阶段：早期原坐标小偏移
验证、自由场点源误差预算、生产 Gaussian 窗口审计、展开坐标正式对比、横向
三方审计以及 Bellhop 坐标极限检查。早期的
`pe_bellhop_flat_surface_cross_validation_complete_report.md` 已明确标为历史报告，
没有覆盖 2026-08 完成的后续验证。

本文将这些证据串成一条当前完整验证链，回答最初的核心问题：

> 在相同介质、相同 Gaussian 发射声场、相同传播距离和平面压力释放反射条件下，
> 当前 PE 与 Bellhop 是否一致；如果不完全一致，差异是否说明 PE 公式有问题？

这里的“完整”指当前**平面海面验证范围**内的完整证据链，不表示已经验证粗糙
海面、随机散射、实际换能器、深度相关声速或海底相互作用。

## 2. 当前结论摘要

当前证据支持以下分级结论：

1. **PE marching 实现通过。** 生产多步 PE 与独立一步精确角谱在完整场和横向
   剖面上达到约 `10^-13` 的复数相对误差，没有证据支持修改平方根传播算子。
2. **Bellhop 的自由场归一化、路径时延和压力释放反射符号通过。** Bellhop 输出
   符合 `1/R` 归一化；直达和一次海面反射路径能被正确识别，未混入海底反射。
3. **轴上 4–8 kHz 相对信道高度一致。** 对生产 Gaussian 源，PE 与 Bellhop 的
   反射/直达比 `Q(f)`、总相对信道 `1+Q(f)`、群时延和相对 PDP 均通过门槛。
4. **展开坐标本身有效。** Bellhop-only 检查表明，当前 97/103 m 镜像展开结果与
   正常坐标的 89°、89.5°近垂直结果在几何归一化后几乎相同。
5. **横向远轴复场尚未完全闭合。** 在 8 kHz、横向偏移 19.53125 m 处，幅度仅差
   `0.0158 dB`，但相位差为 `0.0612 rad`，超过 `0.05 rad` 门槛；复场误差
   `0.0613` 也超过 `0.02` 门槛。
6. **横向差异不来自 PE。** 独立角谱与 PE 达到浮点精度一致，而 Bellhop 与独立
   角谱保留相同的约 `0.061 rad` 差异。27 个 Bellhop 步长、角扇区和 `.sbp`
   采样组合也未消除它。剩余问题更可能位于 Gaussian 源的 2D/3D 映射或
   Bellhop 射线束横向相位表示。

因此，当前最准确的总体判断是：

> **当前 PE 传播公式和轴上平面海面相对信道得到强验证；PE 与 Bellhop 的完整
> 横向复场尚未达到预设 hard threshold，但该差异不构成修改 PE marching 的依据。**

## 3. 为什么验证链需要分层

直接拿一个 PE 结果和一个 Bellhop 结果比较绝对压力，会同时混入以下因素：

- PE marching 是否正确；
- 有限横向窗口和 FFT 周期延拓；
- sponge 是否扰动中心接收场；
- PE 初始面和 Bellhop 点源的幅相参考；
- Bellhop 的 `1/R` 归一化；
- Gaussian 源到 Bellhop `.sbp` 的角谱映射；
- 接近 90°时 Bellhop 的几何表示；
- 压力释放海面反射的 `-1` 相位；
- 2D Bellhop 射线束与 3D 横向 PE 场的定义差异。

为了避免把上述误差错误归因于 PE 公式，当前验证采用以下层次：

| 层级 | 比较 | 主要目的 |
|---|---|---|
| Level 0 | Bellhop–解析自由场 | 确认 Bellhop 归一化、传播相位和时延 |
| Level 1 | 生产 PE–独立一步角谱 | 检查 marching、FFT 顺序、步长累计和载波处理 |
| Level 2 | PE/角谱–解析 Gaussian/Weyl | 确认有限宽 Gaussian 的自由场物理结果 |
| Level 3 | PE–Bellhop 轴上相对信道 | 独立模型交叉验证直达/反射幅相和时延 |
| Level 4 | PE–角谱–Bellhop 横向剖面 | 定位远轴复场差异属于哪一侧 |
| Level 5 | Bellhop 展开–89°–89.5° | 验证旋转/镜像展开不是几何伪差 |

其中 PE–角谱是**实现级 hard check**，因为两者使用相同的完整 Helmholtz 横向
色散关系；真正的外部模型证据来自连续 Weyl/解析结果和 Bellhop。

## 4. 当前正式物理环境与数值设置

最新正式展开坐标验证使用以下统一配置：

| 类别 | 设置 |
|---|---|
| 介质 | 均匀水体，`c=1500 m/s` |
| 原物理坐标 | Tx `(0,0,100 m)`，Rx `(0,0,3 m)` |
| 发射场 | 生产 Gaussian，`sigma_src_m=0.3 m` |
| 频率 | `4–8 kHz`，65 点，间隔 `62.5 Hz` |
| 海面 | 平面压力释放，反射系数 `-1` |
| PE 网格 | `984×984`，实际横向窗口 `192.1875 m` |
| PE 横向采样 | `dx=dy=0.1953125 m` |
| 数值边界 | `sponge off` |
| Bellhop 正式波束数 | `10001` |
| Bellhop 收敛复核 | `5001/10001`，4/6/8 kHz |
| Bellhop 步长 | `0.05 m` |
| Bellhop 角扇区 | `[-30°,30°]`（展开坐标） |
| `.sbp` | 2401 点，低于 `-120 dB` 截断 |
| 关闭功能 | 粗糙面、随机海面、气泡、Doppler、海底路径和通信处理 |

生产 Gaussian 初场为：

\[
\Psi_0(x,y)=\exp\!\left[-\frac{x^2+y^2}{2\sigma^2}\right],
\qquad \sigma=0.3\ \mathrm{m}.
\]

映射到 Bellhop 的轴对称角度指向性为：

\[
D(\theta,f)=\cos\theta\,
\exp\!\left[-\frac{(k\sigma\sin\theta)^2}{2}\right].
\]

该 `.sbp` 在轴向归一化为 0 dB。验证不允许按距离或频率逐点拟合源强。

## 5. 展开坐标与比较量

### 5.1 镜像展开

原始竖直路径被映射到 Bellhop 水平距离轴：

\[
L_d=z_{tx}-z_{rx}=97\ \mathrm{m},
\qquad
L_r=z_{tx}+z_{rx}=103\ \mathrm{m}.
\]

Bellhop 只计算匹配半空间中的自由传播：

\[
H_d=P_{BH}(97\ \mathrm{m}),
\qquad
H_r=-P_{BH}(103\ \mathrm{m}).
\]

反射路径的负号表示压力释放海面。由于没有在 Bellhop 距离终点设置真实边界，
不会引入人为的终点反射或 `r=0` 接收退化问题。

### 5.2 主要验收量

绝对源强只作为诊断。正式宽带硬验收使用：

\[
Q(f)=\frac{H_{reflect}(f)}{H_{direct}(f)},
\qquad
H_{rel}(f)=1+Q(f).
\]

这些量消除共同的固定源强归一化，但保留反射相对直达的真实幅度比、相位差和
额外传播时延。横向剖面使用：

\[
G(\rho,f)=\frac{H(\rho,f)}{H(0,f)}.
\]

项目采用 `exp(-i\omega t)` 合成约定；保存响应按当前公共 physical/direct-DSP
相位接口读取，不再由验证脚本重复恢复未知载波。

## 6. 各阶段工作与结果

### 6.1 早期原坐标准垂直验证：路径和时延先通过

早期案例使用 Tx 深度 80 m、Rx 深度 10 m、3/6/9 m 水平偏移以及平面海面。
它完成了第一轮路径拓扑和到达时间检查：

- PE 和 Bellhop 都识别出直达与一次海面反射；
- PE–Bellhop 六条路径最大时延差 `0.027285 ms`；
- Bellhop–解析像源最大时延差 `4.0135e-06 ms`；
- `H_f=H_direct_f+H_reflect_f` 最大闭合误差 `3.5762e-18`；
- Bellhop `.shd` 与 arrival 合成最大 TL 差 `0.04491 dB`。

该阶段的跨几何绝对幅度严格门槛未全部通过：公共标定跨度 `1.6540 dB`，反射
TL 残差 RMS/最大值为 `1.6295/2.0751 dB`。随后确认这些比较混入小窗口、
sponge、源归一化以及准垂直非共线几何，因此该阶段现作为**历史路径证据**，
不再承担当前同源 Gaussian 主幅度验收。

### 6.2 自由场点源审计：排除 marching 公式错误

反射关闭后，Bellhop normalization audit 确认其输出为 `1/R` 型；转换到
`exp(ikR)/(4\pi R)` 只使用一次固定的 `1/(4\pi)` 系数，不逐点拟合。

生产 PE 多步 marching 与独立一步精确角谱最大相对 L2 误差为
`1.4287e-13`，通过 hard check。但有限空间平面截断的点源初场与解析球面波、
Bellhop 曾出现最大 `18.2012 dB` 差异。后续窗口、Weyl 初始化和 sponge 误差
预算表明，该差异来自有限孔径/FFT 周期像与 sponge 对点源场的作用，而不是
平方根传播算子。因此项目停止用该有限窗口点源案例修改 PE 核心，转向实际生产
Gaussian 源。

### 6.3 生产 Gaussian 的窗口前提：建立可信 PE 参考

生产 Gaussian 源下重新执行窗口与 sponge 审计，并进一步检查完整
`PE上行 → 海面反射 → PE下行` 链路。结果为：

- `50 m + default sponge` 会显著扰动中心接收场，不适合作为严格验证参考；
- 固定海面下，完整反射链路在约 `192.1875 m + no sponge` 时严格收敛；
- 不同海况、15 个固定随机海面 realization 的 4 kHz 检查中，192 m 为
  `15/15` 严格通过；
- 160 m 的接收响应和中心场为 `15/15` 通过，但全过程边缘门槛为 `0/15`。

这一步不直接证明 Bellhop 与 PE 一致，但为正式交叉验证建立了不会被强 sponge
或小窗口主导的 PE 数值前提。

### 6.4 正式展开坐标 Gaussian 宽带对比

正式 `984×984`、4–8 kHz/65 点结果如下：

| 指标 | 实测值 | 门槛 | 状态 |
|---|---:|---:|:---:|
| PE 多步–独立一步角谱 | `4.7395e-13` | `1e-10` | PASS |
| Bellhop 5001/10001 最大 TL 差 | `1.3815e-6 dB` | `0.1 dB` | PASS |
| Bellhop 5001/10001 最大相位差 | `9.6636e-6 rad` | `0.02 rad` | PASS |
| `Q(f)` 最大 TL 差 | `4.9373e-4 dB` | `0.25 dB` | PASS |
| `Q(f)` 最大相位差 | `1.8271e-3 rad` | `0.05 rad` | PASS |
| 内部频点最大群时延差 | `62.4 ns` | `20 us` | PASS |
| 相对 PDP 峰时延差 | `0` | `0.125 ms` | PASS |
| `1+Q` 最大 TL 差 | `0.00808 dB` | `0.25 dB` | PASS |
| `1+Q` 最大相位差 | `0.0296 rad` | `0.05 rad` | PASS |
| 横向剖面最大 TL 差 | `0.015822 dB` | `0.25 dB` | PASS |
| 横向剖面最大相位差 | `0.061204 rad` | `0.05 rad` | **FAIL** |
| 横向剖面最大复场误差 | `0.061277` | `0.02` | **FAIL** |
| `H` 分量闭合 | `3.8790e-18` | `1e-12` | PASS |

轴向反射/直达幅度比约为 `0.94175`。直达和反射路径长度相差 6 m，对应额外
时延：

\[
\Delta\tau=\frac{103-97}{1500}=4\ \mathrm{ms}.
\]

PE 与 Bellhop 的相对 PDP 峰位置完全一致。正式脚本仍保留 `passed=false`，因为
预设规则要求所有 hard checks 同时通过；不能用轴向通过覆盖横向失败。

![展开坐标正式对比](../results/validation/pe_bellhop_unfolded_flat_gaussian/pe_bellhop_unfolded_flat_gaussian.png)

### 6.5 PE–角谱–Bellhop 横向三方审计

为判断 `0.0612 rad` 横向相位差属于哪一侧，验证脚本直接读取正式 MAT，并在
4/6/8 kHz、97/103 m 距离上独立实现一步角谱传播：

| 比较 | 最大 TL 差 | 最大相位差 | 最大复场误差 |
|---|---:|---:|---:|
| PE–AS direct | `3.55e-12 dB` | `3.05e-13 rad` | `5.10e-13` |
| PE–AS reflect | `1.88e-12 dB` | `1.59e-13 rad` | `2.43e-13` |
| Bellhop–AS direct | `0.015822 dB` | `0.061204 rad` | `0.061166` |
| Bellhop–AS reflect | `0.012736 dB` | `0.051673 rad` | `0.051650` |

PE 与独立角谱在完整横向剖面上达到数值精度一致；Bellhop 与角谱保留了和
PE–Bellhop 相同的偏差。因此横向 hard check 的失败已被定位到 Bellhop/源映射
一侧，而不是 PE marching 或 PE 场提取。

![PE–AS–Bellhop 横向审计](../results/validation/pe_as_bellhop_transverse/pe_as_bellhop_transverse.png)

### 6.6 Bellhop-only 数值参数扫描

8 kHz 下固定 10001 波束，扫描：

- Bellhop 步长：`0.1/0.05/0.025 m`；
- 角扇区半宽：`20°/30°/45°`；
- `.sbp` 采样数：`1201/2401/4801`。

27 个案例全部成功运行，但通过数为 `0/27`。最佳案例最大复场误差仍为
`0.061166`，相位仍为 `0.061204 rad`。因此当前差异不能归因于这些已扫描的
Bellhop 步长、角扇区或 `.sbp` 采样密度。

### 6.7 展开坐标与正常近垂直坐标的极限检查

Bellhop-only 脚本比较三种情况：

1. 当前展开坐标：直达 97 m、镜像反射 103 m；
2. 正常坐标 89°：水平偏移 1.693141 m；
3. 正常坐标 89.5°：水平偏移 0.846506 m。

正常坐标使用真实压力释放海面和远离目标路径的匹配海底。去除已知路径长度、
`1/L` 扩展、Gaussian 指向性及传播相位后：

- 最大归一化幅度差：`1.2253e-6 dB`；
- 最大归一化相位差：`9.2934e-5 rad`；
- 最大解析到达时间误差：`3.131 ns`；
- 三个案例都得到直达和一次海面反射，海底反射数均为 0；
- 修正反射系数与 `-1` 的最大复数残差为 `1.23e-4`。

全部检查通过。这证明当前展开方案是均匀介质、平面海面条件下正常 Bellhop
近垂直传播的可靠数值代理。89.5°并未表现出相对 89°严格单调减小的误差，后续
无需继续逼近 89.9°或 89.99°。

![Bellhop 坐标极限检查](../results/validation/bellhop_rotation_limit/bellhop_rotation_limit.png)

## 7. 综合判定矩阵

| 要回答的问题 | 当前证据 | 判定 |
|---|---|:---:|
| PE 平方根 marching 是否实现错误 | PE–独立 AS 为 `10^-13` 量级 | 否 |
| PE 传播方向、FFT 顺序和载波处理是否自洽 | 一步/多步、群时延和公共接口闭合通过 | 是 |
| Bellhop 自由场归一化是否明确 | `|p|R` 稳定，统一 `1/(4π)` 转换 | 是 |
| 直达和一次海面反射路径是否一致 | 早期原坐标和当前极限检查均通过 | 是 |
| 平面压力释放反射的相对相位是否一致 | `Qcorr≈-1`，轴向 `Q(f)` 通过 | 是 |
| 4–8 kHz 轴上相对信道是否一致 | TL、phase、group delay、PDP 全通过 | 是 |
| 展开坐标能否代表近垂直 Bellhop | 与 89°/89.5°归一化结果一致 | 是 |
| 整个横向复场是否达到原 hard threshold | 8 kHz 外侧相位和复场失败 | 否 |
| 横向失败是否来自 PE | PE–AS 通过，Bellhop–AS 保留误差 | 否 |
| 是否应修改 PE marching | 当前没有支持证据 | **不建议** |

## 8. 为什么“正式状态 false”与“核心结论通过”并不矛盾

正式展开验证采用合取规则：12 项检查只要有一项失败，整体 `passed=false`。
目前 10 项通过、2 项失败；失败项都来自远轴横向复场，而不是轴上接收信道：

- 当前实际接收机位于传播轴上；
- 轴上 4–8 kHz `Q(f)`、`1+Q(f)`、群时延和 PDP 已通过；
- PE 与独立角谱的横向结果也已通过；
- Bellhop-only 参数扫描未改变残差。

因此应同时保留两句话：

1. **正式完整横向验收仍为 `passed=false`；**
2. **当前 PE marching 和轴上平面海面相对信道已获得强一致性证据。**

删除其中任何一句都会使结论失真。

## 9. 当前建议

### 9.1 后续 PE–Bellhop 验证方式

- 使用镜像展开作为主定量方案；
- 保留正常坐标 89°作为轻量几何回归；
- 不继续追求 89.9°、89.99°，避免把工作变成 Bellhop 端点角数值测试；
- 继续使用相对量 `Q(f)`、`1+Q(f)` 和轴归一化横向场，禁止逐频/逐点拟合；
- 不修改 PE marching、FFT convention 或 Gaussian 数学定义。

### 9.2 若必须闭合横向 hard check

下一项有价值的工作是研究 Bellhop `.sbp` 对 3D 径向 Gaussian 的 2D 等效映射：

- 检查 2D/3D spreading 与角谱 Jacobian；
- 检查 Bellhop beam pressure 在横向截面中的相位曲率定义；
- 必要时增加独立 2D Helmholtz/角谱参考，而不是继续扫描已经无效的波束数和步长。

这属于源与模型维度解释，不是 PE 核心修复。

### 9.3 结论边界

当前结论不能直接外推到：

- 粗糙海面、Kirchhoff/SSA 漫散射；
- 随机海面 Monte Carlo 或统计信道；
- 深度相关 SSP、折射和焦散；
- 海底反射、吸收或多次边界作用；
- 气泡、Doppler 和通信接收处理；
- 实际换能器。

更换实际换能器后，必须用其测量或建模的复幅相指向性替换 Gaussian `.sbp`，
同时重新确认 PE 横向窗口和边缘能量；不能直接继承 `192.1875 m` 结论。

## 10. 主要脚本、报告与数据索引

| 内容 | 入口/报告 | 主要结果目录 |
|---|---|---|
| 历史原坐标矩阵 | `validate_pe_bellhop_flat_surface_matrix_vertical.m` | `results/validation/pe_bellhop_flat_surface_matrix/` |
| 自由场四层审计 | `validate_pe_bellhop_freefield_vertical.m` | `results/validation/pe_bellhop_freefield/formal/` |
| Gaussian 窗口/sponge | `validate_gaussian_window_convergence_vertical.m` | `results/validation/pe_gaussian_window_sponge/` |
| 完整反射链路窗口 | `validate_reflected_chain_window_convergence_vertical.m` | `results/validation/pe_reflected_chain_window/` |
| 正式展开坐标对比 | `validate_pe_bellhop_unfolded_flat_gaussian_vertical.m` | `results/validation/pe_bellhop_unfolded_flat_gaussian/` |
| 横向 PE–AS–Bellhop | `validate_pe_as_bellhop_transverse_vertical.m` | `results/validation/pe_as_bellhop_transverse/` |
| Bellhop 参数扫描 | `validate_bellhop_transverse_parameter_scan_vertical.m` | `results/validation/pe_as_bellhop_transverse/bellhop_scan/` |
| Bellhop 坐标极限 | `validate_bellhop_rotation_limit_vertical.m` | `results/validation/bellhop_rotation_limit/` |

分阶段详细报告：

- `reports/pe_bellhop_unfolded_flat_gaussian_report.md`；
- `reports/pe_as_bellhop_transverse_audit_report.md`；
- `reports/bellhop_transverse_parameter_scan_report.md`；
- `reports/bellhop_rotation_limit_report.md`；
- `reports/pe_bellhop_freefield_validation_report.md`；
- `reports/pe_reflected_chain_window_validation_report.md`；
- `reports/pe_random_surface_window_robustness_report.md`。

早期 `reports/pe_bellhop_flat_surface_cross_validation_complete_report.md` 保留为历史
证据，不应覆盖本文所汇总的后续结果。

