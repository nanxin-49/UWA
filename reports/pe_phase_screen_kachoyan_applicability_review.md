# PE 粗糙海面相位屏的 Kachoyan--Macaskill 适用性审查

日期：2026-09-16  
范围：只审查当前 4 kHz、fixed-seed PM、1-transverse PE 对比；未重跑 MATLAB/Bellhop，未修改任何模型代码。

## 为什么需要这份报告

当前项目已经分别验证了 PE marching、入射自由场、Bellhop internal-wall 几何、
`Reflect2D`、pressure-release 相位、Gaussian beam `p/q` 状态、receiver selection
和主要数值收敛项。独立 Helmholtz BIE 对光滑正弦面的测试还显示：Bellhop 的
local-specular 结果在这些受控环境中明显比当前 PE Model-0 更接近 Helmholtz
参考。由此，剩余 PE--Bellhop 差异不能再笼统归因于步长、beam count、坐标旋转、
复数约定或接收列选错。

当前 PE Model-0 在表面使用局部乘法算子

```text
psi_ref(x) = -exp(+i 2 k eta(x)) psi_inc(x),
```

它正是近法向 Kirchhoff phase-changing-screen 思想的一种实现。fixed-PM 环境中
已经测得 receiver-line `E_G=0.9583`、phase RMS `1.1077 rad`，远高于已确认的
数值误差预算；受控正弦面中，误差又随高度和表面波数有规律地增长。因此需要
回到该相位屏的理论前提，回答三个问题：

1. 当前表面尺度是否实际违反 Kirchhoff 高频/缓变条件；
2. 当前 Gaussian 入射场是否足够接近论文中的法向入射；
3. 若表面尺度条件形式上成立，为什么局部 `2k eta` 算子仍会偏离 Helmholtz/BIE
   和 Bellhop local-specular 结果。

本报告的目的不是把 Bellhop 当作无条件精确解，也不是为某个 solver 做经验校准；
而是利用论文的渐近条件、全部已完成正弦/PM 测试和独立 BIE，判断当前残差是否
符合 phase-screen 近似超出适用范围后的物理表现，并为下一步是否需要非局部表面
算子提供依据。

## Kachoyan--Macaskill 论文的核心内容

Kachoyan 与 Macaskill 研究压力释放粗糙表面的声散射，并将能够处理任意粗糙度和
多次散射的积分方程结果，与 Kirchhoff/tangent-plane 类近似进行比较。与本项目
直接相关的结论可概括为：

1. **phase-changing screen 是有条件的等效描述。** 在 Kirchhoff 近似有效、
   声学频率较高且接近法向入射时，粗糙面反射可近似为平面上的相位改变；近法向
   情况下主要相位扰动为

   ```text
   Delta phi approximately 2 k eta.
   ```

   这给出了当前 PE Model-0 的理论来源，但不表示该公式对任意入射角谱和任意
   粗糙谱都精确。
2. **关键小参数是表面相对尺度，而不只是绝对高度。** 典型适用要求是
   `kL >> 1` 和 `h/L << 1`：声波长应远小于表面相关尺度，表面高度相对横向
   尺度应较小。论文并不简单要求 `kh << 1`，所以 `2kh` 大于 1 不能单独作为
   否定 Kirchhoff 近似的门限。
3. **近法向条件属于模型本身的一部分。** 对有限入射角，相位应依赖入射和反射
   波矢的法向分量；把所有角谱分量统一写成 `2k eta` 会遗漏 `k_z` 差异。对于
   有宽角谱的波束，“中心射线近法向”并不足以保证整个场都处在近法向极限。
4. **宽尺度随机面比单尺度光滑面更困难。** 论文发现，包含宽广尺度范围的
   modified Pierson--Moskowitz 类随机面，其 Kirchhoff/phase-screen 精度可能比
   光滑 Gaussian 谱表面差，并将其与多尺度结构和较低可微性联系起来。长波部分
   满足条件，不能自动保证短波部分及不同尺度之间的耦合也满足同一局部近似。
5. **这些条件是渐近适用性条件，不是充分误差界。** 即使 `kL` 很大、平均坡度
   很小，局部 screen 仍可能遗漏边界上不同位置之间的波动耦合、边界密度/幅度、
   曲率衍射和多次横向相互作用。因此本报告同时检查尺度条件、入射角谱和 BIE
   结果，而不把单个 `kL` 数字当成最终裁决。

论文原文为 B. J. Kachoyan and C. Macaskill, “Acoustic scattering from an
arbitrarily rough surface,” *JASA* 82(5), 1720--1726 (1987)，见
[DOI 10.1121/1.395165](https://doi.org/10.1121/1.395165)。

## 结论摘要

当前 fixed-PM 工况**没有简单、全面地违反** Kachoyan--Macaskill 的
`kL >> 1` 和 `h/L << 1` 条件：即使用较严格的 `L=1/q`，最短保留尺度仍有
`kL=35.56`，而以 RMS 高度计的最不利 `h/L=0.0878`。因此，不能把现有大偏差
归结为“声波分辨不了海面尺度”或笼统的 Kirchhoff 高频条件完全失效。

更有力的解释是三项条件联合失效或变弱：

1. PE 入射场不是单一法向平面波：`theta_rms=8.144 deg`，95%/99% 角谱范围达到
   `15.97/21.28 deg`；实际 Bellhop wall hits 的最小 `|u.n|=0.831901`，对应局部
   法线入射角可达 `33.71 deg`。把所有角谱分量都乘以 `exp(i 2 k eta)` 不再准确。
2. `2 k sigma_eta=6.2467 rad`，表面造成多弧度、跨多个相位周期的调制。论文并不
   要求 `kh << 1`，所以这不是形式上的禁区；但有限角、局部几何和非局部耦合的
   小误差会被放大成明显的相干场误差。
3. PM 面跨越 12 个 Fourier 尺度，且高波数端虽只占部分高度方差，却主导坡度和
   曲率。`q>=0.2 rad/m` 的分量贡献约 51.7% 高度方差、86.7% 坡度方差和 97.8%
   二阶导数方差。单一相关尺度检验会掩盖这种宽谱导数敏感性。

受控 BIE 结果进一步表明：把 `2k eta` 改成逐入射分量的 `2k_z eta` 能显著改善
低波数正弦面，却对 `K=0.47 rad/m` 的高波数面只改善约 `1.22x`；继续加入局部
坡度修正也没有关闭误差。现有证据因此支持：**法向相位近似是偏差的一部分，但
主要剩余量来自局部乘法 phase screen 未包含的非局部表面耦合/衍射与幅度效应。**

## 重新分析：为什么尺度条件看似满足，结果仍然很差

这里首先需要修正一个容易产生误解的说法：当前实验并不是“完整满足论文的全部
条件后仍然失败”，而是**较好满足了 `kL >> 1` 和 `h/L << 1` 两个表面尺度条件，
但没有完整满足 phase-screen 等效所需的入射场、观测量和算子条件**。论文的两个
无量纲参数说明 Kirchhoff 类近似可能进入有效区，却不是针对任意波束、任意接收
距离和逐点复场误差的充分上界。

### 1. `kL` 和 `h/L` 是渐近参数，不是二元 PASS/FAIL 门限

`kL=35.6` 可以合理称为大于 1，却不是无穷大；相应的 `1/(kL)` 量级仍约为
`2.8%`。对强相干叠加问题，几个百分点的局部边界误差不一定只产生几个百分点的
接收场误差：它可移动干涉零点、焦散和相位绕转位置。PM case 又有
`2k sigma_eta=6.25 rad`，所以局部近似误差经过多周相位调制后可能形成
`O(1)` 的归一化复场差异。

同样，`h/L` 的 RMS 值小并不等价于每个局部片段、每个谱分量和每个相互作用项
都有统一的小参数。fixed PM 的高 q 部分贡献 86.7% slope variance 和 97.8%
curvature variance；用一个 `L_eff` 压缩整个宽谱会丢失这一信息。

### 2. 论文讨论的基准是平面波；当前是有限宽度 Gaussian 波束

论文结论中的 `kL`、`h/L` 明确针对 plane-wave surface scatter。当前入射场为
`sigma=0.3 m` 的 Gaussian，其中心方向虽为法向，但完整角谱并不窄：

```text
theta_rms = 8.144 deg
theta_95  = 15.97 deg
theta_99  = 21.28 deg
```

再叠加最高 `6.65 deg` 的局部表面倾角，fixed-PM Bellhop hits 中相对局部法线的
入射角可达约 `33.7 deg`。这些射线不属于 grazing，但也不能再视为对整个场都
严格近法向。因此当前实验满足“中心射线法向”，不满足“所有重要角谱分量都处于
法向渐近极限”。

G1/G2 已经直接验证这一点：将统一 `2k eta` 改为逐分量 `2k_z eta` 后，弱低 K
case 的 phase error 改善约 25 倍。这说明论文条件与项目结果并不矛盾；项目实际
输入比论文的单一法向平面波更宽。

### 3. 当前 Model-0 只是 phase-only 乘法，不是完整 Kirchhoff 表面积分

当前 validation Model-0 执行的是

```text
incident field on mean plane
-> multiply by -exp(+i 2k eta(x))
-> propagate again from the mean plane.
```

完整 Kirchhoff/tangent-plane 表面积分还涉及真实表面坐标上的 incident field、
Green kernel、法向导数、obliquity/amplitude factor 和弧长 Jacobian。把该积分继续
约化为 mean-plane 上的纯相位乘法，需要比 `kL >> 1`、`h/L << 1` 更具体的
近法向和传播近似。也就是说，论文支持的是某一参数区间内的等效关系，并不能
自动证明本项目这条最简 phase-only 实现已经保留了完整 Kirchhoff 算子的全部
leading-order 项。

Model-2 已加入局部入射角和局部坡度相位，但高 K case 几乎没有继续改善。这表明
缺失量并非另一个可以写成 `c eta` 或局部 slope phase 的标量修正，更可能位于
amplitude/normal-derivative 以及不同表面位置、不同横向波数之间的非局部耦合。

### 4. band-limited `eta` 不代表 phase screen 本身仍是同一带宽

即使 `eta(x)` 只包含到 `q_max=0.471 rad/m`，非线性函数
`exp(i 2k eta(x))` 也会产生卷积级联和高次谐波。对单一正弦面，这对应
Jacobi--Anger/Bessel harmonics；对 PM 面，则对应各 Fourier mode 的多重卷积。
当 `2k sigma_eta=6.25 rad` 时，这些高阶项不会都很小。

所以用原始海面 `q_max` 得到的 `kL_min=35.6` 只能描述几何输入带宽，不能直接
界定 phase-screen 输出场的谱宽、焦散或相干零点。高度 sweep 和 fixed-PM 的
差别也由此变得清楚：`kh` 不必小才能使用 Kirchhoff，但较大的 phase modulation
index 会使任何被遗漏的 amplitude 或 spectral-coupling 项更容易在复场中显现。

### 5. 当前接收量比论文条件本身更苛刻

物理 Rx 位于 `z=3 m`，也就是反射后仅传播约 3 m。当前比较的是该近表面接收线
上的 normalized coherent complex field，而不是只比较总散射能量、平均强度或
远场主瓣。近表面场对高横向波数、局部聚焦、干涉零点和边界附近的非传播/弱传播
分量更敏感；一个在能量意义下仍“相当好”的 Kirchhoff approximation，可以在
逐点 complex-field 指标上表现明显较差。

这与论文的观察相容：论文既讨论近场强度涨落，也指出随入射变斜，视场外的
specular/near-specular surface portions 会越来越重要。当前有限波束和有限接收
窗口进一步加强了这种空间非局部性。

### 6. PM 的大误差不是单纯的指标或 Bellhop 故障

fixed-PM 的 `E_G=0.958` 确实会受到相干零点和归一化的敏感性影响，但
`E_aligned=0.679`、`rho_shape=0.771`、TL RMS `1.274 dB` 和 phase RMS
`1.108 rad` 同时表明：它不只是一个可去除的全局相位常数。另一方面，Bellhop
wall intersection、pressure-release phase、`p/q` rotation、receiver selector 和
收敛检查都通过，因此也没有证据把该残差重新归于 internal-wall bookkeeping。

需要保留一个重要限定：独立 BIE 已经在四个正弦环境中证明 Bellhop 明显更接近
Helmholtz 解，但尚未对完整 fixed-PM realization 运行 BIE。因此“PM 大残差主要
来自 PE phase screen”是由正弦 BIE、Model-1/2 趋势和 PM 谱诊断共同支持的最强
当前判断，而不是已经由 PM 全波解最终证明的事实。

### 7. 与论文随机面结论并不矛盾

论文的结论需要细分：随机调制表面在某些 case 中可能比简单周期正弦面更适合
Kirchhoff approximation，因为误差更局域；但包含很宽尺度范围、具有更强 fractal
特征或较低 correlation differentiability 的 modified-PM 面，又比平滑 Gaussian
spectrum 更困难。不能把论文简化成“随机面一定更好”或“PM 一定失败”。

当前 fixed PM 是有限带宽且几何可微的，并非真正数学分形面；不过其高 q 分量
已经主导几何导数，同时强相位调制和宽入射角谱并存，因此恰好集中体现了论文所
警告的宽尺度敏感性。

## 进一步建议

下一步不建议继续调节 `2k eta` 系数，也不建议用 amplitude/phase calibration
强迫 PE 与 Bellhop 重合。建议按以下顺序缩小物理归因范围：

1. **先建立完整 Kirchhoff surface-integral 中间基准。** 使用已经保存的 incident
   field，在真实表面坐标上保留 Green kernel、normal derivative、obliquity factor
   和 surface Jacobian；不拟合参数。先对现有四个 BIE 正弦 case 比较
   `phase screen -> full Kirchhoff integral -> Helmholtz BIE`。这会直接判断失败发生
   在 Kirchhoff approximation 本身，还是发生在 Kirchhoff 到纯 phase-screen 的
   二次约化。
2. **做入射角谱收窄审计。** 对弱低 K、强高度低 K 和弱高 K 三个代表 case，
   validation-only 地增大 Gaussian `sigma`，把 `theta_rms` 依次降至约
   `4 deg、2 deg、1 deg`，并始终对 BIE 比较。如果误差随角宽明显下降，可定量
   分离“near-normal 不充分”与“surface nonlocality”。这不是修改生产源，而是
   验证论文平面波极限。
3. **做 receiver distance 审计。** 保持 source/surface 不变，把反射后接收距离从
   `3 m` 扩展到 `10/30 m`，同时比较反射系数谱和 receiver field。若近场误差快速
   衰减，而远场角谱闭合，则当前问题主要是 phase screen 的 near-field 重建；若
   不衰减，则是边界算子本身的谱权重错误。
4. **优先比较 outgoing transverse spectrum，而不只比较接收线。** 按论文的
   reflection-coefficient 思路，保存每个入射 `k_x` 到每个出射 `k_x'` 的耦合，
   区分 specular phase、off-specular redistribution 和 amplitude error。当前单条
   receiver line 会把这些机制混合在一起。
5. **最后才进入 PM 全波基准。** 先选同一 fixed realization 的短支持或低
   `K_max` 版本，用现有 BIE/WGF 思路收敛；随后固定低 q coefficients，只逐级加入
   高 q modes。这样可以判断误差从哪个 `K_max` 开始增长，而不是直接面对完整 PM
   的计算成本。必要时再采用 windowed Green function、FMM 或 H-matrix 加速。

建议的最小下一步不是完整 PM BIE，而是第 1 项：在四个已有正弦 BIE case 上加入
**无拟合的完整 Kirchhoff surface integral**。它复用所有现有 source、surface、
receiver 和 BIE artifacts，变量最少，同时最直接回答“论文近似失败”还是“当前
phase-only 实现过度简化”。

## 当前工况和谱范围

- 声学参数：`f=4 kHz`，`c=1500 m/s`，`k=16.7551608 rad/m`，
  `lambda=0.375 m`。
- PE 源：1T Gaussian，`sigma=0.3 m`；fixed-PM 对比的 PE window/numerics 为
  `W=192.1875 m`、`nx=984`、`dx=0.1953125 m`、march step `0.05 m`。
- PM 约定：`U=6 m/s`，seed `260001`，span `160 m`，无 Hs 重归一化、平滑或
  拟合。工程中使用
  `E_1D(q)=alpha/(2q^3) exp[-beta g^2/(U^4 q^2)]`，
  `alpha=8.10e-3`，`beta=0.74`。
- 固定 realization：12 个 Fourier modes，
  `q=0.0392699...0.4712389 rad/m`；对应表面波长 `2pi/q=160...13.333 m`。
  理论 PM 谱峰 `q_peak=0.19140 rad/m`；本 realization 由
  `q_eff=sigma_slope/sigma_eta=0.28956 rad/m` 给出的导数加权尺度为
  `L_eff=1/q_eff=3.454 m`。
- 表面统计：RMS/max height `0.186410/0.499456 m`，RMS/max slope
  `0.053977/0.116674`，RMS/max curvature `0.020907/0.044594 1/m`，最小曲率
  半径 `22.424 m`（约 `59.8 lambda`，`kR_min=375.7`）。

这里同时给出两种尺度定义：`L=1/q` 是更严格的局部变化尺度；若把 `L` 取为
完整表面波长 `2pi/q`，所有 `kL` 会再大 `2pi` 倍。因此以下判断不会通过选择
较宽松的尺度定义人为“通过”。

## 论文条件与当前数值

| 论文条件/隐含前提 | 当前数值 | 是否满足 | 可能后果 |
|---|---:|---|---|
| 高频/大尺度：`kL >> 1` | PM 谱峰：`k/q_peak=87.54`；导数加权：`k/q_eff=57.86`；最短尺度：`k/q_max=35.56`。若用 `L=2pi/q_max`，最小值为 `223.4` | 主体满足；高波数端是“数十”而非极强渐近极限 | 不是当前大偏差的首要单因；但最短尺度的高阶误差会最先显现 |
| 小高度/缓变面：`h/L << 1` | 以 `h=sigma_eta`：谱峰 `0.0357`、有效尺度 `0.0540`、最短尺度 `0.0878`；以全局 `max|eta|` 和最短尺度计为 `0.235` | RMS 意义下满足或中等满足；极值组合仅边缘性小 | 平均几何温和，但局部极值处不再是很强的渐近小参数 |
| 表面局部光滑 | RMS/max slope `0.0540/0.1167`（约 `3.09/6.65 deg`）；`R_min=22.42 m=59.8 lambda` | 几何光滑性本身满足 | 排除由尖角、网格折线或曲率半径接近波长造成的主要失败 |
| 近法向单一入射 | `theta_rms=8.144 deg`；角谱 q90/q95/q99=`13.43/15.97/21.28 deg`；实际 hit 最大局部入射角约 `33.71 deg` | **不均匀满足** | `2k eta` 对宽角谱的相位增量系统性偏大；每一 `k_x` 分量应感受不同 `k_z` |
| phase screen 可用局部乘法表示 | 当前 Model-0 为 `-exp(+i 2k eta(x)) psi_inc(x)` | 高波数/宽谱下证据不足 | 忽略表面不同位置及角谱分量之间的非局部耦合、幅度/Jacobian、曲率衍射等 |
| 相位扰动可被稳定表示 | `2k sigma_eta=6.2467 rad`，`2k max|eta|=16.74 rad` | 不是论文明确禁区，但数值/物理敏感性很高 | 很小的局部相位公式误差即可改变相干叠加、产生相位包络和零点移动 |
| 宽尺度随机面的局部尺度均处于渐近区 | `q_max/q_min=12`；高 q 强烈支配 slope/curvature | **仅整体指标满足，逐尺度不充分** | 单一 `L` 无法证明整个 PM realization 都可由同一局部 Kirchhoff screen 描述 |

## 上一个 Goal 的全部测试环境补充审查

上一个 controlled-comparison / surface-operator Goal 实际覆盖了 flat、弱极限
正弦、高度 sweep、波数/曲率 sweep、四个 Helmholtz-BIE 判别点和 fixed PM。
所有粗糙面 case 都使用相同的 4 kHz、`c=1500 m/s`、`sigma=0.3 m` 入射场，
因此前述 `theta_rms=8.144 deg` 和 q95/q99=`15.97/21.28 deg` 对下表每一行
都成立。表中正弦面的严格尺度取 `L=1/K`、典型高度取 `h_rms=A/sqrt(2)`；
`h/L|max=AK` 同时也是最大坡度。使用完整波长 `2pi/K` 只会令 `kL` 再增大
`2pi`，不会改变判断。

| 环境 | A (m) | K (rad/m) | `kL=k/K` | `h_rms/L` / `h_max/L` | `2k h_rms` / `2kA` (rad) | `R_min/lambda` | 已有 PE--BH `E_G` / phase RMS | 论文尺度条件判断 |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| flat Stage 0 | 0 | 0 | infinity | 0 / 0 | 0 / 0 | infinity | flat BH error `0.005899` / `0.005899 rad` | 粗糙面条件平凡满足；给出 Bellhop 数值 floor |
| weak-limit | 0.0025 | 0.10 | 167.55 | 0.000177 / 0.000250 | 0.0592 / 0.0838 | 106667 | `0.001020 / 0.001011` | **强满足** |
| weak-limit | 0.005 | 0.10 | 167.55 | 0.000354 / 0.000500 | 0.1185 / 0.1676 | 53333 | `0.002041 / 0.002024` | **强满足** |
| weak/height | 0.010 | 0.10 | 167.55 | 0.000707 / 0.001000 | 0.2370 / 0.3351 | 26667 | `0.004093 / 0.004059` | **强满足** |
| height | 0.020 | 0.10 | 167.55 | 0.001414 / 0.002000 | 0.4739 / 0.6702 | 13333 | `0.008234 / 0.008166` | **强满足** |
| height/BIE Region II | 0.050 | 0.10 | 167.55 | 0.003536 / 0.005000 | 1.1848 / 1.6755 | 5333 | `0.021036 / 0.020869` | 尺度强满足；相位已非小量 |
| height | 0.100 | 0.10 | 167.55 | 0.007071 / 0.010000 | 2.3695 / 3.3510 | 2667 | `0.044179 / 0.043877` | 尺度满足；强相位敏感 |
| strong-height/BIE | 0.200 | 0.10 | 167.55 | 0.014142 / 0.020000 | 4.7391 / 6.7021 | 1333 | `0.100485 / 0.100142` | 尺度满足；多周相位敏感 |
| wavenumber | 0.020 | 0.20 | 83.78 | 0.002828 / 0.004000 | 0.4739 / 0.6702 | 3333 | `0.008369 / 0.007579` | **满足** |
| wavenumber | 0.020 | 0.35 | 47.87 | 0.004950 / 0.007000 | 0.4739 / 0.6702 | 1088 | `0.014140 / 0.009428` | 满足，但离极强渐近区更近 |
| high-K/BIE | 0.020 | 0.47 | 35.65 | 0.006647 / 0.009400 | 0.4739 / 0.6702 | 604 | `0.022816 / 0.012906` | 基本满足；高 K 残差已明确 |
| fixed PM | RMS 0.1864 | 0.0393--0.4712 | 35.56--426.67 | 有效 0.0540；最短尺度 0.0878 | RMS 6.2467 | min 59.8 | `0.958269 / 1.107714` | 整体满足、极值/宽谱仅有限渐近 |

注：flat 行的 `0.005899` 是 internal-wall Bellhop 对 flat reference 的数值重建
误差，不是粗糙反射模型误差；弱极限 `A=0.0025/0.005 m` 使用 Stage-1Y 后冻结的
comparison convention，旧 raw-convention FAIL 结果不作为当前物理结论。

### 高度 sweep：尺度条件满足仍不保证相位屏精度

在 `K=0.10 rad/m` 的全部高度点，`kL` 恒为 `167.55`，最大坡度仅从
`0.00025` 增至 `0.02`，最小曲率半径即使在 `A=0.20 m` 时仍为
`500 m=1333 lambda`。从 Kachoyan--Macaskill 的尺度条件看，这组面都非常平缓。
然而 PE--BH phase RMS 从 `0.00101 rad` 单调增至 `0.10014 rad`。

这组结果不支持“`kL` 或 `h/L` 先失效”的解释，反而显示：当 `2kA` 从
`0.0838` 增至 `6.702 rad` 时，宽角谱中每个分量的相位误差及被遗漏的非局部
幅度/耦合效应被逐步放大。论文不要求 `kh << 1`，所以不能把 `2kA>1` 当作
形式上的失效门限；它在本项目中应被解释为**误差敏感度指标**。

独立 BIE 的三个低 K 判别点给出更直接的证据：

| BIE case | Model-0 PE--BIE `E_G` | `k_z`-aware Model-1 | angle+slope Model-2 | Bellhop--BIE |
|---|---:|---:|---:|---:|
| `A=0.01, K=0.10` | 0.00315757 | 0.000514551 | 0.000514113 | 0.000008092 |
| `A=0.05, K=0.10` | 0.0160900 | 0.00398989 | 0.00395441 | 0.000038896 |
| `A=0.20, K=0.10` | 0.0809081 | 0.0507808 | 0.0500519 | 0.000099760 |

在弱面上，逐分量 `2k_z eta` 几乎关闭主要相位误差；高度增大后，即使
`kL`、坡度和曲率条件仍非常好，Model-1/2 仍留下 `E_G≈0.05`。因此强高度点
的剩余量不是局部法线几何没有写对，而是局部 screen 结构本身不充分。

### 波数 sweep：固定相位幅度下的高 K 证据

`A=0.02 m` 的四个波数点具有完全相同的 `2kA=0.6702 rad`。若误差主要由高度
相位量级决定，结果应近似不随 K 变化；实际情况相反：

- `K=0.10 -> 0.47 rad/m` 时，`E_G` 从 `0.00823` 增至 `0.02282`；
- phase RMS 从 `0.00817` 增至 `0.01291 rad`；
- TL RMS 从 `0.00914` 增至 `0.16344 dB`，增大约 17.9 倍；
- 同时 `h/L|max` 仍只有 `0.0094`，`kL` 仍为 `35.65`，曲率半径仍有
  `604 lambda`。

高 K BIE case 中，Model-0/1/2 的 `E_G` 分别为
`0.022906/0.021567/0.021565`，Bellhop--BIE 仅 `1.368e-4`。因此该 case 在论文
两个主尺度条件仍成立时，已经明确排除了“只需修正入射角或局部坡度”的解释。
它是当前非局部边界耦合判断中最干净的受控证据。

### 测试环境综合判定

| 环境族 | `kL` / `h/L` | 近法向性 | 当前偏差主特征 | 判定 |
|---|---|---|---|---|
| flat | 平凡满足 | 同一宽 Gaussian 角谱 | PE--AS 精确；Bellhop 给出有限数值 floor | 不支持传播/约定错误 |
| weak low-K sine | 强满足 | 仅平均近法向，不是单一法向波 | `k_z`-aware 修正非常有效 | 法向 `2k eta` 是主要一阶误差 |
| stronger-height low-K sine | 仍强满足 | 同上 | 相位多周化后 Model-1/2 留下显著残差 | 强相位放大 + 非局部效应 |
| weak high-K sine | 满足 | 同上 | 固定相位高度下 TL/complex error 随 K 增长；局部坡度修正无效 | 高 K 非局部/幅度耦合主导 |
| fixed PM | 主体满足、宽谱极值较弱 | 部分 wall hits 达约 33.7 deg | 高 q 主导导数，且 `2k sigma_eta=6.25` | 三种不利因素叠加，偏差最大 |

## 短尺度/高波数分量诊断

下表直接由冻结的 12-mode coefficient file 分解得到；“curvature”列使用
`q^4 |eta_q|^2` 权重，表示二阶导数方差贡献。

| retained band | height variance | slope variance | curvature variance |
|---|---:|---:|---:|
| `q >= 0.20 rad/m` | 51.66% | 86.74% | 97.84% |
| `q >= 0.30 rad/m` | 37.60% | 75.08% | 92.38% |
| `q >= 0.35 rad/m` | 25.17% | 60.45% | 82.82% |
| `q >= 0.40 rad/m` | 18.61% | 48.93% | 71.52% |
| 单独 `q_max=0.47124 rad/m` | 17.77% | 47.07% | 69.22% |

最高波数 mode 的振幅为 `0.11113 m`，自身 `qA=0.0524`，所以它并不是陡峭或
几何奇异的波；问题在于其导数权重极大。换言之，高 q 分量**没有明显违反**
`kL >> 1`，也没有单独违反小坡度，却把 PM 面推进到“局部相位修正不够、需要
非局部波动耦合”的敏感区。这与论文对宽尺度、低可微随机面的警告一致。

## 入射角与现有误差证据

冻结 Gaussian 入射谱给出：

- `kx_rms=2.3570 rad/m`，`<kz>/k=0.9899505`；
- `|kx|` 的 q95/q99 为 `4.6097/6.0809 rad/m`；
- 相对于 Model-0，遗漏的平均相位系数
  `<2(k-kz)>=0.33676 rad/m`，其谱内 RMS 为 `0.58972 rad/m`。

用 PM 的 `sigma_eta` 估算，这一有限角修正约为 `0.0628 rad`（均值尺度）或
`0.1099 rad`（分量 RMS 尺度）。它足以解释一部分误差，但仍远小于 fixed-PM
实测 receiver-line phase RMS `1.1077 rad`，因此“非严格法向”不是全部原因。

受控正弦/BIE 审计与此一致：

- `k_z`-aware Model-1 对弱低 K、Region-II 低 K、强高度低 K 的误差改善分别约
  `25.0x / 5.11x / 1.61x`；
- 对弱高 K (`A=0.02 m`, `K=0.47 rad/m`) 仅改善 `1.22x`；
- 该高 K case 中 Bellhop 对 BIE 的 `E_G=1.37e-4`，而原 PE Model-0 为
  `2.29e-2`；加入局部角度+坡度耦合仍未通过高 K gate，最终审计结论为
  `NONLOCAL_EFFECT_REQUIRED`。

fixed-PM 本身的 receiver-line 差异为 `E_G=0.9583`、TL RMS `1.274 dB`、phase
RMS `1.108 rad`。与此同时，Bellhop wall residual 为 `7.27e-15 m`、grazing
fraction 为 0、pressure-release phase error 为 `7.11e-15 rad`，说明该差异不能
合理归因于 internal-wall 交点、反射相位或坐标映射失败。

## 最终判断

最可能的偏差来源按证据强度排序为：

1. **局部 multiplicative phase-screen 的非局部性缺失**：高 q 对坡度/曲率的
   支配、BIE high-K 结果和 G3 局部坡度修正失败共同支持这一点。
2. **法向相位 `2k eta` 对有限角谱的过度简化**：量级和误差方向已由 G1/G2
   验证，但只能关闭低 K 问题的一部分。
3. **强多周相位调制放大模型误差**：`2k sigma_eta=6.25 rad` 不是 Kirchhoff
   的形式禁区，却令相干场对相位公式、幅度项和非局部耦合高度敏感。
4. **宽谱 PM 的单一尺度判据不充分**：总体 `kL` 和 RMS `h/L` 看似合格，不能
   推出每个尺度及其相互耦合都满足 phase-screen 等效。

因此，当前较大偏差应描述为：**并非 `kL`/`h/L` 的简单硬失败，而是一个只在
近法向、局部 Kirchhoff 条件下成立的 `2k eta` 相位屏，被用于宽入射角谱、强相位
调制且导数由高 q 主导的 PM 面后，遗漏了显著的有限角和非局部表面波动效应。**
现有数值证据不足以支持继续用局部经验相位/坡度修正；下一步若继续提升 PE
粗糙面模型，应以独立 Helmholtz/BIE 为基准审查非局部边界算子，而不是调节
`2k eta` 的经验系数。

## 依据

- B. J. Kachoyan and C. Macaskill, “Acoustic scattering from an arbitrarily rough
  surface,” *JASA* 82(5), 1720--1726 (1987),
  [DOI 10.1121/1.395165](https://doi.org/10.1121/1.395165)；
  [AIP/JASA issue entry](https://pubs.aip.org/asa/jasa/issue/82/5)。
- `reports/pe_bellhop_pm_canonical_mapper_report.md`
- `reports/pe_bellhop_controlled_comparison_stage0_status.md`
- `reports/pe_bellhop_controlled_comparison_stage1_status.md`（只采用 Stage-1Y 后的
  convention-fixed 结果）
- `reports/pe_bellhop_controlled_comparison_stage2_height_sweep_report.md`
- `reports/pe_bellhop_controlled_comparison_stage3_slope_curvature_report.md`
- `reports/pe_bellhop_controlled_comparison_stage4_validity_map_report.md`
- `reports/pe_bellhop_controlled_comparison_stage5_phase_attribution_report.md`
- `reports/pe_bellhop_controlled_comparison_stage6_fixed_pm_report.md`
- `reports/pe_bellhop_helmholtz_bie_final_report.md`
- `reports/pe_surface_operator_bie_G0_high_K_diagnostic.md`
- `reports/pe_surface_normal_approximation_G1_audit.md`
- `reports/pe_surface_kz_aware_G2_validation.md`
- `reports/pe_surface_angle_slope_G3_validation.md`
- `reports/pe_surface_operator_improvement_final_report.md`
- `results/validation/pe_bellhop_controlled_comparison/stage6_fixed_pm/fixed_pm_fourier_coefficients.csv`
