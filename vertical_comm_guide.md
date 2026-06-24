# 垂直水声通信与海面统计散射模型说明

本文说明当前垂直水声通信模型中的传播、海面反射/散射和通信信道关系。它面向阅读模型的研究或工程人员，重点放在物理含义、数学形式、近似边界和验证结论，而不是脚本接口或运行日志。

当前公开信道入口为 `vertical_channel_model.m`，垂直 WAPE/PE 推进核心为 `vertical_wape_propagator.m`，海面边界模型集中在 `pm_surface_boundary_model.m`。这些文件名只用于定位实现；下文主要使用物理量和公式描述。

## 1. 问题背景

模型描述海底发射端到近海面接收端的上行水声通信链路。总频域信道由直达传播和海面反射/散射贡献叠加：

\[
H(f)=H_{\rm dir}(f)+H_{\rm ref}(f).
\]

在参考频率 \(f_{\rm ref}\) 处，窄带等效信道为

\[
h_{\rm total}=H(f_{\rm ref}),\qquad
h_{\rm dir}=H_{\rm dir}(f_{\rm ref}),\qquad
h_{\rm ref}=H_{\rm ref}(f_{\rm ref}).
\]

后续通信链路只消费最终的 \(H(f)\)。因此海面模型改变的是 \(H_{\rm ref}(f)\) 以及由此导致的总信道变化，不改变调制、噪声、同步、均衡和判决流程。

## 2. 坐标与传播约定

海面定义为

\[
z=0,
\]

深度方向向下为正。发射端深度为 \(z_{\rm tx}\)，接收端深度为 \(z_{\rm rx}\)，上行链路满足

\[
0\le z_{\rm rx}<z_{\rm tx}.
\]

直达项从 \(z_{\rm tx}\) 推进到 \(z_{\rm rx}\)。反射项先从 \(z_{\rm tx}\) 推进到海面，再由海面边界模型给出反射/散射场，最后从海面推进回 \(z_{\rm rx}\)。

## 3. PM 海面高度谱

粗糙海面由 Pierson-Moskowitz, PM, 谱描述。令横向波数

\[
K=\sqrt{K_x^2+K_y^2}.
\]

一维 PM 谱写作

\[
E_{1D}(K)
=
\frac{\alpha_{\rm PM}}{2K^3}
\exp\left(
-\beta_{\rm PM}\frac{g^2}{U^4K^2}
\right),
\qquad K>0.
\]

对应的各向同性二维高度谱为

\[
\Phi_{2D}(K_x,K_y)
=
\frac{E_{1D}(K)}{2\pi K},
\qquad K>0.
\]

实际使用时，谱会按目标有效波高 \(H_s\) 归一化。海面高度标准差取

\[
\sigma_\eta=\frac{H_s}{4},
\]

并要求

\[
\iint W_\eta(K_x,K_y)\,dK_xdK_y=\sigma_\eta^2.
\]

这里 \(W_\eta\) 表示海面高度的波数谱。它说明不同尺度、不同方向的海面起伏各有多少方差贡献，而不是只表示“有没有”某个波纹尺度。

当 \(H_s=0\) 时，

\[
\sigma_\eta=0,\qquad W_\eta=0.
\]

这表示海面没有统计起伏，模型退化为平整自由海面反射。

## 4. Kirchhoff 空间海面模型

`kirchhoff_spatial` 是默认海面模型。它先从 PM 谱生成一次具体随机海面

\[
\xi(x,y),
\]

再用相位屏近似描述海面高度起伏导致的反射相位扰动：

\[
\Psi_{\rm ref}(x,y)
=
R_0\,\Psi_{\rm inc}(x,y)
\exp\left(i\Delta\phi(x,y)\right).
\]

其中 \(R_0\) 是平整海面的反射系数。压力释放自由海面的默认近似为

\[
R_0=-1.
\]

相位扰动写作

\[
\Delta\phi(x,y)=2k_0\Gamma\,\xi(x,y),
\qquad
k_0=\frac{2\pi f}{c_0}.
\]

\(\Gamma\) 是入射和反射方向共同决定的有效垂向相位因子。法向入射近似下，\(\Gamma\) 可理解为固定几何因子；斜入射修正时，它由局部传播方向决定。

该模型的特点是直观，并对应某一次具体海面 realization。它会通过随机相位屏自然产生角谱展宽和 speckle 式起伏，但它不是从海面谱直接给出平均散射功率的闭式统计核。

## 5. PM 谱驱动的 SSA-like 统计散射分支

`ssa_stat_kernel` 分支不生成具体 \(\xi(x,y)\)。它直接由 \(W_\eta(K_x,K_y)\) 构造波数域统计散射功率，再合成海面反射场。

令入射场横向谱为

\[
\Psi_{\rm inc}(K),\qquad
P_{\rm inc}(K)=|\Psi_{\rm inc}(K)|^2.
\]

反射场分为相干镜面项和非相干散射项：

\[
\Psi_{\rm ref}(K)=\Psi_{\rm coh}(K)+\Psi_{\rm sca}(K).
\]
* 即使主传播方向近似垂直，只要海面是粗糙的，反射就不可能只有一个完美镜面分量。
* 粗糙海面让镜面反射变弱，缺失的能量一部分来自非相干散射
### 5.1 相干镜面项

粗糙度会降低镜面相干反射。当前采用

\[
R_{\rm coh}
=
R_0
\exp\left[
-\frac12(\gamma_i+\gamma_s)^2\sigma_\eta^2
\right].
\]

在镜面方向 \(\gamma_s=\gamma_i\)，因此

\[
R_{\rm coh}=R_0\exp[-2\gamma_i^2\sigma_\eta^2].
\]

法向入射时 \(\gamma_i=k_0\)，pressure-release / Dirichlet 自由海面的 \(R_0=-1\)，于是

\[
R_{\rm coh}=-\exp[-2k_0^2\sigma_\eta^2].
\]

Broschat 1993 的 PM 海面 SSA coherent reflection 主要用于核对这个相干镜面反射系数，而不是下面的非相干散射功率 \(P_{\rm sca}\)。也就是说，Broschat 型结果约束的是镜面相干项随粗糙度和频率衰减的趋势；非相干散射功率分布需要由一阶 Dirichlet 散射核单独给出。

当前默认仍采用上述一阶 coherent reflection。为了评估二阶 coherent correction 是否值得纳入，本项目还提供一个可选的 Broschat-style 二阶相干反射系数对比：

\[
R_2
=
R_1
+
2\gamma_i\exp[-2\gamma_i^2\sigma_\eta^2]
\iint
W_\eta(q)\,[\gamma(K_i+q)-\gamma_i]\,dq_xdq_y .
\]

其中

\[
R_1=-\exp[-2\gamma_i^2\sigma_\eta^2],
\qquad
\gamma(K)=\sqrt{k_0^2-|K|^2}.
\]

数值实现中，\(\gamma(K)\) 使用复数平方根分支：

\[
\gamma(K)=\sqrt{\mathrm{complex}(k_0^2-|K|^2,0)}.
\]

因此传播分量给出正实根，倏逝分量给出正虚根。该二阶项只用于比较 coherent reflection coefficient；它不是二阶非相干散射功率 \(P_{\rm sca}\)，也不是完整 SSA2 T-matrix 或 NLSSA 模型。当前主测试采用法向入射 \(K_i=0\)、\(\gamma_i=k_0\)；斜入射时只使用有效二维近似，结论应以法向或近垂直几何为准。

当前 reduced-grid 对比显示，在 \(H_s\le 0.5\) m 的弱到中等海况下，SSA2 coherent loss 与 SSA1 的差异不超过约 \(0.1\) dB；在 \(H_s=1\) m、4000 Hz 的强海况低频点，差异约为 \(0.37\) dB。因此本项目暂把 `ssa2_broschat_coherent` 作为诊断对照，而不是默认模型。是否采用它应取决于目标海况、频段和后续是否实现完整二阶非相干散射。

于是

\[
\Psi_{\rm coh}(x,y)=R_{\rm coh}\Psi_{\rm inc}(x,y).
\]

这个指数项可理解为 SSA coherent reflection 的统计平均。它和 Kirchhoff realization 相位屏中的单次高度相位扰动不是同一个量：Kirchhoff 相位屏先生成具体 \(\xi(x,y)\)，再给每个空间点加相位；SSA coherent reflection 则直接对随机海面统计平均后得到 \(R_{\rm coh}\)。若海面高度近似为零均值高斯随机变量，反射相位扰动的方差越大，不同海面 realization 的镜面相干叠加越容易相互抵消。

当 \(H_s=0\) 时，\(\sigma_\eta=0\)，因此

\[
R_{\rm coh}=R_0.
\]

压力释放平整自由海面下即为 \(R_{\rm coh}=-1\)。

### 5.2 工程基线核：`pm_convolution`

`pm_convolution` 是工程基线核。它把非相干散射理解为：入射方向 \(K'\) 的谱能量，如果海面中存在波数差 \(K-K'\) 对应尺度的起伏，就可以被重新分配到散射方向 \(K\)。

其散射核为

\[
S_{\rm PM}(K,K')\propto C_{\rm sca}W_\eta(K-K').
\]

对应原始散射功率为

\[
P_{\rm sca}^{\rm raw}(K)
=
C_{\rm sca}
\sum_{K'}
W_\eta(K-K')P_{\rm inc}(K')
\Delta K_x\Delta K_y.
\]

这里 \(W_\eta(K-K')\) 不只是判断“有没有合适波纹”，而是给出该波数差处海面高度起伏的统计强度。数值越大，说明该尺度和方向的海面起伏越强，对应的散射贡献越大。

\(\Delta K_x\Delta K_y\) 是离散波数网格上的面积权重。连续形式中散射功率来自积分

\[
\iint W_\eta(K-K')P_{\rm inc}(K')\,dK'_xdK'_y,
\]

离散求和时需要用每个小网格面积把积分近似为求和。

\(C_{\rm sca}\) 是工程归一化常数，用于数值敏感性分析；它不是实验标定的绝对散射截面常数。

### 5.3 一阶 Dirichlet 几何核：`ssa1_geometry`

`ssa1_geometry` 在工程卷积基础上加入压力释放 Dirichlet 自由海面的一阶几何因子。对横向波数 \(K\)，定义传播垂向波数

\[
\gamma(K,f)=\sqrt{\max(k_0^2-|K|^2,0)}.
\]

当前只把 \(|K|\le k_0\) 的传播分量纳入散射功率。直观地说，这些分量可以作为远场传播波携带能量离开海面；\(|K|>k_0\) 的倏逝分量在垂向上不形成传播功率流，当前原型不把它们作为通信反射能量注入。

压力释放 Dirichlet 条件下，一阶几何因子为

\[
G_{\rm SSA1}(K,K';f)=4\gamma(K,f)\gamma(K',f).
\]

于是统计核为

\[
S_{\rm SSA1}(K,K';f)
\propto
G_{\rm SSA1}(K,K';f)W_\eta(K-K').
\]

原始散射功率为

\[
P_{\rm sca}^{\rm raw}(K)
=
C_{\rm sca}
\sum_{K'}
G_{\rm SSA1}(K,K';f)
W_\eta(K-K')
P_{\rm inc}(K')
\Delta K_x\Delta K_y.
\]

等价地，令

\[
A(K')=\gamma(K',f)P_{\rm inc}(K'),
\]

则

\[
P_{\rm sca}^{\rm raw}(K)
=
4C_{\rm sca}\gamma(K,f)
\left[W_\eta*A\right](K)
\Delta K_x\Delta K_y.
\]

这里的卷积表示对所有入射横向波数 \(K'\) 的贡献求和。快速 FFT 形式和显式逐项求和形式在周期卷积假设下是同一个离散公式的两种计算方式。

当前 `ssa1_geometry` 只支持压力释放 Dirichlet 边界，即 \(R_0=-1\)。任意阻抗边界、Neumann 边界或从一般反射系数到 SSA 几何项的映射尚未实现。

## 6. 周期卷积与 zero-padding 线性卷积

当前实现支持两种卷积方式。

`periodic` 是默认兼容路径。它把有限波数网格看作周期延拓，使用 FFT 快速计算周期卷积。`ssa1_debug_dense` 在小网格上逐项显式求和，用来验证这个周期卷积公式和快速 FFT 结果一致。

`zero_padded` 是线性卷积审计路径。它先把谱搬到带符号波数顺序，在扩展网格上做零填充 FFT 线性卷积，再裁剪回原始波数范围。这样可以降低周期折返 aliasing 风险，更适合检查有限波数窗口边界对散射功率的影响。

两种方式的物理含义不同。周期卷积假设超出网格边界的波数会从另一侧折回；zero-padding 则把网格外未表示的谱当作零处理。因此二者不要求逐点相同，但都必须满足能量审计和接口不变量。

## 7. 能量约束

非相干散射项是随机生成的。为避免随机散射导致非物理放大，模型施加能量约束。

定义

\[
E_{\rm inc}=\sum_K|\Psi_{\rm inc}(K)|^2,
\]

\[
E_{\rm coh}=\sum_K|\Psi_{\rm coh}(K)|^2,
\]

\[
E_{\rm sca}^{\rm raw}=\sum_KP_{\rm sca}^{\rm raw}(K).
\]

允许注入的散射能量不超过

\[
E_{\rm sca}^{\rm max}=\max(E_{\rm inc}-E_{\rm coh},0).
\]

因此实际使用的散射功率为

\[
P_{\rm sca}(K)
=
P_{\rm sca}^{\rm raw}(K)
\min\left(
1,
\frac{E_{\rm sca}^{\rm max}}
{\max(E_{\rm sca}^{\rm raw},\varepsilon)}
\right).
\]

从而保证

\[
E_{\rm coh}+E_{\rm sca}\le E_{\rm inc}.
\]

这是数值稳定和物理合理性约束，不代表已经完成绝对散射截面的实验标定。

## 8. 随机散射谱

若启用随机散射，非相干散射场按复高斯随机相位生成：

\[
\Psi_{\rm sca}(K)=\sqrt{P_{\rm sca}(K)}\,z(K),
\qquad
z(K)\sim\mathcal{CN}(0,1).
\]

它存在的原因是：\(P_{\rm sca}(K)\) 只给出每个散射方向上的平均功率预算，而通信仿真需要一条具体的复数信道 realization。随机散射谱把统计功率变成一个可传播、可反变换、可被通信链路消费的具体反射场样本。

最终

\[
\Psi_{\rm ref}(K)=\Psi_{\rm coh}(K)+\Psi_{\rm sca}(K).
\]

若关闭随机散射，模型仍计算 \(P_{\rm sca}\) 和能量审计量，但不把 \(\Psi_{\rm sca}\) 注入反射场。此时通信链路只接收相干镜面反射项。

宽带频率轴上，随机散射 realization 还需要说明跨频率相关性。默认处理仍是各频点独立生成复高斯随机谱，这保持了早期实现的兼容性。当前还提供两个仅用于随机信道生成诊断的工程相关模式：一种是在不同频点复用同一个随机相位样本，另一种是在频率序列上使用一阶自回归相关样本。它们的作用是检查宽带 \(H_f\) 的频域连续性对随机散射 realization 的敏感性；它们不是新的 SSA 物理散射公式，也不是经过实验标定的海面时间频率相关模型。

## 9. 通信链路解释

海面模型只改变频域信道 \(H(f)\)。后续 MPSK 通信仍沿用同一链路：

\[
H(f)\rightarrow H_{\rm baseband}(f)\rightarrow h_{\rm bb}(t),
\]

然后进行符号卷积、噪声注入、同步、均衡和判决。

因此比较 `kirchhoff_spatial`、`pm_convolution` 和 `ssa1_geometry` 时，应主要观察

\[
|H(f)|,\qquad |H_{\rm ref}(f)|,\qquad \arg H(f),
\]

以及同一通信流程下的 BER/SER 趋势。

## 10. 当前验证结论

当前 reduced-grid 数值验证支持以下结论：

1. 当 \(H_s=0\) 时，统计海面模型退化为平整自由海面反射：

   \[
   \sigma_\eta=0,\qquad W_\eta=0,\qquad R_{\rm coh}=R_0,\qquad E_{\rm sca}=0.
   \]

2. `pm_convolution` 和 `ssa1_geometry` 在 \(H_s=0\) 时都与平整 `kirchhoff_spatial` 响应达到舍入误差级一致。

3. 非零海况下能量审计满足

   \[
   E_{\rm coh}+E_{\rm sca}\le E_{\rm inc}.
   \]

4. 相同随机种子下统计散射结果可复现；改变海面随机种子会改变反射散射 realization，但直达项保持稳定。

5. `random_scatter=false` 时，模型仍保留散射功率和能量统计，但不会向通信链路注入随机散射场。

6. `ssa1_debug_dense` 的小网格显式求和与 `ssa1_geometry` 的周期 FFT 卷积结果一致，验证了一阶 Dirichlet 几何核的离散实现。

7. `zero_padded` 线性卷积路径已经可运行，可用于检查周期卷积 aliasing 风险；它与周期卷积不要求逐点相等。

8. 增大 \(H_s\) 时，相干镜面项下降，散射功率预算增强。增大 \(C_{\rm sca}\) 时，原始散射功率单调增强，但最终注入散射能量仍受能量约束限制。

9. SSA 物理趋势验证覆盖 \(H_s=0,0.05,0.2,0.5,1.0\) 和 \(f=4,6,8,10\ {\rm kHz}\)。在法向镜面基准下，`pm_convolution` 与 `ssa1_geometry` 都应满足

   \[
   R_{\rm coh}=R_0\exp(-2k_0^2\sigma_\eta^2).
   \]

   该验证同时检查固定频率下 \(H_s\) 增大时 \(|R_{\rm coh}|\) 与 \(E_{\rm coh}/E_{\rm inc}\) 不增，以及固定非零 \(H_s\) 下频率升高时 \(|R_{\rm coh}|\) 不增、相干指数更负。

10. 新增的 Kirchhoff/SSA1 统计对照不是证明两者逐点等价。`kirchhoff_spatial` 是具体海面 realization 相位屏模型；`ssa1_geometry` 是 PM 谱驱动的一阶 pressure-release Dirichlet 统计反射/散射模型。二者应在统计趋势上互相支持，例如粗糙度增强时 coherent loss 增强、反射角谱展宽增强、反射通道幅相波动增强。

    对照中共同观察

    \[
    |h_{\rm ref}|,\qquad |h_{\rm total}|,\qquad
    \left|\frac{h_{\rm ref}}{h_{\rm dir}}\right|,
    \]

    以及反射角谱的 RMS 横向波数和高波数能量比例。Kirchhoff 的展宽指标来自具体 realization 的反射谱；SSA1 的展宽指标来自统计散射合成后的反射谱。有限 seed 下不要求每个样本或每个相邻海况严格单调。

11. `surface_ssa_scatter_scale` 的标定被作为工程诊断处理。当前 reduced-grid 标定使用 `kirchhoff_spatial` 的多 seed 统计作为 realization-based reference，比较 \(E_{\rm sca}/E_{\rm inc}\)、反射角谱展宽、\(|h_{\rm ref}|\) 和 \(|h_{\rm ref}/h_{\rm dir}|\) 等统计量，给出候选 \(C_{\rm sca}\) 的相对误差。该标定不把 \(C_{\rm sca}\) 解释为实验绝对散射截面。

12. 跨频率相关随机散射验证表明，默认 independent 模式保持原有行为；相关模式只改变随机非相干反射样本，不改变直达项，也不改变
    \[
    H(f)=H_{\rm dir}(f)+H_{\rm ref}(f)
    \]
    的接口语义。它们用于后续宽带随机信道生成研究，而不是替代当前的海面散射核。

13. 新增适用性汇总只读取已有 reduced-grid 结果，汇总 SSA1/SSA2 coherent loss 差异、periodic/zero-padding 差异、Kirchhoff 标定残差和频率相关性诊断。该表是工程使用指南，不是完整物理适用域图。

## 11. 周期卷积与 zero-padding 对比

当前 reduced-grid 对比覆盖 \(H_s=0,0.05,0.2\)，并分别检查 `pm_convolution` 与 `ssa1_geometry`。

主要观察量包括

\[
\frac{E_{\rm sca}^{\rm raw}}{E_{\rm inc}},\qquad
\frac{E_{\rm sca}^{\rm limited}}{E_{\rm inc}},\qquad
\frac{E_{\rm ref}}{E_{\rm inc}},
\]

\[
|h_{\rm ref}|,\qquad |h_{\rm total}|,
\]

以及散射角谱和反射角谱的横向波数 RMS 展宽。

结果显示：

1. 当 \(H_s=0\) 时，`periodic` 与 `zero_padded` 都退化到同一个平整海面结果，与 `kirchhoff_spatial` 的差异约为 \(7.76\times 10^{-17}\)。

2. 非零海况下，`zero_padded` 与 `periodic` 的原始散射能量差异较小。在当前 reduced 设置中，
   `pm_convolution` 的 \(E_{\rm sca}^{\rm raw}/E_{\rm inc}\) 相对差异约为 \(4.90\times 10^{-4}\)，
   `ssa1_geometry` 约为 \(4.64\times 10^{-4}\)。

3. 反射频响差异比散射能量差异更明显，因为随机散射谱经过反射路径传播后会改变相位叠加。在弱海况 \(H_s=0.05\) 下，当前单 seed 对比中 `pm_convolution` 的 \(H_{\rm ref}\) 相对差异约为 \(3.19\times 10^{-3}\)，`ssa1_geometry` 约为 \(3.09\times 10^{-3}\)。

4. 能量审计仍通过。当前对比中总信道不变量保持在舍入误差级，能量守恒误差为 0。

因此，`periodic` 与 `zero_padded` 的判断标准不是逐点相等，而是接口不变量、能量审计和趋势是否可解释。二者差异反映的是有限波数窗口和周期折返假设对统计散射功率分配的影响。

## 12. 当前限制

当前模型仍有明确边界：

- `pm_convolution` 是工程统计基线核，不是严格 SSA 几何核。
- `ssa1_geometry` 只覆盖压力释放 Dirichlet 自由海面的一阶几何因子。
- \(C_{\rm sca}\) 是工程归一化常数，不是实验标定的绝对散射截面。
- 周期卷积是默认兼容路径；zero-padding 线性卷积是 aliasing 审计路径。
- 任意阻抗边界、Neumann 边界和一般反射系数到 SSA 几何项的映射尚未实现。
- second-order SSA coherent correction、NLSSA、高阶多次散射、full T-matrix solver 和实验标定海面散射截面尚未实现。
- 当前不向通信反射场注入倏逝散射分量。
- 当前统计散射 realization 的跨频率相关性只有 independent、共享随机样本和一阶自回归这类工程随机模型；尚未建立物理或经验标定的宽带频率相关模型。因此宽带 \(H_f\) 的频域连续性仍需要后续专门研究。
- 在固定 \(H_s\) 归一化下，改变风速 \(U\) 主要改变 PM 谱形状；不应简单解释为海况强度随风速单调增强。

因此，当前 `ssa_stat_kernel` 应理解为 PM-spectrum-driven first-order pressure-release Dirichlet SSA statistical reflection/scattering model。其中 `pm_convolution` 是工程基线，`ssa1_geometry` 是一阶 Dirichlet 几何核；它还不是完整的海面声散射理论闭环。
