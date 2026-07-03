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
\Delta\phi(x,y)=k_0\Gamma\,\xi(x,y),
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

当前实现中，`pm_convolution` 的相位尺度与 Kirchhoff 相位屏和 SSA coherent reflection 保持同一法向定义：

\[
q_z^{\rm eff}=k_0(\cos\theta_i+\cos\theta_s).
\]

因此工程散射强度中使用 \((q_z^{\rm eff})^2\)，法向时为 \((2k_0)^2\)。旧的 \(2k_0(\cos\theta_i+\cos\theta_s)\) 口径会在法向时变成 \(4k_0\)，相当于重复计算双程相位因子；它现在只作为历史审计量记录，不再作为 `pm_convolution` 主公式。需要强调的是，这个修正只统一工程 baseline 的相位尺度，仍不把 `pm_convolution` 变成严格 SSA 散射截面模型。

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

### 8.1 镜面项与非相干项的接收贡献诊断

为了判断非相干散射项在当前设置下是否可以忽略，模型提供一个默认关闭的诊断分解。它不改变主反射场，只把同一个海面反射场写成

\[
\Psi_{\rm ref}(x,y)
=
\Psi_{\rm coh}(x,y)+\Psi_{\rm sca}(x,y),
\]

其中

\[
\Psi_{\rm coh}(x,y)=R_{\rm coh}\Psi_{\rm inc}(x,y),
\qquad
\Psi_{\rm sca}(x,y)=\Psi_{\rm ref}(x,y)-\Psi_{\rm coh}(x,y).
\]

然后分别把 \(\Psi_{\rm coh}\) 和 \(\Psi_{\rm sca}\) 从海面推进到接收深度，得到

\[
h_{\rm ref}^{\rm coh},\qquad
h_{\rm ref}^{\rm sca},\qquad
h_{\rm ref}^{\rm total}
\approx
h_{\rm ref}^{\rm coh}+h_{\rm ref}^{\rm sca}.
\]

这个诊断回答的是一个工程问题：在给定频率、风速、\(H_s\)、网格、随机种子和 \(C_{\rm sca}\) 下，接收点处的非相干反射分量相对于镜面相干分量有多大。常用判据为

\[
20\log_{10}
\frac{|h_{\rm ref}^{\rm sca}|}
{|h_{\rm ref}^{\rm coh}|}
\]

以及海面谱能量预算中的

\[
10\log_{10}
\frac{E_{\rm sca}}{E_{\rm coh}}.
\]

若二者都低于约 \(-20\) dB，通常可以把非相干项视为当前设置下的弱修正；若接近或高于 \(-10\) dB，则不应轻易忽略。这个阈值是工程判据，不是物理定理。

在固定 \(H_s\) 的风速扫描中，\(\sigma_\eta=H_s/4\) 和

\[
\iint W_\eta(K_x,K_y)\,dK_xdK_y
=
\sigma_\eta^2
\]

保持不变。此时改变风速 \(U\) 主要改变归一化 PM 谱的形状，也就是不同海面波数尺度的能量分布；它不表示海况强度随风速单调增强。

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

14. 镜面/非相干反射贡献对比使用固定 \(H_s\) 的风速-频率扫描。它分别传播相干镜面场和非相干散射场到接收点，比较 \(|h_{\rm ref}^{\rm sca}|/|h_{\rm ref}^{\rm coh}|\) 以及 \(E_{\rm sca}/E_{\rm coh}\)。该对比用于判断当前设置下非相干项是否可忽略，不改变主反射场或通信链路。

    在 \(H_s=0.5\) m、\(C_{\rm sca}=1\)、\(64\times64\) reduced-grid smoke 设置下，\(U=3,8\) m/s 和 \(f=4,8\) kHz 的测试点均显示非相干项远大于相干镜面项。因此这些设置下不应忽略非相干散射。该结论只对应当前粗糙度、网格和工程归一化常数；更弱海况或重新标定 \(C_{\rm sca}\) 后需要重新判断。

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
- 当前只实现了可选的 Broschat-style SSA2 coherent reflection coefficient 诊断；二阶非相干散射功率、完整 SSA2 T-matrix、NLSSA、高阶多次散射和实验标定海面散射截面仍未实现。
- 当前不向通信反射场注入倏逝散射分量。
- 当前统计散射 realization 的跨频率相关性只有 independent、共享随机样本和一阶自回归这类工程随机模型；尚未建立物理或经验标定的宽带频率相关模型。因此宽带 \(H_f\) 的频域连续性仍需要后续专门研究。
- 在固定 \(H_s\) 归一化下，改变风速 \(U\) 主要改变 PM 谱形状；不应简单解释为海况强度随风速单调增强。

因此，当前 `ssa_stat_kernel` 应理解为 PM-spectrum-driven first-order pressure-release Dirichlet SSA statistical reflection/scattering model。其中 `pm_convolution` 是工程基线，`ssa1_geometry` 是一阶 Dirichlet 几何核；它还不是完整的海面声散射理论闭环。

## 13. 入射与海面反射声场的可视化

为了直观比较声波到达海面前后的形状，当前实现可以在参考频点记录两段中心截面：

1. 入射段从发射深度 \(z_{\rm tx}\) 向上传播到海面 \(z=0\)；
2. 反射段从海面 \(z=0\) 向下传播到接收深度 \(z_{\rm rx}\)。

静态传播图显示的是 PE/WAPE 复包络的幅度

\[
20\log_{10}|\Psi(x,z)|,
\]

而不是某个瞬间的真实声压。入射和反射图使用同一个幅度参考，因此可以同时观察波束形状和反射衰减，不能把两幅图分别归一化后再比较强弱。

海面平面图比较

\[
|\Psi_{\rm inc}(x,y)|,\quad \arg\Psi_{\rm inc}(x,y),
\quad
|\Psi_{\rm ref}(x,y)|,\quad \arg\Psi_{\rm ref}(x,y),
\]

用于观察粗糙边界造成的幅度起伏和相位畸变。角谱图比较

\[
|\Psi_{\rm inc}(K_x,K_y)|^2
\quad\text{与}\quad
|\Psi_{\rm ref}(K_x,K_y)|^2,
\]

并给出横向波数 RMS 和 90% 能量半径。角谱变宽表示能量被重新分配到更大的横向波数，但该图仍是离散反射场诊断，不能解释为实验标定的双站散射截面。

为了显示波峰和波谷的运动，可以用名义声速 \(c_0\) 为复包络恢复载波：

\[
p_{\rm inc}(x,z,t)=
\Re\left\{
\Psi_{\rm inc}(x,z)
\exp\left[i k_0(z_{\rm tx}-z)-i\omega t\right]
\right\},
\]

\[
p_{\rm ref}(x,z,t)=
\Re\left\{
\Psi_{\rm ref}(x,z)
\exp\left[i k_0z-i\omega t\right]
\right\}.
\]

这里的动画只是在单频 PE 复包络上恢复一个载波周期，用于展示入射向上、反射向下的传播方向和波前结构；它不是额外的时域声场求解。压力释放自由海面的 \(R_0=-1\) 已包含在 \(\Psi_{\rm ref}\) 中，因此不应再人为增加一次相位反转。

该诊断默认关闭，仅在显式启用时保存参考频点的截面和海面二维复场。它不改变

\[
H(f)=H_{\rm dir}(f)+H_{\rm ref}(f)
\]

及其通信链路接口。

当前 reduced-grid 验证使用 \(64\times64\) 横向网格、\(z_{\rm tx}=20\) m、\(z_{\rm rx}=2\) m，并覆盖 `ssa1_geometry`、`kirchhoff_spatial` 和 \(4,6,8\) kHz 的宽带不变量检查。结果为：

- 开启与关闭该诊断时，\(H_f\)、\(H_{\rm dir}\) 和 \(H_{\rm ref}\) 的最大差异均为 0；
- 当 \(H_s=0\) 时，反射场满足压力释放边界的相位反转，\(\Psi_{\rm ref}\approx-\Psi_{\rm inc}\)，相对误差为 \(2.62\times10^{-16}\)；
- 平整海面反射前后幅度一致性的相对误差为 \(2.08\times10^{-16}\)；
- 宽带总信道代数不变量的最大误差为 \(1.55\times10^{-17}\)；
- 入射场海面端点、反射场海面起点和反射场接收端点与原传播结果的差异均为 0。

在默认展示算例 \(f=6000\) Hz、\(H_s=0.2\) m、`ssa1_geometry` 和 \(128\times128\) 网格下，入射角谱横向波数 RMS 为 \(3.3332\ {\rm rad/m}\)，反射角谱为 \(3.3665\ {\rm rad/m}\)；90% 能量半径分别为 \(5.0613\ {\rm rad/m}\) 和 \(5.1131\ {\rm rad/m}\)。这说明该次 realization 中存在轻微角谱展宽，但这些数值只对应当前网格、海况和随机种子，不能外推为普适散射规律。静态图和一个载波周期的 MP4 动画均已成功生成。
## 附：无 PE 的垂直平面波相干反射诊断

为单独检查海面边界的相干镜面反射强度，可以不调用 PE/WAPE，而令入射波为垂直平面波
\[
p_{\rm inc}=e^{ikz}.
\]
在海面 \(z=0\) 处，其横向复包络为常数：
\[
\Psi_{\rm inc}(x,y)=1.
\]

在这种设置下，SSA1 pressure-release / Dirichlet 相干反射为
\[
R_{\rm SSA1}=-\exp(-2k_0^2\sigma_\eta^2).
\]
如果不对 PM 谱做 \(H_s\) 归一化，则
\[
\sigma_\eta^2=\sum W_\eta(K_x,K_y)\Delta K_x\Delta K_y,
\qquad
H_{s,\rm raw}=4\sigma_\eta .
\]
因此风速 \(U\) 同时改变 PM 谱形状和积分粗糙度方差，不能与固定 \(H_s\) 的风速 sweep 混为一谈。

Kirchhoff realization 对照不是直接给出一个统计平均公式，而是先生成具体海面 \(\eta(x,y)\)，再形成相位屏
\[
G(x,y)=R_0\exp[i\Delta\phi(x,y)].
\]
垂直平面波的相干镜面分量是该相位屏的空间平均：
\[
R_{\rm K,coh}=\langle G(x,y)\rangle_{x,y}.
\]
脚本同时记录两种相位约定：\(\Delta\phi=2k_0\eta\) 的修正后法向双程高度相位，以及 legacy 诊断口径 \(\Delta\phi=2k_0\cdot2\eta\)。前者的 Gaussian coherent expectation 与 SSA1 法向公式一致；后者用于说明旧的额外 2 倍因子会带来更强的相干衰减。

修正后，\(\Delta\phi=2k_0\eta\) 是主 Kirchhoff 公式；\(\Delta\phi=4k_0\eta\) 只作为 `legacy_4k_eta` 对照。

需要区分三个量：\(\langle |G|^2\rangle\) 表示相位屏的总反射表面功率，对纯 pressure-release 相位屏可接近 1；单个 realization 的 \(|\langle G\rangle_{x,y}|\) 会受到有限孔径和随机样本残余影响；多 realization 下的 \(|\mathbb{E}[\langle G\rangle_{x,y}]|\) 才对应严格意义上的 ensemble coherent specular strength。粗糙海面可以在总反射功率基本不变的同时，使镜面相干项因相位抵消而显著下降。

## 附：Raw PM 风速驱动下的 Kirchhoff 两分支统一口径

当前项目有两种海面粗糙度幅度口径：

- `surface_roughness_scale_mode='target_hs'`：默认模式。PM 谱或海面 realization 会按 `sea_hs_target` 缩放，因此改变风速 \(U\) 主要改变 PM 谱形状，不直接表示海况强度随风速增强。
- `surface_roughness_scale_mode='raw_pm'`：风速驱动模式。不按 `sea_hs_target` 缩放，风速 \(U\) 直接决定 PM 谱积分方差和有效波高。

在 raw PM 模式下，当前统一采用项目显式海面使用的离散谱方差作为主口径：
\[
\sigma_{\eta,\rm raw}^2
=
\sum_{K_x,K_y}\Phi_{2D}(K_x,K_y)\Delta K_x\Delta K_y,
\qquad
H_{s,\rm raw}=4\sigma_{\eta,\rm raw}.
\]

因此，`kirchhoff_kdomain` 和 `kirchhoff_kstat` 的 `sigma_eta_raw_m`、`Hs_raw_m` 都应按上式解释。两者的对齐方式不同：

- `kirchhoff_kdomain` 显式生成具体海面 \(\eta(x,y)\)，raw PM 下记录的是该 realization 的 `std(eta(:))` 和对应 \(H_s\)。有限 seed 数下，样本标准差会围绕离散谱目标值波动。
- `kirchhoff_kstat` 不生成具体 \(\eta(x,y)\)，而是用高度谱构造相关函数。它内部使用连续谱归一化
\[
C_\eta(0)=
\frac{1}{(2\pi)^2}
\sum W_\eta(K_x,K_y)\Delta K_x\Delta K_y.
\]
为了让该式与项目离散 PM 方差一致，raw PM 下使用
\[
W_\eta=\Phi_{2D}(2\pi)^2,
\]
从而 \(C_\eta(0)=\sum\Phi_{2D}\Delta K_x\Delta K_y\)。

显式分支中还需要注意一个离散随机场细节：当前海面 realization 由 unconstrained complex spectrum 经 `real(ifft2(...))` 得到。对这种构造，取实部会使方差约损失一半，因此 raw PM 生成时需要乘以 \(\sqrt{2}\)，使显式海面的样本方差回到目标离散谱方差。这个修正只用于保证 raw PM 方差口径正确；在 `target_hs` 模式下，后续按目标 \(H_s\) 的缩放会抵消该差异，默认目标波高结果不应因此改变。

最近一次 reduced-grid 对比使用
`wind_list=[3 5 8 10 12 15]`、`seed_count=32`、`nx=ny=128`、\(f_0=6000\) Hz。结果显示：

- \(H_{s,\rm raw}\) 的两分支相对误差最大约为 \(1.45\times10^{-2}\)，说明海况强度已经基本对齐；
- `kirchhoff_kstat` 全 \(K\) 相位屏能量闭合误差最大约为 \(2.22\times10^{-15}\)；
- 宽带/标量信道代数不变量 \(H_f=H_{\rm direct,f}+H_{\rm reflect,f}\) 保持在约 \(10^{-17}\) 量级。

归一化后，两个分支的主要差异不再是 \(H_s\) 口径，而是模型本身：`kirchhoff_kdomain` 是具体海面相位屏 realization，接收端反射幅度包含有限孔径和 speckle 波动；`kirchhoff_kstat` 是统计相位屏模型，直接生成相干反射项、非相干散射功率谱和可复现随机散射 realization。因此二者不应要求逐点一致，更合理的比较对象是 ensemble 趋势、相干衰减、非相干能量以及接收端统计量。

## 附：Kirchhoff 相位屏法向因子修正

当前 Kirchhoff realization 相位屏采用
\[
\Delta\phi(x,y)=k_0(\cos\theta_i+\cos\theta_r)\eta(x,y).
\]
法向入射和镜面反射时
\[
\cos\theta_i=\cos\theta_r=1,
\]
因此
\[
\Delta\phi=2k_0\eta.
\]

这与垂直平面波的双程高度路径差一致：海面高度扰动为 \(\eta\) 时，入射-反射的等效路径差为 \(2\eta\)，相位扰动为 \(k_0\cdot2\eta\)。旧实现等效使用
\[
\Delta\phi=2k_0(\cos\theta_i+\cos\theta_r)\eta,
\]
在法向时变成 \(4k_0\eta\)，会使 Kirchhoff coherent average 出现类似 \(\exp(-8k_0^2\sigma_\eta^2)\) 的过强相干衰减。该旧口径现在只保留在独立平面波诊断脚本中作为 `legacy_4k_eta` 对照，不再作为主模型公式。
