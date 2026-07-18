# 垂直水声通信与海面统计散射模型说明

本文说明当前垂直水声通信模型中的传播、海面反射/散射和通信信道关系。它面向阅读模型的研究或工程人员，重点放在物理含义、数学形式、近似边界和验证结论，而不是脚本接口或运行日志。

当前公开信道入口为 `vertical_channel_model.m`，垂直 WAPE/PE 推进核心为 `vertical_wape_propagator.m`，海面边界模型集中在 `pm_surface_boundary_model.m`。这些文件名只用于定位实现；下文主要使用物理量和公式描述。

## 当前优化 joint-frequency 验证说明

F=64 的独立验证构建器把解析协方差和伪协方差展开为
`t=C_eta/sigma_eta^2` 的收敛级数，再执行标量空间 FFT、全局增广频率基和
分块 K/-K 因子分解。该变换保留物理 joint 统计定义，并非 shared seed、
shared phase 或 AR(1) 的替代方案。

U=8 m/s 的审计结果选择 150 m/384² PM 网格，并保持与 PE 相同的采样间隔。
该网格给出隐含 `Hs=1.3636 m`、离散/无限域捕获率 `0.99780` 和
`Kpeak/Kmin=2.570`，之后通过中央裁剪映射到 50 m/128² PE 网格，不逐条
realization 重归一化。离散条件库接口只接受已经验证的精确风速节点，不做
静默插值。U=8 的 128/128 接收端验证已经通过：joint reflected-only
PDP/LFM correlation 为 0.9809/0.9911，协方差误差 0.3554，明显优于
independent 的 0.9703，并接近 kdomain split-sample floor 0.3403。
U=8 properness 未被拒绝，默认保留 full rank 42；两节点库已经完成精确节点
抽样和未知风速拒绝测试。完整结果见
`reports/u8_joint_optimization_two_node_library_report.md`。

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

## 3. 发射声场初始条件

当前 PE/WAPE 信道模型中的初始发射声场不是理想数学点源，也不是已经调制好的时域 PSK 波形。代码在发射深度 \(z=z_{\rm tx}\) 的横向平面上给定一个二维高斯复包络：

\[
\Psi_{\rm tx}(x,y;z_{\rm tx})
=
\exp\left[
-\frac{(x-x_{\rm tx})^2+(y-y_{\rm tx})^2}{2\sigma_{\rm src}^2}
\right].
\]

这里 \((x_{\rm tx},y_{\rm tx},z_{\rm tx})\) 是发射端位置，\(\sigma_{\rm src}\) 是横向高斯源宽度。默认参数为

\[
x_{\rm tx}=0,\qquad
y_{\rm tx}=0,\qquad
z_{\rm tx}=100\ {\rm m},\qquad
\sigma_{\rm src}=0.3\ {\rm m}.
\]

该初始包络在中心点幅度为 1，初始相位为 0，随后对每个频率 \(f\) 用

\[
k_0=\frac{2\pi f}{c_0}
\]

进入 PE/WAPE 上行传播。也就是说，声学传播部分处理的是频域/窄带复包络 \(\Psi\)，真实声压可理解为复包络乘以载波相位后取实部。对入射上行场，文档和可视化中采用的约定是

\[
p_{\rm inc}(x,y,z,t)
=
\Re\left\{
\Psi_{\rm inc}(x,y,z)
\exp\left[i k_0(z_{\rm tx}-z)-i\omega t\right]
\right\}.
\]

通信脚本中的 MPSK 符号不是用来重新定义这个空间初始声场，而是在得到信道频响 \(H(f)\) 或基带冲激响应 \(h_{\rm bb}\) 后，通过乘法或卷积作用到符号序列上。因此，“声学初始发射场”对应上面的高斯空间包络；“通信发射信号”对应后续的 PSK 符号流，两者处在不同层级。

## 4. PM 海面高度谱

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

实际使用时有两种粗糙度幅度口径，由 `surface_roughness_scale_mode` 控制。

默认 `target_hs` 模式会按目标有效波高 \(H_s\) 归一化。海面高度标准差取

\[
\sigma_\eta=\frac{H_s}{4},
\]

并要求

\[
\iint W_\eta(K_x,K_y)\,dK_xdK_y=\sigma_\eta^2.
\]

这里 \(W_\eta\) 表示海面高度的波数谱。它说明不同尺度、不同方向的海面起伏各有多少方差贡献，而不是只表示“有没有”某个波纹尺度。

`raw_pm` 模式不按 `sea_hs_target` 缩放，风速 \(U\) 直接决定 PM 谱积分方差。当前项目统一采用显式海面使用的离散谱方差作为 raw PM 主口径：

\[
\sigma_{\eta,\rm raw}^2
=
\sum_{K_x,K_y}\Phi_{2D}(K_x,K_y)\Delta K_x\Delta K_y,
\qquad
H_{s,\rm raw}=4\sigma_{\eta,\rm raw}.
\]

因此固定 \(H_s\) 的风速扫描和 raw PM 风速扫描含义不同：前者主要改变 PM 谱形状，后者同时改变谱形状和海况强度。

当 \(H_s=0\) 时，

\[
\sigma_\eta=0,\qquad W_\eta=0.
\]

这表示海面没有统计起伏，模型退化为平整自由海面反射。

## 5. Kirchhoff 边界模型

当前 Kirchhoff 类分支有两条路径：显式相位屏路径 `kirchhoff_spatial` / `kirchhoff_kdomain`，以及统计相位屏路径 `kirchhoff_kstat`。前者生成具体海面 realization，后者不生成具体海面，而是从 PM 谱直接生成相干项和非相干统计散射谱。

### 5.1 显式相位屏：`kirchhoff_spatial` 与 `kirchhoff_kdomain`

`kirchhoff_spatial` 是默认海面模型。它先从 PM 谱生成一次具体随机海面

\[
\eta(x,y),
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
\Delta\phi(x,y)=k_0\Gamma\,\eta(x,y),
\qquad
k_0=\frac{2\pi f}{c_0}.
\]

\(\Gamma\) 是入射和反射方向共同决定的有效垂向相位因子。当前主公式在法向入射和镜面反射时给出

\[
\Delta\phi=2k_0\eta.
\]

这对应入射-反射双程高度路径差 \(2\eta\)。旧的 \(4k_0\eta\) 口径只保留在历史诊断中，不再作为主模型公式。

`kirchhoff_kdomain` 是同一显式相位屏思想的波数域接口，便于和统计卷积路径比较。显式分支的特点是直观，并对应某一次具体海面 realization。它会通过随机相位屏自然产生角谱展宽和 speckle 式起伏，但单次 realization 不是从海面谱直接给出 ensemble 平均散射功率的闭式统计核。

在 `raw_pm` 模式下，显式海面 realization 记录的是 `std(eta(:))` 和对应 \(H_{s,\rm raw}\)。由于当前 realization 由 unconstrained complex spectrum 经 `real(ifft2(...))` 得到，取实部会使方差约损失一半，因此 raw PM 生成时乘以 \(\sqrt{2}\)，使样本方差回到目标离散谱方差。有限 seed 数下，显式分支的 \(H_{s,\rm raw}\) 仍会围绕离散谱目标值波动。

### 5.2 统计相位屏：`kirchhoff_kstat`

`kirchhoff_kstat` 是 Kirchhoff / 统计相位屏分支。它不显式生成某个具体 \(\eta(x,y)\)，而是由 PM 高度谱构造海面高度相关函数、相干反射项、非相干相位屏谱和可选随机反射场 realization。该分支独立于 `ssa_stat_kernel`，当前作为粗糙海面非相干散射的主统计生成路径。

第一版采用近垂直 Kirchhoff 相位屏

\[
G(\mathbf r)=\exp(i\alpha\eta(\mathbf r)),
\qquad
\alpha=2k_0,
\qquad
R_0=-1.
\]

由高度谱得到相关函数

\[
C_\eta(\boldsymbol\rho)
=
\frac{1}{(2\pi)^2}
\int W_\eta(\mathbf K)
e^{i\mathbf K\cdot\boldsymbol\rho}\,d\mathbf K,
\qquad
\sigma_\eta^2=C_\eta(0).
\]

离散实现采用

\[
C_\eta=\mathrm{ifft2}(W_\eta)\,N_xN_y
\frac{\Delta K_x\Delta K_y}{(2\pi)^2}.
\]

因此在 `raw_pm` 模式下，为了让 \(C_\eta(0)\) 与项目离散 PM 方差一致，K-Stat 使用

\[
W_\eta=\Phi_{2D}(2\pi)^2,
\]

从而

\[
C_\eta(0)=
\sum \Phi_{2D}\Delta K_x\Delta K_y.
\]

相位屏均值和相干反射为

\[
\langle G\rangle
=
\exp\left(-\frac12\alpha^2\sigma_\eta^2\right),
\qquad
R_{\rm coh}=R_0\langle G\rangle.
\]

压力释放自由海面、近法向条件下即

\[
R_{\rm coh}=-\exp(-2k_0^2\sigma_\eta^2).
\]

相干反射场为

\[
\hat\Psi_{\rm coh}(\mathbf K)=R_{\rm coh}\hat\Psi_{\rm inc}(\mathbf K).
\]

非相干涨落相关函数为

\[
C_{\delta G}(\boldsymbol\rho)
=
\exp(-\alpha^2\sigma_\eta^2)
\left[
\exp(\alpha^2C_\eta(\boldsymbol\rho))-1
\right].
\]

对它做傅里叶变换得到非相干相位屏谱 \(S_{\delta G}\)。注意这里使用的是 \(S_{\delta G}\)，不把包含相干 delta 尖峰的总相位屏谱当作非相干散射谱。

非相干散射功率由统计卷积给出：

\[
P_{\rm sca}(\mathbf K_s)
=
\frac{|R_0|^2}{(2\pi)^2}
\int
S_{\delta G}(\mathbf K_s-\mathbf K_i)
|\hat\Psi_{\rm inc}(\mathbf K_i)|^2\,d\mathbf K_i.
\]

若 `surface_kstat_random_scatter=true`，随机反射谱 realization 为

\[
\hat\Psi_{\rm sca}^{(m)}(\mathbf K)
=
\sqrt{P_{\rm sca}(\mathbf K)}Z_m(\mathbf K),
\qquad
Z_m\sim\mathcal{CN}(0,1),
\]

最终

\[
\hat\Psi_{\rm ref}^{(m)}
=
\hat\Psi_{\rm coh}
+
\hat\Psi_{\rm sca}^{(m)}.
\]

随机种子使用 `sea_seed + surface_kstat_seed_offset + frequency_index - 1`，默认 `surface_kstat_seed_offset=200000`。若关闭随机散射，分支仍保留 \(P_{\rm sca}\)、\(S_{\delta G}\) 和能量诊断，但反射场只注入相干项。

K-Stat 记录全波数相位屏能量闭合：

\[
|\langle G\rangle|^2
+
\frac{1}{(2\pi)^2}
\sum S_{\delta G}(\mathbf K)\Delta K_x\Delta K_y
\approx 1.
\]

同时记录传播窗 \(K_h\le k_0\) 以及可选可信角窗内的非相干能量。传播窗能量只是诊断量，不会被强制重新归一化到 \(1-|R_{\rm coh}|^2\)。

## 6. PM 谱驱动的 SSA-like 统计散射分支

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
### 6.1 相干镜面项

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

这个指数项可理解为 SSA coherent reflection 的统计平均。它和 Kirchhoff realization 相位屏中的单次高度相位扰动不是同一个量：Kirchhoff 相位屏先生成具体 \(\eta(x,y)\)，再给每个空间点加相位；SSA coherent reflection 则直接对随机海面统计平均后得到 \(R_{\rm coh}\)。若海面高度近似为零均值高斯随机变量，反射相位扰动的方差越大，不同海面 realization 的镜面相干叠加越容易相互抵消。

当 \(H_s=0\) 时，\(\sigma_\eta=0\)，因此

\[
R_{\rm coh}=R_0.
\]

压力释放平整自由海面下即为 \(R_{\rm coh}=-1\)。

### 6.2 工程基线核：`pm_convolution`

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

### 6.3 一阶 Dirichlet 几何核：`ssa1_geometry`

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

## 7. 周期卷积与 zero-padding 线性卷积

当前实现支持两种卷积方式。

`periodic` 是默认兼容路径。它把有限波数网格看作周期延拓，使用 FFT 快速计算周期卷积。`ssa1_debug_dense` 在小网格上逐项显式求和，用来验证这个周期卷积公式和快速 FFT 结果一致。

`zero_padded` 是线性卷积审计路径。它先把谱搬到带符号波数顺序，在扩展网格上做零填充 FFT 线性卷积，再裁剪回原始波数范围。这样可以降低周期折返 aliasing 风险，更适合检查有限波数窗口边界对散射功率的影响。

两种方式的物理含义不同。周期卷积假设超出网格边界的波数会从另一侧折回；zero-padding 则把网格外未表示的谱当作零处理。因此二者不要求逐点相同，但都必须满足能量审计和接口不变量。

## 8. 能量约束与能量诊断

不同统计分支的能量处理不同。

`ssa_stat_kernel` 的非相干散射项是工程随机生成的。为避免随机散射导致非物理放大，该分支施加能量约束。

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

`kirchhoff_kstat` 不使用 `surface_ssa_scatter_scale` 和上述 SSA 能量限幅公式。它首先检查统计相位屏本身的能量闭合：

\[
|\langle G\rangle|^2
+
\frac{1}{(2\pi)^2}
\sum S_{\delta G}(\mathbf K)\Delta K_x\Delta K_y
\approx 1.
\]

然后单独记录传播窗和可信角窗内的非相干能量。这些窗口能量只是诊断，不会被重新归一化为总非相干能量。

## 9. 随机散射谱

若启用随机散射，`ssa_stat_kernel` 和 `kirchhoff_kstat` 都会把非相干功率谱变成一个具体复谱 realization：

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

宽带频率轴上，`ssa_stat_kernel` 还提供几个仅用于随机信道生成诊断的工程相关模式：默认各频点独立生成复高斯随机谱，也可以在不同频点复用同一个随机相位样本，或在频率序列上使用一阶自回归相关样本。它们的作用是检查宽带 \(H_f\) 的频域连续性对随机散射 realization 的敏感性；它们不是新的 SSA 物理散射公式，也不是经过实验标定的海面时间频率相关模型。`kirchhoff_kstat` 当前使用由 `sea_seed + surface_kstat_seed_offset + frequency_index - 1` 定义的可复现独立频点样本。

### 9.1 镜面项与非相干项的接收贡献诊断

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

这个诊断回答的是一个工程问题：在给定频率、风速、\(H_s\)、网格、随机种子和散射分支配置下，接收点处的非相干反射分量相对于镜面相干分量有多大。常用判据为

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

## 10. 通信链路解释

海面模型只改变频域信道 \(H(f)\)。后续 MPSK 通信仍沿用同一链路：

\[
H(f)\rightarrow H_{\rm baseband}(f)\rightarrow h_{\rm bb}(t),
\]

然后进行符号卷积、噪声注入、同步、均衡和判决。

因此比较 `kirchhoff_spatial`、`kirchhoff_kdomain`、`kirchhoff_kstat`、`pm_convolution` 和 `ssa1_geometry` 时，应主要观察

\[
|H(f)|,\qquad |H_{\rm ref}(f)|,\qquad \arg H(f),
\]

以及同一通信流程下的 BER/SER 趋势。

除通信链路外，当前还提供 LFM 测试信号验证脚本
`validate_kstat_vs_kdomain_lfm_channel_vertical.m`。该脚本不调用
`comm_main_vertical_psk.m`，也不使用 PSK、BER/SER、通信噪声或均衡。它先生成复解析基带 LFM 包络

\[
s_{\rm bb}(t)=w(t)\exp\{i\pi\mu(t-T/2)^2\},
\qquad
\mu=B/T,
\]

再把 `vertical_channel_model` 输出的 \(H_f\) 和 \(H_{\rm ref}(f)\) 插值到 LFM FFT 频率轴，通过频域乘法得到

\[
r_{\rm total}(t)=\mathcal F^{-1}\{S_{\rm LFM}(f)H(f)\},
\qquad
r_{\rm ref}(t)=\mathcal F^{-1}\{S_{\rm LFM}(f)H_{\rm ref}(f)\}.
\]

匹配滤波输出使用

\[
y_{\rm MF}(t)=\mathcal F^{-1}\{R(f)S_{\rm LFM}^*(f)\}.
\]

这个验证是 channel-level waveform validation：它用于直观看发射 LFM、接收总波形、反射波形和匹配滤波输出，并比较 ensemble 接收统计；它不是新的 PE/WAPE 空间源项，也不是通信系统 BER/SER 验证。

为便于直接观察发射/接收变化，报告脚本 `plot_lfm_tx_rx_comparison_vertical.m` 会读取已保存的 LFM 验证结果，额外生成 `lfm_tx_vs_rx_total_waveform_compare.png` 和 `lfm_tx_vs_rx_reflect_waveform_compare.png`。这两张图把发射 LFM、`kirchhoff_kdomain` 接收端和 `kirchhoff_kstat` 接收端画在同一坐标中；实部和包络都做了归一化，因此它们用于比较波形形状和时延/展宽趋势，不用于读取绝对接收幅度。

## 11. 当前验证结论

当前 reduced-grid 数值验证支持以下结论：

1. 默认路径保持兼容。`kirchhoff_spatial` 仍是默认边界模型；新增 `kirchhoff_kdomain`、`kirchhoff_kstat` 和 `raw_pm` 不改变默认 `target_hs` 配置。

2. 当 \(H_s=0\) 时，统计海面模型退化为平整自由海面反射：

   \[
   \sigma_\eta=0,\qquad W_\eta=0,\qquad R_{\rm coh}=R_0,\qquad E_{\rm sca}=0.
   \]

   `ssa_stat_kernel` 与 `kirchhoff_kstat` 均覆盖该退化检查。

3. SSA 物理趋势验证覆盖 \(H_s=0,0.05,0.2,0.5,1.0\) 和 \(f=4,6,8,10\ {\rm kHz}\)。在法向镜面基准下，`pm_convolution` 与 `ssa1_geometry` 都应满足

   \[
   R_{\rm coh}=R_0\exp(-2k_0^2\sigma_\eta^2).
   \]

   `ssa1_debug_dense` 的小网格显式求和与 `ssa1_geometry` 的周期 FFT 卷积结果一致；`zero_padded` 线性卷积路径可用于检查周期 aliasing 风险。

4. `kirchhoff_kstat` 验证覆盖相干项、非相干相位屏谱、随机 seed 可复现性和全 \(K\) 能量闭合。最近一次检查中，相位屏能量闭合误差最大约为 \(2.22\times10^{-15}\)，满足数值精度预期。

5. raw PM 风速驱动的 `kirchhoff_kdomain` / `kirchhoff_kstat` 对比使用
   `wind_list=[3 5 8 10 12 15]`、`seed_count=32`、`nx=ny=128`、\(f_0=6000\) Hz。按离散谱方差统一后，\(H_{s,\rm raw}\) 的两分支相对误差最大约为 \(1.45\times10^{-2}\)，说明海况强度已经基本对齐。

6. 归一化后，两个 Kirchhoff 分支不应要求逐点一致。`kirchhoff_kdomain` 是具体海面相位屏 realization，`kirchhoff_kstat` 是统计相位屏模型。合理比较对象是 ensemble 趋势、相干衰减、非相干能量、接收端 \(|h_{\rm ref}|\) / \(|h_{\rm total}|\) 统计和反射角谱展宽。

7. `validate_kstat_vs_kdomain_phase_screen_vertical.m` 是不接 PE/WAPE 的边界相位屏验证脚本。它用近垂直平面波入射，比较 `kirchhoff_kdomain` 多 realization ensemble 的相干均值、非相干功率谱、径向谱形状和能量闭合是否逼近 `kirchhoff_kstat` 统计相位屏公式。该脚本不要求单个 seed 的复数场逐点一致；若后续接入 PE，只应进一步比较接收端统计量，例如 \(E[|h_{\rm ref}|^2]\)、\(\mathrm{std}(|h_{\rm ref}|)\) 或平均 PDP。

8. `validate_kstat_vs_kdomain_channel_stats_vertical.m` 是接入完整 `vertical_channel_model` 后的接收端统计验证脚本。它比较 `kirchhoff_kdomain` 和 `kirchhoff_kstat` 在 \(E[|h_{\rm ref}|^2]\)、peak-synced reflected/total PDP ensemble mean、\(|h_{\rm total}|\) 分布、direct-path 一致性和 \(H_f=H_{\rm dir}+H_{\rm ref}\) 不变量上的统计一致性。该脚本不运行调制、噪声、均衡或 BER/SER；PDP 结论只表示当前宽带随机实现口径下的接收端统计检查，不表示已经建立了经过物理或实验标定的海面时间频率相关模型。最近一次 smoke run 使用 \(H_s=[0.1,0.2]\) m、8 个 seed、\(64\times64\) 网格和 8 个频点，硬检查通过：最大 \(H_f\) 不变量误差约 \(7.76\times10^{-18}\)，最大 K-Stat 能量闭合误差约 \(9.99\times10^{-16}\)，\(E[|h_{\rm ref}|^2]\) 最大相对误差约 0.097，reflected PDP 相关系数不低于 0.979，\(|h_{\rm total}|\) quantile NRMSE 最大约 0.193；其中 \(H_s=0.2\) m 的 reflected PDP L2 误差约 0.313，超过名义 0.25 统计目标，说明 full reduced 或更多 seed 仍是最终结论所需。

9. `validate_kstat_vs_kdomain_lfm_channel_vertical.m` 是独立的 LFM 测试信号验证脚本。它使用完整 `vertical_channel_model` 得到的 \(H_f\) 和 \(H_{\rm ref}(f)\)，把复解析 LFM 基带信号通过频域乘法送入信道，并比较接收总波形、反射波形、匹配滤波输出、PDP 和 \(E[|h_{\rm ref}|^2]\) 等统计量。最近一次 reflected-focused run 使用 \(H_s=0.2\) m、32 个 seed、\(64\times64\) 网格、32 个频点、4--8 kHz LFM、\(f_s=24\) kHz、\(T=20\) ms。硬检查通过：最大 \(H_f\) 不变量误差约 \(1.39\times10^{-17}\)，最大 K-Stat 能量闭合误差约 \(1.11\times10^{-16}\)，direct path 分支/seed 差异为 0。但 reflected-only 指标没有通过：\(E[|h_{\rm ref}|^2]\) 相对误差约 0.580，reflected LFM 包络相关系数约 0.0279，reflected PDP 相关系数约 0.276，reflected 匹配滤波峰值相对误差约 0.743。total-channel 指标明显更好，例如 total LFM 包络相关系数约 0.982、total PDP 相关系数约 0.930，因此只看 total LFM 会掩盖反射路径统计不一致。

10. 所有相关回归都必须保持

   \[
   H(f)=H_{\rm dir}(f)+H_{\rm ref}(f).
   \]

   最近 raw PM 对比和可视化诊断中的该代数不变量保持在约 \(10^{-17}\) 量级。

## 12. 周期卷积与 zero-padding 对比

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

## 13. 当前限制

当前模型仍有明确边界：

- `kirchhoff_kstat` 第一版是近垂直 Kirchhoff 统计相位屏模型，使用 \(\alpha=2k_0\)，尚未加入一般斜入射的 \(\gamma_i+\gamma_s\) 几何因子。
- `kirchhoff_kstat` 内部能量闭合的是相位屏统计量；传播窗和可信角窗能量只作为诊断记录，不会被强制归一化。
- `kirchhoff_kdomain` / `kirchhoff_spatial` 是具体 realization 模型，有限 seed 下的 \(H_s\)、\(|h_{\rm ref}|\) 和角谱指标会有样本波动。
- `pm_convolution` 是工程统计基线核，不是严格 SSA 几何核。
- `ssa1_geometry` 只覆盖压力释放 Dirichlet 自由海面的一阶几何因子。
- \(C_{\rm sca}\) 是工程归一化常数，不是实验标定的绝对散射截面。
- 周期卷积是默认兼容路径；zero-padding 线性卷积是 aliasing 审计路径。
- 任意阻抗边界、Neumann 边界和一般反射系数到 SSA 几何项的映射尚未实现。
- 当前只实现了可选的 Broschat-style SSA2 coherent reflection coefficient 诊断；二阶非相干散射功率、完整 SSA2 T-matrix、NLSSA、高阶多次散射和实验标定海面散射截面仍未实现。
- 当前不向通信反射场注入倏逝散射分量。
- 当前统计散射 realization 的跨频率相关性只有 independent、共享随机样本和一阶自回归这类工程随机模型；尚未建立物理或经验标定的宽带频率相关模型。因此宽带 \(H_f\) 的频域连续性仍需要后续专门研究。
- 在固定 \(H_s\) 归一化下，改变风速 \(U\) 主要改变 PM 谱形状；只有 `raw_pm` 模式才把风速同时解释为海况强度变化。

因此，当前 `kirchhoff_kstat` 应理解为 Kirchhoff statistical phase-screen model；`ssa_stat_kernel` 应理解为 PM-spectrum-driven first-order pressure-release Dirichlet SSA statistical reflection/scattering reference branch。其中 `pm_convolution` 是工程基线，`ssa1_geometry` 是一阶 Dirichlet 几何核；它们还不是完整的海面声散射理论闭环。

### 13.1 raw PM 网格与跨频 K-Stat 前置验证

2026-07-11 新增的独立验证接口不改变公共默认值。raw PM 审计将两个指标分开：二维离散谱方差与解析无限域方差之比，以及只考虑 `K_min--K_max` 的理想径向支持覆盖。前者还包含矩形网格求积误差，因此单独接近 1 不能证明谱峰已经解析。

固定 `dx=50/128 m` 时，U=5 推荐先在 100 m/256^2 PM 网格上建立海面统计，再以相同 `dx` 裁剪到 50 m/128^2 PE 窗口；32 seed 的平均映射能量误差为 2.04%。U=12 和 U=15 分别需要约 300 m/768^2 与 400 m/1024^2 才能覆盖低波数能量，因此不建议同步扩大完整 PE 网格。高风速更适合低/高 K 分解，但两个高度频带必须先合并到 `C_eta` 或 `eta`，不能在线性层面直接相加两个非线性相位屏散射谱。

跨频 joint kstat 使用

\[
C_{G,ij}(\rho)=R_{0,i}R_{0,j}^*
e^{-(\alpha_i^2+\alpha_j^2)\sigma_\eta^2/2}
\left[e^{\alpha_i\alpha_jC_\eta(\rho)}-1\right]
\]

和

\[
P_{G,ij}(\rho)=R_{0,i}R_{0,j}
e^{-(\alpha_i^2+\alpha_j^2)\sigma_\eta^2/2}
\left[e^{-\alpha_i\alpha_jC_\eta(\rho)}-1\right].
\]

伪协方差在 K 域耦合 `K` 与 `-K`，实现通过增广谱对协方差 EVD 抽样，而不是 shared seed、shared phase 或 AR(1)。U=5、4--8 kHz、F=32、64^2、64 realization 的边界验证中，joint 模式相对同一显式海面的跨频 covariance 误差为 0.67%，PDP 和 LFM matched-filter correlation 为 0.9837 和 0.9842；independent 模式分别为 95.26%、0.4890 和 0.5163。

该 U=5 条件的解析 `||P||_F/||C||_F` 只有 `2.09e-10`，所以此节点近似 proper；这不是删除 `P` 的通用依据。joint 原型尚未接入 PE，下一步必须在接收端重新验证 `C_H`、`P_H`、PDP 和 LFM。详细设置和成本见 `reports/raw_pm_joint_kstat_prerequisite_validation_report.md`。

### 13.2 cached joint-kstat 接收端原型

该独立验证路径不替换公共传播器。固定源、接收机、介质、PE 网格和频率轴后，它缓存海面入射场、直达接收响应、相干反射接收响应、surface-to-receiver 衍射因子、吸收屏以及同 `dx` 的 PM-to-PE 中央裁剪。每个 realization 只重复

\[
\Psi_{\rm ref,sca}^{(m)}\rightarrow H_{\rm ref,sca}^{(m)}(f),
\]

并保持

\[
H_{\rm total}=H_{\rm dir}+H_{\rm ref,coh}+H_{\rm ref,sca}.
\]

在 U=5 m/s、raw PM、4--8 kHz、F=32、PM 100 m/256²、PE 50 m/128²、训练 128 和测试 64 样本下，joint 模式相对显式同一海面 kdomain 的 reflected-only covariance 相对误差为 0.2704，PDP correlation 为 0.9714，LFM correlation 为 0.9917；independent 模式分别为 0.9750、-0.0445 和 0.5721。不能改用总信道重新计算这些指标，因为强直达项会掩盖反射误差。

F=32 时真实时延分辨率为 0.25 ms，最大无模糊时延为 7.75 ms。当前 99% 散射尾部接近该窗口，因此通信 CIR 验证前建议提高到 F=64；零填充只能插值，不能提高物理分辨率。

L=64 时，显式 kdomain 和 joint 的接收端样本 `||P||_F/||C||_F` 都约为 0.26。高维有限样本即使来自 proper 过程也会得到非零样本伪协方差，因此该结果尚不能证明过程 improper。完成 `scripts/validation/analyze_receiver_properness_null_vertical.m` 的 proper-null 校准前必须保留 P。完整设置、时间、内存和验收状态见 `reports/cached_joint_kstat_pe_receiver_validation_report.md`。

### 13.3 U=5、F=64 条件接收端统计生成器

固定条件统计模型将确定性分量与随机散射分量分开保存：

\[
H_{\rm total}^{(m)}=H_{\rm dir}+H_{\rm ref,coh}+H_{\rm ref,sca}^{(m)}.
\]

`estimate_conditional_channel_stats_vertical` 从 F×L 的散射接收样本估计均值、协方差、伪协方差、复EVD和实增广EVD。`sample_conditional_channel_vertical` 支持proper和improper路径；`build_physical_cir_vertical`检查等间隔频率轴并明确区分 `1/B` 真实分辨率、`1/df`无模糊时延和零填充插值间隔。

U=5、raw PM、4--8 kHz、F=64、训练128和测试64的验证中，joint接收端PDP/LFM correlation为0.9525/0.9839；统计生成器full-rank为0.9548/0.9858。F=64给出0.25 ms真实分辨率和15.75 ms无模糊窗；包含99%能量的最短圆周时延区间为4.657 ms。圆周区间用于避免把带限脉冲的峰前旁瓣错误解释为窗口末端长尾，不等同于零填充。

joint properness按central 95% null interval规则不能拒绝。显式kdomain也在central 95%区间内，但单侧p=0.0405，因此仍保留augmented-real路径。当前推荐full数值秩34；99.9%和99%候选秩为7和5，可作为压缩模式。完整数据、性能和限制见 `reports/u5_conditional_channel_generator_f64_report.md`。

## 14. 入射与海面反射声场的可视化

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

## 附：边界模型诊断脚本说明

部分脚本会绕开 PE/WAPE，只在海面边界处使用垂直平面波
\(\Psi_{\rm inc}(x,y)=1\)。这类脚本的目的不是生成通信信道，而是单独检查边界模型的相干反射、相位因子和 raw PM 方差口径。

这些诊断保留三个结论：

1. pressure-release 法向相干反射应满足
   \(R_{\rm coh}=-\exp(-2k_0^2\sigma_\eta^2)\)。`kirchhoff_kstat`、`ssa_stat_kernel` 的相干项，以及显式 Kirchhoff 多 realization 的 ensemble coherent average，都应在适用条件下支持这个趋势。

2. 当前 Kirchhoff 主相位约定是 \(\Delta\phi=2k_0\eta\)。旧的 \(\Delta\phi=4k_0\eta\) 只作为 `legacy_4k_eta` 历史对照，用来说明额外 2 倍因子会造成过强相干衰减。

3. raw PM 风速对比应使用正文第 4 节和第 5.2 节中的离散谱方差口径。附录不再单独维护另一套 raw PM 公式，避免和主实现说明分叉。

## 专题：精确离散伴随接收投影

验证专用伴随原型只适用于当前 cached uniform PE 路径：CPU double、固定发射机和频率/空间网格、单个最近网格点接收机、无气泡且无 Doppler。它不扩展公共传播模式，也不修改 `vertical_channel_model`、`vertical_wape_propagator`、默认 surface model 或 `comm_main_vertical_psk`。

对每个频率，令 `A_i` 表示 cached surface-to-receiver marching operator，`r` 表示 cached 接收网格点的单位脉冲。原型构造

\[
q_i=A_i^Hr,\qquad
a_{{\rm PE},i}=\operatorname{conj}(\psi_{{\rm inc},i})\odot q_i.
\]

这是已实现离散算子的精确共轭转置，不是互易反传、逆传播或损耗补偿。实现反向遍历深度步，共轭 Fresnel 与空间屏，保留 sponge 衰减，并且不引入经验 FFT 归一化。

joint-kstat 随机变量保持为

\[
\delta G_i=R_{0,i}\exp(i\alpha_i\eta)-R_{{\rm coh},i}.
\]

其当前协方差和伪协方差已经包含 `R0`，所以接收权重不能再次乘 `R0`。若 `E` 是现有 PM 到 PE 中央裁剪，则解析权重按 `a_PM=E^H*a_PE` 零嵌入，空间周期协方差始终在 PM 网格上计算：

\[
C_H(i,j)=a_{{\rm PM},i}^H C_{\delta G,ij}a_{{\rm PM},j},\qquad
P_H(i,j)=a_{{\rm PM},i}^H P_{\delta G,ij}a_{{\rm PM},j}^*.
\]

总均值单独按 `H_direct + H_ref_coh + mu_scatter` 构造；当前 joint-kstat 的散射均值为零。FFT 收缩使用未 shift 的 MATLAB 排列，lag 零点位于 `(1,1)`，并已由显式 PM 周期稠密矩阵验证。

完整验收的最大伴随/投影误差为 `7.77e-15`/`5.10e-15`，dense/FFT 的 `C/P` 误差为 `6.18e-16`/`9.35e-16`，同 realization 的 cached-forward/projection 误差在 F=9 和 F=64 下分别为 `2.52e-13` 和 `1.07e-13`。F=9 的 4096 条样本和 F=64 的 512 条样本均使解析/样本 `C/P` 误差落在 `1.25x` split-sample floor 内。F=64、PE 128²、PM 256² 时，主要解析成本来自 `F²` 频率对的 PM FFT/收缩。详细设置、reflected-only 通信派生指标、时间、内存、限制和接入判断见 `reports/adjoint_pe_receiver_projection_feasibility_report.md`。

## 15. 两节点条件统计信道的通信接入与验证

`comm_main_vertical_psk.m` 现在提供默认关闭的外部信道验证入口。仅当设置环境变量 `COMM_EXTERNAL_CHANNEL_FILE` 时，脚本才从 MAT 文件中的 `external_channels` 或 `external_channel` 读取 `H_f/f_axis_hz` 或 `h_t`；未设置该变量时仍执行原来的 PE 信道路径，公共默认配置和输出结构不变。外部信道不包含噪声，AWGN 仍由通信接收端在卷积之后单独注入。

`build_communication_taps_vertical.m` 统一完成物理宽带频响到符号率等效 taps 的转换。4--8 kHz 物理带宽对应真实时延分辨率 0.25 ms；插值到通信 FFT 网格只增加时延采样密度，不能提高这一物理分辨率。正式两节点验证采用最短圆周能量窗保留至少 99.9% 的 CIR 能量，再把该窗口展开为因果 taps，以避免把带限 IFFT 的峰前旁瓣误删为负时延。转换同时检查等间隔频率轴、频带覆盖和 H--IFFT--H 数值闭合。

2026-07-13 的正式验证覆盖 U=5 和 U=8 m/s、raw PM、4--8 kHz、F=64；比较 fresh `kirchhoff_kdomain + cached PE`、fresh `joint-kstat + cached PE`、条件统计 full-rank 和 99.9% 低秩模型。每条曲线使用 32 条独立信道、每信道 8000 个 QPSK 符号，并扫描 0:2:20 dB。统计 full-rank 相对 kdomain 的 BER/SER log10 RMSE 在 U=5 为 0.127/0.122，在 U=8 为 0.194/0.189；所有 Eb/N0 点的信道聚类 bootstrap 置信区间均重叠。四种来源的平均 BER 在 U=8 均高于 U=5，风速趋势一致。

共享潜变量的 128 对 full-rank/99.9% 低秩样本显示：U=8 的所有 Eb/N0 点差值置信区间包含零；U=5 仅 20 dB 出现约 `1.08e-4` 的小幅显著 BER 增量。因此 full-rank 仍是参考模型，99.9% 低秩可作为压缩候选，不能据此自动推广到更激进的 99% 截断。

当前已均衡 BER 在高 Eb/N0 出现约 `2e-3`--`7e-3` 的误码平台，无噪声运行也存在同量级残余误码。其主要来源是现有有限线性卷积块与循环 FFT 均衡器之间没有 CP、保护间隔或 overlap-save，以及固定正则化和有限 taps 截断。因此两节点库可用于同一接收机下的相对算法比较和趋势研究，但现阶段不应把该平台解释为信道物理模型的绝对高信噪比性能极限。完整设置、置信区间和结果路径见 `reports/two_node_statistical_channel_communication_validation_report.md`。
