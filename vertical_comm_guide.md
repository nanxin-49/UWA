# 垂直水声通信与海面统计散射模型说明

本文说明当前 MATLAB 垂直水声通信模型的物理含义、频域信道定义、海面反射/散射模型和验证结论。它面向阅读模型的研究或工程人员，而不是脚本接口说明。

## 1. 问题背景

模型描述一个海底发射端到近海面接收端的上行水声通信链路。发射信号从较深位置向上传播，接收端位于海面以下。总接收信道由两部分组成：

\[
H(f)=H_{\rm dir}(f)+H_{\rm ref}(f),
\]

其中 \(H_{\rm dir}\) 为直达传播贡献，\(H_{\rm ref}\) 为经海面反射或散射后返回接收深度的贡献。通信链路只使用最终频率响应 \(H(f)\)，因此海面模型的变化不要求重写调制、噪声注入、同步或均衡流程。

在参考频率 \(f_{\rm ref}\) 处，窄带等效复信道写为

\[
h_{\rm total}=H(f_{\rm ref}),\qquad
h_{\rm dir}=H_{\rm dir}(f_{\rm ref}),\qquad
h_{\rm ref}=H_{\rm ref}(f_{\rm ref}).
\]

## 2. 坐标与传播约定

海面定义为

\[
z=0,
\]

深度方向向下为正。若发射端深度为 \(z_{\rm tx}\)，接收端深度为 \(z_{\rm rx}\)，则上行通信满足

\[
0\le z_{\rm rx}<z_{\rm tx}.
\]

传播计算在横向平面 \((x,y)\) 上表示复声场，并沿 \(z\) 方向推进。直达项从 \(z_{\rm tx}\) 推进到 \(z_{\rm rx}\)。反射项先从 \(z_{\rm tx}\) 推进到海面，再由海面边界模型给出反射场，最后从海面推进回 \(z_{\rm rx}\)。

## 3. PM 海面高度谱

粗糙海面由 Pierson-Moskowitz（PM）谱描述。令

\[
K=\sqrt{K_x^2+K_y^2},
\]

一维波数谱写为

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

实际使用时，谱会按目标有效波高 \(H_s\) 重新归一化。海面高度标准差取

\[
\sigma_\eta=\frac{H_s}{4},
\]

并要求离散谱对应的连续方差满足

\[
\iint W_\eta(K_x,K_y)\,dK_xdK_y=\sigma_\eta^2.
\]

因此，\(W_\eta\) 表示目标海况下的海面高度波数谱，而不是某一次具体随机海面 realization。

当 \(H_s=0\) 时，

\[
\sigma_\eta=0,\qquad W_\eta=0,
\]

海面退化为平整自由海面。

## 4. Kirchhoff 空间海面模型

`kirchhoff_spatial` 是默认海面模型。它先从 PM 谱生成一个具体随机海面

\[
\xi(x,y),
\]

再用相位屏近似描述海面高度起伏导致的反射相位扰动：

\[
\Psi_{\rm ref}(x,y)
=
R_0\,\Psi_{\rm inc}(x,y)
\exp\left(i\Delta\phi(x,y)\right),
\]

其中 \(R_0\) 为平整海面反射系数。压力释放自由海面的默认近似为

\[
R_0=-1.
\]

相位扰动采用

\[
\Delta\phi(x,y)=2k_0\Gamma\,\xi(x,y),
\qquad
k_0=\frac{2\pi f}{c_0}.
\]

\(\Gamma\) 表示入射与反射方向共同决定的有效垂向相位因子；法向近似下可理解为固定的几何因子，斜入射修正时由局部传播方向决定。

该模型的特点是直观、可与具体海面 realization 对应，但每次随机海面都会改变相位屏。它不是严格的 SSA/NLSSA 散射截面模型，也不提供从海面谱直接生成统计散射功率的闭式通道。

## 5. SSA 统计散射分支

`ssa_stat_kernel` 分支不生成具体的 \(\xi(x,y)\)。它直接从 PM 高度谱 \(W_\eta(K_x,K_y)\) 构造统计散射功率，并在波数域合成海面反射场。

令入射场的横向波数谱为

\[
\Psi_{\rm inc}(K),\qquad
P_{\rm inc}(K)=|\Psi_{\rm inc}(K)|^2.
\]

反射场分为相干镜面项与非相干散射项：

\[
\Psi_{\rm ref}(K)
=
\Psi_{\rm coh}(K)+\Psi_{\rm sca}(K).
\]

### 5.1 相干镜面项

粗糙度会降低镜面相干反射强度。当前模型采用

\[
R_{\rm coh}
=
R_0
\exp\left[
-\frac12(2k_0\Gamma)^2\sigma_\eta^2
\right].
\]

因此

\[
\Psi_{\rm coh}(x,y)=R_{\rm coh}\Psi_{\rm inc}(x,y).
\]

当 \(H_s=0\) 时，\(\sigma_\eta=0\)，于是

\[
R_{\rm coh}=R_0.
\]

在压力释放自由海面默认条件下，

\[
R_{\rm coh}=-1.
\]

### 5.2 工程基线核：`pm_convolution`

`pm_convolution` 是当前统计散射分支的工程基线。它使用海面高度谱与入射功率谱的卷积来分配非相干散射功率：

\[
S_{\rm PM}(K,K')
\propto
C_{\rm sca}W_\eta(K-K').
\]

对应的散射功率可写为

\[
P_{\rm sca}^{\rm raw}(K)
=
C_{\rm sca}
\sum_{K'}
W_\eta(K-K')P_{\rm inc}(K')
\Delta K_x\Delta K_y.
\]

这里 \(C_{\rm sca}\) 是工程归一化常数，用于敏感性分析和数值调节；它不是实验标定的物理散射截面常数。

该核的优点是简单、稳定，并能验证“从 PM 谱直接生成统计散射通道样本”的代码路径。它的限制是没有显式入射/散射方向几何因子，因此不应解释为严格 SSA 一阶角度核。

### 5.3 一阶 Dirichlet 几何核：`ssa1_geometry`

`ssa1_geometry` 加入一阶 SSA 的压力释放 Dirichlet 几何因子。对横向波数 \(K\)，定义传播垂向波数

\[
\gamma(K,f)=\sqrt{\max(k_0^2-|K|^2,0)}.
\]

当前版本只把 \(|K|\le k_0\) 的传播分量纳入散射功率；倏逝分量不参与第一版统计散射能量。

压力释放 Dirichlet 条件下，一阶几何因子为

\[
G_{\rm SSA1}(K,K';f)
=
4\gamma(K,f)\gamma(K',f).
\]

因此统计核写为

\[
S_{\rm SSA1}(K,K';f)
\propto
G_{\rm SSA1}(K,K';f)W_\eta(K-K').
\]

相应的非相干散射功率为

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

等价地，可写成更紧凑的卷积形式：

\[
A(K')=\gamma(K',f)P_{\rm inc}(K'),
\]

\[
P_{\rm sca}^{\rm raw}(K)
=
4C_{\rm sca}\gamma(K,f)
\left[W_\eta * A\right](K)
\Delta K_x\Delta K_y.
\]

这里的卷积为当前网格上的周期卷积。`ssa1_debug_dense` 用小网格显式求和验证该周期卷积形式，与快速卷积结果一致。

`ssa1_geometry` 只对应压力释放 Dirichlet 边界，即 \(R_0=-1\)。任意阻抗边界、Neumann 边界或从一般反射系数到 SSA 几何因子的映射尚未实现。

## 6. 能量约束

非相干散射项是随机生成的。为避免随机散射导致非物理放大，模型对散射能量施加约束。

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
E_{\rm sca}^{\rm max}
=
\max(E_{\rm inc}-E_{\rm coh},0).
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

这是一条数值能量审计约束，不代表已经完成绝对散射截面的实验标定。

## 7. 随机散射谱与确定性

若启用随机散射，非相干散射谱按复高斯随机相位生成：

\[
\Psi_{\rm sca}(K)
=
\sqrt{P_{\rm sca}(K)}\,z(K),
\qquad
z(K)\sim\mathcal{CN}(0,1).
\]

最后得到海面反射场

\[
\Psi_{\rm ref}(K)
=
\Psi_{\rm coh}(K)+\Psi_{\rm sca}(K),
\]

再反变换回横向空间域，并沿原有反射路径传播到接收深度。

为保持可复现性，统计散射随机谱使用独立于 Kirchhoff 随机海面 realization 的种子流。频率索引改变时，随机散射样本随之确定性改变；相同海况、相同种子、相同频率轴会给出相同统计散射结果。

若关闭随机散射，模型仍计算 \(P_{\rm sca}\) 和能量审计量，但不把 \(\Psi_{\rm sca}\) 注入反射场。此时通信链路只接收相干镜面反射项。

## 8. 通信链路解释

海面模型只改变频域信道 \(H(f)\)。后续 MPSK 通信仍沿用同一信道消费方式：

\[
H(f)\rightarrow H_{\rm baseband}(f)\rightarrow h_{\rm bb}(t),
\]

然后进行符号卷积、噪声注入、同步、均衡和判决。也就是说，`kirchhoff_spatial`、`pm_convolution` 与 `ssa1_geometry` 的差异体现在海面反射贡献 \(H_{\rm ref}(f)\)，而不是通信接收机结构。

因此比较不同海面模型时，重点观察

\[
|H(f)|,\qquad |H_{\rm ref}(f)|,\qquad \arg H(f),
\]

以及在相同通信流程下得到的 BER/SER 趋势。

## 9. 验证结论

当前 reduced-grid 验证支持以下结论：

1. 当 \(H_s=0\) 时，统计海面模型退化为平整自由海面反射：

   \[
   \sigma_\eta=0,\qquad
   W_\eta=0,\qquad
   R_{\rm coh}=R_0,\qquad
   E_{\rm sca}=0.
   \]

2. `pm_convolution` 和 `ssa1_geometry` 在 \(H_s=0\) 时都与平整 `kirchhoff_spatial` 响应达到舍入误差级一致。

3. 对非零海况，随机散射能量满足

   \[
   E_{\rm coh}+E_{\rm sca}\le E_{\rm inc}.
   \]

4. 当关闭随机散射时，模型仍保留散射功率和能量统计，但不会向通信链路注入随机散射场。

5. `ssa1_debug_dense` 的小网格显式求和结果与 `ssa1_geometry` 的周期卷积结果一致，验证了一阶 Dirichlet 几何核的离散实现。

6. 增大 \(H_s\) 时，相干镜面项下降；散射功率预算增强。增大 \(C_{\rm sca}\) 时，原始散射功率单调增强，但最终散射能量仍受能量约束限制。

## 10. 当前限制

当前模型仍有明确边界：

- `pm_convolution` 是工程统计核，不是严格 SSA 几何散射核。
- `ssa1_geometry` 只覆盖压力释放 Dirichlet 一阶 SSA 几何因子。
- \(C_{\rm sca}\) 是工程归一化常数，不是实验标定的绝对散射截面。
- 当前卷积采用周期边界；非周期 zero-padding 卷积尚未实现。
- 任意阻抗边界、Neumann 边界和一般反射系数到 SSA 几何项的映射尚未实现。
- NLSSA、高阶多次散射和实验标定散射截面尚未实现。

因此，当前 `ssa_stat_kernel` 应理解为“PM 谱驱动的统计散射信道生成分支”，其中 `ssa1_geometry` 是一阶 Dirichlet SSA 几何核的实现；它还不是完整的海面声散射理论闭环。
