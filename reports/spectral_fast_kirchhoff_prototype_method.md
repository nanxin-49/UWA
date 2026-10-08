# Spectral Fast Kirchhoff prototype：离散结构与近似边界

日期：2026-10-07。本文是 validation-only 算法说明；数值结果和成功/失败判断以
`spectral_fast_kirchhoff_prototype_report.md` 及本轮 CSV/MAT 为准。
物理符号、法向和时间约定直接继承 `solve_full_kirchhoff_2d.m`。

## 1. 完整离散 reference

均匀曲面参数网格 `x_j=x_0+j dx`，真实高度 `eta_j`，接收深度 `z_r>max eta`。
`n_j=(n_xj,n_zj)` 指向水域，`q_j` 为原 solver 的弧长梯形权重。
定义 lag `t=x_i-x_j`、`h=z_r-eta_j`、`R=sqrt(t²+h²)`：

\[
G(t,h)=\frac{i}{4}H_0^{(1)}(kR),\quad
D_z(t,h)=\frac{ik}{4}H_1^{(1)}(kR)\frac{h}{R},\quad
D_x(t,h)=\frac{ik}{4}H_1^{(1)}(kR)\frac{t}{R}.
\]

`D_x,D_z` 是**对源点坐标**的偏导。特别地 `D_z=-partial_h G`，不能把它当作
接收坐标导数并再换一次符号。于是原 reference 的严格离散式为

\[
u_i=-\sum_jq_j\{G(t_{ij},z_r-\eta_j)d_j
 +[D_z(t_{ij},z_r-\eta_j)n_{zj}+D_x(t_{ij},z_r-\eta_j)n_{xj}]u_j\},
\quad d_j=n_{xj}u_{xj}+n_{zj}u_{zj}.
\]

这保持完整有限支撑曲面，而不是接收中心小孔径。
Prototype-A 把 `n=(0,1)`、`q` 改为 `dx` 梯形权重；这是**物理/几何简化**。
Prototype-B 直接恢复上述完整 `n,q,d,D_nG`，没有新的物理模型。

## 2. 入射侧：同一 weighted SVD，有限次 IFFT

继承保存的离散 Gaussian 源 `(co_a,kxa,kza,x_src0,n_original,zTx)`，
`kza=sqrt(k²-kxa²)`，不截到 `theta99_inc`。本数据的全部保存箱都是传播波。

\[
u_{inc}(x,z)=\frac1{n_0}\sum_a co_a\,
 e^{ik_{xa}(x-x_{src0})}e^{ik_{za}z_{Tx}}e^{-ik_{za}z}.
\]

先作**严格因式分解**：

\[
e^{-ik_{za}\eta}=e^{-ik\eta}C(\eta,k_{xa}),\qquad
C=e^{i(k-k_{za})\eta}.
\]

在上一阶段同样的 257 点均匀 `[-A,A]` 高度表上，
`B=C diag(sqrt(w))`、`w=|co|²/sum|co|²`；取 weighted SVD 右子空间 `V_r`。
只重建这些已测案例需要的小表，没有新的 A–K/rank 扫描，也没有 per-receiver 拟合。

\[
C_r(\eta,a)=\sum_{m=1}^{r_{in}}b_m(\eta)a_m(a),\quad
b(\eta)=C(\eta,:)\operatorname{diag}(\sqrt w)V_r,
\quad a_m(a)=\overline{V_{am}}/\sqrt{w_a}.
\]

这是**全局低秩近似**。零能量源箱可以严格移除，非零尾部全部保留。
`b(eta)` 在固定 33 点 Chebyshev-Lobatto 表上计算，通过 barycentric 插值求实际高度值；
这是单独的**数值插值近似**，没有逐位置 SVD。强载波 `exp(-ik eta)` 始终精确求值。

\[
F_m(x)=\frac1{n_0}\sum_a co_a e^{ik_{za}z_{Tx}}a_m(a)
 e^{ik_{xa}(x-x_{src0})},\quad
u_j=e^{-ik\eta_j}\sum_m b_m(\eta_j)F_m(x_j).
\]

每个 `F_m` 是一次 IFFT。计算 `ux,uz` 时分别在相同谱系数上乘 `i kx` 和 `-i kz`，
再使用同一 `b(eta)`；这得到源场的物理坐标偏导，而不是对限制在曲面的迹作切向求导。
沿 `eta(x)` 求导会混入 `eta' uz`，不能替代 `partial_x u_inc`。

源周期 `W=192.1875 m` 与计算 `dx` 必须对齐，`M=W/dx` 为整数。
将整数 harmonic 装入 `mod(harmonic,M)` 箱，MATLAB IFFT 后乘 `M/n_original`，
这是**严格离散 Fourier 级数采样恒等式**。粗 scaling 网格上将多个原 harmonic
累加到同一箱是原级数在粗网格的精确取样；不代表物理频谱已充分解析。
不会通过改变原源 Nyquist 模式来偷偷降低基线成本。

## 3. 出射侧：高度分离后的有限 Toeplitz 核

无限面角谱形式可写为

\[
G(t,h)=\frac1{2\pi}\int_{\mathbb R}\frac{i}{2k_z(\kappa)}
 e^{i\kappa t}e^{ik_z(\kappa)h}\,d\kappa,
\qquad \operatorname{Im} k_z\ge0.
\]

对源点的 `z` 导数对应 `-i kz`；水平导数对应 `-i kappa`。
传播区 `|kappa|<k` 与倏逝区 `|kappa|>k` 均属于该表示。直接对 `1/kz` 在均匀
角谱网格上积分会需要特别处理 `kz=0` 分支、无限 lag 支撑与周期 Green images。
本版以**与 reference 相同的 Hankel 采样**构造其有限 Toeplitz Fourier symbol，
不以未经验证的均匀角谱 quadrature 替换 reference，也不加入分支 epsilon 或经验尺度。

在所有可能用到的 lag `ell=-(N-1),...,N-1` 和固定 33 个高度节点 `eta_q` 上，定义

\[
\widetilde H_\alpha(\ell,\eta)=e^{ik\eta}H_\alpha(\ell dx,z_r-\eta),
\qquad H_\alpha\in\{G,D_z,D_x\}.
\]

依然完整保留强高度相位：`H_alpha=e^{-ik eta} Htilde_alpha`。
将 `[Gtilde; Dztilde/k; Dxtilde/k]` 叠成 **lag × 固定高度**的小列数矩阵，
用同一个右高度子空间作全局 SVD，统一取 `r_out`：

\[
H_{\alpha,r}(\ell,\eta)
=e^{-ik\eta}\sum_{p=1}^{r_{out}}h_{\alpha p}(\ell)b^{out}_p(\eta).
\]

`D/k` 仅为导数量纲匹配，重构时严格乘回 `k`，不是按接收结果拟合的尺度因子。
该步是**出射全局低秩近似**；高度 basis 同样采用 Chebyshev 插值。
它不建立 `receiver × surface` dense kernel，也不建立 `surface × incident spectrum`
矩阵。33 个高度节点只受高度/核平滑性控制，不随 `N` 增长。

首批强高度控制中，`r_out=4` 的误差不够；统一 `r_out=6` 是当前候选。
因此必须单独报告 `r_in/r_out`，不能把上一阶段的入射 rank 2 宣传为整个 FK 的 rank。
完整 outgoing rank=1–6 表保留，未隐藏失败行。

## 4. 两个 Green 项和真实法向的 FFT 组织

令

\[
v_{Gp,j}=q_j e^{-ik\eta_j}b^{out}_p(\eta_j)d_j,\quad
v_{zp,j}=q_j e^{-ik\eta_j}b^{out}_p(\eta_j)n_{zj}u_j,\quad
v_{xp,j}=q_j e^{-ik\eta_j}b^{out}_p(\eta_j)n_{xj}u_j.
\]

那么在选定低秩核下，下式是**严格离散重排**：

\[
u_r=-\sum_p[h_{Gp}*v_{Gp}+h_{zp}*v_{zp}+h_{xp}*v_{xp}].
\]

每个 `*` 是有限线性卷积，通过 FFT(kernel) × FFT(density) → IFFT。
Prototype-A 只计算前两项，Prototype-B 还计算完整水平导数项。
不同高度只出现在 density 乘子 `eta_j,b_out(eta_j)` 中；kernel 只依赖 lag，
所以全局非局部相干作用真正进入 FFT，没有 dense contraction 回退。

## 5. 归一化、支撑、padding、对齐与 evanescent

- MATLAB forward FFT 无 `1/L`、inverse FFT 含 `1/L`。积分权重 `q_j` 已含 `dx/ds`，
  不能再乘经验 `dx`、`L` 或 `2pi` 去校准 FullFK。
- surface physical support 固定 `[-60,60] m`，正弦高度 C2 plateau/support 为 42/50 m。
  外侧平面区域仍属于积分，不缩为零高度区域内的小孔径。
- surface 和 receiver 均为同一均匀网格上的 `N` 个点，receiver 单深度。当前接口拒绝
  非均匀格点或未对齐的源周期；不把任意 receiver 插值误差藏入算法 fidelity。
- kernel lag 从 `-(N-1)` 到 `N-1`，按 `lag mod L` 装入 FFT 长度 `L=pN`。
  `L>=2N-1` 时输出前 N 点严格等于完整 Toeplitz 线性卷积。无需 `3N` 长度，
  因为只求同一 N 点输出，而不是保存全部线性卷积的外侧输出。
- padding=1 将正负 lag 重叠累加，是**有意失效的周期折叠控制**；padding=2/4 不改变
  物理 support、dx、height、source，只改变算术嵌入长度，二者应在 roundoff 内一致。
- 无 sponge，无 padding-dependent taper。receiver 谱图的 Hann 窗只用于谱诊断，
  不影响求解和误差指标。
- 入射保存谱本身均为传播波；所有非零箱保留。出射使用 Hankel 核，没有
  `theta99_inc`、`|kx|<k` 或 evanescent cutoff；采样的 near-field/evanescent 效应
  保留在 kernel Fourier symbol 中。有限网格/有限 lag 的 symbol 不是无限面远场谱。

## 6. 复杂度、计时与不可外推部分

入射：`O(r_in M log M + q r_in N)`；出射 apply：`O(r_out L log L)`，
实际常数包括 2 或 3 个导数通道。setup 中固定 `q=33` 的 Hankel 表及 SVD
为 `O(q N + q²N)`；原源 weighted basis 是固定 257 高度 × 保存源谱的小表。
全局 basis 是 physics-input kernel 压缩，不以 FullFK 结果或每个 receiver 做拟合。

`r_in+r_out` 与导数通道的 FFT 次数都显式披露。没有消失的 dense 工作转移到缓存：
`plan` 含 `N×r`、`(2N-1)×r` 数组，setup 没有 `N×N` 数组。
reference 验证路径仍按原 solver 计算 dense Green blocks；它从 Fast apply 路径隔离。

计时同时报告 Fast setup+apply/cold、warm apply、Full exact incident、Full dense apply、
CPU time 和 wall time；warm 速度不能替代冷启动速度。固定 `N_surface=N_receiver`
的 513/1025/2049/4097 只是四点经验 scaling，不声称精确渐近指数或跨硬件实时速度。
数组 payload 与原 solver block 存储估计不是 peak RSS。各个 mesh 的 fidelity 是对
相同离散 reference，不等于最粗网格已解析连续 Helmholtz 解。

后续若要 PM pilot，须先明确统一高度界、空间导数、出射带宽与 near-field 解析范围。
本任务只给出确定性单频 prototype 的成功/停止建议，不自动运行 PM、Monte Carlo、
SSP、3D、BER 或生产集成。
