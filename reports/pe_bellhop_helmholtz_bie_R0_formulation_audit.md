# PE--Bellhop--Helmholtz BIE R0 formulation audit

日期：2026-09-14  
状态：**R0 PASS / FORMULATION FROZEN WITH EXPLICIT LIMITS**

## 1. Decision

第一轮采用二维 sound-soft rough-surface 的 **Dirichlet half-plane Green
function combined-layer equation**，配合 slow-rise smooth finite section 和高阶
Nyström 离散。它直接求解完整 Helmholtz 边值问题，不使用 PE phase screen、
Bellhop `Reflect2D`、ray 或 local-specular approximation。

第一轮不采用：

- closed-obstacle CFIE：其几何与无限压力释放海面不同；
- locally-rough bounded-contour formulation 作为主解：该方法把 flat reflected
  field 作为已知背景，flat case 的 defect density 恒为零，不能充分执行 R1
  的非平凡 BIE self-validation；
- quasi-periodic Green function：当前 Gaussian angular spectrum 不是单一 Bloch
  入射，第一轮会额外引入 Bloch 分解和 Wood-anomaly 处理；
- free-space kernel 的 abrupt finite cut：不接受为 reference。

主要理论依据：

- Zhang 与 Chandler-Wilde 的无限 rough sound-soft surface 二类积分方程，使用
  combined double/single layer 和 Dirichlet half-plane Green function，并证明对
  所有 wavenumber 唯一可解：
  <https://doi.org/10.1002/mma.361>
- Meier 与 Chandler-Wilde 对 rough-surface finite-section 稳定性和收敛性的分析；
  对一般较大 surface 采用 endpoint flattening：
  <https://doi.org/10.1002/mma.210>
- corrected/slow-rise WGF 说明普通截断与 slow-rise operator windowing 的区别，
  并给出随 window 增大的高阶收敛结构：
  <https://arxiv.org/abs/1507.04445>
- locally rough bounded-contour 方法保留为独立交叉核查候选，而非 R1--R5 主解：
  <https://arxiv.org/abs/1302.0161>

## 2. Frozen physical problem

采用项目物理坐标 `(x,z)`，其中 `z` 向下为正：

```text
water domain D = { (x,z) : z > eta(x) }
surface Gamma  = { (x,eta(x)) }
f              = 4000 Hz
c              = 1500 m/s
k              = 2*pi*f/c
z_tx           = 100 m
z_rx           = 3 m
sigma          = 0.3 m
time convention = exp(-i*omega*t)
```

总场满足

```text
(Delta + k^2) u = 0             in D
u = u_inc + u_ref = 0           on Gamma
```

`u_ref` 满足向 `z -> +infinity` 传播/辐射的 rough-surface radiation
condition。二维 Green function 采用

```text
Phi_k(X,Y) = (i/4) H_0^(1)(k |X-Y|),
```

与 `exp(-i*omega*t)` 的 outgoing convention 一致。

## 3. Incident Gaussian angular spectrum

BIE 入射场不是新的 point source，而是和 1-transverse PE 相同的 upward
Gaussian angular spectrum：

```text
u_inc(x,z) = (1/2*pi) integral A(kx)
             * exp(i*kx*x - i*kz*(z-z_tx)) dkx

A(kx) = sqrt(2*pi)*sigma*exp(-sigma^2*kx^2/2)
kz    = sqrt(k^2-kx^2), Re(kz)>=0, Im(kz)>=0.
```

传播分量和 evanescent 分量使用同一固定 square-root branch；后者从 Tx 平面
向海面衰减。实际离散必须直接复用冻结的 AS sample/grid definition，而不是根据
BIE/Bellhop 误差重新拟合 spectrum。

flat pressure-release 解析反射场为

```text
u_ref,flat(x,z_rx) = -(1/2*pi) integral A(kx)
                     * exp(i*kx*x + i*kz*(z_tx+z_rx)) dkx.
```

它是 R1 的独立 method-of-images/angular-spectrum reference。

## 4. Half-plane Green function and orientation

选择固定辅助直线 `z=h`，满足

```text
h < min_x eta(x),
```

因此该直线严格位于 water domain 外侧。为避免不同 surface 使用不同积分表示，
R1--R5 统一冻结

```text
h = -2.0 m,
```

它低于（物理上高于）本轮 `|eta| <= 0.20 m` 的全部表面；若以后扩大表面高度，
必须先重新检查 `h < min eta`，不能逐 case 静默改变 `h`。

对 `Y=(y,z_y)` 定义其关于 `z=h` 的镜像

```text
Y_h = (y, 2*h-z_y)
```

和 Dirichlet half-plane Green function

```text
G_h(X,Y) = Phi_k(X,Y) - Phi_k(X,Y_h).
```

`G_h=0` on `z=h`，image cancellation 也改善远距离 kernel decay。

法向固定为**指向水体 D 内部**：

```text
nu(y) = (-eta'(y), 1) / sqrt(1+eta'(y)^2).
```

不得在实现中改成 outward normal 后仍沿用下面的 jump sign。

## 5. Frozen combined-layer equation

设置实 coupling parameter

```text
beta = k > 0.
```

反射场表示为

```text
u_ref(X) = integral_Gamma [dG_h(X,Y)/dnu_Y - i*beta*G_h(X,Y)]
                            psi(Y) ds_Y.
```

由于 `nu` 指向 D，且从 D 一侧取 trace，double-layer jump 冻结为

```text
gamma_D D_h = +1/2 I + K_h.
```

因此边界方程为

```text
(1/2 I + K_h - i*beta*S_h) psi = -u_inc|Gamma,
```

等价地

```text
(I + 2*K_h - 2*i*beta*S_h) psi = -2*u_inc|Gamma.       (R0.1)
```

The minus sign is the outgoing combined-field sign for the frozen
`exp(-i*omega*t)` / `H_0^(1)` convention.  An initial R1 smoke run exposed
that the provisional plus sign produced a refinement-dependent near-null
system (at 8 points per wavelength `RCOND` was approximately `2.47e-18`).
Re-deriving the sign from the exterior trace and using `D_h-i*beta*S_h`
removed that pathology; this correction was made before accepting any R1
result and is not a fitted phase or amplitude adjustment.

R1 必须通过 flat image solution 同时验证 `Phi`、normal、jump 和 `-i*beta`
四个符号。禁止在 R1 中尝试多个符号并选择误差最小者；上述定义是唯一候选。

## 6. Smooth finite section

当前 Gaussian incident boundary data 随 `|x|` 衰减，不是非衰减 plane wave。
因此第一版使用 half-plane kernel 上的 slow-rise finite section：

```text
w_L(x) = 1                       |x| <= cL
       = C-infinity slow-rise    cL < |x| < L
       = 0                       |x| >= L
```

冻结 `c=0.7`。离散方程只在 `Gamma_L` 上求解，但所有积分 density 均乘
`w_L`；不允许 rectangular/abrupt cut。reference 值必须由至少三个 `L`
的收敛结果定义，不允许把单个有限窗口称为 Helmholtz reference。

冻结的有限窗口方程及场重建为

```text
[I + 2*K_h*W_L - 2*i*beta*S_h*W_L] psi_L = -2*u_inc  on Gamma_L,

u_ref,L = (D_h - i*beta*S_h)[w_L*psi_L].
```

其中 `W_L psi=w_L psi`，identity jump term 不乘 window。由于本 Goal 的
Gaussian boundary forcing 本身横向衰减，这一 smooth finite-section 与非衰减
plane-wave 问题不同；是否足够准确仍完全由 R1/R2 的三档 `L` 收敛决定。

物理表面与数值窗口严格分开：

```text
eta_bench = chi_phys * eta_raw
supp(chi_phys) lies completely inside {|x| < 0.7 L_min}
```

`chi_phys` 是 PE/Bellhop/BIE 共同物理输入；`w_L` 只属于 BIE 数值求解。
R2 分别记录 physical-support sensitivity 和 `L` convergence。

如果 R1/R2 显示 slow-rise finite section 未达到 gate，则停止；下一候选是把
同一 equation 升级为 analytic-flat-tail corrected WGF，而不是扩大 rough case、
改相位或调 normalization。

## 7. Nyström discretization

第一版采用 smooth-panel Nyström：

1. 以 `gamma(s)=(s,eta_bench(s))` 参数化真实墙面；Jacobian 为
   `sqrt(1+eta'(s)^2)`。
2. 面板使用 Gauss--Legendre nodes；主档每 panel 16 nodes。
3. `S_h` 的 logarithmic self singularity 使用 Kress/Alpert product quadrature；
   不以 `realmin`、删 diagonal 或 receiver offset 规避。
4. `K_h` 使用解析 diagonal limit；near panels 使用加密或专用 close
   quadrature。
5. receiver evaluation 使用同一 layer potential，但 receiver 不参与
   collocation；对 near target 单独做 close-evaluation convergence。
6. 三档 spatial density 必须同时满足 acoustic 与 geometry resolution：

```text
h_node <= min(lambda/N_lambda,
              (2*pi/Kmax)/16,
              0.1/max(abs(kappa)))
N_lambda = 8, 12, 16.
```

7. 初期使用 dense solve，并记录 condition estimate 与 residual；只有 unknown
   count 或内存要求触发后才引入 GMRES/FMM，R1/R2 不提前改变 solver family。

## 8. Receiver evaluation and comparison convention

BIE 只输出 reflected field `H_BIE_ref(x,z_rx)`。定义

```text
G_BIE = H_BIE_ref,rough / H_BIE_ref,flat.
```

三方比较沿用当前冻结 receiver line、M95/M99 mask 和指标：

```text
E_G, phase RMS, TL RMS, rho_shape, phi0, E_aligned.
```

Bellhop 继续使用 Stage-1X/1Y 已冻结的单次 conjugation mapping。BIE 的
`exp(-i*omega*t)` 约定由 R1 image solution 独立固定；禁止 per-case
conjugation selection、complex scalar fitting、phase subtraction 或 amplitude
renormalization。

## 9. Mandatory R1/R2 evidence

R1 flat 必须至少保存：

- collocation residual；
- 位于不同 nodes 的 off-grid boundary residual；
- receiver-line complex L2、phase RMS、TL RMS、`rho_shape`；
- three-level spatial/window convergence；
- normal/jump/time-convention fingerprint；
- NaN/Inf count。

Gate 保持 Goal 原值：

```text
off-grid boundary residual <= 1e-8
receiver complex L2        <= 1e-6
receiver phase RMS         <= 1e-6 rad
NaN/Inf count              = 0
```

R2 weak rough 还必须加入 3 spatial、2 quadrature、3 window levels，并把
最终两档变化的 envelope 定义为 `U_BIE`。不得用 PE/Bellhop difference 定义
`U_BIE`。

## 10. Limits and stop conditions

- R0 formulation 只覆盖 2D、uniform c、pressure-release、Gaussian localized
  incidence 和 graph surface `z=eta(x)`。
- 它不覆盖 3D、layered SSP、impedance surface、overhang 或 point-source
  singular incident field。
- slow-rise finite section 的可接受性必须由 R1/R2 convergence 证明；理论审查
  不能替代数值 gate。
- BIE2D 可作为 kernel/quadrature 对照，但不是可直接调用的完整 rough-surface
  solver；不得把其 closed-curve example 直接当结果。
- fixed PM、6/8 kHz、multi-seed、FMM 均不属于 R1--R5 第一轮。

## 11. R0 gate conclusion

所选 equation、Green function、normal orientation、jump relation、time
convention、radiation direction、window treatment、Nyström discretization 和
receiver evaluation 已明确冻结，且没有把 closed obstacle 或 abrupt truncation
误作 rough-surface reference。

```text
R0 = PASS
NEXT = implement R1 flat self-validation only
```

R1 若不能在 spatial/window refinement 下达到解析 image gate，立即停止，不进入
rough surface。
