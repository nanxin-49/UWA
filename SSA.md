# `ssa1_geometry` 需要插入的 SSA 一阶几何核公式

## 1. 公式来源

`ssa1_geometry` 建议采用 Thorsos & Broschat, 1995, JASA 中针对 **Dirichlet / pressure-release rough surface** 的一阶 SSA 或一阶微扰极限公式。

当前项目中的自由海面可近似为压力释放边界，因此对应 Dirichlet 边界条件：

$$
p=0
$$

该边界下，一阶散射幅度与入射垂向波数和海面高度谱成正比。

---

## 2. 波数定义

声波参考波数为：

$$
k_0=\frac{2\pi f}{c_0}
$$

入射横向波数：

$$
\mathbf K'=(K_x',K_y')
$$

散射横向波数：

$$
\mathbf K=(K_x,K_y)
$$

波数差：

$$
\mathbf q=\mathbf K-\mathbf K'
$$

垂向波数定义为：

$$
\gamma(\mathbf K,f)
=
\sqrt{k_0^2-|\mathbf K|^2}
$$

其中：

$$
|\mathbf K|=\sqrt{K_x^2+K_y^2}
$$

入射垂向波数：

$$
\gamma_i=\gamma(\mathbf K',f)
$$

散射垂向波数：

$$
\gamma_s=\gamma(\mathbf K,f)
$$

第一版建议只保留传播波：

$$
|\mathbf K|<k_0,\qquad |\mathbf K'|<k_0
$$

否则令：

$$
\gamma=0
$$

---

## 3. 当前简化核

当前 `pm_convolution` 使用的是：

$$
S(\mathbf K,\mathbf K';f)
=
C_{\rm sca}
W_\eta(\mathbf K-\mathbf K')
$$

它只保留了 PM 海面谱控制的波数耦合，没有显式包含入射/散射方向的几何因子。

---

## 4. 需要插入的 SSA 一阶几何因子

Dirichlet 边界下一阶散射幅度可写成：

$$
T^{(1)}(\mathbf K,\mathbf K';f)
\propto
2i\gamma_i
\hat{\eta}(\mathbf K-\mathbf K')
$$

从幅度转成功率或能量通量时，引入垂向通量比：

$$
\frac{\gamma_s}{\gamma_i}
$$

因此一阶功率型几何因子为：

$$
G_{\rm SSA1}(\mathbf K,\mathbf K';f)
=
\frac{\gamma_s}{\gamma_i}
\left|2i\gamma_i\right|^2
$$

化简得到：

$$
G_{\rm SSA1}(\mathbf K,\mathbf K';f)
=
4\gamma_i\gamma_s
$$

也就是：

$$
G_{\rm SSA1}(\mathbf K,\mathbf K';f)
=
4
\sqrt{k_0^2-|\mathbf K'|^2}
\sqrt{k_0^2-|\mathbf K|^2}
$$

---

## 5. 新的 `ssa1_geometry` 散射核

将当前核替换为：

$$
S_{\rm SSA1}(\mathbf K,\mathbf K';f)
=
C_{\rm norm}
G_{\rm SSA1}(\mathbf K,\mathbf K';f)
W_\eta(\mathbf K-\mathbf K')
$$

即：

$$
S_{\rm SSA1}(\mathbf K,\mathbf K';f)
=
C_{\rm norm}
4\gamma_i\gamma_s
W_\eta(\mathbf K-\mathbf K')
$$

其中：

- $C_{\rm norm}$ 为数值归一化或能量归一化系数；
- $W_\eta(\mathbf K-\mathbf K')$ 为 PM 海面高度谱；
- $\gamma_i$ 为入射垂向波数；
- $\gamma_s$ 为散射垂向波数。

---

## 6. 对应散射功率

非相干散射功率为：

$$
P_{\rm sca}(\mathbf K,f)
=
\int
S_{\rm SSA1}(\mathbf K,\mathbf K';f)
\left|
\Psi_{\rm inc}(\mathbf K',f)
\right|^2
d\mathbf K'
$$

代入几何核：

$$
P_{\rm sca}(\mathbf K,f)
=
\int
C_{\rm norm}
4\gamma(\mathbf K',f)\gamma(\mathbf K,f)
W_\eta(\mathbf K-\mathbf K')
\left|
\Psi_{\rm inc}(\mathbf K',f)
\right|^2
d\mathbf K'
$$

离散形式：

$$
P_{\rm sca}[\mathbf K]
=
4C_{\rm norm}
\gamma[\mathbf K]
\sum_{\mathbf K'}
W_\eta[\mathbf K-\mathbf K']
\gamma[\mathbf K']
\left|
\Psi_{\rm inc}[\mathbf K']
\right|^2
\Delta k_x\Delta k_y
$$

---

## 7. FFT 加速形式

定义：

$$
A[\mathbf K']
=
\gamma[\mathbf K']
\left|
\Psi_{\rm inc}[\mathbf K']
\right|^2
$$

则：

$$
B[\mathbf K]
=
\left[
W_\eta * A
\right][\mathbf K]
$$

最后：

$$
P_{\rm sca}[\mathbf K]
=
4C_{\rm norm}
\gamma[\mathbf K]
B[\mathbf K]
\Delta k_x\Delta k_y
$$

FFT 实现对应：

$$
B
=
\mathcal F^{-1}
\left\{
\mathcal F[W_\eta]
\mathcal F[A]
\right\}
$$

---

## 8. 与现有 `pm_convolution` 的区别

原核：

$$
S_{\rm pm}(\mathbf K,\mathbf K')
\propto
W_\eta(\mathbf K-\mathbf K')
$$

新核：

$$
S_{\rm SSA1}(\mathbf K,\mathbf K')
\propto
4\gamma_i\gamma_s
W_\eta(\mathbf K-\mathbf K')
$$

也就是说，`ssa1_geometry` 只是在原有 PM 卷积基础上插入：

$$
4\gamma_i\gamma_s
$$

作为一阶 SSA / 一阶微扰极限下的几何因子。

---

## 9. 退化检查

当：

$$
H_s=0
$$

则：

$$
W_\eta=0
$$

所以：

$$
P_{\rm sca}=0
$$

此时反射场应只剩相干镜面反射：

$$
\Psi_{\rm ref}
=
R_0\Psi_{\rm inc}
$$

若：

$$
R_0=-1
$$

则：

$$
\Psi_{\rm ref}
=
-\Psi_{\rm inc}
$$

---

## 10. 需要在 metadata 中记录

建议记录：

```text
kernel_mode = 'ssa1_geometry'
formula_source = 'Thorsos & Broschat 1995 JASA, Dirichlet SSA first-order / perturbation-limit geometry'
G_SSA1_formula = 'G_SSA1(K,K'') = 4*gamma(K'')*gamma(K)'
boundary_condition = 'pressure-release / Dirichlet'
evanescent_included = false
```

## Addendum: Relation to `kirchhoff_kstat`

`surface_boundary_model='kirchhoff_kstat'` is not an SSA branch. It is a Kirchhoff / Gaussian statistical phase-screen branch used to generate coherent plus incoherent reflected spectra without a concrete `eta(x,y)` realization.

The kstat branch uses:

```text
alpha = 2*k0
<G> = exp(-0.5*alpha^2*sigma_eta^2)
R_coh = R0*<G>
C_deltaG(rho) = exp(-alpha^2*sigma_eta^2) * (exp(alpha^2*C_eta(rho)) - 1)
S_deltaG = F{C_deltaG}
P_sca(Ks) = |R0|^2/(2*pi)^2 int S_deltaG(Ks-Ki) |Psi_inc(Ki)|^2 dKi
```

Only `S_deltaG` is used for the incoherent spectrum. The coherent delta spike from the total phase-screen spectrum is not used as scatter power.

SSA1 remains useful as a weak-roughness and small-angle reference because the pressure-release coherent term reduces to the same normal-incidence form:

```text
R_coh = -exp(-2*k0^2*sigma_eta^2)
```

This agreement does not make the kstat branch an SSA/NLSSA/T-matrix implementation, and it does not provide a calibrated scattering cross section.

`surface_roughness_scale_mode='raw_pm'` does not change the SSA formulas above. It only disables the usual `sea_hs_target` rescaling so wind speed directly sets the PM-spectrum roughness used by Kirchhoff branches. Raw-PM Kirchhoff comparisons use the project discrete PM variance as the main convention; the K-Stat branch maps that convention into its continuous correlation formula with `W_eta=Phi2D*(2*pi)^2`.
