# Bellhop 2020 tilted straight internal-wall validation-only POC

## 结论

**几何与 beam-state POC：PASS；native Cartesian 绝对幅度对照：FEASIBLE_WITH_LIMITS。**

本阶段只实现并验证了固定直线

\[
r_w(z)=R_0+a z,
\qquad R_0=100\ \mathrm{m},\quad a=0.005.
\]

9 组 step/beam-count 运行均满足：解析交点、局部标架、镜面反射、一次
pressure-release \(\pi\) 相变、`Amp`、`p/q`、proper rotation 和正 range
增量。没有修改 PE、通信代码、`InfluenceGeoHatCart` 或正式 Bellhop
executable，也没有加入 Kirchhoff、SSA 或其他散射公式。

同一条物理直线的 native Bellhop `.ati` case 在代表性 arrival 的反射点、
方向和时延上与 rotated-wall case 一致。native 原始 range chart 对一个向
负 range 返回的 reflected branch 仍产生稳定约 `0.286 dB` 的 Cartesian
幅度偏移；phase 差约 `0.003 rad`。该偏移被记录为 native backward-range
influence 的诊断，不做拟合、重置或校正，因此不改变物理反射实现。

## 实现

新增文件：

- `scripts/validation/support/bellhop_internal_tilted_wall_poc/Step.f90`
- `scripts/validation/support/bellhop_internal_tilted_wall_poc/bellhop.f90`
- `scripts/validation/support/bellhop_internal_tilted_wall_poc/build_bellhop_internal_tilted_wall_poc.ps1`
- `scripts/validation/support/run_bellhop_internal_tilted_wall_poc_vertical.m`
- `scripts/validation/support/run_bellhop_native_tilted_wall_vertical.m`
- `scripts/validation/validate_bellhop_internal_tilted_wall_poc.m`

构建输出为独立 binary：

`results/validation/bellhop_internal_tilted_wall_poc/bin/bellhop_iwall_tilted_2020.exe`

SHA-256：

`6186A7D2E046E65B8D095F3B6EC965AE2D3AEE043CCDBD511B8B9C627AA0D5E1`

`.iw2` sidecar 的三行依次为 `R0=100`、`a=0.005` 和 mapped receiver
range `103 m`。binary 在运行时拒绝其他 slope、非均匀/有损声速、非零源
横向深度、非 coherent Cartesian geometric-hat 配置和缺失 mapped receiver。

### 解析直线求交

在 `ReduceStep2D` 的 Euler 和 midpoint 两次 step reduction 中使用

\[
F(r,z)=r-R_0-a z,
\qquad
h_w=-\frac{F(x_0)}{u_r-a u_z},
\]

仅在射线从水体侧 \(F<0\) 向外穿越且分母为正时把 \(h_w\) 加入现有
top、bottom、SSP 和 range-segment 限制的最小值。这样 accepted node 就是
真实 wall intersection；不会先穿墙再补反射。与 flat POC 相同，只有真正
被选为最短 wall event 的 step 才绕过通用 infinitesimal-step guard。

### 原生反射与局部标架

直线 wall 使用常量标架

\[
\mathbf t_w=\frac{(a,1)}{\sqrt{1+a^2}},\qquad
\mathbf n_w=\frac{(1,-a)}{\sqrt{1+a^2}},\qquad \kappa_w=0.
\]

其中 \(\mathbf n_w=-J\mathbf t_w\)，sense 与 Bellhop `TOP` 一致。命中后
唯一的物理反射调用仍是：

```fortran
CALL Reflect2D( is, WallHS, 'TOP', WallT, WallN, WallKappa, RTop, NTopPTS )
```

因此镜面方向、vacuum `Phase += pi`、`Amp` 以及原生 curvature `p/q`
更新都来自 Bellhop 2020 `Reflect2D`；直线的 `kappa=0` 不产生额外 kick。

### 反射后 proper rotation

反射节点只执行固定的 determinant-`+1` 变换

\[
\begin{bmatrix}r'\\z'\end{bmatrix}
=
\begin{bmatrix}-1&0\\0&-1\end{bmatrix}
\begin{bmatrix}r\\z\end{bmatrix}
+\begin{bmatrix}2R_0\\0\end{bmatrix},
\qquad \mathbf t'=-\mathbf t.
\]

这是坐标重映射，不再次调用 `Reflect2D`，不改变 `p/q`、`tau`、`Amp` 或
`Phase`。transformed branch 的 chart seam 在 trace 完成后被隔离，只有
post-wall 节点数组传给未修改的 `InfluenceGeoHatCart`。

## 验证设置

- Bellhop 2020 `2020_11_4`，2D；
- 均匀、无损 `c=1500 m/s`；
- source/receiver transverse depth `0 m`；
- 4 kHz；Gaussian `.sbp`：`sigma=0.3 m`，2401 个角度样点；
- angle fan `[-30,30] deg`；
- matched outer domain `z in [-1000,1000] m`；
- mapped receiver `r'=103 m`，对应 native 物理 receiver `(r,z)=(97,0)`；
- steps `0.2, 0.1, 0.05 m`；beams `2001, 5001, 10001`；
- native comparison ATI points `(95,-1000)`, `(100,0)`, `(105,1000)`，
  即同一条精确直线；native field 在 `97 m` 减去同设置 free-field `C`
  run 以隔离一次 top-bounce reflected field，另用 `A` run 检查代表性
  arrival。

## 逐项数值结果

以下为 9 组中最坏值（结果表保存在
`results/validation/bellhop_internal_tilted_wall_poc/tilted_wall_diagnostics_summary.csv`）：

| 检查量 | 最坏值 | 结论 |
|---|---:|---|
| line residual \(|r-R_0-a z|\) | `3.5385e-12 m` | PASS |
| analytic point error | `6.4375e-12 m` | PASS |
| tangent error | `0` | PASS |
| normal error | `0` | PASS |
| specular direction error | `6.6992e-16` | PASS |
| pi-rotation direction error | `0` | PASS |
| `|(Phase_ref-Phase_inc)-pi|` | `7.1054e-15 rad` | PASS |
| reflection Amp jump | `0` | PASS |
| `|kappa|` | `0` | PASS |
| reflection `p` error | `0` | PASS |
| reflection `q` error | `0` | PASS |
| rotation `p` error | `0` | PASS |
| rotation `q` error | `0` | PASS |
| path-to-range-plane error | `8.8676e-12 m` | PASS |
| corresponding time error | `<5.92e-15 s` | PASS |
| target-ray fan-angle error | `3.9021e-05 deg` | PASS |
| target-ray time error (finite fan sampling) | `1.53e-08 s` | diagnostic |
| minimum transformed `dr` | `0.0430491 m` | strictly positive |

The analytic target ray that reaches native `(97,0)` has source angle
`-0.0166876978 deg`, wall point approximately
`(99.99985437,-0.02912549) m`, and path length
`102.9998543725 m`. The finest Bellhop fan has a nearest ray at
`-0.018 deg`; its target-angle discretization explains the small diagnostic
time residual, while the range-plane path residual remains at machine precision.

## Native ATI comparison

The finest native `A` run at `97 m` returned one selected arrival with:

| quantity | native ATI | analytic / rotated check |
|---|---:|---:|
| source angle | `-0.0167267192 deg` | target `-0.0166876978 deg` |
| receiver angle | `179.442093 deg` | analytic `179.443735 deg` |
| top/bottom bounces | `1 / 0` | required `1 / 0` |
| delay | `0.0686666891 s` | analytic `0.0686665696 s` |
| delay path length | `103.00003365 m` | analytic `102.99985437 m` |
| arrival amplitude | `0.010033899` | rotated C field `0.009708729` |

For the 3×3 C scan, native reflected-only pressure versus the transformed wall
pressure had worst-case:

- TL difference: `0.28615 dB`;
- phase difference: `0.00297 rad`;
- complex relative difference: `0.03254`.

The same offset is stable across step and beam scans and is not accompanied by
any direction or delay bias. It is retained as a limitation of evaluating the
native backward branch in Bellhop's original range chart; no amplitude
renormalization was introduced.

The ASCII `A` phase (`180 deg`) is Bellhop's arrival/reflection convention and
does not include the same propagation/source phase carried by the coherent
`C` pressure. The phase comparison above therefore uses the native C-field
reflected component, while the A run is used for independent path, bounce and
delay checks.

## 收敛与下一阶段

The wall residual and beam-state errors are unchanged at machine precision as
step decreases and beam count increases. The native angle/delay result is
consistent with the analytic mirror-source construction, and the transformed
branch is strictly range-increasing for every scanned ray.

**可以进入 sinusoidal-wall 几何验证，但有明确限制：**下一阶段应继续把
native ATI 与 rotated branch 的 reflected direction、intersection、curvature,
phase 和 `p/q` 作为硬门；native backward-range Cartesian absolute amplitude
必须继续单独报告，不能当作已解决的等幅结论。sinusoidal/PM wall、分层 SSP、
3D 和多次反射在本任务中均未实现。
