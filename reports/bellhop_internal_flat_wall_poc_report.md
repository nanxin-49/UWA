# Bellhop 2020 internal flat-wall validation-only POC 报告

## 结论

**结果：PASS。**

flat internal wall POC 已满足本阶段全部停止条件和验收要求：

- ray 在真实 `r=100 m` 交点停止；
- 物理反射直接调用 Bellhop 2020 标准 2D `bellhop.f90` 内部的 `Reflect2D`；
- pressure-release 相位只增加一次 pi；
- flat wall 的 `kappa=0`，没有额外 curvature kick；
- 反射和随后 pi rotation 都未重置、拟合或修正 `p/q`；
- 坐标接缝没有传入 receiver influence；
- `InfluenceGeoHatCart` 源码和公式均未修改；
- 9 组 step/beam-count 组合全部恢复

  \[
  H_{wall}=-P_{BH}(103\ \mathrm m)
  \]

  到 SHD 保存精度。

因此可以安全进入下一阶段 **tilted straight wall** 验证，但本任务没有实现 tilted、sinusoidal 或 PM wall。

## 实现范围

本 POC 被主动限制为：

- Bellhop 2020 `2020_11_4`；
- 2D；
- 均匀、无损 `c=1500 m/s`；
- flat pressure-release internal wall，`r_wall=100 m`；
- 每条 ray 恰好一次 wall reflection；
- mapped receiver range 包含 103 m；
- 当前 `.sbp` Gaussian source angular pattern；
- coherent TL 和 Cartesian geometric-hat influence；
- 无海底、外边界或多次粗糙面反射。

没有修改 PE、Gaussian 源定义、通信代码、正式 Bellhop 2020 源码目录或正式 executable。

## 源码结构

### validation overlay

完整的两个修改后 Fortran 文件保存在：

- `scripts/validation/support/bellhop_internal_flat_wall_poc/Step.f90`；
- `scripts/validation/support/bellhop_internal_flat_wall_poc/bellhop.f90`。

它们以官方 `AcousticsToolbox_2020/Bellhop` 文件为基线。构建脚本先复制官方源码到结果目录，再覆盖这两个文件并生成独立 executable，因此不会覆盖正式安装。

### `Step2D` / `ReduceStep2D`

`Step2D` 新增 validation-only 输入 `WallActive`、`WallRange` 和输出 `WallHit`。`ReduceStep2D` 在 Euler 和 midpoint 两次 reduction 中都计算

\[
h_{wall}=\frac{R_0-r_0}{u_r},
\]

并将它加入现有 top、bottom、SSP 和 range-segment step limits 的最小值。

初次正式扫描发现，当累计步进使最后一段极接近墙面时，Bellhop 通用 small-step guard 会把正确的 wall step 放大并导致越墙。最终实现仅对已经被选为最短 crossing event 的 `h_wall` 绕过该 guard；其他 Bellhop small-step 行为保持不变。这一修正保证了 wall event 的几何位置优先于通用最小步长。

### 原生 `Reflect2D`

flat wall 使用固定有向标架

\[
\mathbf t_w=(0,1),\qquad
\mathbf n_w=(1,0),\qquad
\kappa_w=0,
\]

它满足 Bellhop `TOP` 的标架 sense。命中后直接执行：

```fortran
CALL Reflect2D( is, WallHS, 'TOP', WallT, WallN, WallKappa, RTop, NTopPTS )
```

`WallHS%BC='V'`，所以镜面方向、vacuum pi phase、`Amp`、`p/q` 更新都由原生 `Reflect2D` 完成。没有复制反射公式，也没有附加 Kirchhoff、SSA 或其他散射项。

### 反射后 pi rotation

原生反射完成后，只对 reflected node 执行

\[
r'=2R_0-r,\qquad z'=-z,
\qquad \mathbf t'=-\mathbf t.
\]

`p`、`q`、`tau`、`Amp` 和 `Phase` 不变。该变换是 determinant 为 `+1` 的 proper rotation，所以不反转 ray-normal handedness，也不引入额外 caustic phase。

每条 ray 保存 `WallBranchStart`。trace 完成后执行等价于

```fortran
ray2D(1:NPost) = ray2D(WallBranchStart:Beam%Nsteps)
```

的打包，然后才调用现有 `InfluenceGeoHatCart`。因此 incident/transformed 两个坐标图之间的接缝从未成为真实传播 segment，也没有进入 field accumulator。

## 独立 binary 和输入

独立 executable：

`results/validation/bellhop_internal_flat_wall_poc/bin/bellhop_iwall_flat_2020.exe`

SHA-256：

`9B3DD06C8A6A6D93DF6AB70DE2C835B09D5ED199BB480CA134025A64D51C7835`

正式 Bellhop 2020 executable 仍为：

`E:/MISC/BELLHOP/AcousticsToolbox_2020/windows-bin-20201102/bellhop.exe`

其 SHA-256 仍为：

`7E7809A64C3BF734AFF6D28D0D4D52B1B4BD203D81676E3241FFD3189941B505`

validation binary 要求独立 `.iw2` sidecar：第一行为 `100`，第二行为 mapped receiver range `103`。它还在运行时拒绝非 uniform/lossless 1500 m/s、非零横向源深、错误 field/beam type 和不包含 103 m 的 post-wall receivers。

## 验证设置

- frequency：4 kHz；
- source transverse depth：0 m；
- receiver transverse depth：0 m；
- reference receiver ranges：102、103 m；
- wall：`r=100 m`；
- mapped target receiver：103 m；
- Gaussian source width used by `.sbp`：`sigma=0.3 m`；
- angle fan：`[-30,30] deg`；
- source-pattern samples：2401；
- artificial matched domain：`z in [-1000,1000] m`；
- steps：`0.2, 0.1, 0.05 m`；
- beam counts：`2001, 5001, 10001`。

每个 wall run 都与同 step、同 beam count、同 `.sbp` 的官方 Bellhop free-field run 配对，不做逐频率或逐位置拟合。

## 数值结果

全部 9 组得到相同的 103 m 轴上复压：

\[
P_{BH}(103)=
0.00485436897724867-0.00840801373124123i,
\]

\[
H_{wall}=
-0.00485436897724867+0.00840801373124123i.
\]

因此所有组合的：

- TL difference：`0 dB`；
- phase difference：`0 rad`；
- complex relative error：`0`（SHD 保存精度）。

### 全部 ray / 全部 9 组的最坏值

| 检查量 | 最坏值 | 结果 |
|---|---:|---|
| wall intersection residual | `3.5385028240853e-12 m` | PASS |
| specular direction error | `0` | PASS |
| pi-rotation direction error | `0` | PASS |
| `|(Phase_ref-Phase_inc)-pi|` | `7.105427357601e-15 rad` | PASS |
| reflection `Amp` jump | `0` | PASS |
| flat-wall `|kappa|` | `0` | PASS |
| reflection `p` error | `0` | PASS |
| reflection `q` error | `0` | PASS |
| rotation `p` error | `0` | PASS |
| rotation `q` error | `0` | PASS |
| central-ray travel-time error vs `103/1500` | `9.71445146547012e-16 s` | PASS |
| imaginary travel time | `0 s` | PASS |
| minimum transformed range increment | `0.0433012701892181 m` | strictly positive |

### 轴上最细 case 的状态

对 `step=0.05 m`、`10001 beams`、`alpha=0`：

- hit：`r=99.999999999996461 m`、`z=0`；
- incident unit tangent：`(1,0)`；
- native reflected unit tangent：`(-1,0)`；
- rotated unit tangent：`(1,0)`；
- phase：`0 -> 3.1415926535898002 rad`；
- amplitude：`1 -> 1`；
- `p=(1,0)` before/after；
- `q=(150000,0)` before/after；
- wall time：`0.066666666666669094 s`；
- mapped receiver time：`0.068666666666666626 s`。

## 保存结果

- `results/validation/bellhop_internal_flat_wall_poc/flat_wall_poc_results.mat`；
- `results/validation/bellhop_internal_flat_wall_poc/flat_wall_convergence.csv`；
- `results/validation/bellhop_internal_flat_wall_poc/flat_wall_diagnostics_summary.csv`；
- `results/validation/bellhop_internal_flat_wall_poc/cases/`：配对 `.env/.sbp/.shd/.prt/.iw2/.iwdiag`；
- `results/validation/bellhop_internal_flat_wall_poc/build_manifest.json`；
- `results/validation/bellhop_internal_flat_wall_poc/bin/bellhop_iwall_flat_2020.exe`。

用于本次构建的便携编译器下载和两个 smoke 目录在 formal run 通过后已删除；它们都是可再生成的临时文件。正式 binary、构建 manifest、overlay 源码和全部 3×3 结果均保留。独立 binary 在删除工具链后再次执行 `step=0.05 m, 10001 beams` 成功，确认运行时不依赖临时 compiler 目录。

## 项目回归

虽然 POC 没有改动生产代码，仍按仓库要求执行了两个低成本回归：

1. scalar-frequency channel：`f0=4000 Hz`、`z_tx=99 m`、`z_rx=3 m`、
   `nx=ny=128`、`xw=yw=32 m`、CPU、`rx_only`、surface reflection on、
   `enforce_1_over_R=false`。公开字段均存在，且
   `H_f=H_direct_f+H_reflect_f` 的闭合误差为 `0`。
2. wideband communication：`COMM_NX=COMM_NY=256`、
   `COMM_NF_MIN=COMM_NF_MAX=8`、`COMM_N_SYM=200`、
   `COMM_EBN0_DB_LIST=0,10`、`COMM_ENFORCE_1_OVER_R=0`。`direct_only` 和
   `direct_plus_reflect` 均完成，两个 scenario 均保留 8 点 `H_f` 和现有
   communication result 结构。

首次尝试的 128×128 wideband smoke 被生产输入检查正确拒绝，因为
`dx=dy=50/128 m` 大于 `sigma_src_m=0.3 m`；随后按采样约束改为 256×256，
没有放宽或绕过该检查。

## 是否进入 tilted straight wall

**可以。** flat stage 已证明：

1. wall event 能在真实交点参与 step reduction；
2. 原生 `Reflect2D` 可作为 internal wall 的唯一物理反射；
3. proper pi rotation 保持 `p/q/tau/Amp/Phase`；
4. post-wall branch 可与坐标接缝隔离后直接复用现有 influence；
5. 不需要修改 `InfluenceGeoHatCart` 核心公式。

下一阶段只应增加解析 tilted-line intersection 和对应常量 `t/n`，继续保持 `kappa=0`。在 tilted stage 通过前，不应进入 sinusoidal 或 PM wall。
