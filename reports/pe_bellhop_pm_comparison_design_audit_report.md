# PE rough-PM 与 Bellhop rough-PM 对比方案审查

日期：2026-09-02  
范围：方案与源码可行性审查；本轮没有运行完整 PE--Bellhop PM 对比，没有修改 PE、Bellhop 反射物理或通信代码。

## 1. 结论

最终结论：**FEASIBLE_WITH_MINOR_PREPARATION**。

当前 Bellhop 2020 validation-only parametric internal wall 已具备严格 90 度 PM 对比所需的局部几何、一次原生 `Reflect2D`、curvature `p/q` kick、proper half-turn 和 coherent receiver influence 链。PE 也已能通过 `surface_elevation_override_xy` 接收外部确定海面，并已有 `192.1875 m`、no-sponge 的完整反射链窗口证据。因此不存在需要重写任一核心物理的阻塞。

开始正式运行前仍需三个小型 validation-only preparation：

1. 增加一个**同维桥接量**：优先实现一横向维 PE/AS runner，或等价地从 y-invariant 三维 PE 场提取 `k_y=0` 分量；它不替换生产 PE。生产二维横向 Gaussian、y-invariant 海面的轴上点值保留为第二层 sensitivity check。
2. 增加一个 canonical profile mapper：只读取 seed `260001` 的既有 Fourier coefficients，用同一符号、datum 和系数为 PE 与 Bellhop 求值；不得重新随机、重归一化或单独平滑。它还必须记录 profile hash、插值误差和有限 support/FFT seam 诊断。
3. 增加一个 reflected-only flat/rough extractor，并先做常量高度符号审计和 flat source audit。粗糙 PM 结果不得用于重新拟合 `.sbp`。

六个核心问题的直接回答如下。

| 问题 | 回答 |
|---|---|
| PE 二维粗糙面与 Bellhop 二维 profile 如何映射？ | 第一阶段使用 `eta_PE(x,y)=eta_1D(x)`；Bellhop 使用同一 `eta_1D`。严格同维主量取 PE 的一横向维/`k_y=0` reduction；现有二维横向 PE 的轴上点接收作为次级检查。 |
| exact 90-degree internal wall 是否适合作为主工况？ | **是。** internal wall 已通过 vertical-tangent、fixed-PM density convergence 和有限角度 covariance。89/89.5 度仅作连续性诊断，不应替代名义严格垂直物理问题。 |
| 同一 PM realization 如何复用？ | 读取同一份 master Fourier coefficients，按同一 signed `eta(s)` 求值；PE 复制到所有 y，Bellhop 经严格 90 度刚体变换形成参数墙。两边不再生成、归一化或平滑。 |
| primary metric 是什么？ | `G=H_ref,rough/H_ref,flat`，比较 `20log10|G_PE/G_BH|` 与 `arg(G_PE conj(G_BH))`。`Q_rough/Q_flat` 与它代数等价，但增加一次 direct division，故只作闭合审计。 |
| 最小测试矩阵是什么？ | Stage 0：4 kHz、flat、单轴接收、reflected-only、严格 90 度；Stage 1：同一 fixed PM、同一设置的 flat/rough pair。只在结果可解释后加入 6/8 kHz。 |
| 当前能否开始实现？ | 可以，但应先完成上述三个 preparation；不应直接进入宽带、Monte Carlo 或总信道比较。 |

## 2. 当前代码能力与边界

### 2.1 PE

当前 `vertical_wape_propagator` 是两横向维 `(x,y)` FFT 加竖直 marching 的三维场模型。均匀介质传播核使用 `sqrt(k^2-k_x^2-k_y^2)`；生产源为

```text
psi0(x,y)=exp(-(x^2+y^2)/(2 sigma^2)),  sigma=0.3 m.
```

`vertical_channel_model_impl` 已暴露 `surface_elevation_override_xy`。`pm_surface_boundary_model` 对该数组做尺寸、实数和有限值检查后直接使用，不重新生成或归一化海面。当前 normal-incidence Kirchhoff 路径为

```text
psi_ref(x,y) = -psi_inc(x,y) exp(i 2 k eta(x,y)).
```

因此同一 fixed profile 可以无生产物理修改地注入。已有窗口审计表明，完整 Tx--surface--Rx 反射链必须以实际 `192.1875 m`、`984 x 984`、`dx=dy=0.1953125 m`、sponge off 作为严格基线；50 m/default sponge 不可用于本次比较。

### 2.2 Bellhop

当前保留的 generic PM overlay 接受成对参数点，计算 finite-support segment intersection、单位 tangent/normal 和 signed geometric curvature，在真实 hit 处只调用一次 Bellhop 2020 原生 `Reflect2D`。随后对反射分支作 proper 180 度 rotation，保存 `p/q/tau/Amp/Phase`，并把隔离后的正-range branch 送入未修改的 `InfluenceGeoHatCart`。

当前可可靠取得：

- SHD 中指定 receiver range/depth 的 coherent reflected pressure；
- `.iwdiag` 中每条成功 ray 的 wall hit、tangent/normal/curvature、incident/reflected direction、`RN/RM`、`p/q`、phase、Amp、tau 和 transformed range 状态；
- wall-hit success/reject、grazing、重复命中和非单调 branch 统计。

跨 chart seam 的通用 eigenray/arrival/PDP 语义尚未作为正式接口验证。因此个别 arrival、完整 arrival ordering 和全域 post-wall field 只能作 diagnostic，不能先于 coherent SHD 与局部 ray audit 成为 hard gate。

## 3. 维度对应审查

### 3.1 方案 A：`eta_PE(x,y)=eta_1D(x)`

这是第一阶段最干净的**表面几何映射**。同一条 profile 既是 Bellhop 的二维边界，也是 PE 中沿 y 平移不变的柱面海面。它消除了“PE 与 Bellhop 看见不同 realization”这一不确定性，也不会把二维 PM surface 的任意 transect 误称为完整二维海面。

但是，对当前生产轴上 PE 点值，它还不是严格同维：Gaussian 在 y 方向仍衍射。虽然 y-invariant phase screen 不混合 `k_y`，每个 `k_y` 分量的竖直传播数仍是

```text
k_z=sqrt(k^2-k_x^2-k_y^2),
```

所以轴上场是所有 `k_y` 的相干叠加；Bellhop 2D 对应的是一个传播平面/线源场。rough/flat normalization 会消除固定源强、许多共同的 y 衍射和基线 spreading，但不能严格消除 rough interaction 对 `k_y` 的依赖。`sigma=0.3 m` 的典型角谱宽度约为 `1/(k sigma)=0.199 rad`，即约 11.4 度，不能先验视为纯 `k_y=0`。

结论：方案 A 必须采用，但生产三维 PE 的轴上点值应是**次级物理 sensitivity**，不是唯一主参考。

### 3.2 方案 B：二维 PM surface 的中心 transect

不推荐作为第一阶段。它有三个额外歧义：

1. Bellhop 只保留 `y=0` 一条切线，而 PE 仍受 y 方向不同 roughness 的相干耦合；
2. 中心切线的 RMS height/slope/curvature 一般不等于二维 realization 的整体统计；
3. 任意切线无法把差异明确分配给 dimensionality 还是 reflection physics。

该方法以后可作为“真实二维 surface 的一条观测切片”诊断，但弱于 y-invariant 的 controlled comparison。

### 3.3 推荐的方案 C：显式同维 reduction

主比较采用一横向维 PE/AS validation runner：

```text
psi0(x)=exp(-x^2/(2 sigma^2))
  -> 1-D FFT, sqrt(k^2-k_x^2) propagation
  -> -exp(i 2 k eta_1D(x))
  -> 1-D FFT propagation to Rx.
```

它与 Bellhop 2D 处于同一空间维数，又保持生产 Gaussian 的 x 截面、精确平方根传播和相同 Kirchhoff phase-screen 公式。另一种数值等价实现是对 y-invariant 的完整 PE 场取 `k_y=0` 投影；由于介质与 phase screen 均不依赖 y，`k_y` 守恒。两者应先在 flat/rough 小例上相互核对。

推荐分层为：

```text
Tier 1（主）：1-transverse PE / ky=0  <->  Bellhop 2D
Tier 2（次）：当前 2-transverse PE，eta(x,y)=eta(x)，轴上点值
Tier 3（以后）：真正 eta(x,y) 的 PE；不能再称与 Bellhop 2D 严格同物理
```

这不是要求修改生产 PE，而是增加一个可审计的 validation bridge。

## 4. 严格 90 度坐标映射

第一阶段应直接采用 exact 90-degree internal wall。设物理二维截面坐标为 `(x,z)`，z 正向向下，Tx 为 `(0,z_tx)`，canonical signed profile 的物理边界为

```text
Gamma_P(s) = [s, eta(s)].
```

使用 source-centered proper rotation

```text
[r; zeta] = [0 -1; 1 0] * ([x; z] - [0; z_tx]).
```

则

```text
Tx_B       = [0,0]
wall_B(s)  = [z_tx-eta(s), s]
Rx_direct  = [z_tx-z_rx, 0].
```

反射后对整个 branch 绕名义墙点 `(R0,0)` 作已经验证的 proper half-turn

```text
r'    = 2 R0-r
zeta' = -zeta,
```

flat case 的 receiver image range 为 `z_tx+z_rx`。该变换只重表示反射后 branch，不再次反射、不重算 curvature、不改变 `p/q/tau/Amp/Phase`。

89/89.5 度不是更“正确”的主工况。以 99 m 传播距离估算，偏离竖直 0.5/1 度分别产生约 `0.864/1.728 m` 的横向漂移，可能命中另一段 PM profile；同时会改变局部入射角和 phase factor。有限角度结果可作为 90 度连续性 smoke check，但不值得用这种真实物理偏差换取 native ATI 兼容，因为本任务不需要 native ATI 作为主参考。

现有 flat/unfolded 和 PM internal-wall 证据使用 `z_tx=100 m`、`z_rx=3 m`，对应 flat image range `103 m`。最小方法学 POC 建议继续使用该 validation geometry，以复用现有基线；它没有物理海底。当前项目代表性安装 `z_tx=99 m`、`z_rx=3 m` 应在方法通过后单独重跑，此时 image range 为 `102 m`，不得混用两个 geometry 的压力结果。

## 5. 同一 fixed PM realization 的复用

### 5.1 canonical 数据

第一阶段固定使用：

| 项目 | 数值 |
|---|---:|
| seed | `260001` |
| U | `6 m/s` |
| span / master N | `160 m / 4097` |
| requested / realized Kmax | `0.5 / 0.471238898 rad/m` |
| master spacing | `0.0390625 m` |

权威输入应是 `fixed_pm_fourier_coefficients.csv`，master profile CSV 是其可核对采样。loader 必须：

- 保持所有 Fourier coefficient 与 phase 不变；
- 不减去当前离散均值，不重新设 Hs，不为任一 solver 单独缩放；
- 记录源文件 SHA-256、求值区间、采样点、插值方式和输出 hash；
- 用同一 signed `eta` 数值生成 PE 与 Bellhop 输入。

现有 master 离散均值约 `-1.69e-5 m`，足够接近零但仍不得在一边单独 recenter。端点的 height、slope 和二阶导数一致到约 `1e-16`，符合该 Fourier series 的 160 m 周期性。

### 5.2 PE 映射

在 PE x 网格上直接由 Fourier coefficients 求 `eta_1D(x)`，再构造

```text
eta_PE = repmat(eta_1D(x), ny, 1).
```

这比先从 4097 点 profile 做低阶插值更干净；如必须插值，应同时保存相对 Fourier 求值的 height/slope/curvature 误差。y 方向不得再加入随机项。

### 5.3 Bellhop 映射

在选定参数样本 `s_j` 上由同一 coefficients 求 `eta(s_j)`，然后用

```text
Gamma_B(s_j)=[R0-eta(s_j), s_j]
```

生成 parametric wall。`N=4097` 作为主设置，`N=2049` 只用于离散收敛预算。不能把 PE 网格采样与 Bellhop wall density 混为一项物理带宽变化。

### 5.4 sign 与 datum

本项目当前 PM internal-wall 文件实际采用 `r=R0-eta`，PE normal phase screen 采用 `+2k eta`。第一阶段应把 CSV 中的 `eta` 明确定义为项目内部 signed profile，并在两侧原样使用；不能在一边改成 `-eta`。正式 PM 前增加一个小常量高度 `eta0` 回归，验证 rough/flat phase 的符号和 wall path change 一致。该测试只确定坐标/phasor convention，不是校准。

### 5.5 160 m support 与 192.1875 m PE 窗口

不建议在 `s=+-80 m` 后直接接常量 flat surface：当前端点 slope 约 `0.06458`，会制造非物理 kink 和集中 curvature。首选处理为：

1. 用同一 160 m Fourier series 在整个 PE/Bellhop 所需区间周期求值；不增加新频率，不平滑；
2. 两边使用完全相同的有限区间，并记录所有 Bellhop hit 是否落在 central master support；
3. PE 必须检查 surface-incident、surface-reflected 和 receiver-reflected 的 outer-5% energy，以及 FFT seam 附近的加权能量。

注意，192.1875 m 不是 160 m 的整数倍，周期求值后的两端一般不会在 PE FFT seam 精确相接。若 seam 能量通过既有严格门槛，则它只是未被有效照亮的数值边界；若失败，不得用任意 taper“修好”。下一选择应是使 x window 成为 160 m 整数倍（优先 320 m；y 仍可保留已验证宽度），或生成一条物理带宽相同、span 与 PE 周期窗一致的新 canonical realization 并重新完成 Bellhop density preflight。后者是新 realization，不属于当前最小 Stage 1。

## 6. source matching

现有 Bellhop `.sbp` 为

```text
D(theta,f)=cos(theta) exp(-(k sigma sin(theta))^2/2),  sigma=0.3 m,
```

2401 个角样本并使用 `-120 dB` floor。flat 4--8 kHz 的轴上 `Q=H_ref/H_direct` 已达到最大约 `4.94e-4 dB / 1.83e-3 rad`，5001/10001 beams 的差异约 `1.38e-6 dB / 9.66e-6 rad`。因此第一阶段可以原样复用该 `.sbp`，不应针对 PM 调参。

现有证据只证明轴上相对量很强，并没有证明 Bellhop 2D `.sbp` 严格重现生产三维 Gaussian 的完整横向场；8 kHz 外侧约 `0.061 rad` 的差异仍在。因此 Stage 0 应补充一横向维 PE/`k_y=0` 与该 `.sbp` 的 flat on-axis 和少量 offset audit。若它不通过，应修正“维度对应的解析 Jacobian 定义”，而不是用 PM 结果拟合角谱。

## 7. 比较量与验收层级

### 7.1 primary metric

每个 solver 内先做 flat normalization：

```text
G_PE = Href_PE_rough / Href_PE_flat
G_BH = Href_BH_rough / Href_BH_flat.
```

跨 solver 报告

```text
DeltaTL_model    = 20 log10 |G_PE/G_BH|
DeltaPhase_model = arg(G_PE conj(G_BH))
complex_error    = |G_PE-G_BH| / max(|G_BH|,epsilon).
```

这是最干净的 primary metric。它消除固定 source scale、2D/3D baseline spreading、flat pressure-release coefficient 和大部分共同路径传播，但不会消除角谱维度或两种 rough-reflection physics 的真实差异。

`Q=Href/Hdirect` 时，`Q_rough/Q_flat=Href_rough/Href_flat`，因为同一 solver 内 direct 项代数消去。直接使用 `G` 少一次可能病态的 direct division；`Q` 形式只用于检查代码闭合和与旧 flat 报告衔接。

### 7.2 secondary metrics

- flat-normalized reflected phase 与复数比的网格/beam/profile convergence；
- Bellhop geometric tau/path 与 flat 的变化；多频后再与 PE roughness-induced group delay 比较；
- 少量 normalized receiver-offset profile 或 reflected angular spectrum；
- Tier 1 与 Tier 2 PE 结果之差，用于量化残余 y diffraction/dimensionality。

单频 4 kHz 不能给出数值 group delay；Stage 1 只能报告相位和 Bellhop tau。PE group delay必须等 Stage 2 至少有邻近频点后由一致的 unwrapped phase 定义得到。

### 7.3 diagnostic-only

- 任一 solver 的 absolute reflected TL/phase；
- Bellhop 单 ray amplitude、hit distribution、local incidence、normal、curvature、`RN/RM/p/q`；
- arrival ordering、PDP 和 total channel；
- native backward-range amplitude baseline。

第一阶段只比较 reflected branch。`H_total=Hdirect+Href` 会让强 direct path 掩盖 reflected-model discrepancy，应在 reflected-only 结果可解释后再加入。

## 8. 可比性矩阵

| Quantity | PE definition | Bellhop definition | directly comparable? | normalization needed? | hard/soft/diagnostic |
|---|---|---|---|---|---|
| reflected complex pressure | 连续 phase-screen 后传播到 Rx 的复场 | transformed branch 经 `InfluenceGeoHatCart` 的 SHD coherent pressure | absolute 值否；同维且 flat-normalized 后可比 | 必须 rough/flat；固定相位 convention | absolute diagnostic；normalized soft cross-model |
| rough/flat reflected ratio | `Href_rough/Href_flat` | `Href_rough/Href_flat` | 是，主比较量 | 已内含 | 数值收敛 hard；跨物理模型 soft |
| reflected TL | `-20log10|Href|` 或 `20log10|G|` | 同定义 | absolute 否；rough-induced 是 | rough/flat | rough-induced primary soft |
| reflected phase | `arg(Href)` 或 `arg(G)` | 同定义 | absolute 需 convention；`arg(G)` 可比 | rough/flat + unwrap rule | primary soft；符号/闭合 hard |
| travel time | PE 无单一 ray；可由 phase slope/波包定义 | ray `tau`/path length | 单频不可直接；窄带 aggregate 可软比较 | 减 flat delay | Bellhop diagnostic；多频 secondary soft |
| group delay | `-d arg(Href)/d omega` 或对 `G` 求导 | coherent SHD phase slope；ray tau 另列 | 多频时可比 coherent delay | flat-normalized、统一 unwrap | Stage 2 secondary soft |
| reflection point | phase screen 上所有受照 surface points | 每条 ray 的 wall hit | 否，无一一对应 | 不适用 | diagnostic |
| local incidence angle | normal phase mode固定使用 factor 2；可由入射谱诊断 | 每条 ray 的 `u dot n` | 不是相同内部变量 | 不适用 | preflight/diagnostic |
| surface tangent/normal | 可由同一 `eta` 导数计算，但 phase screen 不使用 local frame | parametric wall 的单位 `t/n`，送入 `Reflect2D` | 几何输入可核对；物理作用不等价 | 同一 profile/sign | geometry hard；field effect diagnostic |
| curvature | 可由 `eta''/(1+eta'^2)^(3/2)` 计算，但当前 phase screen 不使用 | signed kappa，进入 RN 与 `p/q` kick | 几何值可核对；作用不等价 | 同一 samples/continuous reference | geometry hard；physics diagnostic |
| angular spreading | PE receiver/surface FFT spectrum | beam fan 的 coherent field/angle distribution | 仅 normalized aggregate 可比 | 轴上或总能量归一化 | secondary soft |
| PDP / arrivals | 宽带复场 IFFT 得连续/离散窗 PDP | coherent frequency response 或离散 ray arrivals | 只有 aggregate PDP 可软比；单 ray 无 PE 对应 | 同带宽、窗和相位 reference | later soft/diagnostic |
| total channel | `Hdirect+Href` | direct 与 wall-reflected coherent sum | 技术上可算但会掩盖反射差异 | 先分别 flat/rough | Stage 2/3 diagnostic，后续应用量 |

## 9. 物理模型差异与误差预算

PE 当前 rough branch 是 local phase screen 加反射前后 diffraction；它不使用 local tangent、curvature 或 ray-by-ray shadowing。Bellhop 使用真实参数曲面、局部 specular direction、pressure-release `Reflect2D`、curvature `p/q` kick 和 beam coherent accumulation。两边都正确时 `G_PE != G_BH` 仍可能是正常的 model discrepancy。

因此不能预注册任意 `0.1 dB/0.01 rad` 跨模型门槛。应先建立经验数值预算：

```text
E_num <= E_PE(window/grid/projection)
       + E_BH(profile/step/beam)
       + E_flat(source/convention).
```

可用的先验量级是：

- PE 既有 160--192 m：4 kHz reflected TL/phase 约 `0.00154 dB / 0.000612 rad`；3--5 kHz worst 约 `0.0192 dB / 0.00238 rad`。新 fixed PM 仍需自己的 192 m edge check，不能直接借用该数值作为证明。
- Bellhop fixed PM `N=2049 -> 4097`：receiver 约 `0.00677 dB / 0.000311 rad`，complex relative error约 `8.39e-4`；step/beam 项应在 Stage 0/1 重新记录。
- flat PE--Bellhop 轴上 `Q` 先验差约 `4.94e-4 dB / 1.83e-3 rad`（4--8 kHz worst）。

解释规则：

1. 随 window/grid/profile/beam 收敛而下降的差异归为数值误差；
2. 符号、flat ratio、一次 pi phase、receiver selector 或 component closure 错误属于实现失败；
3. 在两边各自收敛、明显超过经验数值预算且对小配置变化稳定的残差，归为 Kirchhoff phase screen 与 local-specular Bellhop 的 model discrepancy，而不是强行校准；
4. 若差异对 PE window seam、Bellhop profile density 或 beam count 不收敛，则暂不作物理解释。

## 10. fixed realization preflight

由当前 master profile 直接计算：

| 量 | 数值 | 判断 |
|---|---:|---|
| acoustic wavelength, 4 kHz | `0.375 m` | 基准 |
| acoustic k | `16.7552 rad/m` | 基准 |
| realized PM Kmax | `0.471239 rad/m` | `Kmax/k=0.02813` |
| shortest surface wavelength | `13.3333 m` | 约 `35.6 lambda` |
| eta RMS / min / max | `0.18639 / -0.40930 / 0.49946 m` | 高度相位不是弱扰动 |
| RMS / max slope | `0.053980 / 0.116674` | 小坡度；最大约 6.65 度 |
| RMS / max geometric curvature | `0.020905 / 0.044594 1/m` | 有限且已被 Bellhop density audit 解析 |
| minimum radius of curvature | `22.424 m = 59.8 lambda` | local curvature 对 Bellhop 可解析 |
| `lambda max|kappa|` | `0.016723` | 远小于 1 |
| nominal vertical min `|u dot n|` from max slope | `0.993262` | 非 grazing |
| existing broad-fan/local-covariance min `|u dot n|` | `0.8022 / about 0.998` | 已有接受 ray 也不接近 grazing |
| PE samples / shortest surface wavelength | 约 `68.3`（dx `0.1953125 m`） | height sampling 充分 |
| BH master samples / shortest surface wavelength | 约 `341` | boundary sampling 充分 |
| RMS phase modulation `2 k sigma_eta` | 约 `6.25 rad` | 不能使用 small-phase 解释 |

该 realization 适合做**受控模型差异比较**：坡度小、曲率半径远大于波长、无 grazing，Bellhop local-specular geometry 数值可解析。它同时并不是“两个模型必然等幅相”的弱粗糙例，因为高度引起的两程相位 RMS 约 6.25 rad，单 realization 的相干结果可能对干涉很敏感。

Kirchhoff/2K 方面，`Kmax/k=0.0281` 表明单个表面谱移对应的小角度约 1.61 度，曲面变化尺度也远大于 acoustic wavelength；这些支持高频局部处理。另一方面，生产 Gaussian 本身具有约 11.4 度角谱宽度，PE normal mode 把 phase factor 固定为 2，而 Bellhop 对每条 ray 使用真实 local normal。现有 2K 验证只能证明近法向替换的趋势和量级，不能把当前 phase screen 宣称为 full-wave truth。此差别正是本对比应量化的物理模型差异之一。

结论：不需要改变 `Kmax`，也不应为“提高一致性”平滑当前 realization。

## 11. 最小执行矩阵

### Stage 0：flat baseline 与 convention regression

| 项目 | 设置 |
|---|---|
| 频率/介质 | `4 kHz`, uniform `c=1500 m/s` |
| geometry | validation baseline `z_tx=100 m`, `z_rx=3 m`; strict 90-degree; reflected-only |
| PE | Tier 1 一横向维或 `k_y=0`; x width `192.1875 m`, dx `0.1953125 m`, sponge off |
| PE secondary | 现有 `984 x 984` 二横向维，sponge off |
| Bellhop | current internal flat wall, step `0.05 m`, current Gaussian `.sbp`, `10001` beams；`5001` 作一次 convergence pair |
| cases | flat `eta=0`；另加小常量 `eta0` 只审计 rough-phase/path sign |
| hard checks | receiver selector、103 m image path、一次 pi、component/ratio closure、constant-height phase sign、Tier-1 flat source audit |

Stage 0 不以 absolute amplitude 为 gate。若 Tier 1 与当前 `.sbp` 的 flat normalized source footprint不一致，应先定位 2D angular Jacobian；不得看 PM 结果后拟合。

### Stage 1：single fixed PM at 4 kHz

| 项目 | 主设置 | 数值审计 |
|---|---|---|
| realization | seed `260001`, same Fourier coefficients | hash/sign/datum一致 |
| PE Tier 1 | `192.1875 m`, dx `0.1953125 m`, no sponge | outer-energy/seam gate；必要时 integer-period x window |
| PE Tier 2 | `984 x 984`, `eta(x,y)=eta(x)`, no sponge | 与 Tier 1 的 dimensionality difference |
| Bellhop wall | `N=4097`, exact 90-degree parametric wall | `N=2049` convergence |
| Bellhop ray settings | step `0.05 m`, `10001` beams | `5001` pair；不做大扫描 |
| receiver | one on-axis receiver | 显式 SHD range/depth selector |
| pair | flat + rough | 两边各自共享全部非 surface 设置 |
| primary output | `G_PE`, `G_BH`, DeltaTL/DeltaPhase/complex error | 与经验 numerical budget 比较 |
| diagnostics | PE edge energy；BH hit/t/n/kappa/RN/RM/p/q/tau/grazing/reject | 不做 PE point-to-ray 配对 |

若 Stage 1 的两边各自收敛且 `G` 差异稳定，结果就是有价值的模型差异，即使不接近零也不构成 internal-wall 失败。

### 后续层级

- Stage 2：只有 Stage 1 可解释后，增加 6/8 kHz；每个频率仍做 matching flat/rough pair。若要 group delay，需要小而一致的频率邻域，不能只用三个稀疏点作高精度导数。
- Stage 3：理解单 realization 后才增加 seed/ensemble；不在本任务范围。
- total channel、通信指标、真正二维 random surface 和 Monte Carlo 均晚于 reflected-only comparison。

## 12. 最小实现边界

下一步实现只需要新增 validation/comparison orchestration 与 helper，不修改：

- `vertical_wape_propagator` 的生产 square-root marching；
- `pm_surface_boundary_model` 的 Kirchhoff 公式；
- Bellhop `Reflect2D`、`InfluenceGeoHatCart` 或 p/q 公式；
- Gaussian 生产源或通信链。

建议新增内容仅为：

1. 一个 canonical fixed-PM loader/mapper；
2. 一个一横向维 PE/`k_y=0` validation reference；
3. 一个 flat/rough reflected-ratio comparison entrypoint及表格/元数据输出。

## 13. 最终判断

**FEASIBLE_WITH_MINOR_PREPARATION**。

exact 90-degree internal-wall、同一 fixed band-limited PM profile、4 kHz、单轴接收、reflected-only、flat/rough normalization 是正确的第一阶段。当前 `U=6 m/s`、`Kmax=0.471 rad/m` realization 的 slope、curvature、grazing 和 Bellhop boundary resolution 均适合该路线；其强 height phase 意味着结果应被解释为 Kirchhoff phase screen 与 local-specular Gaussian-beam 模型的受控 discrepancy audit，而不是预设必须相等的回归。

完成三个小 preparation 后即可实施 Stage 0/1。无需回退到 89/89.5 度、无需 Bellhop3D、无需修改 PE/Bellhop 核心物理，也不应提前进入宽带或 Monte Carlo。
