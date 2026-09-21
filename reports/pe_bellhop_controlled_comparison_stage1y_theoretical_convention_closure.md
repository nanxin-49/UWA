# PE--Bellhop controlled comparison — Stage 1Y theoretical convention closure

日期：2026-09-10  
范围：只读源码、既有审计结果和既有解析闭合数据；本轮没有调用 PE、Bellhop 或修改物理代码。  
最终判定：**CONVENTION_THEORETICALLY_CLOSED**

## 1. 结论摘要

Stage 1Y 已从源码和数学定义独立闭合 PE/Bellhop 的复数约定：

1. 当前 PE helper 使用 `exp(-i*omega*t)` 的时间谐波约定；其 forward
   reduced-envelope 算子严格等于 `exp(+i*(kz-k0)*L)`，恢复载波后为
   `exp(+i*k0*L)`。
2. Bellhop 2020 的 coherent influence 直接累加
   `Amp*exp(-i*(omega*tau-phase))`。因此无额外 beam phase 的自由射线 raw
   SHD 场具有 `exp(-i*k*L)` 空间符号；vacuum 反射在 `Reflect2D` 中只增加
   一次 `+pi`。
3. SHD reader 按交错 real/imag 直接构造 `real+1i*imag`，没有隐藏
   conjugation 或 complex transpose。X-source 的 `ScalePressure` 因子为实数，
   不改变空间相位符号。
4. 因而当前固定的绝对场转换
   `B_abs = conj(B_raw)`（`phase_sign=-1`）是必要且有源码依据的；对粗糙/平面
   比值，正式 comparison field 再定义为
   `G_BH_comparison = conj(G_BH_abs) = G_BH_raw`，以便与 PE 的
   `exp(+i*2*k*eta)` screen 比较。这是固定的、与 case 无关的 convention
   变换，不是按误差择优、拟合或幅相校正。
5. 既有 `eta=+0.05,0,-0.05 m` constant-height audit 在先写理论预测后验证了
   `G_PE=G_BH_raw=exp(+i*2*k*eta)`；Bellhop ratio phase residual 仅约
   `2.56e-5 rad`。既有 free-field/incident audit 也显示 PE--AS 约
   `3.7e-13`，Bellhop--AS 约 `6.1e-3`，且所有预注册 4 kHz gates 通过。

因此：

```text
Stage 1Y = CONVENTION_THEORETICALLY_CLOSED
Stage 1 convention-fixed regression = PASS
Stage 2 = unlocked by theory gate (本轮不执行 Stage 2)
```

## 2. 冻结 provenance

| 项目 | 当前冻结值/来源 |
|---|---|
| 频率/介质 | 4 kHz；均匀 `c=1500 m/s` |
| 几何 | source-to-wall 100 m；wall-to-Rx 3 m；flat image path 103 m |
| PE | 1-transverse helper；`W=192.1875 m`、`nx=984`、`dx=0.1953125 m`、step `0.05 m` |
| source/run | Gaussian sigma `0.3 m`；Bellhop Cartesian line source `X`；coherent `C` |
| `.sbp` | 2401 samples、`[-30,30] deg`、`-120 dB` floor；未重拟合 |
| Bellhop | `E:\\MISC\\BELLHOP\\AcousticsToolbox_2020\\windows-bin-20201102\\bellhop.exe` |
| wall comparison | generic internal wall；single native `Reflect2D`；`T(r,z)=(200-r,-z)` 仅坐标重映射 |
| Stage-0 frozen floors | `T_E=0.0108994`；`T_phi=0.0108992 rad`；`T_TL=0.0107979 dB` |

本轮只读取当前 worktree 的源码/报告以及既有 external project result；没有读取或使用
`cash/`。当前 official source fingerprints（用于复核 provenance）为：

| 文件 | SHA-256 |
|---|---|
| official `influence.f90` | `8188087F79548C4B640DD4C2C1AD547A5DC46DBA729FD42154A491039F544C1C` |
| official `ReflectMod.f90` | `0B1ED8FCD3C156A46C9AA54A935028A25737CCC8AB80E5F90B1FF84DB85141EA` |
| official `Step.f90` | `C569C37F50ADAFD4D881D4817034DF8760FB8FB29DB17F127D8933920A62D23C` |
| official `bellhop.f90` | `A5B58C10B573EC2CB39FE6BC304F12DA08F516EA89CF621F18A6B13A193C968E` |
| PE 1-D helper | `2C653F959108EB1B8E3821461B02063FEF12AD3C8A8D678FD2BA4BDFEF2F3235` |
| SHD unfolded reader | `C9B43B50CBACC995D566EA5B852E784466B219924C560ABB4331AB59C472000C` |

## 3. PE convention audit

### 3.1 源码链路

文件 `scripts/validation/support/run_pe_1d_surface_reflection_validation.m`：

| 位置 | 证据 | 含义 |
|---|---|---|
| lines 26--29 | `kx=...`、`kz=sqrt(complex(k0^2-kx.^2,0))`、`source=exp(-0.5*(x/sigma).^2)` | 离散横向角谱和实 Gaussian 初场 |
| line 30 | `phase_screen = surface_reflect_coeff .* exp(1i*2*k0*eta)` | pressure-release `-1` 与 rough phase `+2*k0*eta` |
| lines 51--58 | `fr=exp(-1i*0.5*ds*(k0-kz))`；两半步得到 `exp(+1i*ds*(kz-k0))` | forward reduced PE propagation |
| lines 29--42 | `surface_incident_field`、`surface_reflected_field` 直接保存 | 没有保存阶段 conjugation |

控制主程序 `scripts/validation/validate_pe_bellhop_controlled_comparison.m`
的 independent AS reference（约 lines 372--379）为：

```matlab
field = -ifft(fft(source) .* exp(1i*(z_tx+z_rx)*(kz-k0)));
```

因此 PE reduced field 的解析符号是：

```text
H_PE,reduced(L) = exp(+i*(kz-k0)*L)
H_PE,physical(L) = exp(+i*k0*L) * H_PE,reduced(L)
```

### 3.2 时间谐波与载波

`scripts/validation/validate_pe_phase_convention_uniform_vertical.m` 的审计将
时间合成为 `exp(-i*omega*t)`，并在 operator identity 中验证：

```text
exp(-i*d*kappa^2/(kz+k0)) = exp(+i*d*(kz-k0))
```

该脚本选择 `carrier_positive=exp(+i*omega*d/c)`，正载波候选通过，
`carrier_negative=conj(carrier_positive)` 不作为 PE forward convention。因此
PE 的正向传播空间相位是 `+kL`，并不是 `-kL`。

### 3.3 FFT、转置和保存

- PE helper 使用 MATLAB `fft/ifft` 配对；其传播乘子已在上表和 operator
  identity 中直接核对，不能把 MATLAB FFT 的内部符号误认为物理 conjugation。
- 当前 `a(:).'`、`x(:).'` 均为 non-conjugate transpose；没有使用 complex
  `'` 把场偷偷共轭。
- `surface_incident_field`、`surface_reflected_field` 和 MAT 中 PE fields
  以当前 complex array 直接保存/读取，没有额外相位变换。

### 3.4 PE pressure-release screen

对平面压力释放面，`surface_reflect_coeff=-1`；对固定 elevation `eta`：

```text
H_PE,rough = - exp(+i*2*k*eta) H_PE,flat
G_PE      = exp(+i*2*k*eta)
```

这给出后文 constant-height closure 的理论预测，不是从观测误差反推。

## 4. Bellhop 2020 convention audit

### 4.1 Coherent beam phase

官方 `E:\\MISC\\BELLHOP\\AcousticsToolbox_2020\\Bellhop\\influence.f90`：

- `InfluenceGeoHatCart` lines 397--479 根据相邻 ray points 构造 ray tangent/normal、
  插值 `q`/`tau`，并调用 `ApplyContribution`。
- `ApplyContribution` lines 625--644 的 coherent 分支为：

```fortran
U = U + CMPLX( Amp * EXP( -i * ( omega * delay - phaseInt ) ) )
```

即 Bellhop raw coherent field 的相位为 `-omega*tau + phaseInt`。均匀介质、
无边界相位时，`tau=L/c`，所以：

```text
B_raw(L) = A(L) exp(-i*kL)
```

这是 Bellhop raw SHD 的空间符号；X-source 的幅度系数不会改变这一结论。

### 4.2 Native `Reflect2D` 相位和 beam state

官方 `Bellhop/ReflectMod.f90`：

| 位置 | 证据 | 结论 |
|---|---|---|
| lines 50--61 | `Th=dot(ray tangent,nBdry)`；`t_ref=t_inc-2*Th*nBdry` | 原生镜面方向 |
| lines 74--86 | `RN=2*kappa/(c^2*Th)`，`RM=Tg/Th`，再合并梯度项 | 原生 curvature/SSP beam kick |
| lines 95--98 | `tau` 不变；`p=p+q*RN`；`q` 不变 | reflection 后 beam state 原生更新 |
| lines 102--108 | vacuum：`Amp` 不变；`Phase=Phase+pi` | pressure-release 只增加一次 `pi` |

当前 validation overlay 只覆盖 wall intersection/segment 和 post-wall branch
隔离；官方 `ReflectMod.f90`、`influence.f90` 未被替换。proper half-turn
`(r,z)->(200-r,-z)` 只映射坐标和方向，保持 `p/q/tau/Amp/Phase`，不再次调用
`Reflect2D`、不重新算 curvature、也不再加 reflection coefficient。

### 4.3 SHD、reader 和 source geometry

- 官方 `Matlab/ReadWrite/read_shd_bin.m` lines 131--132 读取交错 float 并构造
  `temp(odd)+1i*temp(even)`；没有 conjugation。
- validation reader `scripts/validation/support/read_bellhop_shd_unfolded_vertical.m`
  lines 24--29 同样 direct-read；range record 只负责定位。
- `select_bellhop_shd_pressure_at_range_vertical.m` 要求 target range/depth
  在 tolerance 内各唯一匹配，避免 nearest-column 或 hidden receiver choice。
- `write_bellhop_unfolded_gaussian_env_vertical.m` 写入 `C *X`；既有 `.sbp`
  采用 `cos(theta)*exp[-(k*sigma*sin(theta))^2/2]`，source factor 对 X 为实数。
  原历史 point-source `R` 的 range-dependent factor 不进入当前 frozen comparison。

因此 Bellhop 的独立预期是：

```text
B_BH,raw(L) = A(L) exp(-i*kL)
```

而项目为了和 PE 的正向空间相量作绝对场对照，固定采用：

```text
B_BH,abs(L) = conj(B_BH,raw(L)) = A(L) exp(+i*kL)
```

主驱动中的 `local_stage0_convert(raw, phase_sign=-1)` 正是这个固定转换；
它不是按本次 rough error 选择的 transformation。

## 5. Independent analytic closure A — constant-height reflection

### 5.1 先验理论预测

令 flat image path 为 `L0=103 m`，固定平面高度为 `eta0`。图像路径为：

```text
L_eta = L0 - 2*eta0
```

PE pressure-release screen：

```text
H_PE,rough/H_PE,flat = exp(+i*2*k*eta0)
```

Bellhop raw propagation：

```text
G_BH,raw = B_raw(L_eta)/B_raw(L0)
          = exp(-i*k*(L_eta-L0))
          = exp(+i*2*k*eta0)
```

所以 raw ratio（也是固定比较层最终使用的 ratio）与 PE 的理论符号相同；
绝对场若先应用 `B_abs=conj(B_raw)`，则
`G_BH,abs=exp(-i*2*k*eta0)`，主驱动再以固定
`G_BH,comparison=conj(G_BH,abs)` 恢复 raw ratio。所有步骤均由定义决定。

### 5.2 既有数据验证

来源：`reports/pe_bellhop_pm_constant_height_sign_audit.md`；原始逐行数据：
`E:\\MISC\\CARPE3D_matlab\\Explain\\results\\validation\\pe_bellhop_pm_constant_height_sign\\stage0d_constant_height_rows.csv`。
该 audit 在读取数据前已固定上面的理论 prediction。

| `eta0` (m) | 预测 phase (rad) | PE phase error | BH raw ratio phase error | PE--BH phase error | BH TL ratio (dB) | path-time error (s) | 状态 |
|---:|---:|---:|---:|---:|---:|---:|:---|
| +0.05 | +1.6755161 | `1.97e-15` | `-2.55696e-5` | `+2.55696e-5` | `+0.0084364` | `-4.16e-17` | PASS |
| 0 | 0 | 0 | 0 | 0 | 0 | `-4.16e-17` | PASS |
| -0.05 | -1.6755161 | `-9.16e-16` | `+2.55865e-5` | `-2.55865e-5` | `-0.0084295` | `-4.16e-17` | PASS |

CSV 中 `eta0=+0.05` 的 raw Bellhop ratio 为
`-0.1046045845 + 0.9954909966i`，理论值为
`-0.1045284633 + 0.9945218954i`；`eta0=-0.05` 的 raw ratio 为
`-0.1044016472 - 0.9935598643i`，理论值为
`-0.1045284633 - 0.9945218954i`。差异为有限 beam/source 数值误差，符号和
相位方向均符合理论；没有 per-case conjugation 或 phase fitting。

该 constant-height test 是非 sinusoidal 的独立证据，因此排除了仅凭
`A cos(Ks)` 对称性选择 convention 的风险。

## 6. Independent analytic closure B — free-field / incident plane

理论上，同一距离 `L`：

```text
PE physical phasor  = exp(+i*kL)
Bellhop raw phasor  = exp(-i*kL) = conj(PE physical phasor)
Bellhop absolute    = conj(raw) = exp(+i*kL)
```

既有 `reports/pe_bellhop_incident_field_comparison_report.md` 对相同 4 kHz
incident plane 的结果为：

| comparison | normalized complex L2 | phase RMS (rad) | 备注 |
|---|---:|---:|---|
| PE--AS | `3.7313678e-13` | `2.9570732e-13` | PE `exp(+i R(kz-k0))` |
| Bellhop--AS, 5001 beams | `0.0061022184` | `0.0061020139` | fixed reader/convention |
| Bellhop--AS, 10001 beams | `0.0061050369` | `0.0061048329` | 10,001-beam frozen endpoint |
| PE--Bellhop comparison | `0.0061048459` | `0.0061048329` | no data-dependent transform |

该结果同时通过所有 incident hard gates（包括 receiver coordinate、half-dx
shared nodes、symmetry、NaN/Inf）。PE--AS 的双精度闭合和 Bellhop 的稳定有限
source/beam floor，与上面的 `+kL`/`-kL` 源码推导一致。

## 7. 为什么 Stage 0 flat PASS 不能单独证明 rough convention

Stage 0 的 flat reflected field 只包含共同 carrier、共同 pressure-release phase
和共同 source/receiver chain。其 ratio 对应 `eta=0`：

```text
G_PE(0)=G_BH(0)=1
```

因此 axis normalization 与 flat denominator 会消除 common carrier；即使 rough
phase sign 尚未闭合，flat field 也可以通过 Stage-0 numerical/source floor。
这与 Stage 1X 发现的 weak-sinusoid fixed conjugation 不矛盾。Stage 0 PASS
证明的是 common flat chain 和数值地板，不是 `exp(+i2k eta)` 与
`exp(-i2k eta)` 的 rough convention。

## 8. 正式 comparison convention（永久冻结）

为避免把不同层次的场混在一起，后续报告使用以下名称：

```text
B_raw          = SHD direct complex pressure
B_abs          = conj(B_raw)                         (phase_sign=-1)
G_BH_abs       = B_abs,rough / B_abs,flat
G_BH_comparison= conj(G_BH_abs) = G_BH_raw
```

PE `G_PE` 保持不变。以上两个 conjugation 出现在不同层次：第一个是绝对
Bellhop raw SHD 到项目正向空间相量的固定转换；第二个是为了让 rough/flat
ratio 采用与 PE screen 相同的物理 `+2*k*eta` 符号。它们不是重复反射、不是
pressure-release phase、也不是 amplitude calibration。主程序实现位置是
`validate_pe_bellhop_controlled_comparison.m` lines 286--307；原始 `G_BH`
和 raw Stage-1 FAIL 结果仍保留为诊断。

## 9. Stage 1Y gate checklist

| gate | 结果 | 证据 |
|---|---|---|
| PE time-harmonic convention 明确 | PASS | PE phase-convention audit + helper |
| PE propagation sign 明确 | PASS | split-step/AS `exp(+i*(kz-k0)L)` |
| PE screen `-exp(+i2k eta)` 明确 | PASS | 1-D helper line 30 |
| FFT/IFFT 与 transpose 无隐藏共轭 | PASS | helper/main source audit |
| Bellhop coherent phase sign 明确 | PASS | `influence.f90` `ApplyContribution` |
| Bellhop `Reflect2D` vacuum `+pi` exactly once | PASS | `ReflectMod.f90` lines 95--108 |
| SHD real/imag reader无隐藏共轭 | PASS | official reader + local reader |
| X normalization不携带复相位 | PASS | `ScalePressure` X branch + `.sbp` writer |
| proper half-turn仅坐标映射 | PASS | internal-wall overlay; state unchanged |
| independent free-field closure | PASS | incident-field report |
| independent constant-height closure | PASS | constant-height sign audit |
| per-case fitting / source refit | **NOT USED** | fixed code path and frozen `.sbp` |
| fixed conjugation可由理论推出 | PASS | sections 3--6 |

## 10. Remaining limits

- Bellhop--AS incident residual约 `0.0061`、Stage-0 flat internal residual约
  `0.0059` 是现有 X Gaussian beam/source mapping 的有限数值 floor；它们不构成
  convention attribution blocker，但也不应被描述为 PE 与 Bellhop 的严格
  machine-precision 相同。
- Constant-height Bellhop ratio 仍有约 `2.56e-5 rad` phase 和 `0.00843 dB`
  TL residual，符合有限 beam/source 表示的诊断量级。
- Stage 1X 的 centered/even weak sinusoid 中 mirror+conjugation 与 conjugation
  数值退化；正式结论不依赖该退化，而依赖 free-field 与 `eta=±0.05` 的独立
  解析闭合。
- 本报告只闭合 complex convention，不声称 PE Kirchhoff screen 与 Bellhop
  local-specular model 在任意粗糙度、斜率或曲率下物理等价。Stage 2 仍必须按
  Goal 顺序运行并单独判定 roughness model discrepancy。

## 11. 最终决策

```text
CONVENTION_THEORETICALLY_CLOSED
```

理由：PE phasor、Bellhop coherent influence、vacuum phase、SHD reader、X
normalization 和 proper rotation 均已从源码闭合；固定 conjugation 由定义独立
推出；free-field 与非 sinusoidal constant-height 既有数据支持预测；无需
per-case transformation、经验幅相校正或修改 Bellhop/PE 核心。

因此按 revised Goal：

```text
Stage 1Y gate = PASS
Stage 1 convention-fixed regression = PASS
Stage 2 = UNLOCKED
```

本轮没有执行 Stage 2，也没有启动 PM、强正弦、频率扩展或 Monte Carlo。
