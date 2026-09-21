# PE--Bellhop controlled rough-surface comparison — authoritative execution record

状态：**COMPLETED / AUTHORITATIVE EXECUTION RECORD**
冻结日期：2026-09-10  
修订说明：已纳入 Stage 1X phase-convention attribution、Stage 1 convention-only regression、Stage 1Y theoretical convention closure，并统一 10,001-beam 执行覆盖。  
执行记录：按 Stage 顺序、单一验证会话和固定比较约定完成；禁止并行运行和结果拟合。

本文件是当前权威执行记录；Stage 0--7 已全部完成。历史执行前置说明如下：

```text
reports/pe_bellhop_controlled_comparison_goal_audit.md
```

再按本文件顺序推进。

历史执行原则：

- 不并行启动多个 Stage；
- 不得为缩小 PE/Bellhop 差异修改核心物理；
- 不得根据结果临时调整 comparison convention；
- 不得通过 amplitude/phase fitting 人为闭合结果。

---

## 0. 当前已完成状态

### Stage 0

flat receiver-line closure 已 PASS。

Stage-0 已冻结：

```text
T_E   = 0.0108994
T_phi = 0.0108992 rad
T_TL  = 0.0107979 dB
```

后续不得根据 rough-surface 结果修改这些阈值。

### 原始 Stage 1

profile：

```text
eta(s)=A cos(0.1 s)
```

原始 raw comparison 在：

```text
A = 0.0100
A = 0.0050
A = 0.0025 m
```

均未达到 weak-limit complex/phase floor。

该原始结果永久保留：

```text
Stage 1 RAW = FAIL
```

不得删除或改写。

### Stage 1X

Stage 1X 已完成：

```text
CONVENTION_DISCREPANCY_IDENTIFIED
```

固定 conjugation relation 能同时闭合三个 A case。

但 Stage 1X 只证明：

> 固定 conjugation 能解释原始 weak-sinusoid phase discrepancy。

它不单独证明：

> 该 conjugation 在数学上一定是正确的正式 comparison convention。

### Convention-only regression

comparison-only fixed convention 后，Stage 0 重新 PASS。

修正后的 Stage 1：

| A (m) | corrected E_G | corrected phase RMS | corrected TL RMS | rho_shape |
|---:|---:|---:|---:|---:|
| 0.0100 | 0.00409346 | 0.00405930 rad | 0.00456998 dB | 0.99999242 |
| 0.0050 | 0.00204123 | 0.00202407 rad | 0.00228491 dB | 0.99999810 |
| 0.0025 | 0.00101951 | 0.00101096 rad | 0.00114248 dB | 0.99999952 |

因此：

```text
Stage 1 convention-fixed regression = PASS
```

但 **Stage 2 尚不能解锁**。

必须先完成 Stage 1Y。

---

## 1. 最终目标

在统一：

- 2-D Cartesian line source；
- physical receiver；
- surface profile；
- reflected-field definition；
- fixed comparison convention；
- numerical floor；

条件下，通过：

```text
flat
  ->
weak sinusoid
  ->
height sweep
  ->
slope/curvature sweep
  ->
validity map
  ->
phase attribution
  ->
fixed PM
  ->
existing M=50 interpretation
```

回答：

> PE Kirchhoff phase screen 与 Bellhop local-specular `Reflect2D`
> 是否在共同 weak-roughness limit 中趋向同一个 reflected complex field，
> 并确定 roughness 增强后两种 Helmholtz approximation 从哪里开始分离。

若 weak limit 在 convention 理论闭环后仍无法成立，则进入第三方 Helmholtz
BIE/BEM reference adjudication。

---

## 2. 历史执行原则

历史执行采用单一验证会话，按已解锁 Stage 顺序推进，并优先复用
authoritative outputs。不得并行启动多个 Stage，不得通过临时的模型角色、
代理分工、结果拟合或重设比较约定来改变实验结论。

允许使用普通工具、MATLAB、Python、shell、Bellhop executable 和文件搜索。

---

## 3. 全程冻结项

| item | frozen value |
|---|---|
| frequency / medium | `4000 Hz`, uniform `c=1500 m/s` |
| geometry | source-to-wall `100 m`, wall-to-Rx `3 m`, image range `103 m` |
| PE | 1-transverse helper；`W=192.1875 m`, `nx=984`, `dx=0.1953125 m`, step `0.05 m`, sponge off |
| source | unit-peak Gaussian `sigma=0.3 m` |
| Bellhop source/run | Cartesian line source `X`, coherent `C` |
| `.sbp` | existing formula, 2401 points, `[-30,30] deg`, `-120 dB` floor |
| Bellhop wall | generic validation-only parametric internal wall |
| reflection | native `Reflect2D` exactly once |
| post-wall map | `r'=200-r`, `z'=-z` |
| receiver map | PE `x` -> Bellhop `(r',z')=(103,-x)` |
| wall support | `s in [-80,80] m` |
| primary profile N | `4097` |
| density check | `2049 -> 4097` |
| Bellhop step | `0.05 m` |
| Bellhop beams | **固定 10,001 beams** |
| pressure release | `-1`, exactly one pi phase |
| output root | `results/validation/pe_bellhop_controlled_comparison/` |
| final report | `reports/pe_bellhop_controlled_comparison_final_report.md` |

禁止修改：

```text
src/**
PE square-root operator
surface phase-screen formula
Reflect2D
InfluenceGeoHatCart
SHD reader/selector
source normalization
communication code
```

只允许新增 validation-only orchestration / support / diagnostics。

---

## 4. Bellhop beam-count 执行覆盖

所有后续新增 Bellhop case 固定使用：

```text
10,001 beams
```

因此：

- 不再运行 20,001；
- 不再运行 40,001；
- Goal 中任何历史 `10001 -> 20001` / `20001 -> 40001` 要求永久失效；
- 已存在的历史高 beam-count 结果仅作 diagnostic；
- 不得作为新的 hard gate。

如果 10,001-beam 下出现：

- 非有限场；
- receiver extraction 不稳定；
- geometry guard fail；
- 与 Stage-0 frozen floor 明显冲突；

则判：

```text
BLOCKED_BY_COMPARABILITY
```

不得通过增加 beams 绕过用户约束。

---

## 5. 统一场与指标

Stage 0 使用：

```text
F_s(x)=H_s,flat(x)/H_s,flat(0)
```

Stage 1 以后使用：

```text
G_PE(x)=Href_PE,rough(x)/Href_PE,flat(x)
G_BH(x)=Href_BH,rough(x)/Href_BH,flat(x)
```

Stage 1Y 通过后，将永久冻结一个唯一正式 comparison convention。

在 Stage 1Y 完成前，不允许把 `conj()` 视为最终理论上已证明的 convention。

M95/M99、weights 与 valid-phase mask 只由 Stage-0 exact-AS flat field 定义。

M99：

```text
w proportional |F_AS|^2
sum(w)=1
```

phase 只在：

- M95；
- PE/BH flat denominator 均高于 axis `-40 dB`；

的位置统计。

rough result 不允许参与 mask 选择。

### Primary metrics

```text
E_G
=
sqrt(
    sum(w |G_PE-G_BH|^2)
    /
    sum(w |G_BH|^2)
)
```

```text
S = sum(w G_PE conj(G_BH))

c =
S /
sqrt(
    sum(w |G_PE|^2)
    sum(w |G_BH|^2)
)
```

报告：

```text
rho_shape = |c|
rho_raw   = real(c)
phi0      = arg(S)
```

phase-aligned diagnostic：

```text
E_aligned
=
min_phi
||G_PE-exp(i phi)G_BH||_w
/
||G_BH||_w
```

LS/global alignment 只能用于归因，不得替代 raw `E_G`。

同时报告：

```text
Delta TL_RMS
Delta TL_P95
Delta phase_RMS
Delta phase_P95
```

deep-null phase 不作为 hard gate。

---

## 6. Stage 0 — flat receiver-line closure

Stage 0 已完成并 PASS。

后续仅当 Stage 1Y 定位到 comparison/extraction convention bug 会影响 flat comparison 时，
才允许重新运行受影响的 Stage-0 子项。

Stage-0 frozen hard gates保持不变。

---

## 7. Stage 1 — weak sinusoid asymptotic test

profile：

```text
eta(s)=A cos(Ks)
K=0.10 rad/m
Gamma(s)=[100-eta(s),s]
```

原始 raw Stage 1 已完成并 FAIL。

Stage 1X 已发现 fixed conjugation discrepancy。

convention-only regression 已重新获得：

```text
A=0.01   PASS
A=0.005  PASS
A=0.0025 PASS
```

因此 Stage 1 当前状态：

```text
RAW comparison: FAIL
convention-fixed regression: PASS
```

最终是否可以解释为：

```text
COMMON WEAK LIMIT CONFIRMED
```

取决于 Stage 1Y。

---

## 8. Stage 1X — weak-sinusoid phase-convention attribution

Stage 1X 已完成。

结果：

```text
CONVENTION_DISCREPANCY_IDENTIFIED
```

固定候选比较：

```text
G_BH(x) vs G_PE(x)
G_BH(x) vs conj(G_PE(x))
G_BH(x) vs G_PE(-x)
G_BH(x) vs conj(G_PE(-x))
```

结果表明：

- raw PE 不闭合；
- fixed conjugation 同时闭合三个 A；
- mirror 本身不起作用；
- mirror+conjugation 在 centered/even sinusoid 下与 conjugation 数值退化。

weak-response coefficient：

```text
phase_RMS/A ≈ 47.92 rad/m
```

与：

```text
4k/sqrt(2) ≈ 47.39 rad/m
```

高度接近。

Stage 1X 只负责：

```text
发现
```

不负责：

```text
证明
```

因此必须进入 Stage 1Y。

---

## 9. Stage 1Y — theoretical convention closure

### 9.1 目的

这是 Stage 2 前最后一个 hard gate。

唯一目标：

> 从 PE 与 Bellhop 各自的 time-harmonic / propagation / file-output
> complex convention 独立证明，正式比较时为什么必须或不必须进行固定 complex conjugation。

禁止通过：

```text
哪个 transformation 误差最小
```

来决定 convention。

### 9.2 PE convention audit

必须从代码和数学定义追踪：

1. 时间谐波 convention：
   - `exp(-i omega t)` 或
   - `exp(+i omega t)`；
2. PE forward propagation kernel 的符号；
3. plane-wave `H_PE(L)` 应为 `exp(+ikL)` 还是 `exp(-ikL)`；
4. pressure-release `-1` 的 phase；
5. rough phase screen `-exp(i 2 k eta)` 的物理 phase sign；
6. FFT/IFFT convention；
7. MATLAB transpose 是否使用 `'` 或 `.'`；
8. 保存/读取 PE complex field 时是否发生 conjugation。

必须给出：

```text
PE expected complex propagation sign
```

### 9.3 Bellhop convention audit

必须从 Bellhop 2020 source 与当前 MATLAB reader 路径追踪：

1. Bellhop time-harmonic convention；
2. ray / Gaussian beam phase accumulation sign；
3. SHD real/imag pressure convention；
4. SHD reader 是否 direct-read / conjugate / transpose-conjugate；
5. 当前历史 `phase_sign=-1` 的来源；
6. source `X` normalization 是否包含 fixed conjugation；
7. proper half-turn 是否只改变坐标而不改变 complex-field convention。

必须给出：

```text
Bellhop expected complex propagation sign
```

### 9.4 Independent analytic closure

不得使用 weak sinusoid 作为唯一证据。

必须先写理论预测，再读取已有结果验证。

至少完成一个，优先两个都完成。

#### A. Free-field propagation

同一传播距离 `L`，预测：

```text
H_PE(L)
H_BH(L)
```

的 complex phase sign。

判断它们是否天然互为 complex conjugate。

#### B. Constant-height reflection

使用已有：

```text
eta=+h0
eta=-h0
```

constant-height audit。

先预测：

```text
G_PE(+h0)
G_PE(-h0)
G_BH(+h0)
G_BH(-h0)
```

的 phase sign，再读取保存结果。

禁止从数据反推 prediction。

### 9.5 Flat PASS 为什么不能单独证明 convention

必须解释：

- `eta=0` 时 roughness-induced phase 为 0；
- axis normalization / flat denominator 可消除 common carrier；

因此：

```text
Stage 0 flat PASS
```

与：

```text
rough case 发现 fixed conjugation
```

不矛盾。

### 9.6 决策

#### CONVENTION_THEORETICALLY_CLOSED

只有全部满足：

- PE phasor convention 明确；
- Bellhop/SHD convention 明确；
- fixed conjugation 可由数学定义独立推出；
- 至少一个非 sinusoid analytic/existing case 支持预测；
- 不需要 per-case transformation；

才允许。

此时永久冻结：

```text
official comparison convention
```

例如：

```text
G_BH_comparison = conj(G_BH_raw)
```

或源码定义等价形式。

随后：

```text
Stage 1 = PASS
Stage 2 unlocked
```

#### CONVENTION_NOT_PROVEN

若 conjugation 只能由“误差更小”支持，而不能从理论/源码独立推出：

```text
Stage 2 remains locked
```

继续最小 convention/comparability investigation。

若无法进一步解决，则考虑 Conditional Helmholtz reference。

#### BLOCKED_BY_ATTRIBUTION

若缺少必要保存数据，但可通过最小补充诊断解决：

- 只允许最小新增 diagnostic；
- 固定 10,001 beams；
- 禁止 stronger sinusoid / PM / frequency / Monte Carlo。

### 9.7 输出

生成：

```text
reports/pe_bellhop_controlled_comparison_stage1y_theoretical_convention_closure.md
```

必须回答：

1. PE time-harmonic convention；
2. PE propagation sign；
3. Bellhop time-harmonic convention；
4. Bellhop SHD complex convention；
5. reader 是否共轭；
6. 为什么需要或不需要固定 conjugation；
7. independent analytic check；
8. 最终 decision。

---

## 10. Stage 2 — height-phase sweep

只有以下全部成立才允许执行：

```text
Stage 1 convention-fixed regression = PASS
Stage 1Y = CONVENTION_THEORETICALLY_CLOSED
```

固定：

```text
K = 0.10 rad/m
```

按序运行：

```text
A = 0.01
    0.02
    0.05
    0.10
    0.20 m
```

已存在 `A=0.01` 可复用，只要 fingerprint 完全一致。

每个点报告：

```text
2kA
max slope
max/RMS curvature
E_G
TL_RMS
phase_RMS
phi0
rho_raw
rho_shape
E_aligned
```

以及全部 numerical / geometry guards。

model discrepancy 本身不触发停止。

solver / mapping / field finiteness 失败才停止。

---

## 11. Stage 3 — slope/curvature sweep

仅 Stage 2 完成后执行。

固定 A 按以下规则一次性冻结：

1. 优先选择 Stage 2 中仍属于 Region I 的最大 A；
2. 如果所有 `A>=0.05 m` 均已离开 Region I，则固定：
   ```text
   A=0.05 m
   ```
3. A 一旦写入 Stage-3 manifest，不得根据 K-sweep 结果更改。

然后扫描：

```text
K = 0.10
    0.20
    0.35
    0.47 rad/m
```

记录：

```text
max slope = A*K

kappa(s)
=
eta''(s)
/
(1+eta'(s)^2)^(3/2)

max |kappa|
RMS |kappa|
minimum radius
```

同时保存 Bellhop：

- hit-point curvature；
- incidence；
- RN/RM；
- p/q；
- tau；
- wall residual。

不得把 PE/Bellhop curvature physics 人为改成一致。

---

## 12. Stage 4 — validity map

仅使用通过 numerical/mapping guards 的 Stage 1–3 case。

### Region I

```text
E_G <= T_E
E_aligned <= T_E
TL_RMS <= T_TL
phase_RMS <= T_phi
rho_shape >= 0.9995
```

### Region II

不满足 I，但：

```text
E_aligned <= 0.10
TL_RMS <= 0.5 dB
phase_RMS <= 0.20 rad
rho_shape >= 0.99
```

### Region III

其余已经数值收敛的 case。

只报告 sampled validity interval，不外推精确 boundary。

---

## 13. Stage 5 — phase discrepancy attribution

对 Region II / III 代表 case 执行。

### Global-phase dominated

若：

```text
rho_shape >= 0.995
E_aligned <= 0.5 * E_G
```

优先归类为 global/coherent phase dominated。

### Spatial distortion

若：

```text
rho_shape < 0.995
```

或 alignment 后 residual 降低不足 50%，则归类为 spatial distortion。

允许 validation-only diagnostic：

- stationary phase；
- local specular direction；
- phase gradient；
- `2k eta`；
- angle-aware `2k_z eta`。

禁止修改 production PE screen。

---

## 14. Stage 6 — fixed PM realization

只运行：

```text
seed = 260001
f = 4 kHz
```

到达本 Stage 后才复制/核对 canonical Fourier coefficients。

SHA-256：

```text
1f7eda465e4ae85b8ac038310edf053f2d062563a3b943bfd83bd59d07015f67
```

固定：

```text
X source
10001 beams
same receiver line
same flat denominator
same masks
same official comparison convention
```

不得启动新的 20,001/40,001-beam case。

报告：

```text
E_G
TL_RMS
phase_RMS
phi0
rho_raw
rho_shape
E_aligned
```

以及：

- RMS/max height；
- RMS/max slope；
- RMS/max curvature；
- `2k sigma_eta`。

然后将 fixed PM 放入 Stage-4 validity map。

原历史 axis：

```text
0.310371 dB / -2.248201 rad
```

只能作 sanity check。

---

## 15. Stage 7 — existing M=50 interpretation

禁止重跑任何 seed。

只读取 authoritative M=50 CSV/MAT。

必须明确标注：

```text
historical R-source provenance
```

不得：

- 改名为 X；
- 做经验 R->X correction；
- 用 M=50 覆盖新的 X-source controlled comparison。

结合 Stage 4–6 解释：

```text
mean PE/BH reflected power ~ equal
circular mean phase != 0
```

是否能够由：

- convention；
- global phase；
- spatial distortion；
- roughness validity region；

解释。

---

## 16. Conditional Helmholtz reference branch

### 16.1 Stage 1X 后

如果：

```text
Stage 1X = NO_CONVENTION_DISCREPANCY
```

且 weak-limit discrepancy 仍远离 floor：

直接进入本 branch。

### 16.2 Stage 1Y 后

如果：

```text
Stage 1X = CONVENTION_DISCREPANCY_IDENTIFIED
```

但：

```text
Stage 1Y = CONVENTION_NOT_PROVEN
```

且经过最小 comparability investigation 仍无法从 phasor/source/reader 定义证明 fixed conjugation：

允许进入本 branch。

### 16.3 Reference 最小任务

只求：

1. flat；
2. weakest sinusoid；
3. one stronger sinusoid。

建立：

```text
2-D homogeneous Helmholtz
line source
Dirichlet pressure-release rough boundary
BIE/BEM reference
```

必须先验证：

- flat analytic/image solution；
- discretization convergence；
- truncation/domain treatment；
- boundary residual。

然后比较：

```text
BEM <-> PE
BEM <-> Bellhop
```

不得直接扩展到 PM ensemble。

---

## 17. 最终状态

最终只允许：

```text
COMMON_LIMIT_CONFIRMED
PHASE_MECHANISM_IDENTIFIED
REFERENCE_ADJUDICATED
BLOCKED_BY_COMPARABILITY
```

---

## 18. COMMON_LIMIT_CONFIRMED 的要求

至少：

1. Stage 0 flat receiver-line closure PASS；
2. Stage 1Y 理论 convention 闭环；
3. weak sinusoid fixed-convention regression PASS；
4. A 减小时 `E_G/phase_RMS` 向 flat floor 下降；
5. `rho_shape` 接近 1；
6. numerical/geometry guards 全部通过。

这时允许结论：

> PE 与 Bellhop 在当前共同 weak-roughness limit 中具有一致的主导 reflected complex-field behavior。

---

## 19. Final report

生成：

```text
reports/pe_bellhop_controlled_comparison_final_report.md
```

必须回答：

1. flat spatial comparison floor 是多少？
2. PE/Bellhop comparison 为什么需要或不需要 fixed conjugation？
3. 该 convention 是否由理论独立证明？
4. A -> 0 时是否趋于共同 floor？
5. weak sinusoid 的 E_G / TL / phase / rho 是多少？
6. discrepancy 首先随 height、slope 还是 curvature 增长？
7. phase discrepancy 是 global offset 还是 spatial distortion？
8. PE/Bellhop 从什么 roughness interval 开始明显分离？
9. fixed PM 落在哪个 validity region？
10. M=50 的 mean-power agreement 与 phase bias 是否得到解释？
11. 是否需要 BIE/BEM reference？
12. 当前证据能否支持目标 PM 工况下 PE rough-surface model 的有效性？

---

## 20. Stage execution discipline

每个 Stage：

```text
run one unlocked stage
    ->
write result
    ->
PASS / PASS_WITH_LIMITS / FAIL
    ->
evaluate next gate
```

禁止：

- 提前跑后续 Stage；
- Stage fail 后继续堆更复杂 case；
- 用后续结果反向修改前面 frozen threshold；
- 删除旧 FAIL 结果；
- 覆盖原始 raw solver data。

---

## 21. 推荐报告

```text
reports/
    pe_bellhop_controlled_comparison_stage1x_phase_convention_audit.md
    pe_bellhop_controlled_comparison_stage1_convention_fix_report.md
    pe_bellhop_controlled_comparison_stage1y_theoretical_convention_closure.md
    pe_bellhop_controlled_comparison_final_report.md
```

---

## 22. 一句话总目标

> **先从数学和源码层面闭合 PE/Bellhop 的固定 complex convention，再通过 weak sinusoid 证明两种方法在共同弱粗糙极限中趋向同一 reflected complex field；随后逐步增加 height、slope 和 curvature 建立 validity boundary，并最终解释当前 fixed-PM 与历史 M=50 结果。**
