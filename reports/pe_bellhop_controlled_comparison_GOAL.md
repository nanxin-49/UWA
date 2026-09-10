# PE--Bellhop controlled rough-surface comparison — authoritative GOAL

状态：**FROZEN / READY FOR SEQUENTIAL EXECUTION**  
冻结日期：2026-09-10

本文件是后续执行的唯一 authoritative Goal。执行者必须先读
`reports/pe_bellhop_controlled_comparison_goal_audit.md`，再严格按本文件顺序
推进。禁止子 Agent、禁止并行启动多个 Stage、禁止为缩小差异修改核心物理。

## 1. 最终目标

在统一 line source、物理 receiver、surface profile、reflected-field 定义和
数值门槛下，通过

```text
flat -> weak sinusoid -> height sweep -> slope/curvature sweep
     -> validity map -> phase attribution -> fixed PM -> existing M=50 interpretation
```

判断 PE Kirchhoff phase screen 与 Bellhop local-specular `Reflect2D` 是否在
`A -> 0` 时回到共同 flat/numerical floor，并定位 roughness 增强后的分离机制。

## 2. 全程冻结项

| item | frozen value |
|---|---|
| frequency / medium | `4000 Hz`, uniform `c=1500 m/s` |
| validation geometry | source-to-wall `100 m`, wall-to-Rx `3 m`, image range `103 m` |
| PE | 1-transverse helper；`W=192.1875 m`, `nx=984`, `dx=0.1953125 m`, step `0.05 m`, sponge off |
| source | unit-peak Gaussian `sigma=0.3 m` |
| Bellhop source/run | line source `X`, coherent `C` |
| `.sbp` | existing formula, 2401 points, `[-30,30] deg`, `-120 dB` floor |
| Bellhop wall | generic validation-only parametric internal wall; native `Reflect2D` once |
| post-wall map | `r'=200-r`, `z'=-z`; no state/phase/amplitude recalculation |
| receiver map | PE `x` -> Bellhop `(r',z')=(103,-x)`; explicit `[102,103] m` selector |
| wall support / primary N | `s in [-80,80] m`, `N=4097`; density check `N=2049` |
| Bellhop step | `0.05 m` |
| pressure release | coefficient `-1`, exactly one pi phase change |
| output root | `results/validation/pe_bellhop_controlled_comparison/` |
| final report | `reports/pe_bellhop_controlled_comparison_final_report.md` |

禁止修改 `src/**`、PE square-root operator、surface phase-screen formula、
`Reflect2D`、`InfluenceGeoHatCart`、SHD reader/selector、source normalization 或
communication code。仅可新增 `scripts/validation/` orchestration/support 和
对应结果/报告。

## 3. 统一场与指标

Stage 0 使用 axis-normalized flat reflected fields
`F_s(x)=H_s,flat(x)/H_s,flat(0)`。Stage 1 以后使用逐 solver flat ratio：

```text
G_PE(x)=Href_PE,rough(x)/Href_PE,flat(x)
G_BH(x)=Href_BH,rough(x)/Href_BH,flat(x)
```

M95/M99、weights 和 valid-phase mask 只由 Stage-0 exact-AS flat field定义并
冻结。令 M99 weights `w proportional |F_AS|^2` 且和为 1：

```text
E_G = sqrt(sum(w |G_PE-G_BH|^2) / sum(w |G_BH|^2))
S   = sum(w G_PE conj(G_BH))
c   = S / sqrt(sum(w|G_PE|^2) sum(w|G_BH|^2))
phi0 = arg(S)
rho_shape = |c|
rho_raw   = real(c)
E_aligned = min_phi ||G_PE-exp(i phi)G_BH||_w / ||G_BH||_w
          = ||G_PE-exp(i phi0)G_BH||_w / ||G_BH||_w
```

同时报告 M95 `Delta TL_RMS`、circular `Delta phase_RMS`、P95、global phase
`phi0`、`rho_raw`、`rho_shape` 和 `E_aligned`。LS/global alignment 只作归因，
不得替代 raw `E_G` 或参与场校准。

## 4. Preparation P0

只做一次：

1. 验证 official Bellhop 2020 executable；从当前源码构建 validation-only
   flat 和 generic PM internal-wall binaries。
2. 建立新 output tree；不得导入旧 R/X case caches。
3. 写一个新的 controlled-comparison entrypoint，必须可按单一 Stage 运行和
   resume，但一次调用只能推进一个已解锁 Stage。
4. 保存 executable/source hashes、`.sbp` fingerprint、receiver permutation、
   requested/actual receiver coordinates和完整 case fingerprint。
5. 对新增 MATLAB 文件执行静态检查。MATLAB 执行使用现有 MATLAB MCP session；
   session 不可用时停止并报告，不得从 shell 启动 MATLAB。

P0 不产生物理结论。

## 5. Stage 0 — flat receiver-line closure

依次执行：

1. PE flat reflected receiver line与独立一步 exact AS at `103 m`。
2. generic parametric wall输入 `r(s)=100 m`，与专用 flat binary及 official
   X-source free field `-P_X(103 m)` 做 Bellhop-chain closure。
3. Bellhop 5001/10001 beam pair。
4. 10001-beam half-dx receiver grid；只比较与 PE grid 共享的节点。
5. PE--Bellhop axis-normalized complex receiver-line comparison。

Hard gates 完全采用审核报告第 5 节。另要求：全部 rays hit一次、wall residual
`<=1e-9 m`、phase jump error `<=1e-10 rad`、p/q rotation error `<=1e-12`、
`min post dr>0`、NaN/Inf=0。

结果判断：

- `PASS`：记录 `F_E/F_phi/F_TL` 和由其计算的 `T_E/T_phi/T_TL`，解锁 Stage 1。
- `FAIL`：立即停止，最终类别 `BLOCKED_BY_COMPARABILITY`；不得运行 rough case。

## 6. Stage 1 — weak sinusoid asymptotic test

统一 profile：

```text
eta(s)=A cos(Ks)
K=0.10 rad/m
Gamma(s)=[100-eta(s),s]
```

先运行 `A=0.01 m`。该 case 必须先完成：

- Bellhop `10001 -> 20001` beam convergence；
- `N=2049 -> 4097` profile convergence；
- 每个 case 的 intersection/frame/kappa/direction/RN/RM/p/q/tau/pi-phase/
  positive-range/receiver-coordinate guards。

beam/profile convergence沿用 Stage-0 limits。若 10001->20001 beams 不通过，
只增加 40001 beams；若 20001->40001 仍不通过，停止为
`BLOCKED_BY_COMPARABILITY`。收敛 endpoint 成为后续全部新 case 的固定 beam
count。

若 `A=0.01` 满足：

```text
E_G <= T_E
phase_RMS <= T_phi
TL_RMS <= T_TL
rho_shape >= 0.9995
```

则 weak limit 在当前 resolution 下通过，直接解锁 Stage 2。若未满足，按顺序
运行 `A=0.005`、再运行 `A=0.0025 m`。最终通过还要求随 A 降低的 `E_G` 不得
超过上一点加 `max(F_E,0.002)`，且最小 A 满足上述 floor gates。

若 `A=0.0025 m` 仍未回到 floor，而 flat/source/receiver/solver convergence
全部通过：立即停止 Stage 2--7，触发第 12 节 Helmholtz reference branch。

## 7. Stage 2 — height-phase sweep

仅 Stage 1 PASS 后执行。固定 `K=0.10 rad/m`，按序运行并复用已有点：

```text
A = 0.01, 0.02, 0.05, 0.10, 0.20 m
```

每个点报告 `2kA`、`E_G`、TL/phase RMS、`phi0`、`rho_shape`、
`E_aligned` 和全部数值/几何 guards。任何 solver convergence 或 mapping guard
失败都停止；模型差异本身不触发停止。

## 8. Stage 3 — slope/curvature sweep

仅 Stage 2 完成后执行。固定 `A=0.20 m`，复用 `K=0.10` 并依次运行：

```text
K = 0.10, 0.20, 0.35, 0.47 rad/m
```

记录解析与 sampled：

```text
max slope = A K
kappa(s)=eta''(s)/(1+eta'(s)^2)^(3/2)
max |kappa|, RMS |kappa|, minimum radius
```

同时记录 Bellhop hit-point curvature 分布和 common field metrics。不得把
curvature、p/q 或 phase screen人为改成一致。

## 9. Stage 4 — validity map

仅使用通过 numerical/mapping guards 的 Stage 1--3 cases：

- **Region I**：`E_G<=T_E`、`E_aligned<=T_E`、TL/phase RMS 分别不超过
  `T_TL/T_phi`、`rho_shape>=0.9995`。
- **Region II**：不满足 I，但 `E_aligned<=0.10`、TL RMS `<=0.5 dB`、
  phase RMS `<=0.20 rad`、`rho_shape>=0.99`。
- **Region III**：其余已数值收敛 case。

分别对固定 K 的 height axis 和固定 A 的 slope/curvature axis定位首次离开
Region I 的区间；只能报告采样区间，不外推精确 boundary。

## 10. Stage 5 — phase discrepancy attribution

对 Region II/III 代表点执行 validation-only diagnostics：

- global-phase dominated：`rho_shape>=0.995` 且 phase alignment 使 residual
  至少降低 50%；
- spatial distortion：`rho_shape<0.995` 或 alignment 后 residual 降低不足 50%。

允许比较 stationary-phase/local-specular prediction 和 `2k eta` vs
angle-aware `2k_z eta`，但只能输出 diagnostics。不得修改 production PE
screen 或 Bellhop reflection。

## 11. Stage 6--7 — fixed PM and existing M=50

### Stage 6 fixed PM

只运行 seed `260001`、4 kHz。到达本阶段后才复制并核对 canonical coefficients
SHA-256：

```text
1f7eda465e4ae85b8ac038310edf053f2d062563a3b943bfd83bd59d07015f67
```

使用同一 receiver line、X source、flat denominator、masks和 primary beam
count。另做一次 10001/20001（必要时 40001）PM receiver-line convergence。
报告全部 common metrics，并按 RMS/max height、slope、curvature和 `2k sigma`
放入 Stage-4 validity map。原 X-source axis结果
`0.310371 dB / -2.248201 rad` 仅作 sanity check。

### Stage 7 existing M=50

禁止重新运行任何 seed。只读取 authoritative M=50 CSV/MAT并记录 hashes。
必须明确标注这些样本为 historical `R` provenance；不得改名为 X 或做经验
R->X correction。结合 Stage 5--6 解释：mean powers近似一致但 circular phase
mean `-0.37534 rad` 的可能机制，并单列 source-provenance limitation。

## 12. Conditional Helmholtz reference branch

仅在 Stage 1 最小 A 未回到 floor 时触发。停止所有 stronger sinusoid/PM工作，
先做 validation-only BIE/BEM feasibility gate，再只求：flat、weakest sinusoid、
一个 stronger sinusoid。必须使用同一 2-D line-source angular spectrum、Dirichlet
pressure-release boundary和receiver line，并分别报告 BEM--PE、BEM--Bellhop。

如果无法建立经网格/边界截断收敛的参考解，最终状态为
`BLOCKED_BY_COMPARABILITY`，不得猜测误差归属。

## 13. 最终状态与报告

最终只允许：

- `COMMON_LIMIT_CONFIRMED`
- `PHASE_MECHANISM_IDENTIFIED`
- `REFERENCE_ADJUDICATED`
- `BLOCKED_BY_COMPARABILITY`

生成 `reports/pe_bellhop_controlled_comparison_final_report.md`，回答用户列出的
十个问题，并附：case/config/hash manifest、Stage gate table、receiver mapping
audit、flat floor、每个 sinusoid/PM 的 raw/aligned metrics和现有 M=50 provenance。

每完成一个 Stage，先写结果和 PASS/PASS_WITH_LIMITS/FAIL，再决定是否进入下一
Stage；禁止预先并行运行后续 Stage。
