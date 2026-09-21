# PE--Bellhop PM 粗糙海面对比：完整执行汇总

更新日期：2026-09-09
范围：固定 band-limited 1-D PM realization、1-transverse PE bridge、生产二维 PE sensitivity、Bellhop validation-only internal wall；不包含 PE/Bellhop 核心物理重写、通信链、Monte Carlo 大规模统计或 PE 与 Bellhop 的幅相拟合。

## 总结结论

阶段 0A--0E、Stage 1A--1C、Stage 2、Stage 3 和 Stage 4 均已按顺序独立运行并生成报告。当前结论为：

> **PE propagation 与 Bellhop internal-wall 的数值/几何链分别通过；在同一固定 PM 频带 realization 和已量化的维度差异下，PE Kirchhoff phase-screen 与 Bellhop local-specular Reflect2D 的 reflected-only 差异稳定存在，应归类为 reflection-model discrepancy，而不是数值失败。**

Stage 3A 的状态为 **PASS_WITH_LIMITS**：它是 8 个 paired fixed-band
phase-ensemble seed 的统计 smoke。随后补充的 Stage 3B 使用同一 PM
谱密度频带、固定 seed 和独立 Gaussian Fourier coefficient amplitude
抽样，8 个 seed 也通过了全部结构检查；Stage 3B 仍是低成本 first ensemble，
不是最终收敛的 ocean Monte Carlo。

Stage 4 在相同冻结链上完成了 50 个 fixed-4 kHz seed 的模型差异统计；
24→32 的均值/功率/圆相位门通过，但 bootstrap 均值 delta-TL 半宽为
`0.25835 dB`，略高于预设 `0.25 dB`，因此状态为
**PRELIMINARY_MODEL_DISCREPANCY**。该阶段由用户在 seed 260050 后停止，
不包含 260051 及以后样本。

2026-09-09 的 source-geometry 审计又隔离了 Bellhop `R`/`X` 定义：10001-beam
弱正弦墙 native/internal 的 `-0.26065 dB` 偏差来自 point-source `R` 的
97/103 m range normalization，line-source `X` 将其降至 `4.24e-6 dB`。但仅重跑的
4 kHz seed-260001 X-source Tier-1 仍为 `0.310371 dB / -2.248201 rad`，与原 R
结果 `0.310434 dB / -2.248201 rad` 几乎相同；因此 source artifact 已被定位，
reflection-model discrepancy 结论不变。详见
`reports/bellhop_source_geometry_rx_audit_report.md`。

同日已把 `X` 正式写入 active incident/Tier-1/frequency/ensemble/Stage-4
entrypoints，并把 source geometry 纳入 ensemble cache fingerprint。按原 4 kHz
incident-plane 配置重跑后，幅度 P95 从 pre-X 的 `0.15794309 dB` 降至
`0.0019036364 dB`，全部 13 个 gate 通过。既有 Stage-2/3/4 数值报告仍保留
pre-migration R provenance，没有重写或假标为 X；未来重跑会使用独立 X-tagged
case。详见 `reports/pe_bellhop_line_source_integration_report.md`。

## 冻结输入与不变项

| 项目 | 值 |
|---|---|
| medium | uniform `c=1500 m/s` |
| geometry | `z_tx=100 m`, `z_rx=3 m`, direct/image `97/103 m` |
| source | Gaussian `sigma=0.3 m`, existing `.sbp`; active comparator uses Bellhop line-source `X` (existing Stage-2/3/4 reports retain historical `R` provenance) |
| PM realization | seed `260001`, `U=6 m/s`, span `160 m`, master `N=4097` |
| PM spectrum | requested `Kmax=0.5 rad/m`, realized `0.471238898 rad/m` |
| canonical coefficients SHA-256 | `1f7eda465e4ae85b8ac038310edf053f2d062563a3b943bfd83bd59d07015f67` |
| PE Tier 1 | 1-transverse, `W=192.1875 m`, `nx=984`, `step=0.05 m`, sponge off |
| PE Tier 2 | y-invariant `eta(x,y)=eta(x)`, `ny=256/512` sensitivity |
| Bellhop | validation-only internal wall, native `Reflect2D`, 5001/10001 beam checks, `step=0.05 m`, explicit SHD selector |
| Stage 2 sector | predeclared `[-15°,15°]`, 5001 beams; unchanged Gaussian directivity, no fit |

禁止项全部保持冻结：PE square-root marching、FFT convention、Gaussian source、Kirchhoff screen、Bellhop `Reflect2D`/`InfluenceGeoHatCart`、RN/RM 与 p/q 公式、pressure-release phase、proper half-turn、通信代码均未作为本比较的核心物理修改对象。

## 阶段结果

| 阶段 | 结果 | 关键证据 |
|---|---|---|
| 0A canonical mapper | **PASS** | PE/Bellhop 由同一 Fourier 系数求值；profile provenance/hash、符号和 datum 一致 |
| 0B 1-transverse bridge | **PASS** | 1-D PE 与独立 exact angular spectrum、生产 PE `ky=0` bridge 达到双精度量级 |
| 0C flat/source audit | **PASS** | flat path `103 m`；压力释放相位；5001/10001 beams；轴上 TL offset `0.2610 dB` 为稳定 2-D/3-D diagnostic |
| 0D constant-height sign | **PASS** | `eta=±0.05 m` 的 PE phase `±2k eta`、Bellhop path-time 和单次 `pi` phase 均通过 |
| 0E numerical budget | **PASS** | PE window/grid/step、Bellhop profile/beam/step 扫描均在冻结预算内；wall residual约 `8e-15 m` |
| 1A Tier 1 | **PASS_WITH_LIMITS** | `G_PE=-0.086869+0.988069i`, `G_BH=-0.690339-0.662864i`; delta TL `0.3104 dB`, phase `-2.2482 rad` |
| 1B dimensionality | **PASS_WITH_LIMITS** | 1T→2T delta约 `-0.1589 dB / +0.01251 rad`, `ny=256→512` 仅 `4.36e-6 dB / 4.49e-6 rad` |
| 1C interpretation | **PASS_WITH_MODEL_DISCREPANCY** | PE/Bellhop residual 明显大于 dimensionality sensitivity，归类为反射模型差异 |
| 2 frequency extension | **PASS_WITH_MODEL_DISCREPANCY** | 4/6/8 kHz delta TL `0.3104/0.3083/0.3048 dB`；均有限、零 grazing；不由三点估计 group delay |
| 3A paired phase ensemble | **PASS_WITH_LIMITS** | seeds `260001:260008`；mean delta TL `0.6195 dB`，circular phase `-0.5496 rad`；几何/phase/pq/tau checks 全通过 |
| 3B coefficient-amplitude ensemble | **PASS_WITH_LIMITS** | seeds `260001:260008`、统一 5001 beams；mean delta TL `0.3207 dB`，circular phase `-0.6893 rad`；5001/5001 wall hits，几何/phase/pq/tau checks 全通过 |
| 4 model-discrepancy statistics | **PRELIMINARY_MODEL_DISCREPANCY** | seeds `260001:260050`；mean delta TL `-0.00748 dB`，circular phase `-0.37534 rad`；全部数值/几何护栏通过，bootstrap 半宽 `0.25835 dB` |

Stage 2 的权威 Bellhop case roots 使用 `f4000k_B5001_*`、`f6000k_B5001_*` 和
`f8000k_B5001_*`。早期未采用预声明 ±15° sector 的不完整试跑目录只保留作
工作痕迹，不属于当前结果，不应复用。

## Stage 3 统计范围

每个非 canonical seed 保持 canonical 每模态谱幅度和 `Kmax`，只使用确定性 seed phase 形成 paired realization；没有 Hs 重新归一化、重心化、平滑、taper、带宽改变或 source fit。8 个 seed 的 profile RMS height 约 `0.1864 m`，RMS slope 约 `0.0540`，RMS curvature 约 `0.0209 1/m`，最小曲率半径约 `18.8--23.6 m`；Bellhop `min |u·n|` 为 `0.948--0.962`，grazing fraction 为零，post-wall range 始终为正。

Stage 3A/3B 报告同时给出 PE/Bellhop `20log10|G|`、model delta TL、相位
（含 circular mean/std 及 p5/p25/median/p75/p95）和 `|G|^2`。Bellhop
backward-range Cartesian coherent amplitude 的约 `0.26 dB` 基线偏移仍只作
diagnostic，没有 renormalization。3B 的权威结果位于
`results/validation/pe_bellhop_pm_amplitude_ensemble/`，报告为
`reports/pe_bellhop_pm_amplitude_ensemble_report.md`。

Stage 3B 的 flat baseline 与所有 rough cases 现已统一使用 5001 beams。flat 与
canonical rough 只有在 Stage 2 内嵌配置指纹完全匹配且产物晚于 executable/profile
输入时才复用；其余缓存必须同时通过完整 request fingerprint 和 `.shd`/`.iwdiag`
输出 SHA-256 校验。正式运行后又执行了一次相同配置复核：flat 和 canonical rough
来自 `stage2_verified`，其余 7 个 seed 均来自 `cache_verified`，没有混用旧的
1001-beam 数值。8 个 seed 的 mean delta TL 为 `0.320686 dB`，circular mean
phase 为 `-0.689310 rad`。

## 文件索引

### Prerequisite Bellhop-only basis

The PE comparison chain consumes the already frozen Bellhop validation-only
implementation and its comparison design audit:

- `reports/bellhop_internal_wall_implementation_report.md`
- `reports/pe_bellhop_pm_comparison_design_audit_report.md`

These documents establish the accepted flat/tilted/sinusoidal/vertical-tangent
internal-wall geometry, fixed-PM density convergence and finite-angle local
native↔internal `Reflect2D` covariance before the PE stages below. Superseded
one-off Bellhop reports and scripts remain quarantined under `cash/` and are not
inputs to the active chain.

### Validation entrypoints

- `scripts/validation/validate_pe_bellhop_pm_canonical_mapper.m`
- `scripts/validation/validate_pe_1d_validation_bridge.m`
- `scripts/validation/validate_pe_bellhop_pm_stage0_flat_source.m`
- `scripts/validation/validate_pe_bellhop_pm_constant_height_sign.m`
- `scripts/validation/validate_pe_bellhop_pm_numerical_budget.m`
- `scripts/validation/validate_pe_bellhop_pm_stage1_tier1.m`
- `scripts/validation/validate_pe_bellhop_pm_stage1_dimensionality.m`
- `scripts/validation/validate_pe_bellhop_pm_stage1_interpretation.m`
- `scripts/validation/validate_pe_bellhop_pm_frequency_extension.m`
- `scripts/validation/validate_pe_bellhop_pm_ensemble.m`
- `scripts/validation/validate_pe_bellhop_pm_amplitude_ensemble.m`
- `scripts/validation/validate_pe_bellhop_pm_model_discrepancy_statistics.m`
- `scripts/validation/postprocess_pe_bellhop_pm_model_discrepancy_statistics.py`（只读后处理，不调用 PE/Bellhop）

### Reports

- `reports/pe_bellhop_pm_canonical_mapper_report.md`
- `reports/pe_1d_validation_bridge_report.md`
- `reports/pe_bellhop_pm_stage0_flat_source_report.md`
- `reports/pe_bellhop_pm_constant_height_sign_audit.md`
- `reports/pe_bellhop_pm_numerical_error_budget.md`
- `reports/pe_bellhop_pm_stage1_tier1_report.md`
- `reports/pe_bellhop_pm_stage1_dimensionality_report.md`
- `reports/pe_bellhop_fixed_pm_4khz_comparison_report.md`
- `reports/pe_bellhop_fixed_pm_frequency_extension_report.md`
- `reports/pe_bellhop_pm_ensemble_comparison_report.md`
- `reports/pe_bellhop_pm_amplitude_ensemble_report.md`
- `reports/pe_bellhop_pm_model_discrepancy_statistics_report.md`

### Results

All generated artifacts are under `results/validation/` in the matching stage directories. Stage 3 machine-readable outputs are `ensemble_seed_summary.csv`, `ensemble_statistics.csv`, `ensemble_percentiles.csv`, `ensemble_checks.csv`, and `ensemble_comparison.mat`; Stage 4 outputs additionally include `per_seed_results.csv`, `convergence_by_sample_count.csv`, `bootstrap_confidence_intervals.csv`, `discrepancy_correlations.csv`, `outlier_audit.csv`, `roughness_bins.csv`, and `result.mat` under `pe_bellhop_pm_model_discrepancy_statistics/`.

## 下一步

可以进入独立的 PE rough-PM 与 Bellhop rotated-rough 解释阶段，但必须继续使用
`G=H_ref,rough/H_ref,flat`、数值预算和 dimensionality sensitivity 三项护栏。
当前 8-seed amplitude ensemble 仍是 first reduced statistic，不应解读为充分的
ocean Monte Carlo 或最差群时延结论；增加样本前应保持相同 `Sk_m3`、Kmax、profile
provenance 和 source 设置。不得为了减小 PE/Bellhop residual 修改任一核心传播或
反射公式。
