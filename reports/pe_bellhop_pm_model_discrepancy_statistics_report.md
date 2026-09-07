# PE--Bellhop PM model-discrepancy statistical study

状态：**PRELIMINARY_MODEL_DISCREPANCY**

本报告基于固定 4 kHz、Stage 3B independent Gaussian coefficient-amplitude ensemble。用户要求在第 50 个完整 seed 后停止，因此最终样本为 **M=50**（260001--260050）；未运行或生成 260051 及以后样本。未修改 PE/Bellhop 核心物理。Bootstrap repetitions=2000, seed=42032。

## Running-prefix convergence

| M | mean delta TL | std | median | p05 | p95 | PE power | BH power | circular mean phase | circular std | R |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 8 | 0.32068639 | 0.64881613 | 0.15333663 | -0.40329861 | 1.2598687 | 1.030838 | 0.95652236 | -0.68931037 | 0.83581852 | 0.70518415 |
| 16 | -0.00010535314 | 0.96250952 | -0.0094920437 | -1.593703 | 1.3364753 | 1.0029116 | 1.0026701 | -0.50496531 | 0.7043924 | 0.78029414 |
| 24 | 0.013296669 | 0.89692079 | -0.0094920437 | -1.5263707 | 1.2640147 | 1.0044896 | 1.0012788 | -0.42918432 | 0.59333699 | 0.83859753 |
| 32 | 0.0038043155 | 1.0323161 | 0.068362941 | -1.6387229 | 1.3862456 | 1.0090472 | 1.010343 | -0.34566819 | 0.52635912 | 0.87063978 |
| 48 | -0.052198634 | 0.95012642 | -0.0094920437 | -1.5116275 | 1.2598687 | 0.99830936 | 1.011428 | -0.3928919 | 0.60263164 | 0.83394949 |
| 50 | -0.0074785801 | 0.96223247 | 0.068362941 | -1.5042559 | 1.412717 | 1.0045099 | 1.0070072 | -0.37534471 | 0.59470097 | 0.83791835 |

24->32 gates: mean delta-TL change 0.0094923533 dB (PASS); PE/Bellhop power relative changes 0.0045372103/0.0090525997 (PASS); circular phase change 0.083516126 rad (PASS); bootstrap mean-delta-TL half-width 0.25835004 dB (FAIL).

## Bootstrap confidence intervals (95%)

| metric | estimate | low | high | half-width |
|---|---:|---:|---:|---:|
| mean_delta_tl_db | -0.0074785801 | -0.26287501 | 0.25382508 | 0.25835004 |
| median_delta_tl_db | 0.068362941 | -0.20516853 | 0.40588331 | 0.30552592 |
| mean_pe_power | 1.0045099 | 0.97543663 | 1.0361784 | 0.030370894 |
| mean_bh_power | 1.0070072 | 0.97528029 | 1.0430127 | 0.033866206 |
| mean_power_difference | -0.0024973135 | -0.063746419 | 0.059014021 | 0.06138022 |
| circular_mean_phase_rad | -0.37534471 | -0.53860432 | -0.25106816 | 0.14376808 |
| resultant_length | 0.83791835 | 0.74029413 | 0.92572642 | 0.092716143 |

## Final M=50 statistics

Mean delta TL = **-0.0074785801 dB**, std = 0.96223247 dB; circular mean phase = **-0.37534471 rad**, circular std = 0.59470097 rad; PE/Bellhop reflected power means = 1.0045099/1.0070072.

## Numerical/applicability guards

- all_finite: PASS
- all_profile_provenance: PASS
- all_beams_hit: PASS
- all_failed_rays_zero: PASS
- all_rejected_rays_zero: PASS
- all_non_grazing: PASS
- all_wall_residual: PASS
- all_phase: PASS
- all_beam_state: PASS
- all_positive_post_range: PASS
- all_pe_edge: PASS
- all: PASS

## Dimensionality reference

Existing 1T->2T sensitivity: -0.158855462885639 dB and 0.0125052786464214 rad; reported alongside, not subtracted from, model discrepancy.

## Geometry predictors and outliers

The full Pearson/Spearman coefficients, p-values and bootstrap intervals are in `discrepancy_correlations.csv`. The largest absolute correlations with delta TL are listed below (these are exploratory associations, not causal claims).

| predictor | Pearson r | p | Spearman rho | p |
|---|---:|---:|---:|---:|
| profile_rms_curvature_per_m | 0.128346 | 0.374396 | 0.148331 | 0.303934 |
| hit_curvature_p95 | -0.122701 | 0.395938 | -0.0618968 | 0.669368 |
| profile_rms_slope | 0.109056 | 0.450914 | 0.141513 | 0.326948 |
| profile_max_slope | -0.0910249 | 0.529565 | 0.0397119 | 0.78423 |
| surface_rms_eta_m | 0.089537 | 0.536343 | 0.154286 | 0.284713 |

Top outlier rows are retained without deletion in `outlier_audit.csv`; all are classified as legitimate realizations (A) because the numerical/applicability guards pass.

## Interpretation

All 50 completed seeds passed the PE edge/seam, Bellhop wall-hit, non-grazing, wall residual, pressure-release phase, beam-state and positive-post-range guards. The comparison therefore remains a model-discrepancy study between two independently validated approximate models; it does not claim either model is exact. Native backward-range amplitude remains diagnostic only. The 24->32 engineering gates are reported unchanged; the user-requested M=50 truncation is a practical sample-size decision, not a new convergence claim.

## Artifacts

Outputs are in `results/validation/pe_bellhop_pm_model_discrepancy_statistics/`: `per_seed_results.csv`, `ensemble_statistics.csv`, `convergence_by_sample_count.csv`, `bootstrap_confidence_intervals.csv`, `surface_geometry_statistics.csv`, `discrepancy_correlations.csv`, `outlier_audit.csv`, `roughness_bins.csv`, `result.mat`, and `figures/`.
