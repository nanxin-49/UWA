# Adaptive Reduced Kirchhoff aperture: Stage 1 database

Status: **validation-only; production PE, BIE, Full Kirchhoff, and surface model unchanged**.

This first stage reuses the accepted deterministic A sweep (`K=0.10`) and K sweep (`A=0.02`), omitting their duplicate `A=0.02, K=0.10` point. It is not yet a complete 6-by-6 A-K matrix. The operator retains every accepted Full-Kirchhoff term and truncates only source points outside `|x_r-x_s|<L`.

L grid: `4, 8, 16, 32, 64, 128, Inf m`. Gate: `Reduced Kirchhoff-BIE complex L2 < 0.001`. `requires_full_K=true` means no finite tested L passed; `L_min=Inf` is then the complete Full-Kirchhoff reference.

## Adaptive-aperture database

| case | A m | K rad/m | A K | A K^2 1/m | L_min m | Ec at L_min | phase error rad | magnitude error | requires full K | PE-BIE Ec | FK-BIE Ec |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|
| A0.005_K0.1 | 0.005 | 0.1 | 0.0005 | 5e-05 | 128 | 0.000253645 | 0.000180822 | 0.000177868 | false | 0.00157783 | 0.000253645 |
| A0.01_K0.1 | 0.01 | 0.1 | 0.001 | 0.0001 | 128 | 0.00025365 | 0.000180728 | 0.00017797 | false | 0.00315757 | 0.00025365 |
| A0.02_K0.1 | 0.02 | 0.1 | 0.002 | 0.0002 | 128 | 0.000253668 | 0.000180368 | 0.000178362 | false | 0.00633038 | 0.000253668 |
| A0.05_K0.1 | 0.05 | 0.1 | 0.005 | 0.0005 | 128 | 0.000253805 | 0.000178554 | 0.000180372 | false | 0.01609 | 0.000253805 |
| A0.1_K0.1 | 0.1 | 0.1 | 0.01 | 0.001 | 128 | 0.00025436 | 0.000178824 | 0.000180892 | false | 0.0339984 | 0.00025436 |
| A0.2_K0.1 | 0.2 | 0.1 | 0.02 | 0.002 | 128 | 0.000256982 | 0.000186783 | 0.000176529 | false | 0.0809081 | 0.000256982 |
| A0.02_K0.05 | 0.02 | 0.05 | 0.001 | 5e-05 | 128 | 0.000253613 | 0.000180431 | 0.000178226 | false | 0.00941375 | 0.000253613 |
| A0.02_K0.2 | 0.02 | 0.2 | 0.004 | 0.0008 | 128 | 0.000254226 | 0.000180389 | 0.00017913 | false | 0.00820412 | 0.000254226 |
| A0.02_K0.3 | 0.02 | 0.3 | 0.006 | 0.0018 | 128 | 0.000256649 | 0.000183803 | 0.000179103 | false | 0.0116157 | 0.000256649 |
| A0.02_K0.47 | 0.02 | 0.47 | 0.0094 | 0.004418 | 128 | 0.000271269 | 0.000205252 | 0.00017734 | false | 0.0229061 | 0.000271269 |
| A0.02_K0.7 | 0.02 | 0.7 | 0.014 | 0.0098 | 128 | 0.00033158 | 0.000279571 | 0.000178534 | false | 0.0484307 | 0.00033158 |

## Interpretation

- Finite gate passes: `11/11`.
- Finite selected L range: `128` to `128 m`.
- Every `L=Inf` reconstruction was checked directly against its accepted Full-Kirchhoff field; max relative discrepancy is `0`.
- Because this Stage-1 database follows two one-dimensional sweeps rather than a complete matrix, it can diagnose monotonic/censored trends but cannot support an empirical law `L=f(A,K)` or a validated adaptive formula.

The next step should be chosen from the observed gate pattern. If finite L values vary across the one-dimensional sweeps, run only discriminating cells of the A-K matrix to test whether height, slope, or curvature is the better organizer. If the criterion is uniformly censored at 128 m, the current aperture definition is effectively global at this accuracy target and an adaptive-aperture PE boundary operator is not yet justified; low-rank/local-stationary kernel work is then more informative than fitting L.

Artifacts: `results/validation/pe_bie_full_kirchhoff_adaptive_aperture_stage1/pe_bie_full_kirchhoff_adaptive_aperture_stage1.mat`, `results/validation/pe_bie_full_kirchhoff_adaptive_aperture_stage1/adaptive_aperture_database.csv`, `results/validation/pe_bie_full_kirchhoff_adaptive_aperture_stage1/adaptive_aperture_scan.csv`; figures: `results/validation/pe_bie_full_kirchhoff_adaptive_aperture_stage1/figures/`.
