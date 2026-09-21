# PE surface normal-approximation G1 audit

Status: **NORMAL_APPROXIMATION_MECHANISM_SUPPORTED_FOR_G2_TEST**

## Frozen incident angular spectrum

- `kx_rms = 2.35702261 rad/m`; `theta_rms = 8.14400577 deg`.
- `<kz>/k = 0.989950507855`; `<2(k-kz)> = 0.336761714 rad/m`; RMS `= 0.589715494 rad/m`.
- `|kx|` q50/q90/q95/q99: `1.60196 / 3.89047 / 4.60971 / 6.0809 rad/m`.
- `|theta|` q50/q90/q95/q99: `5.48642 / 13.4263 / 15.9693 / 21.2801 deg`.

## Case comparison

The componentwise diagnostic uses `exp(+i 2 kz eta)` before the unchanged receiver march. It is not fitted and is not yet a selectable PE model.

| case | A | K | PE-BIE phase RMS | mean-kz prediction | component prediction | exact angle correction | correlation | projection | angle residual | E_G | improvement |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| weak_low_K | 0.01 | 0.1 | 0.00311764 | 0.00234895 | 0.00411333 | 0.00311506 | 0.9992 | 1 | 0.000124496 | 0.000514551 | 25.04 |
| region_II_low_K | 0.05 | 0.1 | 0.0158945 | 0.0117448 | 0.0205667 | 0.0155753 | 0.98064 | 1.0007 | 0.00311257 | 0.00398989 | 5.107 |
| strong_height_low_K | 0.2 | 0.1 | 0.0803364 | 0.046979 | 0.0822667 | 0.0623047 | 0.78465 | 1.0117 | 0.0498092 | 0.0507808 | 1.613 |
| weak_high_K | 0.02 | 0.47 | 0.0130149 | 0.00476307 | 0.00834078 | 0.0073272 | 0.57519 | 1.0217 | 0.0106476 | 0.0215668 | 1.222 |

## Gate

The frozen mechanism gate requires median low-K error-vector correlation >=0.9, median low-K correction/error RMS >=0.5, and phase-RMS improvement >=1.1 in every case.

- `authoritative_inputs`: true
- `finite`: true
- `low_K_direction`: true
- `low_K_major_fraction`: true
- `multi_case_improvement`: true
- `all`: true

Interpretation: G1 tests whether the omitted finite-angle phase has the correct scale and error-vector direction across frozen cases. High-K residual after this diagnostic is evidence for a separate slope-coupling test, not permission to fit a coefficient.

Artifacts: `results\validation\pe_surface_operator_bie_reference\G1_normal_approximation\G1_normal_approximation_audit.mat`, `results\validation\pe_surface_operator_bie_reference\G1_normal_approximation\G1_normal_approximation_cases.csv`.
