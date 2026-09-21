# PE surface kz-aware Model-1 G2 validation

Status: **NORMAL_APPROXIMATION_CONFIRMED**; high-K limited: **1**.

- Model-0: `R0*exp(+i*2*k*eta)*psi_inc`.
- Model-1: `componentwise R0*exp(+i*2*kz(kx)*eta)*Psi_inc(kx)`.
- No BIE-fitted parameter; production PE is unchanged.
- Flat Model-0/Model-1 relative field difference: `1.83272751e-13`.
- Independent flat BIE error: `1.23551523e-07`; flat Bellhop validation error: `0.00589939216`.

| case | M0 E | M1 E | BH E | M0 phase | M1 phase | BH phase | M0 TL | M1 TL | BH TL | M0 rho | M1 rho | M0 aligned | M1 aligned | E gain | phase gain |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| weak_low_K | 0.00315757 | 0.000514551 | 8.09198e-06 | 0.00311764 | 0.000124496 | 5.88617e-06 | 0.00434838 | 0.00433654 | 4.82302e-05 | 0.99999502 | 0.99999987 | 0.00315588 | 0.000504365 | 6.137 | 25.04 |
| region_II_low_K | 0.01609 | 0.00398989 | 3.88955e-05 | 0.0158945 | 0.00311257 | 2.79147e-05 | 0.0217418 | 0.0216824 | 0.000235259 | 0.99987389 | 0.99999528 | 0.0158817 | 0.00307111 | 4.033 | 5.107 |
| strong_height_low_K | 0.0809081 | 0.0507808 | 9.97601e-05 | 0.0803364 | 0.0498092 | 7.41103e-05 | 0.0869964 | 0.0867584 | 0.000579793 | 0.99757897 | 0.99954045 | 0.0695923 | 0.0303204 | 1.593 | 1.613 |
| weak_high_K | 0.0229061 | 0.0215668 | 0.000136812 | 0.0130149 | 0.0106476 | 9.77021e-05 | 0.163732 | 0.162912 | 0.000831167 | 0.9997761 | 0.99980514 | 0.021162 | 0.0197419 | 1.062 | 1.222 |

## Gates

- `flat_not_degraded`: true
- `legacy_reproduced`: true
- `weak_preserved`: true
- `region_II_significant`: true
- `strong_height_significant`: true
- `consistent_all_cases`: true
- `finite`: true
- `all`: true

If high-K remains limited while G2 passes, the Goal proceeds to G3 local-slope coupling without fitting Model-1.

Artifacts: `results\validation\pe_surface_operator_bie_reference\G2_kz_aware\G2_kz_aware_validation.mat`, `results\validation\pe_surface_operator_bie_reference\G2_kz_aware\G2_kz_aware_cases.csv`.
