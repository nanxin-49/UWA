# PE surface angle--slope Model-2 G3 validation

Status: **NONLOCAL_EFFECT_REQUIRED**.

- Formula: `kzr=((1-s^2)kzi+2*s*kxi)/(1+s^2); R0*exp(+i*(kzi+kzr)*eta) per incident component`.
- No fitted parameter; production PE/marching unchanged.
- Flat Model-0/Model-2 relative difference: `1.83272751e-13`.

| case | M0 E | M1 E | M2 E | BH E | M0 phase | M1 phase | M2 phase | BH phase | M1/M2 E | M1/M2 phase | min kzr | nonreturn |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| weak_low_K | 0.00315757 | 0.000514551 | 0.000514113 | 8.09198e-06 | 0.00311764 | 0.000124496 | 0.000122674 | 5.88617e-06 | 1.001 | 1.015 | 4.61806 | 0 |
| region_II_low_K | 0.01609 | 0.00398989 | 0.00395441 | 3.88955e-05 | 0.0158945 | 0.00311257 | 0.00306696 | 2.79147e-05 | 1.009 | 1.015 | 4.3236 | 0 |
| strong_height_low_K | 0.0809081 | 0.0507808 | 0.0500519 | 9.97601e-05 | 0.0803364 | 0.0498092 | 0.049066 | 7.41103e-05 | 1.015 | 1.015 | 3.20833 | 0 |
| weak_high_K | 0.0229061 | 0.0215668 | 0.0215648 | 0.000136812 | 0.0130149 | 0.0106476 | 0.010644 | 9.77021e-05 | 1 | 1 | 4.38825 | 0 |

## Gates

- `flat_not_degraded`: true
- `low_K_not_degraded`: true
- `high_K_significant`: false
- `high_K_better_than_Model0`: true
- `returning_branch`: true
- `finite`: true
- `all`: false

Artifacts: `C:\Users\ASUS\.codex\worktrees\baf5\Explain\results\validation\pe_surface_operator_bie_reference\G3_angle_slope\G3_angle_slope_validation.mat`, `C:\Users\ASUS\.codex\worktrees\baf5\Explain\results\validation\pe_surface_operator_bie_reference\G3_angle_slope\G3_angle_slope_cases.csv`.
