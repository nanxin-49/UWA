# Strict-normal approximation Gaussian-width validation

Conclusion: **STRICT_NORMAL_APPROXIMATION_PARTIAL**

All cases use 4 kHz, the fixed BIE/Bellhop/receiver convention, and 10,001 beams. Within each surface family, only sigma and its analytic Gaussian SBP change.

The low-K BIE matrix is assembled once and solved for three independent Gaussian right-hand sides. A separate scalar-versus-multiple-RHS audit agreed to `2.84e-16` in receiver field, `3.44e-16` in density, and `1.67e-16` in boundary residual. This is an algebraic batching optimization, not a BIE physics change.

Bellhop retains the established `+/-30 deg` cap and uses the unchanged analytic Gaussian formula. For narrower sources only, the numerical fan stops where directivity reaches `0.006`, above the validation binary termination threshold; omitted angular power is tabulated and no tail is refitted.

## Low-slope A=0.05 m, K=0.10 rad/m

| sigma m | theta rms | theta95 | theta99 | footprint99 m | M0 E | M1 E | BH E | M0 phase | M1 phase | M0 TL | M1 TL | M0-M1 E | Region I |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---:|
| 0.3 | 8.144 | 15.969 | 21.28 | 40.039 | 0.01609 | 0.00398989 | 3.88955e-05 | 0.0158945 | 0.00311257 | 0.0217418 | 0.0216824 | 0.0155767 | false |
| 0.5 | 4.8536 | 9.5468 | 12.509 | 23.047 | 0.00972358 | 0.00500264 | 4.69288e-05 | 0.00912794 | 0.00353931 | 0.029109 | 0.0307095 | 0.00841837 | false |
| 1 | 2.4202 | 4.7007 | 6.2731 | 11.328 | 0.00800711 | 0.00804448 | 0.000502282 | 0.00459875 | 0.00432398 | 0.0569331 | 0.058922 | 0.00158309 | false |
| 2 | 1.2093 | 2.3484 | 3.1319 | 6.6406 | 0.0120502 | 0.0121528 | 0.00422061 | 0.00468893 | 0.0046847 | 0.0964189 | 0.0974015 | 0.000231016 | false |

Numerical support audit:

| sigma m | Bellhop fan half-angle | omitted angular power | outer-5% incident energy | BIE boundary residual |
|---:|---:|---:|---:|---:|
| 0.3 | 30 | 0.00032301 | 3.15584e-06 | 3.30533e-09 |
| 0.5 | 22.2667 | 6.58127e-06 | 9.80912e-15 | 2.01135e-09 |
| 1 | 10.9857 | 6.19236e-06 | 6.92123e-27 | 1.48637e-09 |
| 2 | 5.475 | 6.10758e-06 | 5.86259e-27 | 1.03062e-09 |

## High-K A=0.02 m, K=0.47 rad/m

| sigma m | theta rms | theta95 | theta99 | footprint99 m | M0 E | M1 E | BH E | M0 phase | M1 phase | M0 TL | M1 TL | M0-M1 E | Region I |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---:|
| 0.3 | 8.144 | 15.969 | 21.28 | 40.039 | 0.0229061 | 0.0215668 | 0.000136812 | 0.0130149 | 0.0106476 | 0.163732 | 0.162912 | 0.00732784 | n/a |
| 2 | 1.2093 | 2.3484 | 3.1319 | 6.6406 | 0.0255935 | 0.0255763 | 0.00427532 | 0.0110112 | 0.0110085 | 0.200692 | 0.200538 | 0.000155735 | n/a |

Numerical support audit:

| sigma m | Bellhop fan half-angle | omitted angular power | outer-5% incident energy | BIE boundary residual |
|---:|---:|---:|---:|---:|
| 0.3 | 30 | 0.00032301 | 3.15584e-06 | 4.24555e-09 |
| 2 | 5.475 | 6.10758e-06 | 5.86259e-27 | 1.20111e-09 |

## Strong-height A=0.20 m, K=0.10 rad/m

| sigma m | theta rms | theta95 | theta99 | footprint99 m | M0 E | M1 E | BH E | M0 phase | M1 phase | M0 TL | M1 TL | M0-M1 E | Region I |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---:|
| 0.3 | 8.144 | 15.969 | 21.28 | 40.039 | 0.0809081 | 0.0507808 | 9.97601e-05 | 0.0803364 | 0.0498092 | 0.0869964 | 0.0867584 | 0.0622654 | n/a |
| 2 | 1.2093 | 2.3484 | 3.1319 | 6.6406 | 0.087081 | 0.0872847 | 0.0169177 | 0.0749878 | 0.0749563 | 0.385599 | 0.389547 | 0.00092382 | n/a |

Numerical support audit:

| sigma m | Bellhop fan half-angle | omitted angular power | outer-5% incident energy | BIE boundary residual |
|---:|---:|---:|---:|---:|
| 0.3 | 30 | 0.00032301 | 3.15584e-06 | 1.4619e-08 |
| 2 | 5.475 | 6.10758e-06 | 5.86259e-27 | 4.1207e-09 |

## Angular trend

- Model-0 phase fit versus `theta_rms^2`: intercept `0.00423885463 rad`, slope `0.586123606 rad/rad^2`, R2 `0.988118`.
- Angle-specific excess phase fit: intercept `-0.000217741095 rad`, slope `0.659270105 rad/rad^2`, R2 `0.983925`.
- Model-0 phase fit versus `1-<cos(theta)>`: intercept `0.00423214664 rad`, slope `1.17877908 rad`, R2 `0.988272`.
- Broadest sampled low-K point passing both complex and phase thresholds: `theta_rms=4.85356 deg`, `theta95=9.54683 deg`.

  This is not a monotone total-field guarantee: the narrower sigma=2 m point fails the complex-error threshold because a common amplitude residual grows. The phase-only evidence supports `theta95 <= 9.54683 deg` for this low-K family, while Model-0/Model-1 agreement reaches `E<0.002` by `theta95=4.70073 deg`.

The linear fit is diagnostic. Full validation requires monotone total errors plus all numerical guards; the PARTIAL classification requires phase contraction and Model-0/Model-1 merging while reporting the remaining nonmonotone residual.

## Checks

- `all_numerical`: true
- `model0_E_monotone`: false
- `model0_phase_monotone`: false
- `model_gap_monotone`: true
- `theta2_phase_trend`: true
- `theta2_excess_trend`: true
- `lowK_has_region_I`: false
- `lowK_has_complex_phase_range`: true
- `phase_contraction`: true
- `model_gap_contraction`: true
- `angle_mechanism`: true
- `highK_residual`: true

The sampled angle is an engineering bound for this 4 kHz Gaussian/source/surface family, not a universal deep-ocean theorem. High-K residual is interpreted separately from strict-normal angle error.

## Interpretation

- Narrowing the low-K incident spectrum reduces Model-0 phase RMS from `0.0158945` to a `0.00459875--0.00468893 rad` floor and reduces the Model-0/Model-1 gap from `0.0155767` to `0.000231016`. This validates a finite-angle strict-normal contribution.
- Total complex error and TL do not improve monotonically; no mandatory point passes the complete Region-I gate. Strict-normal angle narrowing alone is therefore insufficient for an absolute field-accuracy claim.
- At high K and `theta95=2.34839 deg`, Model-0 and Model-1 agree to `0.000155735`, yet their BIE errors remain `0.0255935/0.0255763` with phase RMS `0.0110112/0.0110085 rad`. The remaining error is consistent with nonlocal/spectral-coupling physics rather than the `2*k*eta` angle approximation.

- Strong-height check (`A=0.20 m, K=0.10 rad/m`): narrowing from `theta95=15.9693` to `2.34839 deg` changes the PE Model-0/Model-1 gap from `0.0622654` to `0.00092382`, while PE--BIE complex errors remain `0.0809081` and `0.087081`. This separates the finite-angle effect from the strong-height phase-screen residual.

## Optional sigma=4 m diagnostic

PE/BIE remained finite (`theta95=1.22986 deg`, Model-0/Model-1 gap `4.90776e-05`), but Bellhop rough/flat contained `933` NaNs, so this point is excluded from hard gates and trend fits. No Bellhop parameter was tuned.

Artifacts: `C:\Users\ASUS\.codex\worktrees\baf5\Explain\results\validation\pe_strict_normal_sigma_sweep\strict_normal_sigma_sweep_validation.mat`, `C:\Users\ASUS\.codex\worktrees\baf5\Explain\results\validation\pe_strict_normal_sigma_sweep\strict_normal_lowK.csv`, `C:\Users\ASUS\.codex\worktrees\baf5\Explain\results\validation\pe_strict_normal_sigma_sweep\strict_normal_highK.csv`.
