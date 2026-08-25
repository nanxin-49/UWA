# Validation of the near-normal 2k height-phase approximation

## 1. Validation goal

This independent calculation studies only the replacement `(gamma_i + gamma_s) eta -> 2k eta` in the roughness height phase. The `2k-phase` branch is a diagnostic construction, not a new SSA theory. It does not validate the Kirchhoff approximation, SSA as an ocean truth model, PM/TMA/JONSWAP, PE propagation, or KStat. Angles are measured from the surface normal. The primary SSA angular distributions are normalized over 0--80 deg and compared with the same raw spectra normalized over 0--60 deg.

## 2. Analytic geometry result

| theta_i (deg) | max theta_s at 1% (deg) | max theta_s at 3% (deg) | max theta_s at 5% (deg) |
|---:|---:|---:|---:|
| 0.0 | 11.0 | 19.5 | 25.5 |
| 5.0 | 10.0 | 19.0 | 25.0 |
| 10.0 | 5.5 | 17.0 | 23.5 |
| 15.0 | NaN | 13.0 | 20.5 |
| 20.0 | NaN | NaN | 16.0 |

## 3. Gaussian SSA results

The lowest-order two-dimensional isotropic Gaussian SSA shape retains all angle-dependent geometry factors. Only the roughness-statistics argument `a` changes between the full and 2k-phase branches. Cases with `s_rms > 0.25` are skipped. Each Hankel integral uses 4001 points over `rho = 0--6l`; negative finite-quadrature tail values are warned about and clipped to zero before angular normalization.

| kh | kl | s_rms | theta95 full | theta99 full | theta95 2k | TV | weighted phase RMS | F(<0.1 rad) | F(<0.3 rad) |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.25 | 10 | 0.050 | 21.36 | 27.76 | 21.85 | 0.00987 | 0.00830 | 1.0000 | 1.0000 |
| 0.25 | 20 | 0.025 | 10.68 | 13.85 | 10.73 | 0.00242 | 0.00209 | 1.0000 | 1.0000 |
| 0.25 | 40 | 0.013 | 5.34 | 6.92 | 5.34 | 0.00060 | 0.00052 | 1.0000 | 1.0000 |
| 0.5 | 10 | 0.100 | 25.60 | 34.37 | 27.13 | 0.01816 | 0.02352 | 0.9938 | 1.0000 |
| 0.5 | 20 | 0.050 | 13.03 | 17.52 | 13.22 | 0.00457 | 0.00619 | 1.0000 | 1.0000 |
| 0.5 | 40 | 0.025 | 6.55 | 8.80 | 6.57 | 0.00114 | 0.00156 | 1.0000 | 1.0000 |
| 1 | 10 | 0.200 | 38.62 | 49.72 | 44.61 | 0.05111 | 0.10130 | 0.7757 | 0.9811 |
| 1 | 20 | 0.100 | 20.68 | 27.03 | 21.61 | 0.01513 | 0.03015 | 0.9863 | 1.0000 |
| 1 | 40 | 0.050 | 10.51 | 13.76 | 10.62 | 0.00380 | 0.00786 | 1.0000 | 1.0000 |
| 2 | 20 | 0.200 | 37.18 | 46.26 | 43.29 | 0.05841 | 0.19165 | 0.5152 | 0.8885 |
| 2 | 40 | 0.100 | 19.75 | 24.89 | 20.58 | 0.01615 | 0.05534 | 0.9224 | 0.9993 |

Before clipping, the largest negative Hankel value relative to the positive peak is 5.594e-06 for the full branch at `(kh, kl) = (2, 40)` and 5.586e-06 for the 2k branch at `(kh, kl) = (2, 40)`. The weak `kh=0.01, kl=20` sanity case gives 3.316e-07 and 3.750e-07, respectively.

## 4. Scattering-angle limit check

The table below shows the three largest 80-deg TV cases. Deltas are `80-deg result - 60-deg result`, computed from the same raw Hankel spectra.

| kh | kl | theta95 full: 60 | theta95 full: 80 | delta95 | theta99 full: 60 | theta99 full: 80 | delta99 | TV60 | TV80 | delta TV |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 2 | 20 | 37.1267 | 37.1811 | 0.0544 | 46.0374 | 46.2638 | 0.226 | 0.053551 | 0.058409 | 0.00486 |
| 1 | 10 | 38.3752 | 38.6201 | 0.245 | 48.7083 | 49.7223 | 1.01 | 0.045186 | 0.051114 | 0.00593 |
| 0.5 | 10 | 25.5940 | 25.5996 | 0.00568 | 34.3469 | 34.3748 | 0.0279 | 0.017926 | 0.018161 | 0.000235 |

Across all retained cases, the largest absolute full-phase changes are 0.245 deg in theta95 and 1.01 deg in theta99; the largest 2k-phase changes are 1.73 deg and 5.95 deg. The largest absolute TV change is 0.00593. Full results are in `ssa_angle_limit_comparison.csv`.

The worst-TV case remains `(kh, kl) = (2, 20)`. The full-phase percentiles are essentially stable: their largest theta95 and theta99 changes occur at `(1, 10)` and are 0.245 deg and 1.01 deg. The 2k-phase high-angle tail is less converged at 60 deg: `(1, 10)` changes by 1.73 deg in theta95 and 5.95 deg in theta99. TV changes by at most 0.00593 at `(1, 10)`, so the qualitative TV ranking and scale are stable, but the worst-case TV values are not numerically identical.

## 5. Figures

![Analytic phase-coefficient error](../results/validation/ssa_2k_phase/phase_coefficient_error_map.png)

![Representative normalized diffuse spectra](../results/validation/ssa_2k_phase/ssa_spectrum_examples.png)

![SSA 2k error summary](../results/validation/ssa_2k_phase/ssa_2k_error_summary.png)

## 6. Short interpretation

Across the retained Gaussian cases, the full-phase 90%, 95%, and 99% energy angles span 4.62--33.13 deg, 5.34--38.62 deg, and 6.92--49.72 deg, respectively, within the normalized 0--80 deg window.

Replacing only the height-phase argument by 2k gives TV distances from 0.00060 to 0.05841. The largest TV distance occurs at `(kh, kl) = (2, 20)`; its full-phase theta95 is 37.18 deg.

The largest energy-weighted phase RMS is 0.19165 rad at `(kh, kl) = (2, 20)`. This and TV identify the strongest tested diagnostic differences without defining a valid/invalid threshold.

9 of 11 retained cases place at least 90% of their normalized diffuse energy below 0.1 rad phase mismatch; 10 of 11 do so below 0.3 rad. These are diagnostic fractions, not universal applicability limits.

The result addresses only the tested Gaussian covariance and angular window. It does not establish applicability for other sea spectra, and the lowest-order SSA calculation is not an exact scattering result.
