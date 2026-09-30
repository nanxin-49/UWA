# Finite-aperture Reduced Kirchhoff prototype

Status: **validation-only; production PE, BIE, and surface model unchanged**.

The prototype is the accepted Full-Kirchhoff line operator with all surface position, Green function, incident normal derivative, Green normal derivative and exact ds terms retained. Only source contributions with `|x_r-x_s|<L` are kept. `L=Inf` is exactly the accepted Full-Kirchhoff field.

Metrics use the Stage-1 fixed Stage-0 M99 and Full-field -40 dB mask with saved incident-energy weights. Reduction factor is `E_PE,BIE / E_RK,BIE`; it is not a fitted correction.

## Aperture scan relative to BIE

| case | L (m) | complex L2 | magnitude L2 | phase RMS (rad) | magnitude corr. | phase corr. | reduction factor |
|---|---:|---:|---:|---:|---:|---:|---:|
| weak_low_K | 4 | 0.0642782 | 0.0323728 | 0.0555036 | 0.00781327 | 0.99846 | 0.0491235 |
| weak_low_K | 8 | 0.0348447 | 0.0191154 | 0.0291303 | 0.0149452 | 0.999576 | 0.0906185 |
| weak_low_K | 16 | 0.0204315 | 0.0154316 | 0.0133911 | 0.0142377 | 0.99991 | 0.154545 |
| weak_low_K | 32 | 0.0108564 | 0.00768369 | 0.00767134 | 0.03529 | 0.999971 | 0.29085 |
| weak_low_K | 64 | 0.00154372 | 0.00111902 | 0.00106346 | 0.21737 | 0.999999 | 2.04543 |
| weak_low_K | Inf | 0.00025365 | 0.00017797 | 0.000180728 | 0.814035 | 1 | 12.4486 |
| strong_height_low_K | 4 | 0.0644478 | 0.0513857 | 0.0389948 | 0.0966735 | 0.99924 | 1.2554 |
| strong_height_low_K | 8 | 0.0348643 | 0.023585 | 0.025643 | 0.213343 | 0.999671 | 2.32066 |
| strong_height_low_K | 16 | 0.0204458 | 0.0146575 | 0.0142532 | 0.319261 | 0.999898 | 3.95719 |
| strong_height_low_K | 32 | 0.0108713 | 0.00768434 | 0.00769805 | 0.548195 | 0.99997 | 7.44233 |
| strong_height_low_K | 64 | 0.00154192 | 0.00109032 | 0.00109034 | 0.977007 | 0.999999 | 52.4722 |
| strong_height_low_K | Inf | 0.000256982 | 0.000176529 | 0.000186783 | 0.999377 | 1 | 314.84 |
| weak_high_K | 4 | 0.0643385 | 0.0394568 | 0.0508487 | 0.225565 | 0.998708 | 0.356025 |
| weak_high_K | 8 | 0.0348526 | 0.02377 | 0.0254552 | 0.382034 | 0.999676 | 0.657227 |
| weak_high_K | 16 | 0.0204348 | 0.0152265 | 0.0135958 | 0.52581 | 0.999908 | 1.12093 |
| weak_high_K | 32 | 0.0108587 | 0.00769568 | 0.00766451 | 0.773193 | 0.999971 | 2.10946 |
| weak_high_K | 64 | 0.00154673 | 0.00110874 | 0.00107795 | 0.993093 | 0.999999 | 14.8093 |
| weak_high_K | Inf | 0.000271269 | 0.00017734 | 0.000205252 | 0.999821 | 1 | 84.4405 |

## Uniform support test

Requested criterion: every representative case must satisfy `Reduced Kirchhoff-BIE complex L2 < 1e-3`.

- `L=4 m`: maximum case error `0.0644478`; all cases below threshold: **false**.
- `L=8 m`: maximum case error `0.0348643`; all cases below threshold: **false**.
- `L=16 m`: maximum case error `0.0204458`; all cases below threshold: **false**.
- `L=32 m`: maximum case error `0.0108713`; all cases below threshold: **false**.
- `L=64 m`: maximum case error `0.00154673`; all cases below threshold: **false**.
- `L=Inf m`: maximum case error `0.000271269`; all cases below threshold: **true**.

No tested finite support meets the `1e-3` criterion simultaneously for all three cases.

## Interpretation

Finite support is clearly better than the local phase screen before reaching `L=Inf`: the reduction factor is greater than one whenever the retained aperture captures a substantial fraction of the complete coherent integral. However, the required support is case-dependent, consistent with Stage-1: weak/low-K needs the largest support, strong-height/low-K the smallest, and weak/high-K is intermediate.

This prototype is already the simplest physically controlled kernel reduction: it is a truncated Full-Kirchhoff kernel, not a correction formed from `u_FK-u_PS`. A fixed translation-invariant convolution kernel was not claimed, because the retained Green and normal-derivative kernels depend on receiver position, surface height and source/receiver geometry. The next reduction study should test whether those kernels admit a low-rank or locally stationary approximation after the finite-aperture support is selected.

**Recommendation:** use an adaptive aperture as the next prototype candidate, parameterized by the observed height/slope/curvature regime; do not freeze `L=16 m` or `L=32 m` without a stated error gate.

Artifacts: `results/validation/pe_bie_full_kirchhoff_finite_aperture/pe_bie_full_kirchhoff_finite_aperture.mat`, `results/validation/pe_bie_full_kirchhoff_finite_aperture/pe_bie_full_kirchhoff_finite_aperture.csv`; figures: `results/validation/pe_bie_full_kirchhoff_finite_aperture/figures/`.
