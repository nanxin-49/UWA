# Full Kirchhoff 2-D PE--BIE validation

Status: **diagnostic-only; production PE/BIE unchanged**.

Configuration: 4 kHz, c=1500 m/s, one transverse dimension, Gaussian sigma=0.3 m, exp(-i omega t), pressure-release Dirichlet surface. The same saved angular-spectrum source, receiver grid, Stage-0 M99 footprint and -40 dB pairwise threshold are used.

Field convention: the Full-Kirchhoff rough/flat ratio is first formed in the native exp(-i omega t) representation, then conjugated exactly once to match the frozen PE-comparison representation already used by `G_BIE`. No per-case phase choice, scalar fit, or calibration is used.

## Gate-0 flat surface

Refined Full-Kirchhoff versus the frozen flat reference: complex L2 0.000185546, magnitude L2 0.000129482, phase RMS 0.000180854 rad, N=411. Gate-0: **true**.

## Rough-case metrics

| case | method | complex L2 | magnitude L2 | phase RMS (rad) | magnitude corr. | phase corr. | N |
|---|---|---:|---:|---:|---:|---:|---:|
| weak_low_K | Model-0 | 0.00315757 | 0.000500626 | 0.00311764 | 0.979548 | 0.999995 | 411 |
| weak_low_K | kz-aware Model-1 | 0.000514551 | 0.000499263 | 0.000124496 | 0.999988 | 1 | 411 |
| weak_low_K | Full Kirchhoff | 0.00025349 | 0.000177857 | 0.000180615 | 0.814297 | 1 | 411 |
| strong_height_low_K | Model-0 | 0.0809081 | 0.0100156 | 0.0803364 | 0.977756 | 0.997629 | 411 |
| strong_height_low_K | kz-aware Model-1 | 0.0507808 | 0.00998813 | 0.0498092 | 0.998129 | 0.99959 | 411 |
| strong_height_low_K | Full Kirchhoff | 0.000257001 | 0.000176535 | 0.000186802 | 0.999377 | 1 | 411 |
| weak_high_K | Model-0 | 0.0229061 | 0.0188514 | 0.0130149 | 0.998704 | 0.999954 | 411 |
| weak_high_K | kz-aware Model-1 | 0.0215668 | 0.018757 | 0.0106476 | 0.99874 | 0.999981 | 411 |
| weak_high_K | Full Kirchhoff | 0.000271132 | 0.000177236 | 0.000205163 | 0.999822 | 1 | 411 |

## Interpretation

Full Kirchhoff is evaluated on the actual curve using surface position, unit normal, ds, Green function, source-normal Green derivative, incident field and incident normal derivative. It is not a flat-plane phase screen. Nominal/refined integration uses 2049/4097 surface samples.

| case | nominal/refined complex L2 | nominal/refined phase RMS (rad) |
|---|---:|---:|
| weak_low_K | 3.7196e-05 | 2.65017e-05 |
| strong_height_low_K | 3.71886e-05 | 2.66886e-05 |
| weak_high_K | 3.71944e-05 | 2.65905e-05 |

### Direct answers

- **Strong-height, low-K:** after the fixed convention mapping, Full Kirchhoff phase RMS is 0.000186802 rad (Model-1: 0.0498092 rad).
- **Weak, high-K:** after the fixed convention mapping, Full Kirchhoff magnitude L2 is 0.000177236 and phase RMS is 0.000205163 rad.
- **Combined diagnosis:** the earlier rough-case phase anomaly was caused by comparing a native exp(-i omega t) Kirchhoff ratio directly with an already-conjugated PE-comparison BIE ratio. The separate convention audit records the z/normal/Green-derivative checks and the `+/-4*k*eta` diagnostic. Production PE and the BIE/Kirchhoff kernels remain unchanged.

Artifacts: `results/validation/pe_bie_full_kirchhoff/pe_bie_full_kirchhoff_validation.mat`, `results/validation/pe_bie_full_kirchhoff/pe_bie_full_kirchhoff_cases.csv`; figures: `results/validation/pe_bie_full_kirchhoff/figures/`.
