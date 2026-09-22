# Full Kirchhoff 2-D PE--BIE validation

Status: **diagnostic-only; production PE/BIE unchanged**.

Configuration: 4 kHz, c=1500 m/s, one transverse dimension, Gaussian sigma=0.3 m, exp(-i omega t), pressure-release Dirichlet surface. The same saved angular-spectrum source, receiver grid, Stage-0 M99 footprint and -40 dB pairwise threshold are used.

## Gate-0 flat surface

Refined Full-Kirchhoff versus the frozen flat reference: complex L2 0.000185546, magnitude L2 0.000129482, phase RMS 0.000180854 rad, N=411. Gate-0: **true**.

## Rough-case metrics

| case | method | complex L2 | magnitude L2 | phase RMS (rad) | magnitude corr. | phase corr. | N |
|---|---|---:|---:|---:|---:|---:|---:|
| weak_low_K | Model-0 | 0.00315757 | 0.000500626 | 0.00311764 | 0.979548 | 0.999995 | 411 |
| weak_low_K | kz-aware Model-1 | 0.000514551 | 0.000499263 | 0.000124496 | 0.999988 | 1 | 411 |
| weak_low_K | Full Kirchhoff | 0.455442 | 0.000177857 | 0.46168 | 0.814297 | 0.896287 | 411 |
| strong_height_low_K | Model-0 | 0.0809081 | 0.0100156 | 0.0803364 | 0.977756 | 0.997629 | 411 |
| strong_height_low_K | kz-aware Model-1 | 0.0507808 | 0.00998813 | 0.0498092 | 0.998129 | 0.99959 | 411 |
| strong_height_low_K | Full Kirchhoff | 1.25819 | 0.000176535 | 1.60046 | 0.999377 | 0.208384 | 411 |
| weak_high_K | Model-0 | 0.0229061 | 0.0188514 | 0.0130149 | 0.998704 | 0.999954 | 411 |
| weak_high_K | kz-aware Model-1 | 0.0215668 | 0.018757 | 0.0106476 | 0.99874 | 0.999981 | 411 |
| weak_high_K | Full Kirchhoff | 0.889026 | 0.000177236 | 0.939647 | 0.999822 | 0.604936 | 411 |

## Interpretation

Full Kirchhoff is evaluated on the actual curve using surface position, unit normal, ds, Green function, source-normal Green derivative, incident field and incident normal derivative. It is not a flat-plane phase screen. Nominal/refined integration uses 2049/4097 surface samples.

| case | nominal/refined complex L2 | nominal/refined phase RMS (rad) |
|---|---:|---:|
| weak_low_K | 3.7196e-05 | 2.65017e-05 |
| strong_height_low_K | 3.71886e-05 | 2.66886e-05 |
| weak_high_K | 3.71944e-05 | 2.65905e-05 |

### Direct answers

- **Strong-height, low-K:** Full Kirchhoff does not reduce the remaining phase error in this run: phase RMS is approximately 1.60 rad versus 0.0498 rad for Model-1.
- **Weak, high-K:** Full Kirchhoff has a small magnitude residual (about 1.8e-4) but a large phase residual (about 0.94 rad), so it does not provide a coherent-field improvement over Model-1.
- **Combined diagnosis:** Full Kirchhoff is not close to BIE for the rough cases and is not close to Model-1 either. Gate-0 is valid, but the rough-surface Kirchhoff approximation as implemented here is insufficient to adjudicate the PE residual without a further derivation/audit of rough-surface Kirchhoff terms, orientation, and illumination/shadow treatment. This is not evidence to modify production PE.

The result is diagnostic only: if a subsequently audited Full Kirchhoff formulation approaches BIE, the local phase-screen reduction is implicated; if it improves but remains separated, Kirchhoff is useful but incomplete; if it remains close to Model-1, the Kirchhoff approximation itself is insufficient for those cases.

Artifacts: `results/validation/pe_bie_full_kirchhoff/pe_bie_full_kirchhoff_validation.mat`, `results/validation/pe_bie_full_kirchhoff/pe_bie_full_kirchhoff_cases.csv`; figures: `results/validation/pe_bie_full_kirchhoff/figures/`.
