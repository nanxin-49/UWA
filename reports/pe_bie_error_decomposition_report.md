# PE--BIE error decomposition

Status: **independent diagnostic; production physics unchanged**.

Four authoritative cases are reused: weak/low-K (R3), region-II medium height/low-K (R4), strong-height/low-K (R5), and weak/high-K (G0). Only missing pointwise Model-1 fields are recomputed with the frozen 4 kHz, c=1500 m/s PE configuration.

## Mask and metrics

The effective mask is the saved Stage-0 99% incident-energy footprint, finite samples, and samples where each compared field is at least -40 dB relative to the largest field in that pair. Weights are the saved Stage-0 energy weights renormalized on that mask. Complex and magnitude L2 values are relative to BIE magnitude energy; phase is wrapped `angle(a*conj(b))`; phase correlation is circular coherence of unit phasors. Classification uses phase RMS / magnitude relative L2: >1.5 phase-dominated, <2/3 magnitude-dominated, otherwise both significant.

| case | model | complex L2 | magnitude L2 | phase RMS (rad) | magnitude corr. | complex phase corr. | diagnosis | N |
|---|---|---:|---:|---:|---:|---:|---|---:|
| weak_low_K | Model-0 | 0.00315757 | 0.000500626 | 0.00311764 | 0.979548 | 0.999995 | phase-dominated | 411 |
| weak_low_K | kz-aware Model-1 | 0.000514551 | 0.000499263 | 0.000124496 | 0.999988 | 1 | magnitude-dominated | 411 |
| region_II_low_K | Model-0 | 0.01609 | 0.00250311 | 0.0158945 | 0.979353 | 0.999877 | phase-dominated | 411 |
| region_II_low_K | kz-aware Model-1 | 0.00398989 | 0.00249627 | 0.00311257 | 0.99974 | 0.999998 | both-significant | 411 |
| strong_height_low_K | Model-0 | 0.0809081 | 0.0100156 | 0.0803364 | 0.977756 | 0.997629 | phase-dominated | 411 |
| strong_height_low_K | kz-aware Model-1 | 0.0507808 | 0.00998813 | 0.0498092 | 0.998129 | 0.99959 | phase-dominated | 411 |
| weak_high_K | Model-0 | 0.0229061 | 0.0188514 | 0.0130149 | 0.998704 | 0.999954 | both-significant | 411 |
| weak_high_K | kz-aware Model-1 | 0.0215668 | 0.018757 | 0.0106476 | 0.99874 | 0.999981 | magnitude-dominated | 411 |

## Interpretation

- `weak_low_K`: Model-0 **phase-dominated**, Model-1 **magnitude-dominated**; phase/magnitude ratios 6.227 -> 0.2494; kz-aware complex-error improvement 6.137-fold.
- `region_II_low_K`: Model-0 **phase-dominated**, Model-1 **both-significant**; phase/magnitude ratios 6.35 -> 1.247; kz-aware complex-error improvement 4.033-fold.
- `strong_height_low_K`: Model-0 **phase-dominated**, Model-1 **phase-dominated**; phase/magnitude ratios 8.021 -> 4.987; kz-aware complex-error improvement 1.593-fold.
- `weak_high_K`: Model-0 **both-significant**, Model-1 **magnitude-dominated**; phase/magnitude ratios 0.6904 -> 0.5677; kz-aware complex-error improvement 1.062-fold.

## Next-step decision rule

The observed diagnosis is case dependent: Model-0 is phase-dominated for all low-K cases, while the kz-aware correction removes most finite-angle phase error in weak/medium low-K cases. Strong height retains a substantial phase residual after the correction, and high-K retains a non-negligible magnitude residual with only modest complex-error improvement. Therefore amplitude alone is not a sufficient explanation for the PE--BIE discrepancy; the next physics investigation should prioritize complete Kirchhoff/nonlocal phase coupling while separately auditing amplitude/geometric factors.

- If magnitude error is much smaller than phase error, prioritize complete Kirchhoff/nonlocal phase coupling.
- If both remain significant, implement a complete Kirchhoff surface integral and audit obliquity, normal derivative, and surface Jacobian terms.
- If magnitude error dominates, first inspect a simple amplitude/geometric correction; do not immediately introduce a nonlocal operator.

Artifacts: `results/validation/pe_bie_error_decomposition/pe_bie_error_decomposition_validation.mat`, `results/validation/pe_bie_error_decomposition/pe_bie_error_decomposition_cases.csv`.
