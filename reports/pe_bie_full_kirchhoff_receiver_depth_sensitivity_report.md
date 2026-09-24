# Receiver-depth sensitivity: PE Model-1 vs Full Kirchhoff/BIE

Status: validation-only; no production PE, Full Kirchhoff, BIE or surface model was modified.

Setup: 2D, 4 kHz, c=1500 m/s, Tx=100 m, sigma=1 m, receiver depths [1,2,3,5,10,20,30] m. The sigma=1 accepted x grid, M99 footprint and energy weights are held fixed across depth.

## Interpretation

The table and figures report Model-1/BIE, Model-1/FK and FK/BIE using complex L2, magnitude L2, phase RMS, magnitude correlation and phase correlation. The ratio columns are E(z)/E(3 m). z_rx=3 m is not a universal maximum-error point: z=1 m is strongly affected, while the Model-1/BIE residual also grows again at larger depths. The high-K case does not improve monotonically with depth, and the strong-height case remains large away from the surface.

- `weak_low_K`: Model-1/BIE complex error is 0.169491 at 1 m, 0.0102641 at 3 m and 0.0676774 at 30 m (ratios 16.51 and 6.594).
- `strong_height_low_K`: Model-1/BIE complex error is 0.881096 at 1 m, 0.0713033 at 3 m and 0.545856 at 30 m (ratios 12.36 and 7.655).
- `weak_high_K`: Model-1/BIE complex error is 0.0643332 at 1 m, 0.0252671 at 3 m and 0.187208 at 30 m (ratios 2.546 and 7.409).

## Classification

- weak_low_K: **PARTIAL receiver-depth effect**; the residual is smallest near 2--3 m and grows at larger depth, while z=1 m is affected by the near-surface field.
- strong_height_low_K: **STRONG receiver-depth effect**; the residual changes by more than an order of magnitude, but remains large away from the surface.
- weak_high_K: **PARTIAL receiver-depth effect**; moving away from 3 m does not improve the residual monotonically and the high-K residual grows at larger depth.

The z=1 m FK nominal/refined relative difference is 0.172024 for weak_low_K. This is not converged enough for a hard FK/BIE physical conclusion at 1 m; depths >=2 m have stable FK-BIE control errors at the reported scale.

## Visualizations

Each case has metric-vs-depth, receiver-line (1, 3, 10, 30 m), shared-color-scale Model-1/FK/BIE magnitude maps, Model-1--FK magnitude/phase maps, and the three-case normalized-error/FK-BIE plots under `results/validation/pe_bie_full_kirchhoff_receiver_depth_sensitivity/figures/`.

Artifacts: `results/validation/pe_bie_full_kirchhoff_receiver_depth_sensitivity/receiver_depth_sensitivity.mat`, `results/validation/pe_bie_full_kirchhoff_receiver_depth_sensitivity/receiver_depth_sensitivity.csv`.
