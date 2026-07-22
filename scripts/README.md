# Script Directory

Scripts are grouped by role and are meant to run from any MATLAB current folder. Each script starts with `bootstrap_project.m`, which adds the project root and script folders to the MATLAB path and changes the working output directory to `results/<category>`.

## Subdirectories

- `validation/`: smoke tests, invariants, communication-chain checks, and disabled-path regressions.
- `comparisons/`: surface-model, bubble-model, and coherent/incoherent comparison studies.
- `experiments/`: Monte Carlo runs, parameter sweeps, and calibration scripts.
- `reporting/`: figure generation, summary tables, visualizations, and report builders.

## Running Examples

```matlab
run('scripts/validation/validate_surface_wavefield_visualization_vertical.m')
run('scripts/comparisons/compare_specular_incoherent_surface_reflection_vertical.m')
run('scripts/experiments/sweep_monte_carlo_surface_channel_vertical.m')
run('scripts/reporting/plot_c35_core_heatmaps_vertical.m')
```

## Raw-PM and Joint-Frequency K-Stat Prerequisite Validation

- `validation/validate_raw_pm_grid_coverage_vertical.m`: scans PM aperture, FFT-grid resolution, wind-speed coverage, and same-dx large-PM-grid to small-PE-window energy mapping.
- `validation/validate_kstat_joint_frequency_boundary_vertical.m`: boundary-only comparison of one shared explicit Gaussian PM surface, independent-frequency kstat, and covariance+pseudocovariance joint-frequency kstat at U=5 m/s and F=32.
- `validation/audit_pe_caching_vertical.m`: counts cacheable PE marches/FFTs, measures the unchanged public-path baseline, and estimates batched cache memory.

These scripts do not change the public defaults or communication chain. Detailed interpretation is in `reports/raw_pm_joint_kstat_prerequisite_validation_report.md`.

## Cached Joint-Kstat Receiver Validation

- `validation/validate_cached_joint_kstat_pe_receiver_vertical.m`: runs the U=5, F=32, 128/64 train/test reflected-only kdomain/independent/joint comparison through the independent cached PE executor.
- `validation/validate_cached_pe_public_consistency_vertical.m`: checks the cached executor against the unchanged public kdomain path with the exact same double-precision boundary input.
- `validation/analyze_receiver_properness_null_vertical.m`: calibrates finite-sample receiver `||P||_F/||C||_F` under a proper complex-Gaussian null model.

The receiver report is `reports/cached_joint_kstat_pe_receiver_validation_report.md`. The main result MAT and PNG diagnostics are under `results/validation/cached_joint_kstat_pe_receiver/`.

## U=5 F=64 Conditional Channel Generator

- `validation/validate_u5_conditional_channel_generator_vertical.m`: set `U5_CONDITIONAL_MODE=smoke` for 32/16 samples or `full` for 128/64 plus 10,000 generated channels.
- `validation/refresh_u5_f64_properness_vertical.m`: refreshes the F=64 three-mode null distributions and prescribed central-interval decisions.
- `validation/generate_u5_f64_sample_bundle_vertical.m`: saves generated H, total H, physical CIR, delay axis, labels, and timing without PE calls.
- `reporting/plot_u5_f64_conditional_validation_vertical.m`: produces QQ, correlation-matrix, eigenvalue, and LFM diagnostics.

Reusable root interfaces are `estimate_conditional_channel_stats_vertical`, `sample_conditional_channel_vertical`, `build_channel_cir_vertical`, and `properness_null_test_vertical`. `build_physical_cir_vertical` remains only as the legacy common-time-shift wrapper. Full results are under `results/validation/u5_conditional_channel_f64/`; interpretation is in `reports/u5_conditional_channel_generator_f64_report.md`.

## Optimized Joint Builder and U=8 Node

- `validation/validate_streaming_joint_builder_u5_vertical.m`: checks the streaming/series F=64 builder against the original U=5 factors and cached-PE receiver statistics.
- `validation/audit_raw_pm_u8_aperture_vertical.m`: audits 100/256², 150/384², and 200/512² PM grids with 128 central-crop mapping realizations.
- `validation/validate_u8_conditional_channel_generator_vertical.m`: uses `U8_CONDITIONAL_MODE=smoke` for 16/8 samples or `full` for 128/128, U=8 properness, rank comparison, and 10,000 H+CIR samples.
- `validation/build_u5_u8_conditional_library_vertical.m`: creates and smoke-tests the exact-node U=5/U=8 library after both full models exist.

The reusable library interfaces are `build_conditional_channel_library_vertical` and `sample_conditional_channel_library_vertical`. They deliberately reject unsupported wind speeds and do not interpolate. The U=8 full 128/128 validation and two-node smoke test pass; see `reports/u8_joint_optimization_two_node_library_report.md` for measured acceptance and remaining high-wind memory limits.

## Two-Node Communication Validation

- `validation/validate_two_node_communication_vertical.m`: set `TWO_NODE_COMM_MODE=smoke` or `full`; compares fresh kdomain PE, joint cached PE, full-rank statistics, and 99.9% statistics at U=5/U=8.
- The full run uses 32 fresh channels per source, 8000 QPSK symbols per channel, Eb/N0 0:2:20 dB, channel-cluster bootstrap intervals, and a separate 128-pair full/low-rank test.
- `comm_main_vertical_psk` accepts optional external H(f) or h(t) through `COMM_EXTERNAL_CHANNEL_FILE`. Empty/unset preserves the original PE scenarios.

Reusable interfaces are `build_communication_taps_vertical`, `evaluate_mpsk_channel_ensemble_vertical`, and `sample_conditional_channel_rank_pair_vertical`. See `reports/two_node_statistical_channel_communication_validation_report.md`.

## PE Carrier-Phase Release Candidate Workflow

The formal workflow is destructive only in the narrow archival sense: the
preparation step moves registered phase-sensitive results to a timestamped
folder under `results/archive/pe_phase_reference_pre_rc/`. It never deletes
or overwrites them. Run the stages in this order:

```matlab
run('scripts/validation/prepare_pe_phase_reference_release_candidate_vertical.m')
run('scripts/validation/validate_pe_channel_phase_reference_vertical.m')
run('scripts/validation/validate_pe_phase_convention_uniform_vertical.m')

setenv('ADJOINT_PE_VALIDATION_MODE','full')
run('scripts/validation/validate_adjoint_pe_receiver_projection_vertical.m')
setenv('U5_CONDITIONAL_MODE','full')
run('scripts/validation/validate_u5_conditional_channel_generator_vertical.m')
setenv('U8_CONDITIONAL_MODE','full')
run('scripts/validation/validate_u8_conditional_channel_generator_vertical.m')
run('scripts/validation/build_u5_u8_conditional_library_vertical.m')

setenv('TWO_NODE_COMM_MODE','full')
run('scripts/validation/validate_two_node_communication_vertical.m')
run('scripts/validation/validate_public_channel_modes_vertical.m')
run('scripts/validation/validate_cached_pe_public_consistency_vertical.m')

setenv('PE_ATLAS_MODE','full')
run('scripts/reporting/generate_pe_propagation_atlas_vertical.m')
run('scripts/validation/finalize_pe_phase_reference_release_candidate_vertical.m')
```

`audit_phase_reference_artifacts_vertical.m` is the read-only inventory
entrypoint and may be run separately. The active run metadata is stored in
`results/validation/pe_phase_release_candidate/current_run.mat`. A stopped
run may resume only when its code fingerprint is unchanged. If a
fingerprint-covered source changes, rerun the preparation step: it archives
the partial run and creates a new `run_id`.

The formal adjoint configuration is PE `128^2` / PM `256^2`; the U=8
conditional node intentionally uses its separately audited PE `128^2` / PM
`384^2` aperture. F=9 uses 4096 realizations, F=64 uses 512, and these runs
are expensive. The atlas must be last because it rejects validation inputs
whose `run_id` does not match the active release run. The finalizer reports
only `PASS`, `FAIL`, or `INCOMPLETE` and writes
`reports/pe_phase_reference_release_candidate_report.md`.

The completed reference run is `phase_rc_20260722_174945` and is `PASS`.

## Exact Adjoint PE Receiver Projection

Before the adjoint suite, the receiver carrier-reference integration can be
checked independently:

```matlab
run('scripts/validation/validate_pe_channel_phase_reference_vertical.m')
```

This reduced-cost audit uses F=65 for the unaliased 4 ms delay test, exercises
the four public surface branches, checks `legacy_reduced`, compares cached and
adjoint receiver outputs in both reduced/direct-DSP form, validates dense/FFT
`C/P`, and tests schema-1 conditional-model migration. Outputs are under
`results/validation/pe_channel_phase_reference/`.

- `validation/validate_adjoint_pe_receiver_projection_vertical.m`: validates the exact discrete conjugate transpose of the cached uniform surface-to-receiver PE, receiver projection, PM-to-PE embedding, dense/FFT receiver statistics, realization statistics, performance, and public regressions.
- Set `ADJOINT_PE_VALIDATION_MODE=smoke` for the reduced run or `full` for the accepted F=9/4096 and F=64/512 validation.
- Reusable validation interfaces are `apply_forward_surface_to_receiver_vertical`, `apply_adjoint_receiver_to_surface_vertical`, `build_adjoint_receiver_projection_vertical`, `run_adjoint_receiver_projection_vertical`, and `contract_kstat_receiver_stats_vertical`.
- `reporting/plot_adjoint_pe_receiver_projection_validation_vertical.m` regenerates the adjoint error, covariance/pseudo-covariance, PDP, LFM/matched-filter, augmented-eigenvalue, and timing/memory figures from a saved validation result.
- Outputs are written to `results/validation/adjoint_pe_receiver_projection/`; measured results and the integration decision are in `reports/adjoint_pe_receiver_projection_feasibility_report.md`.

This v1 path is limited to uniform sound speed, CPU double, fixed grids and frequency axis, one nearest-grid receiver, no bubbles, and no Doppler. Its PE operator remains a validation path; receiver outputs now use the same central phase-reference layer as the public API. The default surface model is unchanged.

## PE Propagation Visual Atlas

- `reporting/generate_pe_propagation_atlas_vertical.m` is the report-only entrypoint. It assembles the Tx-to-surface and surface-to-Rx center slices, surface and receiver planes, physical boundary branches, exact-adjoint sensitivity, accepted receiver `C/P`, wind-node statistics, and a signal-free carrier-reconstruction animation.
- `reporting/plot_pe_propagation_atlas_vertical.m` renders the fixed 14-item atlas. Spatial panels use a shared amplitude reference and a `[-50,0] dB` scale; phase is masked below `-40 dB`.
- Run `set PE_ATLAS_MODE=smoke` before MATLAB for the PE 64² / PM 128² / F=9 check. With the variable unset, the formal run reads the accepted PE 128² / PM 256² / F=64 adjoint and conditional-model results.
- Outputs are written to `results/visualization/pe_propagation_atlas/`, including a reusable MAT file, source/normalization manifest, numerical-closure summary, 13 PNG figures, and one MP4.

The atlas distinguishes physical boundary models from computational acceleration paths. `q` is an exact receiver-sensitivity kernel, not a physical reverse-propagated pressure field. LFM is applied only after obtaining `H(f)` and is shown without noise, modulation, synchronization, or equalization.

## PE/Bellhop Flat-Surface Cross-Validation

- `validation/validate_pe_bellhop_flat_surface_vertical.m` generates and runs a standard Bellhop ASCII-arrivals case matched to the existing PE model in a uniform medium with a flat pressure-release surface.
- `validation/validate_pe_phase_convention_uniform_vertical.m` independently checks the PE reduced-envelope operator, longitudinal carrier sign, group delay, and validation-local FFT convention against one-step angular-spectrum propagation.
- `validation/validate_pe_bellhop_flat_surface_matrix_vertical.m` runs the authoritative 3/6/9 m analytic--PE--Bellhop timing/amplitude comparison with an open Bellhop fan, plus the C0--C4 PE grid/window/step convergence matrix.
- `reporting/generate_bellhop_flat_surface_visuals_vertical.m` reuses that saved matrix, runs Bellhop `R`/`C`/`I`, parses `.ray`/`.shd` locally, checks 5001/10001-beam convergence, and creates ray, shared-scale TL-field, and receiver-depth slice figures.
- Set `BELLHOP_EXE` when `bellhop.exe` is not already on the MATLAB or system path.
- The validator runs scalar direct-only/direct-plus-reflection regressions and a 65-frequency PE case, compares direct/single-surface arrival times and TL, reconstructs matched PDPs, and checks public PE invariants.
- Outputs are written to `results/validation/pe_bellhop_flat_surface/`; the Markdown report records exact parameters, formulas, thresholds, results, and limitations.

This stage deliberately excludes rough-surface scattering, bottom bounces, stochastic channels, and communication processing. Validators now consume the public `H_*_reduced_f` and `H_*_physical_f` fields instead of applying an independent hidden carrier convention.

Run the phase audit before the matrix validator; the latter refuses to run
unless the saved audit passed with carrier sign `+1`. The current matrix result
passes timing and public-interface invariants but intentionally retains failed
strict checks for cross-geometry source normalization and the doubled-window
C2 phase comparison. Reports are under
`results/validation/pe_phase_convention_uniform/` and
`results/validation/pe_bellhop_flat_surface_matrix/`.

Visualization outputs, including the original Bellhop files, CSV tables, MAT
result, three PNG figures, and Markdown report, are under
`results/visualization/bellhop_flat_surface/`. The display checks pass without
changing the matrix `passed=false` status or any PE/communication source file.

Use environment variables already supported by individual scripts to reduce grid size, seed count, or output file names for quick checks.

## Li et al. (2009) Explicit Rough-Surface Validation

```matlab
run('scripts/validation/validate_li2009_explicit_surface_vertical.m')
```

This independent pure-acoustic workflow uses one shared raw-PM realization
across the full 12 kHz, 6 ms CW frequency synthesis, applies only the explicit
pressure-release Kirchhoff `2*k*eta` screen, and projects the reflected field
through a monostatic bottom--surface--bottom PE path. It excludes noise,
electronics, SSA, kstat, modulation, and BER. Outputs are under
`results/validation/li2009_explicit_surface/`; interpretation and unresolved
transverse-grid sensitivity are documented in
`reports/li2009_explicit_surface_validation_report.md`.
