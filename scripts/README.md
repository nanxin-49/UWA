# 脚本与验证运行索引

脚本按职责分组，入口应能从任意 MATLAB 当前目录运行。推荐先执行根目录 `setup_vertical_project`；它只初始化路径，不改变当前工作目录。仍使用 `scripts/bootstrap_project.m` 的旧脚本会先委托给统一 setup，再按既有规则切换到 `results/<category>`。新入口应优先使用绝对项目根和明确输出目录，不依赖当前目录。

## 目录

- `validation/`: smoke tests, invariants, communication-chain checks, and disabled-path regressions.
- `comparisons/`: surface-model, bubble-model, and coherent/incoherent comparison studies.
- `experiments/`: Monte Carlo runs, parameter sweeps, and calibration scripts.
- `reporting/`: figure generation, summary tables, visualizations, and report builders.
- `validation/support/`: Bellhop/Weyl/Li2009、properness、发布元数据和 artifact 检查等验证专用辅助函数；不属于公共 API。

## 运行方式

```matlab
setup_vertical_project
run('scripts/validation/validate_surface_wavefield_visualization_vertical.m')
run('scripts/comparisons/compare_specular_incoherent_surface_reflection_vertical.m')
run('scripts/experiments/sweep_monte_carlo_surface_channel_vertical.m')
run('scripts/reporting/plot_c35_core_heatmaps_vertical.m')
```

公共可复用实现已移到 `src/`。下文脚本名仍按 `validation/`、`reporting/` 等相对于 `scripts/` 的短路径描述。

## 2k 相位与海谱诊断

```matlab
run('scripts/validation/validate_2k_phase_approx.m')
run('scripts/validation/validate_2k_ocean_spectra.m')
```

- PNG/CSV：`results/validation/ssa_2k_phase/`
- 报告：`reports/validation_2k_phase_report.md`、`reports/validation_2k_ocean_spectra_report.md`

这两项是相位系数近似和最低阶 SSA/海谱趋势诊断，不是完整 PE、KStat、真实海洋散射或高阶 SSA 的验收。

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

Reusable interfaces are initialized by `setup_vertical_project` and implemented under `src/`（`properness_null_test_vertical` is validation support）: `estimate_conditional_channel_stats_vertical`, `sample_conditional_channel_vertical`, `build_channel_cir_vertical`, and `properness_null_test_vertical`. `build_physical_cir_vertical` remains only as the legacy common-time-shift wrapper. Full results are under `results/validation/u5_conditional_channel_f64/`; interpretation is in `reports/u5_conditional_channel_generator_f64_report.md`.

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

`validation/audit_phase_reference_artifacts_vertical.m` is the read-only inventory
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

### Unfolded Gaussian-source comparison

- `validation/validate_pe_bellhop_unfolded_flat_gaussian_vertical.m` is the
  current validation-only entrypoint for the accepted unfolded-coordinate
  construction.  It maps the 100 m to 3 m vertical path to Bellhop ranges
  97 m (direct) and 103 m (image/reflected), then applies the flat
  pressure-release factor `-1` at the image receiver.
- `validation/support/write_bellhop_unfolded_gaussian_env_vertical.m` writes
  the per-frequency `.env` and `.sbp`; the `.sbp` pattern is derived from the
  production Gaussian angular spectrum and is not a pointwise receiver fit.
- `validation/support/run_bellhop_unfolded_gaussian_vertical.m` runs Bellhop
  and reads the standard `.shd`/`.arr` outputs.  The primary acceptance uses
  normalized reflected/direct responses; absolute source pressure is an
  independent Weyl-reference diagnostic.
- Outputs are written to
  `results/validation/pe_bellhop_unfolded_flat_gaussian/` and the report is
  `reports/pe_bellhop_unfolded_flat_gaussian_report.md`.

### Reflection-free four-level audit

- `validation/validate_pe_as_freefield_vertical.m` compares production
  multi-step marching with an independently coded one-step exact discrete
  angular-spectrum propagator. The sponge is exactly off for this hard gate.
- `validation/validate_bellhop_freefield_normalization_vertical.m` reproduces
  the matched-halfspace construction of Bellhop's official free-space point
  source example. It measures `|p|R`, spatial phase sign, the constant source
  phase, beam/step convergence, and zero-bounce arrival delay before applying
  the single global `1/(4*pi)` Green-function conversion.
- `validation/validate_pe_bellhop_freefield_vertical.m` is the orchestration
  entrypoint. It keeps explicit initial-plane and physical-source phase
  references, compares PE/Bellhop/analytic complex fields, saves arrival and
  convergence CSV files, and does not relax a failed tolerance.
- `validation/support/virtual_point_source_initial_field_vertical.m` is enabled only through the
  additive `source_mode='custom_field_fn'` validation hook. Public defaults
  remain `source_mode='gaussian'`.

```matlab
setenv('BELLHOP_EXE','E:/stable/path/to/bellhop.exe')
addpath('scripts/validation')
validation = validate_pe_bellhop_freefield_vertical();
```

Artifacts are under `results/validation/pe_bellhop_freefield/formal/`; the
summary is `reports/pe_bellhop_freefield_validation_report.md`.

### Point-source error-budget follow-up

- `validation/validate_pe_point_source_error_budget_vertical.m` runs the
  no-sponge fixed-`dx` window series, compares spatial truncation with a
  discrete spectral-cell Weyl initializer, validates a separate continuous
  Weyl integral, and then scans the complete window/thickness/strength matrix.
  Its CSV explicitly reports `dA=20log10(|H_sponge|/|H_no_sponge|)` and
  `dTL=-dA`, plus center (`rho<=2 m`), ratio-defined edge-band, and total
  terminal-plane energies.
- `validation/support/weyl_point_source_reference_vertical.m` removes the grazing `1/kz`
  singularity by separate propagating/evanescent substitutions; it is the
  reliable free-space reference.
- `validation/support/weyl_point_source_initial_field_vertical.m` is the discrete FFT initializer
  under test. Its failure to converge within the public 100 m window is
  reported, not hidden.
- `validation/rerun_pe_bellhop_point_source_postbudget_vertical.m` performs the
  final representative Bellhop rerun only after the continuous Weyl gate.

Outputs are in `results/validation/pe_point_source_error_budget/` and
`results/validation/pe_bellhop_freefield/post_error_budget/`.

### Production Gaussian window/sponge audit

- `validation/validate_gaussian_window_convergence_vertical.m` calls the
  production Gaussian path, performs the independent PE--AS hard gate, and
  establishes a no-sponge reference using fixed production sampling.
- `validation/validate_gaussian_sponge_vertical.m` evaluates the complete
  50/80/128 m x 6-ratio x 7-strength matrix with explicit `dA`/`dTL`, center,
  edge, total-energy, and edge-to-center diagnostic metrics.
- `validation/validate_gaussian_sponge_wideband_vertical.m` compares the
  3--5 kHz physical `H(f)`, unwrapped phase, and group delay for large
  no-sponge, production/default, and recommended cases.
- `reporting/generate_gaussian_sponge_validation_figures.m` emits the 15
  requested figures plus an optional noiseless LFM diagnostic and the final
  Q1--Q7 report.

Outputs are under `results/validation/pe_gaussian_window_sponge/`; the report
is `reports/pe_gaussian_window_sponge_validation_report.md`. Extended sponge
ratios and windows are validation-only opt-ins; normal public limits and
defaults are unchanged.

### Full reflected PE -> surface -> PE window audit

- `validation/validate_reflected_chain_window_convergence_vertical.m` creates
  one maximum-domain PM surface, crops it without per-window Hs rescaling,
  runs the adaptive 4 kHz no-sponge window sequence, and saves incident,
  reflected, receiver, center, and edge diagnostics.
- `validation/validate_reflected_chain_wideband_vertical.m` runs the gated
  3--5 kHz/33-point reflected and total-channel comparison. Per-case MAT
  checkpoints make the production-size run resumable.
- `validation/validate_reflected_chain_diagnostics_regression_vertical.m`
  checks seeded/override equality, diagnostics transparency, receiver-center
  consistency, the phase-screen invariant, and channel closure.
- `reporting/generate_reflected_chain_window_validation_figures.m` generates
  the single-frequency field and convergence figures; the wideband validator
  adds reflected and total-channel response figures.

Formal outputs are under
`results/validation/pe_reflected_chain_window/`; the reviewable conclusion is
`reports/pe_reflected_chain_window_validation_report.md`. The validators use
additive validation-only options and do not alter production defaults.

### Random-surface reflected-chain window robustness

- `validation/validate_random_surface_window_robustness_vertical.m` generates
  each 256 m master surface once, takes unscaled 160/192 m central crops, and
  runs the 4 kHz three-Hs-by-five-seed matrix with resumable per-case files.
- `reporting/generate_random_surface_window_robustness_figures.m` creates the
  gate scatter, edge-versus-response, pass-rate, and stage-energy figures.
- `validation/validate_random_surface_window_robustness_wideband_vertical.m`
  is a checkpointed 3--5 kHz follow-up for selected worst cases. The current
  formal run was stopped before all candidate--256 m pairs completed; do not
  use partial checkpoints as group-delay qualification evidence.

Formal 4 kHz outputs are under
`results/validation/pe_random_surface_window_robustness/`; the report is
`reports/pe_random_surface_window_robustness_report.md`. The result is
192.1875 m/no-sponge `15/15` strict passes, while 160.15625 m has `15/15`
response/center passes but `0/15` complete edge passes. No production default
is changed by these validators.

- `validation/validate_pe_bellhop_flat_surface_current_vertical.m` is the
  current formal entrypoint. It requires a stable absolute, non-Temp
  `BELLHOP_EXE`, records the binary/source SHA-256 values, runs the independent
  phase audit, 3/6/9 m paths, small-offset limit, source-aware diagnostic,
  sampling/aperture/sponge matrix, public regressions, Bellhop R/C/I fields,
  and the numbered ten-figure atlas.
- `validation/validate_pe_bellhop_flat_surface_vertical.m` generates and runs a standard Bellhop ASCII-arrivals case matched to the existing PE model in a uniform medium with a flat pressure-release surface.
- `validation/validate_pe_phase_convention_uniform_vertical.m` independently checks the PE reduced-envelope operator, longitudinal carrier sign, group delay, and validation-local FFT convention against one-step angular-spectrum propagation.
- `validation/validate_pe_bellhop_flat_surface_matrix_vertical.m` runs the authoritative 3/6/9 m analytic--PE--Bellhop timing/amplitude comparison with an open Bellhop fan, plus the C0--C4 PE grid/window/step convergence matrix.
- `reporting/generate_bellhop_flat_surface_visuals_vertical.m` reuses that saved matrix, runs Bellhop `R`/`C`/`I`, parses `.ray`/`.shd` locally, checks 5001/10001-beam convergence, and creates ray, shared-scale TL-field, and receiver-depth slice figures.
- Formal mode requires `BELLHOP_EXE` even if another Bellhop is on the MATLAB
  or system path. Temporary paths and binary-hash changes are rejected.
- The validator runs scalar direct-only/direct-plus-reflection regressions and a 65-frequency PE case, compares direct/single-surface arrival times and TL, reconstructs matched PDPs, and checks public PE invariants.
- Outputs are written to `results/validation/pe_bellhop_flat_surface/`; the Markdown report records exact parameters, formulas, thresholds, results, and limitations.

This stage deliberately excludes rough-surface scattering, bottom bounces, stochastic channels, and communication processing. Validators now consume the public `H_*_reduced_f` and `H_*_physical_f` fields instead of applying an independent hidden carrier convention.

Current formal command:

```matlab
setenv('BELLHOP_EXE','E:/stable/path/to/bellhop.exe')
addpath('scripts/validation')
validation = validate_pe_bellhop_flat_surface_current_vertical();
```

The current run `bellhop_current_20260723_rc5` is `FAIL_CORE` with
`amplitude_status=OPEN`: all phase/path/delay/beam/public checks pass, while
the no-sponge aperture-only fields fail the `-40 dB` edge-validity prerequisite.
Outputs are versioned under
`results/validation/pe_bellhop_flat_surface_current/<run_id>/` and
`results/visualization/pe_bellhop_flat_surface_current/<run_id>/`.

The formal atlas contains 10 numbered PNGs plus MAT, CSV, manifest, and text
summary. Older three-figure outputs remain historical evidence and do not
override the current run ID.

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
