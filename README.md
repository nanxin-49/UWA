# Vertical Underwater Acoustic Channel

MATLAB project for a seabed-to-near-surface vertical underwater acoustic
channel. The active path uses PE/WAPE-style propagation, optional rough
sea-surface reflection/scattering, receiver-side statistical models, and an
MPSK communication consumer.

## Start Here

The public channel API is:

```matlab
output = vertical_channel_model(paramsV);
```

Minimal entrypoints:

```matlab
cd('E:/MISC/CARPE3D_matlab/Explain')
explain_main_vertical       % channel demonstration
comm_main_vertical_psk      % end-to-end MPSK demonstration
```

For quick checks, use reduced even grids such as `64^2` or `128^2`, set
`show_figures=false`, and disable the expensive features that are not under
test.

## Current Channel Semantics

The public default is `paramsV.channel_phase_reference='direct_dsp'`.
`H_direct_f`, `H_reflect_f`, and `H_f` therefore share one receiver phase
reference and are ready for the MATLAB-IFFT communication path:

```matlab
H_f = H_direct_f + H_reflect_f;
```

The raw PE envelopes remain available as `H_*_reduced_f`; absolute physical
phasors under `exp(-1i*omega*t)` are available as `H_*_physical_f`.
`phase_reference_meta` records the geometry, nominal delays, convention, and
frequency-dependent phase factors. `legacy_reduced` remains available only
for regression and migration.

For the standard `z_tx=100 m`, `z_rx=3 m`, `c0=1500 m/s` geometry, the
nominal direct/reflected longitudinal spans are 97/103 m and their reference
delay difference is 4 ms. This deterministic reference is not a per-sample
PDP peak alignment.

## Core Code

- `vertical_channel_model.m`: public configuration and output boundary.
- `vertical_wape_propagator.m`: direct and reflected PE/WAPE propagation.
- `pm_surface_boundary_model.m`: explicit Kirchhoff, joint K-stat, and SSA
  research surface boundaries.
- `apply_pe_channel_phase_reference_vertical.m`: central receiver phase
  reference conversion.
- `build_channel_cir_vertical.m`: convention-explicit direct-DSP or physical
  CIR reconstruction with only a common time-origin shift.
- `comm_main_vertical_psk.m`: reference communication consumer of `H_f`.
- `build_cached_joint_kstat_pe_executor_vertical.m`: fixed uniform cached PE
  validation path.
- `build_adjoint_receiver_projection_vertical.m`: exact discrete-adjoint
  single-receiver projection prototype.
- `contract_kstat_receiver_stats_vertical.m`: dense/FFT receiver `C/P`
  contraction on the PM grid.
- `estimate_conditional_channel_stats_vertical.m` and
  `sample_conditional_channel_vertical.m`: receiver statistical generator.

The public default surface model remains `kirchhoff_spatial`. Adjoint and
cached routines are validation/research interfaces; they do not replace the
public PE propagator for unsupported environments.

## Repository Layout

- `scripts/validation/`: deterministic regression and statistical checks.
- `scripts/comparisons/`: model-to-model studies.
- `scripts/experiments/`: sweeps and Monte Carlo experiments.
- `scripts/reporting/`: plots, reports, and propagation atlas generation.
- `reports/`: detailed validation evidence and feasibility conclusions.
- `results/`: generated MAT, text, table, figure, animation, and external-tool
  artifacts; large generated files are normally ignored by Git.
- `old/`: legacy reference code outside the active execution path.

## Documentation Map

| Document | Audience | Purpose |
|---|---|---|
| `README.md` | New repository reader | Quick entrypoint, active code, and navigation |
| `vertical_comm_guide.md` | Researchers and developers | Current physics, phase, statistics, interfaces, and interpretation |
| `PROJECT_CONTEXT.md` | Future coding models and maintainers | Authoritative current index plus append-only project log |
| `scripts/README.md` | Validation operator | Script catalogue, prerequisites, cost, and outputs |
| `reports/` | Technical reviewer | Configuration-specific evidence and decisions |

Do not infer current behavior from an old dated report without checking the
current-state index in `PROJECT_CONTEXT.md` and the implementation.

## Key Validation

Carrier-phase release candidate `phase_rc_20260722_174945` completed with
decision **PASS**. It includes the formal F=65 delay/sign audit, full F=9 and
F=64 adjoint/statistical suites, U=5/U=8 conditional models, two-node
communication, public/cached regression, and a same-run propagation atlas.
See `reports/pe_phase_reference_release_candidate_report.md` for the gate
table and `results/visualization/pe_propagation_atlas/` for the figures.

Run the carrier-phase and branch integration audit:

```matlab
run('scripts/validation/validate_pe_channel_phase_reference_vertical.m')
```

Run the exact-adjoint reduced smoke test:

```matlab
setenv('ADJOINT_PE_VALIDATION_MODE','smoke')
run('scripts/validation/validate_adjoint_pe_receiver_projection_vertical.m')
```

Additional commands, expected runtimes, external Bellhop prerequisites, and
output locations are maintained in `scripts/README.md`.

## Active Limits

- The exact-adjoint/cached v1 path is uniform sound speed, CPU double, fixed
  grids/frequencies, one nearest-grid receiver, and no bubbles or Doppler.
- `kirchhoff_kstat` is a near-vertical Gaussian/Kirchhoff statistical phase
  screen, not a complete rough-surface scattering theory.
- Conditional wind libraries use validated discrete nodes; no silent wind
  interpolation is performed.
- Unknown external `H(f)` without project phase metadata is assumed to be
  DSP-ready and is not silently rotated.
- The F=9 grid `4:0.5:8 kHz` cannot resolve the standard 4 ms reference delay:
  its unambiguous delay is 2 ms. Use the F=65, 62.5 Hz audit grid for this test.
