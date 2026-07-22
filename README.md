# Vertical Underwater Acoustic Channel Project

This repository contains a MATLAB vertical underwater acoustic channel and MPSK communication model. The core path is PE/WAPE-style propagation to obtain channel frequency responses, followed by optional surface reflection/scattering diagnostics and communication metrics.

## Core Files

- `vertical_channel_model.m`: public channel API, `output = vertical_channel_model(paramsV)`.
- `vertical_wape_propagator.m`: upward-marching propagation core and frequency loop.
- `pm_surface_boundary_model.m`: rough sea-surface boundary and statistical scattering model.
- `comm_main_vertical_psk.m`: end-to-end MPSK communication demo using `H_f`.
- `modem_psk.m`, `noise_inject_vertical.m`: communication helpers.
- `bubble_effective_medium.m`, `bubble_environment_vertical.m`, `bubble_hall_spectrum.m`: bubble-layer physics helpers.
- `build_surface_empirical_channel_model_vertical.m`, `sample_surface_empirical_channel_vertical.m`: empirical channel-model helpers.
- `main_vertical.m`: channel-only demo.
- `explain_main_vertical.m`: compatibility alias for the channel-only demo.

Longer implementation notes remain in `PROJECT_CONTEXT.md`, `vertical_comm_guide.md`, `SSA.md`, `BUBBLE_EXTENSION_SPEC.md`, and `README_CODE_STRUCTURE_AND_OUTPUTS.md`.

## Directory Layout

- `scripts/validation/`: reduced-grid regression and diagnostic checks.
- `scripts/comparisons/`: model-to-model comparison studies.
- `scripts/experiments/`: sweeps, Monte Carlo runs, and calibration studies.
- `scripts/reporting/`: figure/report generation scripts.
- `results/`: local generated artifacts, grouped by purpose.
- `reports/`: Markdown research notes and narrative reports.
- `old/`: legacy reference code kept out of the active execution path.

Scripts under `scripts/` call `scripts/bootstrap_project.m`, which adds the project root to the MATLAB path and routes relative output files into the matching `results/` subfolder.

## Quick Runs

From MATLAB:

```matlab
cd('E:/MISC/CARPE3D_matlab/Explain')
explain_main_vertical
comm_main_vertical_psk
run('scripts/validation/validate_surface_wavefield_visualization_vertical.m')
```

Flat-surface PE/Bellhop stage-1 cross-validation (requires the external
Acoustics Toolbox `bellhop.exe`):

```matlab
setenv('BELLHOP_EXE', 'C:/path/to/bellhop.exe')
addpath('scripts/validation')
validate_pe_bellhop_flat_surface_vertical
```

This deterministic validator uses the existing PE API, a uniform SSP, and a
flat pressure-release surface. It writes the Bellhop environment, arrivals,
PDP/TL figure, comparison table, MAT result, and Markdown report under
`results/validation/pe_bellhop_flat_surface/`.

The stricter follow-on audit keeps the main PE code unchanged and writes to
separate result directories:

```matlab
validate_pe_phase_convention_uniform_vertical
validate_pe_bellhop_flat_surface_matrix_vertical
```

The first entry verifies the reduced-envelope carrier sign against an
independent angular-spectrum reference. The second uses an open Bellhop angle
fan at 3/6/9 m offsets and runs the C0--C4 PE grid/window/step convergence
matrix. Its current timing and algebra checks pass, while strict cross-geometry
amplitude and doubled-window phase checks remain flagged; see the generated
matrix report rather than treating the earlier single-geometry timing as an
independent delay validation.

Bellhop flat-surface ray and TL displays can be generated independently from
the saved 3/6/9 m matrix result:

```matlab
setenv('BELLHOP_EXE', 'C:/path/to/bellhop.exe')
addpath('scripts/reporting')
generate_bellhop_flat_surface_visuals_vertical
```

This entrypoint runs Bellhop `R`, coherent `C`, and incoherent `I` modes,
checks 5001/10001-beam convergence, and writes ray geometry, shared-scale TL
fields, and a receiver-depth TL slice under
`results/visualization/bellhop_flat_surface/`. It reuses saved PE receiver
points and does not construct a PE range-depth field or alter the strict
matrix conclusion.

The complete consolidated record of the single-geometry smoke case, phase
audit, 3/6/9 m strict matrix, Bellhop R/C/I visualization, regressions, and
remaining limitations is in
`reports/pe_bellhop_flat_surface_cross_validation_complete_report.md`.

For quick validation, prefer reduced grids (`nx=128` or `nx=256`, `ny=128` or `ny=256`) and `show_figures=false`.

## Surface Boundary Models

`surface_boundary_model` currently supports:

- `kirchhoff_spatial`: default explicit PM sea-surface realization and spatial phase screen.
- `kirchhoff_kdomain`: FFT k-domain interface for the same explicit Kirchhoff phase screen.
- `kirchhoff_kstat`: Kirchhoff statistical phase-screen branch. It does not generate a concrete `eta(x,y)`; it builds `C_eta`, `S_deltaG`, coherent reflection, incoherent scatter power, and optional random reflected spectra from the PM height spectrum.
- `ssa_stat_kernel`: SSA-oriented statistical research branch used as a weak-roughness/small-angle reference.

Roughness amplitude scaling is controlled by `surface_roughness_scale_mode`:

- `target_hs` (default): scale the PM surface spectrum or realization to `sea_hs_target`.
- `raw_pm`: use the PM spectrum amplitude implied directly by `sea_wind_speed`; metadata records `sigma_eta_raw_m`, `Hs_raw_m`, `Hs_target_m`, and `scale_factor`.

Quick kstat validation:

```matlab
run('scripts/validation/validate_kirchhoff_kstat_vertical.m')
```

Raw-PM wind-driven comparison between the explicit `kirchhoff_kdomain` and statistical `kirchhoff_kstat` branches:

```matlab
run('scripts/comparisons/compare_kirchhoff_kdomain_kstat_wind_vertical.m')
```

It writes the MAT/CSV summary and figures under `results/comparisons/`.
See `vertical_comm_guide.md` for the discrete raw-PM variance convention used to align `kirchhoff_kdomain` and `kirchhoff_kstat`.

## Output Policy

Generated `.mat`, `.png`, `.csv`, `.mp4`, `.log`, and similar run artifacts belong under `results/`. Large binary outputs are intentionally ignored by Git; the committed source of truth is the MATLAB code and Markdown documentation.
