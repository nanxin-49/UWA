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
