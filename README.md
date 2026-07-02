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

## Output Policy

Generated `.mat`, `.png`, `.csv`, `.mp4`, `.log`, and similar run artifacts belong under `results/`. Large binary outputs are intentionally ignored by Git; the committed source of truth is the MATLAB code and Markdown documentation.
