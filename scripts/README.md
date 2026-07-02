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

Use environment variables already supported by individual scripts to reduce grid size, seed count, or output file names for quick checks.
