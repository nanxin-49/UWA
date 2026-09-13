# Bellhop internal-wall visualization archive manifest

Date: 2026-09-09

The visualization-only scripts and generated products from the Bellhop
internal-wall work have been moved to the recoverable archive
`cash/bellhop_internal_wall_visualization_archive_20260909/`. They are kept for
reproducibility and historical inspection, but are not active project
entrypoints or current validation inputs.

Archived scripts:

- `scripts/reporting/generate_bellhop_internal_wall_reflection_visuals.py`
- `scripts/reporting/generate_bellhop_internal_wall_visual_enhancement.py`
- `scripts/reporting/generate_bellhop_internal_wall_tl_diagnostic.py`
- `scripts/reporting/generate_bellhop_internal_wall_tl_2d.py`
- `scripts/validation/run_bellhop_internal_wall_tl_grid.py`
- `scripts/validation/validate_bellhop_internal_wall_visual_enhancement.m`

Archived result trees:

- `results/visualization/bellhop_internal_wall_reflection/`
- `results/visualization/bellhop_internal_wall_visual_enhancement/`
- `results/validation/bellhop_internal_wall_tl_grid/`
- `results/validation/bellhop_internal_wall_visual_enhancement/`

The authoritative Bellhop implementation and numerical validation conclusions
remain in `reports/bellhop_internal_wall_implementation_report.md`. This
manifest intentionally does not duplicate archived files or enumerate their
contents.

## 2026-09-10 close-out

The current visualization consolidation keeps only the three final products under
`results/visualization/bellhop_internal_wall_visuals/`:

- `internal_wall_ray_trajectories.png`
- `internal_wall_incident_tl_2d.png`
- `internal_wall_tl_2d.png`

The active renderer and the dense input `dense_wall_fields.mat` remain outside the
archive. Superseded visualization reports, duplicate dense MAT files, and raw case
sidecars were moved to `cash/bellhop_internal_wall_visualization_archive_20260910/`.
The consolidated authoritative report is
`reports/bellhop_internal_wall_visualization_report.md`.

## 2026-09-13 full-fan consolidation

The active visualization was rebaselined to the accepted two-period full-fan
case `r=100-2.4 sin((2pi/30)z)`, with 5001 beams over `[-30,30] deg` and the
receiver eigenray explicitly highlighted. Both TL figures now display the
complete sampled wall support `z=[-75,75] m`; areas outside the computed SHD
receiver grid remain background/NaN and are not extrapolated.

The following superseded or intermediate material was moved to the recoverable
archive `cash/bellhop_internal_wall_visualization_consolidation_20260912/`:

- the previous consolidated, standalone-large-sinusoid and staged full-fan
  visualization reports;
- the solver-free standalone wall-profile plotting entrypoint;
- duplicate visualization trees for the old baseline, standalone wall preview
  and staged full-fan output;
- regenerable Bellhop case sidecars from the full-fan validation.

The active retained products are:

- `reports/bellhop_internal_wall_visualization_report.md`;
- `scripts/reporting/generate_bellhop_internal_wall_visuals_vertical.m`;
- `scripts/validation/validate_bellhop_internal_large_sinusoidal_wall_visual_case.m`;
- the three PNGs and selected-ray CSV under
  `results/visualization/bellhop_internal_wall_visuals/`;
- the parsed full-fan MAT files and geometry/profile CSV files under
  `results/validation/bellhop_internal_wall_two_period_full_fan/`.
