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
