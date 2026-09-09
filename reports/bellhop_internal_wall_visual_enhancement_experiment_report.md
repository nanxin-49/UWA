# Bellhop internal-wall visual-enhancement experiment report

## 1. Status and scope

**Status: PASS (display-only extension).** The stored enhancement validation
summary passes all of its hard-check fields for both cases. This experiment is
an extended diagnostic for making the reflection and chart mapping visible; it
does not replace the existing weak-wall regression, and it does not change the
authoritative conclusions in
`reports/bellhop_internal_wall_implementation_report.md`.

No Bellhop run was repeated during this close-out. The figures are read-only
post-processing of the already complete `.iwdiag`/`.iw3` sidecars and the
stored validation summary.

## 2. Configuration and provenance

The validation entrypoint is
`scripts/validation/validate_bellhop_internal_wall_visual_enhancement.m`.
Both cases use `f=4000 Hz`, `c=1500 m/s`, `step=0.05 m`, `5001` beams, source
depth `0 m`, nominal wall line `R0=100 m`, physical reference receiver
`r=97 m`, mapped receiver `r'=103 m`, and angle fan `[-15,15] deg`.

| case | wall geometry | profile samples | stored sidecars |
|---|---|---:|---|
| tilted | `r=100+0.05 z` | analytic display profile; hit diagnostics from 5001 beams | `.iwdiag`, `.iw2` |
| sinusoidal | `r=100-2 sin(0.04 z)` over `z∈[-32,32] m` | 161 | `.iwdiag`, `.iw3` |

The build and executable provenance is recorded in
`results/validation/bellhop_internal_wall_visual_enhancement/build_manifest.json`.
The executable is a separate validation binary; official Bellhop 2020 sources
and the production PE/Reflect2D/InfluenceGeoHatCart paths were not modified.

## 3. Hard checks and geometry safety

The values below are read from
`results/validation/bellhop_internal_wall_visual_enhancement/visual_enhancement_summary.csv`.
The `step_m` summary column is retained as `NaN` by the existing validator's
case-summary schema; the requested `0.05 m` setting is independently fixed by
the case filenames and validator options above.

| check | tilted | sinusoidal |
|---|---:|---:|
| all finite | PASS | PASS |
| wall intersection residual (m) | `3.53850e-12` | `3.52723e-12` |
| unit frame error | `1.11022e-16` | `2.22045e-16` |
| specular direction error | `6.68443e-16` | `6.68443e-16` |
| rotation direction error | `0` | `0` |
| one-`pi` phase check | PASS | PASS |
| amplitude jump | `0` | `0` |
| dynamic p/q state check | PASS | PASS |
| `max |kappa|` (1/m) | `0` | `2.84018e-3` |
| curvature check | PASS | PASS |
| path-length error (m) | `8.75389e-12` | `8.81073e-12` |
| travel-time error (s) | `5.83592e-15` | `5.87382e-15` |
| imaginary travel time (s) | `0` | `0` |
| minimum post-map range increment (m) | `4.67645e-2` | `4.72099e-2` |
| positive post-map range | PASS | PASS |

The tilted wall has zero curvature as expected. The sinusoidal case supplies
the intended nonzero-curvature display diagnostic: `max |kappa| =
2.84018e-3 1/m` and `max |p_reflect kick| = 0.626241`. These values are
diagnostics of the stored extension and are not a replacement acceptance gate
for the weak case.

## 4. Generated visual products

The read-only reporting entrypoint is
`scripts/reporting/generate_bellhop_internal_wall_visual_enhancement.py`.
It writes to
`results/visualization/bellhop_internal_wall_visual_enhancement/`:

- `01_tilted_wall_end_to_end.png` and `02_sinusoidal_wall_end_to_end.png` are
  independent case figures. Each has a complete Tx-to-wall-to-receiver main
  panel and a wall-neighborhood panel containing all eight selected rays (the
  seven-ray `[-15,-10,-5,0,5,10,15]` fan plus the target beam when it is
  distinct). The physical reconstruction main panels keep `(r,z)` coordinates;
  the wall-neighborhood panel re-orders the display axes to `(z,r)` (horizontal
  `z`, vertical range `r`) so the short wall-normal range window and full fan
  remain readable. This is a plotting convention only; screen slopes in the
  auxiliary panels are not angle measurements.
- The main and zoom panels use smaller receiver, mapped-receiver, and wall-hit
  markers. Small arrows show propagation direction, so the true `rot = -ref`
  180-degree direction relation is not mistaken for a plotting error. Gray
  dashed segments are native physical reflected branches; colored solid
  segments are proper-rotation branches; purple dotted segments are explicitly
  labeled coordinate-map connectors and are not propagation.
- `01_visual_enhancement_end_to_end.png` is retained as a compatibility
  two-panel overlay using the same selected rays and marker/arrow conventions.
- `02_visual_enhancement_delta_r_vs_z.png` remains the true-scale
  `Delta-r = r_wall(z)-R0` versus `z` geometry plot.
- `ray_direction_audit.csv` records every plotted ray's alpha, `norm(rot+ref)`,
  incident/native/rotated displacement-direction residuals, receiver-plane
  residuals, native/rotated segment lengths, absolute/relative length mismatch,
  and PASS/FAIL status. `visualization_manifest.csv` indexes the independent
  figures, compatibility output, geometry plot, and audit CSV.

### 4.1 Plotted-ray direction audit

For every selected ray, the postprocessor reconstructs the native endpoint from
`hit + ((97-hit_r)/ref_ur) ref` and the rotated endpoint from
`mapped_hit + ((103-mapped_hit_r)/rot_ur) rot`. It then checks unit displacement
vectors against the sidecar directions, `norm(rot+ref)`, receiver-plane
endpoint residuals, and equality of native/reflected and rotated branch
lengths. The audit thresholds are `1e-10` for direction norms and endpoint
residuals and `1e-10 m` for the length mismatch.

| case | plotted rays | max `norm(rot+ref)` | max endpoint residual (m) | max length mismatch (m) | audit |
|---|---:|---:|---:|---:|---|
| tilted | 8 | `0` | `1.82609e-12` | `0` | **PASS** |
| sinusoidal | 8 | `0` | `1.84031e-12` | `0` | **PASS** |

Thus the gray native branch and colored proper-rotation branch are parallel as
an unoriented line but oppositely directed as propagation vectors, exactly as
required by the single proper `pi` rotation. The arrows and CSV audit preserve
that distinction.

### 4.2 Physical inverse-rotation reconstruction

The final physical-coordinate reconstruction is written separately as
`03_tilted_wall_physical_reconstruction.png` and
`04_sinusoidal_wall_physical_reconstruction.png`. For every selected ray, the
stored mapped forward segment is transformed with
`T^{-1}(r',z')=(2R0-r',-z')` and plotted in the original physical `(r,z)` chart
over the native `Reflect2D` backward segment. The gray dashed reference is
drawn first and the colored solid inverse-mapped branch second, so coincident
segments are directly visible; both arrow sets point from the wall toward the
physical receiver at `r=97 m`. The mapped positive-range branch is shown only
in the separate computational-coordinate panel and is not treated as a
physical propagation segment.

`inverse_rotation_ray_audit.csv` contains the per-ray direction, hit,
endpoint, receiver-plane, and path-length checks, while
`inverse_rotation_visualization_manifest.csv` indexes the new products. All
16 selected rays pass the `1e-10` direction/position/length thresholds; the
stored audit has zero failed rows and zero length mismatch. The wall-neighborhood
panel is explicitly annotated as a `(z,r)` display-axis reordering. The mapped
forward branch is not drawn as a third physical subplot; its geometry remains
available in the stored sidecars and is used by the inverse-rotation audit.

### 4.3 Stored-data transmission-loss diagnostic

The read-only postprocessor
`scripts/reporting/generate_bellhop_internal_wall_tl_diagnostic.py` writes
`05_tl_diagnostic.png` and `tl_diagnostic.csv`. The upper panels use the
stored direct pressure as the local reference and map the internal-wall
receiver columns `r'=103,102 m` back to physical-equivalent `r=97,98 m`.
The resulting wall-minus-direct relative TL is `0.52168/0.34751 dB` for
tilted and `0.52416/0.34799 dB` for sinusoidal at physical-equivalent
`r=97/98 m`, respectively. These are diagnostic relative levels, not an
absolute source-normalized TL claim.

The lower panels use all 5001 stored rays. Relative reflected-beam level spans
`0--7.65162 dB` in both cases, while the pressure-release reflection jump is
exactly `0 dB` in the stored precision. The native sinusoidal ATI pressure
record is identically zero in the existing MAT file, so it is explicitly
marked unavailable and is not converted to a logarithmic TL value. The tilted
native ATI points are shown only as a diagnostic; their non-flat values are not
a new hard-gate result. The figure therefore supplements, but does not replace,
the existing geometry, phase, delay, or local Reflect2D validation.

### 4.4 Dense two-dimensional coherent TL map

The sparse receiver records were not sufficient for a standard Bellhop-style
TL background, so the accepted validation-only tilted and sinusoidal binaries
were rerun with only the receiver grid changed. The existing `.sbp`, wall
profiles, `Reflect2D` path, beam count (`5001`), step (`0.05 m`), and uniform
`c=1500 m/s` were retained. The dense mapped grid contains `603` ranges
(`100.05--130 m`, including `102` and `103 m`) and `321` transverse depths
(`-32--32 m`). The raw SHD fields and metadata are under
`results/validation/bellhop_internal_wall_tl_grid/`.

`scripts/reporting/generate_bellhop_internal_wall_tl_2d.py` converts the saved
mapped field with `r=2R0-r'`, `z=-z'`, reverses both array axes to restore
increasing physical coordinates, and plots

\[
  TL(r,z)=-20\log_{10}|p(r,z)|
\]

as the main color background in `06_internal_wall_tl_2d.png`. The plotted field
is the coherent **post-wall reflected branch** produced by the current
branch-isolated internal-wall validation binary; it is not a direct-plus-
reflected total field. Zero-pressure cells with no beam support are masked.
White curves show selected reflected ray trajectories, while the wall is shown
in black. The color range `40--70 dB` is a display limit; raw TL extrema and
grid metadata are recorded in `internal_wall_tl_2d_summary.csv`.

## 5. Interpretation and limitations

The enhancement case is intended only to extend visual coverage from the
existing weak tilted/sinusoidal examples. It demonstrates the distinction
between the native physical backward branch and the post-reflection proper
rotation used for a monotone positive-range chart. It does not assert a new
physical model, PE equivalence, rough-surface communication performance, or
full PM/PE receiver-field validation.

The requested MATLAB load/static check was attempted with the installed MATLAB
2025b executable, but MATLAB failed before startup with
`System Error: File system inconsistency`. This is an environment startup
failure, not a validation-data failure; Python syntax/render checks and the
stored MATLAB-produced summary remain independently available.
