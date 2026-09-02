# Bellhop 2020 fixed-seed PM internal-wall validation overlay

This directory contains the modified `bellhop.f90` and `Step.f90` files used
by the validation-only fixed-seed 1D PM rough-wall POC. Their bases are copied from
`AcousticsToolbox_2020/Bellhop` (`2020_11_4`).

The build script applies the two overlays under
`results/validation/bellhop_internal_pm_wall_poc/build` and produces a
separate executable. It never replaces the official Bellhop executable.

The validation executable reads an `.iwpm` sidecar containing `R0`, mapped
receiver range, the fixed profile seed, and the sampled `(range, depth)` wall
profile. The profile samples are generated once in MATLAB and the same points
are written to the rotated internal wall and the native C-ATI case. The
internal wall uses only the finite sampled support; it does not add ATI-style
constant-depth/infinite-range endpoint extensions.

The completed fixed-realization convergence entrypoint (archived under
`cash/bellhop_internal_wall_superseded_20260902/`) stored each case under
`results/validation/bellhop_internal_pm_fixed_realization/` and obtained all
profile densities by interpolation from one master Fourier realization. The
validation-only `.iwdiag` log additionally exports the accepted segment and
local interpolation parameter, plus the native `Tg`, `Th`, `RM` and applied
`RN` values for density-to-density beam-state comparison.

The same isolated binary retains the local covariance audit sidecar
`.iwcov` for reproducibility. Mode `0` instruments a native C-ATI first top
reflection and stops after that `Reflect2D` exit; mode `1` keeps the parametric
internal wall and its existing proper half-turn. The one-off covariance runner
and detailed stage report are archived under
`cash/bellhop_internal_wall_superseded_20260902/`; the accepted conclusion is
merged into the authoritative Bellhop implementation report. No receiver-field
equality is implied by this local audit.

At an accepted wall segment, the unit tangent is obtained from the ordered
parametric polyline, the outward TOP normal is reconstructed from that tangent,
and signed geometric curvature is computed from wrapped tangent turning over
local arc length. The physical reflection calls native 2-D
`Reflect2D` once with pressure-release conditions. It then applies the already
validated proper pi rotation and passes only the transformed, monotonically
increasing branch to the unchanged Cartesian geometric-hat influence routine.

The scope is uniform, lossless `c=1500 m/s`, one source at transverse depth
zero, one wall reflection, pressure release, and the existing Gaussian `.sbp`.
No Kirchhoff, SSA, phase-screen or PE operation is present. This is a
fixed-realization geometry validation, not a PM scattering model.
