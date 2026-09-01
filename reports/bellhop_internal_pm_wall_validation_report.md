# Bellhop 2020 fixed-seed 1-D PM internal-wall validation

## Decision

**Status: NOT_FEASIBLE for the requested native-ATI ↔ rotated-wall PM
comparison (validation stopped before any PE comparison).**

The validation-only overlay is implemented and the isolated Bellhop 2020
binary builds and runs.  The hard native/rotated acceptance cannot be claimed
for the specified zero-clearance native case: a native `C`-ATI top boundary is
a forward-going sea-surface reflector, whereas the rotated internal-wall
branch is a backward-going wall reflector that is only made forward-going by
the post-reflection half-turn.  With `z=eta(r)`, the official native case
therefore does not supply the same returned branch at the frozen native
receiver range.  The PM profile also creates near-vertical sampled segments;
the native C-ATI `Dss` curvature and the `Reflect2D` `1/Th` kick become
sampling/grazing sensitive rather than convergent.

No amplitude or phase fitting, p/q reset, receiver normalization, PE
operation, Kirchhoff screen, SSA, or change to the official Bellhop executable
was used.

## Implemented files

- `scripts/validation/validate_bellhop_internal_pm_wall.m`
- `scripts/validation/support/sample_fixed_pm_profile_vertical.m`
- `scripts/validation/support/run_bellhop_internal_pm_wall_poc_vertical.m`
- `scripts/validation/support/run_bellhop_native_pm_wall_vertical.m`
- `scripts/validation/support/bellhop_internal_pm_wall_poc/bellhop.f90`
- `scripts/validation/support/bellhop_internal_pm_wall_poc/Step.f90`
- `scripts/validation/support/bellhop_internal_pm_wall_poc/build_bellhop_internal_pm_wall_poc.ps1`

The profile generator uses one fixed MATLAB `rng(seed,'twister')` realization
and writes those exact samples to both the native `.ati` and rotated `.iwpm`
files.  The rotated overlay extends the sampled wall in constant-depth pieces,
limits both `ReduceStep2D` trial reductions at the first segment intersection,
calls the unmodified native TOP/vacuum `Reflect2D` once, and applies only the
validated proper half-turn `(r',z')=(2R0-r,-z)` before the unchanged Cartesian
influence path.

## Build and executable smoke

The overlay was compiled from the official OALIB 2020 source tree:

`E:/MISC/BELLHOP/AcousticsToolbox_2020`

using the existing MinGW Fortran toolchain.  The resulting independent
executable is
`results/validation/bellhop_internal_pm_wall_poc/bin/bellhop_iwall_pm_2020.exe`;
the source and executable hashes are recorded in
`results/validation/bellhop_internal_pm_wall_poc/build_manifest.json`.

Direct executable smoke runs with supplied sampled profiles (N=129 and
N=257) completed one wall reflection for the full internal fan.  The logged
invariants were machine-level wall residual, one `+pi` pressure-release phase
jump, zero q jump, zero p/q change under the proper rotation, and positive
post-wall range increments.  These are overlay smoke checks, not a formal
MATLAB fixed-seed acceptance result.

## PM geometry diagnostics

The MATLAB validator records, for every sampled profile and step/beam case:

- RMS height/slope and maximum slope;
- RMS and maximum curvature from the same Fourier realization;
- profile range-turn count and minimum absolute rotated `dr`;
- wall intersection residual, local tangent/normal, signed `kappa`;
- incident/reflected directions, `RN/RM`, p/q before and after `Reflect2D`;
- pressure-release phase and amplitude jump;
- travel time, post-rotation minimum range increment, grazing and vertex hits;
- explicit SHD target-range selection and native total-minus-direct field.

The convergence gate treats geometry, phase, delay, p/q transport and range
monotonicity as hard checks.  The known native backward-range Cartesian
amplitude offset is diagnostic only and is never renormalized.

## Why the specified native comparison stops

For the requested construction the same samples are used as

`native: (r,z)=(R0+s, eta(s))`,

`rotated wall: (r,z)=(R0-eta(s), s)`.

The native object is Bellhop's top boundary.  Its normal reflection leaves the
ray on the positive native range axis; it does not create the negative-range
post-wall segment that the proper half-turn is designed to remap.  Reading a
native 97 m column as a returned wall contribution therefore either gives no
native reflection or mixes a different path.  Changing the receiver column
cannot repair that physical chart mismatch; the frozen explicit selector only
prevents the earlier 102/103 m indexing error.

In addition, a zero-clearance source at `z=0` is not robust for a zero-mean PM
ATI: a positive first extension puts the source above the native boundary, and
upward fan rays can intersect the rough top before the nominal wall range.
The rotated case has matched half-spaces and does not share this native
surface-crossing condition.

The PM spectrum is not artificially smoothed or normalized.  At the tested
sample densities the rotated `r=R0-eta` polyline can contain very small `dr`
segments.  Bellhop's own `Dss` override then amplifies the curvature term
through `Reflect2D`'s grazing-sensitive `RN=2*kappa/(c^2*Th)`.  A change of
profile density consequently changes the p kick and maximum sampled kappa by
orders of magnitude in the standalone smoke, which fails the requested
curvature-kick convergence criterion.

## MATLAB host limitation

The formal MATLAB driver could not be launched on this host: R2025b exits at
startup with `System Error: File system inconsistency` (the same local
DDUX/thread-pool runtime fault documented by the preceding covariance audit).
No launcher workaround was retained.  This is separate from the successful
standalone Fortran binary smoke and is not treated as numerical evidence of a
pass.

## Final answers

1. **Fixed-seed PM wall:** overlay implemented, but the requested native ↔
   rotated hard validation is **not feasible** with native `z=eta(r)` and the
   frozen zero-clearance/backward-range chart; no pass is claimed.
2. **Geometry/curvature/p/q:** the internal overlay uses the shared samples,
   Bellhop TOP frame and native `Dss`/`Reflect2D`; internal smoke preserves the
   expected direction, phase and beam-state invariants.  Native equivalence is
   the failed gate, not a modified curvature formula.
3. **Receiver phase/delay:** no valid native returned branch exists at the
   frozen comparison column, so phase/delay equivalence is not accepted.
4. **Backward-range amplitude:** remains diagnostic only; no calibration was
   applied.
5. **PM slope/curvature suitability:** the unfiltered PM realization is too
   close to a sampled vertical/near-grazing polyline for this rotated-wall
   route to provide a density-convergent native comparison.
6. **Next stage:** do **not** proceed to
   `PE_rough <-> Bellhop_rotated,rough`.  A future attempt would first need a
   separately approved native coordinate construction (and an explicit source
   clearance/receiver geometry), followed by a new validation design; this
   task intentionally does not implement that redesign.
