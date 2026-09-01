# Bellhop 2020 smooth sinusoidal internal-wall POC

## Result

Status: **Historical conclusion superseded by the beam/frame covariance audit**.

The follow-up audit in
`reports/bellhop_curved_wall_beam_frame_audit_report.md` corrected the earlier
coherent-field extraction: the rotated target is the 103 m SHD column, not
the first 102 m column. With the explicit target column and a native
total-minus-direct comparison, the reflected-field phase difference is about
`1.02e-5 rad`; the remaining `-0.2606 dB` is the known backward-range
`ScalePressure` diagnostic. The former `-2.095 rad` value was a receiver
column/indexing artifact, not a Bellhop beam-frame discrepancy.

The validation-only curved-wall mechanism reaches the sampled wall and
preserves the native reflection state. The end-to-end coherent comparison
should be interpreted using the corrected target-column bookkeeping described
above. No PM random wall was implemented.

## Implemented scope

The overlay is built from Bellhop 2020 sources and is independent of the formal
`bellhop.exe`:

- `scripts/validation/support/bellhop_internal_sinusoidal_wall_poc/Step.f90`
- `scripts/validation/support/bellhop_internal_sinusoidal_wall_poc/bellhop.f90`
- `scripts/validation/support/bellhop_internal_sinusoidal_wall_poc/build_bellhop_internal_sinusoidal_wall_poc.ps1`
- `scripts/validation/support/run_bellhop_internal_sinusoidal_wall_poc_vertical.m`
- `scripts/validation/support/run_bellhop_native_sinusoidal_wall_vertical.m`
- `scripts/validation/validate_bellhop_internal_sinusoidal_wall_poc.m`

The internal wall is represented by the same sampled points used for native
`C`-ATI.  The overlay adds the exact piecewise-linear segment intersection to
both `ReduceStep2D` calls.  At the accepted intersection it interpolates node
tangent/normal exactly as `ComputeBdryTangentNormal`/native C-ATI does, uses the
same `Dss` curvature override, and calls the unmodified Bellhop `Reflect2D`
body.  The reflected node is then transformed only by the proper rotation;
the sidecar profile is checked against the requested analytic sinusoid before
the run starts.

\[
  (r',z')=(2R_0-r,-z),\qquad R_0=100\ {m m},
\]

without a second reflection or any p/q, amplitude, phase, or travel-time
reset.  Only that transformed branch is passed to the existing influence
routine.

In plain text, the fixed remapping is `(r',z')=(2R0-r,-z)` with `R0=100 m`.

For a positive-range ATI representation around the unchanged `z=0` datum, the
validation uses the equivalent signed orientation `K=-0.01 m^-1` in
`r=R0-A sin(Kz)`.  This gives `dr/dz>0` over `z∈[-100,100]` and is just the
vertical reflection of the same smooth sinusoid; it does not add a phase-screen
or scattering model.  Weak and stronger amplitudes were `A=0.25 m` and
`A=0.5 m`.

## Numerical checks

The small scan used `A={0.25,0.5} m`, profile counts `{41,81,161}`, step
lengths `{0.2,0.1,0.05} m`, and beam counts `{2001,5001}` at 4 kHz.  The raw
table is [sinusoidal_wall_convergence.csv](../results/validation/bellhop_internal_sinusoidal_wall_poc/sinusoidal_wall_convergence.csv).

| quantity | observed range / maximum |
|---|---:|
| wall intersection residual | `7.1e-15` m typical; `3.54e-12` m maximum |
| analytic tangent/normal error | `3.64e-6` maximum, decreasing with profile density |
| specular direction error | `6.72e-16` maximum |
| proper-rotation direction error | `0` |
| pressure-release phase increment error | `7.11e-15` rad |
| reflection Amp jump | `0` |
| p reflection jump | up to `7.83e-3` (native curvature kick) |
| q reflection jump | `0` |
| p/q change under rotation | `0` |
| sampled signed curvature | `1.40e-5`–`1.46e-5 m^-1` for A=0.25; `2.79e-5`–`2.92e-5 m^-1` for A=0.5 |
| path-to-range-plane error | below `8.8e-12` m |
| transformed minimum range increment | `4.31e-2` m or larger |
| imaginary travel time | `0` |

Thus the wall intersection, local frame, nonzero curvature, `Reflect2D` phase
and p/q kick, and proper rotation all behave as intended.

The native `.ati` and internal `.iw3` files are emitted from the same profile
arrays, and the overlay deliberately uses the same C-ATI node interpolation and
`Dss` expression.  Consequently the native-to-internal sampled tangent,
normal, and signed-curvature comparison is an exact construction check; the
reported `3.64e-6` frame bound is the independent analytic-sinusoid check.

## Native versus rotated coherent field (historical extraction, superseded)

The native case uses the same profile samples in `.ati`, the same `.sbp`,
frequency, source, and numerical settings.  Native reflected pressure is
formed as native total minus the matched free-field run. The earlier
temporary aggregation read the first rotated SHD range (102 m) while
comparing it with the 97 m native reflected field, which produced the
reported `-2.095 rad` value. That comparison is invalid and is retained only
as historical provenance. The corrected target-column result is documented in
the beam/frame audit report; no amplitude or phase renormalization was
applied.

The independent native arrivals run remains geometrically sensible (one top
bounce, zero bottom bounces, delay about `0.0686666 s`, path about `103 m`, and
pressure-release phase `180°`). The follow-up audit finds no unresolved
reflection-frame or Cartesian influence phase discrepancy after explicit
receiver-column pairing.

## Decision

- Native `C`-ATI and internal sampled geometry are consistent to the reported
  tangent/normal/curvature tolerances.
- Bellhop's native `Reflect2D` is actually reused; no duplicate reflection
  formula or extra scattering phase is present.
- The proper π rotation leaves the stored p/q, amplitude, phase, and travel
  time unchanged and keeps the post-wall branch strictly increasing in range.
- The beam/frame audit finds no native↔rotated reflection, single-ray
  contribution, or corrected coherent-field phase mismatch.
- A fixed 1-D PM rough wall is still not implemented in this task; proceeding
  to it remains a separate validation decision.
