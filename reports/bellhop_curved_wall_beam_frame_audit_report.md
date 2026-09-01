# Bellhop curved-wall beam/frame covariance audit

## Scope and case

This is a Bellhop-only, validation-only audit.  It does not modify PE,
communication code, `InfluenceGeoHatCart`'s propagation formula, or the
production Bellhop executable.  The independent Bellhop 2020 overlay logs
the native ATI case and the rotated internal-wall case using the same weak
sinusoid,

\[
r_w(z)=100-0.25\sin(-0.01z)\;\mathrm{m},
\]

at 4 kHz, with 161 profile samples, 0.1 m nominal ray step and 5001 beams.
The wall is pressure-release and each internal ray is allowed one reflection.
The complete raw logs are in
`results/validation/bellhop_curved_wall_beam_frame_audit/cases/`; compact
tables are saved as [reflection_frame_comparison.csv](../results/validation/bellhop_curved_wall_beam_frame_audit/reflection_frame_comparison.csv)
and [receiver_contribution_fan.csv](../results/validation/bellhop_curved_wall_beam_frame_audit/receiver_contribution_fan.csv).

## Three representative rays

The reflection fan spacing is 0.012 degrees, so the requested left/central/
right rays are the nearest available records: -0.036, 0, and +0.024 degrees.
The maximum native-minus-rotated errors over all 5001 paired first-reflection
records are:

| quantity | maximum absolute difference |
|---|---:|
| wall intersection range | 1.00e-6 m |
| wall intersection depth | 5.00e-7 m |
| wall tangent/normal component | 2.51e-11 |
| signed curvature | 3.13e-7 1/m |
| `Tg`, `Th` | 1.67e-14, 2.29e-15 |
| `RN_raw`, `RN_curv`, `RN_final` | 4.17e-10 |
| reflected tangent/normal component | 5.00e-11 |
| pressure-release phase increment | 0 rad error; exactly +pi |
| Amp jump | 0 |
| travel-time component | 6.67e-10 s |
| frame handedness | 0 |

The small center-ray curvature difference is a profile/ATI floating-point
segment-side effect (the native ATI intersection is about 1 micrometre on the
other side of the central sample); it produces only a 4.17e-10 difference in
the final `RN` term.  It is not a systematic frame rotation or curvature-sign
error.  For the off-axis representative rays the signed curvature and all
`RN` terms agree to the printed precision.

For each representative ray the audit stores the complete values before and
after `Reflect2D`: incident ray frame, wall tangent/normal, `kappa`, `Tg`/
`Th`, `RN`/`RM`, `p/q`, reflected ray frame, Amp, Phase and tau.  The raw
three-row comparison is in the CSV above.

## Q1 — `Reflect2D` exit covariance

**No physical mismatch was found at the `Reflect2D` exit.**  Native and
rotated wall tangent/normal, reflected ray tangent/normal, `Tg`/`Th`, `RN` and
`RM` agree within the tolerances above.  The observed absolute `q` difference
at the center ray is 1.50e-3 on a value of 1.50e5 (about 1e-8 relative), caused
by the same micrometre-scale ATI range quantization; `p` differs by at most
6.26e-5 and the curvature term is otherwise identical.  No p/q reset or
fitting is present.

The log is produced by the existing 2-D `Reflect2D` routine, called with
`'TOP'` and `HS%BC='V'`.  Thus the mirror direction, pressure-release phase
change, and Gaussian-beam curvature update are Bellhop-native.

## Q2 — single-ray `InfluenceGeoHatCart` contribution

`InfluenceGeoHatCart` was not changed.  Its validation-only logger records the
segment tangent/normal, interpolated p/q/tau, footprint, KMAH phase, geometric
factors and the complex contribution before the normal field update.

Only the two rays whose beam footprint contains the target receiver produce a
nonzero accepted contribution at this receiver (the neighbouring right ray is
outside the Cartesian hat footprint).  For the accepted fan:

| metric | value |
|---|---:|
| phase difference range | -1.41e-11 to 3.35e-5 rad |
| amplitude-ratio range | 1.0000000000 to 1.0000039705 |
| summed logged-contribution phase difference | 1.02e-5 rad |
| summed logged-contribution TL difference | 1.05e-5 dB |

There is therefore no approximately 2.095 rad single-ray phase difference.
Across the accepted left/central fan samples the phase is a common near-zero
quantity; no curvature-sign-dependent phase jump is observed.  The unaccepted
right representative is absent from this table because its hat footprint does
not contain the target receiver, not because its reflection state was dropped.

## Q3 — coherent accumulation

The coherent field was checked with a separate direct-only native run so that
the reflected field is formed as `native_total - native_direct` at the native
receiver column.  The receiver columns are intentionally different physical
coordinates: rotated target range 103 m, native target range 97 m.

| field | complex pressure |
|---|---|
| rotated internal-wall at 103 m | -0.00485955644 + 0.00840499904i |
| native total at 97 m | 0.000147136801 - 0.000267009280i |
| native direct at 97 m | 0.00515463902 - 0.00892809685i |
| native reflected (difference) | -0.00500750222 + 0.00866108757i |

Rotated versus native reflected gives phase difference **1.02e-5 rad** and
TL difference **-0.26064 dB**.  The absolute complex difference is
`2.9575e-4` (relative complex error `2.9562e-2`), dominated by the stable
backward-range amplitude factor.  This is the previously observed
`ScalePressure` diagnostic, not an irregular curvature-dependent error.  The
phase agrees; coherent accumulation does not create a new discrepancy.

The former approximately -2.095 rad observation was caused by reading the
rotated SHD first range column (102 m) instead of the requested 103 m target
column.  At the correct column, the phase bias disappears.  No phase or
amplitude calibration was applied.

## Proper rotation and beam state

The post-wall mapping remains the fixed proper rotation

\[
(r',z')=(2R_0-r,-z).
\]

Equivalently, `x' = M x + b` with `M=-I`, `b=(2R0,0)`, and
`det(M)=+1`; both tangent and ray-normal vectors are multiplied by `M`.
This is an orientation-preserving half-turn, not a second physical boundary
reflection.

The rotation logger records transformed ray tangent and ray-normal, p/q, Amp,
Phase, tau, handedness and q sign.  For the complete rotated fan,
handedness is 1 to better than 3.2e-10, `q_sign=+1`, and p/q/tau/Amp/Phase
are carried through unchanged.  The only pressure-release phase increment is
the single `Reflect2D` increment of +pi.  The transformed branch has positive
range increments; the minimum recorded post-wall increment is 0.08639 m.
The representative receiver travel time is within 1.8e-7 s of 103/1500 s;
the center ray error is 2.5e-8 s.

## Classification and decision

Three-way physical classification:

* **A:** no — no reflection-frame/p/q mismatch at `Reflect2D` exit.
* **B:** no — accepted single-ray Cartesian contributions are phase-covariant.
* **C:** no in the strict “accumulation creates the error” sense — the
  corrected coherent reflected field is also phase-covariant.

For the requested label, this is **C-like bookkeeping/indexing-only**: the
only 2.095-rad discrepancy was receiver-column/index pairing, not Bellhop
reflection physics, beam-frame covariance, or coherent accumulation.  The
minimum non-physical fix is to select the explicit target range before SHD
comparison and to subtract native direct from native total at the same native
range.  No change to `Reflect2D`, p/q, curvature, receiver normalization or
`InfluenceGeoHatCart` is indicated.

Consequently the curved-wall beam/frame audit passes its localization goal.
There is no final repair to implement in this task.  A tilted/sinusoidal-to-PM
progression is not started here; a fixed-seed PM wall remains a separate,
future validation stage.

## Reproducibility and source

Validation entrypoint:
[validate_bellhop_curved_wall_beam_frame_audit.m](../scripts/validation/validate_bellhop_curved_wall_beam_frame_audit.m)

Runner and support overlay:
[run_bellhop_curved_wall_beam_frame_audit_vertical.m](../scripts/validation/support/run_bellhop_curved_wall_beam_frame_audit_vertical.m)
and the validation-only Bellhop 2020 overlay under
`scripts/validation/support/bellhop_curved_wall_beam_frame_audit/`.

## Engineering freeze

The passing state is frozen with an explicit SHD selector,
`select_bellhop_shd_pressure_at_range_vertical`, which requires a unique
receiver range and depth within `1e-6 m`; no nearest-column fallback remains
in the flat, tilted, sinusoidal, or covariance validators.

The permanent regression
`validate_bellhop_shd_receiver_range_pairing_vertical` was rerun on the saved
Bellhop fields. It selected rotated range 103 m at column 2 and native
total/direct range 97 m at column 1. The correctly paired phase difference is
`1.0219e-5 rad`; deliberately selecting the rotated 102 m column reproduces
`-2.0946 rad`. The TL difference is `-0.26064 dB` and relative complex error
is `0.029562`, consistent with the already documented backward-range
amplitude diagnostic. Its compact result is saved as
`receiver_range_pairing_regression.csv`.

All internal-wall validators and their support runners pass MATLAB static
analysis without warnings. A fresh full binary audit was not accepted as new
evidence during freeze because this host's MATLAB R2025b batch process
intermittently fails to complete external validation-executable launches and
then terminates in the DDUX shutdown library. The standalone validation
binary remains runnable, and the receiver-pairing regression completes and
prints a passing table before the same shutdown fault. No process-launch
workaround was retained in source, and no existing accepted physics result
was replaced.
