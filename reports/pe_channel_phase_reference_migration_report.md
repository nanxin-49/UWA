# PE Receiver Carrier-Phase Reference Migration Report

Date: 2026-07-22

## Decision

The issue existed. Direct and surface-reflected PE paths were stored as
separately reduced envelopes and were added without restoring a common
receiver carrier reference. The correction is implemented as a deterministic
post-propagation phase layer. The PE split-step operator, surface models,
joint-kstat residual, exact adjoint kernel, and PM spatial covariance are not
changed.

The public default is now `direct_dsp`. `legacy_reduced` remains available for
regression and model migration.

## Formula

For surface reference (z_s=0),

\[
\tau_{\rm dir,0}=\frac{z_{\rm tx}-z_{\rm rx}}{c_0},\qquad
\tau_{\rm ref,0}=\frac{z_{\rm tx}+z_{\rm rx}-2z_s}{c_0},
\]

\[
\Delta\tau_0=\tau_{\rm ref,0}-\tau_{\rm dir,0}.
\]

The MATLAB-IFFT receiver representation is

\[
H_{\rm dir}^{\rm dsp}=H_{\rm dir}^{\rm red},\qquad
H_{\rm ref}^{\rm dsp}=e^{-i2\pi f\Delta\tau_0}H_{\rm ref}^{\rm red}.
\]

Under (p(t)=\Re\{P(f)e^{-i2\pi ft}\}), the absolute physical components are

\[
H_{\rm dir}^{\rm phys}=e^{+i2\pi f\tau_{\rm dir,0}}H_{\rm dir}^{\rm red},
\qquad
H_{\rm ref}^{\rm phys}=e^{+i2\pi f\tau_{\rm ref,0}}H_{\rm ref}^{\rm red}.
\]

For the standard 100 m transmitter, 3 m receiver, and 1500 m/s sound speed,
the spans are 97 and 103 m and the nominal relative reference delay is 4 ms.
No per-realization peak alignment is performed.

## Interfaces and compatibility

- `apply_pe_channel_phase_reference_vertical` is the central conversion API.
- `vertical_channel_model` accepts `channel_phase_reference='direct_dsp'`
  (default) or `legacy_reduced`.
- Existing public `H_direct_f/H_reflect_f/H_f` fields remain present and now
  use the configured reference. Explicit reduced and physical fields were
  added together with `phase_reference_meta`.
- Cached and adjoint runners expose all three representations. Their stored
  propagation operators and projection weights remain reduced.
- Receiver `C/P` contraction exposes direct-DSP defaults and `_reduced`
  baselines. With (D=\operatorname{diag}(e^{-i2\pi f\Delta\tau_0})),
  (C=D C_{\rm red}D^H) and (P=D P_{\rm red}D^T).
- Conditional model/library schemas are `2.0.0` and `2.0.0-discrete`.
  `upgrade_conditional_channel_phase_vertical` migrates known schema-1
  project artifacts, including complex and augmented-real eigensystems.
- Unknown external `H(f)` without project metadata is assumed DSP-ready and
  is not silently rotated.
- `build_channel_cir_vertical` separates direct-DSP IFFT from
  absolute-physical FFT/N reconstruction. The old CIR helper is retained as a
  common-time-shift wrapper.

The joint surface variable remains

\[
\delta G_i=R_{0,i}e^{i\alpha_i\eta}-R_{{\rm coh},i}.
\]

Its (C/P) already contains (R_0). This migration does not add another
reflection coefficient or height phase to `deltaG`, `q`, or `a`.

## Validation settings and results

Primary entrypoint:

```matlab
run('scripts/validation/validate_pe_channel_phase_reference_vertical.m')
```

It covers:

- analytic F=65, 4--8 kHz, 62.5 Hz frequency spacing;
- all four public surface branches on a reduced PE grid;
- scalar direct-only and `legacy_reduced` reproduction;
- PM 16^2 / PE 8^2 cached, exact-adjoint, dense, and FFT paths;
- synthetic schema-1 conditional-model migration;
- external-H communication metadata policy.

Measured summary:

| Check | Result |
|---|---:|
| Nominal/measured relative delay | 4 / 4 ms |
| Cached/adjoint direct-DSP error | `3.75814e-16` |
| Dense/FFT covariance error | `3.64406e-16` |
| Dense/FFT pseudo-covariance error | `4.16399e-16` |
| Schema migration covariance error | `3.61301e-16` |
| Schema migration pseudo-covariance error | `4.26365e-16` |
| Reduced/DSP/physical round-trip conversion | pass |
| Communication/CIR phase-reference guards | pass |
| Overall phase audit | pass |

Additional regressions:

- exact-adjoint smoke: adjoint `7.77e-15`, projection `5.10e-15`, runner
  reduced regression zero, dense/FFT approximately `1e-15`, F=9 split-floor
  gate passed;
- public/cached F=32 double consistency: total error `6.9292e-16`;
- minimal synthetic and PE-derived wideband communication validation:
  completed successfully.
- regenerated two-node communication smoke: U=5/U=8 cached PE, joint PE,
  full-rank, and 99.9% conditional paths all completed using explicit
  direct-DSP metadata; legacy total-only caches were not guessed or reused.

F=9 (`4:0.5:8 kHz`) is not a carrier-delay acceptance grid. Its 2 ms
unambiguous window aliases 4 ms to zero phase at every frequency sample.

## Changed areas

The implementation touches the public output packaging, central phase helper,
cached/adjoint runners, receiver-statistics contraction, conditional-model
schema/migration, CIR/communication guards, phase/Bellhop validators, and U=5
/ U=8 conditional validation loaders. The public default surface model and
PE marching formulas are unchanged.

## Remaining work

- The existing expensive F=64/512 ensemble was not regenerated in this task.
  Its covariance conclusions are invariant under unitary deterministic phase
  rotation, but it should be rerun before removing compatibility support.
- Bellhop source normalization and doubled-window/sponge sensitivity remain
  independent open issues.
- Old project MAT files without reliable geometry metadata cannot be migrated
  safely from a signed delay scalar alone and are rejected.
- `legacy_reduced` should remain for at least one compatibility cycle.
