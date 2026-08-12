# PE/Bellhop current flat-surface validation

- run_id: `bellhop_current_20260723_rc5`
- core decision: **FAIL_CORE**
- amplitude status: **OPEN**
- Bellhop SHA-256: `e6f9c1bcfd2b0945bfb59607909af589fa7b3f823df8bd5121425736eaebd796`

## Layered acceptance

| Check | Value | Limit | Pass |
|---|---:|---:|:---:|
| phase_audit | 0 | 0 | 1 |
| phase_conversion | 0 | 1e-12 | 1 |
| component_closure | 3.5762241e-18 | 1e-12 | 1 |
| pe_analytic_delay | 0.027281981 | 0.25 | 1 |
| pe_bellhop_delay | 0.027280961 | 0.25 | 1 |
| bellhop_analytic_delay | 3.9304163e-06 | 0.25 | 1 |
| small_offset_pe_delay | 0.027281981 | 0.25 | 1 |
| small_offset_bellhop_delay | 3.4187815e-06 | 0.25 | 1 |
| path_clusters | 0 | 0 | 1 |
| bellhop_coherent_beams | 0.19041019 | 1 | 1 |
| bellhop_incoherent_beams | 0.1019304 | 0.5 | 1 |
| pe_numerical_convergence | 1 | 0 | 0 |
| public_regression | 3.5762241e-18 | 1e-12 | 1 |
| public_defaults | 0 | 0 | 1 |

## Amplitude diagnostic

Scale spread 1.65401 dB; reflection RMS 1.62948 dB; max 2.0751 dB. Source-aware weighting is diagnostic and never replaces raw amplitudes.

## Interpretation

Physical PE fields are obtained through the public phase-reference layer; no validator adds a second carrier. The 4 kHz Bellhop arrivals are extended across 3--5 kHz only for delay/PDP diagnostics, not as a multifrequency Bellhop amplitude model.

The public PE core and default `kirchhoff_spatial` surface remain unchanged. See the same-run ten-figure atlas under `E:\MISC\CARPE3D_matlab\Explain\results\visualization\pe_bellhop_flat_surface_current\bellhop_current_20260723_rc5`.

## Fixed conclusions

1. Carrier sign and reference semantics close: the validator consumes public
   physical/direct-DSP fields and never adds a second carrier.
2. PE, analytic image-source, and Bellhop direct/surface paths agree in delay;
   the largest PE--Bellhop residual is `0.027281 ms`.
3. Sampling, longitudinal step, and sponge sensitivity pass, but aperture
   convergence is not established. With sponge disabled, the 32/64 m edge
   levels are `-1.087/-14.532 dB`, both above the required `-40 dB`; their
   phase/TL difference is `3.526 rad/-3.896 dB`.
4. Absolute amplitude is not yet source-matched. One global direct scale gives
   `1.654 dB` cross-geometry spread and reflected RMS/max residuals of
   `1.629/2.075 dB`; source-aware weighting remains diagnostic only.
5. The evidence supports keeping the public PE marching core and carrier
   correction unchanged. The next Bellhop work should enlarge/validate the
   no-sponge aperture and establish a common source normalization, not alter
   the phase-reference layer.
