# Stage 0 flat reflected-line closure

状态：**PASS**

- f/c: 4000 Hz / 1500 m/s; source `X`, run `C`; beams [5001,10001].
- AS footprints: M95 29.4922 m (303), M99 40.0391 m (411).

| metric | L2(M99) | phase RMS(M99) |
|---|---:|---:|
| pe_as | 4.2785525e-13 | 2.3912218e-13 |
| flat_internal | 0.0058993922 | 0.0058992479 |
| chain_internal_free | 0 | 0 |
| dedicated_vs_internal | 0 | 0 |
| beam_convergence | 5.854857e-06 | 5.8535549e-06 |

## Checks

- pe_as_complex: PASS
- pe_as_l2: PASS
- pe_outer5: PASS
- receiver_coordinates: PASS
- flat_l2: PASS
- flat_phase: PASS
- flat_phase_p95: PASS
- flat_tl_p95: PASS
- flat_shape: PASS
- flat_global_phase: PASS
- flat_aligned: PASS
- chain_l2: PASS
- chain_phase: PASS
- chain_tl: PASS
- beam_l2: PASS
- beam_phase: PASS
- beam_tl: PASS
- half_dx: PASS
- geometry: PASS
- finite: PASS
- all: PASS

No subsequent stage is unlocked unless this report is reviewed and passes.
