# Stage 0 flat reflected-line closure

状态：**FAIL**

- f/c: 4000 Hz / 1500 m/s; source `X`, run `C`; beams [1001,2001].
- AS footprints: M95 12.0117 m (17), M99 16.5161 m (23).

| metric | L2(M99) | phase RMS(M99) |
|---|---:|---:|
| pe_as | 1.7338745e-13 | 1.2344818e-13 |
| flat_internal | 0.43693585 | 0.22514281 |
| chain_internal_free | 0 | 0 |
| dedicated_vs_internal | 0 | 0 |
| beam_convergence | 0.0001305839 | 0.00014087325 |

## Checks

- pe_as_complex: PASS
- pe_as_l2: PASS
- pe_outer5: FAIL
- receiver_coordinates: PASS
- flat_l2: FAIL
- flat_phase: FAIL
- flat_phase_p95: FAIL
- flat_tl_p95: FAIL
- flat_shape: FAIL
- flat_global_phase: FAIL
- flat_aligned: FAIL
- chain_l2: PASS
- chain_phase: PASS
- chain_tl: PASS
- beam_l2: PASS
- beam_phase: PASS
- beam_tl: PASS
- half_dx: PASS
- geometry: PASS
- finite: PASS
- all: FAIL

No subsequent stage is unlocked unless this report is reviewed and passes.
