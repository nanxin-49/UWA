# Stage 0 flat reflected-line closure

状态：**FAIL**

- f/c: 4000 Hz / 1500 m/s; source `X`, run `C`; beams [101,201].
- AS footprints: M95 6.00586 m (5), M99 12.0117 m (9).

| metric | L2(M99) | phase RMS(M99) |
|---|---:|---:|
| pe_as | 0.04641982 | 0.020546284 |
| flat_internal | 1.2354817 | 0.48227774 |
| chain_internal_free | 0 | 0 |
| dedicated_vs_internal | 0 | 0 |
| beam_convergence | 0.01346778 | 0.010560505 |

## Checks

- pe_as_complex: FAIL
- pe_as_l2: FAIL
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
- beam_l2: FAIL
- beam_phase: FAIL
- beam_tl: PASS
- half_dx: PASS
- geometry: FAIL
- finite: PASS
- all: FAIL

No subsequent stage is unlocked unless this report is reviewed and passes.
