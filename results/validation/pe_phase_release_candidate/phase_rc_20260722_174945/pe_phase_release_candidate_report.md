# PE carrier-phase release candidate

- Decision: **PASS**
- run_id: `phase_rc_20260722_174945`
- revision: `435d31fea80285ae2c255c624b5cb65294159c2d`
- fingerprint: `14ccc63217ae08a11948e9fb489f2abe8148331f079f6ec7a136a2e607712df1`
- created: `2026-07-22T17:49:46+08:00`


## Acceptance checks

| Check | Pass | Value | Limit | Note |
|---|:---:|---:|---:|---|
| `phase_reference_validation` | 1 | 1 | 1 | all phase, migration, and compact research gates |
| `nominal_delay_4ms` | 1 | 3.46945e-15 | 0.25 | F=65 error <= one delay sample |
| `migration_full_closure` | 1 | 4.26365e-16 | 1e-12 | mean, C/P, complex and augmented EVD |
| `uniform_phase_audit` | 1 | 1 | 1 | independent angular-spectrum sign audit |
| `adjoint_full` | 1 | 1 | 1 | exact adjoint, projection, dense/FFT, F9 and F64 |
| `f9_split_floor` | 1 | 0.593078 | 1.25 | analytic/sample versus split floor |
| `f64_split_floor` | 1 | 0.850502 | 1.25 | analytic/sample versus split floor |
| `u5_full` | 1 | 1 | 1 | full conditional validation |
| `u8_full` | 1 | 1 | 1 | full conditional validation |
| `conditional_model_schema` | 1 | 1 | 1 | schema 2.x, direct_dsp, zero common shift |
| `conditional_library` | 1 | 1 | 1 | exact U5/U8 nodes |
| `public_regression` | 1 | 1 | 1 | scalar direct-only/direct-plus-reflect and unchanged default |
| `cached_public_f64` | 1 | 9.88957e-16 | 1e-10 | F64 double consistency |
| `two_node_communication` | 1 | 8.52312e-16 | 1e-10 | four sources per node, direct_dsp, finite, paired latent comparison |
| `propagation_atlas` | 1 | 1 | 1 | 13 PNG, MP4, MAT, manifest, summary; same run_id |
| `matlab_static_parse` | 1 | 0 | 0 | no checkcode parse errors in RC sources |

## Architecture decision

Recommend integration into the main branch.

- Public PE: general propagation entry.
- Cached forward: fixed-path high-confidence regression oracle.
- Adjoint projection: exact receiver realization generation in a fixed environment.
- Analytic FFT C/P: receiver second-order statistics without receiver Monte Carlo.
- Conditional generator: large communication Monte Carlo after physical/statistical validation.
