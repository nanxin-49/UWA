# Stage 1 weak sinusoid

状态：**FAIL**

- A=0.01 m, K=0.1 rad/m; source `X`, run `C`.
- Floors: F_E=0.0058993922, F_phi=0.0058992479, F_TL=0.00079785258; thresholds T_E=0.010899392, T_phi=0.010899248, T_TL=0.010797853.

| comparison | L2(M99) | phase RMS | TL RMS (dB) | rho |
|---|---:|---:|---:|---:|
| profile | 1.2024712e-06 | 1.0348829e-06 | 5.3182868e-06 | 1 |
| beam | 1.543684e-06 | 1.5104209e-06 | 2.7691186e-06 | 1 |
| model | 0.47247164 | 0.47919524 | 0.0045698832 | 0.91922588 |

## Geometry

- wall residual max: 7.102e-15 m; phase jump error max: 7.105e-15 rad; min post dr: 0.0433257 m; center tau error: 4.677e-15 s.

## Checks

- profile_l2: PASS
- profile_phase: PASS
- profile_tl: PASS
- beam_l2: PASS
- beam_phase: PASS
- beam_tl: PASS
- geometry: PASS
- weak_E: FAIL
- weak_phase: FAIL
- weak_tl: PASS
- weak_shape: FAIL
- finite: PASS
- all: FAIL

Stage 2 remains locked unless this weak-limit result passes.
