# Stage 1 weak sinusoid

状态：**FAIL**

- A=0.005 m, K=0.1 rad/m; source `X`, run `C`.
- Floors: F_E=0.0058993922, F_phi=0.0058992479, F_TL=0.00079785258; thresholds T_E=0.010899392, T_phi=0.010899248, T_TL=0.010797853.

| comparison | L2(M99) | phase RMS | TL RMS (dB) | rho |
|---|---:|---:|---:|---:|
| profile | 6.0271858e-07 | 5.1665583e-07 | 2.6958914e-06 | 1 |
| beam | 1.4362984e-06 | 1.4263785e-06 | 1.4641187e-06 | 1 |
| model | 0.23875403 | 0.23959927 | 0.0022849061 | 0.9792833 |

## Geometry

- wall residual max: 7.108e-15 m; phase jump error max: 7.105e-15 rad; min post dr: 0.0433135 m; center tau error: 4.691e-15 s.

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
