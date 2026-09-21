# Stage 1 weak sinusoid

状态：**FAIL**

- A=0.0025 m, K=0.1 rad/m; source `X`, run `C`.
- Floors: F_E=0.0058993922, F_phi=0.0058992479, F_TL=0.00079785258; thresholds T_E=0.010899392, T_phi=0.010899248, T_TL=0.010797853.

| comparison | L2(M99) | phase RMS | TL RMS (dB) | rho |
|---|---:|---:|---:|---:|
| profile | 3.0923133e-07 | 2.5941921e-07 | 1.4618467e-06 | 1 |
| beam | 0 | 0 | 0 | 1 |
| model | 0.11969441 | 0.11980048 | 0.0011424793 | 0.99478804 |

## Geometry

- wall residual max: 7.105e-15 m; phase jump error max: 7.105e-15 rad; min post dr: 0.0433074 m; center tau error: 4.677e-15 s.

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
