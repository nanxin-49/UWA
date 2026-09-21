# Stage 2 height-sweep case

- A=0.1 m, K=0.1 rad/m, `2kA=3.35103216`.
- source tag: `height_A0p1_K0p1_N4097_B10001`

| max slope | max/RMS curvature (1/m) | E_G (M99) | phase RMS (rad) | TL RMS (dB) | phi0 (rad) | rho_raw | rho_shape | E_aligned | geometry | finite |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---:|:---:|
| 0.00999999883 | 0.001 / 0.000700689815 | 0.0441790479 | 0.043877023 | 0.0457022488 | -0.0211065624 | 0.999023639 | 0.999246207 | 0.0388527853 | PASS | PASS |

Geometry guards: wall residual `3.54e-12 m`; pressure-release phase-jump error `7.11e-15 rad`; p/q rotation errors `0/0`; min post-wall dr `0.0435453 m`; center tau error `2.78e-17 s`.

Comparison convention: `G_BH_comparison = conj(G_BH_abs) = G_BH_raw`.
