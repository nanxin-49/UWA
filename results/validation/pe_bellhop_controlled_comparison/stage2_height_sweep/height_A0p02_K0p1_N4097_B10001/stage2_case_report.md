# Stage 2 height-sweep case

- A=0.02 m, K=0.1 rad/m, `2kA=0.670206433`.
- source tag: `height_A0p02_K0p1_N4097_B10001`

| max slope | max/RMS curvature (1/m) | E_G (M99) | phase RMS (rad) | TL RMS (dB) | phi0 (rad) | rho_raw | rho_shape | E_aligned | geometry | finite |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---:|:---:|
| 0.00199999977 | 0.0002 / 0.000140143012 | 0.00823396805 | 0.00816578444 | 0.00914010217 | -0.00270573185 | 0.999966108 | 0.999969768 | 0.0077780338 | PASS | PASS |

Geometry guards: wall residual `7.1e-15 m`; pressure-release phase-jump error `7.11e-15 rad`; p/q rotation errors `0/0`; min post-wall dr `0.0433501 m`; center tau error `4.7e-15 s`.

Comparison convention: `G_BH_comparison = conj(G_BH_abs) = G_BH_raw`.
