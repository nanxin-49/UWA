# Stage 2 height-sweep case

- A=0.2 m, K=0.1 rad/m, `2kA=6.70206433`.
- source tag: `height_A0p2_K0p1_N4097_B10001`

| max slope | max/RMS curvature (1/m) | E_G (M99) | phase RMS (rad) | TL RMS (dB) | phi0 (rad) | rho_raw | rho_shape | E_aligned | geometry | finite |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---:|:---:|
| 0.0199999977 | 0.002 / 0.00140122191 | 0.100485227 | 0.100141661 | 0.0914359332 | -0.0609643815 | 0.994943993 | 0.996795793 | 0.0801428322 | PASS | PASS |

Geometry guards: wall residual `3.52e-12 m`; pressure-release phase-jump error `7.11e-15 rad`; p/q rotation errors `0/0`; min post-wall dr `0.0437893 m`; center tau error `1.39e-17 s`.

Comparison convention: `G_BH_comparison = conj(G_BH_abs) = G_BH_raw`.
