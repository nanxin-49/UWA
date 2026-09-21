# Stage 2 height-sweep case

- A=0.01 m, K=0.1 rad/m, `2kA=0.335103216`.
- source tag: `reused_stage1_convention_fixed`

| max slope | max/RMS curvature (1/m) | E_G (M99) | phase RMS (rad) | TL RMS (dB) | phi0 (rad) | rho_raw | rho_shape | E_aligned | geometry | finite |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---:|:---:|
| 0.000999999883 | 0.0001 / 7.00715847e-05 | 0.00409345672 | 0.00405930291 | 0.00456997571 | -0.00125773333 | 0.999991624 | 0.999992415 | 0.00389576987 | PASS | PASS |

Geometry guards: wall residual `7.1e-15 m`; pressure-release phase-jump error `7.11e-15 rad`; p/q rotation errors `0/0`; min post-wall dr `0.0433257 m`; center tau error `4.68e-15 s`.

Comparison convention: `G_BH_comparison = conj(G_BH_abs) = G_BH_raw`.
