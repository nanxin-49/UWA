# Stage 2 height-sweep case

- A=0.05 m, K=0.1 rad/m, `2kA=1.67551608`.
- source tag: `height_A0p05_K0p1_N4097_B10001`

| max slope | max/RMS curvature (1/m) | E_G (M99) | phase RMS (rad) | TL RMS (dB) | phi0 (rad) | rho_raw | rho_shape | E_aligned | geometry | finite |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---:|:---:|
| 0.00499999941 | 0.0005 / 0.000350354768 | 0.0210357376 | 0.0208690864 | 0.0228496268 | -0.00818966942 | 0.999778741 | 0.99981227 | 0.0193849209 | PASS | PASS |

Geometry guards: wall residual `3.54e-12 m`; pressure-release phase-jump error `7.11e-15 rad`; p/q rotation errors `0/0`; min post-wall dr `0.0434233 m`; center tau error `4.16e-17 s`.

Comparison convention: `G_BH_comparison = conj(G_BH_abs) = G_BH_raw`.
