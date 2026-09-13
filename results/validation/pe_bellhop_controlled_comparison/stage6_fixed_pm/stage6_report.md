# Stage 6 canonical fixed-PM comparison

状态：**PASS**；validity region: **III**。

- seed 260001, 4 kHz, X/C, 10,001 beams, N=4097.
- coefficient SHA-256 `1f7eda465e4ae85b8ac038310edf053f2d062563a3b943bfd83bd59d07015f67`.
- same Stage-0 flat denominator, receiver line, M95/M99 masks, and official comparison convention.

| E_G | TL RMS (dB) | phase RMS (rad) | phi0 | rho_raw | rho_shape | E_aligned |
|---:|---:|---:|---:|---:|---:|---:|
| 0.958268507 | 1.27399043 | 1.10771425 | -0.798545974 | 0.538113006 | 0.771212836 | 0.678503463 |

Profile: RMS/max height 0.186410373/0.499455828 m; RMS/max slope 0.0539769697/0.116674298; RMS/max curvature 0.0209071545/0.0445943842 1/m; minimum radius 22.4243482 m; 2k sigma_eta 6.24667156.

Bellhop: wall residual 7.27e-15 m; min incidence 0.831901249; grazing fraction 0; hit max |kappa| 0.0425258129 1/m; min post dr 0.0396983251 m; phase error 7.11e-15 rad; p/q rotation 0/0; central/tau range 0.0683459103843 / [0.0680368108941,0.0797840904231] s.

Axis sanity: Delta TL 0.310370739 dB and Delta phase -2.24820073 rad; historical X values 0.310371 dB / -2.2482 rad are diagnostic only.

Cross-model discrepancy is not a Stage-6 hard gate. No amplitude/phase fit, PM smoothing, source refit, or core-physics change is used.
