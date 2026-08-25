# Random-surface reflected-chain window robustness

Status after 4 kHz stage: **A**. PE marching and production defaults were not changed.

## Environment

- Uniform c=1500 m/s; Tx=(0,0,100 m); Rx=(0,0,3 m); production Gaussian sigma=0.3 m.
- Kirchhoff spatial, U=5 m/s, Hs=[0.05,0.5,1.0] m, seeds 12345--12349.
- Fixed dx=50/256 m; 160/192 m are exact central crops of each 256 m master and are never renormalized; sponge is off.

## 4 kHz pass rates

| sea | Hs (m) | nominal W (m) | pass | rate | worst TL (dB) | worst phase (rad) | worst receiver L2 | worst edge5 | worst boundary (dB) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| smoke | 0.200 | 40 | 1/1 | 1.000 | 0.000601077 | 0.000173209 | 0.000416625 | 2.18166e-05 | -49.5102 |
| smoke | 0.200 | 48 | 1/1 | 1.000 | 0.000514879 | 0.000274041 | 0.000366324 | 2.78839e-06 | -58.4266 |
| smoke | 0.200 | 64 | 1/1 | 1.000 | 0 | 0 | 0 | 5.49347e-07 | -65.4201 |

Preliminary recommendation: **A** -- 192 m + no sponge passed every 4 kHz realization.

The final Q1--Q6 decision is appended after the selected 3--5 kHz qualification.
