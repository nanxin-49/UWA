# PE--Bellhop--Helmholtz BIE R3 weak three-way closure

日期：2026-09-14  
状态：**R3 PASS**

- Surface: `A=0.01 m`, `K=0.1 rad/m`, fixed C2 physical taper `42 -> 50 m`.
- Convention: `G_BH=conj(G_BH_native) and G_BIE=conj(G_BIE_native), fixed mapping from exp(-iwt) native Helmholtz convention to the PE comparison convention; no fitted scalar`.
- Thresholds (existing weak-limit floor + U_BIE): E `0.0108994522`, phase `0.0108992901 rad`, TL `0.0107982122 dB`.

| pair | E_G | phase RMS rad | TL RMS dB | rho_raw | rho_shape | phi0 rad | E_aligned |
|---|---:|---:|---:|---:|---:|---:|---:|
| PE_BIE | 0.00315757147 | 0.00311763566 | 0.00434837965 | 0.999995015 | 0.99999502 | -0.00010327209 | 0.00315588382 |
| BH_BIE | 8.09197504e-06 | 5.88616577e-06 | 4.82301973e-05 | 1 | 1 | 9.8175707e-08 | 8.09137946e-06 |
| PE_BH | 0.00315854889 | 0.00311863962 | 0.00434761806 | 0.999995012 | 0.999995017 | -0.000103369618 | 0.00315685856 |
| BH_native_BIE_native | 8.09197504e-06 | 5.88616577e-06 | 4.82301973e-05 | 1 | 1 | -9.8175707e-08 | 8.09137946e-06 |
| PE_BIE_native_diagnostic | 0.457769749 | 0.464112028 | 0.00434837965 | 0.895223375 | 0.895223375 | -1.76873651e-06 | 0.457769982 |

## Bellhop geometry

- wall residual max: `3.53850105e-12 m`; phase-jump error max: `7.10542736e-15 rad`.
- p/q rotation errors: `0 / 0`; minimum post-wall dr: `0.0433012702 m`.

## Gates

- `receiver_grid`: PASS
- `bellhop_geometry`: PASS
- `PE_BIE_E`: PASS
- `PE_BIE_phase`: PASS
- `PE_BIE_TL`: PASS
- `PE_BIE_shape`: PASS
- `BH_BIE_E`: PASS
- `BH_BIE_phase`: PASS
- `BH_BIE_TL`: PASS
- `BH_BIE_shape`: PASS
- `finite`: PASS
- `all`: PASS

Artifact: `results\validation\pe_bellhop_helmholtz_bie_reference\R3_weak_three_way\R3_weak_three_way_validation.mat`.
