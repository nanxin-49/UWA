# PE--Bellhop--Helmholtz BIE R1 flat validation

日期：2026-09-14  
状态：**R1 PASS**

- source: `discrete-periodic-Gaussian|N=984|W=192.1875|sigma=0.3|f=4000|c=1500|exp(-iwt)`
- formulation: `Dirichlet half-plane Green D_h-i*k*S_h; smooth finite section`
- convention/normal/jump: `exp(-i*omega*t)`; `into water domain z>eta(x)`; `gamma_D D_h = +1/2 I + K_h`

- M95/M99 radii: `29.4921875 / 40.0390625 m`

| sweep | parameter | nodes | boundary residual | linear residual | image L2 M99 | phase RMS M99 | TL RMS M95 (dB) | rho | successive L2 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| points_per_wavelength_actual | 8 | 2560 | 2.39e-08 | 2.94e-15 | 1.31e-07 | 1.02e-07 | 4.37e-07 | 1 | NaN |
| points_per_wavelength_actual | 10 | 3200 | 4.02e-09 | 3.34e-15 | 1.26e-07 | 9.61e-08 | 4.37e-07 | 1 | 1.15e-08 |
| points_per_wavelength_actual | 12 | 3840 | 9.28e-10 | 3.49e-15 | 1.24e-07 | 9.27e-08 | 4.37e-07 | 1 | 7.69e-09 |
| half_width_m | 60 | 2560 | 2.39e-08 | 2.94e-15 | 1.31e-07 | 1.02e-07 | 4.37e-07 | 1 | NaN |
| half_width_m | 70 | 2992 | 4.13e-08 | 3.18e-15 | 8.43e-08 | 7.23e-08 | 3.46e-07 | 1 | 6.83e-08 |
| half_width_m | 80 | 3416 | 6.04e-08 | 3.3e-15 | 6.89e-08 | 6.35e-08 | 2.25e-07 | 1 | 3.13e-08 |

## Hard gates

- `source_dft`: PASS
- `linear_residual`: PASS
- `boundary_residual`: PASS
- `complex_l2`: PASS
- `phase_rms`: PASS
- `tl_rms_reported`: PASS
- `rho_shape_reported`: PASS
- `spatial_convergence`: PASS
- `window_convergence`: PASS
- `finite`: PASS
- `all`: PASS

Artifacts: `results\validation\pe_bellhop_helmholtz_bie_reference\R1_flat\R1_flat_validation.mat`, `results\validation\pe_bellhop_helmholtz_bie_reference\R1_flat\R1_flat_convergence.csv`.

Rough stages remain locked unless every R1 hard gate passes.
