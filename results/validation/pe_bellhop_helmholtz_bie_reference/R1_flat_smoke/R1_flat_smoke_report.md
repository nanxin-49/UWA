# PE--Bellhop--Helmholtz BIE R1 flat validation

日期：2026-09-14  
状态：**R1 FAIL / ROUGH STAGES LOCKED**

- source: `discrete-periodic-Gaussian|N=984|W=192.1875|sigma=0.3|f=4000|c=1500|exp(-iwt)`
- formulation: `Dirichlet half-plane Green combined layer; smooth finite section`
- convention/normal/jump: `exp(-i*omega*t)`; `into water domain z>eta(x)`; `gamma_D D_h = +1/2 I + K_h`

- M95/M99 radii: `29.4921875 / 40.0390625 m`

| sweep | parameter | nodes | boundary residual | linear residual | image L2 M99 | phase RMS M99 | TL RMS M95 (dB) | rho | successive L2 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| points_per_wavelength_actual | 2 | 640 | 1.07 | 6.52e-15 | 0.617 | 0.456 | 3.31 | 0.853089476 | NaN |
| points_per_wavelength_actual | 3 | 960 | 0.0532 | 1.08e-14 | 0.0127 | 0.00918 | 0.0453 | 0.999919585 | 0.523 |
| points_per_wavelength_actual | 4 | 1280 | 18.9 | 5.9e-11 | 0.767 | 0.219 | 1.33 | 0.79248577 | 0.768 |
| half_width_m | 60 | 640 | 1.07 | 6.52e-15 | 0.617 | 0.456 | 3.31 | 0.853089476 | NaN |
| half_width_m | 65 | 696 | 1.04 | 8.65e-15 | 0.532 | 0.379 | 2.99 | 0.885550576 | 0.581 |
| half_width_m | 70 | 748 | 1.05 | 8.72e-15 | 0.514 | 0.399 | 3.29 | 0.891761698 | 0.576 |

## Hard gates

- `source_dft`: PASS
- `linear_residual`: PASS
- `boundary_residual`: FAIL
- `complex_l2`: FAIL
- `phase_rms`: FAIL
- `tl_rms_reported`: PASS
- `rho_shape_reported`: PASS
- `spatial_convergence`: FAIL
- `window_convergence`: FAIL
- `finite`: PASS
- `all`: FAIL

Artifacts: `C:\Users\ASUS\.codex\worktrees\baf5\Explain\results\validation\pe_bellhop_helmholtz_bie_reference\R1_flat_smoke\R1_flat_validation.mat`, `C:\Users\ASUS\.codex\worktrees\baf5\Explain\results\validation\pe_bellhop_helmholtz_bie_reference\R1_flat_smoke\R1_flat_convergence.csv`.

Rough stages remain locked unless every R1 hard gate passes.
