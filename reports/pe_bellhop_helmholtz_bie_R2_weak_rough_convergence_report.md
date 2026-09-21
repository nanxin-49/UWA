# PE--Bellhop--Helmholtz BIE R2 weak-rough convergence

日期：2026-09-14  
状态：**R2 PASS**

## Frozen benchmark

- `eta_bench(x)=A sin(Kx) chi(x)`, `A=0.01 m`, `K=0.1 rad/m`.
- `chi=1` for `|x|<=42 m`; fixed physical support ends at `|x|=50 m`.
- This physical taper is fixed in every run. The BIE window is a separate numerical window outside that support.
- source: `discrete-periodic-Gaussian|N=984|W=192.1875|sigma=0.3|f=4000|c=1500|exp(-iwt)`

## Numerical uncertainty

- `U_BIE` complex L2 M99: `6.00190522e-08`
- phase RMS M99 uncertainty: `4.21974014e-08 rad`
- TL RMS M95 uncertainty: `3.59585756e-07 dB`

| sweep | parameter | nodes | boundary residual | linear residual | L2 to reference | successive L2 | successive phase | successive TL dB |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| points_per_wavelength_actual | 8 | 2560 | 4.74e-08 | 2.88e-15 | 3.94e-09 | NaN | NaN | NaN |
| points_per_wavelength_actual | 10 | 3200 | 7.2e-09 | 3.07e-15 | 1.67e-09 | 2.48e-09 | 2.38e-09 | 7.26e-09 |
| points_per_wavelength_actual | 12 | 3840 | 1.92e-09 | 3.43e-15 | 0 | 1.67e-09 | 1.59e-09 | 4.71e-09 |
| self_quadrature_order | 64 | 3840 | 1.83e-09 | 3.3e-15 | 3.08e-08 | NaN | NaN | NaN |
| self_quadrature_order | 96 | 3840 | 1.92e-09 | 3.43e-15 | 0 | 3.08e-08 | 3.08e-08 | 3.63e-09 |
| half_width_m | 60 | 3200 | 7.2e-09 | 3.07e-15 | 2.84e-07 | NaN | NaN | NaN |
| half_width_m | 70 | 3736 | 1.19e-08 | 3.43e-15 | 6e-08 | 2.33e-07 | 1.65e-07 | 1.22e-06 |
| half_width_m | 80 | 4272 | 1.54e-08 | 3.68e-15 | 0 | 6e-08 | 4.22e-08 | 3.6e-07 |

## Gates

- `source_dft`: PASS
- `linear_residual`: PASS
- `boundary_residual`: PASS
- `spatial_convergence`: PASS
- `quadrature_convergence`: PASS
- `window_convergence`: PASS
- `finite`: PASS
- `all`: PASS

Artifacts: `results\validation\pe_bellhop_helmholtz_bie_reference\R2_weak_rough\R2_weak_rough_convergence.mat`, `results\validation\pe_bellhop_helmholtz_bie_reference\R2_weak_rough\R2_weak_rough_convergence.csv`.

R3 remains locked unless every R2 gate passes.
