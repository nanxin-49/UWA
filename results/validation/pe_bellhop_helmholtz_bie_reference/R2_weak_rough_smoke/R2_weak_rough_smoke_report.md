# PE--Bellhop--Helmholtz BIE R2 weak-rough convergence

日期：2026-09-14  
状态：**R2 FAIL / R3 LOCKED**

## Frozen benchmark

- `eta_bench(x)=A sin(Kx) chi(x)`, `A=0.01 m`, `K=0.1 rad/m`.
- `chi=1` for `|x|<=42 m`; fixed physical support ends at `|x|=50 m`.
- This physical taper is fixed in every run. The BIE window is a separate numerical window outside that support.
- source: `discrete-periodic-Gaussian|N=984|W=192.1875|sigma=0.3|f=4000|c=1500|exp(-iwt)`

## Numerical uncertainty

- `U_BIE` complex L2 M99: `0.233955219`
- phase RMS M99 uncertainty: `0.182466718 rad`
- TL RMS M95 uncertainty: `1.5453902 dB`

| sweep | parameter | nodes | boundary residual | linear residual | L2 to reference | successive L2 | successive phase | successive TL dB |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| points_per_wavelength_actual | 2 | 640 | 0.185 | 1.55e-15 | 0.281 | NaN | NaN | NaN |
| points_per_wavelength_actual | 3 | 960 | 0.00836 | 2.01e-15 | 0.00528 | 0.269 | 0.204 | 1.63 |
| points_per_wavelength_actual | 4 | 1280 | 0.00333 | 2.12e-15 | 0 | 0.00528 | 0.00327 | 0.0372 |
| self_quadrature_order | 16 | 1280 | 0.00333 | 2.07e-15 | 1.26e-05 | NaN | NaN | NaN |
| self_quadrature_order | 32 | 1280 | 0.00333 | 2.12e-15 | 0 | 1.26e-05 | 1.26e-05 | 6.78e-07 |
| half_width_m | 60 | 640 | 0.185 | 1.55e-15 | 0.232 | NaN | NaN | NaN |
| half_width_m | 65 | 696 | 0.18 | 1.87e-15 | 0.234 | 0.187 | 0.143 | 1.33 |
| half_width_m | 70 | 748 | 0.194 | 1.91e-15 | 0 | 0.234 | 0.182 | 1.55 |

## Gates

- `source_dft`: PASS
- `linear_residual`: PASS
- `boundary_residual`: FAIL
- `spatial_convergence`: FAIL
- `quadrature_convergence`: FAIL
- `window_convergence`: FAIL
- `finite`: PASS
- `all`: FAIL

Artifacts: `C:\Users\ASUS\.codex\worktrees\baf5\Explain\results\validation\pe_bellhop_helmholtz_bie_reference\R2_weak_rough_smoke\R2_weak_rough_convergence.mat`, `C:\Users\ASUS\.codex\worktrees\baf5\Explain\results\validation\pe_bellhop_helmholtz_bie_reference\R2_weak_rough_smoke\R2_weak_rough_convergence.csv`.

R3 remains locked unless every R2 gate passes.
