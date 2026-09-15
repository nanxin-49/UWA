# PE--Bellhop--Helmholtz BIE R4 Region-II adjudication

日期：2026-09-15  
状态：**BELLHOP_CLOSER_TO_HELMHOLTZ_REFERENCE**

- Surface: `A=0.05 m`, `K=0.1 rad/m`, same fixed C2 taper.
- U_BIE: `6.02438379e-08`; required 5U separation: `3.0121919e-07`; observed: `0.0160510649`.

| pair | E_G | phase RMS | TL RMS dB | rho | phi0 | E_aligned |
|---|---:|---:|---:|---:|---:|---:|
| PE_BIE | 0.0160899604 | 0.0158944947 | 0.021741798 | 0.999873888 | -0.00258212044 | 0.0158816508 |
| BH_BIE | 3.88955103e-05 | 2.79146697e-05 | 0.000235258791 | 0.999999999 | 9.77517499e-08 | 3.88953866e-05 |
| PE_BH | 0.0160946927 | 0.0158993331 | 0.0217393071 | 0.999873812 | -0.0025822031 | 0.0158864322 |
| BH_native_BIE_native | 3.88955103e-05 | 2.79146697e-05 | 0.000235258791 | 0.999999999 | -9.77517499e-08 | 3.88953866e-05 |

## Refinement stability

| kind | parameter | E_PE,BIE | E_BH,BIE | separation | winner |
|---|---:|---:|---:|---:|---|
| spatial | 8 | 0.0160899597 | 3.88955559e-05 | 0.0160510641 | BELLHOP |
| spatial | 10 | 0.0160899601 | 3.88956183e-05 | 0.0160510645 | BELLHOP |
| spatial | 12 | 0.0160899604 | 3.88955103e-05 | 0.0160510649 | BELLHOP |
| window | 60 | 0.0160899601 | 3.88956183e-05 | 0.0160510645 | BELLHOP |
| window | 70 | 0.0160899621 | 3.88948896e-05 | 0.0160510672 | BELLHOP |
| window | 80 | 0.0160899625 | 3.8894381e-05 | 0.0160510681 | BELLHOP |

## Gates

- `source_dft`: PASS
- `boundary_residual`: PASS
- `spatial_convergence`: PASS
- `window_convergence`: PASS
- `bellhop_geometry`: PASS
- `ranking_stable`: PASS
- `five_U_separation`: PASS
- `finite`: PASS
- `all`: PASS

Artifacts: `C:\Users\ASUS\.codex\worktrees\baf5\Explain\results\validation\pe_bellhop_helmholtz_bie_reference\R4_region_II\R4_region_II_validation.mat`, `C:\Users\ASUS\.codex\worktrees\baf5\Explain\results\validation\pe_bellhop_helmholtz_bie_reference\R4_region_II\R4_refinement_stability.csv`.
