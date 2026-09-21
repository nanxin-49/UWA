# PE--Bellhop--Helmholtz BIE G0_high_K

日期：2026-09-15  
状态：**BELLHOP_CLOSER_TO_HELMHOLTZ_REFERENCE**

- Surface: `A=0.02 m`, `K=0.47 rad/m`, same fixed C2 taper.
- U_BIE: `6.01391724e-08`; required 5U separation: `3.00695862e-07`; observed: `0.0227692636`.

| pair | E_G | phase RMS | TL RMS dB | rho | phi0 | E_aligned |
|---|---:|---:|---:|---:|---:|---:|
| PE_BIE | 0.0229060754 | 0.0130148657 | 0.163731694 | 0.999776097 | -0.00877050329 | 0.021161953 |
| BH_BIE | 0.000136811849 | 9.770212e-05 | 0.000831167008 | 0.999999991 | -1.27009726e-06 | 0.000136805022 |
| PE_BH | 0.0228902995 | 0.0129918067 | 0.163704445 | 0.99977646 | -0.00877091285 | 0.0211448602 |
| BH_native_BIE_native | 0.000136811849 | 9.770212e-05 | 0.000831167008 | 0.999999991 | 1.27009726e-06 | 0.000136805022 |

## Refinement stability

| kind | parameter | E_PE,BIE | E_BH,BIE | separation | winner |
|---|---:|---:|---:|---:|---|
| spatial | 8 | 0.0229060739 | 0.000136811982 | 0.0227692619 | BELLHOP |
| spatial | 10 | 0.0229060748 | 0.000136812208 | 0.0227692626 | BELLHOP |
| spatial | 12 | 0.0229060754 | 0.000136811849 | 0.0227692636 | BELLHOP |
| window | 60 | 0.0229060748 | 0.000136812208 | 0.0227692626 | BELLHOP |
| window | 70 | 0.0229060748 | 0.000136806941 | 0.0227692679 | BELLHOP |
| window | 80 | 0.0229060749 | 0.000136805177 | 0.0227692698 | BELLHOP |

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

Artifacts: `results\validation\pe_surface_operator_bie_reference\G0_high_K\G0_high_K_validation.mat`, `results\validation\pe_surface_operator_bie_reference\G0_high_K\G0_high_K_refinement_stability.csv`.
