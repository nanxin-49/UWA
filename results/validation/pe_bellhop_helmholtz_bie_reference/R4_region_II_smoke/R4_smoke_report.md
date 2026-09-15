# PE--Bellhop--Helmholtz BIE R4 Region-II adjudication

日期：2026-09-15  
状态：**REFERENCE_NOT_YET_DISCRIMINATING**

- Surface: `A=0.05 m`, `K=0.1 rad/m`, same fixed C2 taper.
- U_BIE: `0.235278479`; required 5U separation: `1.1763924`; observed: `0.0140190749`.

| pair | E_G | phase RMS | TL RMS dB | rho | phi0 | E_aligned |
|---|---:|---:|---:|---:|---:|---:|
| PE_BIE | 0.0162596279 | 0.0159922136 | 0.025493741 | 0.999871357 | -0.00266515401 | 0.0160407733 |
| BH_BIE | 0.00224055295 | 0.00163829083 | 0.0132766748 | 0.999997495 | -8.30421029e-05 | 0.00223912684 |
| PE_BH | 0.0160946927 | 0.0158993331 | 0.0217393071 | 0.999873812 | -0.0025822031 | 0.0158864322 |
| BH_native_BIE_native | 0.00224055295 | 0.00163829083 | 0.0132766748 | 0.999997495 | 8.30421029e-05 | 0.00223912684 |

## Refinement stability

| kind | parameter | E_PE,BIE | E_BH,BIE | separation | winner |
|---|---:|---:|---:|---:|---|
| spatial | 2 | 0.271521346 | 0.271114336 | 0.000407010148 | BELLHOP |
| spatial | 3 | 0.0166992715 | 0.00470767182 | 0.0119915996 | BELLHOP |
| spatial | 4 | 0.0162596279 | 0.00224055295 | 0.0140190749 | BELLHOP |
| window | 60 | 0.271521346 | 0.271114336 | 0.000407010148 | BELLHOP |
| window | 65 | 0.266953467 | 0.266541612 | 0.00041185528 | BELLHOP |
| window | 70 | 0.269243175 | 0.268828786 | 0.00041438909 | BELLHOP |

## Gates

- `source_dft`: PASS
- `boundary_residual`: FAIL
- `spatial_convergence`: FAIL
- `window_convergence`: FAIL
- `bellhop_geometry`: PASS
- `ranking_stable`: PASS
- `five_U_separation`: FAIL
- `finite`: PASS
- `all`: FAIL

Artifacts: `C:\Users\ASUS\.codex\worktrees\baf5\Explain\results\validation\pe_bellhop_helmholtz_bie_reference\R4_region_II_smoke\R4_region_II_validation.mat`, `C:\Users\ASUS\.codex\worktrees\baf5\Explain\results\validation\pe_bellhop_helmholtz_bie_reference\R4_region_II_smoke\R4_refinement_stability.csv`.
