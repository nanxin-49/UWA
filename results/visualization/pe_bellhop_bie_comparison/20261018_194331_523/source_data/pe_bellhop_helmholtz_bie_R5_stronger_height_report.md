# PE--Bellhop--Helmholtz BIE R5_stronger_height

日期：2026-09-15  
状态：**BELLHOP_CLOSER_TO_HELMHOLTZ_REFERENCE**

- Surface: `A=0.2 m`, `K=0.1 rad/m`, same fixed C2 taper.
- U_BIE: `6.17857235e-08`; required 5U separation: `3.08928617e-07`; observed: `0.0808083477`.

- Strong-case boundary cross-check: residual `2.87565501e-08`, main/high-order receiver-field L2 `3.29044637e-08`, acceptance limit `6.17857235e-08`.
- This is an explicit strong-height numerical limit; the R1 analytic boundary gate remains unchanged at `1e-8`.

| pair | E_G | phase RMS | TL RMS dB | rho | phi0 | E_aligned |
|---|---:|---:|---:|---:|---:|---:|
| PE_BIE | 0.0809081077 | 0.0803363527 | 0.0869963866 | 0.997578969 | -0.041348938 | 0.0695923116 |
| BH_BIE | 9.97600752e-05 | 7.41103188e-05 | 0.000579793191 | 0.999999995 | 7.33358644e-08 | 9.97600249e-05 |
| PE_BH | 0.0809227918 | 0.0803511394 | 0.0870017811 | 0.997577772 | -0.0413487689 | 0.0696095237 |
| BH_native_BIE_native | 9.97600752e-05 | 7.41103188e-05 | 0.000579793191 | 0.999999995 | -7.33358644e-08 | 9.97600249e-05 |

## Refinement stability

| kind | parameter | E_PE,BIE | E_BH,BIE | separation | winner |
|---|---:|---:|---:|---:|---|
| spatial | 8 | 0.0809081051 | 9.97602857e-05 | 0.0808083448 | BELLHOP |
| spatial | 10 | 0.0809081066 | 9.97602446e-05 | 0.0808083463 | BELLHOP |
| spatial | 12 | 0.0809081073 | 9.97601785e-05 | 0.0808083471 | BELLHOP |
| spatial | 14 | 0.0809081077 | 9.97600752e-05 | 0.0808083477 | BELLHOP |
| window | 60 | 0.0809081066 | 9.97602446e-05 | 0.0808083463 | BELLHOP |
| window | 70 | 0.0809081047 | 9.97570547e-05 | 0.0808083477 | BELLHOP |
| window | 80 | 0.0809081043 | 9.97566449e-05 | 0.0808083476 | BELLHOP |

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

Artifacts: `results\validation\pe_bellhop_helmholtz_bie_reference\R5_stronger_height\R5_stronger_height_validation.mat`, `results\validation\pe_bellhop_helmholtz_bie_reference\R5_stronger_height\R5_stronger_height_refinement_stability.csv`.
