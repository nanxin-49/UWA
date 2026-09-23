# Full Kirchhoff operator visualization atlas

Status: **post-processing only**. No PE, BIE, Full-Kirchhoff, or sweep computation was rerun.

The plotted scalar density is `I(x_r,x_s)=-[G(x_r,x_s) d_n p_inc(x_s)+d_n_s G(x_r,x_s) p_inc(x_s)] w_s`. It contains the actual surface position, unit normal, source-normal Green derivative, incident normal derivative, and saved arc-length quadrature weight. It is therefore the correct visual object for this implementation; `G` alone is not called the Full-Kirchhoff kernel.

## Reconstruction gate

| case | representative receiver x m | mask samples | field reconstruction relative L2 | selected-row relative error | coherent cancellation ratio |
|---|---:|---:|---:|---:|---:|
| weak_low_K | 19.9219 | 411 | 0 | 0 | 22.8237 |
| strong_height_low_K | 19.9219 | 411 | 0 | 0 | 22.8455 |
| weak_high_K | 19.9219 | 411 | 0 | 0 | 22.8708 |

The reconstruction gate sums the plotted density over the saved surface samples and compares it with the saved native Full-Kirchhoff receiver field before applying the frozen comparison conjugation/normalization.

Figures and summary data: `results/visualization/pe_bie_full_kirchhoff_operator_atlas`.
