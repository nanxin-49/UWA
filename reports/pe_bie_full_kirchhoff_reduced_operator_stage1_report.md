# Full Kirchhoff to reduced-operator Stage-1 analysis

Status: **validation-only diagnostic; production PE/BIE/surface models unchanged**.

Configuration and convention are inherited from the accepted 4 kHz deterministic Full-Kirchhoff sweep. Metrics use a fixed Stage-0 M99 footprint, a Full-Kirchhoff-reference -40 dB mask, and saved incident-energy weights. No scalar fitting or case-dependent phase adjustment is used.

## PE phase screen versus Full Kirchhoff

| case | complex L2 | magnitude L2 | phase RMS rad | magnitude corr. | phase corr. | error r95 m | corr(error,|eta|) | corr(error,|slope|) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| weak_low_K | 0.00316587 | 0.00053118 | 0.003121 | 0.797726 | 0.999995 | 38.6719 | 0.271335 | -0.119443 |
| strong_height_low_K | 0.0808828 | 0.010018 | 0.0803104 | 0.977034 | 0.997631 | 37.8906 | 0.00795687 | 0.10146 |
| weak_high_K | 0.0228908 | 0.0188542 | 0.0129838 | 0.998529 | 0.999954 | 33.0078 | 0.524647 | -0.510576 |

## Physical-term ablations relative to Full Kirchhoff

| case | model | complex L2 | magnitude L2 | phase RMS rad | magnitude corr. | phase corr. |
|---|---|---:|---:|---:|---:|---:|
| weak_low_K | PE phase screen | 0.00316587 | 0.00053118 | 0.003121 | 0.797726 | 0.999995 |
| weak_low_K | flat normal | 2.12279e-06 | 3.12569e-07 | 2.09965e-06 | 1 | 1 |
| weak_low_K | ds=dx | 3.12505e-07 | 3.12505e-07 | 7.27543e-10 | 1 | 1 |
| weak_low_K | flat surface z=0 | 0.230058 | 0.000252264 | 0.230841 | 0.567778 | 0.973537 |
| strong_height_low_K | PE phase screen | 0.0808828 | 0.010018 | 0.0803104 | 0.977034 | 0.997631 |
| strong_height_low_K | flat normal | 0.000131859 | 0.000124992 | 4.19931e-05 | 0.9999 | 1 |
| strong_height_low_K | ds=dx | 0.000124987 | 0.000124987 | 2.91442e-07 | 0.9999 | 1 |
| strong_height_low_K | flat surface z=0 | 1.21726 | 0.00500692 | 1.55247 | 0.00400517 | 0.258933 |
| weak_high_K | PE phase screen | 0.0228908 | 0.0188542 | 0.0129838 | 0.998529 | 0.999954 |
| weak_high_K | flat normal | 9.80551e-05 | 2.7069e-05 | 9.42267e-05 | 0.999999 | 1 |
| weak_high_K | ds=dx | 2.70275e-05 | 2.69988e-05 | 1.24501e-06 | 0.999999 | 1 |
| weak_high_K | flat surface z=0 | 0.463521 | 0.00938391 | 0.469871 | 0.0168444 | 0.892627 |

`flat normal` retains the actual surface position and arc-length measure but consistently uses `n=(0,1)` in both Green and incident normal derivatives. `ds=dx` retains the curved normal and actual surface position but removes only the arc-length Jacobian. `flat surface z=0` removes the physical surface height/path perturbation and recomputes the incident field consistently at z=0.

## Receiver-centered finite-aperture integral

| case | L m | complex L2 | magnitude L2 | phase RMS rad | correction relative error |
|---|---:|---:|---:|---:|---:|
| weak_low_K | 0.5 | 0.493482 | 0.316477 | 0.494674 | 155.875 |
| weak_low_K | 1 | 0.280112 | 0.209487 | 0.199399 | 88.4784 |
| weak_low_K | 2 | 0.131048 | 0.078735 | 0.104136 | 41.3939 |
| weak_low_K | 4 | 0.0642782 | 0.032374 | 0.0555029 | 20.3035 |
| weak_low_K | 8 | 0.0348453 | 0.0191176 | 0.0291297 | 11.0065 |
| weak_low_K | 16 | 0.0204334 | 0.015433 | 0.0133926 | 6.45427 |
| weak_low_K | 32 | 0.010858 | 0.00768511 | 0.0076721 | 3.42971 |
| weak_low_K | 64 | 0.00153326 | 0.00108533 | 0.00108319 | 0.484309 |
| weak_low_K | Inf | 0 | 0 | 0 | 0 |
| strong_height_low_K | 0.5 | 0.483117 | 0.320938 | 0.482086 | 5.97306 |
| strong_height_low_K | 1 | 0.282456 | 0.207716 | 0.206197 | 3.49216 |
| strong_height_low_K | 2 | 0.13169 | 0.0749309 | 0.108492 | 1.62816 |
| strong_height_low_K | 4 | 0.0644481 | 0.0513867 | 0.0389942 | 0.796809 |
| strong_height_low_K | 8 | 0.0348648 | 0.0235841 | 0.0256445 | 0.431054 |
| strong_height_low_K | 16 | 0.0204485 | 0.0146587 | 0.0142559 | 0.252817 |
| strong_height_low_K | 32 | 0.0108708 | 0.00768387 | 0.00769771 | 0.134401 |
| strong_height_low_K | 64 | 0.00153487 | 0.00109041 | 0.00107974 | 0.0189765 |
| strong_height_low_K | Inf | 0 | 0 | 0 | 0 |
| weak_high_K | 0.5 | 0.494546 | 0.317989 | 0.496678 | 21.6045 |
| weak_high_K | 1 | 0.280766 | 0.210254 | 0.199652 | 12.2654 |
| weak_high_K | 2 | 0.131172 | 0.079685 | 0.103468 | 5.73032 |
| weak_high_K | 4 | 0.0643382 | 0.0394571 | 0.0508495 | 2.81065 |
| weak_high_K | 8 | 0.0348525 | 0.0237707 | 0.0254552 | 1.52256 |
| weak_high_K | 16 | 0.0204366 | 0.0152283 | 0.0135963 | 0.892784 |
| weak_high_K | 32 | 0.0108599 | 0.0076968 | 0.00766505 | 0.474424 |
| weak_high_K | 64 | 0.00153339 | 0.00107986 | 0.0010894 | 0.0669871 |
| weak_high_K | Inf | 0 | 0 | 0 | 0 |

The finite-aperture model keeps all accepted Full-Kirchhoff terms but sets contributions with `|x_r-x_s|>L` to zero. Its correction diagnostic compares `(u_L-u_PS)` with the complete `(u_FK-u_PS)`; it is not a fitted convolution kernel.

## Interpretation and ranking

- **weak_low_K:** PS/FK `Ec=0.00316587`; flat-normal `Ec=2.12279e-06`; no-Jacobian `Ec=3.12505e-07`; first aperture below half the PS error: **64 m**.
- **strong_height_low_K:** PS/FK `Ec=0.0808828`; flat-normal `Ec=0.000131859`; no-Jacobian `Ec=0.000124987`; first aperture below half the PS error: **8 m**.
- **weak_high_K:** PS/FK `Ec=0.0228908`; flat-normal `Ec=9.80551e-05`; no-Jacobian `Ec=2.70275e-05`; first aperture below half the PS error: **32 m**.

The physical importance ordering is determined from the tables rather than assumed. A term whose ablation error is far below PS/FK cannot explain the phase-screen failure by itself. The aperture needed to beat the phase screen measures the minimum nonlocal support for a reduced operator under this receiver-centered truncation.

The flat-surface ablation shows that physical surface height/path phase is essential, but it does **not** mean PE omits height phase: PE already contains the local `2*k*eta` proxy. The residual diagnosis is that this local proxy does not reproduce the complete surface-to-receiver path phase and nonlocal coherent integration. Curved-normal and Jacobian corrections are too small to close that residual on their own.

A reduced short-range Kirchhoff integral is the preferred next prototype if a finite `L` consistently beats the phase screen across all three cases. A phase-screen-plus-small-kernel form should only be implemented after extracting a case-independent kernel from this aperture study; directly inserting the measured `u_FK-u_PS` would be tautological. SSA remains out of scope because complete Full Kirchhoff already agrees with BIE in the tested envelope.

Artifacts: `results/validation/pe_bie_full_kirchhoff_reduced_stage1/pe_bie_full_kirchhoff_reduced_stage1.mat`, `results/validation/pe_bie_full_kirchhoff_reduced_stage1/pe_bie_full_kirchhoff_reduced_stage1_summary.csv`, `results/validation/pe_bie_full_kirchhoff_reduced_stage1/pe_bie_full_kirchhoff_reduced_stage1_radius_scan.csv`, `results/validation/pe_bie_full_kirchhoff_reduced_stage1/pe_bie_full_kirchhoff_reduced_stage1_sweep.csv`; figures: `results/validation/pe_bie_full_kirchhoff_reduced_stage1/figures/`.
