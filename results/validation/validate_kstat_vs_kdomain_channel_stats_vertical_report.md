# K-Stat vs K-Domain Channel Statistics Validation

Generated: 2026-07-06 14:49:45

## Setup

- Full `vertical_channel_model` path with PE/WAPE propagation.
- No modulation, noise, demodulation, or BER/SER loop.
- Hs values: `[0.1 0.2]` m; seed count: `8`; grid `64 x 64`; Nf `8`.
- Wideband band: `[4000 8000]` Hz; reference frequency: `6000.000` Hz.
- Roughness mode: `target_hs`; PM shape wind: `5 m/s`.

## Key Results

- Max H invariant error, kdomain: `7.758e-18`.
- Max H invariant error, kstat: `7.758e-18`.
- Max K-Stat phase-screen energy error: `9.992e-16`.
- Max direct-path branch/seed delta: `0`.
- Max E[|h_ref|^2] relative error: `0.09729`.
- Max reflected mean-PDP L2 error: `0.3134`.
- Min reflected mean-PDP correlation: `0.9793`.
- Max |h_total| quantile NRMSE: `0.1932`.
- Hard checks passed: `1`.
- Statistical targets passed: `0`.

## Interpretation

- `kirchhoff_kdomain` is an ensemble of explicit sea-surface phase-screen realizations.
- `kirchhoff_kstat` is an ensemble of statistical phase-screen reflected-field realizations.
- Single-seed `h_ref`, `h_total`, or PDP samples are not expected to match pointwise.
- Valid comparison targets are `E[|h_ref|^2]`, mean PDP shape, `|h_total|` distribution, and implementation invariants.
- PDP is a wideband receiver-side statistic. The current kstat frequency correlation is still an engineering simplification, so PDP agreement is an implementation-level statistical check, not a calibrated sea-surface time-frequency model.

## Figures

- `E:\MISC\CARPE3D_matlab\Explain\results\validation\kstat_kdomain_E_abs_h_reflect2_vs_Hs.png`
- `E:\MISC\CARPE3D_matlab\Explain\results\validation\kstat_kdomain_abs_h_total_distribution.png`
- `E:\MISC\CARPE3D_matlab\Explain\results\validation\kstat_kdomain_mean_pdp_reflect_compare.png`
- `E:\MISC\CARPE3D_matlab\Explain\results\validation\kstat_kdomain_mean_pdp_total_compare.png`
- `E:\MISC\CARPE3D_matlab\Explain\results\validation\kstat_kdomain_pdp_error_vs_Hs.png`

## Summary Table

    Hs_m    seed_count    E_abs_h_reflect2_rel_error    pdp_reflect_l2_error    pdp_reflect_corr    pdp_total_l2_error    pdp_total_corr    abs_h_total_quantile_nrmse    max_kstat_energy_error
    ____    __________    __________________________    ____________________    ________________    __________________    ______________    __________________________    ______________________
    0.1         8                  0.082146                   0.15693               0.99638               0.13749            0.99729                 0.081205                    9.992e-16      
    0.2         8                  0.097289                   0.31338               0.97932              0.094821            0.99878                  0.19316                   2.2204e-16      

