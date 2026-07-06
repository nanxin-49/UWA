# K-Stat vs K-Domain Phase-Screen Validation

Generated: 2026-07-03 12:09:08

## Setup

- Surface-boundary-only validation; no PE/WAPE and no communication chain.
- Frequencies: `[4000 6000 8000]` Hz.
- Hs values: `[0.05 0.1 0.2 0.5]` m.
- M values: `[8 16 32 64]`; grid `128 x 128`; PM shape wind `5 m/s`.

## Key Results

- Max K-Stat energy closure error: `4.663e-15`.
- Max K-Domain phase-screen power closure error: `0`.
- Max sigma_eta relative error: `4.163e-16`.
- Final-M mean |Rcoh| relative error where kstat |Rcoh| > 1.0e-03: `0.07082`.
- Final-M mean |Rcoh| absolute error: `0.003294`.
- Final-M mean radial spectrum L2 error: `0.0213`.
- Final-M mean radial spectrum correlation: `0.9997`.

## Interpretation

- kdomain single realizations are concrete phase screens; kstat is a statistical phase-screen model.
- Single-seed complex fields are not expected to match pointwise.
- Valid comparisons are ensemble coherent mean, incoherent power spectrum, radial spectral shape, energy closure, and statistical trends.
- If PE/WAPE propagation is added later, compare receiver-side statistics such as `E[|h_ref|^2]`, `std(|h_ref|)`, or average PDP, not single-realization `h` equality.

## Figures

- `E:\MISC\CARPE3D_matlab\Explain\results\validation\kstat_kdomain_Rcoh_vs_Hs_f.png`
- `E:\MISC\CARPE3D_matlab\Explain\results\validation\kstat_kdomain_incoh_energy_vs_Hs_f.png`
- `E:\MISC\CARPE3D_matlab\Explain\results\validation\kstat_kdomain_convergence_vs_M.png`
- `E:\MISC\CARPE3D_matlab\Explain\results\validation\kstat_kdomain_radial_spectrum_compare.png`
- `E:\MISC\CARPE3D_matlab\Explain\results\validation\kstat_kdomain_2d_spectrum_typical.png`

## Summary Table Preview

    f_hz    Hs_m    M     Rcoh_kdomain_abs    Rcoh_kstat_abs    incoh_energy_kdomain    incoh_energy_kstat    radial_spectrum_error    energy_closure_error_kstat
    ____    ____    __    ________________    ______________    ____________________    __________________    _____________________    __________________________
    4000    0.05    64          0.9159             0.91601           4.2597e+07             4.3199e+07              0.023768                   2.2204e-16        
    6000    0.05    64         0.82036             0.82087            8.642e+07             8.7557e+07              0.020945                   2.2204e-16        
    8000    0.05    64         0.70263             0.70404           1.3375e+08             1.3538e+08              0.018173                   4.4409e-16        
    4000     0.1    64         0.70263             0.70404           1.3375e+08             1.3538e+08              0.018173                   4.4409e-16        
    6000     0.1    64         0.44916             0.45404           2.1075e+08              2.131e+08              0.016397                   1.7764e-15        
    8000     0.1    64         0.23675             0.24569           2.4923e+08             2.5223e+08               0.02019                   1.3323e-15        
    4000     0.2    64         0.23675             0.24569           2.4923e+08             2.5223e+08               0.02019                   1.3323e-15        
    6000     0.2    64        0.033967            0.042499           2.6397e+08             2.6795e+08               0.02536                   6.6613e-16        
    8000     0.2    64       0.0023746           0.0036439           2.6426e+08             2.6843e+08              0.025216                   4.6629e-15        
    4000     0.5    64        0.002207          0.00015486           2.6419e+08             2.6844e+08              0.025688                   3.6637e-15        
    6000     0.5    64       0.0003513          2.6753e-09           2.6424e+08             2.6844e+08              0.024458                   6.6613e-16        
    8000     0.5    64       0.0011319          5.7514e-16           2.6421e+08             2.6844e+08              0.016995                   1.1102e-16        

