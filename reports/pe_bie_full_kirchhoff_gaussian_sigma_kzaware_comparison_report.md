# Gaussian source angular-spectrum sensitivity with kz-aware PE

Status: validation-only. Production PE, Full Kirchhoff, BIE, surface model, masks, receiver convention and physical parameters were not modified.

Setup: 2D one-transverse, 4 kHz, c=1500 m/s, Tx=100 m, Rx=3 m, sigma=[0.3, 0.5, 1, 2, 5] m. Cases: weak_low_K (A=0.05, K=0.10), strong_height_low_K (A=0.20, K=0.10), weak_high_K (A=0.02, K=0.47). Model-1 is the existing `model1_kz_aware` PE option.

Model-0, Full Kirchhoff and BIE fields are reused from the accepted sigma sensitivity MAT. Existing strict-normal Model-1 MAT files are reused when present; missing Model-1 cases are computed only with `run_pe_1d_surface_reflection_validation`, using the same grid, step, source and surface parameters. Reused=8, computed=7.

Mask: common sigma-specific 99-percent source-energy footprint from the base sigma sweep, followed by a pairwise -40 dB threshold relative to the larger compared-field peak; energy weights are renormalized after masking.

## Master complex-error table

|case|sigma|theta99 (deg)|Model0-BIE Ec|Model1-BIE Ec|FK-BIE Ec|Model0-FK Ec|Model1-FK Ec|M0/M1|M1/M0|
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
|weak_low_K|0.3|21.2801|0.01609|0.00398989|0.000253723|0.0160829|0.00399757|4.03269|0.247974|
|weak_low_K|0.5|12.5086|0.013855|0.00493999|1.10164e-05|0.0138462|0.00494001|2.80466|0.356549|
|weak_low_K|1|6.27314|0.0101692|0.0102641|8.20441e-06|0.0101669|0.0102641|0.990755|1.00933|
|weak_low_K|2|3.13187|0.0169587|0.0172961|5.37393e-06|0.0169587|0.0172961|0.980494|1.01989|
|weak_low_K|5|1.22986|0.0113287|0.0113617|7.24358e-06|0.0113288|0.0113617|0.997103|1.00291|
|strong_height_low_K|0.3|21.2801|0.0809081|0.0507808|0.000257278|0.0808828|0.0507816|1.59328|0.627636|
|strong_height_low_K|0.5|12.5086|0.0728833|0.0510946|4.40203e-05|0.0728567|0.0510946|1.42644|0.701047|
|strong_height_low_K|1|6.27314|0.0713342|0.0713033|3.28167e-05|0.0713291|0.0713035|1.00043|0.999567|
|strong_height_low_K|2|3.13187|0.0959002|0.0968171|2.14855e-05|0.0959004|0.0968176|0.990529|1.00956|
|strong_height_low_K|5|1.22986|0.0763957|0.076478|2.89695e-05|0.076396|0.0764783|0.998924|1.00108|
|weak_high_K|0.3|21.2801|0.0229061|0.0215668|0.000271135|0.0228908|0.0215714|1.0621|0.941531|
|weak_high_K|0.5|12.5086|0.0225558|0.0222385|9.40723e-05|0.0225479|0.0222416|1.01427|0.985933|
|weak_high_K|1|6.27314|0.0251597|0.0252671|9.67738e-05|0.0251589|0.0252701|0.995753|1.00427|
|weak_high_K|2|3.13187|0.0281717|0.028124|9.41935e-05|0.0281746|0.0281274|1.0017|0.998306|
|weak_high_K|5|1.22986|0.0288541|0.028863|9.07573e-05|0.0288575|0.0288663|0.999693|1.00031|

## Phase and magnitude summary

|case|sigma|M0 mag|M1 mag|FK mag|M0 phase (rad)|M1 phase (rad)|FK phase (rad)|
|---|---:|---:|---:|---:|---:|---:|---:|
|weak_low_K|0.3|0.00250311|0.00249627|0.000180263|0.0158945|0.00311257|0.000178551|
|weak_low_K|0.5|0.00351606|0.00389406|5.24037e-07|0.0134017|0.00303976|1.10039e-05|
|weak_low_K|1|0.00888724|0.0095481|3.1008e-08|0.00494182|0.00376657|8.20422e-06|
|weak_low_K|2|0.0163787|0.0167316|6.63201e-08|0.00439746|0.00438321|5.37367e-06|
|weak_low_K|5|0.0106085|0.0106437|4.64779e-08|0.00397518|0.00397502|7.2434e-06|
|strong_height_low_K|0.3|0.0100156|0.00998813|0.00017645|0.0803364|0.0498092|0.000187289|
|strong_height_low_K|0.5|0.0140814|0.0155936|1.89938e-06|0.071541|0.0486811|4.39753e-05|
|strong_height_low_K|1|0.0355494|0.0381909|1.24484e-07|0.0618843|0.0602746|3.28036e-05|
|strong_height_low_K|2|0.0654652|0.0668742|2.65672e-07|0.0702028|0.0701375|2.14858e-05|
|strong_height_low_K|5|0.0424292|0.0425699|1.86444e-07|0.0636058|0.0636107|2.89622e-05|
|weak_high_K|0.3|0.0188514|0.018757|0.000177216|0.0130149|0.0106476|0.000205184|
|weak_high_K|0.5|0.0194504|0.0194777|2.24386e-06|0.0114238|0.0107346|9.40461e-05|
|weak_high_K|1|0.0229507|0.0231157|2.32754e-06|0.0103128|0.0102063|9.6745e-05|
|weak_high_K|2|0.0261509|0.0261014|2.72306e-06|0.0104818|0.0104771|9.41411e-05|
|weak_high_K|5|0.026609|0.0266186|2.62681e-06|0.0111636|0.0111636|9.07234e-05|

## Interpretation

- sigma=0.3 is the broadest tested source spectrum (theta99 about 21.28 deg); sigma=5 is the narrowest (about 1.23 deg).
- The Model-1 column isolates the existing kz-aware correction without changing the rough-surface physics. Improvement ratio >1 means the correction lowers normalized complex error; values near 1 mean the residual is not explained by the source-angle approximation.
- In the generated table, weak_low_K improves by about 4.03x at sigma=0.3 and 2.80x at sigma=0.5, but is essentially unchanged for sigma>=1. strong_height_low_K improves by about 1.59x and 1.43x at sigma=0.3/0.5, then also becomes approximately unchanged. weak_high_K improves by only about 1.06x at sigma=0.3 and is effectively unchanged for the narrower spectra.
- Full Kirchhoff remains the BIE control and should remain near its accepted 1e-4--1e-5 level. Any Model-0/Model-1 residual that persists while FK-BIE stays small is a PE approximation residual, not a BIE or Kirchhoff failure.
- For future vertical-link use, report theta95/theta99 together with sigma; sigma alone is not a universal accuracy guarantee. If Model-1 improvement is concentrated at broad theta99 and disappears for narrow spectra, finite-angle propagation is implicated. If a strong residual remains for narrow spectra or high-K geometry, nonlocal rough-surface coupling remains the dominant explanation.

## Visualizations

- [Angular spectrum](../results/validation/pe_bie_full_kirchhoff_gaussian_sigma_kzaware/figures/angular_spectrum_sigma.png) and [source-beam magnitude](../results/validation/pe_bie_full_kirchhoff_gaussian_sigma_kzaware/figures/source_beam_magnitude_sigma.png).
- Error curves: `model0_bie_ec_vs_sigma.png`, `model1_bie_ec_vs_sigma.png`, `fk_bie_ec_vs_sigma.png`; corresponding `_mag_` and `_phase_` figures are in the same directory.
- [Improvement versus theta99](../results/validation/pe_bie_full_kirchhoff_gaussian_sigma_kzaware/figures/improvement_vs_theta99.png).
- Representative receiver fields: `weak_low_K_representative_fields.png`, `strong_height_low_K_representative_fields.png`, and `weak_high_K_representative_fields.png`.

Artifacts: `results/validation/pe_bie_full_kirchhoff_gaussian_sigma_kzaware/gaussian_sigma_kzaware_comparison.mat`, `results/validation/pe_bie_full_kirchhoff_gaussian_sigma_kzaware/gaussian_sigma_kzaware_master.csv`, `results/validation/pe_bie_full_kirchhoff_gaussian_sigma_kzaware/gaussian_sigma_kzaware_detail.csv`; figures are under `results/validation/pe_bie_full_kirchhoff_gaussian_sigma_kzaware/figures/`.
