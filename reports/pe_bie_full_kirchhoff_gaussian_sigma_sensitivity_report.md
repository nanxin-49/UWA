# Gaussian source angular-spectrum sensitivity study

Status: validation-only; PE, Full Kirchhoff, BIE and surface physics unchanged.

Setup: 2D one-transverse, 4 kHz, c=1500 m/s, Tx=100 m, Rx=3 m, sigma=[0.3 0.5 1 2 5] m. Cases: weak low-K (A=0.05,K=0.10), strong-height low-K (A=0.20,K=0.10), weak high-K (A=0.02,K=0.47).

Mask: sigma-specific 99-percent source-energy footprint (strict-normal footprint when available; source-path footprint for sigma=5) plus pairwise -40 dB threshold, with weights renormalized on that footprint. Source uses the existing discrete Gaussian exp(-i omega t) evaluator; only sigma varies.

## Results

|case|sigma|theta90|theta95|theta99|PE-BIE Ec|PE-BIE mag|PE-BIE phase|FK-BIE Ec|FK-BIE mag|FK-BIE phase|PE-FK Ec|
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
|weak_low_K|0.3|13.4263|15.9693|21.2801|0.01609|0.00250311|0.015816|0.000253723|0.000180263|0.000177669|0.0160829|
|weak_low_K|0.5|7.96317|9.54683|12.5086|0.013855|0.00351606|0.012578|1.10164e-05|5.24037e-07|1.03276e-05|0.0138462|
|weak_low_K|1|4.028|4.70073|6.27314|0.0101692|0.00888724|0.00372039|8.20441e-06|3.1008e-08|6.17645e-06|0.0101669|
|weak_low_K|2|2.01275|2.34839|3.13187|0.0169587|0.0163787|0.00262877|5.37393e-06|6.63201e-08|3.21234e-06|0.0169587|
|weak_low_K|5|0.782601|0.894409|1.22986|0.0113287|0.0106085|0.00397518|7.24358e-06|4.64779e-08|7.2434e-06|0.0113288|
|strong_height_low_K|0.3|13.4263|15.9693|21.2801|0.0809081|0.0100156|0.0799395|0.000257278|0.00017645|0.000186364|0.0808828|
|strong_height_low_K|0.5|7.96317|9.54683|12.5086|0.0728833|0.0140814|0.071541|4.40203e-05|1.89938e-06|4.39753e-05|0.0728567|
|strong_height_low_K|1|4.028|4.70073|6.27314|0.0713342|0.0355494|0.0618843|3.28167e-05|1.24484e-07|3.28036e-05|0.0713291|
|strong_height_low_K|2|2.01275|2.34839|3.13187|0.0959002|0.0654652|0.0419667|2.14855e-05|2.65672e-07|1.2844e-05|0.0959004|
|strong_height_low_K|5|0.782601|0.894409|1.22986|0.0763957|0.0424292|0.0636058|2.89695e-05|1.86444e-07|2.89622e-05|0.076396|
|weak_high_K|0.3|13.4263|15.9693|21.2801|0.0229061|0.0188514|0.0129506|0.000271135|0.000177216|0.000204171|0.0228908|
|weak_high_K|0.5|7.96317|9.54683|12.5086|0.0225558|0.0194504|0.0114238|9.40723e-05|2.24386e-06|9.40461e-05|0.0225479|
|weak_high_K|1|4.028|4.70073|6.27314|0.0251597|0.0229507|0.0103128|9.67738e-05|2.32754e-06|9.6745e-05|0.0251589|
|weak_high_K|2|2.01275|2.34839|3.13187|0.0281717|0.0261509|0.00626596|9.41935e-05|2.72306e-06|5.62768e-05|0.0281746|
|weak_high_K|5|0.782601|0.894409|1.22986|0.0288541|0.026609|0.0111636|9.07573e-05|2.62681e-06|9.07234e-05|0.0288575|

## Direct answers

- sigma=0.3 is the broad-angle endpoint: theta99=21.2801 deg. sigma=5 narrows this to theta99=1.22986 deg.
- weak low-K: PE-BIE Ec changes from 0.01609 (sigma=0.3) to 0.0113287 (sigma=5), while PE-BIE phase RMS changes from 0.015816 to 0.00397518 rad; FK-BIE remains at 2.53723e-4 and 7.24358e-6. This supports a finite-angle contribution but not complete error removal.
- strong-height low-K: PE-BIE Ec changes from 0.0809081 to 0.0763957 and phase RMS from 0.0799395 to 0.0636058 rad. A large residual remains after narrowing, so strong-height path/nonlocal phase effects are not explained by angle width alone.
- weak high-K: PE-BIE Ec changes from 0.0229061 to 0.0288541; phase RMS decreases from 0.0129506 to 0.0111636 rad while magnitude error increases from 0.0188514 to 0.026609. Narrowing therefore does not improve total high-K error.

## Interpretation

The result is **PARTIAL_FINITE_ANGLE_CONTRIBUTION**: sigma=0.3 is demonstrably broad, and narrowing the source reduces the low-K phase-screen component, but PE does not converge monotonically to Full Kirchhoff/BIE across the tested cases. Full Kirchhoff remains close to BIE (approximately 1e-5 to 2.7e-4 complex L2), so the remaining PE residual is attributable primarily to the local phase-screen/nonlocal rough-surface coupling rather than to Full Kirchhoff error. A narrow source should be used for deep-vertical links only when its measured theta95/theta99 lies within the validated case envelope; sigma alone is not a universal accuracy control.

Artifacts: results/validation/pe_bie_full_kirchhoff_gaussian_sigma_sensitivity/gaussian_sigma_sensitivity.mat, results/validation/pe_bie_full_kirchhoff_gaussian_sigma_sensitivity/gaussian_sigma_sensitivity.csv; figures are under results/validation/pe_bie_full_kirchhoff_gaussian_sigma_sensitivity/figures.
