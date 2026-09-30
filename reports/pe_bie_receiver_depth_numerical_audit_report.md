# Receiver-depth numerical audit

Status: validation-only numerical审查; production PE, Full Kirchhoff, BIE and surface physics were not modified. All primary trends use the same sigma=1 M99 support and energy weights at every measured depth; no depth interpolation is used.

## Scope and conventions

The audit reuses `results/validation/pe_bie_full_kirchhoff_receiver_depth_sensitivity/receiver_depth_sensitivity.mat` and reconstructs the flat PE reflected field at the seven stored receiver depths. FK and BIE are transformed using the already accepted conjugation convention; no fitted amplitude or phase scalar is applied. “Raw” fields are also retained, but raw cross-solver comparison is labelled diagnostic because the solvers use different stored field conventions.

## Unified result table

The complete table is `results/validation/pe_bie_receiver_depth_numerical_audit/receiver_depth_numerical_audit.csv`. The columns include raw PE/FK/BIE complex errors, common-convention PE--FK/PE--BIE/FK--BIE errors, magnitude and phase components, and the fixed-mask sample count.

### weak_low_K

- z=1 m: common PE--FK Ec=0.00452328 (mag 0.00431239, phase 0.00136508), PE--BIE Ec=0.169491, FK--BIE Ec=0.170855.
- z=3 m: common PE--FK Ec=0.0102641 (mag 0.00954812, phase 0.00376659), PE--BIE Ec=0.0102641, FK--BIE Ec=8.20441e-06.
- z=10 m: common PE--FK Ec=0.0287294 (mag 0.0261965, phase 0.0117978), PE--BIE Ec=0.0287292, FK--BIE Ec=7.7971e-06.
- z=30 m: common PE--FK Ec=0.0676778 (mag 0.05981, phase 0.0317064), PE--BIE Ec=0.0676774, FK--BIE Ec=6.79542e-06.
### strong_height_low_K

- z=1 m: common PE--FK Ec=0.0278426 (mag 0.0172636, phase 0.0218552), PE--BIE Ec=0.881096, FK--BIE Ec=0.881852.
- z=3 m: common PE--FK Ec=0.0713035 (mag 0.038191, phase 0.0602747), PE--BIE Ec=0.0713033, FK--BIE Ec=3.28165e-05.
- z=10 m: common PE--FK Ec=0.214875 (mag 0.104589, phase 0.188752), PE--BIE Ec=0.214873, FK--BIE Ec=3.11584e-05.
- z=30 m: common PE--FK Ec=0.54586 (mag 0.236355, phase 0.506688), PE--BIE Ec=0.545856, FK--BIE Ec=2.6901e-05.
### weak_high_K

- z=1 m: common PE--FK Ec=0.00878048 (mag 0.00807051, phase 0.00345911), PE--BIE Ec=0.0643332, FK--BIE Ec=0.0670301.
- z=3 m: common PE--FK Ec=0.0252701 (mag 0.023118, phase 0.0102087), PE--BIE Ec=0.0252671, FK--BIE Ec=9.67738e-05.
- z=10 m: common PE--FK Ec=0.079131 (mag 0.0722178, phase 0.0324525), PE--BIE Ec=0.0791218, FK--BIE Ec=9.49493e-05.
- z=30 m: common PE--FK Ec=0.18723 (mag 0.168091, phase 0.0837963), PE--BIE Ec=0.187208, FK--BIE Ec=8.92715e-05.

## Flat control

The flat control compares PE flat reflected fields against the closed-form pressure-release image field on the measured grid. It is a propagation/convention control, not a new rough-surface case. A non-small flat error at a depth would invalidate attributing the corresponding rough residual to the surface model.

The PE envelope and analytic/FK carrier fields are not numerically identical before the accepted rough/flat convention is applied; raw absolute cross-solver complex error is therefore diagnostic only. The primary audit uses the same flat reference, fixed M99 support and weights at every depth, without a fitted depth-dependent scalar.

## z=1 m convergence status

The previous FK nominal/refined z=1 diagnostic is retained in the MAT file. It reported relative change 0.172024 for `weak_low_K`; the audit does not promote this to a converged physical result. Because the stored BIE receiver evaluator used `max(r,1)` for all sub-unit distances, the old z=1 BIE values are not accepted as a hard gate. z=1 is therefore **UNVERIFIED** until a corrected near-boundary BIE evaluation and a demonstrably converged FK quadrature are available.

## Interpretation

For z>=2 m, FK--BIE remains the numerical control while PE--FK and PE--BIE increase again at larger receiver depths in the stored cases. This is evidence of a propagation/field-comparison trend, but not by itself proof of nonlocal scattering: the flat control, common raw/normalized comparison and fixed mask must be checked together. The present audit does not relabel that increase as physical until the flat-control and convention diagnostics are read jointly.

The receiver-depth plots are discrete samples at 1,2,3,5,10,20,30 m; lines connecting markers are visual aids only.

## Answers to the requested checks

1. **30 m:** the increase is reproducible under the fixed 117-sample M99 mask and is not explained by FK--BIE control error (about 1e-5--1e-4 for z>=2 m). It is numerically credible as a normalized comparison trend, but not yet an unconditional physical receiver-depth law.
2. **Normalization/evaluation:** the trend is unchanged when the saved rough/flat ratios are compared directly; no depth-specific mask or interpolation is used.
3. **Error character:** the component balance is case dependent: weak_low_K is amplitude-dominant at 30 m (0.05981 versus phase 0.03171), while strong_height_low_K is phase-dominant (magnitude 0.2364 versus phase 0.5067). Thus the growth is not a pure amplitude effect; both components must be retained.
4. **1 m:** not resolved. The stored FK 4097-to-8193 change is 0.172 and the previous BIE evaluator clipped sub-unit distances; z=1 m is UNVERIFIED.
5. **Depth conclusion:** evidence supports a reproducible trend only for z>=2 m; it is insufficient to attribute the growth to nonlocal scattering without a common flat carrier control.

## Artifacts

Figures are under `results/validation/pe_bie_receiver_depth_numerical_audit/figures`; MAT: `results/validation/pe_bie_receiver_depth_numerical_audit/receiver_depth_numerical_audit.mat`; CSV: `results/validation/pe_bie_receiver_depth_numerical_audit/receiver_depth_numerical_audit.csv`.
