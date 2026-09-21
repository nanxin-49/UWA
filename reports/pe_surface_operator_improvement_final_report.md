# PE rough-surface reflection-operator improvement Goal: final report

Final conclusion: **NONLOCAL_EFFECT_REQUIRED**

## Executed stages

- **G0 PASS** — added `A=0.02 m, K=0.47 rad/m`; BIE uncertainty is
  `6.01e-8`, every numerical/geometry gate passes, and Bellhop remains much
  closer to BIE (`E_G=1.37e-4`) than Model-0 PE (`2.29e-2`).
- **G1 PASS** — the frozen Gaussian spectrum has `kx_rms=2.357 rad/m`,
  `theta_rms=8.144 deg`, and `<kz>/k=0.989951`. The omitted finite-angle phase
  has the correct scale and error-vector direction, especially in the weak and
  Region-II low-K cases.
- **G2 PASS: NORMAL_APPROXIMATION_CONFIRMED** — validation-only Model-1 uses
  componentwise `exp(+i 2 kz eta)`. It leaves flat propagation unchanged and
  improves all four controlled cases without fitting.
- **G3 FAIL / STOP** — adding analytic local-slope specular coupling produces
  negligible additional improvement and fails the frozen high-K gate. All
  reflected spectral branches remain valid and finite.
- **G4 COMPLETE** — the remaining discrepancy requires investigation of a
  nonlocal surface operator. No further local empirical correction is allowed.
- **G5 NOT RUN** — fixed PM remains locked by the G3 stop condition.

## Code scope

Only validation code changed:

- `scripts/validation/support/run_pe_1d_surface_reflection_validation.m`
  retains default `model0_normal` and adds independent validation-only
  `model1_kz_aware` and `model2_angle_slope` switches.
- `scripts/validation/validate_pe_bellhop_helmholtz_bie_reference.m` adds the
  gated G0 high-K entry while reusing the accepted BIE/Bellhop infrastructure.
- `scripts/validation/validate_pe_surface_normal_approximation_audit.m`,
  `validate_pe_surface_kz_aware_operator.m`, and
  `validate_pe_surface_angle_slope_operator.m` implement G1--G3.

Production PE, PE marching, the production surface model, Bellhop,
`Reflect2D`, `InfluenceGeoHatCart`, and the BIE solver are unchanged.

## Authoritative outputs

- `reports/pe_surface_operator_bie_G0_high_K_diagnostic.md`
- `reports/pe_surface_normal_approximation_G1_audit.md`
- `reports/pe_surface_kz_aware_G2_validation.md`
- `reports/pe_surface_angle_slope_G3_validation.md`
- `reports/pe_surface_nonlocal_G4_review.md`
- artifacts under
  `results/validation/pe_surface_operator_bie_reference/G0_high_K/` through
  `G3_angle_slope/`.

The validated Model-1 is a useful controlled-case improvement, but this Goal
does not promote it into production because strong-height and high-K residuals
remain far above the BIE uncertainty and Bellhop reference error.

