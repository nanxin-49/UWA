# PE--Bellhop--Helmholtz BIE third-reference final report

日期：2026-09-15  
最终状态：**BELLHOP_CLOSER_TO_HELMHOLTZ_REFERENCE**

## Scope and result

This validation adds an independent two-dimensional pressure-release Helmholtz
reference.  It does not call the PE phase screen, Bellhop `Reflect2D`,
`InfluenceGeoHatCart`, or any ray/local-specular formula.  R1--R5 were executed
in order; optional R6 was not required because R4 and R5 were already strongly
discriminating and refinement-stable.

| stage | surface | PE--BIE E_G | Bellhop--BIE E_G | PE--Bellhop E_G | U_BIE | result |
|---|---|---:|---:|---:|---:|---|
| R3 | A=0.01 m, K=0.10 rad/m | 0.00315757 | 8.09198e-6 | 0.00315855 | 6.00191e-8 | weak closure PASS |
| R4 | A=0.05 m, K=0.10 rad/m | 0.01608996 | 3.88955e-5 | 0.01609469 | 6.02438e-8 | Bellhop closer |
| R5 | A=0.20 m, K=0.10 rad/m | 0.08090811 | 9.97602e-5 | 0.08092279 | 6.17857e-8 | Bellhop closer |

No complex scalar, amplitude fit, phase subtraction, or per-case conjugation
choice was used.  Bellhop and BIE are both native `exp(-i*omega*t)` Helmholtz
ratios and are mapped by the same fixed conjugation into the established PE
comparison convention.  Native Bellhop--BIE metrics are also saved and agree
with the comparison-convention metrics.

## 1. BIE formulation

The reference uses a sound-soft rough graph in `D={z>eta(x)}`, the outgoing
fundamental solution

```text
Phi(X,Y) = i/4 H_0^(1)(k|X-Y|),   exp(-i*omega*t),
```

a Dirichlet half-plane Green function with fixed auxiliary line `h=-2 m`, and
the inward-to-water normal.  The accepted combined representation is

```text
u_ref = (D_h - i*k*S_h) psi,
(I + 2*K_h*W_L - 2*i*k*S_h*W_L) psi = -2*u_inc.
```

The trace is `gamma_D D_h=+1/2 I+K_h`.  A slow-rise numerical window is kept
outside the fixed physical surface support.  Panel Gauss--Legendre Nyström
collocation, high-order self/adjacent-panel quadrature, an independent off-grid
boundary residual, and direct receiver evaluation are used.  R1 smoke testing
identified and corrected the provisional unstable `+i*k*S_h` sign before any
result was accepted; this was a convention derivation correction, not a fitted
field adjustment.  Full formulation details and literature are in the R0
audit.

## 2. Flat validation accuracy

The accepted R1 12-points-per-wavelength result gives:

- off-grid boundary residual: `9.27624e-10`;
- image-solution complex L2 on M99: `1.23552e-7`;
- energy-weighted phase RMS on M99: `9.26794e-8 rad`;
- TL RMS on M95: `4.3736e-7 dB`;
- spatial last-level change: `7.68644e-9`;
- window last-level change: `3.12997e-8`.

All R1 hard gates pass.

## 3. BIE numerical uncertainty

R2 used one fixed C2-tapered benchmark
`eta=A sin(Kx) chi(x)`, with `A=0.01 m`, `K=0.10 rad/m`, `chi=1` for
`|x|<=42 m`, and physical support ending at `|x|=50 m`.  The physical taper was
unchanged in every run; only BIE spatial density, close-panel quadrature, and
the outer numerical window were varied.

The resulting receiver-line numerical uncertainty is:

- `U_BIE = 6.00191e-8` complex L2 on M99;
- `4.21974e-8 rad` phase RMS;
- `3.59586e-7 dB` TL RMS.

R4 and R5 repeat spatial/window refinement.  Their complex uncertainties are
`6.02438e-8` and `6.17857e-8`, respectively.

## 4. Weak three-way closure

R3 passes every frozen weak-limit gate.  PE--BIE and PE--Bellhop are both about
`3.16e-3`, while Bellhop--BIE is `8.09e-6`.  Bellhop wall intersection,
pressure-release phase, p/q rotation, positive transformed range, and receiver
selection guards all pass.  Thus the independent reference closes with both
approximations inside the existing weak-limit allowance, with Bellhop already
substantially closer to BIE.

## 5. Region-II adjudication

At `A=0.05 m`, PE--BIE is `0.01608996` and Bellhop--BIE is `3.88955e-5`.
The separation is `0.01605106`, compared with the required
`5 U_BIE = 3.01219e-7`.  Every spatial and window level selects Bellhop as the
closer model.  The difference is therefore reference-resolved rather than a
BIE numerical uncertainty.

## 6. Stronger-height trend

At `A=0.20 m`, PE--BIE grows to `0.08090811`, whereas Bellhop--BIE remains
`9.97602e-5`.  The same ranking holds at every refinement level, and the
observed separation `0.08080835` is far above `5 U_BIE = 3.08929e-7`.

One explicit limitation is retained: the strong-height independent boundary
residual reaches a `1.1e-8`--`2.9e-8` plateau rather than the R1 analytic
`1e-8` gate.  A separate order-10 panel calculation agrees with the main
receiver field to `3.29045e-8`, below the declared `U_BIE=6.17857e-8`.
Therefore this does not alter the resolved R5 ranking, but it must remain in
the uncertainty record.  The stricter R1 analytic gate was not changed.

## 7. Is fixed PM now necessary?

No, not for the first-round adjudication objective.  The smooth tapered
R3--R5 sequence already validates the reference, establishes weak closure,
and gives a stable Region-II/strong-height answer.  Fixed PM would answer a
new question about broadband random-surface geometry rather than resolve the
current PE--Bellhop ambiguity.  It should therefore be a separate follow-up,
starting with one fixed, band-limited realization and no Monte Carlo.

## 8. PM feasibility and WGF/FMM

The present dense direct solver is adequate for the tapered sinusoidal cases
but is not the recommended PM production route.  Systems of roughly
`3,800--4,500` unknowns already require minutes per dense solve.  A wider PM
support and higher retained surface wavenumber would increase both unknown
count and close-evaluation difficulty.

Before PM, retain the same physical band limit and realization in PE,
Bellhop, and BIE, and separately converge physical support/taper and BIE
window.  An analytic-flat-tail corrected WGF is the preferred formulation
upgrade if the slow-rise finite section no longer converges.  An iterative
solver with FMM/H-matrix acceleration is likely required for practical PM
refinement; it is not needed to substantiate the present R1--R5 conclusion.

## Authoritative artifacts

- `reports/pe_bellhop_helmholtz_bie_R0_formulation_audit.md`
- `reports/pe_bellhop_helmholtz_bie_R1_flat_validation_report.md`
- `reports/pe_bellhop_helmholtz_bie_R2_weak_rough_convergence_report.md`
- `reports/pe_bellhop_helmholtz_bie_R3_weak_three_way_closure_report.md`
- `reports/pe_bellhop_helmholtz_bie_R4_region_II_adjudication_report.md`
- `reports/pe_bellhop_helmholtz_bie_R5_stronger_height_report.md`
- `results/validation/pe_bellhop_helmholtz_bie_reference/R1_flat/`
- `results/validation/pe_bellhop_helmholtz_bie_reference/R2_weak_rough/`
- `results/validation/pe_bellhop_helmholtz_bie_reference/R3_weak_three_way/`
- `results/validation/pe_bellhop_helmholtz_bie_reference/R4_region_II/`
- `results/validation/pe_bellhop_helmholtz_bie_reference/R5_stronger_height/`

Optional R6 and fixed PM were not run because they are not required by the
completed first-round Goal.
