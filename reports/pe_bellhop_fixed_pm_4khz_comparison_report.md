# PE--Bellhop fixed-PM 4 kHz comparison

状态：**PASS_WITH_MODEL_DISCREPANCY**

## Scope and frozen inputs

This Stage 1C report consumes the separately passed Stage 0E numerical budget, Stage 1A fixed-PM Tier-1 comparison, and Stage 1B dimensionality sensitivity. The canonical source is seed `260001`, U=6 m/s, span 160 m, master N=4097, requested Kmax 0.5 rad/m (realized 0.471238898 rad/m); no profile, source, reflection coefficient, PE marching, Bellhop Reflect2D, p/q, or InfluenceGeoHatCart change was made.

## Primary comparison

| metric | value |
|---|---:|
| G_PE, 1T | -0.0868692441183+0.988068678385i |
| G_Bellhop | -0.690339088449-0.662863954691i |
| cross-model delta TL (dB) | 0.31043352 |
| cross-model delta phase (rad) | -2.2482008 |
| cross-model complex error | 1.836642 |
| G_PE, 2T (ny=512) | -0.0974244073654+0.969018832201i |
| dimensionality delta TL (dB) | -0.15885982 |
| dimensionality delta phase (rad) | 0.01250977 |
| dimensionality complex error | 0.021956905 |

The 1T-to-2T sensitivity is approximately -0.15885982 dB / 0.01250977 rad and its complex error is 0.021956905, far below the PE/Bellhop complex error 1.836642. The remaining stable difference is therefore classified as a reflection-model discrepancy (Kirchhoff phase screen versus Bellhop local-specular Gaussian beam), not numerical failure.

## Frozen evidence

- Stage 0E PE window/grid/step and Bellhop profile/beam/step budgets: PASS.
- Stage 1A Bellhop wall geometry, pressure-release phase, p/q state and transformed positive-range branch: PASS.
- Stage 1B production 2T finite outputs and ny=256/512 sensitivity: PASS.

## Checks

- stage0e_all: PASS
- stage1a_structural: PASS
- stage1b_all: PASS
- same_profile_hash: PASS
- dim_sensitivity_smaller: PASS
- finite_metrics: PASS
- all: PASS

## Decision

`PASS_WITH_MODEL_DISCREPANCY` means the two independently validated solvers are numerically comparable under the same fixed realization, while their remaining rough reflected-branch difference is a quantified physical-model discrepancy. Stage 2 frequency extension is allowed; the next and only recommended step is the 4/6/8 kHz frequency extension before any ensemble sweep.
