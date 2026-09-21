# PE rough-surface operator G4 nonlocality review

Status: **NONLOCAL_EFFECT_REQUIRED**

## Scope and frozen evidence

This review uses only the accepted 4 kHz, 1-transverse, sigma=0.3 m,
10,001-beam X-source Bellhop, and independent Helmholtz BIE artifacts. It does
not modify production PE, PE marching, Bellhop, or the BIE reference. No scalar,
phase, or amplitude parameter is fitted.

G0 added the missing small-height/high-wavenumber case
`A=0.02 m, K=0.47 rad/m`. All BIE refinement and Bellhop geometry gates pass.
Its errors are `E_G=0.0229061` for Model-0 PE and `1.36812e-4` for Bellhop.

## What the local hierarchy established

| case | Model-0 E_G | kz-aware Model-1 E_G | angle+slope Model-2 E_G | Bellhop E_G |
|---|---:|---:|---:|---:|
| A=0.01, K=0.10 | 0.00315757 | 0.000514551 | 0.000514113 | 8.09198e-6 |
| A=0.05, K=0.10 | 0.0160900 | 0.00398989 | 0.00395441 | 3.88955e-5 |
| A=0.20, K=0.10 | 0.0809081 | 0.0507808 | 0.0500519 | 9.97601e-5 |
| A=0.02, K=0.47 | 0.0229061 | 0.0215668 | 0.0215648 | 1.36812e-4 |

Model-1 replaces the strict-normal phase `2 k eta` by the unfitted
componentwise phase `2 kz(kx) eta`. It preserves flat response to
`1.83e-13`, improves every controlled case, and reduces phase RMS by factors
25.0, 5.11, 1.61, and 1.22. Therefore the finite-angle normal approximation
is a real and important Model-0 error source.

Model-2 additionally uses the exact local planar specular relation

```text
kzr = ((1-s^2) kzi + 2 s kxi)/(1+s^2),  s=deta/dx,
phase = (kzi+kzr) eta.
```

It has no non-returning spectral components and all minimum reflected `kz`
values remain positive. Nevertheless it improves Model-1 by only 0.01% to
1.5%; at high K the phase improvement is 0.03%. The G3 failure is therefore
not a grazing, branch-selection, NaN/Inf, or PE marching failure.

## Why another local phase correction is not justified

The exact Dirichlet rough-surface problem couples boundary position, normal
derivative, and the field over the entire illuminated surface. In spectral
language, a nonflat boundary maps one incident transverse wavenumber into a
distribution of outgoing transverse wavenumbers. A local tangent-plane phase
captures the leading eikonal shift but omits the noncommuting propagation and
boundary operators, boundary-density/Jacobian effects, curvature diffraction,
and repeated lateral interaction represented by the Helmholtz boundary
integral equation.

Model-2 already supplies the local incident-angle/local-normal geometry asked
for by G3. Its negligible incremental improvement, while Bellhop remains two
orders of magnitude closer to BIE, leaves no evidence-based coefficient or
additional local phase term to add. Doing so would be case-specific fitting.

## Smallest defensible next model family

The next investigation should be validation-only and nonlocal. Two defensible
routes are:

1. a surface Dirichlet-to-Neumann or coordinate-flattened one-way operator,
   implemented as a pseudodifferential/spectral coupling operator and checked
   term-by-term against BIE; or
2. a controlled Kirchhoff surface integral that maps the saved incident field
   to the receiver field, with its amplitude, normal derivative, and Green
   kernel fixed analytically rather than calibrated.

The first proof should reuse the four sinusoidal cases and the existing BIE
artifacts. It must preserve flat closure, expose its spectral coupling matrix,
and pass spatial/window refinement before any fixed-PM run. A full random PM
test or Monte Carlo is not justified yet.

## Gate decision

- G2 finite-angle mechanism: confirmed.
- G3 local angle--slope mechanism as a sufficient correction: rejected.
- G4 decision: **NONLOCAL_EFFECT_REQUIRED**.
- G5 fixed PM: **LOCKED / NOT RUN**.

