# Analytic exterior identity: saved-array test

For the same frozen angular spectrum, each incident plane mode is a(x)=a0 exp(i*kx*x) at z=0, with d_n u_inc=-i*kz*a. The outgoing Weyl transform is integral G(x-s,z) exp(i*kx*s) ds = i exp(i*kx*x+i*kz*z)/(2*kz). Its source-normal derivative is one half of the same plane mode. Thus each of the two Kirchhoff terms contributes minus one half of the mirrored mode: the infinite flat double-term integral is exactly the negative mirror field, with the outgoing/limiting-absorption interpretation. The frozen modes have nonzero kz.

The C2 physical geometry is exactly flat for abs(x)>=50 m. There the true normal, ds, incident field, incident normal derivative and both Green terms are identical to flat. Therefore the infinite full integral can be evaluated by

u_complete = u_flat_analytic + FK_rough([-L,L]) - FK_flat([-L,L]), L>50 m.

This is analytic evaluation of the unchanged flat exterior; it is not a fitted scalar, rough/flat normalization, per-point calibration, FastFK approximation, or setting local_radius_m finite. The subtraction has compact support because the physical integrands outside the surface support are IDENTICAL, not because a physical contribution is discarded. Finite-grid residuals and exact support checks must still be audited.

The saved-array test records all four previously computed quadrature/window variants and 411-point line candidates in separate files. No acoustic computation ran in this test. Original finite-window main data/metrics and all prior audits remain unchanged. Candidate fields are not promoted until an independent direct integrand subtraction and CSV audit support the identity. This finding alone does not qualify the historical BIE reference or close the Goal.
