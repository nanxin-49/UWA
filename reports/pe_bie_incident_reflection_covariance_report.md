# PE / Full Kirchhoff / BIE incident and reflected field covariance audit

Validation-only audit at 4 kHz, c=1500 m/s, Tx=100 m, sigma=1 m. Cases are weak_low_K (A=0.05, K=0.1), strong_height_low_K (A=0.2, K=0.1), and weak_high_K (A=0.02, K=0.47). Existing receiver-depth fields are reused at z=2 and 3 m. No production physics was changed.

## Incident field check

At the reference plane z=0, PE `surface_incident_field` is compared with the same angular-spectrum source evaluation used by Full Kirchhoff. At the true rough surface, PE is continued from z=0 with the same one-way kz spectrum and compared with the saved FK incident field. The BIE solver does not store an independent incident trace; its source term is the same angular-spectrum construction, so FK/BIE incident agreement is a formulation identity rather than a separately measured receiver field.

- `weak_low_K`: reference-plane Ec=2.62103e-13, true-surface Ec=2.6307e-13; reference-plane and true-surface plots are in `results/validation/pe_bie_incident_reflection_covariance/figures`.
- `strong_height_low_K`: reference-plane Ec=2.62103e-13, true-surface Ec=2.66403e-13; reference-plane and true-surface plots are in `results/validation/pe_bie_incident_reflection_covariance/figures`.
- `weak_high_K`: reference-plane Ec=2.62103e-13, true-surface Ec=2.64463e-13; reference-plane and true-surface plots are in `results/validation/pe_bie_incident_reflection_covariance/figures`.

## Reflected field checks

- `weak_low_K`, z=2 m: PE--FK Ec=0.00740366 (magnitude 0.00694316, phase 0.00257046), PE--BIE Ec=0.00740365, FK--BIE Ec=8.26518e-06.
- `weak_low_K`, z=3 m: PE--FK Ec=0.0102641 (magnitude 0.00954812, phase 0.00376659), PE--BIE Ec=0.0102641, FK--BIE Ec=8.20441e-06.
- `strong_height_low_K`, z=2 m: PE--FK Ec=0.0496138 (magnitude 0.0277787, phase 0.0411373), PE--BIE Ec=0.049614, FK--BIE Ec=3.30609e-05.
- `strong_height_low_K`, z=3 m: PE--FK Ec=0.0713035 (magnitude 0.038191, phase 0.0602747), PE--BIE Ec=0.0713033, FK--BIE Ec=3.28165e-05.
- `weak_high_K`, z=2 m: PE--FK Ec=0.0170753 (magnitude 0.015639, phase 0.00685597), PE--BIE Ec=0.0170731, FK--BIE Ec=9.69427e-05.
- `weak_high_K`, z=3 m: PE--FK Ec=0.0252701 (magnitude 0.023118, phase 0.0102087), PE--BIE Ec=0.0252671, FK--BIE Ec=9.67738e-05.

### kx-region diagnostics

The following unweighted spectral summaries use the common M99/Hann support. They are diagnostics because bins with very small spectral amplitude can have unstable phase. Each tuple is `(mean amplitude difference, RMS phase difference)` for `|kx|/k` bands `[0,.25)`, `[.25,.5)`, `[.5,.8)`, `[.8,1)`.

- `weak_low_K`, z=2 m PE--FK: [0.00,0.25): 3.78e-07 / 0.137; [0.25,0.50): 4.02e-10 / 0.217; [0.50,0.80): 6.02e-11 / 0.195; [0.80,1.00): 1.01e-11 / 0.254; 
- `weak_low_K`, z=3 m PE--FK: [0.00,0.25): 3.67e-07 / 0.167; [0.25,0.50): 4e-10 / 0.22; [0.50,0.80): 5.99e-11 / 0.131; [0.80,1.00): 1e-11 / 0.198; 
- `strong_height_low_K`, z=2 m PE--FK: [0.00,0.25): 1.54e-05 / 0.386; [0.25,0.50): 7.2e-09 / 0.354; [0.50,0.80): 1.05e-09 / 0.398; [0.80,1.00): 1.83e-10 / 0.315; 
- `strong_height_low_K`, z=3 m PE--FK: [0.00,0.25): 1.5e-05 / 0.362; [0.25,0.50): 7.19e-09 / 0.416; [0.50,0.80): 1.05e-09 / 0.381; [0.80,1.00): 1.82e-10 / 0.464; 
- `weak_high_K`, z=2 m PE--FK: [0.00,0.25): 2.86e-06 / 0.366; [0.25,0.50): 7.27e-10 / 0.185; [0.50,0.80): 1.02e-10 / 0.0935; [0.80,1.00): 1.87e-11 / 0.267; 
- `weak_high_K`, z=3 m PE--FK: [0.00,0.25): 2.84e-06 / 0.453; [0.25,0.50): 7.8e-10 / 0.247; [0.50,0.80): 1.1e-10 / 0.122; [0.80,1.00): 1.95e-11 / 0.312; 

## Interpretation

The incident-field test checks the boundary input before reflection. If its fixed-mask complex error is small while the reflected PE--FK/BIE error is larger, the discrepancy begins in the reflection operator or in its spatial/spectral reconstruction. The saved spectra report normalized |P(k_x)| and PE--FK/FK--BIE spectral amplitude differences for the same M99 support; inspect the high-|k_x| tail first. Strong height is expected to emphasize path-phase and nonlocal coherence, while high K can additionally redistribute energy across kx.

For a reduced or fast Kirchhoff operator, retain the true surface path phase and the nonlocal coherent kernel. The incident angular-spectrum representation and its obliquity/normal derivative should remain explicit; a purely local height phase is insufficient when the reflected spectrum changes.

## Artifacts

- MAT: `results/validation/pe_bie_incident_reflection_covariance/incident_reflection_covariance.mat`
- CSV: `results/validation/pe_bie_incident_reflection_covariance/incident_reflection_covariance.csv`
- Figures: `results/validation/pe_bie_incident_reflection_covariance/figures`

The receiver samples are the saved discrete x grid; no depth interpolation was introduced. The audit does not claim that z=1 m is validated.
