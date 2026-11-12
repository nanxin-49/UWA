# History review receipt
Checked 2026-11-11 via local Git and connected GitHub plugin.
Local starting HEAD: a37fd759684d124ce476f36acda5294254f32e79.
GitHub source commit: f3314fded85993832771f48e398f2a79edee582a.
Alternate local/GitHub historical commit: 0ea9bb14c9a98b277c32e160112c05b91a15d6ef.
The two old scripts specify piston_radius_m=0.05 and radius_m=0.05, and evaluate circular 2*J1(u)/u with u=k*a*sin(theta).
Old report, two checkpoints and old piston_frequency_pattern.png were inspected; old beam-map plotting source reads and asserts radius 0.05 m.
The old source and frequency runner SHA values match their historical checkpoint entries; no claim of a fresh full historical artifact audit.
No historical default runner was executed; no retired cash/ artifact was accessed.
Decision: 20 cm total-width 2-D strip vs 10 cm total-width strip as a half-size mapping, with literal 5 cm total-width auxiliary control. This is not dimensional equivalence of circular and strip radiation.
GitHub fixed source: https://github.com/nanxin-49/UWA/blob/f3314fded85993832771f48e398f2a79edee582a/scripts/comparisons/run_piston_frequency_minimal_comparison.m
