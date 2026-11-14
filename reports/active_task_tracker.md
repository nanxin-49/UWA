# Active Research Task Tracker (2-D rough-surface mainline)

Updated: 2026-11-14 (configuration/navigation alignment only; no MATLAB, SSH, propagation or communication run).

This is the short **execution-navigation ledger** for the current 2-D research mainline, not a second authoritative scientific plan or an automatic approval to run. Detailed physics, cost gates, acceptance evidence, and historical failures remain in the applicable reports and `PROJECT_CONTEXT.md`. The current user's task scope overrides this ledger.

## Frozen choices for upcoming tasks

- **Forward main source**: ideal 2-D rigid-baffle uniform-velocity strip, **total width W=0.10 m**, half-width a=0.05 m; peak normal velocity v0=0.001 m/s; source (x,z)=(0,100 m), upward. The historical 3-D circle radius 5 cm is only a nominal size correspondence, not evidence of a measured real transducer aperture.
- **Environment**: c=1500 m/s, rho=1025 kg/m^3, positive z downward; receive depth z=3 m; phase contract v2 and `exp(-i omega t)` native acoustic phasors. Do not normalize away absolute source differences between models.
- **Comparisons**: use BIE as a qualified *numerical reference with uncertainty*, not an exact solution; compare the same reflected complex field on identical true receiver coordinates. Keep PE Model-0, FullFK, and any FastFK algorithmic acceleration distinct.
- **PM input provenance**: begin with the **existing frozen** 4-kHz PM seed 260001, then 260002/260003 only after the pilot; preserve each seed's original Fourier coefficients, 160-m period, 12 nonzero modes, Kmax=0.471238898 rad/m, and the documented C2 embedding (|x| <=42 m plateau, vanishing at |x| >=50 m). Reuse the exact saved surfaces, do not re-sample.
- **Legacy controls**: the 20-cm piston P0/P1 and earlier Gaussian results remain immutable historical controls. Previous P0 covered **both 20 and 10 cm** at 4/6/8 kHz; previous P1 included **both widths at 4 kHz** (free, flat, weak sinusoid). Neither certifies 10-cm PM or broadband production use; the accepted P1 PE used a validation-only continuous-spectrum implementation.

## Current stages

| ID | Task and minimum exit evidence | State | Dependency / next permitted action |
|---|---|---|---|
| B0 | Reuse prior P0 source qualification and 4-kHz P1 flat/weak-sinusoid evidence for W=10 cm, with correct source/phase metadata | COMPLETE (limited) | No blanket rerun; provenance: `results/comparisons/piston_2d_single_frequency/20261112_remote_p0p1_01/report.md` |
| R1 | 4 kHz, W=10 cm, fixed PM seed 260001: PE Model-0 / BIE / full double-term FullFK, paired input and numerical-convergence checks, reflection complex-L2 and spatial residual plots | NEXT - NOT AUTHORIZED | Obtain explicit compute scope, budget and reference gates; reuse fixed PM surface. [Navigation alignment only, 2026-11-14](agents_configuration_alignment_2026-11-14.md); no R1 computation |
| R2 | Extend the same 4-kHz comparison to seeds 260002 and 260003 only if R1 has usable independent reference and cost evidence | PLANNED | Requires R1 result review and explicit run scope |
| R3 | Reflection-model decision from R1/R2: retain PE as baseline or candidate; assess FullFK physics separately from FastFK accuracy/cost. Test any *new physical simplification* independently | WAITING | Do not automatically discard PE or replace full Kirchhoff with a reduced physical formula |
| W1 | W=10 cm, coherent 4-8 kHz direct/reflected/total H(f), frequency-grid, phase/group-delay/holdout controls; conditionally validate FastFK against FullFK and physical anchors | WAITING | Depends on reflection-model decision and bounded physical-reference evidence |
| C1 | Verified CIR/linear-convolution controls followed by fixed-seed QPSK, matched-CSI and mismatched-CSI BER/EVM comparisons | WAITING | Requires coherent W1; common source excitation, power/noise and receiver settings |
| BH1 | Bellhop parallel line: distinguish existing Tx99 m / Rx3 m geometry, qualify X line-source mapping and flat/direct/first-reflection paths before rough-PM comparisons | PARALLEL - NOT AUTHORIZED | Independent of R1-R3; no direct absolute-Pa model ranking before source/phase qualification |

## Update rule for Codex (at every substantive task boundary)

1. Read this table and the exact user-authorized scope before any work. Move only tasks genuinely started to `IN_PROGRESS`; one task can be active in the PM line while Bellhop proceeds independently.
2. After a meaningful task, update the affected row to `COMPLETE (scope-limited)`, `INCONCLUSIVE`, `BLOCKED`, `STOPPED`, or its next honest state, with direct report/checkpoint links and key bounded findings. Do not edit old checkpoints to match new conclusions.
3. Keep a single short last-action record below. Update project-wide `reports/project_progress.json` and `PROJECT_CONTEXT.md` only when that broader status actually changes, under `AGENTS.md`.
4. The next stage remains **not authorized** until the user expressly requests the new scope; a table row is never permission for heavy computation, remote SSH, multi-process parallelism, or a new method.

## Last action

- **2026-11-14 — Configuration alignment complete (scope-limited); Git publication authorized:** Synchronized current context, README navigation and progress source/views with AGENTS: forward W=10 cm; W=20 cm and Gaussian historical controls retained. [Report and static-check receipt](agents_configuration_alignment_2026-11-14.md). The user subsequently requested a scoped commit/push; unrelated pending work stays outside this commit. R1 remains NOT AUTHORIZED (4 kHz, fixed PM seed260001); explicit compute budget, resource and numerical-reference gates remain unresolved. No MATLAB, laboratory SSH or experiment; no research stage advanced.
