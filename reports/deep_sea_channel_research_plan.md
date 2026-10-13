# Research Plan for a Deep-to-Near-Surface Acoustic Communication Channel

Date: 2026-10-09 (Asia/Shanghai).

Status: recommended research plan with limited WP0/WP1 progress below; not an approved engineering operating envelope. Unperformed experiments and acceptance levels remain proposals. Production physics is unchanged.

### Current execution update — 2026-10-12

2026-10-13 matrix-free t=1 mapped Helmholtz minimum pilot completed:
[report](deepsea3d_mapped_direct_report.md), [protocol](deepsea3d_mapped_direct_protocol.md).
Exact source-on-graph Dirichlet, unchanged dealiased operators and outgoing DtN, L0-preconditioned GMRES;
no Taylor expansion. Plane/dual/quarter controls pass, dual increment error1.0935e-5 to frozen reference.
Original PM restart8/12 each capped96 inner steps remain linearly unresolved; best true PDE/DtN
1.329e-3/1.273e-3 and true boundary1.04033e-4 fail original gates. Pressure stability is insufficient.
This improves Taylor32 residuals substantially but does not certify physical error. Stop iteration/grid expansion;
next bounded work should address preconditioning before formal convergence checks, with matched3D BIE fallback.
PM INCONCLUSIVE, G1 OPEN/WP2=false, no support24/FullFK–PM. Kernel wall464.135s, frozen312 files unchanged.

2026-10-13 bounded TFE dealiasing diagnosis completed:
[report](deepsea3d_tfe_dealiasing_report.md), [protocol](deepsea3d_tfe_dealiasing_protocol.md).
Cubic-exact2n padding and matched symmetric Fourier space, independent3n/4n true geometry;
dual/quarter local controls pass. Original PM N32 last coefficient norm improves2.363x,
but true PDE/DtN worsen .410575→.475266/.287852→.326480.99.60% of residual energy is retained.
Dealiasing does not rescue the current Taylor sum; do not stack further orders. Next preferred
bounded reference is a direct matrix-free t=1 mapped Helmholtz solve, with matched3D BIE fallback.
PM INCONCLUSIVE, G1 OPEN, WP2=false; no support24/FullFK–PM. Frozen284 files unchanged.
Kernel wall486.343s includes initial invalid comparator and saved-field independent probe.

2026-10-13 independent 3D TFE minimum prototype completed:
[report](deepsea3d_tfe_prototype_report.md), [protocol](deepsea3d_tfe_prototype_protocol.md).
Mapped Helmholtz with exact outgoing Fourier DtN, independently checked against manufactured modes
and the certified dual case: all dual numerical gates pass, fine increment differs6.308e-7 from frozen FE.
The same frozen PM seed261013/support12 fails at N40; h.125 still has physical PDE/DtN residuals .120/.173.
Quarter-amplitude local diagnostics pass but are not the target or a formally converged reference.
No support24 or FullFK–PM comparison; reference INCONCLUSIVE, full G1 OPEN, WP2=false.
Next bounded reference work may address dealiasing/TFE truncation or a matched3D BIE prototype.
Sixteen rough solves plus one plane control, saved kernel wall278.451s; frozen263 files unchanged.

2026-10-13 fixed-realization 3D PM feasibility completed:
[report](deepsea3d_pm_feasibility_report.md), [protocol](deepsea3d_pm_feasibility_protocol.md).
One isotropic finite-band Hm0=.5m parent (128m, kmin2pi/128, kmax1rad/m), seed261013,
6kHz, Tx1500/Rx5, q=.10; fixed alpha determines finite-band Tp3.638671s. No realized-RMS rescaling.
Nine bounded ordinary-FE probes fail order/source/true-boundary gates at N20; receiver mesh/domain
stability does not make the reference credible. Exact finite-mode exponential trace confirms failure.
Changing C2 support changes the physical sea; cutoff/support pressure differences remain unqualified.
No FullFK comparison or new conditional G1. Single-case INCONCLUSIVE, full G1 OPEN, WP2=false.
Next need matched TFE/3D BIE and justified noncompact support, then original physical/numerical gates;
no automatic WP2 release from the dual/frequency results below. Frozen solvers/production/data unchanged.

2026-10-13 fixed-y 4/8kHz conditional physics checks completed:
[report](deepsea3d_g1_frequency_report.md), [protocol](deepsea3d_g1_frequency_protocol.md).
Geometry/alpha2.5/dual/C2 frozen, 6kHz reused. Independent references require N16 at4kHz,
N20 at8kHz; mesh/domain/true-boundary/M32-to40 and cancellation controls pass.
Initial8kHz Gauss h.05 failed1e-3; one documented h.025 refinement passed without enlarging FE.
All three discrete points CONDITIONAL_G1_PASS. Reflected/increment complex errors decrease
at the primary settings, but order-sensitivity ranges overlap6-to8kHz: full complex-error
monotonicity is unresolved. Phase-error decrease and increment magnitude-error increase6-to8
are separable in observed variants. This does not certify a continuous band or generic trend.
Full G1 target OPEN, WP2=false; frozen solvers, production and historical data unchanged.

2026-10-13 fixed-sea q=.10 azimuth comparison completed:
[report](deepsea3d_g1_azimuth_report.md), [protocol](deepsea3d_g1_azimuth_protocol.md).
Tx changed from(150.5,0,1500) to(0,150.5,1500), Rx unchanged; no sea rotation/translation.
Flat backgrounds match; y-direction independent FE and FK numerical guards passed.
CONDITIONAL_G1_PASS for y: reflected/increment errors2.843e-4/1.375e-4.
Independent reference directional complex changes5.111%/2.467% are distinguished from
model-error reductions1.920%/1.762% and observed numerical/reference variations.
Empirical model-error ranges remain disjoint after saved lower-order pressure sensitivity.
No universal azimuth trend follows. Full target G1 OPEN, WP2=false; frozen data/solvers unchanged.

2026-10-13 q=.10 geometry-only continuation completed:
[report](deepsea3d_g1_q10_report.md), [protocol](deepsea3d_g1_q10_protocol.md).
Tx=(150.5,0,1500), Rx=(0,0,5), fixed alpha2.5/6kHz/dual/C2 and surface center(0,0).
WP0 plane controls were reused without reusing old moved-surface rough cases.
New N12/N16 precheck required N20; independent mesh/domain/true-boundary gates passed.
FullFK numerical guards passed, reflected/increment native errors2.899e-4/1.399e-4:
CONDITIONAL_G1_PASS only at this discrete geometry. Errors decreased2.03%/2.94% versus
frozen q0, with disjoint empirical numerical-variant ranges; total relative error slightly increased.
The complete G1 target domain stays OPEN, WP2=false; no amplitude scan/PM/FastFK/SSP/dynamics.

2026-10-13 alpha2.5 domain and physical check completed:
[report](deepsea3d_alpha25_domain_report.md), [protocol](deepsea3d_alpha25_domain_protocol.md).
Frozen L24-to32/h=.0625/N20 domain change3.242e-7 and true half-cell residual2.977e-7
pass the original1e-4 gates; L48 was not needed. Only after this, matched FullFK
mesh/Gauss/aperture controls passed; reflected/increment native complex errors
2.959e-4/1.442e-4 give CONDITIONAL_G1_PASS at alpha2.5 for this discrete case.
This supersedes the alpha2.5 physical INCONCLUSIVE status in the earlier stage records below.
The complete G1 target domain remains OPEN, WP2=false; no larger-amplitude/PM/FastFK/WP2 work.
Frozen solvers, production and historical acceptance data remain unchanged.

2026-10-13 alpha2.5 fine-reference check completed: [report](deepsea3d_alpha25_fine_report.md).
L24/h=.0625 N16/N20 passed order/grid/true-boundary controls; no observed FE instability.
Next a bounded domain check, then matched FullFK comparison if the reference remains stable.
No domain extension or FK solve was run in this stage; alpha2.5 physical G1 is INCONCLUSIVE,
the complete domain remains OPEN and WP2 is not authorized.

2026-10-13 alpha2.5 bounded precheck completed: [report](deepsea3d_alpha25_preflight_report.md).
N12 is insufficient; N16-to20 change and true half-cell boundary residual pass the existing
coarse precheck. Ordinary FE has no observed catastrophic instability on this configuration.
Next only a controlled finer-grid reference check, then domain/FK work if stable.
No formal alpha2.5 G1 pass or WP2 authorization follows; no full validation was started.

2026-10-13 amplitude-only continuation: [report](deepsea3d_g1_amplitude_report.md).
Frozen alpha1 reused; new .5/1.5/2 each passed independent order/mesh/domain/true-boundary
controls and conditional pointwise G1 at unchanged geometry/spectrum/phase/support.
No credible physical failure or FE instability was found. This certifies four sampled points,
not a continuous interval or the complete G1 domain; WP2 remains blocked by coverage gaps.
If pursuing larger amplitude, first probe alpha2.5 order stability at bounded cost.

Latest continuation: [finite-roughness single-case report](deepsea3d_g1_finite_report.md)
and [protocol](deepsea3d_g1_finite_protocol.md). An independent order16 Dirichlet FE reference
with grid/domain/order and true-boundary residual controls conditionally passes one
Tx1500/Rx5, q0, 6 kHz compact dual-direction surface, amplitudes (.02,.01) m.
FullFK reflected/increment errors are 1.36e-4/1.06e-4; finite-amplitude effects are substantial.
This closes the first-order-only gap for this case, not the full target domain G1.
WP2 remains unauthorized. Next use one controlled physical-domain extension; do not require
a dense BIE or repeat this completed single case as a blanket prerequisite.

[Chinese stage report](deepsea3d_wp0_wp1_report.md) and [frozen protocol](deepsea3d_wp0_wp1_protocol.md) supersede the initial execution order where specified here.

- The 2D PM physical-domain audit is complete (15 controlled inputs, four guarded cases). Do not repeat it as a blanket prerequisite. Its conditional boundary is not a universal cutoff or 3D acceptance criterion.
- Newer 2D single-Rx FastIn/composite prototypes are conditionally frozen only at original Rx(0,3), with cancellation-position failures. They add algorithm fidelity, not BIE physics; no acceleration/rank work is required for this 3D stage.
- WP0/G0 is PASS for Tx1500/Rx5, q=0/.10, unit Green point source, image fields, shifted planes, rotations and native/DSP algebra.
- Minimum WP1 now has the full two-term single-Rx 3D Kirchhoff integral with analytic gradients and exact infinite-flat background subtraction for compact C2 perturbations. Numerical gates pass at 6 kHz, then 4/8 kHz. Millimetric weak/dual-direction controls do not represent the future Hm0=.5/1 m PM pilot.
- One independent first-order Dirichlet shape-reference formulation is implemented: 12 amplitude/geometry comparisons, Gauss grid/support guards and Fourier cross-check. This is a weak-limit check, not finite-roughness BIE certification. Full G1 remains OPEN; these results do not authorize WP2/PM ensembles.
- The earlier next-step finite-amplitude reference is now completed conditionally for the case above. For harder cases test FE stability anew, or use TFE/matched half-space BIE. An open BEM patch or a 1% FK–first-order discrepancy still does not establish full G1.
- Section 9 numerical 1e-3 and physical screening 1% proposals remain distinct from the historical 2D physical 1e-3 gate. No tolerance was relaxed to pass this campaign.

## 1. Recommended direction and intended contribution

**Develop a receiver-focused, physically informed statistical model of the complex broadband response of one deep-to-near-surface link. Use a small, independently checked 3D physical reference to identify the random structure that must be preserved, then build the simplest generator that meets communication tolerances at a defensible total cost.**

My recommended initial physical candidate is an analytic homogeneous-medium incident field combined with a direct, single-receiver 3D rough-surface Kirchhoff integral. This is a candidate data-generation engine, not an already validated truth model. Its first task is to pass analytic and independent physical checks. It avoids requiring production PE integration before the homogeneous problem is understood and avoids computing a dense receiver field for a single-link output.

My preferred statistical hypothesis to test is a small number of random reflected-path parameters plus a jointly frequency-dependent residual. A conditional complex Gaussian model of the entire response is the mandatory simpler comparator. Neither wins in advance. If direct physical evaluation is sufficiently inexpensive, retaining it as the channel generator is also a successful result; a learned model is not required merely to claim acceleration.

The scientific contribution should be a **validated domain of applicability, a justified stochastic representation, and its communication-relevant accuracy/cost tradeoff**. A useful negative finding, such as a physically meaningful regime where few-path or Gaussian models fail, is part of that contribution. Publication novelty must be established against primary literature before claiming a new model; hybrid physical/statistical channel modeling already has precedents. [E4]

The main research questions are:

1. Which rough-surface mechanisms control the single-link broadband response in the proposed deep, near-vertical geometry?
2. Can their effects be represented by a small set of jointly distributed path variables, or is a richer joint-frequency process necessary?
3. Over which geometry, surface spectrum, and receiver-depth ranges does the selected physical/statistical approximation meet the chosen communication tolerance?
4. When do offline construction and condition updates pay for themselves relative to direct physical sampling?
5. Which state variables must be retained to extend the static model to consistent time evolution and Doppler?

## 2. Requirements, assumptions, and scope boundaries

| Category | Current specification | Status |
|---|---|---|
| Application | Deep transmitter to a cable-deployed receiver near a vessel or surface platform; one Tx/Rx link | User requirement |
| Depth | Tx deeper than 1000 m; Rx 3–10 m | User requirement; maximum Tx depth and bottom clearance unknown |
| Geometry | Near vertical, with horizontal displacement permitted | User requirement; operational displacement distribution unknown |
| Spatial physics | 3D propagation and a rough surface with two horizontal coordinates | User requirement; reductions need demonstrated error bounds |
| First stage | Static surface, uniform sound speed | User requirement |
| Frequency | Start with 4–8 kHz | User preference; changing band requires a scientific reason |
| Later capability | Evolving sea, time correlation, and Doppler | User requirement |
| Sea model | PM is a reasonable initial candidate | User preference; directionality, cutoff, and later alternatives open |
| Source | Transducer model and directivity unknown | Explicit uncertainty |
| Deployment | Cable deployment confirmed; hull clearance unspecified | Explicit uncertainty |
| Computing | Limited offline work allowed; count construction, validation, updates, and online sampling | User requirement; numerical budgets unknown |
| Acceptance | Communication performance plus necessary physical/statistical checks | User requirement; proposed tolerances appear in Section 9 |

The initial physical study is a uniform, lossless water half-space with a pressure-release surface, fixed endpoints, and an explicitly normalized point source. It isolates free-surface effects. It is not a calibrated shipboard link budget. Constant sound speed does not imply zero absorption in a real ocean; absorption is omitted only in this first controlled model.

Before making engineering range, transmit-power, or outage claims, add or bound frequency-dependent absorption, source/receiver responses, ambient and platform noise, bottom-related paths, and sound-speed structure as appropriate. Tx depth is not water depth. A source near the bottom may require a bottom-reflected contribution, so bottom clearance is a high-priority future input.

Hull scattering, wakes, bubbles, and motion are not declared negligible. Initially excluding them defines the idealized study domain; later deployment sensitivity determines whether the idealized result transfers. Unknown engineering information should not prevent analytic controls, but it limits the final claim.

## 3. What existing work changes in this plan

| Existing result | Appropriate use now | Consequence |
|---|---|---|
| PE/WAPE, receiver adjoint, cross-frequency statistics, and MPSK tools | Reuse conventions, interfaces, and independently verified numerical components | Do not assume that old local-screen statistics are statistics of a new rough-boundary solver |
| Completed FullFK–BIE deterministic/PM-domain comparisons | Reuse 2D regressions, conditional passes and guarded failures | Do not repeat the completed audit; no automatic physical transfer to 3D |
| Fixed-PM 4–8 kHz FastFK study, ranks 3/12 | Algorithm-fidelity and phase-convention regression | No reason to retain those ranks, sample counts, or surface cutoff in the new model |
| Single-receiver timing reversal | Match computation to the requested output | Suspend dense-field FastFK optimization as the main research activity |
| 2026-10-12 single-Rx conditional freeze | Reuse algorithm contracts and cancellation diagnostics | Not additional BIE physics or a 3D speed forecast |
| Finite-band CIR agreement | Reuse native-to-DSP conversion tests | Re-establish delay support and frequency-grid convergence |
| External AI suggestions | Retain useful questions and candidate methods | Recompute geometry, roughness moments, and budgets; do not inherit unsupported thresholds |

The measured fixed 4 kHz, 2D Tx100/Rx3 case gave 0.98571 ms for FullFK single-point evaluation versus 36.0479 ms for FastFK whole-line evaluation followed by extraction. Including preparation gave 111.55524 versus 214.69579 ms. FullFK was therefore approximately 36.57 and 1.925 times faster, respectively; pointwise complex difference was about 0.0177%. These are workload-specific measurements, not projected 3D timings. [L2]

In that implementation, most FullFK cost was incident-field reconstruction on the actual surface. An analytic 3D point-source field changes this preparation problem substantially, but it also changes the source model from the old 2D Gaussian fixture. Any resulting cost reduction must be reported as a new formulation/workload comparison, not an equivalent acceleration of the old fixture.

Continue to preserve the production WAPE model and all historical references. Research alternatives should be introduced through validation workflows and only integrated into public interfaces after their domain and semantics are established.

## 4. Proposed physical and statistical architecture

### 4.1 Physical sample generation

For a source at position x_t and a surface point s, start from the native outgoing 3D Green function

\[
G_k(s,x_t)=\frac{e^{ikR_t}}{4\pi R_t},\qquad R_t=\lVert s-x_t\rVert,\quad k=2\pi f/c,
\]

with time convention exp(-i omega t). Its analytic spatial derivative is

\[
\nabla_sG_k=G_k\left(ik-\frac{1}{R_t}\right)\frac{s-x_t}{R_t}.
\]

A declared source-strength spectrum multiplies this field. A directional source later requires the corresponding angular weighting and derivatives, not an undocumented beam factor. The same source normalization must be used for the direct and reflected responses.

Compute the receiver-directed surface integral using the actual two-coordinate surface, consistent normals, surface measure, and spherical source/receiver distances. Derive the chosen 3D boundary approximation before implementation. Do not copy the old Hankel kernel or assume that changing a 2D kernel to a 3D kernel certifies the entire formulation.

At fixed quadrature, direct single-point incident construction and integration can be organized with work proportional to the number of surface samples per frequency, rather than receiver-count times surface-count. Oscillatory quadrature, accuracy-driven surface size, and reference-solver complexity remain major unknowns. The proposed first implementation should prioritize transparent convergence and streamed storage over compression.

Use a hierarchy of references:

| Level | Reference | What it can establish |
|---|---|---|
| A | Exact flat pressure-release image solution; uniformly shifted flat surface; coordinate rotations | Normalization, geometry, signs, conventions, and simple limits |
| B | Refinement of the same 3D approximate integral | Its numerical convergence only |
| C | Independently formulated 3D boundary solution on manageable matched cases, plus weak-roughness limiting checks | Physical approximation error under those matched conditions |
| D | Eventually, measured channel responses with known geometry/environment | Transfer to actual deployment |

A converged 3D Kirchhoff implementation is not Level C evidence for itself. Agreement of two Kirchhoff implementations is algorithm fidelity. Existing 2D BIE results remain useful regressions but cannot replace a three-dimensional physical check. If a matched 3D reference is infeasible, report the remaining physical uncertainty and limit the claim; do not label the reference gap as a pass.

### 4.2 Statistical output and candidate generators

All stochastic and communication quantities should use one documented DSP frequency convention. Preserve the native physical response separately. With a declared direct travel-time reference tau_d, use a response such as

\[
\widetilde H(f;\xi)=H_{\rm DSP}(f;\xi)e^{i2\pi f\tau_d},
\]

where xi contains the random sea state realization and, in later stages, motion. The direct reference must be prescribed by geometry; no realization-specific fitted alignment is permitted in physical-solver comparisons.

Compare three mandatory routes:

| Route | Representation | Reason to retain it |
|---|---|---|
| P0: direct physical sampling | Generate a surface and evaluate the single-link response | Accuracy/cost baseline and a possible final solution |
| S1: conditional joint-frequency Gaussian | Estimate mean, covariance C, and, where needed, pseudo-covariance P; sample a regularized joint model | Simple baseline; tests whether second-order structure is enough |
| S2: random paths plus residual | Joint random path delays, amplitudes, and phases, with a conditional joint-frequency residual | Preferred mechanism-based hypothesis if individual arrivals are stable and identifiable |

An illustrative S2 form is

\[
\widetilde H(f;\xi)=A_d(f)+\sum_{p=1}^{P}A_p(f;\xi)e^{-i2\pi f\Delta\tau_p(\xi)}+R(f;\xi).
\]

The amplitudes must have a constrained, smooth frequency representation; otherwise arbitrary complex amplitudes can absorb any delay and make the decomposition unidentifiable. Start with one reflected cluster, then add complexity only if held-out evidence justifies it. Do not require that every rough realization has one specular point, and do not interpret sub-resolution fitted paths as uniquely observed physical arrivals.

Estimate path variables jointly and retain dependencies between delays, amplitudes, and residual statistics. If a phase/delay alignment helps fit a conditional residual, save the transformation and sample its distribution so the final generator restores the original phase and timing. Never use such alignment to hide physical-reference disagreement.

A small ensemble reflected mean does not imply a weak reflected realization. Evaluate reflected-only coherence, reflected power, and direct-to-reflected power separately. A nonzero mean does not require nonzero centered P. Joint proper/improper modeling should follow the actual covariance definitions. [E2]

Add a non-Gaussian mixture, amplitude-dependent residual, or other richer family only if both simpler candidates fail important held-out statistics or communication criteria. Deep learning and extensive surrogate training are deferred until a measured need and a sufficient independent dataset exist.

## 5. Scenario specification and sea-surface design

### 5.1 Parameter priorities

The following are research settings, not user-approved installation limits or published offset quantiles:

| Role | Geometry and sea setting | Purpose |
|---|---|---|
| First analytic control | Tx1500 m, Rx5 m, c1500 m/s, flat surface, q=0 and 0.10 | Establish deep geometry and absolute conventions |
| First numerical roughness controls | Same geometry; compact weak sinusoid and dual-direction perturbations | Executed WP1; amplitudes .001 and (.002,.001) m, not PM |
| Subsequent roughness mechanism pilot | Hm0=0.5 and 1 m crossed with q=0 and 0.10 | Deferred until finite-roughness G1 and surface definition are established |
| Nearby validation | Tx3000 m; Rx3 and 10 m; q up to 0.25; Hm0 up to 2 m | Test transfer one controlled extension at a time, then selected interactions |
| Later stress cases | Tx6000 m, q=0.50, Hm0=3 m | Extend only after the core domain and cost are understood |

Here q=rho/(z_t+z_r), rho is horizontal displacement, and the flat-reflection normal angle is arctan(q). q=0.10 corresponds to 5.71 degrees and rho=150.5 m at Tx1500/Rx5. q=0.25 corresponds to 14.04 degrees. The optional q=0.50 case is a wider-angle stress test, not part of an asserted near-vertical operating specification.

Published deep-vehicle deployments and descriptive sea-state sources motivating these choices are documented in the companion reassessment. They support realistic application context, not the exact proposed grid. [L4]

### 5.2 Define the spectrum before defining the mesh

Start with a fully developed PM candidate and define its spectral measure and units explicitly. PM is a model with assumptions, including wind conventions; a label alone does not specify a realization. [E1]

The surface specification must include Hm0=4 sqrt(m0), peak period, directional distribution, long-wave and short-wave limits, coefficient law, Fourier normalization, and seed. For a fixed PM formulation, height and period cannot be chosen independently without changing the model family. A useful implementation path is to select target Hm0, solve for the corresponding parameter of that declared PM family, and report its implied period; any finite-band adjustment must be explicit.

Use a normalized isotropic directional spectrum as a reproducibility control. Then compare a declared directional spectrum with its mean wave direction parallel and perpendicular to horizontal offset. The spreading function and concentration are proposed sensitivity settings until supported by a target environment; none should be called a measured sea condition.

Separate three lengths:

1. The parent surface domain, large enough to represent long waves and local slope/curvature.
2. The acoustic integration support, determined by contribution and convergence checks.
3. The quadrature spacing, determined by both acoustic oscillation and resolved surface geometry.

Long-wave modes larger than the integration patch can affect its local height and tilt. Do not remove them merely by generating a small periodic patch. Conversely, a large illumination footprint is not automatically the required integration aperture. Use a common parent realization when changing the integration window or numerical mesh.

The PM gravity-wave tail does not provide finite slope/curvature moments without a short-wave limit. Choose and document the physical high-wavenumber treatment before statistical production. Run nested spectral-band checks with common low-frequency coefficients. If omitted bands materially change receiver observables, extend the spectrum or narrow the model's claimed domain. Numerical mesh refinement must not silently add physical roughness.

Record actual height/slope moments and receiver clearance for every realization. If a surface crosses a fixed receiver, classify it as outside the wet-receiver idealization; report its occurrence instead of silently resampling and changing the ensemble.

## 6. Work packages, deliverables, and decision gates

### WP0 — Lock definitions and establish analytic controls

**Question:** Is the proposed single-link response represented consistently before roughness is introduced?

Specify coordinates, source units, native/DSP conversion, direct time reference, source/receiver response assumptions, and the surface spectrum. Obtain the flat 3D solution

\[
H_{\rm native}(f)=\frac{e^{ikL_d}}{4\pi L_d}-\frac{e^{ikL_r}}{4\pi L_r},\quad
L_d=\sqrt{\rho^2+(z_t-z_r)^2},\quad L_r=\sqrt{\rho^2+(z_t+z_r)^2}.
\]

Check direct, reflected, and total fields separately. For vertical geometry, reflected excess delay is 2 z_r/c; Rx5 m at c1500 m/s gives 6.6667 ms and a 150 Hz two-path interference period. This persists when Tx depth changes in the vertical, uniform, flat control.

Include a uniformly shifted flat reflector to show how random phase/delay can suppress the ensemble mean without creating diffuse multipath. Add horizontal coordinate rotation tests so that isotropic-source, isotropic-surface results do not depend on axis naming.

**Deliverable D0:** a versioned problem specification, independent analytic test cases, and a phase/CIR metadata contract.

**Gate G0:** analytic algebra and numerical representations agree within the proposed numerical tolerances; no convention fitting is used. Failure stops rough-surface interpretation.

### WP1 — Build a small receiver-focused 3D physical reference campaign

**Question:** Can the selected rough-surface approximation produce a defensible response at an acceptable cost?

Reuse the completed minimal candidate and controls recorded above rather than restarting WP1. Compact perturbations use exact infinite-flat background plus rough/flat full-density differences, with raw finite-aperture flat errors retained separately. Future noncompact random surfaces require justified boundary/support treatment. Enlarging a C2 embedding changes the physical surface; it is sensitivity, not proof of infinite-ocean convergence.

Apply 4 and 8 kHz only after the first control passes. Organize at most 12 initial independent-reference evaluations across easy, representative, and deliberately difficult cases after a one-case runtime/memory probe. This is a scope ceiling for an exploratory pilot, not a claim that 12 cases establish population validity. Repeat a selected reference with tighter resolution and support to quantify its uncertainty.

Prefer a matched 3D BIE/BEM reference where computationally feasible. Independently derived weak-roughness limits can support small-amplitude checks. Exact flat tests alone do not validate finite roughness. If the independent 3D method is itself uncertain, resolve that uncertainty before using it to rank approximations.

**Deliverable D1:** convergence report, measured cost breakdown, and a case-by-case applicability table including failures and reference limitations.

**Gate G1 and branches:**

- Numerical failure: fix integration/support/conventions before changing physical models.
- Converged Kirchhoff disagreement with an independent reference: try a justified SSA or two-scale alternative for the failing mechanism, or restrict the accepted surface domain. Do not increase FastFK rank to repair physical error.
- Adequate Kirchhoff agreement: provisionally adopt it for the tested domain and begin the stochastic pilot.
- Unaffordable or unavailable independent reference: retain a clearly conditional computational study; do not claim the target physical model is validated. Replan the reference scope before a large dataset.

### WP2 — Run a small paired mechanism experiment

**Question:** Does roughness primarily randomize a small number of path variables, or does it require a richer frequency process?

Use the four pilot conditions in Section 5: Hm0=0.5/1 m crossed with q=0/0.10 at Tx1500/Rx5. Start with eight independent base coefficient seeds. Within each seed, pair the candidate methods and conditions using common standardized coefficients where the spectral construction permits. Save the actual spectra so that coupling is not confused with an identical physical surface.

This gives 32 condition-realization responses. At the initial 4/6/8 kHz screening stage it represents 96 scalar-frequency evaluations per physical candidate, excluding controls, refinements, and independent-reference checks. These counts are budgeting definitions, not a statistically powered final sample. Use the same realizations across frequency; never redraw a surface independently at each frequency.

After the sparse checks, compute a genuinely converged broadband response for selected realizations and then the pilot set if cost permits. A first trial grid is 257 frequencies over 4–8 kHz, with spacing 15.625 Hz and periodic transform window 64 ms. Compare with 513 frequencies at fixed band for a 128 ms window and evaluate intermediate frequencies to test interpolation. Neither grid is accepted a priori. Fixed-band refinement increases the unaliased window, while nominal delay resolution remains about 1/B=0.25 ms; zero-padding only interpolates.

Use geometric path-delay estimates and observed tail energy to enlarge the window when needed. A small energy tail can still affect deep fades, so also test its communication effect. Save transform endpoints, bin ordering, taper, amplitude normalization, reference time, and carrier definition. Do not infer an absent reflected component from its falling outside a periodic CIR window.

Extract reflected power, direct/reflected ratio, reflected-only coherent fraction, phase and delay variation, frequency covariance and pseudo-covariance, band-limited delay structure, and fades. Interpret peak count with the finite bandwidth in mind.

**Deliverable D2:** paired effect plots and an explicit decision on S1 versus S2 feasibility, with uncertainty and counterexamples.

**Gate G2:** proceed only when the chosen observables and frequency representation are numerically stable. If path extraction is unstable or residual structure remains strong, favor S1 or a richer joint model instead of forcing a few-path interpretation.

### WP3 — Fit and validate the conditional static generator

**Question:** Which representation generalizes to unseen surfaces and nearby conditions at the lowest useful cost?

Use the pilot for design and tuning only. A provisional first development batch is 32 new independent base seeds, followed by 32 disjoint validation seeds across the selected core conditions. These are batch sizes for cost and precision assessment, not a guarantee of adequate covariance or tail estimation. Reserve a separate final test set after selecting the model and tolerances.

At hundreds of frequencies, small ensembles cannot identify an unrestricted full covariance. Use shrinkage, a reduced basis, or a structured kernel selected on development data. Verify positive semidefiniteness, pseudo-covariance compatibility, and sensitivity to basis dimension. Do not let low sample rank masquerade as physical low dimensionality.

For S2, regularize path count and spectral amplitude smoothness and preserve joint path/residual dependencies. For S1, assess proper and improper variants; retain P only where its estimated structure is stable and useful. Compare amplitude distributions, joint-frequency functionals, delay statistics, fades, and held-out predictive scores. Paired complex errors apply to methods supplied the same realization; an independently sampled statistical generator should be compared at the distribution/functional level, not by arbitrary sample pairing.

Test interpolation in operating conditions separately from independent-seed generalization. Hold out at least one interior condition, such as an intermediate offset, and one boundary condition. Do not claim extrapolation to Tx6000 m, strong seas, or an unseen directional family from interpolation evidence.

**Deliverable D3:** frozen generator, model card, training/validation manifest, uncertainty report, and model-selection record.

**Gate G3:** the simplest model satisfying the physical-domain restrictions and held-out statistical/communication criteria wins. If a richer model does not show a stable advantage, retain the simpler one. If all statistical candidates fail, keep direct physical sampling and diagnose the missing structure before expanding model complexity.

### WP4 — Validate communication performance and total cost

**Question:** Does replacing physical channel sampling change relevant communication conclusions within the declared tolerance?

Use single-carrier QPSK as the primary static benchmark because the repository already has MPSK infrastructure. Propose fc=6 kHz, symbol rate 2000 symbols/s, root-raised-cosine roll-off 0.25, and uncoded raw rate 4000 bit/s. The nominal occupied band is about 4.75–7.25 kHz, within the preferred band. These are research choices, not specifications for the unknown transducer. Add a wider-band waveform sensitivity test only after the first benchmark is trustworthy.

Establish two receiver modes: a perfect-CSI/timing control and a practical pilot-based channel estimator with a documented equalizer and synchronization procedure. Derive sample-rate, pulse-filter, noise-variance, Es/N0 and Eb/N0 mappings explicitly. Validate AWGN, flat two-path, and noiseless known-channel cases before evaluating new channel models. The old processing-related BER floor must be resolved or the benchmark rebuilt independently; it must not serve as the channel acceptance tolerance.

For a concrete starting configuration, propose 1024 payload QPSK symbols per packet, a 128-symbol known pilot, two complex samples per symbol, and a 10-symbol pulse-filter span. Start with a fractionally spaced linear MMSE equalizer; determine and freeze its length from the accepted channel-memory and pulse-support study, then check a longer filter. Use a regularized pilot-based impulse-response estimate in the practical receiver. Set guard duration from the converged channel support rather than a historical CIR window. If linear equalization is an inadequate benchmark under relevant spectral nulls, evaluate a documented DFE extension and keep its settings identical across generator comparisons. Report net rate including pilot and guard overhead. These packet choices belong to the explicitly static first study; the later dynamic study must test variation within a packet rather than assume it remains frozen.

Use two distinct experiments:

1. **Generator fidelity:** with the same receiver algorithm and noise policy, compare performance distributions under physical and generated channel ensembles.
2. **Receiver-model mismatch:** evaluate receiver designs or priors on the same held-out physical channels, transmitted data, and pilot/noise realizations. Keep a perfect-CSI control to separate estimation losses.

For fade and outage tests, retain channel gain under a fixed transmitted-energy/noise convention. A per-realization normalization experiment may isolate frequency-selective distortion, but must be labeled secondary because it removes part of the fading effect.

Measure surface generation, source construction, scattering, broadband assembly, fitting, validation, condition updates, and sample generation separately and end to end. Include first-call startup and warmed execution as separate categories. No claimed single-link speedup may use an unnecessary dense-field baseline.

**Deliverable D4:** BER/SNR and packet-error comparisons with uncertainty, fixed-energy fade/outage comparisons, and complete cost/break-even curves.

**Gate G4:** acceptance requires the proposed communication tolerances with adequate Monte Carlo precision, plus no unresolved physical/statistical failure in the claimed domain. Speed alone cannot compensate for incorrect tails or communication behavior.

### WP5 — Establish the domain and deployment sensitivity

Expand selected factors after the core model is stable: Tx3000 m, Rx3/10 m, q up to 0.25, Hm0 up to 2 m, and wave direction relative to the offset. Prioritize interactions indicated by WP2 rather than an exhaustive Cartesian product. Use a small space-filling or deliberately selected boundary design once the range is frozen.

Then assess source directivity, receiver pattern, absorption, plausible bottom clearance, and representative sound-speed profiles. These additions may change the preferred propagation tool; PE, rays, or a hybrid can be appropriate here after source/phase consistency is checked. The first-stage analytic homogeneous propagator is a deliberate controlled starting point, not a universal replacement for PE.

If deployment geometry or measured channel data become available, compare trends and distributions under matched conditions. Without them, the strongest defensible deliverable is a simulation-validated idealized model with an explicit domain, not a sea-trial-validated shipboard model.

**Deliverable D5:** applicability map, sensitivity ranking, update-cost map, and clear boundaries between validated, interpolated, stress-tested, and untested conditions.

### WP6 — Extend to temporal evolution and Doppler

Represent an initial linear evolving sea using consistent modal coefficients, for example

\[
\eta(x,y,t)=\sum_m a_m\cos(\boldsymbol K_m\cdot(x,y)-\Omega_m t+\phi_m),
\qquad \Omega_m^2=gK_m
\]

in the deep-water gravity-wave regime. Finite-depth dispersion uses g K tanh(K h); the relevant water depth and spectral regime must be checked. Capillary or other short-wave extensions require their corresponding dispersion model. [E3]

Retain coefficients across time and evolve phases; add source/receiver trajectories separately. Compare deterministic single-wave and constant-velocity controls before a random time series. For slow motion, a path has a Doppler contribution approximately -(f/c) dL/dt under the chosen DSP convention. This is a diagnostic relation, not a proof of a frozen-boundary approximation.

At kilometer depths, acoustic flight time can be appreciable relative to surface evolution. Define whether t denotes transmit, boundary-interaction, or receive time. Test the error of frozen-surface snapshots against a retarded-interaction or suitable moving-boundary reference on simple controls before assuming the whole path sees one instantaneous sea surface. Refining the snapshot interval alone does not establish that approximation.

Validate temporal correlation, path-dependent Doppler, amplitude–Doppler dependence, and communication performance under temporal sampling refinement. Include common Doppler/time scaling from platform motion and residual path variation where relevant. A short coherence time does not justify independent redraws; a small delay spread does not establish WSSUS.

Begin with independent continuous records spanning many dominant wave periods; extend record length based on effective independent information. Split validation by record and seed, not by neighboring snapshots. Fit a temporal stochastic state model only after the evolving physical reference is adequate.

**Deliverable D6:** a continuous-time-compatible generator, temporal validity domain, and a justified update rule.

**Gate G6:** static fidelity must survive temporal extension; correlation and Doppler must agree within predefined uncertainty. If the frozen-boundary or temporal state model fails, retain the static model as the completed first-stage result while revising the dynamic formulation.

## 7. Experimental units, replication, and uncertainty

The independent static unit is the parent surface realization at a specified environmental condition. Frequency bins, receiver samples, repeated timing calls, and repeated packets on one surface are not additional independent sea realizations. For time-varying work, independent records/seeds are the primary units; snapshots within a record are correlated.

Use common random inputs for paired algorithm comparisons and compatible factor changes. Randomize method order within a timing block and block by machine/session. Keep numerical refinement on exactly the same physical surface. Separate seeds for surface coefficients, communication data/noise, and run ordering; record all streams in a manifest. No schedule or new experiment is generated by this planning document.

Pilot batches are exploratory. Estimate uncertainty at the parent-realization level using an appropriate clustered analysis; paired differences are preferable where pairing is real. Repeated packets improve conditional BER estimation, not environmental coverage. Regularized covariance comparisons must acknowledge finite-sample bias.

After the pilot, freeze primary metrics, thresholds, model choices, and a final test size using pilot variability and the available budget. Do not repeatedly inspect one final test set and stop at a favorable result. If adaptive sampling is necessary, predefine the stopping rule or use suitable sequential intervals. Rare-outage claims require enough independent channel events; 32 or 64 surfaces cannot certify extreme tail probabilities. Report unresolved precision as inconclusive, not as a pass.

## 8. Cost model and resource decisions

For one operating condition and M requested channel realizations, compare

\[
T_{\rm physical}(M)=T_{\rm setup,P}+M T_{\rm sample,P},
\]

\[
T_{\rm statistical}(M)=T_{\rm data}+T_{\rm fit}+T_{\rm validation}+T_{\rm update}+M T_{\rm sample,S}.
\]

If sample costs are stable and T_sample,P exceeds T_sample,S, the break-even count is

\[
M_* = \max\left(0,\frac{T_{\rm data}+T_{\rm fit}+T_{\rm validation}+T_{\rm update}-T_{\rm setup,P}}
{T_{\rm sample,P}-T_{\rm sample,S}}\right).
\]

Use matched scopes and avoid double-counting shared validation costs. Report both research-development cost and recurring deployment cost, and how the result changes when conditions require rebuilding. If the denominator is nonpositive, the statistical route has no online speed advantage under those measurements.

Report total time, amortized time, memory peak, batch size, hardware, precision, startup, and uncertainty. Benchmark workloads of 1, 100, and 10,000 generated channels where feasible; these are comparison workloads, not user throughput requirements. A target such as tenfold warmed sampling speedup may be useful as a research aspiration, but total break-even and adequate fidelity determine usefulness.

The original 2D timing cannot budget 3D BIE. Before any expensive batch, run one candidate and one independent-reference case, forecast memory and total cost, and choose a bounded checkpointed campaign. Stop scaling when the reference is unaffordable, accuracy is inadequate, or the intended generation volume falls below break-even. Then reduce the experimental scope or change the formulation explicitly rather than silently relaxing convergence.

## 9. Proposed acceptance criteria

These are initial design tolerances proposed by this plan, not validated universal thresholds or previously confirmed user requirements. Freeze or revise them after the pilot and communication sensitivity study, before final testing. Report all failures as well as passes.

| Layer | Proposed initial criterion | Interpretation and safeguard |
|---|---|---|
| Analytic identities / convention conversion | Approximately 1e-10 normalized discrepancy for well-conditioned double-precision algebraic controls | Does not impose this tolerance on a discretized boundary integral |
| Numerical integration / broadband representation | Initially target 1e-3 change in normalized band response under refinement | Separate direct and reflected components; also track phase/delay and communication changes |
| Independent physical comparison | Use 1% band-relative reflected-field error as a provisional screening level | Not sufficient alone; inspect spectral structure, fades, and communication sensitivity; reference uncertainty should be substantially smaller |
| Statistics | Initial screening: mean reflected power within 0.5 dB and stable covariance/functionals within roughly 10% | Use uncertainty-aware metrics and an absolute floor when a denominator is near zero; tail failures can reject a model despite passing averages |
| Primary communication | Absolute SNR shift no greater than 0.5 dB at uncoded QPSK BER=1e-3 | Proposed criterion; report intervals and require uncertainty appreciably below 0.5 dB; do not extrapolate beyond the measured crossing |
| Packet-error sensitivity | Absolute PER difference no greater than 0.02 near a declared reference PER of 0.10, using a frozen packet size | Secondary proposal; specify overhead/conditioning and cluster uncertainty by channel |
| Cost | Document break-even and useful throughput under measured scopes | No predetermined single-link speedup is treated as achieved |
| Time variation | Tolerances set from Doppler/correlation sensitivity before final time-series testing | No universal Hz, coherence-time, or update-interval value is assumed |

Normalize same-realization complex errors carefully. Report component-relative metrics and an absolute scale based on the declared direct free-field band norm. Never divide only by the total field at deep interference nulls. Conversely, a small total-field error dominated by the direct path must not hide a wrong reflected component.

Phase and delay estimates near zeros require amplitude-aware interpretation; do not discard those frequencies from all tests, because the nulls may matter to communication. Use integrated metrics, raw response plots, and communication checks together. Demonstrate that numerical, reference, and Monte Carlo uncertainty do not consume the claimed tolerance.

## 10. Conditional schedule and the next executable task

The estimates below are effort bands for planning, not measured runtimes or a delivery promise. Reference development is the main uncertainty. Gates, not calendar dates, govern progression.

| Work package | Indicative focused effort | Prerequisite |
|---|---|---|
| WP0 | 3–5 working days | Current information is sufficient for idealized controls |
| WP1 | 2–4 weeks, potentially longer for an independent 3D reference | G0 |
| WP2 | 1–2 weeks after measured cost is known | G1 and a frozen surface definition |
| WP3 | 2–3 weeks plus reference-data cost | G2 |
| WP4 | 1–2 weeks; basic communication sanity checks can begin during WP1 | Reliable physical/statistical inputs and corrected benchmark |
| WP5 | 2–3 weeks for a selected extension set | G3/G4 |
| WP6 | Separate 3–6+ week study after the static result | Static model accepted and dynamics specification available |

The first executable task should implement **WP0 and the minimal WP1 candidate**, not the entire plan:

1. Create the problem/metadata specification and analytic 3D flat controls.
2. Implement the single-receiver 3D integral under a documented boundary approximation, with analytic incident field and derivatives.
3. Validate flat and weak deterministic two-coordinate surfaces, and separate numerical convergence from physical-reference checks.
4. Record one realistic cost/memory breakdown and identify a feasible independent-reference case.
5. Produce a go/no-go report for the 32-response pilot; do not launch that pilot until its physical domain and cost are defensible.

This task should produce a reference-quality research prototype, not a production replacement or a claim of fast statistical channel generation.

## 11. Deliverables and repository organization

The static research is complete only when D0–D5 provide a reproducible physical domain, a selected generator or direct-sampling solution, held-out statistical and communication evidence, and full cost accounting. D6 is a distinct later deliverable.

For future implementation, follow existing repository responsibility boundaries:

- Independent reference solvers and validation-only helpers: `scripts/validation/support/`.
- Invariants, analytic controls, and acceptance workflows: `scripts/validation/`.
- Controlled route comparisons: `scripts/comparisons/`.
- Surface/condition ensembles: `scripts/experiments/`.
- Saved-data figures and summaries: `scripts/reporting/`.
- Validated reusable components, when promoted: appropriate `src/surface/`, `src/receiver/`, `src/statistics/`, or `src/channel/` modules.
- Case-specific data: `results/validation/`, `results/comparisons/`, or `results/experiments/`; conclusions in `reports/`.

Persist raw native complex fields, DSP responses, frequency axes, phase metadata, realized surface coefficients, source definition, reference settings, numerical guards, seed manifests, timing scopes, fitted parameters, and train/test membership. The three essential status labels are algorithm fidelity, physical validity, and statistical/communication adequacy.

Preserve the public `vertical_channel_model` interface and `H_f = H_direct_f + H_reflect_f`. New production integration is a separate stage requiring the repository's scalar-frequency and wideband regressions. Use the shared MATLAB MCP session for execution. Do not recover obsolete evidence from `cash/`.

## 12. Unresolved inputs and how the plan handles them

| Missing input | Current treatment | Needed before |
|---|---|---|
| Source model, calibration, directivity | Ideal normalized point-source controls; explicit sensitivity later | Device-specific received power and range claims |
| Hull clearance and platform geometry | Idealized hull-free domain, with no transfer claim | Shipboard validation |
| Maximum depth and offsets | Proposed nested research settings, not engineering limits | Final applicability envelope |
| Sea-state distribution and short-wave model | Explicit PM candidate, cutoff/direction sensitivity | Population weighting or ocean-wide claims |
| Water depth and source bottom clearance | No bottom in first half-space study | Claims about a bottom-mounted installation |
| Computing and update budgets | Measured cost curves and break-even workloads | Optimization target and deployment choice |
| Measured channel data | Independent numerical reference plus explicit measurement gap | Field-validated model claim |
| Required modem and tolerance | Proposed QPSK benchmark and tolerances above | Final engineering acceptance |

This uncertainty register does not require another clarification round before analytic work. It prevents provisional choices from becoming hidden assumptions in later conclusions.

## 13. Evidence, references, and planning assistance

### Local evidence

- **L1.** [Research background and handoff](deep_sea_channel_research_background_for_new_ai.md): confirmed requirements and full historical evidence register.
- **L2.** [Single-receiver timing report](fast_fk_single_receiver_timing_report.md): fixed 4 kHz timings, source fixture, repetition protocol, and exclusions.
- **L3.** [Multifrequency consistency report](spectral_fast_kirchhoff_multifrequency_consistency_report.md): ranks 3/12, phase conventions, BIE anchors, and diagnostic CIR limitations.
- **L4.** [External-reply reassessment](deep_sea_channel_external_reply_reassessment.md): corrected geometry, PM arithmetic, statistical interpretation, published deployment examples, and scenario provenance.
- **L5.** [Technical guide](../vertical_comm_guide.md), [2D FullFK reference](../scripts/validation/support/solve_full_kirchhoff_2d.m), and [WAPE implementation](../src/propagation/vertical_wape_propagator.m): existing model definitions, numerical conventions, and boundaries.
- **L6.** [FullFK PM validity boundary](full_kirchhoff_pm_validity_boundary_report.md): conditional physical approximation limits; not a universal cutoff or 3D certification.

### Primary and institutional sources

- **E1.** Pierson and Moskowitz (1964), [A proposed spectral form for fully developed wind seas based on the similarity theory of S. A. Kitaigorodskii](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/JZ069i024p05181). Publisher abstract/record checked; supports the fully developed-sea and wind-convention context, not a new numeric cutoff in this plan.
- **E2.** Neeser and Massey (1993), [Proper Complex Random Processes with Applications to Information Theory](https://www.isiweb.ee.ethz.ch/archive/massey_pub/pdf/BI437.pdf). Original author-hosted paper; used for centered covariance/pseudo-covariance semantics, not a claim that the target channel is Gaussian.
- **E3.** MIT Marine Hydrodynamics, [Solution of the Dispersion Relationship](https://web.mit.edu/fluids-modules/www/potential_flows/LecturesHTML/lec19bu/node5.html). Institutional teaching source for the linear gravity-wave dispersion relation and its deep-water limit; not a validation of the proposed moving acoustic boundary model.
- **E4.** Qarabaqi and Stojanovic (2013), [Statistical Characterization and Computationally Efficient Modeling of a Class of Underwater Acoustic Communication Channels](https://doi.org/10.1109/JOE.2013.2278787), IEEE Journal of Oceanic Engineering 38(4), 701–717. The accessible abstract describes a physical/statistical channel formulation. The author-hosted full PDF could not be fetched during this turn, so no detailed theorem, validity region, or numeric performance claim is imported. Read the full paper before claiming novelty or adopting its assumptions for this geometry.

Sources were checked on 2026-10-09. Equations for free-space propagation, image geometry, proposed cost accounting, and planning values are explicitly given above; none is a newly measured project result.

### Planning assistance

The local experimental-design skill informed independent-unit definitions, paired seeds, timing blocks, and separation of pilot/development/final validation. Its provenance is Kassis, T., Agarwal, V., He, Y., Patel, D., and Brueckner, A. M. (2026), [Scientific Agent Skills: A Library of Procedural Knowledge for Research Agents](https://doi.org/10.48550/arXiv.2609.00065). The arXiv author/year record was checked; the current record is revised September 2026. This is a workflow acknowledgment, not evidence of acoustic validity. No acoustic-formula-audit skill workflow was used for this plan.
