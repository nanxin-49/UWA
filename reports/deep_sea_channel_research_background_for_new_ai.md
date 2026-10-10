# Deep-Sea Acoustic Channel Research: Background, Requirements, and Request for an Independent Reassessment

Prepared: 2026-10-09 (Asia/Shanghai).

Follow-on document: [recommended complete research plan](deep_sea_channel_research_plan.md). It provides a preferred route, work packages, proposed experiments and tolerances, and decision branches. It is a planning recommendation, not new validation evidence or a replacement for the confirmed requirements below.

Revised the same day after reviewing an external AI proposal and receiving further user clarification: the receiver is cable-deployed, the source transducer model is unknown, and offset/sea-state limits have no engineering specification yet. The user requested published scenario examples. Section 8 now separates field evidence, provisional test settings, and a suggested decision sequence. The companion [external-reply reassessment](deep_sea_channel_external_reply_reassessment.md) records the detailed corrections; this handoff remains self-contained.

Purpose: provide a self-contained handoff to a new AI researcher who may not have access to the repository. This document records the user's clarified objectives, summarizes existing evidence, and asks for a fresh assessment of the research direction. It does not prescribe PE, FullFK, FastFK, or any previously proposed roadmap.

Provenance: the requirements below were explicitly confirmed by the user in the planning conversation. Experimental statements come from the existing reports and implementation excerpts listed in Section 9. The sources were inspected during preparation; no new acoustic simulation, timing experiment, or statistical validation was performed. This is an internal research handoff, not an independently reviewed publication.

## 1. The actual research objective

**User-confirmed requirement.** Develop a physically credible statistical acoustic channel model and a fast channel generator for communication from a deep transmitter to a receiver near the sea surface. The eventual purpose is communication research, rather than spatial sound-field visualization or accelerating a particular solver for its own sake.

The desired output is the complex channel response of a single Tx/Rx link, its relevant statistical structure, and an efficient way to generate channel realizations. Full spatial fields can be useful for diagnosis, but are not the primary application output.

| Aspect | Confirmed requirement | Not yet fixed |
|---|---|---|
| Transmitter | Extend to transmitter depths exceeding 1000 m | Maximum depth, bottom clearance, source pattern and normalization |
| Receiver | Depth 3–10 m; cable-deployed from a ship or another surface platform | Position and clearance relative to the hull, platform geometry, receiver motion and directivity |
| Link geometry | Near-vertical, with horizontal displacement allowed | Displacement range and resulting angular coverage |
| First research stage | Three-dimensional propagation over a two-dimensional rough sea surface; static single-link realizations | Final numerical representation and approximation level |
| Sound speed | Uniform for the first stage | Representative value and later environmental extensions |
| Later stages | Evolving sea surface, channel time correlation, and Doppler | Motion models, update interval, observation duration, and Doppler range |
| Frequency | 4–8 kHz is the preferred starting band | Another representative band may be proposed with justification |
| Sea surface | Pierson–Moskowitz (PM) is currently considered a reasonable starting model | Sea-state coverage, directional spectrum, spectral limits, and alternatives for time evolution |
| Computational strategy | A limited amount of offline physical simulation or training is acceptable | Offline budget, online generation target, memory budget, and update cost when conditions change |
| Acceptance | Primarily communication-performance tolerances, supported by physical and channel-statistical checks | Numerical tolerances and the communication benchmark |
| Communication system | No mandatory modulation, rate, synchronization, equalization, or coding scheme has been selected | The new study should propose a fair benchmark and its tolerances |
| Method selection | Fully open to reassessment | No requirement to retain PE, FK, a Gaussian generator, or a particular compression method |

Cable deployment defines the intended installation but does not establish that the receiver is directly below the hull or sufficiently far from it. Hull, wake, bubble, noise, and motion effects remain questions for physical prioritization. The source transducer model is currently unknown. An isotropic point source may be a labeled analytic control, but must not be presented as the real hardware specification.

The existing 100 m tests are not the final application. In particular, the project contains a 100 m water-depth Bellhop visualization with Tx at 99 m and Rx at 3 m, while the independent FullFK/FastFK studies use a different Tx100/Rx3 research fixture. Neither demonstrates the required deep-water, three-dimensional application. [R1, R2]

## 2. Quantities of interest and terminology

Depth is positive downward, and the mean sea surface is at zero depth. In the target three-dimensional problem, the rough surface is represented by a height/depth function over two horizontal coordinates. A static realization is one fixed surface; it is not a time history.

For a fixed environment and geometry, write the receiver response as

\[
H(f;\eta)=H_{\mathrm{direct}}(f)+H_{\mathrm{reflect}}(f;\eta),
\]

where \(\eta\) denotes the surface realization. This is the direct-plus-surface-reflection decomposition used by the existing framework, not a decision that all other paths or physical effects are negligible in the future application. Any additional contribution must be justified and incorporated consistently.

Relevant statistical outputs may include the mean, cross-frequency covariance, and pseudo-covariance:

\[
\mu_i=\mathbb E[H(f_i)],\qquad
C_{ij}=\mathbb E[(H(f_i)-\mu_i)(H(f_j)-\mu_j)^*],
\]

\[
P_{ij}=\mathbb E[(H(f_i)-\mu_i)(H(f_j)-\mu_j)].
\]

The pseudo-covariance describes complex statistics not represented by covariance alone; it should not be dropped without evidence. These definitions do not assume that the channel is Gaussian or that second-order statistics are sufficient for communication performance.

Specify what the ensemble conditions on: geometry, source pattern, sea-state parameters, surface directionality, and phase/time reference. Report direct and reflected components separately where possible. A small ensemble mean can result from a strong reflected path whose phase or delay varies between realizations; it does not prove that the path disappears, that many diffuse paths exist, or that a Gaussian model fails. Likewise, a nonzero mean does not force nonzero centered pseudo-covariance. These distinctions matter when comparing random few-path, proper/improper Gaussian, and non-Gaussian generators. [R13]

A channel impulse response (CIR) is obtained from a consistently defined complex frequency response. Its delay coordinate must have a stated reference, such as excess delay relative to the direct path. For future time variation, observation time and propagation delay must be distinguished, for example by \(H(f,t)\) and \(h(\tau,t)\).

Key terms used in the existing work:

- **PE/WAPE:** parabolic-equation / wide-angle parabolic-equation propagation methods used by the existing production framework.
- **Local phase screen:** a surface model that applies a local complex multiplier to an incident field on a reference plane.
- **FullFK:** the project's direct evaluation of a Kirchhoff boundary integral on the actual rough profile. “Full” means retaining the specified geometric and integral terms; it does not mean an exact solution of all scattering physics.
- **FastFK:** the existing low-rank and FFT-based approximation to that FullFK operator, designed to produce a dense receiver-line field.
- **BIE:** an independent boundary-integral-equation reference for the specified two-dimensional Helmholtz problem, subject to numerical convergence checks.
- **Canonical PM:** one designated, fixed PM-derived realization, not a universal or statistically representative ocean surface.

## 3. Existing research framework and why it changed direction

### 3.1 Propagation, statistics, and communication tools already exist

**Existing evidence.** The MATLAB repository contains a PE/WAPE channel framework, explicit rough-surface models, joint cross-frequency surface statistics, cached propagation, an exact discrete receiver adjoint, conditional statistical generators, and MPSK/QPSK communication tools. [R1]

For the supported fixed-environment PE operator, the receiver adjoint can replace repeated surface-to-receiver propagation by a contraction. If \(A_i\) maps a surface-plane field to the receiver plane and \(r\) selects the receiver, then

\[
r^H A_i\psi_{\mathrm{surface},i}
=q_i^H\psi_{\mathrm{surface},i},\qquad q_i=A_i^H r.
\]

The implemented numerical equivalence has been checked under specified conditions: uniform medium, CPU double precision, fixed geometry and grids, nearest-grid-point receiver sampling, and no bubble or Doppler effects in this validation chain. This does not establish equivalence to the independent FullFK rough-boundary operator. [R1]

The existing statistical route includes a local phase-screen residual of the form

\[
\delta G_i=R_{0,i}\exp(i\alpha_i\eta)-R_{\mathrm{coh},i},
\]

and contractions of its cross-frequency covariance and pseudo-covariance to receiver statistics. Conditional generators and communication comparisons have been reported for fixed wind-speed nodes U=5 and U=8 m/s and a 64-frequency axis over 4–8 kHz. These are capabilities within the earlier model and registered conditions, not validation of the new deep-sea target. [R1]

**Interpretation.** There is useful infrastructure to reuse, but numerical agreement within a model does not certify the model's physical accuracy. Replacing its reflection physics may require a new statistical derivation or new reference samples.

### 3.2 Local reflection approximations motivated a different physical reference

**Existing evidence.** Controlled comparisons used common incident data and corrected phasor conventions to compare a local model, a vertical-wavenumber-aware correction, FullFK, and BIE. Some rough-surface cases showed appreciable residual disagreement for the local models, while FullFK was substantially closer to the BIE reference. The cited report concerns two-dimensional, uniform-medium, 4 kHz tests with a specified Gaussian angular-spectrum source. [R3]

The FullFK implementation evaluates

\[
u_{\mathrm{ref}}(\mathbf X)
=-\int_{\Gamma}\left[
G(\mathbf X,\mathbf s)\,\partial_nu_{\mathrm{inc}}(\mathbf s)
+\partial_{n_{\mathbf s}}G(\mathbf X,\mathbf s)\,u_{\mathrm{inc}}(\mathbf s)
\right]\,ds.
\]

Here \(\Gamma\) is the actual rough curve in the present two-dimensional implementation, the normal points into the water, and the native time convention is \(\exp(-i\omega t)\). The operator retains the actual surface position, normal, arc-length measure, incident normal derivative, and both Green-function terms. This is the existing implementation's convention-specific formula, not a proposed universal three-dimensional replacement. [R3, R4]

**Interpretation.** These tests motivated retaining the coherent contribution of the actual rough boundary instead of relying solely on a local phase multiplier. They do not establish that all local approximations fail in every relevant regime, or that FullFK is the best final method.

### 3.3 FullFK itself has a conditional physical validity domain

**Existing evidence.** A later PM study varied spectral content while retaining shared Fourier coefficients and the same source, geometry, receiver metric, and frequency. In a controlled continuation of one realization at 4 kHz:

| Underlying Fourier cutoff | Reported baseline weighted FullFK–BIE error | Guarded classification against the historical threshold |
|---|---:|---|
| \(K_{\max}=0.70685835\) rad/m | \(8.5275345\times10^{-4}\) | Pass |
| \(K_{\max}=0.74612826\) rad/m | \(1.1974915\times10^{-3}\) | Fail |

The associated numerical refinement indicators were approximately \(6.82\times10^{-6}\); the report also lists refined-solve errors separately. The classifications use the historical \(10^{-3}\) criterion in a frozen central receiver region. At the same cutoff, changing the amplitudes of added modes changed acceptance, so cutoff alone was insufficient to define a universal validity rule. [R5]

**Boundary.** This is a single-realization, finite-aperture, two-dimensional conditional result. It is not a population failure probability, a universal PM cutoff, or validation for transmitter depths exceeding 1000 m. Increasing FastFK rank cannot remove the physical approximation error of FullFK.

## 4. What the existing FastFK actually accomplished

### 4.1 Two different compression problems

The incident compression separates the height dependence of the angular-spectrum reconstruction, including the residual factor

\[
\exp\!\left(i(k-k_z)\eta\right),\qquad k_z=\sqrt{k^2-k_x^2},
\]

while retaining the strong height carrier. The outgoing compression concerns the height dependence of Green-function lag kernels and their derivatives. FFT convolution then produces the field on an aligned uniform receiver line. An incident rank of three is not the rank of the entire reflection operator. [R6, R7, R8]

The direct outgoing integral has a receiver-by-surface cost proportional to \(N_rN_s\), where \(N_r\) and \(N_s\) are receiver and surface sample counts. In the dense-line benchmark, both increased together. At one receiver, that same direct outgoing sum instead scales with \(N_s\). These statements concern the outgoing calculation; they exclude incident reconstruction and other preparation.

### 4.2 Accuracy and dense-line acceleration are real but task-specific

**Existing evidence.** The deterministic prototype selected incident/outgoing ranks 3/6 and reported approximately 140.86-fold apply-only acceleration at 4097 surface and 4097 receiver samples; including preparation gave approximately 14.73-fold acceleration in that benchmark. [R7]

A PM pilot used three fixed coefficient realizations and spectral-cutoff controls. Ranks 3/6 failed its algorithm criterion; the tested common choice 3/8 passed the algorithm comparisons. A high-cutoff case still failed the FullFK–BIE physical criterion. The pilot is not a sea-state ensemble validation. [R9]

The subsequent broadband study held one PM realization fixed and used a common rank 3/12 across 65 frequencies over 4–8 kHz. The maximum weighted FastFK–FullFK complex error was \(2.30669959\times10^{-4}\); the maximum whole-grid error was \(2.87485157\times10^{-4}\). BIE was evaluated at three anchor frequencies, 4, 6, and 8 kHz, not at all 65 frequencies. Ranks 3/8 failed from 6 kHz, while 3/10 narrowly failed at 8 kHz. The selected rank is a tested configuration, not a universal minimum. [R10]

### 4.3 The geometry and surface behind these results

The broadband and later single-receiver tests used the following frozen fixture: [R9, R10, R11]

- Two-dimensional propagation with one transverse coordinate, uniform sound speed 1500 m/s, Tx depth 100 m and Rx depth 3 m.
- A saved Gaussian angular-spectrum source with sigma 0.3 m and 984 source samples; this is not automatically a calibrated three-dimensional point source.
- Canonical PM coefficient seed 260001, 12 underlying Fourier modes, and actual underlying cutoff \(K_{\max}=0.471238898\) rad/m.
- A fixed smooth C2 embedding envelope equal to one for \(|x|\le42\) m and zero for \(|x|\ge50\) m, inside the integration interval \([-60,60]\) m. The embedded surface should not be confused with an unmodified, infinite, band-limited PM ocean.
- 4097 surface samples over 4096 intervals, giving \(\Delta x=120/4096=0.029296875\) m. Refining this grid while holding Fourier coefficients fixed changes numerical sampling, not physical sea-state bandwidth.
- A frozen central receiver footprint for the main weighted physical metric. Whole-grid algorithm fidelity is distinct from physical certification of all edge receivers.

### 4.4 Broadband consistency is not complete CIR validation

The native FK phasor convention differs from the IFFT-ready transfer-function convention. The broadband validation used a documented whole-envelope conversion and a common direct-path time origin, without empirical phase fitting or peak alignment. The primary Full/Fast comparison gave a maximum phase difference of approximately \(1.74394\times10^{-4}\) rad, a maximum group-delay difference of approximately 70.19 ns, and a reflected CIR relative complex error of \(7.17177\times10^{-5}\). [R10]

The 65-frequency spacing is 62.5 Hz, corresponding to a 16 ms periodic IFFT window. The report explicitly treats the CIR as a diagnostic: complete late-arrival support and absence of time aliasing have not been established. Agreement of two methods after the same finite-band transform is not sufficient evidence that either has captured the entire physical impulse response. [R10]

## 5. The single-receiver test and why a fresh direction is needed

### 5.1 Direct experimental comparison

**Existing evidence.** The 2026-10-08 test reused the frozen 4 kHz fixture above. FullFK was given exactly one receiver at x=0 m, z=3 m and formed a 1-by-4097 integration kernel. FastFK retained its 4097-receiver computation and then extracted the same point. [R11]

Both methods were warmed up. Nine rounds alternated execution order; each round rebuilt the preparation and timed a batch of ten evaluations. The following are medians of the corresponding measurements:

| Timed stage | FullFK, strictly one receiver | Current FastFK, whole line then extraction |
|---|---:|---:|
| Preparation | 110.4937 ms | 178.9923 ms |
| Evaluation per call | 0.98571 ms | 36.0479 ms |
| Preparation plus one evaluation | 111.55524 ms | 214.69579 ms |

The combined row is the median of per-round totals; it need not equal the sum of the separate medians. The source report defines speedup as Full time divided by Fast time. Its values are 0.02734445 for evaluation and 0.51959677 including preparation: **FullFK is faster by approximately 36.57 and 1.925 times**, respectively. [R11]

The native complex reflected responses were:

```text
FullFK: -0.1017881599952594 + 0.05510467672944143i
FastFK: -0.1017948747443700 + 0.05508535879071448i
```

The pointwise relative complex difference was \(1.76692919\times10^{-4}\), or approximately 0.0177%. Thus, the single-point test confirms fidelity at this point while reversing the performance ranking for this workload. [R11]

### 5.2 What the timer includes and excludes

FullFK preparation reconstructs the incident field and its spatial derivatives on the actual surface and assembles the normal derivative and receiver configuration. FastFK preparation builds its incoming/outgoing bases, lag kernels, interpolation, and plan. These are different preparation procedures. [R11]

The benchmark excludes file I/O, hash verification, reporting, and historical sea-surface generation. Common configuration assembly was separately recorded. It did not execute PE propagation, new BIE calculations, a broadband timing sweep, or a receiver-count sweep. “Preparation plus one evaluation” measures rebuilding plans in an already warmed implementation, not first MATLAB startup or first JIT compilation. [R11]

### 5.3 Interpretation, not an imposed new algorithm

The existing FFT acceleration solves a dense receiver-line problem efficiently. A single-link communication task does not inherently require that output. Previous large dense-line speedups therefore cannot be used as single-link speedups.

Using the separate median preparation and evaluation times, incident preparation accounts for roughly 99% of the timed FullFK two-part workload. This identifies a candidate bottleneck in this implementation, not proof of the bottleneck after real PE input, three-dimensional geometry, broadband execution, or cross-realization reuse.

Possible investigation targets include input reuse, incident reconstruction, receiver-specific evaluation, direct statistical calculation, or offline statistical modeling. None is selected by this handoff. It remains possible that direct physical evaluation is already cheap enough for a useful reference dataset, or that another physical formulation is preferable.

The test does **not** show that every fast single-receiver method must lose to FullFK. It also does not establish a receiver-count crossover, a three-dimensional cost estimate, or a speedup for a new hybrid algorithm.

## 6. Evidence that must not be generalized

| Established within stated conditions | Not established by that evidence |
|---|---|
| Exact discrete PE adjoint and statistical contractions | Physical equivalence to a different rough-boundary model |
| Local-model conditional generators and communication comparisons | Accuracy of the target deep-sea channel |
| FullFK–BIE agreement in selected cases | Universal Kirchhoff validity across sea states, depths, and dimensions |
| FastFK–FullFK agreement | Independent physical correctness of either solver |
| Dense-line speedup | Single-link, broadband, or end-to-end speedup |
| One canonical PM broadband case | Population robustness over independent sea realizations |
| Consistent finite-band CIR diagnostics | Complete, converged, unaliased delay support |
| Static independent channel samples | Correct time correlation, evolving sea realizations, or Doppler |
| Ranks, aperture, and sample counts sufficient for old tests | Appropriate settings for the new three-dimensional deep-water target |

The historical complex-field threshold of \(10^{-3}\) is not a user-mandated future acceptance threshold. The new objective prioritizes communication-performance tolerances, but an insensitive or flawed communication benchmark must not conceal physical or statistical errors.

For example, the existing two-node QPSK study reports a noise-free BER floor of approximately 0.2%–0.5% associated with its communication processing. It is therefore not an unquestioned acceptance oracle for the new study. [R12]

Similarly, the label “PASS” always refers to the stated test and criterion. It should not be interpreted as approval of the entire future research program.

## 7. Open questions for the new research assessment

These are deliberately unresolved, rather than silently filled with historical test settings.

1. **Physical scope:** Which effects materially influence the intended deep-to-near-surface link? Assess source/receiver geometry, finite illumination, rough-surface scattering, relevant additional paths, and platform-related effects. Explain which can be omitted initially and why.
2. **Environmental envelope:** What maximum depth, offset range, sea-state range, directional spectrum, and spectral limits should define the first credible study? PM is a candidate, not an instruction to force all conditions into one spectrum.
3. **Reference hierarchy:** Which analytic controls, independent solvers, convergence tests, and eventual measurements can support three-dimensional validity without making the reference campaign prohibitively expensive?
4. **Output-driven formulation:** Can a receiver response or its statistics be obtained without resolving unnecessary spatial outputs? Does the useful reduction occur in the propagation, surface response, receiver functional, frequency dependence, or stochastic variables?
5. **Statistical sufficiency:** Which statistics and distributions matter for the selected communication benchmark? Are Gaussian and second-order descriptions adequate, including for fades or outage, or are additional structures needed?
6. **Total computational cost:** What is the cost of initial reference generation, model construction, independent validation, condition updates, and online sampling? Which quantities remain reusable when frequency, sea realization, depth, offset, or sea state changes?
7. **Communication acceptance:** Propose modulation, bandwidth/rate, noise normalization, receiver knowledge, synchronization, equalization, and tolerances. Separate errors in the communication implementation from errors in the channel model.
8. **Time variation:** How can the static formulation later support evolving surfaces and platform motion with consistent temporal correlation and Doppler? What extra state and validation would be required?

No measured-data availability, real-time deadline, memory ceiling, or fixed numerical communication tolerance has been supplied. Absence of a supplied requirement should be recorded as an open decision, not interpreted as evidence that the constraint or dataset does not exist.

The external AI reply also left several pitfalls that future work must avoid: uniform-medium wave propagation is not an identity operator; multiplying by a spatial phase screen can redistribute angular energy; and independent azimuthal 2D slices are not automatically a valid solution for a general two-coordinate rough boundary. Physical short-wave truncation, numerical surface sampling, illuminated area, and converged integration support are distinct design choices. [R13]

## 8. Requested output from a new AI researcher

The new assessment should produce a research design rather than immediately modifying the existing code.

- Restate the application and identify the smallest set of assumptions needed to make the first study well posed.
- Compare a small, justified set of candidate physical/statistical routes. For each, explain applicability, three-dimensional scaling, reference requirements, offline cost, update cost, online cost, and main failure modes.
- Distinguish reusable evidence from reusable code: an implementation may be useful even where its old physical conclusions do not transfer.
- Propose a layered validation strategy covering physical response, numerical convergence, statistical fidelity, and communication performance.
- Select the first small experiment that would best discriminate between the leading routes. State what outcomes would support, reject, or redirect each route.
- Provide a staged plan from the static three-dimensional single link to time-varying operation, with decision points and stopping criteria.
- Use verifiable primary literature for any new scientific claims. Do not treat the local project reports as an exhaustive literature review or assume untested methods are already implemented.

### 8.1 Suggested sequence after the external reassessment

This is a proposed research sequence, not a user-mandated solver architecture.

1. **Definitions and analytic controls.** Fix source normalization, spectral conventions, phase/time references, and statistical conditioning. Establish the 3D flat pressure-release image solution before any rough-surface integration. For horizontal displacement rho, the direct and reflected distances are \(L_d=\sqrt{\rho^2+(z_t-z_r)^2}\) and \(L_r=\sqrt{\rho^2+(z_t+z_r)^2}\). Thus a vertical Tx1000/Rx5 control at c=1500 m/s has reflected excess delay 6.6667 ms, not zero or a few microseconds. Its ideal two-path interference period is 150 Hz. These are lossless point-source controls, not measured ocean results.
2. **Small 3D reference and cost pilot.** Test flat, single-direction, perpendicular-direction, and mixed-direction surfaces. Separate quadrature/aperture convergence, missing transverse coupling, and physical boundary-model error. A direct 3D Kirchhoff calculation can check another implementation of that approximation; independent reference evidence is still needed for its physical accuracy. Measure cost before selecting ensemble size. Integrating the production PE with FK is not a prerequisite if another formulation answers the target more directly.
3. **Mechanism and statistical candidate comparison.** Compare a direct physical baseline, a random few-path model with a residual, and conditional statistical candidates. Use reflected power, arrival structure, joint-frequency dependence, pseudo-covariance, fade/outage distributions, and communication sensitivity. Coherent fraction alone must not select the distribution or scattering model.
4. **Generator and communication validation.** Use held-out surfaces and operating conditions, account for initial construction and updates as well as generation, and establish a trustworthy communication baseline. Separate generator fidelity from receiver-model mismatch. For the latter, use the same held-out physical channels and pilot observations, with a perfect-CSI control and a practical estimation experiment.
5. **Time evolution.** Extend using consistent evolving surface modes and platform motion, then validate correlation and Doppler. Static independent samples do not supply a temporal model. Neither WSSUS nor an update interval follows from an assumed separation of time scales alone.

Decision branches should be explicit. If flat or deterministic controls fail, correct the implementation before collecting stochastic data. If a reduced model disagrees with a converged 3D comparator, restrict or abandon the reduction. If a few-path candidate reproduces held-out joint statistics and communication results, prefer its lower cost; otherwise test a richer residual or another family. If Kirchhoff fails an independent physical check, changing its rank cannot repair the failure. The first pilot should resolve these choices before a large parameter sweep.

### 8.2 Published scenarios versus proposed research settings

**Published evidence.** The Nereus trials used a downward-facing ship transducer lowered approximately 2–3 m below the surface, with acoustic baffling, and directional hardware operating over 7–14 kHz. The campaign reached 10,903 m vehicle depth. This provides a concrete deep-to-suspended-surface application analogue, but does not define this project's offset or sea-state envelope. [Singh et al., OCEANS 2009](https://www.researchgate.net/profile/John-Bailey-4/publication/224119665_Acoustic_communication_performance_of_the_WHOI_Micro-Modem_in_sea_trials_of_the_Nereus_vehicle_to_11000_m_depth/links/02e7e53a0cde1c4f40000000/Acoustic-communication-performance-of-the-WHOI-Micro-Modem-in-sea-trials-of-the-Nereus-vehicle-to-11-000-m-depth.pdf).

WMO's descriptive sea-state guidance lists slight, moderate, and rough height ranges of 0.5–1.25, 1.25–2.5, and 2.5–4 m. These are descriptive guidance, not communication-system operating limits or a PM parameter conversion. [WMO Table 3700 via NOAA/NCEI](https://www.nodc.noaa.gov/gtspp/document/codetbls/wmocodes/table3700.html).

**Planning proposal, not confirmed requirements.** Consider Tx depths 1500 and 3000 m initially, with 6000 m later; Rx 5 m initially, followed by 3 and 10 m. Let \(q=\rho/(z_t+z_r)\). Proposed q values 0, 0.05, 0.10, and 0.25 span normal angles 0–14.04 degrees; q=0.50 is an optional 26.57-degree stress case. At Tx1500/Rx5, q=0.10 means rho=150.5 m. These offsets are analyst-selected angular coverage, not published vessel-offset statistics.

Use a flat control and provisional spectral wave heights \(H_{m0}=4\sqrt{\mathbb E\eta^2}\) of 0.5, 1, and 2 m; reserve 3 m for a later stress case. Keep 4–8 kHz. Start with Tx1500/Rx5, q=0/0.10, and Hm0=0.5 m after deterministic checks, rather than running all combinations. Peak period, spreading, cutoff, and sampled moments still need specification; independent height/period choices are not automatically a canonical PM sea. The initial uniform-medium, ideal-source, hull-free control must remain labeled as such. The companion reassessment gives the evidence trail and detailed interpretation. [R13]

## 9. Source register and audit trail

The links below are repository-relative for local navigation. The surrounding document contains the information needed by a reader without repository access. Source IDs identify the evidence behind the claims above; these are project artifacts, not external publications.

| ID | Source | What it supports |
|---|---|---|
| R1 | [Current technical guide](../vertical_comm_guide.md), especially sections 4–9 | Existing local surface models, receiver adjoint, statistical generators, phase/CIR semantics, and validation boundaries |
| R2 | [100 m Bellhop rebaseline plan](bellhop_100m_rebaseline_plan.md), current revision at the beginning | Accepted Tx99/Rx3 shallow test installation; older sections of that file are historical |
| R3 | [Full Kirchhoff 2-D PE–BIE validation](pe_bie_full_kirchhoff_validation_report.md) | Controlled local-model versus FullFK comparisons and fixed convention handling |
| R4 | [Direct FullFK implementation](../scripts/validation/support/solve_full_kirchhoff_2d.m) | Actual curve integral, normal orientation, native time convention, and receiver-specific evaluation |
| R5 | [FullFK PM validity boundary](full_kirchhoff_pm_validity_boundary_report.md), numerical uncertainty and current validity domain sections | Conditional cutoff/amplitude findings and numerical guards |
| R6 | [Spectral rank audit](full_kirchhoff_spectral_rank_audit_report.md) | Incident height/spectrum separation, weighting, and scope of low-rank evidence |
| R7 | [FastFK prototype report](spectral_fast_kirchhoff_prototype_report.md) | Deterministic 3/6 configuration and dense-line timing |
| R8 | [FastFK preparation implementation](../scripts/validation/support/prepare_spectral_fast_kirchhoff_2d.m) | Distinct incoming and outgoing bases and lag-kernel construction |
| R9 | [Fixed PM pilot](spectral_fast_kirchhoff_pm_pilot_report.md) | Canonical surface definition, embedding, physical versus numerical cutoff, and algorithm/physics separation |
| R10 | [Multifrequency consistency report](spectral_fast_kirchhoff_multifrequency_consistency_report.md) | Uniform 3/12 rank, frequency/phase errors, BIE anchors, native-to-DSP bridge, and diagnostic CIR limits |
| R11 | [Single-receiver timing report](fast_fk_single_receiver_timing_report.md) and [test entrypoint](../scripts/validation/validate_fast_fk_single_receiver_timing.m) | New single-point timings, repetition protocol, exact inputs, exclusions, and complex outputs |
| R12 | [Two-node communication validation](two_node_statistical_channel_communication_validation_report.md) | Earlier QPSK comparison and reported processing-related BER floor |
| R13 | [External-reply reassessment](deep_sea_channel_external_reply_reassessment.md) | Corrections to the supplied AI proposal, current deployment clarifications, primary scenario sources, and provisional staged plan; no new simulations |

The requirements in Section 1 originate from direct user clarification on 2026-10-09, not from old simulation defaults. Historical reports can contain then-current recommendations later superseded by subsequent evidence. This handoff records clarified intent without changing the production implementation, its interfaces, or its certified scope. No contents of the retired `cash/` directory were used.

## 10. Ready-to-send request

Copy this document together with the following request into the new AI conversation:

> Please independently reassess the research direction described above. My goal is a physically credible statistical acoustic channel model and fast channel generator for a near-vertical link from a transmitter deeper than 1000 m to a receiver suspended 3–10 m below a ship or another surface platform, with horizontal displacement allowed. The first stage must address three-dimensional propagation and a two-dimensional static rough surface in a uniform-sound-speed medium. Later work should support evolving sea surfaces, temporal channel correlation, and Doppler. The preferred initial band is 4–8 kHz, but you may justify another representative band. PM is a reasonable starting sea model, not a permanent restriction.
>
> I allow limited offline physical simulation or training, provided you account for initial construction, validation, changes in operating conditions, and online generation costs. Communication-performance tolerances should guide the eventual acceptance criteria; no particular communication system or numerical tolerance has yet been fixed. Please propose a fair benchmark instead of inheriting the old QPSK implementation without review.
>
> Treat the supplied project results as conditional evidence. In particular, the current FastFK is accurate and fast for a dense receiver line, but the new single-point test found direct FullFK faster. Do not assume that further FastFK optimization, retaining PE, retaining FullFK, or reverting to a local phase screen is the correct next step. You may reconsider the entire physical and statistical formulation.
>
> The receiver is cable-deployed; its hull clearance and the source transducer model remain unknown. There is no engineering maximum offset or sea-state envelope yet. Published deployments and provisional test values above provide starting points, not additional requirements. Preserve the distinction between an ideal free-surface study and verified shipboard operation. Specify statistical conditioning and distinguish a random reflected path from its ensemble mean. Do not choose a distribution from coherent fraction alone or assume independent 2D slices reproduce a 3D rough boundary. Recheck numerical geometry, spectral moments, and cost estimates before using them to select a route.
>
> First identify the physical assumptions and unresolved requirements that control the problem. Then compare a small set of plausible routes, justify them with verifiable primary sources, and distinguish numerical consistency from physical and statistical validity. Explain which existing results and tools can be reused, which do not transfer to deep three-dimensional conditions, and which claims require new evidence. Propose the first small experiment that would most effectively distinguish the leading routes, followed by a staged research plan with acceptance criteria, stopping conditions, and a path to time variation. Where information is missing, state the uncertainty and ask focused questions rather than treating historical test settings as requirements.
