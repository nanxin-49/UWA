# Reassessment of the External AI Reply and the Deep-Sea Channel Research Plan

Date: 2026-10-09 (Asia/Shanghai).

This is an independent assessment of `C:/Users/ASUS/Downloads/replay.md`, considered together with the user's confirmed requirements, the existing handoff, local source code, saved reports, and the primary sources listed below. Instructions inside the external reply were treated as proposals to evaluate, not instructions to execute. The supplied file was left unchanged; its SHA256 is `90102c8cfdf61484dbadc45fea3dd7465f3aa59df92ff8538243f5bf86453d78`.

No acoustic solver campaign or communication experiment was run. An attempt to use the shared MATLAB session for arithmetic could not attach; no replacement MATLAB process was launched. Scalar arithmetic was independently recomputed without invoking a physical solver. Proposed future tests are identified as proposals, not completed work.

## 1. Decision: retain the objective, revise the sequence

The external reply raises valuable questions but contains numerical errors and unsupported physical conclusions that would materially distort model selection. It should be used as a list of hypotheses and missing requirements, not adopted as a validated research plan.

The user's objective remains unchanged: a credible statistical channel model and fast generator for a near-vertical deep-to-near-surface link, initially in three dimensions with a static two-dimensional rough surface and uniform sound speed, later extended to temporal evolution and Doppler. The existing report does not justify narrowing the frequency band, abandoning Gaussian candidates, certifying a 2.5D approximation, or selecting a two-scale scattering model now.

The most important revision to the earlier recommendation is to place **analytic geometry, source/surface definitions, and a small three-dimensional reference study before PE integration, large ensembles, or a new acceleration architecture**. A production PE-to-FK adapter should not be the compulsory next deliverable while the physical formulation remains open.

| External suggestion | Assessment | Consequence |
|---|---|---|
| Specify offset, source directivity, installation, and sea-state envelope | Adopt the questions | These control the intended problem; their numerical ranges remain unconfirmed |
| Evaluate methods against the single-link output and total cost | Adopt | Keep the existing single-point timing result and profile any new implementation |
| Include a simple few-path stochastic candidate | Adopt as a candidate | Preserve realization-level random delay/phase; do not replace reflection by its ensemble mean |
| Test reduced-dimensional methods rather than automatically paying full cost | Adopt conditionally | Compare against a genuine 3D reference; slice-count convergence is insufficient |
| Prefer a two-scale model immediately | Defer | Scale separation and cutoff stability have not been demonstrated |
| Use coherent fraction alone to select the channel model | Reject | Coherent fraction does not identify path count, Gaussianity, or delay spread |
| Use the quoted coherent-reflection table as a design constant | Reject | The arithmetic is inconsistent with its own inputs |
| Start 144,000 BIE solves as a cheap pilot | Reject | Its assumed solver cost contradicts the saved implementation timings |
| Diagnose the existing noise-free communication error floor | Adopt | Build a trustworthy benchmark; do not assume the repair is trivial |

## 2. Corrections that change the physical interpretation

### 2.1 A vertical link still has a surface-reflected path

For a point source and receiver below an infinite flat pressure-release surface in a uniform, lossless medium, the image construction gives

\[
L_d=\sqrt{\rho^2+(z_t-z_r)^2},\qquad
L_r=\sqrt{\rho^2+(z_t+z_r)^2},\qquad
\Delta\tau=(L_r-L_d)/c.
\]

With native time dependence \(\exp(-i\omega t)\), a normalized point-source baseline is

\[
H_{\rm flat}(f)={e^{ikL_d}\over4\pi L_d}-{e^{ikL_r}\over4\pi L_r}.
\]

This is an analytic control under the stated assumptions, not a model of the ship, seabed, or rough ocean. Directional sources require consistent direction-dependent weights.

For the illustrative values \(z_t=1000\) m, \(z_r=5\) m, and \(c=1500\) m/s:

| Horizontal offset rho | Reflection angle from the mean normal | Reflected-minus-direct delay | Two-path interference period, 1/Delta-tau |
|---:|---:|---:|---:|
| 0 m | 0 degrees | 6.666667 ms | 150.000 Hz |
| 100 m | 5.682 degrees | 6.633580 ms | 150.748 Hz |
| 500 m | 26.451 degrees | 5.962836 ms | 167.705 Hz |
| 2000 m | 63.320 degrees | 2.981418 ms | 335.411 Hz |

The external reply's microsecond-scale direct/specular separations do not agree with this geometry. Its statement that the vertical case has only a direct path is false when surface reflection is included. For this flat geometry, increasing offset decreases the direct/reflected path separation; that result does not predict the rough-surface delay spread.

The interference period is not a universal coherence-bandwidth formula. Null depth depends on relative amplitudes and phases. Nevertheless, the flat example already disproves the claim that small offset necessarily produces a frequency-flat channel across 4 kHz.

Use \(\theta_n=\arctan[\rho/(z_t+z_r)]\) for angle from the surface normal and \(\theta_g=\pi/2-\theta_n\) for grazing angle. The external reply calls an angle that is zero at vertical incidence a grazing angle. Its slope criterion involving \(\tan\theta_0\) cannot be transferred without resolving this convention.

### 2.2 The coherent-reflection numbers do not pass arithmetic checks

At 4 kHz and 1500 m/s, \(k=2\pi f/c=16.755161\) rad/m. Taking the reply's own illustrative RMS height of 0.35 m gives \(k\sigma=5.864306\), not 0.294. At 8 kHz it gives 11.728613, not 0.587.

Inserting those values into the reply's normal-incidence factor \(\exp[-2(k\sigma)^2]\) gives approximately \(10^{-29.871}\) and \(10^{-119.483}\), respectively. These are arithmetic consistency checks using its assumed height, **not certified ocean reflection predictions**.

Likewise, its assumed 0.90 m at 8 kHz gives \(k\sigma=30.159289\) and a base-10 logarithm of the factor of approximately -790.054, not a factor of 0.011. The proposed design constant of 0.011, and the band/sea-state decisions derived from it, must be discarded.

The RMS heights are also inconsistent with the spectrum printed in the reply. For that literal one-sided scalar spectrum,

\[
S(K)=\alpha K^{-3}\exp[-\beta(K_p/K)^2],\qquad
\int_0^\infty S(K)\,dK={\alpha\over2\beta K_p^2}.
\]

Its own \(\alpha=0.0081\), \(\beta=1.25\), and \(K_p=g/U^2\) give RMS heights approximately 0.1451, 0.5802, and 1.3055 m at the stated wind values 5, 10, and 15 m/s, rather than its table. This does not establish that its wind convention or spectral normalization is the appropriate PM model.

The project instead implements a radial spectrum proportional to \(\alpha/(2K^3)\exp[-\beta g^2/(U^4K^2)]\), and a two-dimensional spectral density obtained by dividing by \(2\pi K\) in its isotropic construction. Frequency spectra, radial wavenumber spectra, directional densities, peak definitions, and wind measurement heights must not be mixed. The original PM paper explicitly discusses the anemometer-height issue. [P1; local surface implementation]

For the \(K^{-3}\) tail, the height-variance tail is finite and decays as \(K_{\max}^{-2}\); the slope-variance moment grows logarithmically with the upper cutoff; the curvature-variance moment grows quadratically. The reply confuses these different moments. An unlimited PM gravity-wave tail is not automatically a suitable differentiable boundary for a normal/curvature-based solver. Use a documented physical short-wave model or cutoff, separate it from numerical resolution, and test sensitivity. Keep the old 12-mode fixture as a regression case, not as the final ocean model.

### 2.3 Loss of the ensemble mean does not identify diffuse multipath

For a Gaussian height variable, the characteristic-function identity is

\[
\mathbb E[e^{ia\eta}]=e^{-a^2\sigma^2/2}.
\]

It is exact for that random variable. Turning it into a prediction of a complete bistatic, finite-aperture, near-surface receiver field requires additional physical assumptions. It is not an exact universal BIE acceptance criterion.

There is a simple counterexample to selecting a scattering model from coherent fraction alone. Let

\[
H_r(f)=A(f)e^{-i2\pi f\tau},
\]

where each realization contains one reflected path and only its delay is random. The ensemble mean may be very small because phases cancel, although each realization still has a strong single reflected arrival. A rigidly shifted flat reflector provides a useful conceptual control: it changes path length without creating roughness-induced angular scattering.

If the purpose is to characterize reflection, compute

\[
\gamma_r(f)={|\mathbb E H_r(f)|^2\over\mathbb E|H_r(f)|^2},
\]

and report the reflected power relative to the direct component separately. Applying the same ratio to total H can hide a random reflected field behind a deterministic direct path. Also declare what is conditioned on, the phase/time origin, and whether any realization-specific alignment is used. Unaligned and explicitly conditioned statistics answer different questions and must be reported separately.

A low coherent fraction neither proves non-Gaussianity nor requires nonzero pseudo-covariance. Conversely, a nonzero mean does not force nonzero centered pseudo-covariance: \(H=\mu+Z\) with proper complex Gaussian Z has mean mu and zero centered pseudo-covariance. The primary definition and affine invariance are given by Neeser and Massey. [P2]

Candidate distributions and generators must therefore be tested using joint-frequency statistics, realization-level delay structure, amplitude/fade distributions, and held-out communication performance. The thresholds 0.6, 0.4, and 0.05 proposed in the reply have no demonstrated decision-theoretic basis here.

### 2.4 Several solver claims are too strong or incorrect

**Uniform-medium propagation is not the identity.** In an angular-spectrum representation, nonzero transverse wavenumbers acquire a phase involving \(\sqrt{k^2-k_x^2-k_y^2}\). Removing a reference carrier does not remove this relative phase. The local WAPE code contains this dependence explicitly. Multiple homogeneous propagation steps may be collapsed under appropriate assumptions; that can make a direct homogeneous propagator attractive, but does not imply that PE has no diffraction or that a one-way method cannot contain transverse scattering.

**A phase screen can redistribute angular energy.** Spatial multiplication \(G(x,y)u(x,y)\) becomes a convolution in transverse wavenumber. Its known physical shortcomings should be described through the boundary/propagation approximation and validated error, not the false claim that it can only attenuate a coherent beam. The existing project explicitly constructs a random screen and incoherent statistics, not just its coherent mean.

**A general rough surface does not become independent 2D slices merely because c is uniform.** A boundary \(\eta(x,y)\) couples transverse directions. Exact invariant-coordinate reductions and approximate azimuthal slice sums are different methods. Increasing the number of independent slices can demonstrate convergence of the slice approximation without bounding its missing coupling. The proposed less-than-1-dB claim above 15 degrees is unsupported; any useful reduced model needs comparison with a genuine 3D reference on representative roughness and directionality.

**Illumination, Fresnel contribution, and integration support are different.** A vertical source at 1000 m with a 30-degree half-angle geometrically illuminates a radius near 577 m on a flat surface. This is not automatically the necessary numerical aperture. A much smaller region may dominate a coherent integral, but its size and truncation error require a calculation. At normal incidence the paraxial first-Fresnel radius is of order \(\sqrt{\lambda z_tz_r/(z_t+z_r)}\); that is also not a justified rough-surface truncation rule. At 6 kHz and a 5 m receiver-to-surface distance, R/lambda is 20, while kR is approximately 126. The latter is not a count of wavelengths. Far-field validity depends on the contributing aperture as well as the range.

### 2.5 Cost estimates and proposed experiments need rebuilding

The existing multifrequency report records BIE wall times of 150.72, 284.86, and 481.05 s for 2560, 3840, and 5120 nodes at its three anchors. The external reply assumes approximately 0.5 s near 4000 unknowns without measurements. C-language speed or parallel scaling claims are likewise unverified.

Its proposed grid contains \(3\times5\times3\times2\times200\times8=144000\) BIE solves. Twenty hours follows arithmetically from the unverified 0.5 s assumption. Applying the saved 150.72–481.05 s values mechanically would instead give approximately 251–802 serial days. Neither extrapolation is a validated budget for the new geometry; the discrepancy makes a measured pilot essential.

Conversely, the reply's statement that 50 realizations times 65 frequencies at the measured 2D FullFK cost implies about ten hours is arithmetically inconsistent: \(50\times65\times0.11155524\) s is about 362.55 s, or 6.04 min. This is only a constant-cost illustration, not a measured broadband or changing-surface runtime.

Inspection of the current FastFK apply code shows loops over outgoing rank and three derivative channels, vector operations, repeated FFTs, and lag-kernel embedding. It does not establish a per-receiver MATLAB loop bottleneck. Profiling is required before attributing its 36 ms to a specific mechanism.

### 2.6 Communication and time-variation claims also need qualification

- A pure small phase discrepancy of 0.001 rad at 6 kHz corresponds to about 39.8 micrometers of path length, or 26.5 ns, not 3 mm or 2 microseconds. A general complex relative error does not uniquely determine delay error, which depends on phase slope across frequency.
- Coherence bandwidth requires a stated correlation definition and threshold; it is not universally equal to the inverse RMS delay spread.
- A short delay spread relative to a coherence time does not establish wide-sense stationarity and uncorrelated scattering. Those are separate second-order assumptions. [P3]
- The quoted 2–6 Hz Doppler spread, 0.2–0.5 s coherence time, and 50–100 ms update interval have not been derived for the specified sea and installation. Independent static redraws do not preserve temporal correlation.
- A constant-velocity translation of one surface cannot generally reproduce all dispersive wave modes. For a time-dependent path, a leading-order Doppler relation is \(f_D\approx-(f/c)\,dL/dt\) under the corresponding convention and slow-motion assumptions; both geometry and surface/platform dynamics determine it.
- Sixty-five old frequency samples do not by themselves define a validated OFDM system. Occupied carriers, cyclic prefix, spectral sampling, tail convergence, pilots, synchronization, and noise normalization still need specification.
- Comparing a receiver using true CSI with a different receiver using estimated/model CSI mixes estimation loss with model mismatch. Evaluate model-induced receiver degradation on the same held-out physical channels and pilot observations, and separately evaluate whether generated channel ensembles reproduce communication results.

## 3. What to change in the handoff requirements

Keep all confirmed application requirements. Add the following distinctions rather than importing the reply's assumed numerical values:

1. Separate an **ideal free-surface single-link model** from an **installation-inclusive model**. The user has confirmed cable deployment, but the position and clearance relative to the hull remain unspecified. Cable deployment alone establishes neither universal shadowing nor negligible platform effects.
2. Parameterize geometry with transmitter/receiver depth, horizontal displacement vector, surface-normal incidence angle, and relative wind/wave direction. Offset is important but not a substitute for roughness, directivity, and receiver proximity.
3. Specify a three-dimensional source model and normalization before using the old 2D Gaussian source as an engineering surrogate. An ideal isotropic source may serve as an explicitly labeled analytic control, not as an assumed device specification.
4. Keep both physical truth and statistical conditioning explicit. Distinguish mean response, random coherent/specular path parameters, diffuse residual, and total/direct/reflected components.
5. Define the spectral measure, directional normalization, wind convention or wave-height/period parameterization, and physical short-wave limit before choosing a sea grid.
6. Treat a three-dimensional solution as a required physical capability, not necessarily a requirement to run a full 3D solver for every generated sample. Any reduction needs an independently demonstrated error boundary.

The user has confirmed that the source transducer model is unknown and that no engineering offset or sea-state envelope is available. The user requested typical scenarios from published evidence. The following research settings are provisional proposals, not newly confirmed requirements.

### 3.1 Published deployment examples and what they establish

The Nereus sea trial is a close application analogue: the ship transducer was lowered from the stern, faced downward, and operated approximately 2–3 m below the surface with acoustic baffling against ship noise. Directional transducers operated over 7–14 kHz; the vehicle reached 10,903 m. Reported link performance included approximately 5000 bit/s at 6000 m range and reliable 200 bit/s at 11,000 m. These are results for that equipment and deployment, not predicted performance for this project. [P5, Sections III–IV and abstract]

Webster's field account also describes acoustic monitoring and command updates during Nereus surveys around 5200 m depth in the Cayman Trough. It gives an additional operational example of a deep vehicle communicating with its support vessel. [P6, Section 3.4.4]

These examples support investigating several-kilometer deep-to-surface links and a suspended near-surface receiver. They do not establish a universal horizontal-offset limit, a required source beam pattern, or a sea-state operating envelope. The inspected Nereus communication paper does not supply the offset and sea-state distributions needed to infer those limits. Its directional hardware and actual ocean environment also do not validate the proposed uniform-medium, isotropic-source control. Keep 4–8 kHz as the user's preferred research band; 7–14 kHz is an evidence-based alternative to assess later, not an automatic replacement.

For a transparent description of sea severity, WMO Code Table 3700 gives indicative wave-height ranges of 0.5–1.25 m for slight, 1.25–2.5 m for moderate, and 2.5–4 m for rough seas. Its notes give priority to descriptive sea state and treat heights as guidance; they are not a modem operating specification or a unique PM wind-speed conversion. [P7]

### 3.2 Proposed initial scenario set, pending physical and cost checks

Define \(q=\rho/(z_t+z_r)\), so that the flat-reflection angle from the mean normal is \(\arctan q\). Define the modeled spectral wave height explicitly as \(H_{m0}=4\sqrt{m_0}\), where \(m_0=\mathbb E[\eta^2]\) for a zero-mean surface. This definition avoids silently equating a descriptive sea-state code with a particular spectral normalization.

| Parameter | Proposed starting values | Status and purpose |
|---|---|---|
| Tx depth | 1500 and 3000 m; 6000 m as a later extension | Analyst-selected shallow-deep, central, and deeper research cases; neither the maximum depth nor a field population distribution |
| Rx depth | 5 m initially; then 3 and 10 m | Central value and limits of the user-confirmed interval |
| Horizontal displacement | q = 0, 0.05, 0.10, 0.25; q = 0.50 only as an optional stress case | Normal angles 0, 2.86, 5.71, 14.04, and 26.57 degrees; proposed angular coverage, not measured vessel offsets |
| Modeled wave height | Flat control, then Hm0 = 0.5, 1, and 2 m; 3 m as a later stress case | Chosen to span mild through more demanding roughness; informed by P7's descriptive scale, not an approved operating limit |
| Acoustic band | 4–8 kHz; 4/6/8 kHz for early scalar diagnostics | User-preferred band; sparse diagnostics do not validate a wideband CIR |
| Source | Normalized isotropic point source for controls, then explicit directivity sensitivity | An ideal reference until hardware data exist; no assumed real transducer beamwidth |
| Sea spectrum | A fully specified PM candidate, with directionality and short-wave cutoff documented | Peak period and wave height must respect the chosen formulation; independent height/period tuning is not automatically canonical PM |

For example, at Tx 1500 m and Rx 5 m, q = 0.05/0.10/0.25 gives offsets 75.25/150.5/376.25 m. At Tx 3000 m the same angles give 150.25/300.5/751.25 m. Scaling displacement with depth makes the geometry comparison explicit. These numbers are proposed tests; the literature search has not established them as typical operational offset quantiles.

Start the rough-surface pilot at Tx 1500 m, Rx 5 m, q = 0 and 0.10, and Hm0 = 0.5 m after the flat and deterministic-surface checks. Vary roughness and displacement in a small crossed comparison at this fixed geometry before expanding depth, receiver depth, and directionality. Do not execute the full Cartesian product of the table. Measure cost and uncertainty before increasing samples. For every sea case, report peak period, directional spreading, cutoff, actual sampled height/slope moments, and receiver clearance from the realized surface; Hm0 alone does not define it.

Initially omit the hull only as an explicit idealized-domain assumption. The resulting conclusion would concern a suspended receiver in that idealized free-surface environment, not verified shipboard performance. Platform geometry, directivity, noise, and motion become a separate extension with their own evidence needs.

## 4. Revised research sequence

| Stage | Primary question and bounded work | Deliverable / decision |
|---|---|---|
| 0. Definition and analytic controls | Fix coordinate/phase conventions, source normalization, spectral measure, and conditioning; establish the flat pressure-release image solution | A specification and reproducible analytic controls; no preferred approximate solver yet |
| 1. Small 3D reference and cost pilot | Validate a candidate 3D implementation on flat and simple directional roughness; test integration support and resolution; measure one solve before defining an ensemble | Separate discretization error, dimensional-reduction error, and physical-model error; record actual cost |
| 2. Mechanism and model comparison | Compare a random few-path candidate, a statistical residual model, and a direct physical baseline on a bounded set of shared inputs; add a two-scale candidate only with evidence for its split | Select the simplest adequate model using multiple observables, not coherent fraction alone |
| 3. Statistical generator and communication validation | Use held-out surfaces/conditions; verify frequency dependence, fades, and delay statistics; compare generators and receiver mismatch under controlled experiments | Communication-relevant accuracy and complete offline/update/online cost accounting |
| 4. Time evolution | Introduce consistent modal dynamics and platform motion; test temporal correlations and Doppler against evolving physical responses | An explicit temporal model and validated update interval; no assumed WSSUS shortcut |

### The first experiment should be smaller than the proposed sweep

Begin with the analytic 3D flat-interface response and verify the native-to-DSP conversion. Include zero offset and a nonzero illustrative offset, explicitly labeled as controls rather than approved deployment limits. Check the direct and reflected fields separately, their sum, and the direct-referenced delay.

Then use a small set of deterministic surfaces with one horizontal orientation, the perpendicular orientation, and mixed directional content. This exposes missing transverse coupling more directly than many unvalidated PM realizations. A direct 3D Kirchhoff calculation can serve as a numerical comparator for another implementation of the same approximation; it does not become an exact physical reference by doing so. Physical validation needs an independent method or a justified limiting case.

Do not prescribe 200 realizations per cell in advance. Measure runtime and estimate statistical uncertainty in a pilot; grow the ensemble only for the observables and conditions that need it. Keep tuning cases separate from held-out validation cases. A 4/6/8 kHz diagnostic can precede a denser band, but three frequencies cannot certify a CIR or correlation bandwidth.

Recommended observables include reflected power, direct-to-reflected power ratio, reflected-only coherent fraction, joint-frequency covariance and pseudo-covariance, random arrival-delay structure, amplitude/fade statistics, and communication sensitivity. Do not select a Gaussian, compound, or few-path family solely from one scalar indicator.

### Candidate routes remain open

- **Direct physical sample generator:** useful as a baseline; online feasibility depends on measured three-dimensional cost and the required rate, neither currently fixed.
- **Random few-path model with residual:** test a model such as \(H=H_d+A_s(f,\xi)e^{-i2\pi f\tau_s(\xi)}+H_{\rm residual}(f,\xi)\), retaining joint random path parameters. It is not equivalent to replacing the path by its ensemble mean.
- **Conditional statistical generator:** retain both proper/improper Gaussian and non-Gaussian candidates until held-out distributions and communication tests decide between them.
- **Two-scale or reduced-dimensional physical model:** promising hypotheses where a controlled split or dimensional reduction can be demonstrated, not defaults justified by an unverified slope or angle threshold.

For communication, retain a perfect-CSI control and a practical estimated-CSI experiment. Do not require a 0.5 dB threshold, 65 OFDM carriers, LDPC, or a specific fading distribution merely because the external reply suggested them. Choose provisional tolerances transparently and check that the reference and Monte Carlo uncertainty are smaller than the claimed differences.

## 5. Remaining information and immediate recommendation

Receiver cable deployment is confirmed. The exact hull clearance, source transducer model/directivity, maximum horizontal offset, and target sea-state envelope remain unknown. Section 3 supplies published application examples and explicitly provisional research settings in response to the user's request; it does not convert those settings into requirements. Measured-data availability and offline/online budgets would further constrain the route.

The immediate research deliverable should be **a correct three-dimensional analytic baseline, explicit source/surface definitions, and a measured, limited reference pilot**. This corrects the earlier emphasis on making a PE-to-FK adapter the next mandatory step. Existing FullFK, FastFK, PE, and statistical code should remain available as tools and conditional evidence while the choice is made.

This reassessment changes documentation and planning only. It does not alter the production solver, declare a new physical validity domain, or authorize the external reply's large computation campaign.

## 6. Evidence and reference notes

### Local evidence inspected

- [Single-receiver timing report](fast_fk_single_receiver_timing_report.md): exact timing scope and pointwise fidelity.
- [Multifrequency consistency report](spectral_fast_kirchhoff_multifrequency_consistency_report.md): BIE node counts and measured wall times; diagnostic CIR limitations.
- [Current technical guide](../vertical_comm_guide.md): original statistical-chain definitions and scope.
- [FastFK apply implementation](../scripts/validation/support/solve_spectral_fast_kirchhoff_2d.m): actual loop and FFT organization.
- [PE/WAPE implementation](../src/propagation/vertical_wape_propagator.m): transverse spectral propagation factor.
- [Surface implementation](../src/surface/pm_surface_boundary_model.m): PM spectral conventions, spatial phase multiplier, and coherent/incoherent statistics.
- [Existing handoff, now clarified](deep_sea_channel_research_background_for_new_ai.md): user-confirmed goals and prior evidence register.

### Primary sources checked on 2026-10-09

- **P1. Pierson and Moskowitz (1964).** [A proposed spectral form for fully developed wind seas based on the similarity theory of S. A. Kitaigorodskii](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/JZ069i024p05181). Publisher abstract and bibliographic record inspected. Supports the need to specify fully developed-sea assumptions and wind-measurement conventions; the integrals above are independently derived from the displayed formulas, not quoted from this paper.
- **P2. Neeser and Massey (1993).** [Proper Complex Random Processes with Applications to Information Theory](https://www.isiweb.ee.ethz.ch/archive/massey_pub/pdf/BI437.pdf), IEEE Transactions on Information Theory 39(4), 1293–1302. Author/institution-hosted original paper inspected, especially the centered covariance definitions and affine invariance of properness. Used to correct the claim that a nonzero mean forces nonzero pseudo-covariance.
- **P3. Subramanian and Hajek (2002).** [Broad-Band Fading Channels: Signal Burstiness and Capacity](https://subramanian.engin.umich.edu/wp-content/uploads/sites/383/2018/04/Burstiness_Capacity.pdf), IEEE Transactions on Information Theory 48(4), 809–827. Author-hosted original paper inspected for the separate WSS and uncorrelated-scattering assumptions. It is not evidence that this acoustic link satisfies them.
- **P4. Voronovich.** [Wave Scattering from Rough Surfaces, second edition](https://link.springer.com/book/10.1007/978-3-642-59936-1). Publisher record confirms the second edition is copyright 1999, with the first edition in 1994. The reply's “1994, second edition” citation mixes editions. No uninspected book theorem is used here to certify a two-scale model.
- **P5. Singh et al. (2009).** [Acoustic Communication Performance of the WHOI Micro-Modem in Sea Trials of the Nereus Vehicle to 11,000 m Depth](https://www.researchgate.net/profile/John-Bailey-4/publication/224119665_Acoustic_communication_performance_of_the_WHOI_Micro-Modem_in_sea_trials_of_the_Nereus_vehicle_to_11000_m_depth/links/02e7e53a0cde1c4f40000000/Acoustic-communication-performance-of-the-WHOI-Micro-Modem-in-sea-trials-of-the-Nereus-vehicle-to-11-000-m-depth.pdf), OCEANS 2009, DOI 10.23919/OCEANS.2009.5422316. Original six-page paper in an author-uploaded copy; deployment and reported range/rate results inspected. This is a primary field report, not evidence for arbitrary sea states or the proposed offset grid.
- **P6. Webster (2010).** [Decentralized Single-Beacon Acoustic Navigation: Combined Communication and Navigation for Underwater Vehicles](https://iop.apl.washington.edu/~swebster/SarahWebster_Publications_files/SWebster-thesis-2010.pdf), Johns Hopkins dissertation, author/institution-hosted copy. Chapter 3 describes field deployments; Section 3.4.4 supplies the approximately 5200 m survey example.
- **P7. WMO Code Table 3700, hosted by NOAA/NCEI.** [Sea State](https://www.nodc.noaa.gov/gtspp/document/codetbls/wmocodes/table3700.html). Table and explanatory notes inspected. Used only as descriptive context for selecting research wave-height levels, not as an acoustic operating standard.

The cited Thorsos 1988 and Bello 1963 papers have identifiable bibliographic records, but their original publisher full texts were not accessible through the browsing tool in this session. They remain reading candidates; no specific quantitative validity bound is attributed to their uninspected full texts. The simple geometry, moment integrals, statistical counterexamples, and operation-count arithmetic above are explicitly shown so that the conclusions do not depend on those access limitations.
