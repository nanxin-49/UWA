# Spectral Fast Kirchhoff fixed PM pilot

Date: 2026-10-07. MATLAB 25.2.0.3312555 (R2025b) Update 6. Validation-only, 4 kHz, c=1500 m/s, Tx=100 m, Rx=3 m.

Algorithm pass: **1**; physical-reference pass: **0**; PM pilot PASS: **0**. No production model changes, Monte Carlo or frequency extension.

The unchanged discrete operator and FFT normalization are defined in [prototype method](spectral_fast_kirchhoff_prototype_method.md); its deterministic evidence is in [prototype report](spectral_fast_kirchhoff_prototype_report.md), with [prior separability audit](full_kirchhoff_spectral_rank_audit_report.md). Both exact carrier phases exp(-ik*eta) remain intact. This pilot changes coefficient-defined surface inputs only, retaining full normals, arc weights, incident physical derivatives and both nonlocal Green terms.

## Provenance and physical surface

Three fixed realizations use saved coefficients: canonical 260001 and independent Gaussian coefficient seeds 260002/260003, U=6 m/s, period 160 m. Saved coefficients are copied without renormalization. Canonical SHA256: `1f7eda465e4ae85b8ac038310edf053f2d062563a3b943bfd83bd59d07015f67`. All solvers consume one coefficient-defined surface; none resamples random numbers.

Canonical is the original fixed realization; the other two are independent Gaussian amplitude draws from its saved Sk band. This is three deterministic pilot inputs, not an ensemble estimate. Analytic Fourier derivatives determine slope, curvature and the full water-facing normal.

Raw PM eta=sum(a cos Kx+b sin Kx). Physical cutoff controls the underlying periodic Fourier series. A fixed C2 envelope is 1 for |x|<=42 m and 0 for |x|>=50 m, within [-60,60] m integration support. It preserves central PM heights exactly and makes the same finite surface compatible with the frozen BIE. Its edge transition is a documented geometric embedding, not a PM bandlimited claim or sponge. Raw and embedded samples/statistics are saved separately.

|case|actual Kmax|raw RMS height m|raw RMS/max slope|raw RMS/max curvature 1/m|embedded max height m|
|---|---:|---:|---|---|---:|
|seed_260001|0.4712389|0.1864104|0.05397697 / 0.1166758|0.02090655 / 0.04459474|0.4994562|
|seed_260002|0.4712389|0.1491987|0.04896024 / 0.0979314|0.01788279 / 0.03899057|0.3235663|
|seed_260003|0.4712389|0.1777665|0.04851077 / 0.1323573|0.01593412 / 0.04635763|0.430661|
|cutoff_low|0.2748894|0.1472515|0.0269472 / 0.06911712|0.005785351 / 0.01442481|0.3880059|
|cutoff_high|0.7461283|0.1966725|0.06746681 / 0.1739605|0.03409498 / 0.07336469|0.4889045|

## Physical versus numerical cutoff

Requested Kmax=0.30/0.50/0.75 rad/m; actual last modes 0.27488936 / 0.4712389 / 0.74612826. Main dx=120/4096=0.029296875 m; Nyquist=107.23303 rad/m, independent of Kmax. Grid refinements reevaluate the SAME coefficients; they add no physical roughness modes. Higher cutoff appends modes 13--19 with Gaussian coefficients, mt19937ar seed 26000175, Sk=alpha/(2K^3) exp(-beta*g^2/(U^4 K^2)), alpha=.0081 beta=.74 g=9.81. Every original coefficient is retained exactly; lower cutoff only removes upper modes. No height rescaling.

## Algorithm fidelity and uniform ranks

Frozen baseline 3/6 is evaluated FIRST on every case. Selected common incoming/outgoing ranks: **3/8**; uniform height bound A=0.5 m, fixed 33 height nodes. No per-case best rank, phase adjustment, scalar fit or incoming99 outgoing truncation.

Relative weighted complex L2 = sqrt(sum w|Fast-Full|^2 / sum w|Full|^2), fixed Stage-0 incident M99 weights plus Full-only -40 dB gate. Whole-grid uses all N receivers and ||Fast-Full||2/||Full||2. Magnitude and wrapped phase RMS use that fixed mask. Native exp(-iwt) fields are compared directly without conjugate/ratio calibration.

|case|study|r_in/r_out|weighted L2|whole-grid L2|magnitude L2|phase RMS rad|effective in/out rank 1e-3|
|---|---|---|---:|---:|---:|---:|---|
|seed_260001|baseline|3/6|0.0043814|0.004269419|0.003511881|0.0026285|3/7|
|seed_260002|baseline|3/6|0.004681543|0.004697074|0.004192559|0.002388715|3/7|
|seed_260003|baseline|3/6|0.004582056|0.004449605|0.003610972|0.002869323|3/7|
|cutoff_low|baseline|3/6|0.003839431|0.003491893|0.003082488|0.002145868|3/7|
|cutoff_high|baseline|3/6|0.004577301|0.004538892|0.003193818|0.003305464|3/7|
|seed_260001|selected|3/8|0.0002476376|0.0002966094|0.0001127168|0.0002443287|3/7|
|seed_260002|selected|3/8|0.0001926868|0.000272106|6.589413e-05|0.0002302067|3/7|
|seed_260003|selected|3/8|0.0001874777|0.0002594737|8.014542e-05|0.0002269272|3/7|
|cutoff_low|selected|3/8|0.0001966262|0.0002785341|7.172997e-05|0.0002244208|3/7|
|cutoff_high|selected|3/8|0.0002091716|0.0002746867|6.836048e-05|0.0002250821|3/7|

The already prescribed common 4/8 diagnostic also reaches the ideal 1e-4 across all five embedded cases: maximum weighted 9.354098e-05, whole-grid 9.640783e-05. It is reported as a diagnostic with limited margin; the primary selected 3/8 remains judged against 1e-3, without per-realization ranks. Full 45-row ladder, including failed combinations, is saved in rank_diagnostics.csv.


Effective ranks are relative Frobenius SVD tails of the weighted incoming matrix and stacked carrier-factored outgoing [G;Dz/k;Dx/k] lag-height table, not guaranteed field error bounds. The outgoing basis depends on height bound, receiver depth and lag support; K appears through spatial densities/normal, not directly in this basis. Separate rank diagnostics and incident/outgoing ablations are saved if baseline fails.

Canonical rank3/10 field error 0.0002874895, rank5/6 error 0.004280826. This separates incoming truncation from outgoing truncation without changing the physical model.

## Reflected aperture spectrum

Saved complex A_FK(kx), A_Fast(kx) use dx*FFT with an identical receiver Hann window. They are finite-aperture diagnostics, not infinite-domain far-field flux. Propagating |kx|<=k, near-evanescent k<|kx|<=1.2k, and higher bins are all retained. Fractional energies are spectral norms, not evanescent transported flux.

Additional unmodified PM controls (no C2 envelope) use all true Fourier heights and normals on [-60,60] m with the same selected ranks. These compare Fast-Full on the original finite PM segment, preventing the algorithm conclusion from relying on a favorable surface taper. They do not independently validate an infinite periodic PM boundary against BIE.

|case|untapered weighted L2|untapered whole-grid L2|
|---|---:|---:|
|seed_260001|0.000247638|0.0003108253|
|seed_260002|0.0001927136|0.000285742|
|seed_260003|0.0001874814|0.0002680339|
|cutoff_low|0.0001966497|0.0002888625|
|cutoff_high|0.0002091612|0.0002832788|

|case|theta95/99 out deg|outside incident99 prop energy|prop/near-evan energy|spectrum L2|near-evan relative L2|
|---|---|---:|---|---:|---:|
|seed_260001|17.453 / 20.865|0.005758169|1 / 1.606397e-09|0.0002502763|0.01015153|
|seed_260002|15.215 / 21.826|0.01147803|1 / 1.619756e-09|0.000205826|0.009514325|
|seed_260003|17.266 / 18.961|0.0008127435|1 / 1.603289e-09|0.0001997769|0.006785197|
|cutoff_low|13.92 / 15.772|0.0001170621|1 / 1.616964e-09|0.0002100902|0.004800335|
|cutoff_high|22.405 / 25.54|0.08684941|1 / 1.604063e-09|0.0002161988|0.01874237|

Near-evanescent relative errors can be large when the norm is nearly zero; absolute powers and total spectrum L2 are saved for interpretation. No spectral cutoff is applied to the Hankel-kernel solver.

## Grid, padding and interpolation

Canonical N=2049/4097/8193 use unchanged Fourier coefficients and analytic normals. 4097 versus8193: Full weighted 6.560588e-06, whole-grid 1.629997e-05; Fast weighted 6.562515e-06; nested surface closure 0 m. Padding2/4 L2 4.01e-16. Height-node33/65 L2 2.513433e-06. The latter is a separate interpolation/basis-quadrature check, not physical rank tuning.

|N|dx|Nyquist|matrix incoming/outgoing rank 1e-3|Fast-Full weighted|whole grid|
|---:|---:|---:|---|---:|---:|
|2049|0.05859375|53.61651|3/7|0.0002476414|0.0002966101|
|8193|0.01464844|214.4661|3/7|0.0002476384|0.000296609|

## Separate FullFK--BIE physical reference

Prespecified canonical and high-cutoff cases: frozen Dirichlet halfplane combined layer, h=-2 m, plateau ratio .85, order8/self96, ppw8/10, high-cutoff halfwidth60/70. FullFK reevaluated at the same 513 receivers; same coefficient surface and incident spectrum. Physical acceptable gate is weighted native field L2<1e-3 with numerical uncertainty checks below1e-4. This gate is separate from Fast-Full fidelity.

|case|L|ppw|Full-BIE weighted L2|whole-grid L2|boundary residual|BIE refinement L2|window L2|
|---|---:|---:|---:|---:|---:|---:|---:|
|seed_260001|60|8|0.0008000959|0.009965395|7.690677e-08|5.011037e-08|0|
|seed_260001|60|10|0.0008000959|0.009965395|4.050464e-08|5.011037e-08|0|
|cutoff_high|60|8|0.001197491|0.009990034|1.217007e-07|6.622239e-08|1.427075e-07|
|cutoff_high|60|10|0.001197491|0.009990034|5.660787e-08|6.622239e-08|1.427075e-07|
|cutoff_high|70|10|0.001197493|0.001442555|6.297332e-08|6.622239e-08|1.427075e-07|

High-cutoff Full quadrature4097/8193 at the same513 receivers: weighted convergence 6.557529e-06, whole-grid 1.606541e-05. Refined Full--BIE weighted 0.001198526, whole-grid 0.001449962. The remaining physical deviation is not resolved by Full quadrature or BIE refinement.

Whole-grid Full--BIE errors include finite-window edge receivers outside the validated M99 physical footprint. The high-cutoff L60/70 whole-grid error changes markedly, while the weighted core is stable. Physical accuracy is therefore not claimed on the entire receiver line or for an infinite periodic PM ocean; algorithm whole-grid fidelity is a separate result.


## CPU, memory and speedup

Warm Fast median5 applications; Full dense median2; cold includes Fast preparation or Full exact incident. Windows CPU resolution may give zero for short calls; wall times determine ratios. Memory is named arrays/estimated block payload, NOT peak RSS or FFT/SVD workspace. BIE dense matrix payload/times are in bie_reference_metrics.csv.

|case|Full apply s|Fast setup/apply s|warm/cold speedup|plan/setup/apply MiB|Full block MiB|
|---|---:|---|---|---|---:|
|seed_260001|2.49338|0.202411 / 0.022181|112.41 / 11.656|5.569 / 24.775 / 0.87521|104.03|
|seed_260002|2.3555|0.165448 / 0.0212255|110.98 / 13.224|5.569 / 24.775 / 0.87521|104.03|
|seed_260003|2.29867|0.152113 / 0.0212911|107.96 / 13.886|5.569 / 24.775 / 0.87521|104.03|
|cutoff_low|2.26519|0.171298 / 0.0253904|89.215 / 12.05|5.569 / 24.775 / 0.87521|104.03|
|cutoff_high|2.22808|0.160228 / 0.0204495|108.96 / 12.903|5.569 / 24.775 / 0.87521|104.03|

Fast workflow wall 51.14 s, CPU 358.03 s (excluding separate BIE/render stage). The same global FFT structure remains; finite pilot timings do not prove asymptotic complexity.

The five BIE reference solves total wall 943.708 s and process CPU 8453.98 s. Process CPU can exceed wall on multithreaded linear algebra; these expensive physical checks are not part of Fast setup/apply costs. Provenance hashes are in provenance.json.

## Decision and reproduction

**PM pilot does not meet all gates. Stop before multi-frequency/Monte Carlo.** Algorithm compression succeeds with common3/8 and stable grids/cutoffs. The high-cutoff physical FullFK--BIE discrepancy exceeds the prespecified1e-3 gate. Increasing FastFK rank cannot remove a FullFK physical approximation error. First establish the PM physical validity domain or compare an independent nonlocal SSA/SSA2 reference; block/local-stationary Fast separation is relevant only if compression subsequently fails. No production model is changed here.

This is a narrowly failed strict physical gate, not evidence that low-rank PM compression collapses or that every PM case invalidates FullFK. Canonical passes the weighted1e-3 physical gate; the high-cutoff discrepancy is predominantly phase (about1.2e-3 rad versus magnitude error1.39e-4). Its higher slope/curvature and wider outgoing spectrum suggest a model-validity boundary to investigate, without a fitted causal formula.

Five independent unit tests (PM derivatives/embedding/unmodified series and prior flat Green/FFT/interpolation controls) pass; saved artifact assertions verify coefficient preservation, Gaussian-seed reproduction, geometry, fixed ranks, spectrum norm closure and BIE/Full refinement data. Test_results.csv and artifact_assertions.mat preserve the verification.

Run after setup_vertical_project: `validate_spectral_fast_kirchhoff_pm_pilot('full')`; checkpoint stages `'fast'`, `'bie'`, `'raw'`, `'quadrature'`; `'render'` reuses saved numerical data. Outputs: results/validation/spectral_fast_kirchhoff_pm_pilot. No PE/WAPE, Full FK or BIE reference change; no communication regressions are relevant to these validation-only artifacts.

![pm_profiles](../results/validation/spectral_fast_kirchhoff_pm_pilot/figures/pm_profiles.png)

![ranks_fields_cutoffs](../results/validation/spectral_fast_kirchhoff_pm_pilot/figures/ranks_fields_cutoffs.png)

![reflected_spectra](../results/validation/spectral_fast_kirchhoff_pm_pilot/figures/reflected_spectra.png)

![cutoff_outgoing_bandwidth](../results/validation/spectral_fast_kirchhoff_pm_pilot/figures/cutoff_outgoing_bandwidth.png)

