# Full Kirchhoff PM validity boundary

Date: 2026-10-07. Numerical campaign complete: **1**. This is a controlled single-realization physics audit, not Monte Carlo or a production model.

**Finding:** At the fixed 4 kHz source, geometry and receiver footprint, the guarded master-prefix sequence passes at Kmax 0.70685835 rad/m (E 0.00085275345, U 6.82e-06) and fails at 0.74612826 rad/m (E 0.0011974915, U 6.82e-06). High-K amplitude and illuminated curvature/spectral structure matter; Kmax or height alone does not define a general validity threshold. The discrepancy is mainly native phase variation.

**Decision:** Retain the FullFK/FastFK route within the measured conditional physical domain. If the required spectrum includes rejected structures, next compare SSA/SSA2 on a limited physical-reference set. This audit tests FullFK against BIE on 15 controlled inputs; the separate FastFK fidelity evidence remains the prior pilot evidence. All Fast ranks/helpers and production physics stay frozen.

## Frozen mathematical and numerical definitions

4 kHz, c=1500 m/s, Tx=100 m, Rx=3 m, identical saved discrete Gaussian source, native exp(-iwt), water-facing normal and arc-length weights. [Prior PM pilot](spectral_fast_kirchhoff_pm_pilot_report.md) supplies exact canonical/high-cutoff Fourier coefficients and BIE anchors. FullFK, BIE, C2 embedding and FastFK helpers are unchanged. No FastFK solve or rank optimization is run; ranks3/8 and padding2 remain frozen.

PM eta(x)=sum(a cosKx+b sinKx); slope and second derivative are analytic. Geometric curvature kappa=eta_second/(1+eta_prime^2)^(3/2). C2 embedding equals1 on |x|<=42 m and0 beyond50 m; Full integrates the same [-60,60] m support. Raw periodic statistics use16384 unique points over160 m. Illuminated statistics use the frozen Stage-0 incident weights inside the canonical M99 footprint. Height/slope/curvature of the embedded surface are separately retained.

Frozen incident spectral theta95/theta99: 15.96928/21.28005 deg (coefficient-energy quantiles). Fixed receiver footprint extends to |x|=39.84375 m; all reported maxima are absolute magnitudes. Outgoing comparisons use aperture spectral energy, not an inferred infinite-ocean flux.

Numerical Full grid4097, dx=.029296875 m, Nyquist107.233 rad/m. Physical Kmax is the largest NONZERO Fourier mode, independently of numerical Nyquist; embedding itself is not bandlimited. No whole-surface renormalization, smoothing, phase alignment or complex scalar fitting.

## Design and coefficient provenance

Canonical seed260001, U=6 m/s, period160 m. Canonical SHA256 `1f7eda465e4ae85b8ac038310edf053f2d062563a3b943bfd83bd59d07015f67`. Master modes1--12 are canonical; modes13--19 are the previously saved high-cutoff pilot extension (additional-mode seed26000175). Every nested cutoff uses an exact master prefix, so every common low-K coefficient is identical. Adaptive midpoints are chosen only around the prespecified1e-3 crossing, capped at two extra cutoffs. Coefficients/design.csv are saved before the corresponding BIE response.

Case A: low added band13--15 and high band17--19 each have added periodic height RMS=.03 m. Only that new band is scaled; total periodic RMS height is matched exactly and remains close to canonical. Extra slope-identification pair uses single added modes13 and19: positive amplitude solves the input-moment quadratic to make illuminated total RMS slope^2=canonical+.02^2 exactly. Case B varies only modes13--19 amplitudes .25/.50/.75/1 at the same positive-amplitude Kmax, retaining each modal phase and all canonical coefficients. At amplitude0, the active Kmax drops to canonical by definition. These are controlled Fourier diagnostics, not additional random PM draws.

The slope pair matches the illuminated metric, not every local slope or every other moment. Record residual height/curvature differences before drawing a causal conclusion. Receiver samples are repeated measurements of one realization; they are not independent statistical replicates. No p-values, ensemble probabilities or empirical fitted formulas are reported.

In controlled coefficient CSVs, Sk_m3 is retained as provenance of the original PM sampling law; it is not a recalibrated PSD for the artificial amplitude modifications. All geometry, fields and statistics use the explicitly stored cosine/sine coefficients. The numeric MAT profiles are authoritative for scaled floating-point coefficients.

## Error norm and exact decomposition

All cases use one frozen canonical M99/Full-only-40dB mask and one fixed incident-energy weight vector on513 identical receivers. E_complex^2=sum(w|FK-BIE|^2)/sum(w|BIE|^2); E_mag^2=sum(w(|FK|-|BIE|)^2)/den. Raw phase RMS is sqrt(sum(w arg(FK conjBIE)^2)). The phase contribution obeys the exact identity E_complex^2=E_mag^2+sum(2w|FK||BIE|(1-cos(delta_phi)))/den. Raw phase mean/std are diagnostics only; neither is subtracted from a field. Whole-line unweighted error is separately recorded and is window-sensitive at the edges.

## Cutoff continuation and surface descriptors

|case|Kmax rad/m|k sigma_eta|RMS/max slope|RMS/max curvature 1/m|Kmax/k|max curvature/k|illuminated RMS curvature/k|
|---|---:|---:|---|---|---:|---:|---:|
|cutoff_n07|0.2748894|2.46722|0.0269472 / 0.0691171|0.00578548 / 0.0144248|0.0164063|0.000860918|0.000485227|
|cutoff_n10|0.3926991|2.8178|0.0385738 / 0.0914414|0.0111757 / 0.0256594|0.0234375|0.00153143|0.000643287|
|cutoff_n12|0.4712389|3.12334|0.053977 / 0.116676|0.0209072 / 0.0445947|0.028125|0.00266155|0.00147826|
|cutoff_n14|0.5497787|3.17382|0.0567532 / 0.146221|0.0227913 / 0.0555545|0.0328125|0.00331567|0.00142664|
|cutoff_n16|0.6283185|3.22233|0.0603639 / 0.157759|0.026053 / 0.0624603|0.0375|0.00372783|0.00171326|
|cutoff_n19|0.7461283|3.29528|0.0674668 / 0.17396|0.034096 / 0.0733647|0.0445313|0.00437863|0.00227302|
|height_match_lowband|0.5890486|3.16352|0.056315 / 0.13093|0.0226061 / 0.0517836|0.0351562|0.00309061|0.00143696|
|height_match_highband|0.7461283|3.16352|0.0582756 / 0.12765|0.02636 / 0.0554118|0.0445313|0.00330715|0.00178611|
|slope_match_lowband|0.5105088|3.19973|0.0579823 / 0.14302|0.0234921 / 0.054774|0.0304687|0.00326908|0.0015933|
|slope_match_highband|0.7461283|3.15533|0.0575486 / 0.136438|0.0256285 / 0.0580026|0.0445313|0.00346178|0.00172315|
|amplitude_025|0.7461283|3.13436|0.0549173 / 0.130452|0.0219683 / 0.0475607|0.0445313|0.00283857|0.00150885|
|amplitude_050|0.7461283|3.1672|0.0576461 / 0.144739|0.0248784 / 0.0550072|0.0445313|0.003283|0.0016748|
|amplitude_075|0.7461283|3.22118|0.0619276 / 0.159278|0.0290823 / 0.0641724|0.0445313|0.00383001|0.00194142|
|cutoff_n17|0.6675884|3.2312|0.0611119 / 0.156518|0.026815 / 0.0670861|0.0398438|0.0040039|0.00166386|
|cutoff_n18|0.7068583|3.24046|0.0619779 / 0.153524|0.0277901 / 0.0727306|0.0421875|0.00434079|0.00171844|

## FullFK--BIE field errors and phase structure

|case|complex L2|magnitude L2|phase RMS rad|phase mean/std rad|phase share of complex error energy|U_numeric|classification|
|---|---:|---:|---:|---|---:|---:|---|
|cutoff_n07|0.00032552899|0.0001330339|0.0003067182|-2.98239e-05 / 0.0003052648|0.83299|NaN|baseline_only|
|cutoff_n10|0.00038511836|0.0001332491|0.0003788092|-3.955889e-05 / 0.000376738|0.88029|NaN|baseline_only|
|cutoff_n12|0.00080009587|0.0001364041|0.0007857323|-0.0001277542 / 0.0007752768|0.97093|NaN|baseline_only|
|cutoff_n14|0.0007606848|0.0001414986|0.0007612633|-0.0001201972 / 0.0007517143|0.9654|NaN|baseline_only|
|cutoff_n16|0.0008354997|0.0001255908|0.0008975803|-0.0001748409 / 0.0008803869|0.9774|NaN|baseline_only|
|cutoff_n19|0.0011974915|0.0001388374|0.001198339|-0.0002879767 / 0.001163222|0.98656|6.815e-06|certified_above|
|height_match_lowband|0.00074140117|0.0001368149|0.0007636849|-0.0001237835 / 0.0007535863|0.96595|NaN|baseline_only|
|height_match_highband|0.00093922807|0.0001289718|0.0009401138|-0.0001844495 / 0.0009218418|0.98114|NaN|baseline_only|
|slope_match_lowband|0.0008720317|0.0001288246|0.000851788|-0.0001430939 / 0.0008396826|0.97818|NaN|baseline_only|
|slope_match_highband|0.00091551124|0.0001291275|0.0009081702|-0.000171776 / 0.0008917769|0.98011|6.797e-06|certified_below|
|amplitude_025|0.00079622603|0.0001319238|0.0007997562|-0.0001350772 / 0.0007882666|0.97255|NaN|baseline_only|
|amplitude_050|0.00087424652|0.0001298457|0.0008852453|-0.0001641559 / 0.000869892|0.97794|NaN|baseline_only|
|amplitude_075|0.0010153809|0.0001363355|0.001023987|-0.0002149764 / 0.001001167|0.98197|6.828e-06|certified_above|
|cutoff_n17|0.00080511331|0.0001273368|0.0008796801|-0.0001641597 / 0.0008642272|0.97499|NaN|baseline_only|
|cutoff_n18|0.00085275345|0.0001249502|0.0009031378|-0.0001748963 / 0.0008860413|0.97853|6.818e-06|certified_below|

## Outgoing spectrum and localization of discrepancy

Outgoing theta95/99 and outside-incident99 energy use the original full4097 receiver Hann-aperture FullFK spectrum. FK/BIE spectral differences use one fixed C2 receiver aperture (plateau25 m, support35 m), excluding finite-window edge errors. The BIE anchor data have513 receivers with Nyquist|kx|/k=.8; these core spectra are explicitly limited to that sampled range, not infinite-ocean far-field flux. Full energy above this sampling Nyquist is reported; no physical outgoing cutoff is imposed. Both complex spectra are saved without per-case amplitude/phase normalization.

|case|theta95/99 out deg|outside incident99 reflected energy|core spectrum L2|delta energy center/middle/tail|core phase mean/std rad|Full energy beyond BIE sampling Nyquist|
|---|---|---:|---:|---|---|---:|
|cutoff_n07|13.92 / 15.7719|0.000117062|0.0003250818|0.60702 / 0 / 0.39298|-2.16441e-05 / 0.000211209|4.185e-08|
|cutoff_n10|16.3307 / 18.017|0.000743294|0.0003918111|0.72968 / 0 / 0.27032|-1.6625e-05 / 0.000228476|4.178e-08|
|cutoff_n12|17.4532 / 20.8647|0.00575817|0.0007895641|0.92864 / 0 / 0.07136|-1.75517e-05 / 0.000460286|4.152e-08|
|cutoff_n14|17.4532 / 19.7196|0.00348349|0.0007562428|0.92508 / 0 / 0.074917|-1.46884e-05 / 0.00045978|4.134e-08|
|cutoff_n16|20.2911 / 23.3769|0.0332063|0.0009137754|0.91138 / 0.020514 / 0.068103|-1.61813e-05 / 0.000490744|4.112e-08|
|cutoff_n19|22.4054 / 25.5402|0.0868494|0.001196616|0.91721 / 0.036534 / 0.046258|-1.86239e-05 / 0.000810322|4.153e-08|
|height_match_lowband|18.017 / 21.0564|0.00634845|0.0007607044|0.91993 / 0 / 0.080071|-1.53582e-05 / 0.00045944|4.129e-08|
|height_match_highband|19.1502 / 22.0187|0.0135055|0.0009398507|0.92845 / 0.0085425 / 0.063008|-1.98772e-05 / 0.000604281|4.131e-08|
|slope_match_lowband|18.3939 / 20.4821|0.00366801|0.0008447284|0.93826 / 0 / 0.061744|-1.22025e-05 / 0.000445674|4.127e-08|
|slope_match_highband|18.5826 / 21.8257|0.0114365|0.00090932|0.93121 / 0.0065254 / 0.062264|-1.86864e-05 / 0.000589629|4.14e-08|
|amplitude_025|18.2054 / 21.0564|0.00761456|0.000804631|0.92593 / 0 / 0.07407|-1.79165e-05 / 0.000450612|4.146e-08|
|amplitude_050|19.5296 / 22.2119|0.0150924|0.0008889369|0.92235 / 0.012607 / 0.065039|-1.81796e-05 / 0.000496803|4.142e-08|
|amplitude_075|20.8647 / 23.7675|0.0405555|0.001025508|0.91222 / 0.031072 / 0.056709|-1.82266e-05 / 0.000681344|4.144e-08|
|cutoff_n17|20.1004 / 23.3769|0.0316633|0.0008829321|0.91215 / 0.014058 / 0.073794|-1.70555e-05 / 0.000460725|4.11e-08|
|cutoff_n18|20.4821 / 23.182|0.0338091|0.0009178068|0.9092 / 0.022451 / 0.068345|-1.76291e-05 / 0.000470015|4.125e-08|

Center is |kx|<=incident99; middle is incident99<|kx|<=Full outgoing99 when that interval exists; tail is the remaining observed range. These are a partition of raw spectral DIFFERENCE energy, not a partition of reflected energy. A phase-dominated error need not be a uniform phase offset; compare phase mean to standard deviation.

The center/tail fractions describe the sampled core spectrum. No 1e-4 continuum angular-spectrum convergence claim is made from513 receiver samples; the physical acceptance metric is the native field at the common receiver locations with independent solver/quadrature/window checks.

## Numerical uncertainty and threshold refinements

Baseline BIE: L60, ppw8, panel order8/self96, h=-2 m, plateau .85. Selected guards: Full4097/8193 at exactly the same receivers; BIE ppw8/order8/self96 versus ppw10/order10/self120; BIE L60/L70 at ppw8/order8. U_numeric is the sum of the three weighted receiver differences, a deterministic convergence indicator rather than a rigorous error bound or statistical confidence interval. Certification requires E+U<1e-3 or E-U>1e-3 and U<1e-4; unrefined cases are explicitly baseline-only. Boundary residuals and linear residuals are saved in numerical_runs.csv. Reused anchor runs retain their original settings/costs.

|case|FK grid difference|BIE ppw/order difference|BIE window difference|sum U|refined FK/BIE complex L2|
|---|---:|---:|---:|---:|---:|
|cutoff_n19|6.557529e-06|1.279262e-07|1.300044e-07|6.81546e-06|0.0011985241|
|slope_match_highband|6.559398e-06|1.068169e-07|1.305891e-07|6.796803e-06|0.00091685799|
|amplitude_075|6.559597e-06|1.331365e-07|1.356663e-07|6.8284e-06|0.0010165948|
|cutoff_n18|6.562191e-06|1.238163e-07|1.314966e-07|6.817504e-06|0.00085420013|

## Matched comparisons and identification limits

|pair|raw height ratio high/low|illuminated slope ratio high/low|illuminated curvature ratio high/low|complex error high/low|
|---|---:|---:|---:|---:|
|height_match|1|1.080532|1.242981|1.266828|
|slope_match|0.9861233|1|1.081496|1.04986|

Height matching and slope matching remove only their specified moment confound. Spectral shape/phase and curvature distribution still differ, so a single scalar validity formula would overstate the evidence. The constant-Kmax amplitude family tests whether bandwidth alone suffices; inspect its full error and illuminated-geometry sequence in validity_metrics.csv.

## Current validity domain and next-stage decision

The guarded master-prefix continuation brackets the1e-3 crossing between Kmax 0.70685835 (E 0.00085275345, U 6.82e-06) and 0.74612826 rad/m (E 0.0011974915, U 6.82e-06). This is a sampled, single-realization conditional domain; no interpolation theorem or universal Kmax threshold is asserted.

Equal periodic height RMS: high/low curvature ratio in the illuminated region 1.243, complex error ratio 1.2668; illuminated slope also changes by a factor 1.0805. Equal illuminated slope RMS: curvature ratio 1.0815, error ratio 1.0499, periodic height ratio 0.98612. Height and RMS slope alone therefore do not organize this whole controlled set. Illuminated curvature/k is a useful candidate descriptor, jointly with high-K amplitude/spectral structure and outgoing broadening; it is not a uniquely identified causal parameter or a fitted validity law.

For added-mode amplitude factors 0/.25/.50/.75/1, E_complex is [0.0008000959;0.000796226;0.0008742465;0.001015381;0.001197491]. Every positive factor has identical physical Kmax 0.74612826 rad/m. Consequently Kmax/k alone cannot determine acceptance. The initial canonical-to-mode14 error decreases even though raw height, slope and Kmax grow; illuminated curvature decreases on that step, illustrating why the illumination and modal arrangement matter.

In this finite controlled set, baseline passing cases have illuminated RMS curvature/k no larger than 0.0017861097, and baseline failing cases no smaller than 0.0019414207. This descriptive gap is candidate screening evidence only: most interior controls are baseline-only, the matched-slope error contrast is modest, and no unseen surface or different source/aperture is certified by that scalar interval. The guarded master-prefix endpoints have descriptors 0.001718438 and 0.0022730243 when both are available.

Canonical/high-cutoff phase shares of complex error energy are 97.093%/98.656%. Core spectral difference energy inside incident99 is 92.864%/91.721%. Their core spectral phase standard deviations exceed the magnitudes of the phase means by factors 26.22/43.51. The discrepancy is thus mainly phase variation in the observed core, not merely a spectral-tail effect or a single removable global phase. No phase adjustment is made.

Preserve the successful FastFK algorithm as a compression of the frozen Full model. Retain FullFK/FastFK only in a measured and numerically guarded physical domain. If the intended PM spectrum includes rejected high-K/amplitude structures, initiate a limited SSA/SSA2 physical-reference comparison instead of additional FullFK compression. This one-realization audit cannot establish how often a target sea state exceeds the domain. No SSA, multi-frequency, Monte Carlo, SSP,3D, BER or production integration is run here.

## Reproduction, verification and provenance

`setup_vertical_project; validate_full_kirchhoff_pm_validity_boundary('init')`; repeated `'next'` executes ONE new baseline job; `'adapt'` appends at most two cutoff cases; finish `'next'`, then `'plan'` and repeated `'refine'`; `'render'` reuses saved data. One-job checkpoints avoid restarting long BIE work after an MCP observation timeout. Outputs are in results/validation/full_kirchhoff_pm_validity_boundary. New source/statistics/metric helpers are validation-only; public interfaces and reference solvers are unchanged.

New numerical solves: wall 3676.2 s, MATLAB CPU 29671.2 s (summed per solve; excludes plotting/orchestration). 5 prior BIE runs reused. Largest single complex dense BIE matrix estimate 212.98 MiB; this is not process peak memory. All run-level wall/CPU/node counts, residuals and matrix estimates are in numerical_runs.csv; Full timings are retained in the numeric MAT state.

Verification: artifact_assertions.json checks completion, eight or fewer nested cutoffs, exact coefficient provenance, matched height/slope moments, native phase identity, finite complex fields, spectral difference partition, residuals and all three convergence guards. tests.csv and code_analysis.txt record MATLAB unit tests and Code Analyzer diagnostics. provenance.json records frozen reference hashes and generated artifact hashes.

Study layout used the experimental-design skill to specify matched input-moment controls and avoid treating one realization as replicated sampling. Workflow reference: Kassis T., Agarwal V., He Y., Patel D., Brueckner A.M. (2026), [Scientific Agent Skills: A Library of Procedural Knowledge for Research Agents](https://doi.org/10.48550/arXiv.2609.00065), current arXiv record verified2026-10-07. This is workflow attribution, not acoustic evidence for the validity boundary.

![cutoff_continuation](../results/validation/full_kirchhoff_pm_validity_boundary/figures/cutoff_continuation.png)

![dimensionless_controls](../results/validation/full_kirchhoff_pm_validity_boundary/figures/dimensionless_controls.png)

![matched_fourier_controls](../results/validation/full_kirchhoff_pm_validity_boundary/figures/matched_fourier_controls.png)

![fk_bie_spectral_difference](../results/validation/full_kirchhoff_pm_validity_boundary/figures/fk_bie_spectral_difference.png)

![native_phase_structure](../results/validation/full_kirchhoff_pm_validity_boundary/figures/native_phase_structure.png)

