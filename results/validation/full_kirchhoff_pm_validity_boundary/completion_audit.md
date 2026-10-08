# Completion audit

All requested physics-audit work is complete. Numerical state: 15/15 baseline cases, 8 nested cutoffs, 12/12 guards on 4 selected cases.

- Exact common Fourier coefficients, saved design before collection, and matched periodic-height/illuminated-slope moments: artifact_assertions.json and design.csv.
- Raw, embedded and illuminated surface statistics, physical Kmax independent of numerical Nyquist: validity_metrics.csv and case_data/*_geometry.csv.
- Outgoing theta95/theta99, energy outside incident99, and native complex FK/BIE aperture spectra: validity_metrics.csv and case_data/*_spectrum.csv. Spectrum conclusions are limited to the sampled aperture; field acceptance uses the common receiver locations.
- Native complex, magnitude and phase errors with an exact error-energy decomposition: validity_metrics.csv, raw receiver CSVs and 3 passing MATLAB tests. No phase alignment or complex scalar fit.
- Full 4097/8193 quadrature, BIE ppw8/order8 versus ppw10/order10, and L60/L70 guards: all 12 jobs complete. Summed changes 6.797e-6 to 6.828e-6, below the requested 1e-4 numerical guard. Unrefined cases are labeled baseline_only.
- Guarded master-prefix endpoints: Kmax .706858347 passes (E .000852753450, U .000006817504), .746128255 fails (E .001197491491, U .000006815460). Refined direct errors .000854200132 and .001198524116 respectively.
- Controlled comparisons support illuminated curvature/k together with high-K amplitude/modal structure. Height, slope or Kmax alone do not identify a universal validity formula. No empirical formula is fitted.
- Decision: retain the FullFK/FastFK route within the measured conditional domain; if required spectra include rejected structures, next perform a limited SSA/SSA2 comparison. No failure probability for a PM population is inferred from one realization.
- Scope: prior Fast ranks3/8/padding2 preserved; 9 frozen solver/helper/input hashes match the PM pilot. No production physics, Fast optimization, multi-frequency, Monte Carlo, SSA implementation, BER, SSP or 3D work.

Four new MATLAB files have zero Code Analyzer diagnostics. Five PNG/PDF figure pairs were visually inspected. PROJECT_CONTEXT.md, scripts/README.md and vertical_comm_guide.md link the completed audit. The report includes the experimental-design workflow attribution. provenance.json records all retained artifact and source hashes.
