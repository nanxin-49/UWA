# PE--Bellhop PM paired multi-realization comparison

状态：**PASS_WITH_LIMITS**

> Provenance notice (2026-09-09): these eight numerical cases retain the
> historical Bellhop point-source `R` convention. The active ensemble engine
> now uses line-source `X` and includes source geometry in its cache fingerprint.

固定 4 kHz、uniform c=1500 m/s、W=192.1875 m、nx=984、step=0.05 m、Bellhop 5001 beams、sector ±15°。Seed 260001 is the canonical fixed realization; the remaining paired profiles preserve the canonical per-mode spectral amplitudes and assign deterministic random phases. No Hs renormalization, recentering, smoothing, tapering, bandwidth change, source fit, or core-physics modification is applied.

Canonical coefficient source: `E:\MISC\CARPE3D_matlab\Explain\results\validation\bellhop_internal_pm_fixed_realization\fixed_pm_fourier_coefficients.csv`; SHA-256: `1f7eda465e4ae85b8ac038310edf053f2d062563a3b943bfd83bd59d07015f67`; requested Kmax `0.5 rad/m`, realized Kmax `0.471238898038 rad/m`.

## Per-seed results

| seed | RMS height (m) | RMS slope | max slope | RMS curvature (1/m) | min radius (m) | G_PE | G_Bellhop | delta TL (dB) | delta phase (rad) | complex error | min |u.n| | wall residual (m) | tau (s) |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 260001 | 0.18639076 | 0.053979813 | 0.1166743 | 0.020904738 | 22.424348 | -0.08686924412+0.9880686784i | -0.6903394321-0.6628643101i | 0.31042904 | -2.2482008 | 1.8366415 | 0.95522675 | 7.1435406e-15 | 0.0683459103843 |
| 260002 | 0.1864099 | 0.053970401 | 0.1206493 | 0.020907036 | 20.296859 | -0.3150082321-1.017332307i | -0.04916982396-0.9392397426i | 1.0794629 | -0.24797559 | 0.29459189 | 0.95340823 | 7.129065e-15 | 0.0687406406284 |
| 260003 | 0.18639322 | 0.053977085 | 0.12824989 | 0.020904306 | 21.930321 | -0.9422923845-0.3825176125i | 0.08090181307-0.9460315262i | 0.59643361 | -1.2704861 | 1.2302541 | 0.94968502 | 7.1616774e-15 | 0.0687759473891 |
| 260004 | 0.18638897 | 0.053975931 | 0.16246669 | 0.020894914 | 19.12951 | 1.101716774-0.1259691411i | 0.8516869041+0.3454066459i | 1.6309051 | -0.49913131 | 0.58057222 | 0.9512267 | 7.1637909e-15 | 0.0684184390875 |
| 260005 | 0.18640895 | 0.053985513 | 0.14202974 | 0.020902493 | 22.059197 | -0.573834416-0.8278513819i | -0.2270327121-0.9533663291i | 0.23829758 | -0.37233776 | 0.37633308 | 0.95970171 | 7.1102543e-15 | 0.0685029728458 |
| 260006 | 0.18640646 | 0.053971362 | 0.15792191 | 0.020896771 | 18.775633 | 0.9424547269-0.5089720029i | 0.8406920857-0.4197604158i | 1.1372492 | -0.032081249 | 0.14402061 | 0.96194311 | 7.2064783e-15 | 0.0686862682193 |
| 260007 | 0.18639037 | 0.053972601 | 0.13249897 | 0.020898945 | 21.148727 | 0.8998201823+0.5840299295i | 0.7672351169+0.5468266274i | 1.1273642 | -0.043506019 | 0.14615937 | 0.9480794 | 7.1668142e-15 | 0.0686436288502 |
| 260008 | 0.18642266 | 0.053970566 | 0.12005875 | 0.020907156 | 23.623456 | 0.7734620758-0.5224296992i | 1.047347396-0.2051719974i | -1.1642988 | -0.400607 | 0.39271311 | 0.94986414 | 7.1633815e-15 | 0.0686915805209 |

## Statistics

| statistic | PE TL (dB) | Bellhop TL (dB) | model delta TL (dB) | model delta phase (rad) |
|---|---:|---:|---:|---:|
| mean / circular mean | 0.27383953 | -0.34564083 | 0.61948036 | -0.54961718 |
| std / circular std | 0.4814091 | 0.40067239 | 0.86054857 | 0.69252932 |
| median / p50 | 0.34653449 | -0.48386588 | 0.83794826 | -0.38647238 |
| p5 | -0.41409642 | -0.66571585 | -0.67339008 | -1.9060006 |
| p25 | 0.029583708 | -0.5345836 | 0.29239617 | -0.69197001 |
| p75 | 0.59997077 | -0.32974636 | 1.1298355 | -0.1968582 |
| p95 | 0.79703174 | 0.30614738 | 1.4581255 | -0.036079919 |

Mean reflected roughness power: PE 1.0707181, Bellhop 0.92710726.

## Checks

- reference_band: PASS
- profile_provenance: PASS
- fields_finite: PASS
- bellhop_geometry: PASS
- reflection_phase: PASS
- beam_state: PASS
- sample_count: PASS
- all: PASS

## Limits

This is a paired fixed-band phase-ensemble smoke/first statistical comparison, not a claim of an independently sampled PM amplitude ensemble. The native backward-range Cartesian amplitude limitation remains diagnostic; geometry, pressure-release phase, p/q state and delay are the hard checks. A larger M or a physically sourced coefficient ensemble can be added only after this paired workflow is accepted.
