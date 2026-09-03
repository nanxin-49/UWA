# PE--Bellhop PM coefficient-amplitude ensemble comparison

状态：**PASS_WITH_LIMITS**

固定 4 kHz、uniform c=1500 m/s、W=192.1875 m、nx=984、step=0.05 m、Bellhop 5001 beams、sector ±15°。The canonical PM spectral-density column is fixed and each non-reference seed independently draws Gaussian cosine/sine coefficients with variance S(k)*Delta-k. No Hs renormalization, recentering, smoothing, tapering, bandwidth change, source fit, or core-physics modification is applied.

Canonical coefficient source: `E:\MISC\CARPE3D_matlab\Explain\results\validation\bellhop_internal_pm_fixed_realization\fixed_pm_fourier_coefficients.csv`; SHA-256: `1f7eda465e4ae85b8ac038310edf053f2d062563a3b943bfd83bd59d07015f67`; requested Kmax `0.5 rad/m`, realized Kmax `0.471238898038 rad/m`.

Flat Bellhop baseline: source `stage2_verified`, beams `5001`, fingerprint `7e2112a91a8ccbdf4c98cfefde1b93d24976a13bfd7f79d655e311acb5cf2e74`. Every rough case must match this beam count; cached outputs are accepted only when their request fingerprint and output hashes match.

## Per-seed results

| seed | beams | case source | hits/beams | RMS height (m) | RMS slope | max slope | RMS curvature (1/m) | min radius (m) | G_PE | G_Bellhop | delta TL (dB) | delta phase (rad) | complex error | min |u.n| | wall residual (m) | tau (s) |
|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 260001 | 5001 | stage2_verified | 5001/5001 | 0.18639076 | 0.053979813 | 0.1166743 | 0.020904738 | 22.424348 | -0.08686924412+0.9880686784i | -0.6903394321-0.6628643101i | 0.31042904 | -2.2482008 | 1.8366415 | 0.95522675 | 7.1435406e-15 | 0.0683459103843 |
| 260002 | 5001 | cache_verified | 5001/5001 | 0.14918319 | 0.048954473 | 0.097931396 | 0.017881158 | 25.647225 | -0.7121605959-0.7403220493i | -0.676815429-0.7045531878i | 0.43590887 | -0.00069114618 | 0.051471386 | 0.9449822 | 7.1827844e-15 | 0.0685096447302 |
| 260003 | 5001 | cache_verified | 5001/5001 | 0.17774805 | 0.04852184 | 0.13235598 | 0.015932837 | 21.571424 | 1.005079831+0.3768785367i | -0.5285172913+0.7651355499i | 1.2463943 | -1.8165395 | 1.7011896 | 0.96520298 | 7.1591563e-15 | 0.0683938746638 |
| 260004 | 5001 | cache_verified | 5001/5001 | 0.16815136 | 0.051218303 | 0.11074242 | 0.018491691 | 23.621503 | 0.415403883+1.002263915i | 0.1319649462+0.9283356828i | 1.2671241 | -0.25170803 | 0.31239353 | 0.9548748 | 7.1640202e-15 | 0.0683686254273 |
| 260005 | 5001 | cache_verified | 5001/5001 | 0.12831896 | 0.036133042 | 0.092076077 | 0.014420358 | 32.264829 | 0.9554989654-0.1471339972i | 1.0293544+0.05292855675i | -0.55638403 | -0.20416042 | 0.20690466 | 0.95131608 | 7.1091732e-15 | 0.0686731380655 |
| 260006 | 5001 | cache_verified | 5001/5001 | 0.16459025 | 0.043550807 | 0.11292278 | 0.012858297 | 28.278152 | -0.9850370284+0.1244068051i | -1.005467989+0.04682845866i | -0.11899713 | -0.079091257 | 0.079700918 | 0.93303466 | 7.1700814e-15 | 0.0687967874836 |
| 260007 | 5001 | cache_verified | 5001/5001 | 0.16746519 | 0.055759113 | 0.14769008 | 0.021310064 | 17.65363 | 0.7883854044-0.5866465533i | 0.8259000677+0.5333200784i | -0.0037557789 | -1.2130993 | 1.1398256 | 0.95549413 | 7.1579215e-15 | 0.0689446550319 |
| 260008 | 5001 | cache_verified | 5001/5001 | 0.16782459 | 0.044340598 | 0.14069022 | 0.013994196 | 23.154015 | -0.8169500266+0.5699101892i | -0.9716971914+0.2269245678i | -0.015228309 | -0.37969577 | 0.37709239 | 0.96700965 | 7.1113462e-15 | 0.0685659459986 |

## Statistics

| statistic | PE TL (dB) | Bellhop TL (dB) | model delta TL (dB) | model delta phase (rad) |
|---|---:|---:|---:|---:|
| mean / circular mean | 0.11810376 | -0.20258263 | 0.32068639 | -0.68931037 |
| std / circular std | 0.36689277 | 0.30778939 | 0.64881613 | 0.83581852 |
| median / p50 | -0.04810401 | -0.17507684 | 0.15333663 | -0.3157019 |
| p5 | -0.24389503 | -0.60582434 | -0.40329861 | -2.0971193 |
| p25 | -0.091000838 | -0.42568841 | -0.041170513 | -1.3639593 |
| p75 | 0.32900477 | 0.00012549809 | 0.63853024 | -0.17289313 |
| p95 | 0.67564863 | 0.19066907 | 1.2598687 | -0.028131185 |

Mean reflected roughness power: PE 1.030838, Bellhop 0.95652236.

## Checks

- reference_band: PASS
- profile_provenance: PASS
- case_fingerprints: PASS
- uniform_bellhop_numerics: PASS
- reflection_success: PASS
- fields_finite: PASS
- bellhop_geometry: PASS
- reflection_phase: PASS
- beam_state: PASS
- sample_count: PASS
- all: PASS

## Limits

This is a first reduced independent PM coefficient-amplitude ensemble, not a converged ocean Monte Carlo. The native backward-range Cartesian amplitude limitation remains diagnostic; geometry, pressure-release phase, p/q state and delay are the hard checks. Increase the seed count only after this reduced amplitude ensemble is accepted.
