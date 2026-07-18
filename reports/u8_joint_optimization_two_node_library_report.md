# Joint-Frequency Optimization and U=8 Conditional Model Report

## Status and scope

The optimized F=64 joint-frequency builder, U=8 raw-PM aperture audit, U=8 cached-PE smoke/full validation, and the exact-node U=5/U=8 library are complete. Public defaults, `vertical_channel_model`, `vertical_wape_propagator`, and the communication chain are unchanged.

## 1. Optimized F=64 joint builder

`build_kirchhoff_kstat_joint_model_streaming_vertical` replaces the full `F x F x N_K` covariance and pseudo-covariance tensors by the converged series

\[
C_G(\rho)=\sum_{n\ge1} b_n b_n^H t(\rho)^n,\qquad
P_G(\rho)=\sum_{n\ge1}(-1)^n b_n b_n^T t(\rho)^n,
\]

where `t=C_eta/sigma_eta^2`. It FFTs scalar spatial modes, obtains one global augmented-frequency SVD basis, then performs blockwise small EVDs for K/-K pairs. Factors are stored in single precision. The compact disk representation can be expanded once in RAM for repeated sampling. No effective-K crop is used: at U=5, 65,509 of 65,536 bins are needed for 99.99% spectral energy.

Measured U=5 F=64 results:

| Item | Full builder | Optimized builder |
|---|---:|---:|
| factor build wall time | 631.327 s | 6.502 s (6.569 s outer wall) |
| MATLAB memory snapshot | 6.626 GiB | 1.664 GiB |
| saved factor/cache file | 815 MB | 109.6 MiB |
| series order / basis rank | n/a | 143 / 32 |
| mean local factor rank | n/a | 22 |
| negative/positive eigenvalue energy | baseline | 0 |

Numerical/statistical equivalence:

- small-grid full augmented covariance relative error: `7.80e-8`;
- deterministic sampled augmented-factor error: `6.12e-8`;
- per-frequency incoherent-power relative error: `6.06e-12`;
- adjacent-frequency covariance relative error: `3.62e-12`;
- pseudo-covariance diagonal relative error: `2.28e-5`;
- receiver PDP/LFM correlation against the old builder: `0.998617 / 0.999270`.

The build-time, memory, disk-size, and accuracy targets pass. Runtime expansion costs about 7.79 s and 264 MiB. In the measured short M=8 sample, expanded-factor sampling was `0.948x` the old throughput (compact direct sampling was `0.815x`), so the strict “sampling no slower” goal is not demonstrated; repeated/batched sampling remains the optimization target.

## 2. U=8 raw-PM aperture audit

All grids use the PE spacing `dx=50/128 m`; 128 independent surfaces were centrally cropped to the 50 m/128² PE window.

| PM grid | Kmin (rad/m) | Kpeak/Kmin | discrete variance (m²) | implied Hs (m) | discrete/infinite capture | ideal radial capture | crop energy mean ± std | crop 95% CI |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 100 m / 256² | 0.062832 | 1.7135 | 0.11664 | 1.3661 | 1.0015 | 0.98764 | 1.0176 ± 0.2822 | [0.9687, 1.0664] |
| 150 m / 384² | 0.041888 | 2.5702 | 0.11621 | 1.3636 | 0.99780 | 0.99982 | 0.9579 ± 0.2518 | [0.9143, 1.0015] |
| 200 m / 512² | 0.031416 | 3.4270 | 0.11655 | 1.3656 | 1.0007 | 0.99987 | 1.0076 ± 0.2928 | [0.9569, 1.0583] |

The recommended U=8 grid is **150 m / 384²**. The 100 m aperture has `Kpeak/Kmin<2`; 150 m passes capture and peak resolution, while the 150-to-200 m implied-Hs change is only 0.146%. The PE grid remains 50 m/128² and is decoupled by a same-dx central crop. The crop ratio has appreciable realization scatter, so its recorded confidence interval must remain metadata rather than be replaced by deterministic renormalization.

## 3. U=8 receiver/model implementation

`scripts/validation/validate_u8_conditional_channel_generator_vertical.m` provides:

- smoke `Ltrain=16`, `Ltest=8`, then full `128/128`, with disjoint seed families;
- 4--8 kHz, F=64, PM 150 m/384², PE 50 m/128²;
- explicit same-surface kdomain, independent kstat, and optimized joint kstat through cached PE;
- reflected-scatter-only mean, C, P, correlation, adjacent correlation, phase increments, eigenvalues, numerical rank, and kdomain split-sample floor;
- independent 2,000-realization proper-null tests for U=8;
- physical CIR, shortest circular 99% interval, PDP, LFM, and distribution checks;
- full, 99.9%, and 99% receiver-statistics generators;
- 10,000 H+CIR samples with an explicit zero-PE-call record;
- wall-time, memory, public-path baseline, file-save time, and break-even accounting.

The smoke run (`16/8`) completed as a code-path check but, as expected for only eight held-out samples, gave unstable PDP correlation `0.6795`. It was not used for acceptance. The full `128/128` run passed all scripted gates.

### 3.1 Receiver statistics

| Metric | optimized joint | independent | kdomain split floor |
|---|---:|---:|---:|
| covariance relative error | 0.35543 | 0.97026 | 0.34031 |
| correlation-matrix relative error | 0.37446 | 0.97184 | 0.31392 |
| adjacent-frequency correlation RMSE | 0.01210 | 0.97947 | 0.02503 |
| PDP correlation | 0.98093 | 0.18262 | — |
| LFM correlation | 0.99113 | 0.51431 | — |

Joint covariance is close to the finite-sample kdomain split floor and is decisively better than independent-frequency sampling. The reflected-scatter-only shortest circular 99% interval is `4.5835 ms`, or `0.2910 Tmax`, with `Tmax=15.75 ms` and physical resolution `0.25 ms`.

U=8 properness was tested independently with 2,000 null draws. The observed `||P||F/||C||F` is `0.23512`; central 90/95/99% intervals are `[0.18470,0.27148]`, `[0.17703,0.28022]`, and `[0.16587,0.30013]`. The upper-tail p-value is `0.3318`; properness is not rejected. `P` remains stored and the improper sampler remains available.

### 3.2 Rank and generated channels

| Candidate | Rank | covariance error vs kdomain | PDP | LFM | magnitude KS |
|---|---:|---:|---:|---:|---:|
| full | 42 | 0.35199 | 0.98376 | 0.98681 | 0.0949 |
| 99.9% | 14 | 0.36122 | 0.98205 | 0.98728 | 0.0885 |
| 99% | 11 | 0.35697 | 0.98268 | 0.98707 | 0.0902 |

All candidates pass PDP/LFM, but full rank is retained as the physical default; ranks 14 and 11 are compression candidates only. The 10,000-channel sample time is `0.01365 s`, CIR construction `0.00177 s`, and bundle saving `1.082 s` for a `28.23 MiB` file. The generation stage records zero PE calls. Component algebra error is `1.55e-17`.

### 3.3 U=8 construction and memory

- fixed cached-PE construction: `21.19 s`;
- streaming joint factor build: `44.55 s`;
- series order / global basis rank: `674 / 63`;
- mean/max local factor rank: `45.86 / 46`;
- builder-only peak MATLAB memory: `3.062 GiB`;
- runtime factor expansion: `5.812 s`;
- full runtime/ensemble memory snapshot: `6.858 GiB`;
- U=8 factor storage: `812.63 MiB`, independent diagonal factors `36 MiB`;
- compact cache file: `1.060 GB` (`1,060,424,797` bytes);
- joint training ensemble: `214.24 s`; statistics estimate: `0.0277 s`;
- public F=64 direct realization: `136.75 s`; cached joint realization: `1.674 s`;
- receiver-model build total: `280.00 s`; measured break-even: `2.05` channels;
- receiver condition model: `716.7 KiB`.

The builder itself meets the `<4 GiB` target. Runtime expansion does not: it trades memory for repeated sampling speed. U=8 storage is also larger than the old U=5 815 MB cache because higher roughness raises the series/basis/local ranks and the PM grid has 2.25 times as many K bins. Production should keep the compact representation on disk and consider grouped or streaming runtime sampling instead of unconditional expansion.

## 4. Discrete U=5/U=8 library

`build_conditional_channel_library_vertical` implements the requested schema and stores deterministic components, `mu/C/P`, eigensystem, selected rank, properness result, validation, timing, and the complete one-condition sampling model for every node. `sample_conditional_channel_library_vertical` accepts an exact wind node only and rejects unknown wind speeds; it never silently interpolates or chooses a nearest node.

`scripts/validation/build_u5_u8_conditional_library_vertical.m` saved and smoke-tested `conditional_channel_library_u5_u8_f64.mat` (`2,257,884` bytes). U=5 and U=8 samples both had zero component-sum error. U=6 was explicitly rejected, confirming that no interpolation or nearest-node fallback occurs.

## 5. Regression results

- U=5 F=64 full generator rerun: pass; receiver PDP/LFM `0.9525/0.9839`.
- cached/public same-input: max total-response difference `7.04e-16`, measured `191.1x` speedup.
- joint boundary: joint covariance error `0.0229`, PDP/LFM `0.9950/0.9952`; independent covariance error `0.9524`.
- raw-PM U=5 grid audit rerun: 100 m/256² capture `0.99833`; central-crop mean error `2.04%` over 32 seeds.
- Kirchhoff kstat suite: all 21 checks passed; same-seed error `0`, changed-seed reflection delta `0.011905`.
- public direct-only/direct-plus-reflect: errors `0` and `1.73e-18`; direct path unchanged.

The MATLAB startup reports stale external BELLHOP search-path warnings. They do not affect these project-local tests, but the user MATLAB path should be cleaned separately.

## 6. Current decision

- Optimization: **passed**, except strict batch-sampling throughput parity remains slightly below target in the short benchmark.
- U=8 raw-PM aperture: **passed**, use 150 m/384².
- U=8 condition model: **passed full 128/128 held-out validation**.
- U=5/U=8 library: **built and smoke-tested with exact-node selection**.
- Stable batch generation of H/CIR at U=5 and U=8: **supported**.
- Move to U=10: **physically possible, but first requires a new aperture audit and memory strategy**.

U=5 to U=8 trends are reasonable but modest at this receiver: joint scatter covariance trace rises from `0.10130` to `0.10720`, full numerical rank from `34` to `42`, while kdomain RMS delay stays nearly unchanged (`0.52634` vs `0.52599 ms`) because geometry and propagation remain fixed. The next priority should be communication-performance validation at the two validated nodes. In parallel, runtime factor grouping/streaming should be addressed before U=10; higher wind will further increase low-wavenumber aperture and joint-factor ranks.
