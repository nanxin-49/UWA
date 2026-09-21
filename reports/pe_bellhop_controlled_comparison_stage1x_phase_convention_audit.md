# PE--Bellhop Stage 1X phase/convention audit

状态：**CONVENTION_DISCREPANCY_IDENTIFIED**

本报告只读取三组已完成的 Stage-1 MAT；未调用 PE/Bellhop，未修改任何场或核心公式。

## Fixed-transform diagnostics

| transform | max E_G (M99) | max phase RMS (rad) | max TL RMS (dB) | min rho | max E_aligned | closure |
|---|---:|---:|---:|---:|---:|---|
| PE | 0.47250886 | 0.47919524 | 0.0045698832 | 0.91922588 | 0.40191482 | NO |
| conj(PE) | 0.0040934078 | 0.0040589298 | 0.0045698832 | 0.99999242 | 0.0038955043 | YES |
| PE(-x) | 0.47250886 | 0.47919524 | 0.0045698832 | 0.91922588 | 0.40191482 | NO |
| conj(PE(-x)) | 0.0040934078 | 0.0040589298 | 0.0045698832 | 0.99999242 | 0.0038955043 | YES |

Per-case values are in stage1x_metrics.csv in the Stage1X result directory.

## Per-case four-way metrics

| A (m) | transform | E_G | phase RMS (rad) | TL RMS (dB) | rho_shape | rho_raw | phi0 (rad) | E_aligned | phase P95 | TL P95 |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.01 | PE | 0.47250886 | 0.47919524 | 0.0045698832 | 0.91922588 | 0.88837649 | -0.25980608 | 0.40191482 | 0.66246812 | 0.0096624181 |
| 0.01 | conj(PE) | 0.0040934078 | 0.0040589298 | 0.0045698832 | 0.99999242 | 0.99999163 | -0.0012563958 | 0.0038955043 | 0.010911606 | 0.0096624181 |
| 0.01 | PE(-x) | 0.47250886 | 0.47919524 | 0.0045698832 | 0.91922588 | 0.88837649 | -0.25980608 | 0.40191482 | 0.66246812 | 0.0096624181 |
| 0.01 | conj(PE(-x)) | 0.0040934078 | 0.0040589298 | 0.0045698832 | 0.99999242 | 0.99999163 | -0.0012563958 | 0.0038955043 | 0.010911606 | 0.0096624181 |
| 0.005 | PE | 0.23876341 | 0.23959927 | 0.0022849061 | 0.9792833 | 0.97149714 | -0.12618585 | 0.20354798 | 0.33123391 | 0.0048266453 |
| 0.005 | conj(PE) | 0.0020412258 | 0.0020240691 | 0.0022849061 | 0.9999981 | 0.99999792 | -0.00060381234 | 0.0019497952 | 0.0054392181 | 0.0048266453 |
| 0.005 | PE(-x) | 0.23876341 | 0.23959927 | 0.0022849061 | 0.9792833 | 0.97149714 | -0.12618585 | 0.20354798 | 0.33123391 | 0.0048266453 |
| 0.005 | conj(PE(-x)) | 0.0020412258 | 0.0020240691 | 0.0022849061 | 0.9999981 | 0.99999792 | -0.00060381234 | 0.0019497952 | 0.0054392181 | 0.0048266453 |
| 0.0025 | PE | 0.11969675 | 0.11980048 | 0.0011424793 | 0.99478804 | 0.99283648 | -0.062648514 | 0.10209664 | 0.16561743 | 0.0024133293 |
| 0.0025 | conj(PE) | 0.0010195118 | 0.001010957 | 0.0011424793 | 0.99999952 | 0.99999948 | -0.00029659367 | 0.00097539605 | 0.0027163315 | 0.0024133293 |
| 0.0025 | PE(-x) | 0.11969675 | 0.11980048 | 0.0011424793 | 0.99478804 | 0.99283648 | -0.062648514 | 0.10209664 | 0.16561743 | 0.0024133293 |
| 0.0025 | conj(PE(-x)) | 0.0010195118 | 0.001010957 | 0.0011424793 | 0.99999952 | 0.99999948 | -0.00029659367 | 0.00097539605 | 0.0027163315 | 0.0024133293 |

## Weak-response coefficient audit

At 4 kHz, k=16.7551608 rad/m; 4*k/sqrt(2)=47.3907513 rad/m.

| A (m) | RMS D_PE | RMS D_BH | RMS(D_PE-D_BH) | RMS(D_PE+D_BH) | RMS mirror | RMS negative mirror | phase RMS/A |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.01 | 25.167706 | 24.449702 | 49.614002 | 0.92387392 | 49.614002 | 0.92387392 | 47.919524 |
| 0.005 | 25.167454 | 24.450047 | 49.614108 | 0.9227282 | 49.614108 | 0.9227282 | 47.919853 |
| 0.0025 | 25.167328 | 24.450179 | 49.614116 | 0.92233824 | 49.614116 | 0.92233824 | 47.920192 |

The negative (conjugate) response is the only physically distinguishable fixed relation that closes all three amplitudes; its mirror-conjugate is numerically degenerate for this centered weak profile. This result is diagnostic and does not rewrite the original Stage-1 PASS/FAIL.

## Provenance audit

- Physical transverse coordinate: x_PE; Bellhop receiver depth is z_BH'=-x_PE, written sorted and inverse-permuted back to PE order.
- Wall parameter/profile: s=z_BH'; Gamma(s)=[100-A*cos(0.1*s),s]; post-wall map T(r,z)=(2R0-r,-z).
- Rough and flat fields use the same existing Gaussian .sbp, source geometry X, coherent run type C, and range selector at 103 m (102 m guard column).
- Source-pattern fingerprint: samples=2401, clip=-120 dB, canonical coefficient SHA-256=1f7eda465e4ae85b8ac038310edf053f2d062563a3b943bfd83bd59d07015f67.
- G_BH uses frozen phase_sign=-1; no per-case normalization, fitting, or amplitude calibration is applied.
- Stage-0 flat denominator and AS-defined M95/M99 footprints are reused for every A and transform.
- Stage-1 geometry/finite checks are retained from each authoritative MAT; this audit does not relax them.

## Decision

**CONVENTION_DISCREPANCY_IDENTIFIED** — Fixed transform conj(PE) closes all three A cases; the mirror-conjugate candidate is numerically degenerate for this centered/even receiver realization.

Before any Stage-1 relabeling, the next permitted action is a narrowly scoped convention-only fix and rerun of the affected frozen Stage-0/Stage-1 cases. Stage 2--7 remain locked by the revised Goal.
