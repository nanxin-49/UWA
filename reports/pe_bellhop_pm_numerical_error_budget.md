# PE--Bellhop PM Stage 0E numerical error budget freeze

状态：**PASS**

本阶段冻结后续 PM 对比前的 PE、Bellhop 和 flat/source 数值误差预算；所有 PM profile density 均由同一 seed-260001 Fourier realization 直接采样。

## PE budget

| family | case | W (m) | nx | step (m) | dx (m) | relative TL (dB) | relative phase (rad) | seam jump (m) | outer 5% energy |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| PE | W160 | 160 | 820 | 0.05 | 0.19512195 | 0 | 0 | 0.012492327 | 3.0190542e-05 |
| PE | W192p1875 | 192.1875 | 984 | 0.05 | 0.1953125 | 9.4630024e-05 | -1.6190678e-05 | 0.078096916 | 3.0394899e-06 |
| PE | W320 | 320 | 1638 | 0.05 | 0.1953602 | -6.9437319e-05 | -1.5220385e-05 | 0.021718511 | 4.6168265e-09 |
| PE | grid492 | 192.1875 | 492 | 0.05 | 0.390625 | -0.030229504 | -0.0012228794 | 0.080794571 | 1.0747022e-05 |
| PE | grid984 | 192.1875 | 984 | 0.05 | 0.1953125 | 0 | 0 | 0.078096916 | 3.0394899e-06 |
| PE | step0p1 | 192.1875 | 984 | 0.1 | 0.1953125 | -6.7502923e-13 | 9.3166973e-14 | 0.078096916 | 3.0394899e-06 |
| PE | step0p05 | 192.1875 | 984 | 0.05 | 0.1953125 | 0 | 0 | 0.078096916 | 3.0394899e-06 |
| PE | step0p025 | 192.1875 | 984 | 0.025 | 0.1953125 | 2.3761029e-12 | -5.3720064e-14 | 0.078096916 | 3.0394899e-06 |

## Bellhop budget

| profile N | beams | step (m) | relative TL (dB) | relative phase (rad) | rays | min |u.n| | grazing fraction | kappa max (1/m) | p curvature kick | q residual | wall residual (m) | min post dr (m) |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 2049 | 5001 | 0.1 | 0.0060473606 | -6.7645572e-05 | 5001 | 0.8319092 | 0 | 0.042518405 | 8.4683955 | 0 | 7.2575631e-15 | 0.079398631 |
| 2049 | 5001 | 0.05 | 0.0060473606 | -6.7645572e-05 | 5001 | 0.8319092 | 0 | 0.042518405 | 8.4683955 | 0 | 7.2717039e-15 | 0.039699316 |
| 2049 | 10001 | 0.1 | 0.0060277602 | -0.00010122792 | 10001 | 0.8319092 | 0 | 0.042518405 | 8.4683955 | 0 | 7.2575631e-15 | 0.079398627 |
| 2049 | 10001 | 0.05 | 0.0060277602 | -0.00010122792 | 10001 | 0.8319092 | 0 | 0.042518405 | 8.4683955 | 0 | 7.2879268e-15 | 0.039699313 |
| 4097 | 5001 | 0.1 | 2.0661523e-05 | 3.3561129e-05 | 5001 | 0.83190125 | 0 | 0.042525813 | 8.469353 | 0 | 7.2420764e-15 | 0.0793969 |
| 4097 | 5001 | 0.05 | 2.0661523e-05 | 3.3561129e-05 | 5001 | 0.83190125 | 0 | 0.042525813 | 8.469353 | 0 | 7.2076213e-15 | 0.03969845 |
| 4097 | 10001 | 0.1 | 0 | 0 | 10001 | 0.83190125 | 0 | 0.042525813 | 8.4694773 | 0 | 7.2420764e-15 | 0.07939665 |
| 4097 | 10001 | 0.05 | 0 | 0 | 10001 | 0.83190125 | 0 | 0.042525813 | 8.4694773 | 0 | 7.2076213e-15 | 0.039698325 |

## Flat/source budget

- Stage 0C max axis-Q TL residual: `0.2609997 dB`; phase: `0.0005149472 rad`; normalized offset magnitude residuals: `0.0059687686 / 0.0055544581`.
- Legacy fixed-realization N=2049 -> 4097 convergence CSV is copied as `pm_fixed_realization_budget.csv`; its accepted field residual is not used to recalibrate any Stage 1 result.

## Frozen checks

- pe_window: PASS
- pe_grid: PASS
- pe_step: PASS
- bellhop: PASS
- flat: PASS
- all: PASS

The reported tolerances are frozen validation budgets: PE family 0.1 dB / 0.02 rad, Bellhop family 0.1 dB / 0.02 rad, and the known flat backward-range influence allowance 0.30 dB. They are not tuned from rough-wall cross-model results.
