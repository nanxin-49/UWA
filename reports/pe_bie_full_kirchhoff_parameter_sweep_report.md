# Full Kirchhoff deterministic sinusoidal parameter sweep

Status: **validation-only; production PE, BIE, and surface physics unchanged**.

Configuration: 2-D/one-transverse, 4 kHz, c=1500 m/s, Tx z=100 m, Rx z=3 m, Gaussian sigma=0.3 m, pressure-release Dirichlet surface, C2 taper (42/50 m), PE/BIE/Full-Kirchhoff fixed source/grid/mask convention. Full-Kirchhoff and BIE native fields are conjugated exactly once to the frozen comparison representation.

Mask: Stage-0 99% incident-energy footprint, finite samples, and both compared fields above -40 dB relative to the pair peak; energy weights are renormalized on that mask.

## Height sweep (K=0.10 rad/m)

| A (m) | PE-BIE Ec | FK-BIE Ec | FK magnitude L2 | FK phase RMS (rad) | slope AK | curvature AK^2 |
|---:|---:|---:|---:|---:|---:|---:|
| 0.005 | 0.00157783 | 0.000253645 | 0.000177868 | 0.000180822 | 0.0005 | 5e-05 |
| 0.01 | 0.00315757 | 0.00025365 | 0.00017797 | 0.000180728 | 0.001 | 0.0001 |
| 0.02 | 0.00633038 | 0.000253668 | 0.000178362 | 0.000180368 | 0.002 | 0.0002 |
| 0.05 | 0.01609 | 0.000253805 | 0.000180372 | 0.000178554 | 0.005 | 0.0005 |
| 0.1 | 0.0339984 | 0.00025436 | 0.000180892 | 0.000178824 | 0.01 | 0.001 |
| 0.2 | 0.0809081 | 0.000256982 | 0.000176529 | 0.000186783 | 0.02 | 0.002 |

## Wavenumber sweep (A=0.02 m)

| K (rad/m) | PE-BIE Ec | FK-BIE Ec | FK magnitude L2 | FK phase RMS (rad) | slope AK | curvature AK^2 |
|---:|---:|---:|---:|---:|---:|---:|
| 0.05 | 0.00941375 | 0.000253613 | 0.000178226 | 0.000180431 | 0.001 | 5e-05 |
| 0.1 | 0.00633038 | 0.000253668 | 0.000178362 | 0.000180368 | 0.002 | 0.0002 |
| 0.2 | 0.00820412 | 0.000254226 | 0.00017913 | 0.000180389 | 0.004 | 0.0008 |
| 0.3 | 0.0116157 | 0.000256649 | 0.000179103 | 0.000183803 | 0.006 | 0.0018 |
| 0.47 | 0.0229061 | 0.000271269 | 0.00017734 | 0.000205252 | 0.0094 | 0.004418 |
| 0.7 | 0.0484307 | 0.00033158 | 0.000178534 | 0.000279571 | 0.014 | 0.0098 |

## Interpretation

Across the tested deterministic range, Full Kirchhoff is compared directly with the same Helmholtz BIE geometry and convention. The FK residual is reported without fitting or amplitude renormalization. The PE residual is the local phase-screen result relative to BIE.

- FK complex-error range: `0.000253613` to `0.00033158`; FK phase-RMS range: `0.000178554` to `0.000279571 rad`.
- PE complex-error range: `0.00157783` to `0.0809081`.
- FK 2049--4097 integration convergence: maximum complex L2 `3.71961e-05`, maximum phase RMS `2.66886e-05 rad`.

## Direct answers

1. **Full Kirchhoff remains approximately equal to BIE throughout this sweep.** Its maximum complex L2 is `0.00033158`, more than two orders below the largest PE residual.
2. **No Full-Kirchhoff failure boundary is observed** through `A=0.20 m` at low K and through `K=0.70 rad/m` at `A=0.02 m` (maximum nominal `A*K=0.02`, `A*K^2=0.0098 1/m`). The small high-K increase is phase-led but remains close to the BIE/integration floor.
3. **The PE--Full-Kirchhoff difference is primarily the local phase-screen reduction.** At low K the PE residual grows mainly in phase with height; at high K it develops a large magnitude residual, consistent with missing nonlocal coherent redistribution/spectral coupling rather than a single local amplitude factor.
4. **Next step:** a reduced nonlocal Kirchhoff operator is justified as the immediate PE diagnostic target. SSA is not required by any case in the tested envelope because Full Kirchhoff already closes to BIE. A low-cost local PE correction may help weak/low-K finite-angle phase error, but these results do not support it as a replacement for nonlocal coupling in strong-height or high-K cases.

Artifacts: `results/validation/pe_bie_full_kirchhoff_parameter_sweep/pe_bie_full_kirchhoff_parameter_sweep.mat`, `results/validation/pe_bie_full_kirchhoff_parameter_sweep/pe_bie_full_kirchhoff_parameter_sweep.csv`; figures: `results/validation/pe_bie_full_kirchhoff_parameter_sweep/figures/`.
