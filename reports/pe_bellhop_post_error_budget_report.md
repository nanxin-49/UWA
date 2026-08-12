# Post-error-budget PE--Bellhop rerun

Status: `passed=false`.

The continuous Weyl reference passed independently. The PE input remains the finite spatial-plane source so its error can be separated from sponge error.

Definitions: `dA = 20 log10(|H_sponge|/|H_no_sponge|)` and `dTL = TL_sponge - TL_no_sponge = -dA`. Positive dTL means additional loss.

The complete 108-row window x sponge-ratio x alpha table is available in
[`pe_point_source_error_budget_report.md`](pe_point_source_error_budget_report.md)
and as machine-readable
[`sponge_error_budget.csv`](../results/validation/pe_point_source_error_budget/sponge_error_budget.csv).

- Maximum finite-window TL error: `12.8509 dB`.
- Maximum sponge-only TL change: `16.2231 dB`.
- Maximum final PE--Bellhop TL error: `10.7027 dB`.
- Bellhop--analytic maximum TL error: `7.9244e-07 dB`.

## Receiver-by-receiver decomposition

| Offset (m) | Finite-window dA (dB) | Finite-window dTL (dB) | Sponge dA (dB) | Sponge dTL (dB) | Sponge dphase (rad) | Final PE-BH dTL (dB) |
|---:|---:|---:|---:|---:|---:|---:|
| 0.000 | 12.850856 | -12.850856 | -4.585012 | 4.585012 | -0.980475 | -8.265844 |
| 0.500 | 5.520422 | -5.520422 | -16.223131 | 16.223131 | -0.787373 | 10.702708 |
| 1.000 | -0.613117 | 0.613117 | -6.267654 | 6.267654 | 2.611615 | 6.880771 |
| 2.000 | 9.611861 | -9.611861 | -5.092287 | 5.092287 | -0.804938 | -4.519575 |

## Default sponge terminal-plane energy

Center energy is integrated over `rho <= 2 m`; edge energy is integrated over the ratio-defined sponge band.

| alpha (Np/m) | |H axis| | phase axis (rad) | center energy | edge energy | total energy | dE center (dB) | dE edge (dB) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.000 | 0.00436760095 | 0.012487356 | 0.000107964731 | 0.0304400942 | 0.0690142487 | 0.000000 | 0.000000 |
| 0.150 | 0.00257627569 | -0.967988027 | 6.2021629e-05 | 0.00927474763 | 0.0285068515 | -2.407387 | -5.161439 |
