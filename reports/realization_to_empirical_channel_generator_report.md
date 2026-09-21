# Realization-Based Simulation to Empirical Channel Generator

Date: 2026-06-11

## Executive Summary

The current vertical underwater acoustic project has reached a realization-based simulation platform:

- PE/WAPE propagation still produces the direct and reflected frequency responses.
- The Kirchhoff rough-surface reflection can be evaluated through the original spatial phase screen or the implicit k-domain interface.
- Boundary coupling and incident-weighted redistribution diagnostics quantify how the Kirchhoff phase screen broadens transverse wavenumber content.
- Multi-sea-state Monte Carlo statistics now run over `Hs x wind x sea_seed`.
- The D2 communication chain is fixed at `receive_window_mode='peak_sync'` and `ebn0_reference='rx_clean'`.

The next stage is empirical statistics-based channel generation: use the Monte Carlo sample ensemble to generate fast random channel summaries, wideband `H_f(f)` samples, and tap-level samples without rerunning PE/WAPE.

## Fixed Communication Chain

D2 is treated as fixed for this phase:

- `receive_window_mode='peak_sync'`
- `ebn0_reference='rx_clean'`
- QPSK, `M=4`
- `EbN0_dB_list=0:2:20`
- No changes to the receive-window, equalizer, or noise-reference logic during this stage

This keeps communication Monte Carlo results comparable across sea conditions.

## C3.5 Monte Carlo Sweep

Input grid:

- `sea_hs_target=[0.05, 0.5, 1.0]`
- `sea_wind_speed=[3, 5, 8, 12]`
- `mc_count=16`
- `seed_list=12345:12360`
- `enable_wideband=true`
- `f_band_hz=[4000 8000]`
- `Nf_min=Nf_max=16`
- `f_ref_hz=6000`
- `surface_boundary_model='kirchhoff_kdomain'`
- C2 and C2.5 diagnostics enabled

Validation:

- `run_summary_table` rows: `192`
- sea-state conditions: `12`
- seeds per condition: `16`
- maximum `H_f-(H_direct_f+H_reflect_f)` invariant error: `1.551583845779546e-17`
- maximum direct-path drift: `1.390486644391991e-17`

Core C3.5 condition summary:

| Hs | wind | mean \|H(f_ref)\| | std \|H(f_ref)\| | reflect rms delta k | high-k fraction | tap RMS delay |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.05 | 3 | 0.045284 | 0.022244 | 3.3206 | 0.2462 | 36.334 |
| 0.05 | 5 | 0.029723 | 0.019989 | 2.6774 | 0.13695 | 31.644 |
| 0.05 | 8 | 0.023033 | 0.013071 | 2.5222 | 0.10977 | 34.560 |
| 0.05 | 12 | 0.025277 | 0.015919 | 2.5068 | 0.10661 | 34.779 |
| 0.5 | 3 | 0.042238 | 0.023532 | 6.5736 | 0.82985 | 142.51 |
| 0.5 | 5 | 0.044515 | 0.019683 | 6.4730 | 0.81499 | 99.835 |
| 0.5 | 8 | 0.045544 | 0.018677 | 4.9028 | 0.56899 | 97.101 |
| 0.5 | 12 | 0.051057 | 0.026667 | 4.3088 | 0.46143 | 89.274 |
| 1.0 | 3 | 0.045586 | 0.018760 | 6.5599 | 0.82548 | 125.33 |
| 1.0 | 5 | 0.047975 | 0.020990 | 6.5822 | 0.82880 | 118.31 |
| 1.0 | 8 | 0.045516 | 0.025720 | 6.3567 | 0.80001 | 120.02 |
| 1.0 | 12 | 0.054750 | 0.021414 | 6.0107 | 0.74638 | 109.35 |

Generated core heatmaps:

- `c35_core_abs_H_ref_mean_heatmap.png`
- `c35_core_abs_H_ref_std_heatmap.png`
- `c35_core_reflect_rms_delta_k_mean_heatmap.png`
- `c35_core_reflect_high_k_fraction_mean_heatmap.png`
- `c35_core_tap_rms_delay_mean_heatmap.png`

## Representative Communication Monte Carlo

Representative sea states:

- weak: `Hs=0.05`, `wind=5`
- mid: `Hs=0.5`, `wind=8`
- strong: `Hs=1.0`, `wind=12`

Run settings:

- `mc_count=16`
- scenarios: `direct_only`, `direct_plus_reflect`
- total channel runs: `3 x 2 x 16 = 96`
- maximum channel invariant error: `1.387778780781446e-17`
- maximum direct-only reflected response: `0`

Selected BER/SER results:

| condition | scenario | Eb/N0 | BER mean | BER std | SER mean | SER std |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| weak | direct_only | 0 | 0.11131 | 0.0087823 | 0.15050 | 0.010392 |
| weak | direct_only | 10 | 0 | 0 | 0 | 0 |
| weak | direct_only | 20 | 0 | 0 | 0 | 0 |
| weak | direct_plus_reflect | 0 | 0.11606 | 0.007598 | 0.15887 | 0.0075971 |
| weak | direct_plus_reflect | 10 | 3.125e-05 | 0.000125 | 6.25e-05 | 0.00025 |
| weak | direct_plus_reflect | 20 | 0 | 0 | 0 | 0 |
| mid | direct_only | 0 | 0.11338 | 0.0057315 | 0.15269 | 0.0086465 |
| mid | direct_only | 10 | 0 | 0 | 0 | 0 |
| mid | direct_plus_reflect | 0 | 0.17438 | 0.089988 | 0.24050 | 0.13120 |
| mid | direct_plus_reflect | 10 | 0.045906 | 0.10180 | 0.06550 | 0.14976 |
| mid | direct_plus_reflect | 20 | 0.037469 | 0.099933 | 0.053562 | 0.14439 |
| strong | direct_only | 0 | 0.11062 | 0.007719 | 0.15063 | 0.0090839 |
| strong | direct_only | 10 | 3.125e-05 | 0.000125 | 6.25e-05 | 0.00025 |
| strong | direct_only | 20 | 0 | 0 | 0 | 0 |
| strong | direct_plus_reflect | 0 | 0.15978 | 0.056577 | 0.21863 | 0.080355 |
| strong | direct_plus_reflect | 10 | 0.018125 | 0.032442 | 0.024375 | 0.044027 |
| strong | direct_plus_reflect | 20 | 0.008750 | 0.020451 | 0.012187 | 0.028374 |

The reflected mid and strong channels retain residual BER at high Eb/N0 for some seed realizations, indicating residual ISI/noise-enhancement effects under the current reduced-grid channel and MMSE receiver.

## C4 Upgrade

C4 now supports three sample modes:

- `summary`: low-dimensional empirical bootstrap summaries.
- `wideband_hf`: bootstrap of stored wideband `H_f(f)` samples from C3.5.
- `tap_level`: derives baseband/tap samples from bootstrapped wideband `H_f(f)`.

The upgraded builder detects whether each sea-state condition has stored wideband samples. The completed C3.5 result has `has_wideband_hf=true` for all 12 conditions.

C4 validation:

- input source: `sweep_monte_carlo_surface_channel_vertical_result.mat`
- conditions checked: 3
- generated summary samples per condition: 1000
- output size: `219191` bytes
- wideband/tap mode validation:
  - `H_f_sample_size=[16 50]`
  - `idx_f_ref=8`
  - `tap_count_mean=501.92`
  - `tap_rms_delay_mean=40.0435`
  - `tap_peak_fraction_mean=0.9857`
  - tap metrics finite: true
  - no PE/WAPE propagation run: true

## Interpretation

This phase completes the realization-based platform:

- PE/WAPE generates physical channel realizations.
- Kirchhoff k-domain diagnostics quantify rough-surface angular-spectrum effects.
- C3/C3.5 aggregate realization statistics.
- Communication Monte Carlo evaluates the fixed D2 link under selected sea states.

The next phase should emphasize empirical statistics-based generation:

- improve sea-state interpolation over the C3.5 grid;
- increase `mc_count` for production statistics;
- fit compact conditional distributions only after the empirical sample base is large enough;
- evaluate C4-generated `H_f` and tap samples against held-out PE/WAPE realizations;
- keep claims limited to empirical reproduction of the implemented model, not new scattering theory.

## Limitations

- All C3.5/C4 statistics are tied to the current reduced grid and parameter range.
- `mc_count=16` is better than the previous `mc_count=4`, but still not a converged ocean-channel distribution.
- Wind trends are conditional on PM spectral shape after `sea_hs_target` rescaling.
- C4 is still an empirical bootstrap generator, not a closed-form statistical channel model.
- No T-matrix, SSA/NLSSA, or physical scattering cross-section model has been implemented.
