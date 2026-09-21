# Bellhop Flat-Surface Ray and TL Visualization

Visualization checks passed: `true`. The source PE/Bellhop strict matrix remains `false`; this display does not change that validation result and does not modify the PE or communication chain.

## Environment and Bellhop modes

- Water depth 100.0 m, uniform sound speed 1500.0 m/s, Tx depth 80.0 m, frequency 4000.0 Hz.
- `R`: 51 central rays; `C`: coherent TL; `I`: incoherent TL.
- TL grid 199 x 241, range 0.25--12.00 m, depth 0.50--99.50 m.
- `ZBOX=99.9 m < 100.0 m`; the explicit 0.010 m ray step terminates displayed rays before the seabed.
- Both field plots use 20--80 dB; the 1.0 m source neighborhood is masked.

## Beam-count convergence

| Mode | Beams low/high | TL RMS difference (dB) | |difference| 95th percentile (dB) | Valid points |
|---|---:|---:|---:|---:|
| coherent | 5001/10001 | 0.190410 | 0.000010 | 46374 |
| incoherent | 5001/10001 | 0.101930 | 0.000001 | 46374 |

## Receiver-point TL

| x (m) | C field TL | I field TL | Arrival synthesis TL | C field - arrival | Globally scaled PE total TL |
|---:|---:|---:|---:|---:|---:|
| 3.0 | 33.80482 | 34.85412 | 33.80956 | -0.00474 | 32.85420 |
| 6.0 | 37.07681 | 34.87424 | 37.12172 | -0.04491 | 36.71302 |
| 9.0 | 49.99538 | 34.90819 | 50.03151 | -0.03613 | 43.03165 |

## Automated checks

| Check | Value | Relation | Limit | Passed |
|---|---:|:---:|---:|:---:|
| ray_count | 0 | == | 0 | 1 |
| ray_direct_present | 0 | == | 0 | 1 |
| ray_surface_present | 0 | == | 0 | 1 |
| ray_bottom_absent | 0 | == | 0 | 1 |
| shade_metadata | 0 | <= | 1e-06 | 1 |
| shade_finite | 0 | == | 0 | 1 |
| shade_axes_positive_down | 0 | == | 0 | 1 |
| shared_color_limits | 0 | == | 0 | 1 |
| coherent_tl_rms | 0.19041019 | <= | 1 | 1 |
| coherent_tl_p95 | 1.0183056e-05 | <= | 3 | 1 |
| incoherent_tl_rms | 0.1019304 | <= | 0.5 | 1 |
| incoherent_tl_p95 | 9.0625631e-07 | <= | 1.5 | 1 |
| receiver_field_arrival_tl | 0.044907742 | <= | 0.5 | 1 |
| saved_pe_invariant | 3.5762241e-18 | <= | 1e-10 | 1 |

## Figures

![Bellhop ray geometry](bellhop_ray_geometry.png)

![Bellhop coherent and incoherent TL](bellhop_tl_fields.png)

![Receiver-depth TL slice](bellhop_receiver_depth_tl.png)

## Limitations

Only the three saved PE receiver responses are overlaid on the TL slice; no nonexistent PE 2-D field is synthesized. Bellhop TL is a single-frequency 4 kHz result, not a 3--5 kHz wideband PDP. Coherent TL contains phase interference, while incoherent TL is an energy envelope; they are not interchangeable.
