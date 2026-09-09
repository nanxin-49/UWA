# PE--Bellhop PM Stage 1A Tier-1 reflected-ratio comparison

状态：**PASS_WITH_LIMITS**

> Provenance notice (2026-09-09): this report is the historical point-source
> `R` run. The active Tier-1 code now defaults to line-source `X`; its executed
> seed-260001 result is documented in `bellhop_source_geometry_rx_audit_report.md`.

固定条件：4 kHz、uniform c=1500 m/s、z_tx=100 m、z_rx=3 m、sigma=0.3 m、xw=192.1875 m、nx=984、step=0.05 m、seed=260001、profile N=4097、beam=10001。PE 使用 Kirchhoff phase screen，Bellhop 使用 native Reflect2D local-specular internal wall；只比较 reflected-only rough/flat ratios。

## Primary metrics

| quantity | value |
|---|---:|
| G_PE | -0.08686924411828591+0.9880686783848931i |
| G_Bellhop | -0.6903390884494519-0.6628639546907513i |
| delta TL (dB) | 0.31043352 |
| delta phase (rad) | -2.2482008 |
| complex relative error | 1.836642 |

## Bellhop rough-wall diagnostics

wall residual `7.2076213e-15 m`; min `|u.n|` `0.83190125`; grazing fraction `0`; max `|kappa|` `0.042525813 1/m`; min post-wall `dr` `0.039698325 m`; pressure-release phase error `7.1054274e-15 rad`; q residual `0`; p/q rotation errors `0 / 0`; tau `0.0683459103843 s`.

## Checks

- profile_provenance: PASS
- bellhop_geometry: PASS
- reflection_phase: PASS
- beam_state: PASS
- field_finite: PASS
- all: PASS

The cross-model delta is intentionally not a closeness gate: it quantifies the expected Kirchhoff phase-screen versus local-specular Gaussian-beam model discrepancy after each flat denominator removes source/absolute-amplitude factors.
