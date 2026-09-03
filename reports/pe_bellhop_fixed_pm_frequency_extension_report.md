# PE--Bellhop fixed-PM frequency extension

状态：**PASS_WITH_MODEL_DISCREPANCY**

固定 seed=260001、U=6 m/s、span=160 m、master N=4097、realized Kmax=0.471238898 rad/m；uniform c=1500 m/s、z_tx=100 m、z_rx=3 m、sigma=0.3 m、W=192.1875 m、nx=984、step=0.05 m、beam=5001、Bellhop sector=[-15,15] deg。PE 与 Bellhop 核心均未修改。

## Frequency results

| f (kHz) | G_PE | G_Bellhop | delta TL (dB) | delta phase (rad) | complex error | wall residual (m) | min |u.n| | max |kappa| (1/m) | tau (s) |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 4 | -0.0868692441183+0.988068678385i | -0.690339432108-0.66286431011i | 0.31042904 | -2.2482008 | 1.8366415 | 7.1435406e-15 | 0.95522675 | 0.042525813 | 0.0683459103843 |
| 6 | 0.786733080258-0.602897316542i | -0.872241809293+0.39282002283i | 0.30828546 | 2.9108725 | 2.022601 | 7.1435406e-15 | 0.95522675 | 0.042525813 | 0.0683459103843 |
| 8 | -0.974950225788-0.172791195008i | 0.0388047943136+0.955210495703i | 0.3048223 | 1.7868077 | 1.5864083 | 7.1435406e-15 | 0.95522675 | 0.042525813 | 0.0683459103843 |

## Native numerical diagnostics

- 4 kHz: PE seam jump 0.078096916 m, outer-5% reflected energy 3.0394899e-06; Bellhop grazing fraction 0, phase-jump error 7.1054274e-15 rad, q residual 0, p/q rotation errors 0 / 0, min post-wall dr 0.047105732 m.
- 6 kHz: PE seam jump 0.078096916 m, outer-5% reflected energy 5.5455936e-12; Bellhop grazing fraction 0, phase-jump error 7.1054274e-15 rad, q residual 0, p/q rotation errors 0 / 0, min post-wall dr 0.047105732 m.
- 8 kHz: PE seam jump 0.078096916 m, outer-5% reflected energy 2.7467523e-13; Bellhop grazing fraction 0, phase-jump error 7.1054274e-15 rad, q residual 0, p/q rotation errors 0 / 0, min post-wall dr 0.047105732 m.

## Checks

- profile_provenance: PASS
- pe_finite: PASS
- bellhop_finite: PASS
- bellhop_geometry: PASS
- reflection_phase: PASS
- beam_state: PASS
- all: PASS

The cross-model residual is intentionally not forced toward zero. This stage tests finite, geometrically valid operation at 4/6/8 kHz with the same PM realization; no group delay is inferred from these three frequencies.
