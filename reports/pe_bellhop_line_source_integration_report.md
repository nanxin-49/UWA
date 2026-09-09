# PE--Bellhop line-source integration and cleanup report

Date: 2026-09-09

Status: **PASS_WITH_PROVENANCE_LIMIT**

## What was integrated

The validated amplitude correction is a physical source-geometry selection,
not an empirical gain: active one-transverse PE--Bellhop comparison entrypoints
now request Bellhop `RunType(4)=X` (Cartesian line source). The Gaussian `.sbp`,
PE source, grids, steps, receiver positions, internal-wall geometry,
`Reflect2D`, `p/q`, `InfluenceGeoHatCart`, and rough-surface formulas are
unchanged.

| file | integration role |
|---|---|
| `scripts/validation/support/write_bellhop_unfolded_gaussian_env_vertical.m` | Writes explicit `R` or `X`; default remains `R` for historical callers. |
| `scripts/validation/support/run_bellhop_native_sinusoidal_wall_vertical.m` | Writes the same explicit source-geometry selector for native ATI covariance tests. |
| `scripts/validation/validate_pe_bellhop_incident_field_vertical.m` | Active incident-plane comparison defaults to `X`, forwards it to Bellhop, isolates X output names, and does not apply the historical point-source global constant to X. |
| `scripts/validation/validate_pe_bellhop_pm_stage1_tier1.m` | Fixed-PM Tier-1 now defaults to `X` and uses X-specific case roots. |
| `scripts/validation/validate_pe_bellhop_pm_frequency_extension.m` | Frequency-extension cases now default to `X` and use source-tagged case roots. |
| `scripts/validation/validate_pe_bellhop_pm_ensemble.m` | Ensemble flat/rough cases now default to `X`; source geometry is included in case names and request fingerprints, preventing reuse of R caches. |
| `scripts/validation/validate_pe_bellhop_pm_model_discrepancy_statistics.m` | The current Stage-4 main program explicitly passes `source_geometry=X` to its ensemble engine. |
| `scripts/validation/validate_bellhop_source_geometry_rx_audit.m` | Retained validation-only R/X single-variable audit. |

## Incident-field regression after integration

The 4 kHz incident-plane audit was rerun at the same grid, step, Gaussian
`.sbp`, receivers and 5001/10001 beam counts. Bellhop `.prt` identifies both
active cases as `Line source (Cartesian coordinates)`.

| metric | pre-integration R | integrated X |
|---|---:|---:|
| PE--Bellhop TL P95 over M95 (dB) | 0.15794309 | 0.0019036364 |
| PE--Bellhop normalized L2 over M99 | 0.0098353183 | 0.0061048459 |
| phase RMS over M99 (rad) | 0.0061048392 | 0.0061048329 |
| 5001--10001 beam L2 | 5.2531469e-6 | 5.2520122e-6 |
| hard gates | FAIL (one amplitude gate) | PASS (13/13) |

PE--AS remains at `3.73e-13` normalized L2. No threshold was relaxed and no
data-dependent scalar was fitted. The improvement therefore isolates the
Bellhop point/line-source angular-amplitude convention.

## Provenance limit

The existing 4/6/8 kHz, 8-seed and 50-seed result reports were generated before
this code-default migration and retain their historical R-source provenance.
They were not silently relabeled or numerically rewritten. The already executed
seed-260001 X Tier-1 audit changed delta TL only from `0.3104335 dB` to
`0.3103707 dB` and left phase at `-2.2482007 rad`, so the reflection-model
discrepancy interpretation is unchanged. A future ensemble rerun will use X,
fresh source-tagged case roots, and fingerprints that reject old R caches.

## Archive cleanup

Superseded design material, the pre-X incident FAIL run, one explicit Stage-4
intermediate report, and large raw R/X audit work directories were moved to the
recoverable `cash/` quarantine. Final reports and machine-readable summary
artifacts remain active. See
`reports/pe_bellhop_source_geometry_integration_archive_manifest.md`.
