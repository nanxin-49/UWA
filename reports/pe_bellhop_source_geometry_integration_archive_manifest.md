# PE--Bellhop source-geometry integration archive manifest

Archive date: 2026-09-09

Recoverable archive root: `cash/pe_bellhop_source_geometry_integration_20260909/`

The following items were moved after the Bellhop line-source `X` convention
was integrated into the active PE--Bellhop validation entrypoints:

| archived item | reason |
|---|---|
| `reports/pe_bellhop_incident_field_comparison_design_audit.md` | Design-stage document superseded by the implemented incident-field validator and its result report. |
| `reports/pe_bellhop_incident_field_comparison_report_R_pre_X.md` | Pre-integration point-source `R` incident-field FAIL result; retained only for provenance. |
| `results/validation/pe_bellhop_incident_field_R_pre_X/` | Raw pre-integration `R` incident-field artifacts; replaced by the active `X` run at the original result path. |
| `results/validation/pe_bellhop_pm_model_discrepancy_statistics/stage3_engine_intermediate_report.md` | Internal Stage-3 engine report superseded by the final Stage-4 report and machine-readable outputs. |
| `results/validation/bellhop_source_geometry_rx_audit/cases/` | Large raw R/X Bellhop case workspace; final audit report, CSV summaries, MAT summary and Tier-1 metrics remain active. |
| `results/validation/bellhop_source_geometry_rx_audit/tier1_X/cases/` | Raw conditional Tier-1 Bellhop case workspace; final Tier-1 CSV/MAT/checks remain active. |
| `results/validation/bellhop_source_geometry_rx_audit/tier1_X/tier1_X_internal_report.md` | Internal generated report superseded by `reports/bellhop_source_geometry_rx_audit_report.md`. |

Current authoritative files remain outside `cash/`:

- `reports/pe_bellhop_incident_field_comparison_report.md`
- `reports/bellhop_source_geometry_rx_audit_report.md`
- `reports/pe_bellhop_pm_model_discrepancy_statistics_report.md`
- `results/validation/pe_bellhop_incident_field/`
- `results/validation/bellhop_source_geometry_rx_audit/*.csv`
- `results/validation/bellhop_source_geometry_rx_audit/bellhop_source_geometry_rx_audit.mat`
- `results/validation/bellhop_source_geometry_rx_audit/tier1_X/stage1a_tier1_*.csv`
- `results/validation/bellhop_source_geometry_rx_audit/tier1_X/stage1a_tier1_comparison.mat`

Nothing in the archive is an active cache or dependency. Do not restore or
reuse it automatically.
