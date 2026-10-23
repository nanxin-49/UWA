# S2b evidence package navigation

Status: S2B_DATA_COMPLETE_REFERENCE_LIMITED. Goal remains incomplete; no automatic S3-S6.

Authoritative interpretation: repository reports/two_dimensional_channel_closeout_report.md, appendix S2b. This file is only a package map.

- data/p0_data.mat, line_metrics_both_weights.csv, line_weights.csv, p0_audit.json: 12 spatial primary + 12 legacy-weight rows, unchanged native 411-point fields.
- data/p1.mat, p1_gate.json, frozen_jobs.md, bie_density_inventory.csv, input_hashes.csv: numeric density snapshots, qualification, budgets and pre-first-job freeze.
- checkpoints/: every successful batch and conservative failure charge. No new BIE boundary solve.
- data/s2b_data.mat, fields.csv, coverage_after.csv, point_metrics.csv, array_metrics.csv: 108/108 real fields, 48 reused + 60 new, fixed nine coordinates, equal nine-point weights. Native and physically mapped fields are identical; no calibration applied.
- data/sensitivity_fields.csv (2250 complex samples), numerical.csv, reference.csv: all sensitivity fields, separate numerical changes, BIE known empirical components; missing medium/strong self-quadrature is NaN.
- data/plot_data.mat: self-contained field/metric input for pure plotting.
- figures_v2/: final F01-F07 editable FIG + PNG and captions. figures/ is retained first rendering.
- audit_v2/independent_audit.json: final independent Cartesian CSV and BIE sensitivity arithmetic, hash check, seven saved-MAT redraws with zero pixel difference, seven FIG reopen checks. audit/ is retained first audit.
- data/summary.json: first assembly snapshot; data/qualification.json: final qualification including FullFK window limitation. No physical PASS.
- provenance/code_hashes.csv + code/: final source snapshots. Acoustic execute/preparation/P0 code unchanged since successful corresponding runs; final audit/assembly later received formatting-only changes recorded by final Code Analyzer. Original historical BIE execution source hashes remain unknown, never backfilled.
- provenance/s2a_preservation_hashes.csv: all old S2a artifacts listed; final integrity check passed, no old file overwritten.
- logs/: execution and failure trace, final Code Analyzer, final cwd/data checks. Figure/I-O timing was not separately instrumented.

Resource accounting: new successful physical batches13.9342315s + all persisted S2a runs9.2915103s + conservative unmetered/failure charges150s =173.2257418s of7200s. Every batch had actual Windows physical Available Bytes checked beforehand; estimated incremental peak96MiB, not a measured peak.

Reproduction requires setup_vertical_project in the existing MATLAB MCP session. Pure render/audit functions accept absolute package paths and unused output folders. Reusing saved density for new points IS new physical computation and needs renewed P1 gates. No communications or production API changes; no commit/push.
