# Visualization QA — 2026-10-08

- Runtime: shared MATLAB MCP, R2025b Update 6. No shell MATLAB launch.
- 11 PNG/vector-PDF pairs. Every PNG: 2400 × 1350 pixels (16:9), 150 DPI export;
  every PDF MediaBox: 1152 × 648 points (16 × 9 inches).
- Code Analyzer: no issues. `git diff --check`: passed (line-ending notice only).
- All eleven layouts reviewed visually. Fixed figure padding prevents tight-crop
  aspect changes; shortened field labels and scale ticks prevent crowding.
  Blue/orange solver lines use solid/dashed styles; BIE uses purple dotted lines.
  This is manual presentation QA, not an accessibility certification.
- Recomputed weighted native field errors from plotted MAT data:
  weak 2.8992586989933e-9; strong height 3.20783541734475e-5;
  canonical PM 2.47637565054783e-4 — agree with source tables.
- Rechecked dense maximum error 2.30669959065543e-4;
  maximum relative phase 1.74393972339186e-4 rad;
  maximum delay difference 7.01937218402365e-8 s;
  diagnostic reflected CIR complex L2 7.17177129079795e-5.
- Exact versus saved rank3 Incident whole-grid L2:
  weak 3.67175896717791e-9; PM 2.88025700055554e-4.
  These are input-compression differences, not independent PE propagation errors.
- No new numerical solver runs or timing measurements. Inputs remain in their
  original validation directories; source hashes in `source_manifest.csv`.
- Production physics/interfaces unchanged, so scalar-channel and wideband
  communication regression are outside this presentation-only task.
- PDFs requested with vector ContentType. Dimensions are inspected; text/font
  portability may still depend on the consuming presentation application.
