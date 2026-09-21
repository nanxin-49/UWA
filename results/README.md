# Results Directory

Generated artifacts are organized here so the repository root stays readable.

- `validation/`: outputs from `scripts/validation`.
- `comparisons/`: outputs from `scripts/comparisons`.
- `experiments/`: sweeps, Monte Carlo outputs, calibration outputs, and derived experiment figures.
- `communication/`: end-to-end communication demo outputs such as `psk_comm_result.mat`.
- `visualization/`: channel-demo and wavefield visualization outputs.
- `reports/`: generated report assets and report-specific tables.
- `legacy_output/`: older tracked or historical output files kept for reference.

Most files in this tree are intentionally ignored by Git. Commit source code and Markdown summaries; keep heavy binary outputs local unless a specific review requires adding a small artifact.
