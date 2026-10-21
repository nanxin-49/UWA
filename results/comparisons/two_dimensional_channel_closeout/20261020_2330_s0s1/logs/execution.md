# S0-S1 execution record

Run label: `20261020_2330_s0s1` (explicit unique directory label; not a measurement timestamp).

MATLAB MCP existing session: `25.2.0.3312555 (R2025b) Update 6`. Initial and final MATLAB current directory: `E:\MISC\BELLHOP\AcousticsToolbox_2020\tests\Munk`. The workflow used absolute project paths and did not change that directory.

Commands evaluated through MATLAB MCP after `addpath('E:\MISC\CARPE3D_matlab\Explain'); setup_vertical_project;`:

```matlab
checks = test_two_dimensional_closeout_metrics();
r = run_two_dimensional_closeout_s0_s1('20261020_2330_s0s1');
a = audit_two_dimensional_closeout_s0_s1(r.run_dir);
```

All five new `.m` files were checked through MATLAB MCP Code Analyzer; each returned zero issues. Formula controls passed 9/9. Historical full-precision rebuild passed 9/9. Independent CSV audit passed. Seven figure PNGs were redrawn from `plot_data/plot_data.mat` into `audit/redraw/`; all seven matched the first rendering pixel for pixel. A human visually inspected F01, F03, F04, F05, F05_grid, F06 and F07 at the delivered PNG size. PNGs use `print -dpng -r140` with fixed figure dimensions; PDFs use `exportgraphics(...,'ContentType','vector')`.

Observed postprocessing wall time: 66.2235 s. New PE propagation, angular spectrum field evaluation, BIE boundary solution or density field evaluation, FullFK integration, and Bellhop jobs: all zero. Solver wall time: zero. Peak memory: not measured.

Earlier diagnostic runs `20261020_2100_s0s1` and `20261020_2145_s0s1` are preserved. The first stopped at an overly strict pixel test because `exportgraphics` auto-cropped two PNGs to different dimensions. Fixed-canvas PNG export produced repeatable pixels. A first audit invocation also exposed a numeric CSV validity-index conversion issue; the audit was corrected and passed before this run. Those diagnostic directories are not the final evidence set.

No production solver, public API, or historical result was changed; no commit or push was made. The original Git worktree was clean before this task.
