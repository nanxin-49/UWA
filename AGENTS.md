# AGENTS.md

## Scope
This repository is a MATLAB vertical underwater acoustic channel and MPSK communication project.
Work from code first. Do not assume external documents are authoritative unless explicitly requested.
See `PROJECT_CONTEXT.md` for project background and current status.

## Retired Local Artifacts: cash/
- The repository-root `cash/` directory is a recoverable quarantine for obsolete or superseded experiments, data, figures and reports. The name is deliberately `cash`, not `cache`.
- Its contents are NOT current project context, validation evidence for the active environment, source dependencies, or a result cache to reuse automatically.
- Do not read, search, enumerate, index, summarize, execute or restore anything inside `cash/` unless the user explicitly requests that archive operation. During an explicitly authorized archive move, path/existence checks for the exact target are allowed; do not subsequently inspect archived contents.
- Normal repository searches, recursive listings and automatic result discovery must exclude `cash/`. When bypassing Git ignore rules, explicitly use an exclusion such as `rg -g '!cash/**'` or an equivalent scoped search.
- `.gitignore` excludes `/cash/` after all MATLAB/Markdown allow-rules. Never stage or commit its contents, including with `git add -f`, without an explicit user request to reverse this policy.
- Keep only a small archive manifest and supersession notice OUTSIDE `cash/` when needed for navigation. Do not recover old results to fill missing active results.
- The 500 m rough-surface Bellhop audit and its figures are retired. The user subsequently approved a physical source clearance of 1 m: water depth 100 m, Tx depth 99 m, Rx depth 3 m. Use `generate_bellhop_100m_source1m_visuals_vertical` for the current representative figures. The earlier zero-clearance epsilon studies remain historical diagnostics, not a gate for the newly approved installation. Read `reports/bellhop_100m_rebaseline_plan.md`, not the archive. Old 500 m runners/plotters remain reference code only and must not be run with their historical defaults.

## Repository Layout and New File Placement
- Classify new files by responsibility and workflow stage, not merely by file extension or algorithm name. Inspect the existing neighboring files and the relevant `README.md` before choosing a location.
- Keep the repository root small. Place a new file at the root only when it is a stable public API wrapper, a required compatibility entrypoint, project setup/navigation, or a top-level authoritative document. Do not put ordinary implementations, validation scripts, figures, logs, MAT files, or convenience experiments at the root.
- `vertical_channel_model.m`, `main_vertical.m`, `explain_main_vertical.m`, and `comm_main_vertical_psk.m` are stable root-level compatibility/public entrypoints. Their reusable or substantive implementations belong under `src/` or `examples/` as described below.

### Reusable implementation: `src/`
- `src/` contains reusable production implementation only. Do not place demos, one-off validation scripts, report builders, or generated results there.
- Place new reusable code in the module matching its primary responsibility:
  - `src/channel/`: public channel implementation, channel-component assembly, carrier-phase reference, and CIR representation conversion.
  - `src/propagation/`: PE/WAPE marching, propagation operators, and source initial fields.
  - `src/surface/`: sea-surface spectra, realizations, boundary/reflection/scattering models, and surface-specific operators.
  - `src/receiver/`: forward/cached receiver propagation, adjoint projection, receiver contraction, and receiver-side execution.
  - `src/statistics/`: statistical model construction, estimation, conditional/empirical models, and sampling.
  - `src/communication/`: communication taps, modulation/demodulation, noise, and communication-ensemble evaluation.
  - `src/bubble/`: bubble spectra, bubble environments, and effective-medium models.
- When a reusable feature spans modules, place the implementation in the module that owns its output/abstraction and keep dependencies consistent with `src/README.md`; do not duplicate it across script folders.
- If a genuinely new `src/` module is added, update `src/README.md` and `setup_vertical_project.m` so the module and dependency direction remain explicit. Do not solve path problems by changing the MATLAB current working directory.

### Demos and workflow scripts
- Put the substantive content of user-facing demonstrations in `examples/`. Preserve an existing root compatibility wrapper when its public name must remain stable.
- `scripts/` contains runnable research, validation, and reporting workflow entrypoints. Classify a new script by what the run is intended to do:
  - `scripts/validation/`: smoke tests, invariants, regressions, diagnostics, audits, acceptance checks, and communication-chain validation.
  - `scripts/validation/support/`: helpers used only by validation workflows, including external-tool readers/writers/runners, independent references, artifact assertions, and validation metadata. These helpers are not public APIs.
  - `scripts/comparisons/`: controlled comparisons between models, methods, approximations, or coherent/incoherent alternatives.
  - `scripts/experiments/`: Monte Carlo studies, parameter sweeps, sensitivity studies, and calibration runs.
  - `scripts/reporting/`: figure generation, plotting from saved results, visualization, summary-table generation, and report builders.
- Choose the folder by the primary purpose of the entrypoint. For example, a script that runs an acceptance test and happens to plot diagnostics belongs in `validation/`; a script that only renders accepted saved data belongs in `reporting/`.
- New script entrypoints must be runnable from any MATLAB current directory after `setup_vertical_project`. Resolve an explicit project root and output directory; do not rely on `cd`, the caller's current directory, or outputting beside the script.
- Reusable computational logic discovered while writing a workflow script belongs under the appropriate `src/` module. Keep only orchestration, case configuration, assertions, and workflow-specific presentation in `scripts/`.

### Documents, results, and temporary files
- Put current one-off validation reports, experiment reports, audit conclusions, and project-status reports in `reports/`.
- Keep the current authoritative technical and project navigation documents at their established locations: `README.md`, `PROJECT_CONTEXT.md`, `vertical_comm_guide.md`, `AGENTS.md`, `src/README.md`, and `scripts/README.md`.
- Put non-authoritative historical specifications and superseded explanatory documents in `docs/history/`; do not use them as the sole basis for current APIs, defaults, or validation status.
- Put generated artifacts under a case-specific subdirectory of `results/`, rather than at the repository root or beside source code:
  - `results/validation/<case>/` for validation data, checkpoints, tables, and validation figures.
  - `results/visualization/<case>/` for standalone visualization products and atlases.
  - `results/experiments/<case>/` for experiment, sweep, Monte Carlo, sensitivity, and calibration outputs.
  - `results/comparisons/<case>/` for comparison-study outputs.
  - `results/communication/<case>/` and `results/channel_demo/<case>/` for communication and channel-demo runs when those categories apply.
- Use `results/archive/` only for an explicitly defined, recoverable workflow archive. It is distinct from the retired local quarantine `cash/` and must not be used as an unstructured dumping ground.
- `old/` and `tmp/` are not approved destinations for new persistent project files. Use `tmp/` only for disposable intermediate work when a task explicitly needs it, and move any retained artifact to its proper documented location before completing the task. Do not place new files in `old/` without an explicit repository policy or user instruction.
- Do not create logs, generated MAT/CSV/media files, downloaded references, or ad hoc notes at the repository root. Use the appropriate `results/<category>/<case>/`, `reports/`, or temporary location.
- When adding a new entrypoint, reusable module, report family, or result category, update the relevant navigation document (`README.md`, `src/README.md`, or `scripts/README.md`) so future work can discover its location and invocation.

## Development Standards

### MATLAB style
- Use MATLAB for programming tasks unless the user explicitly requests another language or an existing project integration requires a different language or tool.
- Keep functions deterministic for the same input parameters and RNG seeds.
- Preserve current naming style: `paramsV` for user input, `cfg` for validated runtime config, `output/results/meta` for structured outputs.
- Prefer explicit struct fields over positional arguments for new physics options.
- Keep complex arrays and frequency-domain variables clearly named with `_f`, time/baseband taps with `_bb` or `_taps`.
- Use clear suffixes when helpful: `_xy` for transverse spatial fields and `_k` for transverse wavenumber-domain fields.
- Do not silently change units. Depth is meters, frequency is Hz, sound speed is m/s, z is positive downward.
- Avoid introducing new toolbox dependencies beyond standard MATLAB unless explicitly approved.

### Development workflow
- Prefer small, reviewable changes.
- For new physics, first add disabled-by-default configuration and regression checks.
- Then add the minimal model.
- Then add advanced models only after the disabled/default path is confirmed unchanged.
- After each stage, summarize changed files, changed interfaces, validation results, and documentation updates.
- Before closing the task, verify that any required updates to `PROJECT_CONTEXT.md`, `vertical_comm_guide.md`, README files, and validation reports have been made.

## Coordinate Convention
- Sea surface is `z=0`.
- Positive `z` points downward.
- Upward propagation means marching from larger `z_tx` to smaller `z_rx`.
- Receiver depth must satisfy `0 <= z_rx < z_tx`.

## Files That Require Extra Caution
- `vertical_channel_model.m`: public compatibility API wrapper; validated implementation is `src/channel/vertical_channel_model_impl.m`.
- `src/propagation/vertical_wape_propagator.m`: propagation core and frequency loop.
- `src/surface/pm_surface_boundary_model.m`: rough-surface reflection model.
- `comm_main_vertical_psk.m`: communication-chain reference consumer of `H_f`.
- Markdown method notes or research-summary documents that describe the implemented physics.

Do not rename exported fields in `output` or `results` without updating all entry scripts and downstream consumers.

## Run Entrypoints
- Channel-only demo: `explain_main_vertical`
- End-to-end communication demo: `comm_main_vertical_psk`
- Reusable channel API: `vertical_channel_model(paramsV)`

## Regression Requirements
- Run at least one scalar-frequency channel case and one wideband communication case after nontrivial edits.
- Confirm `direct_only` and `direct_plus_reflect` scenarios both execute when relevant.
- Confirm no interface regression in `output.h_direct`, `output.h_reflect`, `output.h_total`, `output.H_f`, `output.f_axis`, `output.idx_f_ref`.
- For disabled-by-default features, confirm the disabled path matches the previous behavior within documented numerical tolerance.
- For new mathematically equivalent implementations, compare against the previous implementation and report the numerical difference.

## Testing Cost Control
- For quick validation, prefer reduced grids such as `nx=128` or `nx=256`, `ny=128` or `ny=256`, and `show_figures=false`.
- Do not run expensive `1024 x 1024` or large wideband tests repeatedly unless necessary.
- When reporting tests, include the exact parameter overrides used.

## Numerical Consistency
- Keep `H_f = H_direct_f + H_reflect_f` at every frequency bin.
- Keep `h_total = H_f(idx_f_ref)` and matching definitions for direct and reflected components.
- Preserve 1/R validation behavior in uniform medium unless intentionally changing the propagation model.
- Preserve deterministic sea surface and noise behavior for fixed seeds.
- Do not mix GPU/CPU or single/double paths in a way that changes results without documenting the expected tolerance.

## Output Compatibility
- Existing output fields must remain present even if a new feature is disabled.
- New fields should be added under `output.<feature>_meta`, `output.roughness_meta`, or `output.config`, not by replacing existing fields.
- Existing `results(ss).channel` structure in `comm_main_vertical_psk.m` must remain usable.

## Feature Integration Rules
- New physical parameters should enter through `paramsV` and be validated in `local_prepare_config`.
- If a new effect changes propagation physics, update both the direct-path loop and `local_march_field`, or refactor them together.
- If a new effect changes time variation or Doppler, define how it maps to both `fd_hz_used` and `H_f/H_baseband`.
- If a new reflection model is added, preserve compatibility with `enable_surface_reflection`, `surface_reflect_coeff`, and existing result fields where possible.
- Prefer additive extensions over breaking output-field changes.

## Documentation Discipline
- Every nontrivial physics or interface change must update the corresponding Markdown method note or research-summary document.
- Documentation updates must include: changed files, changed interfaces, formulas implemented or affected, assumptions, limitations, validation settings, validation results, and remaining issues.
- At the end of each completed implementation or validation task, check whether the main documentation must be updated before the final response.
- Keep `PROJECT_CONTEXT.md` aligned with the current project state, active models, public interfaces, validation status, and known limitations.
- Keep `vertical_comm_guide.md` aligned with the implemented physical/communication methods, equations, configuration semantics, validation interpretation, and result-reading guidance.
- If a task only changes generated outputs or runs an existing script without changing project behavior, documentation may be left unchanged, but the final response should say why no documentation update was needed.
- Do not claim a new physical model when the code only rewrites an existing model in an equivalent mathematical form.
- Do not claim fast statistical channel generation unless the code actually estimates or uses statistical parameters such as means, variances, covariances, distributions, or a stochastic generator.

## Do Not Do
- Do not rewrite the whole project unless explicitly asked.
- Do not convert scripts into classes or packages unless explicitly asked.
- Do not change public script names or function signatures unless required by the task.
- Do not remove existing plotting or saving behavior unless explicitly asked.
- Do not replace the current WAPE implementation with a different propagation method.
