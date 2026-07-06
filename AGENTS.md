# AGENTS.md

## Scope
This repository is a MATLAB vertical underwater acoustic channel and MPSK communication project.
Work from code first. Do not assume external documents are authoritative unless explicitly requested.

The long-term project context is a seabed-to-near-surface vertical underwater acoustic communication channel. The model uses PE/WAPE-style propagation to obtain channel frequency responses and communication metrics.

## MATLAB Style
- Keep functions deterministic for the same input parameters and RNG seeds.
- Preserve current naming style: `paramsV` for user input, `cfg` for validated runtime config, `output/results/meta` for structured outputs.
- Prefer explicit struct fields over positional arguments for new physics options.
- Keep complex arrays and frequency-domain variables clearly named with `_f`, time/baseband taps with `_bb` or `_taps`.
- Use clear suffixes when helpful: `_xy` for transverse spatial fields and `_k` for transverse wavenumber-domain fields.
- Do not silently change units. Depth is meters, frequency is Hz, sound speed is m/s, z is positive downward.
- Avoid introducing new toolbox dependencies beyond standard MATLAB unless explicitly approved.

## Coordinate Convention
- Sea surface is `z=0`.
- Positive `z` points downward.
- Upward propagation means marching from larger `z_tx` to smaller `z_rx`.
- Receiver depth must satisfy `0 <= z_rx < z_tx`.

## Files That Require Extra Caution
- `vertical_channel_model.m`: public channel API and config validation boundary.
- `vertical_wape_propagator.m`: propagation core and frequency loop.
- `pm_surface_boundary_model.m`: rough-surface reflection model.
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

## Development Workflow
- Prefer small, reviewable changes.
- For new physics, first add disabled-by-default configuration and regression checks.
- Then add the minimal model.
- Then add advanced models only after the disabled/default path is confirmed unchanged.
- After each stage, summarize changed files, changed interfaces, validation results, and documentation updates.
- Before closing the task, verify that any required updates to `PROJECT_CONTEXT.md`, `vertical_comm_guide.md`, README files, and validation reports have been made.

## Do Not Do
- Do not rewrite the whole project unless explicitly asked.
- Do not convert scripts into classes or packages unless explicitly asked.
- Do not change public script names or function signatures unless required by the task.
- Do not remove existing plotting or saving behavior unless explicitly asked.
- Do not replace the current WAPE implementation with a different propagation method.
