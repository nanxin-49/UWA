# PROJECT_CONTEXT.md

## Purpose
This file records the current code structure, execution paths, and data flow of the MATLAB vertical underwater acoustic channel and MPSK communication project.
It is intended as a code-first reference for future maintenance and feature work.

## Repository Role of Each Main MATLAB File

### `CARPE3D_vertical.m`
- Public channel entrypoint: `output = CARPE3D_vertical(paramsV)`.
- Accepts user-facing configuration in `paramsV`.
- Normalizes and validates parameters through `local_prepare_config`.
- Calls `propWAPE_vertical(cfg)` to perform propagation.
- Packages all outputs into a stable `output` struct.
- Enforces 1/R validation in uniform-medium mode when `enforce_1_over_R=true`.

### `propWAPE_vertical.m`
- Core upward-marching WAPE propagation engine.
- Builds transverse grids, source field, spectral operators, and absorption profile.
- Resolves receiver state, frequency axis, optional GPU path, and optional reflection path.
- Produces direct-path, reflected-path, and total complex channel responses:
  - `H_direct_f`
  - `H_reflect_f`
  - `H_f`
- Exposes reference-frequency scalar gains:
  - `h_direct`
  - `h_reflect`
  - `h_total`

### `pm_surface_kirchhoff_module.m`
- Rough-surface reflection submodule.
- Synthesizes a 2D Pierson-Moskowitz rough sea surface.
- Applies Kirchhoff phase distortion to the incident surface field.
- Supports:
  - reflection coefficient control
  - normal-incidence phase mode
  - oblique phase mode using TX/RX geometry
- Returns reflected field and reflection metadata.

### `explain_main_vertical.m`
- Channel-only demonstration script.
- Builds a scalar-frequency `paramsV`.
- Calls `CARPE3D_vertical(paramsV)`.
- Saves the output struct and standard figures.

### `comm_main_vertical_psk.m`
- End-to-end communication demonstration script.
- Builds a wideband channel configuration.
- Calls `CARPE3D_vertical(paramsV)` for each scenario.
- Converts acoustic frequency response to a discrete baseband channel.
- Runs MPSK modulation, channel convolution, noise injection, equalization, and BER/SER statistics.

### `modem_psk.m`
- Utility wrapper for M-PSK modulation, demodulation, and error-rate measurement.
- Supported modes:
  - `modulate`
  - `demodulate`
  - `error_rate`

### `noise_inject_vertical.m`
- Utility wrapper for noise injection at complex baseband.
- Current built-in model: `awgn`.
- Also supports a custom extension hook via `custom_noise_fn`.

## Main Execution Paths

### Path 1: Channel Demo
`explain_main_vertical -> CARPE3D_vertical -> propWAPE_vertical -> output struct + saved figures`

### Path 2: Communication Demo
`comm_main_vertical_psk -> CARPE3D_vertical -> propWAPE_vertical -> H_f/f_axis -> baseband channel -> MPSK link simulation`

## Channel Data Flow

### Step 1: User parameters
- The user-facing configuration starts as `paramsV`.
- `paramsV` may contain both core physics fields and future extension fields.

### Step 2: Runtime configuration
- `CARPE3D_vertical` calls `local_prepare_config(paramsV)`.
- `local_prepare_config`:
  - merges defaults
  - copies user fields into `cfg`
  - keeps backward compatibility for some legacy names
  - validates geometry, grid, bandwidth, flags, and surface settings
  - derives:
    - `lambda0`
    - `dx`
    - `dy`
    - `dz_abs`
    - `path_span`
    - `numstep`
    - `dz_step`

### Step 3: Receiver state resolution
- `propWAPE_vertical` initializes `rx_state_used` from:
  - `x_rx`
  - `y_rx`
  - `z_rx`
- If `cfg.rx_position_fn` is not empty, it overrides the nominal RX location at `t=0`.
- `rx_state_used` then becomes the actual receiver state used by the propagation run.

### Step 4: Frequency axis resolution
- `local_resolve_frequency_axis(cfg, rx_state_used)` defines `f_axis`.
- Cases:
  - explicit sweep: if `cfg.f0` is already a vector
  - wideband auto-grid: if `enable_wideband=true`
  - scalar-frequency mode otherwise
- `idx_f_ref` selects the reference frequency used for scalar outputs and stored field slices.

### Step 5: Direct-path propagation
- For each frequency in `f_axis`, the code:
  - computes `lambda_f`, `k0`, and `dz_step_f`
  - constructs a split-step spectral propagator
  - marches the source field from `z_tx` to `z_rx`
  - samples the received field at `(ix_rx, iy_rx)`
- The result is stored as `H_direct_f(ifq)`.

### Step 6: Reflected-path propagation
- If `enable_surface_reflection=true`, the code performs:
  1. `tx -> surface` propagation via `local_march_field`
  2. rough-surface reflection via `pm_surface_kirchhoff_module`
  3. `surface -> rx` propagation via `local_march_field`
- The result is stored as `H_reflect_f(ifq)`.

### Step 7: Total frequency response
- At each frequency bin:
  - `H_f(ifq) = H_direct_f(ifq) + H_reflect_f(ifq)`
- After the frequency loop:
  - `h_direct = H_direct_f(idx_f_ref)`
  - `h_reflect = H_reflect_f(idx_f_ref)`
  - `h_total = H_f(idx_f_ref)`

### Step 8: Output packaging
- `CARPE3D_vertical` returns a stable `output` struct containing:
  - field snapshots
  - centerline diagnostics
  - rough-surface products
  - frequency-domain responses
  - reference-frequency scalar gains
  - receiver state and Doppler placeholder result
  - config echo and derived scalar diagnostics

## Current Communication Data Flow

### Step 1: Scenario construction
- `comm_main_vertical_psk.m` defines two scenarios:
  - `direct_only`
  - `direct_plus_reflect`
- The only scenario change is `enable_surface_reflection`.

### Step 2: Symbol generation
- Random bits are generated in the script.
- `modem_psk('modulate', bits_tx, M)` maps them to unit-power M-PSK symbols.

### Step 3: Channel acquisition
- Each scenario calls `channel = CARPE3D_vertical(paramsV)`.
- The communication chain mainly consumes:
  - `channel.f_axis`
  - `channel.H_f`
  - `channel.idx_f_ref`
  - `channel.h_direct`
  - `channel.h_reflect`
  - `channel.h_total`

### Step 4: Baseband channel construction
- `local_build_baseband_response`:
  - centers `H_f` around the reference frequency
  - builds a baseband frequency axis
  - interpolates the acoustic response onto that axis
- `local_build_channel_taps`:
  - converts the baseband frequency response to time-domain taps
  - truncates taps to a target cumulative energy ratio

### Step 5: Signal propagation at baseband
- The transmitted symbol stream is convolved with `h_bb` using:
  - `rx_clean = conv(tx_symbols, h_bb, 'same')`

### Step 6: Noise injection
- `noise_inject_vertical` adds noise to `rx_clean`.
- In the current demo:
  - model: `awgn`
  - control variable: `Eb/N0`
  - SNR conversion uses `bits_per_symbol`

### Step 7: Equalization and decisions
- `local_mmse_equalize` performs one-shot frequency-domain MMSE equalization using the known `h_bb`.
- `modem_psk('demodulate', rx_eq, M)` performs coherent hard-decision demodulation.
- `modem_psk('error_rate', bits_tx, bits_rx, M)` computes BER and SER.

### Step 8: Results packaging
- Each scenario stores:
  - physical channel scalars
  - acoustic frequency response
  - derived baseband response
  - BER/SER curves
  - effective SNR metadata
  - the full `channel` struct

## Important Interfaces and Extension Points

### `paramsV` / `cfg`
- New fields can be introduced through `paramsV`.
- `local_prepare_config` copies user fields into `cfg`.
- Existing validated fields must remain compatible with current callers.

### `rx_position_fn`
- Expected form:
  - `rx_xyz = rx_position_fn(t_s, state)`
- Current call site uses `t_s = 0`.
- Must return `[x_rx, y_rx, z_rx]`.

### `doppler_fn`
- Expected form:
  - `fd_hz = doppler_fn(t_s, tx_state, rx_state, env_state)`
- Current call site uses `t_s = 0`.
- Current result is recorded in `fd_hz_used` only.
- No Doppler compensation is currently applied to `H_f` or the communication chain.

### Sound-speed environment
- `env_mode='uniform'`: uses `c0`.
- `env_mode='layered'`:
  - first tries `cz_func(z)`
  - then tries `cz_table_z/cz_table_c`
  - otherwise uses a built-in fallback profile

### Surface reflection controls
- Main controls in `cfg`:
  - `enable_surface_reflection`
  - `sea_wind_speed`
  - `sea_hs_target`
  - `sea_seed`
  - `surface_reflect_coeff`
  - `surface_phase_mode`
  - `surface_oblique_clip`

### Runtime controls
- `save_mode`:
  - `rx_only`
  - `slice`
- `show_figures` toggles plotting and rough-surface visualization.
- `use_gpu` enables the GPU path only when a GPU is available.

## Current Output Structure

### Core field outputs
- `psifinal_xy`
- `psiout`
- `x`
- `y`
- `z_track`
- `Axz`
- `Ayz`

### Validation outputs
- `A_center`
- `R_center`
- `fit_slope`
- `fit_err_rms`
- `pass_1_over_R`
- `fit_mask`

### Reflection outputs
- `surface_elevation`
- `delta_phi`
- `psi_ref`
- `roughness_meta`

### Frequency-response outputs
- `f_axis`
- `H_direct_f`
- `H_reflect_f`
- `H_f`
- `idx_f_ref`

### Reference-frequency scalar outputs
- `h_direct`
- `h_reflect`
- `h_total`

### State and summary outputs
- `rx_state_used`
- `fd_hz_used`
- `rx_amplitude`
- `rx_phase_rad`
- `path_loss_db`
- `direct_to_reflect_db`
- `phase_diff_rad`
- `config`

## Known Structural Sensitivities

### Validation boundary
- `local_prepare_config` is the main contract boundary.
- New models often require new fields here and may also trigger existing constraints.

### Dual propagation logic
- There are two propagation implementations:
  - the main direct-path frequency loop
  - `local_march_field`
- Any physics change that affects marching should be applied consistently to both paths.

### Reference-frequency capture behavior
- Some outputs are only retained for `idx_f_ref`.
- This includes field slices and rough-surface products.
- Do not assume those stored arrays represent all frequencies.

### Communication coupling
- `comm_main_vertical_psk.m` depends on the meaning and consistency of:
  - `H_f`
  - `f_axis`
  - `idx_f_ref`
  - `h_total`
- Breaking these semantics will break the communication demo even if the channel demo still runs.
