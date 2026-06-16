# PROJECT_CONTEXT.md

## Purpose
This file records the current code structure, execution paths, and data flow of the MATLAB vertical underwater acoustic channel and MPSK communication project.
It is intended as a code-first reference for future maintenance and feature work.

## Repository Role of Each Main MATLAB File

### `vertical_channel_model.m`
- Public channel entrypoint: `output = vertical_channel_model(paramsV)`.
- Accepts user-facing configuration in `paramsV`.
- Normalizes and validates parameters through `local_prepare_config`.
- Calls `vertical_wape_propagator(cfg)` to perform propagation.
- Packages all outputs into a stable `output` struct.
- Enforces 1/R validation in uniform-medium mode when `enforce_1_over_R=true`.

### `vertical_wape_propagator.m`
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

### `pm_surface_boundary_model.m`
- Rough-surface reflection submodule.
- Synthesizes a 2D Pierson-Moskowitz rough sea surface.
- Applies Kirchhoff phase distortion to the incident surface field through a selectable boundary interface.
- Supports:
  - reflection coefficient control
  - normal-incidence phase mode
  - oblique phase mode using TX/RX geometry
  - spatial-domain and implicit wavenumber-domain Kirchhoff phase-screen boundary models
- Returns reflected field and reflection metadata.

### `explain_main_vertical.m`
- Channel-only demonstration script.
- Builds a scalar-frequency `paramsV`.
- Calls `vertical_channel_model(paramsV)`.
- Saves the output struct and standard figures.

### `comm_main_vertical_psk.m`
- End-to-end communication demonstration script.
- Builds a wideband channel configuration.
- Calls `vertical_channel_model(paramsV)` for each scenario.
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
`explain_main_vertical -> vertical_channel_model -> vertical_wape_propagator -> output struct + saved figures`

### Path 2: Communication Demo
`comm_main_vertical_psk -> vertical_channel_model -> vertical_wape_propagator -> H_f/f_axis -> baseband channel -> MPSK link simulation`

## Channel Data Flow

### Step 1: User parameters
- The user-facing configuration starts as `paramsV`.
- `paramsV` may contain both core physics fields and future extension fields.

### Step 2: Runtime configuration
- `vertical_channel_model` calls `local_prepare_config(paramsV)`.
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
- `vertical_wape_propagator` initializes `rx_state_used` from:
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
  2. rough-surface reflection via `pm_surface_boundary_model`
  3. `surface -> rx` propagation via `local_march_field`
- The result is stored as `H_reflect_f(ifq)`.
- The selected surface boundary model acts only in step 2; it does not change the two propagation segments.

### Step 7: Total frequency response
- At each frequency bin:
  - `H_f(ifq) = H_direct_f(ifq) + H_reflect_f(ifq)`
- After the frequency loop:
  - `h_direct = H_direct_f(idx_f_ref)`
  - `h_reflect = H_reflect_f(idx_f_ref)`
  - `h_total = H_f(idx_f_ref)`

### Step 8: Output packaging
- `vertical_channel_model` returns a stable `output` struct containing:
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
- Each scenario calls `channel = vertical_channel_model(paramsV)`.
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
  - `surface_boundary_model`
  - `surface_boundary_check_equivalence`
  - `surface_boundary_debug`
  - `surface_boundary_equivalence_tol`
  - `surface_boundary_coupling_diagnostics`
  - `surface_boundary_coupling_debug`
  - `surface_boundary_redistribution_diagnostics`
  - `surface_boundary_redistribution_debug`

### Kirchhoff surface boundary interface
- Default behavior is unchanged:
  - `surface_boundary_model='kirchhoff_spatial'`
  - `G_xy = surface_reflect_coeff * exp(1i*delta_phi)`
  - `psi_ref_xy = G_xy .* psi_inc_xy`
- Optional k-domain interface:
  - `surface_boundary_model='kirchhoff_kdomain'`
  - `Psi_inc_k = fft2(psi_inc_xy)`
  - `Psi_ref_k = fft2(G_xy .* ifft2(Psi_inc_k))`
  - `psi_ref_xy = ifft2(Psi_ref_k)`
- This is the discrete FFT product-convolution form of the same Kirchhoff phase screen. It represents an implicit periodic-grid operator with kernel proportional to `G_hat(K-K')`, but the code does not build a four-dimensional dense `B_xi(K,K')` matrix.
- FFT convention:
  - MATLAB `fft2` is unnormalized.
  - MATLAB `ifft2` includes `1/(nx*ny)`.
  - `kx/ky = [0:N/2-1,-N/2:-1]*2*pi/L`.
  - The implicit convolution is circular on the transverse periodic grid.
- Added `roughness_meta` fields:
  - `boundary_model`
  - `boundary_operator_form`
  - `boundary_dense_matrix_used`
  - `boundary_fft_convention`
  - `boundary_equivalence_error`
  - `boundary_coupling_diagnostics`
  - `boundary_coupling_debug`
  - `boundary_redistribution_diagnostics`
  - `boundary_redistribution_debug`
  - `boundary_debug_stats`
- When `enable_surface_reflection=false`, `roughness_meta.enabled=false` and the boundary metadata is still present with `boundary_operator_form='not_executed'`.
- Current limits:
  - This is not a T-matrix, SSA, or NLSSA rough-surface solver.
  - It is not a fast statistical channel generator.
  - It keeps the PM sea synthesis, phase factor, reflection coefficient, and two-segment reflected path unchanged.

### Boundary screen coupling diagnostics
- `surface_boundary_coupling_diagnostics=false` by default. When enabled, the code only adds metadata and does not change the reflected field or channel outputs.
- Diagnostics are evaluated on the retained reference frequency only. In wideband runs, `roughness_meta.boundary_coupling_diagnostics` therefore describes `f_axis(idx_f_ref)`.
- Boundary screen:
  - `G_xy = surface_reflect_coeff .* exp(1i*delta_phi)`
- Screen Fourier spectrum:
  - `G_hat_k = fft2(G_xy)`
  - `P_k = abs(G_hat_k).^2`
- Energy metrics:
  - `E_total = sum(P_k(:))`
  - `E_zero = P_k(1,1)`
  - `E_nonzero = E_total - E_zero`
  - `nonzero_power_fraction = E_nonzero / max(E_total, eps)`
- Coupling-radius metrics:
  - `DeltaK = sqrt(KX.^2 + KY.^2)`
  - `rms_delta_k = sqrt(sum(DeltaK(:).^2 .* P_k(:)) / max(E_total, eps))`
  - `rms_delta_k_nonzero` uses only nonzero wavenumber bins and normalizes by `E_nonzero`.
  - `energy_radius_50/90/95` are weighted cumulative radii over all spectrum bins.
  - `nonzero_energy_radius_90` is the 90 percent weighted radius after excluding the zero-wavenumber bin.
- Interpretation:
  - Larger `nonzero_power_fraction` and `rms_delta_k` indicate a less uniform Kirchhoff screen and stronger potential off-diagonal `K' -> K` coupling in the implicit `G_hat(K-K')` convolution.
  - These metrics are diagnostic summaries of the Kirchhoff phase-screen spectrum only. They are not strict boundary-operator matrix elements, scattering cross sections, T-matrix terms, SSA/NLSSA terms, or a statistical channel generator.
- `surface_boundary_coupling_debug=true` adds small summaries only:
  - 16-bin radial spectrum energy fractions.
  - top 8 nonzero spectrum peaks with `kx`, `ky`, `DeltaK`, and energy fraction.
  - Full `G_hat_k` and `P_k` arrays are not stored.

### Incident-weighted redistribution diagnostics
- `surface_boundary_redistribution_diagnostics=false` by default. When enabled, the code compares the current incident angular spectrum with the reflected angular spectrum generated by the same Kirchhoff phase screen.
- The diagnostic does not change `psi_ref`, `H_reflect_f`, `H_f`, or communication-chain inputs.
- Diagnostics are evaluated on the retained reference frequency only. The sweep script described below runs separate scalar-frequency cases when a frequency trend is needed.
- Spectra:
  - `Psi_inc_k = fft2(psi_inc_xy)`
  - `Psi_ref_k = fft2(psi_ref_xy)`
  - `Psi_flat_ref_k = fft2(surface_reflect_coeff .* psi_inc_xy)`
  - `P_k = abs(Psi_k).^2`
  - `DeltaK = sqrt(KX.^2 + KY.^2)`
- Spectrum moments:
  - `centroid_kx = sum(KX(:).*P_k(:)) / max(sum(P_k(:)), eps)`
  - `centroid_ky = sum(KY(:).*P_k(:)) / max(sum(P_k(:)), eps)`
  - `rms_delta_k = sqrt(sum(DeltaK(:).^2.*P_k(:)) / max(sum(P_k(:)), eps))`
  - `energy_radius_90` is the weighted cumulative `DeltaK` radius containing 90 percent of `P_k`.
- Reported redistribution metrics:
  - `incident_rms_delta_k_rad_per_m`
  - `reflect_rms_delta_k_rad_per_m`
  - `rms_delta_k_increase_rad_per_m`
  - `centroid_shift_kx/ky/mag_rad_per_m`
  - `reflect_high_k_fraction`, using `incident_energy_radius_90` as the high-k threshold
  - `flat_ref_rms_delta_k_rad_per_m`
  - `rough_vs_flat_rms_delta_k_increase_rad_per_m`
  - `rough_vs_flat_energy_radius_90_increase_rad_per_m`
  - `rough_vs_flat_high_k_fraction_increase`
- `surface_boundary_redistribution_debug=true` adds compact radial energy summaries and top 8 spectral peaks for incident, reflected, and flat-reflected spectra. Full `Psi_inc_k`, `Psi_ref_k`, and `Psi_flat_ref_k` arrays are not stored.
- Interpretation:
  - Larger reflected `rms_delta_k`, high-k fraction, and rough-vs-flat increase indicate stronger angular-spectrum broadening for the current incident field.
  - These are incident-weighted Kirchhoff phase-screen diagnostics only. They are not strict off-diagonal operator entries, T-matrix elements, SSA/NLSSA terms, scattering cross sections, or fast statistical channel parameters.

### Surface redistribution trend script
- `sweep_surface_boundary_redistribution_vertical.m` runs reduced scalar-frequency diagnostics with:
  - `nx=ny=128`, `xw=yw=50`, `z_tx=100`, `z_rx=3`, `sigma_src_m=0.4`
  - `show_figures=false`, `save_mode='rx_only'`, `enforce_1_over_R=false`
  - `surface_boundary_model='kirchhoff_kdomain'`
  - `surface_boundary_redistribution_diagnostics=true`
- It scans:
  - `sea_hs_target = [0.05, 0.5, 1.0]`
  - `sea_wind_speed = [3, 5, 8, 12]`
  - scalar `f0 = [4000, 6000, 8000]`
- It saves `sweep_surface_boundary_redistribution_vertical_result.mat` and lightweight PNG trend plots for RMS broadening and high-k fraction.
- Wind-speed sweep interpretation requires caution because the PM surface is rescaled to fixed `sea_hs_target=0.5`; changing wind speed changes the spectral shape after scaling, not a monotonic wave-height amplitude.

### Surface Monte Carlo channel statistics script
- `monte_carlo_surface_channel_vertical.m` runs reduced wideband Monte Carlo statistics over random Kirchhoff rough-surface realizations.
- The script fixes all physical and numerical parameters except `sea_seed`.
- Default reduced setup:
  - `mc_count=16`, overridable through environment variable `SURFACE_MC_COUNT`.
  - `seed_list=12345+(0:mc_count-1)`.
  - `nx=ny=128`, `xw=yw=50`, `z_tx=100`, `z_rx=3`, `sigma_src_m=0.4`.
  - `enable_wideband=true`, `f_band_hz=[4000 8000]`, `Nf_min=Nf_max=16`, `f_ref_hz=6000`.
  - `surface_boundary_model='kirchhoff_kdomain'`.
  - `surface_boundary_coupling_diagnostics=true`.
  - `surface_boundary_redistribution_diagnostics=true`.
  - `show_figures=false`, `save_mode='rx_only'`, `enforce_1_over_R=false`, `use_gpu=false`.
- Channel ensemble formulas:
  - `H_m(f)` is one Monte Carlo realization.
  - `mean_H(f) = (1/M) * sum_m H_m(f)`.
  - `var_H(f) = (1/M) * sum_m abs(H_m(f)-mean_H(f)).^2`.
  - Amplitude statistics are computed from `abs(H_m(f_ref))`.
  - Principal phase samples use `angle(H_m(f_ref))`.
  - Circular mean phase is `angle(mean(exp(1i*theta_m)))`.
  - Circular phase variance is `1 - abs(mean(exp(1i*theta_m)))`.
  - Reference reflection/direct ratio is `h_reflect_m / max_complex(h_direct_m, eps)`.
- Diagnostic ensemble statistics:
  - Mean, variance, standard deviation, and 5/25/50/75/95 percent quantiles are reported for selected C2/C2.5 metrics.
  - Coupling metrics include `nonzero_power_fraction` and `rms_delta_k_rad_per_m`.
  - Redistribution metrics include `reflect_rms_delta_k_rad_per_m`, `reflect_high_k_fraction`, and `rms_delta_k_increase_rad_per_m`.
- Tap metrics:
  - The script builds a baseband response by interpolating `H_f` around `idx_f_ref`.
  - It computes `h_bb = ifft(ifftshift(H_baseband))`, keeps taps until the requested energy ratio is reached, and reports compact delay/peak metrics.
  - Default tap controls are `symbol_rate_hz=1000`, `tap_fft_len=512`, and `tap_energy_ratio=0.999`.
- Outputs:
  - `monte_carlo_surface_channel_vertical_result.mat` with `params_base`, `seed_list`, `channels_summary`, `summary_table`, and `mc_stats`.
  - Compact numeric sample arrays for `H_f`, `H_direct_f`, `H_reflect_f`, and reference-frequency scalars are stored; full field arrays and full 2-D spectra are not saved.
  - PNG summaries: `H_f` mean magnitude envelope, `abs(H(f_ref))` histogram, phase histogram, reflected RMS delta-k histogram, high-k fraction histogram, and tap RMS delay histogram.
- Scope:
  - This is empirical Monte Carlo over the existing Kirchhoff phase-screen channel. It is not a new scattering formula, a closed-form stochastic channel generator, a BER/SER simulation, a T-matrix method, SSA/NLSSA, or a calibrated scattering cross-section model.

### Multi-sea-state Monte Carlo sweep script
- `sweep_monte_carlo_surface_channel_vertical.m` extends the C3 Monte Carlo workflow from one fixed sea state to an `Hs x wind` grid.
- Default reduced sweep:
  - `sea_hs_target = [0.05, 0.5, 1.0]`.
  - `sea_wind_speed = [3, 5, 8, 12]`.
  - `mc_count=16`, overridable through `SURFACE_MC_SWEEP_COUNT`.
  - `SURFACE_MC_SWEEP_MAX_CONDITIONS` can truncate the condition list for smoke tests.
  - `SURFACE_MC_SWEEP_RESULT_FILE` can redirect the `.mat` output for smoke tests.
  - `SURFACE_MC_SWEEP_FIGURE_PREFIX` can redirect PNG output prefixes for smoke tests.
  - Each condition uses `seed_list=12345+(0:mc_count-1)`.
  - Channel, boundary, diagnostics, and tap settings match `monte_carlo_surface_channel_vertical.m`.
- Sea-state ensemble formulas:
  - `H_{c,m}(f)` is the channel response for sea-state condition `c` and seed sample `m`.
  - `mean_H_c(f) = (1/M) * sum_m H_{c,m}(f)`.
  - `var_H_c(f) = (1/M) * sum_m abs(H_{c,m}(f)-mean_H_c(f)).^2`.
  - Reference-frequency amplitude, phase, reflection/direct ratio, C2/C2.5 diagnostic, and tap statistics are computed independently inside each sea-state condition.
- Outputs:
  - `sweep_monte_carlo_surface_channel_vertical_result.mat` with `base_params`, `sea_hs_values`, `sea_wind_values`, `seed_list`, `condition_specs`, `condition_results`, `run_summary_table`, `condition_summary_table`, and `sweep_stats`.
  - `condition_results(cc).mc_stats.H_f.samples` stores the compact `[Nf x mc_count]` frequency-response sample matrix. Full channel structs, spatial fields, and 2-D spectra are not saved.
  - `run_summary_table` has one row per `(Hs, wind, seed)`.
  - `condition_summary_table` has one row per `(Hs, wind)`.
  - `sweep_stats` stores trend matrices indexed as `[numel(sea_hs_values) x numel(sea_wind_values)]`.
- PNG summaries:
  - heatmaps for `abs_H_ref_mean`, `abs_H_ref_std`, `reflect_rms_delta_k_mean`, `reflect_high_k_fraction_mean`, `rms_delta_k_increase_mean`, and `tap_rms_delay_symbols_mean`.
  - line plots of `reflect_high_k_fraction_mean` and `reflect_rms_delta_k_mean` versus wind speed, grouped by `Hs`.
- Scope:
  - This is still empirical Monte Carlo over the implemented Kirchhoff phase-screen channel. It is not a closed-form statistical channel model, a calibrated sea-surface scattering law, a BER/SER simulation, or a new rough-surface scattering solver.
  - Wind-speed trends under fixed `sea_hs_target` should be interpreted with the same caution as the C2.5 sweep: PM spectral shape changes are rescaled to the requested wave height.

### Communication Monte Carlo PSK statistics script
- `monte_carlo_comm_surface_psk_vertical.m` runs reduced QPSK BER/SER Monte Carlo statistics over the current vertical channel and communication chain.
- Default sea conditions:
  - weak: `sea_hs_target=0.05`, `sea_wind_speed=5`.
  - strong: `sea_hs_target=1.0`, `sea_wind_speed=5`.
- Default scenarios:
  - `direct_only`: `enable_surface_reflection=false`.
  - `direct_plus_reflect`: `enable_surface_reflection=true`, `surface_boundary_model='kirchhoff_kdomain'`.
- Default communication setup:
  - `M=4`, `n_sym=1000`, `EbN0_dB_list=0:2:20`, `symbol_rate_hz=1000`.
  - `mc_count=4`, overridable through `SURFACE_COMM_MC_COUNT`.
  - `seed_list=12345+(0:mc_count-1)`.
  - fixed bit stream with `bits_seed=9000`.
  - varying AWGN seed: `noise_seed_base + condition_index*100000 + scenario_index*10000 + seed_index*100 + ebn0_index`.
- Data flow:
  - each run calls `vertical_channel_model(paramsV)`;
  - builds baseband taps from `H_f` with the same `ifft(ifftshift(...))` logic as `comm_main_vertical_psk.m`;
  - applies the explicit receive-window policy used by `comm_main_vertical_psk.m`;
  - injects AWGN with `noise_inject_vertical`;
  - equalizes with frequency-domain MMSE using known effective taps;
  - demodulates with `modem_psk` and stores BER/SER.
- D2 communication policy:
  - `receive_window_mode='peak_sync'` by default.
  - `peak_sync` finds the largest tap, discards pre-peak taps, applies `conv(...,'full')`, and keeps the first `N` received samples.
  - `ebn0_reference='rx_clean'` by default, so AWGN power is referenced to the received clean waveform instead of transmit-symbol power.
  - `same_legacy` remains available only as a local diagnostic mode in scripts that define it.
- Statistics:
  - `BER_{c,s,m}(e)` and `SER_{c,s,m}(e)` are error rates for sea condition `c`, scenario `s`, seed sample `m`, and Eb/N0 bin `e`.
  - `mean_BER_{c,s}(e) = (1/M) * sum_m BER_{c,s,m}(e)`.
  - `var_BER_{c,s}(e) = (1/M) * sum_m (BER_{c,s,m}(e)-mean_BER_{c,s}(e)).^2`.
  - The same mean, variance, standard deviation, and 5/25/50/75/95 percent quantiles are reported for SER and effective SNR.
- Outputs:
  - `monte_carlo_comm_surface_psk_vertical_result.mat` with `params_base`, `comm_cfg`, `sea_conditions`, `scenarios`, `seed_list`, `EbN0_dB_list`, `run_summary_table`, `curve_summary_table`, `curve_stats`, and compact BER/SER sample arrays.
  - PNG summaries for weak/strong direct-plus-reflect BER/SER and direct-only versus direct-plus-reflect comparisons.
- Scope:
  - This is an empirical reduced communication experiment using the current channel and current PSK chain. It is not a new channel model, BER theory, stochastic generator, or rough-surface scattering solver.
  - Zero BER at high Eb/N0 is stored as zero when it occurs; plot-only floors are used only to make semilog figures readable.

### Minimal communication closed-loop validation script
- `validate_comm_link_minimal_vertical.m` validates the communication chain in layers before interpreting PE channel BER/SER.
- Validation order:
  - ideal unit channel: `h_bb=1`;
  - single complex tap: `h_bb=0.7*exp(1i*0.8)`;
  - known short multipath: `h_bb=[1;0.35*exp(1i*0.5);0.15*exp(-1i*0.7)]`;
  - reduced PE direct-only `h_bb` current path and diagnostic peak-aligned path.
- The script uses:
  - `M=4`, `n_sym=2000`, `EbN0_dB_list=0:2:20`;
  - fixed `bits_seed=9000`;
  - deterministic AWGN seeds from `noise_seed_base=7000`.
- It compares two receive slicing modes:
  - `same`: current `conv(tx_symbols,h_bb,'same')` behavior used by the communication scripts.
  - `causal_head`: diagnostic-only `conv(...,'full')` followed by the first `N` samples, appropriate for causal taps whose main energy starts at tap 1.
  - `peak_sync`: D2 policy that shifts the effective tap origin to the dominant tap before causal-head receive-window selection.
- Outputs:
  - `validate_comm_link_minimal_vertical_result.mat` with `comm_cfg`, `test_cases`, `summary_table`, `case_results`, and `pe_alignment`.
  - PNG curves for simple-channel BER/SER and PE tap magnitudes before/after diagnostic peak alignment.
- Scope:
  - This script is diagnostic only. It does not change `comm_main_vertical_psk.m`, `monte_carlo_comm_surface_psk_vertical.m`, or PE propagation outputs.

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

## 2026-06-10 Kirchhoff k-domain Boundary Validation

Changed files for this interface:
- `vertical_channel_model.m`: adds user-facing `paramsV` fields and validation.
- `vertical_wape_propagator.m`: forwards the surface boundary configuration and preserves disabled-path metadata.
- `pm_surface_boundary_model.m`: implements spatial and implicit k-domain boundary application paths.
- `PROJECT_CONTEXT.md`: records formulas, interface, assumptions, limits, and validation.

Reduced-grid validation settings:
- Common channel setup: `f0=4000`, `c0=1500`, `z_tx=100`, `z_rx=3`, `xw=yw=50`, `nx=ny=128`, `sigma_src_m=0.4`, `stepz_lamb=0.5`, `sponge_ratio=0.12`, `alpha_max_np_per_m=0.15`, `env_mode='uniform'`, `show_figures=false`, `save_mode='rx_only'`, `use_gpu=false`, `sea_seed=12345`, `enforce_1_over_R=false`.
- Single-frequency cases compared default, explicit `kirchhoff_spatial`, explicit `kirchhoff_kdomain`, and direct-only.
- Wideband case used `enable_wideband=true`, `f_band_hz=[4000 8000]`, `Nf_min=Nf_max=16`, `f_ref_hz=6000`.

Validation results:
- Default vs explicit spatial `H_f` relative error: `0`.
- Single-frequency spatial vs kdomain:
  - `psi_ref` relative error: `4.753077448146337e-16`
  - `H_reflect_f` relative error: `6.785387727900559e-15`
  - `H_f` relative error: `5.456834559265761e-16`
  - `roughness_meta.boundary_equivalence_error.rel_l2`: `4.753077448146337e-16`
  - `roughness_meta.boundary_equivalence_error.max_abs`: `2.5033568429418485e-17`
  - equivalence check passed for tolerance `1e-10`.
- Direct-only regression:
  - `max(abs(H_reflect_f))`: `0`
  - `max(abs(H_f-(H_direct_f+H_reflect_f)))`: `0`
  - `roughness_meta.enabled=false`
- Wideband spatial vs kdomain:
  - spatial invariant `max(abs(H_f-(H_direct_f+H_reflect_f)))`: `0`
  - kdomain invariant `max(abs(H_f-(H_direct_f+H_reflect_f)))`: `0`
  - `H_f` relative error: `7.666749014695518e-16`
  - reference-frequency kdomain equivalence `rel_l2`: `4.827102620732236e-16`
  - reference-frequency kdomain equivalence `max_abs`: `2.964296775167217e-17`

Remaining issues:
- No dense boundary matrix is materialized, so matrix-level inspection of individual `B_xi(K,K')` entries is not available.
- The k-domain path intentionally remains mathematically equivalent to the old Kirchhoff phase-screen model and should not be interpreted as a higher-order rough-surface scattering model.

## 2026-06-10 Boundary Coupling Diagnostic Validation

Changed files for this diagnostic:
- `vertical_channel_model.m`: adds default-off diagnostic flags.
- `vertical_wape_propagator.m`: forwards diagnostic flags only for the reference-frequency reflection call and preserves disabled metadata.
- `pm_surface_boundary_model.m`: computes scalar coupling metrics from `G_xy` and optional compact debug summaries.
- `PROJECT_CONTEXT.md`: records formulas, interpretation, limits, and validation.

Reduced-grid validation settings:
- Common channel setup: `f0=4000`, `c0=1500`, `z_tx=100`, `z_rx=3`, `xw=yw=50`, `nx=ny=128`, `sigma_src_m=0.4`, `stepz_lamb=0.5`, `sponge_ratio=0.12`, `alpha_max_np_per_m=0.15`, `env_mode='uniform'`, `show_figures=false`, `save_mode='rx_only'`, `use_gpu=false`, `sea_seed=12345`, `surface_boundary_model='kirchhoff_kdomain'`, `enforce_1_over_R=false`.
- Diagnostic-off vs diagnostic-on used `sea_hs_target=0.5`.
- Weak roughness used `sea_hs_target=0.05`; strong roughness used `sea_hs_target=1.0`.
- Wideband case used `enable_wideband=true`, `f_band_hz=[4000 8000]`, `Nf_min=Nf_max=16`, `f_ref_hz=6000`.

Validation results:
- Diagnostic off/on consistency:
  - off `boundary_coupling_diagnostics.enabled=false`
  - on `boundary_coupling_diagnostics.enabled=true`
  - `H_f` relative error: `0`
  - `H_reflect_f` relative error: `0`
  - `psi_ref` relative error: `0`
- Weak vs strong roughness diagnostics:
  - weak `nonzero_power_fraction`: `0.50513062402654019`
  - strong `nonzero_power_fraction`: `0.99999482166633924`
  - weak `rms_delta_k_rad_per_m`: `0.71355283273946624`
  - strong `rms_delta_k_rad_per_m`: `6.5608130546379275`
  - weak `energy_radius_90_rad_per_m`: `0.95702627363155124`
  - strong `energy_radius_90_rad_per_m`: `8.9856140404778042`
  - weak `nonzero_energy_radius_90_rad_per_m`: `1.3119676091847634`
  - strong `nonzero_energy_radius_90_rad_per_m`: `8.9856140404778042`
- Direct-only regression:
  - `max(abs(H_reflect_f))`: `0`
  - `max(abs(H_f-(H_direct_f+H_reflect_f)))`: `0`
  - `boundary_coupling_diagnostics.enabled=false`
- Wideband regression:
  - `max(abs(H_f-(H_direct_f+H_reflect_f)))`: `0`
  - diagnostic metadata enabled at reference frequency.
  - `idx_f_ref=8`, `f_axis(idx_f_ref)=5866.666666666667 Hz`
- Debug metadata check:
  - radial bin count: `16`
  - top nonzero peak count: `8`
  - full `G_hat_k` stored: `false`
  - full `P_k` stored: `false`

Remaining issues:
- The metrics summarize the Fourier spectrum of the imposed Kirchhoff screen. They do not identify individual physically rigorous scattering mechanisms.
- Strong roughness can push the phase-screen approximation outside a conservative small-perturbation interpretation; treat the diagnostic trend as a numerical screen-spectrum trend, not as validated ocean-surface scattering physics.

## 2026-06-10 Incident-Weighted Redistribution Diagnostic Validation

Changed files for this diagnostic:
- `vertical_channel_model.m`: adds default-off redistribution diagnostic flags.
- `vertical_wape_propagator.m`: forwards redistribution flags only for the reference-frequency reflection call and preserves disabled metadata.
- `pm_surface_boundary_model.m`: compares `Psi_inc_k`, `Psi_ref_k`, and flat-reflected spectra without changing propagation results.
- `sweep_surface_boundary_redistribution_vertical.m`: reduced scalar trend sweep for `Hs`, wind speed, and frequency.
- `PROJECT_CONTEXT.md`: records formulas, interpretation, limits, and validation.

Reduced-grid validation settings:
- Common channel setup: `c0=1500`, `z_tx=100`, `z_rx=3`, `xw=yw=50`, `nx=ny=128`, `sigma_src_m=0.4`, `stepz_lamb=0.5`, `sponge_ratio=0.12`, `alpha_max_np_per_m=0.15`, `env_mode='uniform'`, `show_figures=false`, `save_mode='rx_only'`, `use_gpu=false`, `sea_seed=12345`, `surface_boundary_model='kirchhoff_kdomain'`, `enforce_1_over_R=false`.

Validation results:
- Redistribution off/on consistency at `f0=4000`, `sea_hs_target=0.5`, `sea_wind_speed=5`:
  - off `boundary_redistribution_diagnostics.enabled=false`
  - on `boundary_redistribution_diagnostics.enabled=true`
  - `H_f` relative error: `0`
  - `H_reflect_f` relative error: `0`
  - `psi_ref` relative error: `0`
- Weak vs strong roughness at `f0=4000`, `sea_wind_speed=5`:
  - weak `sea_hs_target=0.05`, `reflect_rms_delta_k=2.4092056625755345`
  - strong `sea_hs_target=1.0`, `reflect_rms_delta_k=6.552320692768185`
  - weak `rms_delta_k_increase=0.12555034900919093`
  - strong `rms_delta_k_increase=4.268665379201842`
  - weak `reflect_high_k_fraction=0.11866935742346503`
  - strong `reflect_high_k_fraction=0.8570406763058447`
  - weak `rough_vs_flat_rms_delta_k_increase=0.12555034900919093`
  - strong `rough_vs_flat_rms_delta_k_increase=4.268665379201842`
- Frequency trend at `sea_hs_target=0.5`, `sea_wind_speed=5`:
  - `f0=[4000,6000,8000] Hz`
  - `rms_delta_k_increase=[3.6659013500300062,4.0216592688428765,4.0046896856288567]`
  - `reflect_high_k_fraction=[0.7832240337445295,0.8228809206147447,0.8185641917866873]`
- Direct-only regression:
  - `max(abs(H_reflect_f))`: `0`
  - `max(abs(H_f-(H_direct_f+H_reflect_f)))`: `0`
  - `boundary_redistribution_diagnostics.enabled=false`
- Wideband regression:
  - `max(abs(H_f-(H_direct_f+H_reflect_f)))`: `0`
  - redistribution metadata enabled at reference frequency.
  - `idx_f_ref=8`, `f_axis(idx_f_ref)=5866.666666666667 Hz`
- Debug metadata check:
  - incident radial bin count: `16`
  - reflected top peak count: `8`
  - full `Psi_inc_k` stored: `false`
  - full `Psi_ref_k` stored: `false`
- Trend script `sweep_surface_boundary_redistribution_vertical` completed successfully. Representative rows:
  - Hs sweep broadening: `0.12555`, `3.6659`, `4.2687` for `Hs=0.05`, `0.5`, `1.0`
  - Wind sweep broadening at fixed `Hs=0.5`: `4.2995`, `3.6659`, `1.1948`, `0.67448` for wind `3`, `5`, `8`, `12 m/s`
  - Frequency sweep broadening: `3.6659`, `4.0217`, `4.0047` for `4000`, `6000`, `8000 Hz`

Remaining issues:
- The wind-speed trend is not monotonic under fixed `Hs_target` scaling; it reflects PM spectral-shape changes after amplitude normalization.
- The diagnostics characterize the simulated incident field and Kirchhoff screen only. They should not be used as calibrated ocean-surface scattering cross sections or strict non-diagonal boundary-operator matrix entries.

## 2026-06-10 Surface Monte Carlo Channel Statistics Validation

Changed files for this script:
- `monte_carlo_surface_channel_vertical.m`: independent reduced-grid Monte Carlo script over `sea_seed`.
- `PROJECT_CONTEXT.md`: records Monte Carlo scope, formulas, outputs, validation settings, and limits.

Reduced-grid Monte Carlo settings:
- Fixed channel setup: `c0=1500`, `z_tx=100`, `z_rx=3`, `xw=yw=50`, `nx=ny=128`, `sigma_src_m=0.4`, `stepz_lamb=0.5`, `sponge_ratio=0.12`, `alpha_max_np_per_m=0.15`, `env_mode='uniform'`, `show_figures=false`, `save_mode='rx_only'`, `use_gpu=false`, `surface_boundary_model='kirchhoff_kdomain'`, `enforce_1_over_R=false`.
- Wideband setup: `enable_wideband=true`, `f_band_hz=[4000 8000]`, `Nf_min=Nf_max=16`, `f_ref_hz=6000`.
- Sea-state setup: `sea_hs_target=0.5`, `sea_wind_speed=5`, `surface_boundary_coupling_diagnostics=true`, `surface_boundary_redistribution_diagnostics=true`.
- Tap setup: `symbol_rate_hz=1000`, `tap_fft_len=512`, `tap_energy_ratio=0.999`.
- Smoke test: `SURFACE_MC_COUNT=4`, `seed_list=12345:12348`.
- Default reduced run: `mc_count=16`, `seed_list=12345:12360`.

Validation results:
- Smoke test with 4 seeds completed successfully and produced finite `H_f`, diagnostic, and tap statistics.
- Default 16-seed run completed successfully.
- All runs shared the same frequency axis and `idx_f_ref`.
- Reference frequency:
  - `idx_f_ref=8`
  - `f_axis(idx_f_ref)=5866.666666666667 Hz`
- Channel invariant:
  - `max(abs(H_f-(H_direct_f+H_reflect_f)))=1.3985787785349405e-17`
- Direct-path consistency across seeds:
  - `direct_drift_max_abs=1.3904866443919908e-17`
- Reference-frequency total channel statistics:
  - `mean(abs(H(f_ref)))=0.04451535630826136`
  - `std(abs(H(f_ref)))=0.019683054200486178`
  - `abs(H(f_ref))` quantiles `[5,25,50,75,95]% = [0.012448453947227933, 0.029477650012892355, 0.047715895428688343, 0.05917630294551255, 0.068142386535915839]`
  - circular mean phase `=-1.5721223443847119 rad`
  - circular phase variance `=0.42130134481217185`
- Reflection/direct ratio statistics at `f_ref`:
  - `mean(abs(h_reflect/h_direct))=0.83366541472645794`
  - `std(abs(h_reflect/h_direct))=0.357222192582143`
  - circular mean phase `=-2.447845029934725 rad`
  - circular phase variance `=0.832299121295223`
- Coupling diagnostic statistics:
  - `nonzero_power_fraction` mean/std `=0.99990160177794407 / 6.5681746413696337e-05`
  - `rms_delta_k_rad_per_m` mean/std `=6.4634066111591792 / 0.052390699662500283`
- Redistribution diagnostic statistics:
  - `reflect_rms_delta_k_rad_per_m` mean/std `=6.4729622076609932 / 0.0676218962238834`
  - `reflect_high_k_fraction` mean/std `=0.81498689692969417 / 0.0092872966143607676`
  - `rms_delta_k_increase_rad_per_m` mean/std `=3.9946260540001335 / 0.0676218962238834`
- Tap metric statistics:
  - `tap_count` mean/std `=512 / 0`
  - `tap_energy_kept` mean `=1`
  - `tap_rms_delay_symbols` mean/std `=99.835328110577819 / 69.550813743973464`
  - `tap_peak_fraction` mean/std `=0.87083037958185672 / 0.18399224840340037`
- Generated result artifacts:
  - `monte_carlo_surface_channel_vertical_result.mat`
  - `monte_carlo_surface_channel_vertical_H_f_mean_magnitude.png`
  - `monte_carlo_surface_channel_vertical_abs_h_ref_hist.png`
  - `monte_carlo_surface_channel_vertical_phase_h_ref_hist.png`
  - `monte_carlo_surface_channel_vertical_reflect_rms_delta_k_hist.png`
  - `monte_carlo_surface_channel_vertical_reflect_high_k_fraction_hist.png`
  - `monte_carlo_surface_channel_vertical_tap_rms_delay_hist.png`

Remaining issues:
- `mc_count=16` is a reduced validation ensemble, not a converged ocean-channel distribution.
- The only Monte Carlo variable is `sea_seed`; sea-state, geometry, frequency grid, and solver settings are fixed.
- The script reports tap metrics only. It does not run BER/SER or modify `comm_main_vertical_psk.m`.
- The statistics are empirical summaries of the existing Kirchhoff phase-screen channel and should not be treated as a closed-form statistical channel generator.

## 2026-06-10 Multi-Sea-State Monte Carlo Sweep Validation

Changed files for this script:
- `sweep_monte_carlo_surface_channel_vertical.m`: independent reduced-grid `Hs x wind x seed` Monte Carlo sweep.
- `PROJECT_CONTEXT.md`: records C3.5 sweep scope, formulas, outputs, validation settings, and limits.

Reduced-grid sweep settings:
- Fixed channel setup: `c0=1500`, `z_tx=100`, `z_rx=3`, `xw=yw=50`, `nx=ny=128`, `sigma_src_m=0.4`, `stepz_lamb=0.5`, `sponge_ratio=0.12`, `alpha_max_np_per_m=0.15`, `env_mode='uniform'`, `show_figures=false`, `save_mode='rx_only'`, `use_gpu=false`, `surface_boundary_model='kirchhoff_kdomain'`, `enforce_1_over_R=false`.
- Wideband setup: `enable_wideband=true`, `f_band_hz=[4000 8000]`, `Nf_min=Nf_max=16`, `f_ref_hz=6000`.
- Diagnostics: `surface_boundary_coupling_diagnostics=true`, `surface_boundary_redistribution_diagnostics=true`.
- Sweep grid: `sea_hs_target=[0.05, 0.5, 1.0]`, `sea_wind_speed=[3, 5, 8, 12]`.
- Tap setup: `symbol_rate_hz=1000`, `tap_fft_len=512`, `tap_energy_ratio=0.999`.
- Smoke test: `SURFACE_MC_SWEEP_COUNT=2`, `SURFACE_MC_SWEEP_MAX_CONDITIONS=2`.
- Default reduced run after the C4 follow-up default-count increase: `12` sea-state conditions, `mc_count=16`, `192` total propagations, `seed_list=12345:12360`.
- The previously recorded validation matrices below came from the earlier cost-control run with `mc_count=4`, `48` total propagations, and `seed_list=12345:12348`.

Validation results:
- Smoke test completed successfully and produced finite `.mat`, table, and PNG outputs.
- Default 48-run sweep completed successfully.
- All condition-level frequency axes and `idx_f_ref` checks passed.
- `sweep_stats` matrix dimensions are `[3 x 4]`, matching `[numel(sea_hs_values) x numel(sea_wind_values)]`.
- Reference frequency:
  - `idx_f_ref=8`
  - `f_axis(idx_f_ref)=5866.666666666667 Hz`
- Channel invariant:
  - `max(abs(H_f-(H_direct_f+H_reflect_f)))=1.3877787807814457e-17`
- Direct-path consistency:
  - maximum condition-level `direct_drift_max_abs=0`
- Summary table finite-field check:
  - all selected channel, phase, ratio, diagnostic, and tap statistics are finite.
- `abs_H_ref_mean` matrix, with rows `Hs=[0.05,0.5,1.0]` and columns `wind=[3,5,8,12]`:
  - `[0.039000343834070653, 0.03775446157444675, 0.021052540344231319, 0.01646692716410688]`
  - `[0.027374350743730395, 0.038604342599260662, 0.055552395529549083, 0.031575982180233053]`
  - `[0.03928118735487425, 0.047635492235929437, 0.058271523033646792, 0.068788019597050512]`
- `reflect_rms_delta_k_mean` matrix:
  - `[3.3081145971986365, 2.7127732766653132, 2.5488843668815639, 2.5248328653472694]`
  - `[6.55751077252848, 6.5011848440695914, 5.027115991544723, 4.4226857324224929]`
  - `[6.5588409967053547, 6.5902779083995151, 6.3846289902553384, 6.039627688146652]`
- `reflect_high_k_fraction_mean` matrix:
  - `[0.24294922660924467, 0.13906324775944112, 0.11172646818271345, 0.10766906380516775]`
  - `[0.82530357342614824, 0.81912102056925717, 0.59471503928237834, 0.48516109208527497]`
  - `[0.82625304890765827, 0.83032699700837842, 0.80360740036150746, 0.75151018098143552]`
- `rms_delta_k_increase_mean` matrix:
  - `[0.82977844353777486, 0.23443712300445163, 0.0705482132207026, 0.046496711686407788]`
  - `[4.0791746188676186, 4.0228486904087308, 2.5487798378838615, 1.9443495787616316]`
  - `[4.0805048430444932, 4.1119417547386528, 3.9062928365944769, 3.5612915344857905]`
- `tap_rms_delay_symbols_mean` matrix:
  - `[37.471572903757632, 35.82348841242306, 35.621289665847, 41.278879008243244]`
  - `[201.14072827568845, 109.35908886719794, 47.846354301354722, 72.62237807072276]`
  - `[134.45187275648524, 103.99139368495908, 95.39888154646296, 78.994400218459091]`
- Generated result artifacts:
  - `sweep_monte_carlo_surface_channel_vertical_result.mat`
  - `sweep_monte_carlo_surface_channel_vertical_abs_H_ref_mean_heatmap.png`
  - `sweep_monte_carlo_surface_channel_vertical_abs_H_ref_std_heatmap.png`
  - `sweep_monte_carlo_surface_channel_vertical_reflect_rms_delta_k_mean_heatmap.png`
  - `sweep_monte_carlo_surface_channel_vertical_reflect_high_k_fraction_mean_heatmap.png`
  - `sweep_monte_carlo_surface_channel_vertical_rms_delta_k_increase_mean_heatmap.png`
  - `sweep_monte_carlo_surface_channel_vertical_tap_rms_delay_mean_heatmap.png`
  - `sweep_monte_carlo_surface_channel_vertical_reflect_high_k_fraction_vs_wind.png`
  - `sweep_monte_carlo_surface_channel_vertical_reflect_rms_delta_k_vs_wind.png`

Remaining issues:
- `mc_count=16` is still a reduced sweep default for cost control, not a converged ensemble estimate.
- The sweep varies only `sea_hs_target`, `sea_wind_speed`, and `sea_seed`; geometry, frequency band, grid, and solver settings stay fixed.
- Wind-speed trends are conditional on PM spectral shape after `sea_hs_target` rescaling and should not be read as an independent wave-height trend.
- No BER/SER or communication-chain Monte Carlo is run in this script.
- The saved statistics remain empirical summaries of the current Kirchhoff phase-screen implementation, not a closed-form statistical channel generator.

## 2026-06-10 C4 Empirical Random Channel Generator Prototype

Changed files for this prototype:
- `build_surface_empirical_channel_model_vertical.m`: reads C3/C3.5 Monte Carlo result files and builds a compact empirical sample model.
- `sample_surface_empirical_channel_vertical.m`: draws nearest-neighbor empirical bootstrap samples from the compact model.
- `validate_surface_empirical_channel_generator_vertical.m`: validates C4 without running PE/WAPE propagation.
- `PROJECT_CONTEXT.md`: records C4 purpose, interface, validation settings, results, and limits.

Purpose and scope:
- C4 converts existing Monte Carlo summary rows into a lightweight empirical random generator for reference-frequency and low-dimensional channel summaries.
- It does not modify `vertical_channel_model`, `vertical_wape_propagator`, `pm_surface_boundary_model`, or communication-chain public outputs.
- It does not run PE/WAPE and does not save full channel structs, spatial fields, or two-dimensional spectra.
- It is an empirical resampler of previously computed C3/C3.5 results, not a new rough-surface scattering theory.

Input and model construction:
- Default input priority:
  - `sweep_monte_carlo_surface_channel_vertical_result.mat`
  - fallback: `monte_carlo_surface_channel_vertical_result.mat`
- The builder reads `run_summary_table` for C3.5 or `summary_table` for C3 and groups rows by sea-state condition.
- Each condition stores `sea_hs_target`, `sea_wind_speed`, `f_ref_hz`, `mc_count`, seed list, compact source samples, and per-metric empirical statistics.
- C4 v1 uses nearest-neighbor matching in normalized `(Hs, wind, f_ref)` space. It does not interpolate between sea states.

Generated sample variables:
- Reference-frequency response:
  - `abs_h_ref`
  - `phase_h_ref_rad`
  - reconstructed `h_ref_complex = abs_h_ref * exp(1i*phase_h_ref_rad)`
- Reflection/direct ratio:
  - `reflect_direct_abs`
  - `reflect_direct_phase_rad`
  - reconstructed `reflect_direct_ratio_complex`
- C2/C2.5 diagnostics:
  - `coupling_nonzero_power_fraction`
  - `coupling_rms_delta_k_rad_per_m`
  - `redistribution_reflect_rms_delta_k_rad_per_m`
  - `redistribution_reflect_high_k_fraction`
  - `redistribution_rms_delta_k_increase_rad_per_m`
- Tap summaries:
  - `tap_count`
  - `tap_rms_delay_symbols`
  - `tap_peak_fraction`

Sampling rule:
- For query sea state `q=(Hs_q, wind_q, f_ref_q)`, select the stored condition
  `c* = argmin_c ||[(Hs_c-Hs_q)/range(Hs), (wind_c-wind_q)/range(wind), (f_ref_c-f_ref_q)/range(f_ref)]||_2`.
- Draw `n_samples` indices with replacement from that condition's Monte Carlo seed rows.
- Copy the corresponding low-dimensional samples into the generated table.
- Metadata records `source_file`, `source_type`, `selection_mode='nearest'`, `sampling_method='empirical_bootstrap_with_replacement'`, matched condition, distance, RNG seed, and sampled source-row indices.

Smoke validation settings:
- Command environment:
  - `SURFACE_EMPIRICAL_SAMPLE_COUNT=200`
  - `SURFACE_EMPIRICAL_MAX_CONDITIONS=2`
- Input file:
  - `sweep_monte_carlo_surface_channel_vertical_result.mat`
- Selected conditions:
  - condition 1: `Hs=0.05`, `wind=3`, `f_ref_hz=5866.666666666667`, `mc_count=4`
  - condition 12: `Hs=1.0`, `wind=12`, `f_ref_hz=5866.666666666667`, `mc_count=4`
- Output file:
  - `surface_empirical_channel_generator_validation_result.mat`
  - size `37379` bytes
- Static check:
  - C4 files do not call `vertical_channel_model`, `vertical_wape_propagator`, or `pm_surface_boundary_model`.

Smoke validation results:
- Condition 1, `Hs=0.05`, `wind=3`:
  - `abs_h_ref` source/generated mean `=0.039000 / 0.038729`, std `=0.032162 / 0.027823`
  - `phase_h_ref_rad` source/generated mean `=-0.85219 / -0.86150`, std `=1.8363 / 1.6075`
  - `reflect_direct_abs` source/generated mean `=1.0942 / 1.0875`, std `=0.35971 / 0.31444`
  - `redistribution_reflect_rms_delta_k_rad_per_m` source/generated mean `=3.3081 / 3.3086`, std `=0.035956 / 0.031058`
  - `redistribution_reflect_high_k_fraction` source/generated mean `=0.24295 / 0.24323`, std `=0.0076396 / 0.0065539`
  - `tap_rms_delay_symbols` source/generated mean `=37.472 / 37.576`, std `=13.653 / 11.866`
- Condition 12, `Hs=1.0`, `wind=12`:
  - `abs_h_ref` source/generated mean `=0.068788 / 0.067981`, std `=0.017683 / 0.014949`
  - `phase_h_ref_rad` source/generated mean `=-1.5262 / -1.4930`, std `=0.66951 / 0.58903`
  - `reflect_direct_abs` source/generated mean `=1.0128 / 1.0014`, std `=0.29010 / 0.24675`
  - `redistribution_reflect_rms_delta_k_rad_per_m` source/generated mean `=6.0396 / 6.0251`, std `=0.40954 / 0.34363`
  - `redistribution_reflect_high_k_fraction` source/generated mean `=0.75151 / 0.74957`, std `=0.055883 / 0.046882`
  - `tap_rms_delay_symbols` source/generated mean `=78.994 / 80.297`, std `=49.800 / 44.629`

Generated C4 validation artifacts:
- `surface_empirical_channel_generator_validation_result.mat`
- `surface_empirical_channel_generator_abs_h_ref_hist.png`
- `surface_empirical_channel_generator_phase_h_ref_rad_hist.png`
- `surface_empirical_channel_generator_reflect_direct_abs_hist.png`
- `surface_empirical_channel_generator_redistribution_reflect_rms_delta_k_rad_per_m_hist.png`
- `surface_empirical_channel_generator_redistribution_reflect_high_k_fraction_hist.png`
- `surface_empirical_channel_generator_tap_rms_delay_symbols_hist.png`

Remaining C4 issues:
- The saved C3.5 source used by this C4 smoke test has only `mc_count=4` per condition, so bootstrap samples replicate a very small empirical support. Rerunning the updated C3.5 script with its new `mc_count=16` default increases the empirical support to 16 seeds per condition.
- Nearest-neighbor matching should not be used to extrapolate outside the scanned `Hs x wind` range.
- C4 v1 does not generate a full frequency response `H_f(f)` and does not replace PE/WAPE propagation.
- Phase samples are bootstrapped from principal phases; no circular distribution is fitted.
- C4 is an empirical reduced generator prototype, not a T matrix, SSA/NLSSA, scattering cross-section model, or closed-form stochastic channel model.

## 2026-06-10 C3.5 Monte Carlo Count Increase

Changed files for this update:
- `sweep_monte_carlo_surface_channel_vertical.m`: increases the default multi-sea-state Monte Carlo count and adds smoke-output redirection controls.
- `PROJECT_CONTEXT.md`: records the new default and smoke validation.

Updated defaults and controls:
- C3.5 default `mc_count` is now `16`.
- Default full reduced sweep size is now `3 Hs x 4 wind x 16 seed = 192` propagations.
- Default seed list is `12345:12360`.
- `SURFACE_MC_SWEEP_COUNT` still overrides the seed count.
- `SURFACE_MC_SWEEP_MAX_CONDITIONS` still truncates conditions for smoke/cost control.
- New optional smoke controls:
  - `SURFACE_MC_SWEEP_RESULT_FILE`
  - `SURFACE_MC_SWEEP_FIGURE_PREFIX`

Smoke validation:
- Command environment:
  - `SURFACE_MC_SWEEP_COUNT=2`
  - `SURFACE_MC_SWEEP_MAX_CONDITIONS=1`
  - `SURFACE_MC_SWEEP_RESULT_FILE=sweep_monte_carlo_surface_channel_vertical_smoke.mat`
  - `SURFACE_MC_SWEEP_FIGURE_PREFIX=sweep_monte_carlo_surface_channel_vertical_smoke_`
- Smoke output:
  - `n_conditions=1`
  - `run_summary_table` rows `=2`
  - `seed_count=2`
  - `condition_summary_table.mc_count=2`
  - `max_invariant_error=6.938893903907228e-18`
  - `direct_drift_max_abs=0`
- The smoke result used a redirected file and PNG prefix, so it did not overwrite the existing C3.5 source file used by C4.

Generated smoke artifacts:
- `sweep_monte_carlo_surface_channel_vertical_smoke.mat`
- `sweep_monte_carlo_surface_channel_vertical_smoke_abs_H_ref_mean_heatmap.png`
- `sweep_monte_carlo_surface_channel_vertical_smoke_abs_H_ref_std_heatmap.png`
- `sweep_monte_carlo_surface_channel_vertical_smoke_reflect_rms_delta_k_mean_heatmap.png`
- `sweep_monte_carlo_surface_channel_vertical_smoke_reflect_high_k_fraction_mean_heatmap.png`
- `sweep_monte_carlo_surface_channel_vertical_smoke_rms_delta_k_increase_mean_heatmap.png`
- `sweep_monte_carlo_surface_channel_vertical_smoke_tap_rms_delay_mean_heatmap.png`
- `sweep_monte_carlo_surface_channel_vertical_smoke_reflect_high_k_fraction_vs_wind.png`
- `sweep_monte_carlo_surface_channel_vertical_smoke_reflect_rms_delta_k_vs_wind.png`

Remaining issue:
- The full 192-propagation C3.5 run was not executed in this update; only the redirected smoke test was run. Existing `sweep_monte_carlo_surface_channel_vertical_result.mat` still reflects the previously saved `mc_count=4` run until the full script is rerun.

## 2026-06-11 Realization-Based Platform and C4 Wideband Bootstrap Update

Changed files for this update:
- `plot_c35_core_heatmaps_vertical.m`: reads the C3.5 result file and regenerates the five core heatmaps.
- `monte_carlo_comm_surface_psk_representative_vertical.m`: runs communication Monte Carlo on weak/mid/strong representative sea states while keeping the D2 link fixed.
- `build_surface_empirical_channel_model_vertical.m`: detects stored wideband `H_f` Monte Carlo samples in C3/C3.5 result files.
- `sample_surface_empirical_channel_vertical.m`: adds `sample_mode='wideband_hf'` and `sample_mode='tap_level'`.
- `validate_surface_empirical_channel_generator_vertical.m`: validates summary, wideband, and tap-level C4 modes.
- `reports/realization_to_empirical_channel_generator_report.md`: records the stage summary and main results.

Fixed communication policy:
- D2 is treated as fixed for this phase.
- `receive_window_mode='peak_sync'`.
- `ebn0_reference='rx_clean'`.
- No changes were made to `comm_main_vertical_psk.m`.

C3.5 full rerun:
- `sea_hs_target=[0.05,0.5,1.0]`.
- `sea_wind_speed=[3,5,8,12]`.
- `mc_count=16`, `seed_list=12345:12360`.
- Total propagation count: `12 x 16 = 192`.
- `run_summary_table` rows: `192`.
- `condition_summary_table` rows: `12`.
- Maximum `H_f-(H_direct_f+H_reflect_f)` invariant error: `1.551583845779546e-17`.
- Maximum direct-path drift: `1.390486644391991e-17`.
- All 12 conditions contain stored wideband `H_f.samples` with `Nf=16`.

C3.5 core heatmaps generated:
- `c35_core_abs_H_ref_mean_heatmap.png`
- `c35_core_abs_H_ref_std_heatmap.png`
- `c35_core_reflect_rms_delta_k_mean_heatmap.png`
- `c35_core_reflect_high_k_fraction_mean_heatmap.png`
- `c35_core_tap_rms_delay_mean_heatmap.png`

C3.5 selected condition results:
- weak sea, `Hs=0.05`, `wind=5`:
  - `abs_H_ref_mean=0.029723`
  - `abs_H_ref_std=0.019989`
  - `reflect_rms_delta_k_mean=2.6774`
  - `reflect_high_k_fraction_mean=0.13695`
  - `tap_rms_delay_symbols_mean=31.644`
- mid representative sea, `Hs=0.5`, `wind=8`:
  - `abs_H_ref_mean=0.045544`
  - `abs_H_ref_std=0.018677`
  - `reflect_rms_delta_k_mean=4.9028`
  - `reflect_high_k_fraction_mean=0.56899`
  - `tap_rms_delay_symbols_mean=97.101`
- strong representative sea, `Hs=1.0`, `wind=12`:
  - `abs_H_ref_mean=0.054750`
  - `abs_H_ref_std=0.021414`
  - `reflect_rms_delta_k_mean=6.0107`
  - `reflect_high_k_fraction_mean=0.74638`
  - `tap_rms_delay_symbols_mean=109.35`

Representative communication Monte Carlo:
- Sea states:
  - weak: `Hs=0.05`, `wind=5`
  - mid: `Hs=0.5`, `wind=8`
  - strong: `Hs=1.0`, `wind=12`
- Scenarios: `direct_only`, `direct_plus_reflect`.
- `mc_count=16`.
- Total channel runs: `3 x 2 x 16 = 96`.
- Maximum channel invariant error: `1.387778780781446e-17`.
- Maximum direct-only reflected response: `0`.
- Selected BER results:
  - weak direct-plus-reflect: `BER_mean` at `[0,10,20] dB = [0.11606, 3.125e-05, 0]`
  - mid direct-plus-reflect: `BER_mean` at `[0,10,20] dB = [0.17438, 0.045906, 0.037469]`
  - strong direct-plus-reflect: `BER_mean` at `[0,10,20] dB = [0.15978, 0.018125, 0.00875]`
- Mid and strong reflected cases retain residual high-Eb/N0 BER for some seeds under the current reduced-grid channel and MMSE receiver.

C4 upgraded generator:
- `summary` mode retains low-dimensional bootstrap behavior.
- `wideband_hf` mode bootstraps stored `H_f(f)` samples from the matched C3.5 condition.
- `tap_level` mode derives baseband/tap samples from bootstrapped `H_f(f)` without rerunning PE/WAPE.
- C4 validation using the new C3.5 source:
  - comparison rows: `18`
  - output file size: `219191` bytes
  - `H_f_sample_size=[16 50]`
  - `idx_f_ref=8`
  - `tap_count_mean=501.92`
  - `tap_rms_delay_mean=40.0435`
  - `tap_peak_fraction_mean=0.9857`
  - tap metrics finite: true
  - no PE/WAPE propagation run: true

Stage conclusion:
- The current platform is a realization-based PE/WAPE simulation and empirical Monte Carlo analysis pipeline.
- The next stage should focus on empirical statistics-based channel generation from the C3.5 sample ensemble.
- C4 remains an empirical bootstrap generator. It is not a T-matrix, SSA/NLSSA, physical scattering cross-section model, or closed-form statistical channel model.

## 2026-06-10 Communication Monte Carlo PSK Validation

Changed files for this script:
- `monte_carlo_comm_surface_psk_vertical.m`: independent reduced-grid QPSK BER/SER Monte Carlo script.
- `PROJECT_CONTEXT.md`: records communication Monte Carlo scope, formulas, random seed policy, outputs, validation settings, and limits.

Reduced communication Monte Carlo settings:
- Fixed channel setup: `c0=1500`, `z_tx=100`, `z_rx=3`, `xw=yw=50`, `nx=ny=128`, `sigma_src_m=0.4`, `stepz_lamb=0.5`, `sponge_ratio=0.12`, `alpha_max_np_per_m=0.15`, `env_mode='uniform'`, `show_figures=false`, `save_mode='rx_only'`, `use_gpu=false`, `surface_boundary_model='kirchhoff_kdomain'`, `enforce_1_over_R=false`.
- Wideband setup: `enable_wideband=true`, `f_band_hz=[4000 8000]`, `Nf_min=Nf_max=16`, `f_ref_hz=6000`.
- Sea conditions: weak `sea_hs_target=0.05`, strong `sea_hs_target=1.0`; both use `sea_wind_speed=5`.
- Scenarios: `direct_only` and `direct_plus_reflect`.
- Communication setup: `M=4`, `n_sym=1000`, `EbN0_dB_list=0:2:20`, `symbol_rate_hz=1000`.
- Seed policy:
  - `bits_seed=9000`, fixed for all runs.
  - `seed_list=12345:12348` for default sea-surface samples.
  - AWGN seed varies with condition, scenario, seed index, and Eb/N0 index.
- Smoke test: `SURFACE_COMM_MC_COUNT=2`.
- Default reduced run: `2` sea conditions, `2` scenarios, `mc_count=4`, `16` total channel runs.

Validation results:
- Smoke test completed successfully and produced finite `.mat`, table, and PNG outputs.
- Default 16-run communication Monte Carlo completed successfully.
- BER/SER sample array dimensions:
  - `ber_samples`: `[11, 4, 2, 2]`
  - `ser_samples`: `[11, 4, 2, 2]`
- Curve summary table rows: `44`, matching `2 sea conditions x 2 scenarios x 11 Eb/N0 values`.
- Channel invariant:
  - `max(abs(H_f-(H_direct_f+H_reflect_f)))=1.3877787807814457e-17`
- Direct-only reflected response:
  - `max(abs(H_reflect_f))=0`
- Finite-field check:
  - all selected BER/SER mean, standard deviation, and quantile fields are finite.
- BER mean curves for `EbN0_dB_list=[0,2,4,6,8,10,12,14,16,18,20]`:
  - weak direct-only: `[0.50025, 0.493375, 0.503, 0.496875, 0.50825, 0.498125, 0.49925, 0.503, 0.49525, 0.49875, 0.4995]`
  - weak direct-plus-reflect: `[0.497, 0.505625, 0.507, 0.500375, 0.51075, 0.50325, 0.498375, 0.512375, 0.5075, 0.498625, 0.507375]`
  - strong direct-only: `[0.504, 0.511125, 0.493875, 0.502125, 0.493875, 0.500875, 0.516625, 0.4955, 0.499875, 0.499375, 0.507375]`
  - strong direct-plus-reflect: `[0.513375, 0.502125, 0.497875, 0.50625, 0.50325, 0.505125, 0.504875, 0.505875, 0.507625, 0.5065, 0.494625]`
- SER mean curves:
  - weak direct-only: `[0.75075, 0.745, 0.7585, 0.75, 0.75975, 0.74775, 0.75925, 0.752, 0.7495, 0.7575, 0.75375]`
  - weak direct-plus-reflect: `[0.747, 0.755, 0.75125, 0.74825, 0.76275, 0.75825, 0.741, 0.76225, 0.754, 0.748, 0.7535]`
  - strong direct-only: `[0.75525, 0.76025, 0.74425, 0.7515, 0.7435, 0.7535, 0.766, 0.74725, 0.75925, 0.75225, 0.7615]`
  - strong direct-plus-reflect: `[0.76, 0.75075, 0.748, 0.7515, 0.75325, 0.751, 0.75575, 0.753, 0.7605, 0.75725, 0.73925]`
- Monotonicity check:
  - BER nonmonotonic increase counts matrix `[scenario x condition] = [[6, 6], [5, 4]]`
  - SER nonmonotonic increase counts matrix `[scenario x condition] = [[4, 6], [4, 4]]`
  - These reduced runs do not show the expected BER/SER decrease with Eb/N0; curves remain close to random hard-decision levels.
- Generated result artifacts:
  - `monte_carlo_comm_surface_psk_vertical_result.mat`
  - `monte_carlo_comm_surface_psk_vertical_BER_weak_strong_direct_plus_reflect.png`
  - `monte_carlo_comm_surface_psk_vertical_SER_weak_strong_direct_plus_reflect.png`
  - `monte_carlo_comm_surface_psk_vertical_BER_weak_direct_vs_reflect.png`
  - `monte_carlo_comm_surface_psk_vertical_BER_strong_direct_vs_reflect.png`
  - `monte_carlo_comm_surface_psk_vertical_SER_weak_direct_vs_reflect.png`
  - `monte_carlo_comm_surface_psk_vertical_SER_strong_direct_vs_reflect.png`

Remaining issues:
- The communication Monte Carlo infrastructure is functional, but current reduced BER/SER curves do not validate a useful demodulation-performance trend.
- The issue appears inherited from the existing baseband tap/MMSE communication chain behavior rather than the Monte Carlo aggregation itself.
- `n_sym=1000` and `mc_count=4` are reduced defaults; they are not sufficient for low-BER tail estimation.
- The script does not modify `comm_main_vertical_psk.m` and does not introduce a new BER/SER theory.

## 2026-06-10 Minimal Communication Closed-Loop Validation

Changed files for this validation:
- `validate_comm_link_minimal_vertical.m`: independent D1 communication-chain validation script.
- `PROJECT_CONTEXT.md`: records D1 validation purpose, cases, results, and diagnosis.

Validation settings:
- Communication: `M=4`, `n_sym=2000`, `EbN0_dB_list=[0,2,4,6,8,10,12,14,16,18,20]`, `bits_seed=9000`, `noise_seed_base=7000`.
- PE reduced direct-only case: `nx=ny=128`, `xw=yw=50`, `z_tx=100`, `z_rx=3`, `sigma_src_m=0.4`, `enable_wideband=true`, `f_band_hz=[4000 8000]`, `Nf_min=Nf_max=16`, `f_ref_hz=6000`, `enable_surface_reflection=false`.

Validation results:
- Unit channel:
  - noiseless `BER=0`, `SER=0`
  - `BER(20 dB)=0`, `SER(20 dB)=0`
  - BER curve `[0.11675, 0.04525, 0.0155, 0.0015, 0.00075, 0, 0, 0, 0, 0, 0]`
  - nonmonotonic counts: BER `0`, SER `0`
- Single complex tap:
  - noiseless `BER=0`, `SER=0`
  - `BER(20 dB)=0`, `SER(20 dB)=0`
  - BER curve `[0.22275, 0.15075, 0.07425, 0.0365, 0.01, 0.002, 0, 0, 0, 0, 0]`
  - nonmonotonic counts: BER `0`, SER `0`
- Known short multipath with current `conv(...,'same')` slicing:
  - `tap_count=3`, `peak_index=1`, `peak_fraction=0.8733624454148472`
  - noiseless `BER=0.506`, `SER=0.758`
  - `BER(20 dB)=0.50575`, `SER(20 dB)=0.7575`
  - This fails even without noise, so it is not an AWGN or PSK mapper issue.
- Same known short multipath with diagnostic `causal_head` slicing:
  - noiseless `BER=0`, `SER=0`
  - `BER(20 dB)=0`, `SER(20 dB)=0`
  - BER curve `[0.14575, 0.07975, 0.03475, 0.00675, 0.0005, 0, 0, 0, 0, 0, 0]`
  - nonmonotonic counts: BER `0`, SER `0`
- PE direct-only current taps with current `conv(...,'same')` slicing:
  - `H_f` invariant error `=0`
  - `max(abs(H_reflect_f))=0`
  - `peak_index_before_alignment=1`
  - `h_full` peak fraction `=0.9976004631466816`
  - current retained `tap_count=1997`
  - noiseless `BER=0.51775`, `SER=0.768`
  - `BER(20 dB)=0.49175`, `SER(20 dB)=0.739`
- PE direct-only peak-aligned taps with current `conv(...,'same')` slicing:
  - peak alignment does not change the result because the dominant tap was already at index 1.
  - noiseless `BER=0.51775`, `SER=0.768`
  - `BER(20 dB)=0.4925`, `SER(20 dB)=0.7415`
- PE direct-only peak-aligned taps with diagnostic `causal_head` slicing:
  - noiseless `BER=0`, `SER=0`
  - BER curve `[0.49725, 0.47675, 0.475, 0.4745, 0.47275, 0.4625, 0.45075, 0.431, 0.42525, 0.3695, 0.342]`
  - SER curve `[0.741, 0.7185, 0.7055, 0.704, 0.701, 0.677, 0.6655, 0.617, 0.6005, 0.5195, 0.479]`
  - nonmonotonic counts: BER `0`, SER `0`
  - effective SNR at nominal `20 dB` Eb/N0 is still `-4.9334090404706892 dB`, because current AWGN injection references transmit-symbol power while the PE channel gain is much smaller.

Diagnosis:
- The PSK mapper/demapper, AWGN injection, and single-tap equalization are functional.
- The current communication scripts fail for multi-tap channels mainly because `conv(tx_symbols,h_bb,'same')` center-crops a causal impulse response. For a long PE `h_bb`, this introduces a large symbol timing offset even when the dominant tap is already at index 1.
- Main-peak alignment alone is not sufficient for PE taps; receive-window selection must also respect causal timing.
- After causal receive-window selection, PE direct-only noiseless BER/SER becomes zero, but noisy BER remains high because the current Eb/N0 noise calibration is referenced to transmitted symbols rather than received/channel-output power.

Generated result artifacts:
- `validate_comm_link_minimal_vertical_result.mat`
- `validate_comm_link_minimal_vertical_simple_channel_BER.png`
- `validate_comm_link_minimal_vertical_simple_channel_SER.png`
- `validate_comm_link_minimal_vertical_pe_h_bb_current_taps.png`
- `validate_comm_link_minimal_vertical_pe_h_bb_aligned_taps.png`

Remaining issues:
- D1 does not change production communication scripts.
- A follow-up fix should replace or parameterize `conv(...,'same')` receive slicing for causal channel taps.
- A separate noise-calibration decision is needed: Eb/N0 can be referenced to transmit-symbol power, received clean-signal power, or post-equalizer noise enhancement, and these produce different BER interpretations.

## 2026-06-10 D2 Communication Receive Window and Eb/N0 Reference Update

Changed files for this update:
- `comm_main_vertical_psk.m`: switches the reference communication demo to explicit receive-window selection and received-signal Eb/N0 reference.
- `monte_carlo_comm_surface_psk_vertical.m`: applies the same communication policy to BER/SER Monte Carlo statistics.
- `validate_comm_link_minimal_vertical.m`: updates D1 diagnostics to include the D2 `peak_sync` policy.
- `PROJECT_CONTEXT.md`: records the D2 communication policy and validation status.

D2 policy:
- Receive window:
  - default `receive_window_mode='peak_sync'`.
  - Find the dominant tap index in `h_bb`.
  - Use `h_eq = h_bb(peak_index:end)` as the equalizer tap vector.
  - Compute `rx_full = conv(tx_symbols,h_eq,'full')`.
  - Use `rx_clean = rx_full(1:n_sym)` as the symbol-aligned receive vector.
  - Record `receive_meta.peak_index_original`, `receive_meta.original_tap_count`, `receive_meta.equalizer_tap_count`, and `receive_meta.discarded_pre_peak_energy_fraction`.
- Eb/N0 reference:
  - default `ebn0_reference='rx_clean'`.
  - `noise_inject_vertical(rx_clean, noise_cfg, rx_clean)` is used for the default communication BER/SER path.
  - local diagnostic support for `ebn0_reference='tx_symbols'` is retained in the scripts, but it is no longer the default.
- Compatibility:
  - `comm_main_vertical_psk.m` still saves the previous result fields and adds `receive_window_mode`, `receive_meta`, `ebn0_reference`, and `h_eq`.
  - `monte_carlo_comm_surface_psk_vertical.m` keeps the same result tables and sample arrays, and adds receive-window metadata in per-run summaries.
  - `comm_main_vertical_psk.m` remains the reference single-run demo; no changes were made to `vertical_channel_model` or propagation outputs.

Expected effect from D1 diagnosis:
- Unit and single-tap channels should remain zero-error in noiseless mode and should show BER/SER reduction with Eb/N0.
- Short causal multipath and PE direct-only channels should no longer fail noiseless validation due to `conv(...,'same')` center cropping.
- Receiver-referenced Eb/N0 should make nominal Eb/N0 correspond to the clean received waveform, so PE path loss no longer drives the effective SNR deeply negative.

Validation status:
- Static diff checks passed for the updated communication scripts.
- `validate_comm_link_minimal_vertical` completed successfully after D2.
- `monte_carlo_comm_surface_psk_vertical` completed successfully with `SURFACE_COMM_MC_COUNT=2`.
- `comm_main_vertical_psk` completed successfully with the reference `256 x 256`, `Nf=32` setup.

D2 minimal validation results:
- Unit channel:
  - noiseless `BER=0`, `SER=0`;
  - `BER(20 dB)=0`, `SER(20 dB)=0`;
  - BER decreases to zero with no nonmonotonic steps.
- Single complex tap:
  - noiseless `BER=0`, `SER=0`;
  - `BER(20 dB)=0`, `SER(20 dB)=0`;
  - BER decreases to zero with no nonmonotonic steps.
- Known short multipath:
  - legacy `same` still fails: noiseless `BER=0.506`, `SER=0.758`;
  - D2 `peak_sync` passes: noiseless `BER=0`, `SER=0`, `BER(20 dB)=0`, `SER(20 dB)=0`.
- PE direct-only reduced D1 case:
  - legacy `same` still fails: noiseless `BER=0.51775`, `SER=0.768`;
  - peak-aligned plus legacy `same` still fails, confirming alignment alone is insufficient;
  - D2 `peak_sync` passes: noiseless `BER=0`, `SER=0`, `BER(20 dB)=0`, `SER(20 dB)=0`.

D2 communication Monte Carlo smoke results:
- Setup: weak/strong sea conditions, direct-only/direct-plus-reflect, `SURFACE_COMM_MC_COUNT=2`, `EbN0_dB_list=[0,2,4,6,8,10,12,14,16,18,20]`.
- Weak direct-only BER mean:
  - `[0.111, 0.06375, 0.017, 0.004, 0, 0, 0, 0, 0, 0, 0]`
- Weak direct-plus-reflect BER mean:
  - `[0.1195, 0.06275, 0.01925, 0.00525, 0.00025, 0, 0, 0, 0, 0, 0]`
- Strong direct-only BER mean:
  - `[0.11825, 0.04925, 0.01525, 0.00275, 0.00025, 0, 0, 0, 0, 0, 0]`
- Strong direct-plus-reflect BER mean:
  - `[0.23825, 0.18975, 0.1535, 0.092, 0.064, 0.04425, 0.03075, 0.02725, 0.019, 0.016, 0.014]`
- The strong reflected case remains harder than direct-only, but now has a monotone decreasing BER trend instead of random-decision behavior.

D2 reference communication demo results:
- `comm_main_vertical_psk` now reports `rx_window=peak_sync`, `ebn0_ref=rx_clean`.
- Direct-only reference run:
  - `BER=[0.118, 0.05325, 0.02425, 0.0035, 0.00025, 0, 0, 0, 0, 0, 0]`
  - `SER=[0.1585, 0.0735, 0.0295, 0.005, 0.0005, 0, 0, 0, 0, 0, 0]`
- Direct-plus-reflect reference run:
  - `BER=[0.151, 0.102, 0.05975, 0.036, 0.01975, 0.0095, 0.004, 0.0015, 0.00075, 0.00025, 0]`
  - `SER=[0.21, 0.139, 0.0765, 0.0485, 0.026, 0.0135, 0.007, 0.0015, 0.001, 0.0005, 0]`
- Effective SNR now tracks nominal Eb/N0 plus QPSK bits-per-symbol scaling and no longer collapses to large negative values under PE path loss.

Remaining D2 issues:
- `peak_sync` is a deterministic dominant-tap synchronization rule; it is not a timing recovery loop.
- `rx_clean` Eb/N0 is a receiver-side performance convention. It is useful for comparing channel distortion at controlled received SNR, but it no longer includes absolute path-loss penalty in the noise power.
- The strong reflected channel can still have residual ISI/noise enhancement, so BER may remain nonzero at 20 dB for small `n_sym`/`mc_count` smoke tests.

## 2026-06-15 SSA-Like Statistical Surface Kernel

Changed files for this update:
- `vertical_channel_model.m`: adds validated `paramsV.surface_ssa_random_scatter`, `paramsV.surface_ssa_scatter_scale`, `paramsV.surface_ssa_seed_offset`, `paramsV.surface_ssa_kernel_mode`, `paramsV.surface_ssa_geometry_source_id`, `paramsV.surface_ssa_kz_branch`, `paramsV.surface_ssa_conv_padding`, and allows `surface_boundary_model='ssa_stat_kernel'`.
- `vertical_wape_propagator.m`: passes the SSA-like random-scatter controls and frequency index into the surface boundary module; keeps the two-segment reflected PE path and `H_f = H_direct_f + H_reflect_f`.
- `pm_surface_boundary_model.m`: adds the SSA-like PM-spectrum statistical boundary branch.
- `vertical_comm_guide.md` and `PROJECT_CONTEXT.md`: document formulas, interfaces, limits, and validation status.

Implemented interface and formulas:
- `ssa_stat_kernel` does not synthesize `xi(x,y)`. It scales the PM height spectrum directly so `sum(W_eta(:))*dkx*dky = sigma_eta^2`, where `sigma_eta = sea_hs_target/4`.
- The coherent term is `R_coh = R0*exp(-0.5*(2*k0*phase_factor_eff)^2*sigma_eta^2)`.
- The raw incoherent power uses the engineering kernel `P_sca_raw = surface_ssa_scatter_scale*abs(R0)^2*(2*k0*phase_factor_eff)^2*dkx*dky*circconv(W_eta,abs(Psi_inc_k).^2)`.
- Energy limiting enforces `E_coh+E_sca <= E_inc`; after random scatter synthesis, the final combined reflected spectrum is also checked against `E_inc`.
- The random scatter seed is `seed_ssa = sea_seed + surface_ssa_seed_offset + frequency_index - 1`, with default `surface_ssa_seed_offset=100000`.
- The current engineering kernel is explicitly identified by `surface_ssa_kernel_mode='pm_convolution'`.
- `ssa1_geometry` and `ssa1_debug_dense` implement the `SSA.md` first-order Dirichlet geometry factor `G_SSA1(K,K';f)=4*gamma(K,f)*gamma(K',f)`.

Metadata:
- `roughness_meta.ssa_stat_kernel_meta` records `sigma_eta_m`, `Hs_target_m`, `R_coh`, `P_sca` stats, `E_inc`, `E_coh`, `E_sca_raw`, `E_sca`, `E_ref`, energy scale factors, `seed_ssa`, formulas, normalization notes, and limitations.
- The metadata now also records `kernel_mode`, `geometry_source_id`, `kz_branch`, `conv_padding`, `E_sca_limited`, `energy_limit_applied`, `energy_conservation_error`, and `propagating_bin_fraction`.
- For this branch, `surface_elevation=[]`, `delta_phi=[]`, and `roughness_meta.surface_realization_generated=false`.

Validation:
- Reduced scalar defaults: `nx=ny=128`, `xw=yw=50`, `z_tx=100`, `z_rx=3`, `sigma_src_m=0.4`, `show_figures=false`, `save_mode='rx_only'`, `enforce_1_over_R=false`.
- Default `kirchhoff_spatial`: invariant error `max(abs(H_f-H_direct_f-H_reflect_f)) = 0`.
- `ssa_stat_kernel`, `sea_hs_target=0`: `R_coh=-1`, `E_sca=0`, and flat-model difference from `kirchhoff_spatial` was `7.76e-17`.
- `ssa_stat_kernel`, `sea_hs_target=0.5`: `E_coh+E_sca = E_inc = 53831.9373`; same seed gave identical `H_f`, changed seed kept `H_direct_f` fixed and changed `H_reflect_f`.
- `surface_ssa_random_scatter=false`: `P_sca.sum=53831.9373` remained available in metadata, while `E_sca=0` and no random scatter was added to `psi_ref`.
- Direct-only with `ssa_stat_kernel`: `H_reflect_f=0` and the public channel fields remained present.
- Reduced wideband PSK smoke: direct-only and direct-plus-reflect both consumed `H_f` and completed the baseband tap/equalizer flow.
- Weak sea comparison with `sea_hs_target=0.05`: `kirchhoff_spatial` gave `|H|=0.074770323`, phase `-1.2009652`, BER `[0.11133 0 0]` at Eb/N0 `[0 10 20]`; `ssa_stat_kernel` gave `|H|=0.027797463`, phase `-1.7734442`, BER `[0.13086 0 0]`. The BER trend is consistent; pointwise channel equality is not expected.

Limitations and remaining issues:
- This is an SSA-like statistical kernel for validating the PM-spectrum-to-random-channel path. It is not a strict SSA/NLSSA scattering solver and does not include calibrated angular geometry factors or scattering cross sections.
- `surface_ssa_scatter_scale` is an engineering normalization knob, not a physical calibration.
- The old Kirchhoff realization models remain the default and should be used for legacy comparison unless the statistical branch is explicitly requested.

## 2026-06-15 SSA Kernel Mode Interface and SSA1 Geometry

Changed behavior:
- `surface_ssa_kernel_mode='pm_convolution'` is the default and reproduces the current engineering PM-spectrum convolution baseline.
- `surface_ssa_kernel_mode='ssa1_geometry'` implements the `SSA.md` first-order Dirichlet / pressure-release power geometry:
  `G_SSA1(K,K';f)=4*gamma(K,f)*gamma(K',f)`.
- The implemented FFT form is `A(K')=gamma(K',f)*abs(Psi_inc(K'))^2`, `B=circconv(W_eta,A)`, and `P_sca_raw(K)=4*C_norm*gamma(K,f)*B(K)*dkx*dky`, with `C_norm=surface_ssa_scatter_scale`.
- `surface_ssa_kernel_mode='ssa1_debug_dense'` computes the same periodic sum explicitly on small grids for FFT-vs-dense validation.
- `surface_ssa_kz_branch='downward_positive_real'` computes `kz(K)=sqrt(max(k0^2-|K|^2,0))`; non-propagating bins are excluded from the propagating fraction metadata but the baseline `pm_convolution` formula remains numerically unchanged.
- `surface_ssa_conv_padding='periodic'` keeps the existing FFT circular convolution. `zero_padded` is accepted at the public config boundary but rejected by the kernel until aliasing-control implementation is added.

Current formula status:
- The first-order Dirichlet SSA geometry factor from `SSA.md` has been implemented.
- No calibrated angular scattering cross section has been implemented.
- No strict boundary-condition mapping from arbitrary `surface_reflect_coeff` to Dirichlet/Neumann/impedance SSA factors has been implemented; `ssa1_geometry` and `ssa1_debug_dense` require `surface_reflect_coeff=-1`.
- The code must not be described as NLSSA, impedance-boundary SSA, or experimentally calibrated rough-surface scattering.

Validation additions:
- `validate_ssa_stat_kernel_vertical.m` checks that `pm_convolution` records the new metadata fields, `ssa1_geometry` degenerates correctly for `Hs=0`, and `ssa1_debug_dense` matches the FFT sum on a small grid through compact `P_sca` metadata.
- Sweep and communication summary tables carry kernel-mode energy audit fields: `E_sca_limited`, `energy_conservation_error`, `energy_limit_applied`, and `propagating_bin_fraction`.
- `sweep_ssa_stat_kernel_surface_channel_vertical.m` and `monte_carlo_comm_ssa_stat_kernel_psk_vertical.m` compare `kirchhoff_spatial`, `ssa_pm_convolution`, and `ssa1_geometry`.
- `plot_ssa_stat_kernel_report_vertical.m` emits `ssa_kernel_mode_*` figures; dense-vs-FFT is available through the validation result rather than default wideband sweeps.

Latest executed validation:
- `validate_ssa_stat_kernel_vertical` completed successfully.
- `Hs=0` flat degeneration:
  - `pm_convolution` vs `kirchhoff_spatial`: max response difference `7.7579e-17`.
  - `ssa1_geometry` vs `kirchhoff_spatial`: max response difference `7.7579e-17`.
- `ssa1_geometry` energy audit: `energy_conservation_error=1.3516e-16`, below the `1e-12` validation tolerance.
- `ssa1_debug_dense` vs FFT `ssa1_geometry`: compact `P_sca_raw` relative sum error `2.5247e-15`.
- `sweep_ssa_stat_kernel_surface_channel_vertical` completed with `SSA_STAT_MC_COUNT=1`, `Hs=[0 0.05 0.2 0.5]`, and model labels `kirchhoff_spatial`, `ssa_pm_convolution`, `ssa1_geometry`.
- In that sweep, both statistical kernels matched the `Hs=0` flat Kirchhoff response within `1.2533e-16`.
- `monte_carlo_comm_ssa_stat_kernel_psk_vertical` completed with `SSA_STAT_COMM_MC_COUNT=1`; all three model labels propagated through the existing `H_f` communication path and produced BER/SER summaries.
- `plot_ssa_stat_kernel_report_vertical` completed and regenerated `ssa_stat_kernel_report_*` plus `ssa_kernel_mode_*` report figures.

## 2026-06-15 SSA Statistical Validation and Report Scripts

Changed files for this validation/report update:
- `validate_ssa_stat_kernel_vertical.m`: new reduced scalar validation script.
- `sweep_ssa_stat_kernel_surface_channel_vertical.m`: new multi-Hs, multi-seed channel statistics script.
- `monte_carlo_comm_ssa_stat_kernel_psk_vertical.m`: new reduced QPSK BER/SER statistics script using the existing `H_f` communication path.
- `sweep_ssa_scatter_scale_sensitivity_vertical.m`: new reduced scalar sensitivity script for `surface_ssa_scatter_scale`.
- `plot_ssa_stat_kernel_report_vertical.m`: new plot-only report script that reads saved `.mat` results.
- `vertical_comm_guide.md` and `PROJECT_CONTEXT.md`: document run order, outputs, validation scope, and limitations.

Validation scope:
- Sea states: `sea_hs_target=[0 0.05 0.2 0.5]`, `sea_wind_speed=5`.
- Surface models: `kirchhoff_spatial` and `ssa_stat_kernel`.
- Channel sweep default seeds: `12345+(0:7)`, override with `SSA_STAT_MC_COUNT`.
- Communication default seeds: `12345+(0:3)`, override with `SSA_STAT_COMM_MC_COUNT`.
- Scatter-scale sweep default seeds: `12345+(0:3)`, override with `SSA_SCALE_SWEEP_MC_COUNT`.
- Scatter-scale sweep defaults: `sea_hs_target=[0.05 0.2 0.5]`, `surface_ssa_scatter_scale=[0 0.25 1 4]`, kernels `ssa_pm_convolution` and `ssa1_geometry`.
- Smoke-test condition limits: `SSA_STAT_SWEEP_MAX_CONDITIONS` and `SSA_STAT_COMM_MAX_CONDITIONS`.
- Scale-sweep smoke-test limit: `SSA_SCALE_SWEEP_MAX_CONDITIONS`.
- Reduced wideband channel defaults: `enable_wideband=true`, `f_band_hz=[4000 8000]`, `Nf_min=Nf_max=16`, `f_ref_hz=6000`, `nx=ny=128`, `show_figures=false`, `enforce_1_over_R=false`.

Acceptance metrics checked or summarized:
- `max(abs(H_f-H_direct_f-H_reflect_f)) <= 1e-10`.
- `h_total == H_f(idx_f_ref)` and `h_reflect == H_reflect_f(idx_f_ref)` within roundoff.
- `direct_only` keeps `H_reflect_f=0` and public channel fields present.
- `Hs=0` gives `sigma_eta_m=0`, `R_coh=surface_reflect_coeff`, `P_sca.sum=0`, `E_sca=0`, and a flat-response match against `kirchhoff_spatial`.
- `E_coh+E_sca <= E_inc + 1e-12*max(E_inc,1)` and `E_ref <= E_inc + 1e-12*max(E_inc,1)` for `ssa_stat_kernel` metadata.
- `surface_ssa_random_scatter=false` leaves finite `P_sca` metadata but sets realized `E_sca=0`, so downstream `H_reflect_f` receives only the coherent term.
- `surface_ssa_scatter_scale=0` gives zero raw scatter energy in both `pm_convolution` and `ssa1_geometry`.
- `E_sca_raw/E_inc` should be nondecreasing as `surface_ssa_scatter_scale` increases for a fixed kernel, sea state, and seed.
- Statistical summaries report mean/std of `|H(f)|`, `|H_reflect(f)|`, `angle(H(f_ref))`, `|h_total|`, `|h_reflect|`, `E_sca/E_inc`, `E_ref/E_inc`, and energy scaling.
- BER/SER summaries compare qualitative trends with paired bit/noise seeds, not pointwise equality between physical models.

Report outputs:
- Channel result: `sweep_ssa_stat_kernel_surface_channel_vertical_result.mat`.
- Communication result: `monte_carlo_comm_ssa_stat_kernel_psk_vertical_result.mat`.
- Validation result: `validate_ssa_stat_kernel_vertical_result.mat`.
- Scale-sensitivity result: `sweep_ssa_scatter_scale_sensitivity_vertical_result.mat`.
- Figures:
  - `ssa_stat_kernel_report_abs_H_f_mean_std.png`
  - `ssa_stat_kernel_report_abs_H_reflect_f_mean_std.png`
  - `ssa_stat_kernel_report_phase_H_f_ref_vs_Hs.png`
  - `ssa_stat_kernel_report_reflect_energy_vs_Hs.png`
  - `ssa_stat_kernel_report_Esca_Einc_vs_Hs.png`
  - `ssa_stat_kernel_report_energy_scale_vs_Hs.png`
  - `ssa_stat_kernel_report_abs_h_total_model_compare.png`
  - `ssa_stat_kernel_report_abs_h_reflect_model_compare.png`
  - `ssa_stat_kernel_report_BER_model_compare.png`
  - `ssa_stat_kernel_report_SER_model_compare.png`
  - `ssa_stat_kernel_report_Hs0_flat_check.png`
  - `ssa_stat_kernel_report_metadata_energy_table.mat`

Limitations:
- These scripts are validation and reporting harnesses around the current implementation. They do not modify `comm_main_vertical_psk.m`, modem/noise/MMSE policy, WAPE marching, or `local_march_field`.
- Full spatial fields, `psi_ref`, `P_sca`, and `W_eta` matrices are intentionally not saved in the statistical result files.
- The current model remains an SSA-like engineering/statistical kernel. Strict SSA angular geometry factors and calibrated scattering cross sections remain future work.

Latest scale-sensitivity smoke result:
- `sweep_ssa_scatter_scale_sensitivity_vertical` completed with `SSA_SCALE_SWEEP_MC_COUNT=1`, `Hs=[0.05 0.2 0.5]`, and scales `[0 0.25 1 4]`.
- For `Hs=0.05`, `E_sca_raw/E_inc` was `[0 0.39478 1.5791 6.3165]` for `ssa_pm_convolution` and `[0 0.09768 0.39072 1.5629]` for `ssa1_geometry`.
- `energy_conservation_error_max=0` for every listed scale-sweep condition.
- The sweep uses `surface_ssa_random_scatter=false`; therefore `E_sca/E_inc=0` in the realized reflected field, while `E_sca_limited/E_inc` records the limited scatter-energy budget for metadata analysis.
