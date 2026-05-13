# BUBBLE_EXTENSION_SPEC.md

## 1. Purpose

This document defines how to add wind-generated bubble-layer effects to the current MATLAB vertical underwater acoustic channel project.

The goal is to extend the existing vertical WAPE channel model so that bubbles affect propagation as a near-surface, frequency-dependent complex medium. Bubbles must not be implemented as receiver-side noise, and must not be implemented as a fixed scalar path loss. The bubble model must produce effective sound speed and physical attenuation fields that enter the WAPE phase screen.

Target coupling:

\[
U_{10}, a, x, y, z, f
\rightarrow
n(a,x,y,z;U_{10})
\rightarrow
k_{\rm eff}(x,y,z,f)
\rightarrow
c_{\rm eff}(x,y,z,f),\alpha_{\rm bub}(x,y,z,f)
\rightarrow
\text{WAPE phase screen}.
\]

Here:

- \(U_{10}\): wind speed at 10 m height, m/s.
- \(a\): bubble radius, m.
- \(z\): depth below sea surface, m. Sea surface is \(z=0\), positive \(z\) points downward.
- \(f\): acoustic frequency, Hz.
- \(n(a,x,y,z;U_{10})\): bubble number density per unit volume per radius interval.
- \(c_{\rm eff}\): effective phase speed in bubbly water, m/s.
- \(\alpha_{\rm bub}\): bubble-induced physical attenuation, Np/m.

The existing numerical sponge absorption \(\alpha_{\rm sponge}(x,y)\) is a boundary treatment. Bubble attenuation is a physical medium loss. They must be added, not substituted.

\[
\alpha_{\rm total}(x,y,z,f)
=
\alpha_{\rm sponge}(x,y)+\alpha_{\rm bub}(x,y,z,f).
\]

---

## 2. Current code context

The current project contains these relevant MATLAB files:

- `CARPE3D_vertical.m`
  - Public channel API.
  - Builds validated runtime config `cfg` from user input `paramsV`.
  - Calls `propWAPE_vertical(cfg)`.
  - Packs channel outputs into `output`.

- `propWAPE_vertical.m`
  - Core vertical WAPE propagation.
  - Builds transverse grid, frequency axis, sponge boundary, direct path, and reflected path.
  - Computes:
    - `H_direct_f`
    - `H_reflect_f`
    - `H_f = H_direct_f + H_reflect_f`
    - `h_direct = H_direct_f(idx_f_ref)`
    - `h_reflect = H_reflect_f(idx_f_ref)`
    - `h_total = H_f(idx_f_ref)`
  - Contains both:
    - the main direct-path propagation loop;
    - `local_march_field`, used by the two reflection-path propagation segments.

- `pm_surface_kirchhoff_module.m`
  - Generates Pierson--Moskowitz rough sea surface.
  - Applies Kirchhoff phase distortion to the incident surface field.
  - Should not be modified for bubble-layer work unless explicitly required later.

- `comm_main_vertical_psk.m`
  - End-to-end communication demo.
  - Calls `CARPE3D_vertical(paramsV)`.
  - Converts `H_f` to baseband response.
  - Builds baseband taps, applies MPSK, channel convolution, AWGN, MMSE equalization, BER/SER.

- `modem_psk.m`
  - MPSK modulation, demodulation, error metrics.
  - Do not modify for bubble-layer propagation work.

- `noise_inject_vertical.m`
  - AWGN or custom receiver-side noise injection.
  - Bubble effects must not be added here.

The key integration point is the WAPE phase screen in `propWAPE_vertical.m`.

Current screen logic is equivalent to:

```matlab
c_local = local_sound_speed(z_curr, cfg);
U_real = (c_local - cfg.c0) / c_local;
screen = exp(-1i * k0 * ds * (U_real - 1i * alpha_xy_work / k0));
```

After adding bubbles, this must become conceptually:

```matlab
c_bg = local_sound_speed(z_curr, cfg);

if cfg.enable_bubbles
    [c_eff_xy, alpha_bub_xy, bubble_step_meta] = ...
        bubble_environment_vertical(x, y, z_curr, f_hz, c_bg, cfg);
else
    c_eff_xy = c_bg * ones(numel(y), numel(x));
    alpha_bub_xy = zeros(numel(y), numel(x));
end

alpha_total_xy = alpha_xy_work + alpha_bub_xy;
U_real_xy = (c_eff_xy - cfg.c0) ./ c_eff_xy;

screen = exp(-1i * k0 * ds * ...
    (U_real_xy - 1i * alpha_total_xy / k0));
```

This change must be applied consistently to both:

1. the main direct-path propagation loop;
2. `local_march_field`, which is used for the reflected path:
   - `z_tx -> 0`
   - `0 -> z_rx`

If only the direct loop is updated, the reflected path will not include bubble effects, which is physically inconsistent.

---

## 3. Physical interpretation

Wind-generated bubbles affect propagation through two main mechanisms:

1. Effective sound-speed change:
\[
c_{\rm bg}(z)\rightarrow c_{\rm eff}(x,y,z,f).
\]

2. Additional physical attenuation:
\[
\alpha_{\rm sponge}(x,y)
\rightarrow
\alpha_{\rm sponge}(x,y)+\alpha_{\rm bub}(x,y,z,f).
\]

In the current vertical geometry:

- Direct path:
  - propagates from \(z_{\rm tx}\) to \(z_{\rm rx}\);
  - usually crosses the near-surface bubble layer once.

- Reflected path:
  - propagates from \(z_{\rm tx}\) to the sea surface;
  - reflects from the rough surface;
  - propagates from the sea surface to \(z_{\rm rx}\);
  - usually crosses the near-surface bubble layer twice.

Therefore, bubble influence is expected to be stronger on the reflected component than on the direct component.

The final channel remains:

\[
H(f)=H_{\rm direct}(f)+H_{\rm reflect}(f),
\]

but both components should already include bubble-layer propagation effects when `enable_bubbles=true`.

---

## 4. WAPE phase-screen formula with bubbles

The current WAPE screen without bubbles is:

\[
S_0(x,y,z,f)
=
\exp\left[
-i k_0 \Delta s
\left(
\frac{c_{\rm bg}(z)-c_0}{c_{\rm bg}(z)}
-i\frac{\alpha_{\rm sponge}(x,y)}{k_0}
\right)
\right].
\]

Here:

\[
k_0=\frac{2\pi f}{c_0}.
\]

With bubbles:

\[
S_{\rm bub}(x,y,z,f)
=
\exp\left[
-i k_0 \Delta s
\left(
\frac{c_{\rm eff}(x,y,z,f)-c_0}{c_{\rm eff}(x,y,z,f)}
-i\frac{\alpha_{\rm sponge}(x,y)+\alpha_{\rm bub}(x,y,z,f)}{k_0}
\right)
\right].
\]

Implementation rule:

- `c_eff_xy` must have size `[numel(y), numel(x)]` or be scalar-expanded to that size.
- `alpha_bub_xy` must have size `[numel(y), numel(x)]`.
- `alpha_bub_xy` must be in Np/m.
- `alpha_total_xy = alpha_xy + alpha_bub_xy`.
- When `enable_bubbles=false`, the numerical result must match the pre-bubble implementation within documented tolerance.

---

## 5. Level 0 empirical bubble model

Level 0 is the first implementation target. It is intentionally simple and is used to validate the code path, sign convention, and direct/reflection consistency before adding Hall--Novarini physics.

### 5.1 Empirical attenuation

Use a configurable near-surface attenuation profile:

\[
\alpha_{\rm bub}(z,f)
=
\alpha_0
\exp\left(-\frac{z}{L_b}\right)
\left(\frac{f}{f_{\rm ref}}\right)^{p_\alpha}.
\]

Where:

- \(\alpha_0\): surface attenuation strength, Np/m.
- \(L_b\): bubble layer e-folding depth, m.
- \(f_{\rm ref}\): reference frequency, Hz.
- \(p_\alpha\): frequency exponent.

For Level 0 one-dimensional background layer:

\[
\alpha_{\rm bub}(x,y,z,f)=\alpha_{\rm bub}(z,f).
\]

Default values should be conservative:

```matlab
enable_bubbles = false
bubble_model = 'level0_empirical'
bubble_alpha0_np_per_m = 0
bubble_layer_decay_m = 0.4
bubble_f_ref_hz = 6000
bubble_alpha_freq_exp = 0
bubble_apply_attenuation = true
```

The default disabled path must have no effect on existing results.

### 5.2 Optional empirical sound-speed perturbation

To test phase effects separately from attenuation:

\[
\Delta c_{\rm bub}(z,f)
=
\Delta c_0
\exp\left(-\frac{z}{L_c}\right)
\left(\frac{f}{f_{\rm ref}}\right)^{p_c}.
\]

Then:

\[
c_{\rm eff}(z,f)
=
c_{\rm bg}(z)+\Delta c_{\rm bub}(z,f).
\]

Recommended defaults:

```matlab
bubble_delta_c0_mps = 0
bubble_sound_speed_decay_m = 0.4
bubble_sound_speed_freq_exp = 0
bubble_apply_sound_speed = true
```

Level 0 may set `bubble_delta_c0_mps = 0` by default so that pure attenuation tests are easier.

### 5.3 Level 0 acceptance checks

With `enable_bubbles=false`:

- `H_f`, `H_direct_f`, `H_reflect_f`, `h_total` must match the old version within documented tolerance.
- `comm_main_vertical_psk.m` must still run.
- `direct_only` and `direct_plus_reflect` scenarios must still run.
- Existing output fields must remain present.

With `enable_bubbles=true`, pure attenuation only:

- set `bubble_alpha0_np_per_m > 0`;
- set `bubble_delta_c0_mps = 0`;
- expect roughly:
\[
|H_{\rm bub}(f)| \le |H_0(f)|
\]
for comparable direct-only tests.
- reflected-path attenuation should usually be stronger than direct-only attenuation because the reflected branch crosses the near-surface layer twice.

With sound-speed perturbation only:

- set `bubble_alpha0_np_per_m = 0`;
- set `bubble_delta_c0_mps ~= 0`;
- expect the dominant change to be phase:
\[
\Delta \phi(f)
\approx
-k_0
\int_{\Gamma}
\frac{c_{\rm eff}(s,f)-c_0}{c_{\rm eff}(s,f)}
ds.
\]

---

## 6. Level 1 Hall-type one-dimensional average bubble layer

Level 1 adds a physically motivated one-dimensional bubble-size distribution depending on radius, depth, and wind speed:

\[
n_{\rm bg}(a,z;U_{10}).
\]

This should still be horizontally uniform:

\[
n(a,x,y,z;U_{10})=n_{\rm bg}(a,z;U_{10}).
\]

### 6.1 Configurable general form

A general separable model may be used:

\[
n_{\rm bg}(a,z;U_{10})
=
A_b(U_{10})F(a)D(z;U_{10}).
\]

Wind-strength factor:

\[
A_b(U_{10})
=
A_{\rm ref}
\left(\frac{U_{10}}{U_{\rm ref}}\right)^{p_U}.
\]

Radius spectrum:

\[
F(a)
=
C_a a^{-p_a},
\qquad
a_{\min}\le a\le a_{\max}.
\]

Depth profile:

\[
D(z;U_{10})
=
\exp\left[-\frac{z}{L_b(U_{10})}\right].
\]

Use this form if a simplified calibrated model is needed.

### 6.2 Hall--Novarini-type empirical form

A closer Hall--Novarini-type model may be implemented as:

\[
\mathcal N(a,z,u_{10})
=
p_0D(z,u_{10})G(a,z)
\left(\frac{u_{10}}{13}\right)^3.
\]

Use:

\[
p_0=1.6\times 10^{10}\ {\rm m^{-4}}.
\]

Depth decay:

\[
D(z,u_{10})
=
\exp\left[-\frac{z}{L(u_{10})}\right],
\]

with:

\[
L(u_{10})=
\begin{cases}
0.4, & u_{10}\le 7.5,\\
0.4+0.115(u_{10}-7.5), & u_{10}>7.5.
\end{cases}
\]

Radius spectrum:

\[
G(a,z)=
\begin{cases}
0, & a<10\ \mu{\rm m},\\
\left(\frac{a_{\rm ref}(z)}{a}\right)^4,
& 10\ \mu{\rm m}\le a\le a_{\rm ref}(z),\\
\left(\frac{a_{\rm ref}(z)}{a}\right)^{\chi(z)},
& a_{\rm ref}(z)<a\le1000\ \mu{\rm m},\\
0, & a>1000\ \mu{\rm m}.
\end{cases}
\]

where:

\[
a_{\rm ref}(z)=54.4+1.984\times10^{-6}z
\quad (\mu{\rm m}),
\]

\[
\chi(z)=4.37+\left(\frac{z}{2.55}\right)^2.
\]

Important implementation note:

- MATLAB internal radius grid should be in meters.
- Convert micrometer thresholds to meters:
  - \(10\ \mu{\rm m}=10^{-5}\ {\rm m}\)
  - \(1000\ \mu{\rm m}=10^{-3}\ {\rm m}\)
- Recommended radius grid:
```matlab
bubble_radius_grid_m = logspace(-5, -3, 80)
```

### 6.3 Void fraction check

The bubble void fraction is:

\[
\beta(z)
=
\int_{a_{\min}}^{a_{\max}}
\frac{4\pi}{3}a^3 n_{\rm bg}(a,z;U_{10})\,da.
\]

For spatially varying plume cases:

\[
\beta(x,y,z)
=
\int_{a_{\min}}^{a_{\max}}
\frac{4\pi}{3}a^3 n(a,x,y,z;U_{10})\,da.
\]

Use `bubble_beta_max` as a safety cap. If computed void fraction exceeds the cap, scale the local bubble number density down so that:

\[
0\le \beta \le \beta_{\max}.
\]

Recommended initial default:

```matlab
bubble_beta_max = 1e-3
```

This value is a numerical safety limit, not a claim that the physical environment always reaches it.

---

## 7. Resonance radius and current 4--8 kHz frequency band

Bubble resonance is frequency-dependent. A simple resonance radius estimate is:

\[
a_{\rm res}(f,z)
=
\frac{1}{2\pi f}
\sqrt{
\frac{3\gamma P_0(z)}{\rho_w}
}.
\]

Hydrostatic pressure:

\[
P_0(z)
=
P_{\rm atm}+\rho_wgz.
\]

Suggested constants:

```matlab
bubble_gamma = 1.4
bubble_rho_w_kg_m3 = 1025
bubble_P_atm_pa = 101325
bubble_g_m_s2 = 9.81
```

For the current \(4\sim8\ {\rm kHz}\) band, the near-surface resonance radius is roughly \(0.4\sim0.9\ {\rm mm}\), which is close to the upper end of the Hall-type \(10\ \mu{\rm m}\) to \(1000\ \mu{\rm m}\) radius range.

Implementation implication:

- Do not use a radius grid that only contains small bubbles.
- Use at least:
\[
10\ \mu{\rm m}\le a\le1000\ \mu{\rm m}.
\]
- If future frequencies go below 4 kHz, the resonance radius may exceed 1 mm, and the model should warn that the Hall upper radius may underpredict resonance loss.

---

## 8. Effective medium conversion

The bubble model must convert the bubble spectrum to either complex sound speed or complex wavenumber.

### 8.1 Complex sound speed form

One usable form is:

\[
\frac{1}{\tilde c_{bb}^2(z,u_{10})}
=
\frac{1}{c_0^2(z)}
+
\frac{1}{\pi f^2}
\int
\frac{
a\mathcal N(a,z,u_{10})
}{
\left(\frac{a_r}{a}\right)^2-1+j d
}
da.
\]

Where:

- \(\tilde c_{bb}\): complex sound speed.
- \(c_0(z)\): background sound speed.
- \(a_r\): resonance radius at \((f,z)\).
- \(d\): damping term.

Complex slowness:

\[
q(z,f)
=
\frac{1}{\tilde c_{bb}(z,f)}
=
q_R(z,f)+jq_I(z,f).
\]

Then:

\[
c_{\rm eff}(z,f)=\frac{1}{q_R(z,f)}.
\]

\[
\alpha_{\rm bub}(z,f)=\omega q_I(z,f),
\qquad
\omega=2\pi f.
\]

If numerical sign conventions produce negative attenuation, use:

\[
\alpha_{\rm bub}(z,f)=\max(\omega q_I(z,f),0).
\]

But this must be validated with a simple attenuation-only test.

### 8.2 Complex wavenumber form

A PE-oriented alternative is:

\[
k_{\rm eff}^2(z,\omega)
=
k_w^2(z,\omega)
+
4\pi
\int_{a_{\min}}^{a_{\max}}
\frac{
a n_{\rm bg}(a,z;U_{10})
}{
\left(\frac{\omega_0^2(a,z)}{\omega^2}-1\right)-i\delta(a,\omega,z)
}
da.
\]

where:

\[
k_w(z,\omega)=\frac{\omega}{c_{\rm bg}(z)}.
\]

After computing:

\[
k_{\rm eff}=k_r+i k_i,
\]

take:

\[
c_{\rm eff}(z,f)=\frac{\omega}{k_r(z,f)},
\]

\[
\alpha_{\rm bub}(z,f)=\max(k_i(z,f),0).
\]

### 8.3 Damping model

For the first Level 1 implementation, use a configurable constant damping term:

```matlab
bubble_damping_model = 'constant'
bubble_delta_const = 0.1
```

Later versions may split damping into thermal, viscous, and radiation components. Do not implement advanced damping until the Level 1 constant-damping model passes validation.

---

## 9. Level 2 horizontally nonuniform plume bubble clouds

Level 2 adds horizontal nonuniformity on top of the one-dimensional background layer:

\[
n(a,x,y,z;U_{10})
=
n_{\rm bg}(a,z;U_{10})M(x,y,z).
\]

Plume modulation:

\[
M(x,y,z)
=
1+
\sum_{m=1}^{N_p}
Q_m
\exp\left[
-\frac{(x-x_m)^2+(y-y_m)^2}{2\sigma_m^2}
\right]
\exp\left[-\frac{z}{L_{p,m}}\right].
\]

Where:

- \(N_p\): number of plume patches.
- \((x_m,y_m)\): plume center.
- \(Q_m\): plume strength.
- \(\sigma_m\): horizontal plume scale.
- \(L_{p,m}\): vertical plume decay scale.

Implementation rules:

- Plume generation must be deterministic for fixed `bubble_seed`.
- Plume centers must lie inside the transverse domain.
- Plume concentration must be clipped or rescaled to satisfy `bubble_beta_max`.
- Plume mode must be disabled by default.
- Do not compute full radius integrals at every `(x,y,z,f)` point if avoidable. Prefer:
  1. compute 1D background effective medium;
  2. use plume strength as a concentration scale;
  3. approximate or table-map concentration scaling to `c_eff` and `alpha_bub`.

Suggested fields:

```matlab
bubble_spatial_mode = '1d'       % 'none', '1d', 'plume'
bubble_plume_count = 0
bubble_plume_seed = 12345
bubble_plume_strength = []
bubble_plume_sigma_m = []
bubble_plume_decay_m = []
bubble_plume_centers_xy = []
```

---

## 10. Proposed new files

Prefer adding new functions instead of inserting all bubble logic into `propWAPE_vertical.m`.

Recommended new files:

### 10.1 `bubble_environment_vertical.m`

Main public bubble environment function.

Suggested signature:

```matlab
function [c_eff_xy, alpha_bub_xy, meta] = ...
    bubble_environment_vertical(x, y, z_curr, f_hz, c_bg, cfg)
```

Responsibilities:

- If disabled:
  - return `c_eff_xy = c_bg * ones(numel(y), numel(x))`;
  - return `alpha_bub_xy = zeros(numel(y), numel(x))`;
  - return `meta.enabled=false`.

- For `level0_empirical`:
  - compute empirical `delta_c` and `alpha_bub`;
  - expand to `[numel(y), numel(x)]`.

- For `hall1d`:
  - call spectrum and effective-medium helpers;
  - expand one-dimensional result to `[numel(y), numel(x)]`.

- For `plume`:
  - apply horizontal plume modulation.

- Must preserve CPU/GPU compatibility:
  - It may return CPU arrays and allow caller to convert to `gpuArray`, or it may detect caller state.
  - Do not mix CPU/GPU arrays in phase-screen multiplication.

### 10.2 `bubble_hall_spectrum.m`

Suggested signature:

```matlab
function [n_a, beta, spec_meta] = ...
    bubble_hall_spectrum(a_grid_m, z_curr, U10, cfg)
```

Responsibilities:

- Compute Hall-type number density on radius grid.
- Compute void fraction.
- Apply `bubble_strength_scale`.
- Enforce or report `bubble_beta_max`.

### 10.3 `bubble_effective_medium.m`

Suggested signature:

```matlab
function [c_eff, alpha_bub, em_meta] = ...
    bubble_effective_medium(a_grid_m, n_a, z_curr, f_hz, c_bg, cfg)
```

Responsibilities:

- Convert `n_a` to complex sound speed or complex wavenumber.
- Return scalar `c_eff` and `alpha_bub` for 1D layer.
- Ensure `alpha_bub >= 0`.
- Return diagnostic quantities:
  - resonance radius;
  - maximum number density;
  - beta;
  - damping model;
  - raw complex quantity.

### 10.4 `bubble_plume_mask.m`

Suggested signature:

```matlab
function [M_xy, plume_meta] = bubble_plume_mask(x, y, z_curr, cfg)
```

Responsibilities:

- Generate or evaluate deterministic plume modulation.
- Return `[numel(y), numel(x)]` field.
- Ensure no negative concentration.

### 10.5 `bubble_field_stats.m`

Suggested signature:

```matlab
function stats = bubble_field_stats(v)
```

Responsibilities:

- Return `min`, `max`, `mean`, `std`, `finite_count`.
- Used for `output.bubble_meta`.

---

## 11. Required changes to existing files

### 11.1 `CARPE3D_vertical.m`

Add defaults to `local_prepare_config`.

Recommended fields:

```matlab
'enable_bubbles', false, ...
'bubble_model', 'off', ...
'bubble_spatial_mode', 'none', ...
'bubble_apply_sound_speed', true, ...
'bubble_apply_attenuation', true, ...
'bubble_alpha0_np_per_m', 0, ...
'bubble_layer_decay_m', 0.4, ...
'bubble_f_ref_hz', 6000, ...
'bubble_alpha_freq_exp', 0, ...
'bubble_delta_c0_mps', 0, ...
'bubble_sound_speed_decay_m', 0.4, ...
'bubble_sound_speed_freq_exp', 0, ...
'bubble_radius_grid_m', logspace(-5, -3, 80), ...
'bubble_strength_scale', 1.0, ...
'bubble_beta_max', 1e-3, ...
'bubble_damping_model', 'constant', ...
'bubble_delta_const', 0.1, ...
'bubble_gamma', 1.4, ...
'bubble_rho_w_kg_m3', 1025, ...
'bubble_P_atm_pa', 101325, ...
'bubble_g_m_s2', 9.81, ...
'bubble_seed', 12345, ...
'bubble_plume_count', 0, ...
'bubble_plume_strength', [], ...
'bubble_plume_sigma_m', [], ...
'bubble_plume_decay_m', [], ...
'bubble_plume_centers_xy', []
```

Validation requirements:

- `enable_bubbles` must be scalar logical.
- If `enable_bubbles=false`, force or allow `bubble_model='off'`.
- `bubble_model` allowed values:
  - `'off'`
  - `'level0_empirical'`
  - `'hall1d'`
  - `'plume'`
- `bubble_spatial_mode` allowed values:
  - `'none'`
  - `'1d'`
  - `'plume'`
- `bubble_radius_grid_m` must be positive, finite, increasing.
- `bubble_beta_max` must be positive.
- `bubble_alpha0_np_per_m` must be nonnegative.
- Decay depths must be positive.
- Damping constants must be positive.
- Plume fields must be validated only when plume mode is enabled.

Add to `output`:

```matlab
output.bubble_meta = bubble_meta;
```

If no bubble metadata is produced, return:

```matlab
output.bubble_meta = struct('enabled', false, 'model', 'off');
```

### 11.2 `propWAPE_vertical.m`

Required changes:

1. Add `bubble_meta` to the output list from `propWAPE_vertical`, or pack it inside an existing returned structure only if interface impact is controlled.
2. In the main direct-path propagation loop:
   - replace scalar `c_local/U_real/screen` logic with bubble-aware screen logic.
3. Refactor or update `local_march_field` so it can also call `bubble_environment_vertical`.
4. Ensure `local_march_field` receives `x`, `y`, and `f_hz`, because the bubble environment needs them.

Suggested `local_march_field` signature:

```matlab
function [psi_end, march_meta] = local_march_field( ...
    psi_start, z_start, z_end, cfg, x, y, alpha_xy, kappa2, f_hz, use_gpu)
```

Inside it, derive:

```matlab
lambda_f = cfg.c0 / f_hz;
k0 = 2*pi / lambda_f;
```

or pass `k0` as before and also pass `f_hz`.

Do not update only one propagation branch.

### 11.3 `comm_main_vertical_psk.m`

Do not change the default communication behavior in the first bubble integration commit.

Later, add optional scenarios:

- `direct_only_no_bubble`
- `direct_only_bubble`
- `direct_plus_reflect_no_bubble`
- `direct_plus_reflect_bubble`

But this should be a later stage, after propagation integration is validated.

### 11.4 `pm_surface_kirchhoff_module.m`

Do not modify for Level 0 or Level 1 unless explicitly necessary.

Bubbles act during propagation to and from the sea surface, not inside the Kirchhoff reflection formula.

### 11.5 `noise_inject_vertical.m`

Do not add bubble effects here.

Bubble attenuation and phase effects are channel propagation effects, not receiver noise.

---

## 12. Metadata requirements

Add diagnostic metadata to help validate and plot results.

Recommended structure:

```matlab
output.bubble_meta.enabled
output.bubble_meta.model
output.bubble_meta.spatial_mode
output.bubble_meta.apply_sound_speed
output.bubble_meta.apply_attenuation
output.bubble_meta.layer_decay_m
output.bubble_meta.radius_grid_m
output.bubble_meta.beta_max
output.bubble_meta.seed
output.bubble_meta.c_eff_stats
output.bubble_meta.alpha_bub_stats
output.bubble_meta.beta_stats
output.bubble_meta.resonance_radius_stats
output.bubble_meta.warning_flags
```

For frequency-dependent results, optionally include compact arrays:

```matlab
output.bubble_meta.f_axis
output.bubble_meta.alpha_ref_f
output.bubble_meta.c_eff_ref_f
output.bubble_meta.beta_ref_f
```

Avoid storing full 3D/4D fields by default, because memory use may become large. Full fields should only be stored when a diagnostic flag is enabled, for example:

```matlab
bubble_save_diagnostics = false
```

---

## 13. Priority and staged implementation plan

### Priority 0: Read project and preserve baseline

Task:

- Read `AGENTS.md`.
- Read the MATLAB files.
- Do not use external documents unless explicitly requested.
- Identify current propagation screen logic and reflection `local_march_field` usage.

Acceptance:

- No code changes.
- Codex reports current integration points.
- Codex identifies that both direct loop and `local_march_field` must be updated.

---

### Priority 1: Add disabled-by-default bubble framework

Task:

- Add bubble config fields to `CARPE3D_vertical.m`.
- Add `bubble_environment_vertical.m`.
- Add output `bubble_meta`.
- Keep default behavior equivalent to old code:
  - `enable_bubbles=false`
  - `bubble_model='off'`
- Wire `bubble_environment_vertical` into both:
  - main direct-path loop;
  - `local_march_field`.
- Before modifying propagation code, run and save a reduced-grid baseline for `enable_bubbles=false` so that later disabled-path regression can be compared numerically.

Acceptance:

- With `enable_bubbles=false`, scalar-frequency direct-only result matches previous baseline within tolerance.
- With `enable_bubbles=false`, direct-plus-reflect result matches previous baseline within tolerance.
- `comm_main_vertical_psk.m` still runs.
- `H_f = H_direct_f + H_reflect_f` remains true.
- Existing output fields remain present.
- `output.bubble_meta.enabled=false`.

Recommended tolerance:

```matlab
relative_error = norm(H_new - H_old) / max(norm(H_old), eps)
```

For CPU deterministic tests, target:

```matlab
relative_error < 1e-10
```

If refactoring changes floating-point order, document any slightly larger tolerance.

---

### Priority 2: Implement Level 0 empirical bubble layer

Task:

- Implement `level0_empirical` in `bubble_environment_vertical`.
- Support:
  - empirical attenuation profile;
  - optional empirical sound-speed perturbation;
  - independent switches for attenuation and sound speed.

Acceptance:

- Pure attenuation test:
  - `bubble_alpha0_np_per_m > 0`
  - `bubble_delta_c0_mps = 0`
  - direct-only `|H_bub| <= |H_0|` in most cases.
- Pure sound-speed test:
  - `bubble_alpha0_np_per_m = 0`
  - `bubble_delta_c0_mps ~= 0`
  - amplitude changes remain small relative to phase changes.
- Reflected path shows stronger bubble influence than direct-only under comparable settings.
- No changes to `noise_inject_vertical.m`.

- During Priority 1 and Priority 2, do not modify `modem_psk.m`, `noise_inject_vertical.m`, or the BER/SER logic in `comm_main_vertical_psk.m`. Bubble effects must first be validated at the channel `H_f` level.
---

### Priority 3: Implement Level 1 Hall-type 1D average bubble layer

Task:

- Add `bubble_hall_spectrum.m`.
- Add `bubble_effective_medium.m`.
- Compute:
  - radius-dependent number density;
  - void fraction;
  - resonance radius;
  - effective sound speed;
  - attenuation.
- Use constant damping first.
- Keep model horizontally uniform.

Acceptance:

- `bubble_model='hall1d'` runs for scalar frequency and wideband frequency axis.
- Radius grid covers \(10^{-5}\) to \(10^{-3}\) m.
- `alpha_bub` is finite and nonnegative.
- `c_eff` is finite and positive.
- `beta` is finite and does not exceed `bubble_beta_max`.
- Frequency dependence is visible in `alpha_bub(f)` or `c_eff(f)`.
- Results remain deterministic for fixed seeds and parameters.
- Direct/reflected channel components remain consistent:
  \[
  H_f = H_{\rm direct,f}+H_{\rm reflect,f}.
  \]

---

### Priority 4: Add diagnostics and comparison scripts

Task:

- Add optional diagnostic plotting or script.
- Compare:
  - no bubble;
  - Level 0 empirical bubble;
  - Hall 1D bubble.
- Output:
  - \(|H(f)|\)
  - \(\angle H(f)\)
  - \(\Delta TL_{\rm bub}(f)\)
  - \(R_H(f)\)
  - BER/SER curves if using communication chain.

Key metrics:

\[
\frac{|H_{\rm bub}(f)|}{|H_0(f)|}.
\]

\[
\Delta TL_{\rm bub}(f)
=
-20\log_{10}
\left(
\frac{|H_{\rm bub}(f)|}{|H_0(f)|}
\right).
\]

\[
R_H(f)=\frac{H_{\rm bub}(f)}{H_0(f)}.
\]

Acceptance:

- Diagnostics run at reduced grid size.
- No large default memory increase.
- Figures and saved results are optional and controlled by flags.

---

### Priority 5: Add plume nonuniformity

Task:

- Add `bubble_plume_mask.m`.
- Implement deterministic horizontal plume modulation.
- Add `bubble_spatial_mode='plume'`.
- Enforce void-fraction cap.

Acceptance:

- `bubble_spatial_mode='1d'` still works unchanged.
- `bubble_spatial_mode='plume'` produces horizontally varying `alpha_bub_xy` and/or `c_eff_xy`.
- Plume result is deterministic for fixed `bubble_seed`.
- No negative bubble densities or negative attenuation.
- Memory use remains controlled.

---

### Priority 6: Integrate bubble comparison into communication workflow

Task:

- Add optional bubble/no-bubble communication scenarios.
- Compare BER/SER and effective baseband taps.
- Do not change default old scenarios unless explicitly requested.

Acceptance:

- Existing `comm_main_vertical_psk.m` scenario structure still works.
- New scenarios can compare no-bubble vs bubble cases.
- BER/SER outputs remain structured and easy to save.
- `results(ss).channel.bubble_meta` is available for each scenario.

---

## 14. Validation checklist

After any nontrivial edit:

1. Run a reduced scalar-frequency channel test:
```matlab
paramsV = struct();
paramsV.f0 = 4000;
paramsV.enable_wideband = false;
paramsV.c0 = 1500;
paramsV.z_max = 100;
paramsV.z_tx = 100;
paramsV.z_rx = 3;
paramsV.xw = 50;
paramsV.yw = 50;
paramsV.nx = 128;
paramsV.ny = 128;
paramsV.show_figures = false;
paramsV.save_mode = 'rx_only';
paramsV.enable_surface_reflection = false;
paramsV.enable_bubbles = false;
out = CARPE3D_vertical(paramsV);
```

2. Run direct-plus-reflect with bubbles disabled:
```matlab
paramsV.enable_surface_reflection = true;
paramsV.enable_bubbles = false;
out = CARPE3D_vertical(paramsV);
```

3. Run direct-only Level 0 attenuation:
```matlab
paramsV.enable_surface_reflection = false;
paramsV.enable_bubbles = true;
paramsV.bubble_model = 'level0_empirical';
paramsV.bubble_alpha0_np_per_m = 0.02;
paramsV.bubble_delta_c0_mps = 0;
out_bub = CARPE3D_vertical(paramsV);
```

4. Run direct-plus-reflect Level 0 attenuation:
```matlab
paramsV.enable_surface_reflection = true;
out_bub_ref = CARPE3D_vertical(paramsV);
```

5. Run reduced wideband communication test only after channel tests pass.

---

## 15. Common mistakes to avoid

- Do not add bubbles in `noise_inject_vertical.m`.
- Do not implement bubbles as a scalar multiplier applied only to `h_total`.
- Do not update only the direct propagation loop.
- Do not forget `local_march_field`, because it is used for reflected path propagation.
- Do not confuse numerical sponge absorption with physical bubble attenuation.
- Do not silently change units between dB/m and Np/m.
- Do not store large full-field bubble diagnostics by default.
- Do not make plume mode the default.
- Do not modify `modem_psk.m` for propagation physics.
- Do not change public output fields without updating all dependent scripts.

---

## 16. Expected final behavior

After all stages are complete, the project should support:

1. Baseline propagation:
```matlab
paramsV.enable_bubbles = false;
```

2. Empirical bubble-layer propagation:
```matlab
paramsV.enable_bubbles = true;
paramsV.bubble_model = 'level0_empirical';
```

3. Hall-type average bubble-layer propagation:
```matlab
paramsV.enable_bubbles = true;
paramsV.bubble_model = 'hall1d';
paramsV.bubble_spatial_mode = '1d';
```

4. Horizontally nonuniform plume propagation:
```matlab
paramsV.enable_bubbles = true;
paramsV.bubble_model = 'plume';
paramsV.bubble_spatial_mode = 'plume';
```

In all enabled modes:

\[
H(f)=H_{\rm direct}(f)+H_{\rm reflect}(f)
\]

must remain valid, and each component should already include the appropriate bubble-layer propagation effect.
