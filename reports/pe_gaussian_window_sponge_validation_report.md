# Production Gaussian finite-window and sponge validation

Status: **default sponge not recommended; strict production recommendation is `W=160 m`, sponge off**. No PE marching, FFT convention, Gaussian definition, surface model, or communication modulation was changed.

## Environment and source audit

- Homogeneous `c=1500 m/s`; direct path only; surface reflection, bubbles, random surface, Doppler, and communication processing off.
- Gaussian: `Psi0=exp(-((x-x_tx)^2+(y-y_tx)^2)/(2 sigma^2))`, `sigma=0.3 m`, unit peak, no additional normalization, frequency independent, centered at `(0,0)`.
- Production geometry: `W=50 m`, `256^2`, `dx=0.1953125 m`, `L=97 m`; single-frequency 4 kHz; wideband 3--5 kHz/33 points.
- `H_f` is direct-DSP referenced; physical phase/group delay use `H_direct_physical_f` and `p(t)=Re(P exp(-i2pift))`.

## Hard check and window reference

PE--AS maximum full/center complex errors are `2.56337e-13 / 1.45326e-13`; the hard gate passed. The 128 m no-sponge field converges to the 160 m field under the preregistered single-frequency limits, so `W_ref=128 m` for the 4 kHz gate while 160 m is retained as the comparison reference.

## Production default sponge, W=50 m

Definitions: `dA=20log10(|Hsp|/|H0|)`, `dTL=-dA`, positive dTL means added loss.

| metric | no sponge | default sponge | change |
|---|---:|---:|---:|
| axis magnitude | 0.0137348876 | 0.0147549456 | dA 0.622250 dB / dTL -0.622250 dB |
| axis phase | -1.322147501 rad | -1.470702805 rad | -0.148555 rad |
| center energy, rho<=2 m | 0.0030909209 | 0.00305726086 | -0.047554 dB |
| edge energy | 0.0738198709 | 0.0389412596 | -2.777633 dB |
| total energy | 0.282743339 | 0.233849674 | -0.824556 dB |

The default changes axis TL by `-0.622250 dB` and phase by `-0.148555 rad`; it is not center-neutral. Edge suppression is only `2.777633 dB`, below the fixed 3 dB requirement.

## Boundary wrap evidence

At 97 m in the production window the outer 5%/10% energy fractions reach `0.113018 / 0.217611` and boundary amplitude reaches `-1.287934 dB` relative to the plane maximum. Thus the no-sponge Gaussian field materially reaches the periodic boundary.

## Wideband results versus 160 m/no-sponge

| case | max |dTL| dB | TL span dB | phase RMS/max rad | group-delay RMS/max us |
|---|---:|---:|---:|---:|
| large_no_sponge | 0.000000 | 0.000000 | 0.000000 / 0.000000 | 0.000000 / 0.000000 |
| production_default | 0.842225 | 1.346365 | 0.131866 / 0.285952 | 232.649327 / 1326.814293 |
| recommended_no_sponge | 0.100293 | 0.168946 | 0.003286 / 0.011962 | 8.635877 / 33.557701 |

The 128 m no-sponge candidate misses the strict wideband `|dTL|<0.1 dB` target by `0.000292612 dB`; it is therefore a near-threshold cost compromise, not the strict recommendation.

## Direct answers

1. **Q1:** 128 m is the first tested 4 kHz window satisfying the single-frequency convergence gate; 160 m is used as the strict wideband reference/recommendation.
2. **Q2:** The production 50 m window is not large enough at L=97 m. No-sponge axis TL/phase errors versus 160 m are `1.074618 dB / 0.232398 rad`.
3. **Q3:** Yes. The final boundary is only `-1.288 dB` below the field maximum.
4. **Q4:** Default single-frequency and wideband errors are tabulated above; wideband max TL, phase RMS, and group-delay RMS are `0.842225 dB`, `0.131866 rad`, and `232.649 us`.
5. **Q5:** No. Default sponge does not provide >=3 dB edge suppression and significantly perturbs axis amplitude/phase.
6. **Q6:** No nonzero scanned sponge meets all fixed center and edge targets. Recommend `W=160 m`, sponge off.
7. **Q7:** Yes: expanding to 160 m with no sponge is sufficient in the tested 3--5 kHz, 97 m direct-path environment.

## Recommendation

- Strict: `window=160 m`, `sponge off (alpha_max=0)`; ratio is inactive. Relative errors are zero by reference definition.
- Cost compromise: `window=128 m`, sponge off; max wideband TL `0.100293 dB`, phase RMS `0.003286 rad`, group-delay RMS `8.635877 us`.
- Do **not** change production defaults automatically yet. The evidence recommends a future reviewed config change from 50 m/default sponge to a larger no-sponge window, followed by communication-cost and full reflected-path qualification.

## Optional LFM probe

Normalized waveform correlations to the 160 m reference are production/default `0.997616` and 128 m/no-sponge `0.999997`; matched-filter peak delays are `[0 0 0] s`. This is engineering interpretation only.

## Figures

![figure 1](../results/validation/pe_gaussian_window_sponge/figures/01_gaussian_initial_field.png)

![figure 2](../results/validation/pe_gaussian_window_sponge/figures/02_no_sponge_center_amplitude_profiles.png)

![figure 3](../results/validation/pe_gaussian_window_sponge/figures/03_no_sponge_center_phase_profiles.png)

![figure 4](../results/validation/pe_gaussian_window_sponge/figures/04_window_convergence_vs_width.png)

![figure 5](../results/validation/pe_gaussian_window_sponge/figures/05_edge_energy_vs_distance.png)

![figure 6](../results/validation/pe_gaussian_window_sponge/figures/06_center_energy_vs_distance.png)

![figure 7](../results/validation/pe_gaussian_window_sponge/figures/07_default_sponge_on_off_terminal_fields.png)

![figure 8](../results/validation/pe_gaussian_window_sponge/figures/08_default_sponge_amplitude_difference.png)

![figure 9](../results/validation/pe_gaussian_window_sponge/figures/09_default_sponge_phase_difference.png)

![figure 10](../results/validation/pe_gaussian_window_sponge/figures/10_axis_tl_change_heatmap.png)

![figure 11](../results/validation/pe_gaussian_window_sponge/figures/11_center_energy_change_heatmap.png)

![figure 12](../results/validation/pe_gaussian_window_sponge/figures/12_axis_phase_change_heatmap.png)

![figure 13](../results/validation/pe_gaussian_window_sponge/figures/13_production_default_vs_reference_Hf.png)

![figure 14](../results/validation/pe_gaussian_window_sponge/figures/14_recommended_vs_reference_Hf.png)

![figure 15](../results/validation/pe_gaussian_window_sponge/figures/15_group_delay_comparison.png)

![optional LFM](../results/validation/pe_gaussian_window_sponge/figures/16_optional_lfm_matched_filter.png)
