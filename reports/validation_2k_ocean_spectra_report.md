# Validation of 2k phase with PM-type, JONSWAP, and TMA spectra

## 1. Goal

Extend the Gaussian 2k-phase benchmark to PM-type, JONSWAP, and TMA ocean-wave spectral shapes. This report studies only the replacement of the roughness height-phase argument `(gamma_i + gamma_s) eta` by `2k eta`. The lowest-order SSA calculation is a diagnostic, not an exact solution. No PE, KStat, directional spreading, random realization, or Monte Carlo calculation is used. The PM/JONSWAP/TMA spectra are treated as omnidirectional radial spectra in this first comparison; directional spreading is not included.

## 2. Surface statistics

| spectrum | Hs (m) | Tp (s) | sigma_eta (m) | s_rms | L_eff (m) | m0 | water depth (m) |
|---|---:|---:|---:|---:|---:|---:|---:|
| PM-type | 0.200 | 4.00 | 0.05000 | 0.03856 | 1.29678 | 0.0025 | 20.0 |
| JONSWAP | 0.200 | 4.00 | 0.05000 | 0.03215 | 1.55508 | 0.0025 | 20.0 |
| TMA | 0.200 | 4.00 | 0.05000 | 0.03219 | 1.55351 | 0.0025 | 20.0 |

The largest combined slope is 0.03856 for PM-type; the smallest effective scale is 1.29678 m for PM-type. `L_eff` is only a convenient statistic defined by sigma_eta/s_rms, not a unique correlation length.

## 3. 2k validation results

| spectrum | frequency (kHz) | kh | s_rms | theta95 full | theta99 full | theta95 2k | TV | weighted phase RMS | F(<0.1 rad) | F(<0.3 rad) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| PM-type | 4.0 | 0.8378 | 0.03856 | 15.619 | 56.218 | 12.160 | 0.011896 | 0.060284 | 0.96718 | 0.98646 |
| PM-type | 6.0 | 1.2566 | 0.03856 | 8.776 | 29.638 | 8.333 | 0.007751 | 0.061487 | 0.98784 | 0.99207 |
| PM-type | 8.0 | 1.6755 | 0.03856 | 8.346 | 15.408 | 8.370 | 0.002479 | 0.037961 | 0.99508 | 0.99873 |
| JONSWAP | 4.0 | 0.8378 | 0.03215 | 6.977 | 20.947 | 6.861 | 0.002491 | 0.027148 | 0.99262 | 0.99740 |
| JONSWAP | 6.0 | 1.2566 | 0.03215 | 6.720 | 13.493 | 6.725 | 0.001847 | 0.025522 | 0.99672 | 0.99871 |
| JONSWAP | 8.0 | 1.6755 | 0.03215 | 6.782 | 12.527 | 6.833 | 0.001568 | 0.017706 | 0.99754 | 0.99982 |
| TMA | 4.0 | 0.8378 | 0.03219 | 7.023 | 22.841 | 6.890 | 0.002665 | 0.028039 | 0.99218 | 0.99724 |
| TMA | 6.0 | 1.2566 | 0.03219 | 6.733 | 13.559 | 6.733 | 0.001878 | 0.026106 | 0.99662 | 0.99863 |
| TMA | 8.0 | 1.6755 | 0.03219 | 6.790 | 12.544 | 6.840 | 0.001574 | 0.018061 | 0.99752 | 0.99980 |

## 4. Comparison with the Gaussian benchmark

The previous Gaussian 80-deg reference covered approximately `kh = 0.25--2` and `s_rms <= 0.2`. In the nine ocean-spectrum cases, 9/9 are inside the kh interval and 9/9 are inside the slope interval. These flags are positional references only, not validity labels. The largest ocean kh is 1.6755, and the largest ocean slope is 0.03856.

For PM/JONSWAP, the deep-water mapping is `K=(2*pi*f)^2/g`. For TMA, `(2*pi*f)^2 = g K tanh(K d)` is solved independently at each frequency. All three frequency spectra are normalized to `m0=(Hs/4)^2`; the radial mapping and covariance checks are reported below.

The frequency grid is 0.03--2.0 Hz with 2000 points; acoustic frequencies are 4, 6, and 8 kHz; each covariance uses `rho = linspace(0,6*L_eff,6001)` and the angular grid is 0--80 deg with 801 points.

The PM-type shape is `f^(-5)*exp(-1.25*(fp/f)^4)` and is explicitly normalized to the prescribed Hs and Tp rather than tied to a wind-speed convention. JONSWAP multiplies that shape by `gamma_j^r(f)` with `gamma_j=3.3`, `sigma=0.07/0.09`. TMA multiplies the JONSWAP raw shape by `Phi(omega_h)=0.5*omega_h^2` for `omega_h<=1`, `1-0.5*(2-omega_h)^2` for `1<omega_h<2`, and 1 otherwise.

Largest frequency-to-wavenumber variance error: 6.939e-16

Largest covariance zero-lag relative error: 6.939e-16

Covariance tail ratios `C_eta(end)/C_eta(1)`: PM-type 6.597e-02, JONSWAP -1.151e-02, TMA -1.243e-02.

Largest pre-clipping negative Hankel value relative to the positive peak: 4.455e-04 (PM-type, 4.0 kHz, full branch).

This exceeds 1e-5 and is retained as a finite-rho-window numerical limitation; the negative tail is clipped for the normalized diagnostic spectra and is not interpreted as an exact positive Hankel transform.

## 5. Figures

![Radial ocean-wave spectra](../results/validation/ssa_2k_phase/ocean_wave_spectra.png)

![8 kHz angular-spectrum examples](../results/validation/ssa_2k_phase/ocean_spectra_2k_examples_8khz.png)

![TV distance summary](../results/validation/ssa_2k_phase/ocean_spectra_2k_summary.png)

## 6. Short conclusion

Across the nine cases, full-phase 95% and 99% energy angles range from 6.72--15.62 deg and 12.53--56.22 deg, respectively, over the 0--80 deg normal-incidence scattering window. The maximum TV distance is 0.011896 for PM-type at 4.0 kHz (`kh=0.8378`).

Within this deliberately small comparison, the spectra do not produce a dramatically broader full-phase tail than the previous Gaussian reference. The 2k replacement changes the normalized angular spectra by TV distances of 0.001568--0.011896, with the largest energy-weighted phase RMS 0.061487 rad. This is consistent with the earlier small-impact trend at the level of these diagnostic cases, while the result is not a validation of real-ocean scattering or a universal statement over wind speed, Hs, or Tp.

The PM/JONSWAP/TMA shapes here are one normalized Hs/Tp condition and an omnidirectional first comparison. Directional spreading, other sea states, and exact higher-order scattering remain outside the scope.
