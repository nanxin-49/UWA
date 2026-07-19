# PE Propagation Visual Atlas

Mode: `full`; reference frequency: `6000 Hz`; PE: `128 x 128`; PM: `256 x 256`.

All spatial magnitude panels use a common `[-50,0] dB` envelope reference. Phase is masked below `-40 dB`. `z` is positive downward.

| File | Meaning | Data source / normalization |
|---|---|---|
| `01_geometry_and_paths.png` | Geometry, directions, windows, sponge, transmitter and receiver. | Formal cache/configuration and public scalar diagnostic. |
| `02_pe_xz_propagation.png` | Direct, Tx-to-surface incident, and surface-to-Rx reflected PE slices. | Public `Axz/psifinal_xy/surface_wavefield_meta` at 6 kHz. |
| `03_surface_boundary_fields.png` | Raw-PM elevation and incident/coherent/reflected/scatter surface fields. | Cached incidence, explicit raw-PM surface, and one joint-kstat realization. |
| `04_surface_angular_spectra.png` | Incident/reflected/scatter angular spectra on a shared scale. | Unshifted stored fields; spectra derived with centered `fft2` only for display. |
| `05_physical_boundary_branches.png` | Flat, explicit Kirchhoff, joint-kstat, and SSA1 illustrative fields. | Same incident field where defined; SSA1 public scalar reference. |
| `06_receiver_plane_components.png` | Direct/coherent/scatter/reflected/total receiver-depth planes. | Shared cached surface-to-receiver forward primitive. |
| `07_receiver_transverse_profiles.png` | Receiver-centred x/y profiles. | Cuts through the same receiver planes and nearest-grid sample. |
| `08_cached_vs_adjoint.png` | Same-input cached PE and exact adjoint projection plus sensitivity kernels. | Cached executor and `q=A^H r`, `a=conj(psi_inc).*q`. |
| `09_receiver_frequency_response.png` | Reflected-only response and group delay. | Same explicit/joint samples plus accepted analytic scatter RMS. |
| `10_receiver_cp_heatmaps.png` | Analytic/sample covariance and pseudo-covariance. | Accepted adjoint F=64 validation; common absolute `|C|` reference. |
| `11_pdp_eigenspectrum_distribution.png` | PDP, covariance modes, IQ and amplitude distribution. | Accepted analytic/projected PDP and deterministic conditional U=5 samples. |
| `12_u5_u8_statistical_contrast.png` | Validated U=5/U=8 receiver statistics; no wind interpolation. | Saved full conditional models and validation results. |
| `13_lfm_and_matched_filter.png` | Noiseless channel-level LFM and matched-filter diagnostics. | Prior public LFM validation plus current accepted joint analytic/sample statistics. |
| `14_pe_carrier_reconstruction.mp4` | Nominal-c0 carrier reconstruction; not time-domain PE. | Stored complex-envelope center slices multiplied by carrier phase. |

Explicit and joint single realizations are illustrative and are not compared pointwise. The adjoint kernel is receiver sensitivity, not reciprocal or inverse propagation. Analytic FFT C/P and conditional models are receiver-statistics paths, not spatial-field generators. Total-channel panels are contextual; validation metrics remain reflected-scatter based.
