# Full Kirchhoff 2-D sign / convention audit

Status: **CONVENTION BUG CONFIRMED AND FIXED IN VALIDATION DRIVER**. Production PE, BIE solver, Full-Kirchhoff kernel, and surface model are unchanged.

## Frozen conventions

- Coordinates: `z` is positive downward; water is `z > eta(x)`; the authoritative physical boundary is `z_s=+eta(x)`.
- Water normal: `n=(-eta'(x),1)/sqrt(1+eta'(x)^2)`. Both Case A (`+eta`) and the diagnostic Case B (`-eta`) recompute position, slope, normal, incident field, and incident normal derivative as one consistent geometry.
- Time/Green function: `exp(-i omega t)` and outgoing `G=(i/4)H_0^(1)(kR)`. The implemented source derivative is consistent with this Green function and the water normal.
- Field definition: BIE `receiver_field` is reflected/scattered field only. The saved `G_BIE` is `conj(G_BIE_native)`. The original Full-Kirchhoff driver formed a native rough/flat ratio but compared it directly with this already-conjugated field.

Flat control: complex L2 0.000185546, magnitude L2 0.000129482, phase RMS 0.000180854 rad; **PASS**.

## Geometry and field-convention audit

| case | geometry | comparison | complex L2 | magnitude L2 | phase RMS (rad) | magnitude corr. | phase corr. |
|---|---|---|---:|---:|---:|---:|---:|
| weak_low_K | A: z_s=+eta | legacy native-vs-comparison | 0.455442 | 0.000177857 | 0.46168 | 0.814297 | 0.896287 |
| weak_low_K | A: z_s=+eta | fixed comparison-vs-comparison | 0.00025349 | 0.000177857 | 0.000180615 | 0.814297 | 1 |
| weak_low_K | A: z_s=+eta | native-vs-native control | 0.00025349 | 0.000177857 | 0.000180615 | 0.814297 | 1 |
| weak_low_K | B: z_s=-eta | fixed comparison-vs-comparison | 0.455442 | 0.000529814 | 0.46168 | 0.813817 | 0.896287 |
| strong_height_low_K | A: z_s=+eta | legacy native-vs-comparison | 1.25819 | 0.000176535 | 1.60046 | 0.999377 | 0.208384 |
| strong_height_low_K | A: z_s=+eta | fixed comparison-vs-comparison | 0.000257001 | 0.000176535 | 0.000186802 | 0.999377 | 1 |
| strong_height_low_K | A: z_s=+eta | native-vs-native control | 0.000257001 | 0.000176535 | 0.000186802 | 0.999377 | 1 |
| strong_height_low_K | B: z_s=-eta | fixed comparison-vs-comparison | 1.25819 | 0.00998604 | 1.59973 | 0.995387 | 0.208382 |
| weak_high_K | A: z_s=+eta | legacy native-vs-comparison | 0.889026 | 0.000177236 | 0.939647 | 0.999822 | 0.604936 |
| weak_high_K | A: z_s=+eta | fixed comparison-vs-comparison | 0.000271132 | 0.000177236 | 0.000205163 | 0.999822 | 1 |
| weak_high_K | A: z_s=+eta | native-vs-native control | 0.000271132 | 0.000177236 | 0.000205163 | 0.999822 | 1 |
| weak_high_K | B: z_s=-eta | fixed comparison-vs-comparison | 0.889129 | 0.0187594 | 0.939586 | 0.998427 | 0.604946 |

Case B is intentionally a different mirrored physical surface and is not an alternative convention for the authoritative BIE case. Its mismatch confirms that the physical boundary coordinate is `z_s=+eta`, not `-eta`.

## `+/-4*k*eta` diagnostic

| case | RMS[legacy delta - (+4keta)] | RMS[legacy delta - (-4keta)] | corrected delta RMS |
|---|---:|---:|---:|
| weak_low_K | 0.928993 | 0.0183667 | 0.000180615 |
| strong_height_low_K | 1.67812 | 0.371275 | 0.000186802 |
| weak_high_K | 1.87875 | 0.18244 | 0.000205163 |

For weak/low-K the old mixed-convention phase difference follows `-4*k*eta` closely. This is exactly the expected signature when a rough/flat ratio is compared with its conjugated counterpart. The fixed comparison removes that signature without fitting or subtracting a phase.

## Direct answers

**A. Cause.** The large rough-case phase error was a field-definition/comparison-convention mismatch: native Full Kirchhoff was compared with conjugated BIE. It was not caused by surface-z, normal orientation, Green derivative, or incident normal derivative.

**B. `4*k*eta` signature.** Yes. The weak/low-K legacy residual is quantitatively close to `-4*k*eta`; the table and wrapped/unwrapped figures record the comparison.

**C. Corrected agreement.** After the single frozen native-to-PE comparison mapping, Full Kirchhoff approaches BIE at the integration-error scale for all three controlled cases. No scalar fit, amplitude normalization, or per-case sign selection is used.

**D. Physics interpretation.** The convention anomaly is excluded. For these controlled cases, Full Kirchhoff being close to BIE indicates that the earlier failure was not evidence that the Kirchhoff approximation itself is insufficient; the residual of the local PE phase-screen reduction remains the relevant model discrepancy. This conclusion is limited to the current 2-D, 4 kHz, smooth deterministic cases.

Artifacts: `results/validation/pe_bie_full_kirchhoff_convention_audit/pe_bie_full_kirchhoff_convention_audit.mat`, `results/validation/pe_bie_full_kirchhoff_convention_audit/pe_bie_full_kirchhoff_convention_audit_summary.csv`, `results/validation/pe_bie_full_kirchhoff_convention_audit/pe_bie_full_kirchhoff_phase_diagnostic.csv`, `results/validation/pe_bie_full_kirchhoff_convention_audit/weak_low_K_phase_trace.csv`, and figures under `results/validation/pe_bie_full_kirchhoff_convention_audit/figures/`.
