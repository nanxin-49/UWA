# Rotated Bellhop field validation adapter

This validation-only module reads an optional `.r25d` sidecar containing:
`point3D_flag fast_Q_flag`, `NRr NRz`, exact receiver ranges in metres, and exact artificial depths in metres. Counts/order and native SHD-coordinate quantization are checked. No field interpolation or gain fit is performed.

The module is copied into an isolated source build by `prepare_bellhop_rotated_field_overlay`. It is not installed into the official toolbox. Point spreading, coordinate handling and range search must pass their fixed-sample gates before any full field.

Native SHD coordinate records remain float32; MATLAB consumers must retain and verify the exact sidecar coordinates used by the influence routine. At zero rotated range the workflow must use another coordinate chart, retaining only the true source singularity.

See `reports/bellhop_rotated_fig12b_stage2_2026-11-14.md` for the preceding accepted scope and unresolved full-field gates.