# PE--Bellhop controlled comparison P0 preparation

状态：**PASS**

本阶段只建立 source-X/fingerprint-safe validation tree并核对依赖；没有调用 PE 或 Bellhop。

## Frozen configuration

- Bellhop：official AcousticsToolbox 2020 executable；run type `C`；source geometry `X`。
- f/c：4000 Hz / 1500 m/s；PE `W=192.1875 m`, `nx=984`, `step=0.05 m`。
- receiver map：`z_BH'=-x_PE`; explicit ranges `[102,103] m`; sorted write/read-back permutation saved in MAT/JSON.
- source-pattern fingerprint：`X|cos(theta)*exp(-0.5*(k*sigma*sin(theta))^2)|N=2401|angle=[-30,30]deg|f=4000Hz|sigma=0.3m`。
- canonical PM coefficient source is deferred until Stage 6; expected SHA-256 `1f7eda465e4ae85b8ac038310edf053f2d062563a3b943bfd83bd59d07015f67`.

## File hashes

| file | SHA-256 |
|---|---|
| official_exe | `7e7809a64c3bf734aff6d28d0d4d52b1b4bd203d81676e3241ffd3189941b505` |
| flat_validation_exe | `1cfaba7ba7bc19738ebffa5361c2d2866cd1abf05c1c3501df6657fad19d5341` |
| parametric_validation_exe | `81de6a512037d63f06ede0537d5a591fde5e41b4ab9a95ff57f53f66917a5360` |
| flat_overlay_bellhop | `ad8ce36e5eb199264b96eb8e1e151cb299fdba3a1c1aa9e2158a60dd86a613bc` |
| flat_overlay_step | `98a1a250b33c2e11309f61e920f3d335d9fadf152398c8847a909ec81f32d43d` |
| parametric_overlay_bellhop | `986010216f03d1ef2156730f076fcec27ad8a39cc1a5eb44096726ee33e0e251` |
| parametric_overlay_step | `1f54c01f35ac9682c80d374abc7cfc83ec9ebf72bd4928feded5e19323c757de` |

## Checks

- official_exe: PASS
- flat_validation_exe: PASS
- parametric_validation_exe: PASS
- flat_overlay_bellhop: PASS
- flat_overlay_step: PASS
- parametric_overlay_bellhop: PASS
- parametric_overlay_step: PASS
- receiver_map_finite: PASS
- receiver_map_within_domain: PASS
- profile_count_valid: PASS
- source_x: PASS
- run_c: PASS
- official_2020_path: PASS
- all: PASS

P0 result is **READY**. Only after this result is reviewed may Stage 0 be unlocked explicitly.
