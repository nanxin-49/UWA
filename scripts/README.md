# 脚本与验证运行索引

脚本按职责分组，入口应能从任意 MATLAB 当前目录运行。推荐先执行根目录 `setup_vertical_project`；它只初始化路径，不改变当前工作目录。仍使用 `scripts/bootstrap_project.m` 的旧脚本会先委托给统一 setup，再按既有规则切换到 `results/<category>`。新入口应优先使用绝对项目根和明确输出目录，不依赖当前目录。

## 目录

- `validation/`: smoke tests, invariants, communication-chain checks, and disabled-path regressions.
- `comparisons/`: surface-model, bubble-model, and coherent/incoherent comparison studies.
- `experiments/`: Monte Carlo runs, parameter sweeps, and calibration scripts.
- `reporting/`: figure generation, summary tables, visualizations, and report builders.
- `validation/support/`: Bellhop/Weyl/Li2009、properness、发布元数据和 artifact 检查等验证专用辅助函数；不属于公共 API。

## PE--Bellhop 展开坐标横向审计

本节中的早期PE–Bellhop、epsilon诊断和历史绘图入口已按用户要求移入
`../cash/pe_bellhop_validation_scripts_2026-09-01/`。完整文件清单及归档原因见
`../reports/pe_bellhop_validation_archive_manifest_2026-09-01.md`；归档内容不应从活动目录调用。

- `../reports/pe_bellhop_complete_validation_overview.md`: 当前 PE--Bellhop
  平面海面验证链的总报告，汇总早期原坐标、自由场、Gaussian 窗口、正式
  展开坐标、横向三方审计、参数扫描和坐标极限结果；旧的“完整汇总”仅作
  历史证据。
- `validation/validate_pe_bellhop_unfolded_flat_gaussian_vertical.m`: 正式的展开坐标平面海面 PE--Bellhop 对比，读取/生成 4--8 kHz Gaussian 直达、镜像反射、波束收敛和 PDP 结果。
- `validation/validate_pe_as_bellhop_transverse_vertical.m`: 只读取正式 MAT，在 4/6/8 kHz 的六个横向采样点上计算独立一步角谱，分离 PE--AS 与 Bellhop--AS 的误差来源，不重新运行完整 PE 或 Bellhop。
- `validation/validate_bellhop_transverse_parameter_scan_vertical.m`: 只运行 Bellhop，在 8 kHz 扫描步长、角扇区和 `.sbp` 采样密度，检查这些数值设置能否解释横向相位差。
- `validation/validate_bellhop_rotation_limit_vertical.m`: Bellhop-only 坐标极限审计，比较当前展开坐标、原坐标 `-89°` 和 `-89.5°` 的同源 Gaussian 直达/海面反射路径。
- `../cash/pe_bellhop_validation_scripts_2026-09-01/scripts/validation/validate_bellhop_rough_surface_near_vertical_limit.m`: **已归档**的历史500 m入口。
- `../cash/pe_bellhop_validation_scripts_2026-09-01/scripts/reporting/plot_bellhop_rough_surface_near_vertical_rays.m`: **已归档**的历史代表性绘图入口。
- `../cash/pe_bellhop_validation_scripts_2026-09-01/scripts/reporting/plot_bellhop_rough_surface_full_paths.m`: **已归档**的历史500 m全路径绘图入口；正常检索不读取`cash/`。
- `../cash/pe_bellhop_validation_scripts_2026-09-01/scripts/validation/validate_bellhop_100m_boundary_source_vertical.m`: **已归档**的100 m水侧边界源epsilon检查。
- `../cash/pe_bellhop_validation_scripts_2026-09-01/scripts/validation/validate_bellhop_100m_arrival_field_audit_vertical.m`: **已归档**的粗糙面arrival/C场诊断。
- `reporting/generate_bellhop_100m_source1m_visuals_vertical.m`: **当前代表绘图入口。** 用户改为实际离底1 m（Tx99 m、水深100 m、Rx3 m），停止epsilon诊断。独立A/E/R显示全部接收贡献及1001条发射声线；接收范围12 m，展示扇区范围120 m，绘图还原物理深度。`render_only=true`仅重绘本轮MAT，不运行Bellhop或读取退役数据。
- `reporting/generate_bellhop_pe_equivalent_ray_paths_vertical.m`: 从当前离底1 m结果中选择A模式唯一的0/0直达和1/0一次海面arrival，以其发射角单独运行两条R声线，并在最接近Rx的位置截断；不显示任何海底或高阶多次反射路径。结果位于当前案例的 `pe_equivalent_paths/` 子目录。
- `../reports/bellhop_100m_rebaseline_plan.md`: 当前 100 m 水深、物理源离底1 m的计划；旧零离底内容仅为历史记录。

正式结果位于 `results/validation/pe_bellhop_unfolded_flat_gaussian/`；横向审计结果位于 `results/validation/pe_as_bellhop_transverse/`，解释报告为 `reports/pe_as_bellhop_transverse_audit_report.md`。
Bellhop-only 扫描结果位于 `results/validation/pe_as_bellhop_transverse/bellhop_scan/`，报告为 `reports/bellhop_transverse_parameter_scan_report.md`。
坐标极限审计结果位于 `results/validation/bellhop_rotation_limit/`，报告为 `reports/bellhop_rotation_limit_report.md`。
旧 500 m 粗糙海面结果已按用户要求移入不参与正常读取的 `cash/` 隔离区，原报告路径仅保留迁移提示。新边界源结果位于 `results/validation/bellhop_100m_rough_surface/boundary_source/explicit_bty/`；父目录中的初次隐式平底诊断不能用作收敛证据。当前状态见 `reports/bellhop_100m_rough_surface_report.md`；归档清单见新计划，不需要打开归档目录。

边界源检查入口已归档（当前活动目录不再提供复现命令）；结果状态、参数和失败原因保留在
`../reports/bellhop_100m_boundary_source_report.md`及归档清单中。

## 运行方式

当前离底1 m代表图的已完成结果及限制见
`../reports/bellhop_100m_source1m_visuals_report.md`：193条ARR、199条E贡献、1001条
扇区；接收范围12 m、展示范围120 m，33条展示轨迹到达内部存储上限。
同一报告也记录只保留PE直达和一次海面反射拓扑的两条Bellhop轨迹图。

```matlab
setup_vertical_project
run('scripts/validation/validate_surface_wavefield_visualization_vertical.m')
run('scripts/comparisons/compare_specular_incoherent_surface_reflection_vertical.m')
run('scripts/experiments/sweep_monte_carlo_surface_channel_vertical.m')
run('scripts/reporting/plot_c35_core_heatmaps_vertical.m')
```

公共可复用实现已移到 `src/`。下文脚本名仍按 `validation/`、`reporting/` 等相对于 `scripts/` 的短路径描述。

## 2k 相位与海谱诊断

```matlab
run('scripts/validation/validate_2k_phase_approx.m')
run('scripts/validation/validate_2k_ocean_spectra.m')
```

- PNG/CSV：`results/validation/ssa_2k_phase/`
- 报告：`reports/validation_2k_phase_report.md`、`reports/validation_2k_ocean_spectra_report.md`

这两项是相位系数近似和最低阶 SSA/海谱趋势诊断，不是完整 PE、KStat、真实海洋散射或高阶 SSA 的验收。

## Raw-PM and Joint-Frequency K-Stat Prerequisite Validation

- `validation/validate_raw_pm_grid_coverage_vertical.m`: scans PM aperture, FFT-grid resolution, wind-speed coverage, and same-dx large-PM-grid to small-PE-window energy mapping.
- `validation/validate_kstat_joint_frequency_boundary_vertical.m`: boundary-only comparison of one shared explicit Gaussian PM surface, independent-frequency kstat, and covariance+pseudocovariance joint-frequency kstat at U=5 m/s and F=32.
- `validation/audit_pe_caching_vertical.m`: counts cacheable PE marches/FFTs, measures the unchanged public-path baseline, and estimates batched cache memory.

These scripts do not change the public defaults or communication chain. Detailed interpretation is in `reports/raw_pm_joint_kstat_prerequisite_validation_report.md`.

## Cached Joint-Kstat Receiver Validation

- `validation/validate_cached_joint_kstat_pe_receiver_vertical.m`: runs the U=5, F=32, 128/64 train/test reflected-only kdomain/independent/joint comparison through the independent cached PE executor.
- `validation/validate_cached_pe_public_consistency_vertical.m`: checks the cached executor against the unchanged public kdomain path with the exact same double-precision boundary input.
- `validation/analyze_receiver_properness_null_vertical.m`: calibrates finite-sample receiver `||P||_F/||C||_F` under a proper complex-Gaussian null model.

The receiver report is `reports/cached_joint_kstat_pe_receiver_validation_report.md`. The main result MAT and PNG diagnostics are under `results/validation/cached_joint_kstat_pe_receiver/`.

## U=5 F=64 Conditional Channel Generator

- `validation/validate_u5_conditional_channel_generator_vertical.m`: set `U5_CONDITIONAL_MODE=smoke` for 32/16 samples or `full` for 128/64 plus 10,000 generated channels.
- `validation/refresh_u5_f64_properness_vertical.m`: refreshes the F=64 three-mode null distributions and prescribed central-interval decisions.
- `validation/generate_u5_f64_sample_bundle_vertical.m`: saves generated H, total H, physical CIR, delay axis, labels, and timing without PE calls.
- `reporting/plot_u5_f64_conditional_validation_vertical.m`: produces QQ, correlation-matrix, eigenvalue, and LFM diagnostics.

Reusable interfaces are initialized by `setup_vertical_project` and implemented under `src/`（`properness_null_test_vertical` is validation support）: `estimate_conditional_channel_stats_vertical`, `sample_conditional_channel_vertical`, `build_channel_cir_vertical`, and `properness_null_test_vertical`. `build_physical_cir_vertical` remains only as the legacy common-time-shift wrapper. Full results are under `results/validation/u5_conditional_channel_f64/`; interpretation is in `reports/u5_conditional_channel_generator_f64_report.md`.

## Optimized Joint Builder and U=8 Node

- `validation/validate_streaming_joint_builder_u5_vertical.m`: checks the streaming/series F=64 builder against the original U=5 factors and cached-PE receiver statistics.
- `validation/audit_raw_pm_u8_aperture_vertical.m`: audits 100/256², 150/384², and 200/512² PM grids with 128 central-crop mapping realizations.
- `validation/validate_u8_conditional_channel_generator_vertical.m`: uses `U8_CONDITIONAL_MODE=smoke` for 16/8 samples or `full` for 128/128, U=8 properness, rank comparison, and 10,000 H+CIR samples.
- `validation/build_u5_u8_conditional_library_vertical.m`: creates and smoke-tests the exact-node U=5/U=8 library after both full models exist.

The reusable library interfaces are `build_conditional_channel_library_vertical` and `sample_conditional_channel_library_vertical`. They deliberately reject unsupported wind speeds and do not interpolate. The U=8 full 128/128 validation and two-node smoke test pass; see `reports/u8_joint_optimization_two_node_library_report.md` for measured acceptance and remaining high-wind memory limits.

## Two-Node Communication Validation

- `validation/validate_two_node_communication_vertical.m`: set `TWO_NODE_COMM_MODE=smoke` or `full`; compares fresh kdomain PE, joint cached PE, full-rank statistics, and 99.9% statistics at U=5/U=8.
- The full run uses 32 fresh channels per source, 8000 QPSK symbols per channel, Eb/N0 0:2:20 dB, channel-cluster bootstrap intervals, and a separate 128-pair full/low-rank test.
- `comm_main_vertical_psk` accepts optional external H(f) or h(t) through `COMM_EXTERNAL_CHANNEL_FILE`. Empty/unset preserves the original PE scenarios.

Reusable interfaces are `build_communication_taps_vertical`, `evaluate_mpsk_channel_ensemble_vertical`, and `sample_conditional_channel_rank_pair_vertical`. See `reports/two_node_statistical_channel_communication_validation_report.md`.

## PE Carrier-Phase Release Candidate Workflow

The formal workflow is destructive only in the narrow archival sense: the
preparation step moves registered phase-sensitive results to a timestamped
folder under `results/archive/pe_phase_reference_pre_rc/`. It never deletes
or overwrites them. Run the stages in this order:

```matlab
run('scripts/validation/prepare_pe_phase_reference_release_candidate_vertical.m')
run('scripts/validation/validate_pe_channel_phase_reference_vertical.m')
run('scripts/validation/validate_pe_phase_convention_uniform_vertical.m')

setenv('ADJOINT_PE_VALIDATION_MODE','full')
run('scripts/validation/validate_adjoint_pe_receiver_projection_vertical.m')
setenv('U5_CONDITIONAL_MODE','full')
run('scripts/validation/validate_u5_conditional_channel_generator_vertical.m')
setenv('U8_CONDITIONAL_MODE','full')
run('scripts/validation/validate_u8_conditional_channel_generator_vertical.m')
run('scripts/validation/build_u5_u8_conditional_library_vertical.m')

setenv('TWO_NODE_COMM_MODE','full')
run('scripts/validation/validate_two_node_communication_vertical.m')
run('scripts/validation/validate_public_channel_modes_vertical.m')
run('scripts/validation/validate_cached_pe_public_consistency_vertical.m')

setenv('PE_ATLAS_MODE','full')
run('scripts/reporting/generate_pe_propagation_atlas_vertical.m')
run('scripts/validation/finalize_pe_phase_reference_release_candidate_vertical.m')
```

`validation/audit_phase_reference_artifacts_vertical.m` is the read-only inventory
entrypoint and may be run separately. The active run metadata is stored in
`results/validation/pe_phase_release_candidate/current_run.mat`. A stopped
run may resume only when its code fingerprint is unchanged. If a
fingerprint-covered source changes, rerun the preparation step: it archives
the partial run and creates a new `run_id`.

The formal adjoint configuration is PE `128^2` / PM `256^2`; the U=8
conditional node intentionally uses its separately audited PE `128^2` / PM
`384^2` aperture. F=9 uses 4096 realizations, F=64 uses 512, and these runs
are expensive. The atlas must be last because it rejects validation inputs
whose `run_id` does not match the active release run. The finalizer reports
only `PASS`, `FAIL`, or `INCOMPLETE` and writes
`reports/pe_phase_reference_release_candidate_report.md`.

The completed reference run is `phase_rc_20260722_174945` and is `PASS`.

## Exact Adjoint PE Receiver Projection

Before the adjoint suite, the receiver carrier-reference integration can be
checked independently:

```matlab
run('scripts/validation/validate_pe_channel_phase_reference_vertical.m')
```

This reduced-cost audit uses F=65 for the unaliased 4 ms delay test, exercises
the four public surface branches, checks `legacy_reduced`, compares cached and
adjoint receiver outputs in both reduced/direct-DSP form, validates dense/FFT
`C/P`, and tests schema-1 conditional-model migration. Outputs are under
`results/validation/pe_channel_phase_reference/`.

- `validation/validate_adjoint_pe_receiver_projection_vertical.m`: validates the exact discrete conjugate transpose of the cached uniform surface-to-receiver PE, receiver projection, PM-to-PE embedding, dense/FFT receiver statistics, realization statistics, performance, and public regressions.
- Set `ADJOINT_PE_VALIDATION_MODE=smoke` for the reduced run or `full` for the accepted F=9/4096 and F=64/512 validation.
- Reusable validation interfaces are `apply_forward_surface_to_receiver_vertical`, `apply_adjoint_receiver_to_surface_vertical`, `build_adjoint_receiver_projection_vertical`, `run_adjoint_receiver_projection_vertical`, and `contract_kstat_receiver_stats_vertical`.
- `reporting/plot_adjoint_pe_receiver_projection_validation_vertical.m` regenerates the adjoint error, covariance/pseudo-covariance, PDP, LFM/matched-filter, augmented-eigenvalue, and timing/memory figures from a saved validation result.
- Outputs are written to `results/validation/adjoint_pe_receiver_projection/`; measured results and the integration decision are in `reports/adjoint_pe_receiver_projection_feasibility_report.md`.

This v1 path is limited to uniform sound speed, CPU double, fixed grids and frequency axis, one nearest-grid receiver, no bubbles, and no Doppler. Its PE operator remains a validation path; receiver outputs now use the same central phase-reference layer as the public API. The default surface model is unchanged.

## PE Propagation Visual Atlas

- `reporting/generate_pe_propagation_atlas_vertical.m` is the report-only entrypoint. It assembles the Tx-to-surface and surface-to-Rx center slices, surface and receiver planes, physical boundary branches, exact-adjoint sensitivity, accepted receiver `C/P`, wind-node statistics, and a signal-free carrier-reconstruction animation.
- `reporting/plot_pe_propagation_atlas_vertical.m` renders the fixed 14-item atlas. Spatial panels use a shared amplitude reference and a `[-50,0] dB` scale; phase is masked below `-40 dB`.
- Run `set PE_ATLAS_MODE=smoke` before MATLAB for the PE 64² / PM 128² / F=9 check. With the variable unset, the formal run reads the accepted PE 128² / PM 256² / F=64 adjoint and conditional-model results.
- Outputs are written to `results/visualization/pe_propagation_atlas/`, including a reusable MAT file, source/normalization manifest, numerical-closure summary, 13 PNG figures, and one MP4.

The atlas distinguishes physical boundary models from computational acceleration paths. `q` is an exact receiver-sensitivity kernel, not a physical reverse-propagated pressure field. LFM is applied only after obtaining `H(f)` and is shown without noise, modulation, synchronization, or equalization.

## PE/Bellhop Flat-Surface Cross-Validation

### Unfolded Gaussian-source comparison

- `validation/validate_pe_bellhop_unfolded_flat_gaussian_vertical.m` is the
  current validation-only entrypoint for the accepted unfolded-coordinate
  construction.  It maps the 100 m to 3 m vertical path to Bellhop ranges
  97 m (direct) and 103 m (image/reflected), then applies the flat
  pressure-release factor `-1` at the image receiver.
- `validation/support/write_bellhop_unfolded_gaussian_env_vertical.m` writes
  the per-frequency `.env` and `.sbp`; the `.sbp` pattern is derived from the
  production Gaussian angular spectrum and is not a pointwise receiver fit.
- `validation/support/run_bellhop_unfolded_gaussian_vertical.m` runs Bellhop
  and reads the standard `.shd`/`.arr` outputs.  The primary acceptance uses
  normalized reflected/direct responses; absolute source pressure is an
  independent Weyl-reference diagnostic.
- Outputs are written to
  `results/validation/pe_bellhop_unfolded_flat_gaussian/` and the report is
  `reports/pe_bellhop_unfolded_flat_gaussian_report.md`.

### Bellhop 2020 parametric internal-wall validation

- The authoritative implementation and status report is
  `../reports/bellhop_internal_wall_implementation_report.md`. It consolidates
  the feasibility, flat, tilted, sinusoidal, beam/frame, old PM, redesign and
  parameterization-fix stages. Superseded reports and obsolete one-off scripts
  are in the recoverable root `cash/` quarantine, not active project context.
- The internal-wall visualization scripts and generated figures/fields described
  in the next bullets are now archived under
  `../cash/bellhop_internal_wall_visualization_archive_20260909/`; they are
  retained for reproducibility only and are not active entrypoints.
- `reporting/generate_bellhop_internal_wall_visuals_vertical.m` is the active
  read-only renderer for the two requested products. It consumes the stored
  dense tilted (`r=100+0.05z`) and sinusoidal (`r=100-2sin(0.04z)`) internal-
  wall fields, draws eight incident/native-backward/inverse-mapped rays with
  arrows plus reordered-axis wall zooms, and writes dense 2-D coherent
  reflected-only TL backgrounds after `T^{-1}(r',z')=(200-r',-z')`. The SHD
  input contains only the transformed post-wall branch, not a total
  direct-plus-reflected field. It also renders a separate dense
  incident-only TL field from the official matched-halfspace `C*X` run.
  Outputs are under
  `results/visualization/bellhop_internal_wall_visuals/`; the dense source
  fields and their validation cases are under
  `results/validation/bellhop_internal_wall_visuals/`.
- `validation/validate_bellhop_internal_flat_wall_poc.m`,
  `validation/validate_bellhop_internal_tilted_wall_poc.m`, and
  `validation/validate_bellhop_internal_sinusoidal_wall_poc.m` are the retained
  flat/straight/curved regression entrypoints. Their isolated source overlays
  and runners remain under `validation/support/`.
- `../cash/bellhop_internal_wall_visualization_archive_20260909/scripts/reporting/generate_bellhop_internal_wall_reflection_visuals.py` is a
  read-only postprocessor for the accepted tilted/sinusoidal `.iwdiag`, `.iw3`,
  and convergence CSV artifacts. It does not launch Bellhop or PE; it renders
  the physical reflection and the isolated proper-rotation chart in separate
  panels under `results/visualization/bellhop_internal_wall_reflection/`.
  The primary `04_internal_wall_end_to_end_overlay.png` overlays the incident
  segment, dashed native backward branch, and solid proper-rotated receiver
  branch in one chart for each wall case.
  Interpretation is documented in
  `../reports/bellhop_internal_wall_reflection_visualization_report.md`.
- `../cash/bellhop_internal_wall_visualization_archive_20260909/scripts/reporting/generate_bellhop_internal_wall_visual_enhancement.py` is the
  read-only postprocessor for the stronger display-only case in
  `results/validation/bellhop_internal_wall_visual_enhancement/`:
  tilted `a=0.05` and sinusoidal `A=2 m, K=0.04 1/m, N=161`, both at
  `step=0.05 m` and `5001` beams. It consumes the stored `.iwdiag`/`.iw3`
  sidecars, writes independent `01_tilted_wall_end_to_end.png` and
  `02_sinusoidal_wall_end_to_end.png` figures (each with complete context and
  a multi-ray wall-neighborhood zoom), the compatibility two-panel overlay,
  the `Delta-r` versus `z` geometry plot, and `ray_direction_audit.csv` to
  `results/visualization/bellhop_internal_wall_visual_enhancement/`, and
  never launches Bellhop or MATLAB. The audit checks `norm(rot+ref)`, plotted
  displacement-direction residuals, and equal native/rotated branch lengths
  for every selected fan/target ray. The experiment is display-only and does
  not replace the weak-wall regression or its conclusion; interpretation,
  audit status, and hard-check status are recorded in
  `../reports/bellhop_internal_wall_visual_enhancement_experiment_report.md`.
  The same postprocessor also writes
  `03_tilted_wall_physical_reconstruction.png`,
  `04_sinusoidal_wall_physical_reconstruction.png`,
  `inverse_rotation_ray_audit.csv`, and
  `inverse_rotation_visualization_manifest.csv`. These products inverse-map
  each stored mapped forward branch with
  `T^-1(r',z')=(2R0-r',-z')` and overlay it on the original physical
  `Reflect2D` backward branch. The mapped forward branch is not drawn as a
  third subplot in these physical reconstruction figures. The physical
  wall-neighborhood panels
  re-order the display axes to `(z,r)` (horizontal `z`, vertical range `r`)
  to match the compact audit view; this is explicitly annotated and screen
  slopes are not angle measurements.
- `../cash/bellhop_internal_wall_visualization_archive_20260909/scripts/reporting/generate_bellhop_internal_wall_tl_diagnostic.py` is a read-only
  postprocessor for the same stored `.iwdiag` and validation MAT data. It
  writes `05_tl_diagnostic.png`, `tl_diagnostic.csv`, and
  `tl_visualization_manifest.csv`. The upper panels show receiver-relative
  TL after mapping `r'=2R0-r`; the lower panels show relative reflected-beam
  level and the pressure-release reflection amplitude jump. Absolute TL and
  native/rotated amplitude remain diagnostic only. A zero stored native
  sinusoidal pressure is reported as unavailable rather than converted to a
  spurious TL value.
- `../cash/bellhop_internal_wall_visualization_archive_20260909/scripts/validation/run_bellhop_internal_wall_tl_grid.py` reruns only the accepted
  validation-only internal-wall binaries with a dense coherent receiver grid
  (`603` mapped ranges from `100.05--130 m`, `321` depths from `-32--32 m`,
  5001 beams, `step=0.05 m`). It reuses the existing `.sbp` and wall profile
  inputs and writes compressed fields under
  `results/validation/bellhop_internal_wall_tl_grid/`; it does not change any
  Bellhop source or the established wall implementation.
- `../cash/bellhop_internal_wall_visualization_archive_20260909/scripts/reporting/generate_bellhop_internal_wall_tl_2d.py` reads those dense SHD
  fields and writes `06_internal_wall_tl_2d.png` plus a grid summary. It
  applies `r=2R0-r'`, `z=-z'` before plotting, so the main product is a
  physical `(r,z)` coherent reflected TL background rather than a sparse
  receiver diagnostic.
- Curved overlays use ordered noncoincident parametric segments, a normalized
  hit tangent, a TOP normal reconstructed from that tangent, signed
  turning-angle/arclength curvature, and finite profile support. Internal-wall
  geometry does not use `dz/dr`, `Dss`, `Delta z/Delta r`, `1/Delta r`, range
  monotonicity or ATI-style infinite endpoint extension.
- `validation/validate_bellhop_internal_vertical_tangent_poc.m` is the focused
  circular-wall regression at exact `dr/ds=0` and non-grazing `u dot n=1`. It
  scans 65/129/257 samples and checks intersection, orthonormal frame, signed
  curvature, `RN`, p/q, one pressure-release phase jump, finite values and
  positive transformed range.
- `validation/support/run_bellhop_internal_pm_wall_poc_vertical.m` and
  `validation/support/bellhop_internal_pm_wall_poc/` are retained as the
  generic parametric-wall runner/build. The fixed-realization convergence
  completed entrypoint is archived at
  `../cash/bellhop_internal_wall_superseded_20260902/scripts/validation/validate_bellhop_internal_pm_fixed_realization_convergence.m`;
  it wrote one master Fourier realization and interpolated `N=513/1025/2049/4097`
  profiles from that same realization. Outputs remain under
  `results/validation/bellhop_internal_pm_fixed_realization/`.
- The invalid strict-90-degree PM comparator and its old profile generator are
  archived. The fixed-PM result and the finite-angle native-ATI/local-wall
  covariance audit are Bellhop-only PASS results; the controlled PE comparison
  and reduced ensemble smoke are complete, while large-scale Monte Carlo remains
  out of scope.
- The completed Bellhop-only fixed-realization local covariance audit is archived at
  `../cash/bellhop_internal_wall_superseded_20260902/scripts/validation/validate_bellhop_pm_local_reflection_covariance.m`;
  its one-off runner is archived at
  `../cash/bellhop_internal_wall_superseded_20260902/scripts/validation/run_bellhop_pm_local_covariance_vertical.m`.
  It reused the same seed-260001 master samples in native C-ATI and a
  source-centered proper rotation at `phi=89` and `89.5` degrees, with
  `N=2049/4097` and three exact paired rays. Results remain under
  `results/validation/bellhop_pm_local_reflection_covariance/`; the detailed
  stage report is archived under
  `../cash/bellhop_internal_wall_superseded_20260902/reports/` and its valid
  conclusion is merged into the authoritative Bellhop implementation report.
- `validation/validate_bellhop_shd_receiver_range_pairing_vertical.m` is the
  permanent receiver-column regression. It requires exact 103 m rotated and
  97 m native matches, verifies that native total/direct use the same column,
  and keeps the 102 m rotated column as a negative control for the former
  approximately `-2.095 rad` indexing error.
- Current Bellhop internal-wall status is **PASS_WITH_LIMITS**: flat, tilted,
  sinusoidal, vertical tangent, fixed-realization PM density convergence, and
  the finite-angle local native↔internal Reflect2D covariance audit pass.
  Complete receiver-field equivalence and large-scale Monte Carlo remain out of
  scope; the controlled PE rough-PM comparison is indexed below.
- The PE rough-PM comparison execution chain is implemented as independent
  Stage 0A--0E, Stage 1A--1C, Stage 2 and Stage 3 validation entrypoints. The
  original scope and comparability matrix remain in
  `../reports/pe_bellhop_pm_comparison_design_audit_report.md`; the complete
  executable index and result interpretation are in
  `../reports/pe_bellhop_pm_complete_execution_summary.md`.
- Stage 0A canonical profile provenance is now implemented by
  `validation/validate_pe_bellhop_pm_canonical_mapper.m` with support helpers
  under `validation/support/`. It is a mapper-only audit (no PE/Bellhop run),
  uses the existing seed-260001 Fourier coefficients, and writes its report to
  `../reports/pe_bellhop_pm_canonical_mapper_report.md` and results to
  `results/validation/pe_bellhop_pm_canonical_mapper/`. Stage 0B and all
  subsequent planned stages are complete; see the complete execution summary.
- Stage 0B one-transverse-dimensional PE bridge is implemented by
  `validation/validate_pe_1d_validation_bridge.m` and
  `validation/support/run_pe_1d_surface_reflection_validation.m`. It passes
  independent 1-D exact-AS and production-PE `k_y=0` checks for flat and weak
  phase-screen cases. Results are under
  `results/validation/pe_1d_validation_bridge/`; report:
  `../reports/pe_1d_validation_bridge_report.md`. Stage 0C and all subsequent
  planned stages are complete; see the complete execution summary.
- Stage 0C flat source/normalization audit is implemented by
  `validation/validate_pe_bellhop_pm_stage0_flat_source.m`. It compares the
  unchanged 1-D PE bridge with the existing Bellhop 2020 Gaussian `.sbp`
  internal-flat case using explicit 97/103 m SHD selectors and 5001/10001
  beams. It passes normalized offset-magnitude and axis-Q checks; the known
  ~0.26 dB backward-range influence offset remains diagnostic. Results are in
  `results/validation/pe_bellhop_pm_stage0_flat_source/`; report:
  `../reports/pe_bellhop_pm_stage0_flat_source_report.md`. Stage 0D and all
  subsequent planned stages are complete; see the complete execution summary.
- The pre-reflection incident-plane audit is implemented by
  `validation/validate_pe_bellhop_incident_field_vertical.m`. It compares the
  complete 1-transverse-dimensional PE incident field at 100 m with Bellhop
  2020 coherent `C *X` line-source output on the same transverse grid, using
  the saved spatial-phase sign and no surface/wall reflection. The default
  4 kHz run checks 5001/10001 beams and a half-dx receiver diagnostic. After
  the X integration all 13 gates pass; PE--Bellhop TL P95 over M95 is
  `0.0019036 dB`, versus `0.157943 dB` in the archived pre-X run. Outputs are
  under `results/validation/pe_bellhop_incident_field/` and the report is
  `../reports/pe_bellhop_incident_field_comparison_report.md`.
- Stage 0D constant-height sign audit is implemented by
  `validation/validate_pe_bellhop_pm_constant_height_sign.m`. It checks
  `eta0=+0.05, 0, -0.05 m` with image span `L=103-2*eta0`, and passes PE,
  Bellhop and cross-model phase-sign, path-time, pressure-release and rotation
  checks. Results are under
  `results/validation/pe_bellhop_pm_constant_height_sign/`; report:
  `../reports/pe_bellhop_pm_constant_height_sign_audit.md`. Stage 0E and all
  subsequent planned stages are complete; see the complete execution summary.
- Stage 0E numerical-budget freeze is implemented by
  `validation/validate_pe_bellhop_pm_numerical_budget.m`. It reuses one
  seed-260001 band-limited realization across PE window/grid/step and Bellhop
  profile/beam/step scans; all frozen numerical gates pass. Results are under
  `results/validation/pe_bellhop_pm_numerical_budget/`; report:
  `../reports/pe_bellhop_pm_numerical_error_budget.md`.
- Stage 1A Tier-1 fixed-PM reflected-only ratio audit is implemented by
  `validation/validate_pe_bellhop_pm_stage1_tier1.m`. Bellhop structural
  geometry/phase/beam-state checks pass. The active code now defaults to
  source geometry X and uses X-specific case roots; the historical report
  remains pre-migration provenance, while the executed X result is recorded
  in the source-geometry audit. The PE/Bellhop ratio difference is diagnostic.
  Results are under
  `results/validation/pe_bellhop_pm_stage1_tier1/`; report:
  `../reports/pe_bellhop_pm_stage1_tier1_report.md`.
- The point-source `R` versus line-source `X` covariance audit is implemented
  by `validation/validate_bellhop_source_geometry_rx_audit.m`. It changes only
  Bellhop `RunType(4)`, checks flat and weak-sinusoidal walls at 5001/10001
  beams, and conditionally runs only the 4 kHz seed-260001 Tier-1 case. At
  10001 beams, `X` removes the native-97 m/internal-103 m `-0.26065 dB`
  range-normalization artifact, while the Tier-1 PE--Bellhop delta remains
  `0.31037 dB / -2.24820 rad`. Results are under
  `results/validation/bellhop_source_geometry_rx_audit/`; report:
  `../reports/bellhop_source_geometry_rx_audit_report.md`.
- Integration and archive details are recorded in
  `../reports/pe_bellhop_line_source_integration_report.md` and
  `../reports/pe_bellhop_source_geometry_integration_archive_manifest.md`.
- Stage 1B production PE dimensionality sensitivity is implemented by
  `validation/validate_pe_bellhop_pm_stage1_dimensionality.m`; it copies the
  same fixed eta(x) across y and compares ny=256/512 against the 1T bridge.
  Results are under `results/validation/pe_bellhop_pm_stage1_dimensionality/`;
  report: `../reports/pe_bellhop_pm_stage1_dimensionality_report.md`.
- Stage 1C interpretation/freeze is implemented by
  `validation/validate_pe_bellhop_pm_stage1_interpretation.m` and classifies
  the 4 kHz result as PASS_WITH_MODEL_DISCREPANCY. The cross-model residual is
  reported separately from the smaller dimensionality sensitivity. Results are
  under `results/validation/pe_bellhop_pm_stage1_interpretation/`; report:
  `../reports/pe_bellhop_fixed_pm_4khz_comparison_report.md`.
- Stage 2 fixed-PM frequency extension is implemented by
  `validation/validate_pe_bellhop_pm_frequency_extension.m`. It runs the
  same realization at 4/6/8 kHz with PE flat/rough and Bellhop flat/internal-
  wall rough pairs; numerical geometry/phase/state checks pass and residuals
  remain model diagnostics. Results are under
  `results/validation/pe_bellhop_pm_frequency_extension/`; report:
  `../reports/pe_bellhop_fixed_pm_frequency_extension_report.md`.
- Stage 3 paired fixed-band PM ensemble smoke is implemented by
  `validation/validate_pe_bellhop_pm_ensemble.m` with seeds 260001:260008.
  It preserves the canonical per-mode spectral amplitudes and records paired
  PE/Bellhop TL, phase, power and native wall diagnostics, including
  machine-readable percentile summaries. Results are under
  `results/validation/pe_bellhop_pm_ensemble/`; report:
  `../reports/pe_bellhop_pm_ensemble_comparison_report.md`. The status is
  PASS_WITH_LIMITS because this is a phase-ensemble smoke, not an independent
  PM amplitude Monte Carlo.
- Stage 3B independent PM coefficient-amplitude ensemble is exposed by
  `validation/validate_pe_bellhop_pm_amplitude_ensemble.m`. It reuses the
  same fixed spectral-density band and canonical seed anchor, then draws
  zero-mean Gaussian cosine/sine coefficients with variance `S(k)*Delta-k`
  for the remaining paired seeds. The reduced 8-seed run uses 5001 beams for
  the flat baseline and every rough case, requires 5001/5001 wall hits, and
  passes all structural checks. Stage-2 reuse requires an exact embedded
  configuration match plus current artifacts; other cache reuse requires a
  full request fingerprint and verified `.shd`/`.iwdiag` output hashes. Its
  minimum accepted seed count is eight. Results belong under
  `results/validation/pe_bellhop_pm_amplitude_ensemble/` and its report is
  `../reports/pe_bellhop_pm_amplitude_ensemble_report.md`.

The complete stage index, frozen parameters, cross-model interpretation and
remaining statistical limitation are summarized in
`../reports/pe_bellhop_pm_complete_execution_summary.md`.

### Reflection-free four-level audit

- `validation/validate_pe_as_freefield_vertical.m` compares production
  multi-step marching with an independently coded one-step exact discrete
  angular-spectrum propagator. The sponge is exactly off for this hard gate.
- `validation/validate_bellhop_freefield_normalization_vertical.m` reproduces
  the matched-halfspace construction of Bellhop's official free-space point
  source example. It measures `|p|R`, spatial phase sign, the constant source
  phase, beam/step convergence, and zero-bounce arrival delay before applying
  the single global `1/(4*pi)` Green-function conversion.
- `validation/validate_pe_bellhop_freefield_vertical.m` is the orchestration
  entrypoint. It keeps explicit initial-plane and physical-source phase
  references, compares PE/Bellhop/analytic complex fields, saves arrival and
  convergence CSV files, and does not relax a failed tolerance.
- `validation/support/virtual_point_source_initial_field_vertical.m` is enabled only through the
  additive `source_mode='custom_field_fn'` validation hook. Public defaults
  remain `source_mode='gaussian'`.

```matlab
setenv('BELLHOP_EXE','E:/stable/path/to/bellhop.exe')
addpath('scripts/validation')
validation = validate_pe_bellhop_freefield_vertical();
```

Artifacts are under `results/validation/pe_bellhop_freefield/formal/`; the
summary is `reports/pe_bellhop_freefield_validation_report.md`.

### Point-source error-budget follow-up

以下段落中的Bellhop重跑入口已移入
`../cash/pe_bellhop_validation_scripts_2026-09-01/`；对应命令仅作历史记录，当前活动目录不再提供这些入口。

- `validation/validate_pe_point_source_error_budget_vertical.m` runs the
  no-sponge fixed-`dx` window series, compares spatial truncation with a
  discrete spectral-cell Weyl initializer, validates a separate continuous
  Weyl integral, and then scans the complete window/thickness/strength matrix.
  Its CSV explicitly reports `dA=20log10(|H_sponge|/|H_no_sponge|)` and
  `dTL=-dA`, plus center (`rho<=2 m`), ratio-defined edge-band, and total
  terminal-plane energies.
- `validation/support/weyl_point_source_reference_vertical.m` removes the grazing `1/kz`
  singularity by separate propagating/evanescent substitutions; it is the
  reliable free-space reference.
- `validation/support/weyl_point_source_initial_field_vertical.m` is the discrete FFT initializer
  under test. Its failure to converge within the public 100 m window is
  reported, not hidden.
- `../cash/pe_bellhop_validation_scripts_2026-09-01/scripts/validation/rerun_pe_bellhop_point_source_postbudget_vertical.m` performs the
  final representative Bellhop rerun only after the continuous Weyl gate.

Outputs are in `results/validation/pe_point_source_error_budget/` and
`results/validation/pe_bellhop_freefield/post_error_budget/`.

### Production Gaussian window/sponge audit

- `validation/validate_gaussian_window_convergence_vertical.m` calls the
  production Gaussian path, performs the independent PE--AS hard gate, and
  establishes a no-sponge reference using fixed production sampling.
- `validation/validate_gaussian_sponge_vertical.m` evaluates the complete
  50/80/128 m x 6-ratio x 7-strength matrix with explicit `dA`/`dTL`, center,
  edge, total-energy, and edge-to-center diagnostic metrics.
- `validation/validate_gaussian_sponge_wideband_vertical.m` compares the
  3--5 kHz physical `H(f)`, unwrapped phase, and group delay for large
  no-sponge, production/default, and recommended cases.
- `reporting/generate_gaussian_sponge_validation_figures.m` emits the 15
  requested figures plus an optional noiseless LFM diagnostic and the final
  Q1--Q7 report.

Outputs are under `results/validation/pe_gaussian_window_sponge/`; the report
is `reports/pe_gaussian_window_sponge_validation_report.md`. Extended sponge
ratios and windows are validation-only opt-ins; normal public limits and
defaults are unchanged.

### Full reflected PE -> surface -> PE window audit

- `validation/validate_reflected_chain_window_convergence_vertical.m` creates
  one maximum-domain PM surface, crops it without per-window Hs rescaling,
  runs the adaptive 4 kHz no-sponge window sequence, and saves incident,
  reflected, receiver, center, and edge diagnostics.
- `validation/validate_reflected_chain_wideband_vertical.m` runs the gated
  3--5 kHz/33-point reflected and total-channel comparison. Per-case MAT
  checkpoints make the production-size run resumable.
- `validation/validate_reflected_chain_diagnostics_regression_vertical.m`
  checks seeded/override equality, diagnostics transparency, receiver-center
  consistency, the phase-screen invariant, and channel closure.
- `reporting/generate_reflected_chain_window_validation_figures.m` generates
  the single-frequency field and convergence figures; the wideband validator
  adds reflected and total-channel response figures.

Formal outputs are under
`results/validation/pe_reflected_chain_window/`; the reviewable conclusion is
`reports/pe_reflected_chain_window_validation_report.md`. The validators use
additive validation-only options and do not alter production defaults.

### Random-surface reflected-chain window robustness

- `validation/validate_random_surface_window_robustness_vertical.m` generates
  each 256 m master surface once, takes unscaled 160/192 m central crops, and
  runs the 4 kHz three-Hs-by-five-seed matrix with resumable per-case files.
- `reporting/generate_random_surface_window_robustness_figures.m` creates the
  gate scatter, edge-versus-response, pass-rate, and stage-energy figures.
- `validation/validate_random_surface_window_robustness_wideband_vertical.m`
  is a checkpointed 3--5 kHz follow-up for selected worst cases. The current
  formal run was stopped before all candidate--256 m pairs completed; do not
  use partial checkpoints as group-delay qualification evidence.

Formal 4 kHz outputs are under
`results/validation/pe_random_surface_window_robustness/`; the report is
`reports/pe_random_surface_window_robustness_report.md`. The result is
192.1875 m/no-sponge `15/15` strict passes, while 160.15625 m has `15/15`
response/center passes but `0/15` complete edge passes. No production default
is changed by these validators.

- `../cash/pe_bellhop_validation_scripts_2026-09-01/scripts/validation/validate_pe_bellhop_flat_surface_current_vertical.m` is the
  current formal entrypoint. It requires a stable absolute, non-Temp
  `BELLHOP_EXE`, records the binary/source SHA-256 values, runs the independent
  phase audit, 3/6/9 m paths, small-offset limit, source-aware diagnostic,
  sampling/aperture/sponge matrix, public regressions, Bellhop R/C/I fields,
  and the numbered ten-figure atlas.
- `../cash/pe_bellhop_validation_scripts_2026-09-01/scripts/validation/validate_pe_bellhop_flat_surface_vertical.m` generates and runs the historical standard Bellhop ASCII-arrivals case matched to the existing PE model in a uniform medium with a flat pressure-release surface.
- `validation/validate_pe_phase_convention_uniform_vertical.m` independently checks the PE reduced-envelope operator, longitudinal carrier sign, group delay, and validation-local FFT convention against one-step angular-spectrum propagation.
- `../cash/pe_bellhop_validation_scripts_2026-09-01/scripts/validation/validate_pe_bellhop_flat_surface_matrix_vertical.m` runs the historical 3/6/9 m analytic--PE--Bellhop timing/amplitude comparison with an open Bellhop fan, plus the C0--C4 PE grid/window/step convergence matrix.
- `../cash/pe_bellhop_validation_scripts_2026-09-01/scripts/reporting/generate_bellhop_flat_surface_visuals_vertical.m` is the archived historical matrix plotting entry.
- Formal mode requires `BELLHOP_EXE` even if another Bellhop is on the MATLAB
  or system path. Temporary paths and binary-hash changes are rejected.
- The current local baseline is OALIB `2020_11_4` at
  `E:/MISC/BELLHOP/AcousticsToolbox_2020/windows-bin-20201102/bellhop.exe`.
  The retained 2017 directory is historical and must be selected explicitly
  only when reproducing a run whose metadata records that binary.
- The validator runs scalar direct-only/direct-plus-reflection regressions and a 65-frequency PE case, compares direct/single-surface arrival times and TL, reconstructs matched PDPs, and checks public PE invariants.
- Outputs are written to `results/validation/pe_bellhop_flat_surface/`; the Markdown report records exact parameters, formulas, thresholds, results, and limitations.

This stage deliberately excludes rough-surface scattering, bottom bounces, stochastic channels, and communication processing. Validators now consume the public `H_*_reduced_f` and `H_*_physical_f` fields instead of applying an independent hidden carrier convention.

Historical formal command (archived; not runnable from the active tree):

```matlab
setenv('BELLHOP_EXE','E:/MISC/BELLHOP/AcousticsToolbox_2020/windows-bin-20201102/bellhop.exe')
addpath('scripts/validation')
% validate_pe_bellhop_flat_surface_current_vertical is archived.
```

The current run `bellhop_current_20260723_rc5` is `FAIL_CORE` with
`amplitude_status=OPEN`: all phase/path/delay/beam/public checks pass, while
the no-sponge aperture-only fields fail the `-40 dB` edge-validity prerequisite.
Outputs are versioned under
`results/validation/pe_bellhop_flat_surface_current/<run_id>/` and
`results/visualization/pe_bellhop_flat_surface_current/<run_id>/`.

The formal atlas contains 10 numbered PNGs plus MAT, CSV, manifest, and text
summary. Older three-figure outputs remain historical evidence and do not
override the current run ID.

Use environment variables already supported by individual scripts to reduce grid size, seed count, or output file names for quick checks.

## Li et al. (2009) Explicit Rough-Surface Validation

```matlab
run('scripts/validation/validate_li2009_explicit_surface_vertical.m')
```

This independent pure-acoustic workflow uses one shared raw-PM realization
across the full 12 kHz, 6 ms CW frequency synthesis, applies only the explicit
pressure-release Kirchhoff `2*k*eta` screen, and projects the reflected field
through a monostatic bottom--surface--bottom PE path. It excludes noise,
electronics, SSA, kstat, modulation, and BER. Outputs are under
`results/validation/li2009_explicit_surface/`; interpretation and unresolved
transverse-grid sensitivity are documented in
`reports/li2009_explicit_surface_validation_report.md`.

## PE--Bellhop PM model-discrepancy statistical study

The Stage-4 fixed-PM study is split into an execution entry and a read-only
postprocessor:

- `validation/validate_pe_bellhop_pm_model_discrepancy_statistics.m` runs the
  frozen Stage-3B 4 kHz seed ensemble and conditionally extends it when the
  prescribed 24->32 gates fail.
- `validation/postprocess_pe_bellhop_pm_model_discrepancy_statistics.py`
  performs no PE or Bellhop calls; it consumes completed per-seed metric CSVs
  and regenerates the M=50 tables, bootstrap intervals, correlations, outlier
  audit, figures, and flat-array `result.mat` after the user stopped the long
  extension at seed 260050.

The authoritative report is
`reports/pe_bellhop_pm_model_discrepancy_statistics_report.md`.  Current state
is `PRELIMINARY_MODEL_DISCREPANCY`: numerical/applicability guards pass, while
the 24->32 bootstrap half-width gate is narrowly above its engineering limit.
