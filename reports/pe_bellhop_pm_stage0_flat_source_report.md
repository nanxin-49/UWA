# PE--Bellhop PM Stage 0C flat source/normalization audit

状态：**PASS**

本阶段只运行 flat reflected/direct geometry；不使用随机 PM、不拟合 SBP、 不比较绝对源强。PE 是 Stage 0B 的一横向维 bridge，Bellhop 使用当前 2020 internal-flat validation binary。

## 配置

- f/c：4000 Hz / 1500 m/s；Tx/Rx：100 / 3 m；direct/image：97 / 103 m
- PE window/grid：192.1875 m / 984；step：0.05 m；sponge off
- Bellhop beams：[5001 10001 ]；step：0.05 m；SBP：2401 points; receiver offsets：11

## Results

| beams | axis Q TL dB | axis Q phase rad | axis Q complex | direct profile mag | reflect profile mag | direct profile complex (diag) | reflect profile complex (diag) | wall/image complex |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 5001 | 0.2609997 | -0.0005149472 | 0.030509198 | 0.0059687552 | 0.0055544581 | 1.7083985 | 1.8110889 | 0 |
| 10001 | 0.2609997 | -0.0005149472 | 0.030509198 | 0.0059687686 | 0.0055544471 | 1.7083973 | 1.8110878 | 0 |

## Native wall invariants

intersection residual max: 3.5385028e-12 m; kappa max: 0 1/m; pressure-release phase error: 7.1054274e-15 rad; Amp jump: 0; p/q rotation errors: 0 / 0; travel-time error: -4.1633363e-17 s; minimum transformed range increment: 0.04330127 m.

## Checks

- path_and_phase: PASS
- wall_geometry: PASS
- pressure_release: PASS
- rotation_state: PASS
- positive_range: PASS
- source_axis: PASS
- source_profiles: PASS
- beam_convergence: PASS
- all: PASS

## Interpretation

The flat axis Q and normalized offset magnitude profiles are the source/directivity audit for the later PM ratio. The complex offset-profile residual is retained as a diagnostic because the PE bridge and Bellhop unfolded chart use opposite transverse phasor orientation; no conjugation, source fitting, or calibration is applied. The reported absolute axis scale ratios are diagnostic only. If this stage passes, Stage 0D constant-height sign audit is the only next step.
