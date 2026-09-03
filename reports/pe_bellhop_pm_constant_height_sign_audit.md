# PE--Bellhop PM Stage 0D constant-height sign audit

状态：**PASS**

本阶段使用 eta0 = +0.05, 0, -0.05 m，保持 Gaussian `.sbp`、uniform c=1500 m/s、pressure-release 和既有 internal-flat POC；不拟合源强或修改 PE/Bellhop 核心。

## Results

| eta0 (m) | expected phase | PE ratio phase err | Bellhop ratio phase err | PE--BH phase err | PE TL ratio (dB) | BH TL ratio (dB) | path time err (s) | wall vs -reference phase |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.05 | 1.6755161 | 1.9706459e-15 | -2.5569633e-05 | 2.5569633e-05 | -5.7859648e-15 | 0.0084364111 | -4.1633363e-17 | 0 |
| 0 | 0 | 0 | 0 | 0 | 0 | 0 | -4.1633363e-17 | 0 |
| -0.05 | -1.6755161 | -9.15934e-16 | 2.558655e-05 | -2.558655e-05 | -6.7502923e-15 | -0.0084295459 | -4.1633363e-17 | 0 |

## Checks

- pe_sign: PASS
- bellhop_sign: PASS
- cross_model_sign: PASS
- phase_magnitude: PASS
- pressure_release_once: PASS
- no_extra_pi: PASS
- path_time: PASS
- geometry: PASS
- rotation_state: PASS
- all: PASS

The physical image span is L=103-2*eta0 m. The PE screen predicts exp(+i 2 k eta0); Bellhop moves the wall to R0-eta0 and reads the transformed branch at the same L. The known Cartesian backward-range amplitude offset is not used as a gate.
