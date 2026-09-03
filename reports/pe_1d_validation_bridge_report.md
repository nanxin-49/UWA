# PE 一横向维 validation bridge Stage 0B 报告

状态：**PASS**

本阶段只验证一横向维 split-step PE、独立一步 exact angular spectrum，以及生产二维横向 PE 的 y-invariant ky=0 投影；没有运行 Bellhop，也没有修改生产 PE。

## 配置

- f/c：4000 Hz / 1500 m/s；Tx/Rx：100 / 3 m；sigma：0.3 m
- PE bridge window/grid：192.1875 m / 984；y cross-check window/grid：50 m / 256；sponge off
- stepz_lamb：0.5；weak test eta amplitude：0.001 m

## Case results

| case | 1D-vs-AS direct | 1D-vs-AS reflected | full-PE ky0 direct | full-PE ky0 reflected | flat sign |
|---|---:|---:|---:|---:|---:|
| flat | 8.96218e-13 | 1.1363e-12 | 7.18587e-13 | 9.45919e-13 | 1 |
| weak_phase_screen | 8.96218e-13 | 1.13801e-12 | 7.18587e-13 | 9.47741e-13 | NaN |

## Checks

- flat_one_d_as: value 8.96218e-13, limit 1e-11, PASS
- weak_one_d_as: value 8.96218e-13, limit 1e-11, PASS
- flat_ky0_direct: value 7.18587e-13, limit 1e-10, PASS
- flat_ky0_reflect: value 9.45919e-13, limit 1e-10, PASS
- weak_ky0_direct: value 7.18587e-13, limit 1e-10, PASS
- weak_ky0_reflect: value 9.47741e-13, limit 1e-10, PASS
- flat_reflection_sign: value 2.64377e-13, limit 1e-12, PASS

## Interpretation

The 1-D helper uses the same reduced square-root propagator and pressure-release phase-screen convention as the production chain. The full PE comparison is made after integrating the receiver plane over y, which selects the conserved ky=0 mode; the expected source factor is the discrete y Gaussian integral and is not fitted. This validates the dimensional bridge, not PE--Bellhop rough-surface agreement.
