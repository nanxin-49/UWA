# PE--Bellhop 反射前入射复声场比较

状态：**PASS**

本报告只比较均匀介质中 `r=100 m` 的反射前 incident plane；不启用 PM、海面反射、internal wall、`Reflect2D` 或 PE 生产二维横向场。

## 配置

- f/c：4000 Hz / 1500 m/s；Gaussian sigma：0.3 m；PE window/grid/step：192.1875 m / 984 / 0.05 m。
- Bellhop：`OALIB AcousticsToolbox 2020_11_4 (2020-11-02 binary)`；run type：`C`；source：`X`；beams：5001 10001；step：0.05 m；angle fan：[-30,30] deg；SHD ranges：[99,100] m。
- executable：`E:\\MISC\\BELLHOP\\AcousticsToolbox_2020\\windows-bin-20201102\\bellhop.exe`；source-pattern fingerprint：`X|cos(theta)*exp(-0.5*(k*sigma*sin(theta))^2)|N=2401|angle=[-30,30]deg|f=4000Hz|sigma=0.3m`。
- footprint：M95 radius 28.7109375 m (295 samples, energy 0.950895936); M99 radius 38.8671875 m (399 samples, energy 0.990139665)。

## 复相位与归一化

PE uses the reduced-envelope one-step convention `exp(i*R*(kz-k0))` for the incident field. Bellhop raw pressure uses the saved normalization-audit spatial sign `-1` (fixed conjugation); no data-dependent conjugation, carrier re-addition, phase fitting, or empirical amplitude correction is applied. Main comparison is axis-normalized at x=0. For line-source X, the historical point-source global source constant is not applied.

## 结果摘要

| comparison | L2 M99 | phase RMS M99 (rad) | phase P95 M95 (rad) | TL P95 M95 (dB) | LS scalar abs | LS scalar phase (rad) |
|---|---:|---:|---:|---:|---:|---:|
| pe_as | 3.7313678e-13 | 2.9570732e-13 | 4.9708963e-13 | 3.9325274e-12 | 1 | 2.4844686e-13 |
| bellhop_as_5001 | 0.0061022184 | 0.0061020139 | 0.01429134 | 0.0019036146 | 0.99999138 | 0.0020724906 |
| bellhop_as_10001 | 0.0061050369 | 0.0061048329 | 0.014293065 | 0.0019036364 | 0.99999144 | 0.0020770223 |
| beam_convergence | 5.2520122e-06 | 5.2507283e-06 | 8.7537556e-06 | 1.948513e-06 | 0.99999994 | -4.5317232e-06 |
| pe_bellhop | 0.0061048459 | 0.0061048329 | 0.014293065 | 0.0019036364 | 1.0000007 | -0.0020770223 |

PE outer-5% incident energy: `3.0394899e-06`; receiver coordinate error max: `0 m`; non-finite count: `0`.

Half-dx receiver diagnostic: `984` shared nodes; maximum coordinate error `0 m`; shared-node normalized L2 M99 `0`; phase RMS M99 `0 rad`.

Symmetry diagnostic (not a hard gate): 983 paired samples for PE/Bellhop/AS; maximum normalized complex pair error `3.76310e-13 / 0 / 8.12532e-15`.

## Hard gates

| check | value | limit | pass |
|---|---:|---:|:---:|
| pe_as_max_complex_full | 5.6156379e-13 | 1e-10 | 1 |
| pe_as_l2_m99 | 3.7313678e-13 | 1e-10 | 1 |
| pe_outer5_energy | 3.0394899e-06 | 1e-05 | 1 |
| receiver_coordinate_error | 0 | 1e-06 | 1 |
| bellhop_beam_l2_m99 | 5.2520122e-06 | 0.002 | 1 |
| bellhop_beam_phase_rms_m95 | 5.2036899e-06 | 0.005 | 1 |
| bellhop_beam_tl_rms_m95 | 8.6875045e-07 | 0.02 | 1 |
| pe_bellhop_l2_m99 | 0.0061048459 | 0.02 | 1 |
| pe_bellhop_phase_rms_m99 | 0.0061048329 | 0.02 | 1 |
| pe_bellhop_phase_p95_m95 | 0.014293065 | 0.05 | 1 |
| pe_bellhop_tl_p95_m95 | 0.0019036364 | 0.1 | 1 |
| receiver_sampling_l2_m99 | 0 | 0.005 | 1 |
| nonfinite_count | 0 | 0 | 1 |

4 kHz incident-plane gates pass. Conditional next step: repeat the same frozen convention and masks at 6/8 kHz.

![incident field comparison](../results/validation/pe_bellhop_incident_field/pe_bellhop_incident_field.png)
