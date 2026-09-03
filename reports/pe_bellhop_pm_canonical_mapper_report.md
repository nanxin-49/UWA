# PE--Bellhop PM canonical mapper Stage 0A 报告

状态：**PASS**

本阶段只验证固定 Fourier PM realization 的跨 solver 数学映射，没有运行 PE/Bellhop propagation，也没有修改任何核心物理。

## 输入与映射

- coefficients：E:\MISC\CARPE3D_matlab\Explain\results\validation\bellhop_internal_pm_fixed_realization\fixed_pm_fourier_coefficients.csv
- coefficients SHA-256：1f7eda465e4ae85b8ac038310edf053f2d062563a3b943bfd83bd59d07015f67
- seed/U/span：260001 / 6 m/s / 160 m
- coefficient count / requested Kmax / realized Kmax：12 / 0.5 / 0.471238898038 rad/m
- series：eta(s)=sum_n(cos_coeff_m(n)*cos(k_n*s)+sin_coeff_m(n)*sin(k_n*s))
- master profile SHA-256：a0ee79059297f16dcabac44b0eec4025b11eaabad94957e533dce4d87105dd33
- PE evaluation interval/count：[-96.0938, 96.0938] m / 984
- Bellhop evaluation interval/count：[-80, 80] m / 4097
- PE/Bellhop sampled profile hashes：38d129e4086286224dd182bd8e1b2717d05867764589efbba76ec6d4c6dc89c1 / 862f2096ab4831b9c190a7a8e388802daa3f776845ae484210b8d6443b5208a8
- sampled hashes intentionally differ because the solver grids differ; the common canonical identity is the coefficient-file SHA-256 above.
- PE mapping：eta_PE(x,y)=eta_1D(x) on the configured x samples.
- Bellhop mapping：Gamma_B(s)=[R0-eta(s),s], R0=100 m; no sign change or renormalization.

## 误差与统计

| quantity | value |
|---|---:|
| master height max error | 8.38218e-15 m |
| master slope max error | 2.01401e-15 |
| master second-derivative max error | 5.30825e-16 1/m |
| Bellhop wall inverse-map error | 7.10543e-15 m |
| RMS height | 0.18639076 m |
| RMS/max slope | 0.0539798126 / 0.116674298 |
| RMS/max curvature | 0.020904738 / 0.0445943842 1/m |
| minimum radius | 22.4243482 m |
| endpoint height/slope/second-derivative differences | -6.27276e-15 / -3.16414e-15 / 5.55112e-16 |

## Hard checks

- master_height: PASS
- master_slope: PASS
- master_second_derivative: PASS
- pe_finite: PASS
- bellhop_finite: PASS
- pe_inverse: PASS
- bellhop_inverse: PASS
- wall_mapping: PASS

## 结论

同一组 Fourier coefficients 在独立 PE/Bellhop 采样点上直接求值；Bellhop wall 的 r=R0-eta 反变换误差为 7.10543e-15 m。本阶段不证明两种传播模型的场一致性，只冻结后续 Stage 0B--1 的 canonical 输入。
