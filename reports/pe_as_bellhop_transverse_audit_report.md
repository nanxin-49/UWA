# PE--AS--Bellhop 横向剖面审计

本审计读取正式 MAT，不重新运行 PE 或 Bellhop。配置为 c=1500 m/s、Gaussian sigma=0.3 m、PE 窗口 192.188 m，频率为 4/6/8 kHz，距离为 97/103 m。

## 最大误差

| 比较 | TL (dB) | phase (rad) | complex |
|---|---:|---:|---:|
| PE-AS direct | 3.5496894e-12 | 3.0494248e-13 | 5.0989432e-13 |
| PE-AS reflect | 1.8804386e-12 | 1.594129e-13 | 2.4339886e-13 |
| Bellhop-AS direct | 0.015822171 | 0.061203913 | 0.061165731 |
| Bellhop-AS reflect | 0.012736042 | 0.051673012 | 0.051650184 |
| PE-Bellhop direct | 0.015822171 | 0.061203913 | 0.061277252 |
| PE-Bellhop reflect | 0.012736042 | 0.051673012 | 0.051725974 |

PE-AS pass: true; Bellhop-AS pass: false; PE-Bellhop pass: false.

PE agrees with the independent angular-spectrum transverse reference; the remaining PE--Bellhop transverse discrepancy is attributed to Bellhop/source mapping or the 2-D ray-beam representation, not the PE marching operator. The audit is diagnostic and does not relax formal thresholds.

下一步：PE-AS 通过时继续 Bellhop 步长、角扇区和 SBP 映射审计；PE-AS 失败时先检查场提取、坐标和相位约定。

![audit](../results/validation/pe_as_bellhop_transverse/pe_as_bellhop_transverse.png)
