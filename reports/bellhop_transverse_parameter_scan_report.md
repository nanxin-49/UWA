# Bellhop 横向参数扫描

仅运行 Bellhop，固定 8000 Hz、10001 beams，并与独立角谱比较；没有重新运行 PE。扫描 27 cases：step=[0.1 0.05 0.025] m，half-angle=[20 30 45] deg，SBP samples=[1201 2401 4801]。Bellhop 横向输出按正式验证的相位约定取共轭，此为固定 convention 转换而非逐点拟合。

## 最优 case

- step=0.1 m，half-angle=30 deg，SBP samples=4801。
- max TL=0.015867013 dB，max phase=0.0612039 rad，max complex error=0.061165714。
- 通过 case 数：0/27。

No scanned Bellhop-only configuration satisfies all transverse limits. The remaining discrepancy is not explained by the tested step, angle, or SBP sampling ranges; continue with the 2-D/3-D source-mapping audit.

详细数据：results/validation/pe_as_bellhop_transverse/bellhop_scan/bellhop_transverse_parameter_scan.csv
