# 二维活塞 P0 验证交付

2026-11-12，run `20261112_remote_p0p1_01`。

P0通过，三频4/6/8 kHz、总宽20/10 cm的同一速度孔径在PE-mode 0、BIE及FullFK所需的压力与梯度完成检查。数值初场1/2 m偏移和phase contract v2载波恢复已核验。Gaussian标量默认回归差0。

P1三环境与独立审计完成。方法、所有冻结门、失败与修正、主要模型差异、20 cm主源结论、交付索引及适用范围见[最终报告](../../../comparisons/piston_2d_single_frequency/20261112_remote_p0p1_01/report.md)。

- `jobs/`保存原始资格与回归MAT、CSV、回执。
- `independent_audit_02/`保存接受的独立数值审计，失败0项。
- `remote_provenance/`保存快照、日志、PID、真实退出码和资源。
- `transfer_integrity.json`核验首批1373项远端/本地SHA。
- `checkpoint_final.json`与`delivery_manifest_sha256.csv`为最终入口。

初期有限窗口梯度失败、周期PE传播失败、许可等待和后续恢复均保留。最终采用验证专用连续角谱及匹配外部平面场表示。原公开WAPE核心与Gaussian生产路径沿用原实现。
