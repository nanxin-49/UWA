# PE--Bellhop Stage 1 convention-only regression

日期：2026-09-10  
依据：reports/pe_bellhop_controlled_comparison_GOAL_revised_stage1x.md

## 变更范围

主驱动 scripts/validation/validate_pe_bellhop_controlled_comparison.m
仅在 Stage 1 rough/flat 比较层增加
G_BH_comparison = conj(G_BH)。Bellhop/PE 场、internal-wall 几何、
Reflect2D、InfluenceGeoHatCart、receiver selector 和 source pattern
均未修改。原始 Stage 1 FAIL MAT 保留在 stage1/。

## 回归结果

Stage 0 重新执行并 PASS。固定比较约定后的 Stage 1 结果为：

| A (m) | corrected E_G | corrected phase RMS (rad) | corrected TL RMS (dB) | rho_shape | 状态 |
|---:|---:|---:|---:|---:|:---|
| 0.0100 | 0.00409346 | 0.00405930 | 0.00456998 | 0.99999242 | PASS |
| 0.0050 | 0.00204123 | 0.00202407 | 0.00228491 | 0.99999810 | PASS |
| 0.0025 | 0.00101951 | 0.00101096 | 0.00114248 | 0.99999952 | PASS |

A=0.01 通过主驱动以 10,001 beams 重新运行；A=0.005 和 A=0.0025
使用已保存的 raw solver MAT 做 comparison-only 重建，未重复启动 Bellhop。
后两组不新增 solver 收敛证据，但所有既有单项 checks 均重新计算并通过。

## 状态

固定共轭关系解释了三组近似线性 phase residual；该结果支持
CONVENTION_DISCREPANCY_IDENTIFIED。原始 Stage 1 FAIL 不删除、不改写，
Stage 2--7 仍需按 revised Goal 另行解锁。
