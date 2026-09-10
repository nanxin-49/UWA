# PE--Bellhop controlled comparison — Stage 1 status

日期：2026-09-10  
权威执行入口：`scripts/validation/validate_pe_bellhop_controlled_comparison.m`  
权威 Goal：`reports/pe_bellhop_controlled_comparison_GOAL.md`

## 已执行配置

4 kHz、均匀 `c=1500 m/s`、X 线源、coherent C、固定弱正弦
`eta(s)=A cos(0.1 s)`。Stage 0 已 PASS。根据用户最新运行约束，后续
Bellhop case 统一使用 10,001 beams；不再启动 20,001-beam。

## Stage 1 结果

| A (m) | profile L2 (M99) | beam L2 (M99) | model L2 (M99) | model phase RMS (rad) | model TL RMS (dB) | rho | 状态 |
|---:|---:|---:|---:|---:|---:|---:|:---|
| 0.0100 | 1.2025e-6 | 1.5437e-6 | 0.47247 | 0.47920 | 0.004570 | 0.91923 | FAIL |
| 0.0050 | 6.0272e-7 | 1.4363e-6 | 0.23875 | 0.23960 | 0.002285 | 0.97928 | FAIL |
| 0.0025 | 3.0923e-7 | 0 (same 10,001-beam endpoint) | 0.11969 | 0.11980 | 0.001142 | 0.99479 | FAIL |

三组 case 的 profile/geometry guards 均通过；wall residual 约 `7.1e-15 m`，
pressure-release phase jump 误差约 `7.1e-15 rad`，最小 transformed `dr`
约 `0.0433 m`，center travel-time 误差约 `4.7e-15 s`，且无 NaN/Inf。

Stage-0 floor 给出的弱限值为 `T_E=0.0108994`、`T_phi=0.0108992 rad`、
`T_TL=0.0107979 dB`。即使 `A=0.0025 m`，`E_G`、phase RMS 和 `rho_shape`
仍未达到 floor gate；TL RMS 已通过。随着 A 减半，模型 L2/phase 近似减半，
数值收敛项保持在 `1e-6` 量级，说明当前阻塞是模型可比性而非 solver
收敛或 internal-wall 几何故障。

## 决策

按唯一权威 Goal：Stage 2--7 锁定，停止 height/slope/PM 扩展，转入
conditional Helmholtz/BEM feasibility branch。A=0.01 与 A=0.005 的历史
20,001-beam endpoint 结果保留为诊断；A=0.0025 遵循最新约束仅使用 10,001
beams，不能宣称原始 10,001→20,001 收敛门已重新执行。

详细逐 case 文件：

- `results/validation/pe_bellhop_controlled_comparison/stage1/A0p01/`
- `results/validation/pe_bellhop_controlled_comparison/stage1/A0p005/`
- `results/validation/pe_bellhop_controlled_comparison/stage1/stage1_report.md`

本轮未修改 PE、Bellhop 核心、`Reflect2D`、`InfluenceGeoHatCart`、SHD
selector 或通信代码。
