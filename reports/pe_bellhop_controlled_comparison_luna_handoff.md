# PE--Bellhop controlled-comparison Luna handoff

本文件是执行交接摘要，不替代唯一权威规范
`reports/pe_bellhop_controlled_comparison_GOAL.md`。

## Current state

- P0：PASS；source `X`、coherent run `C`、official AcousticsToolbox 2020
  executable 和当前 validation-only binaries 的 provenance 已记录在
  `results/validation/pe_bellhop_controlled_comparison/p0_report.md`。
- Stage 0：PASS；authoritative report 为
  `results/validation/pe_bellhop_controlled_comparison/stage0/stage0_report.md`。
  全尺寸配置为 `nx=984`、`W=192.1875 m`、PE/Bellhop step `0.05 m`、
  5001→10001 beams。PE--AS L2 为 `4.2785525e-13`，flat internal-wall
  M99 L2 为 `0.0058993922`，beam-convergence L2 为 `5.854857e-06`，
  几何/phase/tau/positive-range/finite gates 全部通过。
- Stage 1：未完成。`A=0.01 m, K=0.10 rad/m` 的 generic parametric
  internal-wall 高 beam case 在受控运行窗口内未结束，已停止；这不是
  physical PASS/FAIL，也不应解释为 PE--Bellhop model discrepancy。
  当前 stage1 状态见
  `results/validation/pe_bellhop_controlled_comparison/stage1/stage1_report.md`。

## Resume rules

1. 先重新核对 Stage 0 report 的 `passed=true` 和 Goal gates。
2. 只运行 Stage 1；保持 primary `10001→20001` beams、profile `N=2049→4097`、
   source `X`、run `C`、同一 receiver map 和 frozen AS masks。低 beam smoke
   结果不得替代 primary gate。
3. 若运行完成性再次阻塞，记录为 blocked/incomplete 并停止，不启动 Stage 2。
4. 只有 Stage 1 weak-limit gates 全部通过后才可顺序进入 Stage 2--7。

入口：

```matlab
setup_vertical_project
addpath(fullfile(pwd,'scripts','validation'))
result = validate_pe_bellhop_controlled_comparison('stage1');
```

禁止读取或恢复 `cash/`，禁止修改 PE/Bellhop 核心、`Reflect2D`、
`InfluenceGeoHatCart`、SHD selector、通信代码或引入额外散射公式。
