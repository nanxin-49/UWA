# PE--Bellhop controlled comparison — Stage 7 historical M=50 interpretation

日期：2026-09-13  
状态：**PASS / INTERPRETED_WITH_HISTORICAL_R_PROVENANCE**  
solver calls / seeds rerun：`0 / 0`

> **Provenance：historical Bellhop point-source R。** 本阶段没有把数据改名为
> X，没有做经验 R→X correction，也没有用 M=50 覆盖新的 X-source controlled
> comparison。

Stage 7 只读取 authoritative `per_seed_results.csv`、`ensemble_statistics.csv`
和 `result.mat`。50 行 seed 范围严格为 `260001--260050`，已发表统计逐项复现，
全部历史 numerical/applicability guards 仍为 PASS。

| M | mean Delta TL (dB) | std (dB) | PE/BH mean power | mean power difference | circular mean phase (rad) | circular std | resultant length |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 50 | -0.00747858 | 0.962232 | 1.004510 / 1.007007 | -0.00249731 | -0.375345 | 0.594701 | 0.837918 |

## Interpretation

- Convention：Stage 1Y 已关闭 active comparison convention；M=50 的非零圆
  均相位不能继续归因于未处理的统一共轭/载波约定。
- Global phase vs spatial distortion：Stage 5 的代表点显示全局 phase alignment
  只降低少量 `E_G`，主要差异为空间失真。
- Roughness validity：canonical fixed-PM 在 Stage 6 落入 Region III，而全部
  Bellhop/PE 数值与几何护栏通过。
- Power：ensemble mean reflected powers 接近，只表示有符号的 realization-level
  幅度差在平均中抵消；不等价于逐 realization 复场闭合。
- Source provenance：R-source 的历史属性限制它与新 X-source Stage6 的直接统计
  合并，但不能通过未经运行的经验修正予以消除。

因此 M=50 支持“具有历史 source provenance 限制的 reflection-model
discrepancy”解释；它不支持把剩余差异重新归类为数值、receiver selection 或统一
phase convention 错误。

产物位于
`results/validation/pe_bellhop_controlled_comparison/stage7_historical_m50_interpretation/`。
