# 论文与答辩交付材料归档清单（2026-11-09）

按用户要求，将与 MATLAB 声传播/信道/通信项目主线无关的交付包和 PPT 过程材料移至 `cash/2026-11-09_proposal_thesis_collateral/`。本清单记录移动前的源路径与文件数量；`cash/` 内材料不作为当前验证证据。

| 原路径（归档目录保留相同相对路径） | 文件数 | 字节数 |
|---|---:|---:|
| `results/visualization/Chapter4_Final_Review` | 14 | 1973498 |
| `results/visualization/Chapter4_Final_Review_v0.4` | 16 | 1980531 |
| `results/visualization/proposal_research_foundation_20261102` | 28 | 8561221 |
| `reports/研究基础PPT三页图文规划_待确认_20261102.md` | 1 | 12936 |

移动前核对发现：PPT 目录的 `README.md` 声称草稿与成品字节一致，但两份 PPTX 的 SHA-256 不同，内部文件也有差异；此处不判定哪份应作为最终稿。

## 同日续归档

按用户后续要求，将上次列出的其余材料及其对应的本地生成输出继续移入同一目录。下表数量和字节数在移动前记录。

| 原路径（归档目录保留相同相对路径） | 文件数 | 字节数 |
|---|---:|---:|
| `reports/ch4_figure_inventory.md` | 1 | 4061 |
| `reports/ch4_figure_redraw_receipt.md` | 1 | 2793 |
| `reports/ch4_method_audit.md` | 1 | 3950 |
| `reports/ch4_p22_evidence_boundary.md` | 1 | 3797 |
| `reports/ch4_qpsk_error_floor_diagnosis.md` | 1 | 3234 |
| `reports/ch4_terminology_reference_crosscheck.md` | 1 | 4416 |
| `reports/开题报告_v0.3.1_术语图表代码一致性核对_2026-11-02.md` | 1 | 7383 |
| `reports/开题报告补充核验报告.md` | 1 | 29454 |
| `reports/第四章_审校意见_配图规范与Codex核验_v0.3.md` | 1 | 24201 |
| `reports/第四章_已取得的研究基础_修订稿_v0.3.md` | 1 | 20420 |
| `results/visualization/ch4_verified_figures` | 90 | 12778093 |
| `results/visualization/p22_known_csi_flow` | 3 | 164928 |
| `results/visualization/proposal_supplement_20261102` | 7 | 2571783 |
| `results/visualization/target_dsb_near_bottom_schematic` | 9 | 731751 |
| `results/validation/ch4_qpsk_error_floor` | 6 | 8561 |
| `scripts/reporting/render_ch4_verified_figures.m` | 1 | 14523 |
| `scripts/reporting/render_proposal_figure31.js` | 1 | 3228 |
| `scripts/reporting/render_target_dsb_near_bottom_schematic.m` | 1 | 4102 |
| `scripts/validation/diagnose_two_node_qpsk_error_floor_ch4.m` | 1 | 5385 |

`.gitignore` 和 `scripts/README.md` 的未提交版本已分别复制到归档目录的 `tracked_worktree_overrides/` 对应路径；工作区文件随后由 HEAD 字节恢复。两份 v0.2 稿件在本轮开始时已从工作区删除，没有可移动的源文件；本次清理将它们从当前 Git 版本移除，旧版仍可从 Git 历史恢复。

本轮未处理新近修改的 `AGENTS.md`、`PROJECT_CONTEXT.md` 和新增的 SSH MATLAB 交接报告。
