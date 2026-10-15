# 深海研究计划 v3 合并核对（2026-10-14，Asia/Shanghai）

正式依据仍为 [deep_sea_channel_research_plan.md](deep_sea_channel_research_plan.md)。本报告核对本地原计划与用户指定的 Downloads/v3 附件；不把附件中的执行命令视作独立授权。用户要求完成核对后执行双轨 Goal，其范围由正式计划、用户不可变约束及配套执行记录共同限定。

## 备份与来源

- 原文件及完整逐字节备份 [backup](deep_sea_channel_research_plan_backup_2026-10-14.md) 的 SHA256 均为 `37E11B845A9CEA0D9E787D73146B3C35CDA8EB368D1A1354274E8A14FEDF7E22`。
- [v3 版本记录](deep_sea_channel_research_plan_v3_2026-10-14.md) 是 Downloads 附件的原样复制，仅作版本记录，不同步修改。
- 配套 [Goal 执行范围](codex_goal_two_track_2026-10-14.md) 来自 Downloads 中同名文档；已将“采用附件 v3”改为“采用唯一正式计划”。它不是第二份研究计划。
- [完整合并差异](../results/validation/deepsea_two_track_20261014/plan_merge.diff) 保留替换前后的证据。

## 完整性核对

| 项目 | 核对结论 |
|---|---|
| 物理定义 | 三维点源、两坐标海面、z 向下、压力释放、均匀无损半空间、exp(-iωt)、源强归一化、direct/reflected 同源、q 定义及船体/海底/吸收限制保留 |
| 公式 | 原文全部 9 个显示公式逐字符保留（仅允许 CRLF/LF 差异）：Green、梯度、DSP 时延参考、路径+残差、镜像场与 Ld/Lr、演化海面、两项成本和 break-even |
| WP0–WP6 | 全部保留；D0–D6/G0–G6、32 响应/96 标量评价、257/513 频点、开发/验证各 32 seed、QPSK 提议与动态适用限制保留 |
| 门槛 | 第 9 节完整原验收表保留：代数 1e-10、数值 1e-3、物理 1%、功率 .5dB/统计约10%、BER=1e-3 的 .5dB、PER=.10 附近差 .02；仍是提议而非本轮已验收 |
| 历史结果 | 原计划从 2026-10-13 matrix-free 记录至第 1 节前的整个历史块逐字符保留；未删除早期失败、修订或阶段性 INCONCLUSIVE |
| 历史基准与当前状态 | L0/direct restart8/12、96 内迭代、PM PDE/DtN/边界失败值保留；dual 仅离散案例条件通过；G1 OPEN/WP2=false 保留 |
| 统计与资源 | 独立单位、固定父海面、数据切分、禁止测试集调参、成本口径与不确定性限制保留 |

## 差异、冲突与处理

1. **下一步任务范围发生实质收紧。**原 WP1 的至多 12 案例与 WP2 的 32 响应属于未来设计；v3 明确当前仅一次有界预条件干预及确定性通信准备。保留原设计和历史原文，以最新决策覆盖旧“下一步建议”，不重跑旧 pilot，不扩大 PM。
2. **统计路线偏好改变。**原文偏好 S2；v3 将 P0/S1/S2 改为证据驱动选择。模型公式与三个比较路线没有丢失。表中历史“Preferred mechanism-based hypothesis”只表示原假设，不能凌驾于正文“不预选胜者/尚未授权拟合”。
3. **时区元数据不一致。**附件写 Asia/Singapore，原项目与用户环境为 Asia/Shanghai。正式文件统一为 Asia/Shanghai；原始附件记录保留原字符串。二者同为 UTC+8，不改变实验时间或结果。
4. **配套文件原先不在仓库。**找到 Downloads 中 `codex_goal_two_track_2026-10-14.md`，纳入 reports 并修正正式依据引用，解决悬空配套指向。
5. **没有发现物理公式、数值门槛或历史数值的矛盾替换。**这项结论限于新旧计划的文档核对，不是重新认证全部历史实验，也不重新检索外部文献。附录中的历史参考信息按原记录保留。

本轮数值门仍从冻结 mapped-direct 协议读取：线性 1e-8、真实 PDE/Dirichlet/DtN 1e-4；不能用第 9 节一般提议替代。新任务预算单独预声明，历史协议不改。

## 结构、引用与版本检查

检查全部本地 Markdown 链接、代码围栏和显示数学定界符、WP0–WP6 与第 1–13 节、备份和版本来源哈希及正式依据引用。机器结果保存在 [plan_structure_audit.json](../results/validation/deepsea_two_track_20261014/plan_structure_audit.json)。外部 URL 沿用历史来源，本轮未作外网可达性复查；不存在将 v3 记录作为另一正式计划的维护入口。

该合并只更新研究政策和导航；生产公式、公共接口、冻结求解器和历史数据不变。双轨研究独立记录结果，任何基础设施 PASS 均不释放 G1/WP2。
