# 2026-08-20 目录迁移记录

本次迁移采用兼容优先策略：根目录公共函数/脚本名称保留，实际实现移入 `src/` 或 `examples/`。验证专用辅助函数移入 `scripts/validation/support/`，2k 产物移入登记结果目录，旧规格移入 `docs/history/`。

完整逐文件旧路径、新路径和迁移前 SHA-256 见 `repository_layout_migration_20260820.csv`。哈希在移动前、包含当时所有未提交修改的工作树上采集；后续为适应新路径所做的入口改名、setup 调用、fingerprint 更新和文档来源更新会使部分当前哈希不同，这是预期行为。

保护规则：

- 未执行 `git reset`、`git checkout --` 或删除用户修改；
- 未清理历史 MAT、PNG、CSV、checkpoint 或正式报告；
- 根目录公共 API 与三个演示调用名称不变；
- 未改变 PE marching、海面/相位/统计公式、公共默认 surface、窗口或 sponge；
- 历史日志内旧路径保留原文，由 `PROJECT_CONTEXT.md` 顶部当前索引说明新位置。

注意：`weyl_point_source_reference_vertical.m` 的权威迁移前哈希以任务开始时的保护快照为准；若 CSV 与后续自动审计发现不一致，应以保护快照重新核对，不得覆盖源文件。
