# 项目进度记录怎么用、怎么更新

先读[进度文字版](project_progress.md)，或双击本机的[离线看板](../results/visualization/project_progress/index.html)。七条工作线的第一句话是快速回顾，展开后才是详细过程与证据。[原二维工作台](../UWA_2D_research_review_workbench.html)继续保留为专题页。

当前执行顺序先核 [AGENTS.md](../AGENTS.md) 与[任务表](active_task_tracker.md)：后续主声源为总宽10 cm的二维条带，20 cm和Gaussian结果保留历史源条件；4 kHz PM比较 → 反射模型决策 → 相干4–8 kHz/CIR → 通信。看板只反映状态，不授予下一阶段、SSH或计算预算；Bellhop保持独立并行。

## 每次项目任务结束时做什么

1. 每次实质任务结束先更新任务表中受影响行和最后行动；仅当项目进度、接口、阻塞或当前研究导航实际改变时，同步本看板。看板唯一手工源为 `reports/project_progress.json`。修改受影响的工作线，说明“为什么做、做了什么、说明什么、还差什么、下一步是什么”。结果时间写在该工作线的 `result_date`；实际检查过的工作线才更新其 `checked_at`，仅静态导航核对需明确标注，不得当作重做数值审计。顶层 `checked_at` 记录最近一次整理日期，不能用它给未重查的其他工作线盖新日期。`updated_at` 是记录更新时间，页面页脚另有生成时间。
2. 给 `events` 追加一条有唯一 ID 的事件。失败、暂停和停止也要记，旧事件不删除、不改成通过。`workstream_id` 指向已有工作线；事件保持日期顺序。同日可以有多条。
3. 已有 MATLAB MCP 会话可用时，执行以下命令。它会刷新两种视图和 `render_receipt.json`，不会运行声场计算、加载 MAT/FIG、读取归档、恢复脚本或接触远程机器。

```matlab
setup_vertical_project;
output = render_project_progress_dashboard();
```

4. 检查回执中的 `failed_checks` 和 `unexpected_references`。缺失入口若已明确登记，`expected_exists=false`；若原本应存在的引用消失，或缺失入口后来恢复，应提示复核并更新记录。文件存在不等于内容正确或数值通过。

维护源和 MATLAB 入口、HTML 模板纳入 Git。生成的 `reports/project_progress.md` 是方便仓库阅读的快照，也纳入 Git；`results/visualization/project_progress/` 的 HTML 和回执是本地生成物，不需要提交整个结果树。Git checkout 后可重新生成 HTML。不要手工修改生成视图，也不要为保存这个看板恢复其他历史脚本。

## 用户明确禁止MATLAB时的纯文本同步（2026-11-14）

本次用户限定为配置维护、不运行MATLAB、不连接SSH、不启动实验或提交Git。因此不探测MCP，也不推断其当前可用性；按AGENTS的文本产物打包例外，使用同一JSON、HTML模板及模板内既有 `progressMarkdown` 函数刷新HTML/Markdown，不手改生成视图。该步骤只校验JSON、相对路径、明确文件存在性及文本/嵌入数据一致性，禁止调用MATLAB、Git历史数值复查、网络或读取 `cash/`。SSH Skill不因配置维护或MCP故障启用。

每次这种同步须使用独立 `results/validation/project_progress_dashboard/<run_id>/` 留存输入/代码/输出SHA、实际stdout/stderr与退出状态、文本检查及检查点，并明确MATLAB渲染、Code Analyzer和浏览器视觉/交互检查未运行。引用缺失须原样报告，不能移入 `cash/` 查找或将缺失改成数值通过；新配置链接必须存在。当前生成回执标明纯文本模式，不覆盖历史run回执。恢复普通维护后仍按上方MATLAB正式入口执行。[本次修复与回执](agents_configuration_alignment_2026-11-14.md)。

## 最小数据约定

- `schema_version=1`；每条工作线固定包含目的、过程、结论范围、缺口、下一步、技术范围、已有结果、历史审计和当前重跑条件。
- 内部状态只用 `closed/partial/paused/pending/stopped`，对外标签用中文，且写清适用边界。状态由证据判断，不根据文件存在性自动升级。
- `evidence` 使用仓库相对路径，含显示名称、类型、可选历史版本及预期是否存在。禁止绝对路径、越界路径和根目录 `cash/` 引用。当前缺失入口不生成失效链接，改为显示路径和已确认的历史提交。
- 历史版本的存在、输出哈希、审计适用范围需要在相应任务中实际核查；生成工具只检查精确路径存在性，不执行 Git、数值或历史版本审计。
- JSON/HTML 按 UTF-8 保存；模板内嵌 JSON，不使用 CDN、网络字体或请求接口。页面可在 `file://` 下打开。

## 2026-11-10 首版验收与执行限制

基线为 `s5s6`、`3025b35c13aaa34fcb5c40983150eb2c74c664d4`。本次核查当前代码入口、正式报告、保存的检查点/审计摘要和精确引用，未重做历史科学数值审计。特别保留：原始 IFFT 失败、共同窗表示通过、已知 CSI 控制范围、三维参考停止决策、Bellhop 的空间相位/复场未过项、目标 PE 孔径和新机器准备缺口。

MATLAB MCP 的连接探测返回 `failed to attach to MATLAB session`。因此新 `.m` 工具的 MATLAB 运行和 Code Analyzer 检查未执行，不能写成 MATLAB 测试通过。遵守 `AGENTS.md`，没有启动本机 MATLAB 进程或改用 SSH。

为使本次记录可直接阅读，初版用 Node 做一次纯文本打包：从同一 JSON 和 HTML 模板嵌入数据，并使用模板中的 Markdown 生成函数输出文字版。它不执行 `.m` 或任何物理计算，不是 MATLAB 运行方式的替代。该一次性脚本、完整日志、退出状态、来源/输出 SHA 与验收检查点留在 `results/validation/project_progress_dashboard/20261110_01/`，不作为第二套长期维护入口。MATLAB 会话恢复后应运行正式入口，并对照初版文本检查一致性，再更新 `runtime_note`。

浏览器工具仅暴露内置浏览器，安全策略拒绝本地 `file://` URL，并明确禁止绕过。未改用其他浏览器或代理访问本地文件。因此桌面/窄屏的实际显示、真实点击和键盘交互尚未目视验收；只执行离线脚本、结构、内容、链接和 CSS 响应式规则检查，不能把这些检查记成浏览器视觉验收通过。

已完成的初版检查：七条工作线、十条事件、29 个证据引用，3 个预期缺失入口均在指定 Git 历史对象中存在，其中 QPSK/共同窗入口的 SHA 与原检查点一致；页面实际脚本在隔离的内存 DOM 中执行正常，与 Markdown 的叙述内容一致，28 个页面本地链接加维护源链接有效。原生 `details` 和响应式 CSS 已实现，但实际浏览器交互/布局仍待验收。通用 MATLAB 静态模式扫描提示三处预期文件写入和一处用于 SHA-256 的内置 Java 调用，已人工检查用途及边界；这不是 MATLAB Code Analyzer 通过结果。

生成工具只有报告视图接口：`render_project_progress_dashboard(overrides)` 可覆盖绝对 `source_path`、`output_dir`、`markdown_path`；读取源/文字版限定在 `reports/`，HTML 限定在 `results/visualization/`。公共信道 API、生产参数、传播/通信公式和求解器均未修改，所以没有执行标量频率或宽带通信数值回归，也无需修改 `vertical_comm_guide.md` 的实现公式。
