# AGENTS配置与项目导航同步（2026-11-14）

范围：仅修复配置与当前导航。依据 `AGENTS.md` 的10 cm条带、阶段授权、收尾记录、执行模式和Git规则；不运行MATLAB、不连接SSH、不启动实验、不提交或暂存Git，不读取 `cash/`。

## 修复内容

- `PROJECT_CONTEXT.md`：在顶部追加当前决策，明确后续总宽 W=0.10 m、半宽 a=0.05 m；以新条目取代历史“20 cm主源”及旧三频优先顺序，保留全部旧记录。
- `README.md`、`scripts/README.md`：补充当前任务表、源宽和SSH Skill触发条件；根README的二维执行导航改指当前任务表，不恢复已删除的旧计划。
- `reports/project_progress.json`：同步当前摘要、二维工作线和接续步骤，追加中文配置事件；旧事件、结果日期、其他工作线保持不变。仅本次核过的二维导航日期更新，历史数值审计未重跑。
- `reports/project_progress_maintenance.md`：统一任务表与看板分工，明确用户禁止MATLAB时的同源文本打包、回执及未运行检查。`reports/project_progress.md` 和离线HTML由JSON/既有模板生成，不手改。
- `reports/active_task_tracker.md`：R1保留 `NEXT - NOT AUTHORIZED`，仅补配置同步链接并更新最后行动；其他研究行与数值状态不变。

`AGENTS.md`、SSH Skill及 `.gitignore` 无需修改。Skill已在当前会话目录中被识别，本次仅静态核对；远程操作授权后才启用。新增Markdown未被Git忽略；本次不暂存、提交或修改忽略规则。

## 当前研究顺序与限制

1. 候选R1：10 cm、4 kHz、冻结PM seed260001配对PE Model-0/BIE/完整FullFK比较；尚未授权计算范围、预算和参考/资源门。
2. R2额外种子仅在R1有可信参考/成本证据及明确授权后考虑；R3再作反射模型决策。
3. W1相干4–8 kHz频响/CIR、C1受控通信须分别通过前置门并取得范围授权。
4. Bellhop保持独立并行及未授权状态；M1、正式S5/S6的停止和预算不变。

10 cm为理想二维无限条带的名义总宽，不证明实物三维有效口径。20 cm及Gaussian旧结果不改写、不重标。旧P0/P1的有限10 cm控制仅按原证据范围复用；验证专用连续角谱PE不证明生产WAPE活塞接入。本次无接口、参数默认值、传播/相位公式或通信方法变更，因此不修改 `vertical_comm_guide.md`，不运行数值回归。

## 静态验证与证据

独立目录：[20261114_agents_alignment_01](../results/validation/project_progress_dashboard/20261114_agents_alignment_01/)。

该目录保存修改前配置版本/SHA和进度源、一次性文本同步脚本、stdout/stderr、实际退出状态、生成回执及检查点。检查覆盖JSON结构、事件唯一性/日期顺序、旧事件及其他工作线未变、新导航链接、研究状态/历史结果日期保留、JSON嵌入HTML与Markdown文本一致性、模板脚本的内存DOM执行及输入/输出SHA；不复查历史数值证据，不把文件存在当数值PASS。

权威静态结果见 [checkpoint.json](../results/validation/project_progress_dashboard/20261114_agents_alignment_01/checkpoint.json) 和 [文本生成回执](../results/validation/project_progress_dashboard/20261114_agents_alignment_01/render_receipt.json)。引用缺失按回执逐项保留与报告，不自动恢复、不读取归档、不将其记为数值通过。旧运行回执不覆盖。

未运行：MATLAB正式渲染、Code Analyzer、历史MAT/CSV/FIG数值审计、浏览器视觉/真实交互检查。本次未探测MCP，不能推断当前连接是否可用。静态配置同步完成不授予R1或任何后续阶段运行许可。

## 后续Git发布授权（2026-11-14）

用户随后明确要求提交并推送。上文“不暂存/提交”描述此前配置维护步骤；本次发布仅包含AGENTS、任务表、SSH Skill及相关配置导航修复。共享文档按内容拆分暂存，之前未提交的自由场/形状图示记录、脚本与旧计划删除保留在工作区。已提交历史事件、结果日期及数值状态不变，不授予MATLAB、实验或实验室SSH操作许可。

此前本地文本同步回执为8项静态检查通过、60个引用无意外缺失、10项输出SHA复核一致；这些数字对应当时完整工作区视图，未作数值资格升级。回执目录及离线HTML为本地生成物，不随本次提交；提交中的进度文字版仅从本次暂存的JSON与既有模板生成。原检查点保持不可变，追加的Git授权记录不追写其历史SHA。
