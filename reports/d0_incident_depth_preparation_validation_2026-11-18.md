# D0执行包本地验证与dry-run报告

2026-11-18。状态 **PREPARATION PASS（限定软件准备范围）**。正式声场、深海角谱、粗糙边界梯度及采样收敛均未运行，不是声学PASS。SSH、远端MATLAB、新建本地MATLAB进程、PS/BIE/FullFK求解和Git提交推送均为0。

最终包：[preparation_20261118_06](../results/experiments/d0_incident_depth/preparation_20261118_06/)，28项manifest；归档SHA256 `a63dbbbd8f092fd8d85c04ceca8e3f8c65a23d21a6007f07aa5a9515cfd8c3dd`。配置与冻结数值门见[技术方案](d0_incident_depth_plan_2026-11-18.md)，手动操作见[命令](d0_incident_depth_manual_commands_2026-11-18.md)。01–05包、准备失败目录及TEST_ONLY夹具全部保留，不上传。

## 已运行检查

| 检查 | 结果 | 独立证据 |
|---|---|---|
|最终包SHA/固定配置/K0单组检查点引用|28项通过，5深度固定，原源/相位/P0参考代码原样复制|包内manifest/provenance；最终seal回执|
|Python离线基础设施|15项通过；解包重定位、篡改/越界/cash引用/链接拒绝、RAM和时限边界、失败包归档/回收、放宽门拒绝|`results/validation/d0_preparation/infrastructure_20261118_05/receipt.json`|
|PowerShell传输mock|10项通过；DryRun零传输、UploadAndStart仅一次启动、Status只读、退出75不重试、Fetch完整性及失败码保留、任务名注入拒绝|`results/validation/d0_preparation/transport_20261118_04/receipt.json`及scripts.json；所有SSH/scp均替换为TEST_ONLY假实现|
|现有MATLAB MCP准备测试|R2025b，既有PID22340，56.3737486 s，57项计数断言，正式模式拒绝，field_calls=0|`results/validation/d0_preparation/matlab_20261118_07/receipt.json`及matlab_diary.log|
|五个正式MATLAB文件Code Analyzer|0消息|同目录code_analyzer.json；准备测试入口另行Analyzer|
|合成FFT及归一化|Parseval最大3.23692e-14，逆FFT相对差最大4.05942e-16；复中心及正实RMS归一化/中心相位/局部方向检查通过|同MATLAB测试回执及TEST_ONLY保存数据|
|六类图件与可编辑FIG|6组FIG/PNG/重开PNG、60条曲线与轴标签一致；检查了代表PNG布局|同TEST_ONLY_SYNTHETIC/figures/；不是D0研究图|
|PowerShell解析、Bash语法、Python解析；最终DryRun；来源SHA保持|按最终seal的static_checks及sources_unchanged回执记录|`results/validation/d0_preparation/final_20261118_02/`|

MATLAB共享会话起始cwd在另一个工作树；未改变cwd，以绝对路径加入所需代码，测试后恢复路径，未清理用户工作区。四个声学调用哨兵覆盖活塞、原P0角谱、新分段角谱及R1入口。预检只读JSON/CSV/SHA，正式模式在本地、非监督器环境被拒绝。

合成测试特意使用25 Hz、c1000、rho1000、宽0.02 m、虚拟深度2/4/6/8/10 m和解析平面波夹具；不调用实际活塞或深海独立角谱。所有合成图均写TEST_ONLY，不混入正式包或研究结果。图件目视检查发现恒定场自动坐标会夸大舍入噪声，最终06版为幅度设置零起点、角度共用±90°，对其他常数曲线设置最小显示范围；未改物理数组。

## 准备失败与修正记录

- 第一轮哨兵断言因Windows路径分隔符比较失败，未调用声场；第二轮因会话当前目录的公共根入口优先于path而失败。改为路径标准化、只对可正确遮蔽的四个声学入口安装哨兵，并独立静态确认D0没有公共信道入口依赖。原 `matlab_20261118_01/02` 保留。
- 第一轮传输mock的高级PowerShell函数将`scp -o`误判为公共参数缩写，测试停止；改成纯参数透传的测试替身。原 `transport_20261118_01/` 保留；没有实际SSH/scp。
- 01包renderer曾有一条END对齐警告；后续修正，最终五个正式文件Analyzer为0。早期候选与各成功测试版本均保留；最终包按最新代码重新冻结。
- 最终静态封存首次误用不存在的默认Git Bash路径而停止；保留final_20261118_01，在实际Git安装D:/Git中定位bash后重新保存final_20261118_02。不把先前条件跳过的Bash检查算通过。准备测试入口的纯缩进警告另行修正。
- 看板同步首次调用被共享会话另一工作树的当前目录入口遮蔽，误刷新 `C:/Users/ASUS/.codex/worktrees/bellhop_rotated_fig12b/Explain/reports/project_progress.md` 及该树 `results/visualization/project_progress/index.html`、`render_receipt.json` 三份生成产物（16条引用警告）。该树源JSON和生产代码未改；没有刷新前快照，未用Git覆盖回退其潜在未提交产物。此误刷新不作为本任务证据。随后通过MCP明确project_path到本仓库并断言renderer绝对路径；修正PowerShell新JSON对象的字段顺序后，本仓库正式刷新92引用、0意外缺失。初次进度JSON内存构建的不存在字段错误与首次MATLAB schema拒绝均未作为通过。

这些是准备工具/夹具问题，不是D0声学失败；不删除或改写旧失败记录，也不据此放宽任何数值门。

## 尚未运行，正式任务必须满足

1. Ubuntu原生MATLAB许可证、真实后台detach/flock/进程组信号和资源采样。Bash语法及Windows mock不构成这些运行行为的验证。
2. 五深度GL64/128、分段角谱两阶、原P0被积函数适配、真实梯度/粗糙入水法向、phase v2及MAT/CSV检查。
3. Tx100/1490 m的0.025 m细评价、Hann FFT采样、相邻相位门和真实平面波误差分解。
4. 正式六类PNG的可读性与物理解释、本地回收后独立保存数据复核及最终研究结论。

脚本已将这些数值门置于远端顺序流程，失败停止并保留checkpoint，不自动改参数或重试；本轮不会为了测试它们越过禁止正式声场的边界。K0选定面仅原±50 m支撑，D1粗糙面覆盖缺口、R1/R2结论不可直接外推等限制保持。

最终静态回执另确认34项已读来源SHA保持、归档1772544 bytes、正式output不存在。MCP看板正确刷新回执位于 `results/visualization/project_progress/render_receipt.json`，仅验证文件引用，不算声学审计。

新增代码只涉及D0研究入口、研究诊断、独立审计、保存绘图、执行基础设施及准备测试；生产源码、接口和相位公式未变。已同步脚本导航、任务表及项目状态；无需通信回归或改写vertical_comm_guide生产方法。下一步由用户手动启动；Codex停止，不持续监控。
