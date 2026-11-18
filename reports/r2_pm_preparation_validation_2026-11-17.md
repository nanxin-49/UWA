# R2执行准备：本地验证报告

2026-11-17。状态：**PREPARATION_COMPLETE_WITH_RUNTIME_CHECKS_UNRUN，停止等待用户手动启动**。固定0.10 m二维均匀条带、4 kHz、原PM谱/窗口、411真实接收点；仅新增seed260002/260003。用户随后允许两Seed并行，因此交付支持两进程、每Seed单线程与共享30% RAM预留。本轮实际SSH连接0、SCP进程0、MATLAB执行0、正式R2声学求解0，不提交/推送Git。

方案见 [技术方案](r2_pm_multiseed_plan_2026-11-17.md)，命令见 [手动操作说明](r2_pm_manual_commands_2026-11-17.md)。最终执行包只使用 **03版**：

| Seed | 包目录（相对仓库） | 归档字节 | input_package.tar SHA-256 |
|---|---|---:|---|
| 260002 | `results/comparisons/piston_2d_pm/r2_preparation_20261117_03/seed260002/` | 300032 | `6262cedc24cf2babf71c24dd9e1443441b055cbb664126d9cb1257b114451caa` |
| 260003 | `results/comparisons/piston_2d_pm/r2_preparation_20261117_03/seed260003/` | 300032 | `b79b5e9eb347eec58ef41edb0dd41300310f851efb2786d41a1832ce924b8c7a` |

每包manifest绑定44项代码/输入/来源/政策，最终package.sha256绑定49项（含本地检查摘要、manifest与操作说明）。两个包输入互不混用；无common.mat或已求解密度。单Seed正式完成后才会产生自己的结果、检查点与封存链。

## 实际通过的检查

证据主目录：`results/validation/r2_pm_preparation/20261117_03/`。完整本地stdout/stderr与实际子进程退出码文件保存在该目录；PS脚本检查另有JSON回执。各检查的含义限于下表，不能当作声学PASS。

| 检查 | 结果 | 证据 |
|---|---|---|
| 两Seed冻结输入/SHA/拒绝测试 | 各28项通过，0失败；每包44项manifest SHA通过 | `seed260002/receipt.json`、`seed260003/receipt.json`，对应stdout/stderr/exit_code |
| PM来源 | 各Seed系数与明确父文件SHA相同；k/Sk与R1相同，无重新抽样/归一化 | 包内`provenance/sources.json`、`input/*.csv` |
| 接收坐标 | 全411点与原R1点表复制一致；x=(-205:205)×.1953125 m、z=3 m | preflight及28项检查 |
| 几何算术检查 | 260002：深度−.3235664… .3033853 m，最大绝对坡度.0898530；260003：−.4306613… .3817594 m，坡度.1323574 | 各Seed receipt的geometry_only；0.01 m采样，只是几何检查，不是声学或连续极值证明 |
| 原PS公式 | R2 local_ps与R1原函数文本完全相同 | 各Seed28项检查 |
| 复用数值核 | 对存在的R1冻结核逐字节SHA核对；FullFK及入射helper对原R1 models快照核对 | 包内`provenance/sources.json`的code记录 |
| PowerShell/Bash语法 | 2个PS文件和6个Bash文件解析通过 | `syntax_receipt.json` |
| 单Seed状态展示 | 12种本地合成状态通过；包括初始化/运行/失败/真实退出码/陈旧PID等 | `status_stdout.log`及逐状态log |
| Windows原生stdin传输 | 7项通过；引号、括号、中文、大于命令行长度输入、真实退出/错误输出保持 | `native_transport/receipt.json` |
| 管理器动作契约 | 8项本地mock通过；ShowCommands无网络，Upload不启动，Start/Continue只发本Seed启动请求 | `mock_manager/receipt.json`；ssh/scp为本地替身，无实际网络或远程进程 |
| RAM并行准入 | 10项实际Bash算术/非法输入/临界字节检查通过；32GiB拒绝、40GiB同基准通过、40/64GiB旧高基准拒绝 | `parallel/receipt.json` |
| 原runtime序列化 | 29项回归通过，initial/starting/max与非法JSON/缺文件/不可覆盖检查保留 | `../20261117_02/runtime_stdout.log`、`runtime_exit_code.txt`及完整错误预期日志；核源码未变，未重复测试 |
| 最终归档重定位 | 两包分别解包到新目录，完整49项SHA检查和shell --dry-run通过，退出0；无jobs/sessions/control/model_stage/finalization新状态 | `sealed_archives_receipt.json`、`sealed_dryrun_<seed>_*` |
| R1不可变性 | 精确选定回收目录的287项SHA全部符合原清单，0失败，0.6428512 s；未加载MAT、未复算声场 | `r1_integrity_receipt.json`、stdout/stderr/exit_code |

28项输入检查覆盖：源/频率/配方/坐标/PM谱固定、无密度导入、PS原函数同文、输入重定位、缺失/篡改系数与坐标/配置/政策/求解器/门/资源helper拒绝、非法Seed/线程/重试/资源比例/并发数量拒绝、重复与越界manifest路径拒绝，以及即使重算SHA也不能放宽敏感性门或换C4窗口。合成拒绝产物位于明确标注的NONPHYSICAL_TAMPER_FIXTURE，不能当正式结果。

## 未运行项和启动前硬边界

现有MATLAB MCP身份查询返回`failed to attach to MATLAB session`。没有启动本机matlab.exe、没有切换SSH或Octave；**MATLAB dry-run、MATLAB合成恢复/政策测试、Code Analyzer、新R2数值门/独立审计/图件及三Seed汇总全未运行**。因此本轮不能声称R2参考合格、恢复状态机已由MATLAB验证，或新MATLAB编排器已实跑。成本未在新Seed上实测。

Windows Git Bash没有`flock`，本地原生Linux槽锁四项未运行，不能把14项全部标PASS。原样保留02版首次尝试的非零退出与错误日志；03版回执明确只通过10项可用算术检查。已将原生锁测试设为用户手动启动后、MATLAB之前的必需检查；两Seed同时持槽、同Seed重复拒绝、释放后再取得槽须在实际Linux通过。

正式Linux启动会在单线程MATLAB进程中先调用`validate_r2_pm_preparation`：七个声学入口哨兵、SHA/几何dry-run、fresh/orphan/audit-only/skip/前序与产物篡改拒绝、v2非单调诊断和保留原1e-4门的合成测试，以及Code Analyzer。错误会阻止进入`run_r2_pm_seed(...,'run',...)`。该安排保证用户可以手动开始验收流程，但不将未来运行检查写成当前已通过。

实际远端当前内存/版本/依赖、双进程资源曲线、进程信号/恢复链、C0–C4解及收回后的科学审查都待用户启动后获取。12GiB/Seed是继承的准入估计；24GiB合计和共享30%门是保守运行条件，不是实测峰值保证。每Seed模型误差是否分辨、FullFK是否仍较小、相位项/峰位/能量分区如何变化，全部待结果；不预填R1数值为新Seed结论。

## 保留的准备失败与修订

首次Prepare缺少R1 C4包中未携带的旧`audit_r1_bie_manual_final.m`，在生成manifest/封存前停止；原部分目录保存在`../20261117_01/failed_prepare_seed260002_01/`，错误说明在同级`first_prepare_failure.txt`。修复只涉及准备来源：从原R1对应研究源码取审阅副本，不修改R1。01/02候选以及02版缺flock失败都保留；最终03另建新目录冻结和验收，不覆盖旧失败，也不把候选当交付。

## 变更和后续

新增R2研究编排、冻结输入适配、最终独立指标复核、每Seed保存图件/三Seed汇总、包准备和手动管理脚本；复用原源、完整FullFK、BIE及独立审计实现。数值公式与门未变，外部公共声道/通信API、生产默认与R1原始结果均未改。此任务是研究执行准备，未做宽带或通信回归，无需修改`vertical_comm_guide.md`；已维护脚本导航和研究状态记录。每Seed实际运行必须保留完整C0–C4、源/Green/FullFK/MAT-CSV/FIG审计，5倍经验分辨规则不变。

下一步仅为用户按手动命令上传/启动03包（可并行），然后按真实检查点和日志审查。**本轮到此停止，不自动启动、不持续监控、不推进R3或更广范围。**
