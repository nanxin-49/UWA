# K1执行包准备验证

2026-11-18。**PREPARATION_PASS；正式K1未运行。** 最终交付 `results/comparisons/pm_wind_cutoff_k1/preparation_20261118_05/`。296项输入/代码manifest，归档15,535,616 bytes，SHA256 `ddc8abe83a30619852b4f0d6ab0a582dc5b2789ffbf4d29946092ddac57ecf60`；主manifest SHA `815e44af687af278bcd45b5a05ee7be60b638aaf0c97901e2ac117dd082b11ec`。

## 已完成

- **12:21原生补测**：[Ubuntu原生11项检查](../results/validation/k1_preparation/native_20261118_07/retrieved/receipt.json)15.646s、0失败；冻结控制器/预检SHA与05包一致。锁继承/竞争/释放、真实RSS、错误身份拒绝、TERM及15s后KILL、实际CLI无租约拒绝均通过。仅四个有界轻量诊断进程且全部退出；MATLAB启动/正式声学/其他进程信号均0。[传输回执](../results/validation/k1_preparation/native_20261118_07/transport_receipt.json)核35项SHA，归档SHA `c826acdc8f15146e766f7dc0dd11069efdbe30e7181f31c3e3d26de852ea77a6`。独立诊断锁不是正式共享K1锁，未调用正式Start。
- 同次可用RAM38.50–39.02GiB、memory full PSI为0%、io为0.10–0.33%，只读准入当时PASS；Swap近耗尽及并发席位未知，仍采用A。05包296项清单再次通过，归档SHA未变；既有MATLAB/汇总/静态/重定位测试收据SHA复核不变，复用其证据，无重复求解。

- 实际资源快照两次；一次有界非声学MATLAB诊断32.918s、退出0，无pool/声学。选择A：1进程×1线程，原因与预算详见[资源决策](k1_resources_strategy_2026-11-18.md)。并发许可证总上限未知，没有假装验证额外多席。
- 打包前原R1回收287项和K0回收574项SHA核验；严格历史锚点匹配，未改原件。源参数、PM系数、C2、相位和411坐标全部匹配，资格SHA及原最终封存绑定见包内anchor/anchor_proof.json。
- [最终MATLAB准备测试](../results/validation/k1_preparation/matlab_20261118_05/receipt.json)：24项检查、7个声学入口哨兵、六例dry-run、锚点重审、输入漂移/路径/孤立数据/重复调用拒绝、已封存跳过/仅审计恢复与v2原门检查；3.1414294s，0新声学。
- PS局部数学函数及PS/FullFK模型区域与原R2逐字一致；既有传播/BIE/FullFK和独立审计依赖保持原实现。K1仅新增研究适配与编排，不修改生产模型或数值门。
- [最终汇总管线测试](../results/validation/k1_preparation/summary_pipeline_20261118_05/receipt.json)：明确TEST_ONLY的历史R1数组和合成检查点，14行指标、5754逐点、七条件交互/谱表及56条FIG曲线/零像素重导通过；34.2882108s。真实模式拒绝测试包；没有六个新物理声场结果。
- [绘图夹具](../results/validation/k1_preparation/plotting_20261118_02/receipt.json)：单例12条、汇总56条曲线及8次零像素重导，未求解；输出已目视检查，保留PNG/FIG，无PDF。
- [最终基础设施离线测试](../results/validation/k1_preparation/infrastructure_20261118_05/receipt.json)：32项，包括完整派发与汇总、失败/缺审计停止后续派发、完成后零重复、封存后遗漏状态的仅状态恢复、RAM/RSS/压力/累计预算停止、PID重用与身份变化拒绝、未确定fork拒绝、SHA/越界/并发上限篡改和归档链接/越界拒绝。Popen、信号均mock，真实MATLAB/SSH/进程终止0。
- [最终PowerShell传输mock](../results/validation/k1_preparation/transport_20261118_05/receipt.json)：8项，DryRun/ShowCommands零网络，Upload/Start/Status/Fetch命令及内外SHA/保存外部指纹/错误SHA/错误任务名拒绝通过，无真实连接。
- [Code Analyzer](../results/validation/k1_preparation/analyzer_20261118_05.json)：11个MATLAB文件均0消息。PowerShell AST、Python编译和Bash语法检查通过。

## 保留的失败与边界

初始沙箱SSH别名解析失败、首版探针重复关键字语法错误、首个绘图夹具结构体初始化错误、首个controller测试的killpg mock符号错误均保留；已修复并以新目录复验，不回写旧证据。01–04为准备候选，只有05交付。最初测试的空分项R1 C0来自导入计时，历史汇总02已明确其计时边界，不伪填分项。

初次准备时原生锁/信号未实跑、10:41快照未过准入；12:21已用独立非声学诊断补测上述原生分支且当时资源准入通过，旧回执保留。完整正式启动器Start和K1仍未运行，原生/离线检查不能替代声学PASS。用户未来Start重新检查；没有后台K1、持续监控、生产修改、宽带/通信/Bellhop或Git提交推送。

本次Python字节码写入被本地目录权限拒绝，改用内存compile检查，不改变权限。native_06传输驱动因换行转义在远端解析前失败，未创建远端目录或启动测试；修正后使用独立native_07通过。失败日志保留，冻结05包不受影响。

新代码为研究适配/独立审计/报告/基础设施；公共信道与通信接口不变，无需标量/宽带通信回归或修改vertical_comm_guide。[用户命令](k1_manual_commands_2026-11-18.md)。

最终解包重定位及无求解Bash dry-run通过，25项静态依赖均已包含或与包内同名文件SHA一致，所需产品仅MATLAB。Windows Git Bash初次缺python3名称，使用仅该进程内的python3函数映射到已有Python复验；不改远端python3/系统环境。回执在results/validation/k1_preparation/relocated_20261118_05/receipt.json。
