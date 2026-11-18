# PM K0/K1 准备与本地测试回执

2026-11-17。**PREPARATION_PASS；正式K0未执行，K1禁用。** 测试只用已有本地MATLAB MCP R2025b Update6，不启动其他MATLAB进程。SSH连接、远端MATLAB、新声学求解、Git提交/推送全部0。没有生产PM、PS、BIE、FullFK、通信接口或原系数文件的修改。

## 最终交付与证据

- [实验方案](pm_wind_cutoff_k0_k1_plan_2026-11-17.md)、[用户命令](pm_wind_cutoff_manual_commands_2026-11-17.md)。
- 最终包：`results/experiments/pm_wind_cutoff/preparation_20261117_05/`，input_package.tar **70656 bytes**，SHA256 `ac8e3845271f8c04662f76786010bfcc38f59aebab03a67d166a35e131b27970`。
- manifest含15项：10个代码/工具、3个原系数CSV、配置、来源。`provenance/sources.json`固定23个已检查源版本及Git HEAD；工作区已有未提交改动，因此Git HEAD不是完整输入，执行以manifest SHA为准。
- 最终MATLAB预检：[final_preflight_20261117_05.json](../results/validation/pm_wind_cutoff_preparation/final_preflight_20261117_05.json)，manifest SHA `cf48069437260281d34a560ba6c3bcee1047d61acdc18154e8c1d68aeb6fd6c6`，45配置准入通过，仅preflight。
- [解包重定位/原版本SHA回执](../results/validation/pm_wind_cutoff_preparation/relocated_20261117_05/validation_receipt.json)：15项SHA、Bash语法、23个来源版本通过。最终05包的6个MATLAB文件与04包已数值测试文件逐项SHA一致，未无故重复数值测试；05只进一步完善非MATLAB资源基准及回收记录。

## MATLAB测试

入口 `scripts/validation/validate_pm_wind_cutoff_preparation.m`。主证据：[receipt.json](../results/validation/pm_wind_cutoff_preparation/20261117_test_04/receipt.json)，[完整MCP输出](../results/validation/pm_wind_cutoff_preparation/20261117_test_04/mcp_response.json)。成功测试耗时 **25.8637878 s**；这不是全部开发墙钟或正式45组K0耗时。

**287项准备断言通过，0失败。** 检查三个Seed、三个U、五个N的45种**系数配置**，256点仅用于快速数学/配对测试，没有将这45次调用封存为正式K0。原12模态恢复误差2.66454e-15 / 5.32907e-15 / 1.22125e-14；U6原前12项精确相同；跨U标准Gaussian相同、跨N低阶表精确相同；每个U的Sk与生产raw_pm公式差0，全局RNG状态不变。

独立保存数据集成测试仅 **3个TEST_ONLY案例**，均为seed260001：U5/N7、U6/N12、U8/N40。每组8192主采样/4096子网格；保存MAT/CSV、显式日志、SHA、检查点，独立Fourier/C2/导数、谱/DFT、理论积分及采样审计通过：

| TEST_ONLY条件 | 标量重建最大绝对差 | 六RMS最大相对采样差 | 最大坡度采样绝对差 |
|---|---:|---:|---:|
|U5/N7|2.42861e-15|4.82381e-7|6.93889e-18|
|U6/N12|7.96585e-15|9.66382e-8|3.46945e-16|
|U8/N40|2.74780e-14|9.01974e-8|1.11022e-16|

通过门：重建<1e-11；六RMS相对差<1e-4；最大坡度绝对差<1e-4。完整逐项回执在测试包output各audit.json。正式45组K0的其余42个保存/采样审计尚未运行，不能以这三例代替。

恢复测试包含：默认preflight、TEST_ONLY/正式模式隔离、K1参数拒绝、完整完成后SHA核验跳过、删除提交指针的**合成测试**仅补独立审计、通过案例跳过、生成日志SHA不变、已完成产物篡改拒绝、输入漂移拒绝、孤立generation_started标记拒绝自动重生。测试中移动的指针仅属于新建TEST_ONLY夹具，未触碰R1/R2或真实结果。

两次汇总绘图及FIG重开各通过27条统计曲线+12条海面/谱曲线；最终两张PNG均已目视查看，数据、单位、原始/C2线型及TEST_ONLY标签可读。三条件测试图用于软件验收，其U与N同时变化，**不能据此推断主效应或交互**。正式K0才会绘完整15条件×3Seed图。未声称Linux与Windows像素一致，也未增加PDF。

Code Analyzer检查7个MATLAB文件：6个0消息；renderer有1条`find`可提前返回的性能info，无错误/警告。[完整Analyzer回执](../results/validation/pm_wind_cutoff_preparation/analyzer_20261117_04.json)。该提示不影响小型冻结矩阵的正确性，没有为了消除提示再改最终已测源码。

本地MCP diary在初期测试为空，因此最终runner采用显式UTF-8日志持久化每组generation/audit和summary，完整本地stdout另保存在mcp_response.json。远端supervisor直接把真实子进程stdout/stderr写文件，保存其真实退出码。不能把本地MCP函数返回等同于某个新MATLAB进程的OS退出码。

## PowerShell / Python / Bash 离线测试

- PowerShell AST解析：打包、管理、mock测试工具通过；仅ShowCommands不联网。
- [最终传输mock](../results/validation/pm_wind_cutoff_preparation/transport_20261117_05/receipt.json)：7项通过，含Upload命令序列、手动Start/Status、Fetch逐文件SHA、错误tar SHA和越界task名拒绝；ssh/scp由本地函数mock替换，真实连接0。test.log保存全输出。
- [最终supervisor离线测试](../results/validation/pm_wind_cutoff_preparation/infrastructure_20261117_06/receipt.json)：13项通过，含manifest正常/篡改/越界拒绝、非零退出、退出0但缺检查点、RAM停止、墙钟停止、启动异常、全检查点闭合的纯元数据mock。Popen和killpg均mock，MATLAB启动0，不把45个纯元数据占位当作K0结果。
- Python标准库AST解析、最终15项manifest和解包重定位通过；`D:/Git/bin/bash.exe -n`验证启动器语法。Windows Git Bash的语法检查不等于Linux原生flock/nohup执行。

保留准备候选01–04及失败测试：01遇空数据块Java SHA接口异常，修复为空块跳过；02通过旧集成；03完善显式日志/采样/理论审计；04完成最终MATLAB测试；05完善恢复RAM最大基准及归档记录。早期传输mock发现PowerShell辅助函数`Ssh`与`ssh`大小写不敏感重名，发生函数递归，已改为Invoke-K0Ssh；当次测试被停止，未调用真实ssh.exe。另有mock作用域/Windows路径匹配修复，失败日志均保留，不回写成功。最终05是唯一交付启动版本。

## 尚未执行的检查与范围

以下均明确为 **UNRUN**：真实SSH认证/上传/回收；Ubuntu MATLAB R2025b批处理许可和跨平台实际执行；Linux原生flock/nohup/继承锁与真实进程组信号链；正式45组K0的全部保存数据/统计/图件审计；所有K1独立BIE、PS、完整FullFK及411点声学误差。用户手动启动将首先执行本包准入，缺依赖、版本、许可、RAM或SHA不通过即停止；不能把准备PASS升级为远端数值PASS。

原生产/接口/通信公式未变，未运行标量信道或宽带通信回归，无需修改vertical_comm_guide。脚本导航、任务表、PROJECT_CONTEXT及进度源更新为准备完成；下一步仅等待用户手动K0，再审查结果。没有自动监控、远程重试、K1、宽带、通信或Bellhop推进。
