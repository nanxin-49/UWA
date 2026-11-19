# D1 准备验证回执

2026-11-18。状态：**PREPARATION_COMPLETE / MATLAB_NUMERICS_UNRUN / WAITING_USER_AUTHORIZATION**。没有SSH连接、MATLAB进程、声学求解、Git提交；没有改动R1/R2/K0/D0或生产代码。

最终交付：[preparation_20261118_03](../results/comparisons/d1_depth/preparation_20261118_03/)，31项代码/输入SHA，tar 215040 bytes，SHA256：

`dcb3197cab1717ba27528b9b1c65b20a4d7bf36fb4993efcac150b962f7dd42e`

01/02候选保留为准备过程，不用于启动。02增加了源/几何声学前置门、相位RMS/分解审计、深度变化经验分辨和更完整检查点；03进一步固定/恢复图件画布、轴和导出尺寸以规避已知openfig缩图风险，真实图件审计仍未运行。最终包只含代码、原系数、接收坐标、原policy和来源文本；不含旧声学场或任何散射密度。

## 已执行

|检查|结果|证据|
|---|---|---|
|最终包31项SHA/冻结配置/原系数十进制token/源/411坐标|通过|包内manifest.json，d1_package.py预检|
|离线基础设施|26项通过、0失败|[receipt.json](../results/validation/d1_preparation/20261118_01/package03_infrastructure/receipt.json)|
|传输mock|10项通过、实际网络0|[receipt.json](../results/validation/d1_preparation/20261118_01/package03_transport/receipt.json)|
|PowerShell AST|manager/test脚本解析无错误|本轮PowerShell Parser实际调用|
|Bash语法|Git Bash `-n`退出0|本轮D:/Git/bin/bash.exe静态调用；未执行启动器|
|解析几何/稠密成本算术|完成；不计算声场|geometry.csv、dense_cost.csv、checked_sources.json|
|MATLAB静态风险扫描|完成19文件，结果为需人工审查，不是PASS|[扫描JSON](../results/validation/d1_preparation/20261118_01/package03_static_scan.json)|

26项涵盖：归档安全重定位；代码篡改；即便重算manifest仍拒绝降门、第三深度、Rx3、系数漂移；越界/禁止归档路径；RAM/总时间/单深度/进程组RSS的等号与越界；tar符号链接/遍历/重复成员；失败部分结果的SHA和真实退出保留；部分输出不能PASS；Python AST；PS距离/载波、独立Green的Rx6、深源分段求积的静态确认。

10项mock涵盖：DryRun零传输、UploadAndStart明确三次脚本/一次copy、拒绝覆盖远端目录、解包前SHA、单次手动Start、Status命令、失败不重试、Fetch原回执/退出码保留、不安全TaskName在传输前拒绝。mock只在新测试目录替换通信函数，不修改正式包，不调用真实ssh/scp。

静态扫描初次以错误的相对路径调用，退出2，日志保留；改为root内`code`后完成，扫描器按高风险规则退出1：14处MAT读取、3处Java SHA、31处文件写入和15处序列化写入。逐类已审查：MAT仅加载本任务新生成的common/comparison/data，原输入没有MAT；Java限SHA-256；写入限独立任务输出。没有动态eval、shell启动或新增工具箱；**此扫描不是MATLAB语法、Code Analyzer或运行验证**。路径调用错误没有影响正式包或产生MATLAB进程。

## 未执行，不能宣称通过

- MATLAB解析/Code Analyzer、MATLAB级dry-run；用户禁止启动MATLAB，本轮也不探测MCP。
- 深度声源/梯度、PS六层、FullFK四层、十次BIE、线性/边界/窗口/经验分辨及独立保存数据审计。
- MAT/CSV真实声学数值、FIG反序列化曲线和PNG像素重导；这些都由正式任务生成并在退出前检查。
- D1实际Linux的flock继承/身份/信号/组RSS/PSI、许可证和R2025b工作负载；本地Bash语法通过不能替代原生测试。
- 实时可用RAM/磁盘/主机压力。历史K1快照不作D1当前准入证据。
- D1b任何新海面、求解或大规模试验。

正式手动Start将先重新做文件/配置/SHA及资源准入，再在单线程R2025b里执行MATLAB输入dry-run和正式入口。MATLAB语法或数值门失败即保存stdout/stderr/diary和失败检查点并停止，无自动重试。本轮交付的是**可审阅、离线工具检查完成的实验包**，不是已运行合格的声学系统。

## 修改与接口

新增D1专用比较入口、模型适配、门函数、独立审计、图件及基础设施；`scripts/README.md`加入导航，`.gitignore`放行正式工具并排除生成快照/篡改夹具。原源函数、BIE/FullFK、R1/R2/K0/D0文件均只读。生产信道API、输出字段、DSP及通信接口无改动；不需要通信模拟或`vertical_comm_guide.md`公式更新。当前研究导航与进度只登记D1准备状态，数值证据仍待用户手动启动和回收。
