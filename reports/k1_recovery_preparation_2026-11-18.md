# K1诊断执行与恢复04包准备

2026-11-18。**RECOVERY PACKAGE PREPARED，等待用户手动启动。** 完整研究矩阵仍未完成，原N40原生abort具体触发组件仍未定位。本轮只执行获批非声学诊断/保存图件重放、本地代码及测试，未运行新的PS/FullFK声场或BIE；未上传或启动恢复包、未改生产接口或提交推送Git。

## A阶段真实执行结果

|检查|结果|时间与证据|
|---|---|---|
|A1：独立MATLAB启动/基础环境/日志|PASS，R2025b Update1、单计算线程、基础许可证可用|进程16.765s；[A1回执](../results/validation/k1_recovery_preparation/20261118_01/A1_A2_retrieved/A1/receipt.json)|
|A2：N40输入只读＋合成数组|PASS；只读取40行系数/冻结配置，没有新表面生成或入射/散射场计算|进程12.601s；FFT往返4.09465e-16、256阶合成线性系统残差2.56287e-14；[A2回执](../results/validation/k1_recovery_preparation/20261118_01/A1_A2_retrieved/A2/receipt.json)|
|A0 Windows本地图件|四次尝试均FAIL，已停止该分支|原两次及用户追加两次；画布变化或帧像素尺寸不符，失败保留，无降门|
|用户追加的Ubuntu保存图件复核|PASS；12条曲线、画布/轴布局、两图各两次重开PNG像素差0|57.673s监督器墙钟、RSS峰3,561,242,624bytes；[监督器回执](../results/validation/k1_recovery_preparation/20261118_01/A0_Linux_retrieved/receipt.json)、[严格图件回执](../results/validation/k1_recovery_preparation/20261118_01/A0_Linux_retrieved/plot_receipt.json)|

A1/A2在独立目录`/home/misc/cwx/k1_recovery_diag_20261118_01`顺序执行，无parpool、无额外并发席位测试，真实退出均0；日志、身份、资源、prefs隔离，32项回收SHA通过。Ubuntu图件在另一个独立目录`/home/misc/cwx/k1_saved_graphics_diag_20261118_01`单次运行，读取旧任务已审计N7的common/comparison并核SHA；不重算N7，36项回收SHA通过。所有本轮创建的MATLAB进程已退出，不持续监控。

Windows四次候选分别尝试直接print、临时副本print、临时普通图窗默认设置、getframe原生帧捕获，未通过原画布或尺寸门。图窗默认设置在退出时恢复，没有把重置画布或缩放图片用于掩盖门失败。候选代码不进入恢复04包；旧本地FAIL继续保留。用户新增授权Ubuntu保存复核后，以**实际执行平台的原严格图件门**作为恢复包绘图资格；这不是Windows修复成功。正式Start仍先在Ubuntu重放N7保存图，失败则在任何新声学计算之前停止。

非声学诊断通过只能说明这些执行路径本次可用，不能证明N40崩溃已修复，不能归因到某个GTP/FFT/FullFK模块。没有更换MATLAB更新、系统库、线程策略或物理参数。

## 原资格、失败现场与新执行的独立链

最终目录：`results/comparisons/pm_wind_cutoff_k1/recovery_preparation_20261118_04/`。只有04为交付，01–03准备候选/失败保留。

- `lineage/original_inputs/`保存原05完整冻结输入，原296项manifest仍有效。
- `lineage/failed_task/`保存原回收575项证据及两次会话、N40失败标记与旧日志；原575项SHA审查前后保持。
- `cases/U6_N07/`逐字节复制旧N7完整封存链；原case manifest和签名不换，不在新代码manifest下伪装重新签名。`lineage/N7_import.json`及N7本地数值审查形成独立导入证明。控制器明确禁止派发N7 worker。
- 其余五面使用新代码/新case manifest，但input和policy逐字节等于原05；每例有父case/master指纹。N40为一次显式新尝试，原失败仍在lineage，不删除或重置旧调用标记。
- 原R1 U6/N12严格锚点继续按资格SHA复用，绝不把其密度导入新五面。风速仍称U_PM；Kmax、L160m、原42–50m C2、PM系数、数值采样、绝对相位v2、0.10m二维条带/4kHz和411真实坐标保持。

新增`run_k1_recovery_model_job`只为输入、direct、flat、PS1–6、FullFK1–4、common保存及原模型门添加立即关闭文件句柄的BEGIN/END记录；原PS局部函数逐字一致（只统一换行符），整体数学代码去除日志与for行展开后相同。MODELS保存审计事件绑定最终模型seal；MATLAB进程失败后按确切PID收集有限大小的本任务原生转储，公开文本去除会话/许可证标识。原生产PS/BIE/FullFK函数与数值门均未修改。

## 执行顺序、预算与恢复

用户手动Start后：完整双链预检 → 无声学MATLAB准备验收 → Ubuntu N7保存图件验收 → 跳过已封存N7 → N40一次MODELS及独立模型审计 → N40独立C0–C4 → U5/N12、U8/N12、U5/N40、U8/N40逐面原完整层级 → 七条件汇总。N40再次崩溃、模型门/参考门/图件门失败均立即停止后续派发，无自动重试、调参、降门或换主机。

|门|继承/实现|
|---|---|
|原累计时间|旧两会话636s＋A1/A2向上取整30s＋追加Ubuntu保存图复核58s＝**724s**；总10800s不变，启动前剩10076s|
|RAM基准|至少119021486080bytes，和诊断及未来启动可用量取最大，恢复继续取历次最大；预留30%，每次派发扣预留后仍需12GiB|
|并发/RSS|1个MATLAB进程×1计算线程，OMP/MKL/OpenBLAS=1；活动会话RSS≤12GiB，无parpool|
|单案例|1800s跨会话累计；N40继承18.833152669s，后续MODELS/BIE/审计继续累加|
|N40 MODELS|进入该阶段后180s内完成；失败保留唯一未封存调用，不自动再次尝试|
|汇总|≤600s，同时服从累计10800s|
|压力/磁盘/CPU|原memory/io PSI准入5%/10%、运行10%/30%连续10样本、load1≤21、5GiB/5000 inode，保持每秒资源采样及原身份核验|
|BIE计数|本次恢复最多25次新BIE，另有原N7的5次，最终累计30；原R1锚点单列。计数以实际durable invocation marker为准|

仍按55–75min作余五面含审计/汇总的条件规划；当前只有N7成本实测，不保证N40耗时。未来资源以Start时重新检查为准；不干扰D1或其他任务，不自动排队。

已封存阶段按SHA跳过；已求解且有有效solved seal只补审计；最终seal已在而顶层状态遗漏则仅修复元数据。未封存调用/孤立数据/不确定fork拒绝恢复，不能通过直接再Start重算。已全部完成时仅复核并退出，不再次启动MATLAB。配置/会话预算账时不匹配时拒绝，并保持原账本不被724s默认值覆盖。

## 本地验证与边界

|检查|结果|
|---|---|
|实际04包MATLAB准备|37项，9.951s，8个声学入口哨兵、五例新dry-run及原24状态/政策检查，0新求解；[回执](../results/validation/k1_recovery_preparation/20261118_01/preparation_04/receipt.json)|
|离线控制器|25项，Popen/信号均mock：N7禁止派发、5例顺序、失败停派发、180s MODELS、RAM/RSS/PSI/累计/单例停止、预算继承、旧账不覆盖、完成零重复、已封存跳过/仅元数据修复和篡改拒绝；[回执](../results/validation/k1_recovery_preparation/20261118_01/infrastructure_03/receipt.json)|
|PowerShell传输|8项mock，DryRun/Upload/Start/Status/Fetch和内外SHA/错误任务名/错误归档SHA拒绝，无真实恢复上传；[回执](../results/validation/k1_recovery_preparation/20261118_01/transport_01/receipt.json)|
|汇总新来源列|TEST_ONLY原R1数组＋合成seal，真实模式拒绝夹具，14行指标、5754逐点、2行N7导入/2行R1锚点/10行本恢复来源通过；绘图器明确mock，无新图件视觉资格或物理结果；[回执](../results/validation/k1_recovery_preparation/20261118_01/summary_metadata_01/receipt.json)|
|依赖|30项MATLAB依赖全部包内同SHA，所需产品仅MATLAB，无新工具箱|
|静态|PowerShell AST、Bash语法、6个Python文件编译通过；9个MATLAB文件Code Analyzer有2条非致命格式信息（ALIGN/NOCOMMA），无语法错误；[记录](../results/validation/k1_recovery_preparation/20261118_01/preparation_04/analyzer.json)|
|归档重定位|1459项输入/证据SHA、解包后Python预检和Git Bash dry-run通过，0网络/求解；[回执](../results/validation/k1_recovery_preparation/20261118_01/relocation_01/receipt.json)|

初版准备生成器遇到替换次数/UTF8读取错误、01包换行差的重复检查未通过、首版offline mock把已退出进程误标存活导致“完成后重启”测试失败，均已修正并保留旧目录。修正只影响准备工具/测试fixture与记录，不改数值门。37项实际包测试、25项controller和8项传输最终通过。Windows失败代码仍为诊断候选，不用于正式恢复，未宣称修复成功。

归档`input_package.tar`为99,563,520bytes，SHA256：

`77bce9c2b26806ea796758bcdb1662359ecf1f560d7dfa1f3cf3282223dc3f6e`

新恢复包没有实际上传或Start；Linux正式恢复CLI锁/信号复用原已做原生11项诊断的身份安全函数及同一K1租约，本新任务端到端正式执行仍UNRUN。七条件完整结果和25次新BIE仍缺，不从准备PASS升级科学资格。

新增仅研究编排、双链导入、事件日志、诊断/保存图件与准备测试；生产API/物理模型未变，无通信回归或vertical_comm_guide实现公式修改必要。相关入口已更新scripts/README，任务表/PROJECT_CONTEXT/进度只更新受影响K1范围。完成后停止，等待用户使用[手动命令](k1_recovery_manual_commands_2026-11-18.md)。
