# D1 CSV精度修复、保存C0审计与手动恢复包

2026-11-18。**修复及保存C0独立审计完成，恢复包准备通过，尚未上传或启动恢复作业。** 用户授权修复、审计并交付手动启动包；本轮通过既有MATLAB MCP会话审计保存数据，没有启动新的本地/远端MATLAB进程，没有新增BIE求解，没有Git提交。原失败目录、03首次执行包及R1/R2/K0/D0保持不变。

## 1. 恢复来源和原失败现场

原D1用户手动会话`20261118T090204_1101719`，Tx100m C0已求解，停在CSV/MAT坐标严格相等断言，未进入后续BIE和1490m。依据[首次诊断](d1_first_run_diagnosis_2026-11-18.md)，本轮用原03包Fetch工具回收至：

`results/comparisons/d1_depth/d1_failed_03_retrieved_20261118_01/`

外层归档SHA与远端文本回执核验，59项回收SHA在审计前后均通过；远端只新增本任务retrieval归档，未覆盖原输入/结果。原输入manifest SHA为`74b1a47e6bd59ece8589510b7287a3b3de00dc6bf3774142f99f2081be19d81a`。完整传输回执位于同名`_transport/`目录。恢复包包含这份原证据的冻结副本，运行不依赖旧远端目录再做任何写入。

## 2. 修复范围和接口

- 新增验证专用`write_d1_csv(table, absolute_file)`：所有double按`%.17g`写出，支持方法字符串、逻辑掩码和相位NaN；拒绝覆盖已有CSV。
- D1 BIE逐点CSV、模型误差/指标CSV、最终汇总及独立分辨CSV统一使用此写法。**不改严格坐标isequal或原数值容差**，不改物理数组、接收坐标、相位、源或边界算法。
- `run_d1_depth_comparison(root,mode)`增加非求解`audit-import`，正式execute在恢复配置下导入原100m common和C0，重新SHA/同源同网格核验并独立审计；C0分支不调用BIE求解器。100m PS/FullFK也复用原保存场，不重算。随后才允许100m C1–C4及1490m全新五层，剩余恰好9次BIE。
- 新增`audit_d1_recovery_origin`、`verify_d1_saved_import`、`validate_d1_recovery_preparation`；保存审计输出到独立目录，不修改原C0 MAT或原低精度CSV。
- Python包预检新增恢复scope、原59项SHA/原失败回执/输入/审计证明/common和C0的SHA一致性检查；监督器正确区分新9次与导入1次。
- PowerShell管理工具根据恢复配置自动选择新任务名，拒绝拿恢复包去Start原失败任务。

公共生产API、WAPE、源函数、共享BIE/FullFK求解器和共享R1/R2审计均未修改。`src/`和通信接口无变更，不需要生产/通信仿真或`vertical_comm_guide.md`公式更新。新增入口导航在`scripts/README.md`。

## 3. 本轮实际MATLAB审计

现有MATLAB R2025b Update6/PID22340，经MCP执行；C0审计临时单计算线程，结束恢复原线程设置。原远端是R2025b Update1，恢复启动时仍在远端同一R2025b路径重审保存C0，不把本地审计替代远端运行门。

[C0审计回执](../results/validation/d1_recovery/20261118_01/saved_c0_audit/receipt.json)：8.5348716s，0失败，新增BIE调用0。保持原审计函数/门：

|检查|实际值|原门|
|---|--:|--:|
|修复后CSV坐标最大差|0 m|严格isequal|
|CSV复场最大差|0 Pa|1e-11 Pa|
|保存密度411点接收场重建相对差|1.3834216633e-15|1e-11|
|源/RHS独立128阶核对相对差|1.8153545149e-14|1e-11|
|全部340平台中点残差相对范数|6.9711832769e-8|1e-4|
|独立19点GL112边界残差相对范数|8.2895050636e-8|1e-4|

原386个坐标差≤4.9737991503e-14m保留在原CSV，不追改失败。修复只在审计副本重写CSV；C0的data.mat与common.mat字节保持原样。另测试真实坐标/复场、pi、极大/极小有限double、eps、NaN、字符串及逻辑掩码的往返，均通过。

本项只认证**该保存C0的独立审计**，不是C0–C4收敛、完整D1a或双深度模型误差资格；后续门仍必须执行。

## 4. 最终恢复包与验证

最终目录：`results/comparisons/d1_depth/recovery_preparation_20261118_03/`。

- tar大小19169280 bytes（约18.28 MiB），125项输入/代码/原证据/导入数据/审计SHA。
- tar SHA256：`5f79124e177bde005060786c41ed0fa62ff02756e0ac6efd21b96e1e640e194e`。
- [最终MATLAB准备验证](../results/validation/d1_recovery/20261118_01/final_matlab/receipt.json)：1.1120996s，0失败；真实包dry-run/audit-import、坐标和复场逐位往返、重复覆盖拒绝；7个MATLAB文件Code Analyzer消息0，新增BIE0。
- [离线检查](../results/validation/d1_recovery/20261118_01/final_offline/receipt.json)：23项通过，含125项SHA/原来源链、归档重定位、篡改/失败审计拒绝、禁止重置累计时间/旧RAM基准/剩余调用数、C0分支无求解器、部分结果不PASS及失败Fetch链完整。
- [传输mock](../results/validation/d1_recovery/20261118_01/final_transport/receipt.json)：12项通过，真实网络0；含原任务名拒绝、自动恢复任务名、非覆盖上传、一次手动Start、Status及失败Fetch/原退出码保留。
- PowerShell解析与Git Bash语法检查通过；Linux原生恢复进程、实际Start和剩余声学未运行，不以mock推断正式PASS。

准备01/02候选保留；01的验证脚本有一个多余逗号提示，02已修复并运行通过，最终03进一步明确manifest为`D1A_RECOVERY_WITH_AUDITED_C0`，不再沿用“无密度输入”scope。只启动最终03。

## 5. 原预算不重置

继承原1×1、无parpool、OMP/MKL/OpenBLAS=1、12GiB组RSS上限、原始可用RAM30%预留、压力/磁盘/异常停止门及一次性启动保护。

原监督器133.78413781605195s计入原5400s累计预算；100m也保守计入相同耗时，不用较短的MATLAB89.619s来放宽深度预算。恢复总剩余最多约5266.216s（87.77min），100m剩余最多2266.216s（37.77min），1490m仍最多2400s；以上是硬上限，不是实测预计。RAM基准取原47975735296 bytes与恢复Start时可用量的较大者，绝不向下重置。

恢复包首次运行先重新资源准入、SHA及C0审计，再进入任何新BIE；失败即停止，不自动重试或调参。监督器和最终检查点区分导入1层与本次新9层，保留真实累计耗时、完整stdout/stderr/diary、每秒资源和阶段SHA。Codex不持续监控。用户手动命令见[恢复命令](d1_recovery_manual_commands_2026-11-18.md)。

物理仍是0.10m二维条带/4kHz、原U6/N12/260001/L160/C2、Tx100/1490m、Rx6m、±36m/411直接点、均匀无吸收且仅向上源，无海底/向下辐射/多次底面反射/结构散射，不是完整真实深海多途。D1b、其他Seed/频率、K1及宽带通信不由本包推进。
