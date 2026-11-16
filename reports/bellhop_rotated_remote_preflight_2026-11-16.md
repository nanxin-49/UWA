# Bellhop旋转坐标迁移Ubuntu：编译器缺失停止（2026-11-16）

远程SSH和MATLAB健康检查通过，**因缺少gfortran停在环境检查阶段**。未迁移A/B/C数据，未编译Linux overlay，未启动最小物理案例或C剩余9项；C保持30/39。依用户指令先报告缺失依赖，不使用sudo安装，不用普通Bellhop替代validation版本。

## 环境实测

| 检查 | 结果 |
|---|---|
| SSH | `matlab-ubuntu-new`，现有密钥非交互登录成功，用户misc |
| Ubuntu/架构 | 24.04.4 LTS / x86_64，56逻辑CPU |
| gfortran | PATH及gfortran-11至15命令均无，dpkg-query未找到gfortran安装记录 |
| make | GNU Make4.3 |
| MATLAB | `/home/misc/MATLAB/R2025b/bin/matlab`，R2025b Update1 / GLNXA64 |
| MATLAB健康检查 | license测试1，单线程算术断言通过，真实退出0；没有声学求解 |
| MATLAB进程 | wall43.26秒，峰值RSS697976KiB |
| 启动前可用RAM | 53634752KiB，约51.15GiB |
| 磁盘空闲 | 5135307260KiB，约4.78TiB |

首次沙箱SSH因不能读取用户SSH配置而无法解析别名；在已授权的正常用户上下文重试连接成功，没有读取私钥内容或修改SSH配置。最初本地CIM诊断等旧实验与本次远程测量无关。

独立远程目录：`/home/misc/cwx/bellhop_rotated_C_20261116_remote01`。仅此目录有本次探测脚本、stdout/stderr、diary、资源记录、退出码及9文件SHA清单；未枚举/修改共享账户其他目录或进程。返回本地后9/9 SHA一致，stderr为空；仅exit0不作为物理验收。

[迁移checkpoint](../results/validation/bellhop_rotated_fig12b/20261116_remote_C_01/checkpoint.json)、[返回SHA核验](../results/validation/bellhop_rotated_fig12b/20261116_remote_C_01/returned_sha256_checks.csv)。原始日志在同目录remote_evidence，Windows所有已完成分量、REAL8构建和旧失败记录未覆盖。

## 源码、预算及未运行门

当前源工作树 `C:\Users\ASUS\.codex\worktrees\bellhop_rotated_fig12b\Explain`，分支experiment/bellhop-rotated-fig12b、HEAD f2cea1a。冻结C计划SHA为cce53c33f09104675eed7ef6131b6b9d6154259b25baaeffa83570d331947e97；已验证Windows优化REAL8 binary SHA为4d3b99f4365f778eff40b52714794eeee627b1f5b54af51baeba8b0b811ecdf2。这些是现有checkpoint记录，本轮没有重复大规模SHA或求解。

Windows历史账本原剩859.165599秒；本轮环境探测及返回SHA保守计60秒，迁移账本累计2800.834401秒、余 **799.165599秒**。旧Windows账本和checkpoint原样保留，新账本不重置或增加原3600秒额度。健康检查实际wall与保守计账分别记录，后续编译/跨平台最小验证/求解审计须按适用计账及停止门核余额。

Linux编译、REAL(KIND=8)/1000d0/ReceiverRange64有效求值链、优化Step、编译参数、SHD结构与MATLAB读取兼容性均**未验证**；没有“跨平台通过”结论。尚未建立可复用Linux executable，不假设远程装有Bellhop。

以下全部未执行：必要源码/分量允许清单快照与迁移、Linux构建、101束3列坐标/反射计数/交点方向/pq/tau/Amp/Phase/复声压对照、C9、完整C组装、截断和B→C收敛、ABC图册。没有新复场、PNG或FIG，不能把环境MATLAB健康检查当成物理资格。

## 下一步边界

需管理员安装gfortran，或给出已安装工具链的明确路径后重新核验版本；本轮未安装任何依赖。依赖满足后使用同一已授权迁移范围，创建新的独立快照目录并逐项SHA核验，不覆盖既有远程文件。先移植/检查Linux调度与构建，再执行原101束3列跨平台物理门，只有通过才复用C30并运行剩9项。跨平台原生bin不要求SHA相同，但源码来源和数值合同必须验证；原始坐标/网格/模型/声束及门槛不变。

执行仍串行；保留原启动RAM1.5GiB、运行256MiB最低门、单Bellhop120秒、总数值预算及输出门。远程另保留至少启动可用RAM的30%作为共享资源余量，磁盘至少20GiB；实际运行前重测而不引用本次51.15GiB作未来许可。不得终止共享用户其他进程，不D/新FFP、不提交/推送/合并。

已更新任务导航与进度源并保存失败依赖记录。本轮没有生产/物理接口或公式改变，vertical_comm_guide无需修改；没有将C阶段或最终论文复现标为完成。
