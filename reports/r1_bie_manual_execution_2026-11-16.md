# R1-BIE C1–C4用户手动执行包

2026-11-16。**本轮仅准备脚本与无声学dry-run，不连接远程主机、不启动C1–C4。** R1参考仍未完成资格。

## 已准备的内容与政策变更

手动包：`results/comparisons/piston_2d_pm/20261116_r1_bie_manual_ready_03/`。

- MATLAB批处理入口：`scripts/validation/run_r1_bie_manual.m`，缺省模式为`dry-run`；正式launcher显式传`run`。
- Linux启动/状态脚本：`scripts/validation/support/launch_r1_bie_manual.sh`。
- 回收封存脚本：`scripts/validation/support/archive_r1_bie_manual.sh`。
- Windows上传/启动、状态、恢复、回收工具：`scripts/validation/support/manage_r1_bie_manual.ps1`，缺省动作只打印命令。
- 可复现打包工具：`scripts/validation/support/prepare_r1_bie_manual_package.ps1`；只操作本地白名单文件，不连接SSH。
- SHA/原子检查点工具、逐阶段独立重建审计、最终参考分辨率审计均单列在`validation/support/`，已随包冻结。

本包采用用户本轮的新执行方式，**不调用旧成本门控制器**：取消1100 s累计硬停止和730 s预测准入，不继承旧S0/C0及单会话倒计时；不按575 s名义余额安排任务。旧累计525 s仅作来源记录。一个R2025b `-singleCompThread`进程顺序执行，不建parpool，不以旧900 s会话上限截断完整队列；记录真实时间，但不设置新的时间/成本硬停止。保留30%可用RAM、每阶段异常/数值/审计门。此政策仅适用于本用户手动启动的固定C1–C4，不授权Codex自行启动或扩展研究范围。

## Windows启动、状态、恢复、回收命令

在PowerShell中设置一次：

```powershell
$p = 'E:\MISC\CARPE3D_matlab\Explain\results\comparisons\piston_2d_pm\20261116_r1_bie_manual_ready_03'
$ctl = Join-Path $p 'code\manage_r1_bie_manual.ps1'
```

**首次上传并后台启动**（只有用户运行此命令才连接远程并开始计算）：

```powershell
& $ctl -Action UploadAndStart -PackageDir $p
```

工具先核本地输入/源码与tar的SHA，使用`matlab-ubuntu-new`，创建不存在的独立目录`/home/misc/cwx/r1_bie_manual_20261116_01`，上传、解包后后台启动。不会覆盖已有远端任务目录；如果该目录已存在，不要重发UploadAndStart，先查状态。上传中断时保留远端部分目录，需检查后处理，不做删除/强制同步。

**查状态**：

```powershell
& $ctl -Action Status -PackageDir $p
```

**异常停止后，用户确认需要继续时恢复**：

```powershell
& $ctl -Action Resume -PackageDir $p
```

这与重新执行远端`start`相同；独占锁及旧任务进程会话检查拒绝并发重复启动，不会重算已审计通过且SHA完整的阶段。

**停止/完成后封存并取回结果，自动逐文件核验SHA**：

```powershell
& $ctl -Action Fetch -PackageDir $p
```

回收工具拒绝对仍运行的任务制作一致性快照；默认创建包目录旁边一个新的`r1_bie_manual_20261116_01_retrieved_<时间>/`，不会覆盖已有结果。可显式传`-DownloadDir 'E:\...\一个尚不存在的目录'`。全流程不修改生产文件、不提交/推送Git。

直接SSH等价命令（已上传后）：

```powershell
ssh -o BatchMode=yes matlab-ubuntu-new "bash /home/misc/cwx/r1_bie_manual_20261116_01/code/launch_r1_bie_manual.sh start /home/misc/cwx/r1_bie_manual_20261116_01"
ssh -o BatchMode=yes matlab-ubuntu-new "bash /home/misc/cwx/r1_bie_manual_20261116_01/code/launch_r1_bie_manual.sh status /home/misc/cwx/r1_bie_manual_20261116_01"
```

不要求Codex轮询。后台脚本自行记录资源并在RAM门触发时仅停止本任务创建且核实过的独立进程组，不访问/终止其他用户任务。

## 固定任务与C0复用

包内包含本机已回传并核对SHA的`20261116_r1_bie_cache_01` C0原生密度、完整反射场、CSV、等价性回执、独立重建回执和检查点。缓存两个核心源码逐字节绑定已验证版本；**本轮未修改它们**。正式新阶段显式`enable_static_cache=true`，原接口默认仍关闭。

| 阶段 | L / m | ppw | qa | N | 全平台中点 |
|---|---:|---:|---:|---:|---:|
| C0导入，不求解 | 60 | 10 | 96 | 3200 | 340 |
| C1 | 60 | 12 | 96 | 3840 | 408 |
| C2 | 60 | 12 | 64 | 3840 | 408 |
| C3 | 60 | 14 | 96 | 4480 | 476 |
| C4 | 70 | 12 | 96 | 4480 | 476 |

固定10 cm条带、v0=0.001 m/s、Tx100 m、4 kHz、c1500 m/s、rho1025 kg/m³、PM seed260001原12模/160 m周期/C2、h=−2 m、GL8、窗口比例0.85，原411个真实接收位置z=3 m，native相量和phase contract v2。不存在C0求解代码路径，也不重做小网格等价性实验。C0可重新审计保存数据，但不能调用C0求解器。

每个新阶段完整求解全部平台中点残差，线性残差<1e-10，离节点边界残差<1e-4；原四组全线/四分段反射场配对差均<1e-4。保留原计划对非单调空间细化的停止检查，不删检查点/不减少接收点。失败即停止后续阶段，不自动降低门、增加网格或改变模型。

## 每阶段记录、独立审计与恢复

各阶段位于`jobs/bie_coarse、bie_main、bie_self、bie_space、bie_window/`。每阶段保存：

- `data.mat`：密度、修正/完整反射场、入射边界数据和明确数值配置；`points.csv`：411点原生复声压。
- `config.json`：阶段配方、源、相位、输入/代码指纹；`timings.json`：缓存/组装/求解/接收/边界等阶段时间。
- `audit.json`及重建MAT：独立Green核从密度重建411点，MAT/CSV闭合、128阶孔径入射核对、全平台中点支持检查、另选19点GL112边界重建。
- `gates.json`：已可检查的配对差及原门；`solved.json`与`committed.json`：含SHA的求解封存和通过审计的检查点。
- `attempts/attempt_XXXX/`：每次尝试的diary、前后RAM记录及异常全文/此前失败审计；保留失败，不自动重试。

连续资源采样在`control/resources/C*.csv`和`sessions/<id>/resources.csv`；这些随运行追加的日志不放进不可变科学产物commit，避免采样竞争导致假SHA漂移。其完整SHA由会话结束/结果回收时封存。阶段commit绑定密度/场/配置/计时/独立审计/数值门，且包含包指纹和前一通过阶段commit的SHA，不能跨代码、输入或前序状态盲用。

恢复规则：

1. `committed.json`及绑定文件、审计、数值门完整且输入指纹一致：直接`SKIP_AUDITED`，不重复求解或审计。
2. `solved.json`及密度/场/CSV/配置SHA有效，但尚无最终commit：`AUDIT_ONLY`，只补审计与检查点，不重算声场。
3. 求解中断、尚未留下完整`data.mat/points.csv`：下次用户手动恢复时新建attempt，仅重做该未通过阶段，保留旧日志。
4. 已通过阶段被修改、输入/代码漂移，或存在没有封存绑定的孤立`data.mat/points.csv`：停止供审查；不默默覆盖或重算。应恢复原SHA文件，不要删除已通过检查点以强迫重算。
5. 数值/独立审计失败：后续阶段不启动。相同输入下重启只会复核同一已封存数据并再次报告失败，不会自动换参数。

启动器使用独占`flock`、旧任务进程会话检查和新的session目录；意外退出后，先读状态和同一session的完整日志，不先重复启动。仍有任务进程存活时拒绝第二个MATLAB进程。

## 全部阶段通过后的产物

自动在新的`finalization/attempts/attempt_XXXX/work/`工作目录完成最终处理，失败重启只创建新的最终处理attempt，不重算已通过的BIE阶段：

1. 由已保存PS/完整FullFK及新BIE参考更新411点误差表。
2. 独立重建审计验证MAT/CSV、源、BIE密度/边界、保存的PS与FullFK字段；不修改或重新求解PS/FullFK传播模型。
3. 独立计算四组BIE敏感性和逐点经验不确定性，更新`5*norm(bunc+model_unc)/norm(BIE_final)`及分辨标记。
4. 生成最终CSV、MAT、报告与F01–F04 PNG/FIG，重开12条误差曲线并核两图重导出像素；保留各次导出失败。
5. 最终commit绑定交付产物SHA，root `checkpoint.json.final_report`指明报告路径。完整后同样停止等R1审查，**不进入R2**。

整个MATLAB进程的完整stdout/stderr、命令、PID/进程组、`time -v`和**真实退出码**在`sessions/<id>/`。阶段PASS是函数数值/审计状态，不伪造独立阶段OS退出码。MATLAB检查点显示完成而尚未出现`matlab_exit_code.txt`时，进程还未封账；退出0本身也不代替数值与审计PASS。

## 时间、内存与未实跑部分

C0实测求解调用55.48 s。原保守N³估算给C1/C2各约96 s、C3/C4各约152 s；加25%余量及收尾约730 s。该数字仅用于安排等待，**不作派发门或时间上限**。建议预期约10–15分钟，负载、独立审计和图形导出可能使时间更长。

C0峰值RSS约1.23 GiB；N4480一份复双精度矩阵约306 MiB。更细阶段含矩阵/分解/缓存的峰值尚无实测，工程预期数GiB，启动器暂用保守12 GiB增量准入估计（须容纳于当前可用RAM的70%）。实际运行监测至少保留30%启动可用RAM，恢复时采用原任务起始与本次起始可用量中的较大值作基准。RAM不足会明确停止；不以旧“内存充足”替代新预检。

本轮本地检查涵盖：PowerShell解析、Bash语法、Linux风格SHA清单校验、MATLAB Code Analyzer、带7个声学入口哨兵的纯dry-run，以及fresh/audit-only/skip、前序变化/文件篡改/孤立数据拒绝的恢复状态测试。没有调用BIE、PE、FullFK或Hankel声学求值；没有SSH连接。首次清单CRLF问题已在本地dry-run捕获，最终清单固定LF。

本轮不能证明真实C1–C4会通过数值门；远端后台锁/信号/资源守门和真实完整队列未实机运行，Linux依赖/版本/许可在用户启动时检查。任何此类失败均保留日志并停止，不自动安装软件、切换环境或重算已通过阶段。生产WAPE、PS、FullFK、通信默认及接口不变，不运行其他seed/频率/CIR/QPSK/BER/Bellhop，不提交/推送Git。
