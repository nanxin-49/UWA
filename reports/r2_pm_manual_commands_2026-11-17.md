# R2 PowerShell 手动操作说明

2026-11-17。最终交付为 **03 版**，仅 seed260002/260003，固定总宽0.10 m、4 kHz、411真实接收点。**下面的远程命令本轮均未执行。** 两个 Seed 可同时各占一个单线程 MATLAB 进程；每个包不启动另一 Seed。01/02 为本地准备候选，不上传它们。

使用 Windows PowerShell 7、现有 OpenSSH `ssh/scp`、`tar` 和 Python 标准库；仅使用批准别名 `matlab-ubuntu-new`。脚本中的任务目录与单线程 R2025b 路径沿用 R1；本轮没有检查远端当前安装或内存。Python 仅做封存/传输检查，MATLAB 承担全部声学和科学分析。

## 1. 本地变量及无网络检查

```powershell
$repo = 'E:\MISC\CARPE3D_matlab\Explain'
$manager = Join-Path $repo 'scripts\validation\support\manage_r2_pm_seed.ps1'
$bundle = Join-Path $repo 'results\comparisons\piston_2d_pm\r2_preparation_20261117_03'
$p2 = Join-Path $bundle 'seed260002'
$p3 = Join-Path $bundle 'seed260003'

& $manager -Action ShowCommands -PackageDir $p2
& $manager -Action ShowCommands -PackageDir $p3
python (Join-Path $p2 'code\r2_pm_preflight.py') $p2
python (Join-Path $p3 'code\r2_pm_preflight.py') $p3
```

ShowCommands 和 Python preflight 不连接网络、不启动 MATLAB。交付归档各自有 `input_package.tar.sha256`，不能改冻结文件后继续沿用旧包。

## 2. 手动上传，随后启动两个后台进程

**用户决定开始计算时才执行本节。** Upload 只上传/校验/解包，不启动；Start 返回后台会话请求，不等待整个计算结束，因此按下面先后发出两次 Start 可让两 Seed 并行。

```powershell
& $manager -Action Upload -PackageDir $p2
& $manager -Action Upload -PackageDir $p3

& $manager -Action Start -PackageDir $p2
& $manager -Action Start -PackageDir $p3
```

远端目录分别为：

```text
/home/misc/cwx/r2_pm_seed260002_20261117_03
/home/misc/cwx/r2_pm_seed260003_20261117_03
```

两个包共享 `/home/misc/cwx/r2_pm_20261117_03_control/` 的 RAM 基准和各自 Seed 槽锁。每次启动要求当前 `MemAvailable − 30%共享基准 ≥ 24 GiB`；首启约需至少34.286 GiB可用内存。若准入失败则停止，不自动退到串行、减少网格或放宽门。共享基准文件保留，不自行删除以绕过门。

进入声学计算前会依次执行：包 SHA → Linux 原生并行/RAM边界与锁验收 → 资源准入 → 同一个 R2025b 单线程进程内带七个声学哨兵的 MATLAB dry-run/恢复与政策测试/Code Analyzer → 本 Seed 正式入口。任何验收异常都会在声学求解前退出。MCP 当前不可连接，MATLAB 检查本轮未运行；Windows Git Bash 没有 flock，Linux 原生锁测试也本轮未运行。

不建议双击 `.sh` 或手工直接运行 `.m` 正式入口：launcher 负责锁、单线程、环境、RAM watchdog、日志和真实退出码。不要把脚本中的 run 子命令当成第二次 Start。

## 3. 查询状态

```powershell
& $manager -Action Status -PackageDir $p2
& $manager -Action Status -PackageDir $p3
```

Status 只读，不启动、恢复、修改或终止任务。启动前后几十秒可能处于 MATLAB_INITIALIZING；没有退出码不是失败，先看本任务 PID、阶段、stdout/stderr 和 launcher 日志，不要重复启动。

每个会话保留 `sessions/<session>/stdout.log`、`stderr.log`、launcher 日志、`runtime.json`、`shared_runtime.json`、`resources.csv`、`time.txt`、`matlab_exit_code.txt`、`process_elapsed_s.txt` 与输出 SHA。检查点位于包根 `checkpoint.json`，各 C 层有独立 solved/committed 封存。进程退出0不是数值 PASS。

## 4. 失败后有条件手动恢复

先查询确认旧进程已退出，再看失败层和原日志。只有已封存 solved 需要补审计，或已 committed 需要跳过时，才使用：

```powershell
& $manager -Action Continue -PackageDir $p2
# 或只对失败的另一个 Seed：
& $manager -Action Continue -PackageDir $p3
```

Continue 保留同包配置和 SHA，不重算已验收层。若求解调用已开始但未封存、数据孤立、SHA不符或前序变化，会停止；**不能删除 marker、改门或多次 Continue 来请求重求解。** 此类情况需另行诊断和明确恢复授权。

确需人工中止时，先用 Status 确认准确会话，再在远端向该会话所记录的 launcher/controller 发送 TERM，使其只终止自身记录的 MATLAB 进程组；不使用 `pkill MATLAB`、不清理其他作业。没有自动墙钟超时，10–20 min/Seed 只是成本观察窗口；资源/数值异常会自动停止。

## 5. 回收（任务停止后）

```powershell
$d2 = Join-Path $repo 'results\comparisons\piston_2d_pm\r2_seed260002_retrieved_20261117_01'
$d3 = Join-Path $repo 'results\comparisons\piston_2d_pm\r2_seed260003_retrieved_20261117_01'
& $manager -Action Fetch -PackageDir $p2 -DownloadDir $d2
& $manager -Action Fetch -PackageDir $p3 -DownloadDir $d3
```

目标目录必须不存在。Fetch 拒绝正在计算的包；按独立回收清单逐项验证 SHA。失败或停止的包也可回收，仍按真实状态审查，不能因为传输 SHA 通过就宣称参考合格。

两 Seed 都经过保存结果审查后，才运行三 Seed 汇总。入口在现有 MATLAB 会话中调用，不通过 PowerShell 启动本地 MATLAB。下面的 R2 `work` 路径由各自回收的 `finalization/pointer.json` 确认，不能猜测失败 attempt：

```matlab
setup_vertical_project;
r1 = 'E:\MISC\CARPE3D_matlab\Explain\results\comparisons\piston_2d_pm\r1_bie_manual_c4_policy_v2_20261116_04_retrieved_20261116_173609\finalization\attempts\attempt_0001\work';
% r2、r3 分别设为回收目录下 pointer.json 指向且经审查的 work 绝对路径。
summary = summarize_r2_pm_fields({r1,r2,r3}, ...
    'E:\MISC\CARPE3D_matlab\Explain\results\comparisons\piston_2d_pm\r2_three_seed_summary_20261117_01');
```

这一步只读保存场，生成六行指标、逐点/分段能量和排序表；本轮没有执行。完成两个 Seed 或汇总也不授权 R3、宽带或通信计算。

技术定义、停止门与资源估计见 [R2方案](r2_pm_multiseed_plan_2026-11-17.md)，本轮实测检查及未运行项见 [本地验证报告](r2_pm_preparation_validation_2026-11-17.md)。
