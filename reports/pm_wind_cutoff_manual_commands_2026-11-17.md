# K0 用户手动上传、启动、查询和回收

2026-11-17。以下命令**本轮未由Codex执行**。最终交付只用05包；01–04是保留的准备候选/测试版本，不上传。目标SSH别名只用matlab-ubuntu-new，远端独立目录`/home/misc/cwx/pm_wind_cutoff_k0_20261117_05`；不覆盖旧目录。需要PowerShell7、本机Python标准库、tar、ssh/scp；远端需要python3、bash、flock、nohup和已可批处理授权的Ubuntu MATLAB R2025b。

在PowerShell中执行：

```powershell
$package = 'E:\MISC\CARPE3D_matlab\Explain\results\experiments\pm_wind_cutoff\preparation_20261117_05'
$tool = Join-Path $package 'code\manage_pm_wind_cutoff_k0.ps1'
$task = 'pm_wind_cutoff_k0_20261117_05'

# 仅显示命令，不建立网络连接
& $tool -Action ShowCommands -PackageDir $package -TaskName $task

# 上传独立快照；本地/远端均核SHA，不覆盖已有任务目录
& $tool -Action Upload -PackageDir $package -TaskName $task

# 用户手动开始，仅K0表面统计、单线程顺序；不会开启K1
& $tool -Action Start -PackageDir $package -TaskName $task
```

归档SHA256：`ac8e3845271f8c04662f76786010bfcc38f59aebab03a67d166a35e131b27970`。本地再核对：

```powershell
Get-FileHash -LiteralPath (Join-Path $package 'input_package.tar') -Algorithm SHA256
```

用户需要时单次查询（工具无循环，Codex不后台监控）：

```powershell
& $tool -Action Status -PackageDir $package -TaskName $task
```

`NO_EXIT_RECEIPT_YET`只表示还没有退出回执，可能启动中、运行中或中断；不能据此判断PASS，也不要重复启动。槽锁忙时Start/Fetch拒绝操作；若MATLAB许可失败，先阅读该session退出回执和stderr，人工处理许可，不能自动重试或改用旧主机。

用户确认任务停止后回收（默认新建带时间戳的本地目录）：

```powershell
& $tool -Action Fetch -PackageDir $package -TaskName $task
```

Fetch核tar SHA和逐文件SHA；下一步仍需审查`sessions/<id>/exit.json`（matlab_exit_code和supervisor_exit_code）、45组`output/<case>/checkpoint.json`、各audit.json、最终`output/checkpoint.json`、CSV和图件。完整退出0且独立审计/45组资格齐备才能报告正式K0完成。只检查文件存在或传输SHA不够。

已保存generated.json但审计未提交、或前面若干组已通过时，修复外部运行条件后，**用户可以再次手动Start同一任务**：会创建新session，核SHA后跳过已审计案例，只补未通过审计的保存数据。没有定时重试。源码/输入/产物SHA漂移会停止；若仅有generation_started.json而没有generated.json，则本包拒绝自动重生，先保留现场并人工审查，在新任务目录安排未完成项，不删除标记或改写通过检查点。

完成的K0再次Start只核验通过案例，不重复生成；K1入口`run_pm_wind_cutoff_k1()`只输出候选计划，传入execute等参数会报`PM:K1Disabled`。K0结果审查与明确K1范围/资源许可之前，不运行声场或其他阶段。

如需用户自行直连启动，与PowerShell Start等价的单次命令为：

```powershell
ssh -o BatchMode=yes matlab-ubuntu-new 'bash /home/misc/cwx/pm_wind_cutoff_k0_20261117_05/code/launch_pm_wind_cutoff_k0.sh start /home/misc/cwx/pm_wind_cutoff_k0_20261117_05'
```

本轮不执行此命令；不要修改启动器为parpool、多线程或同时派发第二K0包。不需要sudo、安装软件或访问其他用户文件。
