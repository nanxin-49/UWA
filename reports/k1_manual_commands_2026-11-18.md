# K1用户手动操作

2026-11-18。只使用最终05包，Codex未上传、未启动K1。10:41压力未过准入；12:21补查PSI下降、只读准入当时通过，原生锁/身份/真实信号11项通过。Swap仍近耗尽，继续1进程×1线程。Start重新检查，未达标即停止，不自动排队或重试。

PowerShell7：

```powershell
$package = 'E:\MISC\CARPE3D_matlab\Explain\results\comparisons\pm_wind_cutoff_k1\preparation_20261118_05'
$tool = Join-Path $package 'code\manage_k1_pm.ps1'

# 离线检查，不连接SSH，不启动MATLAB
& $tool -Action DryRun -PackageDir $package

# 用户上传到唯一目录，拒绝覆盖既有目录
& $tool -Action Upload -PackageDir $package

# 用户启动一次，远程控制器顺序完成六例及汇总
& $tool -Action Start -PackageDir $package

# 用户需要时单次查询/停止后回收；没有轮询循环
& $tool -Action Status -PackageDir $package
& $tool -Action Fetch -PackageDir $package
```

目标别名matlab-ubuntu-new；固定远端目录`/home/misc/cwx/k1_pm_20261118_05`，不使用旧主机。如果当前PowerShell未读取用户SSH配置，可给上述联网动作追加 `-SshConfig "$env:USERPROFILE\.ssh\config"`；不用提供IP或凭据。

输入归档15,535,616 bytes，SHA256：`ddc8abe83a30619852b4f0d6ab0a582dc5b2789ffbf4d29946092ddac57ecf60`。

Start正常后无需Codex监控。每个案例独立MAT/CSV、日志、计时、SHA和阶段检查点；任何案例失败停止后续派发。Status显示最近session的真实退出回执与全局检查点；NO_EXIT_RECEIPT_YET不等于PASS，也不能据此重复启动。

未通过资源准入时不会启动声学MATLAB。待共享负载缓解后用户可再次手动Start。中断恢复同样用Start：跳过SHA有效封存阶段、已求解未审计只补审计、完成状态遗漏只修复状态；不会自动重复已通过案例。未封存的求解调用/数据会拒绝恢复，保留现场，必须人工审查，不删标记蒙混重算。

Fetch仅在任务槽空闲且无已验证存活子进程/不确定fork时允许，保存远端归档SHA、本地SHA和内部逐文件验证回执。传输PASS不等于声学PASS。完成后重点检查：

- `sessions/<id>/exit.json`与每案例`sessions/<id>/{stdout,stderr}.log`、resource.csv、exit.json。
- `cases/<case>/jobs/*/{solved,committed,audit,gates}.json`及独立finalization封存。
- 顶层checkpoint.json应为K1_COMPLETE_SEVEN_CONDITIONS；summary/attempt_*内有14行指标、5754行逐点误差、交互表、原始/C2谱差表、report.md、PNG/FIG及审计回执。

当前许可环境只有基本MATLAB，包不需要PCT；单进程失败不能自动扩进程、换线程或处理账户登录。不要修改数值门或资源基准。正常完成后停止，等待结果审查，不自动进入其他Seed/频率/通信/Bellhop。
