# K1恢复04：用户手动上传、启动、状态与回收

仅使用`recovery_preparation_20261118_04`。原05及原失败现场不动，01–03为准备候选，不能启动。Codex本轮只启动了三项已结束的非声学/保存图件诊断，没有上传或启动正式恢复任务。

在Windows PowerShell执行：

```powershell
$package = 'E:\MISC\CARPE3D_matlab\Explain\results\comparisons\pm_wind_cutoff_k1\recovery_preparation_20261118_04'
$tool = Join-Path $package 'code\manage_k1_pm_recovery.ps1'

# 离线双链检查，不连接SSH或计算声场
& $tool -Action DryRun -PackageDir $package

# 只上传到全新目录，拒绝覆盖
& $tool -Action Upload -PackageDir $package

# 手动启动一次，远端后台自行执行
& $tool -Action Start -PackageDir $package

# 按需查询；完成或停止后回收
& $tool -Action Status -PackageDir $package
& $tool -Action Fetch -PackageDir $package
```

联网工具使用`matlab-ubuntu-new`，新根为`/home/misc/cwx/k1_pm_recovery_20261118_04`；必要时联网动作追加 `-SshConfig "$env:USERPROFILE\.ssh\config"`。不要复用原05的`$tool`或TaskName。

归档99,563,520bytes，SHA256 `77bce9c2b26806ea796758bcdb1662359ecf1f560d7dfa1f3cf3282223dc3f6e`。Upload校验外层SHA并做双链预检，Start再验证输入与实时资源。

启动后首先无声学准备验收与Ubuntu N7保存图件验收；通过后跳过N7，N40只允许一次带事件日志的MODELS计算，再独立C0–C4，随后其余四面。1进程×1计算线程，继承724s/10800s、119021486080bytes以上RAM基准/30%预留/12GiB与原数值门，最多25次新BIE。N40再次失败或任一原门失败停止后续派发，不自动重试，不删未封存调用。

Start成功返回SESSION后可以关闭本地终端。不要因Status暂无退出回执就再次Start；需要时再单次Status。Codex不持续轮询。

进度重点文件：

- `sessions/<session>/exit.json`：控制器真实状态；`preflight`与`saved_graphics`必须先通过。
- `cases/U6_N40/model_stage/recovery_events.jsonl`：INPUT/DIRECT/FLAT/PS1–6/FullFK1–4/COMMON_SAVE/MODEL_GATES/SAVED_MODEL_AUDIT的BEGIN/END。
- `cases/<case>/sessions/<session>/{stdout.log,stderr.log,resource.csv,exit.json}`：实际退出、阶段/模型事件与每秒资源。
- `control/accounting.json`及`case_accounting.json`：继承累计及单例时间。
- 失败时的`owned_crash_evidence.json`：仅该PID转储公开摘要；没有文件也会记录缺失，不能据此重试。
- 完成时顶层应为`K1_COMPLETE_SEVEN_CONDITIONS`，`new_BIE_invocations=25`、`imported_BIE_invocations=5`、累计30，汇总14行指标/5754逐点及真实来源列、PNG/FIG和原始/C2谱差表。

Fetch只在任务槽空闲、无存活本任务子进程/不确定fork时允许，保留外层SHA及内部逐文件SHA。回收后交由保存数值/图件审查；传输成功不代替数值审计。

Windows本地导出画布问题仍保留FAIL；Ubuntu执行平台原严格图件门已经保存重放通过，正式Start再次核验。原N40崩溃具体原生触发组件仍未定位，此包提供有界新尝试与更完整证据，不保证故障不再发生。详见[准备报告](k1_recovery_preparation_2026-11-18.md)。
