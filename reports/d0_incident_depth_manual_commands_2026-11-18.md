# D0手动上传、启动、查询与回收

2026-11-18。仅使用 **PowerShell 7 (`pwsh`)**。Codex本轮未执行以下网络命令。使用最终06包，不使用01–05候选或TEST_ONLY测试目录。

```powershell
$Package = 'E:\MISC\CARPE3D_matlab\Explain\results\experiments\d0_incident_depth\preparation_20261118_06'
$Tool = Join-Path $Package 'code\manage_d0_incident_depth.ps1'
$TaskName = 'd0_incident_depth_20261118_01'

# 只核冻结配置、输入、代码及tar SHA，无网络、无MATLAB进程
& $Tool -Action DryRun -PackageDir $Package -TaskName $TaskName

# 用户手动执行一次：新目录上传、核SHA，随后后台顺序启动五深度
& $Tool -Action UploadAndStart -PackageDir $Package -TaskName $TaskName

# 按需单次状态查询，不会启动或重启计算
& $Tool -Action Status -PackageDir $Package -TaskName $TaskName

# 确认任务停止后回收；运行中会拒绝。失败也可回收，完整保留失败证据
& $Tool -Action Fetch -PackageDir $Package -TaskName $TaskName
```

上传和启动如需分开，可分别执行 `-Action Upload` 和 `-Action Start`。已经使用UploadAndStart后不要再执行Start；one-shot标记会拒绝重复派发。远端目录已存在时Upload拒绝覆盖，不自动清理；遇传输、许可证或启动失败，先检查原Status及回收日志，不另起同一任务或删标记。

远端任务目录为 `/home/misc/cwx/d0_incident_depth_20261118_01`。后台单R2025b、`-singleCompThread`、无parpool；30%初始可用RAM预留、3600 s会话/900 s单深度上限。超过门即停止，不自动放宽参数或重试。Codex不需要保持在线，也不会持续监控。

Fetch默认创建同级全新 `<TaskName>_retrieved_<时间>` 及 `<该目录>_transport`。后者保留tar与远端archive SHA文本回执；先核传输SHA，再拒绝路径越界/链接/设备成员，再核内部输入及产物SHA。可用 `-DownloadDir 'E:\...\new_unique_directory'` 指定不存在的目的目录。

成功回收只说明完整性通过，不代表D0数值或物理结论通过。把工具输出的 `FETCH_VERIFIED=...` 路径发给Codex，下一轮统一检查真实退出码、五深度数值门、MAT/CSV/FIG/PNG及波前/角谱，形成研究报告。尚未执行D1、PS/BIE/FullFK或宽带通信。

最终归档SHA256：`a63dbbbd8f092fd8d85c04ceca8e3f8c65a23d21a6007f07aa5a9515cfd8c3dd`。包内 `input_package.tar.sha256` 和 `manifest.json` 可独立核验。
