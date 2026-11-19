# D1 恢复包：直接手动启动

2026-11-18。修复和保存C0独立审计已通过，最终恢复03包尚未上传/启动。它在新目录导入已审计的100m C0，不重解C0，不重算原100m PS/FullFK；先重审导入数据，再做剩余9层BIE和最终审计。

执行包：[input_package.tar](../results/comparisons/d1_depth/recovery_preparation_20261118_03/input_package.tar)，19169280 bytes、125项manifest，SHA256 `5f79124e177bde005060786c41ed0fa62ff02756e0ac6efd21b96e1e640e194e`。详细[修复与审计报告](d1_recovery_preparation_2026-11-18.md)。

打开 **PowerShell 7**，完整粘贴：

```powershell
$pkg = 'E:\MISC\CARPE3D_matlab\Explain\results\comparisons\d1_depth\recovery_preparation_20261118_03'
$tool = Join-Path $pkg 'code\manage_d1_depth.ps1'

& $tool -Action DryRun -PackageDir $pkg
& $tool -Action UploadAndStart -PackageDir $pkg
```

工具自动使用新任务名`d1_depth_recovery_20261118_01`，目标`matlab-ubuntu-new:/home/misc/cwx/d1_depth_recovery_20261118_01`，不会覆盖旧失败任务。这个恢复包自含原来源和C0副本，不需要在远端搬改原失败目录。UploadAndStart只执行一次；不要同时再调用Upload/Start，也不要使用旧03首次包或旧TaskName。

后台成功派发后可以关闭本地终端；`SESSION=...`只是派发回执，不是数值PASS。以后重新打开PowerShell时重新定义上述`$pkg/$tool`，按需执行一次：

```powershell
& $tool -Action Status -PackageDir $pkg

# 结束后回收；运行期间会被锁拒绝，失败也可回收。
& $tool -Action Fetch -PackageDir $pkg
```

Status显示`completed_steps`、检查点耗时及退出回执；最终只有`D1A_COMPLETE_WAITING_REVIEW`、MATLAB/监督器退出0/0、导入1层与新9层计数匹配及全部数值/审计/FIG/SHA通过，才进入双深度结果解释。失败先Fetch，不重复Start或换任务名规避停止保护。

单进程单线程、30%RAM预留、12GiB RSS、原累计90min/单深度40min预算保持；原133.784s计入累计与100m预算，剩余上限约87.77min/100m37.77min，1490m40min。Start重新测资源，原RAM基准不降低。Codex不会自动启动或持续查询；D1b/宽带仍不进入。
