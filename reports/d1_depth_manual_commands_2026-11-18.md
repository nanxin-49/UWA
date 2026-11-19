# D1a 远程包及手动命令

2026-11-18。**尚未上传/启动，等待用户授权。** 最终包为03；01/02是被替代的准备候选。仅D1a、Tx100/1490m、Rx6m、±36m/411点、0.10m条带/4kHz、原U6/N12/260001/L160/C2。D1b没有执行入口或自动派发路径。

执行包：[input_package.tar](../results/comparisons/d1_depth/preparation_20261118_03/input_package.tar)，215040 bytes、31项manifest；SHA256 `dcb3197cab1717ba27528b9b1c65b20a4d7bf36fb4993efcac150b962f7dd42e`。完整方案见[D1计划](d1_depth_comparison_plan_2026-11-18.md)，已执行/未执行验证见[准备回执](d1_depth_preparation_validation_2026-11-18.md)。

在 **PowerShell 7** 中使用包内冻结的管理工具：

```powershell
$pkg = 'E:\MISC\CARPE3D_matlab\Explain\results\comparisons\d1_depth\preparation_20261118_03'
$tool = Join-Path $pkg 'code\manage_d1_depth.ps1'
$task = 'd1_depth_20261118_01'

# 仅本地输入和归档SHA核验：不连接SSH，不启动MATLAB。
& $tool -Action DryRun -PackageDir $pkg -TaskName $task

# 以下步骤只在审阅方案并授权D1a预算后由用户手动执行。
& $tool -Action Upload -PackageDir $pkg -TaskName $task
& $tool -Action Start -PackageDir $pkg -TaskName $task

# 需要时手动查询一次；没有自动轮询。
& $tool -Action Status -PackageDir $pkg -TaskName $task

# 任务结束后回收；运行期间会被锁拒绝，失败任务也可回收。
& $tool -Action Fetch -PackageDir $pkg -TaskName $task
```

也可用一次`-Action UploadAndStart`代替分别Upload/Start，**两种方式二选一**。本轮Codex没有执行任何网络命令。目标SSH别名仅`matlab-ubuntu-new`，新目录`/home/misc/cwx/d1_depth_20261118_01`，不覆盖已存在目录；R2025b显式路径，单进程单计算线程。预计25–45min，硬上限90min/单深度40min、12GiB工作估计和RSS上限、启动可用RAM至少18GiB、全程30%启动可用RAM预留及PSI门。

Start返回会话路径只证明后台监督器已派发，不表示声学开始或PASS。Status应保留监督器和MATLAB真实退出码；只有状态`D1A_COMPLETE_WAITING_REVIEW`、双退出0、两深度全部阶段/独立审计/FIG门及SHA均有效，才能进入保存数据解释。任何失败先Status/Fetch检查原目录，不重跑Start、不自动换TaskName规避停止门、不继续D1b。Linux原生与MATLAB运行本轮均未验证，首跑可能在前置门停止。

Fetch会保留远端归档SHA文本回执，下载到唯一新目录并安全解包、校验内部manifest；失败传输目录不覆盖。正式输出包括：每深度common/comparison MAT、五层data MAT/411点CSV/密度/计时/检查点、model_metrics和逐点误差CSV、4组PNG/FIG及重开PNG/审计、总depth_comparison/depth_change CSV、report.md、完整stdout/stderr/diary、资源TSV、身份/真实退出码与SHA清单。

回收后再分析深度是否改变模型误差及其经验可分辨性。Codex不持续监控；D1b新周期/支持/资源与模型适配仍需单独方案和授权。
