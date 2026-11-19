# D1b-v2 BIE 规模扩展试点准备

2026-11-19。**本地准备完成，等待用户审查和手动运行。** 本轮 SSH、远程 MATLAB、4096/8192 正式求解均为 0。原型 02 只用于冻结输入和选择 2-worker 基线。 [原型 02 资源审查](d1b_v2_resource_review_2026-11-19.md)

## 固定设置与顺序

直接从原型 02 已审计封存取同一个 N160、低频关闭的 L640 PM 面，`U_PM=6 m/s`、Seed260001、C2 平台 ±220 m、支持 ±240 m。复用 Tx1490 m、Rx6 m、0.10 m 二维条带、4 kHz、c=1500 m/s、rho=1025 kg/m³、0.001 m/s 峰值振速及 ±200 m 2,279 个直接接收点。BIE 为原 02 缓存核，12 ppw、面板 GL8、自项 GL96、辅助平面 z=−2 m、数值窗平台比例 0.85。隔离生成的 MATLAB 文件只增加分阶段资源记录与串行对照，原核、生产接口及误差容差沿用。

| 顺序 | 节点 | 计算半宽 | 配置 | 下一阶段放行条件 |
|---|---:|---:|---|---|
| 1 | 1,536 | 24 m | 串行 | 保存矩阵、密度、接收场、边界和独立审计通过，SHA 封存 |
| 2 | 1,536 | 24 m | 2 Process workers | 上述完整串行对照与独立审计通过，SHA 封存 |
| 3 | 4,096 | 64 m | 同一 2-worker pool | 先释放活跃稠密矩阵，再从保存矩阵分块重算残差、独立 GL112 与逐 worker 门通过，SHA 封存 |
| 4 | 8,192 | 128 m | 同一 2-worker pool | 重新验证 4096 保存数据、审计、独立残差 SHA 后才允许进入；本档通过后关闭池并结束 |

4096、8192 的单复双精度稠密矩阵分别为 0.25、1.00 GiB，阶段顺序执行。两档都在主进程逐块串行重算全部基础矩阵行，逐 worker 检验。在固定密度上，每 worker 独立选取至多 12 个接收行、12 个平台中点串行重现；19 个 GL112 独立点全部另以串行执行并逐 worker 对照。近项替换和 `A\rhs` 仍由客户端串行完成。完整的两种模式密度和接收复场对照在 1536 节点完成。规模梯度用于资源评价，两个计算半宽均小于目标 ±300 m，后续仍需 C2 边缘、非局域散射、固定海面窗口和 C0–C4 目标参考门。

矩阵、RHS、接收复场一致性 `<1e-11`，密度、中点及独立审计串并行一致性 `<1e-10`，保存系统线性相对残差 `<1e-10`，所有原型平台中点总体及独立 GL112 边界残差 `<1e-4`。两档分别保存矩阵、RHS、密度和审计原始向量。4096 任一门失败即停止。22400 节点仅用于成本算术，正式执行入口的最大配置为 8192。

## 固定资源合同

- `matlab-ubuntu-new`，单 MATLAB R2025b 主进程，PCT checkout，显式 `Processes` profile 仅两个 worker。主进程 `-singleCompThread`，worker `NumThreads=1`，实际登记 PID 与计算线程数。唯一任务目录与原型共享互斥锁，所有规模严格顺序。
- 沿用原型 02 的 12 GiB 任务 RSS 上限、至少 18 GiB 启动可用空间及起始 MemAvailable 的 30% 固定 RAM 储备。设置 1,800 s 总预算，1,770 s 停止计算以预留终止处理；启动至少 20 GiB 磁盘、5,000 inode，运行空闲磁盘至少 2 GiB，输出体积上限 8 GiB。继承启动 load1≤21、memory/io full avg10≤5%/10% 和运行压力门。RSS 或可用 RAM 失败时终止任务自身进程，保存实际退出与已完成档结果，自动重试关闭。
- 监督器约每 250 ms 记录 MATLAB 与 worker、任务进程组和所有归属进程的 RSS、MemAvailable、逐 PID CPU ticks/线程数和 CPU 核当量。MATLAB 以阶段事件保存缓存、组装、边界、独立审计、写盘等耗时和 `ticBytes`/`tocBytes` 的双向逐 worker 流量。监控区分采样峰与真实瞬时峰，阶段时间包含通信和调度。逐阶段文件 SHA、日志和退出状态均可回收。

## 启动与回收

PowerShell 7 可从任意文件系统目录执行，下面不切换当前目录。先 DryRun，通过后再分别手动 Upload、Start。

```powershell
$tool = 'E:\MISC\CARPE3D_matlab\Explain\scripts\validation\support\manage_d1b_scale.ps1'
$package = 'E:\MISC\CARPE3D_matlab\Explain\results\validation\d1b_scale\preparation_20261119_04'
$run = 'd1b_scale_20261119_01'
& $tool -Action DryRun -PackageDir $package -TaskName $run
```

```powershell
& $tool -Action Upload -PackageDir $package -TaskName $run
```

```powershell
& $tool -Action Start -PackageDir $package -TaskName $run
```

```powershell
& $tool -Action Status -PackageDir $package -TaskName $run
```

终态后回收。默认在准备包旁创建带时间戳的目录。已有任务目录或输出时拒绝覆盖。

```powershell
& $tool -Action Fetch -PackageDir $package -TaskName $run
```

`SCALE_PILOT_COMPLETE_REVIEW_REQUIRED` 以及两档独立封存表示可开始保存审查。回收 `measured_report.md`、`resource_summary.json`、`timing_summary.csv`、逐阶段通信与 RSS/CPU 原始样本。`render_d1b_scale_saved(retrieved_dir, new_figure_dir)` 从 Fetch 后的封闭资源表生成可编辑 FIG、PNG 和重开图数据检查。代码把稠密直接解占比同时按整档含写盘/独立审计和数值阶段计时报告，用理想线性强缩放估计 2→4 worker 收益，并分别从 4096/8192 对 22400 作筛查外推。两个尺寸各一次测量；外推范围反映锚点差异。RSS 拟合加 50% 裕量用于初筛，仍需核对大规模 LU/临时内存、30% RAM 储备、各层预算和完整科学门后才可申请 22400。

## 本地证据

冻结包包含 22 个输入清单项和 manifest。PowerShell 验证相对与绝对目录，覆盖 PowerShell 当前目录与 .NET 进程工作目录不一致的情况。MATLAB MCP 已做只读 dry-run、合成 128 节点与冻结原型的矩阵/密度/接收/中点逐位对照和合成保存审计。离线测试覆盖 RAM、RSS、时限、磁盘、输出配额以及缺少、修改或失败的阶段封存。资源报告与曲线只用 `TEST_ONLY` 数据校验生成链。真实时间、RSS、通信、稠密解占比和 22400 实测基线均等待用户手动启动与回收。

[冻结包](../results/validation/d1b_scale/preparation_20261119_04/input_package.tar) SHA-256 `3481159271c545cce01091f6f689e1b3de17cd78db405f98a33f75b3ee646802`。输入清单 SHA-256 `2e7a069a1d63089807a3999ec59b47de1634e5a9f16f227c5dd5393be9ed8aa2`。 [本地总检查点](../results/validation/d1b_scale/local_validation_20261119_01/checkpoint.json)
