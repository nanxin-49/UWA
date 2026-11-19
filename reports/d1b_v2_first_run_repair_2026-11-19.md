# D1b-v2 首跑失败与通信统计修复

2026-11-19。用户手动启动的 `d1b_v2_20261119_01` 在 2-worker 第一轮 L06 的 BIE 计时后读取传输统计时退出。MATLAB 退出 `1`，监督器退出 `78`，耗时 `191.731 s`，采样峰值任务总 RSS `3,773,747,200 bytes`（3.51 GiB），最小可用 RAM `128,528,973,824 bytes`。回执未触发资源停止门。远端 R2025b Update 1 确认 PCT 25.2 已安装、checkout 成功、默认 `Processes` profile，另有 `Threads` profile。已完成海面预检和两轮串行共八个基准阶段，2 个 Process worker 经实际登记且各限制为 1 个计算线程。2-worker 阶段结果未封存，4-worker 未到达。 [MATLAB 日志](../results/validation/d1b_v2/first_failure_20261119_01/matlab.log) · [检查点](../results/validation/d1b_v2/first_failure_20261119_01/checkpoint.json) · [监督器退出回执](../results/validation/d1b_v2/first_failure_20261119_01/exit.json)

原 `local_bytes_end` 将 `tocBytes(pool,token)` 的输出按 table 访问 `Properties`，触发 `MATLAB:structRefFromNonStruct`。R2025b 带输出的调用返回 `numWorkers×2` 数值矩阵，列为发往 worker 与从 worker 收到的字节。修复按该格式保留逐 worker 字节数和列名，校验矩阵形状、非负及有限值，并在 pool 登记时先检查一次。科学计算、配置、PCT worker 数与资源门沿用首包。 [MathWorks R2025b `tocBytes`](https://www.mathworks.com/help/releases/R2025b/parallel-computing/parallel.pool.tocbytes.html) · [修复入口](../scripts/validation/run_d1b_resource_prototype.m)

只读取回六份首跑文件并与远端 SHA 核对。新包 `preparation_20261119_02` 的 19 个清单项中只有 MATLAB 原型入口内容发生变化，原 `01` 包 SHA 复核一致。绝对路径的 PowerShell DryRun 通过。现有本地 MATLAB 验证了模拟的 2×2、4×2 通信输出、串行空值和四种无效统计，Code Analyzer 零消息。当地测试没有开进程池，修复版本的实际远端完整基准仍待用户手动运行。 [离线回执](../results/validation/d1b_v2/repair_validation_20261119_02/offline_receipt.json) · [MATLAB 测试回执](../results/validation/d1b_v2/repair_validation_20261119_02/matlab_receipt.json)

新任务名 `d1b_v2_20261119_02`。首跑目录和冻结包保留。新包按同一次运行采集串行、2、4 worker 时序，便于比较同一会话的性能。完整 D1b 声场仍待准入。 [手动命令](d1b_v2_manual_commands_2026-11-19.md)
