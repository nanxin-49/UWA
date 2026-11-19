# K1首跑中断只读诊断

2026-11-18，北京时间17:48:37检查。仅按`remote-matlab-ssh`技能读取用户手动启动的既有任务，不启动MATLAB、不重试、不发信号、不修改远端文件。源总宽0.10m、4kHz、seed260001、原C2和411点，最终05包。

## 结论与证据

会话`20261118T093555_1124145`实际累计611.424s，控制器退出78，状态`K1_FAILED_STOPPED_NO_FURTHER_DISPATCH`。首例U6/N7退出0，五层BIE及最终独立审计封存链由冻结文件验证器重新核验为`AUDITED_CASE_PASS`。这是远端既有审计/文件SHA链核验；未在本地重新计算MAT/CSV/FIG数值，尚不完成K1科学审查。

第二例U6/N40的MATLAB约18.833s后退出-9。其stderr明确记录：

```text
double free or corruption (out)

MATLAB is exiting because of fatal error
```

这是原生内存管理致命错误的直接证据；单凭-9不能判为OOM、关闭终端、用户误杀或某个特定函数错误。尚无崩溃栈定位具体原生组件，不能认定是PM/N40物理参数引起。`stop_reason=null`、冻结控制器的停止逻辑及资源记录表明本包未记录资源/时限主动停止。最后18次采样的最低MemAvailable=118,994,448,384 bytes（110.822GiB），进程组RSS峰=1,378,869,248 bytes（1.284GiB），memory PSI=0、io最高0.25%；远低于12GiB上限，远高于RAM预留。1秒采样不能排除未记录的瞬时事件或系统限制，未查内核日志，不宣称彻底排除所有系统OOM。

U6/N40最后阶段为MODELS；已有`model_stage/solve_invocation.json`一次调用标记，但无`solved.json`、`committed.json`或`common.mat`，没有开始该海面的五层BIE。四个其余新海面未派发。任务记录的最后子进程身份已不存活，不存在应自动恢复监控的任务。

U6/N7墙钟采样约543.644s，整组RSS峰3,960,299,520 bytes（3.688GiB），最低MemAvailable=117,401,104,384 bytes（109.338GiB）；其成本是正式新案例实测，不能用未完成的U6/N40推断剩余总成本。检查时系统可用约111.734GiB、SwapFree约3.724GiB、PSI为0；当前状态不能追溯覆盖此前采样。

原始只读证据、逐文件SHA、完整采样CSV及日志内容保存在[SSH快照](../results/validation/k1_interruption_review/20261118_01/ssh_stdout.json)，[短回执](../results/validation/k1_interruption_review/20261118_01/receipt.json)。根因状态为**原生致命退出已观察，具体触发原因未定位**。

## 下一步边界

先用既有工具手动Fetch，保留全部原件；不能删除调用标记、改数值门或直接再Start。原代码`local_mark_invocation`会拒绝未封存的旧调用，直接Start不会安全地从MODELS内部续算，反而增加新的失败会话及累计账时。

```powershell
$package = 'E:\MISC\CARPE3D_matlab\Explain\results\comparisons\pm_wind_cutoff_k1\preparation_20261118_05'
$tool = Join-Path $package 'code\manage_k1_pm.ps1'
& $tool -Action Fetch -PackageDir $package
```

回收后可独立检查任务目录内崩溃证据/调用链并设计非声学诊断；任何涉及重算U6/N40未封存MODELS阶段的恢复必须显式保留本次失败证据、重新登记独立尝试和授权，不能冒充自动检查点恢复。已封存U6/N7不应重算；原R1/R2/K0和最终05包不修改。本轮未启动求解、未降低数值门、未改生产接口、未提交推送，不开展其他频率/Seed/通信/Bellhop。
