# K1远程资源评估、执行决策与预算

2026-11-18。**准备完成，未启动K1。** 仅通过matlab-ubuntu-new做了资源快照和一次有界非声学MATLAB诊断；无parpool、无BIE/PS/FullFK新求解，无Git提交推送。

## 12:21补查：原生基础设施通过，资源准入恢复

用户继续授权测试及远程检查后，北京时间2026-11-18 12:21:05–12:21:22在独立`/home/misc/cwx/k1_native_probe_20261118_07`完成有界非声学诊断。可用RAM38.504→39.018GiB，SwapFree仍仅约0.14–0.17MiB；memory full avg10=0%，io full avg10=0.10–0.33%，load1约3.24–3.31，CPU忙约2.65–2.85%，磁盘约4836.97GiB可用。冻结准入函数当时返回PASS；没有启动MATLAB或声学计算。以下10:41快照是历史状态，用户正式Start仍重新检查。

仍选择A：**1进程×1计算线程**。两进程24GiB加30%当前基准约11.7GiB后仅余约3.3GiB，Swap近耗尽且额外并发许可证未知；短时PSI下降不足以支持可靠并行。单例12GiB加约11.7GiB系统预留后尚余约15.3GiB。65–85min条件预算、1800s/例及10800s累计和全部原门保持。

11项原生检查15.646s、0失败：fd9经bash/nohup/exec继承、独立flock竞争拒绝/退出释放、会话RSS、错误身份拒绝、真实TERM及15s后KILL、实际CLI无租约拒绝等。四个最多24秒的轻量诊断子进程全部回收，其他进程未收到信号。使用独立诊断锁，未触碰正式共享K1锁，未调用正式Start。[原生回执](../results/validation/k1_preparation/native_20261118_07/retrieved/receipt.json)、[前快照](../results/validation/k1_preparation/native_20261118_07/retrieved/host_before.json)、[后快照](../results/validation/k1_preparation/native_20261118_07/retrieved/host_after.json)。35项回收SHA与外层归档SHA均通过；05包未变。环境/许可证已有非声学诊断证据未变，本次未重复占用MATLAB许可证。

## 实际主机与许可证检查

检查快照截至北京时间2026-11-18 10:41:45；这是当时状态，手动启动会重新准入，不保证此后资源不变。

| 项目 | 实测 |
|---|---|
|CPU|Intel Xeon Gold5120，2插槽×14物理核=28物理核，56逻辑线程；x86_64；2个NUMA节点；affinity 0–55|
|CPU负载|load约3.4；两个1秒快照忙碌约2.98%/1.08%，CPU不是主要准入问题|
|RAM|125.335GiB总量；MemAvailable约44.94–45.01GiB|
|Swap|总约8GiB，空闲由约0.234MiB降到0|
|压力PSI|memory full avg10=19.45%→27.49%；io full avg10=34.49%→38.24%|
|其他进程|存在6条MATLAB和1条matlab_helper进程记录；ANSYS进程RSS约65GiB。只读资源列，没有访问其任务文件/参数或终止进程|
|磁盘|/home/misc/cwx所在卷可用约4838GiB，空闲inode约2.186亿|
|限制|nofile软1024/硬1048576；nproc477341；地址空间/CPU时间无限；栈软8MiB。cgroup v2 cpu.max/memory.max未暴露，不据此宣称无cgroup限制|
|MATLAB|显式/home/misc/MATLAB/R2025b/bin/matlab；R2025b Update1、GLNXA64；Intel oneMKL2024.1；默认计算线程上限28|
|PCT|ver仅列MATLAB，PCT未安装。license(test,PCT)虽为true，不能等同安装、checkout成功或worker可用；未创建pool，未尝试未安装产品checkout|
|基本许可证|一次额外前台MATLAB非声学诊断成功，32.918s、实际退出0。许可证总席位/额外多进程并发上限未知，未从test或历史并发推断|

证据：[第一次成功资源快照](../results/validation/k1_preparation/resource_20261118_03/host_stdout.json)、[MATLAB诊断](../results/validation/k1_preparation/resource_20261118_04/matlab_diagnostic.json)、[诊断后快照](../results/validation/k1_preparation/resource_20261118_05/host_stdout.json)。初始沙箱SSH别名解析失败和探针语法失败均保留；使用原用户SSH配置后成功，没有改用旧主机。诊断日志写在独立k1_resource_probe目录，不修改共享MATLAB偏好设置。

## 历史BIE成本

从R1/R2既有timings和资源CSV提取15个阶段；没有使用K0 RSS预测BIE。

|阶段|节点数|R1秒|R2 260002秒|R2 260003秒|
|---|---:|---:|---:|---:|
|C0|3200|55.481928*|60.230084|56.596749|
|C1|3840|80.841146|83.857962|79.176863|
|C2|3840|65.133454|68.352375|66.782072|
|C3|4480|108.641842|112.945933|108.395258|
|C4|4480|106.646117|112.621083|108.542554|
|五层合计||416.744487|438.007437|419.493496|

*R1 C0为导入的函数计时；该份timings没有分项，不能把未知分项填零。R1回收会话4:02.95只含C4续行与收尾，不是完整五层墙钟。

R2两个独立单线程任务实际墙钟610.78s/562.07s，历史报告记录重叠约558s；260002原图件失败保留，不把其成本记录当数值PASS。最大4480节点时，R2 C3基础组装约29.5–30.0s、近场替换约6.2–6.5s、稠密求解约8.1–8.9s、边界残差约56.4–59.3s。完整两Seed中组装/近场/边界合计约84–85%，稠密求解仅约7–8%。单独加速LU不能保证显著总收益，其他多线程收益本轮未实测。

BIE阶段采样进程组RSS峰约2.154GiB；含models/绘图收尾的整会话峰约3.606GiB。time-v单进程峰约1.89/1.92GiB，和会话RSS求和口径不同，不能互换。4480²个complex double单矩阵约0.299GiB；另有矩阵副本、LU、求积缓存、场与MATLAB开销。矩阵/组装大致二次、LU大致三次规模关系，不对总程序盲目作三次外推。K1最大节点仍4480，N40改变海面物理带宽，不改变这个节点上限。

采用**12GiB/案例**保守工作估计，约为历史整组峰的3.3倍。[15阶段CSV](../results/validation/k1_preparation/historical_20261118_02/BIE_stage_costs.csv)、[历史资源回执及来源SHA](../results/validation/k1_preparation/historical_20261118_02/receipt.json)。

## A/B/C比较与最终选择

|方案|有利负载下的规划耗时|工作内存|许可证与可重复性|决定|
|---|---|---|---|---|
|A：一个进程、单计算线程、六例顺序|约65–85min，含启动/独立审计/图件余量|估计12GiB；单个活动MATLAB会话|一次额外基本MATLAB已验证；保持R1/R2计算线程路径|**选择**|
|B：两例独立进程并行，每例单线程|理想约35–50min；当前内存/I/O压力下不能保证|估计24GiB，另留30%基准约13.5GiB|历史两进程有效，但当前额外两席未验证，总并发限制未知|不选：压力明显，不能以名义可用RAM推断可靠并行|
|C：一个进程多线程BLAS/PCT|没有合格实测总时间；只改善7–8%的LU部分不能保证总体缩短|仍需至少单例预算，内部线程/临时内存需新评估|BLAS默认28线程可用，但声学等价性未重验；PCT未安装|不选：收益未证、资源和复现路径不占优|

选择是基于本次检查，不是预设单线程。上限固定为**1个MATLAB进程×1计算线程**；设置-singleCompThread及OMP/MKL/OpenBLAS=1，不创建parpool，不并行C0–C4，不进行未经验证的NUMA绑定。

按约45GiB启动可用量，单例12GiB和30%预留约13.5GiB后仍有约19.5GiB未分配余量；两例名义算术虽可容纳，实际PSI/Swap证据不足以支持可靠并行。检查时连正式单例准入的压力门也不通过；包不会强行启动或自动排队等待。

## 冻结预算与安全门

- 单例墙钟1800s；累计10800s（3h，含预检、恢复和汇总）；汇总最多600s，任一先到即停止。65–85min是条件估计，不是保证或K1实测。
- 12GiB案例工作估计；活动MATLAB进程组RSS上限12GiB。小型Python控制器及系统开销由额外余量覆盖；MemAvailable全局门包含所有进程影响。
- 保留>=30%启动可用RAM；恢复使用历次基准最大值，不重置降低。派发前要求MemAvailable减预留仍>=12GiB。
- 启动要求memory full avg10<=5%、io full avg10<=10%、load1<=21、磁盘>=5GiB、空闲inode>=5000。检查时19–27%/34–38%不满足，未来Start重新取样。
- 每秒检查本任务会话RSS、MemAvailable和PSI。RAM/RSS/时间超门立即停止自己的已验证会话；memory full>10%或io full>30%连续10次采样停止。不改其他进程、不清Swap。
- 用Popen创建记录、PID/start_ticks/独立session与进程组验证身份；重用PID或身份不匹配拒绝终止。任一案例失败停止后续派发，无自动调参、重试、降门或转主机。
- 保持原v2非单调诊断和全部数值门；每例PS六层、BIE五层、FullFK四层及411点独立审计，封存后才PASS。模型误差是否超过5倍经验分辨率如实分类，不改阈值迎合结果。

## 矩阵与锚点

seed260001：U6/N12严格历史锚点，另依次U6/N7、U6/N40、U5/N12、U8/N12、U5/N40、U8/N40六个新面。已核源、原系数、C2、phase contract2.0.0、411点及资格SHA；C2差1.11e-16m，R1正式资格SHA为eb6d86a8e845ee2785efab4b93470f91d26026ed0d415bdbd808023152b43265。只有完全匹配才复用锚点，六新面不带入其密度，各新建C0–C4，合计30次BIE层。

原始/C2系数、海面、谱与指标随每案例输入冻结，预检核Parseval与谱差；最终另输出raw_C2_surface_spectral_differences.csv。N40绝不标为完整PM谱已收敛。固定0.10m二维条带、4kHz、原相位、原411真实坐标；不扩其他Seed、宽带、CIR、通信或Bellhop。

最终包与操作见[手动命令](k1_manual_commands_2026-11-18.md)，[测试与SHA](k1_preparation_validation_2026-11-18.md)。原生Linux锁竞争/真实信号已由上述独立诊断补测；完整正式启动器Start及K1数值仍未执行，基础设施PASS不等于声学PASS。
