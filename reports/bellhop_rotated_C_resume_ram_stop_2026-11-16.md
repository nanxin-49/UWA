# C组剩余9项恢复前检查：再次RAM停止（2026-11-16）

本轮恢复授权有效，但恢复前SHA检查后再次触发1.5GiB RAM门，因此按用户明确停止条件结束，没有派发Bellhop。C保持30/39项，尚无完整C步长收敛资格。

## 实际执行及证据

工作树 `C:\Users\ASUS\.codex\worktrees\bellhop_rotated_fig12b\Explain`，分支 `experiment/bellhop-rotated-fig12b`；现有共享MATLAB R2025b Update6、PID32280，显式project_path，无新MATLAB/SSH/网络/GitHub。冻结C配置、REAL8及平面墙优化可执行文件均未改动；没有修改生产PE或官方Bellhop。

[恢复停止回执](../results/validation/bellhop_rotated_fig12b/20261115_C_run_01/resume_20261116_01/checkpoint_stop.json)保存本轮数值、预算、未运行项及错误索引。历史root checkpoint/stage/delivery先复制至该恢复目录的history，再更新当前停止状态；2026-11-15原始SHA清单和旧失败输入保留。原空diary同时保留，完整错误补记在preflight_error.txt；不将空diary作为成功证据。

恢复预检严格完成了A277、B1044、C827文件清单SHA，以及A/B78项和C30项component/checkpoint SHA；冻结计划、优化二进制及A/B阶段检查点SHA一致。第12阶下行旧失败输入无PID、SHD、process记录。所有这些断言先于失败的RAM门；没有调用单项派发器。重新计算A/B/C分量均为0。

## 内存诊断与释放边界

| 观测阶段 | 系统可用RAM | MATLAB报告占用 |
|---|---:|---:|
| 初始共享会话 | 2.293716GiB | 2.098831GiB |
| 停止后的独立诊断 | 1.254295GiB | 2.525875GiB |

系统总可见RAM15.245850GiB。初始可用值达标，但该值不是后续启动许可；SHA复核后的检查低于1.5GiB，触发 `prepare_bellhop_rotated_C_remaining` 第23行 `Resume RAM launch gate`。

失败瞬间的精确可用数值未被该断言持久化，明确记为unknown。后续1.254295GiB是停止后的采样，不冒充失败瞬间数据。CIM查询不可访问，早期shell产生的0占位无效，未用于资源判断。

现有会话未发现可确认属于本任务的大型临时工作区数组，因此没有释放任何不明所有权变量、图窗或函数缓存，没有执行clear all、全局清理、垃圾回收请求或终止进程。MATLAB自身占用报告上升约0.427GiB；shell可见其工作集约0.91→1.24GiB，同时存在浏览器、微信、Codex/ChatGPT及桌面等应用。两种MATLAB占用口径不同，不与系统可用量直接相加解释。

读取实际 `black_file_hash.m` 确认其用1MiB分块流式SHA，而非整文件加载。证据支持共享系统内存压力及应用/运行时驻留增长，尚不能识别单一分配对象或证明内存泄漏。没有可安全释放、且预期能解决缺口的本任务大型数组，故保持共享会话和其他应用不动。停止后即使资源回升也不自动重试。

## 原预算与资源

- 上轮最终累计2580.834401秒，余额1019.165599秒，没有重置。
- 本轮SHA预检/诊断保守计120秒，停止交付预留40秒，含静态检查与看板。不是实测Bellhop耗时；新增Bellhop耗时为0。
- 最终累计2740.834401秒，原预算余额 **859.165599秒**。
- 累计输出10.630438GiB，小于15GiB；磁盘空闲132.446716GiB，大于20GiB。磁盘门通过，RAM门失败，不以预算/磁盘充足替代RAM要求。

预算为保守计账，含未使用收尾余量，不擅自退还。此前C30项实际native270.728181秒保持历史数据。

## 未执行与下一步

仍缺第12阶下行和上下第13–16阶，共9项；无PID、无新声压。本轮未运行完整C组装、N12→N16、B→C复L2/最大复差/包络差/有效相位差，也未生成A/B/C完整场图册、PNG或FIG，不能报告上述门通过。既有A/B图册完整保留。

用户授权的恢复入口已新增：`prepare_bellhop_rotated_C_remaining`只复核与记账；`run_bellhop_rotated_C_remaining_next`计划对原无PID任务使用独立后缀目录，并在启动前核新旧输入SHA完全一致。后者本轮未调用；其存在不等于求解资格。没有更改物理参数、接收坐标、步长、束数或任何数值/资源阈值。

本轮到此停止，等待人工审查。若未来重新授权，必须再次核RAM、剩余859.17秒能否覆盖9项和全部审计、磁盘与SHA；30个已完成项直接复用。不得自动反复重试，不启动D或新FFP，不提交、推送或合并。

已更新受影响tracker、PROJECT_CONTEXT、方法计划、README导航与进度源；看板只读渲染不运行声场。未改变生产信道或通信方法，vertical_comm_guide不需修改。静态检查不替代尚未运行的完整C验证。
