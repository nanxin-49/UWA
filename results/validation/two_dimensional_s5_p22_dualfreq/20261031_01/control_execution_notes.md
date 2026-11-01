# 2026-10-31 执行控制失败记录

- 许可探针首次启动脚本在 RAM 算术门使用 `awk` 的默认数字格式，输出 `2.6217e+10`，Bash 报 `无效的算术运算符`；脚本继续运行了轻量探针。探针 MATLAB 自身退出 0、12 s，出现 `LICENSE_PROBE_OK` 与 R2023b，未出现账号提示。**不能把该次 RAM 算术门记为通过。** 在任何 BIE 前另以整数格式重查物理 RAM：总量 `33223528448` byte、可用 `26218598400` byte、MATLAB 进程数 0，输入清单 SHA 通过。三层 BIE 各自的 `attempt_20261031_01/preflight.txt` 均使用整数格式，独立审计重算 30% RAM、磁盘、预算和预测 SHA 门。
- 首次运行新独立审计脚本在读取 4 kHz CSV 时调用 `complex(str, str)`，报 `TypeError: complex() can't take second arg if first is a string`，尚未创建审计输出目录。仅将 CSV 两列显式转为 `float` 后重跑；第二次审计成功，日志保存在 `audit_dualfreq_attempt2.log`。首次失败的命令和 traceback 保留在本轮工具执行记录中；本文件记录其原因与最小修复，不伪装成原始 stdout 文件。
- 2026-10-30 原 `freeze` 退出 124、两次三层合并尝试（一次无输出、一次登录失败且退出未观测）均在原始远端任务与旧本地镜像中保留，未覆盖。
