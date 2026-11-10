# 新 Ubuntu 最小项目环境与兼容性基线（2026-11-10）

## 结论与边界

本轮状态为 **`NEW_UBUNTU_MINIMAL_PROJECT_BASELINE_PASS`**。结论只覆盖隔离快照中的 CPU double、128×128 小型公共信道案例，以及 P22 已保存共同窗和已知 CSI QPSK 结果的原数值审计。九点频率循环只检查接口与跨平台一致性，不认证 CIR。短程案例关闭 `enforce_1_over_R`，不建立 1/R 或目标深海传播资格；正式 S5/S6 继续暂停。两频点旧最终审计代码来源缺口继续保留。

本地项目保持原工作树、分支和已有文件。新机仅使用 `ssh -o BatchMode=yes matlab-ubuntu-new`，未访问旧别名。远端授权根 `/home/misc/cwx` 不是完整 Git 工作区，因此本轮只新建 `/home/misc/cwx/project_baselines/20261110_121507_new_ubuntu_minimal/project/` 作为隔离快照；没有初始化、整理或覆盖根目录中原有文件。远端运行结果保存在该快照的 `results/validation/new_ubuntu_project_baseline/20261110_121507_new_ubuntu_minimal/`。

## 来源、环境与文件身份

- 本地 HEAD：`747af7ee88abda0bd7f554ddbd940ad3bc5f8cb9`。白名单共 23 个文件：公共入口和 setup、所需 `src/` 依赖、两项 P22 独立审计、本轮新兼容性入口，以及 8 个冻结 P22 MAT/CSV/JSON 输入。逐文件路径、SHA-256、Git 状态和类型见[来源清单](../results/validation/new_ubuntu_project_baseline/20261110_121507_new_ubuntu_minimal/project/results/validation/new_ubuntu_project_baseline/20261110_121507_new_ubuntu_minimal/source_manifest.json)，未提交工作树状态见同目录 `source_worktree_status.txt`。两项恢复审计的历史原文、当前独立输出路由补丁分别见 `audit_common_route.diff`、`audit_qpsk_route.diff`；新入口为 `scripts/validation/run_new_ubuntu_minimal_channel_baseline.m`。没有把远端原有散文件当作输入。
- MATLAB Dependency Analyzer 对三项入口解析出的最小 MATLAB 文件依赖已纳入白名单；所查 MAT 顶层内容为数值、逻辑和结构体。原输入 8/8 SHA 与旧检查点一致。上传后及运行后白名单 SHA 均为 **23/23**，缺失或不符为 0。保留所需目录结构，未通过 `cd` 修复入口路径。
- 本机使用现有 MATLAB MCP 会话：**R2025b Update 6**；远端显式调用 `/home/misc/MATLAB/R2025b/bin/matlab`：**R2025b Update 1**。远端 `-singleCompThread` 实测 `maxNumCompThreads=1`，并设置单计算线程环境变量；无 GPU、parpool、软件安装和系统配置变更。新入口 Code Analyzer 问题为 0。
- 每个远端数值作业启动前检查物理 RAM、负载、目标磁盘与 MATLAB 进程；按至少 3 GiB 增量和保留 30% 物理 RAM 余量判门，磁盘要求留 5 GiB。五个数值作业均过门并串行执行。完整逐作业预检、启动时间、PID、stdout/stderr、真实退出码、`/usr/bin/time -v`、耗时和峰值 RSS 在[回传完整日志目录](../results/validation/new_ubuntu_project_baseline/20261110_121507_new_ubuntu_minimal/project/results/validation/new_ubuntu_project_baseline/20261110_121507_new_ubuntu_minimal/retrieved_full/20261110_121507_new_ubuntu_minimal/jobs/)及[最终检查点](../results/validation/new_ubuntu_project_baseline/20261110_121507_new_ubuntu_minimal/project/results/validation/new_ubuntu_project_baseline/20261110_121507_new_ubuntu_minimal/checkpoint_final.json)。首个单线程探测的包装 PID 写入有时序问题，shell PID、启动时间、退出码、耗时、RSS 与完整输出均已保留；没有重试该探测或重复启动数值作业。

## 数值验收

公共案例固定均匀 `c0=1500 m/s`、`nx=ny=128`、`xw=yw=25 m`、`z_tx=z_max=10 m`、`z_rx=1 m`、`sigma_src_m=0.4`、`stepz_lamb=1`、平面海面 `sea_hs_target=0`、`sea_seed=12345`、`surface_boundary_model='kirchhoff_spatial'`、`channel_phase_reference='direct_dsp'`、无气泡、CPU double、`show_figures=false`、`save_mode='rx_only'`、`enforce_1_over_R=false`。完整 `paramsV` 和原生 `output` 分别保存于每例本机及远端 MAT，JSON 回执记录配置与内部检查。

| 案例 | 远端结果 | 本机与远端逐元素复数比较 | 其他检查 |
| --- | --- | --- | --- |
| 4 kHz、无海面反射 | 退出 0，30.64 s | 14 个分量最大绝对差 0 | 公共字段、形状、有限值、频轴、参考索引与分量恒等式通过 |
| 4 kHz、有海面反射 | 退出 0，18.16 s | 19 个分量最大绝对差 0 | 同上；最大本机分量闭合误差 `2.7756e-17` |
| 有反射、4–8 kHz 九点、参考 6 kHz | 退出 0，20.26 s | 19 个分量最大绝对差 0 | 频轴及索引 5 一致；最大本机分量闭合误差 `6.2063e-17` |

逐元素门为 `abs(remote-local) <= 1e-12 + 1e-9*abs(local)`；未拟合幅相、移峰或放宽容差。配置语义、频率轴及参考索引均一致，`H_f=H_direct_f+H_reflect_f` 和参考频点的 `h_total` 等既有接口恒等式通过。各案例详见结果包中的 `comparison_direct_4k.json`、`comparison_reflect_4k.json` 和 `comparison_reflect_9f.json`。

P22 两项均只读取指定保存结果，审计输出显式指向本轮独立目录；没有调用生成新通信实验的入口或写回原数据包：

| 保存结果审计 | 原门与本轮结果 | 跨平台复核 |
| --- | --- | --- |
| 共同窗 20261030_01，含原宽带与 257 点输入 | 原 `1e-3` 门通过；本机最大相对差 `8.2248e-6`，原始 65/129 点失败门仍保留 | 指标最大差 `3.5406e-19` |
| 最终 QPSK 20261030_03 | 原 `1e-12` 重放及 CSV 门通过；远端重放差 `1.4218e-15`、CSV 差 `3.9080e-14`；原配置、种子、误比特计数保持 | 重放指标差 `1.4218e-15`；CSV 指标差 0 |

两项远端审计分别耗时 50.70 s、35.50 s，真实退出码均为 0；其原数值断言及生成的 JSON/FIG/PNG 已检查。P22 结论仍限于 PE 相位屏、共同窗表示及已知 CSI 控制通信，不能升级为正式多模型 S5/S6 或物理 CIR/BER 资格。

## 资源、回传与接续

本机新增 MATLAB 工作连同进度看板刷新保守计 **110/180 s**。远端全部 7 次 MATLAB 进程（含两个环境探测）累计 **204.02/600 s**，最长单次 **50.70/180 s**，最大峰值 RSS **1,573,464 kB**。输入上传及运行后 SHA 23/23；远端结果整体回传清单自身 SHA 为 `1c18452c9fa515d217ddb48807810703717db532c1ead02d22f6e0e6918f255d`，回传 **92/92** 文件 SHA 相符。跨平台比较的是保存数值数组，不要求新生成 MAT 文件逐字节相同。全量结果、原始日志、比对 JSON、源清单、阶段与最终检查点均在[本轮本地结果包](../results/validation/new_ubuntu_project_baseline/20261110_121507_new_ubuntu_minimal/project/results/validation/new_ubuntu_project_baseline/20261110_121507_new_ubuntu_minimal/)。

新机器现可在这份隔离快照和同等资源门下承担**已验证的小型 CPU double 公共信道运行与 P22 保存结果审计**。尚未检查完整项目在 Linux/R2025b 的所有 toolbox、外部可执行文件、Linux MEX、大小写路径、未纳入白名单的代码与数据、较大网格及长时/并行作业。完整迁移前须另行核定目标工作区归属、完整依赖与输入身份，设计新增案例及资源预算，并按其自身数值基线验收。本轮停止，不启动进一步重计算或完整迁移。
