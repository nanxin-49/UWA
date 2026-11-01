# 二维 S5 P22 双频物理资格先导：授权前检查（2026-10-29）

## 状态与范围

本轮截至阶段 3，**没有同步文件、没有新远端 MATLAB 进程、没有新声学求解**。阶段 0 `PASS`（接续及源身份核验）；阶段 1 `BLOCKED`（远端工作根不是项目检出且共享根分区仅余 5,640,626,176 byte，`df` 显示 100% 使用）；阶段 2 `PASS`（事前冻结试验及预测规则，尚无待测频点预测数值）；阶段 3 `BLOCKED`（环境资格未通过且目录、同步、1800 s 预算未获明确批准）。阶段 4–5 `NOT_RUN`。S4=`CORE_COMPLETE` 仅原 4 kHz 范围；原 S5 A=`PARTIAL`、B=`BLOCKED`、S6 C=`BLOCKED`。

工作树起始为 `3D_try`，HEAD `a30655b921030a35e73613410b65e4f245f1ee78`，相对 `origin/3D_try` ahead 2，起始工作树干净。上一轮阶段 0–4 的限定轻量验收见 [回执](two_dimensional_s5_s6_preflight_2026-10-29.md)；本轮没有重复其九点解析 CIR 或解析 QPSK。4 kHz 参考见 [S2b 正式包](../results/comparisons/two_dimensional_channel_closeout/20261023_s2b_1305/formal_reference_v2/README.md)，其中弱正弦 P22 BIE 原生反射复场为 `0.116633796622219 + 0.032179471122302503i`，平面为 `0.116630666757294 + 0.032190841211023599i`。这是旧证据的读取，不是本轮重新求解。

共享 MATLAB MCP R2025b 中实际重开旧 S2a `numeric_inputs.mat`、S2c `frozen_source.mat` 和 S3a 冻结源，计算 `max|fft(S2a Tx)-frozen.coeff|=0`、`max|S2a kx-frozen.kx|=0`、`max|S2a Tx-ifft(frozen.coeff)|=2.2204460492503131e-16`，S2c/S3a `coeff/kx/kz` 逐项最大差均为 0。旧弱正弦参考及代码记录的是 Gaussian σ=0.3 m、N=984、宽 192.1875 m、Tx 深 100 m、c=1500 m/s、4 kHz、原生 `exp(-iωt)`、A=0.01 m、K=0.1 rad/m、C2 内 42 m/外 50 m。源定义具有可核对的同一性；旧 4 kHz 结果仍仅具其原有经验数值参考范围。

## 远端只读资源回执

经 `ssh -G` 确认 `matlab-ubuntu` 配置用户 `misc`，BatchMode 返回 `SSH_OK` 和 `id -un=misc`。`readlink -f /home/misc/cwx` 仍为 `/home/misc/cwx`；此根目录是多个任务文件夹的集合，不含 `.git`、`setup_vertical_project.m` 或本项目完整检出。已知旧任务 `/home/misc/cwx/uwa2d_s5a_20261027` 有独立 `code/`、`input/`、`jobs/`，其 `input/frozen_source.mat` SHA-256 与本地相同：`8d60d3d9d6aa3aa16110ac33d81e0d193a61526aef2c0c2e9edd531cae531f55`。拟用全新目录 `/home/misc/cwx/uwa2d_s5_p22_dualfreq_20261029_01` 当前不存在。没有读取其他任务内部、清理、覆盖或传输。

远端 `/home/misc/MATLAB2023b/bin/matlab` 可执行；安装元信息为 R2023b 23.2.0.2365128，但**尚未启动 MATLAB 核对运行时版本或工具箱**。物理 RAM 总 33,223,528,448 byte、available 29,858,402,304 byte；12 核，load average 约 `0.03/0.08/0.08`，当时无 MATLAB 进程。`/` 与 `/tmp` 同一分区，`df -B1` 可用 5,640,626,176 byte，显示使用率 100%；inode 空闲 119,673,924。旧任务 `code/` 36,103 byte、`input/` 25,967 byte、`jobs/` 2,790,322 byte，8 kHz 子任务约 3,247,914 byte。旧文件体积不能证明新求解或 MATLAB 临时文件上界。由于共享分区近满且根目录不是完整项目，阶段 1 安全执行资格未通过；必须由用户明确指定/批准独立目录、同步范围并确认足够可用空间后重新只读预检。

当前远端旧任务的七个通用 helper 与本地 SHA 相同；远端旧 BIE helper 为旧版 `76f1a9...`，本地九点接口修复版为 `5ea51e...`。旧固定 6 kHz runner 含硬编码频率及旧平台数值，不能直接作为新双频 runner。未来仅在获批新目录内放置经静态审查的双频验证 runner，不覆盖旧任务。

续行只读复查：拟建目录仍不存在；可用磁盘降至 `5,619,503,104 byte`，available RAM 为 `29,846,220,800 byte`，12 核、load 约 `0.10/0.10/0.09`，仍无 MATLAB 进程。未发现环境门解除。本地新增 [验证专用双频 runner](../scripts/validation/run_two_dimensional_s5_p22_dualfreq_remote.m)，只允许两项冻结频率，在半步 BIE 前要求独立预测文件，逐频保存原生 direct/reflect/total 与共同原点 DSP 表示，且不改生产求解器。共享 MATLAB MCP 的 Code Analyzer 返回 1 条“先前隐藏消息不再产生”的 info、0 error/warning；实际无声学计算的输入拒绝控制 2/2 通过。首次测试 harness 用英文 `Assertion failed` 匹配本地化 MATLAB 错误文本造成 1/2 误报，随后改用 `MATLAB:assertion:failed` 错误 ID 复测通过；原误报在本报告保留。新 runner 的物理分支 **尚未执行**，Code Analyzer 与拒绝控制不构成物理兼容性或数值资格。

再次只读复查：拟建目录仍不存在、根分区可用 `5,605,466,112 byte`，约 29.85 GB 物理 RAM available，MATLAB 进程 0；磁盘门仍未解除。本地另增 [半步预测冻结入口](../scripts/validation/freeze_two_dimensional_s5_p22_halfstep_prediction.m)，只读取经 SHA 核对的 4 kHz CSV 和 4062.5 Hz 五层 BIE 保存结果，要求五层参考敏感性 ≤`1e-4`、实际退出码为 0，并拒绝已有半步 BIE 目录或覆盖已有预测。它将原生复数预测、端点、源/参考/五层 SHA 和公式存成 `prediction.mat/json` 与单独 SHA；双频 runner 在半步 BIE 前复核该 SHA。两个入口 Code Analyzer 均为 0 error/0 warning，分别 3 条和 2 条过时信息提示。

最终只读复查：`/home/misc/cwx` 仍为原任务集合，拟建目录仍不存在；共享根分区可用进一步降至 `5,590,794,240 byte` 且 `df` 仍显示 100% 使用，inode 充足，MATLAB 进程仍为 0。没有来自用户的目录、文件同步或 1800 s 预算批准。资源/授权阻断在连续续行中持续，物理阶段没有启动。

[synthetic fixture](../results/validation/two_dimensional_s5_p22_dualfreq_preflight/20261029_01/synthetic_prediction_fixture/) 使用真实冻结源与旧 4 kHz CSV、**人工构造**的 4062.5 Hz 五层复场，不调用 BIE；预测代数误差为 0、人工五层相对变化 `4e-6`、重复冻结拒绝 1/1，JSON 的独立 SHA 复核一致。其人工 `prediction.json` 只测试软件契约，绝不是新频点物理预测或 BIE 参考；正式半步预测仍须在真实 4062.5 Hz 五层及其独立审计合格后另行冻结。

## 事前冻结的试验与预测规则

物理参数冻结：二维 `Tx=(0,100 m)`、`Rx=P22=(0,3 m)`、`c=1500 m/s`、冻结 984 项 Gaussian 系数、弱正弦 `η(x)=0.01 sin(0.1x)χ_C2(x)`（φ=0，原窗口）、native `exp(-iωt)`；所有分量与频率共用 `tau_ref=96.875/1500 s`。顺序是先 4062.5 Hz 基网格，再 4031.25 Hz 独立半步。每频同源 PE、FullFK 400/800 求积与 BIE 五层：`main96/self64/space14/window60/window70`，L=`60/60/60/60/70` m，ppw=`12/12/14/10/10`，自积分阶数=`96/64/96/96/96`。BIE 经验敏感性须 ≤`1e-4`；模型可分辨须差异 >5 倍经验数值分辨率；不降低门槛。

**预先冻结的 4031.25 Hz BIE 复数预测器：**令 `R(f)` 为同源、原生相位的 P22 弱正弦 BIE *反射*复场，`F(f)=-u_inc(x=0,z=-3;f)` 为相同冻结源的解析无限平面反射复场，`Q(f)=R(f)/F(f)`。当且仅当 4000 Hz 旧 BIE 参考与 4062.5 Hz 新 BIE 参考都合格，先在未读取 4031.25 Hz BIE 结果的情况下计算

`Q_pred(4031.25) = [Q(4000)+Q(4062.5)]/2`, `R_pred(4031.25)=F(4031.25) Q_pred(4031.25)`。

只对复数比值做线性半步预测，不拟合源幅度、整体相位或路径时移；`F` 在各频由相同 984 项源和更新的 `kz` 独立解析计算，须检查 `|F|>0` 且有限。正式执行时先保存预测原生复数、两个端点值、`F/Q`、代码与配置 SHA、时间戳及只读审计，再允许启动 4031.25 Hz BIE。独立参考取得后，分别报告 `|R_pred-R_BIE|/|R_BIE|`（预测）、五层变化/主参考（数值敏感性）、`|R_PE/FK-R_BIE|/|R_BIE|`（模型差异）。若端点资格或源身份失败，不生成伪预测 PASS，也不运行依赖的半步 BIE。

物理原生 direct/reflect/total 必须逐频满足 `H_total=H_direct+H_reflect`。DSP 记录用 `H_DSP(f)=conj(H_native(f))*exp(i 2π f tau_ref)`，并做逆映射误差检查。此批两个频点不足以生成正式 CIR、PDP、BER 或频率重建 `L2≤1e-3` 的全带资格。

## 资源申请，尚未批准

计划新目录：`/home/misc/cwx/uwa2d_s5_p22_dualfreq_20261029_01`，独立 `input/`、`code/`、`jobs/`、`audit/`。拟只同步本地冻结源 `results/comparisons/two_dimensional_channel_closeout/20261026_s3a_seed260003_01/input/frozen_source.mat`，4 kHz `formal_reference_v2/fields_core.csv`，以及下列验证源码的任务快照：`closeout2d_source_at_frequency.m`、`closeout_s2b_fullfk_analytic_exterior.m`、`closeout_s2b_surface.m`、`fk_single_receiver_incident.m`、`run_pe_1d_surface_reflection_validation.m`、`run_pe_physical_surface_model0_validation.m`、`solve_helmholtz_bie_halfplane_vertical.m`、`run_two_dimensional_s5a_weak_bie_remote.m`（仅 6 kHz 兼容性 smoke）、`run_two_dimensional_s5_p22_dualfreq_remote.m`、`freeze_two_dimensional_s5_p22_halfstep_prediction.m`，另有待准备的串行 launcher。已存在的 7 个 helper 共约 22.3 KiB，冻结源 18.8 KiB，旧参考 CSV 49.3 KiB；新增 runner/launcher、日志、MAT 与 MATLAB 临时空间尚需上界核验。**此清单是申请，不是传输授权。**

历史成本投影：两频十层 BIE 核时间约 `226.98+232.06=459.04 s`，若十次启动各计 40 s 则约 `859 s`，另需 PE/FullFK、兼容性、审计和失败余量。预计最大 BIE 增量内存约 `1.460 GB`，仅为公式估计。申请远端新 MATLAB 进程累计 ≤`1800 s`，每进程严格 `<900 s`，串行、零并行；每次作业前至少 30% available 物理 RAM 余量并复查 CPU、磁盘、`/tmp`、已有进程、剩余预算。磁盘目前无法给出安全临时上界，**不以旧结果的数 MB 代替可用空间保障**。共享磁盘空间得到处理且重新测量后，方可提出最终可执行预算。任何源 SHA、数值参考、内存、磁盘、单作业或累计预算门失败即停并保留原输入、日志、实际退出码和部分结果。

## 下一步

请先确认共享远端有足够磁盘空间，且允许将本任务置于上述全新目录；之后重新只读预检并提交最终同步 SHA 清单与 1800 s 计算申请。未获三个独立明确授权（目录、必要同步、物理预算）前，阶段 4 不启动。本轮没有修改求解器、没有覆盖旧结果、没有提交或推送。

## 2026-10-30 远端恢复后复核与正式授权申请

以下是新日期的续行；前述 2026-10-29 磁盘阻断记录保留。用户表示远端已恢复，**尚未批准**目录、同步或物理计算。SSH `BatchMode` 用户为 `misc`，`/home/misc/cwx` canonical 路径不变。新目录 `/home/misc/cwx/uwa2d_s5_p22_dualfreq_20261029_01` 仍不存在；根目录仍是既有任务集合而非 Git 检出，因此只申请在该根下创建一个独立任务目录，不初始化、清理、覆盖或同步到根目录及旧任务。当前本地 `3D_try` HEAD `a30655b921030a35e73613410b65e4f245f1ee78`、ahead 2；此前未提交文件保留。新代码用独立 SHA 清单标识，不冒充已提交版本。

实时只读资源：`/` 与 `/tmp` 同一分区，可用 `623,059,312,640 byte`、已用 67%；inode 可用 `119,662,043`。物理 RAM total `33,223,528,448 byte`、available `29,761,282,048 byte`；12 核、load `0.12/0.10/0.09`、MATLAB 进程 0。`/home/misc/MATLAB2023b/bin/matlab` 可执行，安装元信息 23.2.0.2365128（R2023b）；运行时版本/工具箱及新代码兼容性仍须获批后用小型进程实际检查。旧任务冻结源 SHA `8d60d3d9d6aa3aa16110ac33d81e0d193a61526aef2c0c2e9edd531cae531f55` 与本地一致。

拟同步到新目录 `input/` 和 `code/` 的**全部 13 个文件、目标相对路径、字节数及逐项 SHA**见 [2026-10-30 同步清单](../results/validation/two_dimensional_s5_p22_dualfreq_preflight/20261030_01/sync_manifest.json)，总计 `118,246 byte`，13/13 本地哈希复核通过、目标无重名；[拟用远端 `sha256sum -c` 文件](../results/validation/two_dimensional_s5_p22_dualfreq_preflight/20261030_01/SHA256SUMS_proposed) SHA 为 `0ecbbeedac27539ed149258849487c2da283cbe5d5dd3d6d0c75c6c0606cbfd8`。包含冻结源、4 kHz 参考 CSV、七个通用验证 helper、旧 6 kHz 兼容性 smoke runner、新双频 runner、半步预测冻结入口及[串行 launcher](../scripts/validation/support/launch_two_dimensional_s5_p22_dualfreq.sh)。launcher 在本地 Git Bash `bash -n` 通过；原 WSL `bash.exe` 因宿主访问拒绝未运行，其失败不是脚本语法结果。所有 SHA 必须在正式传输前和传输后重新核对；任何文件变化需更新申请范围。

旧 6/8 kHz 作业目录均只有数 MB，不能作为 MATLAB 临时空间上界。新 launcher 每个作业前要求任务盘和 `/tmp` 各至少 10 GiB 可用、预计增量内存不超过当时 available RAM 的 70%、无其他 MATLAB 进程，并保存 shell 预检、stdout/stderr、真实退出码、进程墙钟及输出 SHA。当前可用盘远大于 10 GiB；预计最大 BIE 增量内存约 1.460 GB，低于当前 70% available 门，但每次仍须重测。launcher 串行、`-singleCompThread`，每作业冻结上限 90–350 s（均严格低于 900 s），累计已用进程秒数加下一作业上限不得超过 1800 s；超时或任一门失败即停、不自动重试。历史双频十层 BIE 估计约 `859 s`（含每层 40 s 启动），只是资源规划，`-singleCompThread` 的实际耗时需首频实测后重新估计。新目录的输出/临时空间持续监测，若不足 10 GiB 或实际耗时推算超预算即停止，保留失败日志。

**本次需要用户明确批准的三项操作：**（1）在 `/home/misc/cwx/uwa2d_s5_p22_dualfreq_20261029_01` 创建全新独立目录；（2）只向该目录同步清单所列 13 个文件并逐项核 SHA；（3）按上述串行门运行最多 1800 s 的远端 MATLAB 双频物理资格先导。批准后先做 R2023b/6 kHz 兼容性 smoke，再 4062.5 Hz 同源 PE/FullFK、五层 BIE；合格后且 4031.25 Hz BIE 尚未运行时冻结复预测和 SHA，再运行半步点。若用户只批准部分事项，只执行相应已批准且不依赖其他批准的操作；绝不把“远端已恢复”视为物理计算批准。当前阶段 1 在**独立任务快照方案**下为只读资源 `PASS`，阶段 3 为 `AWAIT_EXPLICIT_APPROVAL`，阶段 4–5 仍 `NOT_RUN`。原 S4/S5/S6 状态不变。
