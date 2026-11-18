# R2：固定 4 kHz、10 cm 条带的有限跨 Seed 复核与手动执行方案

2026-11-17。**本轮仅准备执行包，不连接 SSH、不启动远端 MATLAB、不运行正式声学求解。** 用户随后明确允许两个 Seed 同时计算；方案因此支持 **两个不同 Seed 各一个单线程进程**，无 parpool、无线程池、无额外 Seed。实际计算由用户手动启动。R1 原始结果、生产代码和已批准数值门均不改。

## 问题与复核单位

以 [R1 最终复核](r1_bie_v2_final_review_2026-11-16.md) 为锚点：seed260001 的 PS/BIE 复场相对 L2 为 0.15923498070661，完整 FullFK/BIE 为 0.000747208248357；相位项分别贡献约 99.7691% 和 99.9437% 的平方复残差。目标是复核新增 seed260002/260003 上的量级、可分辨性、排序和空间分布。单个 Seed 是一个复核单位，411 点是该海面的空间采样，不能算成 411 次独立重复。三 Seed 结果只作描述性对照，不作总体置信区间、显著性检验或海况适用域宣称，也不新设工程合格误差。

固定顺序为每个 Seed 内 models → C0 → C1 → C2 → C3 → C4 → final；两个独立 Seed 可并行，同步不改变公式或采样。固定环境、源码和输入 SHA 限制计算批次差异；两个运行的实际耗时、RSS、可用内存分别记录，不能将并行墙钟直接当作算法加速比。

## 冻结配置与输入来源

- 理想二维无限长刚性障板均匀条带，总宽 0.10 m，半宽 0.05 m，v0=0.001 m/s，源中心 x=0、z=100 m；4000 Hz，c=1500 m/s，rho=1025 kg/m³。
- 原生 `exp(-i omega t)`、absolute physical peak Pa、phase contract 2.0.0。海面 helper 的值是物理深度 z_s；不改 PE/Green/Fourier 符号。
- 160 m 周期、12 非零模、Kmax=0.471238898038469 rad/m，原 Sk；C2 平台 |x|≤42 m、|x|≥50 m 归零。无再抽样、去均值、逐 Seed 归一化或带宽变化。
- 原 411 个真实接收点：x=(-205:205)×0.1953125 m，z=3 m。冻结时从 R1 已保存 C4 点表复制坐标十进制 token，并核对全线；直接逐坐标求值，不插值。
- 260002/260003 系数逐字节复制自 `results/validation/spectral_fast_kirchhoff_pm_pilot/seed_<seed>_coefficients.csv`，并与其明确父来源 `results/validation/pe_bellhop_pm_amplitude_ensemble/seed_<seed>_coefficients.csv` 做 SHA 对照。二者是同 PM 谱的独立 Gaussian 余弦/正弦系数样本，**不是固定模态振幅的纯相位旋转样本**。Sk 和 k 必须与 R1 完全一致，实际几何统计允许随 Seed 自然变化。
- 每包只有本 Seed 系数、R1 系数锚点、坐标、JSON 配置、政策和代码，不含 R1 或另一 Seed 的 BIE 密度、common.mat 或 Gaussian 场。源、源梯度、边界 rhs 与 C0–C4 密度均按本 Seed 重新建立。

冻结目录为 `results/comparisons/piston_2d_pm/r2_preparation_20261117_03/seed260002/` 和 `seed260003/`。每包 `manifest.json` 绑定所有输入/源码，`provenance/sources.json` 记录来源和 SHA；`policy.json` 从 R1 回收包逐字节复制。`input_package.tar.sha256` 绑定交付归档。正式运行后才产生 jobs、会话、检查点和最终结果，准备阶段不伪造最终场或 PASS。

## 原算法与参考资格

PS 继续使用 R1 的精确平面反射背景与屏修正 `(1-exp(-2ik*z_s))*p_inc(x,0)`，通过原平方根 PE 角谱算子及物理载波传播到 z=3 m，含倏逝谱，直接 Fourier 求值到 411 点。PE 是传播算子，PS 是反射近似。`run_r2_piston_pm_job.m` 只适配冻结输入与 Seed，PS 局部函数保留 R1 原文。

完整 FullFK 继续使用 `closeout_s3a_pm_fullfk`，包含入射法向导数项和 Green 法向导数项、真实法向、弧长 Jacobian 及解析平面外部背景；不替换成 FastFK 或单项 Kirchhoff。源仍调用 R1 的 `baffled_strip_piston_field_vertical`，不改任何生产实现。

| 模型层 | 固定配置 | 原数值门 |
|---|---|---|
| PS 六层 | W/dx/dz：128/.05/.1875、256/.05/.1875、512/.05/.1875、1024/.05/.1875、1024/.025/.1875、1024/.025/.09375 | W512→1024 和 dx.05→.025 的全线/四段最坏复 L2 <1e-4；步长差 <1e-10 |
| FullFK 四层 | panels/extent：400/50、800/50、960/60、1440/90 | 三个配对全线/四段最坏复 L2 均 <1e-5；最终模型为 800/50 |
| BIE C0 | L60、ppw10、qa96、3200 节点 | 本 Seed 新求解与独立审计 |
| BIE C1 | L60、ppw12、qa96、3840 节点 | C0→C1 <1e-4 |
| BIE C2 | L60、ppw12、qa64、3840 节点 | C2→C1 <1e-4 |
| BIE C3 | L60、ppw14、qa96、4480 节点 | C1→C3 <1e-4；最终参考始终为 C3 |
| BIE C4 | L70、ppw12、qa96、4480 节点 | C1→C4 <1e-4；仅窗口检查，不替换 C3 |

BIE 固定 GL8 面板、halfplane h=-2 m、窗口平台比 .85，沿用验证过的 opt-in 静态缓存。缓存只缓存本次调用的求积/插值模板及源面几何，不复用跨 Seed 散射解。单调性沿用 `r1-bie-reference-v2`：保存第一/第二空间差、比值、增幅、非单调标记及 legacy v1 状态，**单调性只作诊断**；四项敏感性原门不变。R1 原 v1 FAIL 留在原证据中。

每层必须保留全部平台中点：C0 340、C1/C2 408、C3/C4 476，另有 19 个 GL112 独立边界点。原门：线性残差 <1e-10，完整中点和独立边界残差 <1e-4，密度接收重建/rhs/MAT-CSV/分量闭合 <1e-11，几何 <1e-12，独立源/梯度 <1e-7。正式入口复用 `audit_r1_bie_manual_stage` 与 `audit_piston_pm_r1`，不减点、不只看保存 PASS 或退出码。

最终另独立重算四敏感性、C3 参考选择、mask、误差及不确定性，核 PS 直接 Green、FullFK 双项、源/梯度与 rhs、BIE 密度接收重建、MAT/CSV；图件保留 PNG/FIG，F02 八条和 F04 四条曲线逐条核对并重导 PNG 比像素。R2 去掉仅针对不同历史声源的 Gaussian 配对图，不将它扩展成新研究。

## 指标、分辨率与三 Seed 汇总

对每个 Seed 和模型 m∈{PS,FullFK}，以该 Seed C3 反射场 b 为分母：

`E_m=||u_m-b||2/||b||2`；`A_m=||||u_m|-|b|||2/||b||2`；`scale=||b||2/sqrt(411)`。

逐点输出复残差实/虚部、绝对值 Pa、幅度差 Pa、包裹相位差 rad，以及 `abs(residual)/scale`。固定 `abs(b)>=0.01*scale` 相位显示掩码，报告有效点数；不按 Seed 调阈值，不作逐点局部相对误差主排名。零掩码点的相位 RMS/最大值写 NaN，不能伪报零。

稳定平方分解为 `|u-b|²=(|u|-|b|)²+4|u||b|sin²(Δphi/2)`。保存幅度/相位能量、相位项占总残差能量比例、按掩码内 BIE 功率归一化的相位 RMS、最大相位和幅度差及位置、最大复残差位置、最大 41 点的能量占比。平方分解使用所有点，显示相位和加权 RMS 才使用掩码。BIE 粗糙场相对同源平面场变化单列，不与模型误差混合。

`bunc=|C0-C1|+|C1-C3|+|C2-C1|+|C1-C4|`，`punc` 和 `func` 分别为原 PS 三配对和 FullFK 三配对逐点绝对差之和。原分辨条件为：

`E_m > 5*||bunc+model_unc||2/||b||2`。

四个 BIE 配对均须先通过原门；随后无论模型误差是否被分辨，都如实报告，不用“很小”替代分辨判断。经验不确定性不是严格误差界或统计置信区间。模型误差/排序改变属于有效研究结果，不触发调参或重选种子。

跨 Seed 的 `summarize_r2_pm_fields({R1_work,R2_260002_work,R2_260003_work},new_output_dir)` 仅读保存数据，要求三个 Seed 顺序固定、公共源/坐标/k/Sk/模型配方完全相同、各自综合/最终审计通过，并独立重核 CSV 与复场及指标。输出六行模型表、2466 行逐点表、24 行分段能量、24 行共享端点区间范数、三行排序/误差比、输入 SHA 回执和报告。先分别报告每 Seed，再描述三者范围、中位趋势及空间峰移动；不把接收线点混成总体样本，也不相干平均不同海面的复场。

范数区间保留 R1 `round(linspace(1,411,5))` 的共享端点，收敛门定义不改；能量区间使用互斥 x<−20、[−20,0)、[0,20)、x≥20 m，各点只计一次。每 Seed 的同一归一化尺度用于整条线，跨 Seed 比较采用各自 BIE 范数归一化。

## 成本与并行资源

R1 现存计时只作估计，不是 R2 实测或保证：C0=55.481928 s（缓存完整回放）、C1=80.841146 s、C2=65.133454 s、C3=108.641842 s、C4=106.646117 s，五次合计约 416.7445 s/Seed；原 models 30.715001 s。再加五层独立审计、综合/图件审计和 MATLAB 初始化，建议规划 **10–20 min/Seed 的观察窗口**，并行总墙钟待实测，不承诺加速两倍。十次独立 BIE 新求解是全部 R2 上限，无额外配置/Seed/频率。摘要分析不含新声学求解。

按 R1 手动机制，不恢复已取消的 1100 s 累计门或成本预测硬门，也不新增自动超时重试。本方案没有自动墙钟强杀；观察窗口超出时人工查询实际阶段、PID/资源和日志，再决定是否停止，不能因启动空档重派作业。每会话保存完整 stdout/stderr、diary、`time -v`、1 s 采样的 MemAvailable/进程组 RSS、真实退出码、阶段前后资源和检查点。两个进程计时分别记账；RSS 进程组和 `time -v`单进程口径不混用。

每 Seed 继承原 **12 GiB 工作内存估计**，两者合计 24 GiB。以共享组 `/home/misc/cwx/r2_pm_20261117_03_control/` 的可用内存基准为准，基准只取当前/既有启动基准的最大值，不在恢复时降低。每次启动要求 `MemAvailable - .30*shared_baseline >= 24 GiB`，即首次大致需 ≥34.286 GiB 可用内存；第二次仍预留完整 24 GiB，是有意保守的准入。不能凭总装机容量或多核数宣称准入通过。12 GiB 是估计，不是物理内存上界；实际低于共享基准30%时，各进程守卫/watchdog停止自己的进程组，两个 Seed 都会观察同一资源门。

每 Seed 独占自己的 workflow lock 和固定 Seed 槽锁（仅 260002/260003 两个槽），不同 Seed 可并行；同 Seed 不可重复运行。准入锁只串行化共享 RAM 基准更新，不串行化计算。每个 MATLAB 使用 `-singleCompThread`，OMP/MKL/OpenBLAS=1，不启 parpool。保留原 R1 runtime 的 initial/starting/max 独立记录，另存 shared_runtime.json，避免篡改原 schema。共享组保留为同一批次证据，不自动删除或重置其 RAM 基准。

## 停止、恢复与证据边界

1. 输入/代码/政策/配置 SHA 或几何不符、发行版本不是 R2025b、MCP/远端依赖失败、RAM 准入不通过，均先停止。准备端 MCP 不可用不会触发 SSH 或本机 MATLAB fallback；用户未来手动启动时，**同一个单线程远端进程先做带七个声学哨兵的 dry-run、恢复状态机/政策测试和 Code Analyzer，再开始正式入口**。这些尚未实际执行的 MATLAB 验收不能在本轮记 PASS；若验收异常则正式求解零派发。
2. 每 Seed 的数值失败停止其后续阶段，不扩大窗、加网格、改门或重抽样。另一个已手动启动的独立 Seed 不因某 Seed 的数学失败自动重启/停止；其结果单独判定。共享资源门失败时两者各自停止。失败/停止都不删除产物。
3. 每次模型计算和每个 BIE 层调用前写 `solve_invocation.json`。一个阶段只允许一次新求解调用。已有 committed 且全部 SHA/前序签名匹配则跳过；已有 solved seal 则只审计；存在未封存调用、孤立 MAT/CSV、篡改或前序变化则拒绝自动重算。**Continue 不代表允许重新求解失败层。** 此时需先检查已有日志/真实进程/检查点，另行审查明确恢复范围；不删除 marker 绕过。
4. 签名包含本包 manifest、Seed 及前序提交 SHA；每次阶段审计在新 attempt 中进行，旧失败审计不覆盖。阶段提交绑定密度/场/配置/计时/调用 marker/solved/审计/门及重建产物。最终从已验收阶段复制到新 finalization attempt，只分析保存密度，不再求解；失败可手动重新 finalization，保留旧 attempt。
5. 模型门、四项 BIE 门或独立审计未完成时不能写参考合格；进程退出0也不能单独当 PASS。最终必须同时核 checkpoint COMPLETE、真实退出0、所有层提交、综合/最终/图件审计与 finalization seal；Fetch 只提供传输完整性，还需保存结果审查。收回后复核输入及输出 SHA、MAT/CSV、图件及分辨率；不从 R1 密度填补缺失。
6. 每个 Seed 最终状态均是 `R2_SEED_COMPLETE_WAITING_REVIEW`，不自动启动另一 Seed、汇总、R3、宽带、CIR/QPSK/BER 或 Bellhop。用户可一次手动启动两个包，程序没有后续派发器。

## 文件与接口

研究适配器：`scripts/comparisons/run_r2_piston_pm_job.m`；有界编排：`scripts/validation/run_r2_pm_seed.m`；验证/最终审计：`validate_r2_pm_preparation`、`support/audit_r2_pm_final`；保存图与汇总：`scripts/reporting/render_r2_pm_seed.m`、`summarize_r2_pm_fields.m`。打包/手动管理入口在 `scripts/validation/support/prepare_r2_pm_package.ps1` 和 `manage_r2_pm_seed.ps1`，Linux 启动/回收/状态与无网络 preflight 同目录。复用的 R1 数值核逐字节冻结并核 SHA。

公共 `vertical_channel_model`/传播器/通信接口、公式和默认行为无变更，无需通信回归或修改 `vertical_comm_guide.md`。脚本导航、任务表、当前状态与进度源需记录准备完成和未执行项；本地验收以 [本地验证报告](r2_pm_preparation_validation_2026-11-17.md) 为准，具体用户命令见 [PowerShell 操作说明](r2_pm_manual_commands_2026-11-17.md)。完成准备后停止。
