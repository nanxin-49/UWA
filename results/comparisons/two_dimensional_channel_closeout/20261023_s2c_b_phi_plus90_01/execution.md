# PHI_PLUS90 执行回执（2026-10-23）

## 2026-10-25 图像导出资格补记（当前）

在本地原 R2025b MCP 会话从保存分析 MAT 重绘 `figures_final_complete_v2/`，仅固定 M05c 右上纵轴刻度；没有新增声学评价、传播或求解。[严格审计 v2](audit/final_independent_audit_phi_v2.json)的五张 PNG 独立重绘逐像素一致，数值/来源/公式/支持门也通过。混合平台同网格 BIE 求解差仍未直接定界，状态 `S2C_B_PHI_PLUS90_DATA_COMPLETE_REFERENCE_LIMITED`，共享扣记维持 `3401.4241040/7506.601 s`；下节是图像修复前快照。

## 2026-10-25 保存密度跨平台边界复核（当前）

在远端新建唯一 `boundary_operator_20261025T065734Z` 只读作业，重新核 `main96/self64` 原 SHA，启动前 `/proc/meminfo MemAvailable=27106107392 B`，预计增量 268435456 B，通过 30% 保留门；共享预算即使按单作业 900 s 最坏情形和 S3A 3000 s 也通过。R2023b 对两份 Windows 保存密度独立评价各 408 个离节点边界点，时间 `60.737829/39.485799 s`，未调用求解器的线性方程求解或写新密度。残差向量与 Windows 原件最大绝对差 `1.04e-16/8.09e-17`。[数值回执](audit/remote_saved_density_boundary_20261025.json)及[前后资源/退出码日志](remote_jobs/boundary_operator_20261025T065734Z/)。进程从 `06:57:34` 至 `06:59:29 UTC`、exit 0，按整段 115 s 保守扣记；最新共享累计 `3401.4241040/7506.601 s`、余额 `4105.1768960 s`，留 S3A 3000 s 后余 `1105.1768960 s`。对应 [FIG 对象复核](audit/visual_semantics_phi_v1.json) 全部数值对象一致，M05c 严格像素仍未过。下节是本次额外复核前快照，保留原始记账过程。

## 2026-10-25 SSH 远端续算与本地汇合（当前）

- 现有 `matlab-ubuntu` 别名在受限文件沙箱中不可解析；使用已授权 SSH 访问后，确认远端 `misc`、真实路径 `/home/misc/cwx`、独立且原先不存在的工作区 `/home/misc/cwx/UWA_phi_plus90`、R2023b `23.2.0.2365128`、12 逻辑 CPU、约 26.2 GiB `MemAvailable`、142 GiB 可用磁盘。未接触其它项目或进程，未安装软件。MATLAB 默认仍为本地 MCP；本次是明确授权的 SSH 远端作业。
- 用本地工作树构建 50 文件选择性快照，含冻结源/配置/417 坐标、S2b 正式参考、S2c-A 协议、K_HIGH 最终审计摘要、现有 PHI 两层检查点、全部实际依赖源码与规范。远端逐文件长度/SHA 50/50 一致；`main96/self64` MAT 原哈希分别 `265baf00…b1478a5`、`4b2efc30…b9413481`，JSON 原哈希分别 `1ffc2de2…b8c9d7`、`39e41088…aff08b`。无 `cash/`、无 GitHub 克隆、无历史结果覆盖。[完整快照清单](remote_recovery/snapshot_manifest.json)。远端只对 BIE 作业包装的内存来源文字改为 Linux `/proc/meminfo MemAvailable`；生产 BIE 求解器 SHA `76f1a9d5…349fa9` 和冻结物理配置未变。
- 跨版本预检加载 984 项 Gaussian 原离散源、417 点/两份密度、解析 C2 海面和法向，验证 Bessel/Hankel 单值控制。Linux 对两份 Windows 密度重新评价 417 点场，相对差均为 `0`，用时 `0.282620/0.192037 s`。首次验证包装脚本在评价 main96 后因 R2023b 空结构体赋值差异失败；只修改专用验证包装，随后通过，无 BIE 重解。按两次验证进程墙钟 `12.011566+10.6699636=22.6815296 s` 保守扣记，包含失败和真实评价；非声学版本探针另有 MATLAB 启动，不计入 PHI 声学核。[兼容性回执](remote_recovery_compatibility.json)。同源码/同密度评价一致无法单独界定跨平台边界**求解**差。
- 每次新物理作业先核旧 SHA、前序回执、唯一作业目录、当前 `/proc/meminfo MemAvailable` 和 30% 余量；`space14/window60/window70/pe_phi/fullfk_phi` 即时可用物理字节依次 `28161662976/28178563072/28196216832/28180152320/28198072320`，均远高于本案最大 1,418,723,328 B 估计增量所需门。独立 SSH MATLAB R2023b 进程有 PID、资源前后记录、日志、真实退出码和 900 s timeout，五项均 exit 0；无并行池或重复作业。[逐项日志](remote_jobs/)与[回传清单](remote_recovery/remote_return_manifest.json)。
- 真正远端续算顺序与声学实耗：`space14` L60/14ppw/qa96 `96.959982 s`，`window60` L60/10ppw/qa96 `50.767938 s`，`window70` L70/10ppw/qa96 `67.880350 s`，PE-mode0 主/半步 417 点及冻结 13 探针倍网格 `2.903151 s`，完整 `FullFK_complete_analytic_flat_exterior_GL8_v1` 的 GL8 400/800 和 L60/L90 `2.352785 s`。所有新作业的 MAT/JSON/资源日志回传，69/69 文件 SHA 与长度相同；四份旧 MAT/JSON 与本地相同，65 新文件增量写入，没有同路径异 SHA 覆盖。Linux 原址仍可追溯，本地本 run 为正式汇合位置。
- 五层参考审计得到整线/九点经验相对变化 `2.5896501e-7/1.5783384e-7`，PE 平面同源控制 `4.2238524e-13`，FullFK GL8 400/800 `6.3234807e-14`、L60/L90 `7.5256492e-16`，组件恒等式最大差 `1.96e-17`。本地 R2025b 仅进行保存数据分析、M05a–e 绘图及独立审计，没有再做声学求解。严格逐像素复绘门 M05c 未过，且混合平台求解差未独立定界；结果状态 `S2C_B_PHI_PLUS90_DATA_COMPLETE_REFERENCE_LIMITED`，详见[审计](independent_audit.md)。
- 七项 PHI 声学实耗 `578.0285102 s`，加兼容过程保守 `22.6815296 s`，累计共享 `3286.4241040/7506.601 s`，余额 `4220.1768960 s`；保留 S3A 3000 s 后余 `1220.1768960 s`。旧恢复门中的 `3042.8783684 s` 是远端续算前快照，不追溯改写。

## 2026-10-25 内存释放后恢复与再次资源停止（当前状态）

用户确认内存已释放，随后实测 Windows 可用物理内存 3,163,262,976 B，高于五层中最大 `space14` 所需的 2,026,747,612 B。当前 MATLAB MCP 会话为 R2025b Update 6、项目根目录；首次入口调用因会话路径未初始化而在求解前停止，`jobs/` 仍为空。运行 `setup_vertical_project` 后确认入口可见；没有从 shell 启动 MATLAB。冻结 PHI 配置、源、411 线/九点 SHA 均与协议相符，BIE 求解器 SHA `76f1a9d5…349fa9` 与 K_HIGH 五份正式回执一致；K_HIGH 七份作业回执均 `COMPLETE`。

按序执行 PHI 新边界解 `main96`（3840 节点、qa96、417 点，实耗 203.6039647 s）和 `self64`（同节点、qa64、417 点，实耗 153.5603395 s）。各作业前即时可用物理 RAM 分别为 2,792,235,008 B、2,648,829,952 B，均通过完整剩余层级 30% 余量门；两份 MAT 含完整边界密度、原生接收复场、源和回执，并有 SHA 锚定。两层线性残差分别约 `3.34e−15/3.30e−15`，离网格边界残差 `3.81e−9/2.04e−9`。从保存 MAT 独立核对同节点/法向与 417 同坐标，自求积变化在整线与中心 11 点分别为 `3.1105e−8/3.1153e−8` 的相对范数。见[两层部分审计](audit/P2_recovery_partial_two_level_audit_20261025.json)，**它不构成五层经验分辨率或正式参考资格**。

第三层 `space14` 作业前的真实可用物理 RAM 降至 1,863,331,840 B，低于 2,026,747,612 B 门槛，故未调用 MATLAB 求解器。查看现有 MATLAB 进程占用后仅作一次资源复核，1,879,687,168 B，仍失败；未忙等、重复启动、清理用户进程或降低节点规模。`space14/window60/window70`、PE 与 FullFK 均未运行；M05b–e 与 PHI 中心误差判断仍缺完整正式场。详见[恢复资源门与作业 SHA](P2_recovery_gate_20261025.json)。

本次两项声学实耗合计 **357.1643042 s**；共享核心累计 **3042.8783684/7506.601 s**，余额 **4463.7226316 s**。保持 S3A 3000 s 后，PHI 仍可用 **1463.7226316 s**；按 K_HIGH 余三层 1.2 倍加其余/失败各 150 s 的保守剩余计划 **1252.49592096 s**，成本门仍通过。当前阻塞仅是即时物理 RAM。以下 2026-10-24 与 2026-10-23 记录保留为各自当时状态，不改写。

资源阻塞期间补齐了保存数据专用的 M05b–e 绘图入口 `scripts/reporting/render_two_dimensional_closeout_s2c_b_phi.m` 与最终 CSV/重绘审计入口 `scripts/validation/support/audit_two_dimensional_closeout_s2c_b_phi_final.m`。绘图代码会从完整分析 MAT 生成可编辑 FIG/PNG，并从全部 417 真点保存同深度 ±x 配对 CSV；审计代码将复核五支持、幅相平方、五层数值分类、配对变化、图像像素和 FIG 结构。两文件 MATLAB Code Analyzer 零消息，**因 PHI 三方法正式场未齐而未执行**，未生成任何 M05b–e 或最终 PASS 回执。

随后审计代码再核对冻结源/配置/求解器 SHA、五层参考审计状态及 PE/FullFK 控制，并将 M05a 从已保存的 `geometry.mat` **独立重绘**，再与正式图逐像素比较；M05b–e 则从完整分析 MAT 独立重绘。该修正仍只属于未执行的验收入口，历史 M05a 数据及图像未被重算或覆盖。

## 2026-10-24 预算授权与 RAM 门补记（覆盖下方旧成本阻塞状态）

用户明确将共享 7200 s 上限增加至少 306.601 s，并要求保留 S3A 3000 s 规划储备。按最小增量记账，共享上限为 7506.601 s；前期扣记 2685.7140642 s，PHI 可用子额度 1820.8869358 s。原保守完整计划 1820.88667776 s 因而通过成本门，修订原件见 [P2_budget_amendment_7506_601.json](P2_budget_amendment_7506_601.json)。原 `cost_gate.json` 是授权前快照，不改写。

接着两次读取 Windows 可用**物理**内存，分别为 1,641,025,536 B 和 1,131,163,648 B；冻结最大 BIE 层级预计增加 1,418,723,328 B，按保留 30% 的规则至少需要 2,026,747,612 B。两次均未通过 RAM 门，最新测量见 [P2_ram_gate_20261024.json](P2_ram_gate_20261024.json)。用户随后要求暂时保持资源阻塞。未启动 BIE、PE 或 FullFK，新 PHI 声学耗时仍 0 s；共享核心扣记仍 2685.7140642 s，S3A 3000 s 储备仍完整。状态维持 `S2C_B_PHI_PLUS90_RESOURCE_REFERENCE_BLOCKED`，此轮直接原因已由成本改为 RAM。

已经增设只对 PHI 适用的预算门、五层 BIE、PE、完整 FullFK 双项积分入口；五个 MATLAB 文件的 Code Analyzer 均为零消息。入口尚未执行，不能据此宣称声场、参考资格或 M05b–e 已完成。下一次仅在用户告知已释放内存后重测物理 RAM；通过后仍逐作业重新过成本、RAM 和数值门。

后续只读准备补充了[五层 BIE 保存场参考审计入口](../../../../scripts/validation/support/audit_two_dimensional_closeout_s2c_b_phi_bie_reference.m)（尚未执行，Code Analyzer 零消息），以及[对称性输入审计](P1_symmetry_input_audit.md)：417/417 接收点有同深度 ±x 配对，五层 BIE 节点与 PHI 高度/坡度满足浮点精度内的偶/奇对称；PE 原 984 点 FFT 网格有一个不配对的 Nyquist 端点。此检查没有新增声场，不改变 RAM 阻塞状态。

保存场分析入口 [`analyze_two_dimensional_closeout_s2c_b_phi.m`](../../../../scripts/comparisons/analyze_two_dimensional_closeout_s2c_b_phi.m) 已将 BASE_M 与 PHI 的同点原生复场、平面控制、五支持、幅相平方分解、BIE 物理响应和模型残差分开组织，并检查预算修订 SHA。PHI BIE 参考审计入口另从保存节点核对五层高度和法向。两个入口 Code Analyzer 均为零消息，但缺新 PHI 场而未执行，表格、M05b–e 和独立结果审计仍不存在。

1. 在干净工作区创建唯一新 run；核对当前 S2b、S2c-A、K_HIGH 最终哈希清单，分别为 247/247、203/203、132/132。已有 MATLAB MCP 会话 R2025b Update 6，工作目录项目根。
2. 冻结 PHI 输入 SHA、原源和 417 点接收位置；从 S2c-A 真实解析几何生成 M05a。初版图例把辅助竖线和中心标记列入图例，在 `figures_v2/` 修正；最终又关闭下划线 TeX 解释，正式图在 `figures_final/`。旧图保留为过程证据。没有生成 PHI 声学场。
3. 按 K_HIGH 七作业回执重估五层 BIE，加入 20% 保守增量及两项 150 s 预留，完整计划 1820.8866778 s 大于保持 S3A 储备后的 1514.2859358 s。P2 成本门失败，未进行即时 RAM 许可检查，也未启动 BIE 代表、PE 或 FullFK。新声学耗时 0 s；共享累计与余额不变。
4. K_HIGH M04 小包复制首版清单生成时，PowerShell 在枚举输出目录时同时碰到正在写入的清单自身，产生一次文件占用报错；已核对最终清单恰含八个载荷，各与原件 SHA 相同，未覆盖来源。保存复场十字五点首次手动复核选错两个深度点，按 P12/P21/P22/P23/P32 的真坐标修正并通过；旧正式 CSV 未改。
5. 从 S2c-A/K_HIGH 保存 CSV 复核三档 A 与 K_HIGH 的空间支持和幅相平方项；核对 S3A 三个已有 PM Fourier 输入的系数 SHA、同父谱和旧 MAT 接收坐标，不执行 PM 求解。M05a 从保存 MAT 重新绘制后像素一致，FIG 重开三轴，静态代码检查零消息，`git diff --check` 通过。产物索引见 [M05图册](mechanism_atlas.md)、[独立审计](independent_audit.md)、[成本JSON](cost_gate.json)及[最终审计JSON](final_audit.json)。
6. 下一回合只读补查旧 PM pilot 的 `bie_reference_metrics.csv`：只有 seed260001 两个旧 L60 层级，seed260002/003 没有对应 BIE 回执。按 K_HIGH 五层实耗直接乘三个 seed 已需 3802.217 s，加本轮 20% 保守增量仅 BIE 就约 4562.660 s，略超共享剩余；这是未经 PM 节点/版本资格确认的规划情景，不作为 S3A 正式估价或执行许可。详见 [S3A只读清单](S3A_readonly_feasibility.md)。

阶段分类为 `S2C_B_PHI_PLUS90_RESOURCE_REFERENCE_BLOCKED`，不等于中心过零机制已验证。未改变生产求解器、公共 API、旧正式数据、3000 s S3A 规划储备；未提交或推送。
