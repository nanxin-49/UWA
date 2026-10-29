# 脚本与验证运行索引

**2026-10-29 二维 S5/S6 最终轻量控制：**`validation/test_two_dimensional_s5_bie_receiver_coordinates.m` 在现有 MATLAB MCP R2025b 会话 5/5 通过；`validation/validate_two_dimensional_s5_nine_point_cir_controls(output_dir)` 以新目录保存九点解析 H/CIR、65/129 网格和端点检查（最终 `20261029_04`）；`validation/validate_two_dimensional_s6_analytic_dsp_controls(output_dir)` 保存无噪声与 Gray-QPSK AWGN 真实计数、PER/EVM、事前二项区间（最终 `20261029_04`）；`validation/audit_two_dimensional_s5_s6_analytic_controls(cir_dir,dsp_dir,audit_dir)` 独立重开 MAT/CSV（最终审计 `20261029_07`）。`validation/audit_two_dimensional_s5_cost_projections(single_dir,dual_dir,audit_dir)` 独立核旧 970 行单/双锚点规划表。详见[正式回执最新节](../reports/two_dimensional_s5_s6_preflight_2026-10-29.md)。这些是解析/软件和预算算术资格，不是物理 S5 CIR 或正式 S6 BER。

**2026-10-29 二维 S5/S6 轻量资格续行：**九点解析 CIR 控制已修正 P11–P33 的行优先点号顺序并显式区分 65 含端点频率、129 嵌套加密频率与 64 个 IFFT bin；BIE 接收单测覆盖原九点坐标、逐点/批量、旧标量深度及非法输入。当前 MATLAB MCP 无法附着共享会话，SSH 别名无法解析，本轮 Code Analyzer、单测及解析 CIR/QPSK 均**未运行**；不得把改动或旧 Analyzer 结果写成 PASS。见[本轮回执](../reports/two_dimensional_s5_s6_preflight_2026-10-29.md)与[检查点](../results/validation/two_dimensional_s5_s6_preflight/20261029_01/checkpoint.json)。

**2026-10-27 二维 S5 前置控制：**`validation/validate_two_dimensional_s5_prerequisites(output_dir)` 在新结果目录验证 native↔DSP、解析双径及 FFT 往返，不调用声学求解器。阶段 A 物理源/参考/资源门阻塞，见[中文回执](../reports/two_dimensional_s5_s6_2026-10-27.md)。

**二维 S5 原九点解析 CIR 代数控制（新增、待运行）：**`validation/validate_two_dimensional_s5_nine_point_cir_controls(output_dir)` 对归一化已知点路径与复系数双径，按全局 `tau_ref` 覆盖原九点，保存 H 分量、CIR、MAT/CSV/JSON 回执并检查复误差与到达采样容差；不调用 PE/BIE/FullFK，也不代表冻结 Gaussian 源的物理平面结果。入口尚未由 MATLAB MCP 运行，不能记作 PASS。

**2026-10-27 二维冻结源平面续行（先前步骤）：**`validation/support/closeout2d_source_at_frequency` 按固定 Fourier 源重新计算 `k_z(f)`；`validation/validate_two_dimensional_s5_source_flat(frozen_source_path,output_dir,available_physical_bytes)` 在新目录仅检查 4/6/8 kHz 的 P22 平面 PE 直达/反射与解析镜像。旧 S3A PM FullFK helper 的 4 kHz 行为不变，新源可显式给频率；该步骤当时尚无粗糙面 BIE 资格，[续行数据](../results/validation/two_dimensional_s5a_source_flat/20261027_s5a_flat_01/)。当前状态见下一条。

**2026-10-27 二维 6 kHz 弱正弦 P22 试点：**`validation/validate_two_dimensional_s5a_weak_p22` 与 `validate_two_dimensional_s5a_weak_pe_sensitivity` 在已有 MATLAB MCP 会话计算同源 PE/FullFK 及半步长/同谱双网格；`validation/run_two_dimensional_s5a_weak_bie_remote` 由验证专用 `support/launch_two_dimensional_s5a_weak_bie.sh` 在获批独立 Ubuntu 目录串行做 smoke、配对模型与五层 BIE，并留真实退出码/日志；`validation/audit_two_dimensional_s5a_weak_bie` 重开回传 MAT/CSV 做 1e−4/五倍门。结果限单频 P22，[中文报告](../reports/two_dimensional_s5_s6_2026-10-27.md)与[原数据包](../results/validation/two_dimensional_s5a_weak_bie/20261027_01/)，不视为 S5/CIR 或 S6 通过。

**2026-10-27 二维 8 kHz 弱正弦 P22 试点：**`validation/run_two_dimensional_s5a_weak_bie_8k_remote` 与 `validation/support/launch_two_dimensional_s5a_weak_bie_8k.sh` 在获批任务子目录执行 R2023b smoke、配对 PE/FullFK 和五层 BIE；`validation/audit_two_dimensional_s5a_weak_bie_8k` 独立重开 MAT 并计算敏感性/五倍分辨率门。[原始日志、87 项输出 SHA、审计和检查点](../results/validation/two_dimensional_s5a_weak_bie_8k/20261027_01/)齐全。A 仅新增一个 8 kHz P22 参考锚点，非宽带/CIR；图入口 `reporting/render_two_dimensional_s5a_weak_bie_8k` 已写但本轮 MATLAB MCP 调用被自动审批 404 阻断，没有 `.fig/PNG`。双锚点成本规划 CSV 同样缺 MATLAB MCP 独立重审，见报告中的限定。

**二维 S5 全带首批预算仍待批准：**[预算提案 v2](../results/validation/two_dimensional_s5a_fullband_cost/20261027_02/budget_scope_proposal_v2.json)请求固定平面解析/PE/FullFK 对照及原弱正弦九点、129 频点 BIE/PE/FullFK，累计上限 48 h；逐层 RAM 保留、`<900 s` 作业、零并行与阶段边界停止均为强制条件。PM 的 S4 参考仅覆盖固定 4 kHz，不构成宽带资格，因此本批明确留 PM 缺口；CIR/S6 也不由该预算自动授权。仅规划提案，未获批未运行；范围和证据边界见[中文报告](../reports/two_dimensional_s5_s6_2026-10-27.md)。

全带复用九点的 BIE helper 已接受配对接收深度向量；`validation/test_two_dimensional_s5_bie_receiver_coordinates.m` 覆盖新旧输入契约。MATLAB Code Analyzer 已检查该单测文件（0 项），helper 本身未在这次调用中直接受检；单测数值执行被自动审批服务 404 阻断。首个获批频点要核对九点复场并测量评价耗时后再按预算继续，见[报告](../reports/two_dimensional_s5_s6_2026-10-27.md)。

九点分属三个深度，PE每频按三层输出整条横向场再抽点；全九点及P22初算统一使用预先冻结的全局 `tau_ref=min_r[hypot(x_rx,z_tx-z_rx)/c]` 转 DSP，并保留 `hypot(x_rx,z_tx∓z_rx)/c` 得到的逐点理论到达差，禁止逐点或逐路径移峰。预算提案 v2 已冻结此相位约定，尚待全频数值验证。

**2026-10-27 二维 129 点 A 源/直达检查：**`validation/validate_two_dimensional_s5a_source_band(frozen_source_path,output_dir,available_physical_bytes)` 在原计划加密频点固定 984 项 Tx Fourier 场，仅计算 PE 与同源角谱直达，检查 65 点嵌套、传播谱及 `1e−10` 源/直达门；[MAT/CSV/审计/检查点](../results/validation/two_dimensional_s5a_source_band/20261027_01/)已保存。没有反射场、CIR 或通信资格。

**2026-10-27 三维固定dual 65点网格四频点验收：**`validation/validate_deepsea3d_65grid_four_frequency.m`和`run_deepsea3d_65grid_four_frequency.m`仅在已完成独立目录执行4750/7250/4125/6562.5Hz的原FE/FullFK门；`validation/support/freeze_deepsea3d_65grid_four.py`事前锁五锚点B预测，`audit_deepsea3d_65grid_four_saved.py`只读重算保存复数/门/成本/SHA，`reporting/render_deepsea3d_65grid_four_v3.m`从审计CSV绘制最终两张FIG/PNG，`validation/audit_deepsea3d_65grid_four_figures.m`重开FIG核坐标。四点条件资格和冻结`e_R/e_Δ/a_D`门各4/4通过，总保守进程账341/720s。见[正式报告](../reports/deepsea3d_65grid_four_holdout_2026-10-27.md)与[结果索引](../results/validation/deepsea3d_65grid_four_20261027/README.md)。目录及入口不授权重复运行、扩展频点或CIR/BER。

**2026-10-26 三维固定dual四频点独立盲测：**`validation/validate_deepsea3d_four_holdout_frequency.m`复用旧FE/FullFK资格门、仅改变4500/5500/6500/7500Hz频率，`validation/run_deepsea3d_four_holdout_frequency.m`按单点成本上限执行；`validation/support/freeze_deepsea3d_four_holdout.py`在新FE前锁定A/B预测，`package_deepsea3d_four_holdout.py`打包19项SHA输入，`audit_deepsea3d_four_holdout_saved.py`独立复算回传MAT/CSV与图。五个新远端MATLAB进程总计141/400s，四点各`CONDITIONAL_G1_PASS`，八个反射预测过原`1e−3`筛查。[正式报告](../reports/deepsea3d_four_frequency_holdout_2026-10-26.md)与[结果索引](../results/validation/deepsea3d_four_holdout_20261026/README.md)。仅GO后续网格设计；不自动重跑/扩频/生成CIR。生产API不变。

**2026-10-26 三维 6 kHz Q1 BIE 限额检查点：**`validation/run_deepsea_bie_q1_checkpoint.m`仅在本 Goal 冻结的平面 Gaussian 制造问题下执行连续节点 Q1 密度及局部 Duffy/Gauss 求积；`n16/q8`、`n24/q8`、`n24/q12` 已各执行一次，独立边界/离面场门失败，路线 `STOP`。只读独立审计在[结果包](../results/validation/deepsea_bie_q1_20261026/README.md)；[技术决策](../reports/deepsea3d_bie_q1_decision_2026-10-26.md)。此入口属于历史限定验证，不授权重跑、扩网格、真实 PM 或生产 API 集成。

**2026-10-26 三维固定dual 5/7 kHz：**`validation/check_deepsea3d_ubuntu_compatibility.m`只做一次冻结6kHz同配置跨平台诊断；`validation/validate_deepsea3d_5k7k_frequency.m`在独立case目录沿用旧 FE/FullFK 阶段门，`validation/run_deepsea3d_5k7k_frequency.m`按频点断点续算。`reporting/render_deepsea3d_5k7k_frequency.m`从五个保存的原生复场与旧锁定预测生成CSV、FIG/PNG；`validation/audit_deepsea3d_5k7k_figure.m`重开FIG；`validation/support/audit_deepsea3d_5k7k_saved.py`独立复算保存MAT/CSV/SHA/符号/预算，`finalize_deepsea3d_5k7k_receipts.py`生成检查点。先运行 `setup_vertical_project`，入口不依赖当前目录；远端求解只限本轮已授权的独立目录和预算，不应重用已完成目录。见[三维报告](../reports/deepsea_dual_5k7k_independent_validation_2026-10-26.md)与[结果索引](../results/validation/deepsea3d_5k7k_20261026/README.md)。生产求解器/公共API、二维S4不变。

**2026-10-26 S3A第三seed与三实现收尾（当前）：** `comparisons/run_two_dimensional_closeout_s3a_seed_job.m` 和 `run_two_dimensional_closeout_s3a_pe_grid_job.m` 仅增加对预算中的260003识别；第三seed用独立远端结果包完成五层BIE、PE、完整FullFK及双网格。保存场分析与F01–F03沿用现有入口。新增 `reporting/render_two_dimensional_closeout_s3a_three_seed(root,dest)` 从三份已审计正式CSV生成六行表、几何表和F08 FIG/PNG；`validation/support/audit_two_dimensional_closeout_s3a_three_seed(root,dest)` 独立交叉核CSV、幅相代数、五倍分辨率和F08逐像素。入口可在 `setup_vertical_project` 后用显式绝对root与dest从任意MATLAB目录运行；目标目录须不存在，防止覆盖。[第三seed检查点及三seed产物](../results/comparisons/two_dimensional_channel_closeout/20261026_s3a_seed260003_01/checkpoint.json)和[唯一主报告首节](../reports/two_dimensional_channel_closeout_report.md)为结果入口。S3A已3/3参考限定合格，S3B/S4不自动开始。

**2026-10-25 S3A第二seed续算（历史快照）：** 原有 comparisons/run_two_dimensional_closeout_s3a_seed_job.m、run_two_dimensional_closeout_s3a_pe_grid_job.m、analyze_two_dimensional_closeout_s3a_seed.m、reporting/render_two_dimensional_closeout_s3a_seed.m 与 validation/support/audit_two_dimensional_closeout_s3a_seed.m 仅按冻结 protocol/budget.json 参数化seed、预算和输入目录，不改生产求解器/API。新增 validation/support/audit_two_dimensional_closeout_s3a_figure_content.m(root,run_dir) 从保存MAT/系数核FIG中F01三条线、F02四面板各九个点号与真值、F03七条线；原审计另核MAT/CSV与FIG复导PNG像素。使用 setup_vertical_project 及显式root/run_dir。结果、SSH脚本、预算/输入SHA、全部退出码/日志/检查点和配对表在 results/comparisons/two_dimensional_channel_closeout/20261025_s3a_seed260002_01/；唯一主报告首节给数值与限制。260001保存结果不改，260003未运行；S3A为2/3。

**2026-10-25 S3A单seed PM试点：** `comparisons/run_two_dimensional_closeout_s3a_seed_job.m` 在用户已授权的独立Ubuntu工作包中逐层运行seed260001五层BIE、PE-mode0及完整FullFK；`comparisons/run_two_dimensional_closeout_s3a_pe_grid_job.m`补同一984模源谱的1968网格PE敏感性。`validation/support/closeout_s3a_pm_surface.m`定义保存系数的C2高度/坡度/曲率，`closeout_s3a_pm_fullfk.m`为完整双项/解析平面外场验证专用实现；`audit_two_dimensional_closeout_s3a_inputs.m`核三个旧输入，`comparisons/analyze_two_dimensional_closeout_s3a_seed.m`只读计算真实九点与411线指标，`reporting/render_two_dimensional_closeout_s3a_seed.m`只读生成FIG/PNG，`validation/support/audit_two_dimensional_closeout_s3a_seed.m`独立重算保存场/CSV并重开FIG、复导PNG。先运行根目录`setup_vertical_project`，各入口使用显式root/run_dir，不依赖当前目录；声学入口受冻结预算与物理RAM门约束，不得在无新阶段授权时重跑。结果在`results/comparisons/two_dimensional_channel_closeout/20261025_s3a_seed260001_01/`，结论见[唯一主报告最新首节](../reports/two_dimensional_channel_closeout_report.md)。单seed参考限定合格，S3A总体1/3，后续seed不自动运行。

**2026-10-25 M05 图像门修复（当前）：**`reporting/render_two_dimensional_closeout_s2c_b_phi.m` 只固定 M05c 右上角显示刻度；保存场和数值判据未变。`validation/support/audit_two_dimensional_closeout_s2c_b_phi_final.m(run_dir)` 原严格数值及逐像素逻辑未放宽，新审计输出 `audit/final_independent_audit_phi_v2.json` 通过五图重绘。图册改指 `figures_final_complete_v2/`；下条为修复前快照，混合平台 BIE 求解资格仍受限。

**2026-10-25 附加 PHI 审计：**`validation/support/audit_two_dimensional_closeout_s2c_b_phi_visual_semantics.m(run_dir)` 只读已保存的 M05a–e FIG/PNG 与独立重绘，逐对象比较轴、曲线、真实九点像素单元、标签和线型；输出 `audit/visual_semantics_phi_v1.json`，保留原严格 PNG 像素门的失败状态。远端任务专用 `remote_recovery/audit_phi_saved_density_boundary.m` 在 Linux 上重算 Windows 原密度的 408 个离节点边界残差而不重解密度；源码与资源包装保存在本案结果目录，数值回执见[审计](../results/comparisons/two_dimensional_channel_closeout/20261023_s2c_b_phi_plus90_01/independent_audit.md)。

**2026-10-25 PHI SSH 续算：**现有五个 PHI 声学入口在 `/home/misc/cwx/UWA_phi_plus90` 的 R2023b 上完成余三层 BIE、PE-mode0 与完整 FullFK；前两层 Windows R2025b 检查点不重算。本地 MCP R2025b 执行 `analyze_two_dimensional_closeout_s2c_b_phi`、`render_two_dimensional_closeout_s2c_b_phi` 和 `audit_two_dimensional_closeout_s2c_b_phi_final`。保存场分析及 M05a–e 图已生成；最终独立审计的数值段过关，M05c 逐像素门未过，混合平台边界求解差未定界，所以当前 `S2C_B_PHI_PLUS90_DATA_COMPLETE_REFERENCE_LIMITED`。远端选择性迁移、Linux 内存适配包装、作业资源/退出码与回传哈希见[本案执行记录](../results/comparisons/two_dimensional_channel_closeout/20261023_s2c_b_phi_plus90_01/execution.md)；下方早期“未运行”条目是当时快照。

脚本按职责分组，入口应能从任意 MATLAB 当前目录运行。推荐先执行根目录 `setup_vertical_project`；它只初始化路径，不改变当前工作目录。仍使用 `scripts/bootstrap_project.m` 的旧脚本会先委托给统一 setup，再按既有规则切换到 `results/<category>`。新入口应优先使用绝对项目根和明确输出目录，不依赖当前目录。

## 目录

**2026-10-25 PHI 续跑状态：**下条“新 PHI 场未计算”是 2026-10-24 快照。用户确认内存释放后，`run_two_dimensional_closeout_s2c_b_phi_bie_job` 已通过 MATLAB MCP 串行完成 `main96/self64` 两层并保存完整密度、417 点场与回执。`space14` 前真实物理 RAM 门及一次复核失败，剩余三层、PE、FullFK 和保存场分析/正式审计均未运行；不可由两层宣称参考合格。新 `reporting/render_two_dimensional_closeout_s2c_b_phi.m` 仅从验收后的分析 MAT 绘 M05b–e、复制既有 M05a，并保存 417 点对称性配对 CSV；`validation/support/audit_two_dimensional_closeout_s2c_b_phi_final.m` 从 CSV 独立复核五支持和数值分类、重绘比较五图、重开 FIG。两个入口仅静态检查，尚未执行或生成正式 M05b–e。最新成本、SHA 和继续条件见[恢复回执](../results/comparisons/two_dimensional_channel_closeout/20261023_s2c_b_phi_plus90_01/P2_recovery_gate_20261025.json)。

- **S2c-B PHI_PLUS90 资源受限几何收尾：**`reporting/render_two_dimensional_closeout_s2c_b_phi_geometry.m(root,new_figure_dir[,saved_geometry_mat])`只读 S2c-A 冻结解析数组（第三参数可从保存 MAT 重绘），生成 M05a PNG/可编辑 FIG/几何 CSV/MAT，拒绝覆盖目录；不调用声学核。正式图为 `20261023_s2c_b_phi_plus90_01/figures_final/`，早期 `figures/figures_v2/` 仅修图过程，`audit_redraw/` 与正式图像逐像素一致。用户已增加共享预算，五层保守计划通过成本门；当前真实物理 RAM 未过 30% 余量门，新 PHI 场未计算。`validation/support/closeout_s2c_b_phi_gate.m`只读取预算授权及配置 SHA；`comparisons/run_two_dimensional_closeout_s2c_b_phi_bie_job.m`、`run_two_dimensional_closeout_s2c_b_phi_pe_job.m`、`run_two_dimensional_closeout_s2c_b_phi_fullfk_job.m`及 FullFK 解析外场助手已静态检查，尚未运行。`validation/support/audit_two_dimensional_closeout_s2c_b_phi_bie_reference.m`将从五层保存场核对 BIE 经验分辨率与真实节点法向；`comparisons/analyze_two_dimensional_closeout_s2c_b_phi.m`将生成 BASE_M/PHI 同点五支持的原生复场、物理响应与模型残差表。二者已静态检查，缺真实 PHI 场而未运行。进入声学计算前须按[执行回执](../results/comparisons/two_dimensional_channel_closeout/20261023_s2c_b_phi_plus90_01/execution.md)重新核验 RAM；另见[输入对称性审计](../results/comparisons/two_dimensional_channel_closeout/20261023_s2c_b_phi_plus90_01/P1_symmetry_input_audit.md)、[M05图册](../results/comparisons/two_dimensional_channel_closeout/20261023_s2c_b_phi_plus90_01/mechanism_atlas.md)及[只读S3A清单](../results/comparisons/two_dimensional_channel_closeout/20261023_s2c_b_phi_plus90_01/S3A_readonly_feasibility.md)。

- **S2c-A只读诊断与协议：**`comparisons/run_two_dimensional_closeout_s2c_a.m(root,new_run_id)`只读取正式S2b，计算解析几何和保存复场指标；`reporting/render_two_dimensional_closeout_s2c_a.m(mat_path,new_dir)`只读MAT生成11组M00–M03 FIG/PNG；`validation/support/audit_two_dimensional_closeout_s2c_a.m(out)`从CSV独立复算/重绘/检查图册链接，拒绝覆盖已有审计目录。正式包`20261023_s2c_a_1642`为`S2C_A_DIAGNOSTIC_PROTOCOL_COMPLETE`，见[中文图册](../results/comparisons/two_dimensional_channel_closeout/20261023_s2c_a_1642/mechanism_atlas.md)和[唯一主报告](../reports/two_dimensional_channel_closeout_report.md)。新增声学场0；冻结协议本身不是求解许可。K_HIGH后来经单独授权完成S2c-B，其它候选/S3未授权。

- **S2c-B K_HIGH完整单案例：**`comparisons/prepare_two_dimensional_closeout_s2c_b_khigh.m`核验S2b/S2c-A、417点BASE_M和C2解析几何；`comparisons/run_two_dimensional_closeout_s2c_b_bie_job.m`先完成代表层，用户单独授权2400 s后完成余四层；`run_two_dimensional_closeout_s2c_b_pe_job.m`计算同一离散源、417点及13冻结探针；`run_two_dimensional_closeout_s2c_b_fullfk_job.m`计算完整双项解析平面外场。`validation/support/audit_two_dimensional_closeout_s2c_b_bie_reference.m`合成五层经验参考；`comparisons/analyze_two_dimensional_closeout_s2c_b_khigh.m`从保存场生成逐点/支持集/配对CSV；`reporting/render_two_dimensional_closeout_s2c_b_khigh.m`从保存MAT绘M04a–d并保留FIG；`validation/support/audit_two_dimensional_closeout_s2c_b_final.m`从CSV独立复算并像素重绘。以上运行都要使用已有MATLAB MCP会话与新run/目录，禁止覆盖正式包。先前P0/单层审计和1800 s停止记录是历史快照，当前证据见[图册](../results/comparisons/two_dimensional_channel_closeout/20261023_s2c_b_khigh_02/mechanism_atlas.md)及[主报告](../reports/two_dimensional_channel_closeout_report.md)。

- **S2b参考闭合正式入口与数据：**同一run_id `20261023_s2b_1305/formal_reference_v2`当前为`S2B_COMPLETE_REFERENCE_QUALIFIED`（声明的经验数值分辨率，不是物理PASS）。`comparisons/run_two_dimensional_closeout_s2b_bie_reference.m`在重新采样的真实物理RAM/同一累计7200s预算门下运行成本代表与当前代码BIE64/96再资格；禁止跳过其代表/预算门。`comparisons/complete_two_dimensional_closeout_s2b_reference.m`作完整支持的参考评价及解析平面外场+真实双项Gauss积分，原求解器/公共API不改，source/height/normal不变。`validation/support/closeout_s2b_fullfk_analytic_exterior.m`不是FastFK或局部截断，外部平面场解析精确计入，local_radius_m=Inf；`audit_two_dimensional_closeout_s2b_tail.m`为保存数组恒等式审查。`audit_two_dimensional_closeout_s2b_reference.m`与`audit_two_dimensional_closeout_s2b_resolution.m`为纯保存MAT/CSV/重绘/分辨率独立审计，禁止覆盖已有审计目录。所有数值作业仅MCP现有会话、显式输出路径、无cd。正式图/原始场/指标/成本/资格证据见[唯一主报告最新附录](../reports/two_dimensional_channel_closeout_report.md)，不自动开S3–S6/通信/Bellhop。下条早期受限流程、文件和“未完成”状态保留历史时间含义。

- S2b实际工作流：`comparisons/run_two_dimensional_closeout_s2b_p0.m`仅保存数据重算411点双权重及公式控制；`comparisons/prepare_two_dimensional_closeout_s2b.m`核验源哈希/保存密度/历史预算并冻结P1；`comparisons/execute_two_dimensional_closeout_s2b_job.m`串行执行单批（首批必须`bie3`），调用前从Windows Memory Available Bytes读取真实物理RAM并作为第三参数传入。禁止重复已有检查点、无BIE边界重解、无Bellhop/通信/PM。`reporting/finalize_two_dimensional_closeout_s2b.m`、`export_two_dimensional_closeout_s2b_sensitivity.m`只汇总保存复场；`render_two_dimensional_closeout_s2b.m`从保存MAT生成F01–F07 FIG/PNG；`validation/support/audit_two_dimensional_closeout_s2b.m(out,figure_folder,audit_folder)`从CSV独立复算并保存MAT重绘。所有输出目录须未使用；入口在setup后不改变cwd。实际包`20261023_s2b_1305`覆盖108/108、P0通过，最终图在`figures_v2`、审计在`audit_v2`；状态`S2B_DATA_COMPLETE_REFERENCE_LIMITED`，Goal未标完成，FullFK窗口和BIE证据限制见[唯一主报告附录](../reports/two_dimensional_channel_closeout_report.md)。不把调用纯绘图入口视为新场计算授权。

- `validation/prepare_pe_bellhop_bie_corrected_data.m`：从固定三案例复数图册及强海面
  411×45 原始场包读取保存入射场，在物理海面 `z=eta`、深度向下、`exp(-iωt)`
  约定下以 `exp(-2iketa)` 重算 PE 保存起始场到接收线和 45 个深度；BIE/Bellhop
  保留原生复场。生成新 `pe_bellhop_bie_corrected/<run_id>/` 数据包、逐案例复数
  CSV、二维复数 CSV、本轮/旧值指标及精确幅相平方分解。无新声源传播、BIE 求解或
  Bellhop 作业，不修改生产算子。`reporting/generate_pe_bellhop_bie_corrected_visuals.m`
  只读该包重绘 12 张 PNG/PDF 与中文报告，包内 `source_data/workflow_code/` 保存
  绘图代码快照。`validation/audit_pe_bellhop_bie_corrected.m` 独立复核 MAT/CSV、
  指标分解、源哈希和 12 张像素一致重绘。见[修正图册报告](../reports/pe_bellhop_bie_corrected_visualization_report.md)。

- `validation/audit_pe_bie_height_phase_sign.m`：读取保存入射场，以同一3 m FFT核诊断正负高度相位。
  原正号重现旧PE，负号与BIE原生比值比较；正入射平移平面给出解析符号控制。
  六个保存起始场评价，不启动BIE/Bellhop或上行传播，不修改生产代码。
  输出 `results/validation/pe_bie_height_sign_audit/<run_id>/`；见[镜像原因与旧指标更正](../reports/pe_bie_height_phase_sign_audit.md)。

- `comparisons/run_two_dimensional_closeout_s0_s1.m`：只读 E01/E02/E03 活动保存复数，冻结 3×3 坐标及来源，复核三档接收线并提取可比单点/点阵与原始场；不运行任何新物理场计算。
  `validation/support/two_dimensional_closeout_metrics.m` 实现逐点与幅相分解；`validation/test_two_dimensional_closeout_metrics.m` 用人工复数控制公式；`validation/support/audit_two_dimensional_closeout_s0_s1.m` 独立复核来源 CSV、汇总和逐像素重绘。
  `reporting/render_two_dimensional_closeout_s0_s1.m` 只读 `plot_data.mat` 生成七张 PNG/PDF。入口 `setup_vertical_project; r=run_two_dimensional_closeout_s0_s1('unique_run_id'); audit_two_dimensional_closeout_s0_s1(r.run_dir);`；报告见[二维收尾唯一主报告](../reports/two_dimensional_channel_closeout_report.md)，结果位于 `results/comparisons/two_dimensional_channel_closeout/<run_id>/`。S2b–S6 未执行。

- `comparisons/run_two_dimensional_closeout_s2a.m`：S2a 固定4 kHz、Gaussian源、z=3 m 的实际源/约定/平面控制；验证专用 `validation/support/run_pe_physical_surface_model0_validation.m` 将物理高度一次性适配到旧PE helper，再运行屏、传播与比值，包含`eta_z=0,±0.01 m`正入射独立路径控制。复用S0–S1、E10–E12和修正图册v2；预算内只作必要PE评价及完整双项FullFK平面、弱、中接收线积分，强档复用E12。`validation/support/solve_full_kirchhoff_2d.m`新增分项输出供闭合审计，方法本身不变。无BIE新求解/密度评价、Bellhop作业、整套x-z场或生产API修改。
  入口 `setup_vertical_project; r=run_two_dimensional_closeout_s2a('new_unique_id',available_physical_bytes); a=audit_two_dimensional_closeout_s2a(r.run_dir);`；先从Windows物理内存计数器取得可用字节数，并仅在已有MATLAB MCP会话运行。`reporting/render_two_dimensional_closeout_s2a(plot_data_path,figure_dir)`只读保存绘图数据重绘F02/F04/F05的PNG和可编辑FIG。输出及审计见[主报告](../reports/two_dimensional_channel_closeout_report.md)与`results/comparisons/two_dimensional_channel_closeout/20261023_s2a_final/`；S2b–S6未执行。

- 两套PE/Bellhop/BIE图册的[展示汇总](../reports/pe_bellhop_bie_visualization_summary.md)与
  [数据字段及只读/重绘说明](../reports/pe_bellhop_bie_visualization_data_guide.md)。
  汇总11张已有图，不新增计算；仓库内Obsidian阅读副本位于
  `results/visualization/pe_bellhop_bie_summary/20261019_presentation/`，不含数值数据副本。

- `validation/prepare_pe_bellhop_bie_xz_data.m`：固定强海面A=.20 m、4 kHz的真实二维
  411×45原始反射场。`'init'`保存原始数据和BIE密度；`'stage'`只写一个10001波束
  Bellhop网格输入；`'pebie'`从保存PE起始场传播并评价保存BIE密度，不重解边界方程；
  `'finish'`仅导入已有完成回执的Bellhop作业并核对3 m截线。
  `reporting/generate_pe_bellhop_bie_xz_visuals(run_dir)`独立重绘幅值、相位、幅相合成和
  原始剖面4张PNG/PDF；`support/audit_pe_bellhop_bie_xz`检查CSV、表示、颜色公式及独立重绘。
  输出`results/visualization/pe_bellhop_bie_xz/<run_id>/`。保留各方法原始源标度，
  不用比值冒充原场、不做幅相拟合；[范围与说明](../reports/pe_bellhop_bie_xz_visualization_report.md)。

- `validation/prepare_pe_bellhop_bie_visualization_data.m`：固定4 kHz、Tx100/Rx3、
  Gaussian sigma=.3、K=.10、A=.01/.05/.20 m的Bellhop–BIE/PE三方图册数据准备。
  `'init'`保存活动PE/BIE副本和输入哈希；`'next'`每次运行或接收一个独立Bellhop作业，
  flat/三个rough各5001/10001波束；`'stage'`只生成未开始作业的输入，允许独立Bellhop
  进程并行；`'finalize'`保存自包含复数MAT/CSV和配对指标。不得重跑PE/BIE或读取cash。
  `reporting/generate_pe_bellhop_bie_comparison_visuals(run_dir)`仅从数据包重绘7张PNG/PDF，
  可指定第二参数为新图目录。输出在`results/visualization/pe_bellhop_bie_comparison/<run_id>/`。
  `support/audit_pe_bellhop_bie_visualization(run_dir,rebuild_dir)`检查输入哈希、复数CSV、
  native表示映射、八作业诊断和7张重绘PNG像素；[图册报告](../reports/pe_bellhop_bie_comparison_visualization_report.md)。

- `validation/run_black_channel_basics.m`：仅 BLACK 数据下载、校验、官方解包/8 阵元短段 BPSK 回放及基础分析。
  `setup_vertical_project; r = run_black_channel_basics();`；`'prepare'` 仅下载检查，不加载 MAT；
  `run_black_channel_basics('render', r.run_dir)` 从保存结果重建三张图、CSV 与短报告。
  `support/black_*` 提供 HDF5 类型检查、时延矩、相位斜率、匹配滤波到达检查与解析合成控制；
  `support/prepare_black_channel_assets` 固定 Zenodo 21287414 和官方 MATLAB 提交。
  `reporting/generate_black_channel_basics_visuals` 与 `write_black_channel_basics_report` 仅负责展示。
  报告正文模板为 `reporting/templates/black_channel_basics_report_zh.md`，按“海试场景—数据含义—三图解读—建模边界”组织，
  完整参数/逐阵元指标保留于附录；生成器从保存结果替换数值，新增运行的工作流快照同时保存该模板。
  原始数据、源码快照、校验和和独立运行目录在 `results/validation/black_channel_basics/`；
  [报告](../reports/black_channel_basics_report.md)。仅通过 MATLAB MCP 现有会话运行，不改变当前目录、生产 API 或物理模型。
  Obsidian 独立阅读副本在 `results/visualization/black_channel_basics/20261019_obsidian/`，仅含一份 Markdown 与 `图片/` 下三张原图；不重新计算。

- `reporting/analyze_deepsea_dual_broadband_decision.m`：只读冻结三频MAT/CSV，分离平面—粗糙物理效应、参考敏感性、FK数值/模型差，输出成本情景及几何相位候选的探索性筛选；不调用散射求解器。
  `validation/preflight_deepsea_dual_broadband.m`锁定5/7kHz候选及30%物理RAM门；`validation/support/audit_deepsea_dual_broadband_decision.m`核对保存复数、成本与停止分类。入口拒绝覆盖当前证据，不能当作后续物理扫描授权。
  [物理决策](../reports/deepsea_dual_broadband_decision_2026-10-15.md)、[成本](../reports/deepsea_dual_broadband_cost_2026-10-15.md)、[协议](../reports/deepsea_dual_broadband_protocol_2026-10-15.md)；数据`results/validation/deepsea_dual_broadband_decision_20261015/`。本轮资源门拒绝，新增物理点0，决策B，仅已有4/6/8kHz可信；G1 OPEN/WP2=false。旧高阶BIE/预条件后续建议停止，不自动执行。

- `validation/validate_deepsea_bie_minimum.m`：三维辅助半空间层核、奇异/近奇异、独立Fourier制造逆问题；固定预算、不构建全局稠密阵、不求原PM。
  `support/deepsea_bie_*`为Gauss/Duffy、分块层核、平面trace/独立径向参考和真实图面已知Cauchy制造；`diagnose_deepsea_bie_windows`检查固定h窗口与图面PV。
- `validation/validate_deepsea_frequency_reconstruction.m`：保留旧257/513线性失败，对已知解析路径作时延分离，独立频点/深衰落/CIR对照；`deepsea_delay_separated_reconstruct`拒绝dual/PM认证用途。
  `audit_deepsea_bie_reconstruction`复核正负门和独立PDE/出射诊断；新unit为`test_deepsea_bie_reconstruction`。
  [BIE报告](../reports/deepsea3d_bie_minimum_report.md)、[重建报告](../reports/deepsea_frequency_reconstruction_report.md)、[协议](../reports/deepsea_bie_reconstruction_protocol_2026-10-14.md)。
  数据`results/validation/deepsea_bie_reconstruction_20261014/`。BIE未过参考门停止，B仅已知路径解析重建通过；G1 OPEN/WP2=false。


- `validation/validate_deepsea3d_preconditioner_decision.m`：右Neumann P=2I−B的小规模线性/显式矩阵诊断与70%可用RAM预检。
  `validation/validate_deepsea3d_preconditioner_checkpoint.m`：同一冻结B下flat→dual→quarter→pm的唯一干预，96内步/640 B调用/600s单案/3600s总核，拒绝覆盖。
  `support/solve_deepsea3d_right_neumann`与`deepsea3d_right_neumann_action`只属于验证；旧求解器不改。
- `validation/validate_deepsea_singlelink_readiness.m`：解析平面native/DSP/CIR与QPSK/AWGN/离散dual基础控制；`support/deepsea_readiness_guard`拒绝本轮数据进入统计/训练/PM BER。
  `support/audit_deepsea_two_track`重算新场原始离散残差并核对物理配置，`estimate_deepsea3d_bie_cost`仅非对角算术微基准与成本外推。
  [A报告](../reports/deepsea3d_preconditioner_decision_report.md)、[B报告](../reports/deepsea_singlelink_readiness_report.md)、[协议](../reports/deepsea_two_track_protocol_2026-10-14.md)。
  数据`results/validation/deepsea_two_track_20261014/`；A原PM失败并按规则停止，B确定性控制通过；G1 OPEN/WP2=false。无生产或旧基准更改。


- `validation/validate_deepsea3d_mapped_direct.m`：无Taylor的t=1矩阵自由mapped Helmholtz原型。
  setup后flat→dual→quarter→pm→pm_restart12；拒绝覆盖，PM要求controls已过门。
  `support/deepsea3d_mapped_action`、`deepsea3d_mapped_l0_solve`、`deepsea3d_mapped_discrete_residual`、
  `solve_deepsea3d_mapped_direct`复用冻结去混叠/模态helper；`certify_deepsea3d_mapped_geometry`给连续map界，
  `audit_deepsea3d_mapped_direct`重建三层门，测试`test_deepsea3d_mapped_direct`。
  [报告](../reports/deepsea3d_mapped_direct_report.md)、[协议](../reports/deepsea3d_mapped_direct_protocol.md)，
  数据`results/validation/deepsea3d_mapped_direct/`。controls通过，原PM短restart96内步仍未收敛；
  停止扩算，建议有界预条件改进。PM INCONCLUSIVE/G1 OPEN/WP2=false，无生产或旧基准改动。

- `validation/validate_deepsea3d_tfe_dealiasing.m`：同空间Fourier-Galerkin去混叠诊断，冻结PM不重生成。
  有效job为dual_fg、quarter_fg、pm_fg、pm_collocation_fixed；pm_collocation仅保存初始无效对照。
  support下`deepsea3d_fg_*`为插值/乘积/plan，`solve_deepsea3d_tfe_galerkin`为新验证分支，
  `observe_deepsea3d_fg_residual`独立真几何评价，`probe_deepsea3d_fg_geometry`只读保存场诊断，
  `audit_deepsea3d_fg_artifacts`核对空间/数据/失败分类，测试`test_deepsea3d_fg_products`。
  [报告](../reports/deepsea3d_tfe_dealiasing_report.md)、[协议](../reports/deepsea3d_tfe_dealiasing_protocol.md)，
  数据`results/validation/deepsea3d_tfe_dealiasing/`。原PM仍失败，停止扩阶；建议直接mapped Helmholtz。
  G1 OPEN/WP2=false，无support24/FullFK–PM或生产修改。

- `validation/validate_deepsea3d_tfe_prototype.m`：独立mapped Helmholtz TFE参考原型，冻结PM系数只读加载。
  setup后init，再run带job名：dual_base/dual_order/dual_order20/dual_horizontal/dual_vertical/dual_top/dual_domain；
  pm16/pm24/pm32/pm40/pm_scaled/pm_vertical/pm_top/pm_domain/pm_horizontal；status只读。
  最后`support/audit_deepsea3d_tfe_artifacts`重建数值门和成本，测试`test_deepsea3d_tfe`。
  support中的`solve_deepsea3d_tfe`及tfe_chebyshev/modal_plan/modal_solve/geometry/operators/forcing仅供验证。
  数据`results/validation/deepsea3d_tfe_prototype/`，[报告](../reports/deepsea3d_tfe_prototype_report.md)、
  [协议](../reports/deepsea3d_tfe_prototype_protocol.md)。已认证dual复现通过，原PM仍INCONCLUSIVE；
  不做support24/FullFK–PM，全G1 OPEN/WP2=false，生产/冻结/历史结果不变。

- `validation/validate_deepsea3d_pm_feasibility.m`：一个固定PM realization的有界独立参考可行性预检。
  setup后init、九次next、observe、trace、audit；status只读。support内的
  `build_deepsea3d_pm_target`冻结父面与系数，`observe_deepsea3d_pm_trace`直接指数边界诊断，
  `audit_deepsea3d_pm_feasibility_artifacts`核对数据/相位/失败分类；新测试`test_deepsea3d_pm_target`。
  输出`results/validation/deepsea3d_pm_feasibility/`，[报告](../reports/deepsea3d_pm_feasibility_report.md)、
  [协议](../reports/deepsea3d_pm_feasibility_protocol.md)。普通FE未过门，无FullFK执行，
  谱截止与非紧支撑影响未可信量化；单案例INCONCLUSIVE，全G1 OPEN/WP2=false。

- `validation/validate_deepsea3d_g1_frequency.m`：固定y方位q=.10/alpha2.5的4/8kHz独立验证，6kHz冻结复用。
  init后每频next至evaluation，evaluate检查M32/M40，observe通过才五次fk；
  若初始独立Gauss失败，只允许一次quadrature h.025，然后audit。各stage可传第二参数4000/8000。
  `support/audit_deepsea3d_g1_frequency_artifacts.m`生成原生三频摘要、数值变体范围、相消与分量趋势。
  数据`results/validation/deepsea3d_g1_frequency/`，[报告](../reports/deepsea3d_g1_frequency_report.md)、
  [协议](../reports/deepsea3d_g1_frequency_protocol.md)。4/6/8各点条件性PASS，
  6→8复误差趋势因阶数敏感性未分辨；全G1 OPEN/WP2=false，无生产或旧数据变更。

- `validation/validate_deepsea3d_g1_azimuth.m`：q=.10/alpha2.5固定海面的x→y方位对照。
  setup后init，重复next至observe/done，observe通过才五次fk，最后audit；status只读。
  与x方位相同预检/网格/域/真实边界规则；helpers
  `support/summarize_deepsea3d_azimuth_comparison.m`、`support/audit_deepsea3d_g1_azimuth_artifacts.m`
  分离方向性复场差、模型误差范围与截断阶数敏感性。
  输出`results/validation/deepsea3d_g1_azimuth/`，
  [报告](../reports/deepsea3d_g1_azimuth_report.md)、[协议](../reports/deepsea3d_g1_azimuth_protocol.md)。
  y方位单案例条件性G1 PASS；全G1 OPEN/WP2=false，海面/冻结求解器/旧验收数据不变。

- `validation/validate_deepsea3d_g1_q10.m`：alpha2.5的固定海面q=.10几何验证。
  setup后`'init'`，重复`'next'`至phase为observe/done，`'observe'`通过才五次`'fk'`，最后`'audit'`；
  `'status'`只读，init拒绝覆盖，低成本预检按需N20、细网格及域检查按需L48。
  `support/summarize_deepsea3d_q10_comparison.m`生成原生误差、相消诊断与数值变体范围；
  `support/audit_deepsea3d_g1_q10_artifacts.m`核对固定海面、唯一几何改动及CSV/相位/分类。
  输出`results/validation/deepsea3d_g1_q10/`；[报告](../reports/deepsea3d_g1_q10_report.md)、
  [协议](../reports/deepsea3d_g1_q10_protocol.md)。q=.10单案例条件性G1 PASS，
  全G1 OPEN/WP2=false；无海面平移、振幅扫描、PM/FastFK/SSP/动态或生产变更。

- `validation/validate_deepsea3d_alpha25_domain.m`：冻结alpha2.5的计算域检查与门控FullFK。
  setup后`'init'`→`'next'`→`'observe'`→`'reference'`；仅域门未通过时追加L48并重复三步；
  参考通过才执行五次`'fk'`，最后`'audit'`。`'status'`只读，init拒绝覆盖。
  `support/audit_deepsea3d_alpha25_domain_artifacts.m`核对CSV复数、原生相位恒等式及同一物理面。
  输出`results/validation/deepsea3d_alpha25_domain/`；
  [报告](../reports/deepsea3d_alpha25_domain_report.md)、[协议](../reports/deepsea3d_alpha25_domain_protocol.md)。
  本轮L24→32已通过，无需L48；反射/增量复误差2.959e-4/1.442e-4，
  alpha2.5条件性G1 PASS，仅此离散案例；全G1 OPEN/WP2=false，冻结代码与旧数据未改。

- `validation/validate_deepsea3d_alpha25_fine.m`：冻结alpha2.5的L24/h=.0625 N16/N20参考验证。
  setup后`'init'`、两次`'next'`、两次`'observe'`、`'audit'`；
  输出`results/validation/deepsea3d_alpha25_fine/`。
  [报告](../reports/deepsea3d_alpha25_fine_report.md)、[协议](../reports/deepsea3d_alpha25_fine_protocol.md)。
  细网格参考检查PASS，未扩域或调用FK；alpha2.5物理G1仍INCONCLUSIVE/全域OPEN/WP2=false。

- `validation/validate_deepsea3d_alpha25_preflight.m`：仅alpha2.5低成本FE稳定性预检。
  setup后`'init'`、重复`'next'`（本次3次）、`'cancellation'`、`'audit'`。
  观测helper为`observe_deepsea3d_fe_cancellation`，不修改冻结模型。
  输出`results/validation/deepsea3d_alpha25_preflight/`；
  [报告](../reports/deepsea3d_alpha25_preflight_report.md)、[协议](../reports/deepsea3d_alpha25_preflight_protocol.md)。
  预检支持进入受控正式检查；alpha2.5物理G1仍INCONCLUSIVE，全目标域OPEN/WP2=false。

- `validation/validate_deepsea3d_g1_amplitude.m`：冻结有限振幅三维基准的统一alpha扩展。
  setup后`'init'`，重复`'next'`至complete，再`'audit'`；分类测试保存要求见
  [报告](../reports/deepsea3d_g1_amplitude_report.md)与[协议](../reports/deepsea3d_g1_amplitude_protocol.md)。
  输出`results/validation/deepsea3d_g1_amplitude/`。alpha1复用；.5/1.5/2独立护栏通过，
  全目标域仍OPEN/WP2=false；冻结求解器与源/海面约定不变。

- `validation/validate_deepsea3d_g1_preflight.m`和`validate_deepsea3d_g1_finite.m`：
  一个Tx1500/Rx5、6 kHz、振幅(.02,.01) m紧支撑双方向面，独立高阶Dirichlet FE。
  先preflight，再`'init'`、十次`'next'`、`'offgrid'`、`'compare'`、`'cancellation'`、`'audit'`。
  新测试`test_deepsea3d_finite_reference`，保存规则见[协议](../reports/deepsea3d_g1_finite_protocol.md)。
  输出`results/validation/deepsea3d_g1_finite/`；[中文报告](../reports/deepsea3d_g1_finite_report.md)。
  单案例条件性G1 PASS，不批准WP2或全目标域；原FullFK及二维/生产实现不变。

- `validation/validate_deepsea3d_wp0_wp1.m`：Tx1500/Rx5、q=0/.10 独立三维WP0及
  紧支撑海面单Rx FullFK。setup后依次`'wp0'`、`'pilot'`、`'extend'`、`'reference'`、
  `'spectral'`、`'cancellation'`、`'audit'`；测试`runtests('test_deepsea3d_controls')`。
  helpers为`validation/support/deepsea3d_*`和`solve_deepsea3d_*`。
  输出`results/validation/deepsea3d_wp0_wp1/`；
  [协议](../reports/deepsea3d_wp0_wp1_protocol.md)、[中文报告](../reports/deepsea3d_wp0_wp1_report.md)。
  G0与原型数值门通过，独立参考仅一阶弱扰动，有限粗糙度G1开放。

- `validation/validate_fk_single_receiver_methods.m`：固定Tx100–Rx(0,3)的独立
  FastIn + FullOut与复合偏导核原型，4/6/8 kHz、12个PM输入及5个控制。
  `setup_vertical_project; validate_fk_single_receiver_methods();`；`'render'`
  仅刷新中文报告；`audit_fk_single_receiver_methods_artifacts()`审计保存证据。
  helpers位于`validation/support/fk_single_receiver_*`，生产PE及原基准不变。
  输出`results/validation/fk_single_receiver_methods/`，
  [方法与协议](../reports/fk_single_receiver_methods_protocol.md)、
  [方案对比报告](../reports/fk_single_receiver_methods_report.md)。
  补强后另输出`plan_controls.csv`、`cancellation.csv`、扫描CSV及
  `freeze_manifest.json`；补强前证据保留于同case`pre_hardening/`。
  [有条件冻结审核](../reports/fk_single_receiver_baseline_freeze_report.md)明确
  Rx0可冻结与相消点FastIn误差失败；严格匹配错误ID，不用无关异常过关。

- `validation/validate_bellhop_fig12b_three_attempts.m`：用户限定的三种追加
  原生方案（`CB R/CR R/CB X`），固定15例，`next/full/status/render`分段或
  复用；`reporting/generate_bellhop_fig12b_three_attempt_visuals.m`绘制保存复场。
  50/500/2000束、1.25/.625 m步长及均匀控制均已完成，仍未可信复现，已停止
  额外搜索。输出`results/validation/bellhop_fig12b_three_attempts/`；
  [最终停止报告](../reports/bellhop_fig12b_three_attempts_stop_report.md)。

- `validation/validate_bellhop_fig12b.m`：独立5000 m Munk/Tx4950 m/10 kHz
  图12(b)复现尝试。`setup_vertical_project; validate_bellhop_fig12b('full');`
  完成或复用34例；也可`init/next/status/render`分段。正负坐标均实际求解，
  50/500/2000束及扩展8000/10000束、步长、严格垂直、解析1/R控制均保留。
  `reporting/generate_bellhop_fig12b_visuals.m`只读本轮MAT绘图；
  `support/audit_bellhop_fig12b_artifacts.m`核对原始文件，文件检查通过不等于
  物理复现通过。输出`results/validation/bellhop_fig12b/`，
  [报告](../reports/bellhop_fig12b_reproduction_report.md)。原生点源r=0置零、
  全角声束未收敛、部分细步长触存储上限；未改Bellhop核心或通信链。

- `validation/validate_fast_fk_single_receiver_timing.m`：冻结已验证 canonical PM、
  4 kHz、4097点、rank3/12，严格FullFK单点(x0,z3)对当前FastFK整场取点。
  `setup_vertical_project; validate_fast_fk_single_receiver_timing();` 执行预热及
  9轮×10次计算，分别记录准备/计算；不重测Full整场、不改reference算法。
  原始数据在 `results/validation/fast_fk_single_receiver_timing/`，结论见
  `reports/fast_fk_single_receiver_timing_report.md`。本机计算阶段Full单点快
  36.57倍，含准备快1.925倍；此前整场速度比不能代表单点任务。

- `reporting/generate_fk_optimization_ppt_visuals.m`：从已有 Prototype、PM pilot、
  PM validity boundary 和 multi-frequency CSV/MAT 重绘 11 张 PPT 宽屏图。
  `setup_vertical_project; generate_fk_optimization_ppt_visuals();` 不运行新求解或计时。
  输出 `results/visualization/fk_optimization_ppt/`，讲述顺序与证据边界见
  `reports/fk_optimization_ppt_visualization_report.md`。

- `validation/validate_spectral_fast_kirchhoff_multifrequency_consistency.m`：
  固定 canonical seed260001、C2 海面与 Gaussian 源的 4–8 kHz 跨频验证。
  在 `setup_vertical_project` 后先 `'init'`，五次 `'next'`；只有冻结 3/8
  fidelity 失败时运行 `'rank'` 检查统一频带 rank，`'timing'` 重测选定 rank。
  三次 `'anchor'` 检查低/中/高频 BIE，三次 `'guard'` 检查 padding、插值和
  网格。仅全部 sparse 门通过后，65 次 `'dense'` 逐频 checkpoint；
  `'render'` 生成报告与图，`'status'` 查看状态。原始 3/8 结果和统一 rank
  诊断保留，不做逐频 tuning。输出目录
  `results/validation/spectral_fast_kirchhoff_multifrequency_consistency/`，报告
  `reports/spectral_fast_kirchhoff_multifrequency_consistency_report.md`。
  FullFK/BIE/FastFK reference 与生产代码保持冻结。
  `test_spectral_fast_kirchhoff_multifrequency_phase` 使用独立解析时延与
  固定横向波数控制；`support/audit_fk_multifrequency_phase` 处理原生复向量
  （禁止意外共轭转置），artifact helper 核对CSV复数符号、相位闭合和CIR。
  主H表为`dense_native_physical_H.csv`，主CIR为`reflected_cir.csv`；MAT的
  `phase_reference_meta`区分native physical、native reduced与DSP输入桥接。
  本轮8 kHz BIE超过MCP300秒等待，但共享会话原计算在481秒后完成checkpoint；
  遇到类似等待超时，应先核对checkpoint/会话，避免盲目重复昂贵任务。

- `validation/validate_full_kirchhoff_pm_validity_boundary.m`：固定 canonical PM
  realization 的单频物理适用域审查。先运行 `'init'` 保存受控 Fourier 设计并复用
  原有锚点；重复 `'next'` 每次只执行一个基准案例；六个初始 cutoff 完成后
  `'adapt'` 增加最多两个边界附近 cutoff；完成 `'next'` 后运行 `'plan'`，再重复
  `'refine'` 逐个执行 Full 积分、BIE ppw/order 和窗口护栏；`'render'` 从保存数据
  生成报告与图。入口在 `setup_vertical_project` 后从任意目录运行。
  结果：`results/validation/full_kirchhoff_pm_validity_boundary/`；报告：
  `reports/full_kirchhoff_pm_validity_boundary_report.md`。`'init'` 拒绝覆盖已有审查。
  `test_pm_validity_metrics` 检查原始相位误差、幅相误差恒等式及冻结离散源。
  本流程不运行 Fast FK 优化、multi-frequency、Monte Carlo 或生产集成。

- `validation/validate_spectral_fast_kirchhoff_pm_pilot.m`: 固定 PM 海面的单频
  FastFK–FullFK pilot。复用 canonical 与两个 Gaussian coefficient seed；保持
  共同 Fourier 模式，检查 physical cutoff、网格、padding 与少量独立 BIE。
  先测冻结 3/6 rank，失败时执行统一的有限 rank 诊断。运行 `'full'`，或分段
  `'fast'` / `'bie'` / `'raw'` / `'quadrature'`；`'raw'` 检查无包络 Fourier PM，
  `'quadrature'` 精化 high-K FullFK 的独立 BIE receiver 比较；
  `'render'` 只读取本轮 MAT。输出在
  `results/validation/spectral_fast_kirchhoff_pm_pilot/`，结论见
  `reports/spectral_fast_kirchhoff_pm_pilot_report.md`。新 helper 仅供验证；
  `test_spectral_fast_kirchhoff_pm_surface` 检查解析导数与 C2 嵌入。
  `assert_spectral_fast_kirchhoff_pm_pilot_artifacts` 核对保存数据、共同系数、
  固定 rank 与几何；算法 PASS 不替代整体 pilot 的物理门槛。

- `validation/validate_spectral_fast_kirchhoff_prototype.m`: 真正的 FFT 非局部
  Full FK prototype。全局入射 weighted basis + 出射 lag/height 分离 + 有限
  Toeplitz 线性 FFT 卷积；选定 `r_in=3,r_out=6,padding=2`，不建立 receiver×surface
  dense kernel。执行 `result=validate_spectral_fast_kirchhoff_prototype()`；
  `'render'` 仅从本轮 MAT 重绘。三个主案例、K=.70/sigma=1/Rx=2 控制、padding
  1/2/4、N=513/1025/2049/4097 的 fidelity/scaling 输出到
  `results/validation/spectral_fast_kirchhoff_prototype/`。方法与数值报告在
  `reports/spectral_fast_kirchhoff_prototype_method.md` 和
  `reports/spectral_fast_kirchhoff_prototype_report.md`。
  `validation/test_spectral_fast_kirchhoff_fft_controls.m` 是独立 flat Green、
  坐标导数、padding 与复插值控制；三个新增 support helpers 仅供验证。

- `validation/validate_full_kirchhoff_spectral_rank_audit.m`: 第一阶段 Full FK
  spectral separability audit。执行 `setup_vertical_project` 后调用
  `audit=validate_full_kirchhoff_spectral_rank_audit()`；复用保存参数/源谱，
  比较残余高度耦合的 weighted/unweighted SVD，并在三个真实曲面上做完整
  FK 运算 reconstruction。结果为
  `results/validation/full_kirchhoff_spectral_rank_audit/`，报告为
  `reports/full_kirchhoff_spectral_rank_audit_report.md`。不是生产 Fast FK；
  不运行 PE/BIE 或 PM Monte Carlo，不改变 cwd。

- `validation/`: smoke tests, invariants, communication-chain checks, and disabled-path regressions.
- `comparisons/`: surface-model, bubble-model, and coherent/incoherent comparison studies.
- `experiments/`: Monte Carlo runs, parameter sweeps, and calibration scripts.
- `reporting/`: figure generation, summary tables, visualizations, and report builders.
- `validation/support/`: Bellhop/Weyl/Li2009、properness、发布元数据和 artifact 检查等验证专用辅助函数；不属于公共 API。

## PE--Bellhop 展开坐标横向审计

本节中的早期PE–Bellhop、epsilon诊断和历史绘图入口已按用户要求移入
`../cash/pe_bellhop_validation_scripts_2026-09-01/`。完整文件清单及归档原因见
`../reports/pe_bellhop_validation_archive_manifest_2026-09-01.md`；归档内容不应从活动目录调用。

- `../reports/pe_bellhop_complete_validation_overview.md`: 当前 PE--Bellhop
  平面海面验证链的总报告，汇总早期原坐标、自由场、Gaussian 窗口、正式
  展开坐标、横向三方审计、参数扫描和坐标极限结果；旧的“完整汇总”仅作
  历史证据。
- `validation/validate_pe_bellhop_unfolded_flat_gaussian_vertical.m`: 正式的展开坐标平面海面 PE--Bellhop 对比，读取/生成 4--8 kHz Gaussian 直达、镜像反射、波束收敛和 PDP 结果。
- `validation/validate_pe_as_bellhop_transverse_vertical.m`: 只读取正式 MAT，在 4/6/8 kHz 的六个横向采样点上计算独立一步角谱，分离 PE--AS 与 Bellhop--AS 的误差来源，不重新运行完整 PE 或 Bellhop。
- `validation/validate_bellhop_transverse_parameter_scan_vertical.m`: 只运行 Bellhop，在 8 kHz 扫描步长、角扇区和 `.sbp` 采样密度，检查这些数值设置能否解释横向相位差。
- `validation/validate_bellhop_rotation_limit_vertical.m`: Bellhop-only 坐标极限审计，比较当前展开坐标、原坐标 `-89°` 和 `-89.5°` 的同源 Gaussian 直达/海面反射路径。
- `../cash/pe_bellhop_validation_scripts_2026-09-01/scripts/validation/validate_bellhop_rough_surface_near_vertical_limit.m`: **已归档**的历史500 m入口。
- `../cash/pe_bellhop_validation_scripts_2026-09-01/scripts/reporting/plot_bellhop_rough_surface_near_vertical_rays.m`: **已归档**的历史代表性绘图入口。
- `../cash/pe_bellhop_validation_scripts_2026-09-01/scripts/reporting/plot_bellhop_rough_surface_full_paths.m`: **已归档**的历史500 m全路径绘图入口；正常检索不读取`cash/`。
- `../cash/pe_bellhop_validation_scripts_2026-09-01/scripts/validation/validate_bellhop_100m_boundary_source_vertical.m`: **已归档**的100 m水侧边界源epsilon检查。
- `../cash/pe_bellhop_validation_scripts_2026-09-01/scripts/validation/validate_bellhop_100m_arrival_field_audit_vertical.m`: **已归档**的粗糙面arrival/C场诊断。
- `reporting/generate_bellhop_100m_source1m_visuals_vertical.m`: **当前代表绘图入口。** 用户改为实际离底1 m（Tx99 m、水深100 m、Rx3 m），停止epsilon诊断。独立A/E/R显示全部接收贡献及1001条发射声线；接收范围12 m，展示扇区范围120 m，绘图还原物理深度。`render_only=true`仅重绘本轮MAT，不运行Bellhop或读取退役数据。
- `reporting/generate_bellhop_pe_equivalent_ray_paths_vertical.m`: 从当前离底1 m结果中选择A模式唯一的0/0直达和1/0一次海面arrival，以其发射角单独运行两条R声线，并在最接近Rx的位置截断；不显示任何海底或高阶多次反射路径。结果位于当前案例的 `pe_equivalent_paths/` 子目录。
- `../reports/bellhop_100m_rebaseline_plan.md`: 当前 100 m 水深、物理源离底1 m的计划；旧零离底内容仅为历史记录。

正式结果位于 `results/validation/pe_bellhop_unfolded_flat_gaussian/`；横向审计结果位于 `results/validation/pe_as_bellhop_transverse/`，解释报告为 `reports/pe_as_bellhop_transverse_audit_report.md`。
Bellhop-only 扫描结果位于 `results/validation/pe_as_bellhop_transverse/bellhop_scan/`，报告为 `reports/bellhop_transverse_parameter_scan_report.md`。
坐标极限审计结果位于 `results/validation/bellhop_rotation_limit/`，报告为 `reports/bellhop_rotation_limit_report.md`。
旧 500 m 粗糙海面结果已按用户要求移入不参与正常读取的 `cash/` 隔离区，原报告路径仅保留迁移提示。新边界源结果位于 `results/validation/bellhop_100m_rough_surface/boundary_source/explicit_bty/`；父目录中的初次隐式平底诊断不能用作收敛证据。当前状态见 `reports/bellhop_100m_rough_surface_report.md`；归档清单见新计划，不需要打开归档目录。

边界源检查入口已归档（当前活动目录不再提供复现命令）；结果状态、参数和失败原因保留在
`../reports/bellhop_100m_boundary_source_report.md`及归档清单中。

## 运行方式

当前离底1 m代表图的已完成结果及限制见
`../reports/bellhop_100m_source1m_visuals_report.md`：193条ARR、199条E贡献、1001条
扇区；接收范围12 m、展示范围120 m，33条展示轨迹到达内部存储上限。
同一报告也记录只保留PE直达和一次海面反射拓扑的两条Bellhop轨迹图。

```matlab
setup_vertical_project
run('scripts/validation/validate_surface_wavefield_visualization_vertical.m')
run('scripts/comparisons/compare_specular_incoherent_surface_reflection_vertical.m')
run('scripts/experiments/sweep_monte_carlo_surface_channel_vertical.m')
run('scripts/reporting/plot_c35_core_heatmaps_vertical.m')
```

公共可复用实现已移到 `src/`。下文脚本名仍按 `validation/`、`reporting/` 等相对于 `scripts/` 的短路径描述。

## Kachoyan--Macaskill (1987) rough-surface reproduction

- `validation/validate_kachoyan_macaskill_1987.m`: validation-only implementation of the 1987 Fredholm collocation, Kirchhoff (kernel-disabled), propagating-mode reflection coefficients, energy balance, near-field reconstruction, Gaussian/PM fixed-seed substitutes, center-128 truncation, and Table-II sinusoid controls.
- Support solver: `validation/support/solve_kachoyan_macaskill_1987.m`.
- Results: `results/validation/kachoyan_macaskill_1987/`; report: `reports/kachoyan_macaskill_1987_reproduction_report.md`.
- Run with `validate_kachoyan_macaskill_1987(struct('run_mode','smoke'))` for a reduced check or omit the override for the paper-size run. The paper does not publish the original random arrays/seeds, so random-surface figures are explicitly approximate replacements.

## 2k 相位与海谱诊断

```matlab
run('scripts/validation/validate_2k_phase_approx.m')
run('scripts/validation/validate_2k_ocean_spectra.m')
```

- PNG/CSV：`results/validation/ssa_2k_phase/`
- 报告：`reports/validation_2k_phase_report.md`、`reports/validation_2k_ocean_spectra_report.md`

这两项是相位系数近似和最低阶 SSA/海谱趋势诊断，不是完整 PE、KStat、真实海洋散射或高阶 SSA 的验收。

## Raw-PM and Joint-Frequency K-Stat Prerequisite Validation

- `validation/validate_raw_pm_grid_coverage_vertical.m`: scans PM aperture, FFT-grid resolution, wind-speed coverage, and same-dx large-PM-grid to small-PE-window energy mapping.
- `validation/validate_kstat_joint_frequency_boundary_vertical.m`: boundary-only comparison of one shared explicit Gaussian PM surface, independent-frequency kstat, and covariance+pseudocovariance joint-frequency kstat at U=5 m/s and F=32.
- `validation/audit_pe_caching_vertical.m`: counts cacheable PE marches/FFTs, measures the unchanged public-path baseline, and estimates batched cache memory.

These scripts do not change the public defaults or communication chain. Detailed interpretation is in `reports/raw_pm_joint_kstat_prerequisite_validation_report.md`.

## Cached Joint-Kstat Receiver Validation

- `validation/validate_cached_joint_kstat_pe_receiver_vertical.m`: runs the U=5, F=32, 128/64 train/test reflected-only kdomain/independent/joint comparison through the independent cached PE executor.
- `validation/validate_cached_pe_public_consistency_vertical.m`: checks the cached executor against the unchanged public kdomain path with the exact same double-precision boundary input.
- `validation/analyze_receiver_properness_null_vertical.m`: calibrates finite-sample receiver `||P||_F/||C||_F` under a proper complex-Gaussian null model.

The receiver report is `reports/cached_joint_kstat_pe_receiver_validation_report.md`. The main result MAT and PNG diagnostics are under `results/validation/cached_joint_kstat_pe_receiver/`.

## U=5 F=64 Conditional Channel Generator

- `validation/validate_u5_conditional_channel_generator_vertical.m`: set `U5_CONDITIONAL_MODE=smoke` for 32/16 samples or `full` for 128/64 plus 10,000 generated channels.
- `validation/refresh_u5_f64_properness_vertical.m`: refreshes the F=64 three-mode null distributions and prescribed central-interval decisions.
- `validation/generate_u5_f64_sample_bundle_vertical.m`: saves generated H, total H, physical CIR, delay axis, labels, and timing without PE calls.
- `reporting/plot_u5_f64_conditional_validation_vertical.m`: produces QQ, correlation-matrix, eigenvalue, and LFM diagnostics.

Reusable interfaces are initialized by `setup_vertical_project` and implemented under `src/`（`properness_null_test_vertical` is validation support）: `estimate_conditional_channel_stats_vertical`, `sample_conditional_channel_vertical`, `build_channel_cir_vertical`, and `properness_null_test_vertical`. `build_physical_cir_vertical` remains only as the legacy common-time-shift wrapper. Full results are under `results/validation/u5_conditional_channel_f64/`; interpretation is in `reports/u5_conditional_channel_generator_f64_report.md`.

## Optimized Joint Builder and U=8 Node

- `validation/validate_streaming_joint_builder_u5_vertical.m`: checks the streaming/series F=64 builder against the original U=5 factors and cached-PE receiver statistics.
- `validation/audit_raw_pm_u8_aperture_vertical.m`: audits 100/256², 150/384², and 200/512² PM grids with 128 central-crop mapping realizations.
- `validation/validate_u8_conditional_channel_generator_vertical.m`: uses `U8_CONDITIONAL_MODE=smoke` for 16/8 samples or `full` for 128/128, U=8 properness, rank comparison, and 10,000 H+CIR samples.
- `validation/build_u5_u8_conditional_library_vertical.m`: creates and smoke-tests the exact-node U=5/U=8 library after both full models exist.

The reusable library interfaces are `build_conditional_channel_library_vertical` and `sample_conditional_channel_library_vertical`. They deliberately reject unsupported wind speeds and do not interpolate. The U=8 full 128/128 validation and two-node smoke test pass; see `reports/u8_joint_optimization_two_node_library_report.md` for measured acceptance and remaining high-wind memory limits.

## Two-Node Communication Validation

- `validation/validate_two_node_communication_vertical.m`: set `TWO_NODE_COMM_MODE=smoke` or `full`; compares fresh kdomain PE, joint cached PE, full-rank statistics, and 99.9% statistics at U=5/U=8.
- The full run uses 32 fresh channels per source, 8000 QPSK symbols per channel, Eb/N0 0:2:20 dB, channel-cluster bootstrap intervals, and a separate 128-pair full/low-rank test.
- `comm_main_vertical_psk` accepts optional external H(f) or h(t) through `COMM_EXTERNAL_CHANNEL_FILE`. Empty/unset preserves the original PE scenarios.

Reusable interfaces are `build_communication_taps_vertical`, `evaluate_mpsk_channel_ensemble_vertical`, and `sample_conditional_channel_rank_pair_vertical`. See `reports/two_node_statistical_channel_communication_validation_report.md`.

## PE Carrier-Phase Release Candidate Workflow

The formal workflow is destructive only in the narrow archival sense: the
preparation step moves registered phase-sensitive results to a timestamped
folder under `results/archive/pe_phase_reference_pre_rc/`. It never deletes
or overwrites them. Run the stages in this order:

```matlab
run('scripts/validation/prepare_pe_phase_reference_release_candidate_vertical.m')
run('scripts/validation/validate_pe_channel_phase_reference_vertical.m')
run('scripts/validation/validate_pe_phase_convention_uniform_vertical.m')

setenv('ADJOINT_PE_VALIDATION_MODE','full')
run('scripts/validation/validate_adjoint_pe_receiver_projection_vertical.m')
setenv('U5_CONDITIONAL_MODE','full')
run('scripts/validation/validate_u5_conditional_channel_generator_vertical.m')
setenv('U8_CONDITIONAL_MODE','full')
run('scripts/validation/validate_u8_conditional_channel_generator_vertical.m')
run('scripts/validation/build_u5_u8_conditional_library_vertical.m')

setenv('TWO_NODE_COMM_MODE','full')
run('scripts/validation/validate_two_node_communication_vertical.m')
run('scripts/validation/validate_public_channel_modes_vertical.m')
run('scripts/validation/validate_cached_pe_public_consistency_vertical.m')

setenv('PE_ATLAS_MODE','full')
run('scripts/reporting/generate_pe_propagation_atlas_vertical.m')
run('scripts/validation/finalize_pe_phase_reference_release_candidate_vertical.m')
```

`validation/audit_phase_reference_artifacts_vertical.m` is the read-only inventory
entrypoint and may be run separately. The active run metadata is stored in
`results/validation/pe_phase_release_candidate/current_run.mat`. A stopped
run may resume only when its code fingerprint is unchanged. If a
fingerprint-covered source changes, rerun the preparation step: it archives
the partial run and creates a new `run_id`.

The formal adjoint configuration is PE `128^2` / PM `256^2`; the U=8
conditional node intentionally uses its separately audited PE `128^2` / PM
`384^2` aperture. F=9 uses 4096 realizations, F=64 uses 512, and these runs
are expensive. The atlas must be last because it rejects validation inputs
whose `run_id` does not match the active release run. The finalizer reports
only `PASS`, `FAIL`, or `INCOMPLETE` and writes
`reports/pe_phase_reference_release_candidate_report.md`.

The completed reference run is `phase_rc_20260722_174945` and is `PASS`.

## Exact Adjoint PE Receiver Projection

Before the adjoint suite, the receiver carrier-reference integration can be
checked independently:

```matlab
run('scripts/validation/validate_pe_channel_phase_reference_vertical.m')
```

This reduced-cost audit uses F=65 for the unaliased 4 ms delay test, exercises
the four public surface branches, checks `legacy_reduced`, compares cached and
adjoint receiver outputs in both reduced/direct-DSP form, validates dense/FFT
`C/P`, and tests schema-1 conditional-model migration. Outputs are under
`results/validation/pe_channel_phase_reference/`.

- `validation/validate_adjoint_pe_receiver_projection_vertical.m`: validates the exact discrete conjugate transpose of the cached uniform surface-to-receiver PE, receiver projection, PM-to-PE embedding, dense/FFT receiver statistics, realization statistics, performance, and public regressions.
- Set `ADJOINT_PE_VALIDATION_MODE=smoke` for the reduced run or `full` for the accepted F=9/4096 and F=64/512 validation.
- Reusable validation interfaces are `apply_forward_surface_to_receiver_vertical`, `apply_adjoint_receiver_to_surface_vertical`, `build_adjoint_receiver_projection_vertical`, `run_adjoint_receiver_projection_vertical`, and `contract_kstat_receiver_stats_vertical`.
- `reporting/plot_adjoint_pe_receiver_projection_validation_vertical.m` regenerates the adjoint error, covariance/pseudo-covariance, PDP, LFM/matched-filter, augmented-eigenvalue, and timing/memory figures from a saved validation result.
- Outputs are written to `results/validation/adjoint_pe_receiver_projection/`; measured results and the integration decision are in `reports/adjoint_pe_receiver_projection_feasibility_report.md`.

This v1 path is limited to uniform sound speed, CPU double, fixed grids and frequency axis, one nearest-grid receiver, no bubbles, and no Doppler. Its PE operator remains a validation path; receiver outputs now use the same central phase-reference layer as the public API. The default surface model is unchanged.

## PE Propagation Visual Atlas

- `reporting/generate_pe_propagation_atlas_vertical.m` is the report-only entrypoint. It assembles the Tx-to-surface and surface-to-Rx center slices, surface and receiver planes, physical boundary branches, exact-adjoint sensitivity, accepted receiver `C/P`, wind-node statistics, and a signal-free carrier-reconstruction animation.
- `reporting/plot_pe_propagation_atlas_vertical.m` renders the fixed 14-item atlas. Spatial panels use a shared amplitude reference and a `[-50,0] dB` scale; phase is masked below `-40 dB`.
- Run `set PE_ATLAS_MODE=smoke` before MATLAB for the PE 64² / PM 128² / F=9 check. With the variable unset, the formal run reads the accepted PE 128² / PM 256² / F=64 adjoint and conditional-model results.
- Outputs are written to `results/visualization/pe_propagation_atlas/`, including a reusable MAT file, source/normalization manifest, numerical-closure summary, 13 PNG figures, and one MP4.

The atlas distinguishes physical boundary models from computational acceleration paths. `q` is an exact receiver-sensitivity kernel, not a physical reverse-propagated pressure field. LFM is applied only after obtaining `H(f)` and is shown without noise, modulation, synchronization, or equalization.

## PE/Bellhop Flat-Surface Cross-Validation

### Unfolded Gaussian-source comparison

- `validation/validate_pe_bellhop_unfolded_flat_gaussian_vertical.m` is the
  current validation-only entrypoint for the accepted unfolded-coordinate
  construction.  It maps the 100 m to 3 m vertical path to Bellhop ranges
  97 m (direct) and 103 m (image/reflected), then applies the flat
  pressure-release factor `-1` at the image receiver.
- `validation/support/write_bellhop_unfolded_gaussian_env_vertical.m` writes
  the per-frequency `.env` and `.sbp`; the `.sbp` pattern is derived from the
  production Gaussian angular spectrum and is not a pointwise receiver fit.
- `validation/support/run_bellhop_unfolded_gaussian_vertical.m` runs Bellhop
  and reads the standard `.shd`/`.arr` outputs.  The primary acceptance uses
  normalized reflected/direct responses; absolute source pressure is an
  independent Weyl-reference diagnostic.
- Outputs are written to
  `results/validation/pe_bellhop_unfolded_flat_gaussian/` and the report is
  `reports/pe_bellhop_unfolded_flat_gaussian_report.md`.

### Bellhop 2020 parametric internal-wall validation

- The authoritative implementation and status report is
  `../reports/bellhop_internal_wall_implementation_report.md`. It consolidates
  the feasibility, flat, tilted, sinusoidal, beam/frame, old PM, redesign and
  parameterization-fix stages. Superseded reports and obsolete one-off scripts
  are in the recoverable root `cash/` quarantine, not active project context.
- The internal-wall visualization scripts and generated figures/fields described
  in the next bullets are now archived under
  `../cash/bellhop_internal_wall_visualization_archive_20260909/`; they are
  retained for reproducibility only and are not active entrypoints.
- `reporting/generate_bellhop_internal_wall_visuals_vertical.m` is the active
  read-only renderer for the two requested products. It consumes the stored
  dense tilted (`r=100+0.05z`) and sinusoidal (`r=100-2sin(0.04z)`) internal-
  wall fields, draws eight incident/native-backward/inverse-mapped rays with
  arrows plus equal-aspect reordered-axis wall zooms with explicit wall normals, and writes dense 2-D coherent
  reflected-only TL backgrounds after `T^{-1}(r',z')=(200-r',-z')`. The SHD
  input contains only the transformed post-wall branch, not a total
  direct-plus-reflected field. It also renders a separate dense
  incident-only TL field from the official matched-halfspace `C*X` run.
  Outputs are under
  `results/visualization/bellhop_internal_wall_visuals/`; the dense source
  field retained for regeneration is under
  `results/validation/bellhop_internal_wall_visuals/`; raw case sidecars and
  duplicate dense inputs are archived under
  `../cash/bellhop_internal_wall_visualization_archive_20260910/`.
- `validation/validate_bellhop_internal_flat_wall_poc.m`,
  `validation/validate_bellhop_internal_tilted_wall_poc.m`, and
  `validation/validate_bellhop_internal_sinusoidal_wall_poc.m` are the retained
  flat/straight/curved regression entrypoints. Their isolated source overlays
  and runners remain under `validation/support/`.
- Historical reflection/enhancement/TL renderers and their superseded figures
  are retained only in the recoverable archives listed by
  `../reports/bellhop_internal_wall_visualization_archive_manifest.md`; do not
  use them as active entrypoints. The consolidated report and the active
  renderer above define the current visualization products.
- Curved overlays use ordered noncoincident parametric segments, a normalized
  hit tangent, a TOP normal reconstructed from that tangent, signed
  turning-angle/arclength curvature, and finite profile support. Internal-wall
  geometry does not use `dz/dr`, `Dss`, `Delta z/Delta r`, `1/Delta r`, range
  monotonicity or ATI-style infinite endpoint extension.
- `validation/validate_bellhop_internal_vertical_tangent_poc.m` is the focused
  circular-wall regression at exact `dr/ds=0` and non-grazing `u dot n=1`. It
  scans 65/129/257 samples and checks intersection, orthonormal frame, signed
  curvature, `RN`, p/q, one pressure-release phase jump, finite values and
  positive transformed range.
- `validation/validate_bellhop_internal_large_sinusoidal_wall_visual_case.m`
  runs the corresponding full `[-30,30] deg` Bellhop-only fan after a
  5001-ray analytic geometry preflight. The finite wall support is
  `z=[-75,75]` m, while the wall-neighborhood plot remains `z=[-30,30]` m.
  Numerical outputs are under
  `results/validation/bellhop_internal_wall_two_period_full_fan/`; the three
  authoritative visualization PNGs are written to
  `results/visualization/bellhop_internal_wall_visuals/`. The TL axes include
  the complete sampled wall support, leaving regions outside the SHD grid as
  background rather than extrapolating pressure.
- `validation/support/run_bellhop_internal_pm_wall_poc_vertical.m` and
  `validation/support/bellhop_internal_pm_wall_poc/` are retained as the
  generic parametric-wall runner/build. The fixed-realization convergence
  completed entrypoint is archived at
  `../cash/bellhop_internal_wall_superseded_20260902/scripts/validation/validate_bellhop_internal_pm_fixed_realization_convergence.m`;
  it wrote one master Fourier realization and interpolated `N=513/1025/2049/4097`
  profiles from that same realization. Outputs remain under
  `results/validation/bellhop_internal_pm_fixed_realization/`.
- The invalid strict-90-degree PM comparator and its old profile generator are
  archived. The fixed-PM result and the finite-angle native-ATI/local-wall
  covariance audit are Bellhop-only PASS results; the controlled PE comparison
  and reduced ensemble smoke are complete, while large-scale Monte Carlo remains
  out of scope.
- The completed Bellhop-only fixed-realization local covariance audit is archived at
  `../cash/bellhop_internal_wall_superseded_20260902/scripts/validation/validate_bellhop_pm_local_reflection_covariance.m`;
  its one-off runner is archived at
  `../cash/bellhop_internal_wall_superseded_20260902/scripts/validation/run_bellhop_pm_local_covariance_vertical.m`.
  It reused the same seed-260001 master samples in native C-ATI and a
  source-centered proper rotation at `phi=89` and `89.5` degrees, with
  `N=2049/4097` and three exact paired rays. Results remain under
  `results/validation/bellhop_pm_local_reflection_covariance/`; the detailed
  stage report is archived under
  `../cash/bellhop_internal_wall_superseded_20260902/reports/` and its valid
  conclusion is merged into the authoritative Bellhop implementation report.
- `validation/validate_bellhop_shd_receiver_range_pairing_vertical.m` is the
  permanent receiver-column regression. It requires exact 103 m rotated and
  97 m native matches, verifies that native total/direct use the same column,
  and keeps the 102 m rotated column as a negative control for the former
  approximately `-2.095 rad` indexing error.
- Current Bellhop internal-wall status is **PASS_WITH_LIMITS**: flat, tilted,
  sinusoidal, vertical tangent, fixed-realization PM density convergence, and
  the finite-angle local native↔internal Reflect2D covariance audit pass.
  Complete receiver-field equivalence and large-scale Monte Carlo remain out of
  scope; the controlled PE rough-PM comparison is indexed below.
- The PE rough-PM comparison execution chain is implemented as independent
  Stage 0A--0E, Stage 1A--1C, Stage 2 and Stage 3 validation entrypoints. The
  original scope and comparability matrix remain in
  `../reports/pe_bellhop_pm_comparison_design_audit_report.md`; the complete
  executable index and result interpretation are in
  `../reports/pe_bellhop_pm_complete_execution_summary.md`.
- Stage 0A canonical profile provenance is now implemented by
  `validation/validate_pe_bellhop_pm_canonical_mapper.m` with support helpers
  under `validation/support/`. It is a mapper-only audit (no PE/Bellhop run),
  uses the existing seed-260001 Fourier coefficients, and writes its report to
  `../reports/pe_bellhop_pm_canonical_mapper_report.md` and results to
  `results/validation/pe_bellhop_pm_canonical_mapper/`. Stage 0B and all
  subsequent planned stages are complete; see the complete execution summary.
- Stage 0B one-transverse-dimensional PE bridge is implemented by
  `validation/validate_pe_1d_validation_bridge.m` and
  `validation/support/run_pe_1d_surface_reflection_validation.m`. It passes
  independent 1-D exact-AS and production-PE `k_y=0` checks for flat and weak
  phase-screen cases. Results are under
  `results/validation/pe_1d_validation_bridge/`; report:
  `../reports/pe_1d_validation_bridge_report.md`. Stage 0C and all subsequent
  planned stages are complete; see the complete execution summary.
- Stage 0C flat source/normalization audit is implemented by
  `validation/validate_pe_bellhop_pm_stage0_flat_source.m`. It compares the
  unchanged 1-D PE bridge with the existing Bellhop 2020 Gaussian `.sbp`
  internal-flat case using explicit 97/103 m SHD selectors and 5001/10001
  beams. It passes normalized offset-magnitude and axis-Q checks; the known
  ~0.26 dB backward-range influence offset remains diagnostic. Results are in
  `results/validation/pe_bellhop_pm_stage0_flat_source/`; report:
  `../reports/pe_bellhop_pm_stage0_flat_source_report.md`. Stage 0D and all
  subsequent planned stages are complete; see the complete execution summary.
- The pre-reflection incident-plane audit is implemented by
  `validation/validate_pe_bellhop_incident_field_vertical.m`. It compares the
  complete 1-transverse-dimensional PE incident field at 100 m with Bellhop
  2020 coherent `C *X` line-source output on the same transverse grid, using
  the saved spatial-phase sign and no surface/wall reflection. The default
  4 kHz run checks 5001/10001 beams and a half-dx receiver diagnostic. After
  the X integration all 13 gates pass; PE--Bellhop TL P95 over M95 is
  `0.0019036 dB`, versus `0.157943 dB` in the archived pre-X run. Outputs are
  under `results/validation/pe_bellhop_incident_field/` and the report is
  `../reports/pe_bellhop_incident_field_comparison_report.md`.
- Stage 0D constant-height sign audit is implemented by
  `validation/validate_pe_bellhop_pm_constant_height_sign.m`. It checks
  `eta0=+0.05, 0, -0.05 m` with image span `L=103-2*eta0`, and passes PE,
  Bellhop and cross-model phase-sign, path-time, pressure-release and rotation
  checks. Results are under
  `results/validation/pe_bellhop_pm_constant_height_sign/`; report:
  `../reports/pe_bellhop_pm_constant_height_sign_audit.md`. Stage 0E and all
  subsequent planned stages are complete; see the complete execution summary.
- Stage 0E numerical-budget freeze is implemented by
  `validation/validate_pe_bellhop_pm_numerical_budget.m`. It reuses one
  seed-260001 band-limited realization across PE window/grid/step and Bellhop
  profile/beam/step scans; all frozen numerical gates pass. Results are under
  `results/validation/pe_bellhop_pm_numerical_budget/`; report:
  `../reports/pe_bellhop_pm_numerical_error_budget.md`.
- Stage 1A Tier-1 fixed-PM reflected-only ratio audit is implemented by
  `validation/validate_pe_bellhop_pm_stage1_tier1.m`. Bellhop structural
  geometry/phase/beam-state checks pass. The active code now defaults to
  source geometry X and uses X-specific case roots; the historical report
  remains pre-migration provenance, while the executed X result is recorded
  in the source-geometry audit. The PE/Bellhop ratio difference is diagnostic.
  Results are under
  `results/validation/pe_bellhop_pm_stage1_tier1/`; report:
  `../reports/pe_bellhop_pm_stage1_tier1_report.md`.
- The point-source `R` versus line-source `X` covariance audit is implemented
  by `validation/validate_bellhop_source_geometry_rx_audit.m`. It changes only
  Bellhop `RunType(4)`, checks flat and weak-sinusoidal walls at 5001/10001
  beams, and conditionally runs only the 4 kHz seed-260001 Tier-1 case. At
  10001 beams, `X` removes the native-97 m/internal-103 m `-0.26065 dB`
  range-normalization artifact, while the Tier-1 PE--Bellhop delta remains
  `0.31037 dB / -2.24820 rad`. Results are under
  `results/validation/bellhop_source_geometry_rx_audit/`; report:
  `../reports/bellhop_source_geometry_rx_audit_report.md`.
- Integration and archive details are recorded in
  `../reports/pe_bellhop_line_source_integration_report.md` and
  `../reports/pe_bellhop_source_geometry_integration_archive_manifest.md`.
- Stage 1B production PE dimensionality sensitivity is implemented by
  `validation/validate_pe_bellhop_pm_stage1_dimensionality.m`; it copies the
  same fixed eta(x) across y and compares ny=256/512 against the 1T bridge.
  Results are under `results/validation/pe_bellhop_pm_stage1_dimensionality/`;
  report: `../reports/pe_bellhop_pm_stage1_dimensionality_report.md`.
- Stage 1C interpretation/freeze is implemented by
  `validation/validate_pe_bellhop_pm_stage1_interpretation.m` and classifies
  the 4 kHz result as PASS_WITH_MODEL_DISCREPANCY. The cross-model residual is
  reported separately from the smaller dimensionality sensitivity. Results are
  under `results/validation/pe_bellhop_pm_stage1_interpretation/`; report:
  `../reports/pe_bellhop_fixed_pm_4khz_comparison_report.md`.
- Stage 2 fixed-PM frequency extension is implemented by
  `validation/validate_pe_bellhop_pm_frequency_extension.m`. It runs the
  same realization at 4/6/8 kHz with PE flat/rough and Bellhop flat/internal-
  wall rough pairs; numerical geometry/phase/state checks pass and residuals
  remain model diagnostics. Results are under
  `results/validation/pe_bellhop_pm_frequency_extension/`; report:
  `../reports/pe_bellhop_fixed_pm_frequency_extension_report.md`.
- Stage 3 paired fixed-band PM ensemble smoke is implemented by
  `validation/validate_pe_bellhop_pm_ensemble.m` with seeds 260001:260008.
  It preserves the canonical per-mode spectral amplitudes and records paired
  PE/Bellhop TL, phase, power and native wall diagnostics, including
  machine-readable percentile summaries. Results are under
  `results/validation/pe_bellhop_pm_ensemble/`; report:
  `../reports/pe_bellhop_pm_ensemble_comparison_report.md`. The status is
  PASS_WITH_LIMITS because this is a phase-ensemble smoke, not an independent
  PM amplitude Monte Carlo.
- Stage 3B independent PM coefficient-amplitude ensemble is exposed by
  `validation/validate_pe_bellhop_pm_amplitude_ensemble.m`. It reuses the
  same fixed spectral-density band and canonical seed anchor, then draws
  zero-mean Gaussian cosine/sine coefficients with variance `S(k)*Delta-k`
  for the remaining paired seeds. The reduced 8-seed run uses 5001 beams for
  the flat baseline and every rough case, requires 5001/5001 wall hits, and
  passes all structural checks. Stage-2 reuse requires an exact embedded
  configuration match plus current artifacts; other cache reuse requires a
  full request fingerprint and verified `.shd`/`.iwdiag` output hashes. Its
  minimum accepted seed count is eight. Results belong under
  `results/validation/pe_bellhop_pm_amplitude_ensemble/` and its report is
  `../reports/pe_bellhop_pm_amplitude_ensemble_report.md`.

The complete stage index, frozen parameters, cross-model interpretation and
remaining statistical limitation are summarized in
`../reports/pe_bellhop_pm_complete_execution_summary.md`.

### Reflection-free four-level audit

- `validation/validate_pe_as_freefield_vertical.m` compares production
  multi-step marching with an independently coded one-step exact discrete
  angular-spectrum propagator. The sponge is exactly off for this hard gate.
- `validation/validate_bellhop_freefield_normalization_vertical.m` reproduces
  the matched-halfspace construction of Bellhop's official free-space point
  source example. It measures `|p|R`, spatial phase sign, the constant source
  phase, beam/step convergence, and zero-bounce arrival delay before applying
  the single global `1/(4*pi)` Green-function conversion.
- `validation/validate_pe_bellhop_freefield_vertical.m` is the orchestration
  entrypoint. It keeps explicit initial-plane and physical-source phase
  references, compares PE/Bellhop/analytic complex fields, saves arrival and
  convergence CSV files, and does not relax a failed tolerance.
- `validation/support/virtual_point_source_initial_field_vertical.m` is enabled only through the
  additive `source_mode='custom_field_fn'` validation hook. Public defaults
  remain `source_mode='gaussian'`.

```matlab
setenv('BELLHOP_EXE',fullfile('path','to','bellhop.exe'))
addpath('scripts/validation')
validation = validate_pe_bellhop_freefield_vertical();
```

Bellhop is an external dependency and is never copied into `results/`. Set
`BELLHOP_EXE` to the executable directly, or set `BELLHOP_TOOLBOX_ROOT` to
the external Acoustics Toolbox root; the shared resolver checks
`windows-bin-20201102/bellhop.exe`. If neither variable resolves to a file,
the validator must stop at its dependency check.

Artifacts are under `results/validation/pe_bellhop_freefield/formal/`; the
summary is `reports/pe_bellhop_freefield_validation_report.md`.

### Point-source error-budget follow-up

以下段落中的Bellhop重跑入口已移入
`../cash/pe_bellhop_validation_scripts_2026-09-01/`；对应命令仅作历史记录，当前活动目录不再提供这些入口。

- `validation/validate_pe_point_source_error_budget_vertical.m` runs the
  no-sponge fixed-`dx` window series, compares spatial truncation with a
  discrete spectral-cell Weyl initializer, validates a separate continuous
  Weyl integral, and then scans the complete window/thickness/strength matrix.
  Its CSV explicitly reports `dA=20log10(|H_sponge|/|H_no_sponge|)` and
  `dTL=-dA`, plus center (`rho<=2 m`), ratio-defined edge-band, and total
  terminal-plane energies.
- `validation/support/weyl_point_source_reference_vertical.m` removes the grazing `1/kz`
  singularity by separate propagating/evanescent substitutions; it is the
  reliable free-space reference.
- `validation/support/weyl_point_source_initial_field_vertical.m` is the discrete FFT initializer
  under test. Its failure to converge within the public 100 m window is
  reported, not hidden.
- `../cash/pe_bellhop_validation_scripts_2026-09-01/scripts/validation/rerun_pe_bellhop_point_source_postbudget_vertical.m` performs the
  final representative Bellhop rerun only after the continuous Weyl gate.

Outputs are in `results/validation/pe_point_source_error_budget/` and
`results/validation/pe_bellhop_freefield/post_error_budget/`.

### Production Gaussian window/sponge audit

- `validation/validate_gaussian_window_convergence_vertical.m` calls the
  production Gaussian path, performs the independent PE--AS hard gate, and
  establishes a no-sponge reference using fixed production sampling.
- `validation/validate_gaussian_sponge_vertical.m` evaluates the complete
  50/80/128 m x 6-ratio x 7-strength matrix with explicit `dA`/`dTL`, center,
  edge, total-energy, and edge-to-center diagnostic metrics.
- `validation/validate_gaussian_sponge_wideband_vertical.m` compares the
  3--5 kHz physical `H(f)`, unwrapped phase, and group delay for large
  no-sponge, production/default, and recommended cases.
- `reporting/generate_gaussian_sponge_validation_figures.m` emits the 15
  requested figures plus an optional noiseless LFM diagnostic and the final
  Q1--Q7 report.

Outputs are under `results/validation/pe_gaussian_window_sponge/`; the report
is `reports/pe_gaussian_window_sponge_validation_report.md`. Extended sponge
ratios and windows are validation-only opt-ins; normal public limits and
defaults are unchanged.

### Full reflected PE -> surface -> PE window audit

- `validation/validate_reflected_chain_window_convergence_vertical.m` creates
  one maximum-domain PM surface, crops it without per-window Hs rescaling,
  runs the adaptive 4 kHz no-sponge window sequence, and saves incident,
  reflected, receiver, center, and edge diagnostics.
- `validation/validate_reflected_chain_wideband_vertical.m` runs the gated
  3--5 kHz/33-point reflected and total-channel comparison. Per-case MAT
  checkpoints make the production-size run resumable.
- `validation/validate_reflected_chain_diagnostics_regression_vertical.m`
  checks seeded/override equality, diagnostics transparency, receiver-center
  consistency, the phase-screen invariant, and channel closure.
- `reporting/generate_reflected_chain_window_validation_figures.m` generates
  the single-frequency field and convergence figures; the wideband validator
  adds reflected and total-channel response figures.

Formal outputs are under
`results/validation/pe_reflected_chain_window/`; the reviewable conclusion is
`reports/pe_reflected_chain_window_validation_report.md`. The validators use
additive validation-only options and do not alter production defaults.

### Random-surface reflected-chain window robustness

- `validation/validate_random_surface_window_robustness_vertical.m` generates
  each 256 m master surface once, takes unscaled 160/192 m central crops, and
  runs the 4 kHz three-Hs-by-five-seed matrix with resumable per-case files.
- `reporting/generate_random_surface_window_robustness_figures.m` creates the
  gate scatter, edge-versus-response, pass-rate, and stage-energy figures.
- `validation/validate_random_surface_window_robustness_wideband_vertical.m`
  is a checkpointed 3--5 kHz follow-up for selected worst cases. The current
  formal run was stopped before all candidate--256 m pairs completed; do not
  use partial checkpoints as group-delay qualification evidence.

Formal 4 kHz outputs are under
`results/validation/pe_random_surface_window_robustness/`; the report is
`reports/pe_random_surface_window_robustness_report.md`. The result is
192.1875 m/no-sponge `15/15` strict passes, while 160.15625 m has `15/15`
response/center passes but `0/15` complete edge passes. No production default
is changed by these validators.

- `../cash/pe_bellhop_validation_scripts_2026-09-01/scripts/validation/validate_pe_bellhop_flat_surface_current_vertical.m` is the
  current formal entrypoint. It requires a stable absolute, non-Temp
  `BELLHOP_EXE`, records the binary/source SHA-256 values, runs the independent
  phase audit, 3/6/9 m paths, small-offset limit, source-aware diagnostic,
  sampling/aperture/sponge matrix, public regressions, Bellhop R/C/I fields,
  and the numbered ten-figure atlas.
- `../cash/pe_bellhop_validation_scripts_2026-09-01/scripts/validation/validate_pe_bellhop_flat_surface_vertical.m` generates and runs the historical standard Bellhop ASCII-arrivals case matched to the existing PE model in a uniform medium with a flat pressure-release surface.
- `validation/validate_pe_phase_convention_uniform_vertical.m` independently checks the PE reduced-envelope operator, longitudinal carrier sign, group delay, and validation-local FFT convention against one-step angular-spectrum propagation.
- `../cash/pe_bellhop_validation_scripts_2026-09-01/scripts/validation/validate_pe_bellhop_flat_surface_matrix_vertical.m` runs the historical 3/6/9 m analytic--PE--Bellhop timing/amplitude comparison with an open Bellhop fan, plus the C0--C4 PE grid/window/step convergence matrix.
- `../cash/pe_bellhop_validation_scripts_2026-09-01/scripts/reporting/generate_bellhop_flat_surface_visuals_vertical.m` is the archived historical matrix plotting entry.
- Formal mode requires `BELLHOP_EXE` even if another Bellhop is on the MATLAB
  or system path. Temporary paths and binary-hash changes are rejected.
- The current local baseline is OALIB `2020_11_4` at
  `<external>/AcousticsToolbox_2020/windows-bin-20201102/bellhop.exe`.
  The retained 2017 directory is historical and must be selected explicitly
  only when reproducing a run whose metadata records that binary.
- The validator runs scalar direct-only/direct-plus-reflection regressions and a 65-frequency PE case, compares direct/single-surface arrival times and TL, reconstructs matched PDPs, and checks public PE invariants.
- Outputs are written to `results/validation/pe_bellhop_flat_surface/`; the Markdown report records exact parameters, formulas, thresholds, results, and limitations.

This stage deliberately excludes rough-surface scattering, bottom bounces, stochastic channels, and communication processing. Validators now consume the public `H_*_reduced_f` and `H_*_physical_f` fields instead of applying an independent hidden carrier convention.

Historical formal command (archived; not runnable from the active tree):

```matlab
setenv('BELLHOP_EXE','<external>/AcousticsToolbox_2020/windows-bin-20201102/bellhop.exe')
addpath('scripts/validation')
% validate_pe_bellhop_flat_surface_current_vertical is archived.
```

The current run `bellhop_current_20260723_rc5` is `FAIL_CORE` with
`amplitude_status=OPEN`: all phase/path/delay/beam/public checks pass, while
the no-sponge aperture-only fields fail the `-40 dB` edge-validity prerequisite.
Outputs are versioned under
`results/validation/pe_bellhop_flat_surface_current/<run_id>/` and
`results/visualization/pe_bellhop_flat_surface_current/<run_id>/`.

The formal atlas contains 10 numbered PNGs plus MAT, CSV, manifest, and text
summary. Older three-figure outputs remain historical evidence and do not
override the current run ID.

Use environment variables already supported by individual scripts to reduce grid size, seed count, or output file names for quick checks.

## Li et al. (2009) Explicit Rough-Surface Validation

```matlab
run('scripts/validation/validate_li2009_explicit_surface_vertical.m')
```

This independent pure-acoustic workflow uses one shared raw-PM realization
across the full 12 kHz, 6 ms CW frequency synthesis, applies only the explicit
pressure-release Kirchhoff `2*k*eta` screen, and projects the reflected field
through a monostatic bottom--surface--bottom PE path. It excludes noise,
electronics, SSA, kstat, modulation, and BER. Outputs are under
`results/validation/li2009_explicit_surface/`; interpretation and unresolved
transverse-grid sensitivity are documented in
`reports/li2009_explicit_surface_validation_report.md`.

## PE--Bellhop PM model-discrepancy statistical study

The Stage-4 fixed-PM study is split into an execution entry and a read-only
postprocessor:

- `validation/validate_pe_bellhop_pm_model_discrepancy_statistics.m` runs the
  frozen Stage-3B 4 kHz seed ensemble and conditionally extends it when the
  prescribed 24->32 gates fail.
- `validation/postprocess_pe_bellhop_pm_model_discrepancy_statistics.py`
  performs no PE or Bellhop calls; it consumes completed per-seed metric CSVs
  and regenerates the M=50 tables, bootstrap intervals, correlations, outlier
  audit, figures, and flat-array `result.mat` after the user stopped the long
  extension at seed 260050.

The authoritative report is
`reports/pe_bellhop_pm_model_discrepancy_statistics_report.md`.  Current state
is `PRELIMINARY_MODEL_DISCREPANCY`: numerical/applicability guards pass, while
the 24->32 bootstrap half-width gate is narrowly above its engineering limit.

## Controlled PE--Bellhop rough-surface Goal workflow

- `validation/validate_pe_bellhop_controlled_comparison.m` is the sequential,
  validation-only driver for the frozen controlled-comparison Goal. Run
  `validate_pe_bellhop_controlled_comparison('p0')`, then `'stage0'`; later
  stages remain locked until the preceding stage passes. Stage 0 writes the
  authoritative flat closure under
  `results/validation/pe_bellhop_controlled_comparison/stage0/`.
- The workflow uses source `X`, coherent `C`, explicit receiver-range
  selection, the rebuilt validation-only internal-wall binaries, and the
  existing 1-transverse PE helper. It does not modify PE, Bellhop core,
  `Reflect2D`, or `InfluenceGeoHatCart`.
- P0 and Stage 0 are PASS. The original raw Stage 1 weak-sinusoid results at
  `A=0.01`, `0.005`, and `0.0025 m` remain permanently retained as FAIL
  diagnostics. The convention-fixed regression is PASS, and the read-only
  Stage 1Y theoretical closure is
  `CONVENTION_THEORETICALLY_CLOSED`; Stages 2 and 3 have now PASSed their
  solver/mapping/finite guards for all prescribed height and wavenumber points;
  see `reports/pe_bellhop_controlled_comparison_stage2_height_sweep_report.md`
  and `reports/pe_bellhop_controlled_comparison_stage3_slope_curvature_report.md`.
  Future runs are limited to 10,001 beams by explicit user override; see
  `reports/pe_bellhop_controlled_comparison_stage1_status.md`.
- Stage 1X is implemented by
  `validation/validate_pe_bellhop_controlled_comparison_stage1x.m`; it reads
  completed `A=0.01,0.005,0.0025` MAT files only and writes diagnostics under
  `results/validation/pe_bellhop_controlled_comparison/stage1x/`.
- The Stage 1X convention-only regression uses the explicit comparison field
  `G_BH_comparison=conj(G_BH)` and writes isolated outputs under
  `results/validation/pe_bellhop_controlled_comparison/stage1_convention_fixed/`;
  the original raw Stage-1 results remain unchanged.
- Stage 1Y report:
  `reports/pe_bellhop_controlled_comparison_stage1y_theoretical_convention_closure.md`.
  It closes the PE/Bellhop time-harmonic, coherent influence, SHD reader, and
  constant-height analytic sign convention without new solver calls or any
  per-case fitting. The fixed comparison definition is `B_abs=conj(B_raw)` and
  `G_BH_comparison=conj(G_BH_abs)=G_BH_raw`.
- Stage 2 height sweep is implemented as the `stage2` branch of
  `validation/validate_pe_bellhop_controlled_comparison.m`, with resumable
  case files under
  `results/validation/pe_bellhop_controlled_comparison/stage2_height_sweep/`.
  A=0.01 is reused from the convention-fixed result and A=0.02/0.05/0.10/0.20
  are completed at 10,001 beams.
- Stage 3 slope/curvature sweep is the `stage3` branch of the same driver. Its
  manifest freezes `A=0.02 m` as the largest Stage-2 Region-I height, then runs
  `K=[0.10,0.20,0.35,0.47] rad/m` with 10,001 beams. Per-case diagnostics and
  resumable results are under
  `results/validation/pe_bellhop_controlled_comparison/stage3_slope_curvature_sweep/`.
  Stage 3 is PASS.
- Stage 4 is the read-only `stage4` branch. It makes no PE/Bellhop solver call,
  deduplicates the shared `(A,K)=(0.02,0.10)` case, applies the frozen Region
  rules, and writes the sampled validity map under
  `results/validation/pe_bellhop_controlled_comparison/stage4_validity_map/`.
  Stage 4 is PASS; the reported transitions are sampled brackets only.
- Stage 5 is the read-only `stage5` branch. It selects the first non-I height,
  maximum sampled height, and maximum sampled curvature cases, then applies the
  frozen `rho_shape` and `E_aligned/E_G` attribution rule. It also records
  analytic stationary/specular geometry, phase gradients, `2k eta`, and the
  angle-aware phase diagnostic under
  `results/validation/pe_bellhop_controlled_comparison/stage5_phase_attribution/`.
  Stage 5 is PASS with all representatives classified as spatial distortion;
  it makes zero solver calls and does not change the production PE screen.
- Stage 6 is the `stage6` branch. It byte-copies and verifies the canonical
  seed-260001 Fourier coefficients, runs only the frozen 4 kHz X/C,
  10,001-beam fixed-PM case, reuses the Stage-0 flat denominator and masks, and
  writes `results/validation/pe_bellhop_controlled_comparison/stage6_fixed_pm/`.
  Stage 6 passes every numerical/geometry guard and is Region III; its
  historical axis values are sanity checks only.
- Stage 7 is the zero-solver `stage7` branch. It reads the authoritative M=50
  CSV/MAT without rerunning seeds, preserves its historical point-source R
  provenance, applies no R-to-X correction, and writes the interpretation to
  `results/validation/pe_bellhop_controlled_comparison/stage7_historical_m50_interpretation/`.
- The completed workflow is summarized in
  `reports/pe_bellhop_controlled_comparison_final_report.md`; final status is
  `PHASE_MECHANISM_IDENTIFIED`. The conditional BIE/BEM branch was not triggered.

## Independent Helmholtz BIE third-reference workflow

- `validation/validate_pe_bellhop_helmholtz_bie_reference.m` is the new
  validation-only entrypoint for the independent 2-D pressure-release
  Helmholtz reference. `R1` flat self-validation and `R2` weak-rough BIE
  convergence and `R3` weak three-way closure are complete and pass.
- Run from any MATLAB working directory after project setup with
  `validate_pe_bellhop_helmholtz_bie_reference('r1')`, then `'r2'`, `'r3'`,
  `'r4'`, and `'r5'`. Each stage reads and enforces the preceding
  authoritative gate.
- `validation/support/solve_helmholtz_bie_halfplane_vertical.m` implements the
  R0-frozen Dirichlet-half-plane-Green combined-layer equation, slow-rise
  finite section, panel Nyström quadrature, independent boundary residual,
  and receiver evaluation. It does not call PE or Bellhop physics code.
- R0 theory, signs, normal orientation, restrictions, and stop conditions are
  frozen in `reports/pe_bellhop_helmholtz_bie_R0_formulation_audit.md`.
- R1 writes validation artifacts under
  `results/validation/pe_bellhop_helmholtz_bie_reference/R1_flat/` and the
  authoritative summary to
  `reports/pe_bellhop_helmholtz_bie_R1_flat_validation_report.md`.
- R2 uses one fixed `A=0.01 m`, `K=0.10 rad/m` C2-tapered benchmark surface;
  it keeps the physical support fixed while separately refining space,
  close-panel quadrature, and the BIE window. Artifacts are under
  `results/validation/pe_bellhop_helmholtz_bie_reference/R2_weak_rough/` and
  the report is
  `reports/pe_bellhop_helmholtz_bie_R2_weak_rough_convergence_report.md`.
- R3 constructs that exact same continuous tapered surface for the existing
  1-transverse PE helper, the 10,001-beam validation-only internal-wall
  Bellhop runner, and BIE. It applies only the previously frozen Bellhop
  rough/flat ratio conjugation for both native Helmholtz results and performs
  no scalar fitting. Its report is
  `reports/pe_bellhop_helmholtz_bie_R3_weak_three_way_closure_report.md`.
- R4 (`A=0.05 m`) and R5 (`A=0.20 m`) repeat BIE spatial/window refinement
  and compare the same tapered surface in PE, Bellhop, and BIE. Both are
  complete and select `BELLHOP_CLOSER_TO_HELMHOLTZ_REFERENCE`; see
  `reports/pe_bellhop_helmholtz_bie_final_report.md`. Optional R6 and fixed PM
  were not run because they are unnecessary for the completed first-round
  adjudication.

## PE rough-surface operator improvement workflow

- `validation/validate_pe_bellhop_helmholtz_bie_reference.m` also accepts
  `'g0'` after the authoritative R5 gate. It runs the frozen
  `A=0.02 m, K=0.47 rad/m`, 10,001-beam three-way diagnostic.
- `validation/validate_pe_surface_normal_approximation_audit.m` is the G1
  read-only angular-spectrum mechanism audit.
- `validation/validate_pe_surface_kz_aware_operator.m` runs G2 with the
  validation-only `model1_kz_aware` switch.
- `validation/validate_pe_surface_angle_slope_operator.m` runs G3 with the
  validation-only `model2_angle_slope` switch.
- These switches live only in
  `validation/support/run_pe_1d_surface_reflection_validation.m`; its omitted
  or `model0_normal` setting exactly preserves the previous behavior.
- G0--G2 pass. G3 fails the high-K improvement gate, so fixed PM is locked.
  See `reports/pe_surface_operator_improvement_final_report.md`; the final
  conclusion is `NONLOCAL_EFFECT_REQUIRED`.
- `validation/validate_pe_strict_normal_sigma_sweep.m` performs the subsequent
  4 kHz Gaussian-width applicability audit without changing production PE.
  Run the `lowk` cases at sigma `0.3/0.5/1/2 m`, the `highk` cases at
  `0.3/2 m`, and the `strong` cases at `A=0.20 m, K=0.10 rad/m` with
  sigma `0.3/2 m`, then `aggregate`. It compares legacy `2*k*eta`, validation-only
  kz-aware PE, 10,001-beam X-source Bellhop, and the accepted Helmholtz BIE;
  artifacts are under `results/validation/pe_strict_normal_sigma_sweep/` and
  the report is `reports/pe_strict_normal_sigma_sweep_report.md`. The accepted
  status is `STRICT_NORMAL_APPROXIMATION_PARTIAL`. Optional sigma `4 m` is a
  diagnostic only because the fixed Bellhop receiver influence produces
  sparse-support `0/0` ratios, so it is excluded from hard gates.
- `validation/validate_pe_bie_error_decomposition.m` is a read-only diagnostic
  that reuses the authoritative R3/R4/R5/G0 PE--Helmholtz artifacts, recomputes
  only the missing validation-only kz-aware field, and separates normalized
  complex, magnitude, and wrapped-phase residuals. It writes
  `results/validation/pe_bie_error_decomposition/` and
  `reports/pe_bie_error_decomposition_report.md`; it does not modify
  production PE, BIE, or Bellhop code.
- `validation/validate_full_kirchhoff_bie_vertical.m` runs the independent
  2-D Full Kirchhoff line-integral diagnostic against the frozen R3/R5/G0
  PE--BIE cases. Its solver is isolated in
  `validation/support/solve_full_kirchhoff_2d.m`; outputs are written under
  `results/validation/pe_bie_full_kirchhoff/` and summarized in
  `reports/pe_bie_full_kirchhoff_validation_report.md`. The native
  `exp(-i*omega*t)` rough/flat ratio is conjugated exactly once to match the
  frozen PE-comparison representation already used by `G_BIE`.
- `validation/validate_full_kirchhoff_convention_audit.m` is the independent
  sign/convention audit for that comparison. It reuses authoritative BIE
  fields, compares internally consistent `z_s=+eta` and `z_s=-eta` charts,
  checks native-vs-native and comparison-vs-comparison fields, and plots the
  legacy phase residual against `+/-4*k*eta`. Outputs are under
  `results/validation/pe_bie_full_kirchhoff_convention_audit/`; the report is
  `reports/pe_bie_full_kirchhoff_convention_audit_report.md`.
- `validation/validate_full_kirchhoff_parameter_sweep.m` performs the
  validation-only deterministic sinusoidal height and wavenumber sweeps using
  the same 4 kHz source, PE receiver grid, Helmholtz BIE convention, and
  fixed Full-Kirchhoff native-to-comparison conjugation. It writes checkpointed
  MAT/CSV data and representative geometry, incident-field, receiver-field,
  error, and sweep figures under
  `results/validation/pe_bie_full_kirchhoff_parameter_sweep/`; the summary is
  `reports/pe_bie_full_kirchhoff_parameter_sweep_report.md`. The checkpoints
  permit continuation after a long MATLAB run without changing the numerical
  settings.
- `validation/validate_full_kirchhoff_reduced_operator_stage1.m` is the first
  validation-only reduced-operator audit. It compares the accepted
  Full-Kirchhoff field with the PE phase screen, consistently removes surface
  height, normal, or arc-length terms, and scans receiver-centered finite
  integration radii. It writes diagnostic MAT/CSV tables and figures under
  `results/validation/pe_bie_full_kirchhoff_reduced_stage1/`; the report is
  `reports/pe_bie_full_kirchhoff_reduced_operator_stage1_report.md`.
- `validation/validate_reduced_kirchhoff_adaptive_aperture_stage1.m` builds
  the deterministic finite-aperture and adaptive-aperture database from the
  accepted A and K sweeps. It searches `L=4,8,16,32,64,128,Inf m` against a
  fixed Reduced--BIE complex-L2 gate, verifies the `Inf` regression against
  accepted Full Kirchhoff, and writes CSV/MAT data, sparse A-K maps,
  relationships, and representative field/convergence figures under
  `results/validation/pe_bie_full_kirchhoff_adaptive_aperture_stage1/`.
- `reporting/generate_full_kirchhoff_operator_visuals.m` is a saved-data-only
  visualization workflow. It reconstructs the validated Full-Kirchhoff
  integrand density, asserts its sum reproduces the saved native field, and
  writes kernel/contribution, geometry/support, common-mask receiver-error,
  aperture, and all-metric sweep figures under
  `results/visualization/pe_bie_full_kirchhoff_operator_atlas/`.
- `reporting/regenerate_full_kirchhoff_figure_conventions.m` is a strictly
  saved-MAT post-processing workflow. It does not call PE, BIE, or
  Full-Kirchhoff solvers; it only refreshes the parameter-sweep and Gaussian
  sigma-sweep figures with the corrected depth-axis direction, LaTeX labels,
  and explicit sigma/model line conventions.
- `reporting/plot_gaussian_sigma_angle_sensitivity_atlas.m` is a saved-data-only
  paper/report plotting workflow for the 15 Gaussian-sigma PE Model-0,
  kz-aware Model-1, Full Kirchhoff, and BIE cases. It validates both MAT
  sources, preserves the original pairwise-mask metrics, recomputes a clearly
  labeled four-field common M99/-40 dB mask, and writes PDF/300-dpi PNG
  figures plus the complete CSV/MAT audit table under
  `results/visualization/pe_bie_gaussian_sigma_angle_sensitivity/`.
- `validation/validate_full_kirchhoff_gaussian_sigma_sensitivity.m` compares
  PE, Full Kirchhoff, and Helmholtz BIE for sigma=`0.3, 0.5, 1, 2, 5 m` on
  three accepted deterministic sinusoidal cases. It reports propagating
  angular-spectrum quantiles, receiver errors, source-beam figures, and
  theta99/error trends under
  `results/validation/pe_bie_full_kirchhoff_gaussian_sigma_sensitivity/`.
- `validation/validate_full_kirchhoff_gaussian_sigma_kzaware_comparison.m`
  reuses that accepted sigma sweep and adds the existing PE
  `model1_kz_aware` result. It compares Model-0, Model-1, Full Kirchhoff and
  BIE on the same sigma-specific footprint, computes Model-0/Model-1
  improvement and remaining-fraction diagnostics, and writes the master/detail
  CSV/MAT tables and figures under
  `results/validation/pe_bie_full_kirchhoff_gaussian_sigma_kzaware/`.
- `validation/validate_full_kirchhoff_receiver_depth_sensitivity.m` evaluates
  the accepted sigma=1 m Model-1, Full Kirchhoff and BIE comparison at
  receiver depths `1, 2, 3, 5, 10, 20, 30 m`, holding the x grid, M99
  footprint and weights fixed. It writes depth metrics, normalized E(z)/E(3 m)
  tables, receiver lines, shared-scale field maps, difference maps and the
  z=1 m Full-Kirchhoff sampling check under
  `results/validation/pe_bie_full_kirchhoff_receiver_depth_sensitivity/`.
- `validation/validate_receiver_depth_numerical_audit.m` performs the
  validation-only receiver-depth numerical audit using the saved sensitivity
  fields, fixed M99 mask/weights, flat control, common-convention raw-field
  diagnostics, angular spectra, and the stored z=1 m FK convergence check.
  Outputs are under
  `results/validation/pe_bie_receiver_depth_numerical_audit/`; report:
  `reports/pe_bie_receiver_depth_numerical_audit_report.md`.
- `validation/validate_incident_reflection_field_covariance.m` first checks
  the PE Model-1 incident field against the saved Full-Kirchhoff angular-
  spectrum source at the reference plane and the sampled rough surface, then
  localizes the reflected-field PE/FK/BIE differences at z=2 and 3 m in both
  receiver space and kx space. It reuses the receiver-depth MAT and writes
  `results/validation/pe_bie_incident_reflection_covariance/`; report:
  `reports/pe_bie_incident_reflection_covariance_report.md`.
