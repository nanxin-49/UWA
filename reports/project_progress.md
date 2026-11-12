<!-- GENERATED: edit reports/project_progress.json, then run render_project_progress_dashboard. -->
# 从这里接着做

垂直水声项目 · 仿真与验证进度

记录更新：2026-11-12；本次核查：2026-11-12。

二维高斯单频已收尾；活塞20/10/5 cm三频自由场quicklook完成，P0未执行。既有M1仍为资源门停止，不进入M2。

## 现在在哪里

这个项目想回答的是：海底附近发出的声音，经过直达和海面反射到达接收机后，会怎样影响通信？为此，一边需要可信的声场和信道，一边需要能检验信道估计与均衡方法的通信实验。

你离开前，已经完成二维单频比较，做过三维固定海面的独立频点检查，还在一个接收点上跑通了受控的宽带与四相位调制（QPSK，每个符号承载两比特）通信实验。因此，重新开始时有不少结果可以接着用。

P22 保存结果与新 Ubuntu 最小项目基线已经复核通过。现在优先补齐约 100 m 场景下的生产模型约定、窗口与步长收敛、宽带/CIR转换和粗糙海面适用范围；统计快速生成、时变海面、多普勒及1500 m扩展留到后续。

2026-11-10 M0：生产原生PE到DSP改为完整包络共轭，统计/迁移/外部桥接同步。P22平面波高度、Gaussian同源角谱、C/P、SSA、CIR及小型已知CSI无噪声通信通过；旧失败和历史结果保留。

 本轮追溯原始窗口证据并冻结M1：不重扫失败小窗口，采用192.1875m起点、独立AS参考及逐步边缘检查。

2026-11-11：本轮M1执行在MCP连接前置门停止，未取得窗口、网格、步长或吸收层新资格；没有合格生产推荐配置，不能进入M2。

用户随后希望回到二维单频工作线并改用活塞。当前已查资料、选定20 cm代表尺寸和严格二维条带方案；先核共源与4 kHz试点，再测新Ubuntu成本，最后安排4/6/8 kHz比较。下面接续步骤面向这个新任务，不接管其他窗口的M1。

同日后续M1：既有MCP版本/PID核实后，完成P4/P6和10组独立AS；已完成对照独立审计0失败。P8在RAM跌破启动可用量30%时中止，Z4/Z8和全部吸收层对照未派发，无生产推荐，不进入M2。

2026-11-11后续活塞图示：旧5 cm经GitHub和本地代码确认是圆形半径，本次以二维总宽10 cm映射半尺寸，并加总宽5 cm控制。20/10/5 cm同振速三频自由场已保存复Pa、图、独立检查；这不关闭P0，后续先共源压力/梯度/相位资格，再平面镜像。

先读下面七条工作线的第一句话，再看末尾的三个接续步骤。技术参数和原始证据都可以展开，不必先读懂阶段编号。

## 各条工作线

### 二维单频比较

状态：高斯基线已收尾 · 活塞自由场图示完成，P0未执行；结果日期：2026-10-26（高斯收尾）；2026-11-11（活塞数据）；2026-11-12（图件封存）；本次核查：2026-11-12。

固定4 kHz高斯比较已收尾；总宽20/10/5 cm活塞三频自由场图示已完成，尚无生产共源与P0资格。

**为什么做：**先在较容易控制的二维问题里，确认不同方法比较的是同一声源、同一海面、同一接收位置，再看它们算出的反射波差多少。

**做了什么：**PE、边界积分方程（BIE，独立数值参考）和完整 Kirchhoff 积分近似（FullFK）的同源比较已完成。平面和三档正弦海面的九点结果齐全，三个预定 PM 海面实现也完成了比较；相关表格、图和历史审计均有保存。 新增严格二维速度孔径积分三频自由场，保存同振速绝对复Pa、接收线、近场/整体场和方向图；独立积分/角谱与保存数据检查无失败。

**说明什么：**已有一套可复查的二维对照基础，也知道某些设置下模型差异会随接收位置变化。写研究基础时可以引用，但需要带上频率、声源、海面和接收范围。 20 cm轴向更强且随频率升高明显收窄，适合作待资格主候选；同振速比较并非等功率指向增益。

**还差什么：**历史5 cm是圆形半径，二维总宽10 cm只作半尺寸映射；另有总宽5 cm辅助控制。活塞尚缺梯度、PE/BIE/FullFK共源、生产相位与Gaussian回归，未做平面/粗糙反射、窗口资格或远端成本。既有S4新增PM与跨平台限制继续保留。

**下一步：**本轮停止。后续先P0核共源压力/梯度、单位、初场平面与相位、倏逝/掠射和Gaussian回归；之后才P1平面镜像及4 kHz试点。无自动推进、SSH或M1预算转用。

**技术范围：**固定二维、4 kHz、984 项 Gaussian 发射场；原九点与 411 点接收线分别统计。S4=CORE_COMPLETE，S3A=3/3，S3B 因预算未运行；PHI_PLUS90 的混合平台 BIE 求解差未直接定界。 活塞quicklook：c1500 m/s、rho1025 kg/m³、Tx100 m/Rx3 m、峰值振速1 mm/s，W20/10/5 cm；20 cm三频-3 dB全宽112/67/49°。

**已有结果：**正式结果、表格和图件入口存在。 活塞run 20261111_2d_strip_01：MAT、CSV、8张交付PNG/FIG、日志、SHA和检查点。

**历史审计：**S4 保存证据审查通过并保留限制；本次没有重做数值审计。 本轮仅核旧源代码/检查点中的半径定义与SHA，不重演旧声学资格。

**当前重跑条件：**专题入口和保存结果可查；重新求解仍要按原阶段、参考和资源门审批。

证据入口：

- [活塞代表参数及三频执行计划（待实现）](../reports/piston_source_2d_single_frequency_plan_2026-11-11.md) · `reports/piston_source_2d_single_frequency_plan_2026-11-11.md`
- [正式收尾报告](../reports/two_dimensional_channel_closeout_report.md) · `reports/two_dimensional_channel_closeout_report.md`
- [S4 结果索引](../results/comparisons/two_dimensional_channel_closeout/20261026_s4_final_closeout_01/result_index.md) · `results/comparisons/two_dimensional_channel_closeout/20261026_s4_final_closeout_01/result_index.md`
- [保存证据审计摘要](../results/comparisons/two_dimensional_channel_closeout/20261026_s4_final_closeout_01/s4_saved_evidence_audit.json) · `results/comparisons/two_dimensional_channel_closeout/20261026_s4_final_closeout_01/s4_saved_evidence_audit.json`
- [原二维研究工作台](../UWA_2D_research_review_workbench.html) · `UWA_2D_research_review_workbench.html`
- [二维活塞自由场图示与历史尺寸溯源](../reports/piston_source_quicklook_2026-11-11.md) · `reports/piston_source_quicklook_2026-11-11.md` · 历史版本 `a37fd75 + source snapshots`
- [活塞quicklook最终检查点（非P0资格）](../results/validation/piston_source_quicklook/20261111_2d_strip_01/checkpoint_final.json) · `results/validation/piston_source_quicklook/20261111_2d_strip_01/checkpoint_final.json`

### 二维宽带与多模型参考

状态：M1_INCOMPLETE · P8 RAM门停止；结果日期：2026-11-11（M1部分对照/资源停止）；正式S5/S6仍暂停；本次核查：2026-11-11。

M0封存保持。M1已完成P4/P6和10组AS，限定数值对照通过；P8因RAM门中止，步长与吸收层资格未完成，无生产推荐配置。

**为什么做：**要把声场变成可用于通信的信道冲激响应（CIR，声音在不同延迟上的复响应），频率之间必须使用一致的声源、相位和时间原点，还要知道反射参考是否可靠。几个单频点相近并不足以保证时域结果可信。

**做了什么：**10 月 27 日已有 129 点同源直达检查和弱正弦 P22 的 6/8 kHz 单点参考。10 月 31 日的独立试点完成 4062.5 Hz、4031.25 Hz 的五层 BIE 参考，并检查了事先锁定的一个插值预测。 本轮完成原生reduced/physical/DSP三表示修复，P22平面波eta=0/±.01m与Gaussian同源角谱、统计和CIR独立控制通过。 本轮对37项相关文件绑定当前SHA，读取旧MAT/CSV并核主海面裁剪；17项限定门通过。设计5次生产PE及10组AS以分离空间/步长误差。 2026-11-11静态复核M0产物181/181 SHA一致，所有清单内MATLAB代码一致；确认吸收层单方向三次增长、x/y相加。连接前置门失败已保存，未启动PE或AS。 同日MCP恢复后执行生产P4/P6与10组AS，参考自身边缘、三频空间对照与4/6kHz PE-AS通过；独立收尾81行比较与MAT/CSV/SHA审计0失败。默认关闭的观察接口逐步记录三段传播，64²开关输出差0，已有M0宽带H消费者回放BER/SER=0。 后续只读诊断核41项源码、run03内部209项产物SHA一致。系统可用RAM下降1.724GiB，但P8采样工作集1.149–1.281GiB且生涯峰值不增长，未证实泄漏。

**说明什么：**说明频率表示与局部参考已有实际进展，后续可以据此设计更有针对性的验证。两个新频点只是小范围试点，还不能替代整个频带。 允许单独安排M1；不是全频物理CIR或三维粗糙海面认证。 M1_INCOMPLETE是执行受阻，不表示物理门数值失败；原吸收层既未判合格也未判不合格。 本轮证明指定均匀/平面Gaussian离散锚点的已完成对照，不形成完整M1配置资格。8kHz AS通过不能代替未完成生产PE或步长实测。

**还差什么：**完整 65/129 点反射参考、九点覆盖、高频与时间窗检查仍缺。4031.25 Hz 的 PE 独立步长/网格敏感性未齐；FullFK 与 BIE 的差异低于约定分辨率，不能据此判优。原正式 S5/S6 仍暂停。 最终封存发现两频点检查点所指审计代码 SHA=2c95abb6… 的原文未定位；固定提交的通用审计脚本为 bd865b9c…，不能顶替完整旧最终审计的来源。6 个产物绑定及 287 条数据清单匹配，本轮新 MATLAB 限定数值重核有效。 旧结果缺当次源码SHA，边缘仅33/65个z采样；新平面4–8kHz网格/步长与全步边缘资格待执行。 本轮M1没有新MAT数据、三频误差、逐步边缘或吸收对照；当前会话版本、真实PID与可靠长作业追踪尚未确认。 当前P8仅直达1036/1036、入射977/1068、反射0/33平面；Z4/Z8及SP4/SP6/SP8均未派发。原吸收层既未判合格也未判不合格。 缺运行期系统占用归因及内部数组分配轨迹；优化尚未实施、没有实测节省或P8重跑。部分MAT没有中途psi_k，不能直接从976步续算。

**下一步：**先审阅P8内存方案：优先缩短无用数组寿命、缓存固定环境屏，再做保持原运算链和逐步诊断的等价/资源验证。当前不派发传播，M1仍INCOMPLETE；后续恢复需当前资源、预算、版本/PID和输入证据门。

**技术范围：**二维 Tx=(0,100 m)，P22=(0,3 m)，原弱 C2 正弦 A=0.01 m、K=0.1 rad/m。10-31 状态为 DATA_COMPLETE_REFERENCE_PASS_MODEL_RESOLUTION_LIMITED；全带重建、物理 CIR/PDP/BER 未生成。10-27 的 48 h 提案尚未获批，不能作为后续执行许可。

**已有结果：**10-27 和 10-31 的结果包、检查点与审计摘要存在。 11-10 合并镜像 234 条与增补 53 条 SHA 清单均核验一致。

**历史审计：**历史两个新频点参考通过但模型分辨受限；11-10 独立 MATLAB 保存标量/源数组/配置/相位/退出记录与预测公式重核通过，未重演原全部场合成与边界数组检查。

**当前重跑条件：**两频点原求解入口已从 24a05ba 恢复；BIE helper 少几项历史输入检查，已列差异但未替换/执行。远端布局、预测冻结及调度环境未恢复为可直接运行环境；正式流程仍暂停。 本轮MCP attachment失败，禁止将远程或新MATLAB进程作为替代；当前M1任务已授权，不重复申请每个子步骤。 当前MCP版本和PID已核实，原连接失败成为历史；本轮停止原因更新为P8 RAM门。原run02空diary与300s传输超时保留，经确认终态后以显式stdout捕获在新run03重跑P4且三表示圆盘差0。

证据入口：

- [正式流程暂停与门槛](../reports/two_dimensional_s5_s6_2026-10-27.md) · `reports/two_dimensional_s5_s6_2026-10-27.md`
- [10-31 两频点检查点](../results/validation/two_dimensional_s5_p22_dualfreq/20261031_01/checkpoint.json) · `results/validation/two_dimensional_s5_p22_dualfreq/20261031_01/checkpoint.json` · 历史版本 `24a05ba72b58042084417e54ca437f03adb988dd`
- [两频点独立审计摘要](../results/validation/two_dimensional_s5_p22_dualfreq/20261031_01/audit_dualfreq_complete_01/audit.json) · `results/validation/two_dimensional_s5_p22_dualfreq/20261031_01/audit_dualfreq_complete_01/audit.json`
- [两频点运行入口（已恢复）](../scripts/validation/run_two_dimensional_s5_p22_dualfreq_remote.m) · `scripts/validation/run_two_dimensional_s5_p22_dualfreq_remote.m` · 历史版本 `24a05ba72b58042084417e54ca437f03adb988dd`
- [P22 可复现基线与对应表](../reports/p22_reproducibility_baseline_2026-11-10.md) · `reports/p22_reproducibility_baseline_2026-11-10.md` · 历史版本 `747af7ee88abda0bd7f554ddbd940ad3bc5f8cb9`
- [约100m静态模型收尾计划](../reports/static_channel_model_100m_closeout_plan_2026-11-10.md) · `reports/static_channel_model_100m_closeout_plan_2026-11-10.md`
- [M0合同修复与验收报告](../reports/m0_phase_contract_closeout_2026-11-10.md) · `reports/m0_phase_contract_closeout_2026-11-10.md` · 历史版本 `3f7071d + working tree`
- [M0独立数值与MAT/CSV审计](../results/validation/m0_phase_contract/20261110_m0_native_dsp_v2_01/independent_audit_final/checkpoint.json) · `results/validation/m0_phase_contract/20261110_m0_native_dsp_v2_01/independent_audit_final/checkpoint.json`
- [M0封存与M1窗口/网格/步长预审设计](../reports/m0_seal_m1_preflight_2026-11-10.md) · `reports/m0_seal_m1_preflight_2026-11-10.md` · 历史版本 `M0 source 969127f`
- [M1历史原始数据只读审计](../results/validation/m1_preflight/20261110_saved_evidence_04/audit.json) · `results/validation/m1_preflight/20261110_saved_evidence_04/audit.json`
- [M1执行收尾：连接门停止](../reports/m1_static_convergence_closeout_2026-11-11.md) · `reports/m1_static_convergence_closeout_2026-11-11.md`
- [M1前置检查点](../results/validation/m1_convergence/20261111_m1_execution_preflight_01/checkpoint_preflight.json) · `results/validation/m1_convergence/20261111_m1_execution_preflight_01/checkpoint_preflight.json`
- [M1数值执行收尾与RAM停止](../reports/m1_static_convergence_qualification_2026-11-11.md) · `reports/m1_static_convergence_qualification_2026-11-11.md`
- [已完成对照独立审计](../results/validation/m1_convergence/20261111_m1_execution_03/partial_closeout_audit/checkpoint.json) · `results/validation/m1_convergence/20261111_m1_execution_03/partial_closeout_audit/checkpoint.json`
- [18项矩阵实际状态](../results/validation/m1_convergence/20261111_m1_execution_03/matrix_status.csv) · `results/validation/m1_convergence/20261111_m1_execution_03/matrix_status.csv`
- [P8内存诊断与无物理改动优化方案](../reports/m1_p8_memory_diagnosis_2026-11-11.md) · `reports/m1_p8_memory_diagnosis_2026-11-11.md`
- [P8只读内存诊断回执](../results/validation/m1_p8_memory_diagnosis/20261111_01/checkpoint.json) · `results/validation/m1_p8_memory_diagnosis/20261111_01/checkpoint.json`

### P22 宽带与通信控制实验

状态：局部通过 · 已知信道控制；结果日期：2026-10-30（P22历史）；2026-11-10（M0小型控制）；本次核查：2026-11-10。

历史P22已知CSI控制保留原资格；本轮新增修复后小型无噪声链路门，不能拿旧P22代替新物理真值。

**为什么做：**先把一条受控的“频响 → 时域信道 → 发射 → 加噪 → 均衡 → 误码统计”链路走通，排除表示或通信实现本身的问题。

**做了什么：**最初直接比较 65/129 点 IFFT 没通过。后来把分量按共同时间窗投影到 257 点网格，并与独立计算的 257 点 PE 对照，通过原有门槛，随后完成平面与弱正弦两种海面的 QPSK 实验。 11-10 最小入口已恢复；原检查点输入/输出 SHA 核验，保存频响共同窗重建及原 QPSK 固定配置重放均通过原审计。 M0发现线性卷积截尾与FFT逆的边界不匹配，保留完整尾部后，在相同新生产H上BER从.07421875变为0；解析taps和CIR门通过。

**说明什么：**共同窗表示与已知 CSI 通信链现在有当前可用入口、匹配原记录的数据和本轮独立复核；可以继续用作受控基线。

**还差什么：**接收端事先知道信道，尚未验证从导频估计信道的能力；也没有完成 PE/BIE/FullFK 全频带的通信影响比较。该实验不升级正式 S5/S6，更不能代替目标深海场景验证。 旧carrier-only DSP、CIR、统计和相位敏感通信结果须明确来源后另存迁移/重算，不能直接标为新合同通过。

**下一步：**基线整理已结束；后续如需导频估计或模型失配研究，另定范围、预算和验收门，不自动开始。

**技术范围：**P22=(0,3 m)，发射深度 100 m，平面/原弱 C2 正弦；共同 257 点、64 ms 表示。65/129 共同窗 CIR 差约 7.1e-6，小于 1e-3。QPSK、已知 CSI 块 FFT LMMSE 均衡，SNR=0/6/12/18/24 dB，每点 262144 bit。原始 IFFT 失败记录保留，不能写成原始方法后来通过。

**已有结果：**原共同窗和 QPSK 结果保持原位；新审计、图件、版本/哈希/日志保存在 p22_reproducibility_baseline/20261110_115037。

**历史审计：**原始 IFFT 失败保留。共同窗/QPSK 原审计通过；11-10 保存数据重审通过，均衡符号差 0，CSV 指标最大差 3.9080e-14。

**当前重跑条件：**三个关键入口和两个独立审计已最小恢复；审计支持显式新输出目录，原数值门不变。当前可复核保存数据；新性能实验未授权。

证据入口：

- [最初 CIR 失败检查点](../results/validation/pe_screen_p22_comm/20261030_01/checkpoint.json) · `results/validation/pe_screen_p22_comm/20261030_01/checkpoint.json`
- [共同窗通过检查点](../results/validation/pe_screen_p22_common_window/20261030_01/checkpoint.json) · `results/validation/pe_screen_p22_common_window/20261030_01/checkpoint.json`
- [共同窗独立审计](../results/validation/pe_screen_p22_common_window/20261030_01/audit.json) · `results/validation/pe_screen_p22_common_window/20261030_01/audit.json`
- [最终 QPSK 检查点](../results/communication/pe_screen_p22_qpsk/20261030_03/checkpoint.json) · `results/communication/pe_screen_p22_qpsk/20261030_03/checkpoint.json` · 历史版本 `24a05ba72b58042084417e54ca437f03adb988dd`
- [QPSK 独立审计](../results/communication/pe_screen_p22_qpsk/20261030_03/audit.json) · `results/communication/pe_screen_p22_qpsk/20261030_03/audit.json`
- [QPSK 运行入口（已恢复）](../scripts/validation/run_pe_screen_p22_qpsk.m) · `scripts/validation/run_pe_screen_p22_qpsk.m` · 历史版本 `24a05ba72b58042084417e54ca437f03adb988dd`
- [共同窗验证入口（已恢复）](../scripts/validation/validate_pe_screen_p22_common_window.m) · `scripts/validation/validate_pe_screen_p22_common_window.m` · 历史版本 `24a05ba72b58042084417e54ca437f03adb988dd`
- [P22 可复现基线与对应表](../reports/p22_reproducibility_baseline_2026-11-10.md) · `reports/p22_reproducibility_baseline_2026-11-10.md` · 历史版本 `747af7ee88abda0bd7f554ddbd940ad3bc5f8cb9`
- [M0合同修复与验收报告](../reports/m0_phase_contract_closeout_2026-11-10.md) · `reports/m0_phase_contract_closeout_2026-11-10.md` · 历史版本 `3f7071d + working tree`

### 三维固定海面与深海模型

状态：局部通过 · 目标海况待验证；结果日期：2026-10-27；本次核查：2026-11-10。

固定测试海面在多个独立频点上有条件通过，但还没有证明整个频带或目标随机海况都可靠。

**为什么做：**二维模型不能自动代表真实三维散射，需要对三维点源、二维海面和指定链路另建独立参考。

**做了什么：**固定 compact-dual 海面已有 4/5/6/7/8 kHz 锚点及后续独立频点检查。65 点预测网格又用事先锁定的四个新频点检查，四点都通过约定筛查门。

**说明什么：**固定案例的跨频预测有有限的独立支持，可以继续研究怎样少算一些频点。65 点是预测输出网格，并非 65 个独立参考都已算过。

**还差什么：**连续频带、时域窗和 CIR/BER 尚未认证。目标 PM 随机海面的独立参考仍未解决；Q1 BIE 路线因独立误差门失败已经停止，不应继续沿用早期“下一步修 BIE”的建议。

**下一步：**暂缓1500m深海分支。约100m公共三维模型按新计划独立核约定和数值误差；二维BIE证据不自动认证三维PM，目标PM独立参考缺口保留。

**技术范围：**固定 Tx=(0,150.5,1500 m)、Rx=(0,0,5 m)，C2 compact-dual 测试海面。四个新点为 4750/7250/4125/6562.5 Hz。G1_target_PM=OPEN，WP2_formal=false。此处 FE 指 Dirichlet 场展开，不是有限元。

**已有结果：**正式报告、四点原生结果、预测锁和检查点存在。

**历史审计：**四点筛查和保存数据审计通过，有条件参考的范围不变。

**当前重跑条件：**历史结果来自旧 Ubuntu R2023b；不能直接视为新 Ubuntu R2025b 已通过项目验证。

证据入口：

- [四个新频点正式报告](../reports/deepsea3d_65grid_four_holdout_2026-10-27.md) · `reports/deepsea3d_65grid_four_holdout_2026-10-27.md`
- [结果索引](../results/validation/deepsea3d_65grid_four_20261027/README.md) · `results/validation/deepsea3d_65grid_four_20261027/README.md`
- [Q1 BIE 停止决策](../reports/deepsea3d_bie_q1_decision_2026-10-26.md) · `reports/deepsea3d_bie_q1_decision_2026-10-26.md`
- [三维总研究计划](../reports/deep_sea_channel_research_plan.md) · `reports/deep_sea_channel_research_plan.md`

### 1500 m 目标深海几何与 PE 窗口

状态：暂缓 · 先完成约100m建模；结果日期：2026-11-01；本次核查：2026-11-10。

已经算清楚几条候选路径什么时候到；还没有验证生产 PE 在这个长距离、大偏移场景下能否算准。

**为什么做：**先明确目标场景里直达、海面和海底一次反射的距离与时延，再判断计算窗口够不够大，避免接收点已经落进吸收边缘却把结果当成物理变化。

**做了什么：**水深 1500 m、Tx1490 m、Rx3/6/10 m、偏移 0–200 m 的解析几何表和图已生成。海底一次反射的候选到达也补算过。另对生产 PE 的窗口、网格和海绵层限制进行了代码核查。

**说明什么：**后续传播验证有了独立的几何参照。Rx 深度 10 m 时，平面海面与平坦海底的一次反射几何到达重合；这只说明到达时间重合，声压强弱还不知道。

**还差什么：**400 m 窗口下，200 m 接收位置落入吸收层；保持既有分辨率的更大窗口又受到当前窗口/网格保护限制。实际内存、运行时间与长程复声场收敛未测量。

**下一步：**按用户2026-11-10要求暂不推进1500 m孔径及大偏移计算；原解析几何和限制保留，后续有需要再另立任务。

**技术范围：**等声速 c=1500 m/s 的镜像法，603 行几何采样；没有计算路径复声压、海底反射系数、粗糙散射或通信。600/800 m 只是候选窗口，不是已通过的配置。

**已有结果：**D/S 和 D/S/B 解析报告、CSV 和图件入口存在。

**历史审计：**解析代数、CSV 和图件检查通过；没有目标 PE 传播通过结论。

**当前重跑条件：**解析几何入口在当前工作树；生产 PE 目标数值试验尚未准入。

证据入口：

- [解析几何与 PE 孔径设计](../reports/rc1_target_deepsea_geometry_and_pe_aperture_design_2026-11-01.md) · `reports/rc1_target_deepsea_geometry_and_pe_aperture_design_2026-11-01.md`
- [海底一次反射几何补充](../reports/rc1_flat_bottom_single_reflection_geometry_2026-11-01.md) · `reports/rc1_flat_bottom_single_reflection_geometry_2026-11-01.md`
- [D/S/B 几何入口](../scripts/reporting/render_target_deepsea_dsb_geometry.m) · `scripts/reporting/render_target_deepsea_dsb_geometry.m`

### Bellhop 代表案例与平面对照

状态：图示已完成 · 对照有未过项；结果日期：100 m 代表案例报告未标日期；连通 smoke 为 2026-11-01；本次核查：2026-11-10。

当前代表图采用 100 m 水深、发射机离底 1 m；声线展示完成，但不能把画出来的路径当成幅相已经验证。

**为什么做：**借助独立的射线方法理解多途结构，并在可比较的平面条件下检查 PE 与 Bellhop 的差异。

**做了什么：**Tx99 m、Rx3 m 的 4 kHz 代表案例已有接收贡献图和完整发射扇区图。另有平面镜像展开对照的保存表：波束收敛、宽带幅相和时延等检查通过，空间相位及复场差仍有未通过项。

**说明什么：**已有直观的路径材料，也已经定位出对照尚未闭合的部分。连通 smoke 说明程序能运行；代表图和部分通过项各有自己的用途。

**还差什么：**未完成角度/宽带扫描，也未认证粗糙面全部幅相。平面对照的空间相位约 0.0612 rad，大于 0.05 rad 门；空间复 L2 约 0.0613，大于 0.02 门。

**下一步：**若要继续比较，先限定相同路径拓扑、声源和相位约定，解释已有失败项；不要拿含底高阶多途直接对照当前 PE 的直达加一次海面反射。

**技术范围：**当前安装：水深100 m，Tx99 m，Rx3 m；停止零离底 epsilon 扫描。代表案例193条 ARR、199条 E 贡献记录、1001条展示声线；这些数量不代表真实独立到达路径数。旧500 m材料已经退役，本页不引用归档内容。

**已有结果：**代表案例报告、图件及平面对照 checks.csv 存在。

**历史审计：**代表图检查完成；平面对照 12 项中有 2 项未通过，不能记为整体通过。

**当前重跑条件：**当前绘图入口存在，仍依赖外部 Bellhop 可执行文件；本次没有重跑。

证据入口：

- [当前安装与代表图报告](../reports/bellhop_100m_source1m_visuals_report.md) · `reports/bellhop_100m_source1m_visuals_report.md`
- [当前重建计划顶部修订](../reports/bellhop_100m_rebaseline_plan.md) · `reports/bellhop_100m_rebaseline_plan.md`
- [平面展开逐项检查表](../results/validation/pe_bellhop_unfolded_flat_gaussian/checks.csv) · `results/validation/pe_bellhop_unfolded_flat_gaussian/checks.csv`
- [当前代表图入口](../scripts/reporting/generate_bellhop_100m_source1m_visuals_vertical.m) · `scripts/reporting/generate_bellhop_100m_source1m_visuals_vertical.m`

### 运行环境与重跑条件

状态：新机最小项目基线通过 · 完整迁移未做；结果日期：P22 原结果至 2026-10-31；本机保存数据复核及新机最小基线 2026-11-10；本次核查：2026-11-10。

新机隔离快照已通过小型公共信道的跨平台数值比较和两项 P22 保存结果审计；范围以外的项目组件仍待检查。

**为什么做：**后续计算应能说明用哪份代码、哪台机器和哪套输入得到结果，避免把环境打通当成科研验证通过。

**做了什么：**11-10 先恢复 P22 五个缺失文件并独立复核保存数据；随后在新机隔离快照上传 23 个白名单文件，逐项 SHA 通过。三组 CPU double 小型信道案例与两项 P22 保存结果审计在 Linux/R2025b Update 1 通过，与本机 Update 6 数值对照通过；92 个回传文件 SHA 匹配。 M0新增独立run，保留SHA回执失败、通信边界失败、CSV舍入和hex导入失败，再在新子目录恢复/复核。条件模型schema3与phase_contract_version=2.0.0阻止旧表示静默混入；代码/结果SHA和完整日志可追溯。 M0源码已有969127f，轻量检查点/失败记录/哈希白名单封存；既有47份生成源码快照保留历史，新快照不再自动纳入。旧窗口源码可恢复至8946b8e，但不能伪补原run绑定。

**说明什么：**当前来源和输入身份可追溯；新机已取得限定于 128×128、4 kHz 直达/反射、4–8 kHz 九点接口及两项 P22 保存审计的最小项目环境与兼容性基线。

**还差什么：**/home/misc/cwx 仍非完整 Git 工作区；仅隔离快照有本轮 23 个文件。未查完整项目全部 toolbox、Linux MEX、外部程序、大网格或长时任务；九点不认证 CIR，正式 S5/S6 继续暂停。两频点旧最终审计代码 SHA=2c95abb6… 的原文仍未定位，不能用固定提交中 SHA=bd865b9c… 的通用脚本顶替；旧失败证据保留。

**下一步：**本轮停止。后续如获授权，再核完整项目目标工作区、代码与输入身份，并按拟执行任务设计兼容性、资源和数值验收；不自动迁移或启动科研重计算。

**技术范围：**本机现有 MCP R2025b Update 6，新增工作连同进度看板刷新保守计 110/180 s；新机 SSH BatchMode、R2025b Update 1、单计算线程经实测，七次远端 MATLAB 进程 204.02/600 s，最长 50.70/180 s。五个数值作业逐项过 RAM/负载/磁盘门；未访问旧主机、归档或他人目录。

**已有结果：**本轮隔离结果包保存逐文件来源清单、补丁差异、完整本机/远端输出与日志、配置及原生数值数组、逐作业资源回执、对照 JSON、上传及回传 SHA 和检查点；原 P22 结果包未覆盖。

**历史审计：**此前基础 FFT smoke 仅说明基础环境；本轮最小真实项目案例和 P22 原保存数据审计建立限定基线，不追认旧 R2023b 结果或两频点完整旧最终审计。

**当前重跑条件：**只可按本轮隔离快照、固定配置和资源门重建限定基线；更大项目范围仍需独立授权与资格检查。

证据入口：

- [新机交接与边界](../reports/new_ubuntu_ssh_matlab_handoff_2026-11-09.md) · `reports/new_ubuntu_ssh_matlab_handoff_2026-11-09.md`
- [新机基础 smoke 检查点](../results/validation/new_ubuntu_matlab_smoke/20261109_232540_new_ubuntu/checkpoint.json) · `results/validation/new_ubuntu_matlab_smoke/20261109_232540_new_ubuntu/checkpoint.json`
- [现行公共信道入口](../vertical_channel_model.m) · `vertical_channel_model.m`
- [P22 可复现基线与对应表](../reports/p22_reproducibility_baseline_2026-11-10.md) · `reports/p22_reproducibility_baseline_2026-11-10.md` · 历史版本 `747af7ee88abda0bd7f554ddbd940ad3bc5f8cb9`
- [新机最小项目基线报告](../reports/new_ubuntu_minimal_project_baseline_2026-11-10.md) · `reports/new_ubuntu_minimal_project_baseline_2026-11-10.md`
- [新机最小项目基线最终检查点](../results/validation/new_ubuntu_project_baseline/20261110_121507_new_ubuntu_minimal/project/results/validation/new_ubuntu_project_baseline/20261110_121507_new_ubuntu_minimal/checkpoint_final.json) · `results/validation/new_ubuntu_project_baseline/20261110_121507_new_ubuntu_minimal/project/results/validation/new_ubuntu_project_baseline/20261110_121507_new_ubuntu_minimal/checkpoint_final.json`
- [M0合同修复与验收报告](../reports/m0_phase_contract_closeout_2026-11-10.md) · `reports/m0_phase_contract_closeout_2026-11-10.md` · 历史版本 `3f7071d + working tree`
- [M0封存与M1窗口/网格/步长预审设计](../reports/m0_seal_m1_preflight_2026-11-10.md) · `reports/m0_seal_m1_preflight_2026-11-10.md` · 历史版本 `M0 source 969127f`

## 最近发生了什么

### 2026-10-26 · 二维核心正式收尾

把已有比较收成一个范围明确的结论；三个预定 PM 实现完成，新增五个因预算没有运行。之后可以复用成果，不必从二维起点重做。

### 2026-10-26 · 停止当前三维 Q1 BIE 参考路线

线性求解残差虽小，独立边界和场误差仍未过门。因此停止继续扩大这条路线，旧计划中的继续调参建议不再有效。

### 2026-10-27 · 正式宽带流程停下来等待审阅

直达与少数频点参考有进展，但全频带反射、九点和通信依赖还不齐。48 小时只是一份待批提案，没有开始全带计算。

### 2026-10-27 · 三维预测接受四个新频点检查

固定海面的预测在四个新点通过有限筛查。它支持继续研究稀疏预测，但没有把 65 点网格变成 65 个独立物理参考。

### 2026-10-30 · 先发现 CIR 比较失败，再修正比较表示

原始 65/129 点 IFFT 差异没有过门；共同窗、按分量和时延处理后，与独立 257 点 PE 对照通过原门槛。原失败记录保留。

### 2026-10-30 · P22 已知信道 QPSK 控制实验跑通

通信链获得受控仿真结果和独立审计。这里尚未从导频估计信道，也没有升级正式多模型宽带资格。

### 2026-10-31 · 补齐两个邻近频点的参考试点

4062.5/4031.25 Hz 的 BIE 参考完成，局部预测也接受检查；模型差异的分辨能力和全带覆盖仍是缺口。

### 2026-11-01 · 目标深海路径几何算清楚，窗口问题浮现

有了直达、海面和海底一次反射的几何到达表；代码检查同时发现 200 m 偏移与现有 PE 窗口/吸收层冲突，尚未进入传播试验。

### 2026-11-09 · 新 Ubuntu 通路打通，项目仍待迁移

连接、MATLAB 版本和独立 FFT 检查通过；项目目录还不是 Git 工作区。机器能运行基础计算与项目能复现是两件需要分别验收的事。

### 2026-11-10 · 建立统一进度首页，补上索引漏记的进展

把二维、三维、P22 控制实验、目标几何与环境放在同一条阅读路径中；登记缺失入口和持续更新规则。本次没有新增声场或通信计算。

### 2026-11-10 · P22 最小入口已恢复，保存数据已独立复核

恢复 5 个缺失工作流文件，隔离审计输出；23 条旧检查点 SHA 和 287 条两频点清单匹配，共同窗/QPSK 原门和两频点限定重核通过。首次回执脚本错误保留并计时，原证据不覆盖。仍缺全带多模型、PE 敏感性和新机项目资格；下一步只能审阅本基线，另行授权后再设计后续工作。

### 2026-11-10 · 两频点最终旧审计代码来源仍缺一环

封存时六项产物 SHA 匹配，但检查点未给路径的最终审计代码 SHA 与固定提交的通用审计脚本不同。已保留首次封存失败记录，不替换旧检查点或假装代码一致。本轮保存标量/配置/相位重核可用，完整旧最终审计仍不能宣称已复现；本轮停止，不从归档补齐。

### 2026-11-10 · 新机最小项目兼容性基线通过

在不覆盖 /home/misc/cwx 原文件的隔离快照中上传 23 项白名单并验证 SHA；本机与新 Ubuntu 三组 128×128 CPU double 公共信道数值逐元素通过原容差，两项 P22 保存结果原审计通过，92 个回传文件 SHA 匹配。只建立本轮小案例与保存审计资格；九点不认证 CIR，完整项目依赖与迁移仍缺，正式 S5/S6 暂停及两频点旧审计来源缺口保留。下一步须另行授权并按目标任务核代码、输入、资源与数值门，本轮停止。

### 2026-11-10 · 提交基线成果，优先收尾约100m静态模型

P22恢复入口与新机最小基线已提交为5d2e9f8。按用户要求暂缓1500m，接下来先核生产入口的声源/高度/相位，再做窗口、宽带CIR和粗糙海面范围验证。本次只是计划，不新增数值结论。

### 2026-11-10 · M0相位合同与最小控制闭合，停在M1之前

修复原生PE到DSP的全包络共轭及C/P、增广实统计、迁移、SSA和桥接，物理/统计/CIR/无噪声控制和独立审计通过。保留回执与CSV打包失败；通信截尾边界失败在相同H上修复为BER=0。旧结果未覆盖，未升级正式S5/S6；下一步仅可另行授权M1。用户要求停止main合并，已abort并保留compare原HEAD。

### 2026-11-10 · 封存M0，改用历史192m证据设计M1

确认969127f已包含M0，核SHA并追加轻量失败/验收凭据。只读旧窗口数据，确认192m合格范围、160m边缘失败、33/65深度采样及源码绑定缺口；冻结5个PE/10组AS的最小新矩阵和保守资源门。本轮无新传播，下一步须新任务授权M1。

### 2026-11-11 · M1完整执行在现有MCP连接门停止

已获5个无吸收层PE、10组AS和阶段门后3个吸收层PE的14400s任务授权。核M0产物181项及源码身份、相位v2和吸收层定义后，MCP返回failed to attach to MATLAB session，无法确认版本、真实PID和日志/状态追踪，故18作业均未派发。无新数值资格和生产推荐，M1_INCOMPLETE；保留错误与SHA检查点。下一步仅在现有MCP恢复且资源/版本/PID门通过后从P4继续，不启动M2。

### 2026-11-11 · 活塞比较先定代表参数与严格二维范围

用户暂无实物参数，委托查阅后自行选择。依据厂家半波长口径量级，取20 cm工程参照；首轮严格二维使用20 cm宽匀速条带、向上发射、同一振速，在4/6/8 kHz沿单接收线比较PE/BIE/FullFK。圆形活塞与无限条带分开解释。已保存共源、采样、窗口、参考与新Ubuntu分阶段并行计划；没有实现新源或运行仿真。下一步先共源控制和4 kHz试点，取得远端成本后再定正式范围；原Gaussian收尾和其他窗口M1任务状态与预算保留。本次现有MCP正式渲染成功，未核或恢复M1执行门。

### 2026-11-11 · M1恢复运行后在P8触发RAM预留门

现有MCP确认R2025b Update 6/PID28836，完成P4/P6与10组AS；三频AS空间和4/6kHz生产同格对照及MAT/CSV/SHA独立审计通过。首次P4空diary/300s超时保留并计时，新run显式捕获stdout重跑且数值差0；首版参考审计错误保留后修正，仅重审已有MAT。P8在入射976步因RAM713232384B低于769405747B停止，Z4/Z8与全部SP未派发。故M1_INCOMPLETE，无生产推荐、不启动M2；下一步仅在继续指令和资源复核后恢复未完成M1。活塞计划及其他工作流记录/预算保留。

### 2026-11-11 · P8确认RAM余量门触发，给出无物理改动优化方案

只读核源码/回执/SHA和当前workspace元数据，没有传播或生产修改。系统可用RAM下降1.724GiB，P8工作集采样变化不足以解释全部下降，未证实泄漏、具体系统来源仍缺证据。列出缩短数组生命周期、固定环境屏缓存和分块逐步统计的方案及验证门；节省量只是静态估算。M1仍INCOMPLETE，不自动重跑或进入M2。

### 2026-11-11 · 活塞三频自由场图示完成，先P0再平面

核GitHub/本地历史确认5 cm为圆形半径；严格二维改用总宽20/10/5 cm速度孔径积分，不复用圆形方向图。保存三频绝对复Pa、接收线与归一化方向图；独立角谱/求积/MAT-CSV/FIG重建检查无失败，仅支持本次可视化。20 cm三频全宽约112/67/49°，8 kHz宽角域出现深谷，大离轴可弱于小孔径。缺生产共源压力/梯度、相位、Gaussian回归与平面反射；本轮停止，下次先P0。未运行PE/BIE/FullFK/SSH，不改Gaussian或转用M1预算。

### 2026-11-12 · 活塞图件修正版检查并封存

继续完成前一日quicklook交付：保留试绘标题重叠、首版负增益截断和一次exportgraphics多容器异常；改用单布局专用标题行，8张交付PNG逐图查看、FIG重建36幅CData差0，24条接收曲线匹配，0失败。声压数据未重算，数值源码SHA未变。交付绘图加图审计50.688s；报告、MAT/CSV、源码快照、日志和SHA封存。仍缺P0共源压力/梯度/相位与Gaussian回归，之后才做P1平面；本轮停止，不自动推进。

## 重新开始先做什么

### 1. 下一任务先做P0声源资格

自由场quicklook已完成，20 cm保留为主候选，10/5 cm作控制；先核生产共源压力/梯度、单位/相位与Gaussian回归，再考虑P1平面镜像。

完成标志：独立源定义与梯度、各模型接口及相位一致，P0门通过后才可安排平面反射试点。

执行边界：本轮仅可视化，到此停止；P0/P1需后续任务范围，旧Gaussian窗口、掩膜、M1预算及正式S5/S6状态均不转用。

### 2. 再核新Ubuntu依赖并测实际成本

白名单隔离快照先与本机已接受的小基线对照，再测代表性8 kHz作业内存和耗时。先评估两个独立单线程进程，不直接开全矩阵。

完成标志：版本和输入SHA明确，跨平台差异有数值门，内存、耗时及许可足以支持具体调度方案。

执行边界：仅matlab-ubuntu-new与新任务目录；不得覆盖已有快照。并发总峰值至少保留当时可用RAM的30%；重任务预算须明确，断线先查原作业。

### 3. 最后安排三频单接收线比较

试点合格后，先比较4/6/8 kHz下平面与弱正弦的PE/BIE/FullFK，沿同一物理接收区间解释模型误差，再决定是否扩展到中强正弦及固定PM。

完成标志：六个环境—频率组合具有可信参考、沿线复场误差、成本和独立审计；未分辨差异按未分辨报告。

执行边界：不预设速度提升或模型胜负；三个单频不是宽带/CIR资格。圆形源三维辐射、正式S5/S6、统计、动态和1500m均另列后续。

### 4. 窗口M1停在RAM门，保留已完成证据

P4/P6及全部AS已通过限定审计，P8部分轨迹保留；步长与吸收层尚未取得资格。

完成标志：后续取得完整P8、Z4/Z8及阶段门后的吸收层证据，才可判M1最终资格。

执行边界：当前不自动重跑或进入M2；保留30%启动可用RAM余量、14400s累计与3600s单作业限制。该预算不转用活塞计划。

## 术语随手查

**PE / WAPE：**抛物方程 / 宽角抛物方程：本项目逐步推进声波的传播算法。它是一种近似方法，仍需检查数值和适用范围。

**BIE：**边界积分方程：在本项目的一些受控案例中充当独立数值参考；它也有离散误差，并非实测真值。

**FullFK：**完整 Kirchhoff 积分近似：保留完整积分项的另一种声场计算方法。完整指积分形式，不代表物理上精确。

**CIR / PDP：**信道冲激响应 / 功率时延分布：描述不同延迟上收到的复信号或能量。空间接收线图不能直接当作时域信道。

**已知 CSI：**接收端事先知道信道信息。适合先验证通信链和均衡控制，但不能证明导频估计已经实现。

**QPSK：**四相位调制：每个符号承载两比特，用来测试信道和均衡后的通信表现。

**PM 海面 / compact-dual：**PM 是按指定谱构造的随机海面；compact-dual 是这里固定的测试海面。固定测试海面通过不代表随机海况总体通过。

**S4 / S5 / S6：**二维专项的阶段编号：单频核心收尾 / 宽带与时域 / 通信。P22 的独立控制实验与原正式阶段需要分别记录。

**未分辨：**观测到的模型差小于预先约定的数值分辨能力。只能说本次比较分不清，不能说两种模型严格相同。

## 记录与核查范围

M1前置任务曾从compare/c90d0ca干净基线核M0的181项产物SHA，48项源码/文档42匹配、6项为既有导航漂移，MATLAB源码一致。新增活塞计划本轮只查官方资料、静态核必要接入点和生成视图，不重演该身份审计、不运行声学求解或MAT数值审计。当前工作树另有M1修改，未将其作为冻结活塞版本。 本次M1执行另核真实版本/PID、生产代码与相位v2，补默认关闭的每步验证观察；完成既有数据MAT/IEEE-754 CSV/SHA和参考/生产对照审计。未接管活塞计划或将其预算用于M1。 后续P8内存任务仅查保存资源轨迹、源码和当前元数据，给出未实施优化方案，不改变M1数值状态。 本次活塞quicklook只做二维Hankel孔径积分；独立角谱含传播与倏逝分量、加密求积、MAT/CSV与FIG数据重建。没有生产模型接入、PE/BIE/FullFK比较、SSH或后续阶段执行。

M1前置任务保留连接失败及18作业未派发的历史回执。本次活塞计划更新时，现有MATLAB MCP两次正式渲染成功，45项引用无意外缺失；最后修改接续步骤后的再次渲染一直未返回，已结束编排等待，未终止MATLAB进程，是否稍后完成未知。最终视图依维护说明用同一JSON/模板纯文本包装刷新。未核或恢复M1执行门；无新MATLAB、SSH、活塞仿真、物理/通信回归或浏览器视觉验收。 同日后续M1使用现有本机MCP R2025b Update 6、PID28836、CPU double单线程串行。P4/P6与10组AS完成；P8用时552.6459804s，在入射976/1067步因可用RAM713232384B低于保留769405747B而停止。全部失败、先导重跑和审计计入14400s账；未使用SSH/GPU/新MATLAB。本次看板正式刷新结果以最新render_receipt为准，未做浏览器目视验收。 内存诊断耗时1.529s，无新传播，单独记录且不改M1原封账。 活塞quicklook使用同一MCP R2025b Update6/PID28836，数值计算6.418s、独立数值审计1.123s；图件试绘失败记录保留，交付figures_delivered明确与试绘版区分，未重算声压。 2026-11-12继续完成图件封存：最终渲染与FIG审计50.688s，FIG审计12.202s、CData差0，8张PNG实际查看；一次图导出异常保留，未重算数值场。

本次核查代码基线：`compare / c90d0ca; M0 identity rechecked; M1 execution attachment gate stopped; M1 resumed in attached R2025b Update 6/PID28836, then P8 RAM-reserve stop; piston free-field quicklook based on a37fd75, local source SHA snapshot`。

唯一维护源：`reports/project_progress.json`。

以后影响进度的任务完成、失败或停止时，修改相关工作线，追加事件，再刷新页面和文字版。旧失败和停止原因保留，记录只说明已有证据，不授予计算许可。

只更新与当次任务有关的内容；没有状态变化时不制造进展。结果日期、核查日期和生成时间分别记录。MATLAB MCP 不可用时注明未执行检查，不改用其他 MATLAB 运行模式。

```matlab
setup_vertical_project; render_project_progress_dashboard();
```
