# PROJECT_CONTEXT.md

## Purpose
This file is the append-only project log and long-term model memory for the
MATLAB vertical underwater acoustic channel and MPSK communication project.
Future coding models should read the current-state index first, then consult
the dated log for provenance. Historical entries are retained even when
superseded; a newer explicitly marked entry is authoritative.

## Authoritative Current-State Index / 权威当前状态索引（2026-08-20）

### S3A 第二阶段：seed260002 三方法真实PM对照（2026-10-25；当前）

- 状态 S3A_SEED260002_PILOT_COMPLETE_REFERENCE_QUALIFIED，S3A整体 IN_PROGRESS_2_OF_3；260001原参考限定合格保持，260003未运行。冻结260002的12模同谱Fourier系数、C2窗口、4 kHz Gaussian源和417真接收点在Ubuntu R2023b完成五层BIE、PE-mode0/半步长/同谱双网格、完整双项FullFK/解析平面外场/GL8。九点BIE PM相对同源平面响应1.8201161101，PE–BIE复L2 0.002089934812，FullFK–BIE 0.000267327198；411线分别为1.4739405374/0.093966108956/0.000400121178，九点与线的分母独立。四项模型差均高于5倍经验数值分辨率；这是固定样本的数值资格，不是严格物理精度界或PM总体结论。[唯一主报告首节](reports/two_dimensional_channel_closeout_report.md)含三张真图、与260001配对、公式和限制；[检查点](results/comparisons/two_dimensional_channel_closeout/20261025_s3a_seed260002_01/checkpoint.json)与独立MAT/CSV/FIG/PNG/SHA审计存原证据。
- 八个声学核517.021011 s；全部远端成功/失败/启动进程保守扣记656/1100 s，加260001的667秒后S3A规划余1677/3000 s、共享累计4724.424104/7506.601 s，共享余2782.176896 s。260003原标准1100秒和额外133秒故障余量尚可保留，但必须另过输入、成本、RAM和参考门；未自动执行。G1_target_PM=OPEN、WP2_formal=false，PHI仍为 S2C_B_PHI_PLUS90_DATA_COMPLETE_REFERENCE_LIMITED。下条是第一阶段快照。

### Codex 执行上下文与短回执（2026-10-25；仅流程）

- 后续任务先读本索引顶部最新状态、实际授权阶段的契约与入口；旧日期段落按精确依赖定位，不整篇重复载入。短回执、日志留存、检查点、增量证据和可比较测量方法见 [执行协议](reports/codex_token_execution_protocol.md)。本条没有启动 MATLAB、重算、绘图或改变任何数值资格；当前 PHI 仍为 `S2C_B_PHI_PLUS90_DATA_COMPLETE_REFERENCE_LIMITED`，混合平台同网格边界求解差未独立定界；`G1_target_PM=OPEN`、`WP2_formal=false`，阶段与资源门保持原规定。

### S3A 第一阶段：seed260001 单种子 PM 三方法对照（2026-10-25；历史快照，结果保留）

- 唯一S3A新声学案例为保存`seed260001`；三个原PM输入的12模同谱/不同系数、C2高度/坡度/曲率、S2b冻结源和旧场缺口已逐值核验。远端Ubuntu R2023b在同一417接收支持完成五层BIE、PE-mode0/半步长/同谱1968网格、完整双项FullFK及解析平面外场。真实九点BIE PM对平面物理响应复L2 `1.7664`；PE–BIE `0.1099207929`，FullFK–BIE `0.0005607659`，各自经验数值分辨率 `1.13255e-5/1.20757e-7`，且均过5倍门。411线另列，不以线插值九点。[唯一主报告首节](reports/two_dimensional_channel_closeout_report.md)先展示中文真图并列明源、公式、分母、参考资格与物理限制；[单seed检查点](results/comparisons/two_dimensional_channel_closeout/20261025_s3a_seed260001_01/checkpoint.json)和独立保存数据/FIG像素审计为执行证据。
- 状态`S3A_SEED260001_PILOT_COMPLETE_REFERENCE_QUALIFIED`只覆盖本固定单样本，S3A总体`IN_PROGRESS_1_OF_3`。八个成功声学核`516.355606 s`，包括失败尝试/启动的所有远端进程保守扣记`667 s`；S3A规划余额`2333/3000 s`，共享累计`4068.424104/7506.601 s`、余`3438.176896 s`。260002/260003未开展，建议下一轮仅审议260002的≤1100秒子预算并保留260003≥1100秒和133秒故障余量。`G1_target_PM=OPEN`、`WP2_formal=false`，不自动开展S3B/S4/Bellhop/三维/通信。以下PHI状态仍是其自身有效结论，但其旧“尚未开始S3A”文字是本轮之前的快照。

### S2c-B PHI_PLUS90：M05 严格审计完成（2026-10-25；当前）

- M05c 的旧 PNG 像素差来自右上角自动纵轴刻度；绘图层固定刻度后，新正式 v2 图与独立重绘五张逐像素相等，[严格审计 v2](results/comparisons/two_dimensional_channel_closeout/20261023_s2c_b_phi_plus90_01/audit/final_independent_audit_phi_v2.json)通过数值、来源、真九点和图像门。声场及判据未更改。混合平台同网格边界求解差仍未独立定界，正式状态 `S2C_B_PHI_PLUS90_DATA_COMPLETE_REFERENCE_LIMITED`；成本和 S3A 保留余额不变。下条为图像修复前快照。

### S2c-B PHI_PLUS90：跨平台边界与图形对象补证（2026-10-25；当前）

- R2023b 用未重解的 Windows `main96/self64` 密度独立复算 408 个离节点边界点，残差向量与 Windows 原值最大差 `1.04e-16/8.09e-17`；接收 417 点重评价差仍为 0。M05a–e 的 FIG 全 34 轴、61 曲线、4 个九点图像对象独立相同，M05c 严格 PNG 像素差仍存在。附加评价进程保守扣记 115 s 后共享累计 `3401.4241040/7506.601 s`，余额 `4105.1768960 s`，S3A 3000 s 保留后余 `1105.1768960 s`。没有同网格双平台重解的直接差界，故仍是 `S2C_B_PHI_PLUS90_DATA_COMPLETE_REFERENCE_LIMITED`。[主报告](reports/two_dimensional_channel_closeout_report.md)和[独立审计](results/comparisons/two_dimensional_channel_closeout/20261023_s2c_b_phi_plus90_01/independent_audit.md)为当前证据。

### S2c-B PHI_PLUS90：SSH 续算真场齐全、参考资格受限（2026-10-25；边界复核前快照）

- 见[唯一主报告最新首节](reports/two_dimensional_channel_closeout_report.md)及[中文 M05a–e 图册](results/comparisons/two_dimensional_channel_closeout/20261023_s2c_b_phi_plus90_01/mechanism_atlas.md)。本地 R2025b 的 `main96/self64` 原检查点逐字节保留；远端 R2023b 在 `/home/misc/cwx/UWA_phi_plus90` 完成 `space14/window60/window70`、PE-mode0、完整 FullFK，并把 MAT/日志/SHA 增量回传。七作业 417 点声学实耗 `578.0285102 s`；加兼容包装保守扣记后共享累计 `3286.4241040/7506.601 s`、余额 `4220.1768960 s`，保留 S3A 3000 s 后余 `1220.1768960 s`。PHI 中心 PE–BIE 复 L2 `5.09698e-4`，整线 `0.0203315`，中心仍远小于整线，但中心较 BASE 增大约 16.83 倍。混合平台边界求解差没有独立上界，严格 M05c 像素复绘未过（FIG 数值数据相同）；状态 `S2C_B_PHI_PLUS90_DATA_COMPLETE_REFERENCE_LIMITED`，不自动开 S3A 或其它阶段。

### S2c-B PHI_PLUS90：两层真实 BIE 后 RAM 再阻塞（2026-10-25；历史快照）

- 用户确认内存释放后，现有 MATLAB MCP R2025b 会话通过初始物理 RAM 门；在唯一 PHI run 内按序完成 `main96/self64` 两层新 BIE，保存各自完整密度和同 417 点原生复场，实耗 `203.6039647+153.5603395=357.1643042 s`。第三层 `space14` 前 RAM 为 `1,863,331,840 B`，一次有依据复核为 `1,879,687,168 B`，均低于保留 30% 所需 `2,026,747,612 B`，第三层未启动。共享累计 `3042.8783684/7506.601 s`，保留 S3A 3000 s 后 PHI 剩余额度 `1463.7226316 s`，保守剩余计划 `1252.49592096 s` 仍在预算内。两层自求积变化只能作部分数值诊断，五层经验参考与 PE/FullFK 尚无，M05b–e 和中心误差结论未得出。状态保持 `S2C_B_PHI_PLUS90_RESOURCE_REFERENCE_BLOCKED`；[2026-10-25 执行回执](results/comparisons/two_dimensional_channel_closeout/20261023_s2c_b_phi_plus90_01/execution.md)及[资源/哈希回执](results/comparisons/two_dimensional_channel_closeout/20261023_s2c_b_phi_plus90_01/P2_recovery_gate_20261025.json)为最新，以下旧“声学 0 s”是历史快照。生产求解器与公共 API 未改。

### S2c-B PHI_PLUS90：预算已授权，RAM 门阻塞（2026-10-24；历史快照）

- 用户将共享核心上限至少增加 306.601 s；按最小增量记录新上限 7506.601 s，PHI 子额度 1820.8869358 s，S3A 3000 s 储备不动。PHI 保守完整计划 1820.88667776 s 已通过成本门。Windows 可用物理内存最新实测 1,131,163,648 B，低于冻结最大 BIE 层级 30% 余量门 2,026,747,612 B；用户要求暂时保持资源阻塞。新 PHI 声学计算仍为 0，状态仍 `S2C_B_PHI_PLUS90_RESOURCE_REFERENCE_BLOCKED`，下方成本门结论仅为授权前历史快照。详见[执行回执补记](results/comparisons/two_dimensional_channel_closeout/20261023_s2c_b_phi_plus90_01/execution.md)。K_HIGH 已有结论及 S3A 保留不变。

### S2c-B PHI_PLUS90：完整参考成本门阻塞（2026-10-23；历史快照）

- `S2C_B_PHI_PLUS90_RESOURCE_REFERENCE_BLOCKED`；唯一新 run `20261023_s2c_b_phi_plus90_01`。M05a 冻结 C2 几何 PNG/FIG、SHA/预算协议、K_HIGH M04 轻量发布包、既有复场独立复核和 S3A 三 seed 只读清单完成；新 PHI 声学场 0。计划 `1820.8866778 s`，本轮在保留 S3A 3000 s 后最多 `1514.2859358 s`，缺 `306.6007420 s`，停止新昂贵求解。共享余额仍 `4514.2859358 s`。PHI 中心过零机理、M05b–e 及三方法 417 点参考均未得出；详情见[主报告最新首节](reports/two_dimensional_channel_closeout_report.md)。K_HIGH 既有合格状态保留，S3A 未运行。

### S2c-B K_HIGH：五层参考与单案例实验完成（2026-10-23；前次最新状态）

- `S2C_B_K_HIGH_COMPLETE_REFERENCE_QUALIFIED`，仅覆盖冻结417点和五层经验数值参考，非模型物理PASS。用户只对K_HIGH将子额度从1800 s提高至2400 s，原7200 s总额、S3A≥3000 s预留、单作业≤900 s、每次真实物理RAM保留30%均保持。五次新BIE边界解、PE模式0、完整FullFK双项积分全部完成；七作业实计1277.0539592 s，累计核心2685.7140642/7200 s，S3A规划余额4514.2859358 s。
- [唯一主报告最新首节](reports/two_dimensional_channel_closeout_report.md)与[中文M04图册](results/comparisons/two_dimensional_channel_closeout/20261023_s2c_b_khigh_02/mechanism_atlas.md)给出真实图片、411线/中心11/九点/五点/边缘104的原始场和经验数值分辨率分类。原生CSV独立复算、四PNG像素重绘及四FIG重开通过；完整审计见`20261023_s2c_b_khigh_02/audit/final_independent_audit_v4.json`。FullFK五点十字的两案变化未过5倍门，保留未分辨。PHI_PLUS90、S3–S6、Bellhop、通信未执行；`G1_target_PM=OPEN`、`WP2_formal=false`。

### S2c-B K_HIGH：单层BIE代表已审计，P1成本门停止（2026-10-23；扩额前历史快照）

- 用户只授权K_HIGH一个新物理海面。独立Goal已创建；[主报告S2c-B附录](reports/two_dimensional_channel_closeout_report.md)与`results/comparisons/two_dimensional_channel_closeout/20261023_s2c_b_khigh_02/`保存冻结配置、新run授权、BASE_M正式417点原生复场、解析C2几何、M04a预实验图、五层BIE成本/RAM门和P0独立审计。S2c-A仍为`S2C_A_DIAGNOSTIC_PROTOCOL_COMPLETE`，S2b仍为`S2B_COMPLETE_REFERENCE_QUALIFIED`。
- S2b正式输出247/247与代码24/24、S2c-A正式输出203/203与代码5/5及22/22上游输入SHA再次一致；K_HIGH配置/源/411线/九点哈希与事前冻结协议一致。P0独立审计的BASE_M CSV/MAT最大复差5.02e-16，解析几何差5.03e-17，|x|≥50 m严格平面；M04a保存MAT重绘像素一致、FIG重开。
- 五层BIE事前保守总估计1749.8 s（含其他物理及失败预留），起初RAM低于2.027 GB门，自动审批据此拒绝；后续代表作业前新实测2,101,432,320物理字节，通过完整五层30%保留门。仅执行L60/12 ppw/qa96 BIE代表，3840节点、417同点原生复场，实耗311.2095523 s；MAT/JSON与SHA/同坐标独立审计通过，另从保存MAT独立重建全部3840节点高度/法向及984项源谱，最大差0。用代表实耗重新保守估计五层加必要预留2157.6201 s，超1800 s子额度357.6201 s，故P1停止，余四层及PE/FullFK均未运行。K_HIGH单层复场及与BASE_M配对CSV只作**参考未合格**的候选响应，不用于5倍经验数值可分辨判定或模型误差；M04b–d正式图未生成。状态`S2C_B_K_HIGH_RESOURCE_REFERENCE_BLOCKED`，Goal不标完整完成；PHI_PLUS90、S3–S6未执行。

### 二维收尾 S2c-A：保存正弦空间诊断与试验冻结完成（2026-10-23；当前阶段状态）

- `S2C_A_DIAGNOSTIC_PROTOCOL_COMPLETE`。唯一新run为`results/comparisons/two_dimensional_channel_closeout/20261023_s2c_a_1642/`，[中文图册](results/comparisons/two_dimensional_channel_closeout/20261023_s2c_a_1642/mechanism_atlas.md)按M00→M01→三档M02a/b→三档M03提供11组FIG/PNG。只分析本地正式S2b保存复场与解析含C2几何，新增声学计算为0；S2b合格状态及完整FullFK版本、全部旧场/指标不变。
- 4932线场/108九点MAT与CSV一致；3336模型诊断、1668 BIE物理响应和40分支持指标分别保存。强档PE全线/中心±1 m的11点/原九点复L2为0.0629112453/0.000209815156/0.0000805108583。仅P21/P22/P23在z=3线；六离深度点不投影、不插值。中心η=0但η′非零，且BIE粗糙相对平面相位−0.0389860 rad；不能据过零唯一归因。平台η″=−K²η，高度与曲率相关完全混杂。
- 独立CSV复算最大指标差1.89e-14、幅相平方恒等式余量1.62e-18；11图保存MAT重绘像素差0、11 FIG重开、22相对图件链接有效，逐张实图查看及三文件Code Analyzer通过。247旧S2b归档SHA、24代码SHA及22上游输入SHA核验通过。保存强档411×45只核验z=3截线，其他深度未升级参考资格；无FullFK伪连续场。
- 冻结BASE_M与K_HIGH/PHI_PLUS90及条件对照的解析几何、原离散源、接收支持、配置SHA和待授权作业。固定A改变K会同时改变坡度/曲率/周期；固定C2改变φ也改变中心多个描述量，不等于只变局部高度。建议S2c-B≤1800 s共享子额度，先K_HIGH代表，两个候选完整参考可能超额；规划预留S3A至少3000 s但不保证可行。历史BIE169–290 s仅量级，下一物理作业仍须新采Windows真实物理可用RAM、30%保留、单作业900 s和数值门。
- 共享7200 s仍扣记1408.660105 s、余5791.339895 s。S2c-A结束，不自动开S2c-B、S3–S6/Bellhop/宽带通信；`G1_target_PM=OPEN`、`WP2_formal=false`。仅工作流/图表/解释与协议变更，无生产API或求解器改动，无范围外通信回归或提交推送。

### 二维收尾 S2b：参考闭合完成（2026-10-23；当前阶段状态，覆盖下方早期受限快照）

- `S2B_COMPLETE_REFERENCE_QUALIFIED`，限固定4 kHz离散Gaussian源、平面及三档C2正弦、M99的411点线/原P11–P33及所报告经验数值分辨率；不是模型物理PASS。正式包在同一run_id的`results/comparisons/two_dimensional_channel_closeout/20261023_s2b_1305/formal_reference_v2/`，[唯一主报告](reports/two_dimensional_channel_closeout_report.md)新增参考闭合附录；早期`S2B_DATA_COMPLETE_REFERENCE_LIMITED`及其文件保留为历史。
- P0的S2a原生411点q/w各12行及旧口径原样保留。108/108原九点真场齐全；正式PE/BIE72格保持原数值，FullFK36格另存完整积分数值版本，原有限L60基线不覆盖。强档PE–BIE空间主复L2仍0.0629112453153，旧入射加权0.0439621797406、旧8.0908%/6.2911%的G口径不改写。
- 在新RAM/成本门通过后作7次BIE边界解（2次成本代表，5次参考），未改原求解器SHA；弱64、中96、强原生14ppw96的整密度与对应历史保存密度逐点完全一致。中/强当前64/96同网格接收复场相对变化约3e-8；全417唯一支持的BIE空间/窗口/自求积经验分辨率约1.5e-7（九点）/2.6e-7（线），不冒充严格误差界。旧执行时SHA仍未知，但由当前已知代码及整密度复现完成数值再资格，未倒填历史哈希。
- FullFK有限窗口误差通过解析平面外场恒等式与独立Gauss8双项差积分分离。物理几何在|x|>=50 m完全平，外部平面贡献解析计入，不是局部截断、拟合或FastFK；真实法向/弧长和local_radius_m=Inf保留。400→800面板复场相对变化不超过1.5e-12，与原有限双积分相减加解析尾项最大差8.22e-11。正式FullFK–BIE线弱/中/强复L2约2.17e-6/1.07e-5/4.29e-5；九点弱/中差异在预冻结5倍分辨率规则下未分辨，不能称更高物理精度认证。
- 累计扣记1408.660105/7200 s，其中1258.660105 s实测、150 s保守旧成本/失败扣记；最大单作业289.742 s。所有新大解/评价前均有真实物理RAM门，无页文件替代。18个新MATLAB文件最终Code Analyzer零消息；正式7张FIG/PNG实际查看、FIG重开、MAT重绘像素差0，CSV/MAT/恒等式/24行分辨率独立审计通过。
- S2b Goal按完成性核验结束，不自动进入S3–S6、Bellhop源映射或通信；`G1_target_PM=OPEN`、`WP2_formal=false`不变。生产求解器/公共API及技术指南实现公式未改，无提交推送。历史执行SHA、严格误差界、FullFK连续x-z图和跨场景推广作为明确限制保留，非本阶段伪造补齐项。

### 二维收尾 S2b：点值齐全、参考受限（2026-10-23；早期快照，已由上方补证覆盖）

- [唯一主报告附录S2b](reports/two_dimensional_channel_closeout_report.md)及`results/comparisons/two_dimensional_channel_closeout/20261023_s2b_1305/`记录本次实际执行；状态`S2B_DATA_COMPLETE_REFERENCE_LIMITED`。P0空间主指标与历史入射权重各12行从不变的411点原生复场复算并通过CSV独立审计；强档PE–BIE空间复L2=0.0629112453153，旧入射加权0.0439621797406保留，历史S2a状态不变。
- 核心PE/BIE/FullFK×平面/弱/中/强×P11–P33已108/108，原48值保留、新增60值（PE18/BIE18/FullFK24）。BIE平面6点为同源解析镜像，弱中12点为保存密度评价；无新BIE边界解。新成功批次实计13.9342315 s；三次S2a实计9.2915103 s，加历史120 s及失败尝试30 s保守扣记，总173.2257418 s/7200 s。每批使用Windows真实可用物理RAM门，估计增量96 MiB，单批最大2.4232 s。
- 匹配九点及接收线最差/最小参考探针上复核BIE空间/窗口/弱档自求积层级、PE半步长与保持同一冻结复谱的984→1968网格、FullFK求积及L60→90固定细网距窗口。FullFK窗口相对变化约2.52e-4–2.69e-4，与主差异相当；中强BIE独立自求积和原执行源码SHA仍缺，不能认证物理精度。P1时新边界重解内存门失败，未强启高内存BIE。数据齐全不等于参考资格闭合，Goal未标完成。
- 最终`figures_v2/`含F01–F07 FIG/PNG；F03只用已有强档PE/BIE真二维场，FullFK连续场缺失不插值。`audit_v2/independent_audit.json`从CSV复核数据、权重、恒等式及敏感性，7张PNG保存MAT重绘像素差0，7张FIG可重开；12个新MATLAB文件Code Analyzer零消息。
- 下一单项建议先审FullFK平面有限窗口尾项与解析镜像参考，事前冻结源周期解释和成本；不自动执行S3–S6/Bellhop源映射/通信。`G1_target_PM=OPEN`、`WP2_formal=false`。生产API/求解器和技术指南实现公式未改，无通信回归、提交或推送；所有旧日志条目为历史，不覆盖本索引当前状态。

### 二维收尾 S2a 实际入口、声源资格与平面审查（2026-10-23；当前阶段状态）

- [唯一主报告](reports/two_dimensional_channel_closeout_report.md)与`results/comparisons/two_dimensional_channel_closeout/20261023_s2a_final/`记录已执行的S2a；`S2A_CONVENTION_SOURCE_BASELINE_COMPLETE`只覆盖固定4 kHz、Gaussian源、z=3 m的PE-mode 0/BIE/FullFK输入资格、平面与三档正弦接收线。S0–S1历史状态不变；S2b–S6未执行，`G1_target_PM=OPEN`、`WP2_formal=false`。
- 实际PE适配—屏—传播—比值入口的`eta_z=0,±0.01 m`正入射控制与独立路径公式最大复差`1.75e-13`。PE/BIE平面原生场相对差`4.70e-13`；FullFK/BIE平面411点M99加权复差`1.85546e-4`，不解释为严格误差界。源与表示合同在核心三方法间闭合；Bellhop物理源强/相位/单位映射未闭合，条件扩展不阻塞核心合格子集。
- 保存/复核四条原生`u_flat/u_rough`线及辅助G、δu、分项与交叉项；三图实际检查并从保存绘图数据逐像素重绘，独立MAT/CSV、哈希和中档FullFK直接积分审计通过。核心S2b的108个方法×案例×九点复场格中48个有真实现值、60个缺失（PE18、BIE18、FullFK24）；旧54值下界和历史图册指标不是当前补算表。
- 新声学工作仅有预算内PE评价与FullFK平面/弱/中接收线积分，强档复用经核验的E12；无新BIE求解/保存密度评价、Bellhop作业、完整x-z场、PM、宽带或FastFK。生产接口/物理未改，故未运行不相关的宽带通信回归或改`vertical_comm_guide.md`。先处理S2b参考与成本缺口，未授权自动开算。

### 二维计划改为同源反射复场并纳入完整Kirchhoff（2026-10-23；仅文档）

- [唯一二维主计划](two_dimensional_channel_closeout_plan.md)现以PE-mode 0、BIE、完整Kirchhoff（FullFK）为核心三方法；Bellhop保留为完成源映射资格后加入的条件扩展。覆盖以下旧条目中“FullFK仅诊断”和旧三方补算清单作为当前任务的解释，历史实绩不变。
- 主量为共同二维源/表示下的u_rough，u_flat控制必须并列；G_m及δu辅助。物理源映射、原生场与固定平面复系数拟合分层保存；约0.27%旧拟合残差不是同源证明，本轮未重算。
- FullFK对应`scripts/validation/support/solve_full_kirchhoff_2d.m`，保留双项Green积分、真实几何/入水法向及弧长，核心设置arc_length、local_radius_m=Inf。它仍是Kirchhoff边界近似，不称为精确解，不以FastFK或局部截断替代。
- S0_S1_COMPLETE_WITH_DECLARED_GAPS保留。当前下一步S2a实际入口/声源资格/平面控制，可在预算内增加必要FullFK固定接收线积分；S2b补核心三方法平面与三档正弦九点，S3有限PM，S4有界收尾。旧54个缺值不含新增FullFK，须重盘点；核心总目标为108个点值覆盖格，非新增作业数。
- 此次仅改写计划、历史报告提示及导航；未运行MATLAB/新场/Goal，未改变生产接口或物理实现，因此未运行生产回归，也无需修改vertical_comm_guide.md实现公式。S2a及后续验证尚未完成，G1_target_PM=OPEN、WP2_formal=false不变。

### 二维收尾计划按高度符号审计修订（2026-10-21；仅文档）

- [主计划](two_dimensional_channel_closeout_plan.md)保留S0_S1_COMPLETE_WITH_DECLARED_GAPS，将当前下一步改为S2a实际入口约定控制与修正基线资格核验，再做S2b九点补齐、S3有限PM、S4收尾。旧8.0908%等保持历史复现地位，不再无条件当作纯PE物理误差。
- 修正图册v2已有三条新版接收线、强正弦411×45场及12张图，保存数据/像素重绘审计通过，优先复用。S2a仍须实际反射/传播/比值入口通过eta_z=0及h=±0.01 m正入射平面波控制；已有符号审计的plane_control直接计算解析exp，尚未覆盖该入口。
- PE复用保存入射场，在验证适配层统一向下几何/屏相位/native表示；旧正号粗糙起始场不直接复用。BIE/Bellhop原生证据按原适用范围保留；PM各分支和辅助mode需独立核验约定。
- 若v2通过S2a资格核验，S2b名义剩余为18个平面分母、24个BIE/Bellhop弱中粗糙点值、12个修正PE弱中离深度点值。强正弦修正场不重复算，按其他合格现值继续去重；点值数不等于作业数。
- 本次仅静态核对报告、代码、CSV/JSON并修订计划、旧报告解释提示、历史首轮提示词及README；未执行MATLAB/新场/Goal，未修改生产代码或修正版脚本。物理参数、点阵、G1_target_PM=OPEN及WP2_formal=false不变；无需生产回归或技术指南方法修订。

### 固定三方可视化高度符号重置（2026-10-21；仅验证/报告工作流）

- 新入口 `prepare_pe_bellhop_bie_corrected_data` 读取旧三案例和强海面二维包保存的复数数值，按物理边界 `z=eta`、深度向下、`exp(-iωt)`，用 `exp(-2iketa)` 重新评价 PE 近海面反射场。BIE、Bellhop沿用原生复数；旧图册不覆盖。新数据、图册与精确幅相平方分解见 [修正报告](reports/pe_bellhop_bie_corrected_visualization_report.md)。
- 强海面固定M99下PE–BIE复L2由旧0.080908变为0.062911，幅值相对L2由约0.010016变为0.0010485，相位RMS为0.062936 rad；Bellhop–BIE复L2不变。平面分母和二维z=3 m截线均闭合。仅保存起始场FFT评价，无新声源传播、BIE求解或Bellhop作业。
- 新数据包包含三份源MAT数值副本、哈希、绘图代码快照、接收线/二维复数CSV和12张PNG/PDF。独立重绘12张PNG像素差为0；MAT/CSV、指标平方分解与源哈希复核通过，见 `metrics/corrected_audit.json`。
- 该重置只涉及固定验证图册，不改变生产高度参数、公共API、旧历史文件或S0–S1保存数据的历史数值。剩余6.29%复差不得单因归因；单频图册仍不支持PM、三维或通信信道结论。生产修复须另审全部调用方及回归。

### 固定三方比较的高度符号审计（2026-10-21；覆盖旧PE误差物理解读）

- [审计](reports/pe_bie_height_phase_sign_audit.md)确认同一η同时进入BIE物理边界z=η（深度向下）和PE的+2kη相位屏，旧图册仅共轭BIE/Bellhop比值以匹配主要相位。平移平面解析路径要求exp(−2ikh)；旧单边共轭不能替代一致的高度/场表示转换。
- 同一保存入射场、同一3m FFT核、同一411点M99下，正号重现旧PE最大差6.22e-14；负号诊断的强海面幅值L2由1.00156%降至0.104850%，幅值相关由−0.977756转为+0.979251。短距离幅值输运解释其近似镜像。负号原生复差6.291125%，仍有模型残差，不能称为完成生产修复或新版验收。
- 旧8.0908%公式可重现，但其纯物理模型误差解读被本更正限制；Bellhop/BIE彼此的共同共轭不改变二者复差。S0–S1的保存数据重建结果仍为历史算术证据，不等价于已验证PE高度符号。未执行S2–S6、PM或通信；G1/WP2状态不变。
- 新增审计入口及独立MAT/CSV，六个保存起始场FFT诊断评价；无新上行传播、BIE求解/场评价或Bellhop作业。源报告/数据与冻结求解器保留；展示汇总、数据说明与导航增加更正链接。生产接口与方法未变，因此不运行生产回归或修改技术指南公式。

### 二维收尾 S0–S1 保存数据审查（2026-10-20；最新阶段状态）

- 已完成 [S0–S1 审查报告](reports/two_dimensional_channel_closeout_report.md)与 `results/comparisons/two_dimensional_channel_closeout/20261020_2330_s0s1/`；状态 `S0_S1_COMPLETE_WITH_DECLARED_GAPS`。S2–S6 未执行，新增物理场计算零；未评价保存 BIE 密度的新接收点。
- 三档正弦三方旧接收线 G 指标 9/9 从完整保存复数和原 411 点 M99 权重重建；强正弦 PE–BIE 整线复 L2 为 0.08090810774，其中相位平方项占 98.4676%。旧公式正确，但整线 G 不能代表中心或完整 3×3 阵列、原始压力或通信性能。
- 强正弦中心 P22 的 G 逐点相对复差 0.0779531；同深度三点 `E_subset=0.0780765`，完整九点 G 因平面分母缺失保持未定义。强正弦九点原生 PE/BIE 幅度可比，逐点相对差均值/中位数/最大为 0.000215237/0.000254523/0.000375290；原生复相位约定未统一，未作九点复误差排名。
- 来源/缺失逐点表、图组和保存 MAT、公式控制 9/9、独立 CSV 审计与七张 PNG 像素一致重绘已保存。S2 最小值下界是共享平面复分母 18 值及弱中正弦粗糙复场 36 值，另需匹配点的数值参考不确定性；不把点值数当求解作业数。`G1_target_PM=OPEN`、`WP2_formal=false`。仅增只读工作流与报告/导航，生产方法及接口未改。

### 二维收尾计划审查（2026-10-20；文档修订，阶段未执行）

- [审查记录与首轮Goal提示词](reports/two_dimensional_channel_closeout_plan_review.md)建议首轮只执行S0–S1；已在原计划0.3节补充执行契约，保存密度的新接收场评价也属于禁止的新增物理场计算。
- 明确完整9点原始场不等于完整9点G；缺失平面分母留待S2。补充零候选相位、部分集E_subset、逐指标覆盖、未知参考不确定性及接收线同深度规则。
- Goal完成仅指已有证据复核、图表、独立重绘和缺口报告完成；方法不通过可以是结果，环境未能运行必需检查则不能写完成。本次未创建Goal或执行S0–S6。
- 仅修改计划、审查和导航文档；物理参数、生产代码、公共接口、历史数值及G1/WP2状态不变，无需生产回归或技术指南修订。

### 二维点阵对比收尾计划（2026-10-20；仅规划，未执行）

- 按用户要求在根目录保存[二维专项执行计划](two_dimensional_channel_closeout_plan.md)：固定中心(0,3)m的3×3 x-z接收点阵，详细展示中心与四个边中点，保存9点均值/中位数/最大值、E_array、幅相分解及缺失/近零标记。
- 核心S0–S4依次覆盖来源与预算冻结、保存数据指标重建、正弦点阵最小补算、3个至总计8个PM输入的匹配对照、二维收尾报告；S5宽带/CIR和S6最小通信为另行指定的条件扩展。所有阶段尚未执行。
- 保持Gaussian角谱源与Bellhop X/internal-wall约定，原始复场与G分开；规定必须绘图、保留逐点数据并独立重绘。历史R源50海面不能并入当前X源样本。
- 本文件是深海总计划下的二维专项，不改变G1_target_PM=OPEN、WP2_formal=false，不启动已停止的三维参考修补。仅新增计划与README导航；无代码、接口、物理模型或数值结果变化，因此无需修改技术指南或运行生产回归。

### BLACK 汇报版报告整理（2026-10-19；仅表述与解释）

- [BLACK 报告](reports/black_channel_basics_report.md)已按海试场景、数据含义、三张图的读法与主要特征重写；补充非专业读者可直接使用的口头汇报段落，逐阵元数值与复现参数移至附录。
- 数值仍来自 `results/validation/black_channel_basics/runs/20261016_175145_949/`，保留三张原图、MAT、CSV及运行状态；没有新下载、回放或物理求解。修订前报告与生成器副本位于 `results/validation/black_channel_basics/report_revision_20261019/`。
- 明确区分归一化幅度/绝对传播损失、窗口相对时延/绝对传播时间、人工20 dB测试SNR/海试实测SNR，以及公共频移/残余谱。-25 dB保留能量比例仅针对保存窗口。
- 补充理想垂直、平面边界、c=1500 m/s的几何推断：Rx6 m对应海面反射超额时延约8 ms，超出主峰后约6.2 ms的现有覆盖；Tx离底2 m对应底反射约2.67 ms，仅作候选，不识别实测反射路径。现有BLACK数据不能认证直达＋一次海面反射二径结构。
- 报告生成器与中文模板同步更新，新运行的源快照包括模板；公共接口、数值分析、生产模型未改。只做报告生成、代码静态检查、数值/链接/输入哈希一致性核验，不重跑通信仿真。

### BLACK 实测信道基础分析（2026-10-16；独立数据工作流）

- 仅下载 BLACK/black_noise，Zenodo 21287414 v1.0，官方 MD5 和本地 SHA-256 通过；固定官方 MATLAB 提交 af1d3814af9199f1e420a248ab27edaced5339cf，未安装全局 add-on。
- MATLAB MCP 既有 R2025b Update 6 + Signal Toolbox 实跑：h_hat 为 176×8×4569，fs_delay=25 kHz、fs_time=125 Hz、fc=18 kHz，实际记录 36.552 s。官方解包、8 阵元短段 BPSK 回放和配套噪声通过。
- 固定 seed1994、8192符号、6.25 ksym/s、100 kHz、首尾20 ms、四分之一记录处起始；重复回放相对误差0，整体SNR20 dB，最大主到达误差0.103957个时延样点（约4.16 us）<2样点。已知双径/负频移/有符号Welch控制通过。
- -25 dB门限下 RMS 时延扩展0.0737–0.0826 ms、90%区间宽0.160–0.280 ms；0.1 s相位拟合的公共频移范围−5.7801至+5.4651 Hz。保存主抽头的残余谱集中在零频，不代表完整多普勒扩展；不从归一化数据估计绝对传播损失。
- 三张PNG/PDF、一张CSV、保存MAT与中文报告完成；从保存数据重建三张PNG像素完全相同，CSV校验和相同。工作流Code Analyzer 0消息，MATLAB当前目录保持不变。未运行其他信道、BER或生产物理仿真，未改生产接口。
- 入口 `run_black_channel_basics`，结果 `results/validation/black_channel_basics/runs/20261016_175145_949/`；重绘验证在其 `render_20261016_175421_919/rebuild_verification.json`；[中文报告](reports/black_channel_basics_report.md)。这项数据分析不更改下方既有物理模型的验证状态。

### 固定三维dual宽带决策（2026-10-15；最新，覆盖下方旧下一步建议）

- 最终B：只交付已验证4/6/8kHz离散物理结果。反射幅度约降6.43%，总场相对平面−4.99/+2.38/+12.82dB；物理增量至少约为保存参考经验敏感性的1.1万倍。参考敏感性、FK求积误差、FK模型差分别记录，不作连续频带或PM推论。
- 已知几何高度相位分离后的稀疏候选仅探索性支持：4/8预测6kHz复误差2.1544e−4；6kHz参与候选筛选，不能充当独立认证。5/7kHz已预声明但因内存门拒绝而未计算。当前可用639.199MiB，保留30%后447.439MiB，小于估计增量1520MiB；新物理点0、核0s。不是新频点物理失败或计算机永久不可行。
- 历史完整频点约44.6–57.9s。257/513源点加逐区间中点验证的新增核外推6.31–8.21/12.65–16.44小时；9/17源点的对应情景10.4–13.5/22.3–29.0分钟，均非采样充分性证明。时间未触发否决，未知点数/当前资源/独立重建证据仍限制交付。
- 停止FE/TFE/BIE连续修补，不再自动执行旧“高阶BIE优先”；冻结控制和失败数据保留。G1_target_PM=OPEN、WP2_formal=false。本Goal结束，不自动新建研究Goal。
  [物理决策](reports/deepsea_dual_broadband_decision_2026-10-15.md)、[成本](reports/deepsea_dual_broadband_cost_2026-10-15.md)、[协议](reports/deepsea_dual_broadband_protocol_2026-10-15.md)；数据results/validation/deepsea_dual_broadband_decision_20261015/。无生产/求解器/接口更改。

### BIE最小验证与已知路径重建（2026-10-14；当前决策）

- 当前多项式预条件器保持FAILED_CLOSED，不重复调参。BIE六个制造逆问题线性约1e-10通过，但独立边界/场全失败；6kHz最细1600未知半单元边界3.385%、离面.422%，实际C2图面制造PV3.616e-4仍不过1e-4门。停止扩大，不运行dual散射或原PM BIE。
- 辅助半空间Green、法向跳跃、Duffy自项/近奇异减项及几何FD有局部支持；未知粗糙密度外区闭合未验证。下一单项先测试同6kHz平面制造问题的高阶密度/近场规则，不用大PM赌算力。
- 同257/513频点，已知几何/分量时延分离的独立频点复误差≤6.22e-12、深衰落≤1.73e-9；通用total-only插值与dual中间频率仍未认证，无PM通信结论。
- A含局部检查与诊断约166.528s，B .991s；采样空闲RAM最低1838.63MiB，生命周期working-set峰1313.44MiB只作任务峰保守上界。G1 OPEN/WP2=false。
  [BIE报告](reports/deepsea3d_bie_minimum_report.md)、[重建报告](reports/deepsea_frequency_reconstruction_report.md)、[协议](reports/deepsea_bie_reconstruction_protocol_2026-10-14.md)；数据results/validation/deepsea_bie_reconstruction_20261014/。

### 双轨检查点与正式计划合并（2026-10-14；当前决策）

- 唯一正式依据：[deep_sea_channel_research_plan.md](reports/deep_sea_channel_research_plan.md)，已核对合并v3。原9公式、WP0–WP6、验收表及完整历史块保留；backup/v3仅版本记录。
- A唯一干预为右P=2I−B，冻结原物理B不改。平面/dual/quarter通过；原PM96内步后B8.3065e-3、真PDE3.3601e-3/DtN3.2457e-3、真边界1.040333e-4仍失败且比旧restart12更差。停止本轮调参/扩网格，转匹配三维BIE窗口/奇异积分可行性评估，不实现大型求解器。
- A新物理核410.957s，原PM241 B调用/288.323s，未超3600s/640调用/600s；内存只具事前70%物理空闲门和估计/快照，任务峰值未知。
- B解析native/absolute-DSP/relative-DSP、4–8kHz端点/DFT、flat与离散dual、AWGN/QPSK已知信道基础控制通过，wall2.049s。513点线性插值仍约1.3%误差，未认证通用插值/完整RRC-FSE-pilot接收机或dual宽带；无PM通信结论。
- 冻结算子回归4/4通过；新数据审计通过。保持G1 OPEN/WP2=false，无生产/冻结求解器/旧数据更改，不提交推送。
  [A报告](reports/deepsea3d_preconditioner_decision_report.md)、[B报告](reports/deepsea_singlelink_readiness_report.md)、[协议](reports/deepsea_two_track_protocol_2026-10-14.md)；新入口validate_deepsea3d_preconditioner_decision、validate_deepsea3d_preconditioner_checkpoint、validate_deepsea_singlelink_readiness；数据results/validation/deepsea_two_track_20261014/。

### Deep-sea matrix-free t=1 mapped Helmholtz pilot（2026-10-13；最新）

- 新验证专用直接法：精确图面点源Dirichlet、原去混叠L0+L1+L2/DtN、L0模态预条件GMRES，
  无Taylor求和/源项展开。连续map Jacobian下界.3221已用解析梯度和采样间隙界认证。
- 平面/dual/quarter控制通过；dual对冻结参考增量差1.0935e-5。四项制造/预算测试通过。
  原PM restart8/12各96内迭代仍未线性收敛：B残差1.822e-3/1.359e-3，离散体1.155e-3/9.212e-4。
  restart12真PDE1.329e-3、DtN1.273e-3、边界1.04033e-4未过门；接收差虽3.65e-5不能认证。
- 对旧Taylor N32体/DtN残差改善357.6/256.4倍，但不是物理解误差或速度认证。
  停止扩算；优先有界改进预条件以关闭离散线性门，之后才评估分辨率/正式收敛；失败再匹配3D BIE。
- PM INCONCLUSIVE、G1 OPEN/WP2=false，无FullFK–PM/support24/MC/FastFK。新核wall464.135s，
  工作数组估计663.19MiB、快照非峰值；冻结312文件不变，生产/历史结果未改。
  入口`validate_deepsea3d_mapped_direct`，数据`results/validation/deepsea3d_mapped_direct/`；
  [中文报告](reports/deepsea3d_mapped_direct_report.md)、[协议](reports/deepsea3d_mapped_direct_protocol.md)。

### Deep-sea TFE dealiasing failure diagnosis（2026-10-13）

- L2为三因子乘积；新验证分支2n填充且最后一次投影，DtN乘积一并处理，5项制造测试通过。
  对称Fourier空间显式去Nyquist；同空间对照纠正并保留初始无效数据及成本，不改旧求解器。
- dual与quarter PM局部控制通过；原PM共享同系数/support12/6kHz/q.10，N16/24/32原网格。
  去混叠N32末阶体场系数降2.363倍，但PDE .410575→.475266、DtN .287852→.326480，未获参考收益。
  独立3n/4n算子评价一致，99.60%物理残差能量仍在保留空间；真几何算子差1.1521e-5。
- 当前t=1有限Taylor和仍明显失效，不能只归因混叠或真残差评价网格；不证明连续TFE必然发散。
  停止扩阶/扩网格；建议下一项矩阵自由直接求解t=1 mapped Helmholtz，失败再考虑匹配3D BIE。
- PM INCONCLUSIVE、G1 OPEN/WP2=false。无support24/FullFK–PM/Monte Carlo。
  新核wall486.343s含初始错误及probe，命名数组估计384.75MiB、非峰值；冻结284文件不变。
  入口`validate_deepsea3d_tfe_dealiasing`，数据`results/validation/deepsea3d_tfe_dealiasing/`，
  [中文报告](reports/deepsea3d_tfe_dealiasing_report.md)、[协议](reports/deepsea3d_tfe_dealiasing_protocol.md)。

### Deep-sea 3D TFE independent reference prototype（2026-10-13）

- 新验证专用TFE：z=s+t(1-s/b)eta，独立mapped Helmholtz与精确Fourier DtN；
  Chebyshev M32/40、传播/倏逝/掠射模态及近内部极点保护，不调用FullFK/普通FE求解。
  冻结PM cfg/系数只读加载，没有重生成海面、重跑FE失败或修改生产。
- 平面与实现测试通过；认证dual的阶数/横向/纵向/域/透明平面/真实边界与PDE门均通过，
  最细TFE与冻结独立参考增量差6.308e-7，成功复现该离散案例。
- 原振幅PM support12仍INCONCLUSIVE：N32→40变化1.923，细h.125体PDE≈.120/DtN≈.173，
  纵向M40不能修复；N40体场尾部增长。quarter幅度仅局部残差通过，不是原目标或正式参考PASS。
  不做support24/FullFK–PM；全G1 OPEN/WP2=false，建议有界dealiasing/TFE改进或匹配3D BIE原型。
- 16粗糙解+1平面，保存核wall278.451s；体场74.25MiB、工作数组估计891MiB、快照非峰值。
  新测试3/3和数据审计通过，冻结263文件不变；入口`validate_deepsea3d_tfe_prototype`，
  数据`results/validation/deepsea3d_tfe_prototype/`；
  [中文报告](reports/deepsea3d_tfe_prototype_report.md)、[协议](reports/deepsea3d_tfe_prototype_protocol.md)。

### Deep-sea 3D PM fixed-realization feasibility（2026-10-13）

- 新目标：各向同性有限带PM、Hm0=.5m、seed261013、6kHz、Tx(0,150.5,1500)/Rx5、q=.10。
  父域128m、kmin2pi/128、kmax1rad/m，固定alpha解离散谱质量，Tp3.638671s；
  实际Hm0=.504033，无抽样后RMS缩放。完整母系数保存，长短截止变体共享系数。
- 9个有界FE探针：N20仍不过阶数/真实边界/源项门；粗/细真实边界残差62.26/919.56。
  小接收网格/域变化不能替代参考认证；精确有限模态指数15点见证残差281.9/3637.9。
  加宽C2改变物理面并恶化参考；支撑与谱截止原始压力差均未获物理解释资格。
- 无FullFK运行或新条件性G1；单案例INCONCLUSIVE，全G1 OPEN、WP2=false。
  下一步需匹配TFE/3D BIE参考及非紧支撑处理，不能直接扩大FE规模或因dual通过放行。
- 新入口`validate_deepsea3d_pm_feasibility`、数据`results/validation/deepsea3d_pm_feasibility/`；
  [中文报告](reports/deepsea3d_pm_feasibility_report.md)、[协议](reports/deepsea3d_pm_feasibility_protocol.md)。
  新测试3/3、数据审计通过，冻结246文件不变；保存核约91.223s、谱载荷180MiB非峰值。
  不改生产/冻结求解器/旧数据，不做Monte Carlo、FastFK、宽带统计或通信仿真。

### Deep-sea G1 fixed-y 4/6/8kHz conditional checks（2026-10-13）

- 固定Tx=(0,150.5,1500)、Rx=(0,0,5)、q=.10/alpha2.5/dual/C2，唯一改动频率4/8kHz。
  6kHz只读复用。4kHz N16足够，8kHz按需N20，N16→20变化8.931e-5接近原1e-4门。
  两频独立网格/域/真实边界/消项通过，L32/h=.0625，无需L48；M32/M40逐频评价及尾界通过。
- 8kHz初始Gauss4 h.05增量变化1.200e-3失败，先标FK数值INCONCLUSIVE，未判物理FAIL；
  有界补一次Gauss4 h.025后变化1.683e-9，失败原始解及成本保留。
  **4/6/8kHz各离散案例条件性G1 PASS**：反射复误差4.034e-4/2.843e-4/2.146e-4，
  增量2.235e-4/1.375e-4/1.217e-4；各过原1%门，无拟合。
- 含截断阶数敏感性的范围支持4→6复误差下降，6→8复误差范围重叠、趋势未分辨。
  两项相位绝对误差三点下降可分辨；增量幅度误差6→8上升，不能概括越高频整体越准。
  8kHz参考total敏感性8.122e-5不忽略；非宽带或连续频率规律认证。
- 新FE含M32共11解，FK11解（含初始失败及精化），核wall102.477s，尾界诊断另.345s，
  最大谱载荷320MiB，快照非峰值；201冻结文件未变，分类4/4通过。
  全G1 OPEN，WP2=false；不改生产/二维/旧数据，不做PM MC/FastFK/SSP/动态/WP2。
  入口`validate_deepsea3d_g1_frequency`；数据`results/validation/deepsea3d_g1_frequency/`；
  [中文报告](reports/deepsea3d_g1_frequency_report.md)、[协议](reports/deepsea3d_g1_frequency_protocol.md)。

### Deep-sea G1 q=.10 x/y azimuth check（2026-10-13；最新）

- 只将Tx从(150.5,0,1500)变为(0,150.5,1500)，Rx仍(0,0,5)，海面中心/波矢/phase/C2不变。
  两方位平面直接/反射/total背景完全一致；新y方位按需N20、h=.0625、L32通过独立门，
  网格/域变化8.564e-7/3.037e-7，半网格真实边界残差3.211e-7，无需L48，消项回放匹配。
- 原FullFK网格/Gauss/孔径通过；**y方位单案例CONDITIONAL_G1_PASS**，
  反射/增量复误差2.843e-4/1.375e-4，相比x误差下降1.920%/1.762%。
  含N16/N18压力敏感性的经验范围仍不重叠；不外推通用方位规律。
- 独立参考y−x反射/增量复变化5.111%/2.467%，差幅2.504e-6m^-1；
  含阶数的参考配对最大位移2.051e-10m^-1。区分方向性场差、模型误差与已观察参考数值变化。
  新保存核wall48.866s、最大谱载荷320MiB，快照非峰值；169冻结文件未变，分类4/4通过。
- 全G1目标域OPEN，WP2=false；不修改生产/二维/旧数据，不做PM/FastFK/SSP/动态或其他频率。
  入口`validate_deepsea3d_g1_azimuth`，结果`results/validation/deepsea3d_g1_azimuth/`；
  [中文报告](reports/deepsea3d_g1_azimuth_report.md)、[协议](reports/deepsea3d_g1_azimuth_protocol.md)。
  下条x方位保持冻结；本轮只认证两个离散方位对照。

### Deep-sea 3D G1 q=.10 geometry-only check（2026-10-13；最新）

- 唯一物理改动Tx=(150.5,0,1500)，Rx=(0,0,5)，q=.10；沿用alpha2.5/6kHz/dual/C2，
  海面中心保持(0,0)，不移向(.5,0)镜面点。只读复用WP0平面解析控制，相位闭合8.53e-13。
- 低成本N12/N16不足，按原门追加N20；新几何独立网格/域变化1.400e-6/4.161e-7，
  最终半网格真实边界残差3.013e-7，尾部递减，内部消项回放逐阶匹配；无需L48。
  随后FullFK网格/Gauss变化5.991e-9/6.923e-8，A12/16/24孔径门通过。
- **q=.10离散案例 CONDITIONAL_G1_PASS**：原生反射/增量复误差2.899e-4/1.399e-4。
  对冻结q0分别下降2.03%/2.94%，已检验数值变体范围不重叠；仅离散经验可辨识。
  粗糙total相消指标1.278，total误差略增至1.781e-4，不能由反射误差推断total误差。
- 新增核wall53.887秒，最大谱载荷320MiB、内存快照非峰值；分类4/4通过，Code Analyzer零问题。
  冻结140文件未变；不改生产/二维/旧数据，不做PM/FastFK/SSP/动态/更大振幅/WP2。
  全G1目标域仍OPEN，WP2=false。入口`validate_deepsea3d_g1_q10`；
  数据`results/validation/deepsea3d_g1_q10/`；[中文报告](reports/deepsea3d_g1_q10_report.md)、
  [预设协议](reports/deepsea3d_g1_q10_protocol.md)。下条q0仍是冻结基准。

### Deep-sea alpha2.5 domain and conditional physical G1（2026-10-13；最新，接续下条）

- 同一冻结Tx1500/Rx5、q0、6 kHz、dual振幅(.05,.025) m及C2 plateau8/support12。
  复用L24/h=.0625/N20，新增L32：域增量变化3.242e-7、N16→20变化2.363e-6、
  半网格真实边界残差2.977e-7，均过原1e-4门；尾部递减，内部消项指标3314。
  按预设规则无需L48，未单独重复N16求解；不声称无限域严格误差界。
- 参考通过后五次原FullFK：网格/Gauss变化1.021e-8/1.905e-7；
  独立A12→16→24孔径变化0/2.561e-15。原生反射/增量复误差2.959e-4/1.442e-4。
  **alpha2.5离散案例 CONDITIONAL_G1_PASS**；全G1目标域仍OPEN，WP2=false。
- 新增保存核wall31.053秒，谱载荷320 MiB，内存快照非峰值；分类回归4/4通过，
  Code Analyzer零问题，CSV/复数/相位护栏与冻结哈希通过。冻结求解器、生产、二维、旧数据未改。
  入口`validate_deepsea3d_alpha25_domain`；数据`results/validation/deepsea3d_alpha25_domain/`；
  [中文报告](reports/deepsea3d_alpha25_domain_report.md)、[预设协议](reports/deepsea3d_alpha25_domain_protocol.md)。
  本条更新下条的物理INCONCLUSIVE，不放行更大振幅/PM/FastFK/WP2。

### Deep-sea alpha2.5 fine-reference check（2026-10-13；最新）

- 同一冻结alpha2.5 cfg，L24/h=.0625、N16/N20、M40完成参考验证，不重算粗网格。
  细阶数差2.36e-6、粗→细增量变化1.77e-6、N20半网格残差3.06e-7，各过原1e-4门。
  h归一化系数尾部与粗网格一致，尾部均衰减；内部消项指标N16/N20为655/3314。
- **细网格参考检查PASS**，未见当前配置普通FE失稳；没有域扩展或FK物理对照。
  alpha2.5物理G1仍INCONCLUSIVE，全目标域OPEN/WP2=false。仅建议受控扩域后物理比较。
- 保存计算26.138秒（FE+回放，未计IO/分析），谱载荷最大180 MiB；内存快照非峰值。
  入口`validate_deepsea3d_alpha25_fine`；数据`results/validation/deepsea3d_alpha25_fine/`；
  [中文报告](reports/deepsea3d_alpha25_fine_report.md)。冻结求解器、生产、二维与旧数据不变。

### Deep-sea alpha2.5 reference preflight（2026-10-13；仅预检）

- 固定原单链与dual谱/phase/C2，只将振幅统一放大2.5；L24/h.125，N12/16必要补N20。
  N12边界残差4.07e-4失败；16→20增量变化2.36e-6，N20半网格残差2.87e-6通过预检门。
  系数尾部衰减，内部消项指标最高3314、接收求和指标15.62；诊断回放逐阶匹配冻结FE。
- 值得继续低成本细网格正式参考检查，普通FE粗配置未见灾难性失稳；未跑FK或域扩展。
  alpha2.5物理G1仍INCONCLUSIVE，全目标域OPEN，WP2=false；不默认开启完整验证。
- 保存计算wall5.944秒，进程内存前后快照约3.12 GiB（非峰值），谱载荷最大45 MiB。
  入口`validate_deepsea3d_alpha25_preflight`；数据`results/validation/deepsea3d_alpha25_preflight/`；
  [中文报告](reports/deepsea3d_alpha25_preflight_report.md)。冻结求解器及历史数据不变。

### Deep-sea 3D G1 amplitude-only extension（2026-10-13）

- 固定Tx1500/Rx5、q0、6 kHz、dual谱/phase、C2支撑及原生约定，只乘统一alpha。
  1.0冻结基准复用，.5/1.5/2.0新增点均获得条件性G1 PASS，无可信物理失败或FE失稳。
- 每新增点独立N12/16、网格/域和原网格/半网格真实边界检查；N16实际足够，非默认。
  最大最终细阶数/网格变化1.83e-6/1.51e-6，半网格残差2.44e-7；alpha2粗阶数差4.62e-5。
  四点反射误差6.19e-5/1.36e-4/2.07e-4/2.57e-4，增量最大1.27e-4；误差主要在相位。
- 仅四个离散点，不是连续振幅区间或PM总体认证；全G1目标域OPEN，WP2=false。
  如继续，仅建议2.5低成本阶数稳定性预检；尚未找到物理失败边界。
- 入口`validate_deepsea3d_g1_amplitude`；结果`results/validation/deepsea3d_g1_amplitude/`；
  [中文报告](reports/deepsea3d_g1_amplitude_report.md)。新增保存计算wall137.35秒，1.0零新增求解。
  未修改冻结三维、生产或二维模型；分类回归4项通过，不重复原WP0/WP1/G1。

### Deep-sea 3D finite-roughness G1 single-case（2026-10-12；接续WP1）

- 一个Tx1500/Rx5、q0、6 kHz、C2双方向面，振幅(.02,.01) m、support12 m。
  独立16阶Dirichlet FE保留倏逝模态，完成自身阶数、网格、域和真实边界残差检查。
  最大网格增量变化1.07e-6，N16原网格/半网格边界残差9.37e-12/2.02e-7，三个高阶测试通过。
- FullFK原生反射/粗糙增量误差1.36e-4/1.06e-4；高阶与一阶增量差70.8%。
  6033.0075 Hz局部相消总场误差.364%，独立参考自身总场变化3.74e-5。
- **该有限振幅案例条件性G1 PASS**；不是严格余项定理、band-relative或PM总体认证。
  全目标域G1仍OPEN，WP2不批准。此条接续下条的一阶参考缺口，不重新做二维研究。
- 入口`validate_deepsea3d_g1_preflight`、`validate_deepsea3d_g1_finite`；数据
  `results/validation/deepsea3d_g1_finite/`；[中文报告](reports/deepsea3d_g1_finite_report.md)。
  保存计算wall137.45 s（含探针/参考/比较，不含IO与开发）；无密集BIE、新依赖或生产变更。
  原三维FK、源/海面/相位helper及二维/生产基准保持哈希。下一步仅受控物理域扩展。

### Deep-sea WP0 / minimum 3D WP1（2026-10-12；validation only）

- Tx1500/Rx5、q=0/.10、c1500、单位3D Green点源、均匀无损压力释放半空间。
  G0 PASS：24组平面/移位/旋转/表示控制，最大代数/相位误差7.28e-12/5.26e-12。
- 单Rx三维FullFK保留真实eta(x,y)、法向/面积、解析源导数及两项Green积分。
  解析无限平面加紧支撑密度差。flat、weak(.001 m)、dual(.002/.001 m)的
  6 kHz数值门通过，再通过4/8 kHz网格检查；不是PM海况或3D FastFK。
- 一阶Dirichlet参考完成12个振幅/几何对照、Gauss网格/支持护栏及Fourier交叉检查。
  相消控制最大条件数约271；8项实现测试通过。
- **有限粗糙度物理G1仍OPEN**：没有有限振幅余项界或匹配3D BIE。
  不批准WP2/PM总体/统计/SSP/动态/CIR/BER；允许继续三维独立物理参考补强。
- 入口`validate_deepsea3d_wp0_wp1`；数据`results/validation/deepsea3d_wp0_wp1/`；
  [中文报告](reports/deepsea3d_wp0_wp1_report.md)、[协议](reports/deepsea3d_wp0_wp1_protocol.md)。
  计划修订避免重做已完成2D域审查/单Rx算法。生产PE、公共接口与既有2D基准不变。

### 固定 Tx–单 Rx Kirchhoff 原型与有条件冻结（2026-10-12；validation only）

- 独立完成FastIn + FullOut与固定Tx–Rx复合偏导核，保留真实高度、法向、
  弧长权重及完整两项积分；未改变生产PE/公共接口或FullFK/BIE/FastFK基准。
- Tx100/Rx(0,3)、c1500、4/6/8 kHz、4097点、Gaussian入射，12个PM输入
  加5控制；算法复误差最大7.35e-4 / 3.30e-13，全部数值护栏通过。
- 本轮建议FastIn为默认验证原型，复合核为高算法精度备选；精确入射已提供
  时继续使用原始标量FullFK。跨realization成本与冷启动分开测量，不能使用
  接收线FFT的历史加速倍数。结果不建立新海面物理真实性、3-D或BER结论。
- 验证补强要求异常标识符精确匹配，141条依赖/拒绝/验证器控制；FastIn
  不绑定Rx，复合核不绑定无关rank，海面及法向可更新。只修正plan契约和
  错误标识符，计算公式不变；原Rx0压力逐值对补强前保存证据核对。
- **仅原Rx(0,3)研究fixture可冻结**。已有强正弦控制的低压位置诊断发现
  FastIn rank3相对复误差约.00225/.00575/.01736，均不满足原1e-3；复合核
  通过但新Rx不再普遍是舍入级误差。未提高基线rank或放宽门槛；相消点
  不纳入主域PASS。raw比例、幅度、相位、有限性和零点undefined均明确。
- 补强前MAT/CSV在本case`pre_hardening/`；新plan_controls/cancellation CSV、
  `freeze_manifest.json`与执行代码哈希保存完整冻结边界；
  [冻结审核与所有未关闭问题](reports/fk_single_receiver_baseline_freeze_report.md)。
- 入口`validate_fk_single_receiver_methods`，支持`full/render`；数据在
  `results/validation/fk_single_receiver_methods/`；
  [中文评价](reports/fk_single_receiver_methods_report.md)、
  [公式与试验协议](reports/fk_single_receiver_methods_protocol.md)。

### Deep-sea Fig12(b) native Bellhop reproduction attempt（2026-10-10）

- Final stop:the user capped further work at three strategies.All three completed
  (`CB R/CR R/CB X`,15 additional native runs with beam/step/free-field controls).
  No strategy passed credible whole-figure reproduction;parameter search stopped.
  Cerveny is visually closest but strongly asymmetric and fails amplitude control;
  the line-source alternative remains axis-empty and is too bright.
  [Final Chinese stop report](reports/bellhop_fig12b_three_attempts_stop_report.md);
  fixed entry:`validate_bellhop_fig12b_three_attempts('full')`;
  data:`results/validation/bellhop_fig12b_three_attempts/`.
- Independent validation only:5000 m Munk,Tx4950 m,10 kHz,flat pressure-release
  surface and assumed fluid halfspace.34 native field/ray runs;known paper
  information and added assumptions are separated in the
  [Chinese report](reports/bellhop_fig12b_reproduction_report.md).
- Partial qualitative reproduction only;not a validated full field.Native point
  source forces the entire r=0 pressure column to zero;strict vertical ray works
  while the degenerate one-beam coherent field is zero.Full-angle near-axis
  magnitude fails a homogeneous1/R control;beam convergence remains incomplete.
- Bilateral[-90.06,-89.94]deg,8000 beams,STEP.625 m is the most refined tested
  narrow fan without storage warnings.2000-to8000 field change remains.1129;
  smaller.3125 m steps hit100000-step storage limits for386 narrow rays.
  34 records have no nonfinite complex pressure;artifact PASS is not physics PASS.
- Entry:`validate_bellhop_fig12b('full')`;data:`results/validation/bellhop_fig12b/`.
  Production physics,public interfaces,100 m Bellhop installation and communication
  chain remain unchanged.No CIR,wideband communication or BER is performed.

### Deep-sea single-link research plan（2026-10-09；planning only）

- Current recommended roadmap: [complete research plan](reports/deep_sea_channel_research_plan.md).
  User target is Tx deeper than 1000 m, cable-deployed Rx3–10 m, 3D propagation
  over a 2D rough surface, initially static/uniform and preferably4–8 kHz.
- Proposed sequence: analytic 3D controls and single-receiver rough-surface
  candidate; independent physical/cost pilot; paired mechanism study;
  conditional statistical generator and communication acceptance; domain
  extension; continuous surface dynamics and Doppler. Analytic homogeneous
  illumination is the first candidate, not a production WAPE replacement.
- Scenario levels, sample counts, tolerances and effort bands are proposals,
  not completed experiments or newly confirmed engineering requirements.
  No solver, public interface or validation status changed in this planning task.
- Background and evidence: [handoff](reports/deep_sea_channel_research_background_for_new_ai.md)
  and [external-reply reassessment](reports/deep_sea_channel_external_reply_reassessment.md).

### Full/Fast FK multi-frequency consistency（2026-10-08；最新扩频验证）

- **PASS**, validation only: canonical seed260001, actual Kmax=.471238898 rad/m,
  unchanged 12 Fourier modes, C2 plateau42/support50 m, uniform c1500,
  frozen Gaussian source and Tx100/Rx3 research fixture. Sparse4/5/6/7/8 kHz
  gates preceded65 frequencies at delta-f62.5 Hz. This fixture is distinct
  from the approved Tx99 Bellhop installation visualization.
- Frozen rank3/8 fails from6 kHz; common3/10 narrowly fails at8 kHz. One
  **3/12** rank with padding2 passes the complete band: max weighted complex
  error2.30670e-4, whole-grid2.87485e-4. BIE anchor errors at4/6/8 kHz are
  8.00096e-4 / 5.30379e-4 / 4.08269e-4 in the frozen M99 core. All three
  surface-grid, interpolation and padding guards pass.
- Native exp(-i omega t) physical FK fields require a documented full-envelope
  Fourier-sign bridge for MATLAB IFFT: `H=conj(P)*exp(+i*omega*tau_direct0)`.
  Supplying conjugated native reduced envelopes to the unchanged helper's
  direct_dsp target closes this bridge at2.07e-13. Native reduced fields are
  not silently treated as public PE reduced fields; native physical outputs
  are preserved separately. The unadapted carrier-only import is diagnostic.
- Max Full/Fast relative phase1.74394e-4 rad, delay difference70.19 ns,
  reflected CIR complex error7.17177e-5. Native group delay is about3.595 ms;
  sampled CIR peak3.69231 ms. The16 ms periodic window is diagnostic and does
  not establish unaliased long-path tails or a final converged CIR.
- Minimum measured cold/warm speedups are9.83 / 41.94 across65 samples.
  Memory figures are named arrays and block estimates, not process peaks.
  Five independent phase/surface tests and saved-artifact assertions pass.
  FullFK/BIE/FastFK references and all production/public interfaces remain
  unchanged. Next scope may be fixed-PM broadband plus a small realization
  set with explicit phase metadata; no database/Monte Carlo/BER/SSP/3D or
  production integration is performed here.
- Entry: `validate_spectral_fast_kirchhoff_multifrequency_consistency`;
  report: `reports/spectral_fast_kirchhoff_multifrequency_consistency_report.md`;
  artifacts: `results/validation/spectral_fast_kirchhoff_multifrequency_consistency/`.

### Full Kirchhoff PM validity boundary（2026-10-07；最新物理适用域审查）

- Validation-only continuation of one canonical realization: 15 controlled cases,
  eight master-prefix cutoffs (modes7/10/12/14/16/17/18/19), identical 4 kHz,
  c1500, Tx100/Rx3, discrete Gaussian source and C2 embedding. Common Fourier
  coefficients are exact; no new random realization, phase fit or global scaling.
- In the sampled master-prefix family, mode18/Kmax=.706858 rad/m has native
  weighted FullFK--BIE error8.5275e-4; mode19/.746128 has1.19749e-3. Four
  selected cases have Full4097/8193, BIE ppw8/order8 to ppw10/order10, and
  L60/L70 guards. Summed numerical changes are about6.8e-6, below1e-4.
  The scope is the frozen M99 receiver core, not whole-line window edges or
  infinite-ocean flux. Detailed converged values/classifications are in the report.
- Equal-height and equal-illuminated-slope controls favor illuminated curvature/k
  jointly with high-K amplitude/modal structure as useful descriptors. They do not
  identify a unique causal scalar. At the same positive-amplitude Kmax=.746128,
  factors .25/.50/.75/1 give errors .000796/.000874/.001015/.001197; the guarded
  .75 case also fails. Kmax alone and height alone cannot define a universal gate.
- The error is mainly varying native phase, with magnitude error near1.3e-4.
  Canonical/high phase contributions account for97.1/98.7 percent of complex
  error energy; sampled core spectral discrepancy is concentrated inside incoming99.
  Raw field/spectral means and standard deviations are saved; no alignment is used.
- **Retain FullFK/FastFK in the measured conditional domain.** If the intended
  spectrum must include rejected high-K/amplitude structures, the next physics task
  is a limited SSA/SSA2 comparison. No population failure rate is inferred from one
  realization. Fast ranks3/8/padding2 and all production/reference physics are frozen;
  no frequency extension, Monte Carlo, SSA implementation or production integration.
- Report: `reports/full_kirchhoff_pm_validity_boundary_report.md`; artifacts:
  `results/validation/full_kirchhoff_pm_validity_boundary/`; entrypoint:
  `validate_full_kirchhoff_pm_validity_boundary`. Three metric/source tests, artifact
  assertions, Code Analyzer and nine frozen reference/input hashes are checked.
  This completes the physical-domain investigation proposed by the prior PM pilot;
  that pilot's overall high-cutoff failure remains valid.

### Fixed PM Spectral Fast FK pilot（2026-10-07；最新验证状态）

- Complete validation-only pilot: three saved coefficient realizations
  (canonical seed260001 and Gaussian seeds260002/260003), 4 kHz, U=6 m/s,
  Tx100/Rx3 m. Canonical physical Kmax=.471239 rad/m; nested cutoff controls
  .274889/.746128 preserve all common Fourier coefficients without scaling.
  Numerical Nyquist107.233 rad/m is independently fixed by dx=.029296875 m.
- Frozen rank3/6 fails with field errors about .0035--.0047. Bounded common
  rank diagnosis identifies outgoing truncation as the main baseline issue.
  Common rank3/8 passes all five cases: weighted errors .000187--.000248,
  whole-grid below .0003 for C2-embedded PM and below .000311 for unmodified
  Fourier PM segments. Common4/8 diagnostic reaches1e-4 on embedded cases.
  Incoming/outgoing matrix ranks at1e-3 are3/7 and do not grow on refinement.
- Canonical Full4097/8193 weighted convergence6.56e-6, padding2/4 about4e-16,
  height-node33/65 difference2.51e-6. High-cutoff outgoing theta99=25.54deg;
  8.685 percent of reflected propagating aperture-spectrum energy lies outside
  incoming99. No incident-angle outgoing cutoff is applied.
- **Algorithm compression PASS; overall PM pilot NOT PASS at the declared
  physical-reference gate1e-3.** Native FullFK--BIE weighted error is8.001e-4
  for canonical and1.199e-3 for high cutoff after Full quadrature refinement.
  High-cutoff BIE spatial/window changes are6.62e-8/1.43e-7, so the latter
  deviation is not removed by numerical refinement. Physical comparisons
  concern the fixed M99 core of the common C2-embedded surface; whole-grid
  BIE edge errors remain window-sensitive and are not physically certified.
- N4097 embedded Fast warm speedups89--112x, cold11.7--13.9x. Memory values
  are named-array payload estimates, not peak RSS. The PM full model's
  physical validity is the next issue; increasing Fast rank cannot correct it.
  **Do not enter multi-frequency/Monte Carlo yet.** Establish a PM validity
  boundary or compare an independent SSA/SSA2 reference first; neither is run.
- Reports/data: `reports/spectral_fast_kirchhoff_pm_pilot_report.md` and
  `results/validation/spectral_fast_kirchhoff_pm_pilot/`. Entrypoint:
  `validate_spectral_fast_kirchhoff_pm_pilot`, checkpoint stages fast/bie/raw/
  quadrature/render. Five numerical unit tests and artifact assertions pass.
  No PE/WAPE, Full FK, BIE, local screen, production API or communications
  physics is changed. Previous prototype-only PM recommendations below are
  historical; this dated pilot judgment supersedes them.

### Spectral Fast Kirchhoff FFT prototype（2026-10-07）

- Validation-only Prototype-A (flat normal/dx) and Prototype-B (true normal,
  arc-length weights and both Green derivative directions) are implemented
  under `scripts/validation/support/`, independent of the frozen Full FK,
  PE/WAPE, BIE and local-screen reference implementations. No production API
  or configuration has changed.
- The complete outgoing Hankel/derivative lag kernels are globally compressed
  in height, then applied using FFT linear convolution. No dense receiver×surface
  kernel remains in the Fast implementation. Strong `exp(-ik*eta)` phases are
  retained exactly on both sides. Outgoing rank is separate from incident rank:
  global `r_in=3,r_out=6` is the selected prototype, not an alleged total rank 3.
- At N_surface=N_receiver=4097, padding=2, sigma=.3 m and Rx=3 m, Prototype-B
  FastFK/FullFK weighted native complex L2 is 2.899e-9 (weak .01/.1),
  3.208e-5 (strong .2/.1), 2.115e-8 (high-K .02/.47). K=.70, sigma=1,
  and Rx=2 controls also pass 1e-3. Padding 2/4 differs near 3e-16;
  padding=1 retains an explicit wrap-around failure. No sponge is used.
- One four-grid timing run found dense Full apply exponent 1.937 versus Fast
  warm apply .854; at N=4097, warm speedup is 140.9x and including-setup speedup
  is 14.73x. N=513 cold speedup is .897x: setup can outweigh benefits on small
  grids. These are finite single-host benchmarks, not asymptotic or real-time
  claims; memory values are array payloads, not peak RSS.
- Reports: `reports/spectral_fast_kirchhoff_prototype_report.md` and its method
  note. Results: `results/validation/spectral_fast_kirchhoff_prototype/`.
  Entrypoint: `validate_spectral_fast_kirchhoff_prototype` after setup.
- Recommendation: a fixed PM random-surface pilot is justified as a next
  validation stage. PM/Monte Carlo, wider bands, SSP, 3D, BER and production
  integration remain unstarted/unverified; this task stops at the prototype.
  Incoming theta99 must not be used as an outgoing cutoff: the sigma=1
  control has about 2.76% of reflected propagating aperture-spectrum energy
  outside the incident99 interval, which this full-kernel method retains.

### Full Kirchhoff spectral rank audit（2026-09-30）

- Stage 1 is complete: validation-only saved-data audit in
  `reports/full_kirchhoff_spectral_rank_audit_report.md`, with artifacts in
  `results/validation/full_kirchhoff_spectral_rank_audit/` and entrypoint
  `validate_full_kirchhoff_spectral_rank_audit` after project setup.
- At 4 kHz, deterministic C2 sinusoidal A/K cases and saved Gaussian spectra,
  the residual `exp(i*(k-kz)*eta)` needs weighted ranks 1–2 / 1–2 / 1–3 for
  relative Frobenius errors 1e-2 / 1e-3 / 1e-4 across 15 small audit cases.
  All saved propagating spectral bins are retained for the engineering result.
  Fixed-band sampling refinement has stable ranks. Chi=A*(k-kz(theta99))
  and source energy distribution explain separability better than kA alone;
  K is absent from the uniform-height matrix but still affects full scattering.
- Three complete-support FK contraction checks at sigma=.3 m, Rx=3 m
  (A/K=.01/.1, .2/.1, .02/.47) have weighted rank-2 ratio errors
  1.664e-6, 6.613e-4, 7.177e-6. Exact reconstruction closes to saved FK at
  about 6e-14. No PE, BIE, reference Full FK or production interface changed.
- Recommendation: investigate a Spectral Fast FK prototype next. Full
  outgoing-kernel FFT complexity O(r*N*log(N)), PM/wider-angle behavior,
  evanescent/outgoing bandwidth and performance remain UNVERIFIED. This task
  stops at the audit; no complete fast solver or production integration exists.
- GitHub PR #3 was verified as Draft/open (not merged) during this audit;
  the local main at 02a29d9 already contains the research-branch integration.

### Unified development mainline and project-level skills（2026-09-30）

- `main` is the sole development mainline. The `phase_show` and
  `cl-verfication` research results have both been merged into `main`; the
  integration baseline is commit `e78c24d` (the mainline may have later
  maintenance commits).
- The `cl-verfication` worktree has been removed. Its 1,144 ignored,
  non-Git files (1,389,680,270 bytes) were preserved locally under
  `results/archive/cl-verification-worktree-unmanaged-20260930/`. Original
  `results/` relative directory structure is retained below that archive root.
  The files were checked against a per-file SHA-256 manifest. The archive
  payload is excluded from Git; only its manifest and preservation note are
  tracked at `reports/cl-verification-worktree-archive-20260930/`.
- Start future development, validation, and new worktrees from `main`.
  The project is entering the unified-mainline and project-level Agent Skills
  phase; the existing project instructions remain in `AGENTS.md`. The
  project-local skills are in `.agents/skills/`: `matlab`,
  `scientific-critical-thinking`, `statistical-analysis`,
  `experimental-design`, `scientific-visualization`, `literature-review`,
  and `scientific-writing`. They were copied from the official
  `K-Dense-AI/scientific-agent-skills` repository at commit
  `65d6e786832e2c52832713117bbbf5096b56f77f`; `skills-lock.json` records
  content hashes for reproducibility.

### Bellhop task supersession / 当前任务补充（2026-08-27）

**Latest user decision:** use a physical 1 m source clearance, i.e. water
depth 100 m, Tx depth 99 m, Rx depth 3 m. Stop epsilon-limit investigation.
`generate_bellhop_100m_source1m_visuals_vertical` is the isolated current entry
for 4 kHz / 89.5-degree representative A/E/R figures. Receiver runs use a
12 m range box; the illustrative 1001-ray fan uses 120 m to show returns.
The old zero-clearance results below remain historical and do not certify
or invalidate this changed physical installation. PE defaults are untouched.

The representative run and both PNGs are complete: 193 A arrivals, 199 E
contribution records (all retraced using 167 launch angles), and 1001 fan rays.
All fan rays touch the sea surface; 961 touch the 100 m bottom. The fan has
33 storage-limited trajectories; receiver A/E/R have no such warnings.
E prefixes are not exact point-hit rays (maximum closest-prefix distance
10.0569 m). See `reports/bellhop_100m_source1m_visuals_report.md`; these figures
do not certify rough multipath convergence or equivalence to 90 degrees.

On 2026-09-01 the user archived 12 historical PE–Bellhop validation/plotting
scripts (including the dedicated `pe_bellhop_current_run_meta_vertical` helper)
under the Git-ignored `cash/pe_bellhop_validation_scripts_2026-09-01/` directory.
The exact list and reason are recorded in
`reports/pe_bellhop_validation_archive_manifest_2026-09-01.md`. These files are
not active entrypoints; the archive is not read during normal project searches.

The PE-topology-only display is now provided by
`generate_bellhop_pe_equivalent_ray_paths_vertical`. It selects the unique
0-top/0-bottom and 1-top/0-bottom A arrivals, retraces their launch angles in
R mode, and truncates each stored ray at its closest approach to Rx. At 4 kHz
the delays are 64.0024394 ms and 67.9360554 ms; the closest-point residuals are
approximately 5.61e-9 m and 6.25e-3 m. No bottom-reflected ray is plotted, and
this is a visualization/filtering result rather than a new PE validation.

- The previous zero-clearance normal-coordinate task required mean water depth
  100 m, physical source clearance 0 m and receiver depth 3 m. The proposed
  fluid bottom uses 1800 m/s, 2000 kg/m^3 and 0.5 dB/wavelength attenuation.
- The user accepted water-side epsilon=[1e-3,1e-4,1e-5] m. Read
  `reports/bellhop_100m_rough_surface_report.md` for the current status and
  `reports/bellhop_100m_boundary_source_report.md` for the six 4 kHz trials.
  Direct/one-surface classes passed all four final-pair checks; rough-surface
  higher bottom-bounce classes failed. Angle/frequency expansion is paused,
  not completed. The full channel is NOT approved by this check.
- The installed input reader clamps near-bottom sources using the SSP limit.
  Explicit flat BTY at physical 100 m plus a homogeneous SSP-only extension
  of 0.1 m avoids that clamp without moving the actual bottom. All six source
  preflights passed; requested 1e-5 m is recorded as approximately 8e-6 m.
  The unclamped flat 1e-3 m ARR is byte-identical before/after this input change.
- Final-pair rough total ARR-synthesized response changes by -0.0385556 dB TL
  and +0.0103577 rad. Primary classes passing does not validate coherent totals:
  the flat D+S sum near cancellation changes by +0.0132625 dB/-0.0106057 rad.
  Results are Bellhop-only unit Cartesian line-source diagnostics, not a new
  Gaussian PE comparison. No bottom-inclusive production recommendation yet.
- A follow-up C-vs-ARR audit on rough 89.5 degrees at 4 kHz (epsilon 1e-4 versus
  1e-5) found C and ARR complex-ratio changes within 0.000533 dB TL and
  2.29e-5 rad. Epsilon sensitivity is therefore already present in Bellhop's
  coherent field; ARR single-precision storage/merging is not its sole cause.
  This remains a diagnostic, not full multipath acceptance or a PE comparison.
- The historical 500 m rough-surface audit, full-ray plots and trial outputs
  have been moved to `cash/`. Entries below describing them are historical,
  not active result locations or evidence for the 100 m channel. Do not read
  or auto-recover `cash/`; the external manifest is in the new plan.
- Keep PE/AS, Gaussian-window and unfolded flat-surface evidence. No PE core,
  source definition, communication code or production default was changed.

### Bellhop executable baseline / Bellhop 执行版本（2026-08-31）

- New Bellhop runs use the official OALIB Windows package `2020_11_4`,
  installed at `AcousticsToolbox_2020/`. The executable is
  `windows-bin-20201102/bellhop.exe` (1,273,704 bytes, SHA-256
  `7E7809A64C3BF734AFF6D28D0D4D52B1B4BD203D81676E3241FFD3189941B505`).
- Active validation scripts read the configured `BELLHOP_EXE`; no machine-local
  executable fallback is embedded in the repository.
- The 2017 installation remains in place for reproducibility. Historical
  reports and saved run metadata keep their recorded 2017 path and hashes;
  do not reinterpret those artifacts as 2020 reruns.
- The bundled Munk ray case completed with exit code 0. This installation
  change does not alter PE marching, the Gaussian source, communication code,
  validation thresholds, or previously generated outputs.

### Bellhop 2020 parametric internal-wall validation（2026-09-02）

- The authoritative implementation and validation status is consolidated in
  `reports/bellhop_internal_wall_implementation_report.md`; superseded stage
  reports and obsolete one-off runners have been moved to the recoverable
  repository-root `cash/` quarantine and are not active project context.
- Current status is **PASS_WITH_LIMITS**. Flat, tilted, smooth sinusoidal, a
  circular vertical-tangent wall, the fixed-realization PM density audit and
  the finite-angle fixed-PM local Reflect2D covariance audit pass. Curved overlays use ordered
  non-coincident parametric segments, normalized hit tangents, reconstructed
  TOP normals and signed turning-angle/arclength curvature. They no longer use
  internal-wall `dz/dr`, `Dss`, `Delta z/Delta r`, `1/Delta r`, strict range
  monotonicity or ATI-style infinite endpoint extension.
- The accepted hit calls unchanged native 2D `Reflect2D` once. Native RN/RM,
  p/q, vacuum phase, Amp and tau are retained; the reflected branch then uses
  the proper half-turn `r'=2R0-r, z'=-z` only as a storage-chart map before the
  unchanged `InfluenceGeoHatCart`.
- Flat reproduces `H_wall=-P_BH(103 m)` for the 3-by-3 step/beam scan. Tilted
  geometry is machine-precision. The sinusoidal 36-case scan has maximum phase
  difference `6.71e-5 rad` after explicit SHD receiver pairing. A radius-100 m
  circular wall remains finite and convergent at exact `dr/ds=0`, with
  `u·n=1`, no NaN/Inf and signed curvature converging to `-0.01 1/m`.
- The former `2.095 rad` discrepancy was a wrong SHD receiver-column pairing,
  not a reflection/frame defect. Native backward-range Cartesian amplitude
  remains a `-0.2606` to `-0.286 dB` diagnostic and is not fitted.
- The fixed PM audit uses seed `260001`, `U=6 m/s`, span `160 m`, master `4097`
  and density `513/1025/2049/4097`, all sampled from one realization. Preflight
  has `min |u·n|=0.8022`, `lambda max|kappa|=0.01672`, and no rejected segments;
  hit/tangent/normal/curvature/reflected-direction/RN/path/phase errors decrease
  through `N=2049`, with 2001/2001 rays successful at every density. Details:
  `cash/bellhop_internal_wall_superseded_20260902/reports/bellhop_internal_pm_fixed_realization_convergence_report.md`.
- The old strict-90-degree native ATI versus rotated-returned PM comparison is
  not the same physical problem and is retired. The next allowed stage is a
  PE rough-wall cross-validation using this same fixed, band-limited realization
  and the accepted 89/89.5-degree rigid transforms. Complete receiver-field
  equivalence and Monte Carlo remain unstarted.
- The completed Bellhop-only local reflection covariance validator and its
  fixed-realization convergence validator are archived under
  `cash/bellhop_internal_wall_superseded_20260902/scripts/validation/`.
  Their paired-ray evidence remains under
  `results/validation/bellhop_pm_local_reflection_covariance/` and writes
  the merged authoritative summary in
  `reports/bellhop_internal_wall_implementation_report.md`. Its hard gates stop
  at the native ATI/internal-wall intersection and native `Reflect2D`出口;
  receiver coherent pressure is diagnostic-only. The full stage reports are
  retained in `cash/bellhop_internal_wall_superseded_20260902/reports/`.

### Bellhop internal-wall visualization archive（2026-09-09）

- The internal-wall visualization-only scripts and generated products used in
  the tilted/sinusoidal ray figures and dense two-dimensional TL figure have
  been moved to the recoverable archive
  `cash/bellhop_internal_wall_visualization_archive_20260909/`.
- The archive includes the reflection-overlay renderer, visual-enhancement
  renderer, sparse TL diagnostic, dense TL-grid runner/renderer, the related
  visualization validation entrypoint, and their validation/visualization
  sidecars, fields, PNGs, CSVs and manifests. The exact moved-path list is in
  `reports/bellhop_internal_wall_visualization_archive_manifest.md`.
- This is an artifact-organization change only. It does not modify Bellhop
  internal-wall geometry, native `Reflect2D`, `InfluenceGeoHatCart`, PE,
  communication code, validation binaries, or the PASS_WITH_LIMITS numerical
  conclusions in `reports/bellhop_internal_wall_implementation_report.md`.
  The archived visual products are historical/reproducibility artifacts and
  are not active entrypoints or current validation inputs. Normal searches must
  exclude `cash/`.

### PE rough-PM ↔ Bellhop rough-PM comparison design（2026-09-02）

- The current design authority is
  `reports/pe_bellhop_pm_comparison_design_audit_report.md`; status is
  **FEASIBLE_WITH_MINOR_PREPARATION**. This design supersedes the earlier
  suggestion that the PE comparison itself should use 89/89.5 degrees.
  Those angles remain Bellhop-only covariance/continuity diagnostics; the
  nominal PE comparison should use the now-validated exact-90-degree
  parametric internal wall.
- The first controlled surface is the existing seed-260001, U=6 m/s,
  160 m, band-limited 1-D PM realization. PE uses
  `eta_PE(x,y)=eta_1D(x)` and Bellhop uses the same coefficients/profile;
  neither side may regenerate, recenter, renormalize or separately smooth it.
- Because the production PE is a two-transverse-dimensional Gaussian field
  while Bellhop is 2-D, the primary dimensional bridge must be a
  validation-only one-transverse-dimensional PE/AS result or equivalent
  `k_y=0` projection. The production y-invariant-surface/on-axis PE result is
  a secondary dimensionality sensitivity, not an exact line-source match.
- The primary metric is reflected-only
  `G=H_ref_rough/H_ref_flat`; compare `G_PE/G_BH`. Absolute reflected pressure,
  total channel, individual ray/arrival matching and native backward-range
  amplitude remain diagnostic. The first run is one frequency (4 kHz), one
  on-axis receiver and a flat/rough pair using the strict 192.1875 m
  no-sponge PE baseline; 6/8 kHz and ensembles come only after interpretation.
- Three minor preparations remain: a same-dimension PE validation reference,
  a canonical Fourier-profile mapper with support/seam metadata, and a shared
  flat/rough reflected-ratio extractor including flat/source and constant-
  height sign audits. No PE/Bellhop core physics change is authorized by the
  design review.

### PE--Bellhop PM Stage 0A canonical mapper（2026-09-03）

- Stage 0A is **PASS**. The validation-only entrypoint is
  `scripts/validation/validate_pe_bellhop_pm_canonical_mapper.m`; its support
  functions are `load_fixed_pm_profile_for_pe_bellhop_validation.m` and
  `evaluate_fixed_pm_fourier_profile.m`.
- The mapper reads the existing
  `results/validation/bellhop_internal_pm_fixed_realization/fixed_pm_fourier_coefficients.csv`
  directly. It evaluates the same periodic series and analytic first/second
  derivatives for PE samples and Bellhop wall samples, with no new random
  numbers, recentering, variance scaling, smoothing or tapering.
- The canonical mapping is `eta_PE(x,y)=eta_1D(x)` and
  `Gamma_B(s)=[R0-eta(s),s]` with `R0=100 m`. The master profile cross-check
  errors are `8.38e-15 m` in height, `2.01e-15` in slope and `5.31e-16 1/m`
  in the second derivative; the wall inverse-map error is `7.11e-15 m`.
- The fixed profile statistics from the coefficient evaluator are RMS height
  `0.18639076 m`, RMS/max slope `0.0539798/0.1166743`, RMS/max geometric
  curvature `0.0209047/0.0445944 1/m`, and minimum radius `22.4243 m`.
- Outputs are under `results/validation/pe_bellhop_pm_canonical_mapper/` and
  the stage report is `reports/pe_bellhop_pm_canonical_mapper_report.md`.
  This stage freezes input provenance only; it does not run PE or Bellhop
  propagation and does not authorize skipping Stage 0B--0E.

### PE--Bellhop PM Stage 0B one-transverse-dimensional bridge（2026-09-03）

- Stage 0B is **PASS**. The validation-only entrypoint is
  `scripts/validation/validate_pe_1d_validation_bridge.m`; the 1-D operator is
  implemented in `scripts/validation/support/run_pe_1d_surface_reflection_validation.m`.
- The bridge uses `xw=192.1875 m`, `nx=984`, `dx=0.1953125 m`,
  `yw=50 m`, `ny=256`, sponge off, `f=4 kHz`, `c=1500 m/s`,
  `z_tx=100 m`, `z_rx=3 m`, and `sigma=0.3 m`. It compares the independent
  one-step exact 1-D angular spectrum with the 1-D split-step helper and with
  the `k_y=0` projection of the unchanged production two-transverse PE.
- Flat and a weak deterministic phase-screen case pass. Maximum 1-D versus
  exact-AS complex error is `1.14e-12`; maximum production-PE `k_y=0` versus
  scaled 1-D complex error is `9.48e-13`; flat pressure-release sign error is
  `2.64e-13`. No Bellhop propagation or PE core changes were made.
- Results are under `results/validation/pe_1d_validation_bridge/` and the
  stage report is `reports/pe_1d_validation_bridge_report.md`. Stage 0C flat
  source/normalization audit is the only next permitted stage.

### PE--Bellhop PM Stage 0C flat source/normalization audit（2026-09-03）

- Stage 0C is **PASS**. The validation-only entrypoint is
  `scripts/validation/validate_pe_bellhop_pm_stage0_flat_source.m`; it uses
  the existing Gaussian `.sbp` pattern and the unchanged 1-D PE bridge.
- The audit uses the 97 m direct / 103 m image geometry, explicit SHD range
  selectors, and Bellhop beam counts 5001/10001 at 0.05 m steps. Flat wall
  intersection residual is `3.54e-12 m`, pressure-release phase error is
  `7.11e-15 rad`, and both beam counts produce identical axis-Q values
  (`0.2610 dB` diagnostic TL offset, `-5.15e-4 rad` phase).
- Normalized offset *magnitude* residuals are `5.97e-3` (direct) and
  `5.55e-3` (reflected). Complex offset profiles retain a documented
  transverse phasor-orientation diagnostic (`1.71/1.81`), not a source fit or
  calibration. The known `~0.26 dB` backward-range Cartesian influence offset
  remains diagnostic only.
- A minimal compatibility fix in
  `scripts/validation/support/read_bellhop_shd_unfolded_vertical.m` locates
  receiver depths relative to the detected range record, which is required
  for `.sbp` SHD files. Results are in
  `results/validation/pe_bellhop_pm_stage0_flat_source/` and the report is
  `reports/pe_bellhop_pm_stage0_flat_source_report.md`. Stage 0D and all
  subsequent planned stages are complete below.

### PE--Bellhop PM Stage 0D constant-height sign audit（2026-09-03）

- Stage 0D is **PASS**. The validation-only entrypoint is
  `scripts/validation/validate_pe_bellhop_pm_constant_height_sign.m` and
  tests `eta0=+0.05, 0, -0.05 m` at 4 kHz with the unchanged 1-D PE bridge,
  Gaussian `.sbp`, and the independent internal-flat validation binary.
- The image range is consistently `L=103-2*eta0 m`; PE predicts
  `exp(+i 2 k eta0)`. PE phase errors are below `2.0e-15 rad` and Bellhop
  phase errors are `2.56e-5 rad`; PE--Bellhop ratio phase differences stay
  below `2.56e-5 rad`. Travel-time residuals are `4.16e-17 s`, pressure-
  release phase is applied once, and flat curvature/rotation-state checks pass.
- The validation overlay now accepts translated flat walls subject only to
  positive wall/post-wall ranges; this is not an official Bellhop change.
  Results are under `results/validation/pe_bellhop_pm_constant_height_sign/`
  and the report is `reports/pe_bellhop_pm_constant_height_sign_audit.md`.
  Stage 0E numerical-budget freeze is completed below.

### PE--Bellhop PM Stage 0E numerical-budget freeze（2026-09-03）

- Stage 0E is **PASS**. The validation-only entrypoint is
  `scripts/validation/validate_pe_bellhop_pm_numerical_budget.m`.
- One seed-260001 band-limited Fourier realization is reused across PE
  window/grid/step and Bellhop profile/beam/step scans. PE stays within
  0.1 dB / 0.02 rad and Bellhop stays within 0.1 dB / 0.02 rad, with zero
  grazing hits, wall residual below `8e-15 m`, and zero q-state residual.
  The p-state change is the native curvature kick, not an error.
- The frozen flat/source allowance remains `0.30 dB` for the known
  backward-range Cartesian influence offset. Results are in
  `results/validation/pe_bellhop_pm_numerical_budget/`; report:
  `reports/pe_bellhop_pm_numerical_error_budget.md`.

### PE--Bellhop PM Stage 1A Tier-1（2026-09-03）

- Stage 1A structural audit is **PASS_WITH_LIMITS** via
  `scripts/validation/validate_pe_bellhop_pm_stage1_tier1.m`. It uses the
  fixed seed-260001 realization at 4 kHz and compares reflected-only
  rough/flat ratios without fitting or renormalizing either model.
- Bellhop geometry and native beam-state checks pass (`|u.n|` minimum
  `0.8319`, wall residual `7.21e-15 m`, phase jump error `7.11e-15 rad`,
  q residual and rotation-state errors zero). The PE/Bellhop ratio difference
  (`0.3104 dB`, `-2.2482 rad`) is retained as a cross-model diagnostic, not a
  structural acceptance gate. Results are in
  `results/validation/pe_bellhop_pm_stage1_tier1/`; report:
  `reports/pe_bellhop_pm_stage1_tier1_report.md`.

### Bellhop R/X source-geometry audit（2026-09-09）

- A validation-only single-variable audit now explicitly selects Bellhop
  `RunType(4)=R` (point source) or `X` (line source), leaving the Gaussian
  `.sbp`, geometry, step, receiver positions, Reflect2D, p/q and
  InfluenceGeoHatCart unchanged. Flat internal-wall regression remains exact.
- For the weak sinusoidal native-ATI/internal-wall comparison at 10001 beams,
  `R` gives `-0.260650 dB`, matching the range-normalization prediction
  `10 log10(97/103)=-0.260655 dB`; `X` reduces the discrepancy to
  `4.24e-6 dB / -3.92e-8 rad`. The 5001-beam sinusoidal result has a common
  `3.156 dB` discretization change for both R and X and is not treated as the
  converged endpoint.
- The conditional 4 kHz seed-260001 X-source Tier-1 rerun passes all Bellhop
  structural checks. Its PE--Bellhop delta is `0.310371 dB / -2.248201 rad`,
  versus the prior R result `0.310434 dB / -2.248201 rad`. Thus future
  dimensionally matched 1-transverse PE Tier-1 comparisons should use X, but
  the existing reflection-model discrepancy conclusion is unchanged. Report:
  `reports/bellhop_source_geometry_rx_audit_report.md`; results:
  `results/validation/bellhop_source_geometry_rx_audit/`.
- The active PE--Bellhop comparison entrypoints now request line-source X:
  the incident-plane validator, fixed-PM Tier-1, frequency extension, ensemble
  engine and Stage-4 statistical driver. Ensemble case roots and request
  fingerprints include the source geometry, so historical R caches cannot be
  silently reused as X.
- The incident-plane comparison was rerun after integration with the original
  Gaussian `.sbp`, grids, step, receivers and 5001/10001 beams. All 13 gates
  pass. PE--Bellhop TL P95 over M95 improved from the archived R value
  `0.15794309 dB` to `0.0019036364 dB`; normalized L2 over M99 is
  `0.0061048459`, and the 5001--10001 beam L2 is `5.2520e-6`. The historical
  point-source global constant is not applied to X; the comparison remains
  axis-normalized and uses no fitted amplitude.
- Existing Stage-2/3/4 numerical reports retain their pre-migration R-source
  provenance and were not relabeled. A future rerun will generate new X-tagged
  cases. Integration report:
  `reports/pe_bellhop_line_source_integration_report.md`; archive manifest:
  `reports/pe_bellhop_source_geometry_integration_archive_manifest.md`.

### PE--Bellhop PM Stage 1B/1C（2026-09-03）

- Stage 1B production dimensionality sensitivity is **PASS_WITH_LIMITS** via
  `scripts/validation/validate_pe_bellhop_pm_stage1_dimensionality.m`.
  The same seed-260001 `eta(x)` is copied unchanged over y; `ny=256/512`
  changes by only `4.36e-6 dB / 4.49e-6 rad`. Relative to the 1-transverse
  bridge, the 2-transverse reflected rough/flat ratio differs by about
  `-0.1589 dB / 0.01251 rad` (complex error `0.02196`).
- Stage 1C freezes the interpretation as **PASS_WITH_MODEL_DISCREPANCY** via
  `scripts/validation/validate_pe_bellhop_pm_stage1_interpretation.m` and
  `reports/pe_bellhop_fixed_pm_4khz_comparison_report.md`. The cross-model
  Tier-1 difference is `0.3104 dB`, `-2.2482 rad`, complex error `1.8366`,
  substantially larger than the dimensionality sensitivity. Stage 2 frequency
  extension is therefore permitted; no model fitting or core-physics change is
  implied.

### PE--Bellhop PM Stage 2 frequency extension（2026-09-03）

- Stage 2 is **PASS_WITH_MODEL_DISCREPANCY** via
  `scripts/validation/validate_pe_bellhop_pm_frequency_extension.m`.
  The fixed seed-260001 realization is evaluated at 4/6/8 kHz with PE
  flat/rough and Bellhop flat/internal-PM rough pairs. Each point is finite
  and geometrically valid: wall residual `7.14e-15 m`, zero grazing fraction,
  pressure-release phase error `7.11e-15 rad`, q residual and p/q rotation
  errors zero.
- Cross-model diagnostics are `(0.3104, 0.3083, 0.3048) dB` TL and
  `(-2.2482, 2.9109, 1.7868) rad` phase at 4/6/8 kHz; no three-point group
  delay is inferred. The sweep uses the unchanged Gaussian source formula and
  a predeclared Bellhop ±15° sector with the accepted 5001-beam setting; the
  sector avoids low-amplitude tail termination in the validation overlay and
  is not a source fit. Results are under
  `results/validation/pe_bellhop_pm_frequency_extension/`; report:
  `reports/pe_bellhop_fixed_pm_frequency_extension_report.md`.

### PE--Bellhop PM Stage 3 paired ensemble（2026-09-03）

- Stage 3 is **PASS_WITH_LIMITS** via
  `scripts/validation/validate_pe_bellhop_pm_ensemble.m` with seeds
  `260001:260008` at 4 kHz. Seed 260001 reuses the canonical realization;
  the other paired profiles preserve the canonical per-mode spectral
  amplitudes and use deterministic seed phases, with no Hs renormalization,
  recentering, smoothing, tapering, or bandwidth change.
- All eight Bellhop cases pass finite-field, geometry, positive-range,
  pressure-release phase and p/q-state checks. Wall residuals are
  `7.11--7.16e-15 m`, grazing fraction is zero, and the reported mean
  reflected roughness powers are PE `1.0707` and Bellhop `0.9271`.
- This is a paired fixed-band phase-ensemble statistical smoke, not an
  independently sampled PM amplitude ensemble. The native backward-range
  Cartesian amplitude limitation remains diagnostic. Results are under
  `results/validation/pe_bellhop_pm_ensemble/`; report:
  `reports/pe_bellhop_pm_ensemble_comparison_report.md`.

### PE--Bellhop PM Stage 3B coefficient-amplitude ensemble（2026-09-03）

- Stage 3B is **PASS_WITH_LIMITS** via
  `scripts/validation/validate_pe_bellhop_pm_amplitude_ensemble.m`, using
  seeds `260001:260008`, 4 kHz, and a uniform 5001 Bellhop beams for a reduced
  first independent ensemble.
- Seed 260001 is copied byte-for-byte from the canonical coefficient file;
  the other seeds draw independent zero-mean Gaussian cosine/sine
  coefficients with variance `S(k)*Delta-k` from the same `Sk_m3` band. No
  Hs renormalization, recentering, smoothing, tapering, or bandwidth change
  is applied. All finite-field, geometry, pressure-release phase, beam-state,
  positive-range, request-fingerprint, output-hash, uniform-numerics and
  reflection-success checks pass; every rough case records 5001/5001 wall hits.
- The flat and canonical rough Stage-2 artifacts are reused only after exact
  embedded-configuration matching and input-age checks. Other cached cases
  require a matching full request fingerprint and verified `.shd`/`.iwdiag`
  SHA-256 hashes. A second identical run recovered the seven non-canonical
  cases only as `cache_verified`, with no 1001/5001-beam mixing. The current
  mean model delta TL is `0.320686 dB`, and the circular mean phase difference
  is `-0.689310 rad`.
- The amplitude-ensemble report is
  `reports/pe_bellhop_pm_amplitude_ensemble_report.md`; machine-readable
  results are under `results/validation/pe_bellhop_pm_amplitude_ensemble/`.
  This is a first reduced statistical result, not a converged ocean Monte
  Carlo or an ensemble worst-delay claim.

The complete executable PE--Bellhop PM comparison chain is indexed by
`reports/pe_bellhop_pm_complete_execution_summary.md`. Stages 0A--0E,
1A--1C, 2, 3A and 3B were run independently with the frozen seed-260001
band-limited realization and no PE/Bellhop core-physics or communication
changes. The current rough-wall cross-model result is therefore a controlled
model-discrepancy finding with the Stage-3 fixed-band ensemble limitation,
not a claim of final independent PM Monte Carlo convergence.

Stage 3B 的 validation-only entrypoint
`scripts/validation/validate_pe_bellhop_pm_amplitude_ensemble.m` 已完成
`260001:260008` 的 reduced first ensemble；它使用 canonical `Sk_m3` band
独立抽样 Gaussian cosine/sine coefficient amplitudes，并与 Stage-3A
phase-only smoke 分开保存。结果仍需按 reduced ensemble 限制解读。

### Active core

- Public channel: root `vertical_channel_model.m` compatibility wrapper ->
  `src/channel/vertical_channel_model_impl.m` ->
  `src/propagation/vertical_wape_propagator.m`.
- Surface boundary: `src/surface/pm_surface_boundary_model.m`; public default remains
  `kirchhoff_spatial`.
- Receiver phase reference:
  `src/channel/apply_pe_channel_phase_reference_vertical.m`.
- Fixed-path research acceleration: cached joint-kstat PE, exact discrete
  adjoint receiver projection, PM-grid FFT `C/P`, and conditional receiver
  statistics.
- Communication consumer: root `comm_main_vertical_psk.m` compatibility
  wrapper, `examples/comm_main_vertical_psk_demo.m`, and
  `src/communication/build_communication_taps_vertical.m`.
- `setup_vertical_project.m` is the single idempotent path initializer. Core
  reusable implementations are under `src/`; validation-only helpers are
  under `scripts/validation/support/`.

### Public receiver phase semantics

- Default `paramsV.channel_phase_reference='direct_dsp'`.
- `H_direct_f/H_reflect_f/H_f` now share the direct receiver reference and
  are ready for the MATLAB-IFFT communication path.
- `H_*_reduced_f` preserves the raw PE envelopes; `H_*_physical_f` restores
  absolute longitudinal carriers under `exp(-i*omega*t)`.
- For the standard 100 m Tx / 3 m Rx / 1500 m/s geometry, the nominal direct
  and reflected spans are 97/103 m and the relative reference delay is 4 ms.
  This is not a per-realization PDP peak alignment.
- `legacy_reduced` is retained only for migration and regression.

### Current validation status

- Release candidate `phase_rc_20260722_174945` is **PASS**. The authoritative
  decision and per-gate evidence are in
  `reports/pe_phase_reference_release_candidate_report.md`.
- The independent F=65 phase audit recovered the nominal 4 ms relative delay
  within the delay grid, with operator error `2.2741e-13`, positive-carrier
  phase RMS `3.1007e-11 rad`, and group-delay difference
  `1.2768e-12 ms`.
- The full adjoint suite passed: exact adjoint `7.7684e-15`, receiver
  projection `5.0950e-15`, dense/FFT `C/P`
  `5.7693e-16/1.0052e-15`; F=9 and F=64 analytic/sample ratios to their
  split-sample floors were `0.593078` and `0.850502`.
- U=5 and U=8 full conditional models, the exact-node library, and the
  two-node full communication validation passed with schema 2.x,
  `target_reference='direct_dsp'`, and zero common CIR shift.
- The F=64 public/cached double consistency error was `9.88957e-16`.
- The same-run propagation atlas contains 13 PNG files, one MP4, MAT,
  manifest, and summary. MATLAB static parsing found no RC-source syntax
  errors.
- F=9 (`df=500 Hz`) aliases a 4 ms delay and is not a phase-delay acceptance
  grid. The authoritative delay audit uses F=65 (`df=62.5 Hz`, 16 ms
  unambiguous window).
- Current flat/uniform PE--Bellhop run `bellhop_current_20260723_rc5` is
  **FAIL_CORE / amplitude OPEN**. Phase conversion and component closure were
  at `0` and `3.5762e-18`; PE--analytic/PE--Bellhop maximum delay residuals
  were `0.027282/0.027281 ms`; Bellhop--analytic was
  `3.9304e-06 ms`; coherent/incoherent 5001--10001 beam RMS differences were
  `0.19041/0.10193 dB`. The only core failure was the no-sponge aperture
  validity prerequisite: 32/64 m edge levels were `-1.087/-14.532 dB`, both
  above the required `-40 dB`, and the 32--64 m phase/TL differences were
  `3.526 rad/-3.896 dB`. The current report and ten-figure atlas supersede
  the 2026-07-21 complete report for present-version conclusions.
- Gaussian direct-field aperture study is **CONVERGED** at 160 m/no-sponge.
  The complete reflected chain is strictly **CONVERGED** at actual
  192.1875 m/no-sponge. The current public 50 m/default-sponge configuration
  remains a compatibility default, not the validation recommendation.
- The 4 kHz random-surface aperture set is **PASS** for 192.1875 m in 15/15
  tested realizations. The paired 3--5 kHz random-surface qualification is
  **INCOMPLETE**; no ensemble worst group-delay claim is authorized.
- The 2k phase-coefficient and ocean-spectrum studies are **diagnostic only**.
  They do not validate PE, KStat, realistic ocean scattering, or higher-order
  SSA. Their scripts/results now live under `scripts/validation/` and
  `results/validation/ssa_2k_phase/`.

### Known limits and next checks

- Exact-adjoint/cached v1 remains uniform, CPU double, fixed grid/frequency,
  one nearest-grid receiver, and no bubble/Doppler path.
- The formal adjoint configuration uses PE `128^2` / PM `256^2`; the accepted
  U=8 conditional branch deliberately retains its separately audited
  PE `128^2` / PM `384^2` aperture.
- Bellhop absolute-amplitude source normalization remains open: one global
  direct-path scale gave `1.654 dB` scale spread and reflected TL RMS/max
  `1.629/2.075 dB`; Gaussian source-aware direct-only weighting reduced its
  diagnostic RMS from `0.815` to `0.502 dB` but is not an acceptance path.
- The current 32/64 m no-sponge fields do not meet the `-40 dB` edge-validity
  prerequisite, so they cannot be used to claim PE aperture convergence.
  Sampling/step and sponge location/strength checks themselves passed.
- The atlas emitted one cosmetic MATLAB title-format warning for `\Delta`;
  the saved fields, checks, and artifact counts are unaffected.
- See `vertical_comm_guide.md` for current formulas,
  `scripts/README.md` for commands, and
  `reports/pe_phase_reference_release_candidate_report.md` for the release
  decision. Older combined reduced-envelope phase conclusions are
  **Superseded** by this release-candidate run.
- Bellhop absolute amplitude, point-source infinite-aperture equivalence, and
  source normalization remain **OPEN**. Li2009 remains a partial
  reproduction without completed transverse-grid convergence.
- Repository layout changed on 2026-08-20 without changing public call names,
  formulas, defaults, or result schemas. Current code paths in this index and
  the technical guide supersede root-level paths written in older log entries.

## Detailed Snapshot and Historical Reference

> Entries below this point predate or accompany later changes. When a phase
> or output-semantic statement conflicts with the current-state index above,
> the 2026-08-20 index and latest dated log entries are authoritative.

## Current Snapshot / 当前状态摘要

- Project role: vertical underwater acoustic random-channel generation and communication validation platform.
- Propagation structure: the channel is decomposed into direct propagation, sea-surface reflection/scattering, and the total frequency response:
  - `H_f = H_direct_f + H_reflect_f`
  - `h_total = H_f(idx_f_ref)`
- Current propagation engine: the project still uses the vertical WAPE/PE-style split-step propagation path for acoustic channel generation. Surface boundary models only replace the sea-surface reflection/scattering operator between the upward incident march and the downward reflected march.
- Current external propagation cross-validation: `validate_pe_bellhop_flat_surface_vertical.m` runs a deterministic uniform-SSP, flat pressure-release surface case through the PE public API and the external Bellhop executable. Historical runs restored PE carriers only in the validator; since 2026-07-22 the public API exposes explicit reduced, direct-DSP, and absolute-physical fields, and validators consume those fields directly.
- Bellhop rough-surface near-vertical limit audit (2026-08-26): `validate_bellhop_rough_surface_near_vertical_limit.m` is a validation-only Bellhop entrypoint. It fixes one deterministic three-cosine 1-D altimetry profile (`CVW *`/`.ati`), uses a Cartesian line source, and changes only the Tx--Rx horizontal offset for `85, 87, 88, 89, 89.5, 89.75, 89.875, 89.9375, 89.96875, 89.99` degrees at 4 kHz. Direct and one-top-bounce arrivals plus `.ray` reflection diagnostics are saved under `results/validation/bellhop_rough_surface_near_vertical/`. The selected primary comparison excludes bottom bounces; the unfiltered arrival files do contain bottom/multiple-bounce paths. The near-vertical coherent reflected response is not continuous under the configured `0.1 dB/0.02 rad/10 us` checks and the 89.99-degree case is not captured as a one-top-bounce arrival by the installed Bellhop build. Therefore no finite angle is currently approved as a strict 90-degree substitute for this rough profile; this result does not change PE marching or production defaults.
- Complete Bellhop path display (2026-08-27): `plot_bellhop_rough_surface_full_paths.m` adds all receiver-contributing E records (including bottom/multiple bounces), independently retraced in R mode, and the entire 1001-ray emitted fan in a common 12 m / 500 m box. Historical launch limits are recovered from original `.prt` logs; A replay must reproduce original arrivals. E beams are not exact point-intersecting rays and may be merged in A output. New figures/CSV/MAT live in `bellhop_rough_surface_near_vertical/full_paths/`; original figures, PE code, defaults and validation status are preserved.
- Current communication policy after the D2 update:
  - receive-window mode defaults to `peak_sync`, which aligns the effective equalizer taps to the dominant baseband tap and then keeps the first `N` samples from a full convolution;
  - Eb/N0 noise reference defaults to `rx_clean`, i.e. receiver-side clean waveform power.
- `same_legacy` / `same` receive slicing is retained only as a historical or diagnostic mode. It is not the current default communication flow.
- Current sea-surface model set:
  - `kirchhoff_spatial`: realization-based Kirchhoff phase-screen path and default surface model;
  - `kirchhoff_kdomain`: wavenumber-domain interface for the same explicit Kirchhoff phase-screen realization model;
  - `kirchhoff_kstat`: Kirchhoff statistical phase-screen branch. It does not generate a concrete sea surface; it builds coherent reflection, incoherent scattering power, optional random reflected spectra, and energy diagnostics from the PM height spectrum;
  - `ssa_stat_kernel` with `pm_convolution`: PM-spectrum engineering baseline statistical kernel;
  - `ssa_stat_kernel` with `ssa1_geometry`: first-order pressure-release / Dirichlet geometry kernel.
- Current Kirchhoff K-Stat status:
  - implemented as an independent branch under `pm_surface_boundary_model.m`, selected by `surface_boundary_model='kirchhoff_kstat'`;
  - uses near-vertical phase-screen factor `alpha=2*k0`, coherent reflection `R_coh=-exp(-2*k0^2*sigma_eta^2)`, and the incoherent phase-screen spectrum `S_deltaG` rather than the total phase-screen spectrum with the coherent delta peak included;
  - metadata records `sigma_eta2_m2`, `G_mean`, `R_coh`, `S_deltaG`, `P_sca`, deterministic random seed fields, FFT normalization notes, full-K energy closure, and propagation/trusted-window energies;
  - disabled/default paths keep compatibility and do not replace the existing default `kirchhoff_spatial` behavior.
- Current roughness scaling policy:
  - `surface_roughness_scale_mode='target_hs'` remains the default and scales the PM spectrum or explicit surface realization to `sea_hs_target`;
  - `surface_roughness_scale_mode='raw_pm'` disables target-Hs scaling, so `sea_wind_speed` sets both PM spectral shape and integrated roughness strength;
  - raw-PM comparisons use the project's discrete spectrum variance convention, `sigma_eta_raw^2=sum(Phi2D(:))*dkx*dky` and `Hs_raw=4*sigma_eta_raw`;
  - in raw-PM mode, `kirchhoff_kstat` maps the project PM spectrum to its continuous correlation formula with `W_eta=Phi2D*(2*pi)^2`, while explicit `kirchhoff_kdomain` uses the current `sqrt(2)*real(ifft2(...))` realization convention.
- Current SSA statistical kernel status:
  - PM-spectrum-driven statistical scattering is connected to the existing reflected PE path and communication `H_f` consumer;
  - `periodic` and `zero_padded` convolution paths are implemented for `pm_convolution` and `ssa1_geometry`;
  - energy audit, compact metadata, Hs=0 degeneration, deterministic seeding, scale monotonicity, dense-vs-FFT periodic validation, and reduced-grid smoke tests have passed;
  - optional frequency-correlated random scattering modes are available for wideband stochastic-channel diagnostics, while the default remains independent per-frequency scatter;
  - reduced-grid scatter-scale calibration and compact applicability reporting scripts are available as engineering diagnostics;
  - optional `ssa2_broschat_coherent` is implemented only as a coherent-reflection sensitivity diagnostic; second-order incoherent scattering is not implemented.
- Current K-Stat / K-Domain validation status:
  - `validate_kirchhoff_kstat_vertical.m` validates the standalone kstat branch for flat-surface degeneration, coherent formula, phase-screen energy closure, propagation-window diagnostics, deterministic seeds, explicit Kirchhoff mean comparison, and weak SSA1 coherent-reference consistency;
  - `compare_kirchhoff_kdomain_kstat_wind_vertical.m` compares explicit and statistical Kirchhoff branches under raw-PM wind sweeps. With the discrete-variance alignment, the full reduced run reported maximum `Hs_raw_rel_error` about `0.0145`, maximum kstat energy-closure error about `2.22e-15`, and `H_f=H_direct_f+H_reflect_f` at about the `1e-17` level;
  - `validate_kstat_vs_kdomain_phase_screen_vertical.m` is a surface-boundary-only plane-wave ensemble check. The latest full reduced run (`f=[4000 6000 8000] Hz`, `Hs=[0.05 0.1 0.2 0.5] m`, `M=[8 16 32 64]`, `128x128`) reported maximum kstat energy-closure error `4.663e-15`, maximum `sigma_eta` relative error `4.163e-16`, final-M mean radial-spectrum L2 error `0.0213`, and final-M mean radial-spectrum correlation `0.9997`;
  - `validate_kstat_vs_kdomain_channel_stats_vertical.m` extends the comparison to the full `vertical_channel_model` PE/WAPE path and receiver-side statistics: `E[|h_ref|^2]`, peak-synced reflected/total PDP ensemble means, `|h_total|` distributions, direct-path consistency, and `H_f` invariants. The latest smoke run (`Hs=[0.1 0.2] m`, `seed_count=8`, `64x64`, `Nf=8`) passed hard checks with maximum `H_f` invariant error `7.76e-18`, maximum kstat energy error `9.99e-16`, maximum `E[|h_ref|^2]` relative error `0.0973`, reflected PDP correlations at least `0.979`, and `|h_total|` quantile NRMSE at most `0.193`; reflected PDP L2 exceeded the nominal statistical target at `0.313` for `Hs=0.2`, so larger ensembles are still needed for the full reduced conclusion;
  - `validate_kstat_vs_kdomain_lfm_channel_vertical.m` is an independent channel-level LFM waveform validation script. It does not call `comm_main_vertical_psk.m` and does not use PSK, BER/SER, noise injection, or equalization. It generates an analytic baseband LFM test signal, interpolates `H_f` and `H_reflect_f` onto the LFM FFT grid, computes received total/reflected LFM waveforms and matched-filter outputs, and compares ensemble receiver statistics between `kirchhoff_kdomain` and `kirchhoff_kstat`. The latest reflected-focused run (`Hs=0.2 m`, `seed_count=32`, `64x64`, `Nf=32`, LFM `4-8 kHz`, `fs=24 kHz`, `T=20 ms`) passed hard checks with maximum `H_f` invariant error `1.39e-17`, maximum kstat energy error `1.11e-16`, and direct-path branch/seed delta `0`. However reflected-only agreement did not pass: `E[|h_ref|^2]` relative error was `0.580`, reflected LFM envelope correlation `0.0279`, reflected PDP correlation `0.276`, and reflected matched-filter peak relative error `0.743`. Total-channel metrics were much better (`rx_total_env_corr=0.982`, `pdp_total_corr=0.930`), so total LFM can mask reflected-path disagreement. `plot_lfm_tx_rx_comparison_vertical.m` reads the saved LFM result and generates normalized TX-vs-RX total/reflected waveform overlays for direct visual inspection;
  - single-seed `kirchhoff_kdomain` fields are not expected to match `kirchhoff_kstat` pointwise. The intended comparison is ensemble coherent mean, incoherent power spectrum, radial spectral shape, energy closure, and receiver-side statistics after PE/WAPE if propagation is included.
- Optional SSA component diagnostics can decompose the surface-reflected field into coherent specular and incoherent scatter components, march both components to the receiver, and report whether the scatter-only received contribution is negligible for the current settings. This path is disabled by default and does not alter `H_reflect_f`.
- A surface-only plane-wave coherent-reflection diagnostic is available in `compare_plane_wave_surface_coherent_reflection_vertical.m`. It bypasses PE/WAPE and communication code, uses a vertical plane wave at the sea surface, and compares raw-PM SSA1 coherent reflection with Kirchhoff phase-screen ensemble coherent reflection under wind-speed sweeps.
- Optional reference-frequency surface-wavefield diagnostics can record:
  - the incident center slice from `z_tx` to the sea surface;
  - the reflected center slice from the sea surface to `z_rx`;
  - the incident and reflected complex fields on the sea-surface plane.
  This path is disabled by default and does not change channel values.
- Current physics boundary:
  - this is not a complete SSA/NLSSA, T-matrix, impedance-boundary, multiple-scattering, or experimentally calibrated sea-surface scattering model;
  - `kirchhoff_kstat` is a near-vertical Kirchhoff/Gaussian statistical phase-screen model, not a full rough-surface scattering theory and not a replacement for calibrated sea-surface scattering measurements;
  - the wideband frequency correlation of statistical scattering realizations is still simplified and needs a physically or empirically calibrated model;
  - under fixed `sea_hs_target` normalization, changing wind speed mainly changes PM spectral shape, not a monotonic sea-state intensity by itself;
  - in the plane-wave raw-PM diagnostic, no `sea_hs_target` normalization is applied, so wind speed changes both the PM spectral shape and the integrated roughness variance.

## Repository Role of Each Main MATLAB File

### `vertical_channel_model.m`
- Public channel entrypoint: `output = vertical_channel_model(paramsV)`.
- Accepts user-facing configuration in `paramsV`.
- Normalizes and validates parameters through `local_prepare_config`.
- Calls `vertical_wape_propagator(cfg)` to perform propagation.
- Packages all outputs into a stable `output` struct.
- Enforces 1/R validation in uniform-medium mode when `enforce_1_over_R=true`.

### `vertical_wape_propagator.m`
- Core upward-marching WAPE propagation engine.
- Builds transverse grids, source field, spectral operators, and absorption profile.
- Resolves receiver state, frequency axis, optional GPU path, and optional reflection path.
- Produces direct-path, reflected-path, and total complex channel responses:
  - `H_direct_f`
  - `H_reflect_f`
  - `H_f`
- Exposes reference-frequency scalar gains:
  - `h_direct`
  - `h_reflect`
  - `h_total`

### `pm_surface_boundary_model.m`
- Rough-surface boundary and statistical scattering submodule.
- Supports concrete Pierson-Moskowitz surface realizations for Kirchhoff-style spatial reflection.
- Supports PM-spectrum statistical scattering through `ssa_stat_kernel`, including `pm_convolution`, `ssa1_geometry`, and small-grid `ssa1_debug_dense`.
- Applies the selected surface boundary operator to the incident surface field.
- Supports:
  - reflection coefficient control
  - normal-incidence phase mode
  - oblique phase mode using TX/RX geometry
  - spatial-domain and implicit wavenumber-domain Kirchhoff phase-screen boundary models
  - SSA-like statistical scatter metadata and energy audit fields
- Returns reflected field and reflection/scattering metadata.

### `explain_main_vertical.m`
- Channel-only demonstration script.
- Builds a scalar-frequency `paramsV`.
- Calls `vertical_channel_model(paramsV)`.
- Saves the output struct and standard figures.

### `comm_main_vertical_psk.m`
- End-to-end communication demonstration script.
- Builds a wideband channel configuration.
- Calls `vertical_channel_model(paramsV)` for each scenario.
- Converts acoustic frequency response to a discrete baseband channel.
- Runs MPSK modulation, channel convolution, noise injection, equalization, and BER/SER statistics.

### `modem_psk.m`
- Utility wrapper for M-PSK modulation, demodulation, and error-rate measurement.
- Supported modes:
  - `modulate`
  - `demodulate`
  - `error_rate`

### `noise_inject_vertical.m`
- Utility wrapper for noise injection at complex baseband.
- Current built-in model: `awgn`.
- Also supports a custom extension hook via `custom_noise_fn`.

## Main Execution Paths

### Path 1: Channel Demo
`explain_main_vertical -> vertical_channel_model -> vertical_wape_propagator -> output struct + saved figures`

### Path 2: Communication Demo
`comm_main_vertical_psk -> vertical_channel_model -> vertical_wape_propagator -> H_f/f_axis -> baseband channel -> MPSK link simulation`

## Channel Data Flow

### Step 1: User parameters
- The user-facing configuration starts as `paramsV`.
- `paramsV` may contain both core physics fields and future extension fields.

### Step 2: Runtime configuration
- `vertical_channel_model` calls `local_prepare_config(paramsV)`.
- `local_prepare_config`:
  - merges defaults
  - copies user fields into `cfg`
  - keeps backward compatibility for some legacy names
  - validates geometry, grid, bandwidth, flags, and surface settings
  - derives:
    - `lambda0`
    - `dx`
    - `dy`
    - `dz_abs`
    - `path_span`
    - `numstep`
    - `dz_step`

### Step 3: Receiver state resolution
- `vertical_wape_propagator` initializes `rx_state_used` from:
  - `x_rx`
  - `y_rx`
  - `z_rx`
- If `cfg.rx_position_fn` is not empty, it overrides the nominal RX location at `t=0`.
- `rx_state_used` then becomes the actual receiver state used by the propagation run.

### Step 4: Frequency axis resolution
- `local_resolve_frequency_axis(cfg, rx_state_used)` defines `f_axis`.
- Cases:
  - explicit sweep: if `cfg.f0` is already a vector
  - wideband auto-grid: if `enable_wideband=true`
  - scalar-frequency mode otherwise
- `idx_f_ref` selects the reference frequency used for scalar outputs and stored field slices.

### Step 5: Direct-path propagation
- For each frequency in `f_axis`, the code:
  - computes `lambda_f`, `k0`, and `dz_step_f`
  - constructs a split-step spectral propagator
  - marches the source field from `z_tx` to `z_rx`
  - samples the received field at `(ix_rx, iy_rx)`
- The result is stored as `H_direct_f(ifq)`.

### Step 6: Reflected-path propagation
- If `enable_surface_reflection=true`, the code performs:
  1. `tx -> surface` propagation via `local_march_field`
  2. rough-surface reflection via `pm_surface_boundary_model`
  3. `surface -> rx` propagation via `local_march_field`
- The result is stored as `H_reflect_f(ifq)`.
- The selected surface boundary model acts only in step 2; it does not change the two propagation segments.

### Step 7: Total frequency response
- At each frequency bin:
  - `H_f(ifq) = H_direct_f(ifq) + H_reflect_f(ifq)`
- After the frequency loop:
  - `h_direct = H_direct_f(idx_f_ref)`
  - `h_reflect = H_reflect_f(idx_f_ref)`
  - `h_total = H_f(idx_f_ref)`

### Step 8: Output packaging
- `vertical_channel_model` returns a stable `output` struct containing:
  - field snapshots
  - centerline diagnostics
  - rough-surface products
  - frequency-domain responses
  - reference-frequency scalar gains
  - receiver state and Doppler placeholder result
  - config echo and derived scalar diagnostics

## Current Communication Data Flow

### Step 1: Scenario construction
- `comm_main_vertical_psk.m` defines two scenarios:
  - `direct_only`
  - `direct_plus_reflect`
- The only scenario change is `enable_surface_reflection`.

### Step 2: Symbol generation
- Random bits are generated in the script.
- `modem_psk('modulate', bits_tx, M)` maps them to unit-power M-PSK symbols.

### Step 3: Channel acquisition
- Each scenario calls `channel = vertical_channel_model(paramsV)`.
- The communication chain mainly consumes:
  - `channel.f_axis`
  - `channel.H_f`
  - `channel.idx_f_ref`
  - `channel.h_direct`
  - `channel.h_reflect`
  - `channel.h_total`

### Step 4: Baseband channel construction
- `local_build_baseband_response`:
  - centers `H_f` around the reference frequency
  - builds a baseband frequency axis
  - interpolates the acoustic response onto that axis
- `local_build_channel_taps`:
  - converts the baseband frequency response to time-domain taps
  - truncates taps to a target cumulative energy ratio

### Step 5: Signal propagation at baseband
- The current default receive-window policy is `peak_sync`:
  - find the dominant tap of `h_bb`;
  - discard pre-peak taps for equalizer/channel application;
  - compute a full convolution with the effective taps;
  - keep the first `N` received samples as the symbol-aligned `rx_clean`.
- `same_legacy` remains available only as a diagnostic/historical slicing mode.

### Step 6: Noise injection
- `noise_inject_vertical` adds noise to `rx_clean`.
- In the current demo:
  - model: `awgn`
  - control variable: `Eb/N0`
  - default reference signal: `rx_clean`
  - SNR conversion uses `bits_per_symbol`

### Step 7: Equalization and decisions
- `local_mmse_equalize` performs one-shot frequency-domain MMSE equalization using the known `h_bb`.
- `modem_psk('demodulate', rx_eq, M)` performs coherent hard-decision demodulation.
- `modem_psk('error_rate', bits_tx, bits_rx, M)` computes BER and SER.

### Step 8: Results packaging
- Each scenario stores:
  - physical channel scalars
  - acoustic frequency response
  - derived baseband response
  - BER/SER curves
  - effective SNR metadata
  - the full `channel` struct

## Important Interfaces and Extension Points

### `paramsV` / `cfg`
- New fields can be introduced through `paramsV`.
- `local_prepare_config` copies user fields into `cfg`.
- Existing validated fields must remain compatible with current callers.

### `rx_position_fn`
- Expected form:
  - `rx_xyz = rx_position_fn(t_s, state)`
- Current call site uses `t_s = 0`.
- Must return `[x_rx, y_rx, z_rx]`.

### `doppler_fn`
- Expected form:
  - `fd_hz = doppler_fn(t_s, tx_state, rx_state, env_state)`
- Current call site uses `t_s = 0`.
- Current result is recorded in `fd_hz_used` only.
- No Doppler compensation is currently applied to `H_f` or the communication chain.

### Sound-speed environment
- `env_mode='uniform'`: uses `c0`.
- `env_mode='layered'`:
  - first tries `cz_func(z)`
  - then tries `cz_table_z/cz_table_c`
  - otherwise uses a built-in fallback profile

### Surface reflection controls
- Main controls in `cfg`:
  - `enable_surface_reflection`
  - `sea_wind_speed`
  - `sea_hs_target`
  - `sea_seed`
  - `surface_reflect_coeff`
  - `surface_phase_mode`
  - `surface_oblique_clip`
  - `surface_boundary_model`
  - `surface_boundary_check_equivalence`
  - `surface_boundary_debug`
  - `surface_boundary_equivalence_tol`
  - `surface_boundary_coupling_diagnostics`
  - `surface_boundary_coupling_debug`
  - `surface_boundary_redistribution_diagnostics`
  - `surface_boundary_redistribution_debug`

### Kirchhoff surface boundary interface
- Default surface model remains Kirchhoff, but the rough-surface phase convention has been corrected:
  - `surface_boundary_model='kirchhoff_spatial'`
  - `delta_phi = k0*(cos_i+cos_r)*eta`
  - normal incidence gives `delta_phi=2*k0*eta`
  - `G_xy = surface_reflect_coeff * exp(1i*delta_phi)`
  - `psi_ref_xy = G_xy .* psi_inc_xy`
- Optional k-domain interface:
  - `surface_boundary_model='kirchhoff_kdomain'`
  - `Psi_inc_k = fft2(psi_inc_xy)`
  - `Psi_ref_k = fft2(G_xy .* ifft2(Psi_inc_k))`
  - `psi_ref_xy = ifft2(Psi_ref_k)`
- This is the discrete FFT product-convolution form of the same Kirchhoff phase screen. It represents an implicit periodic-grid operator with kernel proportional to `G_hat(K-K')`, but the code does not build a four-dimensional dense `B_xi(K,K')` matrix.
- FFT convention:
  - MATLAB `fft2` is unnormalized.
  - MATLAB `ifft2` includes `1/(nx*ny)`.
  - `kx/ky = [0:N/2-1,-N/2:-1]*2*pi/L`.
  - The implicit convolution is circular on the transverse periodic grid.
- Added `roughness_meta` fields:
  - `boundary_model`
  - `boundary_operator_form`
  - `boundary_dense_matrix_used`
  - `boundary_fft_convention`
  - `boundary_equivalence_error`
  - `boundary_coupling_diagnostics`
  - `boundary_coupling_debug`
  - `boundary_redistribution_diagnostics`
  - `boundary_redistribution_debug`
  - `boundary_debug_stats`
- When `enable_surface_reflection=false`, `roughness_meta.enabled=false` and the boundary metadata is still present with `boundary_operator_form='not_executed'`.
- Current limits:
  - This is not a T-matrix, SSA, or NLSSA rough-surface solver.
  - It is not a fast statistical channel generator.
  - It keeps the PM sea synthesis, reflection coefficient, and two-segment reflected path unchanged.
  - The pre-fix `delta_phi=2*k0*(cos_i+cos_r)*eta` convention is retained only in standalone diagnostics as a legacy comparison.

### Boundary screen coupling diagnostics
- `surface_boundary_coupling_diagnostics=false` by default. When enabled, the code only adds metadata and does not change the reflected field or channel outputs.
- Diagnostics are evaluated on the retained reference frequency only. In wideband runs, `roughness_meta.boundary_coupling_diagnostics` therefore describes `f_axis(idx_f_ref)`.
- Boundary screen:
  - `G_xy = surface_reflect_coeff .* exp(1i*delta_phi)`
- Screen Fourier spectrum:
  - `G_hat_k = fft2(G_xy)`
  - `P_k = abs(G_hat_k).^2`
- Energy metrics:
  - `E_total = sum(P_k(:))`
  - `E_zero = P_k(1,1)`
  - `E_nonzero = E_total - E_zero`
  - `nonzero_power_fraction = E_nonzero / max(E_total, eps)`
- Coupling-radius metrics:
  - `DeltaK = sqrt(KX.^2 + KY.^2)`
  - `rms_delta_k = sqrt(sum(DeltaK(:).^2 .* P_k(:)) / max(E_total, eps))`
  - `rms_delta_k_nonzero` uses only nonzero wavenumber bins and normalizes by `E_nonzero`.
  - `energy_radius_50/90/95` are weighted cumulative radii over all spectrum bins.
  - `nonzero_energy_radius_90` is the 90 percent weighted radius after excluding the zero-wavenumber bin.
- Interpretation:
  - Larger `nonzero_power_fraction` and `rms_delta_k` indicate a less uniform Kirchhoff screen and stronger potential off-diagonal `K' -> K` coupling in the implicit `G_hat(K-K')` convolution.
  - These metrics are diagnostic summaries of the Kirchhoff phase-screen spectrum only. They are not strict boundary-operator matrix elements, scattering cross sections, T-matrix terms, SSA/NLSSA terms, or a statistical channel generator.
- `surface_boundary_coupling_debug=true` adds small summaries only:
  - 16-bin radial spectrum energy fractions.
  - top 8 nonzero spectrum peaks with `kx`, `ky`, `DeltaK`, and energy fraction.
  - Full `G_hat_k` and `P_k` arrays are not stored.

### Incident-weighted redistribution diagnostics
- `surface_boundary_redistribution_diagnostics=false` by default. When enabled, the code compares the current incident angular spectrum with the reflected angular spectrum generated by the same Kirchhoff phase screen.
- The diagnostic does not change `psi_ref`, `H_reflect_f`, `H_f`, or communication-chain inputs.
- Diagnostics are evaluated on the retained reference frequency only. The sweep script described below runs separate scalar-frequency cases when a frequency trend is needed.
- Spectra:
  - `Psi_inc_k = fft2(psi_inc_xy)`
  - `Psi_ref_k = fft2(psi_ref_xy)`
  - `Psi_flat_ref_k = fft2(surface_reflect_coeff .* psi_inc_xy)`
  - `P_k = abs(Psi_k).^2`
  - `DeltaK = sqrt(KX.^2 + KY.^2)`
- Spectrum moments:
  - `centroid_kx = sum(KX(:).*P_k(:)) / max(sum(P_k(:)), eps)`
  - `centroid_ky = sum(KY(:).*P_k(:)) / max(sum(P_k(:)), eps)`
  - `rms_delta_k = sqrt(sum(DeltaK(:).^2.*P_k(:)) / max(sum(P_k(:)), eps))`
  - `energy_radius_90` is the weighted cumulative `DeltaK` radius containing 90 percent of `P_k`.
- Reported redistribution metrics:
  - `incident_rms_delta_k_rad_per_m`
  - `reflect_rms_delta_k_rad_per_m`
  - `rms_delta_k_increase_rad_per_m`
  - `centroid_shift_kx/ky/mag_rad_per_m`
  - `reflect_high_k_fraction`, using `incident_energy_radius_90` as the high-k threshold
  - `flat_ref_rms_delta_k_rad_per_m`
  - `rough_vs_flat_rms_delta_k_increase_rad_per_m`
  - `rough_vs_flat_energy_radius_90_increase_rad_per_m`
  - `rough_vs_flat_high_k_fraction_increase`
- `surface_boundary_redistribution_debug=true` adds compact radial energy summaries and top 8 spectral peaks for incident, reflected, and flat-reflected spectra. Full `Psi_inc_k`, `Psi_ref_k`, and `Psi_flat_ref_k` arrays are not stored.
- Interpretation:
  - Larger reflected `rms_delta_k`, high-k fraction, and rough-vs-flat increase indicate stronger angular-spectrum broadening for the current incident field.
  - These are incident-weighted Kirchhoff phase-screen diagnostics only. They are not strict off-diagonal operator entries, T-matrix elements, SSA/NLSSA terms, scattering cross sections, or fast statistical channel parameters.

### Surface redistribution trend script
- `sweep_surface_boundary_redistribution_vertical.m` runs reduced scalar-frequency diagnostics with:
  - `nx=ny=128`, `xw=yw=50`, `z_tx=100`, `z_rx=3`, `sigma_src_m=0.4`
  - `show_figures=false`, `save_mode='rx_only'`, `enforce_1_over_R=false`
  - `surface_boundary_model='kirchhoff_kdomain'`
  - `surface_boundary_redistribution_diagnostics=true`
- It scans:
  - `sea_hs_target = [0.05, 0.5, 1.0]`
  - `sea_wind_speed = [3, 5, 8, 12]`
  - scalar `f0 = [4000, 6000, 8000]`
- It saves `sweep_surface_boundary_redistribution_vertical_result.mat` and lightweight PNG trend plots for RMS broadening and high-k fraction.
- Wind-speed sweep interpretation requires caution because the PM surface is rescaled to fixed `sea_hs_target=0.5`; changing wind speed changes the spectral shape after scaling, not a monotonic wave-height amplitude.

### Surface Monte Carlo channel statistics script
- `monte_carlo_surface_channel_vertical.m` runs reduced wideband Monte Carlo statistics over random Kirchhoff rough-surface realizations.
- The script fixes all physical and numerical parameters except `sea_seed`.
- Default reduced setup:
  - `mc_count=16`, overridable through environment variable `SURFACE_MC_COUNT`.
  - `seed_list=12345+(0:mc_count-1)`.
  - `nx=ny=128`, `xw=yw=50`, `z_tx=100`, `z_rx=3`, `sigma_src_m=0.4`.
  - `enable_wideband=true`, `f_band_hz=[4000 8000]`, `Nf_min=Nf_max=16`, `f_ref_hz=6000`.
  - `surface_boundary_model='kirchhoff_kdomain'`.
  - `surface_boundary_coupling_diagnostics=true`.
  - `surface_boundary_redistribution_diagnostics=true`.
  - `show_figures=false`, `save_mode='rx_only'`, `enforce_1_over_R=false`, `use_gpu=false`.
- Channel ensemble formulas:
  - `H_m(f)` is one Monte Carlo realization.
  - `mean_H(f) = (1/M) * sum_m H_m(f)`.
  - `var_H(f) = (1/M) * sum_m abs(H_m(f)-mean_H(f)).^2`.
  - Amplitude statistics are computed from `abs(H_m(f_ref))`.
  - Principal phase samples use `angle(H_m(f_ref))`.
  - Circular mean phase is `angle(mean(exp(1i*theta_m)))`.
  - Circular phase variance is `1 - abs(mean(exp(1i*theta_m)))`.
  - Reference reflection/direct ratio is `h_reflect_m / max_complex(h_direct_m, eps)`.
- Diagnostic ensemble statistics:
  - Mean, variance, standard deviation, and 5/25/50/75/95 percent quantiles are reported for selected C2/C2.5 metrics.
  - Coupling metrics include `nonzero_power_fraction` and `rms_delta_k_rad_per_m`.
  - Redistribution metrics include `reflect_rms_delta_k_rad_per_m`, `reflect_high_k_fraction`, and `rms_delta_k_increase_rad_per_m`.
- Tap metrics:
  - The script builds a baseband response by interpolating `H_f` around `idx_f_ref`.
  - It computes `h_bb = ifft(ifftshift(H_baseband))`, keeps taps until the requested energy ratio is reached, and reports compact delay/peak metrics.
  - Default tap controls are `symbol_rate_hz=1000`, `tap_fft_len=512`, and `tap_energy_ratio=0.999`.
- Outputs:
  - `monte_carlo_surface_channel_vertical_result.mat` with `params_base`, `seed_list`, `channels_summary`, `summary_table`, and `mc_stats`.
  - Compact numeric sample arrays for `H_f`, `H_direct_f`, `H_reflect_f`, and reference-frequency scalars are stored; full field arrays and full 2-D spectra are not saved.
  - PNG summaries: `H_f` mean magnitude envelope, `abs(H(f_ref))` histogram, phase histogram, reflected RMS delta-k histogram, high-k fraction histogram, and tap RMS delay histogram.
- Scope:
  - This is empirical Monte Carlo over the existing Kirchhoff phase-screen channel. It is not a new scattering formula, a closed-form stochastic channel generator, a BER/SER simulation, a T-matrix method, SSA/NLSSA, or a calibrated scattering cross-section model.

### Multi-sea-state Monte Carlo sweep script
- `sweep_monte_carlo_surface_channel_vertical.m` extends the C3 Monte Carlo workflow from one fixed sea state to an `Hs x wind` grid.
- Default reduced sweep:
  - `sea_hs_target = [0.05, 0.5, 1.0]`.
  - `sea_wind_speed = [3, 5, 8, 12]`.
  - `mc_count=16`, overridable through `SURFACE_MC_SWEEP_COUNT`.
  - `SURFACE_MC_SWEEP_MAX_CONDITIONS` can truncate the condition list for smoke tests.
  - `SURFACE_MC_SWEEP_RESULT_FILE` can redirect the `.mat` output for smoke tests.
  - `SURFACE_MC_SWEEP_FIGURE_PREFIX` can redirect PNG output prefixes for smoke tests.
  - Each condition uses `seed_list=12345+(0:mc_count-1)`.
  - Channel, boundary, diagnostics, and tap settings match `monte_carlo_surface_channel_vertical.m`.
- Sea-state ensemble formulas:
  - `H_{c,m}(f)` is the channel response for sea-state condition `c` and seed sample `m`.
  - `mean_H_c(f) = (1/M) * sum_m H_{c,m}(f)`.
  - `var_H_c(f) = (1/M) * sum_m abs(H_{c,m}(f)-mean_H_c(f)).^2`.
  - Reference-frequency amplitude, phase, reflection/direct ratio, C2/C2.5 diagnostic, and tap statistics are computed independently inside each sea-state condition.
- Outputs:
  - `sweep_monte_carlo_surface_channel_vertical_result.mat` with `base_params`, `sea_hs_values`, `sea_wind_values`, `seed_list`, `condition_specs`, `condition_results`, `run_summary_table`, `condition_summary_table`, and `sweep_stats`.
  - `condition_results(cc).mc_stats.H_f.samples` stores the compact `[Nf x mc_count]` frequency-response sample matrix. Full channel structs, spatial fields, and 2-D spectra are not saved.
  - `run_summary_table` has one row per `(Hs, wind, seed)`.
  - `condition_summary_table` has one row per `(Hs, wind)`.
  - `sweep_stats` stores trend matrices indexed as `[numel(sea_hs_values) x numel(sea_wind_values)]`.
- PNG summaries:
  - heatmaps for `abs_H_ref_mean`, `abs_H_ref_std`, `reflect_rms_delta_k_mean`, `reflect_high_k_fraction_mean`, `rms_delta_k_increase_mean`, and `tap_rms_delay_symbols_mean`.
  - line plots of `reflect_high_k_fraction_mean` and `reflect_rms_delta_k_mean` versus wind speed, grouped by `Hs`.
- Scope:
  - This is still empirical Monte Carlo over the implemented Kirchhoff phase-screen channel. It is not a closed-form statistical channel model, a calibrated sea-surface scattering law, a BER/SER simulation, or a new rough-surface scattering solver.
  - Wind-speed trends under fixed `sea_hs_target` should be interpreted with the same caution as the C2.5 sweep: PM spectral shape changes are rescaled to the requested wave height.

### Communication Monte Carlo PSK statistics script
- `monte_carlo_comm_surface_psk_vertical.m` runs reduced QPSK BER/SER Monte Carlo statistics over the current vertical channel and communication chain.
- Default sea conditions:
  - weak: `sea_hs_target=0.05`, `sea_wind_speed=5`.
  - strong: `sea_hs_target=1.0`, `sea_wind_speed=5`.
- Default scenarios:
  - `direct_only`: `enable_surface_reflection=false`.
  - `direct_plus_reflect`: `enable_surface_reflection=true`, `surface_boundary_model='kirchhoff_kdomain'`.
- Default communication setup:
  - `M=4`, `n_sym=1000`, `EbN0_dB_list=0:2:20`, `symbol_rate_hz=1000`.
  - `mc_count=4`, overridable through `SURFACE_COMM_MC_COUNT`.
  - `seed_list=12345+(0:mc_count-1)`.
  - fixed bit stream with `bits_seed=9000`.
  - varying AWGN seed: `noise_seed_base + condition_index*100000 + scenario_index*10000 + seed_index*100 + ebn0_index`.
- Data flow:
  - each run calls `vertical_channel_model(paramsV)`;
  - builds baseband taps from `H_f` with the same `ifft(ifftshift(...))` logic as `comm_main_vertical_psk.m`;
  - applies the explicit receive-window policy used by `comm_main_vertical_psk.m`;
  - injects AWGN with `noise_inject_vertical`;
  - equalizes with frequency-domain MMSE using known effective taps;
  - demodulates with `modem_psk` and stores BER/SER.
- D2 communication policy:
  - `receive_window_mode='peak_sync'` by default.
  - `peak_sync` finds the largest tap, discards pre-peak taps, applies `conv(...,'full')`, and keeps the first `N` received samples.
  - `ebn0_reference='rx_clean'` by default, so AWGN power is referenced to the received clean waveform instead of transmit-symbol power.
  - `same_legacy` remains available only as a local diagnostic mode in scripts that define it.
- Statistics:
  - `BER_{c,s,m}(e)` and `SER_{c,s,m}(e)` are error rates for sea condition `c`, scenario `s`, seed sample `m`, and Eb/N0 bin `e`.
  - `mean_BER_{c,s}(e) = (1/M) * sum_m BER_{c,s,m}(e)`.
  - `var_BER_{c,s}(e) = (1/M) * sum_m (BER_{c,s,m}(e)-mean_BER_{c,s}(e)).^2`.
  - The same mean, variance, standard deviation, and 5/25/50/75/95 percent quantiles are reported for SER and effective SNR.
- Outputs:
  - `monte_carlo_comm_surface_psk_vertical_result.mat` with `params_base`, `comm_cfg`, `sea_conditions`, `scenarios`, `seed_list`, `EbN0_dB_list`, `run_summary_table`, `curve_summary_table`, `curve_stats`, and compact BER/SER sample arrays.
  - PNG summaries for weak/strong direct-plus-reflect BER/SER and direct-only versus direct-plus-reflect comparisons.
- Scope:
  - This is an empirical reduced communication experiment using the current channel and current PSK chain. It is not a new channel model, BER theory, stochastic generator, or rough-surface scattering solver.
  - Zero BER at high Eb/N0 is stored as zero when it occurs; plot-only floors are used only to make semilog figures readable.

### Minimal communication closed-loop validation script
- `validate_comm_link_minimal_vertical.m` validates the communication chain in layers before interpreting PE channel BER/SER.
- Validation order:
  - ideal unit channel: `h_bb=1`;
  - single complex tap: `h_bb=0.7*exp(1i*0.8)`;
  - known short multipath: `h_bb=[1;0.35*exp(1i*0.5);0.15*exp(-1i*0.7)]`;
  - reduced PE direct-only `h_bb` current path and diagnostic peak-aligned path.
- The script uses:
  - `M=4`, `n_sym=2000`, `EbN0_dB_list=0:2:20`;
  - fixed `bits_seed=9000`;
  - deterministic AWGN seeds from `noise_seed_base=7000`.
- It compares receive slicing modes:
  - `same`: legacy `conv(tx_symbols,h_bb,'same')` diagnostic behavior retained for pre-D2 comparisons.
  - `causal_head`: diagnostic-only `conv(...,'full')` followed by the first `N` samples, appropriate for causal taps whose main energy starts at tap 1.
  - `peak_sync`: D2 policy that shifts the effective tap origin to the dominant tap before causal-head receive-window selection.
- Outputs:
  - `validate_comm_link_minimal_vertical_result.mat` with `comm_cfg`, `test_cases`, `summary_table`, `case_results`, and `pe_alignment`.
  - PNG curves for simple-channel BER/SER and PE tap magnitudes before/after diagnostic peak alignment.
- Scope:
  - This script is diagnostic only. It does not change `comm_main_vertical_psk.m`, `monte_carlo_comm_surface_psk_vertical.m`, or PE propagation outputs.

### Runtime controls
- `save_mode`:
  - `rx_only`
  - `slice`
- `show_figures` toggles plotting and rough-surface visualization.
- `use_gpu` enables the GPU path only when a GPU is available.

## Current Output Structure

### Core field outputs
- `psifinal_xy`
- `psiout`
- `x`
- `y`
- `z_track`
- `Axz`
- `Ayz`

### Validation outputs
- `A_center`
- `R_center`
- `fit_slope`
- `fit_err_rms`
- `pass_1_over_R`
- `fit_mask`

### Reflection outputs
- `surface_elevation`
- `delta_phi`
- `psi_ref`
- `roughness_meta`

### Frequency-response outputs
- `f_axis`
- `H_direct_f`
- `H_reflect_f`
- `H_f`
- `idx_f_ref`

### Reference-frequency scalar outputs
- `h_direct`
- `h_reflect`
- `h_total`

### State and summary outputs
- `rx_state_used`
- `fd_hz_used`
- `rx_amplitude`
- `rx_phase_rad`
- `path_loss_db`
- `direct_to_reflect_db`
- `phase_diff_rad`
- `config`

## Known Structural Sensitivities

### Validation boundary
- `local_prepare_config` is the main contract boundary.
- New models often require new fields here and may also trigger existing constraints.

### Dual propagation logic
- There are two propagation implementations:
  - the main direct-path frequency loop
  - `local_march_field`
- Any physics change that affects marching should be applied consistently to both paths.

### Reference-frequency capture behavior
- Some outputs are only retained for `idx_f_ref`.
- This includes field slices and rough-surface products.
- Do not assume those stored arrays represent all frequencies.

### Communication coupling
- `comm_main_vertical_psk.m` depends on the meaning and consistency of:
  - `H_f`
  - `f_axis`
  - `idx_f_ref`
  - `h_total`
- Breaking these semantics will break the communication demo even if the channel demo still runs.

## Chronological Project Log / 按时间追加的项目日志

Historical entries are evidence of what was implemented or concluded at that
time. They are not silently rewritten when later work changes the public
semantics. Superseding entries explicitly identify the affected conclusion.

## 2026-06-10 Kirchhoff k-domain Boundary Validation

Changed files for this interface:
- `vertical_channel_model.m`: adds user-facing `paramsV` fields and validation.
- `vertical_wape_propagator.m`: forwards the surface boundary configuration and preserves disabled-path metadata.
- `pm_surface_boundary_model.m`: implements spatial and implicit k-domain boundary application paths.
- `PROJECT_CONTEXT.md`: records formulas, interface, assumptions, limits, and validation.

Reduced-grid validation settings:
- Common channel setup: `f0=4000`, `c0=1500`, `z_tx=100`, `z_rx=3`, `xw=yw=50`, `nx=ny=128`, `sigma_src_m=0.4`, `stepz_lamb=0.5`, `sponge_ratio=0.12`, `alpha_max_np_per_m=0.15`, `env_mode='uniform'`, `show_figures=false`, `save_mode='rx_only'`, `use_gpu=false`, `sea_seed=12345`, `enforce_1_over_R=false`.
- Single-frequency cases compared default, explicit `kirchhoff_spatial`, explicit `kirchhoff_kdomain`, and direct-only.
- Wideband case used `enable_wideband=true`, `f_band_hz=[4000 8000]`, `Nf_min=Nf_max=16`, `f_ref_hz=6000`.

Validation results:
- Default vs explicit spatial `H_f` relative error: `0`.
- Single-frequency spatial vs kdomain:
  - `psi_ref` relative error: `4.753077448146337e-16`
  - `H_reflect_f` relative error: `6.785387727900559e-15`
  - `H_f` relative error: `5.456834559265761e-16`
  - `roughness_meta.boundary_equivalence_error.rel_l2`: `4.753077448146337e-16`
  - `roughness_meta.boundary_equivalence_error.max_abs`: `2.5033568429418485e-17`
  - equivalence check passed for tolerance `1e-10`.
- Direct-only regression:
  - `max(abs(H_reflect_f))`: `0`
  - `max(abs(H_f-(H_direct_f+H_reflect_f)))`: `0`
  - `roughness_meta.enabled=false`
- Wideband spatial vs kdomain:
  - spatial invariant `max(abs(H_f-(H_direct_f+H_reflect_f)))`: `0`
  - kdomain invariant `max(abs(H_f-(H_direct_f+H_reflect_f)))`: `0`
  - `H_f` relative error: `7.666749014695518e-16`
  - reference-frequency kdomain equivalence `rel_l2`: `4.827102620732236e-16`
  - reference-frequency kdomain equivalence `max_abs`: `2.964296775167217e-17`

Remaining issues:
- No dense boundary matrix is materialized, so matrix-level inspection of individual `B_xi(K,K')` entries is not available.
- The k-domain path intentionally remains mathematically equivalent to the old Kirchhoff phase-screen model and should not be interpreted as a higher-order rough-surface scattering model.

## 2026-06-10 Boundary Coupling Diagnostic Validation

Changed files for this diagnostic:
- `vertical_channel_model.m`: adds default-off diagnostic flags.
- `vertical_wape_propagator.m`: forwards diagnostic flags only for the reference-frequency reflection call and preserves disabled metadata.
- `pm_surface_boundary_model.m`: computes scalar coupling metrics from `G_xy` and optional compact debug summaries.
- `PROJECT_CONTEXT.md`: records formulas, interpretation, limits, and validation.

Reduced-grid validation settings:
- Common channel setup: `f0=4000`, `c0=1500`, `z_tx=100`, `z_rx=3`, `xw=yw=50`, `nx=ny=128`, `sigma_src_m=0.4`, `stepz_lamb=0.5`, `sponge_ratio=0.12`, `alpha_max_np_per_m=0.15`, `env_mode='uniform'`, `show_figures=false`, `save_mode='rx_only'`, `use_gpu=false`, `sea_seed=12345`, `surface_boundary_model='kirchhoff_kdomain'`, `enforce_1_over_R=false`.
- Diagnostic-off vs diagnostic-on used `sea_hs_target=0.5`.
- Weak roughness used `sea_hs_target=0.05`; strong roughness used `sea_hs_target=1.0`.
- Wideband case used `enable_wideband=true`, `f_band_hz=[4000 8000]`, `Nf_min=Nf_max=16`, `f_ref_hz=6000`.

Validation results:
- Diagnostic off/on consistency:
  - off `boundary_coupling_diagnostics.enabled=false`
  - on `boundary_coupling_diagnostics.enabled=true`
  - `H_f` relative error: `0`
  - `H_reflect_f` relative error: `0`
  - `psi_ref` relative error: `0`
- Weak vs strong roughness diagnostics:
  - weak `nonzero_power_fraction`: `0.50513062402654019`
  - strong `nonzero_power_fraction`: `0.99999482166633924`
  - weak `rms_delta_k_rad_per_m`: `0.71355283273946624`
  - strong `rms_delta_k_rad_per_m`: `6.5608130546379275`
  - weak `energy_radius_90_rad_per_m`: `0.95702627363155124`
  - strong `energy_radius_90_rad_per_m`: `8.9856140404778042`
  - weak `nonzero_energy_radius_90_rad_per_m`: `1.3119676091847634`
  - strong `nonzero_energy_radius_90_rad_per_m`: `8.9856140404778042`
- Direct-only regression:
  - `max(abs(H_reflect_f))`: `0`
  - `max(abs(H_f-(H_direct_f+H_reflect_f)))`: `0`
  - `boundary_coupling_diagnostics.enabled=false`
- Wideband regression:
  - `max(abs(H_f-(H_direct_f+H_reflect_f)))`: `0`
  - diagnostic metadata enabled at reference frequency.
  - `idx_f_ref=8`, `f_axis(idx_f_ref)=5866.666666666667 Hz`
- Debug metadata check:
  - radial bin count: `16`
  - top nonzero peak count: `8`
  - full `G_hat_k` stored: `false`
  - full `P_k` stored: `false`

Remaining issues:
- The metrics summarize the Fourier spectrum of the imposed Kirchhoff screen. They do not identify individual physically rigorous scattering mechanisms.
- Strong roughness can push the phase-screen approximation outside a conservative small-perturbation interpretation; treat the diagnostic trend as a numerical screen-spectrum trend, not as validated ocean-surface scattering physics.

## 2026-06-10 Incident-Weighted Redistribution Diagnostic Validation

Changed files for this diagnostic:
- `vertical_channel_model.m`: adds default-off redistribution diagnostic flags.
- `vertical_wape_propagator.m`: forwards redistribution flags only for the reference-frequency reflection call and preserves disabled metadata.
- `pm_surface_boundary_model.m`: compares `Psi_inc_k`, `Psi_ref_k`, and flat-reflected spectra without changing propagation results.
- `sweep_surface_boundary_redistribution_vertical.m`: reduced scalar trend sweep for `Hs`, wind speed, and frequency.
- `PROJECT_CONTEXT.md`: records formulas, interpretation, limits, and validation.

Reduced-grid validation settings:
- Common channel setup: `c0=1500`, `z_tx=100`, `z_rx=3`, `xw=yw=50`, `nx=ny=128`, `sigma_src_m=0.4`, `stepz_lamb=0.5`, `sponge_ratio=0.12`, `alpha_max_np_per_m=0.15`, `env_mode='uniform'`, `show_figures=false`, `save_mode='rx_only'`, `use_gpu=false`, `sea_seed=12345`, `surface_boundary_model='kirchhoff_kdomain'`, `enforce_1_over_R=false`.

Validation results:
- Redistribution off/on consistency at `f0=4000`, `sea_hs_target=0.5`, `sea_wind_speed=5`:
  - off `boundary_redistribution_diagnostics.enabled=false`
  - on `boundary_redistribution_diagnostics.enabled=true`
  - `H_f` relative error: `0`
  - `H_reflect_f` relative error: `0`
  - `psi_ref` relative error: `0`
- Weak vs strong roughness at `f0=4000`, `sea_wind_speed=5`:
  - weak `sea_hs_target=0.05`, `reflect_rms_delta_k=2.4092056625755345`
  - strong `sea_hs_target=1.0`, `reflect_rms_delta_k=6.552320692768185`
  - weak `rms_delta_k_increase=0.12555034900919093`
  - strong `rms_delta_k_increase=4.268665379201842`
  - weak `reflect_high_k_fraction=0.11866935742346503`
  - strong `reflect_high_k_fraction=0.8570406763058447`
  - weak `rough_vs_flat_rms_delta_k_increase=0.12555034900919093`
  - strong `rough_vs_flat_rms_delta_k_increase=4.268665379201842`
- Frequency trend at `sea_hs_target=0.5`, `sea_wind_speed=5`:
  - `f0=[4000,6000,8000] Hz`
  - `rms_delta_k_increase=[3.6659013500300062,4.0216592688428765,4.0046896856288567]`
  - `reflect_high_k_fraction=[0.7832240337445295,0.8228809206147447,0.8185641917866873]`
- Direct-only regression:
  - `max(abs(H_reflect_f))`: `0`
  - `max(abs(H_f-(H_direct_f+H_reflect_f)))`: `0`
  - `boundary_redistribution_diagnostics.enabled=false`
- Wideband regression:
  - `max(abs(H_f-(H_direct_f+H_reflect_f)))`: `0`
  - redistribution metadata enabled at reference frequency.
  - `idx_f_ref=8`, `f_axis(idx_f_ref)=5866.666666666667 Hz`
- Debug metadata check:
  - incident radial bin count: `16`
  - reflected top peak count: `8`
  - full `Psi_inc_k` stored: `false`
  - full `Psi_ref_k` stored: `false`
- Trend script `sweep_surface_boundary_redistribution_vertical` completed successfully. Representative rows:
  - Hs sweep broadening: `0.12555`, `3.6659`, `4.2687` for `Hs=0.05`, `0.5`, `1.0`
  - Wind sweep broadening at fixed `Hs=0.5`: `4.2995`, `3.6659`, `1.1948`, `0.67448` for wind `3`, `5`, `8`, `12 m/s`
  - Frequency sweep broadening: `3.6659`, `4.0217`, `4.0047` for `4000`, `6000`, `8000 Hz`

Remaining issues:
- The wind-speed trend is not monotonic under fixed `Hs_target` scaling; it reflects PM spectral-shape changes after amplitude normalization.
- The diagnostics characterize the simulated incident field and Kirchhoff screen only. They should not be used as calibrated ocean-surface scattering cross sections or strict non-diagonal boundary-operator matrix entries.

## 2026-06-10 Surface Monte Carlo Channel Statistics Validation

Changed files for this script:
- `monte_carlo_surface_channel_vertical.m`: independent reduced-grid Monte Carlo script over `sea_seed`.
- `PROJECT_CONTEXT.md`: records Monte Carlo scope, formulas, outputs, validation settings, and limits.

Reduced-grid Monte Carlo settings:
- Fixed channel setup: `c0=1500`, `z_tx=100`, `z_rx=3`, `xw=yw=50`, `nx=ny=128`, `sigma_src_m=0.4`, `stepz_lamb=0.5`, `sponge_ratio=0.12`, `alpha_max_np_per_m=0.15`, `env_mode='uniform'`, `show_figures=false`, `save_mode='rx_only'`, `use_gpu=false`, `surface_boundary_model='kirchhoff_kdomain'`, `enforce_1_over_R=false`.
- Wideband setup: `enable_wideband=true`, `f_band_hz=[4000 8000]`, `Nf_min=Nf_max=16`, `f_ref_hz=6000`.
- Sea-state setup: `sea_hs_target=0.5`, `sea_wind_speed=5`, `surface_boundary_coupling_diagnostics=true`, `surface_boundary_redistribution_diagnostics=true`.
- Tap setup: `symbol_rate_hz=1000`, `tap_fft_len=512`, `tap_energy_ratio=0.999`.
- Smoke test: `SURFACE_MC_COUNT=4`, `seed_list=12345:12348`.
- Default reduced run: `mc_count=16`, `seed_list=12345:12360`.

Validation results:
- Smoke test with 4 seeds completed successfully and produced finite `H_f`, diagnostic, and tap statistics.
- Default 16-seed run completed successfully.
- All runs shared the same frequency axis and `idx_f_ref`.
- Reference frequency:
  - `idx_f_ref=8`
  - `f_axis(idx_f_ref)=5866.666666666667 Hz`
- Channel invariant:
  - `max(abs(H_f-(H_direct_f+H_reflect_f)))=1.3985787785349405e-17`
- Direct-path consistency across seeds:
  - `direct_drift_max_abs=1.3904866443919908e-17`
- Reference-frequency total channel statistics:
  - `mean(abs(H(f_ref)))=0.04451535630826136`
  - `std(abs(H(f_ref)))=0.019683054200486178`
  - `abs(H(f_ref))` quantiles `[5,25,50,75,95]% = [0.012448453947227933, 0.029477650012892355, 0.047715895428688343, 0.05917630294551255, 0.068142386535915839]`
  - circular mean phase `=-1.5721223443847119 rad`
  - circular phase variance `=0.42130134481217185`
- Reflection/direct ratio statistics at `f_ref`:
  - `mean(abs(h_reflect/h_direct))=0.83366541472645794`
  - `std(abs(h_reflect/h_direct))=0.357222192582143`
  - circular mean phase `=-2.447845029934725 rad`
  - circular phase variance `=0.832299121295223`
- Coupling diagnostic statistics:
  - `nonzero_power_fraction` mean/std `=0.99990160177794407 / 6.5681746413696337e-05`
  - `rms_delta_k_rad_per_m` mean/std `=6.4634066111591792 / 0.052390699662500283`
- Redistribution diagnostic statistics:
  - `reflect_rms_delta_k_rad_per_m` mean/std `=6.4729622076609932 / 0.0676218962238834`
  - `reflect_high_k_fraction` mean/std `=0.81498689692969417 / 0.0092872966143607676`
  - `rms_delta_k_increase_rad_per_m` mean/std `=3.9946260540001335 / 0.0676218962238834`
- Tap metric statistics:
  - `tap_count` mean/std `=512 / 0`
  - `tap_energy_kept` mean `=1`
  - `tap_rms_delay_symbols` mean/std `=99.835328110577819 / 69.550813743973464`
  - `tap_peak_fraction` mean/std `=0.87083037958185672 / 0.18399224840340037`
- Generated result artifacts:
  - `monte_carlo_surface_channel_vertical_result.mat`
  - `monte_carlo_surface_channel_vertical_H_f_mean_magnitude.png`
  - `monte_carlo_surface_channel_vertical_abs_h_ref_hist.png`
  - `monte_carlo_surface_channel_vertical_phase_h_ref_hist.png`
  - `monte_carlo_surface_channel_vertical_reflect_rms_delta_k_hist.png`
  - `monte_carlo_surface_channel_vertical_reflect_high_k_fraction_hist.png`
  - `monte_carlo_surface_channel_vertical_tap_rms_delay_hist.png`

Remaining issues:
- `mc_count=16` is a reduced validation ensemble, not a converged ocean-channel distribution.
- The only Monte Carlo variable is `sea_seed`; sea-state, geometry, frequency grid, and solver settings are fixed.
- The script reports tap metrics only. It does not run BER/SER or modify `comm_main_vertical_psk.m`.
- The statistics are empirical summaries of the existing Kirchhoff phase-screen channel and should not be treated as a closed-form statistical channel generator.

## 2026-06-10 Multi-Sea-State Monte Carlo Sweep Validation

Changed files for this script:
- `sweep_monte_carlo_surface_channel_vertical.m`: independent reduced-grid `Hs x wind x seed` Monte Carlo sweep.
- `PROJECT_CONTEXT.md`: records C3.5 sweep scope, formulas, outputs, validation settings, and limits.

Reduced-grid sweep settings:
- Fixed channel setup: `c0=1500`, `z_tx=100`, `z_rx=3`, `xw=yw=50`, `nx=ny=128`, `sigma_src_m=0.4`, `stepz_lamb=0.5`, `sponge_ratio=0.12`, `alpha_max_np_per_m=0.15`, `env_mode='uniform'`, `show_figures=false`, `save_mode='rx_only'`, `use_gpu=false`, `surface_boundary_model='kirchhoff_kdomain'`, `enforce_1_over_R=false`.
- Wideband setup: `enable_wideband=true`, `f_band_hz=[4000 8000]`, `Nf_min=Nf_max=16`, `f_ref_hz=6000`.
- Diagnostics: `surface_boundary_coupling_diagnostics=true`, `surface_boundary_redistribution_diagnostics=true`.
- Sweep grid: `sea_hs_target=[0.05, 0.5, 1.0]`, `sea_wind_speed=[3, 5, 8, 12]`.
- Tap setup: `symbol_rate_hz=1000`, `tap_fft_len=512`, `tap_energy_ratio=0.999`.
- Smoke test: `SURFACE_MC_SWEEP_COUNT=2`, `SURFACE_MC_SWEEP_MAX_CONDITIONS=2`.
- Default reduced run after the C4 follow-up default-count increase: `12` sea-state conditions, `mc_count=16`, `192` total propagations, `seed_list=12345:12360`.
- The previously recorded validation matrices below came from the earlier cost-control run with `mc_count=4`, `48` total propagations, and `seed_list=12345:12348`.

Validation results:
- Smoke test completed successfully and produced finite `.mat`, table, and PNG outputs.
- Default 48-run sweep completed successfully.
- All condition-level frequency axes and `idx_f_ref` checks passed.
- `sweep_stats` matrix dimensions are `[3 x 4]`, matching `[numel(sea_hs_values) x numel(sea_wind_values)]`.
- Reference frequency:
  - `idx_f_ref=8`
  - `f_axis(idx_f_ref)=5866.666666666667 Hz`
- Channel invariant:
  - `max(abs(H_f-(H_direct_f+H_reflect_f)))=1.3877787807814457e-17`
- Direct-path consistency:
  - maximum condition-level `direct_drift_max_abs=0`
- Summary table finite-field check:
  - all selected channel, phase, ratio, diagnostic, and tap statistics are finite.
- `abs_H_ref_mean` matrix, with rows `Hs=[0.05,0.5,1.0]` and columns `wind=[3,5,8,12]`:
  - `[0.039000343834070653, 0.03775446157444675, 0.021052540344231319, 0.01646692716410688]`
  - `[0.027374350743730395, 0.038604342599260662, 0.055552395529549083, 0.031575982180233053]`
  - `[0.03928118735487425, 0.047635492235929437, 0.058271523033646792, 0.068788019597050512]`
- `reflect_rms_delta_k_mean` matrix:
  - `[3.3081145971986365, 2.7127732766653132, 2.5488843668815639, 2.5248328653472694]`
  - `[6.55751077252848, 6.5011848440695914, 5.027115991544723, 4.4226857324224929]`
  - `[6.5588409967053547, 6.5902779083995151, 6.3846289902553384, 6.039627688146652]`
- `reflect_high_k_fraction_mean` matrix:
  - `[0.24294922660924467, 0.13906324775944112, 0.11172646818271345, 0.10766906380516775]`
  - `[0.82530357342614824, 0.81912102056925717, 0.59471503928237834, 0.48516109208527497]`
  - `[0.82625304890765827, 0.83032699700837842, 0.80360740036150746, 0.75151018098143552]`
- `rms_delta_k_increase_mean` matrix:
  - `[0.82977844353777486, 0.23443712300445163, 0.0705482132207026, 0.046496711686407788]`
  - `[4.0791746188676186, 4.0228486904087308, 2.5487798378838615, 1.9443495787616316]`
  - `[4.0805048430444932, 4.1119417547386528, 3.9062928365944769, 3.5612915344857905]`
- `tap_rms_delay_symbols_mean` matrix:
  - `[37.471572903757632, 35.82348841242306, 35.621289665847, 41.278879008243244]`
  - `[201.14072827568845, 109.35908886719794, 47.846354301354722, 72.62237807072276]`
  - `[134.45187275648524, 103.99139368495908, 95.39888154646296, 78.994400218459091]`
- Generated result artifacts:
  - `sweep_monte_carlo_surface_channel_vertical_result.mat`
  - `sweep_monte_carlo_surface_channel_vertical_abs_H_ref_mean_heatmap.png`
  - `sweep_monte_carlo_surface_channel_vertical_abs_H_ref_std_heatmap.png`
  - `sweep_monte_carlo_surface_channel_vertical_reflect_rms_delta_k_mean_heatmap.png`
  - `sweep_monte_carlo_surface_channel_vertical_reflect_high_k_fraction_mean_heatmap.png`
  - `sweep_monte_carlo_surface_channel_vertical_rms_delta_k_increase_mean_heatmap.png`
  - `sweep_monte_carlo_surface_channel_vertical_tap_rms_delay_mean_heatmap.png`
  - `sweep_monte_carlo_surface_channel_vertical_reflect_high_k_fraction_vs_wind.png`
  - `sweep_monte_carlo_surface_channel_vertical_reflect_rms_delta_k_vs_wind.png`

Remaining issues:
- `mc_count=16` is still a reduced sweep default for cost control, not a converged ensemble estimate.
- The sweep varies only `sea_hs_target`, `sea_wind_speed`, and `sea_seed`; geometry, frequency band, grid, and solver settings stay fixed.
- Wind-speed trends are conditional on PM spectral shape after `sea_hs_target` rescaling and should not be read as an independent wave-height trend.
- No BER/SER or communication-chain Monte Carlo is run in this script.
- The saved statistics remain empirical summaries of the current Kirchhoff phase-screen implementation, not a closed-form statistical channel generator.

## 2026-06-10 C4 Empirical Random Channel Generator Prototype

Changed files for this prototype:
- `build_surface_empirical_channel_model_vertical.m`: reads C3/C3.5 Monte Carlo result files and builds a compact empirical sample model.
- `sample_surface_empirical_channel_vertical.m`: draws nearest-neighbor empirical bootstrap samples from the compact model.
- `validate_surface_empirical_channel_generator_vertical.m`: validates C4 without running PE/WAPE propagation.
- `PROJECT_CONTEXT.md`: records C4 purpose, interface, validation settings, results, and limits.

Purpose and scope:
- C4 converts existing Monte Carlo summary rows into a lightweight empirical random generator for reference-frequency and low-dimensional channel summaries.
- It does not modify `vertical_channel_model`, `vertical_wape_propagator`, `pm_surface_boundary_model`, or communication-chain public outputs.
- It does not run PE/WAPE and does not save full channel structs, spatial fields, or two-dimensional spectra.
- It is an empirical resampler of previously computed C3/C3.5 results, not a new rough-surface scattering theory.

Input and model construction:
- Default input priority:
  - `sweep_monte_carlo_surface_channel_vertical_result.mat`
  - fallback: `monte_carlo_surface_channel_vertical_result.mat`
- The builder reads `run_summary_table` for C3.5 or `summary_table` for C3 and groups rows by sea-state condition.
- Each condition stores `sea_hs_target`, `sea_wind_speed`, `f_ref_hz`, `mc_count`, seed list, compact source samples, and per-metric empirical statistics.
- C4 v1 uses nearest-neighbor matching in normalized `(Hs, wind, f_ref)` space. It does not interpolate between sea states.

Generated sample variables:
- Reference-frequency response:
  - `abs_h_ref`
  - `phase_h_ref_rad`
  - reconstructed `h_ref_complex = abs_h_ref * exp(1i*phase_h_ref_rad)`
- Reflection/direct ratio:
  - `reflect_direct_abs`
  - `reflect_direct_phase_rad`
  - reconstructed `reflect_direct_ratio_complex`
- C2/C2.5 diagnostics:
  - `coupling_nonzero_power_fraction`
  - `coupling_rms_delta_k_rad_per_m`
  - `redistribution_reflect_rms_delta_k_rad_per_m`
  - `redistribution_reflect_high_k_fraction`
  - `redistribution_rms_delta_k_increase_rad_per_m`
- Tap summaries:
  - `tap_count`
  - `tap_rms_delay_symbols`
  - `tap_peak_fraction`

Sampling rule:
- For query sea state `q=(Hs_q, wind_q, f_ref_q)`, select the stored condition
  `c* = argmin_c ||[(Hs_c-Hs_q)/range(Hs), (wind_c-wind_q)/range(wind), (f_ref_c-f_ref_q)/range(f_ref)]||_2`.
- Draw `n_samples` indices with replacement from that condition's Monte Carlo seed rows.
- Copy the corresponding low-dimensional samples into the generated table.
- Metadata records `source_file`, `source_type`, `selection_mode='nearest'`, `sampling_method='empirical_bootstrap_with_replacement'`, matched condition, distance, RNG seed, and sampled source-row indices.

Smoke validation settings:
- Command environment:
  - `SURFACE_EMPIRICAL_SAMPLE_COUNT=200`
  - `SURFACE_EMPIRICAL_MAX_CONDITIONS=2`
- Input file:
  - `sweep_monte_carlo_surface_channel_vertical_result.mat`
- Selected conditions:
  - condition 1: `Hs=0.05`, `wind=3`, `f_ref_hz=5866.666666666667`, `mc_count=4`
  - condition 12: `Hs=1.0`, `wind=12`, `f_ref_hz=5866.666666666667`, `mc_count=4`
- Output file:
  - `surface_empirical_channel_generator_validation_result.mat`
  - size `37379` bytes
- Static check:
  - C4 files do not call `vertical_channel_model`, `vertical_wape_propagator`, or `pm_surface_boundary_model`.

Smoke validation results:
- Condition 1, `Hs=0.05`, `wind=3`:
  - `abs_h_ref` source/generated mean `=0.039000 / 0.038729`, std `=0.032162 / 0.027823`
  - `phase_h_ref_rad` source/generated mean `=-0.85219 / -0.86150`, std `=1.8363 / 1.6075`
  - `reflect_direct_abs` source/generated mean `=1.0942 / 1.0875`, std `=0.35971 / 0.31444`
  - `redistribution_reflect_rms_delta_k_rad_per_m` source/generated mean `=3.3081 / 3.3086`, std `=0.035956 / 0.031058`
  - `redistribution_reflect_high_k_fraction` source/generated mean `=0.24295 / 0.24323`, std `=0.0076396 / 0.0065539`
  - `tap_rms_delay_symbols` source/generated mean `=37.472 / 37.576`, std `=13.653 / 11.866`
- Condition 12, `Hs=1.0`, `wind=12`:
  - `abs_h_ref` source/generated mean `=0.068788 / 0.067981`, std `=0.017683 / 0.014949`
  - `phase_h_ref_rad` source/generated mean `=-1.5262 / -1.4930`, std `=0.66951 / 0.58903`
  - `reflect_direct_abs` source/generated mean `=1.0128 / 1.0014`, std `=0.29010 / 0.24675`
  - `redistribution_reflect_rms_delta_k_rad_per_m` source/generated mean `=6.0396 / 6.0251`, std `=0.40954 / 0.34363`
  - `redistribution_reflect_high_k_fraction` source/generated mean `=0.75151 / 0.74957`, std `=0.055883 / 0.046882`
  - `tap_rms_delay_symbols` source/generated mean `=78.994 / 80.297`, std `=49.800 / 44.629`

Generated C4 validation artifacts:
- `surface_empirical_channel_generator_validation_result.mat`
- `surface_empirical_channel_generator_abs_h_ref_hist.png`
- `surface_empirical_channel_generator_phase_h_ref_rad_hist.png`
- `surface_empirical_channel_generator_reflect_direct_abs_hist.png`
- `surface_empirical_channel_generator_redistribution_reflect_rms_delta_k_rad_per_m_hist.png`
- `surface_empirical_channel_generator_redistribution_reflect_high_k_fraction_hist.png`
- `surface_empirical_channel_generator_tap_rms_delay_symbols_hist.png`

Remaining C4 issues:
- The saved C3.5 source used by this C4 smoke test has only `mc_count=4` per condition, so bootstrap samples replicate a very small empirical support. Rerunning the updated C3.5 script with its new `mc_count=16` default increases the empirical support to 16 seeds per condition.
- Nearest-neighbor matching should not be used to extrapolate outside the scanned `Hs x wind` range.
- C4 v1 does not generate a full frequency response `H_f(f)` and does not replace PE/WAPE propagation.
- Phase samples are bootstrapped from principal phases; no circular distribution is fitted.
- C4 is an empirical reduced generator prototype, not a T matrix, SSA/NLSSA, scattering cross-section model, or closed-form stochastic channel model.

## 2026-06-10 C3.5 Monte Carlo Count Increase

Changed files for this update:
- `sweep_monte_carlo_surface_channel_vertical.m`: increases the default multi-sea-state Monte Carlo count and adds smoke-output redirection controls.
- `PROJECT_CONTEXT.md`: records the new default and smoke validation.

Updated defaults and controls:
- C3.5 default `mc_count` is now `16`.
- Default full reduced sweep size is now `3 Hs x 4 wind x 16 seed = 192` propagations.
- Default seed list is `12345:12360`.
- `SURFACE_MC_SWEEP_COUNT` still overrides the seed count.
- `SURFACE_MC_SWEEP_MAX_CONDITIONS` still truncates conditions for smoke/cost control.
- New optional smoke controls:
  - `SURFACE_MC_SWEEP_RESULT_FILE`
  - `SURFACE_MC_SWEEP_FIGURE_PREFIX`

Smoke validation:
- Command environment:
  - `SURFACE_MC_SWEEP_COUNT=2`
  - `SURFACE_MC_SWEEP_MAX_CONDITIONS=1`
  - `SURFACE_MC_SWEEP_RESULT_FILE=sweep_monte_carlo_surface_channel_vertical_smoke.mat`
  - `SURFACE_MC_SWEEP_FIGURE_PREFIX=sweep_monte_carlo_surface_channel_vertical_smoke_`
- Smoke output:
  - `n_conditions=1`
  - `run_summary_table` rows `=2`
  - `seed_count=2`
  - `condition_summary_table.mc_count=2`
  - `max_invariant_error=6.938893903907228e-18`
  - `direct_drift_max_abs=0`
- The smoke result used a redirected file and PNG prefix, so it did not overwrite the existing C3.5 source file used by C4.

Generated smoke artifacts:
- `sweep_monte_carlo_surface_channel_vertical_smoke.mat`
- `sweep_monte_carlo_surface_channel_vertical_smoke_abs_H_ref_mean_heatmap.png`
- `sweep_monte_carlo_surface_channel_vertical_smoke_abs_H_ref_std_heatmap.png`
- `sweep_monte_carlo_surface_channel_vertical_smoke_reflect_rms_delta_k_mean_heatmap.png`
- `sweep_monte_carlo_surface_channel_vertical_smoke_reflect_high_k_fraction_mean_heatmap.png`
- `sweep_monte_carlo_surface_channel_vertical_smoke_rms_delta_k_increase_mean_heatmap.png`
- `sweep_monte_carlo_surface_channel_vertical_smoke_tap_rms_delay_mean_heatmap.png`
- `sweep_monte_carlo_surface_channel_vertical_smoke_reflect_high_k_fraction_vs_wind.png`
- `sweep_monte_carlo_surface_channel_vertical_smoke_reflect_rms_delta_k_vs_wind.png`

Remaining issue:
- The full 192-propagation C3.5 run was not executed in this update; only the redirected smoke test was run. Existing `sweep_monte_carlo_surface_channel_vertical_result.mat` still reflects the previously saved `mc_count=4` run until the full script is rerun.

## 2026-06-11 Realization-Based Platform and C4 Wideband Bootstrap Update

Changed files for this update:
- `plot_c35_core_heatmaps_vertical.m`: reads the C3.5 result file and regenerates the five core heatmaps.
- `monte_carlo_comm_surface_psk_representative_vertical.m`: runs communication Monte Carlo on weak/mid/strong representative sea states while keeping the D2 link fixed.
- `build_surface_empirical_channel_model_vertical.m`: detects stored wideband `H_f` Monte Carlo samples in C3/C3.5 result files.
- `sample_surface_empirical_channel_vertical.m`: adds `sample_mode='wideband_hf'` and `sample_mode='tap_level'`.
- `validate_surface_empirical_channel_generator_vertical.m`: validates summary, wideband, and tap-level C4 modes.
- `reports/realization_to_empirical_channel_generator_report.md`: records the stage summary and main results.

Fixed communication policy:
- D2 is treated as fixed for this phase.
- `receive_window_mode='peak_sync'`.
- `ebn0_reference='rx_clean'`.
- No changes were made to `comm_main_vertical_psk.m`.

C3.5 full rerun:
- `sea_hs_target=[0.05,0.5,1.0]`.
- `sea_wind_speed=[3,5,8,12]`.
- `mc_count=16`, `seed_list=12345:12360`.
- Total propagation count: `12 x 16 = 192`.
- `run_summary_table` rows: `192`.
- `condition_summary_table` rows: `12`.
- Maximum `H_f-(H_direct_f+H_reflect_f)` invariant error: `1.551583845779546e-17`.
- Maximum direct-path drift: `1.390486644391991e-17`.
- All 12 conditions contain stored wideband `H_f.samples` with `Nf=16`.

C3.5 core heatmaps generated:
- `c35_core_abs_H_ref_mean_heatmap.png`
- `c35_core_abs_H_ref_std_heatmap.png`
- `c35_core_reflect_rms_delta_k_mean_heatmap.png`
- `c35_core_reflect_high_k_fraction_mean_heatmap.png`
- `c35_core_tap_rms_delay_mean_heatmap.png`

C3.5 selected condition results:
- weak sea, `Hs=0.05`, `wind=5`:
  - `abs_H_ref_mean=0.029723`
  - `abs_H_ref_std=0.019989`
  - `reflect_rms_delta_k_mean=2.6774`
  - `reflect_high_k_fraction_mean=0.13695`
  - `tap_rms_delay_symbols_mean=31.644`
- mid representative sea, `Hs=0.5`, `wind=8`:
  - `abs_H_ref_mean=0.045544`
  - `abs_H_ref_std=0.018677`
  - `reflect_rms_delta_k_mean=4.9028`
  - `reflect_high_k_fraction_mean=0.56899`
  - `tap_rms_delay_symbols_mean=97.101`
- strong representative sea, `Hs=1.0`, `wind=12`:
  - `abs_H_ref_mean=0.054750`
  - `abs_H_ref_std=0.021414`
  - `reflect_rms_delta_k_mean=6.0107`
  - `reflect_high_k_fraction_mean=0.74638`
  - `tap_rms_delay_symbols_mean=109.35`

Representative communication Monte Carlo:
- Sea states:
  - weak: `Hs=0.05`, `wind=5`
  - mid: `Hs=0.5`, `wind=8`
  - strong: `Hs=1.0`, `wind=12`
- Scenarios: `direct_only`, `direct_plus_reflect`.
- `mc_count=16`.
- Total channel runs: `3 x 2 x 16 = 96`.
- Maximum channel invariant error: `1.387778780781446e-17`.
- Maximum direct-only reflected response: `0`.
- Selected BER results:
  - weak direct-plus-reflect: `BER_mean` at `[0,10,20] dB = [0.11606, 3.125e-05, 0]`
  - mid direct-plus-reflect: `BER_mean` at `[0,10,20] dB = [0.17438, 0.045906, 0.037469]`
  - strong direct-plus-reflect: `BER_mean` at `[0,10,20] dB = [0.15978, 0.018125, 0.00875]`
- Mid and strong reflected cases retain residual high-Eb/N0 BER for some seeds under the current reduced-grid channel and MMSE receiver.

C4 upgraded generator:
- `summary` mode retains low-dimensional bootstrap behavior.
- `wideband_hf` mode bootstraps stored `H_f(f)` samples from the matched C3.5 condition.
- `tap_level` mode derives baseband/tap samples from bootstrapped `H_f(f)` without rerunning PE/WAPE.
- C4 validation using the new C3.5 source:
  - comparison rows: `18`
  - output file size: `219191` bytes
  - `H_f_sample_size=[16 50]`
  - `idx_f_ref=8`
  - `tap_count_mean=501.92`
  - `tap_rms_delay_mean=40.0435`
  - `tap_peak_fraction_mean=0.9857`
  - tap metrics finite: true
  - no PE/WAPE propagation run: true

Stage conclusion:
- The current platform is a realization-based PE/WAPE simulation and empirical Monte Carlo analysis pipeline.
- The next stage should focus on empirical statistics-based channel generation from the C3.5 sample ensemble.
- C4 remains an empirical bootstrap generator. It is not a T-matrix, SSA/NLSSA, physical scattering cross-section model, or closed-form statistical channel model.

## 2026-06-10 Communication Monte Carlo PSK Validation

Changed files for this script:
- `monte_carlo_comm_surface_psk_vertical.m`: independent reduced-grid QPSK BER/SER Monte Carlo script.
- `PROJECT_CONTEXT.md`: records communication Monte Carlo scope, formulas, random seed policy, outputs, validation settings, and limits.

Reduced communication Monte Carlo settings:
- Fixed channel setup: `c0=1500`, `z_tx=100`, `z_rx=3`, `xw=yw=50`, `nx=ny=128`, `sigma_src_m=0.4`, `stepz_lamb=0.5`, `sponge_ratio=0.12`, `alpha_max_np_per_m=0.15`, `env_mode='uniform'`, `show_figures=false`, `save_mode='rx_only'`, `use_gpu=false`, `surface_boundary_model='kirchhoff_kdomain'`, `enforce_1_over_R=false`.
- Wideband setup: `enable_wideband=true`, `f_band_hz=[4000 8000]`, `Nf_min=Nf_max=16`, `f_ref_hz=6000`.
- Sea conditions: weak `sea_hs_target=0.05`, strong `sea_hs_target=1.0`; both use `sea_wind_speed=5`.
- Scenarios: `direct_only` and `direct_plus_reflect`.
- Communication setup: `M=4`, `n_sym=1000`, `EbN0_dB_list=0:2:20`, `symbol_rate_hz=1000`.
- Seed policy:
  - `bits_seed=9000`, fixed for all runs.
  - `seed_list=12345:12348` for default sea-surface samples.
  - AWGN seed varies with condition, scenario, seed index, and Eb/N0 index.
- Smoke test: `SURFACE_COMM_MC_COUNT=2`.
- Default reduced run: `2` sea conditions, `2` scenarios, `mc_count=4`, `16` total channel runs.

Validation results:
- Smoke test completed successfully and produced finite `.mat`, table, and PNG outputs.
- Default 16-run communication Monte Carlo completed successfully.
- BER/SER sample array dimensions:
  - `ber_samples`: `[11, 4, 2, 2]`
  - `ser_samples`: `[11, 4, 2, 2]`
- Curve summary table rows: `44`, matching `2 sea conditions x 2 scenarios x 11 Eb/N0 values`.
- Channel invariant:
  - `max(abs(H_f-(H_direct_f+H_reflect_f)))=1.3877787807814457e-17`
- Direct-only reflected response:
  - `max(abs(H_reflect_f))=0`
- Finite-field check:
  - all selected BER/SER mean, standard deviation, and quantile fields are finite.
- BER mean curves for `EbN0_dB_list=[0,2,4,6,8,10,12,14,16,18,20]`:
  - weak direct-only: `[0.50025, 0.493375, 0.503, 0.496875, 0.50825, 0.498125, 0.49925, 0.503, 0.49525, 0.49875, 0.4995]`
  - weak direct-plus-reflect: `[0.497, 0.505625, 0.507, 0.500375, 0.51075, 0.50325, 0.498375, 0.512375, 0.5075, 0.498625, 0.507375]`
  - strong direct-only: `[0.504, 0.511125, 0.493875, 0.502125, 0.493875, 0.500875, 0.516625, 0.4955, 0.499875, 0.499375, 0.507375]`
  - strong direct-plus-reflect: `[0.513375, 0.502125, 0.497875, 0.50625, 0.50325, 0.505125, 0.504875, 0.505875, 0.507625, 0.5065, 0.494625]`
- SER mean curves:
  - weak direct-only: `[0.75075, 0.745, 0.7585, 0.75, 0.75975, 0.74775, 0.75925, 0.752, 0.7495, 0.7575, 0.75375]`
  - weak direct-plus-reflect: `[0.747, 0.755, 0.75125, 0.74825, 0.76275, 0.75825, 0.741, 0.76225, 0.754, 0.748, 0.7535]`
  - strong direct-only: `[0.75525, 0.76025, 0.74425, 0.7515, 0.7435, 0.7535, 0.766, 0.74725, 0.75925, 0.75225, 0.7615]`
  - strong direct-plus-reflect: `[0.76, 0.75075, 0.748, 0.7515, 0.75325, 0.751, 0.75575, 0.753, 0.7605, 0.75725, 0.73925]`
- Monotonicity check:
  - BER nonmonotonic increase counts matrix `[scenario x condition] = [[6, 6], [5, 4]]`
  - SER nonmonotonic increase counts matrix `[scenario x condition] = [[4, 6], [4, 4]]`
  - These reduced runs do not show the expected BER/SER decrease with Eb/N0; curves remain close to random hard-decision levels.
- Generated result artifacts:
  - `monte_carlo_comm_surface_psk_vertical_result.mat`
  - `monte_carlo_comm_surface_psk_vertical_BER_weak_strong_direct_plus_reflect.png`
  - `monte_carlo_comm_surface_psk_vertical_SER_weak_strong_direct_plus_reflect.png`
  - `monte_carlo_comm_surface_psk_vertical_BER_weak_direct_vs_reflect.png`
  - `monte_carlo_comm_surface_psk_vertical_BER_strong_direct_vs_reflect.png`
  - `monte_carlo_comm_surface_psk_vertical_SER_weak_direct_vs_reflect.png`
  - `monte_carlo_comm_surface_psk_vertical_SER_strong_direct_vs_reflect.png`

Remaining issues:
- The communication Monte Carlo infrastructure is functional, but current reduced BER/SER curves do not validate a useful demodulation-performance trend.
- The issue appears inherited from the existing baseband tap/MMSE communication chain behavior rather than the Monte Carlo aggregation itself.
- `n_sym=1000` and `mc_count=4` are reduced defaults; they are not sufficient for low-BER tail estimation.
- The script does not modify `comm_main_vertical_psk.m` and does not introduce a new BER/SER theory.

## 2026-06-10 Minimal Communication Closed-Loop Validation

Changed files for this validation:
- `validate_comm_link_minimal_vertical.m`: independent D1 communication-chain validation script.
- `PROJECT_CONTEXT.md`: records D1 validation purpose, cases, results, and diagnosis.

Validation settings:
- Communication: `M=4`, `n_sym=2000`, `EbN0_dB_list=[0,2,4,6,8,10,12,14,16,18,20]`, `bits_seed=9000`, `noise_seed_base=7000`.
- PE reduced direct-only case: `nx=ny=128`, `xw=yw=50`, `z_tx=100`, `z_rx=3`, `sigma_src_m=0.4`, `enable_wideband=true`, `f_band_hz=[4000 8000]`, `Nf_min=Nf_max=16`, `f_ref_hz=6000`, `enable_surface_reflection=false`.

Validation results:
- Unit channel:
  - noiseless `BER=0`, `SER=0`
  - `BER(20 dB)=0`, `SER(20 dB)=0`
  - BER curve `[0.11675, 0.04525, 0.0155, 0.0015, 0.00075, 0, 0, 0, 0, 0, 0]`
  - nonmonotonic counts: BER `0`, SER `0`
- Single complex tap:
  - noiseless `BER=0`, `SER=0`
  - `BER(20 dB)=0`, `SER(20 dB)=0`
  - BER curve `[0.22275, 0.15075, 0.07425, 0.0365, 0.01, 0.002, 0, 0, 0, 0, 0]`
  - nonmonotonic counts: BER `0`, SER `0`
- Known short multipath with pre-D2 legacy `conv(...,'same')` slicing:
  - `tap_count=3`, `peak_index=1`, `peak_fraction=0.8733624454148472`
  - noiseless `BER=0.506`, `SER=0.758`
  - `BER(20 dB)=0.50575`, `SER(20 dB)=0.7575`
  - This fails even without noise, so it is not an AWGN or PSK mapper issue.
- Same known short multipath with diagnostic `causal_head` slicing:
  - noiseless `BER=0`, `SER=0`
  - `BER(20 dB)=0`, `SER(20 dB)=0`
  - BER curve `[0.14575, 0.07975, 0.03475, 0.00675, 0.0005, 0, 0, 0, 0, 0, 0]`
  - nonmonotonic counts: BER `0`, SER `0`
- PE direct-only taps with pre-D2 legacy `conv(...,'same')` slicing:
  - `H_f` invariant error `=0`
  - `max(abs(H_reflect_f))=0`
  - `peak_index_before_alignment=1`
  - `h_full` peak fraction `=0.9976004631466816`
  - retained `tap_count=1997`
  - noiseless `BER=0.51775`, `SER=0.768`
  - `BER(20 dB)=0.49175`, `SER(20 dB)=0.739`
- PE direct-only peak-aligned taps with pre-D2 legacy `conv(...,'same')` slicing:
  - peak alignment does not change the result because the dominant tap was already at index 1.
  - noiseless `BER=0.51775`, `SER=0.768`
  - `BER(20 dB)=0.4925`, `SER(20 dB)=0.7415`
- PE direct-only peak-aligned taps with diagnostic `causal_head` slicing:
  - noiseless `BER=0`, `SER=0`
  - BER curve `[0.49725, 0.47675, 0.475, 0.4745, 0.47275, 0.4625, 0.45075, 0.431, 0.42525, 0.3695, 0.342]`
  - SER curve `[0.741, 0.7185, 0.7055, 0.704, 0.701, 0.677, 0.6655, 0.617, 0.6005, 0.5195, 0.479]`
  - nonmonotonic counts: BER `0`, SER `0`
  - effective SNR at nominal `20 dB` Eb/N0 is still `-4.9334090404706892 dB`, because the pre-D2 AWGN reference used transmit-symbol power while the PE channel gain is much smaller.

Diagnosis:
- The PSK mapper/demapper, AWGN injection, and single-tap equalization are functional.
- The pre-D2 communication slicing failed for multi-tap channels mainly because `conv(tx_symbols,h_bb,'same')` center-cropped a causal impulse response. For a long PE `h_bb`, this introduced a large symbol timing offset even when the dominant tap was already at index 1.
- Main-peak alignment alone is not sufficient for PE taps; receive-window selection must also respect causal timing.
- After causal receive-window selection, PE direct-only noiseless BER/SER becomes zero, but noisy BER remains high under the pre-D2 transmit-symbol Eb/N0 reference.

Generated result artifacts:
- `validate_comm_link_minimal_vertical_result.mat`
- `validate_comm_link_minimal_vertical_simple_channel_BER.png`
- `validate_comm_link_minimal_vertical_simple_channel_SER.png`
- `validate_comm_link_minimal_vertical_pe_h_bb_current_taps.png`
- `validate_comm_link_minimal_vertical_pe_h_bb_aligned_taps.png`

Remaining issues:
- D1 did not change production communication scripts; the later D2 update changed the default receive policy to `peak_sync` and the default Eb/N0 reference to `rx_clean`.
- A separate noise-calibration decision is needed: Eb/N0 can be referenced to transmit-symbol power, received clean-signal power, or post-equalizer noise enhancement, and these produce different BER interpretations.

## 2026-06-10 D2 Communication Receive Window and Eb/N0 Reference Update

Changed files for this update:
- `comm_main_vertical_psk.m`: switches the reference communication demo to explicit receive-window selection and received-signal Eb/N0 reference.
- `monte_carlo_comm_surface_psk_vertical.m`: applies the same communication policy to BER/SER Monte Carlo statistics.
- `validate_comm_link_minimal_vertical.m`: updates D1 diagnostics to include the D2 `peak_sync` policy.
- `PROJECT_CONTEXT.md`: records the D2 communication policy and validation status.

D2 policy:
- Receive window:
  - default `receive_window_mode='peak_sync'`.
  - Find the dominant tap index in `h_bb`.
  - Use `h_eq = h_bb(peak_index:end)` as the equalizer tap vector.
  - Compute `rx_full = conv(tx_symbols,h_eq,'full')`.
  - Use `rx_clean = rx_full(1:n_sym)` as the symbol-aligned receive vector.
  - Record `receive_meta.peak_index_original`, `receive_meta.original_tap_count`, `receive_meta.equalizer_tap_count`, and `receive_meta.discarded_pre_peak_energy_fraction`.
- Eb/N0 reference:
  - default `ebn0_reference='rx_clean'`.
  - `noise_inject_vertical(rx_clean, noise_cfg, rx_clean)` is used for the default communication BER/SER path.
  - local diagnostic support for `ebn0_reference='tx_symbols'` is retained in the scripts, but it is no longer the default.
- Compatibility:
  - `comm_main_vertical_psk.m` still saves the previous result fields and adds `receive_window_mode`, `receive_meta`, `ebn0_reference`, and `h_eq`.
  - `monte_carlo_comm_surface_psk_vertical.m` keeps the same result tables and sample arrays, and adds receive-window metadata in per-run summaries.
  - `comm_main_vertical_psk.m` remains the reference single-run demo; no changes were made to `vertical_channel_model` or propagation outputs.

Expected effect from D1 diagnosis:
- Unit and single-tap channels should remain zero-error in noiseless mode and should show BER/SER reduction with Eb/N0.
- Short causal multipath and PE direct-only channels should no longer fail noiseless validation due to `conv(...,'same')` center cropping.
- Receiver-referenced Eb/N0 should make nominal Eb/N0 correspond to the clean received waveform, so PE path loss no longer drives the effective SNR deeply negative.

Validation status:
- Static diff checks passed for the updated communication scripts.
- `validate_comm_link_minimal_vertical` completed successfully after D2.
- `monte_carlo_comm_surface_psk_vertical` completed successfully with `SURFACE_COMM_MC_COUNT=2`.
- `comm_main_vertical_psk` completed successfully with the reference `256 x 256`, `Nf=32` setup.

D2 minimal validation results:
- Unit channel:
  - noiseless `BER=0`, `SER=0`;
  - `BER(20 dB)=0`, `SER(20 dB)=0`;
  - BER decreases to zero with no nonmonotonic steps.
- Single complex tap:
  - noiseless `BER=0`, `SER=0`;
  - `BER(20 dB)=0`, `SER(20 dB)=0`;
  - BER decreases to zero with no nonmonotonic steps.
- Known short multipath:
  - legacy `same` still fails: noiseless `BER=0.506`, `SER=0.758`;
  - D2 `peak_sync` passes: noiseless `BER=0`, `SER=0`, `BER(20 dB)=0`, `SER(20 dB)=0`.
- PE direct-only reduced D1 case:
  - legacy `same` still fails: noiseless `BER=0.51775`, `SER=0.768`;
  - peak-aligned plus legacy `same` still fails, confirming alignment alone is insufficient;
  - D2 `peak_sync` passes: noiseless `BER=0`, `SER=0`, `BER(20 dB)=0`, `SER(20 dB)=0`.

D2 communication Monte Carlo smoke results:
- Setup: weak/strong sea conditions, direct-only/direct-plus-reflect, `SURFACE_COMM_MC_COUNT=2`, `EbN0_dB_list=[0,2,4,6,8,10,12,14,16,18,20]`.
- Weak direct-only BER mean:
  - `[0.111, 0.06375, 0.017, 0.004, 0, 0, 0, 0, 0, 0, 0]`
- Weak direct-plus-reflect BER mean:
  - `[0.1195, 0.06275, 0.01925, 0.00525, 0.00025, 0, 0, 0, 0, 0, 0]`
- Strong direct-only BER mean:
  - `[0.11825, 0.04925, 0.01525, 0.00275, 0.00025, 0, 0, 0, 0, 0, 0]`
- Strong direct-plus-reflect BER mean:
  - `[0.23825, 0.18975, 0.1535, 0.092, 0.064, 0.04425, 0.03075, 0.02725, 0.019, 0.016, 0.014]`
- The strong reflected case remains harder than direct-only, but now has a monotone decreasing BER trend instead of random-decision behavior.

D2 reference communication demo results:
- `comm_main_vertical_psk` now reports `rx_window=peak_sync`, `ebn0_ref=rx_clean`.
- Direct-only reference run:
  - `BER=[0.118, 0.05325, 0.02425, 0.0035, 0.00025, 0, 0, 0, 0, 0, 0]`
  - `SER=[0.1585, 0.0735, 0.0295, 0.005, 0.0005, 0, 0, 0, 0, 0, 0]`
- Direct-plus-reflect reference run:
  - `BER=[0.151, 0.102, 0.05975, 0.036, 0.01975, 0.0095, 0.004, 0.0015, 0.00075, 0.00025, 0]`
  - `SER=[0.21, 0.139, 0.0765, 0.0485, 0.026, 0.0135, 0.007, 0.0015, 0.001, 0.0005, 0]`
- Effective SNR now tracks nominal Eb/N0 plus QPSK bits-per-symbol scaling and no longer collapses to large negative values under PE path loss.

Remaining D2 issues:
- `peak_sync` is a deterministic dominant-tap synchronization rule; it is not a timing recovery loop.
- `rx_clean` Eb/N0 is a receiver-side performance convention. It is useful for comparing channel distortion at controlled received SNR, but it no longer includes absolute path-loss penalty in the noise power.
- The strong reflected channel can still have residual ISI/noise enhancement, so BER may remain nonzero at 20 dB for small `n_sym`/`mc_count` smoke tests.

## 2026-06-15 SSA-Like Statistical Surface Kernel

Changed files for this update:
- `vertical_channel_model.m`: adds validated `paramsV.surface_ssa_random_scatter`, `paramsV.surface_ssa_scatter_scale`, `paramsV.surface_ssa_seed_offset`, `paramsV.surface_ssa_kernel_mode`, `paramsV.surface_ssa_geometry_source_id`, `paramsV.surface_ssa_kz_branch`, `paramsV.surface_ssa_conv_padding`, and allows `surface_boundary_model='ssa_stat_kernel'`.
- `vertical_wape_propagator.m`: passes the SSA-like random-scatter controls and frequency index into the surface boundary module; keeps the two-segment reflected PE path and `H_f = H_direct_f + H_reflect_f`.
- `pm_surface_boundary_model.m`: adds the SSA-like PM-spectrum statistical boundary branch.
- `vertical_comm_guide.md` and `PROJECT_CONTEXT.md`: document formulas, interfaces, limits, and validation status.

Implemented interface and formulas:
- `ssa_stat_kernel` does not synthesize `xi(x,y)`. It scales the PM height spectrum directly so `sum(W_eta(:))*dkx*dky = sigma_eta^2`, where `sigma_eta = sea_hs_target/4`.
- The coherent term is now interpreted as the pressure-release / Dirichlet SSA coherent reflection coefficient:
  `R_coh = R0*exp(-0.5*(gamma_i+gamma_s)^2*sigma_eta^2)`.
  In normal specular reflection this reduces to `R_coh = -exp(-2*k0^2*sigma_eta^2)` for `R0=-1`.
- The raw incoherent power uses the corrected engineering baseline scale `P_sca_raw = surface_ssa_scatter_scale*abs(R0)^2*(k0*phase_factor_eff)^2*dkx*dky*circconv(W_eta,abs(Psi_inc_k).^2)`, where `phase_factor_eff` represents the effective `cos_i+cos_s` factor. The pre-correction `2*k0*phase_factor_eff` scale is retained only in metadata as a legacy audit value.
- Energy limiting enforces `E_coh+E_sca <= E_inc`; after random scatter synthesis, the final combined reflected spectrum is also checked against `E_inc`.
- The random scatter seed is `seed_ssa = sea_seed + surface_ssa_seed_offset + frequency_index - 1`, with default `surface_ssa_seed_offset=100000`.
- The current engineering kernel is explicitly identified by `surface_ssa_kernel_mode='pm_convolution'`.
- `pm_convolution` remains an engineering PM-spectrum baseline, not a strict SSA geometry kernel or a calibrated scattering cross section.
- `ssa1_geometry` and `ssa1_debug_dense` implement the `SSA.md` first-order Dirichlet geometry factor `G_SSA1(K,K';f)=4*gamma(K,f)*gamma(K',f)`.

Metadata:
- `roughness_meta.ssa_stat_kernel_meta` records `sigma_eta_m`, `Hs_target_m`, `R_coh`, `P_sca` stats, `E_inc`, `E_coh`, `E_sca_raw`, `E_sca`, `E_ref`, energy scale factors, `seed_ssa`, formulas, normalization notes, and limitations.
- The metadata now also records `kernel_mode`, `geometry_source_id`, `kz_branch`, `conv_padding`, `E_sca_limited`, `energy_limit_applied`, `energy_conservation_error`, and `propagating_bin_fraction`.
- For this branch, `surface_elevation=[]`, `delta_phi=[]`, and `roughness_meta.surface_realization_generated=false`.

Validation:
- Reduced scalar defaults: `nx=ny=128`, `xw=yw=50`, `z_tx=100`, `z_rx=3`, `sigma_src_m=0.4`, `show_figures=false`, `save_mode='rx_only'`, `enforce_1_over_R=false`.
- Default `kirchhoff_spatial`: invariant error `max(abs(H_f-H_direct_f-H_reflect_f)) = 0`.
- `ssa_stat_kernel`, `sea_hs_target=0`: `R_coh=-1`, `E_sca=0`, and flat-model difference from `kirchhoff_spatial` was `7.76e-17`.
- `ssa_stat_kernel`, `sea_hs_target=0.5`: `E_coh+E_sca = E_inc = 53831.9373`; same seed gave identical `H_f`, changed seed kept `H_direct_f` fixed and changed `H_reflect_f`.
- `surface_ssa_random_scatter=false`: `P_sca.sum=53831.9373` remained available in metadata, while `E_sca=0` and no random scatter was added to `psi_ref`.
- Direct-only with `ssa_stat_kernel`: `H_reflect_f=0` and the public channel fields remained present.
- Reduced wideband PSK smoke: direct-only and direct-plus-reflect both consumed `H_f` and completed the baseband tap/equalizer flow.
- Weak sea comparison with `sea_hs_target=0.05`: `kirchhoff_spatial` gave `|H|=0.074770323`, phase `-1.2009652`, BER `[0.11133 0 0]` at Eb/N0 `[0 10 20]`; `ssa_stat_kernel` gave `|H|=0.027797463`, phase `-1.7734442`, BER `[0.13086 0 0]`. The BER trend is consistent; pointwise channel equality is not expected.

Limitations and remaining issues:
- This is an SSA-like statistical kernel for validating the PM-spectrum-to-random-channel path. It is not a strict SSA/NLSSA scattering solver and does not include calibrated angular geometry factors or scattering cross sections.
- `surface_ssa_scatter_scale` is an engineering normalization knob, not a physical calibration.
- The old Kirchhoff realization models remain the default and should be used for legacy comparison unless the statistical branch is explicitly requested.

## 2026-06-15 SSA Kernel Mode Interface and SSA1 Geometry

Changed behavior:
- `surface_ssa_kernel_mode='pm_convolution'` is the default and reproduces the current engineering PM-spectrum convolution baseline.
- `surface_ssa_kernel_mode='ssa1_geometry'` implements the `SSA.md` first-order Dirichlet / pressure-release power geometry:
  `G_SSA1(K,K';f)=4*gamma(K,f)*gamma(K',f)`.
- The implemented FFT form is `A(K')=gamma(K',f)*abs(Psi_inc(K'))^2`, `B=circconv(W_eta,A)`, and `P_sca_raw(K)=4*C_norm*gamma(K,f)*B(K)*dkx*dky`, with `C_norm=surface_ssa_scatter_scale`.
- `surface_ssa_kernel_mode='ssa1_debug_dense'` computes the same periodic sum explicitly on small grids for FFT-vs-dense validation.
- `surface_ssa_kz_branch='downward_positive_real'` computes `kz(K)=sqrt(max(k0^2-|K|^2,0))`; non-propagating bins are excluded from the propagating fraction metadata but the baseline `pm_convolution` formula remains numerically unchanged.
- `surface_ssa_conv_padding='periodic'` keeps the existing FFT circular convolution.
- `surface_ssa_conv_padding='zero_padded'` now runs for `pm_convolution` and `ssa1_geometry`: inputs are `fftshift`ed to signed-k order, a full zero-padded linear convolution is computed on a `(2*ny-1,2*nx-1)` grid, and the original signed-k support is cropped back before returning to FFT order.

Current formula status:
- The first-order Dirichlet SSA geometry factor from `SSA.md` has been implemented.
- No calibrated angular scattering cross section has been implemented.
- No strict boundary-condition mapping from arbitrary `surface_reflect_coeff` to Dirichlet/Neumann/impedance SSA factors has been implemented; `ssa1_geometry` and `ssa1_debug_dense` require `surface_reflect_coeff=-1`.
- The code must not be described as NLSSA, impedance-boundary SSA, or experimentally calibrated rough-surface scattering.

Validation additions:
- `validate_ssa_stat_kernel_vertical.m` checks that `pm_convolution` records the new metadata fields, `ssa1_geometry` degenerates correctly for `Hs=0`, `zero_padded` convolution runs for both `pm_convolution` and `ssa1_geometry`, and `ssa1_debug_dense` matches the periodic FFT sum on a small grid through compact `P_sca` metadata.
- Sweep and communication summary tables carry kernel-mode energy audit fields: `E_sca_limited`, `energy_conservation_error`, `energy_limit_applied`, and `propagating_bin_fraction`.
- `sweep_ssa_stat_kernel_surface_channel_vertical.m` and `monte_carlo_comm_ssa_stat_kernel_psk_vertical.m` compare `kirchhoff_spatial`, `ssa_pm_convolution`, and `ssa1_geometry`.
- `plot_ssa_stat_kernel_report_vertical.m` emits `ssa_kernel_mode_*` figures; dense-vs-FFT is available through the validation result rather than default wideband sweeps.

Latest executed validation:
- `validate_ssa_stat_kernel_vertical` completed successfully.
- `Hs=0` flat degeneration:
  - `pm_convolution` vs `kirchhoff_spatial`: max response difference `7.7579e-17`.
  - `ssa1_geometry` vs `kirchhoff_spatial`: max response difference `7.7579e-17`.
- `ssa1_geometry` energy audit: `energy_conservation_error=1.3516e-16`, below the `1e-12` validation tolerance.
- `ssa1_debug_dense` vs FFT `ssa1_geometry`: compact `P_sca_raw` relative sum error `2.5247e-15`.
- `sweep_ssa_stat_kernel_surface_channel_vertical` completed with `SSA_STAT_MC_COUNT=1`, `Hs=[0 0.05 0.2 0.5]`, and model labels `kirchhoff_spatial`, `ssa_pm_convolution`, `ssa1_geometry`.
- In that sweep, both statistical kernels matched the `Hs=0` flat Kirchhoff response within `1.2533e-16`.
- `monte_carlo_comm_ssa_stat_kernel_psk_vertical` completed with `SSA_STAT_COMM_MC_COUNT=1`; all three model labels propagated through the existing `H_f` communication path and produced BER/SER summaries.
- `plot_ssa_stat_kernel_report_vertical` completed and regenerated `ssa_stat_kernel_report_*` plus `ssa_kernel_mode_*` report figures.

## 2026-06-15 SSA Statistical Validation and Report Scripts

Changed files for this validation/report update:
- `validate_ssa_stat_kernel_vertical.m`: new reduced scalar validation script.
- `sweep_ssa_stat_kernel_surface_channel_vertical.m`: new multi-Hs, multi-seed channel statistics script.
- `monte_carlo_comm_ssa_stat_kernel_psk_vertical.m`: new reduced QPSK BER/SER statistics script using the existing `H_f` communication path.
- `sweep_ssa_scatter_scale_sensitivity_vertical.m`: new reduced scalar sensitivity script for `surface_ssa_scatter_scale`.
- `plot_ssa_stat_kernel_report_vertical.m`: new plot-only report script that reads saved `.mat` results.
- `vertical_comm_guide.md` and `PROJECT_CONTEXT.md`: document run order, outputs, validation scope, and limitations.

Validation scope:
- Sea states: `sea_hs_target=[0 0.05 0.2 0.5]`, `sea_wind_speed=5`.
- Surface models: `kirchhoff_spatial` and `ssa_stat_kernel`.
- Channel sweep default seeds: `12345+(0:7)`, override with `SSA_STAT_MC_COUNT`.
- Communication default seeds: `12345+(0:3)`, override with `SSA_STAT_COMM_MC_COUNT`.
- Scatter-scale sweep default seeds: `12345+(0:3)`, override with `SSA_SCALE_SWEEP_MC_COUNT`.
- Scatter-scale sweep defaults: `sea_hs_target=[0.05 0.2 0.5]`, `surface_ssa_scatter_scale=[0 0.25 1 4]`, kernels `ssa_pm_convolution` and `ssa1_geometry`.
- Smoke-test condition limits: `SSA_STAT_SWEEP_MAX_CONDITIONS` and `SSA_STAT_COMM_MAX_CONDITIONS`.
- Scale-sweep smoke-test limit: `SSA_SCALE_SWEEP_MAX_CONDITIONS`.
- Reduced wideband channel defaults: `enable_wideband=true`, `f_band_hz=[4000 8000]`, `Nf_min=Nf_max=16`, `f_ref_hz=6000`, `nx=ny=128`, `show_figures=false`, `enforce_1_over_R=false`.

Acceptance metrics checked or summarized:
- `max(abs(H_f-H_direct_f-H_reflect_f)) <= 1e-10`.
- `h_total == H_f(idx_f_ref)` and `h_reflect == H_reflect_f(idx_f_ref)` within roundoff.
- `direct_only` keeps `H_reflect_f=0` and public channel fields present.
- `Hs=0` gives `sigma_eta_m=0`, `R_coh=surface_reflect_coeff`, `P_sca.sum=0`, `E_sca=0`, and a flat-response match against `kirchhoff_spatial`.
- `E_coh+E_sca <= E_inc + 1e-12*max(E_inc,1)` and `E_ref <= E_inc + 1e-12*max(E_inc,1)` for `ssa_stat_kernel` metadata.
- `surface_ssa_random_scatter=false` leaves finite `P_sca` metadata but sets realized `E_sca=0`, so downstream `H_reflect_f` receives only the coherent term.
- `surface_ssa_scatter_scale=0` gives zero raw scatter energy in both `pm_convolution` and `ssa1_geometry`.
- `E_sca_raw/E_inc` should be nondecreasing as `surface_ssa_scatter_scale` increases for a fixed kernel, sea state, and seed.
- Statistical summaries report mean/std of `|H(f)|`, `|H_reflect(f)|`, `angle(H(f_ref))`, `|h_total|`, `|h_reflect|`, `E_sca/E_inc`, `E_ref/E_inc`, and energy scaling.
- BER/SER summaries compare qualitative trends with paired bit/noise seeds, not pointwise equality between physical models.

Report outputs:
- Channel result: `sweep_ssa_stat_kernel_surface_channel_vertical_result.mat`.
- Communication result: `monte_carlo_comm_ssa_stat_kernel_psk_vertical_result.mat`.
- Validation result: `validate_ssa_stat_kernel_vertical_result.mat`.
- Scale-sensitivity result: `sweep_ssa_scatter_scale_sensitivity_vertical_result.mat`.
- Figures:
  - `ssa_stat_kernel_report_abs_H_f_mean_std.png`
  - `ssa_stat_kernel_report_abs_H_reflect_f_mean_std.png`
  - `ssa_stat_kernel_report_phase_H_f_ref_vs_Hs.png`
  - `ssa_stat_kernel_report_reflect_energy_vs_Hs.png`
  - `ssa_stat_kernel_report_Esca_Einc_vs_Hs.png`
  - `ssa_stat_kernel_report_energy_scale_vs_Hs.png`
  - `ssa_stat_kernel_report_abs_h_total_model_compare.png`
  - `ssa_stat_kernel_report_abs_h_reflect_model_compare.png`
  - `ssa_stat_kernel_report_BER_model_compare.png`
  - `ssa_stat_kernel_report_SER_model_compare.png`
  - `ssa_stat_kernel_report_Hs0_flat_check.png`
  - `ssa_stat_kernel_report_metadata_energy_table.mat`

Limitations:
- These scripts are validation and reporting harnesses around the current implementation. They do not modify `comm_main_vertical_psk.m`, modem/noise/MMSE policy, WAPE marching, or `local_march_field`.
- Full spatial fields, `psi_ref`, `P_sca`, and `W_eta` matrices are intentionally not saved in the statistical result files.
- The current model remains an SSA-like engineering/statistical kernel. Strict SSA angular geometry factors and calibrated scattering cross sections remain future work.

Latest scale-sensitivity smoke result:
- `sweep_ssa_scatter_scale_sensitivity_vertical` completed with `SSA_SCALE_SWEEP_MC_COUNT=1`, `Hs=[0.05 0.2 0.5]`, and scales `[0 0.25 1 4]`.
- For `Hs=0.05`, `E_sca_raw/E_inc` was `[0 0.39478 1.5791 6.3165]` for `ssa_pm_convolution` and `[0 0.09768 0.39072 1.5629]` for `ssa1_geometry`.
- `energy_conservation_error_max=0` for every listed scale-sweep condition.
- The sweep uses `surface_ssa_random_scatter=false`; therefore `E_sca/E_inc=0` in the realized reflected field, while `E_sca_limited/E_inc` records the limited scatter-energy budget for metadata analysis.

## 2026-06-16 Core MATLAB File Rename

Changed core names:
- `CARPE3D_vertical.m` -> `vertical_channel_model.m`.
- `propWAPE_vertical.m` -> `vertical_wape_propagator.m`.
- `pm_surface_kirchhoff_module.m` -> `pm_surface_boundary_model.m`.

Behavioral intent:
- This is a naming refactor only. It does not change `paramsV`, `cfg`, `output`, `H_f`, `H_direct_f`, `H_reflect_f`, `h_total`, or communication-chain result semantics.
- No compatibility wrappers for the old core names are kept; active scripts call the new names directly.
- The `old/` directory remains archival and was not migrated.

Validation after rename:
- Reduced scalar smoke through `vertical_channel_model` passed with `max(abs(H_f-H_direct_f-H_reflect_f))=0`.
- Direct-only smoke passed with `max(abs(H_reflect_f))=0`.
- `validate_ssa_stat_kernel_vertical` passed after updating expected error IDs to `pm_surface_boundary_model:*`.
- `monte_carlo_comm_ssa_stat_kernel_psk_vertical` completed with the renamed entrypoint and existing `H_f` communication path.
- `git diff --check` passed.

## 2026-06-17 SSA-Like Statistical Kernel Increment

Purpose:
- Continue treating `ssa_stat_kernel` as a PM-spectrum-driven SSA-like statistical scattering kernel prototype.
- Preserve `kirchhoff_spatial` as the realization-based phase-screen comparison path.
- Keep `pm_convolution` as an engineering baseline and `ssa1_geometry` as the current first-order Dirichlet research kernel.
- Do not claim complete SSA/NLSSA, T-matrix, impedance-boundary, multiple-scattering, or experimentally calibrated sea-surface scattering physics.

Implementation updates:
- `pm_convolution` and `ssa1_geometry` now support both `surface_ssa_conv_padding='periodic'` and `surface_ssa_conv_padding='zero_padded'`.
- `periodic` remains the compatible FFT circular convolution path.
- `zero_padded` computes an FFT-accelerated linear convolution using `(2*ny-1,2*nx-1)` padding, signed-k ordering via `fftshift`, center cropping back to the original signed-k support, and `ifftshift` back to the solver's FFT order.
- `ssa1_debug_dense` remains a small-grid explicit periodic-sum validator for the periodic FFT path.
- `compare_ssa_conv_padding_vertical.m` adds a reduced scalar-frequency comparison between `periodic` and `zero_padded` for `pm_convolution` and `ssa1_geometry`.

Metadata additions:
- `ssa_stat_kernel_meta` records `conv_operator`, `conv_padding_size`, `conv_crop_start_index`, `conv_crop_end_index`, and `conv_crop_rule`.
- It also records compact angular-spectrum diagnostics for incident, coherent, scatter-power, and reflected spectra: total spectral energy, centroid, RMS transverse wavenumber, and 90% energy radius.
- Existing energy audit fields remain unchanged: `E_inc`, `E_coh`, `E_sca_raw`, `E_sca_limited`, `E_sca`, `E_ref`, `energy_scale_applied`, `energy_limit_applied`, and `energy_conservation_error`.

Validation updates:
- `validate_ssa_stat_kernel_vertical` now expects `zero_padded` to run for both `pm_convolution` and `ssa1_geometry` and records its raw-scatter-energy difference from the periodic path.
- The validation still covers `Hs=0`, direct-only invariants, deterministic seeding, changed-seed reflection sensitivity, `random_scatter=false`, scale monotonicity, non-Dirichlet rejection for `ssa1_geometry`, and periodic dense-vs-FFT consistency.
- Wideband statistics scripts now carry convolution mode/operator and angular-spectrum broadening metrics in their run and condition summaries.

Latest numerical validation:
- `validate_ssa_stat_kernel_vertical` completed successfully after the zero-padding implementation.
- `Hs=0` flat degeneration remained at roundoff level:
  - `pm_convolution` vs flat `kirchhoff_spatial`: max response difference `7.7579e-17`.
  - `ssa1_geometry` vs flat `kirchhoff_spatial`: max response difference `7.7579e-17`.
- Energy audit checks remained within tolerance:
  - `ssa1_geometry` energy conservation error `1.3516e-16`.
  - zero-padded `pm_convolution` energy conservation error `0`.
  - zero-padded `ssa1_geometry` energy conservation error `0`.
- `ssa1_debug_dense` still matches the periodic FFT implementation for compact raw scatter power with relative sum error `2.5247e-15`.
- The reduced zero-padded wideband invariant check passed with `max(abs(H_f-H_direct_f-H_reflect_f))=8.9974e-19` over four frequency bins.
- `sweep_ssa_stat_kernel_surface_channel_vertical` completed a light statistical sweep with one seed over `Hs=[0 0.05 0.2 0.5]` and model labels `kirchhoff_spatial`, `ssa_pm_convolution`, and `ssa1_geometry`.
- In the light sweep, both statistical kernels matched the `Hs=0` flat Kirchhoff response within `1.2533e-16`.
- `compare_ssa_conv_padding_vertical` completed with one seed over `Hs=[0 0.05 0.2]` and kernels `pm_convolution` / `ssa1_geometry`.
- Periodic-vs-zero-padded comparison highlights:
  - `Hs=0`: both convolution paths match the flat `kirchhoff_spatial` reference within `7.7579e-17`.
  - `pm_convolution`, `Hs=0.05`: `E_sca_raw/E_inc` changed from `1.5791` to `1.5784`; `rel_diff_H_reflect_f=0.0031907`; `rel_diff_H_f=0.0032389`.
  - `pm_convolution`, `Hs=0.2`: `E_sca_raw/E_inc` changed from `25.266` to `25.254`; `rel_diff_H_reflect_f=0.0063383`; `rel_diff_H_f=0.0019823`.
  - `ssa1_geometry`, `Hs=0.05`: `E_sca_raw/E_inc` changed from `0.39072` to `0.39054`; `rel_diff_H_reflect_f=0.0030922`; `rel_diff_H_f=0.0031332`.
  - `ssa1_geometry`, `Hs=0.2`: `E_sca_raw/E_inc` changed from `6.2515` to `6.2486`; `rel_diff_H_reflect_f=0.0061289`; `rel_diff_H_f=0.0019182`.
  - All periodic-vs-zero-padded comparison rows kept `energy_conservation_error=0`; maximum invariant error stayed below `3.6e-18`.
- Comparison artifacts:
  - `compare_ssa_conv_padding_vertical_result.mat`
  - `ssa_conv_padding_compare_raw_scatter_energy.png`
  - `ssa_conv_padding_compare_relative_difference.png`

Documentation update:
- `vertical_comm_guide.md` was rewritten as a human-readable mathematical model note.
- The guide separates `kirchhoff_spatial`, `pm_convolution`, and `ssa1_geometry` by physical role.
- It states that `periodic` is the default compatible convolution path and `zero_padded` is a linear-convolution aliasing audit path.
- It keeps the current scope conservative: PM-spectrum-driven SSA-like statistical scattering prototype, not complete SSA/NLSSA, T-matrix, impedance-boundary, or experimentally calibrated rough-surface scattering physics.
- It documents that periodic and zero-padded paths are not expected to match pointwise; acceptance is based on channel invariants, energy audit, and interpretable trends.
- It now states two additional limitations:
  - statistical scattering realization frequency correlation is still simplified, so wideband `H_f` continuity needs a dedicated model;
  - under fixed `sea_hs_target`, wind-speed sweeps primarily change PM spectral shape and should not be read as monotonic sea-state intensity sweeps.

Stage conclusion:
- Current status: PM-spectrum-driven SSA-like statistical scattering has engineering integration, `periodic`/`zero_padded` convolution paths, energy audit, compact metadata, core reduced-grid validation, and an initial periodic-vs-zero-padded comparison report.
- The module can serve as a data-generation component for random-channel and communication-algorithm evaluation.
- It remains outside the scope of complete physical SSA/NLSSA sea-surface scattering theory, calibrated absolute scattering cross sections, impedance-boundary SSA, and higher-order multiple scattering.

Recommended next stage, not implemented here:
- frequency-correlated random scattering across wideband bins;
- larger multi-seed convergence statistics;
- numerical or experimental calibration of `surface_ssa_scatter_scale`;
- larger-scale linkage between scatter statistics and BER/SER trends.

## 2026-06-17 SSA Coherent Reflection Formula Correction

Purpose:
- Correct the `ssa_stat_kernel` coherent reflection coefficient so it matches the pressure-release / Dirichlet first-order SSA coherent reflection interpretation.
- At this stage the Kirchhoff realization phase-screen path was not changed; a later correction updates its phase convention to `delta_phi=k0*(cos_i+cos_r)*eta`.
- Keep `pm_convolution` as an engineering PM-spectrum convolution baseline and `ssa1_geometry` as the first-order pressure-release / Dirichlet scattering geometry.

Formula correction:
- Previous SSA coherent implementation used:
  `R_coh = R0*exp(-0.5*(2*k0*phase_factor_eff)^2*sigma_eta^2)`.
- Because `phase_factor_eff=2` in normal mode, that produced the unintended normal exponent `-8*k0^2*sigma_eta^2`.
- The corrected SSA coherent implementation uses:
  `R_coh = R0*exp(-0.5*(gamma_i+gamma_s)^2*sigma_eta^2)`.
- For normal specular reflection, `gamma_i=gamma_s=k0`, so with `R0=-1`:
  `R_coh = -exp(-2*k0^2*sigma_eta^2)`.
- Broschat 1993 is used as the coherent reflection coefficient reference for PM rough surfaces; it is not treated as the source of the noncoherent `P_sca` distribution.

Implementation notes:
- `phase_factor_eff` is still retained as the effective `cos_i+cos_s` geometry factor.
- `coherent_gamma_sum_eff_rad_per_m = k0*phase_factor_eff` is now the quantity used in the SSA coherent exponent.
- The SSA correction did not change Kirchhoff at that time. The current Kirchhoff phase screen now uses `delta_phi=k0*(cos_i+cos_r)*eta`, giving `2*k0*eta` for normal incidence.
- The `pm_convolution` noncoherent baseline remains an engineering scalar normalization; after the Kirchhoff phase audit its phase scale is `k0*phase_factor_eff`, with the older `2*k0*phase_factor_eff` value kept only as legacy metadata.
- `ssa1_geometry` still computes the first-order Dirichlet scattering power through `G_SSA1=4*gamma_s*gamma_i` and still rejects non-Dirichlet `surface_reflect_coeff`.

Metadata additions:
- `R_coh_raw`, `R_coh`.
- `coherent_exponent`.
- `coherent_gamma_sum_eff_rad_per_m`.
- `coherent_vertical_factor_eff`.
- `coherent_reflection_formula`.
- `phase_factor_eff_legacy_note`.
- `W_eta_variance_target`, `W_eta_variance_discrete`, `W_eta_variance_rel_error`.

Validation results:
- `validate_ssa_stat_kernel_vertical` completed successfully.
- Normal-mode coherent reflection matched `R0*exp(-2*k0^2*sigma_eta^2)` within roundoff.
- The validation explicitly rejected the old unintended `R0*exp(-8*k0^2*sigma_eta^2)` behavior; the checked difference was `2.6753e-09`.
- `coherent_gamma_sum_eff_rad_per_m` equaled `2*k0` in normal mode.
- PM spectrum normalization passed: `W_eta_variance_rel_error <= 8.6736e-15` across the Hs trend cases.
- `Hs=[0 0.05 0.2 0.5]` gave monotone decreasing `|R_coh|` and nondecreasing `E_sca_limit`.
- `ssa1_geometry` metadata recorded `G_SSA1` and `pressure-release / Dirichlet`; non-Dirichlet reflection coefficient rejection still passed.
- `ssa1_debug_dense` vs periodic FFT raw scatter relative sum error remained `2.5247e-15`.
- `compare_ssa_conv_padding_vertical` completed after the coherent correction, with energy conservation error `0` for all compared rows.

Recommended naming:
- In papers or method notes, describe the module as:
  `PM-spectrum-driven first-order pressure-release / Dirichlet SSA statistical scattering branch`.
- Avoid naming it complete SSA, NLSSA, T-matrix, impedance-boundary scattering, multiple-scattering, or experimentally calibrated sea-surface scattering.

## 2026-06-17 SSA Physical Trend Validation Closeout

Purpose:
- Add a reduced-grid physical trend validation for the sea-surface SSA branch without running or modifying the communication chain.
- Close the wording around the current SSA module as a PM-spectrum-driven first-order pressure-release Dirichlet SSA statistical reflection/scattering model.
- Keep `kirchhoff_spatial` as the default and preserve the existing `H_f = H_direct_f + H_reflect_f` channel semantics.

New validation script:
- `validate_ssa_physical_trends_vertical.m` scans:
  - `Hs=[0 0.05 0.2 0.5 1.0]`;
  - `f=[4000 6000 8000 10000]` Hz;
  - `sea_wind_speed=5`;
  - reduced grid `nx=ny=128`;
  - `surface_ssa_random_scatter=false`;
  - `surface_ssa_conv_padding='periodic'`.
- Compared kernels:
  - `pm_convolution`, retained as an engineering PM-spectrum convolution baseline;
  - `ssa1_geometry`, retained as the first-order pressure-release / Dirichlet geometry kernel.
- The script saves compact metadata tables only; it does not save full 2-D fields or full spectra.

Recorded quantities:
- Frequency, target sea state, `sigma_eta`, `k0`, `R_coh`, `abs(R_coh)`, coherent exponent, effective vertical factor, and PM spectrum variance audit.
- Energy audit fields: `E_inc`, `E_coh`, `E_sca_raw`, `E_sca_limited`, `E_sca`, `E_sca_limit`, `E_ref`, and `energy_conservation_error`.
- Kernel metadata: `kernel_mode`, `conv_padding`, `kernel_formula`, `G_SSA1_formula`, `boundary_condition`, `formula_source`, and limitations.

Trend checks:
- `Hs=0`: `R_coh=R0`, `W_eta_variance_discrete=0`, `E_sca_raw=0`, and `E_sca=0`.
- Normal coherent formula at every scanned frequency and sea state:
  `R_coh = R0*exp(-2*k0^2*sigma_eta^2)`.
- Fixed frequency and increasing `Hs`: `abs(R_coh)` and `E_coh/E_inc` must not increase; `E_sca_limit/E_inc` must not decrease.
- Fixed nonzero `Hs` and increasing frequency: `abs(R_coh)` must not increase and the coherent exponent must become more negative.
- `ssa1_geometry` must record `G_SSA1(Ks,Ki;f)=4*gamma(Ks)*gamma(Ki)` and `pressure-release / Dirichlet`, and it must reject non-Dirichlet `surface_reflect_coeff`.
- `pm_convolution` must remain clearly labeled as an engineering baseline and must not be described as strict SSA.

Expected artifacts:
- `validate_ssa_physical_trends_vertical_result.mat`, with `summary_table`, `trend_table`, and `validation_report`.
- `validate_ssa_abs_R_coh_vs_Hs.png`.
- `validate_ssa_abs_R_coh_vs_frequency.png`.
- `validate_ssa_coherent_energy_fraction_vs_Hs.png`.
- `validate_ssa_scatter_budget_fraction_vs_Hs.png`.
- `validate_ssa_W_eta_variance_check.png`.

Documentation update:
- `vertical_comm_guide.md` now states that Broschat 1993 is used for coherent reflection coefficient interpretation, not as the source of `P_sca`.
- It separates the Kirchhoff realization phase screen from the SSA coherent statistical average exponent.
- It records the first-order Dirichlet noncoherent power kernel:
  `P_sca(Ks)=4*C_norm*gamma(Ks)*[W_eta*(gamma(Ki)*abs(Psi_inc(Ki))^2)](Ks)*DeltaKx*DeltaKy`.
- It states that `C_norm` / `surface_ssa_scatter_scale` is an engineering normalization parameter, not an experimentally calibrated absolute scattering cross-section constant.

Remaining physical boundaries:
- At this stage there was no second-order SSA coherent correction. This historical limitation was later superseded by the coherent-only `ssa2_broschat_coherent` diagnostic documented in the 2026-06-23 section.
- No NLSSA.
- No full T-matrix solver.
- No impedance-boundary or Neumann-boundary SSA factor.
- No multiple scattering.
- No evanescent scattering injection into the communication reflected field.
- No experimentally calibrated bistatic scattering cross section.

## 2026-06-18 SSA1 vs Kirchhoff Surface Statistics Comparison

Purpose:
- Compare the PM-spectrum-driven first-order pressure-release Dirichlet SSA statistical reflection/scattering model against the default Kirchhoff realization phase-screen model.
- Validate statistical trend compatibility, not sample-wise equivalence.
- Keep communication scripts, BER/SER, WAPE marching, default `kirchhoff_spatial`, and `H_f = H_direct_f + H_reflect_f` semantics unchanged.

New comparison script:
- `compare_ssa1_kirchhoff_surface_statistics_vertical.m`.
- Default scan:
  - `Hs=[0 0.05 0.2 0.5 1.0]`;
  - `f=[4000 6000 8000 10000]` Hz;
  - `sea_seed=12345+(0:15)`;
  - `sea_wind_speed=5`;
  - reduced grid `nx=ny=64`;
  - `show_figures=false`, `save_mode='rx_only'`, `enforce_1_over_R=false`;
  - `surface_boundary_redistribution_diagnostics=true`.
- Default model set:
  - `kirchhoff_spatial`;
  - `ssa_stat_kernel` with `surface_ssa_kernel_mode='ssa1_geometry'`.
- Optional engineering baseline:
  - `ssa_stat_kernel` with `surface_ssa_kernel_mode='pm_convolution'`;
  - enabled only by environment flag, and not described as strict SSA.

Outputs:
- `compare_ssa1_kirchhoff_surface_statistics_vertical_result.mat`, containing:
  - `run_table`;
  - `summary_table`;
  - `trend_table`;
  - `validation_report`.
- Report figures:
  - `compare_ssa1_kirchhoff_abs_h_reflect_vs_Hs.png`;
  - `compare_ssa1_kirchhoff_abs_R_coh_vs_Hs.png`;
  - `compare_ssa1_kirchhoff_coherent_fraction_vs_Hs.png`;
  - `compare_ssa1_kirchhoff_scatter_budget_vs_Hs.png`;
  - `compare_ssa1_kirchhoff_rms_delta_k_vs_Hs.png`;
  - `compare_ssa1_kirchhoff_phase_variance_vs_Hs.png`.

Recorded run-level quantities:
- Channel statistics: `h_reflect`, `abs_h_reflect`, `phase_h_reflect`, `h_total`, `abs_h_total`, and `h_reflect/h_direct` magnitude/phase.
- SSA metadata where applicable: `R_coh`, `abs_R_coh`, `E_inc`, `E_coh`, `E_sca`, `E_sca_limit`, `E_ref`, and `energy_conservation_error`.
- Angular-spectrum spread:
  - SSA rows use `ssa_stat_kernel_meta.reflected_spectrum_stats`;
  - Kirchhoff rows use `boundary_redistribution_diagnostics.reflect_rms_delta_k_rad_per_m`, `reflect_high_k_fraction`, and `centroid_shift_mag_rad_per_m`;
  - unavailable fields remain `NaN` rather than fabricated.

Validation interpretation:
- Hard checks:
  - all rows satisfy `max(abs(H_f-H_direct_f-H_reflect_f)) <= 1e-10`;
  - `Hs=0` `ssa1_geometry` and `kirchhoff_spatial` reflected/total channel samples match to roundoff;
  - `ssa1_geometry` keeps `energy_conservation_error <= 1e-12`;
  - `ssa1_geometry` metadata records `G_SSA1=4*gamma(Ks)*gamma(Ki)` and `pressure-release / Dirichlet`.
- Compatibility observations:
  - as `Hs` increases, SSA1 `abs(R_coh)` and `E_coh/E_inc` decrease while `E_sca_limit/E_inc` increases;
  - Kirchhoff reflected spectrum should broaden from flat to rough sea through RMS delta-k or high-k fraction;
  - frequency trends are recorded as compatibility evidence, not as proof of realization-by-realization equivalence.

Full 16-seed reduced-grid result:
- The default comparison was completed with the default 16 seeds, 5 sea states, 4 frequencies, 2 models, and `nx=ny=64`.
- Total reduced scalar cases: `640`.
- Result tables:
  - `run_table`: 640 rows;
  - `summary_table`: 40 rows;
  - `trend_table`: 156 rows.
- Hard checks:
  - `hard_checks_passed=true`;
  - `hard_fail_count=0`;
  - maximum channel invariant error `max(abs(H_f-H_direct_f-H_reflect_f))=1.5516e-17`;
  - maximum `Hs=0` SSA1/Kirchhoff reflected/total channel difference `2.2205e-16`;
  - maximum SSA1 `energy_conservation_error=9.9322e-16`.
- Compatibility checks:
  - `compatibility_observed_fraction=1.000`;
  - `compat_fail_count=0`.
- Kirchhoff reflected-spectrum broadening from `Hs=0` to `Hs=1.0` was observed at every tested frequency through the RMS delta-k metric:
  - 4000 Hz: increase `3.4532 rad/m`;
  - 6000 Hz: increase `3.2032 rad/m`;
  - 8000 Hz: increase `3.1521 rad/m`;
  - 10000 Hz: increase `3.1419 rad/m`.
- Kirchhoff frequency-variability compatibility metrics were positive for every nonzero tested sea state:
  - `Hs=0.05`: `0.80958`;
  - `Hs=0.2`: `1.8987`;
  - `Hs=0.5`: `0.23620`;
  - `Hs=1.0`: `0.058270`.
- The six `compare_ssa1_kirchhoff_*.png` figures were regenerated from the 16-seed result and copied into `ssa1_geometry_report_figures/`, replacing the earlier smoke-summary figures.

Scope boundary:
- This comparison does not prove Kirchhoff and SSA1 are pointwise equivalent.
- Kirchhoff remains a concrete rough-surface realization phase-screen model.
- SSA1 remains a PM-spectrum-driven first-order pressure-release Dirichlet statistical reflection/scattering model.
- At the time of this comparison, the project did not include second-order SSA. The later `ssa2_broschat_coherent` addition covers only the coherent reflection coefficient; second-order incoherent scattering, NLSSA, T-matrix, impedance boundary, multiple scattering, evanescent scattering injection, and experimentally calibrated scattering cross sections remain outside the current model.

## 2026-06-23 SSA1 vs Broschat-Style SSA2 Coherent Reflection Comparison

Purpose:
- Add an optional coherent reflection order inside `ssa_stat_kernel` for comparing the default first-order coherent coefficient against a Broschat 1993-style second-order coherent correction.
- Keep default `surface_ssa_coherent_order='ssa1'`.
- Keep `kirchhoff_spatial`, WAPE marching, communication scripts, BER/SER logic, random scatter synthesis, and the first-order `P_sca` formula unchanged.

Implemented interface:
- New validated parameter:
  - `surface_ssa_coherent_order='ssa1'` by default;
  - optional value `surface_ssa_coherent_order='ssa2_broschat_coherent'`.
- The parameter is passed from `vertical_channel_model.m` through `vertical_wape_propagator.m` to `pm_surface_boundary_model.m`.
- `ssa2_broschat_coherent` is accepted only for pressure-release / Dirichlet reflection, `surface_reflect_coeff=-1`; otherwise the code raises:
  `pm_surface_boundary_model:Ssa2DirichletReflectCoeffRequired`.

Implemented formula:
- SSA1 coherent reflection remains:
  `R1 = R0*exp(-0.5*(gamma_i+gamma_s)^2*sigma_eta^2)`.
- In normal pressure-release reflection:
  `R1 = -exp(-2*k0^2*sigma_eta^2)`.
- The optional Broschat-style coherent correction uses the user-specified 2-D extension:
  `R2 = R1 + 2*gamma_i*exp(-2*gamma_i^2*sigma_eta^2)*sum(W_eta(q)*(gamma(K_i+q)-gamma_i))*dkx*dky`.
- The square-root branch is:
  `gamma(K)=sqrt(complex(k0^2-|K|^2,0))`.
  This keeps the full complex coherent integral and records propagating-only and evanescent contributions separately.
- For normal tests, `K_i=[0,0]` and `gamma_i=k0`.
- For oblique mode, the implementation records an effective 2-D approximation using an assumed `+kx` incident azimuth; physical conclusions should rely on normal or near-vertical cases.

Metadata additions:
- `coherent_order`.
- `R_coh_ssa1`, `R_coh_ssa2`, and selected `R_coh`.
- `abs_R_coh_ssa1`, `abs_R_coh_ssa2`.
- `coherent_loss_ssa1_db`, `coherent_loss_ssa2_db`, `coherent_loss_delta_db`.
- `delta_R_abs`, `delta_R_rel`.
- `ssa2_correction_integral`, real/imag components, `ssa2_gamma_i_eff_rad_per_m`, `ssa2_K_i_eff_rad_per_m`.
- `ssa2_sqrt_branch`, `ssa2_formula`, `ssa2_formula_source`, `ssa2_limitations`, `ssa2_evanescent_handling`.
- `ssa2_propagating_integral`, `ssa2_evanescent_integral`, and `ssa2_evanescent_fraction_abs`.

New validation script:
- `compare_ssa1_ssa2_coherent_reflection_vertical.m`.
- Default scan:
  - `Hs=[0 0.05 0.2 0.5 1.0]`;
  - `f=[4000 6000 8000 10000]` Hz;
  - `surface_phase_mode='normal'`;
  - `surface_boundary_model='ssa_stat_kernel'`;
  - `surface_ssa_kernel_mode='ssa1_geometry'`;
  - `surface_ssa_random_scatter=false`;
  - reduced grid `nx=ny=128`.
- The script saves:
  - `compare_ssa1_ssa2_coherent_reflection_vertical_result.mat`;
  - `run_table`, `summary_table`, and `validation_report`.
- Figures:
  - `compare_ssa1_ssa2_abs_R_coh_vs_Hs.png`;
  - `compare_ssa1_ssa2_coherent_loss_vs_Hs.png`;
  - `compare_ssa1_ssa2_loss_delta_vs_Hs.png`;
  - `compare_ssa1_ssa2_loss_delta_vs_frequency.png`;
  - `compare_ssa1_ssa2_negligible_region_heatmap.png`;
  - `compare_ssa1_ssa2_evanescent_fraction_heatmap.png`.

Executed validation:
- `compare_ssa1_ssa2_coherent_reflection_vertical` completed successfully.
- All checks passed:
  - channel invariant maximum `4.0492e-18`;
  - maximum energy conservation error `0`;
  - maximum `W_eta` variance relative error `8.6736e-15`;
  - `Hs=0` gives `R1=R2=R0`;
  - normal SSA1 formula matched exactly within reported roundoff;
  - `ssa2_broschat_coherent` rejected non-Dirichlet `surface_reflect_coeff=-0.8`.
- Maximum absolute coherent-loss delta over the default grid was `0.370648 dB`.
- With the `0.1 dB` negligible threshold:
  - weak to moderate cases through `Hs<=0.5 m` stayed at or below about `0.1 dB`;
  - the strongest low-frequency point `Hs=1 m`, `f=4000 Hz` exceeded the threshold at about `0.37065 dB`;
  - high-frequency high-roughness cases can saturate the loss calculation near machine epsilon and should not be overinterpreted as accurate absolute loss.
- `validate_ssa_stat_kernel_vertical` completed successfully after the new coherent-order metadata additions.

Interpretation:
- `ssa2_broschat_coherent` is a coherent-coefficient diagnostic, not a complete second-order SSA scattering model.
- The first-order `ssa1_geometry` noncoherent scattering power remains:
  `G_SSA1=4*gamma_s*gamma_i`.
- The current recommendation is not to make SSA2 the default. It can be used as a sensitivity check, especially for strong roughness and lower frequencies.
- A true second-order rough-surface model would still require second-order incoherent scattering, consistent normalization, and a broader validation basis.

## 2026-06-23 SSA Scatter-Scale Calibration and Frequency-Correlation Diagnostics

Purpose:
- Keep the core WAPE propagation, communication scripts, BER/SER logic, and default `kirchhoff_spatial` behavior unchanged.
- Add reduced-grid engineering diagnostics for two practical SSA statistical-kernel questions:
  - how to choose the numerical normalization `surface_ssa_scatter_scale` against a Kirchhoff realization-based reference;
  - how sensitive wideband random reflected channels are to the assumed cross-frequency correlation of the random scatter spectrum.

Checkpoint status:
- A checkpoint commit was created before these follow-up diagnostics:
  - `ca913a5 feat: add SSA coherent diagnostics and validation scripts`.
- Push to `origin/codex/rename-core-matlab-files` was attempted from the current environment but was blocked by the environment's external remote export policy. The local commit remains available.

New random-scatter frequency-correlation interface:
- `vertical_channel_model.m` now validates:
  - `surface_ssa_frequency_correlation_mode`, default `'independent'`;
  - allowed values: `'independent'`, `'shared_seed_phase'`, `'ar1_frequency'`;
  - `surface_ssa_frequency_correlation_rho`, default `0.8`, constrained to `0 <= rho < 1`.
- `vertical_wape_propagator.m` passes these fields into the surface-boundary configuration.
- `pm_surface_boundary_model.m` uses them only when `surface_boundary_model='ssa_stat_kernel'` and `surface_ssa_random_scatter=true`.
- Default `'independent'` preserves the previous seed rule and random spectrum generation.
- `'shared_seed_phase'` reuses the same complex Gaussian random scatter spectrum across frequency bins.
- `'ar1_frequency'` creates a deterministic first-order autoregressive complex Gaussian sequence along the frequency index.
- Metadata records:
  - `frequency_correlation_mode`;
  - `frequency_correlation_rho`;
  - random spectrum base and innovation seeds;
  - frequency index and generation rule.

Important scope boundary:
- These frequency-correlation options are stochastic-channel generation diagnostics.
- They are not new SSA physics, not a sea-surface time-evolution model, and not an experimentally calibrated cross-frequency coherence law.
- They do not alter the direct path, coherent reflection formula, `P_sca` formula, or `H_f = H_direct_f + H_reflect_f` semantics.

New scatter-scale calibration script:
- `calibrate_ssa_scatter_scale_vertical.m`.
- It uses `kirchhoff_spatial` multi-seed reduced-grid statistics as a realization-based reference and compares `ssa_stat_kernel + ssa1_geometry` across candidate `surface_ssa_scatter_scale` values.
- Calibration metrics include:
  - `abs(h_reflect)` mean/std;
  - `abs(h_reflect/h_direct)` mean/std;
  - reflected spectrum RMS delta-k;
  - SSA scatter budget terms such as `E_sca_limit/E_inc`.
- Outputs:
  - `calibrate_ssa_scatter_scale_vertical_result.mat`;
  - `calibrate_ssa_scatter_scale_summary.png`;
  - `run_table`, `summary_table`, `calibration_table`, `candidate_table`, `global_scale_table`, and `validation_report`.
- Interpretation:
  - `surface_ssa_scatter_scale` remains an engineering/numerical normalization parameter;
  - it is not an experimentally calibrated absolute scattering cross-section constant.

Scatter-scale smoke validation:
- Command overrides used:
  - `SSA_SCALE_CAL_HS_LIST='0.05 0.2'`;
  - `SSA_SCALE_CAL_F_LIST_HZ='4000 8000'`;
  - `SSA_SCALE_CAL_SEED_COUNT='2'`;
  - `SSA_SCALE_CAL_SCALE_VALUES='0 1 4'`;
  - `SSA_SCALE_CAL_GRID_N='64'`.
- The script completed successfully.
- Recommended condition-wise scales in the smoke run:
  - `Hs=0.05`, `4000 Hz`: scale `0`;
  - `Hs=0.05`, `8000 Hz`: scale `1`;
  - `Hs=0.2`, `4000 Hz`: scale `1`;
  - `Hs=0.2`, `8000 Hz`: scale `1`.
- Global candidate objective in the smoke run:
  - scale `0`: mean objective `0.31582`;
  - scale `1`: mean objective `0.085728`;
  - scale `4`: mean objective `0.085728`.
- Smoke-run global recommendation: scale `1`.
- Because this was a two-seed smoke run, it should not be treated as a final physical calibration.

New frequency-correlation validation script:
- `validate_ssa_frequency_correlation_vertical.m`.
- It compares:
  - `'independent'`;
  - `'shared_seed_phase'`;
  - `'ar1_frequency'`.
- It computes random reflected-channel residuals by subtracting a `surface_ssa_random_scatter=false` coherent reference from random-scatter runs.
- Validation metrics include:
  - adjacent-frequency random-reflection coherence;
  - repeatability with fixed seed;
  - direct-path stability;
  - `H_f = H_direct_f + H_reflect_f` invariant error.
- Outputs:
  - `validate_ssa_frequency_correlation_vertical_result.mat`;
  - `validate_ssa_frequency_correlation_summary.png`;
  - `run_table`, `summary_table`, `trend_table`, and `validation_report`.

Frequency-correlation smoke validation:
- Command overrides used:
  - `SSA_FREQ_CORR_SEED_COUNT='2'`;
  - `SSA_FREQ_CORR_GRID_N='64'`;
  - `SSA_FREQ_CORR_RHO='0.85'`.
- The script completed successfully.
- Adjacent-frequency random-reflection coherence means:
  - independent: `0.24155`;
  - shared seed phase: `0.99993`;
  - AR(1): `0.78367`.
- Maximum channel invariant error: `7.1524e-18`.
- Repeatability error for fixed seed: `0`.
- Direct-path difference across correlation modes: `0`.
- Interpretation:
  - correlated modes increase random reflected-frequency continuity relative to the independent baseline;
  - only the random reflected component changes;
  - the channel interface and direct path remain stable.

New applicability summary script:
- `summarize_ssa_applicability_vertical.m`.
- It only reads existing compact result files; it does not run PE/WAPE or communication simulations.
- It combines available diagnostics from:
  - SSA1/SSA2 coherent comparison;
  - periodic vs zero-padded convolution comparison;
  - scatter-scale calibration;
  - frequency-correlation validation.
- Outputs:
  - `summarize_ssa_applicability_vertical_result.mat`;
  - `summarize_ssa_applicability_vertical_table.csv`.
- The table is an engineering guide with columns for:
  - recommended kernel;
  - recommended scale when available;
  - SSA2 coherent-loss attention flag;
  - padding-difference attention flag;
  - wideband frequency-correlation risk note.
- It is not an experimental validity map and not a full SSA/NLSSA applicability proof.

Current recommendation:
- Keep `kirchhoff_spatial` as the default realization-based surface model.
- Keep `ssa_stat_kernel + ssa1_geometry` as the current PM-spectrum-driven first-order pressure-release / Dirichlet statistical reflection/scattering research branch.
- Keep `surface_ssa_coherent_order='ssa1'` by default; use `ssa2_broschat_coherent` as a coherent-loss sensitivity diagnostic.
- Use `surface_ssa_frequency_correlation_mode='independent'` by default for backward compatibility.
- Treat `shared_seed_phase` and `ar1_frequency` as optional wideband random-channel diagnostics until a physical or empirical frequency-correlation model is established.

## 2026-06-24 Incident and Surface-Reflected Wavefield Visualization

Purpose:
- Show the shape of the incident wave approaching the sea surface and the reflected/scattered wave leaving the surface.
- Keep the diagnostic separate from channel generation and communication processing.

Interface:
- New disabled-by-default configuration:
  - `surface_wavefield_diagnostics=false`;
  - `surface_wavefield_slice_axis='x'`;
  - `surface_wavefield_max_z_samples=256`.
- New stable output field:
  - `output.surface_wavefield_meta`.
- When disabled, the field remains present with `enabled=false` and empty arrays.
- When enabled, only the reference-frequency center slice is sampled. Full three-dimensional propagation volumes are not saved.

Recorded data:
- incident complex-envelope slice for `z_tx -> 0`;
- reflected complex-envelope slice for `0 -> z_rx`;
- sea-surface incident and reflected complex fields;
- transverse and depth coordinates;
- reference frequency, `k0`, boundary model, normalization amplitude;
- endpoint consistency errors against the existing surface and receiver fields;
- carrier-reconstruction formulas and their visualization-only limitation.

Core `output.surface_wavefield_meta` fields:
- state and coordinates:
  - `enabled`, `reference_frequency_hz`, `k0_rad_per_m`, `omega_rad_per_s`, `c0_m_per_s`;
  - `slice_axis`, `fixed_coordinate_m`, `transverse_coordinate_m`;
  - `incident_z_m`, `reflected_z_m`, `x_m`, `y_m`;
- complex fields:
  - `incident_field_slice`, `reflected_field_slice`;
  - `surface_incident_xy`, `surface_reflected_xy`;
- interpretation and audit:
  - `surface_boundary_model`, `ssa_kernel_mode`, `phase_mode`;
  - `normalization_amplitude`;
  - `endpoint_consistency`;
  - `carrier_reconstruction`.

Visualization script:
- `visualize_surface_reflection_wavefield_vertical.m`.
- Default demonstration:
  - `ssa_stat_kernel + ssa1_geometry`;
  - `f=6000 Hz`;
  - `Hs=0.2 m`;
  - `128 x 128` transverse grid;
  - `zero_padded` SSA convolution.
- The model can be changed with `SURFACE_WAVEFIELD_MODEL`, including `kirchhoff_spatial`.
- Outputs:
  - `surface_wavefield_xz_comparison.png`;
  - `surface_wavefield_xy_comparison.png`;
  - `surface_wavefield_spectrum_comparison.png`;
  - `surface_wavefield_instantaneous_pressure.mp4`;
  - compact `surface_wavefield_visualization_result.mat`.

Physical interpretation:
- Static `x-z` plots show PE/WAPE complex-envelope magnitude, not instantaneous pressure.
- Surface-plane phase plots show the complex phase of the incident and reflected envelopes.
- Spectrum plots show discrete angular-spectrum redistribution. They are not calibrated bistatic scattering cross sections.
- The animation reconstructs one nominal carrier period using:
  - incident: `real(Psi_inc*exp(i*k0*(z_tx-z)-i*omega*t))`;
  - reflected: `real(Psi_ref*exp(i*k0*z-i*omega*t))`.
- The pressure-release phase reversal is already included in `Psi_ref`.
- The carrier restoration uses nominal `c0` and is not a separate time-domain propagation solution.

Reduced validation:
- Script: `validate_surface_wavefield_visualization_vertical.m`.
- Settings:
  - `nx=ny=64`;
  - `z_tx=20 m`, `z_rx=2 m`;
  - `f=6000 Hz`, plus a `[4000,6000,8000] Hz` wideband invariant case;
  - `ssa1_geometry` and `kirchhoff_spatial`.
- Results:
  - diagnostics-on versus diagnostics-off differences in `H_f`, `H_direct_f`, and `H_reflect_f`: `0`;
  - flat pressure-release phase-reversal relative error: `2.6158e-16`;
  - flat reflected-magnitude relative error: `2.0828e-16`;
  - maximum wideband channel invariant error: `1.5516e-17`;
  - incident-surface, reflected-surface, and reflected-receiver endpoint errors: `0`;
  - all validation checks passed.

Generated-example diagnostics:
- Example settings:
  - `ssa_stat_kernel + ssa1_geometry`;
  - `f=6000 Hz`, `Hs=0.2 m`;
  - `nx=ny=128`;
  - `surface_ssa_conv_padding='zero_padded'`;
  - `sea_seed=12345`.
- Angular-spectrum results:
  - incident RMS transverse wavenumber: `3.3332019367144 rad/m`;
  - reflected RMS transverse wavenumber: `3.36650468954104 rad/m`;
  - RMS increase: `0.0333027528266383 rad/m`;
  - incident 90% energy radius: `5.06133446181766 rad/m`;
  - reflected 90% energy radius: `5.11307005973656 rad/m`;
  - 90% energy-radius increase: `0.0517355979189071 rad/m`.
- These values describe one reduced-grid realization. They are visualization diagnostics, not universal rough-surface scattering coefficients.
- Static PNG generation completed and `surface_wavefield_instantaneous_pressure.mp4` was written successfully.

Regression coverage after integration:
- Existing `validate_ssa_stat_kernel_vertical` completed successfully:
  - all Hs=0, energy, seed, scatter-scale, Dirichlet-boundary, zero-padding, and dense-vs-FFT checks passed;
  - maximum `W_eta` variance relative error: `8.6736e-15`;
  - dense periodic sum versus FFT raw-scatter relative difference: `2.5247e-15`;
  - maximum reduced wideband invariant error: `2.1346e-18`;
  - reported energy-conservation errors remained `0`.
- A one-condition reduced wideband communication smoke run completed through the existing `H_f -> H_baseband -> h_bb -> peak_sync/MMSE` consumer with diagnostics left at their default disabled state.
- `git diff --check` passed.

## 2026-06-30 Specular vs Incoherent SSA Reflection Contribution Diagnostics

Purpose:
- Quantify whether the SSA1 incoherent scatter component can be ignored under selected sea-surface and frequency settings.
- Keep the default `kirchhoff_spatial` model, WAPE marching, communication scripts, BER/SER logic, and `H_f = H_direct_f + H_reflect_f` semantics unchanged.

New disabled-by-default interface:
- `vertical_channel_model.m` validates `surface_ssa_component_diagnostics=false` by default.
- When enabled with `surface_boundary_model='ssa_stat_kernel'`, `vertical_wape_propagator.m` decomposes the sea-surface reflected field as:
  - `Psi_coh = R_coh * Psi_inc`;
  - `Psi_sca = Psi_ref - Psi_coh`.
- It then marches `Psi_coh` and `Psi_sca` separately from the sea surface to the receiver and stores the result in `output.surface_ssa_component_meta`.
- The original reflected field and public channel fields are unchanged:
  - `H_reflect_f` remains the result of the full reflected field;
  - `H_f = H_direct_f + H_reflect_f` remains the only channel composition used by downstream scripts.

Core diagnostic fields:
- Receiver components:
  - `h_reflect_coh_f`;
  - `h_reflect_sca_f`;
  - `h_reflect_total_f`;
  - `scatter_to_coherent_rx_db_f`;
  - `component_sum_error_rel_f`.
- Surface spectral energy audit:
  - `surface_E_inc_f`;
  - `surface_E_coh_f`;
  - `surface_E_sca_f`;
  - `surface_E_sca_limited_f`;
  - `surface_E_ref_f`;
  - `scatter_energy_fraction_db_f`;
  - `energy_conservation_error_f`;
  - `W_eta_variance_rel_error_f`.

New comparison script:
- `compare_specular_incoherent_surface_reflection_vertical.m`.
- Default reduced scan:
  - wind speed `[3,5,8,12] m/s`;
  - frequency `[4000,6000,8000,10000] Hz`;
  - fixed `Hs_target=0.5 m`;
  - `seed=12345+(0:7)`;
  - `nx=ny=128`;
  - `ssa_stat_kernel + ssa1_geometry`;
  - `surface_ssa_random_scatter=true`.
- Environment variables can reduce the run for smoke tests:
  - `SPEC_INCOH_WIND_LIST`;
  - `SPEC_INCOH_F_LIST_HZ`;
  - `SPEC_INCOH_SEED_COUNT`;
  - `SPEC_INCOH_GRID_N`;
  - `SPEC_INCOH_HS_TARGET`;
  - `SPEC_INCOH_SCATTER_SCALE`;
  - `SPEC_INCOH_CONV_PADDING`;
  - `SPEC_INCOH_RESULT_FILE`.

Outputs:
- `compare_specular_incoherent_surface_reflection_vertical_result.mat`;
- `compare_specular_incoherent_surface_reflection_vertical_summary.csv`;
- `run_table`, `summary_table`, `decision_table`, and `validation_report`;
- PNG figures:
  - `specular_incoherent_Ecoh_Esca_vs_frequency.png`;
  - `specular_incoherent_rx_ratio_vs_frequency.png`;
  - `specular_incoherent_rx_abs_vs_frequency.png`;
  - `specular_incoherent_decision_heatmap.png`;
  - `specular_incoherent_energy_limit_heatmap.png`.

Decision rule:
- `negligible`: both the receiver scatter/coherent ratio and surface scatter/coherent energy ratio are below `-20 dB`.
- `borderline`: between `-20 dB` and `-10 dB`.
- `not_negligible`: above `-10 dB`, or the scatter/coherent receiver ratio has large seed-to-seed spread.
- This is an engineering diagnostic for the current grid, `Hs_target`, `surface_ssa_scatter_scale`, seed set, and SSA1 Dirichlet model. It is not an experimentally calibrated scattering threshold.

Interpretation boundary:
- The default sweep fixes `Hs_target`, so changing `sea_wind_speed` changes the normalized PM spectral shape, not total roughness variance.
- `surface_ssa_scatter_scale` remains an engineering normalization parameter, not an absolute scattering cross section.
- Scatter-only receiver propagation is a diagnostic decomposition; it is not a new channel model and does not replace the full reflected channel.

Smoke validation:
- Command overrides used:
  - `SPEC_INCOH_WIND_LIST='3 8'`;
  - `SPEC_INCOH_F_LIST_HZ='4000 8000'`;
  - `SPEC_INCOH_SEED_COUNT='2'`;
  - `SPEC_INCOH_GRID_N='64'`;
  - `SPEC_INCOH_RESULT_FILE='compare_specular_incoherent_surface_reflection_vertical_smoke_result.mat'`.
- The script completed successfully and generated all five PNG figures.
- Diagnostics-disabled versus diagnostics-enabled differences:
  - `H_f`: `0`;
  - `H_direct_f`: `0`;
  - `H_reflect_f`: `0`.
- Validation maxima:
  - channel invariant error: `6.9929e-18`;
  - receiver component sum relative error: `1.3409e-15`;
  - energy-conservation error: `1.6554e-16`;
  - `W_eta` variance relative error: `2.2204e-15`.
- Smoke decision table:
  - `U=3 m/s`, `4000 Hz`: scatter/coherent receiver ratio `+74.55 dB`, scatter/coherent energy ratio `+76.20 dB`, decision `not_negligible`;
  - `U=3 m/s`, `8000 Hz`: scatter/coherent receiver ratio `+282.63 dB`, scatter/coherent energy ratio `+197.63 dB`, decision `not_negligible`;
  - `U=8 m/s`, `4000 Hz`: scatter/coherent receiver ratio `+75.10 dB`, scatter/coherent energy ratio `+76.20 dB`, decision `not_negligible`;
  - `U=8 m/s`, `8000 Hz`: scatter/coherent receiver ratio `+283.82 dB`, scatter/coherent energy ratio `+197.63 dB`, decision `not_negligible`.
- Interpretation:
  - For this smoke case, `Hs_target=0.5 m` and `surface_ssa_scatter_scale=1` cause very strong SSA coherent mirror loss, so the coherent receiver component is extremely small and the scatter-only reflected contribution dominates.
  - This result supports retaining the incoherent term for the tested roughness/frequency settings; weaker `Hs_target` and calibrated smaller `surface_ssa_scatter_scale` should be checked before making a broader simplification.

## Surface-Only Plane-Wave Coherent Reflection Diagnostic

Purpose:
- `compare_plane_wave_surface_coherent_reflection_vertical.m` isolates the sea-surface boundary response from PE/WAPE propagation.
- The incident field is a vertical plane wave; on the sea-surface plane this is represented by a constant transverse envelope:
  - `Psi_inc(x,y)=1`.
- The diagnostic compares coherent specular reflection strength, not full channel gain, communication BER/SER, or receiver propagation.

Raw PM roughness:
- This script intentionally does not apply `Hs_target` normalization.
- For each wind speed, the PM spectrum is integrated directly:
  - `sigma_eta_raw^2 = sum(Phi2D(:))*dkx*dky`;
  - `Hs_raw = 4*sigma_eta_raw`.
- Therefore wind speed changes both PM spectral shape and integrated roughness variance in this diagnostic.

Compared coherent-reflection quantities:
- SSA1 pressure-release normal-incidence coherent reflection:
  - `R_SSA1 = -exp(-2*k0^2*sigma_eta_raw^2)`.
- Kirchhoff realization coherent reflection:
  - for each random surface realization, the script forms a phase screen `G(x,y)=R0*exp(i*delta_phi(x,y))`;
  - the coherent mirror component for one finite surface is the spatial mean `mean(G(:))`;
  - the ensemble coherent estimate is `abs(mean(R_sample))` across realizations, where `R_sample=mean(G(:))`;
  - `mean(abs(R_sample))` is also reported as a finite-aperture residual, but it is not the strict ensemble coherent coefficient.
- Two Kirchhoff phase conventions are reported:
  - `kirchhoff_corrected`: `delta_phi=2*k0*eta` for normal incidence, the corrected two-way height phase convention whose Gaussian coherent expectation matches SSA1;
  - `kirchhoff_legacy_4k_eta`: `delta_phi=2*k0*2*eta`, the pre-fix diagnostic convention that produces stronger `exp(-8*k0^2*sigma_eta^2)`-type coherent attenuation.
- The script also reports `mean(abs(G(:)).^2)`, which should stay near 1 for a pure pressure-release phase screen. This distinguishes total reflected surface power from coherent specular strength.
- The script also records realization variance `std(eta)^2` against the PM spectrum variance and writes a controlled i.i.d. Gaussian validation table. The controlled validation checks that `abs(E[-exp(i*2*k0*eta)])` follows `exp(-2*k0^2*sigma_eta^2)` within Monte Carlo tolerance.

Outputs:
- `compare_plane_wave_surface_coherent_reflection_vertical_result.mat`;
- `compare_plane_wave_surface_coherent_reflection_vertical_summary.csv`;
- `compare_plane_wave_surface_coherent_reflection_vertical_controlled_gaussian.csv`;
- `plane_wave_coherent_absR_vs_wind.png`;
- `plane_wave_coherent_power_vs_wind.png`;
- `plane_wave_coherent_loss_vs_wind.png`;
- `plane_wave_kirchhoff_vs_ssa_absR_heatmap.png`;
- `plane_wave_raw_pm_sigma_vs_wind.png`.

Validation checks:
- flat synthetic case gives `|R_SSA1|=1` and `|R_Kirchhoff|=1`;
- PM variance identity is checked from the discrete spectrum;
- Kirchhoff phase-screen total reflected power is checked against 1;
- all checks are recorded in `validation_report`.

## 2026-07-02 Kirchhoff K-Stat Statistical Phase-Screen Branch

Changed files:
- `vertical_channel_model.m`: accepts `surface_boundary_model='kirchhoff_kstat'` and validates `surface_kstat_random_scatter`, `surface_kstat_seed_offset`, `surface_kstat_conv_padding`, and `surface_kstat_trusted_angle_deg`.
- `vertical_wape_propagator.m`: passes the kstat configuration to `pm_surface_boundary_model.m` and preserves disabled `roughness_meta.kirchhoff_kstat_meta` when surface reflection is off.
- `pm_surface_boundary_model.m`: implements `local_apply_kirchhoff_kstat` without generating a concrete `eta(x,y)`.
- `scripts/validation/validate_kirchhoff_kstat_vertical.m`: validates flat-surface degeneration, coherent formula, phase-screen energy closure, propagation-window diagnostics, deterministic seeds, explicit Kirchhoff phase-screen mean comparison, and weak SSA1 coherent-reference consistency.

Interface:
- Select the branch with `paramsV.surface_boundary_model='kirchhoff_kstat'`.
- `surface_kstat_random_scatter=true` synthesizes `Psi_sca_k=sqrt(P_sca).*Z`, `Z~CN(0,1)`.
- `surface_kstat_random_scatter=false` returns the coherent reflected field only, while still recording `P_sca` and all energy diagnostics.
- `surface_kstat_seed_offset=200000` gives `seed_kstat=sea_seed+surface_kstat_seed_offset+frequency_index-1`.
- `surface_kstat_conv_padding` supports `periodic` and `zero_padded`, reusing the existing FFT convolution helper.
- `surface_kstat_trusted_angle_deg=NaN` disables the optional trusted-angle window; finite values record `K_h<=k0*sind(angle)` energy without renormalizing it.

Implemented formulas:
- PM height spectrum is converted to the continuous convention used by this branch so that
  `C_eta(0)=sum(W_eta(:))*dkx*dky/(2*pi)^2=(Hs_target/4)^2`.
- `C_eta_xy = real(ifft2(W_eta))*numel(W_eta)*dkx*dky/(2*pi)^2`.
- Near-vertical first version uses `alpha=2*k0`.
- `<G>=exp(-0.5*alpha^2*sigma_eta^2)` and `R_coh=R0*<G>`. For pressure-release `R0=-1`, this is `-exp(-2*k0^2*sigma_eta^2)`.
- `C_deltaG(rho)=exp(-alpha^2*sigma_eta^2)*(exp(alpha^2*C_eta(rho))-1)`, implemented in the algebraically equivalent stable form `exp(alpha^2*(C_eta-sigma_eta^2))-exp(-alpha^2*sigma_eta^2)`.
- `S_deltaG=fft2(C_deltaG)*dx*dy`. Only `S_deltaG` is used for incoherent scatter; the coherent delta spike from total `S_G` is not used.
- `P_sca(Ks)=|R0|^2*circconv(S_deltaG,abs(Psi_inc_k).^2)*dkx*dky/(2*pi)^2`.

Metadata:
- `output.roughness_meta.kirchhoff_kstat_meta` records `sigma_eta2_m2`, `alpha_rad_per_m`, `G_mean`, `R_coh`, `S_deltaG`, `P_sca`, seed fields, FFT normalization notes, negative-spectrum clipping diagnostics, and spectrum moment summaries.
- Full-grid phase-screen energy is checked as `abs(G_mean)^2 + sum(S_deltaG(:))*dkx*dky/(2*pi)^2`.
- Propagating-window energy `K_h<=k0` and optional trusted-angle energy are recorded separately and are not normalized to `1-abs(G_mean)^2`.

Validation summary:
- Reduced kstat validation with `KSTAT_VALIDATE_ENSEMBLE_COUNT=2` completed successfully.
- `Hs=0` flat kstat versus `kirchhoff_spatial` maximum `H_f` difference: `7.7579e-17`.
- Coherent formula check `abs(R_coh - R0*exp(-2*k0^2*sigma_eta^2))`: `0`.
- Full-grid phase-screen energy closure error for the reduced rough case: `6.6613e-16`.
- `H_f = H_direct_f + H_reflect_f` invariant for the rough kstat case: `2.6026e-18`.
- Fixed seed repeat difference in `H_f`: `0`; changed seed changed reflected response by `0.011905` while direct path drift stayed `0`.
- Explicit Kirchhoff phase-screen mean comparison relative error with 2 seeds: `0.26709`, recorded as a finite-Monte-Carlo trend check.
- Weak roughness coherent reference against `ssa_stat_kernel + ssa1_geometry`: difference `0`.

Interpretation boundary:
- `kirchhoff_kstat` is a Kirchhoff / Gaussian statistical phase-screen model, not SSA/NLSSA, not a T-matrix, and not a calibrated scattering cross-section model.

## 2026-07-02 Raw-PM Kirchhoff Wind Comparison

Changed files:

- `vertical_channel_model.m`: adds `paramsV.surface_roughness_scale_mode`.
- `vertical_wape_propagator.m`: forwards the roughness scale mode into `pm_surface_boundary_model.m`.
- `pm_surface_boundary_model.m`: supports `roughness_scale_mode='target_hs'` and `roughness_scale_mode='raw_pm'`.
- `scripts/comparisons/compare_kirchhoff_kdomain_kstat_wind_vertical.m`: compares `kirchhoff_kdomain` and `kirchhoff_kstat` across wind speeds using raw PM roughness.

Interface:

- `surface_roughness_scale_mode='target_hs'` remains the default and preserves prior behavior: PM spectra or realizations are scaled to `sea_hs_target`.
- `surface_roughness_scale_mode='raw_pm'` disables `sea_hs_target` amplitude scaling. `sea_hs_target` is retained in config/metadata but does not set the surface variance.
- Roughness metadata records `roughness_scale_mode`, `sigma_eta_raw_m`, `Hs_raw_m`, `Hs_target_m`, `scale_factor`, and raw PM variance audits.

Formula and convention notes:

- Explicit `kirchhoff_kdomain` uses the existing discrete PM realization convention, where the raw realization is generated from `Phi2D*dkx*dky` and `Hs_raw=4*std(eta(:))`.
- Raw-PM wind comparisons use the project discrete PM variance as the main convention: `sigma_eta_raw^2=sum(Phi2D(:))*dkx*dky`.
- `kirchhoff_kstat` still evaluates `C_eta(0)=sum(W_eta(:))*dkx*dky/(2*pi)^2` internally, so in `raw_pm` mode it sets `W_eta=Phi2D*(2*pi)^2` and records `scale_factor=2*pi`.
- The metadata includes both `pm_variance_raw_discrete_m2` and `pm_variance_raw_continuous_m2`; `sigma_eta_raw_m` and `Hs_raw_m` use the discrete convention for the two Kirchhoff branches.

Comparison script defaults:

- `wind_list=[3,5,8,10,12,15]` m/s.
- `seed_count=32`, `seed_list=12345+(0:31)`.
- `f0=6000` Hz, `nx=ny=128`, `save_mode='rx_only'`, `show_figures=false`.
- Environment overrides: `KIRCH_WIND_LIST`, `KIRCH_SEED_COUNT`, `KIRCH_SEED_LIST`, `KIRCH_F0_HZ`, `KIRCH_GRID_N`, `KIRCH_NX`, `KIRCH_NY`, `KIRCH_XW_M`, and `KIRCH_YW_M`.

Outputs:

- `results/comparisons/compare_kirchhoff_kdomain_kstat_wind_vertical_result.mat`.
- `results/comparisons/compare_kirchhoff_kdomain_kstat_wind_vertical_summary.csv`.
- `kirchhoff_raw_pm_Hs_vs_wind.png`, `kirchhoff_coherent_R_vs_wind.png`, `kirchhoff_abs_h_reflect_vs_wind.png`, `kirchhoff_incoherent_energy_vs_wind.png`, and `kirchhoff_propagating_energy_vs_wind.png`.

Validation intent:

- The script checks `H_f = H_direct_f + H_reflect_f` for both branches.
- For `kirchhoff_kstat`, it records phase-screen full-K energy closure and propagating-window incoherent energy without renormalizing the propagating window.
- The comparison is a raw-PM wind-speed diagnostic. It should not be mixed with older sweeps where wind speed changed PM spectral shape after fixed-`Hs_target` rescaling.
- Its full phase-screen energy is internally conserved by the coherent plus `S_deltaG` split.
- SSA1 is retained as a weak-roughness/small-angle reference check; it is not the main generation formula for this branch.
- The first version uses near-vertical `alpha=2*k0`; oblique `gamma_i+gamma_s` phase statistics remain future work.

## 2026-07-11 Raw-PM Grid and Joint-Frequency K-Stat Prerequisite Validation

New independent prototype interfaces:

- `raw_pm_spectrum_grid_vertical.m` reproduces the project PM spectrum and records analytic infinite-domain variance, `K_min`, `K_peak`, `K_nyquist`, implied `Hs`, discrete/infinite variance ratio, and an idealized radial-support ratio.
- `sample_raw_pm_surface_vertical.m` reproduces the explicit raw-PM surface convention for mapping tests.
- `sample_kirchhoff_kstat_joint_frequency_vertical.m` generates zero-mean joint-frequency Kirchhoff phase screens from the analytic cross-frequency covariance and pseudo-covariance. It uses augmented `K,-K` spectral pairs; `mode='independent'` is retained only as an engineering comparison.
- These interfaces do not change `vertical_channel_model`, public defaults, the existing independent kstat branch, or the communication chain.

Validation scripts and results:

- `validate_raw_pm_grid_coverage_vertical.m` scans wind, aperture, and resolution. At fixed `dx=50/128 m`, U=5 is well resolved by a 100 m/256^2 PM grid (`variance ratio=0.9983`). U=15 requires approximately 400 m/1024^2 (`variance ratio=0.9913`, idealized support ratio `0.9967`).
- A same-dx central crop from U=5, 100 m/256^2 PM to 50 m/128^2 PE coordinates retained ensemble mean energy within `2.04%` over 32 seeds. This supports PM/PE grid separation for the prototype, but larger seed confidence intervals are still required.
- `validate_kstat_joint_frequency_boundary_vertical.m` compares one shared explicit Gaussian surface, independent kstat, and covariance+pseudocovariance joint kstat at U=5, 4--8 kHz, F=32, 64^2, and 64 realizations.
- Joint kstat achieved cross-frequency covariance relative error `0.0067`, adjacent-correlation RMSE `2.15e-4`, boundary PDP correlation `0.9837`, and boundary LFM matched-filter correlation `0.9842`. Independent kstat gave `0.9526`, `0.9899`, `0.4890`, and `0.5163`, respectively.
- The analytic U=5 zero-lag pseudo/covariance Frobenius ratio was `2.09e-10`; finite explicit/joint sample ratios were `0.0027/0.0034`. Pseudo-covariance remains implemented because this near-proper result is condition-specific.
- `audit_pe_caching_vertical.m` shows that fixed source/environment/frequency incident fields, direct response, coherent reflection, propagation kernels, and medium terms are cacheable. For F=32 and L=64, PE step and FFT counts are projected to fall by about `42.9x` and `41.1x` when only random surface-to-receiver propagation is repeated.

Current boundary:

- The joint model is a boundary-only prototype and is not yet connected to PE receiver outputs.
- Full joint storage scales as `F^2*Nxy`; large PM grids need frequency low rank, K-domain blocking, or a low/high-wavenumber decomposition before production use.
- The next stage is one U=5 cached joint-kstat + PE receiver validation, not a multi-wind channel library.
- Detailed report: `reports/raw_pm_joint_kstat_prerequisite_validation_report.md`.

## 2026-07-11 Cached Joint-Kstat + PE Receiver Validation

The validation-only cached executor now stores the fixed incident surface
field, direct receiver response, coherent reflected receiver response,
surface-to-receiver PE factors/screens, and the same-dx PM-to-PE crop. The
public propagator and defaults are unchanged.

At U=5 m/s, raw PM, 4--8 kHz, F=32, PM 100 m/256^2, PE 50 m/128^2,
Ltrain=128, and Ltest=64, reflected-only joint kstat achieved covariance
relative error 0.2704 versus 0.9750 for independent kstat. PDP and LFM
correlations against explicit same-surface kdomain were 0.9714 and 0.9917.
The explicit kdomain train/test covariance split itself differed by 0.2610,
so the joint result is near the present sampling floor.

Double same-input cached/public validation gave a maximum total-response
difference of 7.04e-16. Public F=32 time was 43.71 s versus 0.233 s for one
cached surface-to-receiver propagation. Joint generation plus cached PE used
74.14 s for L=128. Cache/model build times were 10.92/32.85 s; stored cache
and factor arrays were 20/360 MiB, with a 3.018 GiB MATLAB memory snapshot.

Seven of eight receiver-entry criteria pass. The unresolved item is whether
receiver pseudo-covariance can be ignored: finite L=64 sample P/C ratios are
about 0.26 for both kdomain and joint, and require a proper-null finite-sample
calibration before interpretation. F=32 also places the 99% scatter tail near
the 7.75 ms ambiguity window. Detailed report:
`reports/cached_joint_kstat_pe_receiver_validation_report.md`.

## 2026-07-12 U=5 Conditional Receiver Generator, F=64

Independent modules now estimate one-condition receiver `mu/C/P`, sample
proper or augmented-real improper channels, and build physical CIRs without
changing the public propagator or communication chain. The validated model is
U=5 m/s raw PM, 4--8 kHz, F=64, PM 100 m/256^2, PE 50 m/128^2,
Ltrain=128, and held-out Ltest=64.

Joint receiver PDP/LFM correlations against explicit kdomain were
0.9525/0.9839; independent-frequency results were 0.0048/0.4104. The shortest
circular 99% energy interval was 4.657 ms versus Tmax=15.75 ms. Joint
properness was not rejected (observed P/C=0.1911, central 95% null interval
[0.1364,0.3476], p_upper=0.7326). Explicit kdomain stayed inside its central
95% interval but has a one-sided p=0.0405 warning; the improper sampler remains
available.

Held-out selection chose full numerical rank 34. Full/99.9%/99% ranks were
34/7/5 with covariance errors 0.3922/0.4000/0.3992 and PDP correlations
0.9548/0.9501/0.9533. Ten thousand H+CIR samples required 0.178 s without PE;
MAT saving required another 1.007 s. The receiver model is about 694 KB.

The expensive training frontend remains the main risk: F=64 joint factors
took 631 s, total model construction 840 s, peak MATLAB memory 6.626 GiB, and
the reusable spatial cache is about 815 MB. Detailed report:
`reports/u5_conditional_channel_generator_f64_report.md`.

## 2026-07-12 Streaming Joint Builder and U=8 Aperture

The F=64 joint builder now uses a converged covariance/pseudo-covariance
series, scalar spatial FFT modes, a global augmented-frequency basis, and
blockwise local EVDs. It never allocates an `F x F x N_K` tensor. At U=5 its
factor build fell from 631.3 s to 6.50 s, the MATLAB memory snapshot from
6.626 to 1.664 GiB, and disk storage from about 815 MB to 109.6 MiB. Receiver
PDP/LFM correlations against the original builder were 0.9986/0.9993.
Short-batch expanded-factor sampling measured 0.948x old throughput, so
sampling parity remains a small engineering issue even though construction
targets pass.

The U=8 raw-PM aperture audit used 128 seeds at common PE spacing. The chosen
grid is 150 m/384^2: implied Hs=1.3636 m, discrete/infinite capture=0.99780,
Kpeak/Kmin=2.570, and central-crop energy ratio 0.9579 with 95% CI
[0.9143,1.0015]. The 100 m grid undersamples below the PM peak; moving from
150 to 200 m changes implied Hs by only 0.146%.

The U=8 full 128/128 cached receiver validation now passes. Joint covariance
error is 0.3554, close to the kdomain split floor 0.3403 and far below the
independent result 0.9703. Joint reflected-only PDP/LFM correlations are
0.9809/0.9911, and T99/Tmax is 0.2910. U=8 properness is not rejected
(observed P/C=0.2351, central 95% null interval [0.1770,0.2802], upper-tail
p=0.3318). Full receiver rank 42 remains the physical default; ranks 14/11
are compression candidates.

Ten thousand U=8 H+CIR samples take 0.0154 s before file saving and make zero
PE calls. The exact-node U=5/U=8 library is built and smoke-tested; unsupported
U=6 is rejected and component sums are exact. U=8 builder-only memory is
3.062 GiB, but expanded runtime sampling reaches a 6.858 GiB snapshot and the
384^2 cache is 1.060 GB, so runtime factor streaming/grouping is required
before higher-wind production. See
`reports/u8_joint_optimization_two_node_library_report.md`.

## 2026-07-13 Two-Node Communication Validation

The existing MPSK entry now has an optional, disabled-by-default external
H(f)/h(t) MAT-file input. The default direct-only/direct-plus-reflect PE path
and receiver-side noise placement are unchanged. Independent reusable code
converts H(f) to symbol taps and evaluates channel ensembles with common bits,
noise seeds, peak synchronization, known-channel MMSE, Wilson intervals, and
channel-cluster bootstrap intervals.

At U=5 and U=8, 32 fresh channels per source and 8000 QPSK symbols per
channel were tested over 0:2:20 dB. Statistical-full versus kdomain BER-curve
log10 RMSE was 0.127/0.194 decades; joint versus kdomain was 0.242/0.086.
All BER and SER cluster-95% intervals overlapped at every Eb/N0 point. Mean
BER increased from U=5 to U=8 by 1.93x for kdomain, 1.35x for joint, and 1.69x
for statistical-full channels, so the wind trend is consistent.

The H-to-tap audit found that the legacy unshifted IFFT prefix retained all
2048 zero-padded samples and could discard 61.8% pre-peak energy. The
validation path now unwraps the shortest circular interval containing 99.9%
energy. It preserves at least 0.999 energy, has about 1e-15 FFT closure, and
keeps the distinction between 0.25 ms physical resolution and 1 ms symbol
tap spacing.

A 128-pair shared-latent test found no significant U=8 BER change between
full rank 42 and 99.9% rank 14. U=5 rank 7 differed only at 20 dB by
1.08e-4 BER. Full rank remains the reference; 99.9% is an acceptable
compression candidate.

The present frequency-domain MMSE has a noiseless BER floor of roughly
0.2--0.5%, caused by finite linear-convolution boundaries, circular FFT
equalization, tap truncation, and fixed regularization. The two-node library
is approved for comparative communication research, but the next priority is
a zero-error noiseless block/equalizer baseline before adding U=10. Detailed
report: `reports/two_node_statistical_channel_communication_validation_report.md`.

## 2026-07-16 Exact Discrete Adjoint PE Receiver Projection

An independent validation-only prototype now implements the exact discrete conjugate transpose of the cached uniform surface-to-receiver PE. The public `vertical_channel_model`, `vertical_wape_propagator`, default surface model, and communication entrypoint remain unchanged. The v1 scope is CPU double, uniform sound speed, fixed transmitter/grid/frequency axis, one nearest-grid receiver, no bubbles, and no Doppler.

The shared forward primitive preserves the existing `fr -> ifft2 -> screen -> fft2 -> fr` ordering. The adjoint reverses depth steps and uses `conj(fr)` and `conj(screen)` without inverse sponge gain or empirical FFT scaling. The receiver weight is `a_PE=conj(psi_inc).*q`. The PM weight is the exact central-crop adjoint `a_PM=E^H*a_PE`, and all covariance operations remain on the PM periodic grid. The current joint variable remains `deltaG=R0*exp(1i*alpha*eta)-Rcoh`; `C_deltaG/P_deltaG` already contain `R0`, so the projection weight does not multiply it again.

Full validation passed. The maximum adjoint and receiver-projection errors were `7.77e-15` and `5.10e-15`; the refactored cached runner matched its original loop exactly. PM 16^2 / PE 8^2 dense-versus-FFT errors were `6.18e-16` for covariance and `9.35e-16` for pseudo-covariance. F=9 with 4096 realizations and F=64 with 512 realizations both placed analytic/sample covariance and pseudo-covariance errors below `1.25x` their split-sample floors. Same-realization cached-forward/projection errors were `2.52e-13` and `1.07e-13`.

At F=64, PE 128^2, and PM 256^2, kernel construction took `0.497 s`, analytic `C/P` took `7.10 s`, and batched projection was about `35.2x` faster per realization than cached forward PE. The PM embedded/FFT weights require about `128 MiB`; the implementation streams frequency pairs and never stores spatial F^2 covariance blocks. The main analytic bottleneck is the F^2 frequency-pair PM FFT/contraction, not kernel construction.

Recommended roles are: adjoint projection for exact fixed-environment receiver realizations; analytic FFT contraction for direct second-order `C/P`; cached forward PE as the physical regression oracle; and the conditional statistical generator for large communication Monte Carlo after analytic/sample validation. The prototype is suitable for main-branch inclusion as an isolated research module, but not yet for the public default propagation path. See `reports/adjoint_pe_receiver_projection_feasibility_report.md`.

## 2026-07-18 PE Propagation Visual Atlas

A report-only atlas now connects the existing PE and rough-surface diagnostics into one reproducible visual narrative. `scripts/reporting/generate_pe_propagation_atlas_vertical.m` reads the accepted F=64 adjoint `C/P`, U=5/U=8 conditional models, and prior LFM validation, while directly collecting one deterministic U=5 spatial realization at 6 kHz. It stores only center x-z slices, selected x-y planes, receiver-frequency quantities, and plot metadata; no full 3-D wavefield volume is retained.

The 14 outputs cover geometry, direct and reflected PE marching, surface components and angular spectra, flat/explicit Kirchhoff/joint-kstat/SSA1 boundary fields, receiver-plane components and profiles, cached-versus-adjoint projection, reflected-only frequency response and `C/P`, PDP/eigen/distribution views, U=5/U=8 receiver statistics, an optional noiseless LFM probe, and a carrier-reconstruction MP4. Spatial comparisons share coordinates, reference amplitude, and color limits. The animation is a complex-envelope carrier reconstruction, not a time-domain PE solve.

Comparison semantics are intentionally restricted. Flat and explicit Kirchhoff fields may be compared for the same incident field; cached PE and adjoint projection use the same `deltaG` and are compared numerically. Independent explicit and joint-kstat realizations are compared only at ensemble/statistical level. Analytic FFT `C/P` and the conditional generator produce receiver statistics or channels, not spatial propagation fields. Acceptance metrics remain reflected-scatter based; total-channel plots are contextual only.

This atlas does not modify `vertical_channel_model`, `vertical_wape_propagator`, the default surface model, or `comm_main_vertical_psk`. Its scope remains uniform sound speed, CPU double, fixed grids/frequencies, one nearest-grid receiver, and no bubbles or Doppler. Outputs and provenance are under `results/visualization/pe_propagation_atlas/`.

## 2026-07-19 Stage-1 PE/Bellhop Flat-Surface Cross-Validation

`scripts/validation/validate_pe_bellhop_flat_surface_vertical.m` adds an external-model validation layer without changing the PE propagation core or public output fields. It locates Bellhop through the optional function override, `BELLHOP_EXE`, the MATLAB path, or the operating-system path; writes a standard Bellhop ASCII-arrivals `.env`; runs the executable in the result directory; parses the `.arr`; and compares it with the existing `vertical_channel_model` direct and reflected components.

The shared deterministic environment is a 100 m isovelocity water column with (c=1500\ \mathrm{m/s}), a flat pressure-release sea surface, transmitter at ([0,0,80]\ \mathrm{m}), receiver at ([5,0,10]\ \mathrm{m}), and center frequency 4 kHz. The 5 m horizontal offset avoids the degenerate zero-range ray geometry while retaining near-vertical propagation. PE uses CPU double, a (128\times128) grid over (32\times32\ \mathrm{m}), `stepz_lamb=0.5`, `sea_hs_target=0`, `surface_reflect_coeff=-1`, `surface_boundary_model='kirchhoff_spatial'`, bubbles off, and 65 equally spaced frequencies from 3--5 kHz. Bellhop uses 5001 geometric hat beams in an angle fan restricted to the direct and one-surface-bounce paths. Bottom-bounce arrivals are excluded from the stage-1 comparison.

The PE outputs are reduced envelopes, so the validation layer reconstructs physical path phase separately:

\[
H_{\rm PE}^{\rm phys}(f)=
H_{\rm dir}(f)e^{-i2\pi fL_{\rm dir}/c}
+H_{\rm ref}(f)e^{-i2\pi fL_{\rm surf}/c},
\]

with (L_{\rm dir}=\sqrt{r^2+(z_{\rm tx}-z_{\rm rx})^2}) and (L_{\rm surf}=\sqrt{r^2+(z_{\rm tx}+z_{\rm rx})^2}). A Hann-windowed, 8-times-zero-padded CIR is used to extract component/PDP peak times. Zero padding interpolates the delay grid but does not improve the 0.5 ms physical resolution set by the 2 kHz bandwidth. The reconstruction exists only in the validator; the public PE arrays and `build_physical_cir_vertical` convention are unchanged.

The accepted run found two main paths in both models. Direct arrival was 46.785563 ms in both to reported precision; the single surface arrival was 60.016332 ms for PE and 60.092516 ms for Bellhop, a 0.076184 ms difference. Maximum reconstructed-PDP peak difference was 0.061538 ms. Raw PE TL was 33.9582/35.7698 dB for direct/reflected paths, while Bellhop was 36.9241/39.0982 dB. The approximately 3 dB common offset is consistent with the finite-width PE Gaussian initial-field normalization versus Bellhop's unit point source. After a single direct-path amplitude calibration, the reflected-path TL residual was 0.3625 dB; this calibration is reported but never written back to PE outputs.

Regression checks also ran a scalar-frequency direct-only case, a scalar flat-reflection case, and the 65-frequency reflected case. `H_f=H_direct_f+H_reflect_f` held to (1.73\times10^{-18}), disabling reflection produced exactly zero reflected response and no direct-path change, and both scalar and wideband 1/R checks passed. Acceptance limits were two paths exactly, arrival/PDP peak error no more than 0.75 ms, raw TL error no more than 6 dB, relative-path and direct-calibrated reflected TL error no more than 3 dB, and channel invariant error no more than (10^{-10}).

Artifacts are written to `results/validation/pe_bellhop_flat_surface/`: Bellhop input/raw arrivals/print output, comparison CSV, MAT result, PDP/TL PNG, and `pe_bellhop_flat_surface_report.md`. This first stage does not validate rough surfaces, Kirchhoff/SSA scattering, stochastic channels, bottom interaction, modulation, or communication performance. The remaining visible PE PDP sidelobe structure is attributable to finite-band reconstruction and frequency-dependent PE envelope amplitude/phase and must not be counted as additional Bellhop eigenrays.

## 2026-07-20 PE Phase Audit and Bellhop Geometry/Convergence Matrix

**Superseded in implementation status by the 2026-07-22 carrier-phase
unification entry below.** The carrier-sign and Bellhop timing measurements
remain valid historical evidence, but statements that carrier restoration is
validation-only or that the public/CIR/communication path is unchanged no
longer describe the current code.

Two validation-only entrypoints now qualify the 2026-07-19 result without
changing `vertical_channel_model`, `vertical_wape_propagator`,
`build_physical_cir_vertical`, or the communication chain. The phase audit
`scripts/validation/validate_pe_phase_convention_uniform_vertical.m` derives
the reduced PE operator identity

\[
\exp\!\left[-i d\frac{\kappa^2}{\sqrt{k^2-\kappa^2}+k}\right]
=\exp[i d(k_z-k)]
\]

and compares the marched direct envelope with an independent one-step
discrete angular-spectrum propagation of the same Gaussian source. Under the
validation synthesis convention `exp(-i*omega*t)`, the correct physical
carrier is `exp(+i*k*d)`. The audit suppresses sponge influence with the
smallest validation-compatible positive coefficient (`1e-14`); it does not
relax the public API validation. At receiver offsets 3/6/9 m, the maximum
operator error was `2.2741e-13`, positive-carrier phase RMS was
`3.1021e-11 rad`, PE/reference group-delay difference was
`1.2768e-12 ms`, and validation-local CIR peak differences were exactly zero
on the sampled delay grid. No-carrier and negative-carrier candidates were
worse by many orders of magnitude.

The follow-on
`scripts/validation/validate_pe_bellhop_flat_surface_matrix_vertical.m` uses
the unchanged flat pressure-release PE path with the same 100 m isovelocity
environment and an open Bellhop launch fan of -89.5 to -60 degrees with
10001 beams. Bellhop arrivals are classified after the run by bottom/top
bounce count and 0.1 ms delay clustering; the validator does not preselect
exactly two launch angles. It tests offsets 3/6/9 m and restores only the
nominal longitudinal carriers: `exp(+ik(z_tx-z_rx))` for the direct component
and `exp(+ik(z_tx+z_rx))` for the image-source component. This makes the
remaining transverse phase, and therefore the extracted path delay, a PE
prediction rather than injected full-path geometry.

Timing and algebra passed the prescribed strict checks. Across all six paths,
the maximum PE-versus-analytic/Bellhop delay difference was `0.027285 ms`;
Bellhop-versus-analytic error was `4.0135e-06 ms`. Each expected Bellhop path
formed one classified cluster. `H_f=H_direct_f+H_reflect_f` held within
`3.5762e-18`, the scalar direct-only response had exactly zero reflection,
and all 1/R checks passed.

Amplitude and one convergence criterion did not meet the deliberately strict
targets, so the matrix result is recorded as `passed=false` rather than being
retuned after inspection. One global direct-path amplitude scale left a
`1.6540 dB` scale spread across geometry (limit 1 dB); reflected-path TL
residual RMS/max were `1.6295/2.0751 dB` (limits 1/2 dB). These trends are
consistent with comparing the finite-width PE Gaussian initial field against
Bellhop's unit point source, especially at the 9 m off-axis receiver, but they
remain a measured model mismatch rather than a proven normalization constant.

The x=6 m PE convergence matrix used C0=128^2/32 m/0.5 lambda,
C1=256^2/32 m/0.5 lambda, C2=256^2/64 m/0.5 lambda,
C3=128^2/32 m/0.25 lambda, and C4=256^2/32 m/0.25 lambda. All reconstructed
arrival changes were zero on the interpolated delay grid; the maximum
direct/reflected relative-TL change was `0.14445 dB`. C1, C3, and C4 phase
RMS values were at most `0.00201 rad`, while the doubled-window C2 reflected
phase RMS was `0.22559 rad`, exceeding the 0.15 rad target. This isolates the
remaining numerical concern to transverse-window/sponge sensitivity rather
than longitudinal step size or transverse sampling density.

Reports and machine-readable outputs are under
`results/validation/pe_phase_convention_uniform/` and
`results/validation/pe_bellhop_flat_surface_matrix/`. The earlier
single-geometry 2026-07-19 arrival times used a full geometric carrier and a
restricted launch fan; they remain useful historical smoke results but must
not be cited as an independent PE delay validation. The 2026-07-20 matrix is
the authoritative phase/timing result. The appropriate next work is a
validation-only study of source normalization/directivity and transverse
window/sponge placement, not a change to the main PE propagation framework.

## 2026-07-21 Bellhop Flat-Surface Visualization Layer

`scripts/reporting/generate_bellhop_flat_surface_visuals_vertical.m` adds an
independent Bellhop display and consistency-check entrypoint. It reads the
saved 3/6/9 m PE/Bellhop matrix, writes and runs standard Bellhop `R`, `C`, and
`I` environments, and parses ASCII `.ray` plus binary two-dimensional
rectilinear `.shd` files without adding the Acoustic Toolbox MATLAB plotting
directories to the path. No propagation-core, public-channel, CIR-builder, or
communication-chain file is modified.

The representative environment remains the 100 m isovelocity (1500 m/s)
water column, 80 m transmitter depth, 10 m receiver depth, and 4 kHz
frequency. The ray fan uses 51 beams over -89.5 to -60 degrees. The coherent
and incoherent TL fields use 5001 and 10001 beams on a 199-depth by 241-range
grid spanning 0.5--99.5 m and 0.25--12 m. The formal figures use 10001 beams,
a shared 20--80 dB color scale, and a 1 m source mask. `ZBOX=99.9 m`, together
with a 0.01 m ray-display step, prevents displayed seabed reflections while
retaining direct and pressure-release surface-reflected rays.

All visualization-layer checks passed. The 5001-to-10001-beam coherent and
incoherent TL RMS differences were 0.19041 and 0.10193 dB; their 95th
percentile absolute differences were 1.0183e-5 and 9.0626e-7 dB. At the
3/6/9 m receiver points, the maximum coherent `.shd` versus saved Bellhop
arrival-synthesis TL difference was 0.04491 dB. The saved PE invariant
`H_f=H_direct_f+H_reflect_f` remained within 3.5762e-18. These successful
display checks do not override the strict cross-geometry matrix result, which
remains `passed=false` because of its documented PE/Bellhop amplitude and C2
window-phase criteria.

Outputs are under `results/visualization/bellhop_flat_surface/`: three PNG
figures, all `.env/.ray/.shd/.prt` files, convergence and receiver CSV tables,
a reusable MAT file, and `bellhop_flat_surface_visual_report.md`. The PE
overlay is restricted to the three saved, globally scaled receiver responses;
the visualization does not invent a PE two-dimensional field.

## 2026-07-22 PE Carrier-Phase Reference Unification

> **Superseded validation status:** the implementation description below
> remains valid, but its reduced-grid-only validation status and stated need
> to rerun F=64 were superseded later on 2026-07-22 by the release-candidate
> closure entry and `reports/pe_phase_reference_release_candidate_report.md`.

### Objective

The direct and surface-reflected PE paths previously removed their own
nominal longitudinal carriers and then added the two reduced envelopes at the
receiver. For the standard geometry the nominal spans are 97 and 103 m, so
the missing relative reference delay is 4 ms. This entry supersedes all prior
statements that public `H_direct_f` and `H_reflect_f` remain unaligned reduced
envelopes or that carrier restoration occurs only inside validators.

### Changed implementation

- Added `apply_pe_channel_phase_reference_vertical.m` as the single
  post-propagation conversion layer.
- Added `paramsV.channel_phase_reference`, default `direct_dsp`, with
  `legacy_reduced` retained for regression.
- Public output now includes selected `H_direct_f/H_reflect_f/H_f`, explicit
  `H_*_reduced_f`, explicit `H_*_physical_f`, reference-frequency scalars,
  and `phase_reference_meta`.
- Cached forward and adjoint runners return reduced, direct-DSP, and physical
  receiver components without changing their PE operators or `q/a` weights.
- Analytic receiver statistics retain reduced `C/P` and expose direct-DSP
  defaults using `C=D*C_red*D^H` and `P=D*P_red*D^T`.
- Conditional models/libraries use schema `2.0.0` / `2.0.0-discrete`.
  `upgrade_conditional_channel_phase_vertical.m` rotates legacy project
  means, `C/P`, complex eigenvectors, and augmented-real structures.
- Added `build_channel_cir_vertical.m`; the old
  `build_physical_cir_vertical.m` remains a common-shift compatibility
  wrapper. Unknown external `H(f)` without project metadata is assumed
  DSP-ready rather than silently rotated.

The joint-kstat definition remains
`deltaG=R0*exp(1i*alpha*eta)-R_coh`; no additional `R0` or deterministic
path carrier is inserted into `deltaG`, the PM lag covariance, `q`, or `a`.

### Validation

`scripts/validation/validate_pe_channel_phase_reference_vertical.m` covers
the F=65 analytic delay audit, all four public surface branches, scalar
direct-only, legacy reproduction, cached/adjoint projection, PM 16^2 dense
versus FFT statistics, schema-1 migration, and unknown-external-H policy.
The generated summary reports:

- nominal/measured relative delay: `4/4 ms`;
- cached/adjoint direct-DSP error: `3.75814e-16`;
- dense/FFT `C/P`: `3.64406e-16/4.16399e-16`;
- migration `C/P`: `3.61301e-16/4.26365e-16`;
- all gates passed.

The existing adjoint smoke also passed its mathematical gates: exact
adjoint `7.77e-15`, projection `5.10e-15`, dense/FFT approximately `1e-15`,
and F=9 split-floor statistics accepted. The reporting helper path was added
to that validator so the smoke run now completes. The F=32 public/cached
double regression passed with `6.9292e-16` total error, and the minimal
wideband communication validation completed. The two-node communication
smoke was then regenerated with explicit direct-DSP metadata for U=5/U=8
cached PE, joint PE, full-rank, and 99.9% conditional sources; total-only
legacy ensemble caches are rejected for reuse because their direct and
reflected components cannot be migrated separately.

### Outputs, documentation, and remaining work

Generated validation data is under
`results/validation/pe_channel_phase_reference/`. The implementation and
compatibility decision are summarized in
`reports/pe_channel_phase_reference_migration_report.md`. `README.md` is now
the concise repository entrypoint; `vertical_comm_guide.md` is the current
technical reference; this file remains the full append-only model log.

The expensive F=64/512 ensemble was not regenerated during this migration.
Its second-order conclusions are invariant under the unitary deterministic
phase rotation, but a future release candidate should rerun the full suite
before removing the `legacy_reduced` compatibility option. Bellhop amplitude
normalization and transverse-window sensitivity remain separate open issues.

## 2026-07-22 Li et al. (2009) Explicit Rough-Surface Validation

An independent pure-acoustic validation path was added without changing the
public PE, surface-boundary defaults, statistical channel generators, or
communication entrypoint. `li2009_u10_to_u195.m` solves the paper's Charnock
Eqs. (12)--(13). `li2009_explicit_surface_validation.m` generates one raw-PM
surface per seed, shares that exact array across every frequency of a 12 kHz,
6 ms CW pulse, applies only `G=-exp(i*2*k*eta)`, and evaluates a monostatic
bottom--surface--bottom uniform split-step PE through an exact discrete
receiver projection. The driver is
`scripts/validation/validate_li2009_explicit_surface_vertical.m`.

The 128-realization reduced baseline uses H=256 m, PE 100 m/128^2, PM
200 m/256^2, dz=4 m, and 65 frequencies from 11.5--12.5 kHz. U10=5/10 m/s
maps to U19.5=5.2684/10.6203 m/s. PM discrete/infinite variance coverage is
0.9953/0.9911, and sample-mean Hs is 0.5908/2.3966 m. Flat pressure-release,
local-amplitude, and `2*k*eta` phase checks are 0, 6.94e-18, and 2.78e-16 rad.

The fixed 20% first-threshold travel-time width increases from 0.3208 to
0.7238 ms, and its shifted-Rayleigh b increases from 0.5531 to 1.6110 ms.
Peak-time width changes only from 2.0713 to 2.0801 ms and fails the 5 percent
broadening gate. Sample convergence passes at 128 for both first-threshold
conditions and the U10=5 peak, but not for the U10=10 peak. Step size,
frequency count, sponge, and PE window comparisons are stable; transverse
grid refinement from dx=0.78125 to 0.52083/0.390625 m remains nonconverged.
The 1024 m/200 m/256^2 coarse-grid check gives the correct flat time but is
not statistically accepted.

Public direct-only/direct-plus-reflect regression passed with component error
3.47e-18 and no direct-path drift. The minimal peak-synced QPSK regression
also completed with zero noiseless and 20 dB BER/SER for the accepted cases.
The full audit, figures, convergence table, interpretation boundaries, and
remaining work are in
`reports/li2009_explicit_surface_validation_report.md` and
`results/validation/li2009_explicit_surface/`.

## 2026-07-22 PE Carrier-Phase Release Candidate Closure

### Goal and scope

The carrier-reference fix was taken through a formal release-candidate
closure without changing the PE marching operator, public surface default,
or unsupported physics. The active run is `phase_rc_20260722_174945`, with
source fingerprint
`14ccc63217ae08a11948e9fb489f2abe8148331f079f6ec7a136a2e607712df1`.
The final decision is **PASS: recommend main-branch integration**.

### Interfaces and artifact policy

- `pe_phase_release_run_meta_vertical.m` creates a reproducible run identity
  containing revision/dirty state, MATLAB version, geometry, frequency axis,
  seed definition, and SHA-256 source fingerprint. A partial run may resume
  only when this fingerprint is unchanged.
- `audit_phase_reference_artifacts_vertical.m` classifies registered MAT
  files without modifying them. `prepare_pe_phase_reference_release_candidate_vertical.m`
  archives phase-sensitive outputs by moving them under
  `results/archive/pe_phase_reference_pre_rc/<timestamp>/`, preserving names
  and relative paths, then creates migration-test copies where valid.
- Required current MAT outputs carry `schema_version`,
  `phase_reference_meta`, and `validation_run_meta`. The atlas refuses inputs
  from a different run ID.
- `finalize_pe_phase_reference_release_candidate_vertical.m` emits only
  `PASS`, `FAIL`, or `INCOMPLETE` and writes the authoritative report at
  `reports/pe_phase_reference_release_candidate_report.md`.

Two archive batches were retained. The first prepared run was invalidated
when the existing U=8 aperture contract was found to be PM `384^2`, not the
adjoint suite's PM `256^2`; its partial results were archived before creating
the final run. No old result was deleted or overwritten. Schema-1 migration
copies were used only for closure tests, not as replacements for the newly
trained schema-2 models.

### Formula and convention

The direct-DSP receiver representation remains

\[
H_{\rm dir}^{\rm dsp}=H_{\rm dir}^{\rm red},\qquad
H_{\rm ref}^{\rm dsp}=e^{-i2\pi f\Delta\tau_0}H_{\rm ref}^{\rm red}.
\]

Receiver statistics are rotated after reduced-domain spatial contraction:

\[
\mu_{\rm dsp}=D\mu_{\rm red},\qquad
C_{\rm dsp}=DC_{\rm red}D^H,\qquad
P_{\rm dsp}=DP_{\rm red}D^T.
\]

No carrier factor enters the marching operator, sea-surface model,
`deltaG`, adjoint kernel `q`, receiver weight `a`, or PM spatial covariance.
The 4 ms standard-geometry delay is a deterministic longitudinal reference,
not an alignment of realization-dependent peaks.

### Formal configurations and results

- F=65 phase audit: 4--8 kHz at 62.5 Hz spacing. Nominal/measured relative
  delay closed within one sample; phase representation round trips and
  migration closure passed.
- Exact adjoint full: F=9 used 4096 realizations in batches of 16; F=64 used
  512 realizations in batches of 8 on PE `128^2` / PM `256^2`. Exact-adjoint,
  projection, dense/FFT `C/P`, and same-realization cached/projection checks
  all passed. F=9 and F=64 analytic/sample covariance ratios to split floors
  were `0.593078` and `0.850502`, both below `1.25`.
- Conditional nodes: U=5 used 128/64 train/test; U=8 used 128/128 and its
  audited PE `128^2` / PM `384^2` aperture. Full, 99.9%, and 99% variants and
  10,000-sample generator checks passed; models are schema 2.x direct-DSP
  with zero common CIR shift.
- Two-node communication: all four channel sources completed at U=5 and U=8
  with 32 channels per source and 8000 QPSK symbols. Full/low-rank paired
  comparisons used the same latent variables and all inputs were finite,
  direct-DSP tagged, and IFFT-consistent.
- Public scalar/direct-only/direct-plus-reflect and F=64 cached/public
  regressions passed. The F=64 total consistency error was
  `9.88957e-16`; the public default remained `kirchhoff_spatial`.
- The rebuilt same-run atlas passed its artifact checks and contains 13 PNG,
  one MP4, MAT, manifest, and summary. Spatial panels remain reduced complex
  envelopes; carrier reference is applied only to receiver combinations.

### Architecture decision and remaining limitations

The approved roles are: public PE for general propagation; cached forward as
the fixed-path regression oracle; adjoint projection for exact fixed-environment
receiver realizations; analytic FFT `C/P` for receiver second-order
statistics; and the conditional generator for large communication Monte
Carlo after physical/statistical validation. This does not extend the exact
adjoint v1 beyond uniform CPU-double, fixed grids/frequencies, one
nearest-grid receiver, or into bubbles/Doppler/layered/GPU modes.

Bellhop source normalization and doubled-window sensitivity remain separate
open studies. The atlas produced one harmless title-format warning involving
`\Delta`; it did not affect saved data, validation checks, or deliverables.

## 2026-07-23 Current-Version PE/Bellhop Flat-Surface Revalidation

### Task and implementation

This task replaced the old mixed-version Bellhop interpretation with an
independent, run-ID-scoped validation batch. New interfaces are
`pe_bellhop_current_run_meta_vertical.m`,
`run_bellhop_flat_surface_arrivals_vertical.m`,
`scripts/validation/validate_pe_bellhop_flat_surface_current_vertical.m`, and
`scripts/reporting/generate_pe_bellhop_comparison_atlas_vertical.m`. The
existing phase validator now accepts an explicit output directory and shared
run metadata; the old Bellhop smoke validator explicitly requests
`direct_dsp` and documents the public physical-field sign. No public PE
marching, default surface model, communication entry, or carrier formula was
changed.

Formal mode requires a stable absolute non-Temp `BELLHOP_EXE`. The accepted
external executable was
the historical AcousticsToolbox 2017 Bellhop executable (external dependency), size 683637 bytes,
SHA-256
`e6f9c1bcfd2b0945bfb59607909af589fa7b3f823df8bd5121425736eaebd796`.
The wrapper writes standard uniform/flat pressure-release `.env` cases,
parses `.arr`, and classifies direct and one-surface paths by top/bottom bounce
count and 0.1 ms delay clusters. It rejects case roots with filename
extensions to avoid Bellhop root-name truncation.

### Configuration, formulas, and semantics

The final run is `bellhop_current_20260723_rc5`: Tx `[0,0,80] m`, Rx depth
10 m, 3/6/9 m baselines, 4 kHz, PE 3--5 kHz/33 frequencies, `128^2` over
32 m, `stepz_lamb=0.5`, uniform 1500 m/s, flat pressure-release surface,
CPU double. The independent phase audit used F=65. Small-offset cases were
0/0.25/0.5/1/3 m; Bellhop was not called at strict zero range. Bellhop fields
used R/C/I and 5001/10001 beams. The PE window audit separated sampling,
no-sponge 32/64 m aperture, sponge ratio 0.10/0.12/0.15, and sponge strength
0.05/0.15/0.30. Window cases used `save_mode='slice'` so edge amplitude was
measured from an actual receiver plane.

The validator never restores a carrier manually. Reduced fields are used
only for conversion closure, public `H_*_physical_f` for physical-delay
comparisons, and direct-DSP fields for public/IFFT regression. Bellhop's
4 kHz arrivals are extended over the PE frequency band only for delay/PDP
diagnostics, not presented as frequency-resolved Bellhop amplitudes.

### Results and decision

The final layered decision is **FAIL_CORE**, with
`amplitude_status=OPEN`. Passed evidence includes: independent phase audit;
reduced/physical/direct-DSP conversion; `H_f` component closure; direct and
surface path classification; baseline and small-offset delays; Bellhop beam
convergence; PE sampling/step; sponge location/strength; scalar direct-only,
direct-plus-reflect, wideband closure; and unchanged public defaults.

Maximum PE--analytic and PE--Bellhop delay errors were `0.027281981 ms` and
`0.027280961 ms`; Bellhop--analytic was `3.930416e-06 ms`. Coherent and
incoherent beam-count RMS differences were `0.190410/0.101930 dB`. Phase
conversion error was zero and the maximum public component closure error was
`3.576224e-18`.

The failed core item is specifically aperture validity. With sponge disabled,
the 32/64 m receiver-plane edge levels were `-1.087/-14.532 dB`, not below the
required `-40 dB`; the 32--64 m phase RMS and relative TL difference were
`3.526186 rad` and `-3.896489 dB`. These results cannot serve as an aperture
convergence reference. This failure does not overturn the independently
closed phase, path, delay, or Bellhop beam results and does not motivate a PE
core or carrier-reference change.

Raw amplitude remained open: one global direct scale produced `1.65401 dB`
cross-geometry spread and reflected RMS/max residuals of
`1.62948/2.07510 dB`. The Gaussian angular-spectrum direct-only diagnostic
reduced RMS from `0.815` to `0.502 dB`, but it never replaces raw Bellhop point
source values or changes acceptance.

### Artifacts, documentation, and remaining work

The final MAT/CSV/Bellhop artifacts are under
`results/validation/pe_bellhop_flat_surface_current/bellhop_current_20260723_rc5/`.
The same-run visualization directory contains 10 numbered PNGs, an atlas MAT
with comparison tables/checks, manifest, summary, and raw R/C/I field evidence.
The current report is
`reports/pe_bellhop_flat_surface_current_validation_report.md`; the 2026-07-21
complete report is retained and marked Superseded. README, technical guide,
and script index were synchronized.

Remaining Bellhop work is deliberately separate from this task: enlarge the
no-sponge aperture until the edge prerequisite is met, then repeat the
aperture comparison; establish a common point/finite-width source amplitude
normalization. Rough surfaces, joint-kstat, bottom interaction, stochastic
channels, and communication performance were not part of this validation.

## 2026-08-12 Reflection-Free PE/Bellhop Free-Space Audit

The validation layer now has an additive custom initial-field hook while the
public default remains Gaussian. `source_mode='custom_field_fn'` is used only
by the free-space validators to initialize a virtual point source on the PE
plane. `alpha_max_np_per_m=0` is an explicit, valid sponge-off setting. The
production marching equations, exported channel fields, reflection branches,
and communication consumers are unchanged.

The evidence chain is separated into four levels: Bellhop normalization,
production multi-step PE versus an independently coded one-step exact
discrete angular spectrum, source-referenced PE versus
`exp(i*k*R)/(4*pi*R)`, and source-referenced PE versus converted Bellhop.
The validation stores both `H_plane` and `H_source=exp(i*k*s0)*H_plane` so the
virtual-source carrier is restored only for physical comparisons.

Bellhop uses the matched upper/lower halfspace construction from its official
free-space point-source example. The audit measures `|p|R`, spatial phase
sign, one constant source phase, beam/step convergence, and zero-bounce delay.
The `1/(4*pi)` conversion is fixed once for every distance and frequency; no
pointwise calibration is permitted. Outputs are under
`results/validation/pe_bellhop_freefield/formal/`, with the decision report in
`reports/pe_bellhop_freefield_validation_report.md`.

The formal run used 3/4/5 kHz, `s0=[5,10,20] m`, `L=[20,40,70,100] m`,
offsets `[0,0.5,1,2] m`, and a 256-square PE grid over 64 m. PE--AS passed at
`1.42870e-13` maximum relative L2 error. Bellhop normalization and analytic
comparison passed (`1.99322e-6 dB` maximum TL error), but PE--Bellhop failed
with `18.2012 dB` maximum TL error, `1.38356 rad` phase RMS, and a
`3.76572 ms` representative group-delay residual. The 32-to-64 m field study
improved but did not converge to the infinite-aperture spherical field. The
status is therefore `passed=false`: the production propagation operator is
verified, while validation-local finite-plane point-source reconstruction is
not. No PE marching formula or threshold was changed to force acceptance.

## 2026-08-12 Point-Source Window/Sponge Error Budget

The PE square-root marching operator and public Gaussian source semantics were
explicitly frozen. A validation-only fixed-`dx=0.25 m` series used 24, 32, 48,
64, 80, and 100 m transverse windows with sponge disabled. Production PE and
independent one-step AS continued to agree within `8.65585e-13`, but neither
the spatially truncated spherical plane nor the discrete FFT-Weyl plane
converged monotonically to the analytic point source. At the public 100 m
limit the spatial-plane amplitude error remained `+12.8509 dB`, equivalently
`TL_PE - TL_analytic = -12.8509 dB`; the last two discrete
Weyl windows differed by `3.23331 dB`.

The raw Weyl `1/kz` spectrum has a grazing-circle integrable singularity, so a
separate continuous axisymmetric Weyl reference was implemented with
propagating `kz` and evanescent `q` substitutions. It matched
`exp(i*k*R)/(4*pi*R)` to `2.41716e-13` relative error and is the reliable
free-space baseline. The discrete FFT-Weyl initializer remains an explicitly
failed diagnostic rather than a replacement production source.

The sponge audit now covers the complete 6-window x 3-ratio x 6-strength
matrix (108 rows). It defines `dA=20*log10(|H_sponge|/|H_no_sponge|)` and
`dTL=TL_sponge-TL_no_sponge=-dA`. Across the complete matrix the maximum
on-axis additional loss is `35.6811 dB` (24 m, ratio 0.15, 0.3 Np/m); in the
100 m window alone the maximum is `5.99769 dB` (ratio 0.15, 0.3 Np/m).

In the final 4 kHz, `s0=10 m`, `L=70 m`, offset 0/0.5/1/2 m Bellhop rerun,
finite-window amplitude errors are `+12.851/+5.520/-0.613/+9.612 dB`, hence
the correctly signed TL errors are `-12.851/-5.520/+0.613/-9.612 dB`. For
the selected ratio 0.12/strength 0.15 sponge, the amplitude changes are
`-4.585/-16.223/-6.268/-5.092 dB`, equivalently additional TL values of
`+4.585/+16.223/+6.268/+5.092 dB`. The axis phase change is `-0.980475 rad`.
At the terminal plane, center (`rho<=2 m`), edge (ratio-defined sponge band),
and total energies change by `-2.40739/-5.16144/-3.83990 dB`, respectively.
Final PE--Bellhop TL errors are `-8.266/+10.703/+6.881/-4.520 dB`.
Bellhop--analytic remained below `8e-7 dB`.
The post-budget comparison therefore remains `passed=false`; its disagreement
is quantitatively assigned to finite-plane/periodic-image and sponge effects,
not to the already verified PE marching operator.

## 2026-08-12 Production Gaussian Window/Sponge Audit

This follow-up stops point-source/Bellhop expansion and freezes the square-root
marching operator, FFT convention, Gaussian mathematics, surface model, and
communication modulation. The exact production Gaussian expression was moved
without algebraic change into `gaussian_source_initial_field_vertical.m` so
production and validation call the same implementation; direct comparison to
the old expression is exactly zero. Validation-only opt-ins permit 0.05--0.20
sponge ratios and windows up to 200 m while normal public defaults and limits
remain unchanged.

The homogeneous direct-only audit uses 4 kHz, `sigma=0.3 m`, fixed
`dx=50/256 m`, path lengths 20/40/70/97/100 m, and no bubbles, surface,
Doppler, randomness, or communication processing. Gaussian PE--independent AS
passed with maximum full/center complex L2 `2.56337e-13/1.45326e-13`.
No-sponge 128 m converged to 160 m at 4 kHz, but 128 m had a 3--5 kHz maximum
TL error of `0.100293 dB`, narrowly exceeding the preregistered 0.1 dB target.

At the production 50 m/97 m case the no-sponge field reaches the FFT boundary:
outer 5%/10% energy fractions are `0.1130/0.2176` and boundary amplitude is
only `-1.2879 dB` relative to the plane maximum. The current ratio 0.12,
alpha 0.15 sponge gives `dA=+0.62225 dB`, `dTL=-0.62225 dB`, axis phase
change `-0.148555 rad`, center (`rho<=2 m`) energy change `-0.047554 dB`,
edge change `-2.77763 dB`, and total change `-0.824556 dB`. None of the 126
window/ratio/strength rows met all fixed center targets plus 3 dB edge
suppression.

Against 160 m/no-sponge over 3--5 kHz, production/default has maximum TL error
`0.842225 dB`, phase RMS/max `0.131866/0.285952 rad`, and group-delay RMS/max
`232.649/1326.814 us`. The 128 m/no-sponge compromise gives
`0.100293 dB`, `0.003286/0.011962 rad`, and `8.636/33.558 us`.
The strict recommendation is therefore 160 m/no-sponge; 128 m/no-sponge is a
near-threshold cost option. No production default was changed. Full report and
16 figures are under `reports/pe_gaussian_window_sponge_validation_report.md`
and `results/validation/pe_gaussian_window_sponge/`.

## 2026-08-13 Full Reflected-Chain Window Convergence

The complete production Gaussian PE -> Kirchhoff surface -> PE chain was
validated without changing the PE marching operator or defaults. A single
seed-12345, U=5 m/s PM surface was generated on the maximum domain, normalized
once to Hs=0.5 m, and centrally cropped at fixed `dx=50/256 m`; crops were
never renormalized. The 4 kHz no-sponge sequence tested actual widths
50.0000, 128.1250, 160.1563, and 192.1875 m.

The 160-to-192 m comparison gave reflected TL/phase differences of
`-0.0015395 dB / 0.0006122 rad` and receiver-center complex L2 `4.546e-4`.
Nevertheless, 160 m had outer-5% energy fraction `1.1769e-4`, narrowly above
the fixed `1e-4` gate. The 192.1875 m field passed with boundary amplitude
`-52.29 dB` and outer-5% energy `1.7699e-5`; it is therefore the strict
minimum qualified window. Edge tracking attributes the main contamination to
the 100 m upward march, with only modest additional spreading after the
Kirchhoff screen and 3 m downward march.

Across 3--5 kHz/33 points, 160.1563 m/no-sponge versus 192.1875 m/no-sponge
has maximum reflected TL, phase, and internal group-delay differences of
`0.019205 dB`, `0.0023817 rad`, and `4.0507 us`. The 50 m/default-sponge case
fails at `1.3203 dB`, `0.17782 rad`, and `70.943 us`; its total-channel TL and
phase differences reach `8.4134 dB` and `0.64120 rad`. All complex component
closures are below `3.9e-18`. The strict recommendation is therefore
192.1875 m/no-sponge; 160 m is a response-accurate cost option but does not
pass the preregistered edge gate. Production defaults were not changed.

Implementation adds `surface_elevation_override_xy` and low-memory
`surface_wavefield_meta` energy/receiver diagnostics as disabled-by-default
interfaces. Formal evidence is in
`reports/pe_reflected_chain_window_validation_report.md` and
`results/validation/pe_reflected_chain_window/`.

## 2026-08-13 Random-Surface Window Robustness

The fixed-realization conclusion was extended at 4 kHz to three target sea
surface strengths (`Hs=0.05/0.5/1.0 m`) and five seeds per strength. Each
realization was generated and normalized once on a 255.8594 m master grid;
160.1563 and 192.1875 m cases are exact central crops at fixed
`dx=50/256 m` and are never renormalized. All paths use sponge off and the
unchanged production Gaussian PE -> Kirchhoff spatial -> PE chain.

Relative to the matching 256 m reference, 192 m passed all response,
center-field, boundary, and outer-5% energy gates in `15/15` realizations.
Its worst 4 kHz errors were `0.0185269 dB` TL, `0.000698941 rad` phase, and
`8.02062e-4` receiver-center complex L2; worst outer-5% energy was
`5.94925e-5`. The 160 m cases passed the response and center-field gates in
`15/15`, but passed the full edge gate in `0/15`; worst outer-5% energy was
`1.89970e-4`. Stronger Hs increased lateral spreading, but did not require a
window larger than 192 m in this sample.

The selected random-realization 3--5 kHz qualification was stopped before a
complete set of candidate--256 m pairs was obtained, so no random-ensemble
worst group-delay claim is made. The production recommendation remains
192.1875 m with sponge off; 160 m is only a monitored cost option. Production
defaults were not changed. See
`reports/pe_random_surface_window_robustness_report.md` and
`results/validation/pe_random_surface_window_robustness/`.

## 2026-08-20 Repository Layout Migration and Documentation Consolidation

### Goal and status

Reorganized the repository with compatibility-first semantics after multiple
validation windows had added files at the root. Status: **IMPLEMENTED** for
layout and documentation; lightweight post-migration regression results are
recorded in `reports/current_project_status_20260820.md`. No expensive
F=64/512 ensemble, complete Bellhop matrix, or 192/256 m formal aperture run
was repeated for this organizational change.

### Changed paths and interfaces

- Root public names remain `vertical_channel_model`, `main_vertical`,
  `explain_main_vertical`, and `comm_main_vertical_psk`.
- `vertical_channel_model.m` is now a short compatibility wrapper around
  `src/channel/vertical_channel_model_impl.m`.
- The three demonstration entrypoints run the corresponding scripts in
  `examples/` and preserve script-workspace behavior.
- Reusable code moved into `src/channel`, `src/propagation`, `src/surface`,
  `src/receiver`, `src/statistics`, `src/communication`, and `src/bubble`.
- Bellhop/Weyl/Li2009/properness/release-metadata helpers moved into
  `scripts/validation/support/`.
- `setup_vertical_project.m` is the idempotent initializer and does not
  change the current directory. `scripts/bootstrap_project.m` delegates path
  setup to it before applying the existing result-directory behavior.
- The root 2k scripts moved to `scripts/validation/`; their PNG/CSV outputs
  moved to `results/validation/ssa_2k_phase/`, and their reports moved to
  `reports/`.
- `SSA.md`, `BUBBLE_EXTENSION_SPEC.md`, and
  `README_CODE_STRUCTURE_AND_OUTPUTS.md` moved verbatim to `docs/history/`.
  They are explicitly non-authoritative. Active SSA metadata now refers to
  `vertical_comm_guide.md` and the cited formal literature.

The migration mapping and pre-migration SHA-256 evidence are recorded in
`reports/repository_layout_migration_20260820.csv` and its Markdown summary.
Dirty files were moved using their current contents; no reset, checkout,
deletion of user work, or overwrite of prior validation results was used.

### Physics and compatibility

This change deliberately does not modify the PE split-step order, Gaussian
source, sponge formula, surface reflection formula, joint-kstat `deltaG`,
adjoint `q/a`, PM covariance contraction, carrier-phase conversion, public
surface default, window default, or bubble default. Existing output field
names and `H_f=H_direct_f+H_reflect_f` semantics remain required.

### Documentation

- `README.md` is now the concise Chinese repository entry and navigation.
- `vertical_comm_guide.md` is the current Chinese technical authority for
  coordinates, PE envelope, surfaces, adjoint/statistics, CIR/LFM,
  communication outputs, validation boundaries, and recommendations.
- `PROJECT_CONTEXT.md` remains append-only long-term model memory: history is
  preserved and this top index points to current paths and conclusions.
- `scripts/README.md`, `src/README.md`, `docs/history/README.md`, and the
  dated status report define non-overlapping operational roles.

### Post-migration validation

- Root/external-directory setup and single-definition `which -all` checks
  passed. Code Analyzer scanned 155 MATLAB files with zero parse errors.
- Reduced scalar direct-only/direct-plus-reflect closure passed; wrapper vs
  implementation difference was zero, total component closure was
  `3.47e-18`, and seeded/override/diagnostic transparency passed.
- PE 32^2 / PM 64^2, 4/6/8 kHz cached/adjoint smoke gave inner-product error
  `3.35e-15` and reduced/direct-DSP projection errors
  `1.04e-15/1.05e-15`.
- The F=9, PE 64^2 QPSK/AWGN smoke had IFFT closure `2.96e-16`, channel
  closure `5.55e-17`, and no NaN/Inf.
- Both relocated 2k scripts completed; all five key CSV SHA-256 values were
  identical to the pre-migration files.
- Several MATLAB batches (including cold setup and FFT-heavy checks) printed
  passing assertions and then hit a local R2025b DDUX/threadpool
  `std::terminate` during process shutdown. Other batches exited normally;
  the shutdown issue is recorded rather than reported as a clean process
  exit.

### Remaining issues

- Public 50 m/default-sponge settings remain unchanged although the strict
  complete-reflection recommendation is 192.1875 m/no-sponge.
- Random-surface paired wideband convergence is incomplete.
- Bellhop absolute amplitude/source normalization and infinite-aperture point
  source equivalence remain open.
- Li2009/SSA transverse convergence and higher-order validation remain open.
- The exact-adjoint/statistics v1 scope is still uniform CPU double, fixed
  grid/frequency, one nearest-grid receiver, and no bubble/Doppler.

## 2026-08-20 Unfolded-coordinate Gaussian PE--Bellhop validator

Added the validation-only entry
`scripts/validation/validate_pe_bellhop_unfolded_flat_gaussian_vertical.m`.
It maps the vertical 100 m-to-3 m geometry to a Bellhop horizontal free-space
case: the direct path is 97 m, while the pressure-release image path is 103 m
and is multiplied by `-1` at the comparison stage. Bellhop uses matched upper
and lower halfspaces, so the artificial depth boundaries do not add rays.

The entry calls the existing production Gaussian source through
`vertical_channel_model` (sigma `0.3 m`), keeps the converged PE aperture
(`192.1875 m`, `N=984`) and sponge off by default, writes a frequency-dependent
`.sbp` directionality from the Gaussian angular spectrum, and compares PE,
independent one-step continuous Weyl/Hankel axis references, and Bellhop using
relative quantities `G=H/H(0)`, `Q=H_reflect/H_direct`, and `1+Q`. It includes
5001/10001 beam convergence checks, 4--8 kHz relative PDP/group-delay outputs,
spatial profiles, the separate PE unfolded identity, and the existing small
offset native-coordinate arrival regression. Bellhop SHD files with source
beam patterns are parsed by the dedicated
`scripts/validation/support/read_bellhop_shd_unfolded_vertical.m`; no external
Acoustic Toolbox plotting path is retained.

The 128x128, 32 m smoke run completed and correctly passed beam convergence,
the PE unfolded identity (`2.51e-13`) and component closure
(`3.47e-18`). Its spatial/wideband checks fail as expected because that small
aperture is intentionally contaminated by periodic-window effects. The full
984x984, 65-frequency run subsequently completed and saved 65 formal Bellhop
shade files, 3 spatial-profile files, 6 beam-audit files, MAT/CSV/PNG output,
and the generated report. Beam convergence, wideband TL/phase, group delay,
PDP, the independent PE--AS identity (`4.74e-13`), native auxiliary arrivals,
and `H_f=H_direct_f+H_reflect_f` closure (`3.88e-18`) passed. The only failed
hard checks were the transverse profile phase (`0.0612 rad` versus `0.05`) and
transverse profile complex L2 (`0.0613` versus `0.02`); axial/wideband relative
channel agreement passed, so the validator remains `passed=false` and no
threshold was relaxed. The point-source normalization layer used the prior
validated saved audit because the installed legacy Bellhop executable crashed
on a fresh ±180° point-source audit; this external limitation is recorded in
the report. A 256x256, 8/16-frequency direct-plus-reflect communication smoke
and the reduced reflected-chain diagnostic regression also completed
successfully. No PE marching, Gaussian source, public window, or sponge
defaults were changed.

The detailed report now recommends a low-cost transverse three-way audit as
the next action: compare the saved 4/6/8 kHz PE slices against an independent
one-step exact angular-spectrum slice and then against Bellhop. This should
separate PE extraction errors from Bellhop source-mapping/2D-ray-beam phase
effects without rerunning the 65-frequency production PE batch. Until that
audit explains the `0.0612 rad` outer-profile discrepancy, the formal status
remains failed even though the axial 4--8 kHz relative channel passed.

## 2026-08-25 Independent PE--AS transverse audit and Bellhop scan

The follow-up entry `scripts/validation/validate_pe_as_bellhop_transverse_vertical.m`
now reads the saved formal result and independently computes one-step discrete
angular-spectrum propagation at 4/6/8 kHz for the 97 m direct and 103 m image
paths. PE--AS agrees to numerical precision: the maximum complex errors are
`5.10e-13` (direct) and `2.43e-13` (reflected). Bellhop--AS retains the formal
outer-profile discrepancy: `0.01582 dB / 0.06120 rad / 0.06117` direct and
`0.01274 dB / 0.05167 rad / 0.05165` reflected (TL/phase/complex). This
independently rules out the PE marching operator and field extraction as the
source of the failed transverse hard check.

The validation-only entry
`scripts/validation/validate_bellhop_transverse_parameter_scan_vertical.m`
then ran 27 Bellhop-only cases at 8 kHz, varying step (`0.1/0.05/0.025 m`),
angle half-width (`20/30/45 deg`), and SBP sampling (`1201/2401/4801`) with
10001 beams. All 27 `.env/.sbp/.shd/.prt` cases parsed successfully, but no
case passed the transverse limits. The best maximum complex error was
`0.061166`, with phase still `0.061204 rad`; numerical step, angle, and SBP
sampling therefore do not explain the discrepancy. Bellhop outputs were
conjugated once to match the formally selected spatial phase convention; no
frequency- or distance-dependent fit was used.

Artifacts are under `results/validation/pe_as_bellhop_transverse/` and the
reports are `reports/pe_as_bellhop_transverse_audit_report.md` and
`reports/bellhop_transverse_parameter_scan_report.md`. The formal unfolded
validator remains `passed=false` solely because its transverse phase/complex
hard checks fail. The next investigation is the 2-D/3-D Gaussian source
mapping (`.sbp` amplitude/Jacobian/phase interpretation), not PE core changes,
window changes, sponge changes, or threshold relaxation.

## 2026-08-25 Bellhop unfolded-coordinate rotation-limit audit

Added the Bellhop-only entry
`scripts/validation/validate_bellhop_rotation_limit_vertical.m`. It compares
the current unfolded geometry (97/103 m direct/image ranges) with native flat
pressure-release-surface cases whose direct rays are `-89` and `-89.5` degrees.
All three cases use the same 4 kHz, c=1500 m/s, Tx=100 m, Rx=3 m Gaussian
`.sbp` pattern; native horizontal offsets are 1.69314 m and 0.846506 m.
The native bottom is placed far away and no bottom-bounce arrival is accepted.

The run found two required paths in every case. Maximum analytic delay error was
`3.13 ns`; after removing the known geometric path length and rotated Gaussian
directivity, native-vs-unfolded amplitude differences were below `1.23e-6 dB`
and phase differences below `9.30e-5 rad`. The corrected reflection coefficient
was within `1.23e-4` of `-1`. All five structural/consistency checks passed.
This supports the unfolded construction as a near-vertical limit for the
uniform flat-surface Bellhop case only; it does not extend to rough surfaces,
depth-dependent SSPs, or actual-transducer source patterns.

## 2026-08-26 Current PE--Bellhop validation overview

`reports/pe_bellhop_complete_validation_overview.md` is now the current
cross-stage summary for the uniform, flat pressure-release-surface comparison.
It links the historical native-coordinate path matrix, Bellhop free-space
normalization, independent PE--AS checks, production-Gaussian aperture work,
the 4--8 kHz unfolded comparison, the transverse three-way audit, the 27-case
Bellhop numerical scan, and the 89/89.5-degree coordinate-limit result.

The summary deliberately retains both parts of the current decision: axial
relative `Q(f)`, `1+Q(f)`, delay, group delay, and PDP checks pass, and PE--AS
agrees at about `1e-13`; however, the formal unfolded validator remains
`passed=false` because the 8 kHz outer transverse phase (`0.061204 rad`) and
complex error (`0.061277`) exceed their fixed hard limits. The three-way audit
assigns that residual to the Bellhop/source-mapping or 2-D/3-D representation
side rather than the PE marching operator. The overview supersedes the old
2026-07 “complete” report only as a status summary; historical files and raw
results remain unchanged.

## 2026-09-03 PE--Bellhop PM model-discrepancy statistical study (M=50)

The fixed-4 kHz Stage-4 PM ensemble used the frozen Stage-3B independent
Gaussian coefficient-amplitude chain.  The user stopped the conditional
extension after 50 complete seeds (`260001:260050`); no PE or Bellhop core
physics, source normalization, receiver selector, or influence formula was
changed.  The active report is
`reports/pe_bellhop_pm_model_discrepancy_statistics_report.md`; seed-level
tables and figures are under
`results/validation/pe_bellhop_pm_model_discrepancy_statistics/`.

At M=50 the mean delta TL is `-0.00748 dB` (sample standard deviation
`0.96223 dB`) and the circular mean phase difference is `-0.37534 rad`.
The 24-to-32 engineering gates pass for mean delta-TL, PE/Bellhop power
change, and circular phase change; the bootstrap mean-delta-TL half-width is
`0.25835 dB`, just above the prescribed `0.25 dB` gate.  Classification is
therefore `PRELIMINARY_MODEL_DISCREPANCY`, not statistical establishment.
All 50 rows pass finite-value, shared-profile, beam-hit, non-grazing,
wall-residual, pressure-release phase, beam-state, positive-range, and PE
edge/seam guards.  Native backward-range amplitude remains diagnostic only,
and the dimensionality sensitivity (`-0.15886 dB`, `+0.01251 rad`) is
reported separately rather than subtracted.

## 2026-09-10 Controlled PE--Bellhop Goal: Stage 0/1 status

The authoritative sequential driver is
`scripts/validation/validate_pe_bellhop_controlled_comparison.m`, governed by
`reports/pe_bellhop_controlled_comparison_GOAL_revised_stage1x.md`.  P0 and Stage 0 passed,
including PE--AS closure, flat parametric internal-wall closure, receiver
mapping, phase/tau/positive-range guards, and 5001/10001-beam checks.

Stage 1 weak-sinusoid runs at `A=0.01`, `0.005`, and `0.0025 m` all passed
profile/geometry/finite guards but did not return to the Stage-0 numerical
floor: model M99 L2 was `0.47247`, `0.23875`, and `0.11969`, respectively;
phase RMS was `0.47920`, `0.23960`, and `0.11980 rad`.  The controlled
comparison is therefore blocked by model comparability at the current weak
limit, not by PE marching, Bellhop geometry, or solver convergence.  The
read-only Stage 1X audit
`reports/pe_bellhop_controlled_comparison_stage1x_phase_convention_audit.md`
identifies a fixed PE conjugation convention that closes all three amplitudes
(the mirror-conjugate is degenerate for the centered/even profile).  This does
not rewrite the original Stage-1 PASS/FAIL: Stage 2--7 remain locked pending a
narrow convention-only rerun.  Per the latest user execution override, future
Bellhop runs use 10,001 beams only; no 20,001-beam endpoint is to be started.

The convention-only regression is recorded under
`results/validation/pe_bellhop_controlled_comparison/stage1_convention_fixed/`.
It uses `G_BH_comparison=conj(G_BH)` only in the rough/flat comparison layer:
Stage 0 was rerun and passed, A=0.01 was rerun through the main driver, and the
A=0.005/A=0.0025 outputs are comparison-only reconstructions from preserved raw
solver MAT files.  The original Stage-1 FAIL data remain under `stage1/`;
Stage 2--7 are not automatically unlocked.

## 2026-09-10 Controlled PE--Bellhop Goal: Stage 1Y convention closure

Stage 1Y is complete and is documented in
`reports/pe_bellhop_controlled_comparison_stage1y_theoretical_convention_closure.md`.
Without new PE/Bellhop runs, the source audit independently closes the PE
`exp(-i*omega*t)`/`exp(+i*(kz-k0)L)` convention, Bellhop 2020 coherent
`exp(-i*(omega*tau-phase))` influence convention, the direct real/imag SHD
reader, the real X-source normalization, and the single vacuum `+pi` phase in
`Reflect2D`. Existing free-field and non-sinusoidal constant-height audits
support the analytic prediction; no per-case conjugation, fitting, or core
physics change is used.

The permanently frozen comparison fields are
`B_abs=conj(B_raw)` and `G_BH_comparison=conj(G_BH_abs)=G_BH_raw`. The original
raw Stage-1 FAIL outputs remain unchanged, while the convention-fixed
regression remains PASS. Stage 1Y is
`CONVENTION_THEORETICALLY_CLOSED`, so Stage 2 is unlocked by the Goal gate;
Stage 2 was subsequently completed and is summarized in the Stage-2 status
entry below. The remaining Bellhop--AS/flat
residual (~`0.006` normalized complex L2) is retained as a finite
source/beam mapping floor, not treated as a convention blocker.

## 2026-09-11 Controlled PE--Bellhop Goal: Stage 2 execution status

Stage 2 height-sweep orchestration is present in
`scripts/validation/validate_pe_bellhop_controlled_comparison.m`. It fixes
`K=0.10 rad/m`, X/C source settings, `N=4097`, and the user-mandated 10,001
beams, and supports case-granular resume without rerunning completed cases.
The prescribed A=`0.01/0.02/0.05/0.10/0.20 m` cases are complete; A=0.01
reused the convention-fixed result and the other four points ran one Bellhop
case each. All solver/mapping/finite guards pass. Metrics rise with height:
`E_G`=`0.004093/0.008234/0.021036/0.044179/0.100485` and phase RMS
`0.004059/0.008166/0.020869/0.043877/0.100142 rad` in that A order.
The authoritative report is
`reports/pe_bellhop_controlled_comparison_stage2_height_sweep_report.md`:
Stage 2 is `PASS`. Its summary was later augmented with the Goal-required
`rho_raw`, maximum slope, and maximum/RMS curvature columns without rerunning
the solvers.

## 2026-09-12 Controlled PE--Bellhop Goal: Stage 3 execution status

Stage 3 is implemented as the `stage3` branch of
`scripts/validation/validate_pe_bellhop_controlled_comparison.m`. The frozen
manifest selects `A=0.02 m`, the largest Stage-2 point satisfying every Region-I
gate, and scans `K=0.10/0.20/0.35/0.47 rad/m` at X/C, `N=4097`, and 10,001
beams. All four solver/mapping/geometry/finite guards pass. The classifications
are `I/II/II/II`; corresponding `E_G` values are
`0.008234/0.008369/0.014140/0.022816`, and phase RMS values are
`0.008166/0.007579/0.009428/0.012906 rad`. The largest sampled curvature is
`0.004418 1/m`; Bellhop hit-curvature, native `Reflect2D` state, wall residual,
and positive transformed range remain finite and consistent. The authoritative
report is
`reports/pe_bellhop_controlled_comparison_stage3_slope_curvature_report.md`.
Stage 4 subsequently completed as a read-only validity-map aggregation with
zero solver calls. At `K=0.10 rad/m`, the sampled Region-I to non-I transition
is bracketed by `A=(0.02,0.05] m`; at `A=0.02 m`, it is bracketed by
`K=(0.10,0.20] rad/m`. These are not interpolated or extrapolated boundaries.
See `reports/pe_bellhop_controlled_comparison_stage4_validity_map_report.md`.

## 2026-09-12 Controlled PE--Bellhop Goal: Stage 5 attribution

Stage 5 completed as a zero-solver-call diagnostic on three Region-II
representatives: `(A,K)=(0.05,0.10)`, `(0.20,0.10)`, and `(0.02,0.47)` in
meter/rad-per-meter units. Global phase alignment reduces `E_G` by only
`7.85%/20.24%/7.94%`, so every case is classified as `spatial distortion`
under the frozen Goal rule rather than global/coherent-phase dominated.
Bellhop differs from the analytic stationary-path phase by only
`1.7e-6/6.2e-6/1.3e-5 rad` RMS, whereas the corresponding PE residual is
`0.02087/0.10014/0.01290 rad`. These quantities are validation diagnostics;
the production PE phase screen and Bellhop physics remain unchanged. See
`reports/pe_bellhop_controlled_comparison_stage5_phase_attribution_report.md`.

## 2026-09-13 Controlled PE--Bellhop Goal: Stages 6--7

Stage 6 copied and verified the canonical seed-260001 coefficients with
SHA-256 `1f7eda465e4ae85b8ac038310edf053f2d062563a3b943bfd83bd59d07015f67`,
then ran the sole allowed 4 kHz X/C, 10,001-beam fixed-PM case. All Bellhop
geometry, incidence, reflection-phase, beam-state, positive-range and finite
guards pass. The axis sanity result is `0.310370739 dB / -2.24820073 rad`,
matching the historical X audit, while the complete M99 line has
`E_G=0.958269`, phase RMS `1.107714 rad`, TL RMS `1.273990 dB`,
`rho_shape=0.771213`, and is Region III. See
`reports/pe_bellhop_controlled_comparison_stage6_fixed_pm_report.md`.

Stage 7 made zero solver calls and reproduced the existing M=50 published
statistics directly from the authoritative CSV: mean delta TL `-0.00747858 dB`,
PE/BH mean powers `1.004510/1.007007`, and circular mean phase
`-0.375345 rad`. It explicitly retains historical point-source R provenance;
the data are neither relabeled X nor empirically corrected. Together with
Stages 4--6, near-equal ensemble power is compatible with realization-level
spatial reflection-model discrepancy, and the nonzero circular phase is not a
remaining convention error. See
`reports/pe_bellhop_controlled_comparison_stage7_historical_m50_interpretation_report.md`.

The controlled Goal is now complete. Its final report is
`reports/pe_bellhop_controlled_comparison_final_report.md`, with status
`PHASE_MECHANISM_IDENTIFIED`. All requirements for `COMMON_LIMIT_CONFIRMED` are
also met, but the more informative final label is used because Stage 5 locates
the residual as spatial reflection-model distortion. The conditional BIE/BEM
branch was not triggered for convention adjudication; it remains necessary only
for a future absolute-accuracy claim in the Region-III PM regime.

## 2026-09-14 Helmholtz BIE third-reference Goal: R0--R2

The R0 formulation audit is complete and frozen in
`reports/pe_bellhop_helmholtz_bie_R0_formulation_audit.md`. The selected
validation-only reference is the two-dimensional sound-soft rough-surface
combined-layer equation built with a Dirichlet half-plane Green function,
`exp(-i*omega*t)`, an inward-to-water normal, and a slow-rise smooth finite
section. The outgoing combined layer is `D_h-i*k*S_h`; an initial R1 smoke
run exposed that the provisional plus sign was refinement-dependent and nearly
singular, so the sign was corrected from the frozen time/radiation/normal
conventions before accepting any result. It is specifically limited to the localized Gaussian angular-spectrum
incidence and graph surfaces required by this Goal. Closed-obstacle BEM,
quasi-periodic/Bloch formulations, and abrupt free-space truncation were
rejected for R1--R5. R1 flat self-validation passes: the accepted 12-ppw
result has boundary residual `9.276e-10`, image-field complex L2 `1.236e-7`,
and phase RMS `9.268e-8 rad`. R2 also passes for the fixed C2-tapered
`A=0.01 m`, `K=0.10 rad/m` benchmark. Its receiver-line BIE uncertainty is
`6.002e-8` complex L2, `4.220e-8 rad` phase RMS, and `3.596e-7 dB` TL RMS.
The physical surface support remains fixed while spatial, close-panel
quadrature, and BIE-window convergence are varied independently. R3 weak
PE--Bellhop--BIE closure also passes on that exact benchmark: PE--BIE has
`E_G=0.003158` and phase RMS `0.003118 rad`, Bellhop--BIE has
`E_G=8.092e-6` and phase RMS `5.886e-6 rad`, and PE--Bellhop has
`E_G=0.003159`. Both native Helmholtz ratios (Bellhop and BIE) are mapped by
the same fixed conjugation into the PE comparison convention; no scalar is
fitted. Bellhop geometry and receiver guards pass. R4 (`A=0.05 m`) gives
PE--BIE/Bellhop--BIE complex errors `0.016090/3.890e-5`; R5 (`A=0.20 m`)
gives `0.080908/9.976e-5`. The ranking is unchanged at every spatial/window
level and both separations exceed `5 U_BIE` by many orders of magnitude.
The final first-round status is `BELLHOP_CLOSER_TO_HELMHOLTZ_REFERENCE`; see
`reports/pe_bellhop_helmholtz_bie_final_report.md`. The R5 strong-height
boundary residual has an explicit `1.1e-8`--`2.9e-8` plateau, while an
independent higher-order receiver field agrees within `3.291e-8`, below its
declared `U_BIE=6.179e-8`. Optional R6 and fixed PM are not needed for this
adjudication and were not run.

## 2026-09-15 PE rough-surface reflection-operator improvement Goal

The BIE workflow was extended with the frozen high-K G0 case
`A=0.02 m, K=0.47 rad/m`. All numerical and Bellhop geometry gates pass;
Model-0 PE/BIE `E_G=0.0229061`, while Bellhop/BIE `E_G=1.36812e-4`.

G1 measures the actual sigma=0.3 m Gaussian angular spectrum:
`kx_rms=2.35702 rad/m`, `theta_rms=8.14401 deg`, and
`<kz>/k=0.9899505`. G2 adds a validation-only `model1_kz_aware` branch to the
1T helper. It uses componentwise `exp(+i*2*kz*eta)`, changes the flat field by
only `1.83e-13`, reproduces legacy Model-0 exactly when selected, and reduces
controlled-case phase RMS by factors `25.0/5.11/1.61/1.22` in
weak/Region-II/strong-height/high-K order. This establishes
`NORMAL_APPROXIMATION_CONFIRMED`.

G3 adds analytic local-normal specular `kzr` and `(kzi+kzr)*eta` phase without
fitted coefficients. It produces only `0.01%--1.5%` further improvement and
fails the high-K gate, despite positive reflected `kz`, zero non-returning
components, and finite fields. The Goal therefore stops before fixed PM with
`NONLOCAL_EFFECT_REQUIRED`. Production PE, PE marching, the production surface
model, Bellhop, and the BIE solver are unchanged. See
`reports/pe_surface_operator_improvement_final_report.md`.

## 2026-09-20 Strict-normal Gaussian-width applicability audit

The validation-only sigma sweep keeps the 4 kHz controlled benchmark, accepted
Helmholtz BIE, X-source Bellhop convention, receiver map, and surface geometry
fixed while changing only the analytic Gaussian width. For the low-slope
`A=0.05 m, K=0.10 rad/m` family, sigma `0.3/0.5/1/2 m` corresponds to
theta95 `15.97/9.55/4.70/2.35 deg`. Legacy Model-0 phase RMS contracts from
`0.01589` to a `0.00460--0.00469 rad` floor, its theta-rms-squared fit has
`R2=0.9881`, and the Model-0/Model-1 complex gap contracts from `0.01558` to
`0.000231`. Total complex/TL error is not monotone and no point passes the
complete all-metric Region-I gate. At high K (`A=0.02 m, K=0.47 rad/m`) and
theta95 `2.35 deg`, both PE variants agree to `1.56e-4` but remain about
`0.0256` from BIE with phase RMS `0.0110 rad`, confirming a separate
nonlocal/spectral-coupling residual. The final status is
`STRICT_NORMAL_APPROXIMATION_PARTIAL`; for this smooth low-K family,
theta95 around `9.5 deg` is supported for the phase criterion, not as a
universal total-field bound. Optional sigma `4 m` is excluded because fixed
Bellhop receiver influence yields sparse `0/0` ratios. Production PE remains
unchanged. See `reports/pe_strict_normal_sigma_sweep_report.md`.

The same audit now includes a strong-height discriminator
`A=0.20 m, K=0.10 rad/m` at sigma `0.3/2 m`. Its Model-0/Model-1 gap falls
from `0.06227` to `0.000924` as theta95 narrows from `15.97` to `2.35 deg`,
while PE--BIE remains `0.08091/0.08708` and phase RMS remains about
`0.075 rad`. Thus the finite-angle normal approximation is removed by a narrow
source, but strong-height phase-screen/nonlocal error remains. This result is
part of the same `STRICT_NORMAL_APPROXIMATION_PARTIAL` report.

## 2026-09-22 Full Kirchhoff sign/convention audit

The independent 2-D Full Kirchhoff validation originally compared its native
`exp(-i*omega*t)` rough/flat ratio directly with `G_BIE`, although the frozen
PE-comparison `G_BIE` had already been mapped as `conj(G_BIE_native)`. Flat
Gate-0 did not expose this mixed representation because both flat quantities
were native. The validation driver now applies the same fixed, single
conjugation to `G_FK_native`; neither the Kirchhoff kernel nor production PE or
BIE physics was changed.

The dedicated `validate_full_kirchhoff_convention_audit.m` confirms that the
authoritative geometry is `z_s=+eta` for the downward-positive coordinate,
with water normal `(-eta',1)/sqrt(1+eta'^2)`, outgoing
`(i/4)H_0^(1)`, and reflected-only BIE output. For weak/low-K, the old mixed
phase residual follows `-4*k*eta` with `0.01837 rad` weighted RMS, while the
corrected phase RMS is `1.806e-4 rad`. Corrected Full-Kirchhoff/BIE complex L2
is `2.535e-4`, `2.570e-4`, and `2.711e-4` for weak/low-K,
strong-height/low-K, and weak/high-K respectively; corresponding phase RMS is
`1.806e-4`, `1.868e-4`, and `2.052e-4 rad`. The mirrored diagnostic geometry
`z_s=-eta`, with its slope, normal, and incident evaluation all recomputed,
retains the large phase mismatch and is not the same physical case. Thus the
previous Full-Kirchhoff rough-case anomaly was a comparison-convention bug,
not evidence that the Kirchhoff approximation fails in these controlled
cases. See `reports/pe_bie_full_kirchhoff_convention_audit_report.md` and the
corrected `reports/pe_bie_full_kirchhoff_validation_report.md`.

## 2026-09-22 Full Kirchhoff deterministic parameter sweep

The validation-only height sweep (`K=0.10 rad/m`, `A=0.005--0.20 m`) and
wavenumber sweep (`A=0.02 m`, `K=0.05--0.70 rad/m`) use the same 4 kHz
Gaussian source, Stage-0 receiver grid/mask, analytic C2-tapered surface,
12-points-per-wavelength Helmholtz BIE, and corrected Full-Kirchhoff field
convention. Full-Kirchhoff/BIE complex L2 remains between `2.536e-4` and
`3.316e-4`; its phase RMS remains between `1.786e-4` and `2.796e-4 rad`.
The 2049--4097 Full-Kirchhoff integration change is at most `3.720e-5` in
complex L2 and `2.669e-5 rad` in phase RMS.

No Full-Kirchhoff failure boundary is observed through nominal `A*K=0.02`
and `A*K^2=0.0098 1/m`. By contrast, PE/BIE complex L2 grows from
`0.00158` to `0.08091` over the height sweep and reaches `0.04843` at
`K=0.70 rad/m`; the low-K growth is phase-led, while high-K develops a large
magnitude residual. This supports the local phase-screen reduction, especially
its missing nonlocal coherent redistribution, as the discrepancy source in
this deterministic envelope. A reduced nonlocal Kirchhoff operator is the
next diagnostic target; SSA is not required by these cases because Full
Kirchhoff already closes to BIE. Production PE, BIE, and surface physics are
unchanged. See `reports/pe_bie_full_kirchhoff_parameter_sweep_report.md`.

## 2026-09-22 Full Kirchhoff reduced-operator Stage 1 diagnostic

The validation-only `validate_full_kirchhoff_reduced_operator_stage1.m`
compares the accepted Full-Kirchhoff field against the PE phase screen for
weak/low-K, strong-height/low-K, and weak/high-K cases. Removing physical
surface height/path phase produces complex L2 errors `0.230`, `1.217`, and
`0.464`; removing the curved normal produces only `2.12e-6`, `1.32e-4`, and
`9.81e-5`; removing the arc-length Jacobian produces `3.13e-7`, `1.25e-4`,
and `2.70e-5`. Physical height/path phase is therefore essential, but PE
already contains its local `2*k*eta` proxy. The remaining discrepancy is the
failure of that local proxy to reproduce complete surface-to-receiver path
phase and nonlocal coherent integration; normal and Jacobian effects are
secondary at the tested slopes.

Receiver-centered truncation of the complete integral requires approximately
`64 m`, `8 m`, and `32 m` respectively to reduce the PE--Full-Kirchhoff error
by half. This support is case-dependent, so it is evidence for a nonlocal
short-range operator rather than a universal tiny convolution kernel. No
reduced production operator was implemented; PE/BIE/Full-Kirchhoff physics
remain unchanged. See
`reports/pe_bie_full_kirchhoff_reduced_operator_stage1_report.md`.

## 2026-09-22 Finite-aperture Reduced Kirchhoff prototype

The validation-only finite-aperture prototype retains the complete accepted
Full-Kirchhoff kernel and truncates only source contributions outside
`|x_r-x_s|<L`. Its three representative-case results are retained in
`reports/pe_bie_full_kirchhoff_finite_aperture_report.md`; the runnable
finite- and adaptive-aperture workflow is now
`validate_reduced_kirchhoff_adaptive_aperture_stage1.m`. No empirical
correction, SSA, or production-model change was made.

## 2026-09-22 Adaptive-aperture Reduced Kirchhoff Stage 1

The validation-only adaptive-aperture database reuses the accepted A sweep
at `K=0.10 rad/m` and K sweep at `A=0.02 m` (eleven unique cases). For each
case it searches `L=4,8,16,32,64,128,Inf m` for the first finite aperture
whose Reduced-Kirchhoff/BIE complex L2 is below `1e-3`. All eleven first pass
at `L=128 m`, where the aperture already covers the entire relevant accepted
surface support on the fixed M99 receiver mask; all `L=Inf` results reproduce
their accepted Full-Kirchhoff fields exactly.

Consequently this strict accuracy gate produces no observed
`L_min(A,K)` variation and does not support an empirical adaptive-aperture
rule or demonstrate an advantage over the fixed `128 m` aperture. This is a
negative but useful result: at this accuracy level the hard spatial truncation
is effectively global. A complete A-K matrix would not resolve that censoring
without a different error target or a different reduction basis. The next
better reduction investigation is a low-rank or locally-stationary kernel
approximation, rather than fitting aperture to height/slope/curvature. See
`reports/pe_bie_full_kirchhoff_adaptive_aperture_stage1_report.md`.

## 2026-09-23 Gaussian source angular-spectrum sensitivity

The validation-only `validate_full_kirchhoff_gaussian_sigma_sensitivity.m`
holds the accepted 4 kHz one-transverse geometry and deterministic sinusoidal
surfaces fixed while sweeping Gaussian source sigma over `0.3, 0.5, 1, 2,
5 m`. It measures propagating theta90/theta95/theta99 quantiles and compares
PE, Full Kirchhoff, and Helmholtz BIE using the frozen Stage-0 M99 and -40 dB
receiver mask. It does not modify PE, BIE, Full Kirchhoff, or the surface
model. The generated report and artifacts are under
`reports/pe_bie_full_kirchhoff_gaussian_sigma_sensitivity_report.md` and
`results/validation/pe_bie_full_kirchhoff_gaussian_sigma_sensitivity/`.

The validation-only `validate_full_kirchhoff_gaussian_sigma_kzaware_comparison.m`
reuses the accepted sigma-sensitivity PE/Full-Kirchhoff/BIE fields and adds the
existing `model1_kz_aware` PE fields. It uses the same sigma-specific M99 plus
pairwise -40 dB mask for all comparisons, and does not change production PE,
Full Kirchhoff, BIE, or surface physics. Missing Model-1 sigma/case files are
computed only with the existing PE validation helper. Artifacts are under
`results/validation/pe_bie_full_kirchhoff_gaussian_sigma_kzaware/`; the report
is `reports/pe_bie_full_kirchhoff_gaussian_sigma_kzaware_comparison_report.md`.

The validation-only `validate_full_kirchhoff_receiver_depth_sensitivity.m`
holds the accepted 4 kHz, sigma=1 m setup fixed and scans receiver depths
1--30 m for the three deterministic sinusoidal cases. It reuses one BIE
boundary density per case for the depth evaluations, holds the accepted
sigma=1 M99 footprint and weights fixed, and performs a nominal/refined
Full-Kirchhoff sampling check at z_rx=1 m. The artifacts are under
`results/validation/pe_bie_full_kirchhoff_receiver_depth_sensitivity/`; the
report is `reports/pe_bie_full_kirchhoff_receiver_depth_sensitivity_report.md`.

The independent numerical audit `validate_receiver_depth_numerical_audit.m`
reuses those saved fields, applies a fixed sigma=1 M99 mask and energy weights,
and adds flat-control, raw/common-convention, angular-spectrum, and z=1 m
convergence diagnostics. Its artifacts are under
`results/validation/pe_bie_receiver_depth_numerical_audit/` and its report is
`reports/pe_bie_receiver_depth_numerical_audit_report.md`. The audit does not
modify production PE, Full Kirchhoff, BIE, or surface physics.
The z_rx=1 m Full-Kirchhoff check is not converged at the nominal sampling
(relative difference about 0.172), so FK/BIE at that depth is reported as a
numerical warning; depths >=2 m have stable FK/BIE control errors.

## 2026-09-23 Full-Kirchhoff operator visualization atlas

The reporting-only `generate_full_kirchhoff_operator_visuals.m` reconstructs
the saved Full-Kirchhoff integrand density for the weak/low-K,
strong-height/low-K, and weak/high-K cases. It does not rerun PE, BIE, or a
parameter sweep. The reconstruction sums to the saved native Full-Kirchhoff
field with zero relative discrepancy at the stored precision. The atlas adds
integrand magnitude/phase maps, single-receiver coherent contributions,
explicit geometry/support/aperture figures, common-mask receiver error plots,
and all-metric height/K sweep plots under
`results/visualization/pe_bie_full_kirchhoff_operator_atlas/`; the interpretation
is recorded in `reports/pe_bie_full_kirchhoff_operator_visualization_report.md`.

## 2026-09-24 Incident/reflection field covariance audit

The validation-only `validate_incident_reflection_field_covariance.m` reuses
the receiver-depth fields at z_rx=2 and 3 m for the three deterministic
sinusoidal cases. It compares the existing PE Model-1 surface incident field
with the same angular-spectrum source used by Full Kirchhoff/BIE, after
applying the fixed PE envelope-to-carrier factor. It then compares reflected
PE/FK/BIE fields in receiver space and transverse kx spectra using the saved
M99 support and weights. Artifacts are under
`results/validation/pe_bie_incident_reflection_covariance/`; the report is
`reports/pe_bie_incident_reflection_covariance_report.md`. No production PE,
Full Kirchhoff, BIE, or surface model was changed.

## 2026-10-08 Full/Fast FK small multi-frequency consistency

Completed the sparse-to-dense frozen canonical PM audit described in the
current-state index. New validation files are the orchestration entrypoint,
`test_spectral_fast_kirchhoff_multifrequency_phase`,
`support/audit_fk_multifrequency_phase`, and
`support/assert_spectral_fast_kirchhoff_multifrequency_artifacts`. The report
and case-specific artifacts preserve the original3/8 failures and every
uniform-rank trial. No reference implementation, public configuration,
exported channel field or communication consumer changed. README navigation
and the technical guide document the native-to-IFFT phase bridge and the
conditional scope of PASS. The source/geometry and C2 embedding are frozen;
no new model, empirical phase fit or per-frequency height normalization was
introduced. The next data-generation stage is permitted under these tested
conditions but was not run.

## 2026-10-08 — FK optimization PPT visualization (presentation only)

- Main research progression paused for the requested presentation task. Render-only
  entry `scripts/reporting/generate_fk_optimization_ppt_visuals.m` reads saved
  prototype/PM pilot/PM boundary/multifrequency evidence; no new scattering solve,
  benchmark, random realization, production physics or public interface change.
- Eleven 16:9 PNG/vector-PDF figures under
  `results/visualization/fk_optimization_ppt/figures/`, including exact saved-source
  Incident reconstruction and native Reflected receiver profiles. Geometry uses
  Tx100/Rx3 research fixture, not the Tx99 Bellhop installation.
- Algorithm fidelity and physical validity remain separate. Fixed canonical
  broadband consistency does not certify a PM population; 16 ms CIR remains
  diagnostic. Report and recommended slide order:
  `reports/fk_optimization_ppt_visualization_report.md`.

## 2026-10-08 — FullFK single receiver versus current FastFK

- Added one validation-only timing entry, `validate_fast_fk_single_receiver_timing`.
  Frozen canonical seed260001, 4 kHz, N4097, Tx100/Rx3, geometry B,
  rank3/12, padding2; the only receiver is existing index2049 (x0,z3).
  FullFK receives scalar coordinates and computes only one integral row;
  unchanged FastFK computes the full grid then extracts that point.
- Shared MATLAB R2025b Update6, one warmup and9 alternating rounds of10 calls:
  median Full single apply0.00098571 s versus Fast grid0.0360479 s;
  Full is36.57 times faster. Preparation plus one application medians are
  0.11155524 /0.21469579 s, so Full is1.925 times faster including setup.
  Relative complex error1.76692919e-4; both outputs reproduce saved fields
  within1e-12 relative tolerance. Frozen reference hashes pass; new entry has
  clean Code Analyzer output. No Full grid solve was run for this benchmark.
- Current whole-grid FastFK advantage does not carry to this single-buoy task;
  prefer direct FullFK single-point evaluation with reusable incident input.
  Historical full-grid timings are documented separately, not mixed into the
  new speed ratios. No physics, public API, communication method or guide
  equation changed; `vertical_comm_guide.md` needs no method revision.
  Report: `reports/fast_fk_single_receiver_timing_report.md`; raw9-round
  CSV/MAT: `results/validation/fast_fk_single_receiver_timing/`.

## 2026-10-18 — Bellhop–BIE / PE–Bellhop–BIE comparison atlas

- Added validation-only `prepare_pe_bellhop_bie_visualization_data`, pure renderer
  `generate_pe_bellhop_bie_comparison_visuals`, and saved-artifact audit helper.
  Fixed 2-D 4 kHz/c1500, Tx100/Rx3, Gaussian sigma .3, matched X line source;
  A .01/.05/.20 m, K .10 rad/m, C2 plateau42/support50. No PM/3-D/Tx99 mixing.
- PE Model-0/BIE reused from active derivative MAT. Missing active R3-R5 originals
  were not recreated; numeric-only crosscheck against results/archive snapshot
  is pointwise identical. Recovered exact internal-wall executable hash matches
  active build manifest and overlay hashes; no compiler installation or cash access.
- Eight new Bellhop jobs (flat plus three rough, 5001/10001 beams), all checkpointed
  and diagnostics passed. Main10001 flat reproduces saved Stage-0 field exactly.
  Main Bellhop–BIE complex L2: 8.091975e-6 /3.889551e-5 /9.976008e-5;
  PE–BIE: .0031575715 /.016089960 /.080908108. Fixed Stage-0 M99 weights,411points.
  Beam-ratio sensitivity:3.51134e-6 /4.89397e-6 /1.11217e-5, not strict error bounds.
- Outside M99,375 zero-flat Bellhop samples give undefined ratios, explicitly NaN
  plus valid masks; raw SHD retained. No interpolation, fitted phase or scaling.
  All411 comparison points valid. G denotes reflected rough/flat ratio, not TL.
- Seven2400x1350 PNG/vector PDF pairs and Chinese per-figure report completed.
  Portable MAT plus packaged renderer rebuild reproduces all7 PNG pixels exactly;
  metric closure0, CSV roundtrip4.8937e-15, native mapping0, all source/binary hashes pass.
  PNG and Poppler-rendered PDF visual QA passed; three new MATLAB files analyzer-clean.
- Result: `results/visualization/pe_bellhop_bie_comparison/20261018_194331_523/`;
  report: `reports/pe_bellhop_bie_comparison_visualization_report.md`.
  README/scripts navigation updated. No production API, frozen solver or physical
  method change; no PE/BIE rerun or communication regression needed. Historical BIE
  refinement checks remain historical; this atlas does not alter PM/3-D validity status.

## 2026-10-19 — Raw PE/Bellhop/BIE x-z field visualization

- User explicitly requested true x-z raw fields, then approved the minimal strong
  A=.20/K=.10/4kHz fixture: existing M99 x411 points, z=.5:.125:6 m (45 levels),
  18495 actual samples per method. No image-row replication or spatial interpolation.
- Added prepare_pe_bellhop_bie_xz_data, evaluate_saved_halfplane_bie_field,
  generate_pe_bellhop_bie_xz_visuals and audit_pe_bellhop_bie_xz (validation/reporting).
  PE evaluates saved surface starting field at new depths; BIE reuses saved4480-node
  density and the frozen off-surface Green formula. New BIE boundary solves0.
- Exactly one new matched internal-wall Bellhop10001-beam grid job completed,
  wall686.2287s, geometry diagnostics passed and coordinate error0. New z3 line
  matches old Bellhop exactly; PE/BIE saved-line closures1.4786e-14/5.2563e-15.
  Saved BIE12/14ppw density probe change1.4232e-8 on44 points, not full-domain certification.
- Four raw-field PNG/PDF pairs: amplitude, phase, phase-hue/amplitude-brightness,
  and raw1/3/6m magnitude profiles. All2400x1350 PNGs, single-page1152x648pt PDFs.
  PDFs embed reviewed PNGs to avoid R2025b tiled/vector export clipping; not vector PDFs.
- Preserve PE envelope plus analytic exp(+ik0(100+z)) carrier restoration, raw BH SHD
  plus its established conjugated-native field, and BIE native field. No fitted
  scaling/phase. PE-restored/BIE flat closure3.855e-13; native flat BH/BIE center
  scale5.77379 and phase .777896rad are diagnostics only, never applied as corrections.
  PE/BIE share amplitude limits; BH has explicit independent native-source limits.
  Old G conjugation does not define a complete absolute-field phase conversion;
  raw phase colors do not constitute a new cross-model phase-error ranking.
- Portable MAT plus packaged script rebuild reproduces4/4 PNG pixels. CSV closure
  7.1089e-16, RGB and carrier formulas0. Four new files analyzer-clean; visual PNG/PDF QA
  passed. Existing three-case atlas preserved. Production APIs/frozen solvers unchanged;
  no communication regression or vertical_comm_guide method update applicable.
- Result: results/visualization/pe_bellhop_bie_xz/20261019_115214_348/;
  report: reports/pe_bellhop_bie_xz_visualization_report.md. README/scripts/current
  comparison report navigation updated. No PM/3-D/absolute pressure or TL claims.

## 2026-10-19 — Presentation summary of the two reflection-field atlases

- Added reports/pe_bellhop_bie_visualization_summary.md and the companion
  pe_bellhop_bie_visualization_data_guide.md. All eleven existing figures are
  explained in plain Chinese, with a suggested presentation sequence, source
  settings, quantitative tables, phase explanation, and an optional spoken summary.
- Keep the z=3m normalized G metrics separate from raw x-z fields. Strong-case
  PE-BIE and BH-BIE weighted complex differences are 8.0908% and 0.0099760%;
  neither is claimed as a full-domain raw-field error. Native source scales and
  phase conventions remain explicit. MAT fields, CSV semantics, masks/weights,
  and non-overwriting saved-data redraw commands are documented.
- Added a repository-local Obsidian reading package at
  results/visualization/pe_bellhop_bie_summary/20261019_presentation/ containing
  the summary, data guide and eleven unchanged PNG copies. No numerical inputs,
  MAT/CSV, executable or source-job outputs are duplicated; original packages
  remain authoritative. No external vault write, solver run or image re-render.
- README and script navigation updated. Production code, physics, APIs, frozen
  solvers and the two original reports/data packages are unchanged; no production
  regression or vertical_comm_guide method update is applicable to this prose task.

## 2026-10-20 — Two-dimensional closeout S0–S1 saved-data audit

- Created and executed a Goal limited to S0–S1 under the plan's section 0.3 and
  the review prompt. The active MATLAB MCP R2025b session ran saved-data-only
  extraction, metrics, artificial-array formula tests, rendering and independent
  redraw; new physical-field computations were zero, including receiver
  evaluation from saved BIE densities. No S2–S6 or production solver work ran.
- Rebuilt all nine historical three-sine/three-pair 411-point M99 G metrics to
  full precision and independently checked source CSV against MAT. The old
  formula is correct for its z=3 m weighted receiver-line G target, but does
  not measure a center point, a complete 3x3 G array or native raw fields.
- Extracted all available center, same-depth three-point G and strong-case
  3x3 native raw samples. Kept missing off-depth flat denominators and
  phase-convention limitations explicit. Seven reviewed figures were rebuilt
  pixel-identically; nine synthetic formula controls passed. The exact source
  hashes, usable-point counts, CSV/MAT, audit, figures and S2 minimal-value
  lower bound are under results/comparisons/two_dimensional_channel_closeout/
  20261020_2330_s0s1/ and explained in
  reports/two_dimensional_channel_closeout_report.md.
- Status S0_S1_COMPLETE_WITH_DECLARED_GAPS. S2–S6 are unexecuted;
  G1_target_PM=OPEN and WP2_formal=false. No public API or physics change,
  commit or push. README, scripts/README and the canonical plan status were
  updated; vertical_comm_guide is unchanged because its methods did not change.

## 2026-10-23 — Two-dimensional closeout S2a actual-entry and source audit

- Executed the S2a-only Goal in the existing MATLAB MCP R2025b session after the
  user restored it. The fixed PE-mode 0/BIE/FullFK source, complex angular
  spectrum, carrier, geometry, normal and reflected-component contracts were
  checked; Bellhop remains source-unresolved and conditional.
- Ran the actual PE physical-height adapter, screen, propagation and ratio
  path for eta_z=0,+0.01,-0.01 m. Independent plane-wave path control differs
  by at most 1.75e-13; saved Gaussian PE/BIE plane fields differ by 4.70e-13.
  Reused verified E12 strong FullFK and computed only bounded flat, weak and
  medium complete two-term FullFK receiver lines with arc_length and infinite
  local radius. No BIE solve/density evaluation, Bellhop job, whole x-z field,
  PM, broadband, FastFK, production API edit, commit or push occurred.
- Saved 4,932 raw line complex-field rows, 12 weighted metrics, G/delta-u,
  increment cross terms, source/version manifest, three inspected PNG/FIG,
  and the 108-cell S2b coverage/cost table. Independent audit checks saved
  fields, 19 current hashes, three direct medium-FullFK reintegrations and
  pixel-identical redraw. Core S2b has 48 available and 60 missing point
  values; BIE reference uncertainty remains unknown. The result is an S2a
  convention/source/baseline closure on declared support, not physics PASS.
- Evidence and reproducible commands: reports/two_dimensional_channel_closeout_report.md
  and results/comparisons/two_dimensional_channel_closeout/20261023_s2a_final/.
  README, scripts/README and the unique plan's status were updated. No method
  change requires a vertical_comm_guide update or production communication
  regression; S2b–S6 remain unexecuted, G1_target_PM=OPEN, WP2_formal=false.

## 2026-10-23 — S2c-B K_HIGH P0 and resource-gated stop

- Created an independent Goal for only K_HIGH. Revalidated 247 formal S2b
  output hashes, 24 code hashes, 203 S2c-A outputs, five S2c-A code hashes,
  22 upstream inputs and the frozen K_HIGH configuration/source/receiver hashes.
  Existing R2025b MATLAB MCP session and repository cwd were confirmed.
- A first zero-physics P0 attempt stopped on MATLAB table column names. The
  preserved final P0 in 20261023_s2c_b_khigh_02 independently reconstructed
  BASE_M and K_HIGH C2 height/slope/curvature, copied no old density, and
  checked BASE_M native 417x3 fields against the formal data. M04a is labeled
  pre-experiment geometry only; saved-MAT redraw is pixel identical and FIG
  reopens. K_HIGH has no acoustic field observation.
- Five required BIE solves were projected conservatively at 1449.8 s; with
  150 s for other physics and 150 s contingency, 1749.8 s fits the 1800 s
  subbudget and leaves about 4041.6 s planning balance. This is an estimate,
  not charged acoustic runtime. A 4480-node job requires at least about
  2.027 GB currently available physical RAM under the 30% reserve rule.
  Initial 1.975 GB and later 1.210 GB readings failed. Auto-review rejected
  a potentially high-cost BIE entry call while the gate failed. No K_HIGH
  PE propagation, BIE density evaluation/solve, FullFK integral, M04b-d or
  other candidate ran. The Goal remains unfinished pending a valid resource
  gate; no accuracy threshold or reference layer was removed.
- Main report, canonical plan and navigation are updated. This validation-only
  P0/preflight work changes no production solver/API, communication method,
  source conventions, S2b/S2c-A artifacts or vertical_comm_guide. No commit,
  push or cash access; S3-S6 remain unexecuted, G1_target_PM=OPEN,
  WP2_formal=false.

## 2026-10-23 — S2c-B K_HIGH representative BIE and cost stop (current)

- This entry supersedes the earlier P0-only RAM state above. A fresh Windows
  physical RAM reading of 2,101,432,320 bytes passed the 30% reserve gate for
  all five planned BIE levels. The existing MATLAB MCP session ran exactly one
  K_HIGH L60/12 ppw/qa96 boundary solve and all 417 frozen receiver positions.
  The saved MAT/JSON receipt records 311.2095523 s, 3840 nodes, linear
  residual 3.40e-15 and off-grid boundary residual 4.13e-9. A tool reply timed
  out near 300 s while MATLAB finished and saved the result; it was not rerun.
- The preregistered five-level cost was recomputed from this actual timing:
  BIE 1857.6201 s plus 150 s other physics and 150 s contingency totals
  2157.6201 s, exceeding the 1800 s S2c-B cap by 357.6201 s. The remaining
  four BIE solves and all K_HIGH PE/FullFK jobs were stopped. No reference
  check or fivefold empirical resolution rule was removed. Shared core charge
  is now 1719.8696573/7200 s; the estimated S3A balance after a hypothetical
  complete K_HIGH remains above 3000 s but does not relax the local cap.
- Independent P1 audit reread saved MAT/JSON, checked all 417 finite complex
  values against exact BASE_M coordinates, solver/config hashes, receipt
  timing and stop arithmetic. Provisional BIE response CSV/support statistics
  are explicitly UNQUALIFIED: new K_HIGH spatial/window/self quadrature
  resolution is absent, so no model-error, mechanism or 5x qualification is
  claimed. M04a is the only figure; M04b-d have no formal input. Status is
  S2C_B_K_HIGH_RESOURCE_REFERENCE_BLOCKED, Goal not marked complete. No
  PHI_PLUS90, S3-S6, Bellhop, communication work, production API edit,
  commit, push or cash access.
