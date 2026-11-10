# M0：100 m 静态信道相位合同修复与验收

任务日期：2026-11-10。基线：`compare / 3f7071dbfb8e7bc7cf7ff22c1962506bd3a7388d`，本轮改动未提交、未推送。新结果 run ID：`20261110_m0_native_dsp_v2_01`。

本报告的数值验收状态为 **M0_COMPLETE_M1_READY**；最终版本、SHA、资源账、文档与看板验收以[最终检查点](../results/validation/m0_phase_contract/20261110_m0_native_dsp_v2_01/checkpoint_final.json)为闭合凭据。只允许下一任务设计和执行 M1，不表示 M1 已通过。

## 数学合同与实际修改

原生 PE 算子仍是 `exp(i L(kz-k0))`，没有修改 marching、伴随传播、海面 deltaG 或 PM 空间统计。名义直达/反射跨度为97/103 m；`tau_D=97/1500`、`tau_R=103/1500`、`D=diag(exp(-i*2*pi*f*.004))`。

| 表示 | 直达 | 反射 |
|---|---|---|
| legacy_reduced | d | r |
| absolute_physical，时间因子 exp(-i omega t) | d exp(i omega tau_D) | r exp(i omega tau_R) |
| direct_dsp，MATLAB IFFT | conj(d) | D conj(r) |

统一关系为 `H_DSP=conj(P_physical)*exp(i*omega*tau_D)`。总频响及参考频点标量字段名不变：`H_f=H_direct_f+H_reflect_f`，`h_j=H_j(idx_f_ref)`；标量不冒充CIR。

- `src/channel/apply_pe_channel_phase_reference_vertical.m` 保留四参数调用，明定输入为 `native_pe_reduced`；新增可选第五参数。已预共轭的外部 reduced 必须显式指定 `input_representation='dsp_adapted_reduced'`，只允许 DSP 目标。DSP-ready 输入不能再送入此 helper。
- 新 `src/channel/assert_channel_dsp_metadata_vertical.m` 验证 `phase_contract_version='2.0.0'` 及输出表示；公共、cached/adjoint 输出都获得该元数据。CIR/taps/外部通信入口拒绝旧项目合同。裸外部数组继续执行原有 DSP-ready 输入合同，直接使用，不按数值猜测或转换。
- `vertical_wape_propagator` 同步 SSA coherent/scattered/total 诊断字段；物理波场重构中的反射载波补齐 Tx→surface 名义载波，两个可视化消费者同步。原生场数组与传播核不变。
- `contract_kstat_receiver_stats_vertical` 使用 `mu=D*conj(mu_red)`、`C=D*conj(C_red)*D^H`、`P=D*conj(P_red)*D^T`。C/P 的共轭是逐元素共轭。
- 条件模型/库为 schema 3；旧1/2须同时提供明确 `source_representation` 和几何。原生输入用 `D*conj(z)`，历史 carrier-only DSP 输入用 `D^2*conj(z_old)`。复特征向量同步转换；增广实变换为 `[Re A,Im A;Im A,-Re A]`，其中A是D或D²。没有来源或只有总场时拒绝猜测迁移。
- conditional生成、配对秩采样和库采样通过迁移/元数据门；empirical模型增加来源元数据，旧未标记数据不能静默转成taps。MC生产脚本保存新元数据；旧U5/U8/streaming重用入口拒绝不明表示。发布候选检查同步schema3与合同指纹，但本轮没有重跑或覆盖旧RC。
- Full/Fast FK的预共轭桥接和其审计显式声明 adapted 输入，维持单次转换；保留原 carrier-only 负控制。直接构造 `conj(P)*exp(i omega tau_D)` 的二维P22/BIE入口不改相位公式。

Gaussian源保持 `exp(-rho^2/(2*sigma^2))`、连续单位峰值、同一发射面跨频复用，没有Pa或总功率标定。海面升高eta向上为正，物理深度是`-eta`；原生屏为`R0*exp(+2i*k0*eta)`，轴向独立镜像路径为`103+2eta`。非零横向波数下，正常入射屏仍是近似，不声称等价于精确`2kz eta`平移。

## 最小验证与独立审计

仅使用现有本机 MATLAB MCP，R2025b Update 6，逐作业单线程、串行。每作业估计额外峰值200 MiB，不超过启动时可用物理RAM的70%，并检查结束时余量；线程设置退出作业后恢复。没有SSH、新MATLAB进程、BIE/FullFK求解或M1扫描。

| 门 | 配置/独立参考 | 结果 |
|---|---|---|
| A 非实代数 | 4/4062.5/6/8 kHz；非实d/r、多列散射；三个分支、物理↔DSP及防重复转换 | 最大物理↔DSP误差4.286e-13，门1e-12通过 |
| B 高度 | P22，256²/50 m，stepz_lamb=.5，constant plane wave，4062.5 Hz，eta=0/±.01 m | 最大DSP复相对差9.948e-14；高度时延误差约3.90e-18 s；门1e-10/1e-12通过 |
| C Gaussian | sigma=.3 m，P22，256²/50 m，4062.5 Hz；同一离散源独立完整角谱/镜像 | 最大分量物理复相对差3.160e-13；direct_only亦通过 |
| 源跨频 | 独立表达式、单位峰值及4/6/8 kHz初场 | 差0 |
| D 接收/统计 | 64²/25 m、3频点；cached/adjoint/public；另PM16²→PE8² dense/FFT；160个非实、不正规复样本 | public/cached总差1.615e-15；dense/FFT C/P差3.074e-16/6.272e-16；样本C/P差8.983e-15/2.083e-13 |
| D 迁移/采样/SSA | 原生及旧carrier-only两种明确来源；复/增广实EVD、proper/improper及库/配对采样；小型生产SSA | 迁移P最大4.166e-13；SSA闭合3.879e-17，全部通过原门 |
| E CIR | 64个独立DFT bin，62.5 Hz；已知双径，共同原点；FFT/N与IFFT关系 | 4 ms到达误差0；FFT往返4.966e-16；共同移动1.035e-15 |
| E 通信 | 解析复taps及新生产3频点小信道；QPSK128符号、已知CSI、无噪声、不移峰 | BER/SER=0；生产BER=0 |
| F 外部/经验 | 既有四个Fourier/桥接单测；新构造empirical MAT输入、wideband/taps及旧表示拒绝门 | 全部通过；未运行需要独立历史PM fixture的非相位测试 |

所有传播制造控制显式设`alpha_max_np_per_m=0`、uniform、bubble=false、show_figures=false、use_gpu=false、save_mode=rx_only、enforce_1_over_R=false，以分离窗口吸收与源/相位合同；平面波不适用1/R检查。它们只认证匹配离散源的周期域公式，不认证无限域窗口或生产吸收层。

[独立审计](../scripts/validation/audit_m0_static_phase_vertical.m)重开保存的MAT/CSV，不调用传播器或相位helper；直接按镜像几何和显式逆DFT重建物理接收场，重算非实样本C/P和旧DSP迁移、直接DFT重建CIR，并核对相同生产频响上的通信修复。物理重建直达标度误差最大约2.50e-13；独立C/P仍满足1e-12。新IEEE-754 CSV与MAT逐位往返差0。最终审计目录见最终检查点，不用脚本自报PASS代替独立检查。

结果内保存`payload.mat`、`metrics.csv`、stdout/stderr、真实MATLAB异常状态、逐作业检查点及输入/代码/产物SHA。共享MCP不退出MATLAB进程，故`exit_code`明确是本作业异常状态，不冒充操作系统进程退出码。早期四个作业的源绑定见不可变初始清单；后续检查点逐项计算源SHA。最终代码快照与哈希闭合见`code_snapshot/`及`artifacts_sha256.csv`。

## 保留的失败与修复

1. Gaussian数值全部通过后，回执SHA对空日志触发Java重载错误。原payload、失败检查点保留；22 s保守计入预算。修复空字节处理，从原payload恢复回执，没有重复传播。
2. 首次生产通信BER=0.07421875。独立定位为线性卷积截尾后却按循环卷积均衡的边界不匹配。`evaluate_mpsk_channel_ensemble_vertical` 增加并默认使用`receive_window_mode='full_linear'`，保留完整尾部；旧行为可显式选`legacy_head`。修复前后生产H逐元素差0、门不变，修复后BER=0。M0未移峰、拟合相位或删除误码。
3. 首次独立审计要求普通CSV文本与MAT零差，发现writetable显示精度造成最大4.544e-28舍入差。原CSV和失败保留；另存十六进制IEEE-754的无损CSV，仍要求零差。首次hex读取的自动类型推断错误也保留，随后显式string导入通过。没有放宽任何物理/统计门。

## 历史结果和准入边界

- 原生PE传播核、经身份校验的原生缓存和空间统计可继续作为原表示数据；旧物理场可保留其原资格。
- 旧carrier-only direct-DSP、条件模型/统计、CIR/PDP、相位敏感比较和通信结果不能换标签当成修复后结果。分量及来源齐全时可显式迁移并另存；只有总场、无来源或未知schema时停止，不猜测。
- 已经按`conj(P)*exp(i omega tau_D)`构造的外部/P22/BIE DSP数据不再共轭；其既有资格仍限原声源、几何、频率和审计范围。本轮没有将旧P22当作新物理真值。
- 旧报告、MAT、CSV和检查点不覆盖；本轮全部新数值与审计均隔离在新run目录。用户明确撤销main合并要求后，合并已abort，compare原HEAD不变；这不是已同步main的交付。

M1建议**起点候选**：P22、50 m窗口、256²、sigma=.3 m、stepz_lamb=.5、CPU double。M0中关闭吸收仅为制造控制；恢复原sponge=.12/alpha=.15时应先独立评估吸收影响，再按既定M1计划检查窗口、横向网格及步长。该配置尚不是收敛后的推荐生产配置。M1须有新的执行范围/资源账；不自动扩展到宽带传播、粗糙海面物理资格、1500 m或正式S5/S6。

主文档、AGENTS、模块/脚本索引和进度JSON已按此合同更新；看板由正式MATLAB渲染器生成。最终累计MATLAB执行时间、各作业耗时、静态分析、SHA/文档与看板门以最终检查点为准（本轮上限600 s、每作业180 s，包含失败与恢复）。

本轮保守封账 **120/600 s**（包含预检、失败与小额未计时包装余量）；最长单作业保守22 s，未接近180 s门。各作业原始计时在[job_receipts.csv](../results/validation/m0_phase_contract/20261110_m0_native_dsp_v2_01/job_receipts.csv)，数值及最终独立审计门均无未解决失败。历史失败不会从流水中删除。
完整修改清单见[changed_files.txt](../results/validation/m0_phase_contract/20261110_m0_native_dsp_v2_01/changed_files.txt)，可复核源快照、Git差异和产物哈希分别在`code_snapshot/`、`implementation.patch`与`artifacts_sha256.csv`。
