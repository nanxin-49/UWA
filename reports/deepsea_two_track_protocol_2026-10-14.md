# 双轨本轮预声明协议（2026-10-14）

依据 [唯一正式计划](deep_sea_channel_research_plan.md)，保持 G1 OPEN/WP2=false。结果独立保存到 `results/validation/deepsea_two_track_20261014/`；入口拒绝覆盖已有 MAT/JSON。

## A：唯一干预与资源门

仅选择右侧二项 Neumann 多项式预条件：对冻结的 L0 预条件系统 Bx=b，令 P=2I−B，解 BP y=b，x=P y。BP=2B−B²。物理 B、边界、DtN、去混叠、精确来源不改；P 不属于物理补偿。若 B 的困难特征值远离 1，该固定多项式可能恶化，不假定必然成功，不切换其它路线。

先 n=8、Chebyshev degree12、L4、h1、b.8、600Hz、.01m 制造输入，仅检查线性、BP 代数、已知向量恢复、左/右残差与原始 PDE 残差关系；不是目标 PM 或 dual 物理认证。仅做一次小规模显式矩阵对照，不改变生产实现。

原 PM 必须匹配旧 pilot L24/h.25/M32/b.8/6kHz/support12、seed261013 原系数。拟定 restart12、最多8周期/96内步，所有 B 调用（含 P 与最终新鲜残差）≤640，单案核≤600s、A 总核≤3600s；原线性门1e-8/真 PDE、边界、DtN门1e-4保留。预条件目标1e-10。若物理残差/线性门没有自洽改善，不扩网格/域。

新增预条件存储以旧 restart12 663.19MiB + 4×18.5625MiB = 737.44MiB 做保守命名数组预估。必须小于实时物理空闲的70%，且不能依赖swap；该估计不是实测峰值。初始空闲650.23MiB，正式求解尚不获资源放行。小规模诊断后再检查一次，若不满足立即 `STOP_MEMORY_PREFLIGHT`，原 PM、dual、quarter 新物理求解均不执行。既有 controls 只读；这种停止不等于“预条件器对 PM 无效”。不等待用户关闭应用，不擅自清理共享 MATLAB 工作区。

若资源门停止，完成匹配三维 BIE/BEM 技术与成本估计，标明没有实测 BIE kernel 时间；不实现大型 BIE，不重跑旧 pilot。峰值可得时记录进程生命周期 PeakWorkingSet64（保守口径）；否则仅记录快照，标为峰值未知。

## B：确定性控制

使用基础 MATLAB、固定数据seed2610141/noise seed2610142；解析 Tx1500/Rx5/c1500，q0/.10，单位点Green。保留绝对 native、conj(native) DSP、几何参考 exp(i2πfτd) 后的 relative DSP，三者不混名。

解析4–8kHz 257/513端点包含频率；DFT 分别取256/512个独立bin（舍弃右端点），中心fc6000，Δf15.625/7.8125Hz，时窗64/128ms，Δt=.25ms。保存端点处理、矩形带限、carrier相位、DFT归一化；零填充单独标为插值，不增加带宽/物理分辨率。验证分量roundtrip/Parseval、6.6667ms/150Hz（q0）、q.10重新计算、中间独立频点与线性插值误差、谱零附近条件数及反射相位延迟。

dual 仅从冻结 frequency_summary.csv 读取已通过的4/6/8kHz，几何Tx(0,150.5,1500)/Rx(0,0,5)、dual振幅(.05,.025)m、wavevectors(.4,0)/(.2,.6)、phase(.2,.7)、C2 plateau8/support12。只做每个离散频点的已知复增益QPSK无噪声回归；禁止把三点插值为“认证宽带dual”。

modem 本轮采用显式Gray QPSK符号、完美CSI定时、解析两路径循环卷积/频域已知信道均衡及足够CP控制。AWGN以Es=1、Eb=.5、E|n|²=1/(2γb)，Eb/N0=0/4/8dB、每点500000bit，与 .5erfc(sqrt(γb)) 作5标准差二项筛查。基础符号级控制不是完整RRC/FSE/pilot接收机验收；这些WP4正式工作保留延期。平面两路径noiseless每独立频bin送QPSK，保留几何增益，不在零点除法；q0平面仍有非零路径幅度差，记录最小增益和条件数。

schema 对 `PM_unvalidated_diagnostic` 一律拒绝 training/statistical_acceptance/PM_BER；本轮所有数据只有 readiness 用途。数据契约存频轴、三类复场、时延、seed、版本哈希、表面定义、证据scope和成本。可选2D→3D限延期，不占本轮预算，不增加G1门。

## 审计

保存基线代码/全部既有src、scripts与results/validation文件哈希，最终逐一比对；新增文件单独记录，不改历史数据。运行新入口、适用冻结小规模回归、Code Analyzer和引用/结构/版本检查，输出两份中文报告与状态。A资源阻塞如实标注，不能写成已进行原PM残差比较。

## 正式检查点前资源复查
小规模诊断后物理空闲升至约1517.7MiB，70%上限1062.39MiB，大于737.44MiB预估。允许按 flat→dual→quarter→原PM 顺序执行新右预条件器控制，逐案仍检查70%物理空闲门。统一restart12/8周期/96内步/640 B调用，求解回调560s，预留40s独立残差评价，单案600s、总3600s；不新增其它预条件器或参数扫。旧控制/旧失败只读。具体新代码哈希已写 run_manifest.json。
