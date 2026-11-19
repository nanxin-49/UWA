# K1首批回收诊断与N7/N12保存场分析

2026-11-18。**N7独立保存数值审查PASS；N40原生致命退出，K1整体STOPPED且矩阵不完整。** 本轮无新声场求解、无恢复或重启；原回收、R1/R2/K0、冻结05包与生产代码不改。

## 来源与可信范围

回收目录为`results/comparisons/pm_wind_cutoff_k1/k1_pm_20261118_05_retrieved_20261118_175259/`，唯一会话`20261118T093555_1124145`。外层归档SHA `c3903d28c47830c290f16dce69a2a73cd7a65c82c9a3506fc2919876cb9e2318`与用户Fetch保留的远端指纹相同；575项内部SHA、296项冻结输入和原05 manifest全部通过。N7 MODELS→C0–C4→finalization双链/检查点闭合，N40仅MODELS未封存调用，其他四个新面没有计算。见[完整性回执](../results/validation/k1_retrieval_review/20261118_01/integrity_02/receipt.json)。

本次使用已有本地MATLAB MCP R2025b Update6/PID22340，临时单计算线程并恢复原设置，30%本地RAM预留及180秒阶段墙钟检查，七个求解入口哨兵阻止新求解。仅在新目录复制数据并验证SHA后调用冻结独立审计器，不在回收原件内运行或覆盖审计。原R1 N12仅复用包内已通过严格输入/资格匹配的锚点，不将其散射密度用于N7。源为总宽0.10m理想二维条带，4kHz、Tx100m/Rx3m、seed260001、phase contract v2、411个真实位置；U始终标为U_PM。

## N40崩溃：已定位到原生abort，尚未定位具体组件

N40 stderr明确为`double free or corruption (out)`及MATLAB fatal error。只读查询本任务目录与确切失败PID对应的文件，在`/home/misc/matlab_crash_dump.1136071-1`找到2835字节转储；SHA为`6ae51c91f254d21d874b6364c22a7e20428fc72a1fa0409c4ab7d1034e8d94a1`。没有遍历其他用户目录，没有读取其他任务日志，未改远端文件或启动MATLAB。公开摘要已去除Session Key/许可证标识。

转储时间为北京时间17:46:05，PID1136071，Ubuntu24.04.4、glibc2.39、MATLAB R2025b Update1；`Abnormal termination: abort()`，当前线程为`GTP_1`，栈只有glibc的pthread_kill/gsignal/abort和未符号化帧。转储没有M文件栈、PS/FullFK函数名或明确的MATLAB原生模块帧，不能从线程名推断parpool或并行工具箱参与，更不能断言N40海面物理失效或某个Bessel/FFT函数存在缺陷。

**证据支持原生内存管理致命退出；具体触发链仍未确定。** 转储记录abort而监督器最终观察-9，两者都保留，不猜测后续信号来源。终端断开不解释该原生错误，记录也不支持本包资源主动停止：`stop_reason=null`，N40墙钟18.833s，峰RSS1.284GiB，最低可用110.822GiB，memory PSI=0、io最高0.25%，远未触发12GiB/30%RAM/30分钟门。1秒采样不是排除所有瞬时OOM或外部系统事件的证明；本轮未读取系统内核日志。

MODELS调用标记17:45:59写入，约6秒后转储记录abort；该差值不是精确函数执行时间。MODELS源码顺序为源/平面场→六层PS→四层FullFK→common保存；现有空diary/stdout及缺失common无法进一步确定失败在哪一项。N40没有独立BIE计算和可用复场。原Fetch工具只回收白名单目录，未包含该home根目录转储；本轮按确切PID补读，未修改冻结归档器。见[崩溃转储摘要](../results/validation/k1_retrieval_review/20261118_01/owned_crash_summary.json)。

## N7原数值门与独立审计

五层保存密度分别重建接收场、几何、入射rhs、MAT/CSV和组件闭合；所有340/408/408/476/476个平台中点与每层额外19个GL112独立点均按原门检查。综合审计另重建源/解析梯度、BIE层势、FullFK完整双项和PS直接Green表示，0数值失败。原始误差FIG12条曲线逐条对应411点场，2张远端PNG已目视检查。

| 检查 | 本地结果 | 原门 |
|---|---:|---:|
| BIE四项层级敏感性 | 8.3151e-8 / 9.8889e-8 / 8.1932e-8 / 2.5934e-14 | 每项<1e-4 |
| 综合BIE密度场重建相对差 | 2.0540e-14 | <1e-11 |
| 综合独立边界残差 | 6.3284e-8 | <1e-4 |
| FullFK双项重建相对差 | 1.1236e-14 | <1e-11 |
| PS独立Green重建相对差 | 7.8345e-9 | <1e-4 |
| MAT/CSV复场最大差 | 6.7874e-14 Pa | <1e-11 Pa |

PS原六层的窗口/网格/步长差为2.7263e-8、2.3274e-9、2.7719e-15；FullFK原四层最坏差2.8450e-11，均过原门。BIE相邻网格差出现18.9278%相对上升，但两项绝对值约1e-7；按既有r1-bie-reference-v2只作非单调诊断，保留旧失败标识，不改变政策或宣称严格真值。[N7数值回执](../results/validation/k1_retrieval_review/20261118_01/numeric_01/numeric_receipt.json)。

## 固定U_PM=6的N7与历史N12有限比较

比较的源、原低阶系数、相位、C2和411点支持一致；N7为本次独立BIE，N12为严格匹配R1历史锚点。下表是**反射复场相对各自BIE的L2误差**，不是两海面之间的场差。

| N | Kmax (rad/m) | PS复L2 | FullFK复L2 | PS误差峰x (m) | FullFK误差峰x (m) |
|---:|---:|---:|---:|---:|---:|
| 7 | 0.274889 | 10.25099% | 0.0208328% | −40.0391 | 3.1250 |
| 12 | 0.471239 | 15.92350% | 0.0747208% | 21.0938 | −18.5547 |

N7相对N12，PS误差低35.62%，FullFK误差低72.12%；N7的FullFK误差约为PS的1/492。两模型均超过原5倍经验分辨率：N7 PS/FullFK阈值分别1.1405e-6、1.0261e-6（无量纲），分辨成立不等于达到工程容许误差。N7相位误差项占总复误差能量99.9192%/99.9918%，全部411点满足原1% BIE线RMS相位掩码。幅度误差小不代表复场误差小。

按`norm(p_N7-p_N12)/norm(p_N12)`计算，两海面的BIE复场差为1.37894（137.894%）；PS/FullFK对应1.37362/1.37897。相干场差可超过100%，它不表示模型误差大于100%，也不是能量增益。截止改变同时改变表面形状和相干相位，不能把两场视为同一物理问题上的网格收敛。

仅两个嵌套截止、一个Seed，不能证明误差随N单调变化，更不能推断N40结果、风速效应或交互作用；411点不是411个独立样本。见[配对误差表](../results/validation/k1_retrieval_review/20261118_01/numeric_01/paired_N7_N12_metrics.csv)、[两海面场差](../results/validation/k1_retrieval_review/20261118_01/numeric_01/paired_field_changes.csv)和[N7空间/相位分析](../results/validation/k1_retrieval_review/20261118_01/numeric_01/N7_spatial/seed_metrics.csv)。

## 谱带宽、C2与计算成本

N7/N12理论连续PM方差覆盖率分别48.3264%/78.0791%；它们都不是完整PM谱。原始高度RMS为0.14725/0.18641m、坡度RMS为0.02695/0.05398；C2后二阶导数RMS为0.007184/0.017355 m⁻¹。低截止虽更平滑，嵌入也会改变谱：N7/N12嵌入后Kmax以上功率分别8.4736e-4/3.2934e-4m²；该嵌入带外功率不能冒称原始Fourier模式或恢复的PM尾谱。没有RMS重归一化。K0的8192点/dx0.01953125m为表面统计采样，与物理Kmax及PS/BIE求解网格分别记录。

N7 MODELS实测26.709s；五层BIE调用59.648/78.619/62.078/102.049/101.773s，合计404.167s，整案例约544s、RSS峰3.688GiB。最大节点C3基础组装29.267s、近场替换5.958s、LU7.995s、边界残差52.874s，仍符合原预算量级。该已完成案例不能证明N40成本；未授权据此并行或增加资源门。

## 本地图件门失败保留与恢复边界

本次数值/配对分析完成后，本地新绘图的严格零像素/布局检查失败在`r2figure:Layout`：FIG重开数据通过，F02两次PNG像素差均为0，但exportgraphics把figure.Position从1800×1400改成1437×785，轴位置保持。第二次独立布局探针重复观察到相同收缩；显式normal窗口仍复现。像素通过不能覆盖布局门失败；本地候选PNG/FIG保留为失败证据，不宣称通过完整新图件封存。该问题发生在Windows R2025b Update6的保存重导阶段，和Ubuntu N40 MODELS的原生abort是两件事，不能混因。[失败回执](../results/validation/k1_retrieval_review/20261118_01/numeric_01/new_figure_audit/failure.json)、[布局探针](../results/validation/k1_retrieval_review/20261118_01/numeric_01/layout_probe.json)。

首版Python MAT变量盘点遇到scipy无法解析已有MATLAB表序列化，保留失败，改为MAT文件头检查及MATLAB原生whos-file盘点；没有转换或修补原数据。所有N7数值检查仅执行一次；随后只汇总已完成审计回执，没有重复五层审计或求解。Code Analyzer最终0消息。

**当前停止，不直接Start。** N40未封存调用按原防重算规则拒绝再次进入；不得删标记或在原目录冒充恢复。后续应先在不求声场的独立诊断中验证运行环境/原生日志捕获，并设计显式新尝试恢复：保留本次失败、跳过已审计N7、绑定原六面输入和数值门，仅在用户授权重算未封存MODELS后执行。暂不指定更换MATLAB更新、线程、参数或关闭某原生组件为修复，因为没有证据证明其有效。整矩阵尚缺五面独立参考及风速/交互分析；本轮不准备或启动新的正式重算。

本轮新增仅保存审计/完整性诊断入口与报告，没有生产物理或通信接口变更，无需通信回归或修改vertical_comm_guide实现公式；任务表、PROJECT_CONTEXT和进度记录只更新受影响K1范围，其他任务保持。无Git提交推送，无持续监控。
