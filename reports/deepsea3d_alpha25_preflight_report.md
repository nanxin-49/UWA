# alpha2.5 高阶Dirichlet参考稳定性预检

**决策：值得进入受控正式G1参考验证；本轮预检通过，alpha2.5物理G1仍INCONCLUSIVE。
全G1目标域OPEN，WP2=false。没有运行FullFK物理比较或正式网格/域扩展。**

固定Tx1500/Rx5、q0、6 kHz、c1500和原dual波矢/phase/C2支撑，仅将(.02,.01) m
放大至(.05,.025) m，最大高度.075 m、k高度界1.885。冻结求解器不变。
低成本配置均为L24 m、h=.125 m、384²点、边界评价M40，包含原网格及半网格真实海面。
alpha2.0仅读取保存结果，未重算既有验证。规则见[预设协议](deepsea3d_alpha25_preflight_protocol.md)。

## 阶数、非配点残差及稳定性

|N|配点真实边界相对残差|半网格残差|wall s|保存复谱载荷MiB|
|---:|---:|---:|---:|---:|
|12|4.070e-4|4.070e-4|1.039|27|
|16|1.244e-6|3.124e-6|1.454|36|
|20|1.554e-9|2.868e-6|1.910|45|

N12未过原1e-4门，因此不能只凭小网格结果选它作参考。12→16的增量变化约8.33e-4，
故补至冻结允许的最高20阶；16→20降为**2.363e-6**，18→20为**9.248e-8**。
N20配点与半网格残差都低于原门槛。半网格残差并未随阶数降到配点尺度，
说明不能仅凭配点闭合或阶数差声称连续边界已充分收敛，必须再检查横向网格。

边界系数L2范数在n3达到.0187704，随后衰减；n20为9.772e-11，仅为峰值的5.206e-9。
最后五阶的系数范数及接收系数幅度均递减，没有观察到高阶尾部重新增长。
源Taylor系数对真实源值的最终误差5.48e-12，未通过拟合消除相位。

接收点阶数求和指标`sum|term_n|/|sum term_n|`为**15.617**。
内部递推各贡献范数和/结果范数的最高指标为**3313.85**：存在显著高阶消项，
但目前双精度计算没有表现出灾难性失稳。二者是数值消项指标，不是连续PDE条件数
或严格误差界。观测使用独立诊断回放；它与冻结FE的逐阶原生接收系数及边界系数
范数相对差均为0。回放仅观测原递推，**不是另一个独立物理oracle**。

证据：[runs.csv](../results/validation/deepsea3d_alpha25_preflight/runs.csv)、
[逐阶复系数](../results/validation/deepsea3d_alpha25_preflight/terms.csv)、
[内部消项指标](../results/validation/deepsea3d_alpha25_preflight/internal_cancellation.csv)、
[预检护栏](../results/validation/deepsea3d_alpha25_preflight/guards.csv)。

## 与alpha2.0比较

相同L24/h.125/M40下，alpha2.0的N12/N16配点残差分别2.560e-5/3.193e-8；
alpha2.5升至4.070e-4/1.244e-6，约增大15.9/39.0倍。
接收求和消项指标从约7.60增至15.62，约翻倍。
这支持“需要更强阶数检查”，不是FullFK物理失效证据，也不证明N16本身失稳：
追加N20后，N16接收增量与N20已相差很小。

alpha2.0的原记录没有同粗网格半网格检查，也没有进程内存快照；不能据此声称
两者峰值内存或完整性能倍率已比较。详细原记录对照在
[alpha2_comparison.csv](../results/validation/deepsea3d_alpha25_preflight/alpha2_comparison.csv)。

## 实际成本、复现与下一步

三个冻结FE预检wall合计**4.403 s**；消项回放**1.542 s**，合计**5.944 s**，
CPU合计**39.297 s**。不含IO、分析及报告。谱载荷最大45 MiB，另有forcing及临时数组。
MATLAB进程前后实际内存快照均约3196.64 MiB（3.12 GiB）；这是既有会话快照，
不是计算峰值，也不表示过程零分配。[成本记录](../results/validation/deepsea3d_alpha25_preflight/cost_summary.json)。

复现：setup后`validate_deepsea3d_alpha25_preflight('init')`，重复`'next'`
至complete（本次三次），然后`'cancellation'`、`'audit'`。
init拒绝覆盖已有checkpoint；next逐次保存，audit只审核/生成表。
新入口和观测helper均通过Code Analyzer；数值回放断言核对全部逐阶系数，
输入断言保证除统一振幅外所有cfg字段与冻结基准一致。原代码及历史数据保留哈希。
MAT仅保存数值/字符串struct；来源、结果hash和验收JSON在本case目录。
工作流沿用MATLAB技能，归属文献见[冻结基准报告](deepsea3d_g1_finite_report.md)末尾。

**建议下一阶段先做L24、h=.0625的N16/N20和半网格检查**，优先观察更高横向
波数对消项与残差的影响；稳定后才继续域扩展和原FullFK反射/增量比较。
本轮没有正式网格/域证据，不能给alpha2.5条件性G1 PASS，1%物理门也尚未检验。
普通FE在当前粗配置下仍可用，不必现在强制转TFE；若细网格出现尾部增长、残差
不降或20阶仍不足，应停止并保留INCONCLUSIVE，再评估TFE或其他独立参考，
不能直接继续加阶、滤波或放宽门槛。
