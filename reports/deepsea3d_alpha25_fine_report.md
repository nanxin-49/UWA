# alpha2.5 细网格高阶Dirichlet参考验证

**结论：L24、h=.0625 的N16/N20细网格参考检查通过，普通FE在该配置下未见失稳。
建议下一阶段受控扩域，稳定后再作FullFK物理对照；本轮不执行。
alpha2.5物理G1仍INCONCLUSIVE，全目标域OPEN，WP2=false。**

沿用冻结预检的Tx1500/Rx5、q0、6 kHz、c1500、dual波矢/phase/C2支撑和原生约定；
振幅固定(.05,.025) m。只计算h=.0625、N16/N20、M40，旧h=.125结果只读。
原FE、诊断回放helper及生产/二维/历史基准不变，没有滤波、拟合或放宽门限。
[预设协议](deepsea3d_alpha25_fine_protocol.md)保留参考1e-4门；未检验FullFK的1%物理门。

## 数值收敛与真实边界

|N|粗网格配点残差|细网格配点残差|粗网格半网格残差|细网格半网格残差|粗→细复增量变化|
|---:|---:|---:|---:|---:|---:|
|16|1.244e-6|1.244e-6|3.124e-6|1.281e-6|1.773e-6|
|20|1.554e-9|1.554e-9|2.868e-6|3.060e-7|1.773e-6|

细网格N16→N20接收复增量变化**2.363e-6**；N18→N20部分和变化**9.248e-8**。
N20增量为`9.95271228e-5 - i 1.51696165e-5`（原生单位Green响应m^-1）。
半网格检查在独立采样的真实海面上评价同一Fourier场，并非递推配点闭合。
N20半网格残差比粗网格降低约9.38倍，全部已有门槛通过；没有仅凭点压力相近判定稳定。
N16配点残差不随加密明显下降，说明该量主要受阶数截断限制，N20补强仍有意义。

两个N的最后五阶接收系数幅度及边界系数范数均递减。N20原始边界系数范数
由粗网格9.772e-11增至细网格1.954e-10，**这是采样数增加造成的范数变化**。
在二维平面上按h转换近似连续L2后，两者为1.22149080e-11/1.22149078e-11，
比例.999999987，未显示尾部增长。范数一致不等于整个复边界场逐值误差已证明。

源Taylor系数误差最大2.02e-11。M40有限Fourier指数评价尾项绝对上界约3.13e-20 m^-1，
只控制指数求值截断，不含浮点误差，**不是FE或连续PDE误差界**。
没有超出20阶，也没有认证更细网格、更大振幅或普通FE一般收敛半径。

证据：[runs.csv](../results/validation/deepsea3d_alpha25_fine/runs.csv)、
[guards.csv](../results/validation/deepsea3d_alpha25_fine/guards.csv)、
[逐阶复系数及h归一化范数](../results/validation/deepsea3d_alpha25_fine/terms.csv)、
[评价检查](../results/validation/deepsea3d_alpha25_fine/evaluation_checks.csv)。

## 消项、内存与实际成本

|N|内部贡献范数和/结果范数最大指标|接收点阶数求和指标|FE wall s|诊断回放wall s|保存谱载荷MiB|
|---:|---:|---:|---:|---:|---:|
|16|654.57|15.6173|7.069|4.015|144|
|20|3313.85|15.6173|8.414|6.640|180|

N20内部指标与粗网格约3313.85一致；存在高阶消项，但未观察到加密后的灾难性放大。
回放的逐阶原生压力和边界系数范数与冻结FE相对差均为0。它仅观测原递推，
不是第二个物理oracle；消项指标也不等于连续问题条件数或严格浮点误差界。

两次FE wall合计15.484 s，两次回放10.655 s，总保存计算**26.138 s**；CPU合计152.828 s。
多线程CPU可大于wall。成本不含IO、分析和文档；旧粗网格及历史求解不重复计入。
最大谱载荷180 MiB，比对应粗网格45 MiB增加4倍；还有forcing/临时数组，非进程峰值。
FE进程前后快照约3193.63–3193.67 MiB，回放快照约3190.62–3193.68 MiB；
这是既有MATLAB会话内存，不表示计算零分配，也不能用于宣称峰值内存已测量。
详见[cost_summary.json](../results/validation/deepsea3d_alpha25_fine/cost_summary.json)、
[observers.csv](../results/validation/deepsea3d_alpha25_fine/observers.csv)、
[内部逐阶消项](../results/validation/deepsea3d_alpha25_fine/internal_cancellation.csv)。

## 复现和下一阶段

setup后依次执行：

```matlab
validate_deepsea3d_alpha25_fine('init');
validate_deepsea3d_alpha25_fine('next'); % N16
validate_deepsea3d_alpha25_fine('next'); % N20
validate_deepsea3d_alpha25_fine('observe'); % N16回放
validate_deepsea3d_alpha25_fine('observe'); % N20回放
validate_deepsea3d_alpha25_fine('audit');
```

init拒绝覆盖已有状态，next/observe各一次运行并checkpoint，audit只核查和生成表。
新入口通过Code Analyzer；回放断言逐阶一致，cfg断言全部物理字段等同冻结预检。
结果在`results/validation/deepsea3d_alpha25_fine/`，纯数值MAT、CSV、验收JSON及哈希保留。
工作流沿用MATLAB技能，归属文献见[冻结基准报告](deepsea3d_g1_finite_report.md)末尾。

**建议继续，但先验证域扩展，不能直接把本轮细网格PASS改为物理G1 PASS。**
可逐个检查L32/L48的N20、h=.0625、真实半网格边界及成本；参考稳定后，再和
同一alpha2.5原FullFK作反射场/粗糙增量原生比较。若扩域出现尾部增长或残差失败，
应保留INCONCLUSIVE并评估稳定参考，不能加阶超过上限或滤波制造通过。
目前没有必要强制切换TFE；全目标域G1仍OPEN，WP2不放行。
