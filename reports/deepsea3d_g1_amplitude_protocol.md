# 三维G1单变量振幅扩展预设协议

基准为冻结的`deepsea3d_g1_finite`；比例1.0直接复用其原始数值状态及验收证据，
不重跑WP0、旧WP1或原G1。只将(.02,.01) m乘共同alpha；Tx1500/Rx5、q0、
6 kHz、c1500、dual波矢/phase、C2 plateau8/support12与所有相位/源约定不变。
先alpha=.5/1.5；两者参考稳定且单次参考wall<=60 s才安排alpha2。
60 s是本轮单次成本护栏，不是物理/精度门；已完成较慢运行不重启。

冻结原FullFK、Dirichlet FE、递推、source/flat/surface/phase helpers，只增加编排入口。
每个新增点先L24、h=.125、N12/16；不满足阶数变化/真实边界门时增加N20。
当前冻结参考最多N20；N20仍不足则INCONCLUSIVE，不能修改冻结参考或提高门限。
所选阶数各自检查L24 h=.0625、L48 h=.125和L48 h=.0625；后者是联合精化主参考。
实际粗糙面原网格和半网格残差均检查；M32/40评价变化及有限Fourier指数尾界单列。
最终点还比较同一次FE输出的N与N-2部分和（同网格严格等价的阶数截断）。
若细网格护栏不满足，可追加N20同配置复查，最多一个补强对；不无限调参。

沿用原协议：阶数/网格/域/真实边界残差各<1e-4；FullFK自身网格/独立求积<1e-3；
可信参考下反射复误差和非零增量复误差各<1%。原始相位差另列。
加权/概率/宽带物理门不在本轮范围；不以total或平面背景掩盖散射增量误差。
未过参考门一律INCONCLUSIVE，即使FK对它偏差很大也不判物理FAIL。
通过参考但FK越过1%则条件性PHYSICAL_FAIL；无幅相拟合或经验校准。

网格、域变化是经验收敛证据，不是严格PDE余项界。真实海面统计在同一固定网格
计算RMS高度、RMS/最大坡度；相位不对齐，系数相位不重新采样。
每个可信新增振幅只补做其FullFK h=.05/.025及Gauss4确认，未复跑旧面数值验证。
全部G1目标域OPEN，WP2=false；不得从四个离散振幅点声明连续振幅区间定理。

复现：`setup_vertical_project; validate_deepsea3d_g1_amplitude('init')`，
重复`'next'`，每次最多一个参考解或一个FK检查块，保存checkpoint；入口显示剩余状态。
全部安排完成后`'audit'`。结果在`results/validation/deepsea3d_g1_amplitude/`。
init拒绝覆盖已有状态；next可从已保存状态续跑。真实测量成本与基准历史成本分开。
