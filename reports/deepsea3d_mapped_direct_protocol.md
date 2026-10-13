# 三维t=1展平Helmholtz直接法预设协议

固定原PM系数/seed261013、Hm0=.5父面、6kHz、Tx1500/Rx5、q=.10、support12及原生约定。
加载冻结prototype.mat，不生成海面，不调用Taylor源项/展开，也不修改既有求解器。

连续方程采用已核查的z=s+(1-s/b)eta、a=1-eta/b，a²(Delta+k²)=L0+L1+L2。
粗糙修正V底边=-u_flat(eta)精确点源函数；物理上边v_z=T v，展平Vs−TV=-eta TV/b。
源在薄层外，无体源；Rx由出射固定平面传播。a>0时连续方程等价；有限Galerkin/周期域仍是数值近似。

离散同一对称Fourier空间、2n三因子去混叠与M32 Chebyshev；L24/h.25/b.8。
记S0为平面Robin模态逆，求B V=rhs，B V=V−S0(−[L1+L2]V,0,−eta TV/b)，
rhs=S0(0,g,0)。矩阵自由，不建稠密矩阵。与TFE相同算子，但不对t做Taylor求和。

GMRES restart8、目标预条件相对残差1e-10；平面/dual/quarter先行，控制必须线性及真残差过门才原PM。
controls最多8周期/64内迭代、100次算子调用、Krylov180s；PM最多12周期/96内迭代、160调用、300s。
逐周期保存最后有效解；中途预算终止保留上一周期解，独立评价另计成本。不追加更大restart/迭代。
所有新核预算900s，开始下个job前审查剩余预算；命名工作数组上限900MiB、物理可用RAM的80%，不靠交换区扩算。
基网格restart8加独立4倍评价估计588.94MiB，一向量18.56MiB；估计/快照均非峰值。

三层证据分别验收：fresh B残差；独立重算离散体/Dirichlet/DtN；独立3n/4n真连续几何PDE/边界/DtN。
线性门1e-8，物理门1e-4不变。真边界使用精确点源函数，不以预条件GMRES残差取代物理认证。
保存原生接收复压力并与冻结dual独立参考交叉检查；quarter仅缩放控制。
若原PM初始过门，再评估收敛验证的价值和预算，不自动扩算；失败定位并停止。
无FullFK–PM/support24/MC/FastFK/WP2；即使初检通过，PM INCONCLUSIVE/G1 OPEN/WP2=false。

入口validate_deepsea3d_mapped_direct('flat'/'dual'/'quarter'/'pm')，setup后从任何cwd可运行，拒绝覆盖。
结果results/validation/deepsea3d_mapped_direct/，MAT完整场、周期/压力CSV、协议与审计。

## 原PM失败后的有界restart诊断

restart8在96内迭代后B残差仍1.82e-3，但连续下降且真实PDE显著改善。
仅追加pm_restart12：restart12、8周期，仍最多96内迭代、160调用、Krylov300s。
命名数组估计663.19MiB，仍受900MiB及80%物理可用RAM限制；总核预算900s不变。
这是restart信息损失/预条件效果对照，不扩迭代总数、物理域或数值网格，也不放宽门。
完成后无论结果如何均停止原PM扩算，再按实际残差决定下一轮参考路线。

## 最终停止状态

controls通过；原PM两restart在96内步仍未关闭线性/物理门，不进行网格/纵向/域精化。
最后audit_deepsea3d_mapped_direct重建三层分类并给出continuous_map_certificate.json。
下一轮若继续应是单独有界的预条件器开发，不在当前结果目录继续加周期。
