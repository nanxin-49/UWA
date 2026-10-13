# 三维 PM 矩阵自由展平 Helmholtz 直接法（2026-10-13）

**原型通过平面、认证dual及1/4振幅PM控制；原振幅PM在既定96内迭代上限内仍未线性收敛。直接法明显改善了当前Taylor有限和，但尚不能提供可信PM独立参考。PM INCONCLUSIVE、G1 OPEN、WP2=false。**

## 1. 连续方程、源场与实际离散

只读加载冻结`deepsea3d_tfe_prototype/prototype.mat`：seed261013、父面Hm0=.5m、6kHz、
Tx=(0,150.5,1500)、Rx=(0,0,5)、q=.10、C2 plateau8/support12、c1500、单位三维Green与exp(-iωt)不变。
没有重新生成海面。quarter仅是同一系数乘.25的控制；平面/dual沿用已建立的控制。

对粗糙修正v=u−u_flat，令z=s+(1−s/b)eta、a=1−eta/b，s∈[0,b]。精确链式法则给出

\[
a^2(\Delta_{physical}+k^2)V=(L_0+L_1+L_2)V=0,
\quad V(0)=-u_{flat}(\eta),\quad V_s-TV=-\eta TV/b\ (s=b),
\]

L1/L2按[已验证去混叠定义](deepsea3d_tfe_dealiasing_report.md)使用，不做t展开。
T=i√(k²−kx²−ky²)，取Im√≥0，物理v_z=Tv转换为上式；eta在T左侧相乘。
粗糙修正中点源奇性与u_flat的同一点源抵消；薄层无体源，固定平面以外的修正场出射。
底边直接计算真实海面上的点源/镜像函数，**没有Taylor源项系数或幅相拟合**。

映射等价要求a>0。本轮用精确C2导数界
|∇eta|≤Σ|A||kappa|+1.875Σ|A|/(support−plateau)，以及h.125几何采样的最近点距离界，
得到原PM连续高度上界≈.5423m、连续Jacobian下界≈.3221>0（b=.8m）。
这是含采样间隙的连续上界，不只是采样最大值；Rx/Tx位于薄层外，背景Tx镜像也在水域外。
该几何证据不构成PDE误差界。

实际仍用L24、h.25、Chebyshev M32、b.8和同一对称Fourier空间、2n三因子去混叠。
连续方程的等价性不等于有限Galerkin/周期域已收敛；纵向DtN也不能认证横向周期重复或无限PM支撑。
本轮没有更改物理支撑、网格或周期域。

## 2. 矩阵自由求解与预算

记S0为已有平面L0 Robin模态逆、Q=L1+L2。求解

\[
B V=V-S_0(-QV,0,-\eta TV/b)=S_0(0,g,0),\quad g=-u_{flat}(\eta).
\]

这是原离散方程的预条件形式，不是Taylor求和。每次应用B复用去混叠算子与L0模态解，
不组装稠密体矩阵。GMRES按[MATLAB官方定义](https://www.mathworks.com/help/matlab/ref/gmres.html)
每周期restart后保存有效解，重算新鲜B残差；同时单独验收未预条件体/Dirichlet/DtN方程。

预设restart8；controls最多64内迭代/100算子调用/Krylov180s，原PM最多96/160/300s。
目标B相对残差1e-10，最终线性门1e-8、独立物理门1e-4不变。
原PM失败后只补**restart12、仍96内迭代**的有界对照，没有增加迭代总数。
命名工作数组估计restart8/12为588.94/663.19MiB，受900MiB及80%物理空闲RAM限制；
所有新核预算900s，不用交换区扩算。周期中断通过MATLAB异常cause链识别预算原因，保留上一有效周期解。

## 3. 控制、原PM及三层残差

四项测试通过：平面预条件算子恒等性；预条件形式与原始体/边界残差逆的一致性；
矩阵自由GMRES恢复制造解；预算回调中断及有效解保存。测试包含非平面算子，平面零修正不是唯一依据。

平面精确复现镜像背景。dual完成4周期，quarter完成3周期；二者线性与3n/4n独立真残差均通过。
dual原生接收增量与冻结认证Dirichlet参考差**1.0935e-5**，过原1e-4门；不是新FullFK–PM认证。

| 案例 | 新鲜B相对残差 | 未预条件离散体残差 | 独立4倍网格体PDE | 真实边界 | 真实DtN |
|---|---:|---:|---:|---:|---:|
| dual |5.095e-11|3.167e-11|5.706e-6|2.224e-5|1.043e-6|
| quarter PM |8.039e-11|3.723e-11|6.032e-6|2.314e-5|1.143e-6|
| 原PM restart8 |1.822e-3|1.155e-3|1.901e-3|1.04033e-4|1.690e-3|
| 原PM restart12 |1.359e-3|9.212e-4|1.329e-3|1.04033e-4|1.273e-3|

原PM两个分支均达到96内迭代后退出（MAX_CYCLES，GMRES flag1），没有达到1e-10目标或1e-8线性门。
restart12仅有限改善，没有关闭Krylov瓶颈。3n→4n体PDE变化约1.32e-9，失败不是评价网格太粗。
restart12投影回保留空间的物理残差约1.3288e-3，几何算子差约1.41e-5；
未完成的离散线性求解是当前主要障碍，不能先将物理误差全部归因于分辨率。

同时存在粗网格Dirichlet表示问题：离散底边误差≈9e-15，真实连续边界仍1.04033e-4，
略高于门。这是采样/投影后的源边界与连续图面的差异，增加Krylov精度本身不能解决。
纵向残差仍在Chebyshev节点检查；正式横向/纵向/域收敛尚未执行。

原生接收增量（m⁻¹、单位Green归一化）为：

- restart8：4.8999620666e-5 − 6.1494355588e-5 i；
- restart12：4.8998849101e-5 − 6.1491588047e-5 i。

两者仅差3.6541e-5相对量，但体/透明边界仍失败，**接收压力稳定不能使参考合格**。
反射/总场原生复数、周期收敛曲线及完整场均保留。

## 4. 是否优于Taylor，以及下一步

相较同网格去混叠Taylor N32，restart12的独立体残差降低约**357.6倍**、DtN残差降低约**256.4倍**。
这是残差指标改善，不是相对于真实解的误差认证，也不是公平运行时加速倍数；两轮检查点数量不同。
直接法避免了当前有限Taylor和的主要问题，具有继续开发价值，但本轮仍没有合格原PM离散解。

**不启动正式网格/纵向/域收敛验证，不继续提高当前GMRES迭代预算。**
下一步优先一次有明确预算的预条件器诊断：改进块/粗空间预条件或更合适的模态预条件，
保持原物理算子无损、频率/海面不变，先在同网格关闭1e-8线性门；不得把预条件中的近似改成物理补偿。
restart8→12效果有限，单纯扩大restart/周期数缺乏足够成本依据。

暂不直接放弃这条展平直接法转大型BIE；若改进预条件仍无法在预算内获得合格离散解，
再转匹配三维BIE最小原型。即使线性解通过，也需处理真实边界表示及独立数值收敛，
然后才能考虑物理对照。最小BIE仍须先核查像源位置、近奇异求积及远区闭合，不能用开放面片替代无限界面。

## 5. 复现、成本与状态

[协议](deepsea3d_mapped_direct_protocol.md)。`setup_vertical_project`后执行
`validate_deepsea3d_mapped_direct`的flat、dual、quarter、pm、pm_restart12五个job；
最后`audit_deepsea3d_mapped_direct`。新job拒绝覆盖；从任何已setup的cwd按显式根路径解析输出。

- 新验证代码：`deepsea3d_mapped_action`、`deepsea3d_mapped_l0_solve`、
  `deepsea3d_mapped_discrete_residual`、`solve_deepsea3d_mapped_direct`；
  另有连续map证书、数据审计和测试。旧去混叠/模态/独立评价代码原样复用，不修改冻结求解器。
- [数据目录](../results/validation/deepsea3d_mapped_direct/)保存完整W、原生压力、周期/摘要CSV、几何证书、成本/审计JSON及冻结哈希。
- 5个保存案例核wall **464.135s**、CPU **2654.547s**；含plan、源场、GMRES及离散/独立残差评价，
  不含IO、测试、审计和报告。CPU是多线程累计。平面/dual/quarter分别3.340/63.291/65.030s，原PM两分支168.631/163.843s。
- 命名工作数组估计最大663.19MiB，进程前后快照最大2675.62MiB；均非峰值监测。
  运行前物理空闲约1025.55MiB；未扩网格或依赖交换区。
- 最终4/4测试、数据审计及静态检查通过；冻结312文件SHA256不变。
  项目/方法/脚本/计划导航已更新，没有修改生产接口或历史基准，未提交或推送。

无FullFK–PM、support24、MC、FastFK或WP2。**PM INCONCLUSIVE、G1 OPEN、WP2=false。**
科学证据审查流程沿用 Kassis et al. (2026), [Scientific Agent Skills](https://doi.org/10.48550/arXiv.2609.00065)，仅作工作流出处。
