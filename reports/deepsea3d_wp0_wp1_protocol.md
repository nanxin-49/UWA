# WP0 / WP1 最小三维验证协议（2026-10-12）

本协议先于结果冻结。validation only；不改变生产 PE、公共接口或二维基准。
Tx=(0,0,1500) m，Rx=(1505q,0,5) m，q=0/.10，c=1500 m/s。
z向下，水域z>eta；压力释放、均匀无损、静态三维单位Green点源。
`(Delta+k^2)u=-delta`，`u=exp(ikR)/(4piR)`，native采用exp(-i omega t)。
源强单位由此方程定义，非Pa标定或发射功率。DSP=`conj(P)*exp(i omega Ld/c)`，
各分量共享同一直达时刻；没有逐频拟合或分别对齐路径。

G0冻结：解析镜像、源归一化、平面移位、旋转和native/DSP闭合；良态代数1e-10。
梯度有限差分作为独立实现检查，2e-6；不能将有限差分容差当积分物理门限。

三维FullFK保留`-int(G d_n u_inc + d_n G u_inc)dS`，真实法向和面积。
初始海面为紧支撑C2扰动，径向plateau8/support12 m，中心为平面镜面点。
weak：A=.001 m、K=(.4,0) rad/m、phase=.2。
dual：A=(.002,.001) m、K=(.4,0),(.2,.6) rad/m、phase=(.2,.7)。
flat：A=0。它们不代表Hm0=.5/1 m PM总体，后者属于未来阶段。

针对无限平面：解析平面背景加粗糙/平面完整积分密度之差；差密度有紧支撑。
此解析减法不改变Kirchhoff近似，不是经验幅相校准。保存原始有限平面积分误差，
避免用小孔径直接截断无限点源照明。L=8/12/16/20 m中的8 m有意截掉扰动；
只以L>=12比较同一表面。support改成18/24是不同嵌入几何的敏感性，不叫收敛。

6 kHz先做h=.1/.05/.025 m、梯形及独立Gauss4积分、90度旋转、孔径检查。
G1数值子门：反射及非零扰动相对网格变化<1e-3，孔径/Gauss反射变化<1e-3，
旋转<1e-10；同时报告相消处total相对变化，不以total替代reflected验收。
通过后扩展4/8 kHz的h=.05/.025配对检查。这里沿用计划建议的数值1e-3，
不是旧二维FullFK–BIE物理门槛；新计划物理1%仅是筛查建议。

独立弱扰动参考：由Dirichlet条件展开得到b=-2 eta partial_z u_inc，
使用平面半空间Rayleigh核`-2 G'(R)(z_rx-h0)/R`积分传播b。
采用Gauss4 h=.05/.025、L16/20对照，振幅scale1/.5/.25，比较扰动及反射总量，
另以正负振幅中心差检查FK的一阶切线。参考验证范围只有一阶形状导数，
有限振幅截断误差不能伪装为BIE不确定度。无匹配有限粗糙度3D参考则G1保持开放。

运行：`setup_vertical_project; validate_deepsea3d_wp0_wp1('wp0')`，然后依次
`'pilot'`, `'extend'`, `'reference'`, `'spectral'`, `'cancellation'`, `'audit'`。
`runtests('test_deepsea3d_controls')`。结果在`results/validation/deepsea3d_wp0_wp1/`，
每个数值块保存checkpoint；不从MAT加载函数句柄。不进入PM MC、FastFK、统计、SSP或BER。
