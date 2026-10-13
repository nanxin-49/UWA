# 三维压力释放 TFE 最小原型协议

冻结输入只读加载PM feasibility.mat中的cfg、sea系数，不调用海面生成器；
dual加载azimuth.mat中已认证参考与cfg。6kHz、原生exp(-iωt)、z向下、单位点Green不变。
原型独立解 mapped Helmholtz，不调用FullFK/普通FE求解器。

对粗糙修正v=u-u_flat，z=s+t(1-s/b)eta，a=1-t eta/b，s∈[0,b]；
v在点源之外满足Helmholtz，v(s=0)=-u_flat(t eta)。
上边界物理v_z=T v，T符号 i sqrt(k²-kx²-ky²)，Im sqrt≥0，保留倏逝与掠射模态。
展平后v_s-Tv=-t eta/b T v。横向周期仅用于紧支撑修正，不将点源周期化或用周期背景代替无限平面镜像；
有限数值域内周期修正的重复影响仍须扩域验收，纵向DtN不证明横向无限域可靠。

链式法则得到 a²(Δ+k²)=L0+t L1+t²L2；每阶解同一Fourier-Chebyshev Robin ODE。
采用底面Dirichlet消元和内部特征分解；接近内部Dirichlet极点的列用完整Robin矩阵直接解。
基底/线性残差/递推消项指标保留。预算命名数组估计上限900MiB，非进程峰值。
纵向薄层首选b=.8m（在PM support12最高海面以下、水中），Rx5位于层外，接收通过DtN场传播。

实现验证：平面镜像；传播/倏逝/掠射垂直模态；非零体源对照完整耦合ODE；
解析平面波展平系数验证L0/L1/L2；已G1 PASS dual的冻结高阶独立参考交叉验证。
dual预设粗L24/h.25/M32/N12，随后h.125、M40、b1.2、L32；
可按明确失败追加一次h.0625小域诊断，但不可将其替代正式扩域。

PM先support12，N16/24/32，横向h.25/纵向M32；失败可用振幅.25作标记诊断，原海面系数不变。
仅在原振幅数值初门通过后做h.125、M40、L32、b1.2验证；最高N40。
support24只在support12参考可信后进入，不默认执行FullFK–PM。

原参考门1e-4：阶数/横向/纵向/计算域/透明平面位置，以及真实边界与半网格、
独立重算t=1物理Helmholtz残差与半网格残差、真实物理DtN残差；
ODE归一残差<1e-8，源项真实边界误差不能被接收点平稳替代。
邻近低阶压力、体系数尾部和消项同列，边界按真实图面精确点源函数评价。
平面背景是解析已知控制，不单独证明非平面方程实现。
如原参考门或内存预算无法满足则INCONCLUSIVE，建议匹配3D BIE小原型，不增大稠密矩阵。
全目标G1 OPEN、WP2=false，无幅相拟合，无生产/冻结解算器/历史结果修改。

## 结果后有界失效诊断补充

原振幅PM的N32/40体PDE与TBC失败，M32→40不改变接收压力，b=.8→1.2仍失败；
quarter幅度诊断通过局部门。为了区分横向截断/周期重复与阶数失效，追加仅一次
原振幅N32的h.125和一次L32/h.25。900MiB命名数组预算不变，仍不做support24。
这两项是失败原因诊断，不能将小压力差补成参考PASS；不追加N>40或更大域。
dual原N12→16变化未过阶数门，另补一次同粗网格N20检验16→20及尾部；
其余已保存的独立横向/纵向/域/透明平面变体保持N16，不删除最初N12失败。

## 复现命令

setup_vertical_project；validate_deepsea3d_tfe_prototype('init')；
按scripts/README中列出的16个job逐一执行validate_deepsea3d_tfe_prototype('run',job)；
最后audit_deepsea3d_tfe_artifacts。'status'只读，init和重复job拒绝覆盖。
数据存于results/validation/deepsea3d_tfe_prototype/，所有路径从脚本显式根目录解析，不依赖cd。
