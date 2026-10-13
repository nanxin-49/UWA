# 三维 FullFK 单案例有限粗糙度 G1 审查

2026-10-12，MATLAB R2025b Update 6，共享MCP会话。

**结果：当前单案例条件性 G1 PASS；整个目标域的G1仍OPEN，不放行WP2。**
采用独立16阶Dirichlet场展开，补上了先前只有一阶参考的有限振幅缺口。
这个结果不覆盖PM总体、宽带、其他几何或更强/更高K海面。

## 方案与预检

|候选|独立性和主要限制|本轮决定|
|---|---|---|
|匹配3D半空间BIE/BEM|独立解Helmholtz边界；需要无限边界处理、奇异求积和收敛实现；旧32m面lambda/8密矩阵类比约16TiB|不启动稠密求解|
|高阶Dirichlet FE|直接由真实压力释放条件展开，不使用Kirchhoff局部近似；可能有高阶消项和导数放大|先小网格预检，然后逐阶/网格/域精化|
|TFE/稳定高阶算法|独立PDE参考，能改善普通FE稳定性；实现和验证成本更高|若FE护栏失败再用|
|单独一阶形状导数|忽略显著有限振幅项|不作为本轮oracle，不重复旧实验|

半空间BIE数学路线见[Chandler-Wilde等](https://doi.org/10.1137/050635262)；
[Bempp文档](https://bempp.com/handbook/api/boundary_operators.html)提供层算子/奇异求积/FMM，
不自动提供本问题的无限粗糙压力释放边界。普通FE的消项风险及稳定替代见
[Nicholls、Reitich (2004)](https://homepages.math.uic.edu/~nicholls/papers/Final/shapedefia.josa.pdf)。
本轮实现经典FE，没有将它标为TFE或声称自动继承稳定TFE的定理。

预检L24 m、h=.25/.125 m、N6，单次wall .734/.370 s，保存谱载荷3.38/13.5 MiB。
结果网格变化小，但4→6阶增量变化3.30%，边界残差8.25e-4，**N6不能作最终参考**。
据此先检查8/10/12/14阶，再实施正式加密；未先启动大规模BIE。
证据：[preflight.csv](../results/validation/deepsea3d_g1_finite/preflight.csv)。

## 同一真实海面及独立方程

Tx=(0,0,1500)、Rx=(0,0,5) m，c1500，6 kHz。海面沿用原dual的波矢
(.4,0)/(.2,.6) rad/m、phase(.2,.7)、C2 plateau8/support12 m，中心(0,0)，
仅将振幅固定为(.02,.01) m，理论高度界.03 m、k高度界.754。
这是有两个水平坐标、非共线波矢的有限振幅紧支撑海面，不是无限PM海面。

保持原native exp(-i omega t)、单位`exp(ikR)/(4piR)`点源和水域z>eta。
令u0为点源加平面镜像总场，粗糙增量满足
`delta_u(x,y,t eta)=-u0(x,y,t eta)`，每阶在平面向水侧辐射：

\[
b_n=-[t^n]u_0(t\eta)-\sum_{j=1}^{n-1}\frac{\eta^j}{j!}D_z^j b_{n-j},\qquad
D_z=\mathcal F^{-1}(i k_z)\mathcal F,
\quad k_z=\sqrt{k^2-k_x^2-k_y^2}.
\]

保留全部采样倏逝模态，kz选正实/正虚支，无滤波、相位对齐、幅度拟合。
点源/镜像Taylor系数由sqrt、倒数、exp幂级数解析计算，不用FullFK核或差分源导数。
每阶传播严格满足离散横向Fourier模式的Helmholtz方程及辐射支；FFT域是数值
周期域，通过扩域诊断镜像污染，不能直接当实际周期海面。

实际粗糙面Dirichlet残差用另外32/40阶传播指数求值，覆盖正/负eta；不是只检查
由递推自动闭合的N阶边界多项式。另有两个有限移位平面波精确高阶解测试及点源
Taylor检查，3/3通过；这些是新增有限振幅检查，不重复WP0或一阶验收。

## 参考可信程度

|自身检查|观测|
|---|---:|
|N8/10/12/14真实边界相对残差|2.26e-5 / 4.18e-7 / 5.54e-9 / 5.60e-11|
|N14→16，L48、h=.0625，增量变化|8.39e-11|
|h=.125→.0625，L48、N14，增量变化|1.071e-6|
|L32→48，h=.125；L24→48，h=.0625|1.72e-7；1.22e-7|
|边界传播求值M32→40|本配置记录差为0，不是严格零余项证明|
|最终N16真实边界残差|9.37e-12|
|半网格偏移独立真实海面残差|2.018e-7|
|解析点源Taylor系数对真实源值误差|9.37e-12|

半网格检查使用同一Fourier场在新采样位置求值、重新计算真实海面和精确点源，
不复用递推collocation点；它比原网格残差大，仍通过预设1e-4门。
传播指数求值的有限Fourier场Taylor尾项另有`sum||b_n||*exp(T)*T^(M+1)/(M+1)!`
绝对上界，T为高度界乘最大|kz|。这只控制M阶指数截断，不含浮点误差，
不是N阶FE或连续PDE误差界。详细护栏见
[reference_guards.csv](../results/validation/deepsea3d_g1_finite/reference_guards.csv)及
[offgrid.csv](../results/validation/deepsea3d_g1_finite/offgrid.csv)。

参考数值变化比下述FullFK差异小约两个数量级。阶数、网格、域及实际边界残差
共同支持**本案例的高可信数值参考**；没有证明无限精度余项界、连续算子逆范数
或任意C2海面的FE收敛半径。因此这是条件性数值物理验收，不能宣传为严格数学认证。

证据：[reference_runs.csv](../results/validation/deepsea3d_g1_finite/reference_runs.csv)、
[逐阶复系数](../results/validation/deepsea3d_g1_finite/reference_terms.csv)、
[tests.csv](../results/validation/deepsea3d_g1_finite/tests.csv)。

## 原生复声压比较与有限振幅效应

原FullFK helper不变。本轮只对新有限振幅面补做h=.05/.025及Gauss检查：
增量网格变化6.85e-9、不同求积变化1.27e-7；没有重跑旧解析或旧弱面研究。
参考反射声压为`-9.98218465e-6 + i 5.03040149e-5`；
FullFK为`-9.97526115e-6 + i 5.03048978e-5`（单位Green响应m^-1）。

|量|复相对误差|幅度相对误差|原始相位差rad|
|---|---:|---:|---:|
|反射场P_ref|1.36094e-4|9.3823e-6|-1.35771e-4|
|粗糙增量DeltaP=P_ref-P_flat|1.05578e-4|7.8116e-5|-7.10217e-5|

绝对声压差6.98e-9。高阶增量与其一阶项差异相对完整增量为**70.8%**，说明
有限振幅效应显著；不是以原一阶模型冒充独立有限粗糙度参考。
两项误差都通过[事先协议](deepsea3d_g1_finite_protocol.md)固定的单频1%条件筛查。
这里新增增量1%规则防止平面背景掩盖散射误差；并未宣称完成计划中的band-relative门。

同一海面的一处附加局部相消诊断：在6000–6075 Hz内由FK定位到6033.0075 Hz，
两种参考配置重新比较，增量变化1.094e-6，总场变化3.738e-5。
参考总场条件数53.78；FK–参考反射/增量误差1.379e-4/1.065e-4，总场误差放大至
**.364%**。保留这一放大，没有把total用于替代反射或散射验收，也不称为全带最差点。
没有出现零分母；将来精确零点须标记相对误差undefined并保留绝对误差。

证据：[comparison.csv](../results/validation/deepsea3d_g1_finite/comparison.csv)、
[cancellation.csv](../results/validation/deepsea3d_g1_finite/cancellation.csv)。

## 成本、复现和下一阶段

正式10个参考运行wall合计70.10 s；两次相消参考23.80 s。
包括已保存小预检、阶数/护栏试探、半网格检查、FK比较和相消搜索的总计算wall **137.45 s**。
CPU和分项在[cost_summary.json](../results/validation/deepsea3d_g1_finite/cost_summary.json)。
成本不含IO、测试、分析、文档及未保存开发；不是整个任务耗时。
最大参考网格1536²；16个复谱数组载荷576 MiB，另外还有forcing/临时数组，
不是进程峰值内存。未组装密集BIE、未安装新依赖。

复现：setup后运行`validate_deepsea3d_g1_preflight`，
`validate_deepsea3d_g1_finite('init')`，重复`'next'`十次（每次一个checkpoint），
然后`'offgrid'`、`'compare'`、`'cancellation'`、`'audit'`；
`runtests('test_deepsea3d_finite_reference')`须在audit前通过并写入tests.csv，详见协议。
保存`finite.mat`（纯数值/struct）、CSV、验收JSON、来源hash。五个新增MATLAB文件
静态分析通过；生产PE、二维基准和原3D FK/源/海面/phase helpers哈希保持。
MATLAB工作流归属：[Kassis等 (2026)](https://doi.org/10.48550/arXiv.2609.00065)，
不作为水声有效性证据。

**下一步：保留本案例作为有限振幅独立参考回归。全目标域G1仍OPEN，WP2=false。**
可继续一个明确的几何/谱/振幅扩展，并为新案例重新核查FE稳定性和参考收敛；
出现消项或残差失败则采用TFE或匹配半空间BIE。无需先重做二维审查，也无需为已
通过的本案例强制再做一次稠密BIE。本轮未开展任何被排除的加速、随机总体或通信研究。
