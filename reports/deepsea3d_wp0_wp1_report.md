# 深海 WP0 / WP1 最小三维单接收点阶段报告

日期：2026-10-12（Asia/Shanghai）。MATLAB R2025b Update 6，共享 MCP 会话。

**结论：WP0/G0通过；最小三维 FullFK 数值原型及一阶独立参考对照完成。
有限粗糙度三维物理 G1 仍开放，不允许进入 WP2/PM Monte Carlo或统计信道生成。**
允许继续一个有明确余项或边界收敛证据的有限粗糙度物理参考小案例。

## 1. 计划核查与修订

已有[二维 PM 物理域审查](full_kirchhoff_pm_validity_boundary_report.md)完成15个受控输入、
四个关键案例的数值护栏；它不再是待执行的普遍前置任务。通过/失败仅限单realization、
源、几何与M99核心孔径，不能继承为三维Kmax门槛。
[单Rx冻结](fk_single_receiver_baseline_freeze_report.md)只覆盖原Rx(0,3)算法域，并有相消点
失败；12个海面算法验证不能视为新增BIE物理验证。本轮没有重做这些工作或开展FastFK。

[研究计划](deep_sea_channel_research_plan.md)新增当前执行记录，区分已完成WP0、最小WP1、
开放的有限粗糙度G1；原Hm0=.5/1 m随机海面pilot移至G1后。先用毫米级弱扰动是为了
检查形状导数，不是以容易案例替代未来随机海况验收。1e-3数值变化门、1%新物理筛查
建议与旧二维1e-3物理门分别保留，没有提高容差或经验幅相拟合。

## 2. WP0：定义与解析基准

Tx=(0,0,1500) m，Rx=(1505q,0,5) m，q=0/.10（rho=0/150.5 m），c=1500 m/s。
z向下，水域z>eta；静态压力释放、均匀无损无限半空间，不含海底、平台、吸收或指向性。
单位Green点源满足`(Delta+k^2)u=-delta`：

\[
G(R)=\frac{e^{ikR}}{4\pi R},\qquad
\nabla G=G(ik-1/R)\frac{\mathbf r-\mathbf r_t}{R}.
\]

Green响应单位为m^-1；源方程系数为1，不是校准Pa声压或指定声功率。native采用
exp(-i omega t)。平面高度h0时镜像源为(0,0,2h0-1500)，直达/反射/总场分别为
G(Ld)、-G(Lr)、二者之和。DSP为`conj(P_native)*exp(i omega Ld/c)`，三分量共享
直达时刻，不分别对齐；反射DSP相位为`-exp(-i omega (Lr-Ld)/c)/(4pi Lr)`。

24组控制覆盖q、4/6/8 kHz及6037.5 Hz、h0=-.2/0/.2 m，另有源导数、法向、完整
积分闭合和旋转测试。最大代数误差7.276e-12、相位表示误差5.258e-12，低于冻结
1e-10门。8项测试全部通过。垂直h0=0的超额时延为6.6666667 ms；不是CIR仿真。
平面任意角度旋转检查为.713 rad，粗糙面实际90度旋转同步旋转Rx、中心和波矢。

证据：[wp0.csv](../results/validation/deepsea3d_wp0_wp1/wp0.csv)、
[tests.csv](../results/validation/deepsea3d_wp0_wp1/tests.csv)。

## 3. WP1：完整积分、数值证据与失败控制

独立helper `solve_deepsea3d_fullfk`计算

\[
P_r=-\int_{\Gamma}\left(G\,\partial_n u_{inc}+u_{inc}\,\partial_nG\right)dS,
\quad \mathbf n=(-\eta_x,-\eta_y,1)/J,\quad J=\sqrt{1+\eta_x^2+\eta_y^2}.
\]

实现按解析`n*dS=(-eta_x,-eta_y,1)dxdy`收缩；这是面积/法向的严格抵消，不是平法向
或dxdy近似。G对源积分点的导数和入射场导数都使用解析式，两项分别保留。
真实二维海面含非共线波矢；未用一维截面或2.5D替代。

海面为径向C2紧支撑扰动，plateau8/support12 m，中心位于平面镜面点(q=.1时x=150 m)。
weak：A=.001 m、K=(.4,0) rad/m、phase=.2；dual：A=(.002,.001) m、
K=(.4,0),(.2,.6)、phase=(.2,.7)。这些参数不是PM或海洋波高标定。

对于这样的紧支撑表面，用解析无限平面背景加粗糙/平面完整积分密度之差。
这是同一个Kirchhoff积分的解析背景减法，不是校准或新物理模型；平面情形的精确
匹配由背景构造保证，不能单独证明粗糙积分正确。另以Gauss求积、一阶参考及Fourier
传播检查粗糙响应。原始平面截断结果保留，不能以构造闭合隐藏截断失败。

|检查|实际设置|结果|
|---|---|---|
|6 kHz网格|q=0/.1，flat/weak/dual，h=.1/.05/.025 m，L16 m|最后两网格最大反射变化9.569e-10、扰动变化8.921e-9|
|共同海面孔径|h=.025，L8/12/16/20 m|L>=12覆盖紧支撑，变化约舍入级；L8为有意截断诊断|
|不同求积|Gauss4与梯形，L16|最大反射差约1.40e-8|
|坐标旋转|实际90度同步旋转|小于1e-10|
|4/8 kHz扩展|先通过6 kHz门；h=.05/.025，两个q、三个面|全部反射/非零扰动网格变化<1e-3|
|相消敏感性|四个局部低压位置，h=.05/.025|条件数205–271，总场相对变化最大1.074e-7|

原始有限平面积分并不通过：q=0、L8/12/20的误差为8.75%/6.21%/4.19%；q=.1时
L8偶然较小(.262%)，L12却为6.07%，不能据此选一个“好看”孔径。解析背景减法
解决当前紧支撑模型的平面尾部处理，不提供无限非紧支撑随机海面验收。

support18/24控制保留相同未嵌入正弦母面，但改变嵌入几何，故仅是边界敏感性。
本案例扰动复声压变化很小；两种嵌入不能叫同一物理海面的孔径加密。
相消搜索只在平面预测干涉零点附近±5 Hz局部最小化，没有幅相拟合，也不证明全带
最坏相消点。若将来出现零分母，应保留绝对变化并标记相对误差undefined。

证据：[runs.csv](../results/validation/deepsea3d_wp0_wp1/runs.csv)、
[guards.csv](../results/validation/deepsea3d_wp0_wp1/guards.csv)、
[cancellation.csv](../results/validation/deepsea3d_wp0_wp1/cancellation.csv)。

## 4. 独立弱粗糙度参考与物理缺口

直接展开Dirichlet边界`u(x,y,eta)=0`，取平面总场u0，则一阶边界值
`b=-eta*partial_z u0=-2 eta*partial_z u_inc`。水侧传播由半空间Rayleigh核给出：

\[
\delta u(\mathbf r)=\int b(x,y)\left[-2G'(R)\frac{z_r-h_0}{R}\right]dxdy.
\]

这是从PDE边界条件独立得到的线性参考，不调用FullFK核、法向或粗糙积分。
Gauss4 h=.05/.025及L16/20护栏的扰动相对变化最大1.444e-7/4.08e-15。
另实现Fourier半空间传播`exp(i kz(z_r-h0))`，保留正虚部kz衰减支；q0 weak的
FFT域L32/48/64、h=.0625对Rayleigh差异5.51e-8/1.07e-8/1.36e-8。
这排除了当前参考的明显传播符号或有限周期域错误，但不认证任意周期域。

两个q、两个面、scale1/.5/.25共12个对照。weak的正负振幅中心差对参考偏差从约
3.5e-4下降至2.6e-5–4.4e-5；dual从约3e-3下降至1.0e-4–1.13e-4。
有限正振幅扰动项偏差随振幅减半约减半：weak约2.45%→1.23%→.61%，
dual约6.8%→3.4%→1.7%。这是与一阶余项趋势相容的有限集证据，没有拟合余项界。

最大反射总量差约.931%，虽然低于计划建议1%，仍不能认证有限粗糙度物理：
平面背景稀释扰动差，而且一阶参考自身遗漏O(eta^2)及以上项。不能把它当收敛BIE。
本轮只建立独立弱极限检查；G1有限粗糙度PASS=false。

证据：[shape_reference.csv](../results/validation/deepsea3d_wp0_wp1/shape_reference.csv)、
[spectral_reference_check.csv](../results/validation/deepsea3d_wp0_wp1/spectral_reference_check.csv)。

## 5. 匹配三维 BIE/BEM 可行性

数学路线存在：sound-soft三维粗糙半空间可用适当半空间Green的组合层BIE，但无限
表面的慢衰减和算子性质需要专门处理；不能将平面圆盘/方片当封闭障碍物直接套BEM。
参考：[Chandler-Wilde、Heinemeyer、Potthast (2006)](https://doi.org/10.1137/050635262)。
紧支撑三维扰动与点源问题也有专门半空间理论，见
[Zhao、Hu、Yan](https://arxiv.org/abs/1812.09009)。这些理论不提供本fixture的数值PASS。

WGF提供一种无限界面截断路线，但层状介质论文不能直接当本Dirichlet单界面的已配置
实现：[Bruno、Pérez-Arancibia](https://arxiv.org/abs/1703.01034)。
[Bempp官方文档](https://bempp.com/handbook/api/boundary_operators.html)有Helmholtz层算子、
奇异求积和FMM接口；它们不自动给出匹配的无限压力释放粗糙海面。未安装或调用新依赖。

成本下限示例：6 kHz波长.25 m，L16的32×32 m面若按lambda/8类均匀节点采样，
约1,048,576未知量，单个complex-double密矩阵即约16 TiB；即便24×24 m也约5.06 TiB。
这是网格类比估算，不是对高阶BEM最优未知量或实际内存的测量。高阶面元、FMM、
层次矩阵、WGF可改变成本，但须先匹配边界公式、奇异/近奇异求积、残差及接收场
收敛，不能从旧二维时长外推三维成本。本轮未尝试这一稠密规模。

**判断：可负担的一阶独立参考已经实现；匹配有限粗糙度BIE在数学上可行，当前没有
已验证且可直接运行的匹配实现，朴素密集方案不合适。其可用实现及实测成本仍开放。**

## 6. 成本、复现与交付

|记录范围|wall s|MATLAB CPU s|
|---|---:|---:|
|106个保存FullFK运行|60.864|349.375|
|12组Rayleigh参考，每组3个护栏运行|11.924|72.359|
|3个Fourier交叉检查|2.029|9.078|
|4次局部相消搜索（42次函数调用）|5.634|37.125|
|以上合计|80.451|467.938|

多线程CPU可大于wall。以上是保留结果的计算成本，不含IO、测试、代码分析、文档、
失败测试/开发或替换前旋转诊断，因此不是整个任务耗时。最大流式FullFK网格约502万
节点，估计块数组载荷32.83 MiB；最大FFT复矩阵64 MiB，都不是进程峰值内存。

入口[validate_deepsea3d_wp0_wp1.m](../scripts/validation/validate_deepsea3d_wp0_wp1.m)，
按[预设协议](deepsea3d_wp0_wp1_protocol.md)执行七个阶段；可从任意当前目录在setup后
调用。各阶段可重跑，pilot重置后续状态；中间保存用于检查，未承诺阶段内部自动续跑。
MAT仅保存数值/字符串struct，无函数句柄。每个相同case确定性复现，计时受会话影响。

9个新MATLAB文件Code Analyzer无诊断。冻结旧代码/输入9项当前SHA256均相同。
保存CSV、state.mat、artifact_assertions.json、cost_summary.json及provenance.json在
`results/validation/deepsea3d_wp0_wp1/`。验收JSON明确记录WP2_allowed=false。
更新了研究计划、PROJECT_CONTEXT、技术指南与scripts运行索引；生产PE、公共API、
二维基准与通信算法未改，用户排除通信/BER，本轮不重复无关生产通信回归。

后续最小任务是一个匹配有限幅度3D参考及其精化，或经验证的高阶展开/余项界。
在此之前，不把数值收敛、一阶比较或小总反射误差升级为真实海况物理有效性。

工作流使用MATLAB技能组织确定性代码、测试和来源记录；归属
Kassis T., Agarwal V., He Y., Patel D., Brueckner A. M. (2026),
[Scientific Agent Skills](https://doi.org/10.48550/arXiv.2609.00065)（2026-10-12核对当前记录）。
此归属不是水声物理有效性的证据。
