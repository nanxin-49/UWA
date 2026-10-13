# 垂直水声 PE/WAPE 信道技术指南

### 三维TFE去混叠与失效诊断（2026-10-13；最新）

验证分支对L1二因子、L2三因子乘积采用2n填充、一次投影及显式Nyquist空间，DtN乘积同样去混叠。
同一冻结PM support12不变；独立3n/4n真实连续几何评价保留完整与投影PDE残差。
dual/quarter控制通过，但原幅度N32去混叠体PDE .4753/DtN .3265，均未过原1e-4门。
末阶系数减小2.363倍不等于物理收敛；99.60%残差能量在保留空间，几何算子差仅1.15e-5。
停在N32/既有网格，不做support24/FullFK–PM；PM INCONCLUSIVE、G1 OPEN/WP2=false。
建议下一轮直接求解t=1展平Helmholtz、矩阵自由且受预算约束，而非继续堆Taylor阶数。
生产PE、通信公式与冻结求解器均不变。[报告](reports/deepsea3d_tfe_dealiasing_report.md)。

### 三维TFE独立参考最小原型（2026-10-13；最新）

验证专用展平z=s+t(1-s/b)eta、a=1-t eta/b，粗糙修正的算子满足a²(Delta+k²)=L0+tL1+t²L2。
精确物理DtN为v_z=T v，T=i sqrt(k²-kx²-ky²)、Im sqrt≥0；展平上边界Vs−TV=−t eta/b TV。
点源在薄层外，Rx由固定平面传播，背景仍是精确无限平面镜像；横向周期修正须独立扩域。
平面与制造解通过，认证dual完整数值检查通过，最细增量对冻结Dirichlet参考差6.308e-7。
同一冻结PM seed261013/support12在N40仍无收敛；h.125体PDE约.120/DtN约.173。
quarter幅度仅诊断局部门通过，不是Hm0=.5目标认证；未做support24或FullFK–PM。
PM参考INCONCLUSIVE，全G1 OPEN/WP2=false；没有改动生产PE、通信公式、冻结解算器或历史数据。
[报告](reports/deepsea3d_tfe_prototype_report.md)、[协议](reports/deepsea3d_tfe_prototype_protocol.md)。

### 三维PM固定realization最小可行性（2026-10-13；最新）

各向同性有限带PM，P(kx,ky)=alpha/(4pi*kappa^4)*exp[-1.25*(kp/kappa)^2]，
alpha=.0081，父域128m、kmin2pi/128、kmax1rad/m，解离散谱质量Hm0=.5得到Tp3.638671s。
mt19937ar seed261013、固定复高斯Fourier系数、共轭实面；实际Hm0=.504033，不作样本RMS拟合。
6kHz、Tx(0,150.5,1500)/Rx(0,0,5)、q=.10。C2支撑12/24m仅物理嵌入诊断，不能冒充无限PM。
普通FE最高N20未过阶数/源项/真实边界门，直接有限模态指数评价也给出失败见证；
接收压力在网格/域下变化小不足以认证参考。谱截止、远区和支撑影响仍缺可信压力证据。
无FullFK物理比较或新G1，单案例INCONCLUSIVE、全G1 OPEN/WP2=false；建议先评估匹配TFE/3D BIE。
生产PE、通信接口、既有G1及相位公式不变，无通信运行。
[报告](reports/deepsea3d_pm_feasibility_report.md)、[协议](reports/deepsea3d_pm_feasibility_protocol.md)。

### 三维G1固定y方位4/6/8kHz检查（2026-10-13；最新）

固定q=.10/alpha2.5/dual/C2/Tx(0,150.5,1500)/Rx(0,0,5)，只改频率4/8kHz，6kHz冻结复用。
4kHz选N16，8kHz预检补N20，阶数变化8.931e-5接近原1e-4门；网格/域/真实边界及消项通过。
各新频点M32/M40真实边界对照与有限Fourier指数尾界通过，不默认M40或旧网格充分。
8kHz Gauss4 h.05曾失败（1.200e-3）；保留失败并精化一次h.025，变化1.683e-9通过原1e-3门。
4/6/8各点条件性G1 PASS，反射复误差4.034e-4/2.843e-4/2.146e-4，增量2.235e-4/1.375e-4/1.217e-4。
含阶数敏感性仅支持4→6复误差可分辨下降，6→8范围重叠；相位误差下降、增量幅误差6→8上升。
三离散点不认证整个4–8kHz带宽。原FE/FullFK、原生相位、生产公共接口和通信公式不变。
全G1 OPEN/WP2=false。[报告](reports/deepsea3d_g1_frequency_report.md)、[协议](reports/deepsea3d_g1_frequency_protocol.md)。

### 三维G1 q=.10水平方位角对照（2026-10-13；最新）

固定alpha2.5/6kHz/dual/C2海面，只改Tx方位x→y：(150.5,0,1500)→(0,150.5,1500)，
Rx=(0,0,5)，不旋转或移动海面。两方位平面解析背景一致；y独立FE按需N20，
h=.0625/L32参考门与真实半网格残差3.211e-7通过，原FullFK网格/Gauss/孔径门通过。
**y方位条件性G1 PASS**，反射/增量复误差2.843e-4/1.375e-4，均过原1%门，无拟合。
参考预测y−x反射/增量复变化5.111%/2.467%；与约2e-10m^-1已观察参考变体位移分开。
模型误差相比x下降1.920%/1.762%，含截断阶数敏感性的经验范围不重叠。
低阶部分和不是额外独立边界认证；经验范围非置信区间或严格PDE界，不外推普遍方向规律。
原FE、完整两项FullFK、生产公共接口及通信公式均未改；全G1 OPEN，WP2=false。
[报告](reports/deepsea3d_g1_azimuth_report.md)、[协议](reports/deepsea3d_g1_azimuth_protocol.md)。

### 三维G1 q=.10固定海面几何检查（2026-10-13；最新）

仅将冻结alpha2.5/6kHz/dual/C2案例的Tx改为(150.5,0,1500)，Rx仍(0,0,5)，
q=.10；海面中心保持(0,0)，不向镜面点(.5,0)移动。WP0平面镜像及相位约定通过。
重新低成本检查N12/N16后按需N20；h=.125→.0625及L24→32变化1.400e-6/4.161e-7，
最终真实半网格残差3.013e-7，参考过原1e-4门。FullFK网格、Gauss及独立孔径过1e-3门。
反射/增量原生复误差2.899e-4/1.399e-4，各过原1%门，**q=.10单案例条件性PASS**。
相比q0误差下降2.03%/2.94%，已检验数值变体范围不重叠；不是连续q或PM总体定理。
粗糙total相消指标1.278，total误差略增；无相消频率搜索，不以total替代反射/增量验收。
原FE递推、完整两项FullFK及原生/DSP公式、生产公共接口和通信链均未改。
全G1 OPEN，WP2=false。[报告](reports/deepsea3d_g1_q10_report.md)、[协议](reports/deepsea3d_g1_q10_protocol.md)。

### alpha2.5计算域与条件性物理G1（2026-10-13；最新，接续下述阶段记录）

冻结Tx1500/Rx5、q0、6 kHz、dual振幅(.05,.025) m与C2 plateau8/support12。
L仅扩展Fourier计算域，不缩放海面；L24→32/h=.0625/N20增量变化3.242e-7，
N16→20变化2.363e-6，真实半网格残差2.977e-7，过原1e-4门且高阶尾部递减。
无需追加L48。参考通过后，FullFK独立A12/16/24紧支撑密度差孔径检查通过，
网格/Gauss变化1.021e-8/1.905e-7；反射/增量原生复误差2.959e-4/1.442e-4，
均过原1%物理门，无幅相拟合。**该alpha2.5离散案例条件性PASS**，全G1仍OPEN/WP2=false。
沿用原FE递推与完整两项FullFK/解析无限平面背景，无新物理公式、公共接口或通信方法。
经验域稳定不等于严格PDE余项界；不覆盖连续振幅、PM总体或宽带。
[报告](reports/deepsea3d_alpha25_domain_report.md)、[协议](reports/deepsea3d_alpha25_domain_protocol.md)。

### alpha2.5细网格参考检查（2026-10-13；最新）

L24/h=.0625 N16/N20检查PASS：阶数差2.36e-6，粗→细增量变化1.77e-6，
N20半网格残差3.06e-7；系数尾部与h归一化粗网格范数一致，内部消项未灾难性放大。
此处PASS仅是普通FE当前配置参考稳定性，不是alpha2.5物理G1或域收敛认证。
未调用FullFK或扩域；alpha2.5物理G1 INCONCLUSIVE，全域OPEN/WP2=false。
[报告](reports/deepsea3d_alpha25_fine_report.md)、[协议](reports/deepsea3d_alpha25_fine_protocol.md)。

### alpha2.5低成本参考预检（2026-10-13；最新）

同一单链/dual/C2，alpha2.5的粗配置N12不足；追加N20后阶数变化2.36e-6、
半网格真实边界残差2.87e-6。普通FE在该配置未见失稳，但高阶消项已明显。
只建议受控细网格正式参考检查；没有新FK物理对照、域扩展或alpha2.5 G1 PASS。
alpha2.5物理状态INCONCLUSIVE，全目标域OPEN/WP2=false；冻结公式和求解器未改。
[报告](reports/deepsea3d_alpha25_preflight_report.md)。

### 三维G1统一振幅扩展（2026-10-13；最新）

冻结Tx1500/Rx5、q0、6 kHz、dual波矢/phase/C2支撑，仅缩放(.02,.01) m。
alpha1复用；.5/1.5/2新增点各自参考阶数/网格/域/真实边界检查通过，并均条件性G1 PASS。
最大反射/增量复误差2.57e-4/1.27e-4，最大半网格残差2.44e-7；未见可信物理失败。
原1%物理和参考1e-4门不变；参考不稳定须INCONCLUSIVE，不判FullFK失败。
这只覆盖四个离散点，全目标域G1仍OPEN/WP2=false；冻结模型和相位约定均未改。
[报告](reports/deepsea3d_g1_amplitude_report.md)、[协议](reports/deepsea3d_g1_amplitude_protocol.md)。

### 三维有限振幅单案例物理接续（2026-10-12；最新）

独立高阶Dirichlet FE对Tx1500/Rx5、q0、6 kHz、振幅(.02,.01) m双方向C2海面
提供有限振幅参考。自身阶数/网格/域/真实边界残差护栏通过，原FullFK反射及粗糙增量
误差分别1.36e-4/1.06e-4；同一海面局部相消总场误差.364%。
本单案例条件性G1 PASS，新增单频反射/增量1%筛查不等于原二维物理1e-3或宽带G1。
没有严格余项界、全目标域/PM总体认证，WP2仍不批准。下节WP1的OPEN为此前阶段状态。
原生场/源/海面/相位约定和生产代码未改变。
[报告](reports/deepsea3d_g1_finite_report.md)、[协议](reports/deepsea3d_g1_finite_protocol.md)。

### 三维深海研究原型边界（2026-10-12）

`validate_deepsea3d_wp0_wp1`独立验证Tx1500/Rx5、q=0/.10、均匀压力释放半空间。
点源`G=exp(ikR)/(4piR)`，native为exp(-i omega t)，DSP=`conj(P)*exp(i omega Ld/c)`，
直达/反射共享时间原点。3D FullFK完整保留法向、面积和两项积分；紧支撑C2海面
以精确无限平面背景加粗糙/平面密度差处理。它不是生产PE边界、PM统计或公共H_f。
WP0解析门通过，6 kHz最小原型数值门及4/8 kHz网格检查通过。
独立Dirichlet形状导数和Fourier检查只支持一阶弱扰动；有限粗糙度G1仍开放。
1e-3是本阶段数值变化门，计划1%是物理筛查建议；不等同于旧二维物理1e-3。
不建立收敛CIR、3D总体可靠性或BER。[阶段报告](reports/deepsea3d_wp0_wp1_report.md)
记录公式、参数、成本、失败的原始平面截断及下一阶段门。

本文件是当前实现的中文权威技术参考，说明物理约定、代码语义、公共字段、验证证据和已知限制。开发时间线保存在 `PROJECT_CONTEXT.md`，具体数值证据保存在 `reports/`；`docs/history/` 中的旧规格不是当前接口依据。

状态日期：2026-09-03。

补充验证（2026-10-12）：固定Tx100–单Rx(0,3)的两种独立Kirchhoff原型见
[中文方案对比](reports/fk_single_receiver_methods_report.md)与
[公式/协议](reports/fk_single_receiver_methods_protocol.md)。FastIn只压缩入射，
出射保留标量FullFK；复合核使用`P=G*u_inc`的物理偏导表，并按真实海面高度
插值，法向及弧长每次重新计算。两者仅在验证目录，不改变本指南公共信道接口、
生产PE或通信方法，不提供新的统计信道、CIR、BER或3-D物理证据。
同日验证补强及[冻结审核](reports/fk_single_receiver_baseline_freeze_report.md)
仅冻结Rx(0,3)原fixture。FastIn可复用入射plan到新Rx，但已有相消点诊断
不满足rank3的1e-3复误差门槛；复合核新Rx必须重建。零点相位及raw相对误差
无定义时明确标记，不用分母保护值宣称准确；这些研究状态不改变公共接口。

## 1. 坐标、时间与载波相位参考

- 海面为 `z=0`，深度向下为正。
- 发射机满足 `z_tx>z_rx>=0`，向上传播对应 z 减小。
- 物理相量约定为 `exp(-1i*omega*t)`；通信合成使用 MATLAB `ifft`。
- PE 推进的是相对参考载波的复包络，不是逐周期求解的瞬时声压。

直达和海面反射的名义参考时延为

\[
\tau_{\mathrm{dir},0}=\frac{z_{\mathrm{tx}}-z_{\mathrm{rx}}}{c_0},\qquad
\tau_{\mathrm{ref},0}=\frac{z_{\mathrm{tx}}+z_{\mathrm{rx}}-2z_s}{c_0},
\]

\[
\Delta\tau_0=\tau_{\mathrm{ref},0}-\tau_{\mathrm{dir},0}.
\]

标准几何 `z_tx=100 m`、`z_rx=3 m`、`z_s=0`、`c0=1500 m/s` 给出 `Delta tau0=4 ms`。这只是两条 PE 约化参考之间被消去的确定性载波延迟，不等于粗糙海面反射峰值的全部实际时延，也不允许分别把直达和反射峰“对齐”。

公共默认 `paramsV.channel_phase_reference='direct_dsp'`：

\[
H_{\mathrm{direct}}=H_{\mathrm{direct}}^{\mathrm{red}},\qquad
H_{\mathrm{reflect}}=e^{-i2\pi f\Delta\tau_0}H_{\mathrm{reflect}}^{\mathrm{red}}.
\]

`legacy_reduced` 仅用于复现旧结果；它把分别约化的两条路径直接相加，不是统一物理相位的总信道。绝对物理相量另外保存为

\[
H_{\mathrm{dir}}^{\mathrm{phys}}=e^{+i2\pi f\tau_{\mathrm{dir},0}}H_{\mathrm{dir}}^{\mathrm{red}},\quad
H_{\mathrm{ref}}^{\mathrm{phys}}=e^{+i2\pi f\tau_{\mathrm{ref},0}}H_{\mathrm{ref}}^{\mathrm{red}}.
\]

载波相位只在各传播分量完成后转换，不进入 PE marching、海面边界、`deltaG`、伴随核或 PM 空间协方差。载波相位发布候选已通过，详见 `reports/pe_phase_reference_release_candidate_report.md`。

### 验证专用 Full/Fast FK 原生物理场的输入桥接（2026-10-08）

独立 FullFK/BIE/FastFK reference 保存的是 `exp(-iωt)` 原生物理相量
`P_native`，冻结表面坐标为 `z=eta(x)`。不能只除去其载波，就自动把得到的
native reduced envelope 当作上述公共 PE reduced 字段。对正序声学频率轴，
MATLAB IFFT 的原生物理对应传递函数在同一个直达时刻参考下是

\[
H_{\mathrm{native,DSP}}(f)
=\operatorname{conj}(P_{\mathrm{native}}(f))
 e^{+i2\pi f\tau_{\mathrm{dir},0}}.
\]

反射路径等价于先构造
`R_native=P_ref*exp(-i2*pi*f*tau_ref0)`，再把 `conj(R_native)` 送入现有
helper 的 `direct_dsp` target。直达包络对应
`conj(P_dir*exp(-i2*pi*f*tau_dir0))`。这是一次数学定义的全包络 Fourier
符号转换，所有频点和分量共享同一个 time origin，不是逐频 phase alignment
或 complex scalar fitting。原生物理 `P_native` 独立保留；不得用 helper 的
`absolute_physical` target 从 DSP-adapted envelope 反推它。本轮不改公共
PE 字段语义或生产 API。

固定非零横向波数的解析控制表明，只旋转载波的未转换导入会给出错误的
原生物理群时延解释，而上述桥接与独立构造的 IFFT transfer 复误差约
`2.07e-13`。因此新报告主结果使用明确标记的 native physical direct-DSP
表示，未转换导入只作诊断。报告见
`reports/spectral_fast_kirchhoff_multifrequency_consistency_report.md`。

## 2. PE/WAPE complex envelope

传播核心位于 `src/propagation/vertical_wape_propagator.m`。均匀介质的一步 cached transverse operator 保持以下顺序：

```matlab
psi_k = fr .* fft2(screen .* ifft2(fr .* psi_k));
```

`fr` 表示半步横向谱传播，`screen` 包含该步介质相位和 sponge 衰减。输出空间图中的颜色通常表示复包络相对幅度，例如

\[
20\log_{10}(|\psi|/|\psi|_{\mathrm{ref}}).
\]

它描述载波振幅和相位随传播的慢变结构；若展示 `real(psi*exp(-i2*pi*f*t))`，那只是单频载波重构动画，不是时域 PE 求解。

公共实现逐频计算直达路径与“Tx→Surface→Rx”反射路径，并始终保持

\[
H_f=H_{\mathrm{direct},f}+H_{\mathrm{reflect},f}.
\]

## 3. Gaussian 源、横向窗口与 sponge

生产 Gaussian 初场由 `src/propagation/gaussian_source_initial_field_vertical.m` 构造。横向 FFT 网格是周期数值域，边缘能量会绕回中心；sponge 用于衰减边缘，但也可能改变目标接收响应，不能把 sponge 当成无限孔径的替代品。

当前公共默认仍为 `xw=yw=50 m`、`sponge_ratio=0.12`、`alpha_max_np_per_m=0.15`，本次目录整理没有改变默认值。近期独立验证给出的严格建议是：

- Gaussian 直达路径：`160 m / no-sponge`；
- 完整反射链：实际 `192.1875 m / no-sponge`；
- 4 kHz 随机海面集合中，192.1875 m 对 3 个海况、每个 5 个 seed 为 `15/15` 严格通过；
- 160.15625 m 的接收频响近似收敛，但外围能量仍超门限，不能标记为严格通过。

这些建议尚未自动成为生产默认。随机海面 3–5 kHz 的完整候选—参考成对集合没有跑完，因此不能声称已经获得 ensemble 最差群时延。证据见 `reports/pe_reflected_chain_window_validation_report.md` 和 `reports/pe_random_surface_window_robustness_report.md`。

### 展开坐标 PE--Bellhop 验证

平面、均匀介质的独立 Bellhop 对比入口为
`scripts/validation/validate_pe_bellhop_unfolded_flat_gaussian_vertical.m`。
它把原竖直方向展开成 Bellhop 的水平距离：`z_tx=100 m`、`z_rx=3 m`
分别对应 97 m 直达和 103 m 镜像反射路径，后者在展开接收端施加
压力释放系数 `-1`。Bellhop 的 `.sbp` 由现有 unit-peak Gaussian 的角谱
生成；主验收比较归一化横向场、`H_reflect/H_direct`、相位和群时延，
不把 Bellhop 点源绝对幅度逐点拟合到 PE。该结论只覆盖均匀声速、平面
海面和当前 Gaussian 源；更换换能器后必须重新生成源指向性并重新确认
横向窗口。结果目录为
`results/validation/pe_bellhop_unfolded_flat_gaussian/`。

Bellhop 2020 validation-only internal-wall 的当前权威状态统一见
`reports/bellhop_internal_wall_implementation_report.md`。实现链为：参数化墙
`Gamma(s)=[r(s),z(s)]` 的有限段求交与 `ReduceStep2D` 截断，命中点单位切向、
由切向构造的 TOP 法向、有符号转角/弧长曲率，原生 2D `Reflect2D` 一次，随后
仅用 `r'=2R0-r, z'=-z` proper half-turn 重映射反射分支，最后进入未修改的
`InfluenceGeoHatCart`。curved adapter 不再依赖 `dz/dr`、`Dss`、
`Delta z/Delta r`、`1/Delta r`、range 单调性或 ATI 无限端点延拓。

flat、tilted、smooth sinusoidal 与半径 100 m 的 vertical-tangent circle 均已
回归通过。flat 的 3×3 step/beam scan 恢复 `H_wall=-P_BH(103 m)`；sinusoidal
正确配对 receiver column 后 phase difference 不超过 `6.71e-5 rad`；circle 在
严格 `dr/ds=0`、非 grazing 的 `u·n=1` 下无 NaN/Inf，signed curvature 随
65/129/257 samples 收敛到 `-0.01 1/m`。曾出现的约 `2.095 rad` 偏差是读取了
错误 SHD range column，不是 reflection/frame 缺陷。native backward-range
Cartesian amplitude 的 `-0.2606` 至 `-0.286 dB` 偏差仍只作 diagnostic，未拟合。

固定 PM realization 的 internal-wall density convergence 已完成，报告为
阶段性报告已归档于 `cash/bellhop_internal_wall_superseded_20260902/reports/bellhop_internal_pm_fixed_realization_convergence_report.md`，结果位于
`results/validation/bellhop_internal_pm_fixed_realization/`。seed `260001`、
`U=6 m/s`、span `160 m`、`N_master=4097` 的同一 band-limited realization
生成 `N=513/1025/2049/4097` 四档 profile。preflight 通过（最小
`|u·n|=0.8022`、`lambda max|kappa|=0.01672`），2001/2001 rays 全部成功，
交点、frame、曲率、反射方向、RN、路径和 phase 随密度收敛。该阶段仍为
Bellhop-only；旧 strict-90-degree native ATI 对照不是等价物理问题并已退役。
随后已完成 Bellhop-only 的 fixed-PM local native↔internal Reflect2D
coordinate-covariance：同一 seed-260001 master realization 在 `phi=89` 和
`89.5` 度、`N=2049/4097`、三条精确配对 ray 上通过交点、frame、共同曲率
参考、反射方向、RN/p/q、pressure-release phase 和 internal half-turn 检查。
详见权威汇总 `reports/bellhop_internal_wall_implementation_report.md`（完整阶段报告已归档于 `cash/bellhop_internal_wall_superseded_20260902/reports/`）；硬门截止
Reflect2D 出口，receiver coherent field 仍是 diagnostic-only。随后完成了
受控 PE rough-PM 对比链及两个低成本 Stage 3 ensemble smoke；结果与限制统一见
`reports/pe_bellhop_pm_complete_execution_summary.md`。

PE rough-PM ↔ Bellhop rough-PM 方案曾在
`reports/pe_bellhop_pm_comparison_design_audit_report.md` 中冻结，现已按
Stage 0A--0E、Stage 1A--1C、Stage 2 和 Stage 3 执行完成。主比较量仍是
reflected-only `G=H_ref_rough/H_ref_flat`；一横向维 PE bridge、生产二维
dimensionality sensitivity、4/6/8 kHz 扩展以及 8-seed phase/amplitude
ensemble 的状态和限制见 `reports/pe_bellhop_pm_complete_execution_summary.md`。

Stage 0A canonical mapper 已完成并通过：
`scripts/validation/validate_pe_bellhop_pm_canonical_mapper.m` 直接读取既有
`fixed_pm_fourier_coefficients.csv`，以同一周期 Fourier series 为 PE 与 Bellhop
采样；输出见 `results/validation/pe_bellhop_pm_canonical_mapper/`，报告见
`reports/pe_bellhop_pm_canonical_mapper_report.md`。本阶段只冻结 profile provenance
和 `Gamma_B(s)=[R0-eta(s),s]` 映射，不运行 propagation，也不能替代后续的一横向维
PE bridge、flat/source、constant-height sign 和 numerical-budget audits。

Stage 0B 一横向维 PE bridge 已通过：
`scripts/validation/validate_pe_1d_validation_bridge.m` 用独立一维 split-step、
一步 exact angular spectrum 和完整生产 PE 的 `k_y=0` 投影交叉核对 flat 与弱
phase-screen case。最大复误差分别约为 `1.14e-12` 和 `9.48e-13`，flat
pressure-release sign 误差约 `2.64e-13`；结果在
`results/validation/pe_1d_validation_bridge/`，报告为
`reports/pe_1d_validation_bridge_report.md`。这只证明同维 PE bridge，不证明
PE--Bellhop 粗糙面场一致；Stage 0C 及后续阶段均已完成，详见完整汇总报告。

Stage 0C flat source/normalization audit 已通过：
`scripts/validation/validate_pe_bellhop_pm_stage0_flat_source.m` 使用现有 Gaussian
`.sbp`、97/103 m 显式 SHD range selector、5001/10001 beams 和 0.05 m step，
核对一横向维 PE 与 Bellhop 2020 internal-flat 的相对传播定义。轴上 Q 的 phase
残差为 `5.15e-4 rad`，normalized offset magnitude residual 为 `5.97e-3`
（direct）和 `5.55e-3`（reflected）；约 `0.261 dB` 的 backward-range amplitude
偏移保留为 diagnostic，不做拟合或归一化。带 `.sbp` 的 SHD 深度记录读取也已按
detected range record 做最小兼容修复。结果见
`results/validation/pe_bellhop_pm_stage0_flat_source/` 和
`reports/pe_bellhop_pm_stage0_flat_source_report.md`；Stage 0D 及后续阶段均已
完成，详见完整汇总报告。

Stage 0D 常量高度符号审计已通过：
`scripts/validation/validate_pe_bellhop_pm_constant_height_sign.m` 使用
`eta0=+0.05, 0, -0.05 m`，并将物理 image span 固定为 `L=103-2*eta0 m`。
PE 的 `exp(+i2k eta0)` 相位误差低于 `2.0e-15 rad`，Bellhop 低于
`2.56e-5 rad`，两者相对相位误差低于 `2.56e-5 rad`；路径时间误差为
`4.16e-17 s`，压力释放相位仅施加一次。结果见
`results/validation/pe_bellhop_pm_constant_height_sign/` 和
`reports/pe_bellhop_pm_constant_height_sign_audit.md`。validation-only
overlay 仅放宽了移墙输入范围检查，未改动 Reflect2D、p/q 或 influence。

Stage 0E 数值误差预算已冻结并通过：
`scripts/validation/validate_pe_bellhop_pm_numerical_budget.m` 对同一
seed-260001 band-limited realization 做 PE 窗口/网格/步长及 Bellhop
profile/beam/step 扫描；PE 与 Bellhop 扫描均落在 0.1 dB、0.02 rad 预算内，
Bellhop 墙残差低于 `8e-15 m` 且 q-state residual 为零。结果见
`reports/pe_bellhop_pm_numerical_error_budget.md`；已冻结的平面源基线仍保留
0.30 dB backward-range influence diagnostic allowance。

Stage 1A Tier-1 固定 PM 比较已完成，状态为 PASS_WITH_LIMITS：
`scripts/validation/validate_pe_bellhop_pm_stage1_tier1.m` 的结构几何、
压力释放相位和 Bellhop p/q 状态检查全部通过。PE Kirchhoff phase-screen 与
Bellhop local-specular Gaussian-beam 的 reflected-only ratio 差异（0.3104 dB、
−2.2482 rad）只作模型差异诊断，不作结构 gate。详情见
`reports/pe_bellhop_pm_stage1_tier1_report.md`。

后续 Bellhop source-geometry 单变量审计明确区分了 point-source `R` 与
line-source `X`。在不改变 `.sbp`、几何、step、beams、receiver、Reflect2D、p/q
或 influence 公式的条件下，10001-beam 弱正弦墙 native/internal 对比从
`R` 的 `-0.260650 dB` 降至 `X` 的 `4.24e-6 dB`；前者与
`10 log10(97/103)=-0.260655 dB` 一致。5001-beam 正弦场有两种 source 共有的
`3.156 dB` 未收敛离散变化，不能作为端点。对 4 kHz、seed=260001 唯一一次
X-source Tier-1 重跑后，跨模型差仅从 `0.310434 dB / -2.248201 rad` 变为
`0.310371 dB / -2.248201 rad`，因此一横向维 PE 的正式 Bellhop comparator 应用
`X`，但 Kirchhoff phase-screen 与 local-specular Reflect2D 的模型差异结论不变。
详见 `reports/bellhop_source_geometry_rx_audit_report.md`。

当前 PE--Bellhop validation 主链已把 line-source `X` 写入 incident、Tier-1、
frequency-extension、ensemble 和 Stage-4 入口，并把 source geometry 纳入 ensemble
case fingerprint，旧 `R` cache 不能冒充 `X`。在相同 4 kHz incident-plane 配置下
重跑后，PE--Bellhop 的 M95 幅度 P95 从旧 R 的 `0.15794309 dB` 降至
`0.0019036364 dB`，13/13 gates 全部通过；M99 phase RMS 仍为
`0.00610483 rad`，PE--AS L2 仍约 `3.73e-13`。X 不使用旧 point-source 全局常数，
且整个比较仍按轴上归一化，不包含经验拟合。既有 Stage-2/3/4 报告保留其 R-source
历史 provenance；只有将来显式重跑后才会生成 X-tagged 统计结果。整合与归档见
`reports/pe_bellhop_line_source_integration_report.md`。

Stage 1B/1C 已冻结。生产二维 PE 在 y 方向复制同一个固定 `eta(x)`，
`ny=256/512` 的变化仅为 `4.36e-6 dB / 4.49e-6 rad`；相对一横向维桥接，
dimensionality sensitivity 约 `−0.1589 dB / 0.01251 rad`（复误差 `0.02196`）。
4 kHz fixed-PM 总结判定为 `PASS_WITH_MODEL_DISCREPANCY`：PE/Bellhop ratio
差异为 `0.3104 dB`、`−2.2482 rad`，明显高于维度误差，归类为 Kirchhoff
phase-screen 与 Bellhop local-specular 模型差异。详情见
`reports/pe_bellhop_fixed_pm_4khz_comparison_report.md`；随后已完成 4/6/8 kHz
频率扩展。

Stage 2 频率扩展已通过，状态为 `PASS_WITH_MODEL_DISCREPANCY`：固定
seed-260001 band-limited PM 在 4/6/8 kHz 分别运行一横向 PE 与 Bellhop
internal-wall flat/rough 配对。三点均满足有限场、wall residual `7.14e-15 m`、
零 grazing、压力释放 phase error `7.11e-15 rad`、q residual 与 p/q rotation
误差为零；跨模型 TL 差异为 `0.3104/0.3083/0.3048 dB`。频率趋势只作
模型 diagnostic，不由三点估计 group delay。Bellhop 使用预声明 ±15° sector
和 5001 beams，以避开 validation overlay 的低振幅尾束终止，不涉及 source
拟合。详情见 `reports/pe_bellhop_fixed_pm_frequency_extension_report.md`；随后已完成
Stage 3A 配对多 realization smoke 与 Stage 3B 独立系数幅度 ensemble。

Stage 3 已完成 8 个 paired fixed-band seed 的 4 kHz 统计 smoke，状态为
`PASS_WITH_LIMITS`。每个 profile 保持 canonical PM 每模态谱幅度和 Kmax，
只由确定性 seed 改变相位；Bellhop wall residual 均约 `7.1e-15 m`、零
grazing、压力释放 phase 与 p/q state 检查全部通过。均值反射粗糙功率为 PE
`1.0707`、Bellhop `0.9271`。这不是独立 PM 幅度抽样的最终 Monte Carlo，
backward-range amplitude 仍只作 diagnostic；详情见
`reports/pe_bellhop_pm_ensemble_comparison_report.md`。随后补充的 Stage 3B 使用
同一 `Sk_m3` 频带独立抽样 Fourier 系数幅度，8 个 seed、统一 5001 beams 也通过
全部结构检查；每个 rough case 均记录 5001/5001 wall hits。flat、canonical rough
和其他缓存只在配置指纹、输入新旧关系及输出哈希满足对应校验时复用，避免旧
1001-beam case 混入；结果见下方独立报告。

完整阶段索引和冻结项汇总见
`reports/pe_bellhop_pm_complete_execution_summary.md`。后续若进入正式粗糙面
统计，仍需保持同一 profile provenance、`G=H_ref,rough/H_ref,flat` 主量和既定
数值/维度护栏，不得用 source 或幅相 renormalization 消除 PE--Bellhop 模型差异。
若需要独立 PM 幅度 ensemble，应运行 validation-only
`validate_pe_bellhop_pm_amplitude_ensemble`；它固定同一 `Sk_m3` 频带，
已完成 `seed=260001:260008`、5001-beam 的 reduced first ensemble，
只改变确定性 Gaussian Fourier 系数抽样，不改变 PE/Bellhop 核心物理。
当前 mean delta TL 为 `0.320686 dB`，circular mean phase difference 为
`-0.689310 rad`；所有 case-fingerprint、uniform-numerics 和 reflection-success
检查通过。
其结果仍不是充分收敛的 ocean Monte Carlo，报告见
`reports/pe_bellhop_pm_amplitude_ensemble_report.md`。

2026-08-27 起，正常坐标粗糙海面 Bellhop 任务按用户要求重建为平均水深
100 m、物理源离底 0 m、Rx 深度 3 m。海底测试模型暂定均匀流体半空间：
1800 m/s、2000 kg/m³、0.5 dB/波长；这是可复现测试参数，不是实测底质。
原 500 m 粗糙海面试验和完整声线图已退役，不能作为 100 m 全多途信道的证据。
旧运行与绘图代码保留，但不得按其 500 m 默认设置继续当前任务。

2026-10-10新增独立的5000 m Munk、Tx4950 m、10 kHz图12(b)复现尝试：
`validate_bellhop_fig12b('full')`，见[假设、34例结果及限制](reports/bellhop_fig12b_reproduction_report.md)。
原生点声源把零水平距离整列置零；严格垂直声线成功不代表相干幅相已验证。
全角近轴幅值未通过解析1/R控制，窄角扇声束收敛仍不足，当前复场不可直接
当作已验证通信H_f。本任务不构造CIR、不做调制/BER，也不改变下述100 m安装。

随后按用户限定追加三种原生束/源方案，共15次声束、步长及均匀控制运行；
仍未可信复现，已停止参数搜索。外观最接近的Cerveny方案左右不对称且幅值
失败；线源替换也未解决轴列。见[最终停止报告](reports/bellhop_fig12b_three_attempts_stop_report.md)。
这些新增诊断同样不构成可直接使用的通信信道响应。

**最新决定覆盖上述贴底安装：用户已接受实际离底1 m，Tx改为99 m，Rx仍为3 m。**
以下epsilon诊断仅保留为历史，不再继续。新入口
`generate_bellhop_100m_source1m_visuals_vertical`仅用于Bellhop代表案例及声线图，
不改变PE默认配置。接收A/E与重建R使用12 m水平范围；展示发射扇区使用120 m，
两者不能混称为同范围多途全集。完整扇区保留触面后的回程、底反射及已保存多次
反射；Bellhop终止规则及轨迹存储上限仍适用，不延长或补造射线。

用户已接受以 epsilon=1e-3/1e-4/1e-5 m 检查水侧极限，物理离底字段仍为 0。
安装版本不仅不支持精确边界源，还会按 SSP 下界自动上移过近的源。本轮使用
显式平底 BTY 固定物理海底 100 m，只延伸均匀 SSP 0.1 m，绕过输入层限制；
不是加深海底。整体输入坐标下移 1 m，海面、海底、Tx、Rx 一起移动，绘图应
再还原至实际 0–100 m。请求最小 epsilon 的 ARR 记录值约为 8e-6 m。

4 kHz、89.5°、10001 波束、步长 0.02 m 的平面/固定粗糙面六项检查已完成。
1e-4→1e-5 m 的四个直达/一次海面类别均满足 0.01 dB、0.002 rad、1 us；
粗糙面高阶含底类别未通过，所以没有继续角度/频率扩展。粗糙面全部 ARR
合成响应变化为 -0.0385556 dB TL、+0.0103577 rad；不能宣称全信道已收敛。
平面 D+S 因强相消，合成变化也大于单条路径误差，详见报告。

随后对同一粗糙面两个 epsilon 运行 Bellhop C 模式并与 ARR 重建比较，C/ARR
复比值的差异仅为 0.000533 dB TL、2.29e-5 rad。这说明高阶路径的 epsilon
敏感性已存在于相干场本身，而不是 ARR 后处理单独制造；该结果仍不批准
完整含底多途。报告见 `reports/bellhop_100m_arrival_field_audit_report.md`。

ARR 按安装版本公式重建 U_native=sum A exp(i phi-i omega tau)，再取共轭
得到项目 exp(-i omega t) convention。Delta TL=TL_new-TL_old，幅度变化
为其相反数；本轮 arrival 加权时延不是宽带群时延。该实验的 X 线源不等于
PE Gaussian 源，绝对输出级不能混比。
报告见 `reports/bellhop_100m_boundary_source_report.md` 和
`reports/bellhop_100m_rough_surface_report.md`；计划及归档清单见
`reports/bellhop_100m_rebaseline_plan.md`。`cash/` 不参与正常读取或自动恢复。

后续图仍区分 A 合并后的 arrival 与 E 接收贡献波束，E 轨迹用独立 R 重建，
不把射线端点人为接到 Rx。比较时分别报告直达、一次海面和含海底的多途；
Bellhop 含海底的全响应不能直接等同于目前 PE 的直达加海面反射响应。

## 4. 海面边界、SSA 与 bubble

海面实现位于 `src/surface/`。`paramsV.surface_boundary_model` 当前允许：

- `kirchhoff_spatial`：生成显式 PM 海面 `eta(x,y)`，在空间域施加 Kirchhoff 相位反射；这是当前公共默认。
- `kirchhoff_kdomain`：同一类显式海面反射在 k 域组织，适合与 joint 统计模型对照。
- `kirchhoff_kstat`：不逐次传播显式海面，而按同一 PM/Kirchhoff 相位屏定义建立随机残差及其跨频统计。
- `ssa_stat_kernel`：第一阶 Dirichlet SSA/微扰极限的统计核诊断分支，不等于 NLSSA 或完整高阶 SSA。

joint-kstat 的随机变量是

\[
\delta G_i=R_{0,i}e^{i\alpha_i\eta}-R_{\mathrm{coh},i}.
\]

它仍包含海面高度引起的相位，但“joint”表示同一随机海面在多个频率上的联合统计：代码构造或采样跨频 `C_deltaG` 与 `P_deltaG`，并非只做单频确定性相位计算。当前 `C_deltaG/P_deltaG` 已包含 `R0`，接收权重不得再次乘 `R0`。

SSA1 使用压力释放/Dirichlet 几何因子

\[
G_{\mathrm{SSA1}}(K,K';f)=4\gamma(K,f)\gamma(K',f),
\]

并只保留已实现的一阶统计散射、传播波支和已声明的周期或零填充卷积。它不包括 NLSSA、多次散射、阻抗边界或实验标定。2k 诊断位于 `scripts/validation/validate_2k_phase_approx.m` 和 `scripts/validation/validate_2k_ocean_spectra.m`，只验证相位系数近似与最低阶 SSA 趋势，不等于 PE、KStat、真实海洋散射或高阶 SSA 验证。

气泡实现位于 `src/bubble/`，支持 `off`、`level0_empirical`、`hall1d` 和已有 plume 配置接口。公共默认 `enable_bubbles=false`。Hall1D/Li2009 相关结果仍是部分复现：横向网格收敛、独立绝对幅度和完整外部交叉验证尚未全部闭环。

## 5. cached forward、精确离散伴随与解析 C/P

`cached forward` 是把固定频率、固定网格、固定均匀环境的 Surface→Rx PE 步进因子预先缓存，再对很多海面 realization 重复执行同一个前向算子。它不是近似传播器，而是固定路径的高可信回归 oracle。

对应离散伴随位于 `src/receiver/`：

```matlab
q_k = conj(fr) .* fft2(conj(screen) .* ifft2(conj(fr) .* q_k));
```

它从接收平面单位源开始，反向遍历深度步，得到

\[
q_i=A_i^Hr,\qquad
a_i=\operatorname{conj}(\psi_{\mathrm{inc},i})\odot q_i.
\]

`q` 是接收灵敏度核，不是真实反向声压场，也不是逆传播。对固定单接收机，可用 `q_i^H psi_surface` 精确替代每条 realization 的 Surface→Rx PE。当前 full 验证中 exact adjoint、投影、cached/adjoint、F=9/F=64 均通过；总体和单样本差异处于双精度量级。范围限定为 uniform、CPU double、固定 Tx/Rx、最近网格点采样、无 bubble/Doppler、固定 PE/PM 网格。

PM 大网格到 PE 小网格使用现有中央裁剪 `E`，统计收缩必须把 PE 权重通过 `E^H` 零嵌回 PM 周期网格。接收统计为

\[
C_H(i,j)=\widetilde a_i^H C_{\delta G,ij}\widetilde a_j,\qquad
P_H(i,j)=\widetilde a_i^H P_{\delta G,ij}\widetilde a_j^*.
\]

`dense` 仅在 PM 不超过 32² 时作基准；正式方法逐频率对流式计算 FFT 收缩，不保存空间 `F^2` 协方差块。direct-DSP 统计由 reduced 统计经 `D=diag(exp(-i2*pi*f*Delta tau0))` 转换：`C_dsp=D*C_red*D'`，`P_dsp=D*P_red*D.'`。

## 6. 条件统计生成器

条件模型用已验证接收端样本或解析统计估计复均值、协方差 `C`、伪协方差 `P`，支持 full 与低秩采样。schema 2.x 明确记录 `target_reference='direct_dsp'`、几何和频率轴。旧 schema 1.x 只能在有可靠几何时迁移；只有旧总信道或缺少几何时不得猜测拆分或旋转。

U=5/U=8、F=64 条件模型及 two-node 通信验证已纳入 2026-07-22 的 PASS 发布候选。其角色是经物理和统计验证后承担大规模通信 Monte Carlo，不替代 PE 物理回归，也不能插值成未验证风速的连续模型。

## 7. CIR、F=64/F=65、LFM 与通信链

F 表示频率采样点数。当前 F=64 主轴覆盖 4–8 kHz，共 64 个等间隔频点，用于宽带统计、条件模型和通信验证。它并不表示 64 条独立信道。

相位验收另用 F=65、间隔 62.5 Hz，其无模糊时延窗为 16 ms，可显示标准几何的 4 ms 相对反射时延。F=9 smoke 的间隔为 0.5 kHz，无模糊窗仅 2 ms，4 ms 会混叠为零相位，不能用于相位验收。

`src/channel/build_channel_cir_vertical.m` 区分 direct-DSP `H(f)` 的 IFFT 和绝对物理相量的展示性重构。时间窗移动只能作用于完整总信道，禁止分别对齐直达与反射。

LFM 是获得 `H(f)` 后的无噪声线性探针：频域输入乘以信道频响，再变换到时域。LFM 输出功率描述信道后扫频脉冲能量随时间的分布；匹配滤波输出功率描述接收 LFM 与已知发射模板相关后的压缩峰和旁瓣。二者可比较 reflected-only 的时延扩展与统计功率，但不构成新的 PE 空间源，也不自动包含同步、均衡或噪声。

通信模块位于 `src/communication/`。外部没有项目 metadata 的 `H(f)` 默认视为已经 DSP-ready，不静默旋转。

## 8. 公共输入输出语义

公共入口为 `vertical_channel_model(paramsV)`，实现位于 `src/channel/vertical_channel_model_impl.m`。原有字段保持：

- `output.H_direct_f`、`H_reflect_f`、`H_f`：按 `channel_phase_reference` 选择后的频响，默认 direct-DSP；
- `output.h_direct`、`h_reflect`、`h_total`：上述数组在 `idx_f_ref` 的标量；
- `output.f_axis`、`idx_f_ref`：频率轴及参考频点；
- `output.H_direct_reduced_f`、`H_reflect_reduced_f`、`H_total_reduced_f`：约化 PE 分量，用于兼容/诊断；
- `output.H_direct_physical_f`、`H_reflect_physical_f`、`H_physical_f`：绝对物理相量展示；
- `output.phase_reference_meta`、`roughness_meta`、`bubble_meta`、`config`：语义和配置元数据。

无论选择何种公共表示，都必须保持频域和参考频点的分量闭合。不得删除或重命名现有输出字段。

## 9. 验证证据边界

已通过：

- direct-DSP 载波相位发布候选、F=65 名义 4 ms 时延和旧条件模型迁移闭合；
- uniform CPU-double 的精确离散伴随、接收投影、PM 零嵌、dense/FFT `C/P`；
- F=9/4096 与 F=64/512 的解析—样本统计在 split-sample floor 内；
- U=5/U=8 条件模型、two-node 通信和已有 PE 图册；
- 固定海面完整反射链 192.1875 m/no-sponge，以及 4 kHz 随机集合 15/15。

开放或不完整：

- 公共默认窗口仍是兼容值，尚未切换到严格推荐窗口；
- 随机海面宽带成对收敛没有完成，不能给出 ensemble 最差群时延；
- Bellhop 的绝对幅度、点源无限孔径和源归一化仍开放；
- Li2009/SSA 是部分复现，横向网格与高阶物理验证不完整；
- 伴随统计原型不覆盖 layered、bubble、Doppler、GPU、多接收机和插值接收。

判断新结果时应把 `PASS`（满足预设门限）、`CONVERGED`（在已测数值配置收敛）、`OPEN`、`INCOMPLETE` 和“仅诊断”分开。不能用 total channel 掩盖 reflected-only 差异，也不能把 explicit Kirchhoff 与 joint-kstat 的独立 realization 做像素级误差验收。

## 10. 推荐角色与配置

- 公共 PE：通用传播入口。
- cached forward：固定路径高可信数值 oracle。
- 伴随投影：固定环境快速、精确地产生单接收端 realization。
- 解析 FFT `C/P`：无需接收端 Monte Carlo 直接获得二阶统计。
- 条件统计生成器：验证通过后的大规模通信抽样。

研究级严格窗口当前优先采用 direct `160 m/no-sponge`、完整反射 `192.1875 m/no-sponge`；若使用公共 50 m/default sponge，必须明确它是兼容默认而非近期窗口收敛验证的推荐配置。

## PE--Bellhop PM 模型差异统计（Stage 4）

Stage 4 只在固定 4 kHz、冻结数值设置下，对 Stage 3B 的独立 Gaussian
系数振幅 realization 做 PE--Bellhop 统计，不改变任一侧的核心传播或
接收公式。当前保留 50 个完整 seed（`260001:260050`），报告为
`reports/pe_bellhop_pm_model_discrepancy_statistics_report.md`，逐 seed 与
bootstrap/相关性表位于
`results/validation/pe_bellhop_pm_model_discrepancy_statistics/`。

M=50 的平均 `delta TL=-0.00748 dB`、圆均值相位差 `-0.37534 rad`；所有
几何、非 grazing、beam-state、正 range、PE seam/edge 防护均通过。24→32
的均值/功率/相位工程门通过，但 bootstrap 均值 delta-TL 半宽为
`0.25835 dB`，略超 `0.25 dB`，故状态仍记为
`PRELIMINARY_MODEL_DISCREPANCY`。该结果用于判断模型差异的稳定性，不用于
校准或强行消除 native backward-range 幅度偏差。

## Full Kirchhoff spectral separability audit — validation boundary (2026-09-30)

### 小规模跨频一致性（2026-10-08；最新扩频结果）

固定 canonical seed260001、12个 Fourier 模态、Kmax=.471238898 rad/m、
C2嵌入和 Tx100/Rx3 的研究夹具，已完成4/5/6/7/8 kHz sparse gates，再完成
4–8 kHz、F65、delta-f62.5 Hz 的 dense diagnostic，结论为 **PASS**。
冻结3/8从6 kHz起不够；3/10在8 kHz仍略超1e-3，因此采用整频带统一3/12，
padding2不变，没有逐频 tuning。最大 Fast–Full 加权/whole-grid复误差为
2.30670e-4 / 2.87485e-4；三个 Full–BIE 物理锚点为8.00096e-4、
5.30379e-4、4.08269e-4，物理范围仍限定在冻结M99核心孔径。

原生物理对应的 direct-DSP 相位与群时延经上节桥接保持一致：最大相对相位
1.74394e-4 rad，群时延差70.19 ns，reflected-only CIR复误差7.17177e-5。
观测群时延约3.595 ms，采样CIR峰约3.69231 ms，二者受有限带宽与时间bin
影响，不能逐路径对齐。1 kHz sparse spacing只有1 ms窗，不能认证时延；
F65的16 ms窗覆盖名义4 ms超额时延，但有限积分边缘允许约53.129 ms的
几何长路径，尚未证明弱散射尾部无混叠。这仍是phase-consistency diagnostic，
不是最终支持收敛的CIR。

三个频率的padding2/4、height nodes33/65和surface4097/8193护栏通过；
65点最小cold/warm加速约9.83/41.94，内存按named arrays和块估计记录。
Reference和生产PE/通信代码不变。允许下一阶段在此条件和显式phase metadata
下开展固定PM宽带及小规模realization数据；未外推到PM总体、SSP、3D、
time-varying或生产API集成。输出与完整限制见对应多频报告。

### PM 物理适用域的受控审查（2026-10-07；最新）

在冻结 4 kHz、c=1500 m/s、Tx=100 m、Rx=3 m、Gaussian 离散源和 C2 嵌入后，
已围绕同一个 canonical realization 完成八个嵌套 cutoff 与七个受控 Fourier
对照。共有系数完全相同；高度匹配、照明 slope 匹配及固定 Kmax 的新增模态
幅度序列分别检查混杂因素，不把接收采样点视为独立 PM realization。

当前 master-prefix 家族在 Kmax=.706858 rad/m 的原始复误差为8.5275e-4，
.746128 rad/m 为1.19749e-3。边界端点、等 slope 高 K 及75%新增模态幅度点
使用 Full4097/8193、BIE ppw/order8→10 与 L60/L70 护栏，合计变化约6.8e-6。
这些是离散、单 realization、固定 M99 核心孔径下的数值护栏结论，不能外推为
所有 PM 海面或整个接收线的统一 Kmax 门槛。

相同 Kmax 的幅度序列从7.96e-4 增至1.20e-3；等高度与等照明 slope 对照
支持把照明加权曲率/k 和 high-K 幅度/谱结构共同作为候选适用域描述量。
尚未识别唯一控制参数，也未拟合经验公式。幅度误差约1.3e-4，主要差异来自
原始相位起伏；没有相位对齐或 complex scalar fitting。共同核心复谱受513个
BIE接收采样的 Nyquist 限制，保留全部复数据并明确区分孔径谱与远场通量。

**保留有限适用域内的 Full FK/Fast FK 路线。** 若目标必须覆盖当前超标的
高 K/幅度结构，下一阶段先做有限 SSA/SSA2 物理模型比较。Fast rank3/8 与
padding2 的算法成果保持冻结，本轮没有扩频、Monte Carlo、BER、SSP、3D
或生产集成，paramsV/output/H_f 和通信方法保持既有语义。详见
`reports/full_kirchhoff_pm_validity_boundary_report.md`；前一阶段 PM pilot 的
高 cutoff 物理失败结论仍成立，其“适用域尚待研究”建议已由本轮受控审查接续。

### 固定 PM pilot 的最新边界（2026-10-07）

验证专用 pilot 已完成三个固定 seed 和一组嵌套物理 cutoff，完整保留高度、
解析法向、ds、强相位与非局部 Green 作用。冻结 3/6 rank 在 PM 高度界
A=.5 m 下失败；统一 3/8 在全部五个嵌入及无包络 PM 案例上通过 FastFK–FullFK
的 1e-3 算法 fidelity。物理 Kmax=.275/.471/.746 rad/m，与 dx 决定的
Nyquist=107.23 rad/m 明确分开。没有每个 seed 单独拟合 rank 或相位。

**整体 PM pilot 未通过预设物理门槛。** FullFK–BIE 的固定 M99 加权误差为
canonical 8.00e-4、high-cutoff 1.199e-3；BIE 与 Full 积分精化已将数值误差压至
明显更小的水平。BIE 使用共同 C2 嵌入海面，其物理结论不覆盖整个接收线
的窗口边缘或无限周期海面。算法精度与模型物理精度仍分别报告。

下一步应先明确 PM 下 Full FK 的适用边界或比较独立 SSA/SSA2 reference，
本轮不进入 multi-frequency、Monte Carlo、生产 PE 或通信/BER。详见
`reports/spectral_fast_kirchhoff_pm_pilot_report.md`。生产 paramsV/output/H_f
语义与本指南通信模型保持不变；下面的 prototype-only PM 建议为上一阶段记录。

### 后续 FFT prototype 验证边界（2026-10-07）

验证专用 `solve_spectral_fast_kirchhoff_2d` 已用全局入射分离与出射高度分离，
把完整 lag Green/法向导数作用转为 FFT 线性卷积；Prototype-B 保留真实海面、
法向、弧长和两项 Green 积分。它没有替换本指南描述的生产信道/通信算法。
选定入射/出射 rank 分别为 3/6；单独的入射 rank 2 不是完整核 rank。
4 kHz、三个确定性正弦面、sigma=.3 m、Rx=3 m 的原生 FastFK–FullFK
复误差均小于 1e-4。padding=2/4 对应同一完整线性卷积，padding=1 会回绕；
验证中没有 sponge 或逐案例经验相位/尺度拟合。

参见 `reports/spectral_fast_kirchhoff_prototype_method.md` 的严格/低秩/几何简化
分类，以及 `reports/spectral_fast_kirchhoff_prototype_report.md` 的 rank、padding、
出射谱、CPU/内存和冷/热启动速度证据。FullFK–BIE 的物理模型有效性与
FastFK–FullFK 的算法 fidelity 仍分开。固定 PM pilot 仅为下一阶段建议；
宽带、随机海面、3D 与 BER 尚未验证，公开 paramsV/output/H_f 不变。

### 上阶段审计记录（2026-09-30，后续 FFT 证据见上节）

The independent 2-D validation in
`reports/full_kirchhoff_spectral_rank_audit_report.md` factors
`exp(-i*kz*eta)=exp(-i*k*eta)*exp(i*(k-kz)*eta)` and preserves the strong
carrier phase exactly. Its weighted SVD uses the square root of saved incident
spectral energy as a column weight and relative Frobenius matrix error.
At 4 kHz, the saved deterministic near-vertical cases need rank 1–2 at 1e-3;
three complete Full FK contractions at Rx=3 m also pass 1e-3 at weighted rank 2.
The 99% angular interval describes effective bandwidth; full saved-spectrum
tails are retained for the strict engineering reconstruction checks.

This is evidence for a future spectral prototype, not a replacement of the
production surface boundary or a communication-channel validation. Outgoing
Green-kernel separability, O(r*N*log(N)) performance, PM surfaces, wider angles
and wideband communication remain unverified. All public paramsV/output/H_f
and results.channel semantics remain unchanged; no new configuration option
or physical model is introduced. The existing PE/local screen, BIE and Full FK
reference implementations were not edited. Run the audit after
`setup_vertical_project` via `validate_full_kirchhoff_spectral_rank_audit()`.
