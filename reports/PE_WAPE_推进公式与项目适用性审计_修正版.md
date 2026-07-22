# PE/WAPE 推进公式与项目适用性审计（修正版）

审计日期：2026-07-22  
审计基准：当前工作区源码优先；项目内 Markdown 和既有结果仅作为实现说明与验证记录，不覆盖源码事实。  
重点文件：`vertical_wape_propagator.m`、`vertical_channel_model.m`、`pm_surface_boundary_model.m`、`comm_main_vertical_psk.m`、`build_physical_cir_vertical.m`、`build_communication_taps_vertical.m`、气泡介质模块及现有验证脚本/报告。

## 1. 审批结论

原审计稿的核心判断大体正确，但不应原样批准。修正后的结论如下。

1. **自由推进因子的代数形式正确。** 代码中的两个半步因子合并后等于
$\exp\{i d[\sqrt{k_0^2-\kappa^2}-k_0]\}$。因此在 $c=c_0$、无介质/吸收屏、采用当前离散 FFT 网格时，它就是去除名义轴向载波后的单向角谱推进，不是窄角 Fresnel 展开。对 $\kappa>k_0$ 的主值复平方根也给出正确的倏逝衰减。

2. **“精确”只能修饰自由因子或均匀参考介质案例。** 实际每步还会乘海绵、声速和气泡屏；一旦 $c\ne c_0$ 或屏随推进坐标变化，整个算法仍是 Thomson–Chapman/Feit–Fleck 型分步近似，而不是一般非均匀 Helmholtz 方程的精确解。它还是单向模型，不生成由介质梯度或界面引起的反向传播。

3. **当前公共宽带输出存在真实的相位参考问题，但应准确表述。** `H_direct_f` 是相对直达段名义轴向距离 $d_{\mathrm{dir}}=z_{\mathrm{tx}}-z_{\mathrm{rx}}$ 去载波后的包络；`H_reflect_f` 是相对两段名义轴向距离 $d_{\mathrm{ref}}=z_{\mathrm{tx}}+z_{\mathrm{rx}}$ 去载波后的包络。两者使用不同的局部载波参考。当前
   `H_f = H_direct_f + H_reflect_f`
   是通过的软件代数不变量，但在补齐共同相位参考以前，不是可直接作物理相干叠加、再直接 IFFT 的总传递函数。

4. **主通信脚本确实直接消费上述 `H_f`。** `comm_main_vertical_psk.m` 将其移到参考频率附近、插值并执行 `ifft(ifftshift(...))`，没有先把直达和反射分量转换到统一的声学或信号处理相位约定。因此当前 BER/SER 可用于既有 reduced-envelope 算法链内部的相对比较，但不能据此声称已经验证了真实直达—海面反射几何时延和绝对物理多径性能。

5. **项目已有独立 CIR/tap 工具，但它们没有自动修复公共 PE 输出语义。** `build_physical_cir_vertical.m` 假设输入频响采用“正时延对应负相位斜率”的 IFFT 约定；`build_communication_taps_vertical.m` 复现通信 tap 转换并改进了圆周能量窗。二者都要求调用方先提供与其约定一致的频响。它们的存在不能证明当前 `vertical_channel_model().H_f` 已经是物理频响。

6. **纵向非均匀介质仍缺少独立验证。** 主循环和 `local_march_field` 都在步终点 `z_curr` 计算介质/气泡屏。对随推进坐标变化的屏，这不是非自治 Strang 分裂所需的中点取样；对仅随深度变化的标量屏，主要表现为终点矩形积分误差。气泡功能默认关闭，因此该问题不是默认均匀水体回归失败，而是启用 layered/bubble 后的验证缺口。

7. **默认粗糙海面模型是具体 realization 的 Kirchhoff 相位屏，不是严格起伏压力释放边界解。** `kirchhoff_spatial` 在法向模式下执行
   $\Psi_{\mathrm{ref}}=R_0\exp(i2k_0\eta)\Psi_{\mathrm{inc}}$，默认 $R_0=-1$。平面极限正确，但模型忽略非局域边界耦合、遮蔽和多次散射。项目另有 `kirchhoff_kdomain`、`kirchhoff_kstat` 和 `ssa_stat_kernel/ssa1_geometry` 等可选分支；这些分支扩大了诊断与统计生成能力，却仍不是完整 SSA2/NLSSA、积分方程或实验标定的粗糙面求解器。

8. **均匀介质、平面海面和软件不变量已有较强证据；开放边界、绝对幅度、粗糙面、气泡和端到端物理通信尚未全部通过。** 现有严格 Bellhop/收敛矩阵的总体状态仍为 `false`，失败项来自跨几何幅度门槛和改变窗口及海绵布局后的反射相位差，不应改写为“全部验证通过”。

建议将主传播器描述为：

> 经均匀参考介质离散角谱与平面压力释放海面案例验证的、恒密度、单向、近轴向 Thomson–Chapman 型 WAPE 复包络传播器；默认外接一次 Kirchhoff 粗糙海面相位屏，并提供若干统计散射与快速条件信道研究分支。当前公共直达/反射频响仍处在不同局部载波参考下，尚未形成由主通信入口直接消费且端到端闭合的物理宽带总传递函数。

## 2. 对原稿的主要修正

| 原稿表述 | 审批结果 | 修正 |
|---|---|---|
| “均匀介质中推进核严格等于单向角谱传播” | 有条件通过 | 只对自由因子，或 $c=c_0$ 且无介质/吸收屏时成立；默认正海绵会使完整步进算子不再等于纯自由角谱传播。 |
| “公共 `H_f` 不合格” | 核心问题成立，措辞需收窄 | `H_direct_f`、`H_reflect_f` 作为各自局部参考下的包络是有定义的；有问题的是把二者直接当成同一相位参考下的物理总频响并送入 IFFT。 |
| “当前公共 CIR/通信工具采用负斜率与 IFFT” | 通过，但需补充 | 项目同时存在验证专用的 `exp(-i\omega t)`、正斜率+FFT 约定。两套约定都可自洽，当前缺陷是转换关系没有在公共 API 中固化。 |
| “默认网格在 8 kHz 约 3.84 点/波长” | 表述含混 | 这是 `vertical_channel_model` 的 `50 m/1024` 公共默认值；`comm_main_vertical_psk` 使用 `50 m/256`，只有约 0.96 个全波长网格点。PE 横向网格应按所需 $\kappa$ 支撑而非仅按全波长判断，但粗糙面产生的高横向波数更需要独立检查。 |
| “默认粗糙面分支是 Kirchhoff 相位屏” | 通过 | 只适用于默认 `kirchhoff_spatial`；项目另有显式 k-domain、Kirchhoff 统计相位屏和 SSA1 统计几何核，不应把所有分支统称为同一个点乘模型。 |
| “项目尚无物理 CIR/快速统计信道能力” | 不成立或不完整 | 已有物理 CIR/tap 工具、cached/joint-kstat 和条件 $\mu/C/P$ 生成器；但这些能力不等于公共 PE `H_f` 的相位语义已修复，也不等于粗糙面物理已实验标定。 |
| 原稿中的 20 m→3 m IFFT 峰值表 | 不纳入正式证据 | 当前仓库没有与该表绑定的保存脚本或结果文件。保留解析路径时延结论，不把未归档的临时运行数字当成项目验证结果。 |
| 建议“让公共 `H_f` 改成物理频响” | 方向正确但可能破坏兼容性 | 不应静默改写既有字段。应新增明确命名的 common-reference/physical/signal-convention 字段或转换函数，并保留原字段直至下游迁移完成。 |

## 3. 当前推进公式与代码对应

令推进距离为 $s$，横向坐标为 $(x,y)$，$\kappa^2=k_x^2+k_y^2$，$k_0=2\pi f/c_0$。代码保存谱域复包络 $\widehat\Psi$，每步执行

$$
\widehat\Psi_{j+1}
=D_f\,\mathcal F\!\left[
D_s\,\mathcal F^{-1}(D_f\widehat\Psi_j)
\right],
$$

其中半步自由因子为

$$
D_f(\kappa)=\exp\!\left[-\frac{i d}{2}
\frac{\kappa^2}{\sqrt{k_0^2-\kappa^2}+k_0}\right].
$$

利用

$$
-\frac{\kappa^2}{\sqrt{k_0^2-\kappa^2}+k_0}
=\sqrt{k_0^2-\kappa^2}-k_0,
$$

可得

$$
D_f^2=\exp\{i d[k_z(\kappa)-k_0]\},\qquad
k_z=\sqrt{k_0^2-\kappa^2}.
$$

MATLAB 主值复平方根使 $\kappa>k_0$ 时 $k_z=+i\gamma$，所以倏逝分量随正推进距离衰减。该结论与 `vertical_wape_propagator.m` 的实现一致。

介质与吸收屏为

$$
D_s=\exp\left[-ik_0d\left(
\frac{c-c_0}{c}-i\frac{\alpha}{k_0}
\right)\right]
=\exp[i k_0d(n-1)]\exp(-\alpha d),
\qquad n=\frac{c_0}{c}.
$$

因此：

- `U_real_xy=(c_eff_xy-c0)./c_eff_xy` 等于 $1-n$，相位符号正确；
- 正的幅度衰减系数 $\alpha\,[\mathrm{Np/m}]$ 产生 $\exp(-\alpha d)$；
- 当前 $\alpha$ 由横向海绵和可选气泡衰减组成，不包含常规海水体吸收；
- `old/propWAPErev11.m` 也使用 `(cin-c0)./cin`；WHOI 2006 报告正文写的是 $(c-c_0)/c_0$。两者在小声速扰动下一阶相同，但不能在文档中当作完全相同的定义。按当前屏要实现 $n-1$ 的推导，代码中的除以 $c$ 更直接一致。

## 4. WAPE 近似的适用边界

对恒密度标量 Helmholtz 方程

$$
(\partial_s^2+\nabla_\perp^2+k_0^2n^2)p=0,
$$

在 $p=\Psi e^{ik_0s}$、时间因子 $e^{-i\omega t}$ 下，形式上的单向算子为

$$
\partial_s\Psi=ik_0\left[
\sqrt{n^2+k_0^{-2}\nabla_\perp^2}-1
\right]\Psi.
$$

当前分步法对应近似

$$
\sqrt{n^2+X}\approx \sqrt{1+X}+(n-1).
$$

它在 $n=1$ 时保留完整自由空间平方根，也在轴向 $X=0$ 时给出正确局部波数；混合误差随折射率偏差与传播角共同增长，量级可写成 $O[(n-1)\sin^2\theta]$。所以“自由因子宽角精确”不能替代对分层声速、气泡相速异常和较大偏轴角的验证。

对冻结且不对易的自由算子 $A$ 与屏算子 $B$，`半步 A—整步 B—半步 A` 是 Strang 结构。若 $B=B(s)$ 随推进坐标变化，要保持非自治问题的二阶性质，通常应使用步中点 $B(s+d/2)$ 或等价对称离散。当前两条传播路径都在 `z_curr=z_start+j*dz_step` 取屏：

- 仅随深度变化、横向均匀的声速/衰减屏与自由算子可交换，此时主要误差是终点矩形积分；
- 横向非均匀气泡羽流等情形还会同时出现非对易分裂误差；
- 上、下行分别重新确定步数，不能假设终点误差必然抵消。

以指数尺度 $L=0.4\,\mathrm{m}$ 和名义最大步长 $d=0.5\lambda$ 为示意，单个上行步对 $e^{-z/L}$ 的终点矩形积分相对精确积分之比为

$$
\frac{d/L}{1-e^{-d/L}}.
$$

在 4、6、8 kHz 时名义值约为 1.253、1.164、1.122。实际代码会用 `ceil` 后重新均分整段距离，所以真实每步误差略依赖具体段长；这些数字是局部量级估计，不是当前完整气泡案例的实测误差。

## 5. 宽带相位与 `H_f` 语义

### 5.1 当前字段实际代表什么

自由因子只推进去掉 $e^{ik_0s}$ 后的包络。对当前平均海面 $z=0$：

$$
d_{\mathrm{dir}}=z_{\mathrm{tx}}-z_{\mathrm{rx}},\qquad
d_{\mathrm{ref}}=z_{\mathrm{tx}}+z_{\mathrm{rx}}.
$$

在项目验证采用的声学时间约定 $e^{-i\omega t}$ 下，恢复名义轴向载波应为

$$
H_{\mathrm{dir}}^{(a)}(f)
=H_{\mathrm{dir}}^{\mathrm{red}}(f)e^{+i2\pi f d_{\mathrm{dir}}/c_0},
$$

$$
H_{\mathrm{ref}}^{(a)}(f)
=H_{\mathrm{ref}}^{\mathrm{red}}(f)e^{+i2\pi f d_{\mathrm{ref}}/c_0}.
$$

这里不应把完整斜距再次写入载波；横向超额路径相位已经由 PE 包络携带。粗糙面高度相位也已经包含在反射包络中。

若要继续使用通信工程常见的“正时延对应 $e^{-i2\pi f\tau}$、IFFT 得到正延迟”约定，则必须显式完成时间因子转换。最直接的实现选择是把完整声学复频响转换为共轭的信号处理频响，并同时一致处理复反射系数、复介质响应和负频率对称性。不能只给现有 `H_f` 机械乘一个符号未定义的载波。

### 5.2 为什么当前总和不是物理相干总和

当前代码在各频点执行

$$
H_f=H_{\mathrm{dir}}^{\mathrm{red}}+H_{\mathrm{ref}}^{\mathrm{red}}.
$$

但两项分别去掉了不同距离的载波。若希望保留以直达轴向距离为公共参考的声学 reduced response，至少应形成

$$
H_{\mathrm{common}}^{\mathrm{red}}
=H_{\mathrm{dir}}^{\mathrm{red}}
+H_{\mathrm{ref}}^{\mathrm{red}}
e^{+i2\pi f(d_{\mathrm{ref}}-d_{\mathrm{dir}})/c_0}.
$$

再乘直达公共载波即可得到声学物理频响。默认几何 $z_{\mathrm{tx}}=100\,\mathrm{m}$、$z_{\mathrm{rx}}=3\,\mathrm{m}$ 的解析轴向时延为

| 路径 | 距离 | 时延 |
|---|---:|---:|
| 直达 | 97 m | 64.6667 ms |
| 一次平均海面反射 | 103 m | 68.6667 ms |
| 相对时延 | 6 m | 4.0000 ms |

当前直接相加会丢掉这项 $4\,\mathrm{ms}$ 的载波相对相位，只剩各自包络中的横向/介质/边界残余相位。因此，`H_f=H_direct_f+H_reflect_f` 只能证明代码内代数闭合，不能证明物理多径闭合。

### 5.3 项目内部现有约定冲突

- `validate_pe_phase_convention_uniform_vertical.m` 和 `vertical_comm_guide.md` 的 2026-07-20 相位审计采用 $e^{-i\omega t}$，正相位斜率表示正时延，并以 FFT 提取验证 CIR；这与 WAPE 载波重构公式一致。
- `build_physical_cir_vertical.m` 的接口说明采用负相位斜率表示正时延，并使用 IFFT；这是另一套信号处理约定。
- `comm_main_vertical_psk.m` 和 `build_communication_taps_vertical.m` 直接 IFFT 输入 `H_f`，但没有把当前 PE reduced envelope 明确转换到后一约定。
- `vertical_comm_guide.md` 的早期 Bellhop小节保留过负号载波公式，紧随其后的 2026-07-20 小节已说明正确的验证层正号。文档中这两段历史公式仍容易被误读。

这些约定并非谁天然“错误”；错误在于字段没有携带明确 convention/reference metadata，且消费者跨约定使用时没有显式转换。

## 6. 网格、源与边界的适用性

### 6.1 高斯初场

初场为固定物理宽度 $\sigma=0.3\,\mathrm{m}$ 的二维高斯。连续无限域近似下，角谱能量的 90% 半径满足

$$
\kappa_{90}=\frac{\sqrt{\ln 10}}{\sigma}\approx5.06\,\mathrm{rad/m}.
$$

对应角度约为 4 kHz 时 17.6°、8 kHz 时 8.7°。这说明主能量近轴向，但该估计不是离散网格、海绵和粗糙散射后的严格谱界。

固定 $\sigma$ 使波束角随频率变窄，可解释为固定孔径换能器的简化模型；若目标是固定角宽，源宽应随波长调整。初场中心幅度固定为 1，没有源级、声阻抗、1 m 参考声压或实际换能器方向图标定，所以 `path_loss_db=-20log10(abs(h_total))` 不是可直接对实验的绝对 TL。

### 6.2 横向采样

- 公共 API 默认 `xw=yw=50 m`、`nx=ny=1024`，$\Delta x=0.04883\,\mathrm{m}$，在 8 kHz 为 3.84 个全波长网格点。
- 主通信 demo 使用 `nx=ny=256`，$\Delta x=0.19531\,\mathrm{m}$，在 8 kHz 约为 0.96 个全波长网格点。

PE 的横向 FFT 网格应按需要表示的 $\kappa$ 范围判断。通信 demo 的 $\kappa_{\mathrm{Nyq}}=\pi/\Delta x\approx16.08\,\mathrm{rad/m}$，仍覆盖初始高斯的主要角谱；但它不能覆盖全部 8 kHz 传播角，并且粗糙面可能把能量重分布到更高 $\kappa$。因此不能仅凭源波束较窄就宣布粗糙反射后的采样充分，需检查传播角窗能量、aliasing 和网格收敛。

### 6.3 海绵与窗口

现有 C0--C4 严格矩阵中，固定 32 m 窗口时加密网格和把 `stepz_lamb` 从 0.5 减到 0.25 的相位变化很小；但 C2 同时扩大到 64 m 窗口并移动了海绵物理布局，反射相位 RMS 相对 C0 为 `0.22558887 rad`，超过 `0.15 rad` 门槛。

这个结果只能说明“窗口+海绵组合敏感”，不能单独归因于窗口、海绵宽度、吸收强度或采样间距。后续必须固定其余物理量逐项扫描。

## 7. 海面模型审批

### 7.1 默认 `kirchhoff_spatial`

代码生成 PM 海面 realization，并执行

$$
\Delta\phi=k_0(\cos\theta_i+\cos\theta_r)\eta,
\qquad
\Psi_{\mathrm{ref}}=R_0e^{i\Delta\phi}\Psi_{\mathrm{inc}}.
$$

法向模式给出 $\Delta\phi=2k_0\eta$，默认压力释放系数 $R_0=-1$。其优点是：

- $H_s=0$ 时严格退化为平面压力释放反射；
- 能以低成本把具体海面 realization 引入相位与角谱展宽；
- `kirchhoff_spatial` 与 `kirchhoff_kdomain` 只是同一离散相位屏的空间域/隐式卷积实现，项目已验证到 roundoff 级等价。

其限制是：

- 不在真实起伏曲面上解非局域 Dirichlet 边界条件；
- 不含遮蔽、再入射、多次散射或完整斜入射耦合；
- 纯相位屏的全平面幅度能量保持不等于传播角窗内能量、通量或散射截面已经物理正确；
- `sea_hs_target=0.5 m` 时 $\sigma_\eta=H_s/4=0.125 m$，故 $k_0\sigma_\eta$ 在 4/8 kHz 约为 2.09/4.19，相干镜面项会很弱。此时更需要独立粗糙面参考，而不能只靠 flat limit。

### 7.2 可选统计分支

当前项目还实现：

- `kirchhoff_kstat`：近法向 Gaussian/Kirchhoff 统计相位屏；
- `ssa_stat_kernel/pm_convolution`：工程统计基线；
- `ssa_stat_kernel/ssa1_geometry`：压力释放 Dirichlet 一阶几何核；
- 可选 SSA2 coherent correction：仅修正相干反射诊断，不含完整二阶非相干散射。

这些分支已有能量闭合、平面退化、seed 可重复性、显式/统计 ensemble 和接收端统计等验证。项目还实现了 joint-frequency/cached PE、条件 $\mu/C/P$ 建模和快速抽样。应准确称其为“当前假设下的统计建模与数值验证能力”，不能据此声称已经得到完整 SSA/NLSSA、绝对散射截面或实验标定的时变海面信道。

## 8. 现有验证可以证明什么

### 8.1 已证明或有较强支持

- 自由因子代数恒等误差约 `2.2741e-13`。
- 对相同离散高斯源和无实际海绵的一步角谱参考，正确正载波的相位 RMS 最大约 `3.1021e-11 rad`，群时延差最大约 `1.2768e-12 ms`。
- 平面压力释放海面下，`H_s=0` 退化、反射符号、直达禁用回归和 `H_f=H_direct_f+H_reflect_f` 软件不变量通过。
- Bellhop 多几何矩阵中，PE 相对解析/Bellhop 的路径峰时延最大误差约 `0.027285 ms`；该比较在验证层恢复了载波，不能反推公共 `H_f` 已携带绝对时延。
- 固定 32 m 窗口下的网格加密和纵向步长减半结果稳定。
- 显式 Kirchhoff、k-stat、joint-frequency、cached receiver 和条件统计生成路径已有多项数值闭合与 held-out 统计验证；它们是分支内部和相互比较的证据。

### 8.2 尚不能推出

- 不能推出当前公共 `H_f` 可直接 IFFT 得到物理因果 CIR。
- 不能推出主通信 demo 已保留直达与一次海面反射的解析相对时延。
- 不能推出分层声速、气泡等效介质、粗糙海面或统计跨频相关已通过独立实验/高保真物理验证。
- 不能推出绝对 TL、源级或换能器方向图正确。早期单几何 Bellhop run 约有 3 dB 公共幅度差；严格多几何矩阵只允许一个公共标定后，反射 TL 残差 RMS/最大约为 `1.6295/2.0751 dB`，并未通过全部预设幅度门槛。
- 不能推出横向开放边界已经收敛；C2 仍显示窗口/海绵组合敏感。
- 不能推出一般宽角、强声速异常、强回散、多次表面往返、海底耦合或时变 Doppler 适用。
- `fd_hz_used` 当前主要是函数返回值；没有把 Doppler 一致映射到 `H_f/H_baseband`。固定 seed 的静态海面也不等于随时间演化的随机信道。

## 9. 优先整改与验证建议

### P0：物理宽带与通信解释

1. **定义唯一、机器可读的频响语义。** 至少记录时间因子、FFT/IFFT 符号、参考距离、载波是否已恢复、频率轴和负频率规则。
2. **保留兼容性，采用新增字段或显式转换函数。** 不要静默改写现有 `H_direct_f`、`H_reflect_f`、`H_f`。可新增 `H_*_common_reduced_f`、`H_*_acoustic_f`、`H_*_signal_f` 等明确字段，或提供一个强制选择 convention 的转换器。
3. **先统一分量参考，再相加。** 直达和反射必须先转换到同一公共载波参考；随后才可构造物理总频响和 CIR。
4. **统一 `build_physical_cir_vertical`、`build_communication_taps_vertical` 与主通信入口。** 每个入口应拒绝缺少 convention/reference metadata 的 PE 频响，或要求调用方显式声明转换。
5. **新增端到端两径验收。** 覆盖 on-axis/off-axis、多组深度、平面海面和多个频带；公共转换后的频响直接进入目标 FFT/IFFT 后，直达与反射峰应落在解析像源时延，并与独立角谱/Bellhop 复频响闭合。验证器不得再次偷偷补入答案相位。
6. **在 P0 完成前限制 BER/SER 的物理表述。** 现有结果可继续用于相同 reduced 链、相同同步/均衡器下的相对比较，但不能称为已验证的绝对海洋多径性能。

### P1：传播、介质与边界可信度

7. 把随推进坐标变化的屏改为中点取样或明确的对称端点离散，同时保持主直达循环与 `local_march_field` 一致。
8. 对 `stepz_lamb=[0.5,0.25,0.125,0.0625]` 做复场误差与收敛阶检查；分别验证纯相速、纯衰减和二者同时存在。
9. 对只随深度变化的屏，与 $\exp[i\int(k(z)-k_0)ds-\int\alpha(z)ds]$ 或局部角谱/WKB 参考比较；再用至少一个独立 Helmholtz/PE 实现检查横向非均匀案例。
10. 独立扫描横向物理窗口、固定物理海绵宽度、最大吸收和网格间距；记录进入海绵前的角谱/能量比例和回绕反射。
11. 对粗糙面执行 flat limit、单正弦面、Gaussian/PM ensemble、相干反射、非相干谱、传播角窗通量和互易检查；以独立 SSA/积分方程小规模结果或实验数据作为外部参考。

### P2：绝对信道与时变通信

12. 通过 Green/self-starter 或明确的换能器孔径/方向图标定源，定义 1 m 参考和绝对 TL。
13. 加入或明确忽略海水体吸收，并给出逐项 dB 损耗账本。
14. 为动态海面、气泡和平台运动定义时间相关、跨频相关和 Doppler 映射，不把独立静态 realization 当作时变信道。
15. 先建立无噪声、零误码的块传输/均衡基线，再解释高 Eb/N0 BER；当前两节点通信报告已记录约 0.2%--0.5% 的无噪声误码平台。
16. 最终使用水池或海试的到达时延、PDP、相干带宽、相干反射、散射统计和 Doppler 谱作外部验收。

## 10. 文献与代码依据

### 项目内依据

- `vertical_wape_propagator.m`：自由半步、终点屏取样、直达/反射分段推进及 `H_f` 组装。
- `vertical_channel_model.m`：公共配置、默认值与输出字段。
- `pm_surface_boundary_model.m`：`kirchhoff_spatial`、`kirchhoff_kdomain`、`kirchhoff_kstat`、`ssa_stat_kernel` 的公式与 metadata。
- `comm_main_vertical_psk.m`：当前 `H_f -> interp1 -> ifft` 通信路径。
- `build_physical_cir_vertical.m`：负相位斜率/IFFT 约定。
- `build_communication_taps_vertical.m`：通信等效 tap 与圆周能量窗。
- `scripts/validation/validate_pe_phase_convention_uniform_vertical.m`：正载波、正斜率/FFT 的声学验证约定。
- `results/validation/pe_phase_convention_uniform/pe_phase_convention_audit_report.md`：均匀介质相位审计结果。
- `results/validation/pe_bellhop_flat_surface_matrix/pe_bellhop_flat_surface_matrix_report.md`：多几何 Bellhop 与 C0--C4 严格矩阵。
- `PROJECT_CONTEXT.md`、`vertical_comm_guide.md`：当前分支状态、统计信道能力、验证范围和已知限制。

### 外部一手资料

- D. J. Thomson and N. R. Chapman, “A wide-angle split-step algorithm for the parabolic equation,” JASA 74, 1848–1854 (1983), [DOI](https://doi.org/10.1121/1.390272).
- M. D. Feit and J. A. Fleck Jr., “Light propagation in graded-index optical fibers,” Applied Optics 17, 3990–3998 (1978), [DOI](https://doi.org/10.1364/AO.17.003990).
- M. D. Collins, “A split-step Padé solution for the parabolic equation method,” JASA 93, 1736–1742 (1993), [DOI](https://doi.org/10.1121/1.406739).
- T. F. Duda, *Initial Results from a Cartesian Three-Dimensional Parabolic Equation Acoustical Propagation Code*, WHOI-2006-14, [技术报告](https://oalib-acoustics.org/website_resources/PE/CARPE3D/2006-14_TR_final.pdf). 该报告明确把模型称为单色、单向 SSF，正文写出 $U=(c-c_0)/c_0$，并指出源定义、网格和侧边吸收仍需研究。
- A. G. Voronovich, “Small-slope approximation in wave scattering by rough surfaces,” Sov. Phys. JETP 62, 65–70 (1985), [原论文](https://www.jetp.ras.ru/cgi-bin/dn/e_062_01_0065.pdf).
- E. I. Thorsos and S. L. Broschat, “An investigation of the small slope approximation for scattering from rough surfaces. Part I: Theory,” JASA 97 (1995), [DOI](https://doi.org/10.1121/1.412001).

## 11. 审计边界

本次只审批和修正文档，没有修改 MATLAB 源码、公共字段或既有验证结果，也没有把临时数值试验写成新的项目证据。公式核对以当前源码、已保存的验证报告和可追溯的一手资料为依据。

没有重新推导 `pm_surface_boundary_model.m` 中所有 SSA/k-stat 高阶统计细节，也没有重新运行昂贵的 MATLAB/Bellhop 矩阵。因此本报告对这些分支给出的是代码语义、已有验证范围和剩余风险，不宣称完成了新的独立物理认证。
