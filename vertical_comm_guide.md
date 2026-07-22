# 垂直水声 PE、海面统计与通信接口指南

本文是当前实现的技术参考，面向需要审阅物理公式、代码数据流和结果语义的研究与开发人员。运行命令集中在 `scripts/README.md`，项目演进记录集中在 `PROJECT_CONTEXT.md`，具体实验数字和图表集中在 `reports/` 与 `results/`。

## 1. 当前核心方案

公共入口为：

```matlab
output = vertical_channel_model(paramsV);
```

核心数据流是：

1. 在发射深度生成二维高斯复包络；
2. PE/WAPE 向上传播到接收深度，得到直达 reduced envelope；
3. 另一路传播到海面，应用海面边界，再传播到接收深度；
4. 在接收端统一两条路径的载波相位参考；
5. 形成宽带 (H(f))，供 CIR、LFM 或通信链消费。

公共默认海面分支仍为 `kirchhoff_spatial`。cached PE、精确伴随投影、解析 FFT `C/P` 和条件统计生成器是研究/加速路径，不是新的海面物理模型。

## 2. 坐标、时间与载波约定

海面为 (z_s=0)，深度向下为正，且

\[
0\le z_{\rm rx}<z_{\rm tx}.
\]

物理复相量采用

\[
p(t)=\Re\{P(f)e^{-i2\pi ft}\}.
\]

正传播时延在该物理相量约定下对应正频率相位斜率。MATLAB `ifft` 的 DSP 频响则采用相反符号：正相对时延对应

\[
e^{-i2\pi f\tau}.
\]

这两个符号分别服务于物理相量和 DSP 合成，不能混用。

### 2.1 PE reduced envelope

当前 split-step PE 的自由传播因子等价于

\[
e^{id(k_z-k_0)},
\]

因此 PE 数组保存的是去掉名义纵向载波 ($e^{ik_0d}$) 后的复包络。直达路径和海面反射路径具有不同的纵向跨度，二者的 raw PE 结果不能在没有统一参考的情况下直接解释为物理总信道。

定义

\[
\tau_{\rm dir,0}=\frac{z_{\rm tx}-z_{\rm rx}}{c_0},
\qquad
\tau_{\rm ref,0}=\frac{z_{\rm tx}+z_{\rm rx}-2z_s}{c_0},
\]

\[
\Delta\tau_0=\tau_{\rm ref,0}-\tau_{\rm dir,0}.
\]

公共默认使用直达参考的 DSP 表示：

\[
H_{\rm dir}^{\rm dsp}=H_{\rm dir}^{\rm red},
\qquad
H_{\rm ref}^{\rm dsp}
=e^{-i2\pi f\Delta\tau_0}H_{\rm ref}^{\rm red}.
\]

绝对物理相量为

\[
H_{\rm dir}^{\rm phys}
=e^{+i2\pi f\tau_{\rm dir,0}}H_{\rm dir}^{\rm red},
\]

\[
H_{\rm ref}^{\rm phys}
=e^{+i2\pi f\tau_{\rm ref,0}}H_{\rm ref}^{\rm red}.
\]

标准几何 (z_{\rm tx}=100\rm,m)、(z_{\rm rx}=3\rm,m)、(c_0=1500\rm,m/s) 对应 97 m 直达跨度、103 m 海面反射跨度和 4 ms 名义相对参考延迟。

4 ms 只恢复被两个 reduced PE 路径分别去除的确定性载波参考。PE 衍射、海面随机相位和有限带宽仍会改变实际 CIR/PDP 的峰值、展宽和群时延。程序不对每条 realization 的直达和反射峰分别对齐。

## 3. 公共相位接口与输出

输入项：

```matlab
paramsV.channel_phase_reference = 'direct_dsp';   % default
paramsV.channel_phase_reference = 'legacy_reduced';
```

默认公共字段为 direct-DSP：

- `H_direct_f`、`H_reflect_f`、`H_f`；
- `h_direct`、`h_reflect`、`h_total`；
- `H_f = H_direct_f + H_reflect_f`；
- `h_total = H_f(idx_f_ref)`。

附加的明确表示为：

- `H_direct_reduced_f/H_reflect_reduced_f/H_total_reduced_f`；
- `H_direct_physical_f/H_reflect_physical_f/H_physical_f`；
- 对应参考频率标量；
- `phase_reference_meta`。

`legacy_reduced` 用于回归旧结果。`H_total_reduced_f` 是历史代数组合，不应再称为相位已经统一的物理总信道。

## 4. PE/WAPE 推进

发射平面采用二维高斯复包络

\[
\Psi_{\rm tx}(x,y)=
\exp\left[-\frac{(x-x_{\rm tx})^2+(y-y_{\rm tx})^2}
{2\sigma_{\rm src}^2}\right].
\]

一步 cached uniform 核心为：

```matlab
psi_k = fr .* fft2(screen .* ifft2(fr .* psi_k));
```

`fr` 包含横向衍射的半步因子，`screen` 包含介质相位和 sponge 衰减。公共传播器还支持现有 layered/bubble/Doppler 配置；精确伴随原型仅覆盖已验证的 uniform CPU-double 固定路径。

## 5. 海面边界模型

### 5.1 显式 Kirchhoff

`kirchhoff_spatial` 和 `kirchhoff_kdomain` 对具体海面 realization
\(\eta(x,y)\) 应用近垂直相位屏：

\[
G_i(x,y)=R_{0,i}e^{i\alpha_i\eta(x,y)},
\qquad
\alpha_i=\frac{4\pi f_i}{c_0}=2k_{0,i}.
\]

压力释放面通常取 (R_0=-1)。`kirchhoff_kdomain` 是同一显式相位屏的波数域接口，不是另一套独立海面理论。

### 5.2 joint-frequency Kirchhoff K-stat

统计分支不需要显式海面 realization。其随机残差严格定义为

\[
\delta G_i=R_{0,i}e^{i\alpha_i\eta}-R_{{\rm coh},i},
\]

\[
R_{{\rm coh},i}=R_{0,i}
\exp\left(-\frac12\alpha_i^2\sigma_\eta^2\right).
\]

当前 $C_{\delta G,ij}$ 与 $P_{\delta G,ij}$ 已包含 $R_0$，接收权重中不得再次乘  $R_0$。joint model 保留跨频 covariance 和 pseudo-covariance；独立频点随机相位不是它的等价替代。

### 5.3 SSA research branch

`ssa_stat_kernel` 提供 `pm_convolution` 工程基线与 `ssa1_geometry` 一阶 Dirichlet 几何核。它们用于物理趋势和模型敏感性研究，不代表已经实现完整 SSA2、多次散射或实验标定海面散射。

## 6. cached forward 与精确离散伴随

固定 uniform surface-to-receiver 算子记为 (A_i)。接收点单位源为 (r)，则

\[
q_i=A_i^Hr.
\]

离散伴随反向遍历传播步，使用 `conj(fr)`、`conj(screen)`，sponge 保持相同衰减幅度，不取倒数，也不经验性补偿 FFT 尺度。

定义

\[
a_i=\operatorname{conj}(\psi_{{\rm inc},i})\odot q_i,
\]

则 reduced scatter 接收值为

\[
H_{{\rm sca},i}^{\rm red}=a_i^H\delta G_i.
\]

`q` 是接收灵敏度核，不是真实逆传播或反向声压场。相位参考转换在得到接收标量后统一施加，因此不改变精确伴随内积关系。

当前 v1 限制为：uniform sound speed、CPU double、固定 Tx、单个最近网格点 Rx、固定 PE/PM 网格与频率轴、无 bubbles/Doppler。

## 7. PM→PE 映射与解析接收统计

PM 大网格到 PE 小网格使用已有中央裁剪 (E)。解析收缩先把 PE 权重零嵌入 PM 网格：

\[
\widetilde a_i=E^Ha_i.
\]

不得在 PE 小网格上直接假设循环平稳协方差。PM 周期网格上的 reduced 统计为

\[
C_H^{\rm red}(i,j)=\widetilde a_i^H
C_{\delta G,ij}^{\rm PM}\widetilde a_j,
\]

\[
P_H^{\rm red}(i,j)=\widetilde a_i^H
P_{\delta G,ij}^{\rm PM}\widetilde a_j^*.
\]

令

\[
D=\operatorname{diag}
\left(e^{-i2\pi f\Delta\tau_0}\right),
\]

则 direct-DSP 接收统计为

\[
\mu_{\rm dsp}=D\mu_{\rm red},\qquad
C_{\rm dsp}=DC_{\rm red}D^H,
\qquad
P_{\rm dsp}=DP_{\rm red}D^T.
\]

空间 dense/FFT 收缩仍在 reduced 层完成；确定性频率旋转不改变特征值、秩、每频功率或 properness 比率。

## 8. 条件统计生成器

schema 2.x 条件模型保存 direct-DSP 的确定性分量、均值、(C/P)、复 EVD 和增广实数 EVD，以及 `phase_reference_meta`。

旧项目 schema 1.x 被视为 `legacy_reduced`。迁移时除了旋转均值和 (C/P)，还要执行

\[
U_{\rm dsp}=DU_{\rm red},
\]

以及增广实数变换

\[
T=
\begin{bmatrix}
\Re D&-\Im D\\
\Im D&\Re D
\end{bmatrix}.
\]

迁移不会从旧的 signed `reference_delay_s` 猜测几何。已知项目模型从模型/library/cache metadata 或已登记的 U=5/U=8 固定验证几何取得参数；无法确认时明确报错。

条件库只允许已经验证的离散风速节点，不做静默最近邻或插值。

旧 MAT 数据按“先审计、后归档、再重建”处理。相位中立的 PE cache、PM spectrum、中央裁剪映射和 joint factor 只有在频率、网格、几何、uniform CPU-double 约束通过核验后才可复用。分别保存直达/反射分量或完整统计且几何可靠的数据可生成 schema-2 迁移副本；只有旧总信道或几何不足的数据必须重建。迁移副本只服务闭合审计，不替代本次重新训练的规范模型。归档目录为 `results/archive/pe_phase_reference_pre_rc/<timestamp>/`，采用移动且保留原相对路径，不删除、不覆盖。所有正式 MAT 必须包含 `schema_version`、`phase_reference_meta` 和带 `run_id`、代码 revision/fingerprint、几何、频率轴、seed 定义的 `validation_run_meta`。

## 9. CIR、LFM 与通信链

`build_channel_cir_vertical` 要求显式声明输入表示：

- `direct_dsp`：使用 MATLAB `ifft`；
- `absolute_physical`：在 $e^{-i\omega t}$ 约定下使用 `fft/N` 展示物理时延；
- `legacy_reduced`：拒绝把未统一的多路径总和直接转换为 CIR。

`build_physical_cir_vertical` 仅作为旧接口包装器保留，其 `reference_delay_s` 现在只表示施加到完整信道的共同时间原点移动。

未知外部 (H(f)) 若没有项目 `phase_reference_meta`，通信辅助函数将其视为已经 DSP-ready，不静默旋转；若 metadata 明确为 `legacy_reduced` 或 `absolute_physical`，通信 IFFT 路径会拒绝输入。

LFM 是在得到 (H(f)) 后施加的线性探针，不是 PE 空间源。同步或接收窗口可以整体移动总信道，但不能分别移动直达与反射分量。

## 10. 频率网格与时延解释

等间隔频率轴的无模糊时延为

\[
T_{\rm amb}=\frac1{\Delta f},
\]

物理分辨率约为 (1/B)。零填充只细化绘图采样，不提高 (1/B) 分辨率。

F=9 的 `4:0.5:8 kHz` 有 (Delta f=500\rm,Hz)，无模糊时延只有 2 ms。标准 4 ms 延迟在这些频点上恰好绕回单位相位，因此 F=9 不能用来验收该问题。

专项相位测试使用 F=65、`4:0.0625:8 kHz`，无模糊时延为 16 ms，并包含准确 6 kHz。

## 11. 验证解释

主要验收层次为：

1. reduced PE 离散伴随内积；
2. cached forward 与 adjoint projection 的同 realization 一致性；
3. PM 零嵌入后的 dense/FFT (C/P) 一致性；
4. direct-DSP 相位旋转的解析/样本一致性；
5. 公共四种海面分支、direct-only、宽带通信回归；
6. F=65 解析两径的 4 ms CIR 验收。

explicit Kirchhoff 与 joint-kstat 的独立 realization 不应做像素点对点验收；应比较 ensemble 均值、(C/P)、功率、PDP、LFM 和分布。所有主要统计指标应优先使用 reflected scatter，避免直达波掩盖反射差异。

正式发布候选 `phase_rc_20260722_174945` 已通过全部门限。独立相位审计的算子误差为 (2.2741\times10^{-13})，正载波相位 RMS 为 (3.1007\times10^{-11}\,\mathrm{rad})，群时延差为 (1.2768\times10^{-12}\,\mathrm{ms})。完整伴随测试的离散内积误差为 (7.7684\times10^{-15})，接收投影误差为 (5.0950\times10^{-15})，dense/FFT 的 (C/P) 误差为 (5.7693\times10^{-16}) 和 (1.0052\times10^{-15})。F=9 与 F=64 的解析—样本 covariance 误差相对 split-sample floor 分别为 0.593078 和 0.850502，均小于 1.25。F=64 public/cached double consistency 为 (9.88957\times10^{-16})。正式判定和全部 gate 见 `reports/pe_phase_reference_release_candidate_report.md`；早期 reduced-grid 报告仅保留为迁移过程证据。

正式 adjoint F=64 验证使用 PE (128^2)/PM (256^2)。U=8 条件节点沿用其独立 aperture audit 通过的 PE (128^2)/PM (384^2)，不能把两者误写为同一 PM 尺寸。正式图册由同一 `run_id` 的新结果生成：空间声场仍表示 reduced complex envelope；频率响应同时显示 reduced、direct-DSP 与理论 (-2\pi f\Delta\tau_0)；F=65 CIR 用于显示 4 ms，F=9 只作为传播 smoke；6 kHz 下 4 ms 等于 24 个载波周期，因此单频空间图不能发现旧相位问题。

## 12. 当前限制

- 公共 PE 仍不是全波、全角度或多次海面散射求解器。
- joint-kstat 是近垂直 Gaussian/Kirchhoff phase-screen 统计模型。
- 精确伴随 v1 不覆盖 layered、bubble、Doppler、GPU、多接收机或插值采样。
- absolute-physical 相量用于物理审计和展示；通信默认消费 direct-DSP。
- Bellhop 严格矩阵仍保留其源归一化和横向窗口敏感性未通过项；载波相位修复不自动解决这些独立问题。
- 任何旧 MAT 若缺少可靠几何 metadata，都不得依据模糊 signed delay 静默迁移。

## 13. Li et al. (2009) 独立显式海面验证

`li2009_explicit_surface_validation.m` 和
`scripts/validation/validate_li2009_explicit_surface_vertical.m` 是独立的纯声学验证路径，
不接入通信主线或统计信道生成器。它采用论文式 (12)--(13) 将 U10 转为 U19.5，
复用 raw-PM 海面约定，并保证同一 6 ms CW pulse 的所有频率显式共享一次生成的
`eta(x,y)`。边界只有

\[
\Psi_{\rm ref}=-\Psi_{\rm inc}\exp(i2k\eta).
\]

海底同址收发几何不满足公共 API 的 `z_rx<z_tx` 约束，因此该验证在独立模块中复用
同一 uniform split-step 核，并以精确离散伴随构造 surface-to-bottom 接收投影。
该投影是计算加速，不是 kstat 或新的统计物理模型。

128 样本缩减案例中，U10=5→10 m/s 的 20% 首阈值标准差由
0.321 ms 增至 0.724 ms，shifted-Rayleigh `b` 由 0.553 ms 增至
1.611 ms；峰值宽度只增加 0.4%，未达到 5% 趋势门槛。推进步长、频点、
sponge 和窗口对照稳定，但横向网格尚未收敛，因此当前只能声明首阈值展宽趋势的
部分复现，不能声明式 (23)--(25) 或论文完整 PIES 系统的定量复现。完整参数审计、
结果和限制见 `reports/li2009_explicit_surface_validation_report.md`。
