# Raw-PM / Kirchhoff-KStat 条件统计信道生成可行性审查

日期：2026-07-11

范围：代码检查、已有验证结果复核、理论可行性分析与执行方案设计。未修改核心传播、海面边界或通信逻辑。

## 执行结论

核心问题的答案是：**部分可行，且接收端统计生成器本身只需小范围新增模块；但在把结果称为“物理合理的宽带 `raw_pm + kirchhoff_kstat` 信道生成器”之前，必须先解决跨频共同随机源和 raw-PM 低波数能量覆盖两个前置问题。**

分类结论如下：

| 能力 | 判断 | 依据 |
| --- | --- | --- |
| `kirchhoff_kstat` 接入公共信道入口 | 可直接实施 | 已由 `vertical_channel_model -> vertical_wape_propagator -> pm_surface_boundary_model` 逐频调用 |
| 相干项/随机项在海面处构造 | 可直接实施 | `Psi_ref_k = Psi_coh_k + Psi_sca_k` 已实现 |
| 接收端同时输出 `H_ref,coh(f)` 与 `H_ref,sca(f)` | 需要小范围修改 | 当前公共输出只有合并后的 `H_reflect_f`；现有分量诊断只服务 `ssa_stat_kernel` |
| 独立 seed realization | 可直接实施 | `sea_seed + surface_kstat_seed_offset + frequency_index - 1` 可复现 |
| 宽带逐频运行 | 可直接实施 | `f_axis` 逐频 PE 已实现且可等间隔 |
| 物理合理的跨频随机相关 | 当前理论上不充分 | kstat 每个频点使用独立复白噪声；不能代表同一海面的跨频响应 |
| `raw_pm` 风速驱动 | 需要小范围修改后实施 | 公式与离散积分已实现，但缺少风速高度定义和 capture ratio；默认 50 m 孔径在高风速严重漏掉低波数能量 |
| `mu_H,C_H,P_H` 建模和低秩抽样 | 需要新增模块 | 当前 C4 仅做已有样本 bootstrap，不估计/生成协方差模型 |
| 当前 `h_bb` 作为物理冲激响应 | 当前理论上不充分 | 它是经插值、限带后的符号率等效 taps，不等同于直接对 4--8 kHz PE 频响作物理 IFFT |
| 方案 C 的低秩算子传播 | 需要较大重构 | 单频 PE 是线性的，但当前没有跨频联合输入协方差或显式传播算子接口 |

推荐排序：**B > A > C**。A 是保留的独立测试基准；B 是主工程路线；C 仅作为跨频物理模型成熟后的研究性加速路线。

---

## Formula Inventory

| id | location | formula | variables_and_units | implementation_notes |
| --- | --- | --- | --- | --- |
| F-001 | `pm_surface_boundary_model.m:489-518` | `E(K)=alpha/(2K^3) exp[-beta g^2/(U^4K^2)]` | `K [rad/m]`, `U [m/s]`, `g [m/s^2]` | 深水色散下由 PM 频率谱变换到径向波数谱的形式在量纲和 Jacobian 上一致 |
| F-002 | `pm_surface_boundary_model.m:503-504` | `Phi_2D(Kx,Ky)=E(K)/(2*pi*K)` | `Phi_2D [m^4]`（按项目的 `dKx dKy` 约定） | 假设二维各向同性；角向积分恢复 `E(K)` |
| F-003 | `pm_surface_boundary_model.m:270-280,337-371` | `sigma_eta^2=sum(Phi_2D) dKx dKy`; kstat 内部 `C_eta(0)=sum(W_eta)dKx dKy/(2pi)^2` | `sigma_eta [m]` | `raw_pm` 为对齐项目离散口径令 `W_eta=(2pi)^2 Phi_2D` |
| F-004 | `pm_surface_boundary_model.m:765-779` | `C_eta(rho)=IFFT(W_eta) N dKx dKy/(2pi)^2` | `C_eta [m^2]` | 周期离散横向网格 |
| F-005 | `pm_surface_boundary_model.m:776-789` | `C_deltaG=exp(alpha^2(C_eta-sigma^2))-exp(-alpha^2 sigma^2)` | `alpha=2k0 [rad/m]` | 近垂直、高斯海面、Kirchhoff 相位屏假设 |
| F-006 | `pm_surface_boundary_model.m:780-789` | `S_deltaG=FFT(C_deltaG) dx dy` | 相位屏非相干谱 | 数值负谱被裁剪到零并记录诊断 |
| F-007 | `pm_surface_boundary_model.m:798-801` | `P_sca=|R0|^2 (S_deltaG * |Psi_inc|^2) dKx dKy/(2pi)^2` | K 域功率 | 支持周期或 zero-padded 卷积 |
| F-008 | `pm_surface_boundary_model.m:803-814` | `Psi_sca=sqrt(P_sca) Z`, `Z~CN(0,I)` | 复随机反射谱 | kstat 当前各频点独立；每个 K bin 也是独立 proper complex Gaussian |
| F-009 | `pm_surface_boundary_model.m:769-775` | `R_coh=R0 exp[-0.5(2k0)^2 sigma_eta^2]` | 复反射系数 | `R0=-1` 时为压力释放近垂直相干反射 |
| F-010 | `vertical_wape_propagator.m:115-320` | `H_f=H_direct_f+H_reflect_f` | 复频响 | 每个频率分别运行直达和两段反射传播 |
| F-011 | `vertical_wape_propagator.m:1019-1049` | `f_axis=linspace(fmin,fmax,Nf)` | `f [Hz]` | 自动轴等间隔；显式 `f0` 向量只做 `unique`，不强制等间隔 |
| F-012 | `comm_main_vertical_psk.m:255-288` | `h_bb=IFFT(ifftshift(interp1(H_f)))` | 符号率等效 taps | 不是未经参考时延处理的全 4--8 kHz 物理 CIR |
| F-013 | 拟议新增 | `mu=E[H]`, `C=E[(H-mu)(H-mu)^H]`, `P=E[(H-mu)(H-mu)^T]` | `F x 1`, `F x F` | 当前代码没有实现 `C` 或 `P` |
| F-014 | 拟议新增 | `H=mu+U_r Lambda_r^(1/2) z` | `z~CN(0,I)` | 仅在 proper complex Gaussian 假设经 held-out PE 验证后采用 |

## Findings by Severity

### Critical

| id | location | formula | issue_type | confidence | source | recommended_fix | validation_test |
| --- | --- | --- | --- | --- | --- | --- | --- |
| C-001 | `pm_surface_boundary_model.m:803-814`; `vertical_wape_propagator.m:237-241` | `seed_f=sea_seed+offset+frequency_index-1` | `cross_frequency_physics_missing` | `1.00` | 当前代码；`vertical_comm_guide.md:597`；现有 LFM 报告 | 在 kstat 中引入同一 realization 跨频共享的潜在随机对象；第一候选是与同一海面 realization 对齐的联合 K/f 随机源，而不是简单独立、共享相位或未标定 AR(1) | 同一 held-out seed 比较 kdomain 与 kstat 的复频率相关矩阵、反射 PDP、LFM 包络和 matched-filter；不得只看总信道 |
| C-002 | `pm_surface_boundary_model.m:489-518`; 默认 `xw=yw=50 m` | `capture=sum(Phi)dKx dKy / [alpha U^4/(4 beta g^2)]` | `finite_grid_low_k_truncation` | `0.99` | 当前公式与本次 MATLAB 数值审计；Pierson--Moskowitz 1964 DOI | 新增解析/高分辨率参考总方差、capture ratio、K_peak/K_min 元数据；在多风速建库前扩大物理孔径或限定通过 capture 门槛的风速范围 | 对每个 U 做孔径/网格收敛；建议 `capture>=0.95` 且加倍孔径后 `Hs_implied` 变化 `<2%` |

### Major

| id | location | formula | issue_type | confidence | source | recommended_fix | validation_test |
| --- | --- | --- | --- | --- | --- | --- | --- |
| M-001 | `vertical_wape_propagator.m:259-278,341-426` | `H_ref=H_ref,coh+H_ref,sca` | `receiver_component_output_missing` | `1.00` | 当前代码 | 将已有 SSA 分量诊断推广到 `kirchhoff_kstat`，逐频传播并输出接收端相干/随机分量；保持原 `H_reflect_f` 不变 | 检查每频 `H_reflect_f-H_ref_coh_f-H_ref_sca_f` 最大绝对误差 `<=1e-10` |
| M-002 | `pm_surface_boundary_model.m:489-518`; 公共参数 `sea_wind_speed` | PM 的 `U` | `wind_definition_missing_conditions` | `0.95` | 当前代码；Pierson & Moskowitz (1964) | 明确 `U` 的测量高度、平均时长和“充分成长海况”适用条件；若实际输入是 `U10`，必须显式定义转换或改用对应谱参数化 | 固定一个已声明的风速约定，与独立 PM 实现比较峰值波数、总方差和 `Hs`；记录适用/外推标志 |
| M-003 | `results/validation/validate_kstat_vs_kdomain_lfm_channel_vertical_report.md` | 反射链宽带统计 | `existing_validation_failed` | `1.00` | 项目已有报告 | 把反射通道作为跨频模型验收主对象；修复前不把 total-channel 高相关当作通过 | 至少 32/64 频点、独立 train/test seed；反射 PDP、复相关和 LFM 指标均达到预注册阈值 |
| M-004 | `comm_main_vertical_psk.m:255-288` | `h_bb=IFFT(interpolated H)` | `delay_semantics_and_aliasing` | `0.99` | 当前代码与 DFT 采样关系 | 新增独立物理 CIR 构建函数：检查等间隔、去参考时延、可选窗、明确 FFT scaling 和 delay axis；通信 taps 由物理/等效 CIR 明确转换 | 已知两径解析频响应恢复正确相对时延；窗/零填充只改变旁瓣/插值，不改变分辨率判据 |
| M-005 | `build_surface_empirical_channel_model_vertical.m`; `sample_surface_empirical_channel_vertical.m` | bootstrap samples | `statistical_generator_not_implemented` | `1.00` | 当前代码 | 新增条件统计估计器与抽样器；保留 C4 bootstrap 作为非参数基准 | 用未参与训练的 PE 集比较均值、协方差、伪协方差、PDP、分布和 BER/SER |
| M-006 | 拟议 `L>=2F` | `rank(C_hat)<=min(F,L-1)` | `covariance_estimation_uncertainty` | `1.00` | 线性代数；Ledoit--Wolf 2004 | 原型至少从 `L_train=4F` 起步并画收敛曲线；使用收缩、Hermitian 化、负特征值裁剪和 held-out 选秩 | 对 `L=2F,4F,8F` 比较 held-out covariance error 和 eigen-spectrum 稳定性 |

### Minor

| id | location | formula | issue_type | confidence | source | recommended_fix | validation_test |
| --- | --- | --- | --- | --- | --- | --- | --- |
| m-001 | `vertical_wape_propagator.m:280-292` | 仅 `idx_f_ref` 保存 `roughness_meta` | `wideband_metadata_incomplete` | `1.00` | 当前代码 | 将紧凑 kstat 元数据按频率保存，避免保存完整二维谱 | 检查所有频点的 seed、`R_coh`、`sigma_eta`、能量闭合和 capture 元数据齐全 |
| m-002 | `vertical_channel_model.m:135-162` | 默认 `target_hs + kirchhoff_spatial` | `proposed_default_not_active` | `1.00` | 当前代码与文档 | 原型脚本显式覆盖；在物理验收前不要改变公共默认值 | 默认回归与当前结果数值一致；原型配置显式记录两项覆盖 |
| m-003 | `pm_surface_boundary_model.m:429-487` | raw PM 元数据 | `capture_metric_missing` | `1.00` | 当前代码 | 新增 `pm_variance_infinite_reference_m2`, `spectral_capture_ratio`, `K_peak`, `K_min`, `K_nyquist` | 与独立数值积分和解析积分比较，相对误差 `<1e-3`（充分网格时） |

### Advisory

| id | location | formula | issue_type | confidence | source | recommended_fix | validation_test |
| --- | --- | --- | --- | --- | --- | --- | --- |
| A-001 | 拟议统计模型 | `P_H=E[(H-mu)(H-mu)^T]` | `complex_properness_unknown_at_receiver` | `0.90` | Schreier & Scharf 2003；当前 kstat 构造 | 第一版始终估计并保存 `P_H`；只有在 held-out 检验支持 proper 时才用仅 `C_H` 的生成器 | 报告 `||P||F/||C||F`、bootstrap CI 和实/虚联合 QQ；必要时用 2F 维实增广协方差抽样 |
| A-002 | 多环境节点 | `C(theta)` | `psd_interpolation_risk` | `0.95` | PSD 锥性质 | v1 只用离散节点；以后可对同一频轴、时延对齐后的协方差作非负权重凸组合，天然保持 PSD；低秩 SPD 情况再研究 log-Euclidean | 任意插值点最小特征值 `>=-tol`，并做留一环境节点验证 |
| A-003 | kstat 相位屏 | proper Gaussian input + complex-linear PE | `higher_order_loss` | `0.95` | 当前代码；复高斯二阶理论 | kstat 分支在固定频点理论上保持 proper Gaussian，但显式 kdomain 相位屏不必高斯；保存 held-out 高阶统计并保留 bootstrap 基准 | 比较偏度、峰度、KS/energy distance、K factor、关键频点联合分布 |

## Reference Mapping

| formula_or_assumption | adopted_source | citation | applicability_conditions | notes |
| --- | --- | --- | --- | --- |
| PM fully developed sea | Pierson & Moskowitz | [JGR 69, 5181--5190 (1964)](https://doi.org/10.1029/JZ069i024p05181) | 充分成长风浪；原始数据风速约 10.29--20.58 m/s，测风高度问题在原文摘要中即被指出 | 计划中的 3--8 m/s 是模型外推，必须标注，不应默认为已验证 |
| improper complex vector needs covariance and complementary covariance | Schreier & Scharf | [IEEE TSP 51(3), 714--725 (2003)](https://doi.org/10.1109/TSP.2002.808085) | 复随机向量二阶描述 | 支持第一版估计 `P_H` 的要求 |
| high-dimensional covariance shrinkage | Ledoit & Wolf | [JMVA 88(2), 365--411 (2004)](https://doi.org/10.1016/S0047-259X(03)00096-4) | F 与 L 同量级时 | 不意味着机械套用；收缩强度由训练/验证确定 |
| nearest PSD/correlation repair | Higham | [IMA J. Numer. Anal. 22(3), 329--343 (2002)](https://doi.org/10.1093/imanum/22.3.329) | 数值误差或插值产生非 PSD 时 | 首选 Hermitian 化和特征值裁剪；复杂约束时再用最近相关矩阵算法 |
| 当前 kstat 公式与边界 | 项目代码 | `pm_surface_boundary_model.m:746-931` | 近垂直、高斯海面、相位屏、无多次散射 | 代码是本次实现状态的权威来源 |
| 当前接收端验证 | 项目验证报告 | `results/validation/validate_kstat_vs_kdomain_channel_stats_vertical_report.md`; `...lfm...report.md` | reduced grid, target_hs | 不能外推为 raw-PM 多风速生产验证 |

### raw PM 数值审计

对代码中已实现的无限波数域径向谱，解析总方差为

\[
\sigma_{\eta,\infty}^2
=\int_0^\infty E(K)\,dK
=\frac{\alpha U^4}{4\beta g^2}.
\]

本次只读 MATLAB 审计使用公共默认横向设置 `xw=yw=50 m`, `nx=ny=256`，按代码同样的二维离散积分计算：

| U (m/s) | 离散方差 (m²) | 无限域参考方差 (m²) | capture ratio | 离散 `sigma_eta` (m) | 离散 `Hs_implied` (m) | `K_peak` / `K_min` (rad/m) |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 3 | 0.0022969 | 0.0023032 | 0.9972 | 0.04793 | 0.1917 | 0.7656 / 0.1257 |
| 5 | 0.0176226 | 0.0177720 | 0.9916 | 0.13275 | 0.5310 | 0.2756 / 0.1257 |
| 8 | 0.1133169 | 0.1164703 | 0.9729 | 0.33663 | 1.3465 | 0.1077 / 0.1257 |
| 10 | 0.1756922 | 0.2843513 | 0.6179 | 0.41916 | 1.6766 | 0.06890 / 0.1257 |
| 12 | 0.2084982 | 0.5896310 | 0.3536 | 0.45662 | 1.8265 | 0.04785 / 0.1257 |
| 15 | 0.2297316 | 1.4395287 | 0.1596 | 0.47930 | 1.9172 | 0.03062 / 0.1257 |

这里的 `K_peak=sqrt(2 beta/3) g/U^2` 是径向 `E(K)` 的峰值位置。表中的 capture ratio 是“当前二维 Riemann 和 / 解析无限域方差”；它同时反映低 K 截断和粗波数采样误差，并不等于理想锐截止积分的纯支持覆盖率。U=8 时峰值已经低于首个非零网格点而方差比仍接近 1，正说明首圈粗离散可能用偶然的积分补偿掩盖欠分辨，必须做孔径与 `dK` 联合收敛。

解释：高风速时 PM 峰值落到网格最小非零波数以下，增加 `nx` 而保持 `xw=50 m` 不能解决该低波数缺失；必须扩大物理孔径。用连续径向近似估算，若希望 U=15 m/s 捕获约 95% 的无限域方差，孔径量级约需 290 m；实际二维离散网格仍需收敛验证。若保持当前约 0.195 m 的横向采样间隔，300 m 孔径会把单轴网格推到约 1536 点，PE 成本显著增加。因此第一版可先在 U=5 m/s 做统计生成器原型，但多风速生产库必须先解决“低 K 孔径”和“PE 横向计算域”耦合问题，或把 PM 统计网格与 PE 网格解耦并定义保能映射。

### 当前频率轴与时延能力

通信 demo 的声学轴为 4--8 kHz、默认 `Nf_min=32`, `Nf_max=64`。当前几何中直达长度约 97 m，镜像反射长度约 103 m，相对几何时延约 4 ms。

| F | B | `Delta f=B/(F-1)` | `Delta tau≈1/B` | `Tmax=1/Delta f` | 判断 |
| ---: | ---: | ---: | ---: | ---: | --- |
| 32 | 4 kHz | 129.03 Hz | 0.25 ms | 7.75 ms | 能分开理想两径，但无模糊窗对散射尾部偏紧；主要是演示配置 |
| 64 | 4 kHz | 63.49 Hz | 0.25 ms | 15.75 ms | 分辨率不变，时延窗更合理；仍需由 held-out PDP 决定是否足够 |

必须区分两类时域量：

1. 物理带通信道：直接对统一等间隔的声学 `H(f)` 做参考时延处理和 IFFT，delay spacing 由总带宽决定，周期由频率间隔决定。
2. 当前通信 `h_bb`：先围绕参考频率取相对频率，再插值到 `fs=symbol_rate=1000 Hz` 的 FFT 网格，因此它是 1 ms 采样的符号率等效基带 taps。`n_sym=2000` 主要增加 IFFT 点数和循环时窗，不增加原始 PE 频率信息。

IFFT 前建议：先移除公共参考传播时延 `exp(-j2*pi*f*tau_ref)`，保留相对多径时延；用明确的矩形/Hann/Tukey 频域窗做旁瓣敏感性分析；窗会改变主瓣宽度。零填充只做时延轴插值，**不能提高真实时延分辨率**。

---

## Patch Guidance

### 1. 推荐数据分解

对固定环境节点，建议保存

\[
H_m(f)=H_{\rm dir}(f)+H_{\rm ref,coh}(f)+H_{\rm ref,sca}^{(m)}(f).
\]

用户任务文本中的分隔横线应理解为分项列举，不应实现成减法。当前代码的相位符号已经包含在各复分量自身，组合仍必须是加法。

主统计对象应为随机散射分量：

\[
\mu_{\rm sca}=E[H_{\rm ref,sca}],\quad
C_{\rm sca}=\operatorname{Cov}(H_{\rm ref,sca}).
\]

固定环境下理论 `mu_sca=0`，而

\[
\mu_H=H_{\rm dir}+H_{\rm ref,coh}+\mu_{\rm sca},\qquad C_H=C_{\rm sca}.
\]

推荐同时保存总信道样本用于端到端验证，但不要直接用 total-channel 的优秀指标掩盖反射随机项错误。接收端统计最合理的存储方式是“确定性直达 + 确定性相干反射 + 随机散射统计 + 派生总信道”，而不是只保存总信道或只保存反射总和。

### 2. 最小代码改动范围

新增模块（名称为建议，不在本阶段创建）：

- `estimate_conditional_channel_stats_vertical.m`：读取 `F x L` 复样本，估计 `mu,C,P`、收缩协方差、特征谱、低秩模型和不确定度。
- `sample_conditional_channel_vertical.m`：proper 或 augmented-real improper 抽样，返回 `H_f`, 物理 `h_tau` 和 metadata。
- `build_physical_cir_vertical.m`：统一频率轴检查、参考时延、窗、IFFT scaling、delay axis 和可选零填充。
- `scripts/experiments/build_kstat_conditional_library_vertical.m`：按离散环境节点建库。
- `scripts/validation/validate_kstat_conditional_generator_vertical.m`：独立 train/test 验证。
- `scripts/validation/validate_raw_pm_spectral_capture_vertical.m`：孔径/网格/U 收敛审计。

小范围修改现有接口：

- `vertical_channel_model.m`：增加 disabled-by-default 的 kstat 接收分量诊断开关及紧凑 metadata 校验；不改已有输出字段。
- `vertical_wape_propagator.m`：复用现有 SSA 分量传播结构，覆盖 kstat，并新增 `H_reflect_coh_f`, `H_reflect_sca_f`（优先放在新的 `output.kstat_channel_meta` 下，避免公共顶层膨胀）。
- `pm_surface_boundary_model.m`：增加 kstat 跨频随机源接口与 raw PM capture metadata；保留现有 independent 模式作为工程诊断兼容项。

可复用模块：

- 公共 `vertical_channel_model(paramsV)` 入口和配置验证框架。
- `vertical_wape_propagator` 的逐频轴、直达/反射 PE、`H_f` 不变量。
- kstat 的 `C_eta -> C_deltaG -> S_deltaG -> P_sca` 公式和单频 proper Gaussian 谱生成。
- 现有 Monte Carlo 样本聚合结构、C4 bootstrap 和通信 BER/SER 链。
- `kirchhoff_kdomain` 作为同一 seed/同一显式海面跨频参考。

不应修改的稳定模块：

- WAPE 步进算子和 `local_march_field` 的物理核心。
- `H_f=H_direct_f+H_reflect_f` 与 reference-frequency 标量语义。
- `comm_main_vertical_psk` 的噪声注入位置；噪声继续在接收端单独加入。
- 现有 `kirchhoff_spatial`, `kirchhoff_kdomain`, `target_hs` 和平整海面分支。

### 3. 建议 metadata 与结果字段

条件库根字段：

- `library.kind`, `schema_version`, `created_at`, `code_revision`。
- `library.fixed_config`：几何、声速、气泡、多普勒、PE 网格、频率轴、窗和参考时延。
- `library.conditions(j).theta`：至少 `sea_wind_speed_mps`，并声明 wind convention。
- `raw_pm_meta`：`sigma_eta_m`, `Hs_implied_m`, `pm_variance_discrete_m2`, `pm_variance_infinite_reference_m2`, `spectral_capture_ratio`, `K_peak`, `K_min`, `K_nyquist`, `aperture_m`, `grid_size`。
- `samples`：`seed_train`, `seed_test`，可选紧凑 `H_direct_f`, `H_reflect_coh_f`, `H_reflect_sca_f`, `H_total_f`。
- `stats.mu_scatter_f`, `stats.mu_total_f`, `stats.C_scatter_f`, `stats.P_scatter_f`。
- `stats.covariance_estimator`, `shrinkage`, `diagonal_loading`, `eigenvalues`, `rank_retained`, `variance_retained`。
- `cir_meta.delta_f_hz`, `bandwidth_hz`, `delay_resolution_s`, `max_unambiguous_delay_s`, `reference_delay_s`, `window`, `zero_padding_factor`。
- `validation.frequency`, `validation.time`, `validation.distribution`, `validation.communication`, `validation.timing`。

### 4. 复高斯、伪协方差和分解选择

当前 kstat 在每个频率上直接生成 proper complex Gaussian K 域随机向量；固定环境下 PE 是复线性传播，因此单频接收随机散射项仍是 proper complex Gaussian。这个结论适用于“当前 kstat 数学模型”，不自动适用于显式 `exp(i2k eta)` 的 kdomain realization，也不证明真实海洋信道是高斯。

第一版仍应估计伪协方差。若 `P` 明显不为零，不能只用 `C` 与 circular `CN(0,I)`；应对 `[Re(H); Im(H)]` 的 `2F x 2F` 实协方差分解并抽样。

分解优先级：

1. EVD/SVD：首选；直接支持 PSD、特征值裁剪和低秩生成。
2. Cholesky：仅在收缩/加载后严格正定且希望全秩生成时使用；不适合作为低秩主接口。
3. 低秩阈值：同时报告累计方差 99%、99.9% 和 held-out 指标；最终 r 由 held-out covariance/PDP/BER 选择，而非只按能量。

数值处理顺序：去均值；`C=(C+C')/2`；估计收缩；EVD；小负特征值裁零；记录被裁负能量；按 held-out 选择秩。对角加载只用于数值稳定，不能用来掩盖样本量不足。

### 5. 环境参数条件化

第一版应只使用离散节点，不做连续插值。每个 `(U,z_rx,bubble config,c(z),frequency axis,PE grid)` 组合都是独立模型版本。

参数影响分类：

| 参数 | 主要影响 | 是否重建统计模型 |
| --- | --- | --- |
| `U`（raw PM） | 海面谱形状、方差、kstat 相干/散射统计 | 是 |
| `z_rx` | 反射/直达传播距离、相位、接收采样位置 | 是 |
| `c(z)` / 环境模式 | 完整 PE 传播算子 | 是 |
| 气泡参数 | 频率相关声速/衰减屏与传播算子 | 是 |
| 发射位置/深度 | 入射场与两段传播 | 是 |
| 多普勒/运动 | 频率/时间映射，超出静态库 | 必须建立新的时变模型，不与 v1 混合 |
| 噪声 | 接收机观测 | 否；不得写入 `H` 或 `h` |

以后若做插值：均值先移除参考时延后对复实/虚部插值；协方差最简单的 PSD 保持方法是同一频轴上的凸组合 `C=sum w_j C_j`, `w_j>=0`, `sum w_j=1`。若节点低秩或接近奇异，直接 matrix-log 插值并不稳健。对 improper 模型应插值增广实协方差而不是分别随意插值 `C` 和 `P`。

### 6. 方案 A/B/C 比较

| 排名 | 方案 | 结论 | 作用 |
| ---: | --- | --- | --- |
| 1 | B：少量 kstat+PE 后建接收端统计 | 推荐主线 | 对现有接口改动有限；生成成本低；能保留样本中已有的跨频协方差，但前提是输入样本的跨频随机源先物理化 |
| 2 | A：每条信道 kstat+PE | 必须保留 | 物理/实现基准和 held-out 测试集；不适合大规模生产 |
| 3 | C：传播反射场统计 | 后续研究 | 单频线性算子允许低秩输入模态逐个传播；但当前缺少跨频联合输入协方差、算子接口，完整海面协方差维度过高 |

方案 C 的可行低秩版本不是显式构造 `A C_Psi A^H`，而是：先求海面随机场/随机反射谱的主模态 `q_r`，对每个模态调用现有 PE 得到接收端响应 `A q_r`，再累加外积。它只有在跨频潜变量定义明确、输入模态秩远小于网格维数时才有收益；不应作为 v1 默认。

### 7. 风险清单

| 风险 | 当前等级 | 后果 | 控制措施 |
| --- | --- | --- | --- |
| raw PM 归一化/风速约定 | 高 | U 标签与真实海况含义不清 | 声明 wind convention；独立公式对照；保留 `target_hs` 控制模式 |
| 有限网格谱能量截断 | 极高（U>=10） | `Hs_implied` 人为饱和 | capture ratio 门槛；孔径收敛；必要时谱网格/PE 网格解耦 |
| 跨频独立随机源 | 极高 | H(f) 不连续、虚假长时延/PDP | 共同潜变量；以 kdomain 同一海面为参考；反射链验收 |
| 样本协方差秩/噪声 | 高 | 虚假特征模态、生成失真 | `L>=4F` 起步、收缩、held-out 选秩和收敛曲线 |
| 高阶非高斯损失 | 中 | 尾部衰落、K factor、BER outage 偏差 | 与 bootstrap/kdomain 比较；必要时混合/椭圆/非参数残差模型 |
| IFFT 时延混叠 | 高 | 路径折返、RMS delay 偏差 | 调小 `Delta f`、参考时延、明确窗和 delay axis |
| PE 计算成本 | 高 | 建库不可承受 | 先单 U/128²；计时；缓存确定性分量；以后研究方案 C |
| 环境参数插值 | 中 | 非 PSD、相位抵消 | v1 离散节点；延迟对齐；PSD 凸组合；留一节点验证 |
| 气泡接入 | 高 | 传播算子改变，旧库失效 | 气泡配置纳入条件键和 schema hash；重新建库 |
| 默认切换过早 | 高 | 破坏历史结果可比性 | 原型显式配置；通过 capture/跨频/held-out 后再讨论默认值 |

---

## Validation Checklist

### 分阶段执行方案

#### 阶段 0：现状确认（本报告已完成代码部分）

- 已确认公共入口、默认配置、seed、频率轴和输出字段。
- 已确认 kstat 公式链完整到随机反射场，并进入 PE。
- 已确认 public receiver output 尚未分离 kstat coherent/scatter。
- 已确认当前跨频独立。
- 已确认 raw PM 默认 50 m 孔径在高 U 严重漏谱。
- 尚需用新增脚本把 capture 审计固化为可回归结果。

#### 阶段 1：单一风速条件原型

- 固定 `U=5 m/s`, `raw_pm`, `kirchhoff_kstat`, bubbles off, Doppler off。
- 先修复/增加跨频共同随机源；保留 `independent` 为诊断对照。
- 接收端输出 `H_dir`, `H_ref,coh`, `H_ref,sca`, `H_total`。
- 建立 `mu,C,P` 和 proper/augmented-real 两种抽样路径。
- 使用独立 train/test PE 样本。

#### 阶段 2：多风速条件统计库

- 在通过 capture 门槛的网格上建立 `U={3,5,8,10,12,15}` 离散节点。
- 不插值；每节点独立版本化。
- 高 U 若无法满足 capture，不得静默建库，应标为 invalid 或缩小风速范围。

#### 阶段 3：大规模数据生成

- 批量产生 `H_f`, `h_tau`，标签含 `U`, `Hs_implied`, capture ratio 和 model version。
- 生成 100/1000/10000 条计时、内存和文件大小报告。

#### 阶段 4：通信性能验证

- 相同符号、调制、均衡器和噪声 seed，对比 held-out 直接 PE 与统计抽样。
- 比较 BER/SER 曲线与风速趋势，而非逐 realization 相等。

#### 阶段 5：扩展环境参数

- 先加入离散 `z_rx`；再加入气泡；最后单独设计多普勒/时变模型。
- 每个会改变传播算子的参数均触发重建，不能只插值旧 U 库。

### 验证指标与建议阈值

用户给出的均值/协方差公式排版中 `|.|*2`、`|.|*F` 应修正为二范数和 Frobenius 范数：

\[
\epsilon_\mu=\frac{\|\mu_{gen}-\mu_{PE}\|_2}{\|\mu_{PE}\|_2},\qquad
\epsilon_C=\frac{\|C_{gen}-C_{PE}\|_F}{\|C_{PE}\|_F}.
\]

建议 v1 预注册阈值（原型阈值，不代表最终物理标定）：

| test_id | setup | metric | pass_threshold | result |
| --- | --- | --- | --- | --- |
| V-001 | raw PM U=5，孔径/网格加倍 | capture 与 `Hs_implied` | capture `>=0.95`; Hs 变化 `<2%` | 当前 50 m capture=0.9916；加倍收敛待测 |
| V-002 | flat surface | `H_ref_sca`, component sum | scatter `<=1e-12`; sum error `<=1e-10` | 待新增分量输出 |
| V-003 | same/different seed | repeatability/independence | same seed bitwise/roundoff equal；different seed scatter changes，direct/coh fixed | 单频已有验证，宽带共同源待测 |
| V-004 | train/test | `epsilon_mu`, `epsilon_C` | 原型目标 `<0.10`, `<0.20`，并报告 CI | 待测 |
| V-005 | train/test | `P` properness | `||P||F/||C||F` 与 bootstrap CI；不预设为零 | 待测 |
| V-006 | 反射频响 | 邻频复相关、相位连续性 | kstat 与 kdomain/held-out 误差小于预注册界限 | 当前未通过 |
| V-007 | 物理 CIR | 主径/反射径、PDP、mean/RMS delay | 几何相对时延误差 `<1` 个真实 delay bin；PDP correlation `>=0.9` | 当前 reflected LFM 0.276，未通过 |
| V-008 | 分布 | KS/QQ/偏度/峰度/K factor | held-out 与生成的差异有 CI；阈值由 pilot 预注册 | 待测 |
| V-009 | 通信 | BER/SER vs Eb/N0 | 曲线 CI 重叠或预定义最大绝对差；风速排序一致 | 待测 |
| V-010 | 计时 | A vs B | `M_break-even=T_build/(T_PE-T_sample)` | 无现有可靠计时，必须实测 |

### 具体下一阶段 Codex 任务

建议下一任务标题：**“实现 U=5 m/s 的 raw-PM/KStat 条件统计信道原型及跨频共同随机源验证”**。

第一批功能：

1. 新增 raw PM capture 审计及 metadata，不改变默认模式。
2. 为 kstat 增加 disabled-by-default 的接收端 coherent/scatter 分量输出。
3. 设计一个物理可解释的跨频共同潜变量版本；保留当前 independent 模式做对照，且不得称 AR(1)/shared phase 为物理标定模型。
4. 新增 `mu,C,P` 估计、收缩 EVD、proper/augmented-real 抽样。
5. 新增独立物理 CIR 构建函数和 train/test 验证脚本。

推荐原型参数：

- `U=5 m/s`, `raw_pm`, `kirchhoff_kstat`。
- `f_band=[4000,8000] Hz`, 软件原型先 `F=32`；物理 PDP 验收改为 `F=64`，若主要能量超过 15 ms 则进一步令 `Delta f<=50 Hz`（4 kHz 带宽至少约 81 点）。
- `nx=ny=128`, `xw=yw=50 m`, `show_figures=false`, bubbles/doppler off。
- `L_train=128 (=4F)`，`L_test=64`；先以 `L=32,64` 做流水线 smoke，再运行完整原型。
- seed 集严格不重叠，并保存列表。

预计计算量：

- F=32 完整原型约 `(128+64+1)*32 = 6176` 个“频点级 channel solve”等价工作量；每个反射频点内部包含发射到海面与海面到接收端的传播，且当前还重复计算直达项。
- 六风速按同样配置约 37,056 个频点级 solve，尚未计入孔径扩大带来的单次成本上升，因此不应在阶段 1 直接全跑。
- 统计抽样本身约为 `O(Fr)` 每条，预计远小于 PE；真实 break-even 必须通过 `tic/toc` 实测，现有产物没有可靠计时字段。

验收标准：

- raw PM capture 通过且元数据完整。
- component sum、`H_f` 不变量和 disabled-path 回归通过。
- reflected-only 跨频相关/PDP/LFM 指标显著优于当前 independent 基线。
- held-out `epsilon_mu<0.10`, `epsilon_C<0.20` 作为首轮工程目标，并报告置信区间。
- properness、非高斯分布和 BER/SER 均有 held-out 结果；未通过时自动回退到增广协方差或 bootstrap 基准，不夸大结论。

是否先调整频率轴：**是，但分两步。** F=32 可用于接口和统计流水线原型；在声称得到有意义的物理 `h(tau)` 前至少使用 F=64，并根据 held-out 反射 PDP 将 `Tmax` 扩展到主要能量窗之外。单纯把通信 IFFT 零填充到更长长度不算调整频率轴。

## 最终判断

当前项目已经具备方案 B 所需的大部分“上游”和“下游”骨架：raw PM 公式、kstat 单频统计散射、可复现随机 realization、完整 PE、稳定的 `H_direct_f/H_reflect_f/H_f/f_axis` 输出，以及通信消费链。因此，**接收端 `mu_H(U),C_H(U)` 生成器的软件主体可用较小改动实现**。

但当前最主要阻碍按优先级为：

1. **跨频相关**：当前逐频独立，直接阻断物理宽带 `h(tau)` 结论。
2. **海面统计网格**：默认 50 m 孔径对 U=10--15 m/s raw PM 低波数能量捕获严重不足。
3. **接收端分量接口**：内部有 coherent/scatter，接收端公共输出尚未分离；这是小改动。
4. **IFFT 时延表示**：当前通信 taps 是符号率等效信道，需新增物理 CIR 定义。
5. **统计估计**：`C/P`、收缩、低秩与 held-out 验证尚未实现，但技术上直接。
6. **计算成本**：可通过 B 摊销，但建库本身仍昂贵，需先单 U 计时与收敛。

因此，不应立即把公共默认从 `target_hs + kirchhoff_spatial` 改成 `raw_pm + kirchhoff_kstat`。正确顺序是先完成 U=5 的跨频/接收端统计原型，再解决高风速 PM 孔径覆盖，最后才讨论 raw PM 作为正式默认模式。
