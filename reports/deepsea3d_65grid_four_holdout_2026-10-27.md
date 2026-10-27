# 三维固定 compact-dual：65 点预测网格的四个独立频点验收

2026-10-27。状态：**`COMPACT_DUAL_65GRID_HOLDOUT_SCREENING_PASS`**，仅指固定点源、固定 dual 海面与五锚点预测器在**四个新增离散频点**的条件性筛查。连续 4–8 kHz 复响应、无混叠 CIR、BER/PER、随机 PM 海况均未认证。唯一[三维总计划](deep_sea_channel_research_plan.md)、[网格设计](deepsea_compact_dual_broadband_grid_design.md)、[原频率门](deepsea3d_g1_frequency_protocol.md)约束本次资格；原九点、二维 S4 和 Q1 BIE 失败证据未改。

## 冻结与执行顺序

五个已合格物理锚点 `4000/5000/6000/7000/8000 Hz`；`η*=0.06251088653876558 m`、单位三维点 Green、Tx=(0,150.5,1500)m、Rx=(0,0,5)m、c=1500m/s、压力释放、C2 plateau8/support12、native `exp(−iωt)` 不变。固定 `Q=R/B*`，只在相邻 1 kHz 锚点间复线性插值，`Rpred=B*Qpred`，`Δpred=Rpred−R0`，`Hpred=D+Rpred`。常规点 4750/7250 Hz 与相消压力点 4125/6562.5 Hz 在同一[预测锁](../results/validation/deepsea3d_65grid_four_20261027/prediction_lock.json)和[原生预测 CSV](../results/validation/deepsea3d_65grid_four_20261027/predictions_pre_FE.csv)中冻结，后者 SHA-256=`badfe577751812388fc519c6872453f1de4dd7b2150799176bb4f0756d368324`；锁 SHA-256=`81f03e3973830c9a02d5441515c299c96c5424268399cb975465322eedcbe8d3`。锁生成 UTC 03:10:28，独立[锁审计](../results/validation/deepsea3d_65grid_four_20261027/prediction_lock_audit.json)于 03:10:52 PASS；首个新 MATLAB 进程 UTC 03:14:14 因目录设置失败、**首个新 FE** UTC 03:17:10。没有看到新 FE 后改点、改模型或旋相。

冻结的逐点预测门：`e_R=|Rpred−RFE|/|RFE|<1e−3`（原门）、非零 `e_Δ=|Δpred−ΔFE|/|ΔFE|<1e−3` 和 `a_D=|Hpred−HFE|/|D|<1e−3`（本轮事前新增设计门）。`e_H=|Hpred−HFE|/|HFE|` 与幅相误差逐点另列；不删近零总场，不用平均值替代四个判定。新频点运行顺序严格为 4750→7250→4125→6562.5 Hz。

## 四个条件性物理参考

Ubuntu MATLAB R2023b 在 `/home/misc/cwx/deepsea3d_65grid_four_20261027/` 的独立目录执行。FE 是 Dirichlet 场展开，并非有限元。四点最终 `L=32 m/h=0.0625 m`，FE 真实边界、离节点、阶数、粗细网格/域、M32/M40 真实边界和有限 Fourier 尾界、出射/forcing/回放与解析平面控制均过旧门；参考 ready 后才各运行 FullFK 双项及原五配置。参考与 FK 数值主门分别严格 `<1e−4`、`<1e−3`，FullFK–FE 反射/非零增量各 `<1%`。未补算、更改或放宽 7.5 kHz 旧门。完整变体在[原生 MAT](../results/validation/deepsea3d_65grid_four_20261027/65grid_four.mat)及每频 `reference_runs.csv`、`reference_guards.csv`、`fk_runs.csv`、`fk_guards.csv` 中。

| 频率/Hz | 组 | FE N | 真实边界 / 半网格 | 阶数 N−4 / 网格 / 域 | FK 粗 Gauss 增量差 | FK–FE 反射 / 增量复误差 | 资格 |
|---:|---|---:|---:|---:|---:|---:|---|
| 4750 | 常规 | 16 | 2.461e−8 / 3.085e−7 | 1.771e−5 / 1.632e−7 / 2.484e−7 | 1.328e−8 | 3.684e−4 / 1.871e−4 | 条件通过 |
| 7250 | 常规 | 20 | 6.983e−8 / 3.527e−7 | 1.637e−5 / 2.859e−6 / 3.197e−8 | 5.225e−4 | 2.343e−4 / 1.206e−4 | 条件通过 |
| 4125 | 相消压力 | 16 | 2.458e−9 / 2.912e−7 | 3.332e−6 / 2.997e−7 / 4.073e−7 | 5.062e−9 | 3.989e−4 / 2.171e−4 | 条件通过 |
| 6562.5 | 相消压力 | 20 | 8.895e−9 / 3.315e−7 | 3.108e−6 / 5.301e−7 / 4.043e−7 | 1.821e−5 | 2.466e−4 / 1.209e−4 | 条件通过 |

7250 Hz 的粗 Gauss **增量**差 `5.225e−4` 虽过原门，数值裕量小于其他三点；同一变体的反射相对变化 `1.015e−3` 是另列诊断，不能偷换成增量门的通过证据，后续外推须保留这一风险。N12 仅是部分频点的预检层，7250/6562.5 Hz 的早期 forcing 分别约 `1.04e−6/2.94e−7`，未冒充最终参考；最终入选参考的 forcing、真实边界和其余资格门通过。所有 FE 敏感性是保存的粗/细/域 `N,N−2,N−4` 部分和相对最终参考的最大**经验绝对复差**，不是严格连续 Helmholtz 误差界。

## 四点原生复数与锁定预测

以下 `R_FE,H_FE` 均为原生复声压，表中实部和虚部乘 `10⁵ m⁻¹`；完整 `D/R0/R/Δ/H` 真值和 `Rpred/Δpred/Hpred` 逐分量实虚部及幅相误差在[四点盲测表](../results/validation/deepsea3d_65grid_four_20261027/four_point_blind_metrics.csv)。`ρ=|Rpred−RFE|/ε_FE`，`C=(|D|+|RFE|)/|HFE|`。

| f/Hz | `R_FE` ×10⁵ | `H_FE` ×10⁵ | `e_R` | `e_Δ` | `e_H` | `a_D` | ρ | C |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 4750 | −1.449605−4.704948i | +2.935299−1.734812i | 5.932e−5 | 3.013e−5 | 8.566e−5 | 5.515e−5 | **34.35** | 2.997 |
| 7250 | −2.509641+4.235407i | −5.759703+8.417056i | 2.753e−5 | 1.418e−5 | 1.329e−5 | 2.559e−5 | 3.94 | 1.002 |
| 4125 | −4.733728−1.352174i | +0.470735−0.371023i | 3.232e−5 | 1.759e−5 | **2.654e−4** | 3.004e−5 | **42.21** | **17.05** |
| 6562.5 | +2.607980+4.175447i | +0.460241−0.665659i | 1.312e−5 | 6.432e−6 | 7.979e−5 | 1.219e−5 | **13.20** | **12.63** |

常规两点、相消两点各自 **2/2 通过**全部三个冻结预测门；四点最大 `e_R/e_Δ/a_D=5.932e−5/3.013e−5/5.515e−5`。相消区的 `e_H` 明显受较小 `|HFE|` 放大，最大在 4125 Hz 为 `2.654e−4`，故并列报告其直达标度 `a_D=3.004e−5`，没有剔除该点。全四点最大包裹反射/总场相位差分别约 `4.809e−5/1.924e−4 rad`；幅度误差逐点留于盲测表。**4750、4125、6562.5 Hz 的预测差高于五倍 FE 经验敏感性**，因而误差在这些保存变体之上可分辨，不能描述为 FE 数值噪声以下；7250 Hz 为 `3.94 ε_FE`，不作更细精度或严格优劣判断。

解析平面/点 Green 与 native→DSP 共同直达时间原点保持旧合同，四点平面直反几何差均为 `6.633361620 ms`，平面/相位控制最大误差 `1.82e−12`。预测与参考未作逐频相位或时间配准。粗糙反射仅四个稀疏新点，不可靠地展开其群时延分支；没有通过事后选支构造时延 PASS，也未生成正式 CIR。

![四点原生反射复平面与锁定残差](../results/validation/deepsea3d_65grid_four_20261027/four_point_reflection_complex_v3.png)

![四点预测复差与 FE 经验敏感性](../results/validation/deepsea3d_65grid_four_20261027/prediction_vs_FE_sensitivity_v3.png)

图只显示四个离散 FE 真点和锁定预测，没有连线推断带内连续值。对应可编辑 [复平面 FIG](../results/validation/deepsea3d_65grid_four_20261027/four_point_reflection_complex_v3.fig)和[误差 FIG](../results/validation/deepsea3d_65grid_four_20261027/prediction_vs_FE_sensitivity_v3.fig)已重开，逐点坐标和三组离散误差标记均与保存 CSV 一致；[FIG 对象审计](../results/validation/deepsea3d_65grid_four_20261027/figure_object_audit.json)为 PASS。早期图形导入失败和图形布局演进的日志/旧图原样保留。

## 审计、成本与范围决定

[独立保存数据审计](../results/validation/deepsea3d_65grid_four_20261027/independent_audit.json)从原生 MAT、逐频 CSV、冻结预测和解析 Green 独立复算四组复数、门槛、`D+R=H`、`R−R0=Δ`、幅相、FE 九变体敏感性、解析 DSP 控制、锁先后顺序、资源与进程；133 个首批回传文件逐 SHA 一致，最终两图另有逐文件 SHA。原 MATLAB 工件断言 PASS；[检查点](../results/validation/deepsea3d_65grid_four_20261027/checkpoint.json)保存输入/代码/输出 SHA、失败进程、资源和下一步门。5/7 和旧四半千赫点只读复用，新增物理频率恰为四个。

本 Goal 的新远端 MATLAB 完整进程墙钟：版本探测保守扣 `95 s`；首次隔离目录缺失在 FE 之前失败 `56 s`；四个成功物理频点 `29/34/27/35 s`；MATLAB 工件审计 `11 s`；首版图源导入失败 `11 s`、成功图与版式修订 `15/14 s`、最终 FIG 对象审计 `14 s`。**总保守扣记 `341/720 s`**，单个新物理频点最大 `35<150 s`；进程最大 RSS `2,720,064 kB`，启动前最低 `MemAvailable=29,187,828 kB`，各作业均留足 3 GiB 增量及至少 30% 物理 RAM，磁盘预留 2 GiB 已核。版本探测的实际输出确认 R2023b，但其完整 stdout/stderr 当时未持久保存；本报告将其 95 s 保守计费并明确这一执行凭据限制，而不伪造日志。无 swap 放行、无并行池或其他用户进程干预。

**决策：**65 点、62.5 Hz 的**预测输出网格**获得固定案例的四个新频点有限筛查支持，最多标记 `COMPACT_DUAL_65GRID_HOLDOUT_SCREENING_PASS`。四点无法约束所有带内极值、未知粗糙散射长尾、频率窗口混叠或时域窗收敛。下一项若获单独授权，应先检查时域窗口/带内极值风险和所需真正独立保留点，再决定是否扩展；本次不自动增加 5187.5/7937.5 Hz、不转 129 点、不扫 65 个 FE、不生成认证 CIR/BER。Q1 BIE `STOP`、二维 S4 `CORE_COMPLETE`、`G1_target_PM=OPEN`、`WP2_formal=false`；生产求解器/API 和 `vertical_comm_guide.md` 实现公式未改。
