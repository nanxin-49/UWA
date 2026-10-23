# K_HIGH：S2c-B 受控波数实验图册 M04（2026-10-23）

**先看 M04b：K 改变了 BIE 参考反射场什么。** 本次唯一新物理海面为 `A=.05 m, K=.2 rad/m, φ=0`，与 `BASE_M` 的区别只是输入 K（基线 K=.1 rad/m）。固定4 kHz、c=1500 m/s、Tx=(0,100) m、984点σ=.3 m离散Gaussian源、原411点 z=3 m线和真实P11–P33九点。`exp(−iωt)` 原生反射场单位是**单位源响应**。

## M04b：BIE 参考物理响应（最先看）

![两海面的BIE原生反射复场及K造成的物理增量](figures_final_v2/M04b_BIE_physical_response.png)

[可编辑FIG](figures_final_v2/M04b_BIE_physical_response.fig) · [逐点原始变化](analysis/BIE_physical_pair.csv)

**看到了什么：**幅度曲线总体接近，但相位和复场明显变化。411点最大 `|Δ_K u_BIE|=0.146457684` 位于 x=18.9453125 m；同点幅度变化最大仅0.132020037 dB（x=−39.0625 m），wrapped 相位变化最大2.8820422 rad（x=−22.65625 m）。五支持集的参考物理变化 L2：整线1.014765，中心11点0.099829，九点0.026930，五点十字0.021366，边缘104点1.074617。

**可以支持什么：**五层新 BIE 边界解已完成，自求积、空间离散、窗口经验检查下，上述物理响应均超过同支持5倍经验分辨率。它是海面参数变化的**参考场响应**，不是 PE/FullFK 的模型误差。相位断点断线，无二维插值。

**还不能证明什么：**BIE 是数值参考而非实验真值；五层差值不是严格误差上界。接收点 x 不是唯一局部反射点，不能据此把变化唯一归于当地几何。

## M04a：完整 C2 海面几何

![BASE_M与K_HIGH的高度坡度曲率](figures_final_v2/M04a_geometry_only.png)

[可编辑FIG](figures_final_v2/M04a_geometry_only.fig) · [解析几何CSV](P0/geometry.csv) · [几何概要](P0/geometry_summary.csv)

**看到了什么：**平台周期62.832→31.416 m、中心坡度0.005→0.010。|x|=42 m为平台边界、|x|=50 m为严格平坦外缘。平台理想 `AK²` 为0.0005/0.0020 m⁻¹；含C2过渡区的实际最大 `|η″|` 却为0.0050572/0.0045973 m⁻¹。

**可以支持什么：**唯一主动改变的输入是K；|x|≥50 m的高度、坡度、曲率严格为0，满足完整FullFK解析平面外场的几何前提。

**还不能证明什么：**K同时改变周期、坡度、曲率、峰谷排列和照射；不能独立识别坡度或曲率的因果效应，也不能以 `AK²` 代替整个海面的最大曲率。

## M04c：两案例模型残差与幅相平方项

![PE和FullFK相对各自BIE的残差与分解](figures_final_v2/M04c_model_error_pair.png)

[可编辑FIG](figures_final_v2/M04c_model_error_pair.fig) · [逐点诊断](analysis/model_point_diagnostics.csv) · [支持集指标](analysis/support_metrics.csv)

**看到了什么：**每列一个方法，两案例共享该列纵轴；列间纵轴不同且明示。展示幅度差dB、wrapped相位差、绝对复残差与精确幅相平方分解。411点复L2，PE由0.0155868变0.0178365，FullFK由1.07328e−5变4.34208e−5。整线残差均主要由相位平方项构成；平方项为代数分解，不是独立机制份额。残差/平方项使用对数轴，零值仅在绘图显示时裁至正下限，CSV原值未改。

**可以支持什么：**在固定同源、同点、同权重下，K_HIGH两种模型的整线误差均增大，变化超过5倍经验门。BIE物理场变化和 `u_method−u_BIE` 模型残差分开计算。

**还不能证明什么：**FullFK残差较小不是实验真值认证。PE横向倍网格只测试冻结13个探针；探针最大差并非严格全线收敛界。

## M04d：整线、中心与真实九点

![411线中心11点与原九点误差](figures_final_v2/M04d_line_center_nine.png)

[可编辑FIG](figures_final_v2/M04d_line_center_nine.fig) · [成对变化与5倍分类](analysis/paired_changes.csv)

**看到了什么：**线图使用真实接收点；九点为真实3×3单元格，无插值。同一方法两案九点共色标。PE九点复L2 `1.27999e−5→3.01212e−5`；FullFK九点 `2.66685e−7→1.98012e−6`。BASE_M FullFK九点标为 `UNRESOLVED`，K_HIGH超过5倍经验门。

**可以支持什么：**中心11点PE `3.0294e−5→6.4999e−5`、边缘104点 `0.0294387→0.0502521`；FullFK中心 `9.0137e−7→7.1490e−6`、边缘 `5.3893e−6→3.6462e−5`。除 **FullFK五点十字的成对误差变化**外，各支持集的两案误差变化在冻结5倍经验门下可分辨；该十字仅报告观察值。

**还不能证明什么：**九点是同一海面的九个接收位置，不是九个独立海面。中心/边缘差异可能共同受入射照射、传播和干涉影响；图不是连续x–z场，空间相关不能单独确定因果机制。

## 数值资格、成本和重绘

五层K_HIGH新边界解依次为L60/12ppw/qa96、同网格qa64、L60/14ppw/qa96、L60/10ppw/qa96、L70/10ppw/qa96；全部417同点有限，线性残差<1e−8。BIE经验分辨率：411线2.5851e−7、九点1.5040e−7，[逐支持数据](P1/BIE_resolution_support.csv)。FullFK为 `FullFK_complete_analytic_flat_exterior_GL8_v1`，真实法向/弧长、`local_radius_m=Inf`、rough-minus-flat双项和解析平面外场；GL8 400/800变化9.7649e−14，L60/90对照7.3729e−16。PE平面控制3.8653e−13。资格只限冻结接收点与经验检查。

用户于代表作业后明确授权[仅K_HIGH子额度2400 s](protocol/P1_budget_amendment_2400.json)，原1800 s停止记录保留。七次新声学作业实计1277.0539592 s，核心累计2685.7140642/7200 s，S3A规划余额4514.2859358 s。每次作业留有即时物理RAM/30%门，最长作业352.726 s<900 s。

[独立最终审计v4](audit/final_independent_audit_v4.json)从CSV复算2502原始场、1668逐点、20支持指标、10配对行、共同坐标/源SHA和幅相恒等式，最大逐点差8.01e−14、恒等式余量2.59e−20；M04a另由P0几何MAT独立重绘；M04b–d由[保存MAT](analysis/analysis_data.mat)经[重绘入口](../../../../scripts/reporting/render_two_dimensional_closeout_s2c_b_khigh.m)重绘PNG像素一致，四FIG重开；[独立重绘目录v4](audit/independent_redraw_v4/)保留。M04a图题保留P0阶段的“PRE-EXPERIMENT GEOMETRY ONLY”，仅表示它是解析几何图，不否定本包随后完成的声场计算。初版最终审计及原先单层BIE/1800 s受阻的[P1审计](audit/P1_independent_audit.json)、[几何初版](figures_v2/M04a_geometry_only.png)仅作历史快照。

**下一阶段建议：**PHI_PLUS90值得作为另行授权的独立单案例评估，以补充本次K变化混合多几何因素的观察。应先以本次实耗和当时RAM重过门，保留同一417点及五层参考；本轮不启动它。单频空间结果不推出PM、三维、宽带或通信结论。`G1_target_PM=OPEN`、`WP2_formal=false`。

图形程序性参考：Kassis et al. (2026), *Scientific Agent Skills: A Library of Procedural Knowledge for Research Agents*, [arXiv:2609.00065](https://doi.org/10.48550/arXiv.2609.00065)。此文不作为声学数值或物理证据。
