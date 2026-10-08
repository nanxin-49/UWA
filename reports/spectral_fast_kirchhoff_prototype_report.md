# Spectral Fast Kirchhoff prototype

日期：2026-10-07。Validation-only。

Numerical accuracy pass: **1**. Combined prototype success gate: **1**.

完整离散推导、近似分类、FFT 归一化与边界定义见 [method note](spectral_fast_kirchhoff_prototype_method.md)。

## 主案例与 rank=1–4

入射 rank 分别为 1–4；出射统一 rank=6，绝不把入射 rank 冒称全核 rank。Rx=3 m、sigma=.3 m、4 kHz。Prototype-A 的几何简化和 Prototype-B 的完整法向/测度分别记录。

|case|geometry|r_in|r_out|complex L2|magnitude L2|phase RMS rad|
|---|---|---:|---:|---:|---:|---:|
|weak|A|1|6|0.00070402285|1.8399789e-05|0.00096178093|
|weak|A|2|6|2.6566713e-06|1.7910103e-06|2.0995888e-06|
|weak|A|3|6|2.0758107e-06|6.7050913e-07|2.1021799e-06|
|weak|A|4|6|2.0737863e-06|6.7050685e-07|2.0996507e-06|
|weak|B|1|6|0.00070337196|1.8386075e-05|0.00096190539|
|weak|B|2|6|1.3852639e-06|1.3828646e-06|9.860981e-08|
|weak|B|3|6|2.8992587e-09|4.2974259e-10|3.4725316e-09|
|weak|B|4|6|7.4000095e-12|6.7988092e-12|3.4030973e-12|
|strong_height|A|1|6|0.014071034|0.00065913613|0.019243584|
|strong_height|A|2|6|0.00071676951|0.00071541696|5.0120205e-05|
|strong_height|A|3|6|0.00027678809|0.00027008467|6.636127e-05|
|strong_height|A|4|6|0.00027282006|0.00026964327|4.3575093e-05|
|strong_height|B|1|6|0.014056227|0.00058338948|0.019246218|
|strong_height|B|2|6|0.00055388907|0.00055247421|4.8488731e-05|
|strong_height|B|3|6|3.2078354e-05|1.5075789e-05|3.1645934e-05|
|strong_height|B|4|6|1.9785813e-05|1.511702e-05|1.2096964e-05|
|high_K|A|1|6|0.0020952712|6.00637e-05|0.0029209591|
|high_K|A|2|6|0.00010860454|5.4768233e-05|9.4260168e-05|
|high_K|A|3|6|0.00010823866|5.4089058e-05|9.4250332e-05|
|high_K|A|4|6|0.00010824025|5.408904e-05|9.4249959e-05|
|high_K|B|1|6|0.0020368155|2.6528049e-05|0.0029153951|
|high_K|B|2|6|6.0023588e-06|5.9959008e-06|3.6412529e-07|
|high_K|B|3|6|2.1155262e-08|3.9414188e-09|2.6008137e-08|
|high_K|B|4|6|1.1425187e-10|9.8106903e-11|7.1638694e-11|
|K070|A|1|6|0.0021706843|0.0001276037|0.0029101901|
|K070|A|2|6|0.00024051152|0.00012011381|0.00020907339|
|K070|A|3|6|0.00024023476|0.00011963423|0.00020906782|
|K070|A|4|6|0.00024023634|0.0001196342|0.00020906744|
|K070|B|1|6|0.002035623|4.3227587e-05|0.0028928656|
|K070|B|2|6|6.007765e-06|5.9644658e-06|9.061657e-07|
|K070|B|3|6|2.1446793e-08|7.0031607e-09|2.5406299e-08|
|K070|B|4|6|1.0608119e-10|8.6778451e-11|7.8821367e-11|
|sigma1|A|1|6|0.0010227931|0.00048439231|0.0020036663|
|sigma1|A|2|6|0.0003824667|0.00038172053|2.993399e-05|
|sigma1|A|3|6|0.00038071747|0.00037996636|2.9922827e-05|
|sigma1|A|4|6|0.00038071773|0.0003799664|2.9927124e-05|
|sigma1|B|1|6|0.00095547853|0.00030371688|0.0020173165|
|sigma1|B|2|6|2.4977358e-05|2.3202126e-05|1.142122e-05|
|sigma1|B|3|6|2.2968075e-05|2.1031904e-05|1.1356867e-05|
|sigma1|B|4|6|2.296921e-05|2.1031404e-05|1.1361108e-05|
|Rx2|A|1|6|0.014135797|0.00067225831|0.019706639|
|Rx2|A|2|6|0.00071795175|0.00071655844|5.0717553e-05|
|Rx2|A|3|6|0.00027673614|0.0002698939|6.6836184e-05|
|Rx2|A|4|6|0.00027272124|0.0002694566|4.396859e-05|
|Rx2|B|1|6|0.014121244|0.00059771453|0.019709741|
|Rx2|B|2|6|0.000556064|0.00055459493|4.9134606e-05|
|Rx2|B|3|6|3.3702052e-05|1.6478861e-05|3.253731e-05|
|Rx2|B|4|6|2.1891148e-05|1.6571336e-05|1.3518631e-05|

比较对象是同一离散曲面、同一入射场定义上的 FullFK 原生场，采用独立固定能量权重与 FullFK-only -40 dB 掩膜。主案例权重来自保存 Stage-0 M99；sigma=1 控制使用冻结源参考面能量的固定 99% footprint。没有逐 rank 掩膜调整、拟合幅度或定相位。这里不以 FastFK–BIE 作为快速算法保真指标。

选定 Prototype-B 3/6 的**不加权全 4097 点原生场 L2**也全部通过，不依赖 M99/幅度门：

|case|whole-grid complex L2|
|---|---:|
|weak|3.6733694e-09|
|strong_height|3.5360728e-05|
|high_K|2.9034529e-08|
|K070|2.995891e-08|
|sigma1|2.2461573e-05|
|Rx2|3.6568069e-05|

该全网格检查加入主 accuracy gate。

## 出射 rank 与 padding

|case|study|r_in|r_out|padding|complex L2|
|---|---|---:|---:|---:|---:|
|weak|padding|3|6|1|0.0093436353|
|weak|padding|3|6|2|2.8992587e-09|
|weak|padding|3|6|4|2.8992587e-09|
|weak|outgoing_rank|4|1|2|0.069842563|
|weak|outgoing_rank|4|2|2|0.00075773149|
|weak|outgoing_rank|4|3|2|4.9570178e-06|
|weak|outgoing_rank|4|4|2|2.4140213e-08|
|weak|outgoing_rank|4|5|2|9.3949911e-11|
|weak|outgoing_rank|4|6|2|7.4000095e-12|
|strong_height|padding|3|6|1|0.0093472412|
|strong_height|padding|3|6|2|3.2078354e-05|
|strong_height|padding|3|6|4|3.2078354e-05|
|strong_height|outgoing_rank|4|1|2|1.0156087|
|strong_height|outgoing_rank|4|2|2|0.34765264|
|strong_height|outgoing_rank|4|3|2|0.043613604|
|strong_height|outgoing_rank|4|4|2|0.0040459861|
|strong_height|outgoing_rank|4|5|2|0.00030583256|
|strong_height|outgoing_rank|4|6|2|1.9785813e-05|
|high_K|padding|3|6|1|0.0093440868|
|high_K|padding|3|6|2|2.1155262e-08|
|high_K|padding|3|6|4|2.1155262e-08|
|high_K|outgoing_rank|4|1|2|0.15166823|
|high_K|outgoing_rank|4|2|2|0.0029739264|
|high_K|outgoing_rank|4|3|2|3.7432309e-05|
|high_K|outgoing_rank|4|4|2|3.3547813e-07|
|high_K|outgoing_rank|4|5|2|2.3566984e-09|
|high_K|outgoing_rank|4|6|2|1.1425187e-10|
|K070|padding|3|6|1|0.009344646|
|K070|padding|3|6|2|2.1446793e-08|
|K070|padding|3|6|4|2.1446793e-08|
|sigma1|padding|3|6|1|0.0092896614|
|sigma1|padding|3|6|2|2.2968075e-05|
|sigma1|padding|3|6|4|2.2968075e-05|
|Rx2|padding|3|6|1|0.0092417533|
|Rx2|padding|3|6|2|3.3702052e-05|
|Rx2|padding|3|6|4|3.3702052e-05|

Padding=1 是有限 lag 序列的循环折叠，故意作为 wrap-around 失败控制；2/4 的长度均 >=2N-1，保持同一线性 Toeplitz 卷积，差异应接近机器精度。没有 sponge、padding-dependent 物理支撑改变或 cutoff 角度。

## 时间、内存与 scaling

MATLAB 25.2.0.3312555 (R2025b) Update 6；FastFK warm apply 五次中位数，FullFK dense apply 两次中位数。cold 总时间分别包含 FastFK 全局 basis/kernel setup 与 FullFK exact incident construction。

|N|FFT size|Full apply s|Fast setup s|Fast apply s|warm speedup|cold speedup|Fast plan MiB|Fast setup arrays MiB|Fast apply arrays MiB|
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
|513|1026|0.0402989|0.0569074|0.0026287|15.33|0.89702|0.52967|6.3021|0.10959|
|1025|2050|0.14471|0.0742725|0.0049341|29.329|2.1784|1.0531|8.6367|0.21896|
|2049|4098|0.553946|0.0962307|0.0080801|68.557|5.8111|2.1|13.306|0.43771|
|4097|8194|2.25483|0.143809|0.0160073|140.86|14.73|4.1937|22.644|0.87521|

log-log empirical exponents: Full apply 1.937, Fast apply 0.8537, Fast setup+apply 0.4675。四个有限样本上的经验斜率不是渐近复杂度证明。Fast 路径中没有 N_rx×N_surface 或 N_surface×N_k 的矩阵；固定 q=33 的 lag×height setup 是 O(q²N)，apply 是 O(r_out N log N)。入射采样 IFFT 及 height-table 插值为 O(r_in N log N + q r_in N)。表中数组存储和 FullFK block estimates 不是实测 peak RSS；SVD/FFT/Bessel workspace 未作峰值监控。CPU 秒、mask、所有单次 case costs、Full incident costs、grid error 都在 CSV。

完整 numerical workflow（不含随后纯重绘/报告刷新）wall 27.904 s、进程 CPU 195.11 s；多线程 CPU 可以大于 wall。本阶段不改变生产信道，因此没有运行被任务排除的宽带/BER/通信回归。独立 test_spectral_fast_kirchhoff_fft_controls 的 flat Hankel/物理源导数/padding 和复多项式插值控制已通过 2 tests；最终保存验证结果在 test_results.csv。

## 验证控制与出射谱

Independent odd/even complex Toeplitz FFT control p=1/2/4 errors: 0.3761 / 2.304e-16 / 2.332e-16。Chebyshev interpolation analytic exponential error: 2.716e-16。

保存谱 CSV 及 field MAT 包含所有 receiver-grid 点。出射谱使用同一 [-60,60] m receiver-plane Hann 窗、同一 FullFK 能量归一；它是有限接收孔径诊断，不能冒称无限面远场散射谱。propagating 与超出 k 的箱均记录，不以 theta99_inc 截断出射。Hankel 核通过有限 lag FFT 保留采样后的传播与近场/倏逝效应。

|case|incident theta95/theta99 deg|FK reflected theta95/theta99 deg|Fast reflected theta99 deg|FK energy outside incident99|spectrum complex L2|
|---|---|---|---:|---:|---:|
|weak|15.9693 / 21.2801|12.2659 / 15.7719|15.7719|0.000319524|3.12169e-09|
|strong_height|15.9693 / 21.2801|12.0827 / 16.5173|16.5173|0.000811755|3.30533e-05|
|high_K|15.9693 / 21.2801|12.4491 / 16.1443|16.1443|0.000323914|2.22608e-08|
|K070|15.9693 / 21.2801|12.4491 / 16.1443|16.1443|0.000450993|2.25778e-08|
|sigma1|4.70073 / 6.27314|5.73776 / 6.81829|6.81829|0.0275885|2.25087e-05|
|Rx2|15.9693 / 21.2801|12.0827 / 16.7041|16.7041|0.000870337|3.45554e-05|

强高度 sigma=1 控制的出射能量在 incident99 区间外约 2.76%，证明直接套用入射角宽会丢失有效散射。本实现使用全 lag 核而没有该截断，high-K 及 K=.70 的对应谱仍与 FullFK 一致。低入射 rank **不能独自保证**出射带宽；本轮经过独立 outgoing rank=6 与全核 FFT 检查后，才支持这几个案例。

入射复谱与归一能量保存在各 case 的 incident_spectrum.csv；angular_summary.csv 保存本表及超出 |kx|=k 的有限孔径谱能量诊断。CPU 中小于计时分辨率的 apply 可记录为 0，尤其 Windows cputime 的约 15.6 ms 粒度；它不代表零计算成本。cold/warm 所有速度比使用实际 wall time。

## 图与复现

![fields_and_rank](../results/validation/spectral_fast_kirchhoff_prototype/figures/fields_and_rank.png)

![padding_audit](../results/validation/spectral_fast_kirchhoff_prototype/figures/padding_audit.png)

![reflected_spectra](../results/validation/spectral_fast_kirchhoff_prototype/figures/reflected_spectra.png)

![cpu_scaling](../results/validation/spectral_fast_kirchhoff_prototype/figures/cpu_scaling.png)

`setup_vertical_project; result=validate_spectral_fast_kirchhoff_prototype()`；纯重绘：`validate_spectral_fast_kirchhoff_prototype('render')`。入口不改变 cwd。输出在 results/validation/spectral_fast_kirchhoff_prototype；spectral_fast_fk_prototype.mat、all_metrics.csv、scaling.csv、各 case field MAT/spectrum CSV/metrics CSV 和 PNG/PDF。

## 判断

**Spectral Fast FK prototype 成功。下一阶段建议进入固定 PM random-surface pilot。** 本任务在该建议处停止，不运行 PM、Monte Carlo、宽带或生产集成。选定 global ranks 是 incoming=3/outgoing=6；相同输入强相位完整保留、真实高度全局非局部相干耦合完整覆盖；速度和精度结论只适用于本报告的二维确定性案例。
