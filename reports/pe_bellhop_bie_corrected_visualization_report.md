# PE–Bellhop–BIE 修正图册：同一海面约定下看接收场

日期：2026-10-21。当前展示请优先使用[本轮图册与完整数据](../results/visualization/pe_bellhop_bie_corrected/20261021_native_reset_v2/corrected_report.md)。[旧七图包](../results/visualization/pe_bellhop_bie_comparison/20261018_194331_523/comparison_report.md)和[旧二维包](../results/visualization/pe_bellhop_bie_xz/20261019_115214_348/xz_field_report.md)保留为历史结果；旧指标不再代表统一高度符号后的 PE 模型误差。

## 一眼看懂修正

![三档海面上旧PE、修正PE与BIE的接收线比较](../results/visualization/pe_bellhop_bie_corrected/20261021_native_reset_v2/figures/R01_height_sign_reset.png)

**图 R01。** 上排直接画 |G|−1：零线代表粗糙面与平面反射的幅值相同；紫色是 BIE，橙色是修正 PE，灰色是旧 PE。灰色在紫色上升处下降，显示旧图中提到的“围绕 1 镜像”；橙色与紫色的起伏转为同向。下排画的是各自约定下的有符号相位差，橙色和灰色对应的参照表示不同，不能把两条曲线直接相减来估算修正量。图中无需幅度缩放或逐案例定相位。

![修正后的三方法接收线反射场幅值与相位](../results/visualization/pe_bellhop_bie_corrected/20261021_native_reset_v2/figures/three_way/TR01_fields.png)

**图 TR01。** 三列为 A=0.01、0.05、0.20 m 的同波数海面；上排是粗糙/平面复数比 G 的幅值，下排是主值相位。现在三条幅值曲线的空间趋势一致，尤其强海面的起伏方向不再对立。相位曲线肉眼重合仍不能证明复数相等；强海面 PE–BIE 相位 RMS 约 62.94 mrad。相位跳到色轴另一端是主值绕回，不是声场消失。

![强海面的真实二维原始反射场幅值](../results/visualization/pe_bellhop_bie_corrected/20261021_native_reset_v2/figures/XZ01_raw_magnitude.png)

**图 XZ01。** 显示同一强海面 x–z 网格上直接计算的原始复反射场幅值，横向 411 点、深度 45 层（0.5–6 m）。PE/BIE 共享一个幅值色标；Bellhop 因源单位不同使用自己的色标。三列共同展示声束包络，不能仅凭彩图颜色相似断言复相位或跨源绝对声压一致。其余 [XZ02 相位](../results/visualization/pe_bellhop_bie_corrected/20261021_native_reset_v2/figures/XZ02_raw_phase.png)、[XZ03 幅相合成](../results/visualization/pe_bellhop_bie_corrected/20261021_native_reset_v2/figures/XZ03_complex_field.png) 和 [XZ04 深度剖面](../results/visualization/pe_bellhop_bie_corrected/20261021_native_reset_v2/figures/XZ04_raw_profiles.png) 保留原始场读法，PE 数值已换成修正后的保存起始场传播结果。

## 改了什么、没改什么

固定环境为二维均匀水体 c=1500 m/s、4 kHz、Tx 深 100 m、Rx 深 3 m，Gaussian 初始场 σ=0.3 m，Bellhop 匹配 X 线源。三个海面共用 K=0.10 rad/m，η=A sin(Kx)w(x)，|x|≤42 m 维持原形、|x|≥50 m 归零。接收线指标继续使用 Stage-0 固定 411 点 M99 入射能量权重。

旧固定比较把同一个 η 用作 BIE 向下的物理边界 z=η，PE 却使用 exp(+2ikη) 相位屏，再只将 BIE/Bellhop 比值共轭。按时间 exp(−iωt)、深度向下的约定，平面整体下移 h 的反射路径缩短 2h，粗糙/平面比是 exp(−2ikh)。因此本轮在**保存的同一海面入射场**上使用 exp(−2ikη) 重算 PE，并用同一 FFT 传播核评价接收线及真实二维网格。BIE、Bellhop 原生复场没有重新求解或拟合；三方法的 G 都在原生表示比较。生产 API、旧图及旧数值文件均未改动。

|强海面、固定 M99|旧固定比较|本轮统一约定|
|---|---:|---:|
|PE–BIE 加权复 L2|8.0908%|6.2911%|
|PE–BIE 幅值相对 L2|约 1.0016%|0.10485%|
|PE–BIE 相位 RMS|80.336 mrad|62.936 mrad|
|PE–BIE dB 幅值比 RMS|0.086996 dB|0.0091097 dB|
|Bellhop–BIE 加权复 L2|0.009976%|0.009976%|

复差可逐点精确分为幅值平方项 (|A|−|B|)² 与相位平方项 2|A||B|(1−cos Δφ)。强海面 PE–BIE 的归一化两项分别为 1.0994×10⁻⁶ 与 3.9567×10⁻³，相位项占复差平方的 99.9722%。因此“新图幅值几乎贴合”和“复 L2 仍有 6.29%”可以同时成立。剩余差异没有被本轮归因给单一物理机制；不能声称 PE 误差已彻底解决。

## 数据与复现

数据包：[完整复数接收线 MAT](../results/visualization/pe_bellhop_bie_corrected/20261021_native_reset_v2/plot_data/comparison_data.mat)、[完整复数二维 MAT](../results/visualization/pe_bellhop_bie_corrected/20261021_native_reset_v2/plot_data/xz_data.mat)、[本轮与旧值逐行表](../results/visualization/pe_bellhop_bie_corrected/20261021_native_reset_v2/metrics/corrected_vs_legacy.csv)、[复差平方分解](../results/visualization/pe_bellhop_bie_corrected/20261021_native_reset_v2/metrics/complex_error_decomposition.csv)、[源副本与哈希](../results/visualization/pe_bellhop_bie_corrected/20261021_native_reset_v2/source_data/source_hashes.csv)、[独立审计](../results/visualization/pe_bellhop_bie_corrected/20261021_native_reset_v2/metrics/corrected_audit.json)。逐案例接收线和二维场另有保存实部/虚部的 CSV；不是只保存彩图像素。

在 MATLAB 已运行 `setup_vertical_project` 后，可从任意当前目录执行：

```matlab
r = prepare_pe_bellhop_bie_corrected_data('新的唯一运行名');
generate_pe_bellhop_bie_corrected_visuals(r.run_dir);
generate_pe_bellhop_bie_corrected_visuals(r.run_dir, fullfile(r.run_dir,'rebuilt_figures'));
audit_pe_bellhop_bie_corrected(r.run_dir);
```

第一个入口只对保存的入射场作 3 个接收线与 45 个深度的 FFT 场评价；新增声源传播、BIE 边界求解、Bellhop 作业均为零。第二个入口只读取本次 MAT，调用包内保存的绘图源码生成 12 张 2400×1350 PNG/PDF。平面 PE 分母相对闭合为 1.82×10⁻¹⁴，强海面二维 z=3 m 与接收线闭合为 1.04×10⁻¹⁶。本轮不改变公共信道函数或通信链，因此无需生产标量/宽带回归；上述检查只验证固定可视化数据链。

独立审计结果为 `PASS_SAVED_DATA_AND_PIXEL_REDRAW`：12 张 PNG 的像素差为 0，接收线 CSV 与 MAT 最大复数差 4.9×10⁻¹⁵，二维 CSV 与 MAT 最大复数差 7.1×10⁻¹⁶，三份包内源 MAT 的 SHA256 全部匹配。

结论仅适用于这三个确定性海面、单频、共同接收线和强海面近海面二维网格。BIE 使用保存的派生副本，本轮没有重新检查边界收敛；5001/10001 波束差只表示敏感性，不是严格误差上界。不能从单频空间图推断时域 CIR、PM 随机海面、三维场或通信性能。
