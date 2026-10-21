# 两套可视化图册的数据说明

> **2026-10-21：** [高度相位符号审计](pe_bie_height_phase_sign_audit.md)发现旧 PE 高度符号与 BIE 几何约定的衔接问题。下文数据字段和历史重绘说明仍有效，但旧单边共轭的 G 比较不能直接当作已统一物理表示的误差。原数据保留；当前修正包的字段、数值和重绘入口见[新图册报告](pe_bellhop_bie_corrected_visualization_report.md)。

更新日期：2026-10-19。配套[展示汇总报告](pe_bellhop_bie_visualization_summary.md)。本说明依据保存 MAT 的实际字段、CSV 表头、绘图代码和审计文件整理；不新增声学计算。

## 1. 先选对数据包

|你要做什么|应使用的数据|不要混淆的地方|
|---|---|---|
|比较弱、中、强海面的三方法差异|第一包 comparison_data.mat|这里的主要数组是粗糙/平面复数比 G，不是原始声压|
|画近海面真实二维幅值、相位|第二包 xz_data.mat|只计算了强海面，不能把它扩充成三档二维结果|
|改图例、字体、布局|对应 MAT 与纯绘图入口|不需要重新运行 PE、BIE 或 Bellhop|
|在表格软件检查数值|逐案例 CSV 或 xz_fields.csv|必须保留实部与虚部；只有幅值不足以重建相位|
|追溯程序、输入或原始输出|各包 source_data/ 和 package_manifest.csv|绘图数据是便于使用的副本，原始作业仍保留|

第一包目录（从仓库根目录起）：

```text
results/visualization/pe_bellhop_bie_comparison/20261018_194331_523/
```

第二包目录：

```text
results/visualization/pe_bellhop_bie_xz/20261019_115214_348/
```

## 2. 第一包：三档海面、3 m 接收线

### 核心文件

|相对包目录的路径|内容|
|---|---|
|plot_data/comparison_data.mat|自包含绘图数据，顶层变量为 data|
|plot_data/weak_low_K.csv|弱海面 A=0.01 m，每行一个接收位置|
|plot_data/region_II_low_K.csv|中海面 A=0.05 m|
|plot_data/strong_height_low_K.csv|强海面 A=0.20 m|
|metrics/current_metrics.csv|本轮从保存 G 数组重算的 9 行配对指标|
|metrics/beam_sensitivity.csv|5001/10001 波束敏感性|
|metrics/historical_report_values.csv|历史报告值；与 current_metrics.csv 按案例和配对核对|
|metrics/artifact_audit.json|数据、表示转换及原有重绘检查|
|figures/|7 张 PNG/PDF，Bellhop–BIE 和三方比较分别归档|

### MAT 中常用字段

|字段|含义、单位与使用方法|
|---|---|
|data.x_m|984×1 接收横坐标，m；深度固定为 3 m|
|data.mask|984×1 逻辑掩膜，固定 Stage-0 M99 支持，共 411 个点|
|data.weights|固定入射能量权重；做指标时在所用固定支持内归一化|
|data.cases(j).name / amplitude_m|案例名、海面高度参数 A，m；j=1/2/3 对应弱/中/强|
|data.cases(j).eta_m|接收横向网格对应的海面高度，m|
|data.cases(j).G_PE / G_BH / G_BIE|用于统一定量比较的复数 G，无量纲；984 个采样位置|
|data.cases(j).G_BH_native / G_BIE_native|保留原生表示的比值；名称中的 native 不代表原始声压 u|
|data.cases(j).G_BH_5001|较少波束数的比较比值，用于敏感性检查|
|data.cases(j).plot_mask|绘图及该案例既定指标使用的共同支持|
|data.cases(j).valid_mask / BH5001_valid_mask|比值有效性逻辑掩膜，不是幅值权重|
|data.cases(j).residuals / transforms|配对残差、幅相及相位断点等绘图变换|
|data.metrics / beam_sensitivity / historical|本轮指标、波束敏感性、历史记录|
|data.config / convention / definitions|配置、复数表示规则和指标定义|

M99 是此前固定的入射能量支持，不是“选取误差最小的 99% 点”。本次比较沿用这 411 个点，约 x=−40.039～40.039 m，间距 0.1953125 m；没有为追求接近而重新挑点。

Bellhop 有 375 个支持外位置的平面场为零，G 无法定义，保存为复数 NaN。CSV 中 ratio_valid、BH5001_ratio_valid 标出可用比值；这些无效点均在 M99 外。PE/BIE 的支持外存储值也不能不看掩膜直接当成有效比较数据。**不要把 NaN 填零后计算全域指标。**

### CSV 与指标读法

CSV 中 PE_real/PE_imag、BH_real/BH_imag、BIE_real/BIE_imag 表示比较用 G 的实部/虚部。BH_native_real/BH_native_imag 在第一包中也是比值的原生表示，不是第二包的 BH_native 声场。incident_weight、M99_mask、ratio_valid 负责权重和有效范围。

配对 A–B 的三个指标为：

$$
E_{AB}=\sqrt{\frac{\sum_i w_i|G_A-G_B|^2}{\sum_i w_i|G_B|^2}},\qquad
\delta\phi_i=\arg(G_A\overline{G_B}),\qquad
\delta a_i=20\log_{10}\frac{|G_A|}{|G_B|}.
$$

相位 RMS 为 √Σwᵢδφᵢ²，单位 rad；幅度比差值 RMS 为 √Σwᵢδaᵢ²，单位 dB。w 在固定 M99 上归一化为和 1。配对中的第二种方法是相对复差分母，因此 A–B 与 B–A 的相对复差不必完全相同。

current_metrics.csv 的 complex_l2 是小数，展示成百分数需乘 100；phase_rms_rad 转 mrad 乘 1000，转度乘 180/π。magnitude_ratio_rms_db 已经是 dB，不再取一次对数。sample_count=411 是位置数，不是随机试验次数。三档确定性海面也不是统计重复实验。

## 3. 第二包：强海面的真实二维原始场

### 核心文件及字段

|相对包目录的路径|内容|
|---|---|
|plot_data/xz_data.mat|自包含二维绘图数据，顶层变量 data|
|plot_data/xz_fields.csv|18495 行，包含 x、z、复数实虚部、幅值和相位|
|plot_data/pe_bie_fields.mat|补齐 Bellhop 前的 PE/BIE 场评价检查点，通常无需使用|
|figures/plot_transforms.mat|幅相合成所用 RGB 数组、映射公式等|
|metrics/field_checks.json / artifact_audit.json|3 m 截线一致性、坐标、CSV、原有独立重绘检查|
|source_data/numeric_inputs.mat|保存的 PE 起始场、BIE 边界密度等数值输入|
|source_data/source_conventions.json|源标度与相位参考说明|

|MAT 字段|含义与维度|
|---|---|
|data.x_m / z_m|横坐标 m、深度 m；x 有 411 点，z 有 45 层|
|data.PE_envelope|原始约化 PE 包络，411×45 复数数组|
|data.PE_carrier_restored|PE 包络恢复解析载波后用于波前显示的场，411×45|
|data.BIE_native|保存边界密度所评价的原生反射场，411×45|
|data.BH_raw_SHD|从 Bellhop SHD 读出的原值，411×45|
|data.BH_native|按既定约定 conj(BH_raw_SHD) 得到的原生场，411×45|
|data.magnitude / phase_rad|各方法的 abs(u)、arg(u)；相位为主值，单位 rad|
|data.conventions / source_diagnostic|表示规则和源幅相诊断，诊断值没有作为拟合系数应用|
|data.checks / geometry|保存截线核对、坐标和反射诊断|

**数组第一维是 x，第二维是 z。**z=0.5:0.125:6 m；CSV 每行用 x_m、z_m 明确标出位置。按 MATLAB 列展开顺序，同一深度的 411 个 x 点先连续存储，再进入下一个深度。CSV 与 MAT 的这些场量满足 u=real+1i×imag、幅值=abs(u)、相位=angle(u)。

PE 显示场使用：

$$
u_{\mathrm{PE,显示}}(x,z)=\psi_{\mathrm{PE}}(x,z)\exp\{ik_0(100+z)\},\qquad k_0=2\pi\times4000/1500.
$$

该因子恢复载波相位，不改变幅值，也没有更改原 Model-0 海面相位屏。第一包针对 G 的固定共轭约定不能直接当作第二包原始 u 的完整跨方法相位转换。

### 色标和比较边界

PE/BIE 幅值上限共用约 0.12248371；Bellhop 原生源上限单列约 0.70721499。它们是源单位，未标定为 Pa，不能据此宣称绝对声压或 TL 相同。原生 Bellhop 源标度及参考相位与 Gaussian 场有差别，没有拟合系数将其强行重合。

相位色条为 [−π,π] 循环。合成图的色相 h=(arg(u)+π)/(2π)、饱和度 1、亮度 v=abs(u)/对应上限，线性映射，无 gamma 增强。深度轴仅在版面上展开；数值网格没有插值或复制行。

## 4. 数据如何得到、检查到了哪一步

第一包复用保存的 PE/BIE 复数曲线，新增 8 个 Bellhop 作业：平面加三档粗糙面，各 5001/10001 波束。第二包从保存 PE 反射起始场传播到 45 个深度，用保存的 4480 节点 BIE 密度评价场，并新增 1 个 Bellhop 二维网格作业。BIE 边界方程没有重新求解。

|既有检查|三案例图册|二维场图册|
|---|---|---|
|独立保存数据重绘|7/7 PNG 逐像素一致|4/4 PNG 逐像素一致|
|CSV 最大绝对往返差|约 4.89×10⁻¹⁵|约 7.11×10⁻¹⁶|
|与旧 3 m 截线核对|保留原有接收线|PE 约 1.48×10⁻¹⁴；BIE 约 5.26×10⁻¹⁵；Bellhop 为 0|
|PNG / PDF|2400×1350；单页矢量 PDF|2400×1350；审核 PNG 整页封装 PDF，非矢量|

这些核对说明数据导出、保存和绘图相互一致，不等价于证明物理模型在所有海况下准确。BIE 12/14 ppw 保存密度在新域 44 个探针处的场变化约 1.42×10⁻⁸，只是敏感性检查，不是整个二维域的严格误差界。原 BIE 目录缺失的来源问题及派生副本交叉核对记录保留在第一包原报告中。

## 5. 只改版式时怎样操作

把 project_root 换成当前仓库路径。以下命令从已有数据重绘，不运行物理求解；为保留原交付，将两包副本复制到一个**尚不存在**的工作目录后操作。

```matlab
project_root = 'E:\MISC\CARPE3D_matlab\Explain';
addpath(project_root);
setup_vertical_project;
base = fullfile(project_root,'results','visualization');
work = fullfile(base,'pe_bellhop_bie_presentation_work');
assert(~isfolder(work),'请换一个新的工作目录，保留已有结果');
mkdir(work);

% 只复制重绘所需 MAT，不重新计算声场。
line_copy = fullfile(work,'line');
xz_copy = fullfile(work,'xz');
mkdir(fullfile(line_copy,'plot_data'));
mkdir(fullfile(xz_copy,'plot_data'));
copyfile(fullfile(base,'pe_bellhop_bie_comparison','20261018_194331_523', ...
    'plot_data','comparison_data.mat'),fullfile(line_copy,'plot_data'));
copyfile(fullfile(base,'pe_bellhop_bie_xz','20261019_115214_348', ...
    'plot_data','xz_data.mat'),fullfile(xz_copy,'plot_data'));
generate_pe_bellhop_bie_comparison_visuals(line_copy);
generate_pe_bellhop_bie_xz_visuals(xz_copy);
```

如果仅需读值，load 对应 MAT 的 data 或 readtable 对应 CSV 即可。不要为改报告文字调用 prepare_* 数据准备入口。离开仓库后重绘，需要将原包 source_data/workflow_code 中相应绘图函数加入路径；绘图数据包本身不要求 Bellhop 可执行程序。

## 6. 本次汇总交付包含什么

汇总报告引用两套既有数据包，另提供仓库内的 Obsidian 阅读副本：报告、数据说明和 11 张原 PNG。阅读副本不复制 MAT、CSV、程序或原始作业输出；这些继续保存在上文两个可追溯的源目录。只把阅读副本移入 Obsidian 时，图文可以独立阅读，数值数据需另行取用。没有写入未指定的外部 vault。
