# Bellhop–BIE 与 PE–Bellhop–BIE 对比图册

固定二维夹具：4 kHz、c=1500 m/s、Tx100/Rx3、Gaussian sigma=.3 m、Bellhop X线源。海面K=.10 rad/m，A=.01/.05/.20 m；C2 plateau42/support50 m。比较量为接收端反射场的粗糙/平面复数比G，不是绝对声压或传播损失。

|参数|本轮统一设置|
|---|---|
|频率 / 声速|4000 Hz / 1500 m/s|
|物理坐标|海面z=0，深度向下；Tx=(0,100 m)，Rx为z=3 m横向线|
|声源|Gaussian sigma=0.3 m；Bellhop X线源使用Stage-0匹配指向性|
|海面|eta=A sin(0.10x) w(x)，A=0.01/0.05/0.20 m|
|C2包络|abs(x)<=42 m保留，abs(x)>=50 m归零，中间五次多项式过渡|
|PE|保存的原始Model-0相位屏|
|Bellhop|压力释放内部壁，[-80,80] m / 4097点，步长0.05 m，角度[-30,30] deg|
|波束|10001主结果；5001作敏感性|
|接收线|984点，周期网格总长192.1875 m，间隔0.1953125 m|
|指标支持|Stage-0 M99共411点，x=-40.03906至40.03906 m；权重为固定入射能量|

指标定义：E_AB=sqrt(sum(w |G_A-G_B|^2)/sum(w |G_B|^2))；相位RMS=sqrt(sum(w angle(G_A conj(G_B))^2))；幅度比差值RMS=sqrt(sum(w [20log10(|G_A|/|G_B|)]^2))。所有sum均在同一M99上，w在该支持内归一化为和1。

## 数据和物理范围

PE Model-0与BIE来自活动目录的误差分解MAT派生副本；原R3–R5活动目录当前缺失，没有重算PE或BIE。复用程序时另发现results/archive旧工作树快照保留原曲线，已用纯数值HDF5读取核对三案例PE/BIE与派生副本逐点相同（相对差0），只作来源交叉核对。Bellhop补算flat及三个rough案例、各5001/10001波束，共8个作业。

Bellhop专用程序从results/archive旧工作树快照复用，SHA-256与活动构建清单完全一致，两份overlay源码哈希也一致；没有访问cash，没有安装编译器或更换普通Bellhop算法。source_data/binary_origin.json保留来源。

共享984点接收网格与Stage-0 M99入射能量权重；主图使用10001波束。Bellhop内部壁剖面为[-80,80] m、4097点、步长.05 m。仅固定表示转换：SHD到native为conj(raw)，native粗糙/平面比再conj进入PE比较表示。无逐案例定相位、幅度拟合或重生成声源。

Bellhop有限发射角覆盖外，平面场为零的375个接收点无法定义粗糙/平面比，显式保存为复NaN并给出ratio_valid掩膜；原始SHD仍保留。它们均在固定M99之外，411个比较点全部有效，比较支持没有缩减；没有插值或补造曲线。

历史BIE数值变化指标仅引用原报告，约6e-8，未在本轮重新认证；不得把历史数值当作本轮新Bellhop的严格误差界。波束变化也只作敏感性。

本轮10001波束平面控制与保存Stage-0平面Bellhop原生场的固定M99复L2为0；记录在metrics/flat_control.csv。八个作业的进程实际耗时单独保存，不能把读取已完成输出的时间当作求解耗时。

## 本轮指标

|案例|配对|复L2|相位RMS / mrad|幅度比差值RMS / dB|
|---|---|---:|---:|---:|
|weak_low_K|PE_BIE|0.0031575715|3.1176357|0.0043483796|
|weak_low_K|BH_BIE|8.091975e-06|0.0058861658|4.8230197e-05|
|weak_low_K|PE_BH|0.0031585489|3.1186396|0.0043476181|
|region_II_low_K|PE_BIE|0.01608996|15.894495|0.021741798|
|region_II_low_K|BH_BIE|3.889551e-05|0.02791467|0.00023525879|
|region_II_low_K|PE_BH|0.016094693|15.899333|0.021739307|
|strong_height_low_K|PE_BIE|0.080908108|80.336353|0.086996387|
|strong_height_low_K|BH_BIE|9.9760075e-05|0.074110319|0.00057979319|
|strong_height_low_K|PE_BH|0.080922792|80.351139|0.087001781|

## 图册

### 共同环境与三个海面

![共同环境与三个海面](figures/01_environment.png)

同一声源和几何下，只改变海面起伏高度。

[矢量PDF](figures/01_environment.pdf)

左图明确Tx与Rx中心竖直对齐，接收线位于3 m深度，没有海底反射。右图比较同一波数、同一边缘处理下的三档高度。Bellhop用变换坐标中的内部壁表示相同海面，不使用Tx99实际水深展示环境。

### Bellhop–BIE 归一反射场

![Bellhop–BIE 归一反射场](figures/bellhop_bie/BB01_fields.png)

两种方法在相同接收坐标上的归一反射幅度和相位直接比较。

[矢量PDF](figures/bellhop_bie/BB01_fields.pdf)

每列对应弱、中、强海面；上排是|G|，下排是相对于平面反射的主值相位。幅度和相位采用跨列共同刻度，主值跨越±pi时留空，不用连线伪造相位连续性。比较中的细小差异应结合对应差值图判断，不能只凭曲线肉眼重合下结论。

本轮Bellhop与BIE的加权复差为8.092e-06至9.976e-05，两条主曲线因此在此尺度下几乎重合。强海面的主值相位多次跨越±pi，图中的断线表示相位绕回，不是接收声场消失。

### Bellhop–BIE 幅相差值

![Bellhop–BIE 幅相差值](figures/bellhop_bie/BB02_residuals.png)

差值图将主曲线难以辨认的幅度与相位差单独展开。

[矢量PDF](figures/bellhop_bie/BB02_residuals.pdf)

上排为20log10(|G_BH|/|G_BIE|)，正值表示Bellhop归一幅度更大；下排为arg(G_BH conj(G_BIE))，正值表示比较表示下Bellhop相位更大。各列独立纵轴，1 mrad=1e-3 rad；没有去掉任何全局相位。

weak_low_K：M99内最大绝对幅度比差值0.000166224 dB，最大绝对相对相位0.024704 mrad。数值对应差值本身，不能用|G|接近1替代误差判断。

region_II_low_K：M99内最大绝对幅度比差值0.00088818 dB，最大绝对相对相位0.134962 mrad。数值对应差值本身，不能用|G|接近1替代误差判断。

strong_height_low_K：M99内最大绝对幅度比差值0.0046381 dB，最大绝对相对相位0.606372 mrad。数值对应差值本身，不能用|G|接近1替代误差判断。

### Bellhop–BIE 指标

![Bellhop–BIE 指标](figures/bellhop_bie/BB03_metrics.png)

加权指标概括三个固定海面的差异，并保留波束敏感性。

[矢量PDF](figures/bellhop_bie/BB03_metrics.pdf)

从原始复数数组在固定M99上重新计算，权重为归一化入射能量。复L2以配对中的第二个方法为分母；相位RMS使用wrapped相位；dB项是幅度比差值，不是绝对TL。对数坐标用于显示不同数量级，没有改变原始数值。

### PE–Bellhop–BIE 归一反射场

![PE–Bellhop–BIE 归一反射场](figures/three_way/TR01_fields.png)

三条曲线使用同一环境与表示约定，不作幅相拟合。

[矢量PDF](figures/three_way/TR01_fields.pdf)

每列对应弱、中、强海面；上排是|G|，下排是相对于平面反射的主值相位。幅度和相位采用跨列共同刻度，主值跨越±pi时留空，不用连线伪造相位连续性。比较中的细小差异应结合对应差值图判断，不能只凭曲线肉眼重合下结论。

本轮Bellhop与BIE的加权复差为8.092e-06至9.976e-05，两条主曲线因此在此尺度下几乎重合。强海面的主值相位多次跨越±pi，图中的断线表示相位绕回，不是接收声场消失。

加入PE后，强海面的幅度空间起伏已可与Bellhop/BIE区分；三种方法的相位在主图尺度下仍看起来接近。PE–BIE复差由弱海面的0.003158增至强海面的0.08091，必须结合TR02与TR03识别，不能把主图相似解读为复数一致。

### 三组配对残差

![三组配对残差](figures/three_way/TR02_residuals.png)

分别观察 PE–BIE、Bellhop–BIE 和 PE–Bellhop 的空间残差。

[矢量PDF](figures/three_way/TR02_residuals.pdf)

上排为无量纲复数差的绝对值，下排为有符号相对相位。它们描述接收坐标上的残差结构，不是海面各点贡献，也不能解释为能量差或独立多径。不同配对使用相同M99坐标。

橙色PE–BIE与灰色PE–Bellhop残差在当前纵轴上接近重合，蓝色Bellhop–BIE位于零附近。这与BB02放大的微小蓝色残差一致；蓝线接近零不表示严格相等。残差在接收线边缘较大，指标仍采用原入射能量权重，未删去这些点。

### 三方指标总结

![三方指标总结](figures/three_way/TR03_metrics.png)

三个配对指标共同说明这组固定实验中的关系。

[矢量PDF](figures/three_way/TR03_metrics.pdf)

从原始复数数组在固定M99上重新计算，权重为归一化入射能量。复L2以配对中的第二个方法为分母；相位RMS使用wrapped相位；dB项是幅度比差值，不是绝对TL。对数坐标用于显示不同数量级，没有改变原始数值。

weak_low_K：PE–BIE复误差为0.00315757，Bellhop–BIE为8.09198e-06，两者比值为390.2。该比值只描述当前固定海面上相对于保存BIE参考的差异，不代表普适模型优劣。

region_II_low_K：PE–BIE复误差为0.01609，Bellhop–BIE为3.88955e-05，两者比值为413.7。该比值只描述当前固定海面上相对于保存BIE参考的差异，不代表普适模型优劣。

strong_height_low_K：PE–BIE复误差为0.0809081，Bellhop–BIE为9.97601e-05，两者比值为811。该比值只描述当前固定海面上相对于保存BIE参考的差异，不代表普适模型优劣。

## 历史核对与波束敏感性

|案例|配对|本轮复L2|历史报告复L2|
|---|---|---:|---:|
|weak_low_K|PE_BIE|0.0031575715|0.00315757|
|weak_low_K|BH_BIE|8.091975e-06|8.09198e-06|
|weak_low_K|PE_BH|0.0031585489|0.00315855|
|region_II_low_K|PE_BIE|0.01608996|0.01608996|
|region_II_low_K|BH_BIE|3.889551e-05|3.88955e-05|
|region_II_low_K|PE_BH|0.016094693|0.01609469|
|strong_height_low_K|PE_BIE|0.080908108|0.08090811|
|strong_height_low_K|BH_BIE|9.9760075e-05|9.97602e-05|
|strong_height_low_K|PE_BH|0.080922792|0.08092279|

历史数字引自source_data中的原报告；差异不通过调参数追平。

|案例|BH5001–BH10001复L2|相位RMS / mrad|
|---|---:|---:|
|weak_low_K|3.5113444e-06|0.0034758894|
|region_II_low_K|4.8939664e-06|0.0041957246|
|strong_height_low_K|1.1121732e-05|0.0048596737|

## 可复现数据

`plot_data/comparison_data.mat`包含原始复数G、native表示、坐标、海面、权重、掩膜和残差；逐案例CSV保留实部/虚部及全部绘图变换。source_data保存输入快照、哈希、Bellhop原始文件、专用程序和作业检查点。metrics分别保存本轮与历史指标。

在setup_vertical_project后运行：`generate_pe_bellhop_bie_comparison_visuals(run_dir)`。指定第二参数可重绘到新目录；无需Bellhop或源工作树路径。

复现示例（run_dir指向本数据包所在文件夹）：

```matlab
setup_vertical_project;
generate_pe_bellhop_bie_comparison_visuals(run_dir, fullfile(run_dir,'rebuilt_figures'));
audit_pe_bellhop_bie_visualization(run_dir, fullfile(run_dir,'rebuilt_figures'));
```

本轮仅新增验证/可视化流程，生产物理与公共接口未改；不运行通信或生产信道回归。图组仅适用于三个固定二维确定性海面，不外推PM、三维或通信性能。

工作流采用matlab与scientific-visualization技能。软件流程归因：Kassis T., Agarwal V., He Y., Patel D., Brueckner A.M. (2026), [Scientific Agent Skills](https://doi.org/10.48550/arXiv.2609.00065)，当前arXiv记录已核对；此引用不是声学证据。
