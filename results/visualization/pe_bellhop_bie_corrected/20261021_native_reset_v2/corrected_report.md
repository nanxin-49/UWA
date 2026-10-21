# PE–Bellhop–BIE 高度符号修正后的可视化

本包与原图册分开保存。固定4 kHz、c=1500 m/s、Tx深100 m、Rx深3 m、Gaussian σ=0.3 m、三个正弦海面A=0.01/0.05/0.20 m，K=0.10 rad/m。海面在|x|≤42 m保持原形，50 m以外归零。

海面几何采用z=η（深度向下），时间约定exp(−iωt)。PE相位屏由旧的exp(+2ikη)改为exp(−2ikη)，并从同一保存入射场FFT传播到接收线与真实x–z网格。BIE、Bellhop沿用保存原生复场；三方G都用原生表示。没有重新传播声源、求解BIE或运行Bellhop，也没有改变生产API。

## 先看直观图

### R01_height_sign_reset

![同一接收线的旧PE、修正PE与原生BIE；上排是|G|−1，下排是两种各自配对的相位差。](figures/R01_height_sign_reset.png)

同一接收线的旧PE、修正PE与原生BIE；上排是|G|−1，下排是两种各自配对的相位差。

### three_way/TR01_fields

![三方法直接比较G的幅值和主值相位；修正PE幅值调制方向与BIE一致。](figures/three_way/TR01_fields.png)

三方法直接比较G的幅值和主值相位；修正PE幅值调制方向与BIE一致。

### three_way/TR02_residuals

![三组配对的点位复差和有符号相位差；展示差异出现的位置，纵轴各案例独立。](figures/three_way/TR02_residuals.png)

三组配对的点位复差和有符号相位差；展示差异出现的位置，纵轴各案例独立。

### three_way/TR03_metrics

![固定M99权重下的复L2、相位RMS及dB幅值比RMS；曲线趋势限于这三个海面。](figures/three_way/TR03_metrics.png)

固定M99权重下的复L2、相位RMS及dB幅值比RMS；曲线趋势限于这三个海面。

### XZ01_raw_magnitude

![强海面真实x–z原始反射场幅值；PE/BIE同色标，Bellhop独立源标度。](figures/XZ01_raw_magnitude.png)

强海面真实x–z原始反射场幅值；PE/BIE同色标，Bellhop独立源标度。

### XZ02_raw_phase

![同一x–z复场的主值相位；色环首尾相接，色差不能单独当作误差。](figures/XZ02_raw_phase.png)

同一x–z复场的主值相位；色环首尾相接，色差不能单独当作误差。

### XZ03_complex_field

![色相表示相位、明暗表示幅值；暗区相位虽存储，但不宜凭颜色解读。](figures/XZ03_complex_field.png)

色相表示相位、明暗表示幅值；暗区相位虽存储，但不宜凭颜色解读。

### XZ04_raw_profiles

![z=1、3、6 m原始幅值剖面；上排PE/BIE共标度，下排Bellhop原生单位。](figures/XZ04_raw_profiles.png)

z=1、3、6 m原始幅值剖面；上排PE/BIE共标度，下排Bellhop原生单位。

Bellhop–BIE子图BB01–BB03和共同环境图01仍在figures/目录；该配对同时转为原生表示，其复差数值不变。彩图画的是原始u(x,z)，其Bellhop源标度与Gaussian峰值尚未标定，不能把三列像素亮度直接当作跨方法误差。接收线G=粗糙反射场/平面反射场，是无量纲比值，不是绝对声压或传播损失。

## 固定接收线指标（本轮与旧值分列）

|海面|配对|本轮复L2|本轮幅值相对L2|本轮相位RMS (mrad)|本轮dB幅值比RMS|旧复L2|
|---|---|---:|---:|---:|---:|---:|
|weak_low_K|PE_BIE|0.0031155506|5.1029695e-05|3.1151365|0.00044323865|0.0031575715|
|weak_low_K|BH_BIE|8.091975e-06|5.5527397e-06|0.0058861658|4.8230197e-05|8.091975e-06|
|weak_low_K|PE_BH|0.0031165556|5.1213159e-05|3.1161386|0.00044483219|0.0031585489|
|region_II_low_K|PE_BIE|0.015586811|0.00025556827|15.585204|0.0022198764|0.01608996|
|region_II_low_K|BH_BIE|3.889551e-05|2.7085074e-05|0.02791467|0.00023525879|3.889551e-05|
|region_II_low_K|PE_BH|0.01559176|0.00025642184|15.590142|0.0022272949|0.016094693|
|strong_height_low_K|PE_BIE|0.062911245|0.0010485013|62.935961|0.0091096633|0.080908108|
|strong_height_low_K|BH_BIE|9.9760075e-05|6.6770078e-05|0.074110319|0.00057979319|9.9760075e-05|
|strong_height_low_K|PE_BH|0.062930454|0.0010484603|62.955235|0.0091093333|0.080922792|

复差精确拆分：E² = 加权(|A|−|B|)²/加权|B|² + 加权[2|A||B|(1−cos Δφ)]/加权|B|²。强海面PE–BIE的两项为1.0993551e-06与0.0039567254，相位项占99.97222%；旧8.0908%包含先前高度/表示失配，不能继续当作修正后的模型误差。

修正后强海面PE–BIE：复L2 6.29112%，幅值相对L2 0.10485%，相位RMS 62.93596 mrad。仍有残余相位差；本次不把它归因给单一机制，也不宣称模型普遍精度。Bellhop–BIE复L2 0.009976%。

## 检查与复现

平面PE分母复现相对L2 1.82e-14；强海面二维z=3 m截线闭合 1.04e-16。每种方法的x–z实数/虚数保存在plot_data/xz_data.mat和xz_fields.csv；三案例G、掩膜、权重和旧PE保存在comparison_data.mat与逐案例CSV。source_data/source_hashes.csv记录三个原始MAT的路径及SHA256；旧图册、原BIE/Bellhop输出保持不变。

在MATLAB执行 `setup_vertical_project; generate_pe_bellhop_bie_corrected_visuals('本包绝对目录','新图目录');`，只读取本包MAT和随包保存的绘图源码；不运行声学求解。生成12对2400×1350 PNG/PDF。

适用边界：固定单频、共同M99接收区、三档确定性海面；BIE为历史派生副本，没有本轮收敛复核；5001/10001波束差是敏感性而非误差上界；二维彩图仅强海面、近海面z=0.5–6 m。
