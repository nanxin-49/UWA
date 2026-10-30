# 三类发射声源最小归一化对照（2026-10-30）

## 范围与现有实现

本次仅比较 4 kHz、均匀自由场、向上前向半空间的归一化角度响应。没有运行或修改 PE、BIE、FullFK、Bellhop、粗糙海面或通信链。现有 `src/propagation/gaussian_source_initial_field_vertical.m` 给出发射平面二维横向初始场

`psi_0(x,y)=exp(-[(x-x_tx)^2+(y-y_tx)^2]/(2 sigma^2))`，峰值 1、无频率参数、复数虚部为零；它不是按换能器表面速度/声压边界条件定义的物理指向性。接口为 `[psi0_xy,meta]=gaussian_source_initial_field_vertical(X,Y,cfg)`，其中 `cfg` 需要 `x_tx,y_tx,sigma_src_m`。项目已有一横向维的离散 FFT 角谱用于 Gaussian sigma 灵敏度研究，见 `scripts/validation/validate_full_kirchhoff_gaussian_sigma_sensitivity.m` 和 `reports/pe_bie_full_kirchhoff_gaussian_sigma_sensitivity_report.md`；该报告的 theta99 是离散谱能量分位角，不是实际发射器的实测波束宽度。本实验用上述可复用入口数值检查一维切片，二维可分离高斯的另一个横向积分在轴向归一化后抵消。

项目 Bellhop 适配脚本 `scripts/validation/validate_pe_bellhop_incident_field_vertical.m` 的 `local_pattern` 使用 `cos(theta)*exp[-(k sigma sin(theta))^2/2]`，作为 Bellhop 方向图输入。该 `cos(theta)` 是额外的斜投影/角度约定因子，并非上述发射平面场的傅里叶变换；现有代码没有证明它对点源或活塞也成立。本次统一采用**纯归一化压力角度权重**，不把该因子加进 Gaussian 或其他两类。若日后做同源物理比较，需要先从相同辐射积分和源面边界条件推导或校准角谱到远场压力的映射，不能混用这些定义。

## 模型与计算

水深 1500 m，`z` 向下，Tx=1490 m，Rx=6 m，竖直距离 `L=1484 m`，水平偏移 `x=0:1:200 m`；`c=1500 m/s`，`f=4000 Hz`，`k=2πf/c=16.75516082 rad/m`，主轴竖直向上。采用 `exp(-iωt)` 与出射 `exp(+ikR)` 约定；本轮只画幅度，公共相位在轴向归一化后不影响结果。

角度 `theta=atan(x/L)`，`R=sqrt(L²+x²)`，最大角度 7.675555°。归一化角度压力权重为：

- Gaussian 初始场的二维 Fourier 角谱：`D_G=exp[-(k sigma sin(theta))²/2]`，`sigma=0.3 m`；这是初始场谱权重，不是完整物理换能器方向图。
- 理想全向点源：`D_O=1`。
- 理想均匀圆形活塞的远场方向图：`D_P=2 J_1(u)/u`，`u=k a sin(theta)`，`a=0.05 m`；`u=0` 单独取 1。

三者只在轴向压力幅度设为 1。接收对照统一取 `A_s(x)=D_s(theta)L/R`，即共同三维远场 `1/R` 幅度相对于 x=0 的值；`20log10(A_s)` 为各自轴上归一化结果，`20log10(A_s/A_O)=20log10(D_s)` 为相对全向的差值。此作图约定没有声源级、辐射功率、活塞近场或源面到远场转换系数的绝对一致性，也不把二维线源的 `1/sqrt(R)` 衰减误称为三维点源结果。

## 结果

| x (m) | theta (°) | Gaussian/轴上 | 全向/轴上 | 活塞/轴上 | Gaussian−全向 (dB) | 活塞−全向 (dB) |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 0 | 1 | 1 | 1 | 0 | 0 |
| 50 | 1.9297 | 0.98522 | 0.99943 | 0.99933 | -0.12442 | -0.000864 |
| 100 | 3.8551 | 0.94236 | 0.99774 | 0.99734 | -0.49601 | -0.003445 |
| 150 | 5.7717 | 0.87560 | 0.99493 | 0.99405 | -1.1097 | -0.007708 |
| 200 | 7.6756 | 0.79107 | 0.99104 | 0.98949 | -1.9575 | -0.013597 |

在此偏移范围内，Gaussian 谱权重影响最明显；5 cm 活塞在 4 kHz、`ka=0.83776`、最大角度约 7.7° 的范围内与全向声源非常接近。Gaussian 可以继续作为可控数值入射场及基线，但**不宜未经换能器参数、边界条件和声源级校准就作为后续物理发射器主模型**。本结论只针对本实验的方向响应与几何，不判断现有 PE/BIE/FullFK 验证资格。

若接入现有传播框架：Gaussian 要固定二维/三维横向初始场、谱归一化及相位/前向谱截断；全向点源需要匹配三维点源 Green 函数或二维线源 Green 函数，处理奇异源与源面注入，并定义同一参考幅相；圆形活塞需要以活塞表面速度/压力和孔径积分建立源面场，处理有限孔径、近远场和前向半空间边界，再转换到相同接口。三者还须统一频率依赖、功率或轴上标度与接收算子，才可称同源求解器验证。

## 检查与复现

入口：`setup_vertical_project; r=run_source_model_minimal_comparison('new_unique_run_id');`。共享 MATLAB R2025b Update 6，实际 run `run_20261030_02`，6.245 s，6 项检查通过、失败数 0：全向恒为 1，活塞轴上有限且为 1，共同几何，三者轴上同标度，有限角范围内幅度次序，Gaussian 数值傅里叶积分与解析角谱一致。该积分最大绝对差 `4.96e-11`（容差 `1e-8`）。从保存 CSV 独立重算 Gaussian 接收曲线最大绝对差 `5.55e-16`；重开 FIG 检出三条方向图曲线。无随机数。

原始数据、MAT、CSV、图 1/2 的 `.fig` 与 `.png`、stdout、审计及哈希检查点见 `results/comparisons/source_model_minimal/run_20261030_02/`。图 1 显示下限为 -30 dB，更低角度权重仍保存在 MAT 中；该截断不影响本场景的 0–7.7° 部分。`run_20261030_01` 是首次绘图检查版本，图例多了自动项；保留而不覆盖，最终使用 `run_20261030_02`。本次没有新增生产模块或默认参数；`vertical_comm_guide.md` 的物理/通信实现未变化，无需修订。

图形制作流程参考：Kassis, T., Agarwal, V., He, Y., Patel, D., & Brueckner, A. M. (2026). *Scientific Agent Skills: A Library of Procedural Knowledge for Research Agents*. arXiv:2609.00065. https://doi.org/10.48550/arXiv.2609.00065
