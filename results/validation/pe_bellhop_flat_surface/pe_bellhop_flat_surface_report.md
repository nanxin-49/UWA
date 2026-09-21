# PE 与 Bellhop 第一阶段平面海面交叉验证报告

## 结论

该精简案例通过全部自动检查。PE 与 Bellhop 均识别出直达和一次海面反射两条主要路径；到达时间、相对路径强度及 PDP 峰位一致。绝对幅度存在近似常数偏差，主要来自当前 PE 高斯初场与 Bellhop 单位点源归一化不同，而不是新增的传播路径或粗糙面效应。

## 共同环境与运行设置

- 水深：`100 m`；均匀声速：`1500 m/s`。
- 发射位置：`[0, 0, 80] m`；接收位置：`[5, 0, 10] m`。
- 中心频率：`4000 Hz`；PE 宽带：`3000--5000 Hz`，`65` 个等间隔频点。
- 海面：平面压力释放边界，`surface_reflect_coeff=-1`、`sea_hs_target=0`；未启用 SSA、随机海面、Kirchhoff 粗糙相位或气泡。
- PE 网格：`128 x 128`，横向窗口 `32 m x 32 m`，`stepz_lamb=0.5`。
- Bellhop：标准 ASCII arrivals 模式，`5001` 条射线，只覆盖直达/一次海面反射所需角扇区。
- Bellhop 依赖：外部 Acoustic Toolbox `bellhop.exe`；本次使用官方托管的 [Windows atWin.zip](https://oalib-acoustics.org/website_resources/AcousticsToolbox/versions/atWin.zip) 分发包，第三方二进制未复制进项目。
- 几何路径长度：直达 `70.1783 m`，海面镜像 `90.1388 m`，理论超时延 `13.307 ms`。

## 到达与幅度结果

| 路径 | PE 到达 (ms) | Bellhop 到达 (ms) | 差值 (ms) | PE TL (dB) | Bellhop TL (dB) | 直达校准后 TL 差 (dB) |
|---|---:|---:|---:|---:|---:|---:|
| direct | 46.785563 | 46.785563 | -0.000000 | 33.9582 | 36.9241 | +0.0000 |
| surface_reflection | 60.016332 | 60.092516 | -0.076184 | 35.7698 | 39.0982 | -0.3625 |

PE 原始幅度相对 Bellhop 的直达标定系数为 `0.71073496`。标定只用于区分源归一化常数与路径相对衰减，不回写 PE 公共输出，也不用于掩盖路径间差异。

## 自动检查

| 检查 | 数值 | 条件 | 阈值 | 通过 |
|---|---:|:---:|---:|:---:|
| path_count | 2 | == | 2 | 1 |
| arrival_time | 0.076184144 | <= | 0.75 | 1 |
| raw_tl | 3.3283895 | <= | 6 | 1 |
| relative_path_tl | 0.36254308 | <= | 3 | 1 |
| direct_calibrated_reflection_tl | 0.36254308 | <= | 3 | 1 |
| pdp_peak_time | 0.061538462 | <= | 0.75 | 1 |
| pe_frequency_sum_invariant | 1.7347235e-18 | <= | 1e-10 | 1 |
| pe_direct_only_invariant | 0 | <= | 1e-10 | 1 |
| pe_1_over_r | 0 | == | 0 | 1 |

![PE 与 Bellhop PDP/TL 对比](pe_bellhop_pdp_comparison.png)

## 差异解释与限制

1. 当前 PE 以有限宽度高斯场启动，Bellhop arrivals 使用单位点源的几何扩展归一化；因此原始 TL 可出现近似常数偏置。更有辨识力的是直达校准后的反射残差及反射/直达相对 TL。
2. PE 保存的是去载波复包络。本验证仅在比较层按直达长度和镜像路径长度恢复`exp(-i 2 pi f tau)`，并从分路径 CIR 峰位提取到达时间；没有修改 `vertical_channel_model` 的公共输出语义。
3. Bellhop 为射线/高斯波束模型，PE 为波动模型；有限频带、有限网格、波束插值和绕射会产生小的幅相差异，不要求逐点完全一致。
4. 为避开 Bellhop 在零水平距离的退化几何，接收机设置 5 m 水平偏移；该传播仍为近垂直上行。阶段一不包含海底反射、粗糙海面、SSA、随机信道或通信调制。

5. 总 PDP 峰间的 PE 旁瓣来自有限带宽重构和 PE 包络随频率的幅相变化；主要路径数由分离的直达/反射分量和 Bellhop bounce count 判定，不能把旁瓣计为额外本征声线。

## 生成文件

- `flat_surface_case.env` / `.arr` / `.prt`：Bellhop 输入及原始输出。
- `arrival_time_comparison.csv`：路径对比表。
- `pe_bellhop_pdp_comparison.png`：PDP 和路径 TL 图。
- `pe_bellhop_flat_surface_validation.mat`：完整可复查结果。

运行入口：`validate_pe_bellhop_flat_surface_vertical`。
