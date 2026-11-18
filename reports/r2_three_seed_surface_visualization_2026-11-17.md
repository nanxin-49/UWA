# R1/R2三个Seed的海面几何补图

2026-11-17。只从已封存的seed260001/260002/260003 PM系数和原C2嵌入重建几何，不抽样新海面、不归一化振幅、不加载求解入口、不做声学计算或SSH。

![三种子海面对照](../results/visualization/r2_three_seed_surfaces/20261117_01/F06_three_seed_surfaces.png)

[可编辑FIG](../results/visualization/r2_three_seed_surfaces/20261117_01/F06_three_seed_surfaces.fig) · [几何CSV](../results/visualization/r2_three_seed_surfaces/20261117_01/surface_profiles.csv) · [审计回执](../results/visualization/r2_three_seed_surfaces/20261117_01/receipt.json)

上图叠加三种子，下三图各自展示，使用相同坐标范围和纵轴尺度。绘制的是实际用于R1/R2比较的C2嵌入海面；海面高程采用`eta=-z_s`，正值向上，原物理深度仍正向下。灰色区域为42<|x|<50 m的C2过渡带；|x|≤42 m保留原Fourier实现，|x|≥50 m归零。纵向为了辨认作了显示放大，不能将图中坡度外观当作等比例几何。

底层160 m周期、12模、PM谱、Fourier系数均不改变。显示区间x∈[−60,60] m，显示采样0.01 m，只用于解析几何绘图，不构成新声场网格或参考资格。411接收线的x支持仍为−40.0390625至40.0390625 m。

| Seed | 显示采样最小/最大高程 (m) | 接收线水平支持内高程RMS (m) |
|---|---:|---:|
| 260001 | −0.4994565 / 0.4092995 | 0.2114029 |
| 260002 | −0.3033853 / 0.3235664 | 0.1575142 |
| 260003 | −0.3817594 / 0.4306613 | 0.2217070 |

统计为所声明显示采样的几何描述，不是连续极值保证、显著性检验或RMS重归一化目标。不同实现的几何幅值自然不同，未为了可比性缩放曲线。

新入口`scripts/reporting/render_r2_three_seed_surfaces.m`可在`setup_vertical_project`后接受显式绝对的新输出目录，从任意MATLAB当前目录运行；输出位于`results/visualization/r2_three_seed_surfaces/20261117_01/`。代码校验原R2总封存SHA、三份common.mat来源和PM谱/k一致性；在原保存几何采样位置核对差<1e-12 m。独立逐模Fourier/C2重建最大差2.77556e-16 m，FIG六条曲线逐点一致，两次FIG重开PNG导出通道差均0。现有MATLAB MCP R2025b Update6运行28.207272 s，PNG已目视查看，标题/单位/图例/过渡带和高程符号一致。保留PNG/FIG、MAT/CSV、来源SHA及工具回执，未输出PDF。

本轮仅新增几何可视化入口、补图/审计与导航。原R2封存目录和[最终报告](r2_saved_closeout_final_2026-11-17.md)不修改，数值结论、研究阶段、接口与生产实现不变，因此不更新PROJECT_CONTEXT/项目进度JSON/通信方法指南、不运行通信回归；仅在任务表记录补图交付。新声学求解0、SSH0、Git提交推送0，不进入R3或宽带。
