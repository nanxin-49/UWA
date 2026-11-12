# 二维活塞自由场 quicklook

计算2026-11-11，交付封存2026-11-12；`20261111_2d_strip_01`。**仅图示完成，P0未执行。**

历史5 cm已确认是**圆活塞半径**：本地与[GitHub固定提交f3314fd](https://github.com/nanxin-49/UWA/blob/f3314fded85993832771f48e398f2a79edee582a/scripts/comparisons/run_piston_frequency_minimal_comparison.m#L12)均写`radius_m=0.05`。因此另算总宽10 cm条带作为半尺寸映射，并保留题目字面的总宽5 cm对照。圆形与条带的辐射模型不等价，未运行旧默认工况。

本次为y向无限的刚性障板匀速条带，源(0,100 m)、向上，c=1500 m/s、rho=1025 kg/m³，同峰值振速0.001 m/s，三频4/6/8 kHz，无海面/海底。使用完整二维速度孔径积分：

`p=(rho*omega*v0/2)*integral[-W/2,W/2] H0^(1)(k*hypot(x-xi,100-z)) dxi`，时间相量`exp(-i*omega*t)`；近场不用远场近似。绝对值保留Pa，归一化另列；差值为`20log10(|p20|/|pW|)`。

| 频率 | 轴上峰值Pa：总宽20/10/5 cm | 20 cm −3 dB全宽 |
|---|---|---|
| 4 kHz | 50.985 / 25.493 / 12.746 | 112.03° |
| 6 kHz | 62.444 / 31.222 / 15.611 | 67.12° |
| 8 kHz | 72.104 / 36.052 / 18.026 | 48.99° |

20 cm轴向约比总宽5 cm强12.04 dB、比历史映射10 cm强6.02 dB，同时波束收窄，8 kHz最明显，远场第一sinc零点约±69.64°；大离轴并非更强。8 kHz接收线x=120 m处20/10/5 cm为11.194/21.155/13.302 Pa。

20 cm值得作**待P0资格的主候选**，小孔径保留作控制。**下一步先P0，再P1平面海面**：本轮未做压力梯度、生产共源/相位、Gaussian回归或反射。没有正式PE/BIE/FullFK比较、SSH、CIR或通信研究，当前停止。

交付版为[figures_delivered/](figures_delivered/)内8张PNG及同名FIG；根目录三个`field_*kHz.mat`保存复Pa、整体场/近场、接收线和等距圆弧；三个`receiver_*kHz.csv`含复声压、绝对峰值Pa、各自轴向归一化与比值。详见[汇总指标](summary_metrics.csv)、[独立数值审计](audit_metrics.csv)、[完整方法/局限报告](../../../../reports/piston_source_quicklook_2026-11-11.md)。旧试绘版保留，勿作为交付版。

现有本机MCP R2025b Update6/PID28836：数值计算6.418 s、独立数值审计1.123 s，9组检查0失败；积分/完整角谱最大复相对差1.01e-14，保存场加密求积差9.30e-11。FIG/MAT重建与逐图记录见`figure_audit_delivered/`、`visual_review.md`；SHA与限制见`checkpoint_final.json`。

局限：无限条带/障板、理想给定振速，非等功率、非真实圆形换能器；省略源前1 cm及振面极限。整体幅度展示网格0.5 m不解析相位波长，不能直接拿作生产初场。quicklook检查不构成活塞数值资格。

![接收线：绝对Pa、轴向归一化与孔径比值](figures_delivered/receiver_line_comparison.png)

![方向图：三频同半径归一化，线为完整积分、圆点为远场sinc](figures_delivered/directivity_comparison.png)
