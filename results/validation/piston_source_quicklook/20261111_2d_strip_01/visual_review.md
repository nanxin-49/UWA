# PNG实际查看与图件数据重建记录

2026-11-12交付封存。交付文件唯一目录：`figures_delivered/`。全部8张PNG已实际查看。

| PNG | 查看结论 |
|---|---|
| field_4kHz_20vs05cm | 总宽20/5 cm、整体与近场、相同绝对Pa色标；标题无遮挡 |
| field_4kHz_20vs10cm | 10 cm为历史半尺寸映射；坐标/色标与5 cm版一致 |
| field_6kHz_20vs05cm | 相比4 kHz束形收窄，近场放大与绝对幅度明确 |
| field_6kHz_20vs10cm | 同一总宽20 cm场复用，不重新归一化，两种幅度与比值标注清楚 |
| field_8kHz_20vs05cm | 宽角域双侧深谷可见；差值保留负值；没有把深谷抹平 |
| field_8kHz_20vs10cm | 总宽20/10 cm的宽角域增益可变负，原始复数仍保留 |
| receiver_line_comparison | 首排峰值Pa；第二排各自除轴向值；第三排比值含负离轴增益；三频共限，|x|=40 m有灰线 |
| directivity_comparison | 等半径97 m，三频三总宽，线型/颜色双重编码；积分线与解析sinc圆点一致；8 kHz零点谷底保留 |

绝对场色标固定[0,70] dB re 1 Pa peak，差值固定[-24,24] dB，均仅显示截断；未修改底层数据。全局/近场的横纵物理尺度保持1:1，放大图独立标坐标；源前s<1 cm省略。该图示采用一般科学报告格式，没有声称期刊或无障碍认证。

独立重开6张FIG，由MAT重新计算36幅CData并逐一匹配，最大差0 dB；公共CLim与z向下检查通过。24条接收曲线逐点匹配，8 kHz负离轴比值完整位于显示范围；代表FIG重新导出位于`figure_audit_delivered/`。方向图数值由独立`audit_metrics.csv`中等半径场/解析sinc比较支持。

保留的试绘问题：根目录首版图有总标题/小标题重叠，接收线纵轴未覆盖8 kHz负比值；`figures_final/`与`figures/`为后续试绘版本，仍有场图标题拥挤。一次uipanel布局因exportgraphics多容器限制失败，异常/部分FIG/源码留在`render_failed_panel/`与`source_snapshot_panel_failed/`；最终改为同一tiledlayout内专用标题行。这里只重绘保存数据，没有重算声学场、改门或删除失败记录。
