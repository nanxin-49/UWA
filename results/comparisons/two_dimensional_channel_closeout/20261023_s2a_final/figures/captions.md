# S2a 图注与实际检查

2026-10-23。三张图只画 z=3 m 原 M99 内 411 个真实接收线样本。纵轴为二维 Gaussian 单位源响应，不是 Pa 或绝对 TL；横轴 x 为米。PE/BIE/FullFK 均在同一原生 `exp(-i omega t)` 表示。蓝色 PE、红色 BIE、绿色 FullFK；F05 蓝色为 PE–BIE、绿色为 FullFK–BIE。每张图有同名可编辑 `.fig` 与 PNG，原始复数在 `../data/line_fields.csv`，无插值/平滑或新 x-z 网格。

- **F02 平面控制**：上排 `|u_flat|` 和主值相位，下排相对 BIE 的原生复残差与 wrapped 相位差。相位跨 ±π 的连线留空；PE 与 BIE 基线重合不表示 FullFK 的误差上界为零。实际检查：三方法平面幅度同源标度可读；小残差轴独立注明数量级；相位绕回未连续误连。
- **F04 三档正弦原生反射场**：各行弱、中、强正弦；左 `|u_rough|`，右原生主值相位。三方法复数来自真实接收线；绿色强档沿用 E12 保存积分，蓝色 PE 经实际适配入口核验后复用 v2。实际检查：三档和各方法均显示；共同幅度轴约 0.02–0.14，原生相位有跳变留白；没有把大尺度幅度贴合说成复数相等。
- **F05 残差与辅助量**：逐行左为同源原生 `|u_rough−u_BIE|`，中为各方法自身平面分母所成 `|G−G_BIE|`，右为 `|δu−δu_BIE|`。三列单位与对象不同，分别标在标题/纵轴；中心图不能替代左侧主量。实际检查：绿色 FullFK 小差仍可辨；蓝色 PE 边缘值未裁去；三档共用同一 x 支持。

从 `data/plot_data.mat` 调用 `render_two_dimensional_closeout_s2a(plot_file,new_dir)` 独立重绘。`audit/independent_audit.json` 记录三张 PNG 逐像素一致、三张 `.fig` 重开成功。PDF 未要求，未额外生成。
