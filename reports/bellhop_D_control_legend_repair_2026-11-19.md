# 直达控制历史PNG图例修复与图件审计

2026-11-19。仅复用保存的201/401/801数据，Bellhop/FFP/新声束计算0。

另存三组PNG/FIG：原201诊断图、三档复数误差图、角间隔与二阶预测图。手动固定图例位置，显式绑定数据句柄，关闭AutoUpdate，去除data1。所有原PNG/FIG保持；6项原图件SHA在修复前后相同。

24条曲线从新FIG重开后与保存数值逐项精确一致；4个图例都在画布内，无多余data条目。三张实际PNG均已查看，裁切问题修复。

新图及审计：[输出目录](../results/visualization/bellhop_D_control_legend_repair/20261119_01/)、[数据审计](../results/visualization/bellhop_D_control_legend_repair/20261119_01/repair_receipt.json)、[目视审计](../results/visualization/bellhop_D_control_legend_repair/20261119_01/visual_audit.json)。脚本：scripts/reporting/repair_bellhop_D_control_legends.m。

本地MCP R2025b Update6函数耗时14.332s，返回图件审计通过；同次出现许可证警告，原文保存于结果目录。未关闭、清理或重新启动共享MATLAB，不用图件修复升级原201/401失败或追改旧801目视待确认回执。

后续仅由本轮明确授权的独立S任务执行，图件修复本身没有声学数值变化，不需要生产通信回归。