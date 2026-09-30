---
name: obsidian-visualization-report
description: Package an existing technical visualization report and its referenced figures into a self-contained Obsidian-ready folder. Use when the user asks to prepare, copy, or organize a report with figures for an Obsidian vault; do not use to generate new scientific results or figures unless separately requested.
---

# Obsidian visualization report packaging

把已有的分析报告及其引用图片整理成可直接放入 Obsidian 的独立目录。默认只整理和转换报告副本，不改源报告、不重跑计算、不重绘图。

## 开始前

1. 确认要打包的权威报告文件。若存在多个版本，依据用户指定、项目当前报告索引或明确的 superseded 标记确认；不把历史/废弃报告当作当前结论。
2. 确认 Obsidian vault 的绝对目标路径及子文件夹名。用户本轮未提供时，先询问并暂停外部写入。不要把曾经给过的 vault 路径视作对其他任务永久有效的授权。
3. 检查报告引用的图文件是否在项目中存在。只收集报告实际引用且属于本报告的图，不遍历归档目录或把整个 results 树复制过去。
4. 先检查目标目录与同名文件。不得静默覆盖用户文件；若目标已有内容，先选择无冲突的新子目录或询问用户。

## 整理步骤

1. 创建独立的报告目录，目录内放报告副本和 `图片/` 子目录。
2. 将权威报告复制为该目录中的 Markdown 文件。保留分析结论、数值、表格和方法说明；不擅自改写技术内容。
3. 从报告中逐项提取本报告引用的图片路径，解析仓库相对路径及报告中明确指向的图路径，去重后复制到 `图片/`。
4. 在报告副本中把图片引用统一转换为 Obsidian wiki embed，并使用相对路径，例如：

   `![[图片/weak_lowK_incident_xz.png]]`

5. 图片文件名保持原样。文件名中的下划线不转义：`strong_height_low_K_ablations.png` 正确；不要写成 `strong\_height\_low\_K...`。不要把 embed 写成 `!图片/...`、普通 Markdown `![...](...)` 或缺少双方 `[[ ]]` 的形式。
6. 删除或改写报告副本中会失效的本机绝对路径链接，使其指向本目录相对图片或项目相对源文件；不改源报告。不要保留机器专属路径作为图片链接。
7. 除非用户明确要求，不生成新图片、不修改图像内容、不复制 CSV/MAT/源码，也不改变任何验证结论。必要的来源说明可保留为项目相对路径或普通文字。

## 完成前验证

- 报告副本中的每个 `![[图片/文件名]]` 都对应 `图片/文件名` 中实际存在的文件。
- 被收集的报告图片均已复制；引用和文件名逐字一致，包括大小写、空格和下划线。
- 图片 embed 使用相对路径，格式为 `![[图片/filename.png]]`，没有多余反斜杠、代码反引号或路径转义。
- 报告副本没有遗留指向原机器绝对路径的图片引用；源报告未被修改。
- 目标目录只包含该报告及其需要的附件，且未覆盖已有文件。

最后向用户给出目标目录、报告副本路径、图片数量和验证结果；若有缺失引用，逐项列出并说明未伪造或替代该图片。
