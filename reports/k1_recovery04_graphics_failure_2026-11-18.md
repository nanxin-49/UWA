# K1恢复04首跑：保存图件像素门失败

2026-11-18。**STOPPED_BEFORE_NEW_ACOUSTICS**。只读诊断会话`20261118T143953_1551655`，没有重新渲染、启动MATLAB作业、发信号或写远端文件。原04包/旧失败与原N7资格均不改。

预检MATLAB退出0；下一步saved_graphics退出1，监督器退出78、65.945894449s，`stop_reason=null`，新案例派发0、新BIE0。N40没有开始，原N40累计18.833152669s不变。累计账为724＋ceil(65.9459)=**790s**，后续必须继承，不能重置。

具体断言为`audit_r2_final_figures:58`的`r2figure:Pixels`：F02第一次重导有**399701个像素通道值不一致**。该图8条误差曲线已与保存场逐条一致，画布/轴布局门没有在这一遍触发失败。`pixel_differences=[399701,0]`的第二个0是预分配初值，因为第一次断言已终止，不能当作第二次PASS；F04审计也未执行。

只读下载了已生成的两张PNG，独立核SHA，未修改图片或重新绘图：

- 初次导出SHA `3b603b1aa386b5ac23d5942ce154f44486d329ee4188f65edeb08b3aeaf581a6`。
- FIG重开重导SHA `4ae18eed062edf7c45f2401547b876748726102f0f550a2ebcbcec74e4c5231f`。

两图目视可见图例框尺寸和坐标文字/曲线栅格呈现差异。审计记录FIG刚重开后尺寸为1024×731，随后原审计器要求恢复1800×1400；与图窗屏幕适配/序列化后重排有关的可能性需要后续FIG属性对照验证，当前不把观察当成已定位的唯一原因。没有原生abort或资源门停止证据，此次是明确的MATLAB审计assert退出，和此前N40的-9/allocator corruption分开。

前次Ubuntu保存图件诊断的零像素PASS仍是那一次运行的真实结果，不能追改；本次表明同一导出链跨会话仍不稳定，04包不能因此继续放行。数值曲线检查通过也不能代替原零像素门。不得跳过saved_graphics、调大像素容限或直接反复Start。

证据位于`results/validation/k1_recovery_failure_review/20261118_01/`：[只读日志及SHA](../results/validation/k1_recovery_failure_review/20261118_01/stdout.json)、[两张PNG传输SHA](../results/validation/k1_recovery_failure_review/20261118_01/image_transport.json)。本地图像读取脚本初用`inspect.py`命名，与Python标准库冲突；已更名并改为标准库纯传输，期间仅重复只读查询，没有额外计算或远端修改。该辅助问题与远端图件失败独立。

下一步先由用户对当前04任务执行Fetch保留全目录，再仅修复保存图件导出/重开链。修复需同时验证原12条曲线、1800×1400画布/轴/图例、两次重开零像素差，并跨独立MATLAB进程复核稳定性；不能只以同一会话单次PASS代替。所有额外原生诊断均应明确有界并计入790s之后的累计预算。N7数组、原PS/BIE/FullFK和科学数值门保持，不重算任何海面；具体修复和新手动包尚未执行。

```powershell
& $tool -Action Fetch -PackageDir $package
```

此处`$package`和`$tool`仍指恢复04，不能换回原05。回收后再准备保留本次失败的新图件修复包，不能在04原目录删改失败记录。本轮未提交推送Git。
