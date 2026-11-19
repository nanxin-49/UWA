# K1目标执行与保存收尾验证

2026-11-19。结论：**限定科学结果/当前新图件PASS，远端退出环境仍有原生异常，旧FAIL保持。**

|检查|实际结果|
|---|---|
|冻结输入准备|37断言，12.8416714s，8个声学入口哨兵，5余例dry-run，0新求解；MODELS剥离日志后的原数学一致。|
|并发编排mock|5场景PASS：四并发完成及零重复、N40 MODELS失败停止派发、并行案例失败停止新派发、RAM只停已验证本任务进程、汇总失败保留案例。未模拟本次封存后原生退出，真实出现的缺口和KeyError保留。|
|实际基础MATLAB许可|4个额外独立进程实际重叠运行并正常退出；无PCT/parpool，不推断许可证总席位。|
|传输/SHA|恢复04失败1575项、首批3876项、最终4040项内部SHA通过。最终tar SHA 8210f7b5bc50e09b8557c9a84edf935d877788223e9404034dc79eed939b38a1。|
|真实原始层级|六新海面各PS六层、BIE五层、完整FullFK四层；导入N7 5个BIE标记+本目标25个新标记=30，每层invocations=1。所有MODELS/C0–C4 seal签名、数值门及独立边界审计保持。|
|本地保存审查|N7复用既有数值资格；U6/N40、U5/N12、U8/N12、U5/N40、U8/N40分别79.377、50.999、50.067、78.880、79.265s，全部0失败/0新求解。余五面临时单计算线程、九入口哨兵、180s及30%本地可用RAM门，N7复用其原哨兵审查记录，结束恢复线程/路径。|
|最终新图件|六例各12曲线/8FIG，两次独立重开；7条件汇总56曲线/8FIG。合计128曲线、56FIG、112次零像素比较，像素容差0，显式v3布局审计。8个汇总PNG及案例各类型预览已实际查看。|
|真实七条件汇总|73.016886s，14模型行、5754逐点行、0未分辨，local-closeout模式，不是TEST_ONLY；原锚点/各自BIE/相位掩码/误差分解/MAT-CSV均独立检查。|
|Code Analyzer|5个当前MATLAB文件无error/warning；汇总入口仅1个info（find可用逻辑索引提高性能），其余4个0消息。|
|生产及阶段边界|git diff对src与四根公开入口无变更；无通信回归必要。未扩Seed/频率/宽带/Bellhop，未提交推送；最后只读确认五个本任务案例进程会话全部退出，未发送信号。|

## 明确保留的失败

原05 N40 MODELS原生abort、恢复04声学前像素FAIL、Windows旧布局失败、两进程大画布当时PASS但后续跨会话FAIL、小画布v1 U5/N12归一化图112750通道差、v2图例边界过近均保留。v3固定1000×550整幅PaperPosition print和内移图例属于显式图件版本更新，不把原1800×1400布局FAIL改成通过。所有4种误差和12条数据曲线不变，未修改原声学数组或数值容差。

一次MCP批量等待在300s超时，随后检查每例及N7的独立图件receipt均已完成，没有重算数值。U8/N40完整CASE_COMPLETE/最终seal后约5秒发生原生Segmentation violation，MATLAB−9；控制器旧错误分类继而KeyError，原exit78保持。只补充只读封存验证、精确本任务PID脱敏崩溃证据和独立保存审查，不重启。科学结果通过保存审查后接受，不认定MATLAB退出问题已经解决。

## 代码及可重现入口

当前科学后处理快照及SHA在results/validation/k1_target_execution/20261118_01/closeout_code_final/：

- reporting/render_k1_recovery_saved：每案例固定单面板输出；validation/audit_k1_recovery_figures：独立数组、显式布局与两个独立重开零像素检查。
- validation/review_k1_partial_saved：保留原两参数入口，追加可选case_id和numeric模式；只审查保存数据，哨兵禁止新求解。
- reporting/summarize_k1_recovery_cases：追加local-closeout及独立资格清单参数；原MODELS/C0–C4链照常验证，另核本地数值/新图件SHA，不伪造失败旧FINAL。
- reporting/render_k1_pm_summary_fields：保存七条件原复场、8个独立面板、56条曲线和两次独立重开零像素检查。

远端准备阶段沿用现有recovery preflight/controller并增加本任务四并发分支，未新建用户上传/启动脚本。真实最后例通过SSH内联命令调用原冻结worker，命令完整保存在last_case_command.json；退出分类缺陷由独立补充证据记录，不回写旧controller结果。生产PM、PS、BIE、FullFK、DSP/通信代码及公共输出接口不改。

最终本地科学资格根：results/validation/k1_target_execution/20261118_01/；主研究报告：[K1最终报告](pm_wind_cutoff_k1_final_review_2026-11-19.md)。最终SHA、检查点和图件视觉回执应与旧远端exit分别阅读。没有待用户手动启动的K1计算。
