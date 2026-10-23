# S2b执行记录

- 用户授权创建并执行一个S2b Goal，未设置token预算。MATLAB已有R2025b Update6会话，经MCP evaluate运行，始终未使用project_path、cd或run_matlab_file改变cwd；未shell启动MATLAB。未提交推送或读取cash。
- 先核验AGENTS/PROJECT_CONTEXT顶索引/主计划/唯一主报告和S2a本地包。开始时git status仅主计划已有修改；保留该修改，仅追加执行状态说明。
- P0实际命令：setup_vertical_project; r2b=run_two_dimensional_closeout_s2b_p0('20261023_s2b_1305'); 24行双权重/4932复场CSV及MAT核对，9个人工公式控制通过，无物理场计算。
- P1先从Windows Get-Counter Memory Available Bytes取得1,351,811,072 bytes，再prepare_two_dimensional_closeout_s2b(r2b.run_dir,1351811072)。历史三个cost.json均扣除实计9.2915103s；未计时旧尝试/旧审计保守预留120s。保存密度评价通过，新BIE重解当时RAM门失败且耗时未测，未授权启动。冻结清单在data/frozen_jobs.md。
- P2每批先重新读取Windows物理可用计数器，再execute_two_dimensional_closeout_s2b_job(out,job_id,bytes)，顺序bie3、flat、bie1、bie2、pe1至pe4、fk1至fk4。代表bie3成功0.1216s后才继续。所有调用串行，低内存批次估计96MiB，至少保留当时30%物理RAM，实际计时/内存/时间戳见checkpoints/*.json和data/cost.csv。13.9342315s为12个成功批次总核内实计，共新增60核心值，并有匹配点/线探针的数值敏感性评价。
- 首次bie3在所有密度层评价后、结果封装时失败：compose返回cell使struct扩成数组。仅将receiver_ids转string后重试一次，无声学输入变更。未取得该失败的核内仪器时间，保守扣记30s并单列bie3_failed.mat/json，不伪称实测。
- 第一次纯保存汇总因receipt.job_id为不同长度char而不能vertcat，转string后重跑成功；该错误没有新物理场。修复记录和最终MAT都保留在当前工作流；未修改或重跑S2a。
- P3保存108格、辅助G/delta/cross、2250敏感性复数和全部指标。F01–F07初版生成、CSV独立审计/保存MAT重绘通过。逐张查看后修复F01 Tx标记贴边与F05长标题挤压，另存figures_v2，audit_v2再次复核并新增BIE敏感性CSV独立算术核验；初版目录保留。
- 12个新MATLAB文件初次Code Analyzer仅4条样式/性能消息，修正后零消息。cwd初次assert比较Windows反斜杠pwd与正斜杠root字符串而失败，未改变实际目录；最终按分隔符归一化核验。
- 所有生产求解器和公共API未改，无通信回归；数据齐全但BIE历史执行代码SHA及中强独立自求积层级缺失，FullFK窗口敏感性与方法差异相当。状态S2B_DATA_COMPLETE_REFERENCE_LIMITED，Goal未标complete；不进入S3–S6/Bellhop源映射/通信，不设置自动继续。
- 图形/审计/I/O未纳入声学核内预算；未单独仪器计时，不报告虚构图形耗时。内存为保守增量估计和批前真实物理可用量，不宣称实测峰值。
Final tabular import check initially received available as text true, then explicitly parsed its string value; all108 finite values and all9/411 supports verified. No acoustic calculation repeated.
