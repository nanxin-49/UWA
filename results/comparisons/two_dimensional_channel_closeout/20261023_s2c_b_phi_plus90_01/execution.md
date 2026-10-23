# PHI_PLUS90 执行回执（2026-10-23）

## 2026-10-24 预算授权与 RAM 门补记（覆盖下方旧成本阻塞状态）

用户明确将共享 7200 s 上限增加至少 306.601 s，并要求保留 S3A 3000 s 规划储备。按最小增量记账，共享上限为 7506.601 s；前期扣记 2685.7140642 s，PHI 可用子额度 1820.8869358 s。原保守完整计划 1820.88667776 s 因而通过成本门，修订原件见 [P2_budget_amendment_7506_601.json](P2_budget_amendment_7506_601.json)。原 `cost_gate.json` 是授权前快照，不改写。

接着两次读取 Windows 可用**物理**内存，分别为 1,641,025,536 B 和 1,131,163,648 B；冻结最大 BIE 层级预计增加 1,418,723,328 B，按保留 30% 的规则至少需要 2,026,747,612 B。两次均未通过 RAM 门，最新测量见 [P2_ram_gate_20261024.json](P2_ram_gate_20261024.json)。用户随后要求暂时保持资源阻塞。未启动 BIE、PE 或 FullFK，新 PHI 声学耗时仍 0 s；共享核心扣记仍 2685.7140642 s，S3A 3000 s 储备仍完整。状态维持 `S2C_B_PHI_PLUS90_RESOURCE_REFERENCE_BLOCKED`，此轮直接原因已由成本改为 RAM。

已经增设只对 PHI 适用的预算门、五层 BIE、PE、完整 FullFK 双项积分入口；五个 MATLAB 文件的 Code Analyzer 均为零消息。入口尚未执行，不能据此宣称声场、参考资格或 M05b–e 已完成。下一次仅在用户告知已释放内存后重测物理 RAM；通过后仍逐作业重新过成本、RAM 和数值门。

后续只读准备补充了[五层 BIE 保存场参考审计入口](../../../../scripts/validation/support/audit_two_dimensional_closeout_s2c_b_phi_bie_reference.m)（尚未执行，Code Analyzer 零消息），以及[对称性输入审计](P1_symmetry_input_audit.md)：417/417 接收点有同深度 ±x 配对，五层 BIE 节点与 PHI 高度/坡度满足浮点精度内的偶/奇对称；PE 原 984 点 FFT 网格有一个不配对的 Nyquist 端点。此检查没有新增声场，不改变 RAM 阻塞状态。

保存场分析入口 [`analyze_two_dimensional_closeout_s2c_b_phi.m`](../../../../scripts/comparisons/analyze_two_dimensional_closeout_s2c_b_phi.m) 已将 BASE_M 与 PHI 的同点原生复场、平面控制、五支持、幅相平方分解、BIE 物理响应和模型残差分开组织，并检查预算修订 SHA。PHI BIE 参考审计入口另从保存节点核对五层高度和法向。两个入口 Code Analyzer 均为零消息，但缺新 PHI 场而未执行，表格、M05b–e 和独立结果审计仍不存在。

1. 在干净工作区创建唯一新 run；核对当前 S2b、S2c-A、K_HIGH 最终哈希清单，分别为 247/247、203/203、132/132。已有 MATLAB MCP 会话 R2025b Update 6，工作目录项目根。
2. 冻结 PHI 输入 SHA、原源和 417 点接收位置；从 S2c-A 真实解析几何生成 M05a。初版图例把辅助竖线和中心标记列入图例，在 `figures_v2/` 修正；最终又关闭下划线 TeX 解释，正式图在 `figures_final/`。旧图保留为过程证据。没有生成 PHI 声学场。
3. 按 K_HIGH 七作业回执重估五层 BIE，加入 20% 保守增量及两项 150 s 预留，完整计划 1820.8866778 s 大于保持 S3A 储备后的 1514.2859358 s。P2 成本门失败，未进行即时 RAM 许可检查，也未启动 BIE 代表、PE 或 FullFK。新声学耗时 0 s；共享累计与余额不变。
4. K_HIGH M04 小包复制首版清单生成时，PowerShell 在枚举输出目录时同时碰到正在写入的清单自身，产生一次文件占用报错；已核对最终清单恰含八个载荷，各与原件 SHA 相同，未覆盖来源。保存复场十字五点首次手动复核选错两个深度点，按 P12/P21/P22/P23/P32 的真坐标修正并通过；旧正式 CSV 未改。
5. 从 S2c-A/K_HIGH 保存 CSV 复核三档 A 与 K_HIGH 的空间支持和幅相平方项；核对 S3A 三个已有 PM Fourier 输入的系数 SHA、同父谱和旧 MAT 接收坐标，不执行 PM 求解。M05a 从保存 MAT 重新绘制后像素一致，FIG 重开三轴，静态代码检查零消息，`git diff --check` 通过。产物索引见 [M05图册](mechanism_atlas.md)、[独立审计](independent_audit.md)、[成本JSON](cost_gate.json)及[最终审计JSON](final_audit.json)。
6. 下一回合只读补查旧 PM pilot 的 `bie_reference_metrics.csv`：只有 seed260001 两个旧 L60 层级，seed260002/003 没有对应 BIE 回执。按 K_HIGH 五层实耗直接乘三个 seed 已需 3802.217 s，加本轮 20% 保守增量仅 BIE 就约 4562.660 s，略超共享剩余；这是未经 PM 节点/版本资格确认的规划情景，不作为 S3A 正式估价或执行许可。详见 [S3A只读清单](S3A_readonly_feasibility.md)。

阶段分类为 `S2C_B_PHI_PLUS90_RESOURCE_REFERENCE_BLOCKED`，不等于中心过零机制已验证。未改变生产求解器、公共 API、旧正式数据、3000 s S3A 规划储备；未提交或推送。
