# 固定dual宽带研究决策数据（2026-10-15）

决策B：当前4/6/8kHz离散交付；没有新物理频点、没有宽带扫描或CIR，G1_target_PM=OPEN、WP2_formal=false。

- [物理报告](../../../reports/deepsea_dual_broadband_decision_2026-10-15.md)、[成本](../../../reports/deepsea_dual_broadband_cost_2026-10-15.md)、[预声明协议](../../../reports/deepsea_dual_broadband_protocol_2026-10-15.md)。
- `saved_evidence.mat/json`：从冻结frequency.mat及baseline6导出的原生参考、FK及数值变体、物理表、成本和解析几何。
- `physical_effects.csv`：物理效应、参考经验敏感性、FK数值误差和模型差；`native_pressures.csv`：三频平面/粗糙/FK分量的原生复声压。
- `historical_costs.csv`：已有核实测；`cost_candidates.csv`：网格成本外推及明确假设，非网格充分性认证。
- `locked_candidate_and_gate.mat`和`stage_gate.json`：预声明5/7kHz候选与预测（predicted_not_observed），以及物理内存门拒绝。没有任何5/7kHz真实参考数据。
- `numerical_audit.json`：MAT/CSV复数和恒等式、计数及停止范围；`final_audit.json`：文件结构、前后哈希、静态检查与状态一致性。
- `frozen_before.json`及`frozen_after.json`保护本轮前5882个现有文件，仅正式计划及四个导航文件允许变化；历史结果与求解器不改。`source_manifest.json`给出报告证据及哈希，`execution_ledger.json`记录实际执行。

只读复核可加载本目录MAT/JSON。三个新入口均从mfilename解析仓库根，不依赖cwd；先setup_vertical_project。执行顺序是analyze_deepsea_dual_broadband_decision、preflight_deepsea_dual_broadband、audit_deepsea_dual_broadband_decision。它们拒绝覆盖当前输出，不应重新运行本轮资源门以改变停止决定。未来若另获任务授权，应使用独立case目录及新协议，保留本轮结果。

本轮新增代码全为报告/验证编排，未改变生产或参考算子，故不运行无关生产标量/宽带通信回归。本轮QA不作为物理突破计数。资源门通过不等于物理通过；本轮资源门未通过更不等于新频点物理失败。
