# 6 kHz Q1 BIE 最小检查点证据索引

状态：`STOP_REFERENCE_GATE_FAILED`。仅运行预先冻结的平面制造问题；没有运行 compact-dual 或原 PM。

- [执行前冻结协议](protocol.md)、[8 项输入 SHA 清单](input_manifest.json)、[SHA-256](input_manifest.sha256)、[远端启动器](launch_q1_ubuntu.sh)。
- `n16_q8/`, `n24_q8/`, `n24_q12/`：各自 `assembled.mat` 含原生复矩阵；`solved.mat` 含原生复密度和 rhs；`result.mat` 含原生复非配置边界、离面场及 Fourier–Bessel 参考；`receipt.json` 含原数值与敏感性。
- 对应 `run_n16q8_01/`, `run_n24q8_01/`, `run_n24q12_01/`：完整 stdout/stderr、真实退出码、进程墙钟、峰值 RSS、事前/事后物理 RAM 与负载。
- [49 项远端结果哈希](output.sha256)、[只读独立审计程序](audit_saved.py)、[独立审计结果](audit.json)、[最终检查点](checkpoint.json)。审计程序不启动 MATLAB 或 BIE 求解。
- [技术决策报告](../../../reports/deepsea3d_bie_q1_decision_2026-10-26.md)。远端原始结果还保存在 `/home/misc/cwx/deepsea3d_bie_q1_20261026/results/validation/deepsea_bie_q1_20261026/`。

本轮新增远端 MATLAB 进程保守墙钟 53/600 s，正式作业 11/13/15 s，所有正式退出码 0，49/49 回传 SHA 一致。版本探测单独保守计 14 s；其远端完整日志未保存，是记录局限。独立复算边界/离面误差分别在最佳 n24/Gauss12 为 `0.0250587454203/0.0409758501290`，均不过 `1e−4`。
