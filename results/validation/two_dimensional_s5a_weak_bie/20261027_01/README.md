# 二维阶段 A：6 kHz 弱正弦 P22 五层参考试点

**状态：单点经验 BIE 参考通过，PE/FullFK 相对 BIE 的差异均未分辨；A 整体 PARTIAL，B/C 未进入。** 冻结源、海面与门槛见[中文报告](../../../../reports/two_dimensional_s5_s6_2026-10-27.md)。

- `remote/jobs/{smoke,model,main96,self64,space14,window60,window70}/`：各作业原 MAT、真实退出码、进程时间、完整 stdout/stderr、预检和 SHA。
- `remote/output_manifest.txt`：远端对 56 个任务文件生成的 SHA；回传本机复核 `56/56`。
- `remote/SHA256SUMS`、`remote/final_sha_check.log`：输入 10 项最终 SHA 复核。
- `remote/launcher_smoke.log`：首次仅启动脚本数值格式拒绝，没有 MATLAB 求解；`launcher_smoke_v2.log` 及其余 launcher 日志保留。
- `audit_v2/independent_audit.mat`、`bie_layers.csv`、`p22_one_frequency_response.csv`、`independent_audit.json`：重开 MAT 后独立复算的五层参考、经验分辨率、模型差和共同 DSP 分量。
- `independent_crosscheck.json`：另用 PowerShell 从保存 CSV 与哈希清单复算，并核对退出/预算。
- `checkpoint.json` 与 `sha256_local.json`：范围、成本、版本、门槛、后续可做事项及本地 SHA。

本目录只有单个 6 kHz 频点，不能用其做 IFFT/CIR、九点结论、PM 总体、QPSK 或 BER 损失推断。
