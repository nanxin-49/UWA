# 二维 S5 P22 远端 MATLAB 非交互许可证故障与恢复门（2026-10-30）

## 当前结论

本项目在 `3D_try` 分支使用远端 Ubuntu 的 `/home/misc/MATLAB2023b/bin/matlab -batch`。**已观察到新进程间歇性进入 MathWorks 账号交互登录，不能可靠地作为无人值守求解器。** 这表明当前运行方式未稳定取得非交互许可；日志本身不能证明是“没有许可证文件”，也不能单凭登录提示区分授权模式、账号关联、网络许可服务或会话状态的问题。不得据此改动声学模型、数值门槛或已保存的频点结果。

用户要求暂时停止本 Goal。原 1800 s 预算按未观测作业的上限保守占用 `1780 s`，剩余 20 s 不足以启动 BIE。4031.25 Hz 仅有 PE/FullFK、`main96` 和 `self64`；`space14/window60/window70` 无结果，五层经验参考、双频资格和 CIR/BER 均未完成。停止期间不启动 MATLAB、不清理共享账号的授权文件、不重试登录、不触碰其他任务。

## 可复核的故障时间线

| 事件 | 原始证据 | 观察与结论 |
|---|---|---|
| 首次预测冻结 | [原 `freeze` stdout](../results/validation/two_dimensional_s5_p22_dualfreq/20261030_01/remote_mirror/jobs/freeze/stdout.log)及 `exit_code.txt` | 出现 MathWorks 邮箱提示，75 s 超时，退出码 124；预测当时未生成。 |
| 受控预测重试 | [续算镜像 `freeze_retry1`](../results/validation/two_dimensional_s5_p22_dualfreq/20261030_01/continuation_remote_mirror_01/jobs/freeze_retry1/) | 新 MATLAB 进程退出 0、12 s；预测 JSON SHA `be3a54b96584b642729407947874d1df4f0bcce1341ba9c9d17a60e53ac5feee`。一次外层 `sha256sum -c` 格式误用另存，不影响 MATLAB 真实退出码。**这一次成功不证明后续新进程也能稳定获许可。** |
| 半步已完成作业 | [4031.25 Hz 作业镜像](../results/validation/two_dimensional_s5_p22_dualfreq/20261030_01/continuation_remote_mirror_01/jobs/f4031p25/) | `model/main96/self64` 各退出 0，用时 54/155/166 s，MAT/JSON SHA 通过；物理审计仅覆盖两层 BIE。 |
| 首次余三层合并尝试 | 同上 `group_remaining3/` | 进程退出 0、54 s，但无任何层结果；外层缺产物门失败。属启动包装/执行问题，不是 BIE 数值结果，原失败保留。 |
| 第二次余三层合并尝试 | [`group_remaining3_retry1/stdout.log`](../results/validation/two_dimensional_s5_p22_dualfreq/20261030_01/continuation_remote_mirror_01/jobs/f4031p25/group_remaining3_retry1/stdout.log)及 [未观测退出说明](../results/validation/two_dimensional_s5_p22_dualfreq/20261030_01/continuation_remote_mirror_01/jobs/f4031p25/group_remaining3_retry1/exit_unobserved.txt) | 出现邮箱、密码提示及“登录失败”；没有 `GROUP_START` 或 BIE 层结果。包装进程未留下真实退出码和耗时，不能补造；按预设 510 s 上限计账。 |

[独立部分审计](../results/validation/two_dimensional_s5_p22_dualfreq/20261030_01/audit_4031_partial_01/audit.json)重开 MAT/JSON 并核 SHA；原镜像 96/96、续算镜像 69/69。完整资源和状态见[续行检查点](../results/validation/two_dimensional_s5_p22_dualfreq/20261030_01/checkpoint_continuation_20261030_02.json)。

## 由许可证所有者或管理员执行的恢复步骤

1. 先确认共享账号 `misc` 对 **MATLAB R2023b** 实际应使用的授权类型和可用许可，由许可证所有者或学校/实验室管理员核对账号关联、有效期、激活方式及必要的许可服务器可达性。MathWorks 的[许可文件与授权类型说明](https://www.mathworks.com/help/install/ug/locate-matlab-license-files.html)说明不同许可类型走不同路径；这里不能根据交互提示擅自假定某一种。
2. 按该授权类型由有权限的人恢复安装的合法激活或许可连接。MathWorks 的[手动激活说明](https://www.mathworks.com/help/install/ug/activate-matlab-installation-manually.html)覆盖 R2022b 及更新版本，包含 Linux 的激活客户端；组织许可的详情应由管理员提供。不要在聊天、仓库、作业日志或临时脚本中提供/保存邮箱、密码、激活密钥或许可证文件内容。Codex 不修改 `/home/misc/.ssh`、全局 shell、MATLAB 安装、用户授权缓存或他人会话。
3. 修复后由用户告知。重新以 `ssh -o BatchMode=yes matlab-ubuntu` 和同一 R2023b 可执行文件做**新的短 `-batch` 进程**，只运行确定性标记与版本输出，保存 stdout、stderr、真实退出码和耗时；要求无账号/密码提示、标记实际出现、退出码 0、进程完全退出。MathWorks [Linux 启动文档](https://www.mathworks.com/help/matlab/ref/matlablinux.html)将 `-batch` 定义为非交互执行，并规定成功自动退出 0、失败非 0；此前“MATLAB 退出 0 但无目标输出”的案例说明还必须核标记和预期产物，不能只看退出码。
4. 因本次故障是**间歇性的**，单次短作业成功后仍须在下一独立新进程启动前重查资源、预测 SHA 和许可证提示；若再次要求交互登录，立即停止，不输入凭据、不把旧成功当成本次资格。
5. 原预算不能续算。完成许可证恢复通知后，按[下一轮申请](two_dimensional_s5_p22_next_goal_prompt_2026-10-30.md)另行取得**最多 1200 s** 的明确物理计算授权，才可在原任务目录用新作业日志补算三层。各进程 <900 s、串行、逐层保留至少 30% 可用物理 RAM，任务盘与 `/tmp` 各至少 10 GiB，输入/预测 SHA 不变；保留所有失败目录。五层 BIE 经验敏感性 `≤1e-4` 与模型差超过五倍经验分辨率的门槛不变。

**恢复验收状态：`NOT_RUN`。** 本报告是故障定位和操作方案，没有声称许可证已经修复，也没有请求或实施新的远端计算。
