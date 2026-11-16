# R1-BIE v2 Status 查询修复（不启动/重启计算）

2026-11-16，Asia/Shanghai。用户报告首次Status显示“没有MATLAB检查点/退出回执”，随后在状态工具更新过程中又遇到`bash -c`第5行的引号语法错误。本轮只修复客户端只读Status及验证，不修改已封存执行包/数值配置/门槛，不启动、Resume、终止或重算MATLAB，不提交/推送Git。

## 查明的状态与根因

指定任务：`/home/misc/cwx/r1_bie_manual_c4_policy_v2_20261116_04/`，session `20261116T091043_470983242_3891066`。

1. **首次提示不是求解错误。** 通过`matlab-ubuntu-new`只读检查，在17:13:07观察到检查点为`RUNNING / C4_SOLVING`，MATLAB PID3891114仍运行，stderr为空；C4求解标记计数1。启动session为17:10:43，首个C4检查点到17:11:34才生成，约51秒的MATLAB启动/输入/SHA核验空档内，旧Status只给出了模糊的“尚无检查点/退出回执”。它不构成数值失败，也不应据此重启。
2. **后来确有状态工具的引号错误。** 我修改客户端Status时把带Python字符串的Bash脚本嵌入`bash -c`；Windows PowerShell/OpenSSH的多层引号未保留，用户得到`python3 -c '''...`语法错误。它发生在只读状态脚本解析阶段，不会进入求解分支。
3. 中间尝试的Base64命令行方案在本地原生命令链测试中暴露长参数截断问题，未作为最终方案。最终使用.NET `ProcessStartInfo.ArgumentList`及UTF-8标准输入，脚本文本不放入命令行，不依赖shell引号拼接，也不依赖长命令行；同步捕获stdout/stderr及真实退出码。

## 已实施修改

- `scripts/validation/support/manage_r1_bie_c4_continuation.ps1`：仅修改`Status`分支；仍使用固定别名`matlab-ubuntu-new`，远端命令为短命令`bash -s -- <已验证的任务路径>`，通过标准输入发送只读脚本。
- 新`invoke_r1_bie_status_stdin.ps1`：PowerShell 7原生进程标准输入传输；UTF-8无BOM、LF、参数数组、并行读取stdout/stderr、真实exit code、隐藏控制台窗口。没有启动MATLAB的代码。
- 新`status_r1_bie_c4_session.sh`：只读取指定任务的session、检查点、记录PID所属SID和匹配本任务入口的进程、资源采样、退出码及必要错误尾部；不写文件、不持锁、不发信号、不启动计算。区分`MATLAB_INITIALIZING`、`RUNNING`、wrapper尚在、真实失败、状态不明及退出0待核验等情况；旧session检查点不充当新进度，PID复用/不匹配进程不被误认。
- 新`test_r1_bie_c4_status.ps1`及`test_r1_bie_c4_status.sh`：本地合成数据和原生标准输入回归；`.gitignore`显式保留这些新工具。
- `prepare_r1_bie_c4_package.ps1`：未来新包显式包含两个Status依赖；本轮没有重新打包、上传或改写04包。

查询输出中的`SESSION_AGE_SECONDS`表示距启动的时间；存在真实运行时长文件时另列`ACTUAL_PROCESS_ELAPSED_SECONDS`，避免把完成后等待查询的时间计作计算成本。退出0仅说明进程退出，不由Status直接宣布参考审计PASS。

## 验证与当前结果

- PowerShell本地7项通过：三个文件语法检查、单/双引号/括号/中文真实往返、超过20KB标准输入、stdout/stderr与退出7的正确保留，以及传输选择核对；无SSH/MATLAB调用。
- Bash本地12种状态通过：未启动、启动初始化、C4运行、旧检查点、非法JSON、wrapper、PID复用、无进程无退出码、退出0、退出失败、非法退出码、非法session。运行时长字段同时核对。
- 修复后通过真正的Windows OpenSSH路径进行一次只读远端Status调用，**成功退出0**。原session当前返回`CHECKPOINT_STATE=COMPLETE`、`CHECKPOINT_STATUS=R1_REFERENCE_COMPLETE_WAITING_REVIEW`、`REAL_EXIT_CODE=0`，C4调用数1。
- 本轮未回收并重建完整MAT/CSV/FIG数值审计，因此不凭该状态输出升级为最终R1科学审查通过。未重复启动、没有持续监控，已停止。

测试及远端输出保存在`results/validation/r1_bie_c4_status_fix/`，失败的早期transport测试也保留。当前04包及其归档SHA `901cd359adfab5b8aa36f247d5828de45803a21f017519e729973fbf7c9d15f7`不变，远端运行代码/清单/失败证据未修改。新版Status从本地传输只读脚本，不要求上传、替换或重新封存正在使用的04包；旧包内复制的旧manager不会自动升级，使用仓库当前脚本。

## 用户命令

```powershell
$project = 'E:\MISC\CARPE3D_matlab\Explain'
$manage = "$project\scripts\validation\support\manage_r1_bie_c4_continuation.ps1"
$pkg = "$project\results\comparisons\piston_2d_pm\20261116_r1_bie_c4_policy_v2_04"
& $manage -Action Status -PackageDir $pkg
```

**不要重新UploadAndStart或Continue来修复Status。** 当前会话已退出0；下一步是用户按既有命令Fetch回收并另行审查完整数值结果，不自动进入R2。本轮未执行Fetch。

此修复不改变物理/通信接口，不需要MATLAB声学回归或修改`vertical_comm_guide.md`；看板与任务表只登记已观察到的退出状态及“待回收审查”的边界。
