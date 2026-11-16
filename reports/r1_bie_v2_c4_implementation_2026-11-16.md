# R1-BIE v2 迁移与仅 C4 手动执行包

2026-11-16。**实施准备完成；未连接SSH、未启动C4、未执行远程Resume、未提交或推送Git。** 批准来源为用户本轮明确指令及[已审查方案](r1_bie_nonmonotonic_policy_c4_plan_2026-11-16.md)。R1参考资格仍未完成，R2未授权。

## 正式交付与变更范围

正式包：`results/comparisons/piston_2d_pm/20261116_r1_bie_c4_policy_v2_04/`。

远端任务目录（由用户手动上传时创建，当前未连接）：`/home/misc/cwx/r1_bie_manual_c4_policy_v2_20261116_04/`，只用`matlab-ubuntu-new`。

- 批处理：`scripts/validation/run_r1_bie_c4_continuation.m`。缺省`dry-run`，正式模式限Ubuntu R2025b、独占launcher会话；只有C4固定配置的一个BIE密度求解调用，没有C0–C3求解分支。
- 新增v2 IO/来源迁移/阶段门/独立最终审计：`scripts/validation/support/r1_bie_continuation_io.m`、`migrate_r1_bie_checkpoints.m`、`r1_bie_continuation_gates.m`、`audit_r1_bie_continuation_final.m`。
- 打包、手动管理、启动与无数值准入：同目录`prepare_r1_bie_c4_package.ps1`、`manage_r1_bie_c4_continuation.ps1`、`launch_r1_bie_c4_continuation.sh`、`r1_bie_c4_preflight.py`。
- 验收：`scripts/validation/validate_r1_bie_c4_preparation.m`及`support/test_r1_bie_c4_package.py`。旧runtime测试/序列化组件、阶段独立审计、R1分析/综合审计/绘图器原字节复用。

旧执行包、原失败结果、BIE原/缓存求解器、PS、FullFK、生产WAPE均未修改。公共声道/通信接口和物理公式不变，因此无需改`vertical_comm_guide.md`，也不运行通信回归。

## 政策与原数值门

新`policy.json`版本`r1-bie-reference-v2`，独立manifest schema 2 / scope `R1_BIE_C4_CONTINUATION`。

`first=1.3977571843886662e-7`，`second=1.5121123234962964e-7`，增长`8.181330805153086%`，`grid_nonmonotonic=true`、`legacy_policy_status=FAIL`均完整保留。旧`r1manual:Nonmonotonic`及原退出1/304 s不改。新阶段门和最终审计核对同一policy SHA；单调性只输出诊断，不抛停止异常。

原四配对（空间10→12、12→14，自求积64→96，窗口60→70）的全411点/四分段最坏复L2均须严格`<1e-4`，等于门槛仍失败。线性`<1e-10`、全部边界与独立19点`<1e-4`、重建/CSV/rhs/闭合`<1e-11`以及几何、源梯度、PS/FK既有门均原样保留。非单调项照常进入经验不确定性，不被删去或归零。

## 本地迁移与验证的实际结果

正式包manifest SHA：`36b0dfef6ddae0b1e08fe96011781f75e04f238ab53f9fb131a8da330c721171`。

Policy SHA：`07ce1a9d14bee14a4e14490bbd282f07af1c98909e7a6f36816cb5726e75b940`。

来源为原回收run `r1_bie_manual_runtime_fix_20261116_02_retrieved_20261116_144235`。原110项回收SHA、36项包SHA、C0–C2提交链、C3 solved及原审计均通过。原件未写入；包内`origin_evidence/`保留清单绑定的110文件及原回收SHA清单。新`migration_sources.json`锚定原manifest/清单，旧新链并存，未篡改原config中的manifest引用。

C0–C2继承已验证的审计；C3只在`migration_audits/C3/`重建审计后封存。新链为C0→C1→C2→C3的`migration_commit.json`，末端由`migration_complete.json`封存。新审计由现有本地MATLAB R2025b Update 6完成，**本轮新增密度求解0**，C3审计约11.07 s。

| C3新工作区审计 | 结果 | 原门 |
|---|---:|---:|
| 保存密度到411点场重建相对差 | 1.40108e-15 | <1e-11 |
| MAT/CSV最大绝对差 / Pa | 6.66785e-14 | <1e-11 |
| 入射rhs相对差 | 1.84534e-14 | <1e-11 |
| 全476平台中点相对残差 | 7.74972e-8 | <1e-4 |
| 独立19点GL112相对残差 | 9.67546e-8 | <1e-4 |

新包的`local_validation.json`与`local_validation_commit.json`绑定本地验收；启动器在MATLAB前同时验证包SHA、迁移链与本地准入。缺验证、缺C3审计、改政策/前序/数据/代码等均拒绝。C3在不同平台的独立边界审计值与远端历史值略有差别，均满足原门，不改历史审计，不声称逐位相等。

| 验收 | 实际完成 |
|---|---|
| MATLAB | 31项检查通过，约32.14 s；6个新增MATLAB文件Code Analyzer 0消息；政策下降/相等/上升、恰等于1e-4/超限、非有限/缺支持/改门拒绝，C4状态与SHA恢复、迁移/跳过、最终审计政策一致和失败路径均覆盖 |
| 求解禁止 | 4个BIE/PS/FullFK求解入口由哨兵遮蔽，确认哨兵解析位置，无调用日志；正式包没有C4数据或求解标记 |
| Python证据与准入 | 25项通过，篡改/缺原MAT、CSV、C3审计、原/新前序链、缺迁移完成/验证封存、旧scope、改政策/批准、越界路径、改数值门均拒绝 |
| Bash/runtime | Bash语法检查及无MATLAB启动dry-run通过；原29项runtime测试通过，覆盖available大/等/小于initial、空值、非法值、文件缺失/读取失败及JSON解析 |
| PowerShell | 两个新工具解析通过；ShowCommands实际执行，无网络连接 |
| SHA空文件 | MATLAB对空文件的SHA已与标准已知值核对通过；新helper跳过空读取，不修改旧helper |

最终审计测试使用`FINAL_TEST_ONLY_*`目录中的**模拟C4窗口场**（复用C1场），只验证保存数据后处理/政策/拒绝逻辑，没有计算窗口传播或C4密度。该测试目录和模拟产物不进入正式包，也不是C4收敛证据；报告中的参考完成标记仅属于测试夹具，不能用于研究结论。

试验历史保留：候选01在旧SHA helper遇到空日志的Java重载错误时失败，随后仅新增v2 SHA实现修正；候选02通过27项，候选03增加最终审计测试后通过31项；核对到本地会话默认8线程后，候选04将迁移/验收临时限制为单线程，31项复验通过并恢复会话原8线程设置，正式交付仅04。Git Bash缺`python3`命令的首次dry-run失败保留；成功的本地测试仅在测试shell中将python3映射到python.exe，正式Ubuntu仍明确依赖python3。所有旧包/失败记录均不覆盖。

验收目录：`results/validation/r1_bie_c4_preparation/20261116_candidate_04/`、`20261116_integrity_04/`、`20261116_shell_04/`及按相同依赖SHA复用的`20261116_runtime_03/`。正式包内保存核心验收回执和shell/runtime日志；原始MATLAB diary与完整模拟夹具留在本地验收目录。完整远程Linux后台锁/信号/RAM杀停生命周期及真实C4/五阶段最终审计仍未运行，不能声明这些实机验收通过。

## 用户手动上传与启动

以下在Windows PowerShell 7执行；这些命令**本轮没有执行**。先核对包内`input_package.tar.sha256`；管理工具会验证本地准入、归档SHA、远端归档SHA及任务目录唯一性，然后单独后台启动一个MATLAB进程。

```powershell
$project = 'E:\MISC\CARPE3D_matlab\Explain'
$pkg = "$project\results\comparisons\piston_2d_pm\20261116_r1_bie_c4_policy_v2_04"
$manage = "$project\scripts\validation\support\manage_r1_bie_c4_continuation.ps1"

# 本地展示命令，不连接远程
& $manage -Action ShowCommands -PackageDir $pkg

# 用户决定启动时才执行：上传、校验并后台启动
& $manage -Action UploadAndStart -PackageDir $pkg
```

远端等价启动入口（包已上传并解包后，仅用户执行）：

```bash
bash /home/misc/cwx/r1_bie_manual_c4_policy_v2_20261116_04/code/launch_r1_bie_c4_continuation.sh start /home/misc/cwx/r1_bie_manual_c4_policy_v2_20261116_04
```

不需要Codex持续轮询。用户可以自行等待，任意时刻用以下命令查询；没有退出码不代表数值通过。

```powershell
& $manage -Action Status -PackageDir $pkg
```

状态返回checkpoint、当前session、stdout/stderr路径及存在时的真实`matlab_exit_code.txt`。完成依据是阶段/最终审计与封存和退出码共同有效，不是单独退出0。

任务停止后回收，禁止覆盖既有下载目录：

```powershell
& $manage -Action Fetch -PackageDir $pkg
# 或指定一个尚不存在的目录
& $manage -Action Fetch -PackageDir $pkg -DownloadDir "$project\results\comparisons\piston_2d_pm\r1_bie_v2_c4_retrieved_01"
```

Fetch拒绝对仍运行的任务快照，回收后验证全部清单SHA。不要重新用UploadAndStart覆盖已存在的远端目录。

## C4配置、输出与恢复边界

C4固定L70 m、ppw12、qa96、GL8、4480节点/560面板，全部476平台中点和411真实接收点，原缓存true；只有此阶段允许新密度求解。`solve_invocation.json`在进入求解器前写入；中断后无有效solved封存时停止，不能自动再调用求解器。

- C4有效`solved.json`：只补审计；有效`committed.json`：跳过求解与阶段审计，只继续最终化。
- 仅求解标记/孤立数据、SHA漂移或损坏封存：停止人工审查，不自动删除或重算。
- 最终化失败：旧attempt保留，下一次可新建最终化attempt，不重算C4；有效最终封存则直接验证。
- 审查失败原因后，如属于允许的“仅补审计/最终化”情况，用户可执行`& $manage -Action Continue -PackageDir $pkg`；它保留原目录并重新启动相同v2状态机，绝不调用旧通用Resume。若C4求解已调用但未封存，Continue会拒绝重求解，须另行明确授权处理。

运行环境仍为Ubuntu MATLAB R2025b，一个`-singleCompThread`进程，无parpool。取消的1100 s累计、730 s预测及旧900 s会话门不恢复。RAM初始值从原已核验文件继承，当前starting单独读取，baseline=max(initial,starting)，保留30%基准可用RAM和12 GiB增量准入估计；空值、非法数据、内存读取/JSON失败均在MATLAB前停止。C4函数粗略2–3分钟，全流程可能5–10分钟或更久；这不是承诺或新的时间硬门，内存12 GiB也不是实测峰值上界。

最终比较仍以C3为BIE参考，C4只增加窗口项：

```
bunc = |u0-u1| + |u1-u3| + |u2-u1| + |u1-u4|
R_m  = 5 * norm(bunc + model_unc_m) / norm(u3)
E_m  = norm(u_m-u3) / norm(u3)
resolved_m = E_m > R_m
```

PS/FK的model_unc、0.01 RMS相位掩码、全线和逐点误差口径不变。完成后保存各阶段MAT/CSV/配置/计时/资源/日志/SHA、独立审计、完整参考不确定性、411点误差、PNG/FIG及曲线/像素重建验证。`resolved=false`表示当前经验分辨率下不能分辨，不代表模型严格一致；经验不确定性不是严格误差上界。

完成后程序停止等待R1审查。不得自行扩展种子、频率、R2、CIR/通信或Bellhop。
