# R1-BIE 非单调政策调整、证据迁移与 C4 续行方案

日期：2026-11-16。状态：**静态方案完成，等待用户批准实施；不是已实现的续行包，也不是运行许可。**

本轮不连接 SSH、不调用 MATLAB、不执行 Resume、不打包或运行新声学任务；不改旧源码、原数值门、原结果和失败记录。原始诊断见 [Nonmonotonic 诊断](r1_bie_nonmonotonic_diagnosis_2026-11-16.md)。当前 Git 基线为 `compare / a8837d1`。

## 1. 建议与证据边界

建议在**新的、有版本与批准记录的续行包**中，将严格空间单调性改为必须保存的诊断；四组敏感性的原 `<1e-4` 门、残差、独立审计及经验分辨规则全部保留。只允许新增 C4 一次 BIE 密度求解；C0–C3 只导入、核验和必要的保存密度重建，不重新组装求解密度。旧包和旧 `FAIL` 永久保持原样。

已有证据：C0 为导入的资格回放，C1/C2 有完整 `committed.json`；C3 的 `solved.json` 和独立审计通过，但 `gates.json=FAIL`、没有 `committed.json`。原会话的新求解数 3 对应 C1/C2/C3。C4 未运行，R1 参考资格仍未完成，R2 未授权。

| 配对 | 保存的全线/四分段最坏复 L2 | 状态 |
|---|---:|---|
| C0→C1，空间 10→12 | 1.3977571843886673e-7 | 原 1e-4 门通过 |
| C1→C3，空间 12→14 | 1.5121123234962967e-7 | 原 1e-4 门通过，严格单调门失败 |
| C2→C1，自求积 64→96 | 1.0420813026046337e-7 | 原 1e-4 门通过 |
| C1→C4，窗口 60→70 m | 未知 | 尚未执行 |

空间邻层差增加 1.14355139e-8（8.18133%），最坏支持同为 Q4，x=20.1171875…40.0390625 m。238/411 点的后一次绝对差增大，不是一个点异常。另一方面，C0→C3 最坏差约 1.2251993e-7，小于 C0→C1。相邻差是两个近似解的差，并不是对未知真解的误差；例如误差向量间的抵消就能改变邻层差的大小。因此严格单调不是这里已经证明的必要条件。

这是调整控制政策的理由，不是证明误差消失、达到舍入地板或最终收敛的理由。空间/自求积敏感性同为约 1e-7，窗口项仍缺；不作 Richardson 外推，不假设已知收敛阶，也不从现有摆动推算 C4 必过。

## 2. 代码审查：两处硬门及不能顺手改变的行为

本轮核对以下 8 个当前文件与回收包 `code/` 的同名文件 SHA，全部一致：主入口、阶段审计、最终审计、原/缓存 BIE 求解器、R1 分析入口、综合独立审计、绘图入口。精确 SHA 与检查回执保存在 `results/validation/r1_bie_c4_continuation_plan/20261116_static_01/`。

| 位置 | 原行为 | 新包拟采用的行为 |
|---|---|---|
| `scripts/validation/run_r1_bie_manual.m:193–218` | 每个可用配对 `<1e-4`；stage≥3 时以 `second>first` 抛出 `r1manual:Nonmonotonic` | 旧入口不改。新入口调用 v2 阶段门，保留四配对门，将单调判断写为诊断 |
| `scripts/validation/support/audit_r1_bie_manual_final.m:7` | `all(values<1e-4) && values(2)<=values(1)` | 旧审计不改。新增 v2 最终审计，保留前项，后项变为必须核对的诊断记录 |
| `scripts/comparisons/run_piston_pm_r1_job.m:69–109` | 完整分析以 C3 为 BIE 参考；四项敏感性相加；计算 5 倍分辨标记 | 原文件、公式、参考选择和输出含义不改，复用其 `analyze` 分支 |
| `scripts/validation/support/r1_bie_manual_io.m:22–67` | manifest schema/scope 和签名固定，旧 commit 必须对应旧指纹 | 不让旧 helper 接受被改签名；新 helper 用独立 v2 schema/scope 验证迁移链 |

注意顺序：旧阶段门保存配对依次为“空间10→12、自求积、空间12→14、窗口”；最终分析依次为“空间10→12、空间12→14、自求积、窗口”。新政策用配对名称映射，独立复算后核对这两种既有输出顺序，不能误把第二项自求积当空间细化。

拟议政策文件 `policy.json`：版本 `r1-bie-reference-v2`，`grid_monotonicity_mode=diagnostic`；包含原政策说明、批准记录路径、四个配对名称/方向/支持/分母、全部原门和经验分辨规则，纳入新 manifest SHA。旧版严格政策不被覆盖。阶段评估、恢复、最终审计必须核对同一 policy SHA/版本；版本不符停止。

诊断至少包含 first/second、差值、比值和百分比、`grid_nonmonotonic=true/false`、各分段/全线指标、最坏支持、逐点残差、峰位置、旧失败标识与会话路径。first=0 时比值为 JSON null 并注明原因，不生成 Inf/NaN JSON；场或必须指标非有限、配对缺失、坐标不符仍停止。非单调本身不改变 `numeric_gates_pass`，用独立的 warning 和 `legacy_policy_status=FAIL` 记录，不能写成“严格单调已通过”。

## 3. 保持原样的验收清单

以下严格不等号、支持、数据方向与分母必须以冻结源码为准；迁移不是放宽容差。

| 类别 | 原条件 / 口径 |
|---|---|
| 固定物理输入 | 总宽0.10 m、半宽0.05 m，v0=0.001 m/s；4 kHz；Tx100 m，c1500 m/s，rho1025 kg/m³；seed260001，160 m/12模/Kmax0.471238898038469；C2平台42 m/归零50 m；phase contract 2.0.0；411真实位置、z=3 m |
| BIE配置 | GL8、h=-2 m、窗口平台比例0.85，真实法向/弧长/原核/入射驱动；原已验证缓存显式 true；双精度单线程 |
| 求解残差 | `boundary.finite`；线性 `<1e-10`，全部保存离节点边界残差 `<1e-4` |
| 完整边界支持 | C0/C1/C2/C3/C4分别340/408/408/476/476个平台中点；坐标与全部应有中点精确一致；19个GL112点是额外独立检查，不替代完整支持 |
| 阶段独立审计 | 接收重建相对差、入射rhs相对差各 `<1e-11`；MAT/CSV及分量闭合绝对差各 `<1e-11 Pa`；保存全中点与独立19点边界相对残差各 `<1e-4`；节点/法向几何差 `<1e-12` |
| 四个BIE配对 | 全411点与四个 `round(linspace(1,411,5))` 分段（共享端点）最坏复L2，各 `<1e-4`；分母是该配对第二场的对应范数，保持原 `max(norm,realmin)` 定义 |
| PS既有数值门 | `ps_values(1:2)<1e-4`，步长项 `<1e-10`；只核验已保存结果，不重跑PE/PS |
| FullFK既有数值门 | 所有 `fk_values<1e-5`；复用完整双项与入射梯度，不重跑FullFK收敛扫参 |
| 最终综合独立审计 | CSV与逐点指标 `<1e-11`、全线L2重建 `<1e-12`；独立角谱源压力/梯度相对差各 `<1e-7`；rhs/BIE密度重建 `<1e-11`；几何 `<1e-12`；独立边界 `<1e-4`；FullFK两项重建 `<1e-11`；PS直接Green重建 `<1e-4`；已保存PS/FK数值门全部保留 |
| 最终分辨审计 | `bunc/punc/func`逐点复算绝对差 `<1e-11`，四配对及E/分辨阈值复算差 `<1e-12`；C3参考、RMS、0.01 RMS相位掩码及其有效点数不改 |
| 图件 | 原12条误差曲线与保存数据一致；FIG重开、两组PNG同尺寸重导出及像素差0；失败保留，不以重画掩盖原失败 |
| 资源与退出 | R2025b一个`-singleCompThread`进程、无parpool；可用RAM≥0.30×RAM基准；保留原12 GiB增量准入估计、完整日志与真实退出码；不恢复已取消的1100 s累计/730 s预测/旧900 s会话门 |

门的边界测试也保持：数值恰等于1e-4应失败；不能将 `<` 改成 `<=`。迁移中的必需文件/配对完整性检查属于证据完整性要求，不能用“原函数曾跳过缺文件”作为缺数据的放行理由。

## 4. 不覆盖旧证据的迁移与重新封存

来源固定为：`results/comparisons/piston_2d_pm/r1_bie_manual_runtime_fix_20261116_02_retrieved_20261116_144235/`。未来新目录采用独立 run_id，例如 `20261116_r1_bie_c4_policy_v2_01`，对应独立远端任务名 `r1_bie_manual_c4_policy_v2_<id>`；这里只是命名方案，尚未创建执行包或远端目录。

### 4.1 只读验证和来源冻结

1. 重新验证原回收清单全部110项SHA及原manifest全部36项绑定；检查路径不可越界、清单不可缺失或重复、文件集与原回收清单对应。将原清单本身的SHA锚定到迁移来源记录，而不把自洽哈希称为外部数字签名。
2. 以旧manifest SHA为起点，依次验证C0/C1/C2 solved/committed的原签名、绑定文件、前序commit SHA，以及audit PASS/failed_checks=0和原gates PASS。C3验证原solved的4项绑定和以C2 commit为前序的签名；额外用回收清单绑定C3 audit、audit_reconstruction、gates及失败attempt，确认其只有非单调硬门失败，而不是掩盖其他失败。
3. 核对common、frozen_source/reference_data、PM系数、PS/FK保存结果、Gaussian条件配对输入及所有求解/源/海面/阶段重建代码。8个本轮核对文件的SHA见静态记录；后续实施须按包依赖完整核验。任何漂移、缺件或额外数值失败均停止迁移，不能自动重算来补齐。
4. 复制原证据到新包 `origin_evidence/`，保留原目录相对布局、session退出1、304 s、原日志/attempt/failed gates/SOLVED及全部旧commit。逐文件复制后再核SHA；不移动、硬链接、改写或重新签署原件。新执行MATLAB路径只加入新包`code/`，绝不将`origin_evidence/code/`加入路径。

### 4.2 新链与原链并存

新包使用独立 schema/scope（拟为2 / `R1_BIE_C4_CONTINUATION`）。先写入不可变 `migration_sources.json`（旧manifest、110项清单、原检查点/session/阶段封存与审核文件SHA）、`approval.json`和`policy.json`，连同代码/输入及旧证据SHA清单纳入**新manifest**。新manifest是新链起点，不能把旧commit签名字符串替换为新指纹后继续称其为原commit。

在新 `jobs/` 放置逐字节相同的C0–C3 data/points/config/timings；保留原config中旧manifest指针，另用 `effective_config.json` 描述新包、来源、政策和审计复用关系，不改原numeric_config。原solve时间在独立的historical字段保留；migration/re-audit时间另计。

| 阶段 | 新包允许的动作 | 新封存的含义 |
|---|---|---|
| C0/C1/C2 | 来源链和SHA有效时复用原独立审计；只重新评估保存场的适用配对/政策 | `IMPORTED_AUDITED_PASS`；不是本轮求解或重新完成原审计 |
| C3 | 导入solved密度和场；在新隔离审计工作区调用未改的阶段独立审计，再复算三个配对与非单调诊断 | `IMPORTED_AUDITED_PASS_WITH_NONMONOTONIC_DIAGNOSTIC`；注明旧政策FAIL，新政策已获批准，不伪造旧committed.json |
| C4 | 唯一可到达新密度求解分支的阶段 | 求解封存→独立审计→四项配对/政策核验→`AUDITED_PASS` |

新 `migration_commit.json`（C0→C1→C2→C3）逐级绑定新manifest、policy、原封存文件SHA、原/新审计及gates_v2、数值载荷SHA和前一新commit SHA；新C4 commit接在新C3之后。新helper必须识别这些新状态，旧helper不参与新链判断。所有封存只创建一次，不覆盖已存在封存；中断的未封存数据标为待审查，不自动当新任务重算。

C3重建审计只评价保存密度，不运行BIE线性求解。现有审计函数会写audit文件，因此输入须复制至新审计工作区，输出保存在 `migration_audits/C3/`，不得指向原C3或覆盖其副本中的legacy audit。C0–C2来源审核不满足时停止；不擅自把全套重新审计/求解当恢复后备路径。

迁移结束保存 `migration_checkpoint.json`：原包状态、原始失败、策略批准和SHA、各导入阶段的新状态、`migration_new_acoustic_solves=0`、历史已完成求解3、C4唯一允许动作。审核与策略迁移通过不等于R1数值完成。

## 5. 仅执行 C4 的新任务状态机

未来准备入口缺省dry-run；实施批准与手动启动仍分别记录。上传和启动继续使用`matlab-ubuntu-new`，但实施/传输前另按远程Skill执行授权流程。本轮不加载远程Skill、不连接，也不提供可误用的现成旧Resume命令。

1. 新launcher验证manifest、origin证据和迁移链、任务目录/独占锁/旧会话PID、版本及运行配置；失败在进入C4前停止。新包拒绝旧C1–C4通用run/Resume模式。
2. 新任务的`control/initial_available_bytes.txt`从已验证原初始文件显式复制并记录来源；运行`starting_available_bytes`单独读取当次MemAvailable，`RAM_baseline_bytes=max(initial,starting)`。两字段绝不互相替代；恢复仍采用同一initial和新starting。Bash整数验证、Python JSON序列化及立即解析验证、12 GiB准入、30%RAM看门均原样保留；缺件/非法值/系统读取失败在MATLAB前停止。
3. MATLAB仅核验/跳过新链中C0–C3的有效迁移commit，**代码中不提供C0–C3求解分支**。C4固定L70/ppw12/qa96，GL8，4480节点、560面板、476个平台中点、411真实接收点，原缓存true。
4. C4保存data.mat（密度/完整反射场/入射/数值配置）、points.csv、config、timings、资源、diary、solved封存。调用原阶段审计，四个敏感性均<1e-4；记录既有非单调warning，不由warning抛错；通过后封存C4 commit。
5. C4已存在且有效commit则只最终化；仅有效solved则只审计；孤立/损坏数据停止。密度求解若失败不自动重试，不创建C5、不改ppw/qa/L；未形成有效solved的中断是否重算须另行授权，不能把“恢复”变为隐藏第二次求解。
6. 最终化使用新的attempt/work目录，调用未改的R1 `analyze`、综合独立审计、v2最终审计和原PNG/FIG重建；所有门通过才封存最终产物和完成状态。最终化失败可在新attempt重做后处理，禁止重求解C4。

计数分开：原会话历史`new_acoustic_solves=3`不改；新会话初始0，首次成功C4保存后1，审计/最终化恢复会话可为0；新任务另存累计C4 solve-invocations和completed-solves，绑定attempt，禁止成功阶段重复调用。最终历史完成密度求解合计4（C1–C4；C0复用），不把C0旧资格回放记成新求解。

成本只作等待预期，不新增硬时间门：已测C3同为4480节点，函数108.64 s；C4面板数相同但物理范围改变，不能承诺相同耗时。粗略预留C4函数约2–3分钟，含MATLAB启动、迁移审计、C4独立审计、最终核/图件约5–10分钟或更久，负载和I/O未知。内存沿用12 GiB增量准入估计与持续RAM门，不将其称为已测峰值上界。

## 6. C4 后的经验不确定性与411点误差更新

令 `u0…u4` 对应C0…C4。最终参考仍为 `u3`（L60/ppw14/qa96），C4是固定ppw12的窗口检查，**不切换为C4参考**。

```
bunc = abs(u0-u1) + abs(u1-u3) + abs(u2-u1) + abs(u1-u4)
punc = abs(PS3-PS4) + abs(PS4-PS5) + abs(PS5-PS6)
func = abs(FK1-FK2) + abs(FK2-FK3) + abs(FK3-FK4)
R_m  = 5 * norm(bunc + model_unc_m) / norm(u3)
E_m  = norm(u_m-u3) / norm(u3)
resolved_m = E_m > R_m
```

`u_PS=PS6`，`u_FullFK=FK2`，`model_unc`分别为punc/func；与现有代码完全一致。bunc/punc/func单位为Pa，E与R无量纲。这里的“完整”是预定四个敏感性贡献全部齐备，不是严格误差上界、统计置信区间或所有系统误差。非单调项照常以绝对值进入bunc，不删除、不取差、不因低于1e-4归零。

输出comparison.mat、receiver_fields_and_errors.csv、model_metrics.csv、independent_reference_resolution.csv、逐点参考不确定性表和policy诊断表。411点都保存复残差实/虚部、绝对值Pa、幅度差Pa、相位差rad和按C3全线RMS归一化残差；相位仍用`abs(u3)>=0.01*norm(u3)/sqrt(411)`掩码，报告有效点数。报告全线L2、峰位置与分段/沿线误差集中区域；保存PNG与FIG及原绘图审计。

最终审计独立从五组保存MAT复算四配对、bunc、PS/FK贡献、E/R及掩码，与CSV/MAT核对；政策SHA必须与阶段一致。数值公式保持原样，独立复算不只信任阶段门自报。新最终报告补充：旧Nonmonotonic失败仍有效、政策调整获批时间、迁移来源、哪些审计继承/重做、C4新增求解数及未分辨项。

现有`resolved`是判别模型差是否超过5倍经验分辨率的**结果标记**，并无`assert(resolved)`；保留该规则。`resolved=false`不自动改为程序失败，也不写成“模型一致”；只能报告当前经验分辨率下不能分辨。`reference_qualified`亦须以全部数值/独立审计实际通过为前提，不能凭`analyze`模式变量先置true就宣布资格完成。历史Gaussian仅作同表面/坐标/原生相位、各自声源BIE尺度的条件配对，不作绝对Pa跨源排名。

## 7. 拟修改/新增文件清单（本轮尚未实施）

| 文件 | 拟议职责 |
|---|---|
| 新 `scripts/validation/run_r1_bie_c4_continuation.m` | 只导入C0–C3、唯一C4求解分支、恢复/最终化与新状态回执 |
| 新 `scripts/validation/support/migrate_r1_bie_checkpoints.m` | 只读核验旧链、复制/新审计工作区、创建新迁移链；无求解调用 |
| 新 `scripts/validation/support/r1_bie_continuation_io.m` | 独立schema/scope、来源/policy/前序SHA绑定与新恢复状态；旧IO不改 |
| 新 `scripts/validation/support/r1_bie_continuation_gates.m` | v2配对门与完整单调诊断，保留全部数值门和配对定义 |
| 新 `scripts/validation/support/audit_r1_bie_continuation_final.m` | 复用原最终审计公式，独立复算并核对v2政策；替代新包中的旧严格单调最终审计 |
| 新 `scripts/validation/support/prepare_r1_bie_c4_package.ps1` | 旧证据核验/复制、新manifest/policy/批准记录/归档；不自动上传启动 |
| 新 `scripts/validation/support/launch_r1_bie_c4_continuation.sh` | 只允许新包入口，保持锁、RAM、stdout/stderr、真实退出码；复用修复后的runtime组件 |
| 新 `scripts/validation/support/manage_r1_bie_c4_continuation.ps1` | 新包上传/手动启动/状态/回收命令；不调用旧通用Resume；可复用原archive工具 |
| 新 `scripts/validation/validate_r1_bie_c4_preparation.m` | 保存场/合成状态测试、声学哨兵和零求解迁移dry-run；不启动正式C4 |
| `scripts/README.md`及新实施报告 | 实施获批后登记入口、命令、政策变化、恢复限制和审计口径 |
| `reports/active_task_tracker.md`、`PROJECT_CONTEXT.md`、进度JSON/生成视图 | 各实际边界如实记录：方案/准备/迁移/C4/最终审计分开；不预报完成 |
| `.gitignore`/`.gitattributes`（必要时） | 忽略新原始包/回收副本，仅显式保留小型回执；保留SHA证据字节，不改旧文件 |

不修改PS、FullFK、生产WAPE、原/缓存BIE求解器、原阶段独立审计、原R1分析/综合审计/绘图公式；旧manual入口、旧最终审计和旧包继续代表旧政策。新政策共用不可变描述及语义；独立最终审计仍单独计算指标，以免共用一个错误计算函数形成循环自证。

## 8. 实施后的准入测试与停止条件

以下是**批准实施后才执行的验收设计，本轮未运行**：

- 合成场覆盖单调下降、相等、实际小幅上升，确认v2在四配对合格时均记录诊断而不中止；其中上升例在旧v1仍按原规则失败。另测任一配对=1e-4或更大、非有限、缺配对/错坐标仍失败，不通过“改warning”吞掉真正数值失败。
- 阶段与独立最终评估使用相同policy SHA，四配对按名称/方向对应；PS/FK/C3参考、bunc和5倍分辨公式/严格比较保持一致。既测试`resolved=true`，也测试`false`和信息未齐全，后者必须保持unknown/unrun。
- 用保存C0–C3完成本地无求解迁移：C3只重建审计；求解器哨兵应记录0调用。测试篡改旧MAT/CSV/审计/manifest、前序链变化、缺C3审计、未获批准、旧scope冒用、新policy漂移、半途迁移及重复封存，均停止并保留旧目录SHA不变。
- 独立测试新C4状态机fresh/SOLVED仅审计/AUDITED仅最终化/孤立数据停止；测试锁重复启动、退出码与资源门。哨兵用合成产物验证最多一次C4分支，不把mock称为真实数值PASS。
- PowerShell/Bash静态检查、JSON往返及非法内存输入测试，含available大于/等于/小于initial；核实RAM字段独立、序列化可靠。测试固定依赖SHA、归档/上传回读SHA与代码路径，审查后再交用户启动。

主要风险：后验调整政策可能产生选择性验收；用显式批准、v1失败留存、四配对及完整不确定性防止隐藏。新旧签名混用可误认历史结果，故用双链和独立schema；C3审计原来不在solved四文件绑定中，故必须锚定回收清单并在新工作区重审计。C4可能出现窗口不收敛、RAM/审计/图件失败；必须停止，不扩网格/窗口/种子或频率，也不将未执行检查改PASS。

审批对象是本方案中的**v2诊断政策、零求解迁移与仅C4续行包实施**。即使实施准备通过，也不自动取得SSH传输/正式启动或R2许可；后续明确用户启动方式和动作范围。当前停止，等待用户审查。
