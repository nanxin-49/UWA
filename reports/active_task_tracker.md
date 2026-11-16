# 当前任务接续（本独立Bellhop分支）

此文件只导航本工作树的当前任务，不授予计算许可，不覆盖主线0.10m二维条带顺序。

| 工作线 | 日期/设置 | 当前状态与证据 | 下一步及停止门 |
|---|---|---|---|
| Bellhop C工况远程声线补齐 | 2026-11-16；5000/4950m、10kHz、点源宽N/A、平面海面/流体底 | 6类18条真实轨迹/1299038节点；18映射、102FIG、298返回SHA及4发布SHA通过，接收声压求值0。[报告](bellhop_C_rays_remote_complete_2026-11-16.md)；[CP](../results/validation/bellhop_rotated_fig12b/20261116_ray_remote_03/checkpoint_delivery.json)；[图册](../results/visualization/bellhop_C_vs_FFP/20261116_rays_complete_01/index.html) | 补齐后停止待审。独立R导出副本保留完整历史，原C ELF不变；原3600+本次120=3720s，余7.860s。ABC/FFP重算0，无D/提交/推送/合并。下列早期停止按历史保留。 |
| Bellhop旋转Fig12b接收距离精度 | 2026-11-14；10kHz；5000m；平坦海面/流体底；全向三维点源，源宽N/A | REAL8读取/1000d0和ReceiverRange64实际链最小验证通过；10981原点，13阶101束3列；25+31+60项检查通过。[中文报告](bellhop_rotated_receiver64_minimal_2026-11-14.md)；[运行](../results/validation/bellhop_rotated_fig12b/20261114_receiver64_minimal_01/audit01/checkpoint.json)；[独立审计](../results/validation/bellhop_rotated_fig12b/20261114_receiver64_minimal_01/saved_audit01/checkpoint.json) | STOPPED_FOR_REVIEW。完整A、第14–16阶、B/C/D不恢复；原25分量和旧FATAL保留；完整场/收敛/FFP资格未取得。 |
| Bellhop旋转模型性能审查 | 2026-11-15；10kHz；5000m；平面墙2段；全向三维点源、源宽N/A | 独立保守排除优化，13阶101束3列；1313次反射及最终声压逐字节相同，234项保存审计通过；一次无采样配对2.4268952→2.1920121s（1.107倍）。[中文审查](bellhop_rotated_performance_review_2026-11-15.md)；[独立审计](../results/validation/bellhop_rotated_fig12b/20261115_performance_review_01/saved_audit01/checkpoint.json) | TIMING_LIMITED_RAM_STOP，重复计时未补齐；细粒度占比不作精确归因。停止待审；不恢复A/14–16/B/C/D，不提交/推送。 |
| Bellhop A组恢复与阶段验收 | 2026-11-15；5000m/4950m/10kHz；全向3D点源、源宽N/A；平面边界 | 39/39正式分量，31复用+8新增；原25 SHA保持。61门原始重建、5区N12→16和海面残压门通过；初步3组PNG/FIG，独立FIG69门通过。[报告](bellhop_rotated_A_complete_2026-11-15.md)；[阶段CP](../results/validation/bellhop_rotated_fig12b/20261115_A_resume_01/checkpoint_stage.json) | A_COMPLETE_STOPPED_FOR_REVIEW。B/C/D预算估计已给，原剩约2687s，上界可能不足；不自动启动，不提交/推送/合并。最终FFP及论文复现未取得。 |
| Bellhop B组Q加密及A/B验收 | 2026-11-15；5000m/4950m/10kHz；全向3D点源、源宽N/A；Q.625m | 39项完成，A39及原25 SHA保持；727组件/61原始重建、B截断及A→B五区均通过；79图审通过。[报告](bellhop_rotated_B_complete_2026-11-15.md)；[阶段CP](../results/validation/bellhop_rotated_fig12b/20261115_B_run_02/checkpoint_stage.json) | B_COMPLETE_STOPPED_FOR_REVIEW。余约1813s/6.58GiB；C/D/FFP未启动，不提交/推送/合并。可选文件记账与图册路径修正保留记录，无原生重跑。 |
| Bellhop C组步长减半 | 2026-11-16；5000m/4950m/10kHz；全向3D点源、源宽N/A；Q.625m、h.5m | C39/39，Windows30复用+Ubuntu9新增，native204.976s；18资格/727组件/195输入/61重建，截断及B→C五区、海面/101图审通过，543返回SHA一致。[报告](bellhop_rotated_C_complete_linux_2026-11-16.md)；[CP](../results/validation/bellhop_rotated_fig12b/20261116_linux_C_01/checkpoint_delivery.json)；[ABC图册](../results/visualization/bellhop_rotated_fig12b/20261116_ABC_linux_01/index.html) | C_COMPLETE_STOPPED_FOR_REVIEW。A/B78及原25保持，原Windows停止记录未覆盖；余155.83s。D/独立FFP和论文最终资格仍缺，本轮不启动、不提交/推送/合并。 |
| Bellhop C远程迁移 | 2026-11-16；gfortran13.3、MATLAB R2025b Update1 | 复用2736文件已校验副本，不重复迁移；独立Linux REAL8+平面优化编译12s，最小101束3列18门通过且复声压零差。首次辅助函数与Java8 PID接口失败保留。[C完成报告](bellhop_rotated_C_complete_linux_2026-11-16.md)；[原拷贝报告](bellhop_rotated_remote_transfer_2026-11-16.md) | Linux仅任务独立目录；C阶段已完成，保留旧复制/CR警告及原停止CP。停止待审，不D/新FFP或Git提交。 |
| Bellhop C与已有FFP共点比较 | 2026-11-16；10kHz、5000m/4950m、全向3D点源N/A；原10981×201网格 | 16可比性/19FIG数据/39独立检查通过；全复L2=2.9902e-5，max复差/E=8.7220e-4，max有效相位.0121599rad；6组PNG/FIG/HTML。[报告](bellhop_C_FFP_common_point_2026-11-16.md)；[CP](../results/comparisons/bellhop_C_vs_FFP/20261116_post02/checkpoint_delivery.json)；[图册](../results/visualization/bellhop_C_vs_FFP/20261116_post02/index.html) | 比较完成，完整反射声线RAM_STOP：R-only首项启动前<1.5GiB，PID0，不重试。只复用真实直达.ray与13阶命中标记，不冒充完整轨迹。余60.27s；不D/FFP/ABC重算或提交，等待人工审查。 |

最后动作：2026-11-14，完成唯一授权最小求解（新进程1个，2.5999938s）及无求解保存字节/原分量SHA审计，更新中文报告与进度看板，等待人工审查。此前RAM停止、输入FATAL及准备阶段失败均保留。
版本封存：2026-11-14按用户要求提交REAL8最小验证源码、中文报告及精简审计记录；早先未采用的标签修补草稿留在工作区，原始MAT/SHD/构建树仍本地保留。本次Git操作不改变人工审查停止门。

最后动作（2026-11-15）：完成性能审查、6个最小案例（累计14.823099s）及无求解234项独立审计；重复计时按1.5GiB RAM启动门停止，下一例无PID。保留REAL8修复、原25分量与全部旧失败，更新报告和看板后停止。尚无完整A恢复许可，未提交/推送。

最后动作（2026-11-15 A完成）：用户批准A-only恢复后完成8项正式计算，累计新native100.3129s；39项组装及原门验收通过，原25 SHA不变，图件检查器假失败另存并由69门独立审计关闭。保存B/C/D规划后停止人工审查，不提交/推送/合并。此前停算说明按各自日期保留。

最后动作（2026-11-15 B完成）：用户仅授权B后完成39项串行，native457.43s；A39与277项清单及原25 SHA一致；B截断/A→B各五区通过，图册迁回目标树后79项独立图审通过，保存收尾问题记录。停止待人工审查，无C/D/FFP或提交。

最后动作（2026-11-15 C停止）：按C-only授权串行完成30项，native270.728s；down_order12启动前RAM门拒绝，无PID。封存复核A/B全部清单及78分量和原25 SHA保持，C413绑定SHA通过。未执行9项、完整C/收敛/ABC图，停止待审；保留输入保护点和封存CSV假失败记录。


最后动作（2026-11-16 C恢复前再次停止）：用户授权剩9项后，先核A/B/C SHA；随后RAM启动门再次失败，未派发Bellhop，C30/39保持。共享会话未见可安全清除的大型任务临时变量，未清理或终止用户会话/进程；保守计账后余859.17s，保存新停止证据及旧CP，等待人工审查。


最后动作（2026-11-16 远程迁移预检）：SSH和MATLAB健康检查通过，9返回文件SHA一致；远端缺gfortran，按用户指令报告并停止安装/迁移/编译/求解。Windows结果和旧CP未覆盖，C30/39，迁移账本余799.17s，不D/FFP或提交。


最后动作（2026-11-16 拷贝完成）：按用户要求帮拷至bellhop_rotated_C_transfer_20261116_01/Explain，2736文件SHA全过，路径表保存；未启动计算。用户要求尽快简核，未重复全量审计；收尾CR与可选环境日志缺失如实记录。



最后动作（2026-11-16 Linux C完成）：原副本复用，编译/跨平台18门通过后串行完成缺失9项，C39项组装与全部原收敛/海面门通过；101图审、543返回SHA、A/B78及原25复核通过，轴图例副本数据不变。最终原预算余155.83s，停止人工审查，不D/FFP，不提交/推送/合并。

最后动作（2026-11-16 C/FFP比较）：已有J07与C原网格和源/底参数核对一致，完成无拟合共点指标和统一色标全场/中央/近底/轴图；近底一致、上部局部残差较大，均不支持贯通中央亮带。准备完整反射R轨迹因RAM门启动前停止，0PID；历史直达与反射命中点明确区分。16/19/39检查通过，原预算余60.27s；不自动重试、不D/新FFP、不提交。


最后动作（2026-11-16 声线补齐）：远程独立R导出版本完成6类18轨迹；复用直达、新5项反射，18映射/102FIG/298返回SHA全过，两PNG目视检查。预算停止与裁断导出失败保留，原C ELF一致，无接收场求值；累计3712.140/3720s，余7.860s，停止人工审查。
最后动作（2026-11-16 Git封存）：用户新授权提交与推送，按性能/元数据、ABC验收、共点与完整声线三部分封存；当前仅准备版本记录，实际提交与推送结果见[封存说明](bellhop_rotated_publication_2026-11-16.md)及对应push_receipt.json。科学状态和资源账本不变，不授予新计算许可。
