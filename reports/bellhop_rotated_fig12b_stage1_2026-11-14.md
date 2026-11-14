# Bellhop 旋转 Munk Fig.12(b)：S0/S1 阶段记录（2026-11-14）

**当前仅直接传播先导与直接 Q overlay 获得限定证据；最终目标保持 active。没有反射场、完整总场或论文复现结论。**

独立工作树 `C:/Users/ASUS/.codex/worktrees/bellhop_rotated_fig12b/Explain`，分支 `experiment/bellhop-rotated-fig12b`，源基线 `834d4a3375f64b2a8024dcd7c61bb2796d9e8c26`。原 FFP 工作树的未提交方向性研究、a28b 分支、生产 PE、官方 Bellhop 和历史结果均保留。未查询科研 GitHub 仓库或执行 fetch/pull/push。用户随后明确授权获取并安装本地工具链，这是仅本地约束的工具链专用例外。

[阶段计划](bellhop_rotated_fig12b_plan.md) · [阶段图册](../results/visualization/bellhop_rotated_fig12b/20261114_s1_02/index.html) · [数值审计](../results/validation/bellhop_rotated_fig12b/20261114_s1_audit_01/checkpoint.json) · [图件审计](../results/validation/bellhop_rotated_fig12b/20261114_s1_figure_audit_01/checkpoint.json)

## 模型与源码合同

物理参数直接来自 `make_munk_5000m_ffp_config`：5000 m、源4950 m、10 kHz；Munk c0=1500、轴/尺度1300 m、epsilon=.00737；示例底 cp1600 m/s、rho1.8 g/cm³、.8dB/波长。底质、源模型、吸收和归一化不是核实的论文参数。最终参考是三维全向点源切面，不能用二维线源替代。

向上 `(r,q)=(4950-z,q0+x)`；向下 `(r,q)=(z-4950,q0-x)`。两者行列式均+1。原 Munk c(z) 写为 Q 的 c(r)，人工横向深度两行完全相同，.ssp 距离单位 km。人工上下边界距源1000 m；本次所有射线距边界不足0.00438倍，未触边。人工边界不是实际海面/海底。

源码检查发现原 PM overlay 强制均匀1500 m/s；Step 已具备 Q 分段限步/梯度跳跃处理。新生成器只在显式 mode2 开启**直接 Q**；原 mode0/1 保留均匀介质门，反射Q明确拒绝。未以删除旧门冒充反射介质变换。生成的 Step 与旧 overlay SHA 完全相同。官方库仅复制到隔离构建目录，未覆盖。

## 已执行数值证据

| 检查 | 设置/结果 | 范围 |
|---|---|---|
| 轴向联合先导 | 上/下各3配置；Q10/5/2.5m，步2/1/.5m，101束，4个接收距离 | 24检查通过；不代表完整场束数收敛 |
| 最细轴时延 | 向上最大1.06787148e-7s；向下7.45851395e-10s | 原门 .05/(2π·10000)=7.95774715e-7s |
| 独立改变Q/步长 | [5,2]、[2.5,2]、[2.5,1]、[2.5,.5]m，上下各4例 | 32检查通过；ARR存在单精度时延输出精度限制 |
| 非轴射线 | 上±.05°、下±5°及轴，各3条；Q2.5m、步.5m | 连续Munk/Snell积分最大路径差3.69012927e-8m、1.68044245e-10m |
| 横向慢度守恒 | 上5.48835435e-14、下3.16175715e-15 s/m | 初始门较宽2e-6；不将此门包装成严格通信相位资格 |
| 独立保存数据审计 | 实际Q线性段用log1p解析积分；独立连续Munk表达式；ARR/CSV/MAT及SHA | 98检查通过，14个直接到达案例 |
| 新binary对官方 | 2组C场（均匀/Munk）复差1.0837e-17/0；上下射线路径与时延差0 | 仅直接Q；不提升旧Bellhop资格 |
| 源相量桥 | 均匀Green复差3.55277246e-8；Munk相位差.00175381198rad | 5项独立保存数据检查通过，见下节 |
| 阶段图件 | 2组PNG/FIG，离线HTML、绘图源MAT；FIG原始数据重核，PNG逐像素重导差0 | 仅几何/轴时延；不是最终声场图册 |

初次 Q 运行在终止域边缘多迈约1e-5m，Bellhop Quad 报 fatal，但进程退出码为0。现有日志门正确拦截；失败保留在 `20261114_s1_direct_01`。修复只在 Beam.Box.r 外增加常值数值保护段，没有改变物理接收点或借外域续传播。另一独立配置轮在RAM约0.39GiB时未派发；`20261114_s1_independent_01`保留，门未降低。RAM恢复后使用新目录02继续。

## 相量与三维扩散：不可混淆的两件事

源码 `InfluenceGeoHatCart -> ApplyContribution` 使用 `exp(-iωτ)`；`ScalePressure` 对几何点源取 `const=-1`，再除 `sqrt(r)`。因此均匀直接输出是 `-exp(-ikR)/R`。FFP 的 `exp(-iωt)` 与正Green合同对应固定变换：

`pressure_ffp_convention = -conj(pressure_Bellhop_geometric_point)`。

这是从源码和解析均匀1/R独立控制确定的源/相量合同，不是按 Munk 或论文结果拟合增益/相位。最初对比遗漏负号导致两项约π相位失败；`20261114_overlay_validation_01`及其失败CSV保留，独立 `20261114_phase_bridge_01`在保存原始复场上审计修正后的合同，不重算或改写原失败。生产 PE 的phase contract v2未改。

对旋转后只依赖r的介质，轴向两条独立横向Jacobi长度均与 `I=integral c(r)dr` 相关。平面分层三维几何输运预测幅度 `sqrt(c_source*c_receiver)/I`；原生R扩散给 `sqrt(c_receiver/(r*I))`。二者比例 `sqrt(c_source*r/I)` 在所测Munk轴点为1.00013658/1.00531807/1.01006273/1.01030184。当前输出**未乘此比例**，不将原生R场称作已与FFP共源归一。还需独立审查非轴及反射后的三维扩散和beam width，不能把1%差当可忽略的任意常数拟合。

## 便携工具链与执行

共享 MATLAB MCP R2025b Update6，先附着失败、用户恢复后连接成功。没有第二个MATLAB、SSH或parpool。安装使用用户后续授权：[MSYS2官方Fortran包](https://packages.msys2.org/packages/mingw-w64-ucrt-x86_64-gcc-fortran)及运行依赖，22包共89221591字节，逐包官方SHA-256匹配；GNU Fortran16.2.0。位置 `C:/Users/ASUS/.codex/tools/msys2-ucrt-bellhop-20261114/installation/ucrt64`，仅构建进程临时PATH。先尝试的[WinLibs](https://winlibs.com/)整包在300s下载超时（退出28），未解压/执行，部分包和日志保留；安装成功来自MSYS2，不冒称WinLibs安装成功。

独立构建14.2918s；禁用原构建脚本的fast-math，保留NaN/Inf语义，只编译隔离副本。官方exe SHA `7e7809a64c3bf734aff6d28d0d4d52b1b4bd203d81676e3241ffd3189941b505`；新直接Q exe SHA `925302d94c50ce6c5a9710f1cd7bd372cec197109e32e246fd0a33fe135d79c5`。日志、真实退出码和构建清单在 `20261114_build_direct_01`。编译成功本身未被记为物理PASS。

S1轻量数值作业/审计/构建目前有计时记录约47s；图件v2渲染8.8651s、独立图件审计8.3912s另列。旧900s FFP预算不转用；原图v1和代码准备早期异常未取得完整计时，完整任务总耗时不作伪精确声明。接下来只允许必要小型局部开发；大型场/远端必须另列预算并获授权。

## 文件/接口与待办

新增验证/支持入口、阶段渲染与报告，导航已更新；未改任何 `src/`、公共信道API、PE或通信方法，所以本阶段没有重新执行生产标量/宽带通信回归。现行限制完全保留。

复现顺序：`setup_vertical_project`；`validate_bellhop_rotated_munk_direct`；`validate_bellhop_rotated_munk_rays`；`audit_bellhop_rotated_munk_direct`；`prepare_bellhop_rotated_q_overlay`生成隔离源码/PS构建入口；按生成脚本显式传OfficialToolboxRoot、CompilerRoot、ProjectRoot和新的OutputDir；`validate_bellhop_rotated_q_overlay`；`audit_bellhop_rotated_phase_bridge`；`render_bellhop_rotated_munk_preflight`；`audit_bellhop_rotated_preflight_figures`。所有入口均接收绝对路径，不依赖cd。当前准备器是直接Q第一步，不是完整反射实现。

尚缺：反射处单侧Munk梯度和p/q、流体底、相量统一后的三维扩散、反射阶数/束数/全复场收敛、FFP共点数据与来源复核、全域/近底各分量、最终图册及六问结论。最终GO/CONDITIONAL_GO/NO_GO未授予；目标保持active，不合并。阶段图件遵循本地 MATLAB/scientific-visualization 技能的数据保真、FIG保存与目视核查工作流；不将工具审计等同物理资格。
进度看板在本独立工作树渲染为 GENERATED_WITH_REFERENCE_WARNINGS：57项引用中16项既有未跟踪结果未随Git工作树提供。均为其他工作线/旧Bellhop保存结果，不是本轮数值失败；完整缺失清单保留于看板回执，不改expected_exists、不从cash补齐、不重跑无关实验。本轮新增报告/检查入口可用。
