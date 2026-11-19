<!-- GENERATED: edit reports/project_progress.json, then run render_project_progress_dashboard. -->
# 从这里接着做

垂直水声项目 · 仿真与验证进度

记录更新：2026-11-19；本次核查：2026-11-19。

K1七条件保存审查限定完成：30层BIE、风速/截止交互齐全；原生退出异常保留

## 现在在哪里

这个项目想回答的是：海底附近发出的声音，经过直达和海面反射到达接收机后，会怎样影响通信？为此，一边需要可信的声场和信道，一边需要能检验信道估计与均衡方法的通信实验。

你离开前，已经完成二维单频比较，做过三维固定海面的独立频点检查，还在一个接收点上跑通了受控的宽带与四相位调制（QPSK，每个符号承载两比特）通信实验。因此，重新开始时有不少结果可以接着用。

P22保存结果与新Ubuntu最小项目基线保留其历史复核范围；M1仍停在RAM门，正式S5/S6仍暂停。当前前向研究顺序以AGENTS.md和reports/active_task_tracker.md为准：10 cm条带4 kHz PM → 反射模型决策 → 相干4–8 kHz频响/CIR → 受控通信；统计、时变、多普勒及1500 m扩展留待另定范围。

2026-11-10 M0：生产原生PE到DSP改为完整包络共轭，统计/迁移/外部桥接同步。P22平面波高度、Gaussian同源角谱、C/P、SSA、CIR及小型已知CSI无噪声通信通过；旧失败和历史结果保留。

 本轮追溯原始窗口证据并冻结M1：不重扫失败小窗口，采用192.1875m起点、独立AS参考及逐步边缘检查。

2026-11-11：本轮M1执行在MCP连接前置门停止，未取得窗口、网格、步长或吸收层新资格；没有合格生产推荐配置，不能进入M2。

历史方案（2026-11-11）曾选20 cm代表尺寸，安排共源、4 kHz试点及三频比较；这些原始源条件和结果继续保留。2026-11-14起后续主源总宽改为10 cm（半宽5 cm）的理想二维无限条带，不能当作已测实物三维有效口径；不重标20 cm结果，也不接管M1预算。

同日后续M1：既有MCP版本/PID核实后，完成P4/P6和10组独立AS；已完成对照独立审计0失败。P8在RAM跌破启动可用量30%时中止，Z4/Z8和全部吸收层对照未派发，无生产推荐，不进入M2。

2026-11-11后续活塞图示：旧5 cm经GitHub和本地代码确认是圆形半径，本次以二维总宽10 cm映射半尺寸，并加总宽5 cm控制。20/10/5 cm同振速三频自由场已保存复Pa、图、独立检查；这不关闭P0，后续先共源压力/梯度/相位资格，再平面镜像。

2026-11-15 R1：10 cm、4 kHz、PM seed260001、411点反射场已保存；PS/FullFK数值检查通过，BIE首层残差通过，但剩余收敛层按本机成本预计超预算而未派发。暂定复L2为15.9235%和0.0747208%，参考未合格；R1停止，R2及以后不启动。

2026-11-16实现验证专用默认关闭缓存并执行S0/C0：4组小规模及3200节点历史密度/复场差全部0。C0调用55.48 s，原N³成本门预测余层730.29 s高于582.80 s余额，C1–C4未派发；累计525/1100 s，R1参考仍未合格。

2026-11-16本次仅准备用户手动C1–C4批处理与新执行政策：取消旧累计/成本预测停止，单个R2025b单线程、无parpool、30%RAM与原数值门保留；已通过本地静态与声学入口哨兵dry-run。C0不重算、已审计阶段按SHA跳过，真实C1–C4尚未运行。

2026-11-17已准备有限R2的seed260002/260003独立执行包；支持两Seed各一单线程进程和共享30% RAM预留，实际结果尚未计算。当前停止等待用户手动启动，不能将准备完成当跨Seed结论通过。

2026-11-17新增PM联合敏感性设计：三U×五N×三Seed，原12模Gaussian可恢复且U6原值保留；新扩展流明确记录。正式45组K0尚未运行，K1声学入口硬禁用，等待用户手动K0，不持续监控。

2026-11-18：用户回传K0正式45组；574项SHA和独立本地MAT/CSV/原FIG117曲线通过。固定风速扩大物理截止时高度变化小也可能伴随显著导数变化；K1保留seed260001七条件，尚未执行声学。

2026-11-18 D0：不同发射深度入射复场06包准备完成，正式声学未运行；用户手动启动/Fetch，1500m背景不引入海底反射，生产/R1/R2/K0不改。

2026-11-18 D0正式回收与保存分析：五深度4kHz/10cm，150SHA/165远端门/127本地检查通过；1490m x200m幅降.6514%但相位224rad，原100m反射误差不可直接外推。D1及其他新求解未自动推进。

2026-11-18：D1a固定10cm/4kHz、原U6/N12/260001/L160/C2、Tx100/1490m、Rx6m，预注册±36m/411点，每深度独立五层BIE；03包准备完成。D1b仅评估新周期/支持与稠密成本，不平铺旧面，待D1a结果后决策。

2026-11-18 D1：原失败数据已回收，59项SHA保持；17位CSV修复后保存C0独立审计通过。恢复03包导入原100m模型/C0、保持旧时间/RAM基准，等待用户手动新任务启动剩余9层；本轮0新BIE/远端MATLAB，其他工作线独立。

2026-11-18 D1a：双深度各五层独立BIE及PS/完整FullFK、保存数据和本地图件审查完成。PS复L2相对降低78.115%，但1490m仍约FullFK误差50倍。只限0.10m/4kHz/原单PM/C2/Rx6m/±36m；0–200m长海面D1b只有可行性方案。

2026-11-19：K1按用户直接目标授权完成六新面独立参考（原N7五层导入＋25新BIE），严格R1锚点单列。固定0.10m/4kHz/U_PM/seed260001/原C2/411点，4040SHA及本地完整保存审查0失败，新128图曲线/56FIG/112次零像素和实际预览通过，14模型行/5754点。PS10.251–53.845%、完整FullFK0.02083–0.19553%，风速/物理截止交互经验可分辨；N40不是完整PM。四并发实测峰9.96GiB、累计2656s/10800s，30%RAM保持。最后一例封存后原生退出−9和控制器78保留，独立保存审查通过后接受结果，不重算、不称退出环境稳定。任务案例进程已退出，停止于K1有限范围。

最新D1a已经完成：此前“待启动/停止/图件待审”是历史状态，失败保留。当前依据最终报告和542SHA/本地独立审计；其他工作线按各自范围。 最新D1已完成精度修复和保存C0独立审计；原首跑失败保留，恢复作业尚未启动，双深度参考仍未完成。此前“待修复/全部未运行”只适用于其历史时点。 最新D1已由用户手动运行并STOPPED；不再是全部声学UNRUN。100m模型自身门通过、C0已求解但审计未封存，1490m未运行。以下准备文字保留其历史时点。 最新D1只准备未运行：所有MATLAB、声学、FIG和Linux原生门仍未执行，等待用户审阅授权后手动Start。以下D0/K1准备文字保留其历史时点；D0正式结论不被D1准备升级。 D0正式入射场已回收审查，准备时UNRUN记录保留为历史。当前研究结论与8张新图见D0正式报告；K1等工作线状态不由D0推进。 D0本轮仅准备完成，正式入射场尚未计算；下一步用户手动启动/Fetch。 K0远端运行日期11-17，正式回收审查日期11-18。当前K0限定完成，K1方案选定但声学执行禁用；物理Kmax与采样dx分别审计，采样通过不等于物理带宽收敛。R1/R2历史资格与失败保留。

## 各条工作线

### 二维单频比较

状态：K1七条件与D1a限定完成；退出异常保留；结果日期：2026-11-18（K1准备；正式声学未运行）；2026-11-18 D0远端计算及本地保存审查完成；2026-11-18 D1仅准备/成本评估，无新声学结果；2026-11-18 用户D1首跑STOPPED；2026-11-18 D1保存C0独立补审通过/恢复包准备完成；2026-11-18 D1a正式限定完成；2026-11-19 K1七条件保存科学资格完成，原进程退出异常不改判；本次核查：2026-11-19（K1完整保存数值/图件；其他研究线未重审）。

D1a固定0.10m/4kHz/原U6/N12/260001/L160/C2、Tx100/1490m、Rx6m/±36m/411点完成。PS/BIE复L2=16.1917%→3.54362%，完整FullFK=0.0765349%→0.0709202%；变化均分辨。全10次BIE无重复，542SHA、80FIG曲线/零像素通过。K1等工作线状态与失败独立保留。 K1按用户直接目标授权完成六新面独立参考（原N7五层导入＋25新BIE），严格R1锚点单列。固定0.10m/4kHz/U_PM/seed260001/原C2/411点，4040SHA及本地完整保存审查0失败，新128图曲线/56FIG/112次零像素和实际预览通过，14模型行/5754点。PS10.251–53.845%、完整FullFK0.02083–0.19553%，风速/物理截止交互经验可分辨；N40不是完整PM。四并发实测峰9.96GiB、累计2656s/10800s，30%RAM保持。最后一例封存后原生退出−9和控制器78保留，独立保存审查通过后接受结果，不重算、不称退出环境稳定。任务案例进程已退出，停止于K1有限范围。

**为什么做：**先在较容易控制的二维问题里，确认不同方法比较的是同一声源、同一海面、同一接收位置，再看它们算出的反射波差多少。

**做了什么：**以下模型比较分别保留其Gaussian或20/10 cm原始声源条件，不能合并当作10 cm PM资格。PE、边界积分方程（BIE，独立数值参考）和完整 Kirchhoff 积分近似（FullFK）的同源比较已完成。平面和三档正弦海面的九点结果齐全，三个预定 PM 海面实现也完成了比较；相关表格、图和历史审计均有保存。 新增严格二维速度孔径积分三频自由场，保存同振速绝对复Pa、接收线、近场/整体场和方向图；独立积分/角谱与保存数据检查无失败。 远端R2025b Update1及实际双进程已验证，源解析梯度对独立角谱/差分/速度/相位通过；BIE/FullFK输入差0，Gaussian标量默认差0；123项部分产物回传SHA匹配。 许可恢复后P0六组、P1自由/平面/弱C2正弦实际完成；最终独立审计0失败，PNG/FIG/MAT/CSV/日志/SHA回传。 2026-11-15完成本轮R1范围内的PS窗口/网格/步长、FullFK求积/范围、10 cm BIE L60/ppw10/qa96一层及10类已执行保存数据独立检查；四张最终PNG/FIG、逐点复场/幅相/残差CSV、MAT、SHA、日志和停止检查点齐备。 2026-11-16新增默认关闭缓存、阶段计时和阶段/累计时间及RAM中止门；4组小对照和C0完整回放通过，独立C0密度接收重建与MAT/CSV通过，65个回传文件SHA匹配。 最新准备完整MATLAB入口、Linux后台启动/状态/回收脚本、Windows管理及打包工具，包35项SHA、PowerShell/Bash检查和带7个声学入口哨兵的dry-run通过；恢复状态机覆盖fresh、audit-only、skip、前序变化、篡改、孤立数据拒绝。 已静态确认回收C1/C2通过提交与C3已求解审计，参考非单调门停止；未改旧证据。 已完成非单调门调整与仅C4有条件续行的静态设计，明确两处政策一致、原数值门不变、双链封存与C3参考不换；尚未实施。 v2独立代码/政策/来源双链迁移完成，正式04包31项MATLAB、25项完整性/准入、29项runtime及归档重定位dry-run通过；新密度求解0。 修复Status的Windows引号传输及初始化误导；7项原生stdin与12种状态回归和实际SSH只读查询通过，未启动/重启求解。 回收最终287项SHA/旧新C4及最终封存已核；新工作区重做C4 476+19点、源/密度/Green/FullFK双项/MAT-CSV/分辨审计0失败，12条FIG及两组PNG对一致、四图实际查看。 2026-11-17新增R2研究编排/冻结适配/独立最终指标/保存图件及三Seed汇总入口，各28输入拒绝检查、29原runtime、12状态、7 stdin、8本地mock、10 RAM边界和两归档重定位dry-run通过；R1指定回收287 SHA保持。 本次核两份229/239回收SHA和每Seed独立C0–C4链，定位260002最终F02像素门失败；本地36.28s保存密度/边界/源/Green/FullFK/CSV/分辨全部通过，12FIG线完全一致；导出探针及原renderer保存重放完成，原件未改。 本轮R2两Seed各五层保存密度及综合/分辨0失败；新固定布局恢复FIG原画布，5图10次PNG比较差0，60条原/新/汇总FIG曲线一致且5图实际查看。各125/84文件新封存，三Seed正式六行、快照六行、2466逐点与24分段表齐备。 本轮检查风速高度/平均时长无定义，统一U_PM；三个原12模Gaussian恢复、U6原值精确保留、13–40独立固定扩展流，跨U重算Sk与跨N嵌套不重归一化。K0 MATLAB入口/独立保存审计/PNG-FIG-CSV/手动上传启动状态回收齐备，05包70656bytes、15项manifest、23项原来源SHA保持。287项本地准备断言与3个TEST_ONLY保存/采样审计、7传输mock、13离线supervisor、重定位通过。 2026-11-18核K0唯一远端会话、实际退出0/0及全部45组生成/独立审计/提交链；574项SHA前后保持。用户恢复MCP后，本地45组MAT/CSV/Fourier/C2/谱/理论/采样和额外曲率指标审计0失败，原FIG117曲线核实；新图例外置六面板PNG/FIG54曲线重开通过，108截止/36风速/18交互描述表及K1七条件连续净距表齐备。 2026-11-18 D0复用原10cm/4kHz源和phase v2，冻结Tx100/300/600/1000/1490m、平面±220m与K0原±50m粗糙辅助；手动UploadAndStart/Status/Fetch、单线程顺序/30%RAM/独立检查点/异常不重试齐备。34项已读来源SHA保持，全部准备失败保留。 K1资源/许可证检查及15阶段历史成本提取；05包冻结六新面原PS/BIE/FullFK/411点层级，24项准备/32离线/8传输mock与汇总14/5754行及56曲线通过，11文件Analyzer0。 补充2026-11-18 12:21原生诊断：11项锁/身份/真实信号通过，四个轻量子进程全部回收；资源快照可用约39GiB、PSI下降，当时只读准入通过。仍选1×1，正式Start和声学未运行。 2026-11-18 D0正式回收完成：R2025b单线程顺序五深度152.710s/真实退出0/0，最低可用36.114GiB（95.55%基准）。本地44.707s/127保存数组检查0失败，原60FIG线匹配；另存8图70线/重开PNG35.345s，原150SHA前后保持，全部图已目视。 2026-11-18 D1：新增隔离入口、Rx3→6m的PS/独立Green适配、深源分段Weyl、每深度五层独立BIE编排、误差/经验分辨/FIG/SHA及手动工具；静态和离线测试通过，正式未启动。 2026-11-18用户D1首跑及本轮只读排查：两次SSH查原会话、日志、输入/C0 CSV，定位精度断言而非RAM停止；原包/失败现场不改。 K1 N7完成本地原五层独立保存审计与有限N7/N12分析，原回收SHA前后不变。 D1本轮：Fetch原59项SHA保持，D1 CSV统一17位往返，原C0六类独立审计通过（密度1.383e-15、340中点6.971e-8、19点8.290e-8），复场/坐标差0；最终恢复包125SHA、真实MATLAB准备/7 Analyzer0/23离线/12mock通过。 K1三项远端非声学/保存图件诊断结束，恢复04实现双链/N7导入/724秒与单例旧时间继承及逐PS/FullFK事件日志；实际包无求解测试通过。 D1目标模式完成：1490m五层新求解，100m按SHA恢复不重算；542SHA、保存审查113.497s/0失败、80FIG曲线/零像素通过。所有IO/图形失败、原v2非单调诊断和累计时间均保留。 K1按用户直接目标授权完成六新面独立参考（原N7五层导入＋25新BIE），严格R1锚点单列。固定0.10m/4kHz/U_PM/seed260001/原C2/411点，4040SHA及本地完整保存审查0失败，新128图曲线/56FIG/112次零像素和实际预览通过，14模型行/5754点。PS10.251–53.845%、完整FullFK0.02083–0.19553%，风速/物理截止交互经验可分辨；N40不是完整PM。四并发实测峰9.96GiB、累计2656s/10800s，30%RAM保持。最后一例封存后原生退出−9和控制器78保留，独立保存审查通过后接受结果，不重算、不称退出环境稳定。任务案例进程已退出，停止于K1有限范围。

**说明什么：**有限三Seed支持PS复场差明显且相位项主导、完整FullFK误差远小于PS但仍可分辨。空间分布不固定；原图件失败不等于声学失败。只关闭本固定二维单频范围，不宣称总体误差规律、真实3-D或宽带通信资格。 K0准备证据只证明所测软件/配对/审计路径；不能升级为正式45组海面统计或声学敏感性结论。 K0证明本45组表面统计与数值采样证据可信，显示高度对高截止不敏感时坡度/二阶导数仍可显著变化；不证明PS/FullFK声学精度或总体统计规律。 D0当前只证明所测软件准备路径，尚无发射深度声学变化结论；有限窗FFT带外功率不等于物理倏逝能量，Tx100m反射误差不向1490m直接外推。 D0证实该理想化均匀水体中深源在固定窗口幅度更平坦、入射更竖直、曲率下降；全线传播角谱模量不变，窗谱更窄。1490m不能全域平面化，局部近似依赖patch；R1/R2反射误差不可直接外推。 D1当前只证明离线工具路径及方案已具体化；不能判断深度是否改善PS/FullFK，也不能认证无限PM或真实深海环境。 本次失败是CSV序列化精度与严格相等检查的兼容问题，不是已证实的物理或资源失败，也不提供完整D1模型误差资格。 当前仅C0保存审计合格，不等于C0–C4收敛或完整D1误差结论。精度修复不改物理数组或科学门。 D1在本固定有限海面内证实深源显著改善PS，但深源PS仍3.54%、约FullFK误差50倍，不可据幅度平坦认为模型等价。一个配对实现与411空间点不支持总体海况统计推断。 K1在同一标准化随机块和固定数值门下表明风速与原始谱带宽都显著影响有限海面上的近似误差，且模型误差有可分辨交互；只属于本单Seed/二维/4kHz/411离散点，不是总体显著性或完整PM收敛。

**还差什么：**D1a不覆盖0–200m或无限PM，D1b仍需新周期/低频/窗口和资源资格。K1已独立完成当前七条件，但只限单Seed、原有限C2、4kHz；N40未穷尽PM，退出原生崩溃根因及控制器异常分类缺口保持。生产集成、FastFK、宽带/CIR/通信及新Seed不由这些结果自动推进。

**下一步：**本轮K1完成后停止，读取final_summary_20261119_01和最终科学检查点；不要用旧恢复包或原FAILED Status再次Start。后续R3、关键Seed复核或D1b须另定明确范围、资源和数值门；不自动扩宽带或通信。

**技术范围：**固定二维、4 kHz、984 项 Gaussian 发射场；原九点与 411 点接收线分别统计。S4=CORE_COMPLETE，S3A=3/3，S3B 因预算未运行；PHI_PLUS90 的混合平台 BIE 求解差未直接定界。 活塞quicklook：c1500 m/s、rho1025 kg/m³、Tx100 m/Rx3 m、峰值振速1 mm/s，W20/10/5 cm；20 cm三频-3 dB全宽112/67/49°。

**已有结果：**正式结果、表格和图件入口存在。 活塞run 20261111_2d_strip_01：MAT、CSV、8张交付PNG/FIG、日志、SHA和检查点。 本轮20261112_remote_p0p1_01保存复MAT、逐点与粗细接收线CSV、五组最终PNG/FIG、日志、资源和SHA。 R1结果位于results/comparisons/piston_2d_pm/20261115_r1_seed260001_01/，137个远端交付文件SHA匹配，最终四图已实际查看。 缓存新run results/comparisons/piston_2d_pm/20261116_r1_bie_cache_01/保存4组对照MAT、C0原生密度/411点场与CSV、8阶段计时、成本决策、独立审计、真实日志/退出码及SHA。本轮不新增重复物理图。 手动包results/comparisons/piston_2d_pm/20261116_r1_bie_manual_ready_03/；本地无声学验收results/validation/r1_bie_manual_preparation/20261116_ready_03_dry_run/。 正式v2包results/comparisons/piston_2d_pm/20261116_r1_bie_c4_policy_v2_04/；原证据、新迁移/审计/本地准入及归档SHA齐全，C4数据不存在。 最终审查产物results/comparisons/piston_2d_pm/20261116_r1_v2_review_01/，含MAT/CSV、回执、原图PNG/FIG按SHA复制及独立重审计，报告reports/r1_bie_v2_final_review_2026-11-16.md。 R2最终包results/comparisons/piston_2d_pm/r2_preparation_20261117_03/seed260002、seed260003/，各300032 bytes；本地检查results/validation/r2_pm_preparation/20261117_03/。01/02候选及准备失败保留，不上传。 新回收r2_seed260002/260003_retrieved_20261117_01；诊断results/validation/r2_seed260002_review/20261117_01/含完整JSON/CSV/日志/PNG/FIG与新工作区重审。 最终限定R2新封存及图表在results/comparisons/piston_2d_pm/r2_saved_closeout_20261117_01/，最终报告reports/r2_saved_closeout_final_2026-11-17.md；两Seed正式回执/源SHA/代码/失败/PNG/FIG/CSV齐全。 K0正式交付results/experiments/pm_wind_cutoff/preparation_20261117_05/；准备测试results/validation/pm_wind_cutoff_preparation/，TEST_ONLY不混作正式K0。 K0原回收pm_wind_cutoff_k0_20261117_05_retrieved_20261117_232303只读保持；本地审查results/validation/pm_wind_cutoff_retrieval/numeric_20261118_01/，清晰图/派生表results/visualization/pm_wind_cutoff_review/20261118_01/。 D0最终results/experiments/d0_incident_depth/preparation_20261118_06/，本地验证results/validation/d0_preparation/，三份D0方案/测试/命令报告在reports/。TEST_ONLY不作研究数据。 K1原生补测与前后资源、35项回收SHA在results/validation/k1_preparation/native_20261118_07/，正式05归档SHA不变。 D0原回收d0_incident_depth_20261118_01_retrieved_20261118_152730只读保持；审查results/validation/d0_retrieval/，最终8图results/visualization/d0_incident_depth/20261118_review_01/；报告reports/d0_incident_depth_review_2026-11-18.md。 D1准备03包位于results/comparisons/d1_depth/preparation_20261118_03/；离线/几何/成本/源码SHA见results/validation/d1_preparation/20261118_01/；三份D1报告在reports/。 K1只读中断证据results/validation/k1_interruption_review/20261118_01/及reports/k1_interruption_diagnosis_2026-11-18.md；不升级U6/N7为本地数值审查通过。 D1失败诊断及原日志/坐标/SHA快照在results/validation/d1_diagnosis/20261118_01/；报告reports/d1_first_run_diagnosis_2026-11-18.md。 K1保存审查与有限对比表/崩溃公开摘要在results/validation/k1_retrieval_review/20261118_01/，报告reports/k1_retrieval_review_2026-11-18.md。 D1最终恢复results/comparisons/d1_depth/recovery_preparation_20261118_03/，125项SHA/19169280bytes；保存审计/准备测试results/validation/d1_recovery/20261118_01/；修复与手动命令两份报告在reports/。 K1恢复04包results/comparisons/pm_wind_cutoff_k1/recovery_preparation_20261118_04/；A诊断/四次Windows失败/Ubuntu新PNG-FIG及准备测试results/validation/k1_recovery_preparation/20261118_01/。 D1正式回收results/comparisons/d1_depth/d1_goal_completed_retrieved_20261118_04/；本地results/validation/d1_goal_closeout/20261118_01/；最终对照PNG/FIG仅results/visualization/d1_depth_summary/20261118_04/；报告reports/d1_depth_comparison_final_2026-11-18.md。 K1当前正式结果results/comparisons/pm_wind_cutoff_k1/final_summary_20261119_01，独立科学封存results/validation/k1_target_execution/20261118_01/final_checkpoint.json，旧exit和失败归档保留。

**历史审计：**S4 保存证据审查通过并保留限制；本次没有重做数值审计。 本轮仅核旧源代码/检查点中的半径定义与SHA，不重演旧声学资格。 本次活塞保存数据独立重算、边界密度重建、FIG逐曲线及传输审计已执行。 2026-11-14本次仅核配置与导航，未重做MAT/CSV、SHA数值链或FIG审计；result_date及历史审计结论不变。 R1本轮10类已执行检查0失败，另外4项参考资格检查未执行；源/密度/Green核/FullFK双项/MAT-CSV重建通过。最终12条误差线一致，两图重导出像素差0。BIE参考未完成资格，不以小残差或退出0替代收敛。 2026-11-16本轮小矩阵/右端/密度/复场/边界对照差0，C0对历史差0，独立密度重建7.69e-15、MAT/CSV最大6.91e-14 Pa；65文件SHA相同。C1–C4与参考经验分辨资格未执行，原STOP仍保留。 本次只对脚本/输入/恢复状态机作本地检查，声学调用与SSH连接0。最终dry-run Code Analyzer 0消息、七个声学入口哨兵有效；C1–C4真实数值验收尚未执行。 最新仅静态核回收110项SHA、包36项SHA、C0–C2提交链及C3 solved绑定；保存审计回执均PASS，未复跑MAT或Green核。 v2本轮现有MCP单线程重建C3：476平台点+19额外GL112点通过，密度接收重建1.40e-15、边界9.68e-8；原110回收SHA不变。31项MATLAB/25项完整性/29项runtime通过，模拟C4最终化测试不作窗口资格。 最新R1固定Seed保存数据重审计42.01s通过，BIE密度重建1.89e-14、PS直接Green2.16e-8、FullFK双项1.10e-14，所有原门保持；原v1 Nonmonotonic失败不追改，v2获批后诊断保留。 本次R2只作本地静态/输入/无声学shell检查，R1 287 SHA未变；没有MATLAB/声学/远程运行，不升级多Seed资格。 本次现有MCP保存数值重审0失败：BIE密度重建1.81e-14、PS Green8.03e-9、FullFK1.16e-14、独立边界6.81e-8；原两F02有144170像素差，本地12FIG核对及原renderer保存重放通过。 本次原287/229/239源SHA前后不变；两Seed完整五层/综合/411点/5倍审计全通过，原FIG24+新FIG24+汇总12曲线一致，10次重开PNG通道差0。本地新封存不追改原远端失败。 本轮23个已读来源版本SHA保持；本地Gaussian恢复最大1.23e-14、谱式差0，TEST_ONLY重建最大2.75e-14、RMS采样差最大4.83e-7。原R1/R2密度和封存不修改、不重算。 本次574项原回收前后SHA不变，45组本地独立保存审计0失败，重建2.75e-14、RMS采样6.75e-7、坡度采样3.35e-5；原117及新54条FIG数据一致。初次MCP不可用快照保留，用户恢复会话后已补齐本次必需审计。 D0仅原源/P0/v2/选定K0来源SHA及静态审查，正式深海场0；既有MCP合成57计数断言、FIG60线核验通过，实际源调用与SSH0。 D0本次150回收SHA/28冻结输入及全阶段链、165远端门通过；127本地保存检查与16阶压力差分（最坏1.70e-8<1e-6）、原60/新70FIG曲线通过；原SHA保持，0新源/传播/反射求解及SSH。首次本地struct初始化失败保留。 D1本轮仅新静态/工具检查：26离线、10传输mock；未读取或改写旧密度/旧MAT，不继承声学PASS。 最新只读诊断核原会话和CSV差、快照SHA；未执行MAT重建或补审计，原C0失败未抹去。 D1本轮原59SHA前后有效，原C0 MAT和common字节保持；新增独立审计通过但不改原失败现场或旧审计状态。 D1当前完整独立保存审查通过；原失败和旧数值门未删改。最终远端数值与本地图件分阶段，必须本地图件通过才合并D1a资格。 2026-11-19 K1：最终4040SHA，原N7独立数值资格复用；五面现有MCP保存五层/综合/分辨全部0失败/0新求解。新六例72曲线＋汇总56曲线/112次零像素通过及实际预览；最后原生退出−9/控制器78未改为PASS，独立验证原seal并审查保存结果后接受科学资格。

**当前重跑条件：**D1a已完成，无需重跑；读取封存MAT/CSV/FIG及最终审查即可。任何新PM支持/频率/Seed/生产集成仍需独立授权，其他工作线原重跑边界保持。 K1已限定完成，无需再Start；原远端Status继续保留退出FAIL，当前科学结论读本地最终封存。不要因退出FAIL自动重复已封存声学计算。

证据入口：

- [当前二维任务表（导航，不授予运行许可）](../reports/active_task_tracker.md) · `reports/active_task_tracker.md`
- [AGENTS执行与源宽规则](../AGENTS.md) · `AGENTS.md`
- [10 cm配置与导航同步记录](../reports/agents_configuration_alignment_2026-11-14.md) · `reports/agents_configuration_alignment_2026-11-14.md`
- [历史20 cm代表参数及三频计划（后续顺序以当前任务表为准）](../reports/piston_source_2d_single_frequency_plan_2026-11-11.md) · `reports/piston_source_2d_single_frequency_plan_2026-11-11.md`
- [正式收尾报告](../reports/two_dimensional_channel_closeout_report.md) · `reports/two_dimensional_channel_closeout_report.md`
- [S4 结果索引](../results/comparisons/two_dimensional_channel_closeout/20261026_s4_final_closeout_01/result_index.md) · `results/comparisons/two_dimensional_channel_closeout/20261026_s4_final_closeout_01/result_index.md`
- [保存证据审计摘要](../results/comparisons/two_dimensional_channel_closeout/20261026_s4_final_closeout_01/s4_saved_evidence_audit.json) · `results/comparisons/two_dimensional_channel_closeout/20261026_s4_final_closeout_01/s4_saved_evidence_audit.json`
- [原二维研究工作台](../UWA_2D_research_review_workbench.html) · `UWA_2D_research_review_workbench.html`
- [二维活塞自由场图示与历史尺寸溯源](../reports/piston_source_quicklook_2026-11-11.md) · `reports/piston_source_quicklook_2026-11-11.md` · 历史版本 `a37fd75 + source snapshots`
- [活塞quicklook最终检查点（非P0资格）](../results/validation/piston_source_quicklook/20261111_2d_strip_01/checkpoint_final.json) · `results/validation/piston_source_quicklook/20261111_2d_strip_01/checkpoint_final.json`
- [活塞P0/P1执行与恢复记录](../reports/piston_p0p1_remote_execution_2026-11-12.md) · `reports/piston_p0p1_remote_execution_2026-11-12.md` · 历史版本 `compare/6c53632 + immutable snapshots`
- [此前许可等待检查点（历史）](../results/validation/piston_source/20261112_remote_p0p1_01/checkpoint_license_wait.json) · `results/validation/piston_source/20261112_remote_p0p1_01/checkpoint_license_wait.json`
- [活塞P0/P1最终报告](../results/comparisons/piston_2d_single_frequency/20261112_remote_p0p1_01/report.md) · `results/comparisons/piston_2d_single_frequency/20261112_remote_p0p1_01/report.md` · 历史版本 `20261112_remote_p0p1_01`
- [活塞最终检查点](../results/validation/piston_source/20261112_remote_p0p1_01/checkpoint_final.json) · `results/validation/piston_source/20261112_remote_p0p1_01/checkpoint_final.json` · 历史版本 `20261112_remote_p0p1_01`
- [10cm形状补图：SSH等待与理论解释](../results/validation/piston_source_shape_comparison/20261114_shape_01/report.md) · `results/validation/piston_source_shape_comparison/20261114_shape_01/report.md`
- [R1固定10 cm PM结果与成本停止报告](../results/comparisons/piston_2d_pm/20261115_r1_seed260001_01/report.md) · `results/comparisons/piston_2d_pm/20261115_r1_seed260001_01/report.md`
- [R1阶段检查点（参考未合格）](../results/comparisons/piston_2d_pm/20261115_r1_seed260001_01/checkpoint.json) · `results/comparisons/piston_2d_pm/20261115_r1_seed260001_01/checkpoint.json`
- [默认关闭BIE缓存实现、C0实测阶段成本与STOP](../reports/r1_bie_cache_implementation_2026-11-16.md) · `reports/r1_bie_cache_implementation_2026-11-16.md`
- [缓存续行最终检查点（C1–C4未派发）](../results/comparisons/piston_2d_pm/20261116_r1_bie_cache_01/checkpoint.json) · `results/comparisons/piston_2d_pm/20261116_r1_bie_cache_01/checkpoint.json`
- [R1-BIE手动启动、状态、恢复、回收与新执行政策](../reports/r1_bie_manual_execution_2026-11-16.md) · `reports/r1_bie_manual_execution_2026-11-16.md`
- [手动包无声学dry-run回执（不是C1–C4数值PASS）](../results/validation/r1_bie_manual_preparation/20261116_ready_03_dry_run/receipt.json) · `results/validation/r1_bie_manual_preparation/20261116_ready_03_dry_run/receipt.json`
- [runtime配置根因修复、新目录与测试](../reports/r1_bie_runtime_json_fix_2026-11-16.md) · `reports/r1_bie_runtime_json_fix_2026-11-16.md`
- [R1-BIE非单调停止只读诊断](../reports/r1_bie_nonmonotonic_diagnosis_2026-11-16.md) · `reports/r1_bie_nonmonotonic_diagnosis_2026-11-16.md`
- [R1-BIE诊断政策、SHA迁移与仅C4续行方案（待批准）](../reports/r1_bie_nonmonotonic_policy_c4_plan_2026-11-16.md) · `reports/r1_bie_nonmonotonic_policy_c4_plan_2026-11-16.md`
- [R1-BIE v2仅C4正式包、验证与手动操作命令（未启动）](../reports/r1_bie_v2_c4_implementation_2026-11-16.md) · `reports/r1_bie_v2_c4_implementation_2026-11-16.md`
- [R1-BIE v2 Status修复与用户会话状态（待数值回收审查）](../reports/r1_bie_v2_status_fix_2026-11-16.md) · `reports/r1_bie_v2_status_fix_2026-11-16.md`
- [R1固定Seed最终回收审查：参考资格、误差位置与R2候选](../reports/r1_bie_v2_final_review_2026-11-16.md) · `reports/r1_bie_v2_final_review_2026-11-16.md`
- [R2有限跨Seed技术方案](../reports/r2_pm_multiseed_plan_2026-11-17.md) · `reports/r2_pm_multiseed_plan_2026-11-17.md`
- [R2本地检查与未运行项](../reports/r2_pm_preparation_validation_2026-11-17.md) · `reports/r2_pm_preparation_validation_2026-11-17.md`
- [R2 PowerShell手动命令](../reports/r2_pm_manual_commands_2026-11-17.md) · `reports/r2_pm_manual_commands_2026-11-17.md`
- [R2逐Seed研究入口（默认dry-run）](../scripts/validation/run_r2_pm_seed.m) · `scripts/validation/run_r2_pm_seed.m`
- [R2 seed260002回收诊断](../reports/r2_seed260002_retrieval_diagnosis_2026-11-17.md) · `reports/r2_seed260002_retrieval_diagnosis_2026-11-17.md`
- [R2保存数据诊断入口](../scripts/validation/audit_r2_seed260002_saved.m) · `scripts/validation/audit_r2_seed260002_saved.m`
- [R2最终保存审查/封存报告](../reports/r2_saved_closeout_final_2026-11-17.md) · `reports/r2_saved_closeout_final_2026-11-17.md`
- [R2限定三Seed正式指标](../results/comparisons/piston_2d_pm/r2_saved_closeout_20261117_01/three_seed_summary/qualified_metrics.csv) · `results/comparisons/piston_2d_pm/r2_saved_closeout_20261117_01/three_seed_summary/qualified_metrics.csv`
- [R2保存收尾入口](../scripts/validation/closeout_r2_saved_seed.m) · `scripts/validation/closeout_r2_saved_seed.m`
- [pm_wind_cutoff_k0_k1_plan_2026-11-17.md](../reports/pm_wind_cutoff_k0_k1_plan_2026-11-17.md) · `reports/pm_wind_cutoff_k0_k1_plan_2026-11-17.md` · 历史版本 `2026-11-17 K0 preparation only`
- [pm_wind_cutoff_preparation_validation_2026-11-17.md](../reports/pm_wind_cutoff_preparation_validation_2026-11-17.md) · `reports/pm_wind_cutoff_preparation_validation_2026-11-17.md` · 历史版本 `2026-11-17 K0 preparation only`
- [pm_wind_cutoff_manual_commands_2026-11-17.md](../reports/pm_wind_cutoff_manual_commands_2026-11-17.md) · `reports/pm_wind_cutoff_manual_commands_2026-11-17.md` · 历史版本 `2026-11-17 K0 preparation only`
- [manifest.json](../results/experiments/pm_wind_cutoff/preparation_20261117_05/manifest.json) · `results/experiments/pm_wind_cutoff/preparation_20261117_05/manifest.json` · 历史版本 `2026-11-17 K0 preparation only`
- [K0正式审查与K1矩阵](../reports/pm_wind_cutoff_k0_review_2026-11-18.md) · `reports/pm_wind_cutoff_k0_review_2026-11-18.md` · 历史版本 `2026-11-18 K0 saved review; no acoustic solves`
- [K0本地45组独立保存审计](../results/validation/pm_wind_cutoff_retrieval/numeric_20261118_01/receipt.json) · `results/validation/pm_wind_cutoff_retrieval/numeric_20261118_01/receipt.json` · 历史版本 `2026-11-18 K0 saved review; no acoustic solves`
- [K0新图与派生表回执](../results/visualization/pm_wind_cutoff_review/20261118_01/receipt.json) · `results/visualization/pm_wind_cutoff_review/20261118_01/receipt.json` · 历史版本 `2026-11-18 K0 saved review; no acoustic solves`
- [K1选定七条件及连续几何界](../results/visualization/pm_wind_cutoff_review/20261118_01/k1_selected_geometry.csv) · `results/visualization/pm_wind_cutoff_review/20261118_01/k1_selected_geometry.csv` · 历史版本 `2026-11-18 K0 saved review; no acoustic solves`
- [D0五深度入射复场技术方案](../reports/d0_incident_depth_plan_2026-11-18.md) · `reports/d0_incident_depth_plan_2026-11-18.md`
- [D0准备验证与未运行门](../reports/d0_incident_depth_preparation_validation_2026-11-18.md) · `reports/d0_incident_depth_preparation_validation_2026-11-18.md`
- [D0用户手动命令](../reports/d0_incident_depth_manual_commands_2026-11-18.md) · `reports/d0_incident_depth_manual_commands_2026-11-18.md`
- [k1_resources_strategy_2026-11-18.md](../reports/k1_resources_strategy_2026-11-18.md) · `reports/k1_resources_strategy_2026-11-18.md` · 历史版本 `K1 preparation only`
- [k1_preparation_validation_2026-11-18.md](../reports/k1_preparation_validation_2026-11-18.md) · `reports/k1_preparation_validation_2026-11-18.md` · 历史版本 `K1 preparation only`
- [k1_manual_commands_2026-11-18.md](../reports/k1_manual_commands_2026-11-18.md) · `reports/k1_manual_commands_2026-11-18.md` · 历史版本 `K1 preparation only`
- [K1原生非声学11项诊断](../results/validation/k1_preparation/native_20261118_07/retrieved/receipt.json) · `results/validation/k1_preparation/native_20261118_07/retrieved/receipt.json` · 历史版本 `native diagnostics only; no acoustics`
- [D0五深度正式研究报告](../reports/d0_incident_depth_review_2026-11-18.md) · `reports/d0_incident_depth_review_2026-11-18.md`
- [D0本地保存数组审查127项](../results/validation/d0_retrieval/saved_review_20261118_02/receipt.json) · `results/validation/d0_retrieval/saved_review_20261118_02/receipt.json`
- [D1a/D1b方案与成本](../reports/d1_depth_comparison_plan_2026-11-18.md) · `reports/d1_depth_comparison_plan_2026-11-18.md` · 历史版本 `2026-11-18 D1 preparation only`
- [D1准备验证与未执行门](../reports/d1_depth_preparation_validation_2026-11-18.md) · `reports/d1_depth_preparation_validation_2026-11-18.md` · 历史版本 `2026-11-18 D1 preparation only`
- [D1手动命令与最终03包](../reports/d1_depth_manual_commands_2026-11-18.md) · `reports/d1_depth_manual_commands_2026-11-18.md` · 历史版本 `2026-11-18 D1 preparation only`
- [K1首跑中断诊断](../reports/k1_interruption_diagnosis_2026-11-18.md) · `reports/k1_interruption_diagnosis_2026-11-18.md` · 历史版本 `read-only failure diagnosis; no retry`
- [D1首跑CSV精度失败诊断](../reports/d1_first_run_diagnosis_2026-11-18.md) · `reports/d1_first_run_diagnosis_2026-11-18.md` · 历史版本 `2026-11-18 read-only diagnosis`
- [K1回收保存审查与N7/N12分析](../reports/k1_retrieval_review_2026-11-18.md) · `reports/k1_retrieval_review_2026-11-18.md` · 历史版本 `N7 numeric PASS; K1 STOPPED; local layout FAIL`
- [D1修复、保存C0审计及恢复准备](../reports/d1_recovery_preparation_2026-11-18.md) · `reports/d1_recovery_preparation_2026-11-18.md` · 历史版本 `2026-11-18 recovery preparation`
- [D1恢复包直接手动命令](../reports/d1_recovery_manual_commands_2026-11-18.md) · `reports/d1_recovery_manual_commands_2026-11-18.md` · 历史版本 `2026-11-18 recovery preparation`
- [D1保存C0独立审计回执](../results/validation/d1_recovery/20261118_01/saved_c0_audit/receipt.json) · `results/validation/d1_recovery/20261118_01/saved_c0_audit/receipt.json` · 历史版本 `2026-11-18 recovery preparation`
- [D1最终恢复包MATLAB准备回执](../results/validation/d1_recovery/20261118_01/final_matlab/receipt.json) · `results/validation/d1_recovery/20261118_01/final_matlab/receipt.json` · 历史版本 `2026-11-18 recovery preparation`
- [K1诊断与保留证据恢复方案（仅规划）](../reports/k1_diagnosis_recovery_plan_2026-11-18.md) · `reports/k1_diagnosis_recovery_plan_2026-11-18.md` · 历史版本 `PLAN_ONLY; no new package or execution`
- [K1恢复04准备与诊断结果](../reports/k1_recovery_preparation_2026-11-18.md) · `reports/k1_recovery_preparation_2026-11-18.md` · 历史版本 `preparation only; official recovery unrun`
- [K1恢复04手动命令](../reports/k1_recovery_manual_commands_2026-11-18.md) · `reports/k1_recovery_manual_commands_2026-11-18.md` · 历史版本 `preparation only; official recovery unrun`
- [K1恢复04首跑图件失败诊断](../reports/k1_recovery04_graphics_failure_2026-11-18.md) · `reports/k1_recovery04_graphics_failure_2026-11-18.md` · 历史版本 `read-only observed failure; no retry`
- [D1双深度最终研究报告](../reports/d1_depth_comparison_final_2026-11-18.md) · `reports/d1_depth_comparison_final_2026-11-18.md` · 历史版本 `2026-11-18 D1 final`
- [D1最终保存独立审查](../results/validation/d1_goal_closeout/20261118_01/receipt.json) · `results/validation/d1_goal_closeout/20261118_01/receipt.json` · 历史版本 `2026-11-18 D1 final`
- [D1最终双深度图件审计](../results/visualization/d1_depth_summary/20261118_04/receipt.json) · `results/visualization/d1_depth_summary/20261118_04/receipt.json` · 历史版本 `2026-11-18 D1 final`
- [D1目标模式执行与失败记录](../reports/d1_goal_execution_2026-11-18.md) · `reports/d1_goal_execution_2026-11-18.md` · 历史版本 `2026-11-18 D1 final`
- [K1目标模式执行记录](../reports/k1_target_execution_2026-11-18.md) · `reports/k1_target_execution_2026-11-18.md` · 历史版本 `2026-11-19 final saved science qualification; old exit failures retained`
- [K1七条件最终研究报告](../reports/pm_wind_cutoff_k1_final_review_2026-11-19.md) · `reports/pm_wind_cutoff_k1_final_review_2026-11-19.md` · 历史版本 `2026-11-19 K1 limited saved closeout`
- [K1最终科学封存](../results/validation/k1_target_execution/20261118_01/final_checkpoint.json) · `results/validation/k1_target_execution/20261118_01/final_checkpoint.json` · 历史版本 `2026-11-19 K1 limited saved closeout`
- [K1汇总14行/5754点审计](../results/comparisons/pm_wind_cutoff_k1/final_summary_20261119_01/receipt.json) · `results/comparisons/pm_wind_cutoff_k1/final_summary_20261119_01/receipt.json` · 历史版本 `2026-11-19 K1 limited saved closeout`
- [K1本地测试与退出限制](../reports/k1_target_closeout_validation_2026-11-19.md) · `reports/k1_target_closeout_validation_2026-11-19.md` · 历史版本 `2026-11-19 K1 limited saved closeout`

### 二维宽带与多模型参考

状态：M1_INCOMPLETE · P8 RAM门停止；结果日期：2026-11-11（M1部分对照/资源停止）；正式S5/S6仍暂停；本次核查：2026-11-11。

M0封存保持。M1已完成P4/P6和10组AS，限定数值对照通过；P8因RAM门中止，步长与吸收层资格未完成，无生产推荐配置。

**为什么做：**要把声场变成可用于通信的信道冲激响应（CIR，声音在不同延迟上的复响应），频率之间必须使用一致的声源、相位和时间原点，还要知道反射参考是否可靠。几个单频点相近并不足以保证时域结果可信。

**做了什么：**10 月 27 日已有 129 点同源直达检查和弱正弦 P22 的 6/8 kHz 单点参考。10 月 31 日的独立试点完成 4062.5 Hz、4031.25 Hz 的五层 BIE 参考，并检查了事先锁定的一个插值预测。 本轮完成原生reduced/physical/DSP三表示修复，P22平面波eta=0/±.01m与Gaussian同源角谱、统计和CIR独立控制通过。 本轮对37项相关文件绑定当前SHA，读取旧MAT/CSV并核主海面裁剪；17项限定门通过。设计5次生产PE及10组AS以分离空间/步长误差。 2026-11-11静态复核M0产物181/181 SHA一致，所有清单内MATLAB代码一致；确认吸收层单方向三次增长、x/y相加。连接前置门失败已保存，未启动PE或AS。 同日MCP恢复后执行生产P4/P6与10组AS，参考自身边缘、三频空间对照与4/6kHz PE-AS通过；独立收尾81行比较与MAT/CSV/SHA审计0失败。默认关闭的观察接口逐步记录三段传播，64²开关输出差0，已有M0宽带H消费者回放BER/SER=0。 后续只读诊断核41项源码、run03内部209项产物SHA一致。系统可用RAM下降1.724GiB，但P8采样工作集1.149–1.281GiB且生涯峰值不增长，未证实泄漏。

**说明什么：**说明频率表示与局部参考已有实际进展，后续可以据此设计更有针对性的验证。两个新频点只是小范围试点，还不能替代整个频带。 允许单独安排M1；不是全频物理CIR或三维粗糙海面认证。 M1_INCOMPLETE是执行受阻，不表示物理门数值失败；原吸收层既未判合格也未判不合格。 本轮证明指定均匀/平面Gaussian离散锚点的已完成对照，不形成完整M1配置资格。8kHz AS通过不能代替未完成生产PE或步长实测。

**还差什么：**完整 65/129 点反射参考、九点覆盖、高频与时间窗检查仍缺。4031.25 Hz 的 PE 独立步长/网格敏感性未齐；FullFK 与 BIE 的差异低于约定分辨率，不能据此判优。原正式 S5/S6 仍暂停。 最终封存发现两频点检查点所指审计代码 SHA=2c95abb6… 的原文未定位；固定提交的通用审计脚本为 bd865b9c…，不能顶替完整旧最终审计的来源。6 个产物绑定及 287 条数据清单匹配，本轮新 MATLAB 限定数值重核有效。 旧结果缺当次源码SHA，边缘仅33/65个z采样；新平面4–8kHz网格/步长与全步边缘资格待执行。 本轮M1没有新MAT数据、三频误差、逐步边缘或吸收对照；当前会话版本、真实PID与可靠长作业追踪尚未确认。 当前P8仅直达1036/1036、入射977/1068、反射0/33平面；Z4/Z8及SP4/SP6/SP8均未派发。原吸收层既未判合格也未判不合格。 缺运行期系统占用归因及内部数组分配轨迹；优化尚未实施、没有实测节省或P8重跑。部分MAT没有中途psi_k，不能直接从976步续算。

**下一步：**先审阅P8内存方案：优先缩短无用数组寿命、缓存固定环境屏，再做保持原运算链和逐步诊断的等价/资源验证。当前不派发传播，M1仍INCOMPLETE；后续恢复需当前资源、预算、版本/PID和输入证据门。

**技术范围：**二维 Tx=(0,100 m)，P22=(0,3 m)，原弱 C2 正弦 A=0.01 m、K=0.1 rad/m。10-31 状态为 DATA_COMPLETE_REFERENCE_PASS_MODEL_RESOLUTION_LIMITED；全带重建、物理 CIR/PDP/BER 未生成。10-27 的 48 h 提案尚未获批，不能作为后续执行许可。

**已有结果：**10-27 和 10-31 的结果包、检查点与审计摘要存在。 11-10 合并镜像 234 条与增补 53 条 SHA 清单均核验一致。

**历史审计：**历史两个新频点参考通过但模型分辨受限；11-10 独立 MATLAB 保存标量/源数组/配置/相位/退出记录与预测公式重核通过，未重演原全部场合成与边界数组检查。

**当前重跑条件：**两频点原求解入口已从 24a05ba 恢复；BIE helper 少几项历史输入检查，已列差异但未替换/执行。远端布局、预测冻结及调度环境未恢复为可直接运行环境；正式流程仍暂停。 本轮MCP attachment失败，禁止将远程或新MATLAB进程作为替代；当前M1任务已授权，不重复申请每个子步骤。 当前MCP版本和PID已核实，原连接失败成为历史；本轮停止原因更新为P8 RAM门。原run02空diary与300s传输超时保留，经确认终态后以显式stdout捕获在新run03重跑P4且三表示圆盘差0。

证据入口：

- [正式流程暂停与门槛](../reports/two_dimensional_s5_s6_2026-10-27.md) · `reports/two_dimensional_s5_s6_2026-10-27.md`
- [10-31 两频点检查点](../results/validation/two_dimensional_s5_p22_dualfreq/20261031_01/checkpoint.json) · `results/validation/two_dimensional_s5_p22_dualfreq/20261031_01/checkpoint.json` · 历史版本 `24a05ba72b58042084417e54ca437f03adb988dd`
- [两频点独立审计摘要](../results/validation/two_dimensional_s5_p22_dualfreq/20261031_01/audit_dualfreq_complete_01/audit.json) · `results/validation/two_dimensional_s5_p22_dualfreq/20261031_01/audit_dualfreq_complete_01/audit.json`
- [两频点运行入口（已恢复）](../scripts/validation/run_two_dimensional_s5_p22_dualfreq_remote.m) · `scripts/validation/run_two_dimensional_s5_p22_dualfreq_remote.m` · 历史版本 `24a05ba72b58042084417e54ca437f03adb988dd`
- [P22 可复现基线与对应表](../reports/p22_reproducibility_baseline_2026-11-10.md) · `reports/p22_reproducibility_baseline_2026-11-10.md` · 历史版本 `747af7ee88abda0bd7f554ddbd940ad3bc5f8cb9`
- [约100m静态模型收尾计划](../reports/static_channel_model_100m_closeout_plan_2026-11-10.md) · `reports/static_channel_model_100m_closeout_plan_2026-11-10.md`
- [M0合同修复与验收报告](../reports/m0_phase_contract_closeout_2026-11-10.md) · `reports/m0_phase_contract_closeout_2026-11-10.md` · 历史版本 `3f7071d + working tree`
- [M0独立数值与MAT/CSV审计](../results/validation/m0_phase_contract/20261110_m0_native_dsp_v2_01/independent_audit_final/checkpoint.json) · `results/validation/m0_phase_contract/20261110_m0_native_dsp_v2_01/independent_audit_final/checkpoint.json`
- [M0封存与M1窗口/网格/步长预审设计](../reports/m0_seal_m1_preflight_2026-11-10.md) · `reports/m0_seal_m1_preflight_2026-11-10.md` · 历史版本 `M0 source 969127f`
- [M1历史原始数据只读审计](../results/validation/m1_preflight/20261110_saved_evidence_04/audit.json) · `results/validation/m1_preflight/20261110_saved_evidence_04/audit.json`
- [M1执行收尾：连接门停止](../reports/m1_static_convergence_closeout_2026-11-11.md) · `reports/m1_static_convergence_closeout_2026-11-11.md`
- [M1前置检查点](../results/validation/m1_convergence/20261111_m1_execution_preflight_01/checkpoint_preflight.json) · `results/validation/m1_convergence/20261111_m1_execution_preflight_01/checkpoint_preflight.json`
- [M1数值执行收尾与RAM停止](../reports/m1_static_convergence_qualification_2026-11-11.md) · `reports/m1_static_convergence_qualification_2026-11-11.md`
- [已完成对照独立审计](../results/validation/m1_convergence/20261111_m1_execution_03/partial_closeout_audit/checkpoint.json) · `results/validation/m1_convergence/20261111_m1_execution_03/partial_closeout_audit/checkpoint.json`
- [18项矩阵实际状态](../results/validation/m1_convergence/20261111_m1_execution_03/matrix_status.csv) · `results/validation/m1_convergence/20261111_m1_execution_03/matrix_status.csv`
- [P8内存诊断与无物理改动优化方案](../reports/m1_p8_memory_diagnosis_2026-11-11.md) · `reports/m1_p8_memory_diagnosis_2026-11-11.md`
- [P8只读内存诊断回执](../results/validation/m1_p8_memory_diagnosis/20261111_01/checkpoint.json) · `results/validation/m1_p8_memory_diagnosis/20261111_01/checkpoint.json`

### P22 宽带与通信控制实验

状态：局部通过 · 已知信道控制；结果日期：2026-10-30（P22历史）；2026-11-10（M0小型控制）；本次核查：2026-11-10。

历史P22已知CSI控制保留原资格；本轮新增修复后小型无噪声链路门，不能拿旧P22代替新物理真值。

**为什么做：**先把一条受控的“频响 → 时域信道 → 发射 → 加噪 → 均衡 → 误码统计”链路走通，排除表示或通信实现本身的问题。

**做了什么：**最初直接比较 65/129 点 IFFT 没通过。后来把分量按共同时间窗投影到 257 点网格，并与独立计算的 257 点 PE 对照，通过原有门槛，随后完成平面与弱正弦两种海面的 QPSK 实验。 11-10 最小入口已恢复；原检查点输入/输出 SHA 核验，保存频响共同窗重建及原 QPSK 固定配置重放均通过原审计。 M0发现线性卷积截尾与FFT逆的边界不匹配，保留完整尾部后，在相同新生产H上BER从.07421875变为0；解析taps和CIR门通过。

**说明什么：**共同窗表示与已知 CSI 通信链现在有当前可用入口、匹配原记录的数据和本轮独立复核；可以继续用作受控基线。

**还差什么：**接收端事先知道信道，尚未验证从导频估计信道的能力；也没有完成 PE/BIE/FullFK 全频带的通信影响比较。该实验不升级正式 S5/S6，更不能代替目标深海场景验证。 旧carrier-only DSP、CIR、统计和相位敏感通信结果须明确来源后另存迁移/重算，不能直接标为新合同通过。

**下一步：**基线整理已结束；后续如需导频估计或模型失配研究，另定范围、预算和验收门，不自动开始。

**技术范围：**P22=(0,3 m)，发射深度 100 m，平面/原弱 C2 正弦；共同 257 点、64 ms 表示。65/129 共同窗 CIR 差约 7.1e-6，小于 1e-3。QPSK、已知 CSI 块 FFT LMMSE 均衡，SNR=0/6/12/18/24 dB，每点 262144 bit。原始 IFFT 失败记录保留，不能写成原始方法后来通过。

**已有结果：**原共同窗和 QPSK 结果保持原位；新审计、图件、版本/哈希/日志保存在 p22_reproducibility_baseline/20261110_115037。

**历史审计：**原始 IFFT 失败保留。共同窗/QPSK 原审计通过；11-10 保存数据重审通过，均衡符号差 0，CSV 指标最大差 3.9080e-14。

**当前重跑条件：**三个关键入口和两个独立审计已最小恢复；审计支持显式新输出目录，原数值门不变。当前可复核保存数据；新性能实验未授权。

证据入口：

- [最初 CIR 失败检查点](../results/validation/pe_screen_p22_comm/20261030_01/checkpoint.json) · `results/validation/pe_screen_p22_comm/20261030_01/checkpoint.json`
- [共同窗通过检查点](../results/validation/pe_screen_p22_common_window/20261030_01/checkpoint.json) · `results/validation/pe_screen_p22_common_window/20261030_01/checkpoint.json`
- [共同窗独立审计](../results/validation/pe_screen_p22_common_window/20261030_01/audit.json) · `results/validation/pe_screen_p22_common_window/20261030_01/audit.json`
- [最终 QPSK 检查点](../results/communication/pe_screen_p22_qpsk/20261030_03/checkpoint.json) · `results/communication/pe_screen_p22_qpsk/20261030_03/checkpoint.json` · 历史版本 `24a05ba72b58042084417e54ca437f03adb988dd`
- [QPSK 独立审计](../results/communication/pe_screen_p22_qpsk/20261030_03/audit.json) · `results/communication/pe_screen_p22_qpsk/20261030_03/audit.json`
- [QPSK 运行入口（已恢复）](../scripts/validation/run_pe_screen_p22_qpsk.m) · `scripts/validation/run_pe_screen_p22_qpsk.m` · 历史版本 `24a05ba72b58042084417e54ca437f03adb988dd`
- [共同窗验证入口（已恢复）](../scripts/validation/validate_pe_screen_p22_common_window.m) · `scripts/validation/validate_pe_screen_p22_common_window.m` · 历史版本 `24a05ba72b58042084417e54ca437f03adb988dd`
- [P22 可复现基线与对应表](../reports/p22_reproducibility_baseline_2026-11-10.md) · `reports/p22_reproducibility_baseline_2026-11-10.md` · 历史版本 `747af7ee88abda0bd7f554ddbd940ad3bc5f8cb9`
- [M0合同修复与验收报告](../reports/m0_phase_contract_closeout_2026-11-10.md) · `reports/m0_phase_contract_closeout_2026-11-10.md` · 历史版本 `3f7071d + working tree`

### 三维固定海面与深海模型

状态：局部通过 · 目标海况待验证；结果日期：2026-10-27；本次核查：2026-11-10。

固定测试海面在多个独立频点上有条件通过，但还没有证明整个频带或目标随机海况都可靠。

**为什么做：**二维模型不能自动代表真实三维散射，需要对三维点源、二维海面和指定链路另建独立参考。

**做了什么：**固定 compact-dual 海面已有 4/5/6/7/8 kHz 锚点及后续独立频点检查。65 点预测网格又用事先锁定的四个新频点检查，四点都通过约定筛查门。

**说明什么：**固定案例的跨频预测有有限的独立支持，可以继续研究怎样少算一些频点。65 点是预测输出网格，并非 65 个独立参考都已算过。

**还差什么：**连续频带、时域窗和 CIR/BER 尚未认证。目标 PM 随机海面的独立参考仍未解决；Q1 BIE 路线因独立误差门失败已经停止，不应继续沿用早期“下一步修 BIE”的建议。

**下一步：**暂缓1500m深海分支。约100m公共三维模型按新计划独立核约定和数值误差；二维BIE证据不自动认证三维PM，目标PM独立参考缺口保留。

**技术范围：**固定 Tx=(0,150.5,1500 m)、Rx=(0,0,5 m)，C2 compact-dual 测试海面。四个新点为 4750/7250/4125/6562.5 Hz。G1_target_PM=OPEN，WP2_formal=false。此处 FE 指 Dirichlet 场展开，不是有限元。

**已有结果：**正式报告、四点原生结果、预测锁和检查点存在。

**历史审计：**四点筛查和保存数据审计通过，有条件参考的范围不变。

**当前重跑条件：**历史结果来自旧 Ubuntu R2023b；不能直接视为新 Ubuntu R2025b 已通过项目验证。

证据入口：

- [四个新频点正式报告](../reports/deepsea3d_65grid_four_holdout_2026-10-27.md) · `reports/deepsea3d_65grid_four_holdout_2026-10-27.md`
- [结果索引](../results/validation/deepsea3d_65grid_four_20261027/README.md) · `results/validation/deepsea3d_65grid_four_20261027/README.md`
- [Q1 BIE 停止决策](../reports/deepsea3d_bie_q1_decision_2026-10-26.md) · `reports/deepsea3d_bie_q1_decision_2026-10-26.md`
- [三维总研究计划](../reports/deep_sea_channel_research_plan.md) · `reports/deep_sea_channel_research_plan.md`

### 1500 m 目标深海几何与 PE 窗口

状态：0–200m目标未验证 · D1a有限范围完成；结果日期：2026-11-01；本次核查：2026-11-18（D1a有限2D反射范围，D1b仅成本/设计）。

已经算清楚几条候选路径什么时候到；还没有验证生产 PE 在这个长距离、大偏移场景下能否算准。

**为什么做：**先明确目标场景里直达、海面和海底一次反射的距离与时延，再判断计算窗口够不够大，避免接收点已经落进吸收边缘却把结果当成物理变化。

**做了什么：**水深 1500 m、Tx1490 m、Rx3/6/10 m、偏移 0–200 m 的解析几何表和图已生成。海底一次反射的候选到达也补算过。另对生产 PE 的窗口、网格和海绵层限制进行了代码核查。 D1a已实际完成1490m源、Rx6m、±36m的固定2D反射五层参考及模型误差比较；该结果不恢复整个0–200m目标。

**说明什么：**后续传播验证有了独立的几何参照。Rx 深度 10 m 时，平面海面与平坦海底的一次反射几何到达重合；这只说明到达时间重合，声压强弱还不知道。

**还差什么：**400 m 窗口下，200 m 接收位置落入吸收层；保持既有分辨率的更大窗口又受到当前窗口/网格保护限制。实际内存、运行时间与长程复声场收敛未测量。 D1b新周期/粗糙支持与稠密资源仍未认证，未模拟海底/多次反射/结构。

**下一步：**按用户2026-11-10要求暂不推进1500 m孔径及大偏移计算；原解析几何和限制保留，后续有需要再另立任务。

**技术范围：**等声速 c=1500 m/s 的镜像法，603 行几何采样；没有计算路径复声压、海底反射系数、粗糙散射或通信。600/800 m 只是候选窗口，不是已通过的配置。

**已有结果：**D/S 和 D/S/B 解析报告、CSV 和图件入口存在。

**历史审计：**解析代数、CSV 和图件检查通过；没有目标 PE 传播通过结论。

**当前重跑条件：**解析几何入口在当前工作树；生产 PE 目标数值试验尚未准入。

证据入口：

- [解析几何与 PE 孔径设计](../reports/rc1_target_deepsea_geometry_and_pe_aperture_design_2026-11-01.md) · `reports/rc1_target_deepsea_geometry_and_pe_aperture_design_2026-11-01.md`
- [海底一次反射几何补充](../reports/rc1_flat_bottom_single_reflection_geometry_2026-11-01.md) · `reports/rc1_flat_bottom_single_reflection_geometry_2026-11-01.md`
- [D/S/B 几何入口](../scripts/reporting/render_target_deepsea_dsb_geometry.m) · `scripts/reporting/render_target_deepsea_dsb_geometry.m`
- [D1a限定结果与D1b可行性](../reports/d1_depth_comparison_final_2026-11-18.md) · `reports/d1_depth_comparison_final_2026-11-18.md` · 历史版本 `2026-11-18 D1 final`

### Bellhop 代表案例与平面对照

状态：图示已完成 · 对照有未过项；结果日期：100 m 代表案例报告未标日期；连通 smoke 为 2026-11-01；本次核查：2026-11-10。

当前代表图采用 100 m 水深、发射机离底 1 m；声线展示完成，但不能把画出来的路径当成幅相已经验证。

**为什么做：**借助独立的射线方法理解多途结构，并在可比较的平面条件下检查 PE 与 Bellhop 的差异。

**做了什么：**Tx99 m、Rx3 m 的 4 kHz 代表案例已有接收贡献图和完整发射扇区图。另有平面镜像展开对照的保存表：波束收敛、宽带幅相和时延等检查通过，空间相位及复场差仍有未通过项。

**说明什么：**已有直观的路径材料，也已经定位出对照尚未闭合的部分。连通 smoke 说明程序能运行；代表图和部分通过项各有自己的用途。

**还差什么：**未完成角度/宽带扫描，也未认证粗糙面全部幅相。平面对照的空间相位约 0.0612 rad，大于 0.05 rad 门；空间复 L2 约 0.0613，大于 0.02 门。

**下一步：**若要继续比较，先限定相同路径拓扑、声源和相位约定，解释已有失败项；不要拿含底高阶多途直接对照当前 PE 的直达加一次海面反射。

**技术范围：**当前安装：水深100 m，Tx99 m，Rx3 m；停止零离底 epsilon 扫描。代表案例193条 ARR、199条 E 贡献记录、1001条展示声线；这些数量不代表真实独立到达路径数。旧500 m材料已经退役，本页不引用归档内容。

**已有结果：**代表案例报告、图件及平面对照 checks.csv 存在。

**历史审计：**代表图检查完成；平面对照 12 项中有 2 项未通过，不能记为整体通过。

**当前重跑条件：**当前绘图入口存在，仍依赖外部 Bellhop 可执行文件；本次没有重跑。

证据入口：

- [当前安装与代表图报告](../reports/bellhop_100m_source1m_visuals_report.md) · `reports/bellhop_100m_source1m_visuals_report.md`
- [当前重建计划顶部修订](../reports/bellhop_100m_rebaseline_plan.md) · `reports/bellhop_100m_rebaseline_plan.md`
- [平面展开逐项检查表](../results/validation/pe_bellhop_unfolded_flat_gaussian/checks.csv) · `results/validation/pe_bellhop_unfolded_flat_gaussian/checks.csv`
- [当前代表图入口](../scripts/reporting/generate_bellhop_100m_source1m_visuals_vertical.m) · `scripts/reporting/generate_bellhop_100m_source1m_visuals_vertical.m`

### 运行环境与重跑条件

状态：新机最小项目基线通过 · 完整迁移未做；结果日期：P22 原结果至 2026-10-31；本机保存数据复核及新机最小基线 2026-11-10；本次核查：2026-11-10。

新机隔离快照已通过小型公共信道的跨平台数值比较和两项 P22 保存结果审计；范围以外的项目组件仍待检查。

**为什么做：**后续计算应能说明用哪份代码、哪台机器和哪套输入得到结果，避免把环境打通当成科研验证通过。

**做了什么：**11-10 先恢复 P22 五个缺失文件并独立复核保存数据；随后在新机隔离快照上传 23 个白名单文件，逐项 SHA 通过。三组 CPU double 小型信道案例与两项 P22 保存结果审计在 Linux/R2025b Update 1 通过，与本机 Update 6 数值对照通过；92 个回传文件 SHA 匹配。 M0新增独立run，保留SHA回执失败、通信边界失败、CSV舍入和hex导入失败，再在新子目录恢复/复核。条件模型schema3与phase_contract_version=2.0.0阻止旧表示静默混入；代码/结果SHA和完整日志可追溯。 M0源码已有969127f，轻量检查点/失败记录/哈希白名单封存；既有47份生成源码快照保留历史，新快照不再自动纳入。旧窗口源码可恢复至8946b8e，但不能伪补原run绑定。

**说明什么：**当前来源和输入身份可追溯；新机已取得限定于 128×128、4 kHz 直达/反射、4–8 kHz 九点接口及两项 P22 保存审计的最小项目环境与兼容性基线。

**还差什么：**/home/misc/cwx 仍非完整 Git 工作区；仅隔离快照有本轮 23 个文件。未查完整项目全部 toolbox、Linux MEX、外部程序、大网格或长时任务；九点不认证 CIR，正式 S5/S6 继续暂停。两频点旧最终审计代码 SHA=2c95abb6… 的原文仍未定位，不能用固定提交中 SHA=bd865b9c… 的通用脚本顶替；旧失败证据保留。

**下一步：**本轮停止。后续如获授权，再核完整项目目标工作区、代码与输入身份，并按拟执行任务设计兼容性、资源和数值验收；不自动迁移或启动科研重计算。

**技术范围：**本机现有 MCP R2025b Update 6，新增工作连同进度看板刷新保守计 110/180 s；新机 SSH BatchMode、R2025b Update 1、单计算线程经实测，七次远端 MATLAB 进程 204.02/600 s，最长 50.70/180 s。五个数值作业逐项过 RAM/负载/磁盘门；未访问旧主机、归档或他人目录。

**已有结果：**本轮隔离结果包保存逐文件来源清单、补丁差异、完整本机/远端输出与日志、配置及原生数值数组、逐作业资源回执、对照 JSON、上传及回传 SHA 和检查点；原 P22 结果包未覆盖。

**历史审计：**此前基础 FFT smoke 仅说明基础环境；本轮最小真实项目案例和 P22 原保存数据审计建立限定基线，不追认旧 R2023b 结果或两频点完整旧最终审计。

**当前重跑条件：**只可按本轮隔离快照、固定配置和资源门重建限定基线；更大项目范围仍需独立授权与资格检查。

证据入口：

- [新机交接与边界](../reports/new_ubuntu_ssh_matlab_handoff_2026-11-09.md) · `reports/new_ubuntu_ssh_matlab_handoff_2026-11-09.md`
- [新机基础 smoke 检查点](../results/validation/new_ubuntu_matlab_smoke/20261109_232540_new_ubuntu/checkpoint.json) · `results/validation/new_ubuntu_matlab_smoke/20261109_232540_new_ubuntu/checkpoint.json`
- [现行公共信道入口](../vertical_channel_model.m) · `vertical_channel_model.m`
- [P22 可复现基线与对应表](../reports/p22_reproducibility_baseline_2026-11-10.md) · `reports/p22_reproducibility_baseline_2026-11-10.md` · 历史版本 `747af7ee88abda0bd7f554ddbd940ad3bc5f8cb9`
- [新机最小项目基线报告](../reports/new_ubuntu_minimal_project_baseline_2026-11-10.md) · `reports/new_ubuntu_minimal_project_baseline_2026-11-10.md`
- [新机最小项目基线最终检查点](../results/validation/new_ubuntu_project_baseline/20261110_121507_new_ubuntu_minimal/project/results/validation/new_ubuntu_project_baseline/20261110_121507_new_ubuntu_minimal/checkpoint_final.json) · `results/validation/new_ubuntu_project_baseline/20261110_121507_new_ubuntu_minimal/project/results/validation/new_ubuntu_project_baseline/20261110_121507_new_ubuntu_minimal/checkpoint_final.json`
- [M0合同修复与验收报告](../reports/m0_phase_contract_closeout_2026-11-10.md) · `reports/m0_phase_contract_closeout_2026-11-10.md` · 历史版本 `3f7071d + working tree`
- [M0封存与M1窗口/网格/步长预审设计](../reports/m0_seal_m1_preflight_2026-11-10.md) · `reports/m0_seal_m1_preflight_2026-11-10.md` · 历史版本 `M0 source 969127f`

## 最近发生了什么

### 2026-10-26 · 二维核心正式收尾

把已有比较收成一个范围明确的结论；三个预定 PM 实现完成，新增五个因预算没有运行。之后可以复用成果，不必从二维起点重做。

### 2026-10-26 · 停止当前三维 Q1 BIE 参考路线

线性求解残差虽小，独立边界和场误差仍未过门。因此停止继续扩大这条路线，旧计划中的继续调参建议不再有效。

### 2026-10-27 · 正式宽带流程停下来等待审阅

直达与少数频点参考有进展，但全频带反射、九点和通信依赖还不齐。48 小时只是一份待批提案，没有开始全带计算。

### 2026-10-27 · 三维预测接受四个新频点检查

固定海面的预测在四个新点通过有限筛查。它支持继续研究稀疏预测，但没有把 65 点网格变成 65 个独立物理参考。

### 2026-10-30 · 先发现 CIR 比较失败，再修正比较表示

原始 65/129 点 IFFT 差异没有过门；共同窗、按分量和时延处理后，与独立 257 点 PE 对照通过原门槛。原失败记录保留。

### 2026-10-30 · P22 已知信道 QPSK 控制实验跑通

通信链获得受控仿真结果和独立审计。这里尚未从导频估计信道，也没有升级正式多模型宽带资格。

### 2026-10-31 · 补齐两个邻近频点的参考试点

4062.5/4031.25 Hz 的 BIE 参考完成，局部预测也接受检查；模型差异的分辨能力和全带覆盖仍是缺口。

### 2026-11-01 · 目标深海路径几何算清楚，窗口问题浮现

有了直达、海面和海底一次反射的几何到达表；代码检查同时发现 200 m 偏移与现有 PE 窗口/吸收层冲突，尚未进入传播试验。

### 2026-11-09 · 新 Ubuntu 通路打通，项目仍待迁移

连接、MATLAB 版本和独立 FFT 检查通过；项目目录还不是 Git 工作区。机器能运行基础计算与项目能复现是两件需要分别验收的事。

### 2026-11-10 · 建立统一进度首页，补上索引漏记的进展

把二维、三维、P22 控制实验、目标几何与环境放在同一条阅读路径中；登记缺失入口和持续更新规则。本次没有新增声场或通信计算。

### 2026-11-10 · P22 最小入口已恢复，保存数据已独立复核

恢复 5 个缺失工作流文件，隔离审计输出；23 条旧检查点 SHA 和 287 条两频点清单匹配，共同窗/QPSK 原门和两频点限定重核通过。首次回执脚本错误保留并计时，原证据不覆盖。仍缺全带多模型、PE 敏感性和新机项目资格；下一步只能审阅本基线，另行授权后再设计后续工作。

### 2026-11-10 · 两频点最终旧审计代码来源仍缺一环

封存时六项产物 SHA 匹配，但检查点未给路径的最终审计代码 SHA 与固定提交的通用审计脚本不同。已保留首次封存失败记录，不替换旧检查点或假装代码一致。本轮保存标量/配置/相位重核可用，完整旧最终审计仍不能宣称已复现；本轮停止，不从归档补齐。

### 2026-11-10 · 新机最小项目兼容性基线通过

在不覆盖 /home/misc/cwx 原文件的隔离快照中上传 23 项白名单并验证 SHA；本机与新 Ubuntu 三组 128×128 CPU double 公共信道数值逐元素通过原容差，两项 P22 保存结果原审计通过，92 个回传文件 SHA 匹配。只建立本轮小案例与保存审计资格；九点不认证 CIR，完整项目依赖与迁移仍缺，正式 S5/S6 暂停及两频点旧审计来源缺口保留。下一步须另行授权并按目标任务核代码、输入、资源与数值门，本轮停止。

### 2026-11-10 · 提交基线成果，优先收尾约100m静态模型

P22恢复入口与新机最小基线已提交为5d2e9f8。按用户要求暂缓1500m，接下来先核生产入口的声源/高度/相位，再做窗口、宽带CIR和粗糙海面范围验证。本次只是计划，不新增数值结论。

### 2026-11-10 · M0相位合同与最小控制闭合，停在M1之前

修复原生PE到DSP的全包络共轭及C/P、增广实统计、迁移、SSA和桥接，物理/统计/CIR/无噪声控制和独立审计通过。保留回执与CSV打包失败；通信截尾边界失败在相同H上修复为BER=0。旧结果未覆盖，未升级正式S5/S6；下一步仅可另行授权M1。用户要求停止main合并，已abort并保留compare原HEAD。

### 2026-11-10 · 封存M0，改用历史192m证据设计M1

确认969127f已包含M0，核SHA并追加轻量失败/验收凭据。只读旧窗口数据，确认192m合格范围、160m边缘失败、33/65深度采样及源码绑定缺口；冻结5个PE/10组AS的最小新矩阵和保守资源门。本轮无新传播，下一步须新任务授权M1。

### 2026-11-11 · M1完整执行在现有MCP连接门停止

已获5个无吸收层PE、10组AS和阶段门后3个吸收层PE的14400s任务授权。核M0产物181项及源码身份、相位v2和吸收层定义后，MCP返回failed to attach to MATLAB session，无法确认版本、真实PID和日志/状态追踪，故18作业均未派发。无新数值资格和生产推荐，M1_INCOMPLETE；保留错误与SHA检查点。下一步仅在现有MCP恢复且资源/版本/PID门通过后从P4继续，不启动M2。

### 2026-11-11 · 活塞比较先定代表参数与严格二维范围

用户暂无实物参数，委托查阅后自行选择。依据厂家半波长口径量级，取20 cm工程参照；首轮严格二维使用20 cm宽匀速条带、向上发射、同一振速，在4/6/8 kHz沿单接收线比较PE/BIE/FullFK。圆形活塞与无限条带分开解释。已保存共源、采样、窗口、参考与新Ubuntu分阶段并行计划；没有实现新源或运行仿真。下一步先共源控制和4 kHz试点，取得远端成本后再定正式范围；原Gaussian收尾和其他窗口M1任务状态与预算保留。本次现有MCP正式渲染成功，未核或恢复M1执行门。

### 2026-11-11 · M1恢复运行后在P8触发RAM预留门

现有MCP确认R2025b Update 6/PID28836，完成P4/P6与10组AS；三频AS空间和4/6kHz生产同格对照及MAT/CSV/SHA独立审计通过。首次P4空diary/300s超时保留并计时，新run显式捕获stdout重跑且数值差0；首版参考审计错误保留后修正，仅重审已有MAT。P8在入射976步因RAM713232384B低于769405747B停止，Z4/Z8与全部SP未派发。故M1_INCOMPLETE，无生产推荐、不启动M2；下一步仅在继续指令和资源复核后恢复未完成M1。活塞计划及其他工作流记录/预算保留。

### 2026-11-11 · P8确认RAM余量门触发，给出无物理改动优化方案

只读核源码/回执/SHA和当前workspace元数据，没有传播或生产修改。系统可用RAM下降1.724GiB，P8工作集采样变化不足以解释全部下降，未证实泄漏、具体系统来源仍缺证据。列出缩短数组生命周期、固定环境屏缓存和分块逐步统计的方案及验证门；节省量只是静态估算。M1仍INCOMPLETE，不自动重跑或进入M2。

### 2026-11-11 · 活塞三频自由场图示完成，先P0再平面

核GitHub/本地历史确认5 cm为圆形半径；严格二维改用总宽20/10/5 cm速度孔径积分，不复用圆形方向图。保存三频绝对复Pa、接收线与归一化方向图；独立角谱/求积/MAT-CSV/FIG重建检查无失败，仅支持本次可视化。20 cm三频全宽约112/67/49°，8 kHz宽角域出现深谷，大离轴可弱于小孔径。缺生产共源压力/梯度、相位、Gaussian回归与平面反射；本轮停止，下次先P0。未运行PE/BIE/FullFK/SSH，不改Gaussian或转用M1预算。

### 2026-11-12 · 活塞图件修正版检查并封存

继续完成前一日quicklook交付：保留试绘标题重叠、首版负增益截断和一次exportgraphics多容器异常；改用单布局专用标题行，8张交付PNG逐图查看、FIG重建36幅CData差0，24条接收曲线匹配，0失败。声压数据未重算，数值源码SHA未变。交付绘图加图审计50.688s；报告、MAT/CSV、源码快照、日志和SHA封存。仍缺P0共源压力/梯度/相位与Gaussian回归，之后才做P1平面；本轮停止，不自动推进。

### 2026-11-12 · P0部分数值完成，PE门失败后加密进程等待许可

在新Ubuntu独立快照验证56项SHA并实跑最多两个单线程MATLAB。源/梯度/速度/相位、BIE/FullFK输入和Gaussian标量默认回归完成；PE512m压力对齐但隐含梯度未过原1e-4门。保留失败并派发8192m加密，发现新进程停在MathWorks账号提示、尚无加密数值结果；已请用户恢复远端非交互许可，不代填账号、不重发或改用本机。123项部分证据回传SHA一致。P0完整资格与P1均未完成，先核原作业状态再续，旧M1/S5S6不升级。

### 2026-11-12 · 远端P0/P1完成并回传独立审计

2026-11-12最终封存。活塞P0与4 kHz P1自由场、平面、弱C2正弦均通过，独立MAT/CSV/密度重建审计失败0项，80条FIG曲线核验一致，五组最终PNG已目视检查。弱面20 cm反射复场PE Model-0/BIE相对L2为0.0048932994，FullFK/BIE为2.0947621e-6，均超过五倍经验数值分辨率门。周期截窗初场失败保留，接受结果采用验证专用连续角谱平方根PE与精确平面背景加紧支撑相屏扰动。20 cm可作为该范围主源。后续生产入口与更广工况需要单独检查，本轮止于P1。 许可等待、512/8192m梯度失败、周期传播失败、native启动崩溃和首次审计pi变量错误均保留，恢复时使用新job且先检查原进程。

### 2026-11-14 · 活塞形状补图待SSH恢复

已复用既有报告和源码定义并实现MATLAB入口；理论求根表明10cm两种孔径在三频前向均无零点。SSH公钥失败且MCP无法附着，PNG/FIG未运行；下一步仅在恢复SSH后执行本轮自由场补图。进度生成视图尚未刷新。

### 2026-11-14 · 10cm孔径形状与远视野图完成

复用480800个旧场值与97m圆弧，补算外侧自由场及200/500m圆弧；五张PNG/FIG已检查，18文件回传SHA匹配，孔径求积最大差6.94e-13。10cm条带/圆形三频前向均无零点，8kHz圆形主瓣更宽。仅自由场图示，不新增传播资格，停止本轮；MCP不可用，生成进度视图未刷新。

### 2026-11-14 · 10 cm后续主线与AGENTS配置导航统一（无求解）

仅修复配置与当前导航：后续主源总宽10 cm，20 cm及Gaussian结果和所有历史失败保持原条件；顺序统一为4 kHz固定PM、反射模型决策、相干4–8 kHz/CIR及受控通信。R1仍未授权，额外种子和各后续阶段保留独立范围/资源/参考门，Bellhop独立并行。仅静态检查并按维护说明纯文本打包看板，无MATLAB、SSH、实验或Git提交，未重做历史数值审计。记录见reports/agents_configuration_alignment_2026-11-14.md；下一步仅可讨论或明确授权有界R1范围，不自动执行。

### 2026-11-15 · R1三模型411点场已算，BIE后续参考层在成本门停止

10 cm、4 kHz、固定PM seed260001的PS/完整FullFK数值检查通过，BIE一层线性及离节点残差通过；暂定复L2为15.9235%和0.0747208%。首层265 s后，余下参考按实测预计2130 s超过773 s余额，未派发；全部进程439/1100 s。已执行独立保存数据检查0失败，BIE网格/自求积/窗口收敛与参考分辨资格仍缺。最终四图、MAT/CSV/日志/SHA和停止检查点齐备，首次图件审计失败保留；下一步等待R1审查，不自动R2。

### 2026-11-16 · 默认关闭BIE缓存与C0等价性通过；剩余参考成本STOP

按用户明确指令实施缓存并完成4组N64/N128×qa64/96原/默认/关/开及C0原3200节点回放，矩阵、右端、密度、411复场和边界差0。C0调用55.48 s，边界残差26.93 s、基础组装15.43 s；余层加收尾按原门预测730.29 s超过决策时582.80 s余额，未派发C1–C4。新增86 s、累计525/1100 s，单线程串行、RAM预留门通过，65文件SHA一致、独立重建通过。缓存不改变物理模型，R1参考未合格；下一步等待成本与参考续行审查，不自动R2。

### 2026-11-16 · 手动C1–C4执行包及断点恢复准备完成（未求解）

按用户指令准备批处理、后台启动/状态/回收工具、阶段密度/411复场/计时/日志/资源/SHA与独立审计检查点；C0仅导入，C1–C4原顺序和1e-4门保留。取消旧1100秒累计及730秒预测停止，仍单进程单线程、无parpool、30%RAM和异常停止。包35项SHA、PowerShell/Bash与7声学入口哨兵dry-run通过，恢复/篡改拒绝测试通过；本轮SSH连接、声学计算均0。等待用户手动启动，真实C1–C4及Linux完整链未实跑，不升级参考资格或R2。

### 2026-11-16 · 手动包runtime.json启动配置修复（未求解）

2026-11-16 — runtime.json故障修复：移除有printf重定向歧义的AWK最大值，验证整数后改Bash比较；保留starting/initial/baseline独立字段，Python标准库序列化并立即解析校验，异常在MATLAB前退出78。29项本地测试与新包无声学dry-run通过；新包/新目录和命令见 reports/r1_bie_runtime_json_fix_2026-11-16.md。未连接SSH/未运行C1–C4，旧失败session/结果不改。

### 2026-11-16 · C3严格非单调门STOP的静态诊断

只读回收CSV/JSON/源码/SHA，C0–C2已提交，C3求解/审计通过但门失败，C4未运行；三次新求解计数合理。空间最坏差在同一Q4增8.18%，绝对均约1e-7且低于1e-4，原严格单调硬门触发，并非单点异常或已确认求解器错误。不修改阈值/求解器/原结果，不启动MATLAB/SSH/Resume；建议先审查门与显式迁移方案。

### 2026-11-16 · 非单调诊断政策与仅C4续行方案完成（未实施）

只做静态审查和方案设计：建议新包将严格空间单调作为完整诊断，原四项1e-4及残差/全部中点/独立审计/5倍分辨规则保持。设计旧新双链SHA迁移、C0–C3不重求解、C3新工作区重建审计和唯一C4求解。110/36项SHA、原封存链及8个代码版本核对一致；没有修改旧包/求解器/门或执行MATLAB/SSH/Resume，等待用户批准实施，R1仍STOPPED。

### 2026-11-16 · v2零求解迁移与仅C4手动包准备通过

按用户批准实现独立v2 policy/manifest/旧新双链，严格单调改为保留旧FAIL的诊断，原四项1e-4及全部数值/边界/独立审计/5倍分辨规则不变。C0–C3零求解迁移，C3本地单线程476+19点独立审计PASS；31项MATLAB、25项完整性/准入、29项runtime及归档重定位dry-run通过。正式04包封存，原110项SHA未变，新密度求解0、SSH0，等待用户手动启动C4；实际最终参考仍未合格，不进入R2，不提交推送Git。

### 2026-11-16 · Status传输修复，用户原C4会话退出0

首次无检查点提示发生在启动空档，SSH只读见C4运行、调用1次。客户端改进引入的bash -c引号错误已用.NET UTF-8标准输入修复，7项传输与12种状态回归及实际只读调用通过。用户同一session现退出0、检查点报告完成；本轮不启动/Resume/终止/改包、不Fetch完整结果，不据状态宣布数值PASS，等待回收审查。

### 2026-11-16 · R1单Seed经验参考及保存数据审查通过

两份回收各287项SHA和证据一致，原/迁移/C4/最终封存有效。现有MCP单线程42.01s重审C4 476+19点、综合源/密度/PS Green/完整FullFK/MAT-CSV/分辨率0失败，12条FIG及两组PNG对一致，四图已查看。PS复L2=15.9235%且相位项99.7691%，FullFK=0.0747208%仍超阈值约478倍；原非单调FAIL保留，原门不变。R1只在固定10cm/4kHz/seed260001/411点范围完成，本轮新增密度求解0/SSH0，不启动R2、不提交推送Git。

### 2026-11-17 · R2有限跨Seed执行准备完成，等待手动启动

仅10 cm、4 kHz、原PM谱和411点，冻结260002/260003既有系数及各自独立C0–C4流程，原v2与数值门/5倍规则不变。用户允许两Seed各单线程并行，共享30%RAM/24GiB准入。03包各44项manifest/49项交付SHA、各28输入测试、29 runtime、12状态、7 stdin、8无网络mock、10 RAM边界及两归档重定位dry-run通过，R1 287 SHA不变。MCP不可连接，MATLAB验收/Analyzer及Linux原生flock测试未运行，列为声学前强制验收；本轮SSH/MATLAB/正式求解均0。保留失败与候选，停止等待用户手动启动，不提交推送Git，不推进R3/宽带/通信。

### 2026-11-17 · R2 seed260002数值复核通过，原PNG像素FAIL保留

两份回收229/239 SHA和输入/阶段链有效；260002原退出1停在F02 PNG像素零差门，五层BIE均提交。MCP保存数值复核36.28s/0失败，12FIG线一致，原renderer新目录保存重放15.97s像素差0。PS15.7052%、FullFK0.0490634%；双进程重叠558s，RAM最低约90.95%基准，非资源停止。260003原退出0/最终封存有效但本轮未完整数值重审。不SSH/Resume/新求解/改门/改原件/提交推送，当前停止于分析。

### 2026-11-17 · R2两Seed保存审查/新封存及正式三Seed配对完成

两Seed各五层BIE保存密度及综合PS/FullFK/411点/5倍分辨0失败；固定轴/刻度/原画布重开，5图10次PNG差0，60条FIG曲线一致。各125/84文件新封存，原287/229/239回收SHA前后不变，260002原FAILED_STOPPED和候选失败保留。三SeedPS15.9235/15.7052/18.0283%、FullFK0.0747208/0.0490634/0.0556842%，均分辨；正式表与原快照资格分开。新求解/SSH/远端MATLAB/Git均0，限定R2完成后停止，不自动R3/宽带通信。

### 2026-11-17 · PM风速—谱截止联合设计与K0手动包准备完成

为区分风速、物理截止和二者交互，建立3×5×3配对嵌套设计；风速定义统一U_PM，原Gaussian可恢复且U6原值保留，新扩展流显式记录。05包15项SHA/23来源/重定位，287项MCP测试与3个TEST_ONLY保存审计、7传输mock、13离线supervisor通过。只证明准备路径，正式45组K0和原生远端链未执行；K1七条件候选须先审查K0并另定资源，入口硬禁用。原R1/R2与生产不改、SSH/新声学求解/远端MATLAB/Git0；停止等待用户手动K0，不持续监控。

### 2026-11-18 · K0正式45组审查完成，K1矩阵选定但不执行

为决定风速与物理截止的最小声场对照，核回收574项SHA/45组检查点与远端退出0/0；MCP恢复后45组独立保存审计0失败，117条原FIG及54条新FIG曲线一致。K0显示N32→40高度比仅1.001–1.020而二阶导数比1.200–1.711，不能把采样通过当截止收敛。连续几何界排除辅助平面/接收相交，但不认证声学；保留seed260001七条件。K1新面独立BIE/PS/FullFK/411点与资源包仍待明确，无新海面生成/声学求解/SSH/Git，下一步仅可准备有界手动方案。

### 2026-11-18 · D0五深度入射场手动包准备，正式声场未运行

为区分Tx100m历史基准与1490m深海目标，复用10cm/4kHz原源、梯度、独立角谱和phase v2，冻结五深度/±220m场及K0原±50m粗糙辅助。06包28项SHA，57项MCP合成准备断言、60条FIG线、15基础设施/10传输mock通过；只证明准备软件路径。真实深海角谱/梯度/采样与Ubuntu运行仍未执行，D1粗糙面范围不足；用户手动启动和Fetch后再分析。不连接SSH、不计算正式场或反射、不推进K1/宽带/Git。

### 2026-11-18 · K1资源实查与最终05手动包准备完成

实际主机28核/56线程、可用约45GiB但Swap耗尽/PSI高，PCT未安装且额外并发席位未知，因此选1×1；结合15阶段历史BIE估65–85min、12GiB/例、30%预留与3h累计。严格锚点/296项SHA及24项MAT、32离线、8传输mock、TEST_ONLY汇总/FIG通过，只证明准备路径；正式声学未运行、生产和数值门不改。下一步用户资源达标后手动Start，当前不排队、不监控、不提交推送。

### 2026-11-18 · K1原生补测与资源补查完成

补充2026-11-18 12:21原生诊断：11项锁/身份/真实信号通过，四个轻量子进程全部回收；资源快照可用约39GiB、PSI下降，当时只读准入通过。仍选1×1，正式Start和声学未运行。 35项回收SHA及05包归档SHA通过；不占用新MATLAB许可证、不调整数值门、不提交推送。下一步仅用户手动Start，Codex停止监控。

### 2026-11-18 · D0五深度回收审查完成：幅度平坦仍不等于平面波

为核实100m基准到1490m目标的入射场变化，核150项回收SHA与165远端数值门，本地127保存检查/原60条FIG及新8图70线通过。1490m x200m幅降.6514%但相位223.8975rad；轴上1/5/10/40m切平面复误差.069/1.603/6.346/89.643%，相位曲率主导。全线传播谱模量不变，窗谱集中不能当源波束变窄。D0只限均匀二维4kHz入射场；D1大范围粗糙支持与深海反射误差仍缺，下一步须独立定范围。本轮0新声源/反射求解、SSH/Git0。

### 2026-11-18 · D1a深度配对包准备，D1b成本评估；正式声学未运行

为隔离发射深度对反射模型误差的影响，冻结0.10m/4kHz/原U6/N12/260001/L160/C2，设Tx100/1490m、Rx6m、共同候选±36m/411点；每深度独立C0–C4，修正旧PS/Green中写死的3m。03包31SHA、26离线+10传输mock通过，MATLAB/原生Linux/声学/FIG均未运行；本轮SSH/MATLAB/Git提交0。D1b需要新周期/支持，禁止平铺，稠密候选22400节点单矩阵7.48GiB但整工作内存无准入资格。下一步仅用户审阅授权后手动D1a，不监控、不推进D1b或宽带。

### 2026-11-18 · K1首跑中断：N7封存、N40原生致命退出

K1用户手动05包首跑STOPPED：U6/N7五层及封存SHA链复核通过；U6/N40在MODELS约18.833s报double free or corruption后-9退出，无BIE层，后四例未派发。资源记录未触发本包停止门；根因尚待原生栈定位。不重启、不删调用标记，先Fetch保留现场。 本轮仅既有文件/资源只读核验，未新启动MATLAB或求解，未降低门或提交推送。

### 2026-11-18 · D1用户首跑C0审计失败，定位CSV坐标精度；未重启

用户03包运行133.784s、退出1/1，最低可用43.70GiB/RSS峰1.66GiB，无资源停止原因。100m PS/FullFK自身门通过并新算C0，但writetable将循环小数坐标写短，386点差≤4.974e-14m触发严格isequal；1490m未运行。本轮两次只读SSH，保留原日志/输入与C0 CSV及SHA，无MATLAB启动/修复/重算/Git。下一步另行授权D1专用往返精度修复和保存C0审计恢复，不自动重复求解或推进D1b。

### 2026-11-18 · K1回收审查：N7数值通过，N40原生根因未全定位

K1回收575项SHA及296冻结输入通过；N7本地五层保存密度/独立边界/PS/FullFK/411点与原5倍分辨0数值失败，PS/FullFK复L2=10.251%/0.0208328%，N12锚点15.9235%/0.0747208%。N40确切PID转储为GTP_1线程abort/glibc栈，具体原生组件未定位；未进入BIE。本地候选图12条数据一致但导出缩画布触发原布局门，保留FAIL。整体仍STOPPED，不重启或删除未封存调用。 未新求解/重启/降门/改生产/Git；只允许后续独立诊断和明确恢复设计，不自动推进矩阵。

### 2026-11-18 · D1精度修复及原C0补审通过，恢复包待手动启动

为继续原D1而不重复已算C0，回收原失败包并核59SHA，17位CSV保留严格isequal和科学门；现有MCP保存C0独立审计8.535s/0失败，340中点及19点重建通过、0新BIE。恢复03包125SHA、7文件Analyzer0/23离线/12传输mock和真实MATLAB准备通过。新任务先重审导入C0/common再剩余9层，旧133.784s及RAM基准保留，不改原失败/生产/R1/R2/K0/D0或Git；尚未上传启动，不推进D1b/K1/宽带。

### 2026-11-18 · K1诊断恢复方案完成，尚未实施

K1诊断恢复方案已制定、尚未实施：先本地保存修图＋两次不同的有界非声学诊断，再由用户另行手动启动一次N40 MODELS尝试；通过后原层级顺序继续余五面，最多25次新BIE。N7数值不重算、原失败保留，新旧资格双链；继承两次旧会话636s/10800s和119021486080bytes RAM基准，原30%/12GiB/1×1/数值门不改。 只读本地账时发现较早压力停止会话，科学资格/执行状态不升级；下一步先授权实施非声学诊断及具体新包，不能把计划当执行许可。

### 2026-11-18 · K1诊断结束，恢复04包准备完成

K1恢复04包准备完成：A1/A2非声学PASS，用户追加Ubuntu保存图12线/布局/零像素门PASS；四次Windows修图仍FAIL保留，候选代码不进恢复包。继承旧636s及新诊断88s=724s，1459SHA/原575及296链保持，37项MAT准备/25离线/8传输mock/重定位通过。N7原资格导入且不重算，N40一次带日志MODELS后其余五面最多25新BIE，正式恢复未上传/启动。原生崩溃具体组件仍未定位。 下一步仅用户手动启动，Codex不持续监控；非声学PASS不保证N40崩溃已修复，不升级矩阵数值资格。

### 2026-11-18 · K1恢复04在声学前图件像素门停止

K1恢复04用户首跑在saved_graphics原PNG零像素门STOPPED：预检退出0，F02第一次重导399701通道值差、8条曲线数据一致；第二次/F04未审计。新worker和BIE均0，N40未开始，N7旧资格保持；累计账790s，原N40单例18.833s不变。仅只读日志/已有PNG诊断，无重画/重启/降门。 后续需跨独立进程验证保存图件稳定性，不能把过去单次PASS或第二个未执行的零值当作通过。

### 2026-11-18 · D1a双深度限定完成；PS显著改善但未接近完整FullFK

为回答发射深度是否改变反射模型误差，在原0.10m/4kHz/单PM/C2和Rx6m/±36m/411点完成两深度各五层独立BIE。PS复L2由16.1917%降到3.54362%，FullFK由0.0765349%降到0.0709202%，变化均超过经验分辨。全D1仅10次BIE，无重复100m；542SHA、本地保存113.497s/0失败、80FIG曲线/零像素并目视通过。CSV精度、原生图例和导入冲突失败及累计1507.537s全部保留，数值与本地图件阶段分别审计。只限当前二维单频有限面；D1b新周期/0–200m/内存成本只有方案，下一步另定授权，不自动推进宽带、通信或K1。

### 2026-11-18 · K1目标模式直接执行：最多四案例并行

K1按最新授权由Codex直接执行，4个基础MATLAB并发许可与100GiB资源通过；图件保存FIG重开基准修复通过双进程原像素门。目标任务最多4案例、各单计算线程，N40 MODELS及首批4个C0已封存，后续仍运行；N7不重算。旧累计/失败/RAM基准继承，科学门不变，不再要求用户手动启动。

### 2026-11-19 · K1七条件限定完成；原生退出异常保留

K1按用户直接目标授权完成六新面独立参考（原N7五层导入＋25新BIE），严格R1锚点单列。固定0.10m/4kHz/U_PM/seed260001/原C2/411点，4040SHA及本地完整保存审查0失败，新128图曲线/56FIG/112次零像素和实际预览通过，14模型行/5754点。PS10.251–53.845%、完整FullFK0.02083–0.19553%，风速/物理截止交互经验可分辨；N40不是完整PM。四并发实测峰9.96GiB、累计2656s/10800s，30%RAM保持。最后一例封存后原生退出−9和控制器78保留，独立保存审查通过后接受结果，不重算、不称退出环境稳定。任务案例进程已退出，停止于K1有限范围。

## 重新开始先做什么

### 1. D1a已完成：读取最终报告，另定D1b或R3范围

深源PS约3.54%、FullFK约0.071%；只限当前有限C2面。D1b需要新周期/配对、低频和稠密成本方案。

完成标志：D1现有任务已完成，无需用户执行脚本；后续工作以新的具体授权为准。

执行边界：不要再Start旧恢复包；不自动扩面/Seed/频率、生产替换、宽带或通信。

### 2. K1最终05包：用户手动上传/启动

1×1顺序六新面及独立审计/汇总，检查时压力未达准入，启动时重新测量。

完成标志：用户手动运行并回收实际退出码、全部数值/独立审计/图/SHA后再判K1资格。

执行边界：准备不是数值PASS；不自动启动/重试/扩Seed、宽带或通信。

### 3. R2限定完成，等待用户明确R3范围

R1/R2已提供三个固定4kHz/10cm/411点海面的正式配对证据；可讨论PS适用性、FullFK物理近似和另行界定FastFK算法误差。

完成标志：用户另行明确模型决策问题、证据范围和任何计算/资源门。

执行边界：本轮已停止，不自动执行R3，不将FullFK误差小写成严格等价，不删除PS或扩大范围。

### 4. 后续频响/CIR/通信仍按独立阶段

三Seed单频空间误差不是相干频响或CIR；只有模型决策及新的频率/相干性资格后，才讨论宽带与通信。

完成标志：同源激励和接收定义下独立通过频率、相位/群延迟、CIR/线性卷积和通信验证。

执行边界：W1/C1等待明确授权；不得因为R2完成自动扩频率、CIR/QPSK/BER或Bellhop。

### 5. W1/C1候选：模型决策后才进入相干频响、CIR和通信

反射模型决策后，再另定10 cm相干4–8 kHz直达/反射/总H(f)及频率网格、相位/群延迟和留出频点控制；通过CIR/线性卷积检查后，才做固定种子的受控QPSK与CSI比较。

完成标志：同一源激励、功率/噪声及接收设置下，频响相干性、CIR和通信链各有独立验证与范围限定。

执行边界：W1和C1仍为WAITING；每阶段单独明确范围、预算及参考门。三个离散单频不等于宽带/CIR资格，正式S5/S6不升级。

### 6. BH1候选：Bellhop独立并行（未授权）

沿既有水深100 m、Tx99 m、Rx3 m及1 m物理安装净空，先资格化X线源映射、平面直达/首次反射，再考虑粗糙PM；不依赖也不阻塞R1–R3。

完成标志：源映射、相位和绝对幅度资格在各自明确范围成立后，才可讨论绝对Pa模型排序。

执行边界：BH1仍为PARALLEL - NOT AUTHORIZED；不运行旧500 m默认入口、不读取cash/，不自动扩大Bellhop范围。

### 7. 窗口M1停在RAM门，保留已完成证据

P4/P6及全部AS已通过限定审计，P8部分轨迹保留；步长与吸收层尚未取得资格。

完成标志：后续取得完整P8、Z4/Z8及阶段门后的吸收层证据，才可判M1最终资格。

执行边界：当前不自动重跑或进入M2；保留30%启动可用RAM余量、14400s累计与3600s单作业限制。该预算不转用活塞计划。

## 术语随手查

**PE / WAPE：**抛物方程 / 宽角抛物方程：本项目逐步推进声波的传播算法。它是一种近似方法，仍需检查数值和适用范围。

**BIE：**边界积分方程：在本项目的一些受控案例中充当独立数值参考；它也有离散误差，并非实测真值。

**FullFK：**完整 Kirchhoff 积分近似：保留完整积分项的另一种声场计算方法。完整指积分形式，不代表物理上精确。

**CIR / PDP：**信道冲激响应 / 功率时延分布：描述不同延迟上收到的复信号或能量。空间接收线图不能直接当作时域信道。

**已知 CSI：**接收端事先知道信道信息。适合先验证通信链和均衡控制，但不能证明导频估计已经实现。

**QPSK：**四相位调制：每个符号承载两比特，用来测试信道和均衡后的通信表现。

**PM 海面 / compact-dual：**PM 是按指定谱构造的随机海面；compact-dual 是这里固定的测试海面。固定测试海面通过不代表随机海况总体通过。

**S4 / S5 / S6：**二维专项的阶段编号：单频核心收尾 / 宽带与时域 / 通信。P22 的独立控制实验与原正式阶段需要分别记录。

**未分辨：**观测到的模型差小于预先约定的数值分辨能力。只能说本次比较分不清，不能说两种模型严格相同。

## 记录与核查范围

最新D1目标：远端1490m五层、双深度独立数值和回收保存审查、本地图件零像素；0–200m/D1b/宽带不包含。 最新D1：原数据Fetch/59SHA、现有MCP保存C0独立审计和恢复包准备，0新BIE/远端MATLAB/Git。 最新：D1两次只读SSH核失败日志/检查点与坐标CSV，不启动MATLAB、不修改远端或重算。 2026-11-18 D1：代码/几何/成本/包与离线工具检查，SSH/MATLAB/声学/Git提交0；未重做R1/R2/K0/D0数值审计。此前记录：本轮D0回收完整性/远端保存门、本地保存数组/CSV/FFT/几何/平面波分解/FIG及新图，0新声源/传播/反射求解、SSH0。 K1只读SSH/非声学环境诊断、R1/R2成本分析、严格锚点和完整手动包准备/无求解测试；正式K1与后台实验0。

D1最终远端650.819s/0/0，累计所有尝试1507.537s；本地保存113.497s/0失败、80条重开FIG/零像素，542SHA保持。原图形失败保留并明确本地图件阶段，不降低门。 D1 C0补审8.535s/0失败/0新BIE；最终包MATLAB准备1.112s/7文件Analyzer0，23离线+12传输mock通过；恢复实际未运行。 D1用户首跑133.784s/退出1/1，RAM最低43.70GiB、owned组RSS峰1.66GiB，stop_reason=null；386个CSV坐标差≤4.974e-14m导致严格相等失败。本轮0新MATLAB/求解，文本模板刷新。 本轮D1：26离线+10传输mock、PowerShell解析/Bash语法通过；MATLAB/原生Linux/声学/FIG均UNRUN。看板按维护规范纯文本模式从同一JSON/模板生成，未调用MCP。历史运行记录：D0远端152.710s、退出0/0、最低可用36.114GiB（95.55%基准），本地保存44.707s/127检查、最终新图35.345s/70线；原150SHA保持。 远端非声学MATLAB诊断32.918s、退出0，未建pool；本地24项准备、32离线/8传输mock和TEST_ONLY汇总通过，无新声学。

本次核查代码基线：`compare / c90d0ca; M0 identity rechecked; M1 execution attachment gate stopped; M1 resumed in attached R2025b Update 6/PID28836, then P8 RAM-reserve stop; piston free-field quicklook based on a37fd75, local source SHA snapshot; piston P0/P1 remote 20261112_remote_p0p1_01 audited; 2026-11-14 AGENTS navigation alignment only; no numerical re-audit; 2026-11-17 R2 uncommitted research preparation, final03 manifests bind code; no production changes; 2026-11-17 K0/K1 preparation05, manifest-bound sources; no production changes; 2026-11-18 K0 retrieved574 SHA and 45 local saved audits, no new generation/solves; 2026-11-18 D0 saved review, 150 SHA/165 remote gates/127 local checks, no new solve; 2026-11-18 D1a preparation03, 31 SHA/26 offline/10 transport mock; no MATLAB/SSH/solve; 2026-11-18 D1 failed session read-only diagnosis; coordinate CSV round-trip failure, no restart; 2026-11-18 D1 recovery03, saved C0 audit PASS (0 new solves), 125 SHA, MATLAB prep/7 Analyzer/23 offline/12 transport mocks; 2026-11-18 D1 goal complete, 542 retrieval SHA/80 FIG curves/zero pixels; 10 BIE total, no repeated Tx100 solve; 2026-11-19 K1 saved seven-condition closeout,30 BIE/4040 SHA/128 new curves/zero pixels; native shutdown failure preserved`。

唯一维护源：`reports/project_progress.json`。

以后影响进度的任务完成、失败或停止时，修改相关工作线，追加事件，再刷新页面和文字版。旧失败和停止原因保留，记录只说明已有证据，不授予计算许可。

只更新与当次任务有关的内容；没有状态变化时不制造进展。结果日期、核查日期和生成时间分别记录。MATLAB MCP 不可用时注明未执行检查，不改用其他 MATLAB 运行模式。

```matlab
setup_vertical_project; render_project_progress_dashboard();
```
