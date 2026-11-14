# Bellhop第13阶接收距离REAL8最小修复验证（2026-11-14）

结论：本轮限定的输入精度修复与最小对照通过，停止等待人工审查。没有恢复完整A组、运行第14–16阶或B/C/D组。此结论不是完整相干场、收敛或独立FFP资格。

## 范围与冻结条件

工作树为 `C:/Users/ASUS/.codex/worktrees/bellhop_rotated_fig12b/Explain`，分支 `experiment/bellhop-rotated-fig12b`。独立产物目录：[20261114_receiver64_minimal_01](../results/validation/bellhop_rotated_fig12b/20261114_receiver64_minimal_01/)。只用现有MATLAB MCP R2025b Update6和本地隔离gfortran16.2.0。没有新MATLAB、SSH、并行、网络或生产代码变更。

保持水深5000m、Tx4950m、10kHz、Munk c0=1500m/s、轴/尺度1300m、epsilon=.00737、平坦压释放海面及示例流体底（cp1600m/s、rho1.8g/cm³、损耗.8dB/波长）。源为全向三维点源，源宽不适用；不涉及主线0.10m二维条带。

唯一新求解为原冻结A计划的up_order13：Q间隔1.25m、步长1m、原角度扇区±0.0060944284213613285度，101束、原生q列[-1,0,1]m，对应物理x=[1,0,-1]m。10981个接收距离、原始顺序及ENV距离文本逐字保留；只将原201列缩成上述3个原列、401束缩成101束。没有排序、删点、坐标偏移、插值或重新生成物理网格。

## 修复前后证据

| 项目 | 修复前 | 本轮隔离版本 |
|---|---|---|
| MATLAB double原始范围 | 64950–69950m，10981点，有限且严格递增 | 原值、顺序、位模式不变 |
| 最小有效点间隔 | 0.00999999999476131m | 相同 |
| ENV接收距离读取 | REAL4 km，再乘REAL4的1000 | REAL(KIND=8) km，乘1000d0 |
| 原生单精度米标签 | 166对重复，0处逆序 | 0对重复，0处逆序 |
| 首个异常对 | 第990/991点，69850.089999999997和69850.100000000006m，旧标签均69850.09375m | 两个有效坐标逐位保留，标签严格递增 |
| ENV距离文本 | 原始%.17g km文本 | 与原失败输入同一距离行 |
| PRT/进程 | ReadRcvrRanges拒绝，尚未追踪射线；旧exit0不能视为成功 | 无FATAL，stdout/stderr均空，exit0，101束追踪结束 |
| 第13阶反射 | 原失败例未计算 | 101束全部13次：7次海面、6次海底 |

第12阶旧输入没有单精度米重复；第13阶上述坐标所在的高距离量级使两阶段舍入产生重复。这一故障不是原MATLAB网格非单调，也不是ENV文本精度不足。旧失败ENV SHA为 `e4746cb9b98c2d5b72df02fb1faac70bf0fa4a8d9408dad49e0f954d23d12527`，原.r25d SHA为 `1845078c5423b885cc0f3f19fe642bda86deac459050ba71fcb04685851bc6f0`；本轮均未改变。

## 双精度数据链

1. 隔离 `SourceReceiverPositions.f90` 的exact-input分支直接读REAL8数组，检查有限性、原始严格递增，再以1000d0换算。该分支不调用Sort/SubTab；旧原生严格单调性检查仍保留，没有放宽验收门。
2. `ReadEnvironmentBell.f90` 在读取接收距离前识别原.r25d，启用此分支。`ReceiverRangeInput64`是REAL8，精确路径禁止从REAL4坐标回退。
3. `RotatedFieldState.f90` 中 `ReceiverRange64`保持REAL8。原.r25d中的米制double作为有效求值坐标直接读入该数组。ENV km→m本身有最大7.275957614183426e-12m的换算舍入，所以不能把换算后的值冒称与原始MATLAB位模式完全相同；原.r25d保留的米制位模式才是实际求值依据。新增ENV64与侧文件的2个double ULP一致性门，原有2个float32 ULP元数据门不变。
4. 本次锁定的C/R、Cartesian geometric-hat路径中，范围搜索、区间判定、`x_rcvr`、`s/n`投影、面外输运因子和距离扩散缩放全部读取REAL8 `ReceiverRange64`；中途没有single回转。`ScalePressure`保留旧接口参数，但exact分支不使用其REAL4距离实参。
5. `.rr64input`记录输入/有效坐标及KIND；`.rr64eval`在真实影响函数投影之前逐调用断言坐标位模式，每对坐标首次写盘；`.rr64scale`对每个距离的实际缩放输入逐位断言。日志覆盖32943个有效坐标对和10981个缩放距离。
6. 原生SHD头仍是float32标签，最大量化误差0.00375000000349246m。这些标签不用于本次实际求值。原生累积/输出复声压仍为complex32，本修复没有宣称全求解器双精度化。其他波束类型、ARR、bellhop3D和无侧文件路径不在本次数值资格范围内。

有效坐标的证明来自求解前、实际执行时、保存字节独立重读三个阶段，而非仅比较标签或声压。保存的物理映射使用冻结的深度/横向索引及原深度数组，正向r=64950+z逐位一致；不通过r-64950相消运算重造原始深度。

## 数值审计

| 检查 | 结果 |
|---|---:|
| [无求解前置门](../results/validation/bellhop_rotated_fig12b/20261114_receiver64_minimal_01/case02/preflight_checks.csv) | 25/25通过；此时新Bellhop启动数0 |
| [单例运行审计](../results/validation/bellhop_rotated_fig12b/20261114_receiver64_minimal_01/audit01/checks.csv) | 31/31通过 |
| [保存数据独立审计](../results/validation/bellhop_rotated_fig12b/20261114_receiver64_minimal_01/saved_audit01/checks.csv) | 60/60通过，未启动求解 |
| 实际求值r/q位模式不符 | 0/32943对 |
| 实际缩放距离位模式不符 | 0/10981点 |
| 原始SHD复声压对既有最小基准 | complex32位模式完全相同，最大复差0、相对L2差0、最大相位差0rad |
| 原始SHD对保存MAT | 最大复差0 |
| 原始mwdiag/r25diag对既有基准 | 两表全部数值差0 |
| 声压幅度范围 | 1.0575421691600425e-8至1.1377109571072535e-8，全部有限 |
| 射线最大步数 | 120743，小于原570000安全界 |
| 末端单位切向量最大残差 | 1.34559030584569e-13，小于原1e-8门 |
| 原25个完成分量 | 全部与原full_01检查点SHA-256一致，未重算 |
| 本轮保护清单 | 35个原计划、后续分量/失败文件前后哈希一致 |

复声压基准只读复用本轮之前已存在的 `20261114_range_metadata_validation_01/labels1/case`。该例已用相同实际.r25d坐标、101束/3列；新旧.ssp/.iwcov/.rotq/.r25d字节相同，ENV只差原生距离标签行。没有为取得对照再次运行旧版本，也没有把原第13阶失败例冒充有可比较声压。该对照证明读取修复没有改变已有同坐标最小场，不等同独立物理参考验证。

## 版本、失败保留与成本

新binary SHA-256：`b30ecec9f6960e024d612320406a2114ee8a9a0995de4c61b2bdab02c46881b3`。旧binary仍为 `93ab5754e45525ab0bff0e5490035c74ffedac83472e9c298b3690eb9928c107`。构建只写build01，官方工具箱与已接受overlay不变；编译/链接工具返回0。编译独立耗时未计时，未保存完整编译stdout，记为未取得；现有build_manifest.json的purpose字段继承旧文案，本报告与阶段checkpoint说明本次实际REAL8/13阶用途，源码副本和SHA清单可复核。

本轮新Bellhop进程仅1个，PID35640，终态exit0；外部等待实测2.5999938s、进程CPU采样2.328125s、PRT自报2.17s，三者口径不同。20ms采样峰值工作集22.10MiB（不保证覆盖未采样的尾段峰值）；启动可用RAM2.80GiB，运行最低观测2.76GiB，保持1.5GiB启动/256MiB运行门。运行加审计入口9.1773136s，独立保存审计6.2230517s，前置门0.5799115s；均未扩展完整场预算。新目录当前约116MiB，包括构建与失败准备副本；最终精确大小与输入/输出SHA见阶段checkpoint。

准备阶段两次非数值失败也保留：overlay初稿的生成换行问题未编译、未执行；case01静态扫描截取了错误的ELSE，导致1个假失败，全部坐标检查已通过。修正生成器/扫描器后形成overlay02/case02，未修改数值门槛；旧准备文件、旧FATAL和所有已完成分量均未覆盖。

## 交付和恢复边界

新增实现/工作流仅四个入口：`support/prepare_bellhop_rotated_receiver64_overlay.m`、`prepare_bellhop_rotated_receiver64_case.m`、`validate_bellhop_rotated_receiver64_minimal.m`、`audit_bellhop_rotated_receiver64_saved.m`，均在scripts/validation内；公开信道API、生产PE、DSP/通信公式和外部底质假设未改变。已更新相关方法/进度/导航；无需变更vertical_comm_guide的生产公式，也无需运行通信回归。本轮无新科学图件需求，没有启动绘图或完整场重建。

本轮提交的是最小输入精度修复的限定PASS。完整A组、第14–16阶、B/C/D、完整相干总场、收敛、FFP共点和最终图册均继续暂停。等待人工审查本报告、实际坐标位模式和SHA证据后再决定后续范围。
代码审查补充：[隔离Fortran源码前后diff](../results/validation/bellhop_rotated_fig12b/20261114_receiver64_minimal_01/receiver64_source_review.patch)。MATLAB Code Analyzer对case准备与独立审计无问题；生成器有4项string写法提示，运行入口有动态列表提示及1项未使用的展示变量警告，没有语法错误。看板已用现有MCP刷新，仍有此前工作树16项缺失引用警告，与本轮116项数值检查分开记录；未恢复任何归档。此前未审查的距离标签修补候选未用于本轮案例，也未接入本轮有效坐标修复。
