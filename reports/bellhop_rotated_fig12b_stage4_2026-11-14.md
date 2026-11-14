# Bellhop旋转Fig12b：显式多次反射稀疏资格（2026-11-14）

当前仅有稀疏复分量资格。完整网格、反射阶数截断收敛、独立FFP比较和最终图册均未执行，不授予最终GO。

## 实现及合同

在独立validation overlay新增mode5，显式请求N次平面反射，按源初始向上/向下分别交替压力释放海面和有损流体底。第j墙的展开距离为`R_first+(j-1)D`。每次调用native Reflect2D一次，再绕该墙半转，变换位置和方向；p/q及面外积分不因坐标半转改变。Q表对应`z=D-abs(mod(zs+sg*r,2D)-D)`，横向接收坐标为`q=q0-sg*(-1)^N*x`。精确接收坐标、源c和积分合同沿用S3；物理复场固定使用`-conj(native)`，以非共轭转置及已保存排列逆映射，无插值/增益拟合。

底波数、密度及损耗与FFP假设相同；实际复反射系数逐墙审计。默认公共API/生产PE未改。mode5取消native累计Amp<.005停止门，保留每次反射本身的原生门；这是防止正常高阶贡献被静默截断，不是补声场。MaxN只在独立副本由100000增至600000，需逐次检查实际步数。人工上下边界仍非实际海面/海底。

新增：`prepare_bellhop_rotated_multifield_overlay`、`make_bellhop_rotated_multifield_config`、`write_bellhop_rotated_multifield_q`、`bellhop_rotated_multifield_reference_factor`、`validate_bellhop_rotated_multifield`、`audit_bellhop_rotated_multifield`。通用Q writer增加显式折叠阶数能力；旧默认ENV/SSP逐字节回归通过。外部进程wrapper增加peak working set与CPU时间字段，后续预算先导使用；此前作业未测得的峰值不追填。

## 失败与局部修复

原10阶、Q20m/步10m先导只命中第一墙，PRT fatal而进程退出0，日志门正确拒绝。`20261114_multifield_sparse_01`和`multifield_miss_diag_01`保留。新增`.mwstep`记录显示第二墙尚差2.9895e-8m时，原生最小步长规则将剩余步扩大成1e-4m并越墙。原因是Q节点与墙重合，两种候选限步只差浮点舍入，`h==hWall`未满足。修复只在RotQEnabled时阻止放大已经位于墙前的有限wall-clipped步；命中门1e-8m和snap门未放宽，未移动接收点或填补压力。

修复后Q20m/步10m的10/12阶命中数正确，但相对展开无反射Q场乘解析反射率的复差为9.1303e-4/1.1209e-3，未过3e-4门。该失败保留在`sparse_02`，未替换为PASS。仅把长路径步长减为2m，得到2.0536e-6/2.5697e-6，门不变；其Q20m仍只是反射计数/截断控制，不是连续Munk相位资格。

## 已执行证据

| 项目 | 结果与设置 |
|---|---|
| 1–3阶，上/下 | z=[0,1000,3000,4900,5000]m、x=[-1,0,1]m；Q5m、步2m、301束；展开对照复差4.03e-8…8.61e-8 |
| 上10/下12阶 | 同15个点，Q20m、步2m、101束；复差2.05e-6/2.57e-6；最高实际步数30201/33004（容量600000） |
| 低振幅保留 | 上10阶最大累计Amp=.002418083；下12阶=.000724712，不因旧.005门丢弃 |
| 原始数据独立审计 | 112门全通过：逐墙位置/计数、复系数、镜反射、proper rotation、q连续、原SHD/MAT与物理排列、实际线性Q时延和面外积分 |
| 独立审计最差值 | 所选合格配置：线性Q相位残差5.8983e-4rad（门.01）；积分相对差3.35e-10；逐墙复反射率差7.97e-12 |
| 关闭回归 | 默认ENV/SSP SHA逐字节相同，关闭扩展后与冻结官方直接Q复场在原1e-4门内通过 |

原始依据：`20261114_multifield_sparse_02`的案例1–6、`sparse_03`的案例7–8、`multifield_audit_02`和`multifield_disabled_02`。独立审计01故意包含粗步长失败长路径时记录4失败，未删除。关闭回归01因选错S1数据容器在加载时停止；02使用已冻结的overlay_validation_01重建并通过，未把前者记作数值失败或PASS。

编译输入`20261114_overlay_multifield_04`；binary为`20261114_build_multifield_03/bin/bellhop_rotated_multifield_2020.exe`，SHA256 `93ab5754e45525ab0bff0e5490035c74ffedac83472e9c298b3690eb9928c107`。历史build_manifest的purpose文字仍沿用旧单墙模板，不能当资格说明；本报告、源码SHA和阶段checkpoint界定实际范围，生成器已更正后续构建的purpose。

本阶段数值均本机现有MCP/串行，不涉及SSH、第二MATLAB或parpool。稀疏轮01/02/03约11.68/9.95/3.79s，关闭02约3.34s，独立审计01/02约2.43/2.36s；未测得任务总峰值内存或准确全任务总时间，不推算。一次RAM观察约239MiB，停止派发；清理本任务已创建的加载变量后回到约1.03GiB，未关闭用户程序。

## FFP读取与剩余门

只读检查本地FFP J07，SHA256 `8e123af0beff9904289a5419fdb82f74473beda075792b7b2f9819bd133ce3be`与原artifact清单一致。头部/配置确认10981×201接收网格、5000/4950m、10kHz、相同Munk/流体底和pressure_1m定义；保留小型reference metadata。尚未进行复场对照或独立重跑FFP，文件SHA符合不等于Bellhop已通过FFP。

下一步只允许小型容量/实际成本先导与完整网格程序准备，再提出可审阅的本地完整场及收敛总预算并等待用户明确批准。须独立检查束数、步长、Q采样和反射阶数，不以当前粗Q高阶计数控制代替它们；源点只允许真实奇异掩膜，源平面其余点由独立左右图表求解。不得复用旧FFP900s预算。没有生产物理/API变化，因此未重跑无关的标量/宽带通信链；本阶段必要的validation关闭回归已执行。