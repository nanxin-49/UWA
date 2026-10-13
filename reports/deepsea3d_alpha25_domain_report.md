# 三维 G1 α=2.5 计算域收敛与 FullFK 条件性物理对照

**结论：α=2.5 的本离散案例获得 CONDITIONAL_G1_PASS。**
独立高阶 Dirichlet FE 的 L24→L32 计算域检查通过，再执行 FullFK 自身网格、
独立求积及孔径检查，最后两项原生复误差均通过既定 1% 物理门。
全 G1 目标域仍为 **OPEN**，**WP2=false**；不认证连续振幅区间、PM总体或宽带。
验证日期 2026-10-13，共享 MATLAB R2025b Update 6。

## 冻结配置与方法

本轮只增加验证编排 `scripts/validation/validate_deepsea3d_alpha25_domain.m`，
未修改冻结 FE/FullFK、生产 PE、二维基准、既有验收数据或公共输出接口。
沿用通过的 `deepsea3d_alpha25_fine/fine.mat`，旧 L24 参考只读复用。
运行前保存 114 个文件哈希，运行后逐一核对；亦逐一核对既有 fine provenance
记录的 71 个代码/历史/细网格产物哈希，均一致。
预设顺序和停止规则见[协议](deepsea3d_alpha25_domain_protocol.md)。

Tx=(0,0,1500) m，Rx=(0,0,5) m，q=0，f=6000 Hz，c=1500 m/s；
dual 振幅(.05,.025) m、波矢(.4,0)/(.2,.6) rad/m、phase(.2,.7) rad；
C2 径向 plateau=8 m、support=12 m，海面 z=0、z 正向水下、压力释放。
物理海面是相同的紧支撑 η(x,y)，不随计算域放缩。源为原单位三维点 Green 源，
原生 exp(-iωt) 约定不变，无幅相对齐、拟合、滤波或模态删减。

独立参考沿用同伦 tη 的 Dirichlet 场展开：

`b_n = -[t^n]u_flat(tη) - Σ_(j=1)^(n-1) η^j/j! D_z^j b_(n-j)`，
`D_z = ifft2(i kz fft2)`，保留全部采样倏逝模态。

FE 的 L 是 Fourier 周期计算域半宽；h=.0625 m，N20，真实边界评价 M40。
FullFK 的 A 是积分正方形半宽，两者单独检查，不把域变化计入物理模型误差。
FullFK 保留完整两项 `−∫(G ∂n u_inc + u_inc ∂n G)dS`、真实高度和法向，
采用解析无限平面背景加紧支撑完整密度差的原实现。

## 阶段一：独立参考计算域通过

|检查|实际结果|原门槛|
|---|---:|---:|
|L24→L32 接收复增量变化|3.2420e-7|<1e-4|
|沿用 L24 粗→细网格变化 N20|1.7728e-6|<1e-4|
|L32 同网格 N16→N20|2.3628e-6|<1e-4|
|L32 同网格 N18→N20|9.2483e-8|<1e-4|
|L32 真实海面配点残差|1.5542e-9|<1e-4|
|L32 双方向半网格真实海面残差|2.9773e-7|<1e-4|
|源 Taylor 系数相对误差|5.5563e-12|<1e-8|
|诊断回放逐阶压力/系数范数相对差|0 / 0|各 <1e-12|

残差是实际 η 上对同一 Fourier 场作独立 M40 传播指数评价，半网格同时平移
x/y 各 h/2，不以自动闭合的递推恒等式代替真实边界检查。

|域|ΔP 原生复增量（m^-1）|P_ref 原生反射压力（m^-1）|半网格残差|
|---|---|---|---:|
|L24，复用|9.95271227953e-5 − i1.51696164851e-5|4.66517264192e-5 − i1.51696164849e-5|3.0597e-7|
|L32，本轮|9.95271239115e-5 − i1.51695838646e-5|4.66517275355e-5 − i1.51695838644e-5|2.9773e-7|

两域 N20 最后五阶压力项和边界系数范数均逐阶递减。h 归一化的第20阶边界
系数 L2 近似为 L24 **1.2214907816623e-11**、L32 **1.2214907816643e-11**，
比例 1.0000000000016；尾部没有随扩域放大。第20阶压力项幅度分别
1.6451149492563e-12 / 1.6451149492774e-12。
完整逐阶复数和范数见[reference_terms.csv](../results/validation/deepsea3d_alpha25_domain/reference_terms.csv)。

内部贡献范数和/结果范数最大指标 **3313.8522**，接收阶数求和指标 **15.6173**，
与 L24 一致。内部存在消项；该诊断量不是数学条件数，回放也不是第二个物理 oracle。
逐阶压力和系数范数与实际冻结求解完全匹配，结合尾部与真实边界门判定稳定。
M40 有限 Fourier 指数求值尾项绝对上界 **3.1265e-20 m^-1**；
它只控制指数评价截断，不含浮点误差，也不是 FE 或连续 PDE 余项界。

L24→L32 已过预设域门，因此未追加 L48。N16 部分和来自相同递推的严格低阶
截断，其接收阶数门已过，未追加独立 N16 边界求解。这里认证的是已检验的
L24/L32 离散经验域稳定性，不能外推为任意更大域的严格收敛定理。
数据：[reference_runs.csv](../results/validation/deepsea3d_alpha25_domain/reference_runs.csv)、
[reference_guards.csv](../results/validation/deepsea3d_alpha25_domain/reference_guards.csv)、
[消项](../results/validation/deepsea3d_alpha25_domain/internal_cancellation.csv)、
[指数评价检查](../results/validation/deepsea3d_alpha25_domain/evaluation_checks.csv)。

## 阶段二：FullFK 数值护栏与物理误差

只有 reference_ready=true 后调用 FullFK，共五次：A16/h=.05/trap、
A16/h=.025/trap（主 FK）、A16/h=.05/gauss4、A12/h=.025/trap、A24/h=.025/trap。
各 A 都覆盖半径12 m 的物理支撑；外围 η 及坡度为零，密度差为零。

|FullFK 自身检查，相对粗糙复增量|实际变化|原门槛|
|---|---:|---:|
|h=.05→.025，A16|1.0214e-8|<1e-3|
|Gauss4 与主 trapezoid|1.9052e-7|<1e-3|
|A12→A16，h=.025|0|<1e-3|
|A16→A24，h=.025|2.5613e-15|<1e-3|

孔径不敏感来自紧支撑密度差已完整覆盖，不能用于声称有限原始粗糙积分自身
近似无限面。`raw_truncated_reflect` 随 A 仍改变，保存作诊断，不用于物理认证。
主 FK 的完整密度差积分与 raw−rawflat 恒等差为5.0248e-19 m^-1，
五次最大同量级；积分贡献绝对值和/增量幅度约81.697。
证据：[fk_runs.csv](../results/validation/deepsea3d_alpha25_domain/fk_runs.csv)、
[fk_guards.csv](../results/validation/deepsea3d_alpha25_domain/fk_guards.csv)。

主 FK ΔP=`9.95232656254347e-5 − i1.51835760691515e-5` m^-1；
P_ref=`4.66478692493897e-5 − i1.51835760689729e-5` m^-1。
相对独立 L32/N20 参考，原生误差如下：

|量|复相对误差|幅度相对误差|原始相位差 rad|
|---|---:|---:|---:|
|反射 P_ref|2.9587e-4（0.02959%）|1.3449e-5|−2.9557e-4|
|粗糙增量 ΔP|1.4417e-4（0.01442%）|1.6935e-5|−1.4317e-4|

复误差定义 `abs(FK−REF)/abs(REF)`；幅度误差
`abs(abs(FK)−abs(REF))/abs(REF)`；相位差 `angle(FK*conj(REF))`。
误差主要在相位，未对齐。两项复误差各小于原 **1%** 门，得到该案例条件性 PASS；
没有以平面背景或总场掩盖粗糙增量。参考域变化3.24e-7与 FK 数值变化1.91e-7
均明显小于当前模型差异；这些是独立经验检查，不合成为严格误差上界。
完整精度见[physical_metrics.csv](../results/validation/deepsea3d_alpha25_domain/physical_metrics.csv)。

## 实际成本与数值护栏交付

|新增计算|wall s|CPU s|命名载荷|
|---|---:|---:|---|
|L32 N20 FE，1024²点|15.3692|97.9688|保存谱320 MiB（单复数组16 MiB）|
|L32 消项回放|11.8140|75.2188|重放相同递推；快照非峰值|
|五个 FullFK 核合计|3.8695|18.5313|最大估计块载荷28.140 MiB|
|本轮合计|31.0526|191.7188|不把上述载荷相加当峰值|

L24 旧参考谱载荷180 MiB；L32 为1.7778倍，FE 本次 wall 约为旧 N20 的1.8265倍，
不是经过重复统计的性能基准。FE 内存前后快照3188.55/3188.64 MiB，
回放3189.55/3189.55 MiB，FK 快照范围3189.55–3196.17 MiB。
全部是共享 MATLAB 会话前后值，**没有测量峰值内存**。CPU 可因多线程大于 wall。
成本只统计本轮保存的 FE、回放、FK 核，不包含 IO、测试、开发和报告，
旧基准成本不重复计入。详见[cost_summary.json](../results/validation/deepsea3d_alpha25_domain/cost_summary.json)。

新入口 Code Analyzer 零问题；既有分类回归4/4通过，覆盖参考未定、FK数值失败、
不能用背景掩盖增量、条件性PASS。本轮 audit 重新从复数压力计算门槛和分类，
核对 cfg、有限性、逐阶回放匹配、尾部、参考先行及全域OPEN/WP2=false。
额外保存 CSV 复数/相位恒等式与固定海面检查，冻结文件前后逐个哈希一致。
验证编排不改生产物理或通信接口，因此未重复无关生产通信、二维基准和旧验收求解。

## 复现及适用边界

在含已通过 L24 细网格证据的独立工作副本中：

```matlab
setup_vertical_project;
validate_deepsea3d_alpha25_domain('init'); % 拒绝覆盖已有本轮 checkpoint
validate_deepsea3d_alpha25_domain('next'); % L32/N20
validate_deepsea3d_alpha25_domain('observe');
state=validate_deepsea3d_alpha25_domain('reference');
% 仅当 reference 将 domains 扩为 [32 48] 时追加 next/observe/reference。
if state.reference_ready
    for job=1:5, validate_deepsea3d_alpha25_domain('fk'); end
end
validate_deepsea3d_alpha25_domain('audit');
```

实际本轮不触发 L48。已完成目录可用 `'status'` 只读、`'audit'` 重建表和核查；
不得覆盖旧证据。`domain.mat` 包含完整 cfg、复用基准、复压力逐阶项、部分和、
范数、残差、回放指标、FK、成本与分类。`provenance.json` 记录新脚本/协议/报告
及产物哈希，`freeze_checks.json`/`prior_provenance_checks.json`记录冻结完整性。

有限振幅独立 FE 和 Kirchhoff 近似在此固定案例的经验对照已完成，但没有严格
连续PDE余项定理，也没有证明普通FE的一般收敛半径。全目标域 G1 OPEN，
WP2 不自动放行；本轮未做更大振幅扫描、PM、FastFK、SSP、动态、多频或通信。
工作流沿用 MATLAB skill；技能归属文献沿用[冻结基准报告](deepsea3d_g1_finite_report.md)末尾。
