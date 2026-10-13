# 三维G1：固定y方位4/6/8kHz物理验证

**4/6/8kHz三个离散案例均为 CONDITIONAL_G1_PASS，6kHz冻结复用。**
复误差4→6kHz下降可分辨，6→8kHz下降仍受参考阶数敏感性限制，不能认定
完整三频的复误差单调规律。相位绝对误差的离散下降及6→8kHz增量幅度误差
上升有数值变体支持。全G1目标域仍 **OPEN**，**WP2=false**。
2026-10-13，共享MATLAB R2025b Update6；[预设协议及求积修订](deepsea3d_g1_frequency_protocol.md)。

## 固定配置及逐频平面控制

Tx=(0,150.5,1500)m、Rx=(0,0,5)m、q=.10、α2.5、c1500；
dual振幅(.05,.025)m、波矢(.4,0)/(.2,.6)rad/m、phase(.2,.7)，
C2 plateau8/support12m、海面中心(0,0)不变。没有海面旋转、平移或再生成。
仅将频率改为4000/8000Hz，6kHz只读复用[方位角验证](deepsea3d_g1_azimuth_report.md)。
原三维Green为exp(ikR)/(4πR)，z向下、压力释放、原生exp(-iωt)不变，无幅相拟合。

每频独立路径及点Green/镜像计算，核对同频WP0平面记录：
4/8kHz最大native/DSP相位闭合误差1.85e-12/3.69e-12，各过原1e-10门；
1/(4πR)幅度归一化误差至多1.11e-16，直接/反射场通过。
Ld/Lr=1502.55623854816/1512.50628097869m、镜像延迟6.633361620ms对频率不变，
复相位随频率正常改变。完整原生声压在各频目录flat_controls.csv。

## 独立参考：不沿用默认充分性

参考仍为冻结Dirichlet高阶场展开，保留全部采样倏逝模态，不改求解器或滤波。
每频从L24/h=.125/N12,N16预检，真实海面配点和双向半网格均检查。

|参考检查|4kHz|8kHz|原门|
|---|---:|---:|---:|
|N12配点残差|2.543e-6|1.448e-2，失败|<1e-4|
|N12→N16增量变化|2.324e-6，通过|1.322e-2，失败|<1e-4|
|选定阶数|N16|按需N20|最多N20|
|最终同网格N−4变化|2.324e-6|8.931e-5|<1e-4|
|最终N−2变化|6.213e-8|5.563e-6|<1e-4|
|h=.125→.0625变化|1.706e-7|4.341e-6|<1e-4|
|L24→32域变化|7.146e-7|1.668e-7|<1e-4|
|L32配点/半网格残差|1.490e-9 / 2.920e-7|5.405e-7 / 6.508e-7|各<1e-4|

两频最终L32/h=.0625，无需L48。8kHz的N16真实边界残差约1.386e-4，
而N20阶数变化8.931e-5已接近门限，**只通过本轮有界配置，不能称其裕量充足**。
最后五阶系数范数与接收项幅度递减；内部消项最大指标4/8kHz为654.602/3318.951，
逐阶回放压力/范数相对差均0/0。消项指标不是数学条件数或第二个物理oracle。

在每频候选的同L/h/N上另算M32，独立比较M40真实边界评价：

|M评价检查|4kHz|8kHz|
|---|---:|---:|
|配点/半网格残差变化|0 / 0|0 / 0|
|接收压力及逐阶项/范数重现差|0|0|
|M32尾界/真实u_flat范数|2.589e-11|4.232e-12|
|M40尾界/真实u_flat范数|3.491e-18|2.624e-19|

均过原1e-4评价门和1e-12重现门，未默认M40充分。
尾界为有限Fourier模态指数评价的Σ||b_n||exp(T)T^(M+1)/(M+1)!，
T=Σ|A|max|kz|；不含舍入、不证明FE或连续PDE余项。主参考参数及全部失败保留。

## FullFK独立求积失败与有界补强

每频参考全部通过才执行原FullFK初始五解：A16/h.05/.025/trap、
A16/h.05/gauss4、A12/24/h.025/trap。A独立于FE域L，全部覆盖12m物理支撑；
沿用解析无限平面+紧支撑完整粗/平密度差及完整两项法向Green积分。

|增量数值护栏|4kHz|8kHz|原门|
|---|---:|---:|---:|
|梯形h.05→.025|4.069e-9|1.471e-8|<1e-3|
|初始Gauss4 h.05|1.017e-9|1.200e-3，失败|<1e-3|
|按需Gauss4 h.025|未追加|1.683e-9，通过|<1e-3|
|A12→16 / A16→24|0 / 5.498e-16|0 / 2.849e-15|各<1e-3|

8kHz初始状态为**FK数值未定INCONCLUSIVE**，没有误判物理FAIL。
结果后明确修订协议，只精化一次Gauss至已有主积分分辨率h=.025，不改FE规模
或门限。原失败解、精化解、成本都保存；精化后才给最终物理状态。

## 三频原生物理误差与相消

|kHz|条件性G1状态|反射复误差|增量复误差|反射幅度误差|增量幅度误差|反射相位差rad|增量相位差rad|
|---:|---|---:|---:|---:|---:|---:|---:|
|4|PASS，新算|4.03405e-4|2.23494e-4|1.18894e-5|1.07097e-4|−4.03233e-4|−1.96153e-4|
|6|PASS，冻结复用|2.84305e-4|1.37459e-4|4.19888e-6|4.87045e-6|−2.84275e-4|−1.37373e-4|
|8|PASS，求积补强后|2.14579e-4|1.21678e-4|3.75655e-6|6.75351e-5|−2.14547e-4|−1.01219e-4|

复误差abs(FK−REF)/abs(REF)，幅度误差abs(abs(FK)−abs(REF))/abs(REF)，
相位angle(FK*conj(REF))。三频反射/非零增量各过原1%，无校准或背景掩盖。
4/8kHz反射误差 **0.04034%/0.02146%**，增量误差 **0.02235%/0.01217%**。

|kHz|平面total相消指标|实际粗糙total指标|total复误差|参考total变体变化|
|---:|---:|---:|---:|---:|
|4|1.005|1.729|3.360e-4|4.398e-6|
|6|1.703|1.254|1.717e-4|1.093e-6|
|8|4.793|1.061|1.096e-4|8.122e-5|

实际粗糙场未见强相消；8kHz的较大参考阶数敏感性仍不可忽略，不能把较小
total主误差当独立参考误差已消失。FK积分贡献求和指标45.107/79.974/144.587，
接收阶数求和指标4.182/11.329/38.885。无额外相消频率搜索或宽带总场认证。

## 哪些频率趋势可信

每频N/N−2/N−4的接收部分和在coarse/fine/domain上构成9个压力变体，
配5个通过的FK数值变体，共45个误差组合。8kHz失败粗Gauss不混入可信物理
范围，使用其精化Gauss；失败原始结果仍保留。低阶部分和是截断敏感性，
不另行认证其边界；范围不是置信区间或严格PDE误差界。

|含阶数敏感性的复误差范围|4kHz|6kHz|8kHz|
|---|---|---|---|
|反射|[3.98553e-4,4.03518e-4]|[2.82151e-4,2.85349e-4]|[2.11971e-4,2.82457e-4]|
|增量|[2.20806e-4,2.23557e-4]|[1.36418e-4,1.37964e-4]|[1.20199e-4,1.60182e-4]|

**4→6kHz两项复误差下降可分辨，6→8kHz两项范围重叠，完整三频复误差
单调下降尚未证实。** 两项相位绝对误差在相邻两段均有不重叠的下降范围；
增量幅度误差4→6下降、6→8上升也可分辨。应分量报告，不概括为“越高频越准确”。
这只是固定海面、几何及三个离散频点的经验结果，不推断一般频率规律或连续4–8kHz带宽。

## 成本、冻结与复现

|新增频点|FE含M32求解数|FK求解数|核wall秒|核CPU秒|
|---:|---:|---:|---:|---:|
|4kHz|5|5|44.554|273.453|
|8kHz|6|6，含失败与精化|57.923|342.188|
|合计|11|11|**102.477**|**615.641**|

M尾界诊断另耗wall **.345s**、CPU1.406s；核成本含FE、M32复算、内部回放及FK，
不含IO、测试、开发/报告。6kHz新增FE/FK求解0。
4/8kHz最大FE谱数组载荷256/320MiB，FK估计块载荷最大28.140MiB；
进程前后快照约3210.08–3226.65MiB，**没有测峰值内存**。

冻结201文件前后哈希一致；新入口/审计helper通过Code Analyzer，分类回归4/4通过。
两频实际FK提前调用均由ReferenceGate拒绝；保存CSV原生复数、幅相恒等式、
M评价及分类审计通过。只新增验证编排，无公共接口、新物理公式或生产通信变更，
未重复无关生产/二维求解。项目状态、技术指南、计划及脚本索引已同步。

```matlab
setup_vertical_project;
validate_deepsea3d_g1_frequency('init'); % 拒绝覆盖已有checkpoint
for f=[4000 8000]
    state=validate_deepsea3d_g1_frequency('status');
    c=state.cases([state.cases.frequency_hz]==f);
    while ismember(c.phase,["preflight","upgrade","fine","domain"])
        state=validate_deepsea3d_g1_frequency('next',f);
        c=state.cases([state.cases.frequency_hz]==f);
    end
    if c.phase=="evaluation",state=validate_deepsea3d_g1_frequency('evaluate',f);c=state.cases([state.cases.frequency_hz]==f);end
    if c.phase=="observe",state=validate_deepsea3d_g1_frequency('observe',f);c=state.cases([state.cases.frequency_hz]==f);end
    if c.reference_ready
        for job=1:5,state=validate_deepsea3d_g1_frequency('fk',f);end
        c=state.cases([state.cases.frequency_hz]==f);
        if c.fk_guards.gauss>=1e-3,validate_deepsea3d_g1_frequency('quadrature',f);end
    end
end
validate_deepsea3d_g1_frequency('audit');
audit_deepsea3d_g1_frequency_artifacts;
```

[完整状态](../results/validation/deepsea3d_g1_frequency/frequency.mat)、
[三频摘要](../results/validation/deepsea3d_g1_frequency/frequency_summary.csv)、
[频率趋势](../results/validation/deepsea3d_g1_frequency/frequency_trends.csv)、
[误差范围](../results/validation/deepsea3d_g1_frequency/numerical_error_envelopes.csv)及两频子目录
保留原始复数、全部失败、评价/消项/护栏和成本；provenance记录代码、协议、报告及产物哈希。
不开展PM MC、FastFK、SSP、动态或WP2；全G1 OPEN/WP2=false。
工作流沿用MATLAB skill，归属文献沿用冻结有限振幅基准报告。
