# 单链路确定性通信基础验证（2026-10-14）

**解析平面宽带、已认证dual离散复增益、native/DSP/CIR与QPSK基础控制通过；通用频率插值与完整实用接收机尚未验收。**这些结果是基础设施准备，不是PM BER/PER、D2–D5或WP4正式验收。**G1 OPEN、WP2=false**。

[正式计划](deep_sea_channel_research_plan.md)、[本轮协议](deepsea_two_track_protocol_2026-10-14.md)、[原始MAT/JSON/CSV](../results/validation/deepsea_two_track_20261014/)。新增入口 `validate_deepsea_singlelink_readiness`；数据/噪声seed2610141/2610142，运行后恢复共享RNG状态，不修改生产modem或PE。

## 复场与时延契约

Tx1500/Rx5/c1500、q0/.10，direct/reflected/total采用单位exp(ikR)/(4πR)。`absolute_DSP=conj(native)`；`relative_DSP=absolute_DSP*exp(i2πfτd)`，τd=Ld/c仅由几何指定。原生复场、绝对DSP、相对DSP均完整保存，没有逐realization幅相配准。

垂直额外延迟6.666666667ms、干涉周期150Hz；q=.10额外延迟6.633361620ms、周期150.753125976Hz。反射相位差独立回算延迟，最大误差7.81e-17s。近零点不剔除：q0最小总增益3.5368e-7、谱条件比300；q.10细网格谱条件比约300.16。分量及direct尺度用于误差，不仅除以深衰落total。

## 宽带与CIR检查

4–8kHz 257/513含端点解析频点，DFT分别取256/512个独立bin，舍弃右端点；fc6000、Δf15.625/7.8125Hz，时窗64/128ms，Δt固定.25ms。直接路径置相对t=0；carrier相位保留，ifft(ifftshift)及FFT符号/归一化明确。

最大分量往返误差3.26e-16、Parseval约5.94e-16；时域补零后的原频点复现差≤2.08e-16，补零只增加频率插值点。固定带宽增加真实采样频点延长时窗，不提升约.25ms带限时延分辨率。矩形带限CIR的长旁瓣不是额外物理路径；没有裁掉尾部以制造验收。

**保留失败方法：简单线性插值。**在独立中间频点，以direct幅度尺度计，q0的257/513源网格误差5.27%/1.33%，q.10为5.22%/1.31%。较细网格仍不达到计划1e-3数值筛查提议；当前精确解析频点控制通过，不能由DFT往返或平滑曲线推断通用插值通过。后续真实宽带信道必须独立检验中间频点/窗口/尾部及通信影响，不将当前网格自动移交PM。

## QPSK与噪声控制

Gray QPSK Es=1、Eb=.5，AWGN复噪声功率1/(2γb)，每分量方差1/(4γb)，Es/N0=Eb/N0+3.0103dB。每点500000bit，与0.5erfc(sqrt(γb))作5标准差二项筛查；这是错误检测控制而非高精度BER曲线或置信区间认证。

| Eb/N0 dB | 错误bit | BER | 理论BER |
|---|---:|---:|---:|
|0|39514|.079028|.07864960|
|4|6281|.012562|.01250082|
|8|87|.000174|.00019091|

已知复增益无噪声控制零误码。解析平面两路径单载波QPSK离散控制Fs4000/Rs2000、每符号两采样，完美CSI/定时、频域均衡，q0/.10及两频网格四案均零误码。CP取完整N−1个带限离散CIR样本，线性卷积与循环模型最大差4.18e-16；这个保守CP只用于代数控制，不是实际包格式或净速率建议。

冻结dual的几何Tx(0,150.5,1500)/Rx(0,0,5)、振幅(.05,.025)m、K(.4,0)/(.2,.6)、phase(.2,.7)、C2 plateau8/support12，仅使用已认证4/6/8kHz各自的独立参考total复增益。三点的native/DSP和无噪声已知增益QPSK均通过；**没有dual宽带CIR、dual多径宽带BER或连续频段认证**。

## Provenance与未完成项

输出schema包括physics_tier、证据scope、源/接收几何、海面/support、三类复场、频轴/时延、seed、版本manifest、G0/G1/WP2与成本。`deepsea_readiness_guard` 对本轮artifact只允许readiness；训练、统计验收、PM_BER、PM诊断伪装readiness及非认证dual5000Hz共5项拒绝检查通过。它约束本轮新验证管线；没有改写或宣称审计了所有既有训练消费者。

未完成：RRC pulse及FSE/MMSE长度、实用pilot估计/同步、noisy频选两路径接收机、包/净速率、真实环境分布、PM统计与通信。它们仍是正式WP4或后续域认证任务；本轮基础控制不替代。可选2D→3D线源极限延期，没有新增G1阻塞门。

## 成本与验证

B整个确定性控制核wall2.049s，包含解析频响/CIR、1500000bit AWGN及flat/dual代数控制，未作暖运行重复benchmark。进程前/后快照2443.45/2505.55MiB，开始物理空闲2211.78MiB；任务峰值未知，不把快照差当峰值。环境MATLAB R2025b Update6，当前cwd保持Bellhop Munk。

数值断言及数据审计通过；冻结算子回归4/4通过。首次`runtests`因当前目录同名Bellhop脚本遮蔽而失败，原始输出保留；改用`matlab.unittest.TestSuite.fromFile`/`TestRunner.withTextOutput`通过，不改变cwd或重置path。新文件Code Analyzer和冻结哈希核对见 [final_audit.json](../results/validation/deepsea_two_track_20261014/final_audit.json)。生产接口未变；已执行6kHz标量参考及解析4–8kHz宽带通信控制，未重跑可能覆盖历史产物的旧通信demo。

启动 `setup_vertical_project` 后从任何cwd调用入口；输出拒绝覆盖。原始MAT保存频响/CIR与输入bit，AWGN计数、噪声功率和seed可重放。没有提交或推送。
