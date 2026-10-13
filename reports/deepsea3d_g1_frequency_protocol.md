# 三维G1：固定y方位4/8kHz受控验证协议

6kHz只读复用 `deepsea3d_g1_azimuth/azimuth.mat` 的y方位条件性PASS。
仅改变频率为4000/8000Hz，Tx=(0,150.5,1500)、Rx=(0,0,5)、q=.10、
α2.5、dual振幅(.05,.025)m、波矢(.4,0)/(.2,.6)、phase(.2,.7)、C2
plateau8/support12、中心(0,0)、c1500、点Green及exp(-iωt)压力释放约定不变。

1. 每频独立三维路径、Green、镜像及native/DSP控制，复用WP0同频平面记录，
   原1e-10代数/相位门；检查abs(P)=1/(4πR)。不要求不同频率复压力相同。
2. 每频L24/h=.125/N12,N16、M40真实配点及双向半网格低成本预检。
   相邻N、N−2、边界与半网格各<1e-4、最后五阶压力和系数范数递减；
   必要时补N20。未通过最大N20则INCONCLUSIVE，不直接放行FK。
3. 所选N独立L24/h=.0625、同网格N−4/N−2、粗→细；N16不足可一次
   受控升级同阶粗/细N20。随后L32同h/N，只有域变化失败而其他门通过时补L48。
   最大N20/L48，最细h=.0625，若仍不足停止，不强行扩展计算规模。
4. 在收敛候选的同L/h/N上另做M32真实配点及半网格评价，对照M40。
   两个残差变化及接收压力变化<原1e-4；原生逐阶项与范数重现<1e-12。
   用全部有限Fourier模态的max|kz|、高度界Σ|A|及Σ||b_n||给出
   exp评价尾界 `Σ||b_n|| exp(T) T^(M+1)/(M+1)!`，T=Σ|A|max|kz|。
   以真实粗糙面上的||u_flat||归一化，M32/M40尾界各<1e-4。
   这是有限Fourier指数求值界，不是FE、舍入或连续PDE余项界。
   M评价/尾界未定则INCONCLUSIVE，不默认M40充分，不追加无限评价阶数。
5. 收敛候选内部消项回放逐阶压力及范数差<1e-12，forcing<1e-8；
   指标不是数学条件数。到此参考ready后才做原FullFK五次：
   A16/h.05/trap、A16/h.025/trap主解、A16/h.05/gauss4、
   A12/h.025/trap、A24/h.025/trap。增量网格/Gauss/两个孔径各<1e-3。
6. 原生反射与非零增量复误差各<原1%为CONDITIONAL_G1_PASS，超门则
   CONDITIONAL_PHYSICAL_FAIL；参考或FK自身数值未定始终INCONCLUSIVE。
   幅度和原始相位另列，无幅相拟合，不用total或平面背景掩盖增量。
7. 结合冻结6kHz分别分析主误差和各自的数值变体范围：所选N、N−2、N−4
   部分和在3个coarse/fine/domain配置，与5个FK共45个误差组合。
   低阶部分和只作截断敏感性，不额外认证其真实边界；未定频点无物理比较范围。
   离散相邻频率范围不重叠才支持该区间内可辨识的误差升/降，重叠则未分辨。
   范围非置信区间或严格误差界，三频点不能认证连续4–8kHz宽带或一般频率规律。
   固定频点保存flat/rough total相消指标、total误差及参考自身变化，不搜索新频率。

入口 `validate_deepsea3d_g1_frequency(stage,frequency_hz)`，省略频率选首个未完案例。
init→每频重复next至evaluation→evaluate→observe→五次fk→audit；
next自动N/网格/域编排，失败转done。status只读，init拒绝覆盖。
最后 `audit_deepsea3d_g1_frequency_artifacts` 核对保存复数、分类、M评价、
唯一频率改动及冻结6kHz，并生成三频摘要/范围/趋势。
结果 `results/validation/deepsea3d_g1_frequency/`；逐个计算checkpoint。
成本含新FE、M评价、回放、FK核，评价尾界诊断另计；旧6kHz零新增求解。
不修改冻结求解器/生产/二维/旧数据，不做PM MC、FastFK、SSP、动态或WP2。
全G1目标域OPEN，WP2=false。

## 数值结果后修订：有界独立求积精化

8kHz初始Gauss4 h=.05护栏失败（增量变化1.200e-3），不能将其判成物理FAIL。
保留全部初始五解，仅允许一次同A16的Gauss4 h=.025，使用已有主梯形网格
分辨率，不扩参考域/阶数/网格或更改门限。仍不通过则INCONCLUSIVE。
新增stage为`'quadrature',8000`；初始失败及精化成本均保留。
8kHz趋势范围使用四个已通过的梯形变体及精化Gauss，共五个FK；
失败的粗Gauss不纳入可信物理范围，单列其失败证据。
