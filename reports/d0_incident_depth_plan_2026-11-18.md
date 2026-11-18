# D0：发射深度与海面入射复声场技术方案

2026-11-18。本轮状态为 **PREPARED / WAITING FOR USER MANUAL START**，不是D0声学计算通过。仅完成代码、冻结输入、远程执行包、静态检查与无正式声场dry-run。SSH连接、远端启动、正式声场及反射求解均为0；不提交或推送Git。

## 1. 输入证据与复用边界

按当前源码核查 `AGENTS.md`、任务表与 `PROJECT_CONTEXT.md` 顶部索引，随后定向读取：

- [K0最新审查](pm_wind_cutoff_k0_review_2026-11-18.md)：正式45组表面统计及保存审查完成；没有声学结论。
- [R1审查](r1_bie_v2_final_review_2026-11-16.md)、[R2最终报告](r2_saved_closeout_final_2026-11-17.md)：限定Tx100 m、0.10 m二维条带、4 kHz、原PM/C2、411点反射场；不向1490 m外推误差数字。
- [P0/P1资格报告](../results/comparisons/piston_2d_single_frequency/20261112_remote_p0p1_01/report.md)：10/20 cm两孔径P0及4 kHz指定P1资格。连续角谱传播资格不等于生产WAPE集成资格。

原样复制、逐字节SHA冻结 `baffled_strip_piston_field_vertical`、`piston_independent_angular_spectrum_vertical`、`apply_pe_channel_phase_reference_vertical` 和现有 `pm_sensitivity_io`。不修改这些实现、R1/R2/K0、公共信道或通信接口。来源代码、原输入、已读报告版本、Git HEAD与已有tracked工作区状态记录在执行包 `provenance/sources.json`；该包不是干净Git提交的替代物。

## 2. 固定物理模型

二维沿y无限长、总宽W=0.10 m的理想刚性障板条带，均匀同相向上法向峰值振速v0=0.001 m/s。x源中心0，z向下，发射深度顺序100/300/600/1000/1490 m。f=4000 Hz，c=1500 m/s，rho=1025 kg/m³；exp(-i omega t)，phase contract 2.0.0，输出为 `absolute_physical` 峰值Pa及Pa/m。

1500 m只作为目标水深背景；1490 m对应距海底10 m的安装意图。**不建立海底边界，不含向下辐射、海底或海面反射、镜像源、多次反射、外壳、背向辐射或结构散射。** 水体均匀无吸收，不含SSP、海流。此理想化模型隔离发射深度影响，不代表完整真实深海传播环境。z=0仅为评价面，不施加压力释放边界。

令d=z_tx-z>0，R=sqrt((x-xi)^2+d^2)，C=rho*omega*v0/2。直接使用已资格实现：

```
p  = C int[-W/2,W/2] H0^(1)(k R) dxi
px = -C k int H1^(1)(k R) (x-xi)/R dxi
pz = +C k int H1^(1)(k R) d/R dxi
```

不再添加载波，不对该物理场作共轭。v2检查只对显式移除纵向载波后的测试表示作往返恢复。这里没有DSP频响、H(f)扫频或CIR。

## 3. 平面评价、相位与归一化

所有深度使用x=-220:0.05:220 m，共8801个真实评价点。100和1490 m另计算x=-220:0.025:220，共17601点。主求积GL64，独立检查GL128；不使用R1/R2接收线替代440 m评价范围，不从稀疏接收点生成复场。

同时保存绝对p/px/pz、p/p(0)的复中心归一化、p/sqrt(mean(abs(p)^2))的正实整线RMS归一化。RMS指**空间RMS**，不是声压相量的时间RMS。RMS归一化保留物理中心相位；相对相位图仅用预先声明的中心参考，不拟合逐深度相位。

局部相位波数为kx_loc=Im(px/p)、kz_loc=Im(pz/p)。几何角theta_g=atan2(x,z_tx)，局部相位方向theta_loc=atan2(kx_loc,-kz_loc)，均相对向上竖直方向，带x符号；有限孔径下两者分别保存。|p|<1e-4全线最大幅度的点掩码。相位曲率为d(kx_loc)/dx，单位rad/m²；轴上等效相位前沿半径k/phi''单独输出，不能把它与海面几何曲率混用。

## 4. 独立角谱及窗口解释

使用P0的Weyl积分及同一速度孔径sinc权重，传播分支kx=k sin(theta)、q=k cos(theta)，倏逝分支kx=k cosh(t)、q=i k sinh(t)。新增 `d0_incident_spectrum_reference` 只把原被积函数作复合GL数值分段，以避免1490 m振荡使自适应积分默认区间上限失效；没有新声源物理模型。

传播theta分段数为ceil(k*(max(d)+max|x|)*(pi/2)/phase_span)，GL16/span8与GL32/span4交叉；倏逝t分8段，各用同阶GL，保留exp(-45)截止及压力尾项绝对界。每深度先与未修改P0自适应参考在同源近面三点核对，再在13个平面点、13个粗糙面点与孔径压力/解析梯度交叉。远场倏逝部分不能直接删除：近掠射传播端点与倏逝分支需相干相加。

连续全线角谱密度定义为

```
P(kx;d) = rho*omega*v0*W*sinc(kx*W/2)*exp(i*q*d)/q
q = sqrt(k^2-kx^2), Im(q)>=0
p(x) = int P(kx) exp(i*kx*x) dkx/(2*pi)
```

单位Pa m。传播分支模量在该均匀无吸收模型中不因深度变化，传播相位变化；固定有限观察窗所截获的场及其FFT会变化。这两种角谱不得混为一谈。q=0的掠射密度奇异点不赋有限值，不用全线密度平方积分宣称有限“能量”；图用theta端点距离1e-4 rad，积分仍用正则变换。

窗口FFT分别保存±200/±220 m的矩形与Hann，去除右重复端点，并恢复起始x的Fourier相位因子。保存FFT复值、相位掩码、传播带/带外标签、Parseval及逆变换检查。**有限窗|kx|>k功率含截断泄漏，不能称为物理倏逝能量。** 倏逝物理由连续分支相干压力/梯度及尾项记录说明；分支范数比只用于列明审计点，非能量比例。

## 5. 平面波问题的可检验量

全线轴向平面波使用实际p(0)，不拟合幅度或相位。另在x0=0/100/200 m和宽1/5/10/40 m的局部patch，使用实际p(x0)及kx_loc(x0)：p_pw=p(x0)*exp(i*kx_loc(x0)*(x-x0))。保存复L2、最大包裹及未包裹相位差、幅度和相位误差分量。

平方误差按恒等式拆分：|p-p_pw|²=(|p|-|p_pw|)²+4|p||p_pw|sin²(delta_phi/2)。独立检查分解闭合，不用包裹相位小误认为累积曲率小。是否可采用平面波须结合patch、偏移及后续模型允许误差判断，不单以传播距离大作为资格。

物理预判（待D0数据检验）：孔径自身进入远场并不保证440 m评价窗为平面波。几何相位项k*(sqrt(d²+x²)-d)及局部二阶项k*x²/(2d)控制横向波前弯曲；二维扩散与有限条带sinc方向性控制幅度。100→1490 m预计改变有限窗内角度分布与曲率，但D0不能由此预报PS/BIE/FullFK反射误差。

## 6. 冻结K0辅助粗糙边界

固定正式回收 `pm_wind_cutoff_k0_20261117_05_retrieved_20261117_232303/output/seed260001_U6_N12`。U_PM=6 m/s（测量高度和平均时长未知），Seed260001、12模、L160 m、Kmax=12*2*pi/160=0.471238898038469 rad/m；原Fourier系数、原42–50 m C2嵌入与全部所需检查点产物SHA验证后复制。选择理由是与R1锚点直接衔接，不按后续声学结果挑面。

只读取 `surface.mat` 中原有abs(x)<=50 m的5121个样点，不重新生成海面，不插值、不平铺、不扩展周期。这里高程eta向上，实际深度z_s=-eta，用户记号p(x,eta(x))落实为p(x,z_actual=-eta)。保存p/px/pz和入水法向n=(-z_s',+1)/sqrt(1+z_s'^2)对应dn p；向空气的外法向梯度为其相反数。各深度使用同一原几何。

用冻结系数独立标量Fourier+C2重建几何；13点检查数值法向梯度和独立角谱。全目标±200 m粗糙海面范围**不足，登记D1独立待定工作**，不以C2归零区或原周期复制代替实际扩大海面。

## 7. 正式计算前冻结的检查门

| 检查 | 严格门 | 来源/含义 |
|---|---:|---|
|GL64/128压力与梯度；孔径/独立角谱；两阶角谱；原P0积分适配核对|相对L2 <1e-7|沿用P0源门，压力和梯度分列|
|五点差分x/z梯度、粗糙入水法向梯度|相对L2 <1e-6|沿用P0；h=1e-4及5e-5 m平面检查|
|v2物理表示往返|最大绝对差 <1e-10 Pa|沿用P0相位门|
|MAT/CSV|相对L2 <1e-12，NaN位置必须一致|采用既有quicklook精度要求|
|粗细共同点压力/梯度|相对L2 <1e-12|复用接收交错点精度|
|粗线幅度/展开相位分别spline预测细点；Hann FFT同波数带内粗细差|相对L2 <1e-4|新D0采样协议，数值门不比既有采样门宽；不把插值当独立源|
|相邻复场相位增量|最大值 <pi/4 rad|新增保守相位采样门|
|FFT Parseval/逆变换；平面波平方误差分解|相对或列明绝对残量 <1e-12|新增代数审计|
|FFT对11个直接DFT波数点|误差/全FFT范数 <1e-10|新增符号/窗口相位审计|
|冻结Fourier/C2几何|最大绝对差 <1e-11 m；坡度同数值门|复用K0重建精度|

矩形窗宽度/采样效应是保留的物理表示诊断，不强行要求不同窗口角谱相等。FIG全部数据线及轴标签重开后必须一致；另从重开FIG导出PNG，不要求跨平台像素逐字相等。正式PNG布局、色彩、掩码与解释的目视验收在Fetch后完成。

所有检查值先保存 `audit.json/audit_checks.csv`，数值门不通过停止后续深度；异常详情写diary及全局失败检查点。不能从进程退出0推断PASS。尚未执行的数值门目前均为UNRUN。

## 8. 执行、资源和停止机制

最终包 `results/experiments/d0_incident_depth/preparation_20261118_06/`，28项输入/代码SHA，tar SHA为 `a63dbbbd8f092fd8d85c04ceca8e3f8c65a23d21a6007f07aa5a9515cfd8c3dd`。

用户手动 `UploadAndStart` 才连接 `matlab-ubuntu-new`，只在 `/home/misc/cwx/d0_incident_depth_<unique-id>` 新目录内工作。执行 `/home/misc/MATLAB/R2025b/bin/matlab -singleCompThread -batch ...`；入口核R2025b、Unix、单计算线程和监督器manifest绑定。无parpool、无后台并行深度。

冻结资源政策：单进程；启动MemAvailable≥4 GiB、磁盘≥2 GiB；保守工作内存估计2 GiB并需可纳入70%可用量。每秒记录可用RAM和owned MATLAB PID RSS，保留启动可用RAM的30%。整个会话含MATLAB启动上限3600 s，每深度从started标记起900 s；这是停止上限，不是成本实测。内存、时限、信号或数值异常即停止，只终止监督器创建且身份已记录的进程组，保留真实退出码。RSS说明仅代表记录PID，不宣称进程组峰值。

Linux锁确保同D0工作线单会话；one-shot启动标记拒绝同包重复启动，已有output拒绝正式入口重跑。**无自动恢复或重试**。每深度独立封存，已通过结果保留；中断后只能先Status/Fetch，任何后续恢复必须另行针对未完成部分设计，不重复已通过深度。后台任务不依赖Codex在线，不需要持续轮询。原K1及其他阶段未被本包授权。

## 9. 产物与本轮交付

每深度：`field.mat`（绝对/归一化复场、梯度、连续及有限窗谱、角度、相位/曲率）、`points.csv`、四组窗口复谱CSV、传播/倏逝密度CSV、`plane_wave.csv`、粗糙面MAT/CSV、代表深度fine.mat、独立分支MAT、指标JSON、数值审计、started及checkpoint。全局：完整physical_config.json、MATLAB diary、stdout/stderr、PID/命令/资源/真实退出码、总检查点与产物SHA。

六类PNG/FIG：绝对Pa、两种归一化幅度、相对相位、复角谱幅相、几何/局部方向、梯度及相位曲率；`depth_summary.csv`含压力、梯度、曲率及几何量。每张FIG重开核全部曲线和轴标签，另存重开PNG。远端生成 `D0_run_report.md` 为数据/审计回执；**中文最终物理解读报告在用户Fetch后统一形成**，本轮不编造深度规律的实验结论。

本轮准备验证见[测试报告](d0_incident_depth_preparation_validation_2026-11-18.md)，用户操作见[手动命令](d0_incident_depth_manual_commands_2026-11-18.md)。生产实现及接口未变，故无需标量/宽带通信回归或修改vertical_comm_guide的生产公式。更新脚本导航、任务表、项目索引与看板；保留K0、R1/R2历史资格及全部准备失败。完成后停止，等待用户手动启动。
