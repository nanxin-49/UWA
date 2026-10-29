# 面向深海近垂直水声链路的信道估计与自适应均衡方法研究

**硕士学位论文选题报告**

| 项目 | 内容 |
|---|---|
| 培养单位 | 待填写 |
| 学科专业 | 待填写 |
| 研究生及学号 | 待填写 |
| 指导教师 | 待填写 |
| 报告日期 | 2026 年 10 月 |

## 目录

1. 选题背景与研究意义
2. 国内外研究现状综述
3. 课题研究内容、技术方案与预期创新
4. 前期研究基础与可行性分析
5. 论文工作计划与进度安排
6. 参考文献

## 插图清单

- 图 3-1 课题研究技术路线
- 图 4-1 典型粗糙海面样本的几何形态与坡度
- 图 4-2 二维粗糙海面接收线的声场及复残差
- 图 4-3 三维固定粗糙海面离散频点的反射复场与预测残差

## 附表清单

- 表 3-1 三级递进信道验证数据集及研究用途
- 表 3-2 单载波接收机拟采用的初始实验参数
- 表 4-1 二维固定粗糙海面样本的同源复场误差
- 表 4-2 三维固定粗糙海面四个独立频点的重构误差
- 表 4-3 解析通信链路的 QPSK 噪声控制结果
- 表 5-1 论文工作总体日程安排

## 符号和缩略语说明

| 符号或缩略语 | 含义 |
|---|---|
| $H(f)$、$h_{\mathrm{bb}}[n]$ | 复频率响应、通信基带等效离散冲激响应 |
| PM | Pierson–Moskowitz 海浪谱 |
| PE、FullFK、BIE | 抛物方程、完整 Kirchhoff 积分、边界积分方法 |
| CIR | 信道冲激响应 |
| LS、LMMSE | 最小二乘估计、线性最小均方误差估计 |
| NLMS、RLS | 归一化最小均方算法、递归最小二乘算法 |
| FS-LE、DFE | 分数间隔线性均衡器、判决反馈均衡器 |
| BER、PER | 比特误码率、数据包错误率 |

## 第1章 选题背景与研究意义

### 1.1 选题背景

深海观测设备向近海面平台回传数据时，声信号经历长距离传播，并在海面反射、散射后与直达波叠加。对发射深度超过 1000 m、接收深度约 3—10 m 的近垂直单链路，直达波与海面路径的传播距离相近但不相等，二者的复相位差随频率变化。相长和相消干涉可形成频率选择性衰落；粗糙海面又使反射能量在角度与时延上展宽。平台运动、海面演化及收发端时钟偏差还会使等效信道随时间变化。这些效应共同影响相干调制中的同步、导频估计和符号判决。[1-3]

本课题以深海向近海面的单载波 QPSK 通信接收机为研究对象。代表性数值几何采用发射深度 1500 m、接收深度 5 m、水平偏移约 150.5 m；第一阶段采用均匀声速和压力释放海面，以突出直达波与海面散射对接收算法的影响。4—8 kHz 作为候选研究频带，用于建立声场、基带信道和导频接收机之间的统一描述。上述参数是受控研究设置，后续将通过接收深度、横向偏移及海面参数变化考察算法的适用范围。

现有相干水声通信研究表明，有限带宽、长多径和相位变化要求信道估计与均衡协同工作。[3,16] 对近垂直链路而言，几何传播关系可提供有用的时延线索，但粗糙海面散射使真实响应不必严格稀疏。研究因此需要同时解决两个问题：如何从复声场构造相位与时延一致的通信基带信道，以及如何在有限导频和动态扰动下利用结构信息抑制估计误差与均衡器的错误传播。

### 1.2 研究意义

第一，建立物理声场到通信基带的跨层转换机制。声学数值模型通常输出固定频率处的原生复声压，而接收机使用的是经过发射成形、接收滤波和时间参考定义后的离散信道。明确复相位约定、公共时延原点、频率采样和无混叠时窗，能够把声场误差传递到估计误差与通信性能，形成可复核的算法研究输入。

第二，在有限导频开销下利用近垂直几何结构。直达与主要海面反射的传播时延可以由几何关系约束，为导频估计提供先验范围。将几何范围与开发数据中的能量分布结合，采用柔性延迟支撑而非硬截断，有望在减少估计方差的同时保留粗糙散射形成的低能量长尾。这一问题对导频长度、净数据率和接收机复杂度具有直接影响。[4,10-11]

第三，改进深衰落与时变条件下的自适应均衡。判决引导更新能够在导频结束后继续跟踪信道，但错误判决可能反过来污染均衡器参数。把导频阶段得到的估计不确定度、接收信号的预测残差和 QPSK 判决距离联合用于更新控制，可为低复杂度接收机提供抑制突发错误传播的途径。[3,6,17-20]

### 1.3 研究目标

本课题拟完成以下工作：构建声学原生复场到单载波通信基带信道的相位与时延一致转换；建立解析、固定粗糙海面和受控统计扩展的三级递进验证数据集；设计几何柔性延迟支撑导频估计算法；设计融合估计置信度与判决可靠性的自适应均衡方法；在统一导频、噪声和测试样本条件下，以信道估计误差、BER、PER、净速率、恢复时间及计算成本评价算法。

## 第2章 国内外研究现状综述

### 2.1 水声信道物理传播建模与粗糙海面散射特性研究现状

水声传播兼有频率相关损耗、多径和时间变化。Stojanovic 与 Preisig 归纳了传播模型和统计刻画方法；Qarabaqi 与 Stojanovic 将物理传播规律和随机局部位移相结合，为通信系统提供计算效率较高的信道模型。[1-2] 面向通信波形与接收机评估的开放模型研究进一步强调模型输入、场景和评价协议的可比较性。[7] 这些研究使物理模型与通信性能评价建立联系，但声场复响应在不同源归一化、相位约定及频率采样条件下转为离散基带抽头时，仍需逐项定义转换规则。

在海面建模方面，Pierson–Moskowitz 谱为充分成长风浪提供了经典谱形式。[8] 对有限粗糙海面，可通过谱截断和随机相位生成受控实现，并用 PE、Kirchhoff 积分或边界积分方法分析声场。Kachoyan 与 Macaskill 对任意粗糙压力释放边界的散射问题进行了研究。[9] 不同数值方法的源定义、边界处理及空间采样会影响复声场比较。面向通信算法时，还需要在统一发射条件下保存宽带复相位、时延参考与未见频点误差，使传播计算能够成为接收机实验的数据基础。

### 2.2 水声通信导频辅助与结构化信道估计研究现状

导频辅助估计通常把已知符号构成的观测矩阵与未知信道抽头联系起来。LS 方法不依赖统计先验，但在导频不足或矩阵病态时方差较大；LMMSE 可利用协方差降低误差，其效果则依赖训练分布与应用环境的一致性。水声信道中存在少量显著到达时，稀疏估计可降低有效未知维数。Berger 等在多载波水声系统中比较了子空间、压缩感知与 LS 估计；Stojanovic 及 Carrascosa 与 Stojanovic 分别研究了 MIMO-OFDM 系统的自适应估计和联合数据检测。[4-5,10]

模型驱动估计进一步关注稀疏程度不确定的问题。[11] 国内研究中，魏昕等利用帧头估计初始化均衡器并选择有效抽头，为训练段与跟踪段衔接提供了实例。[17] 现有方法多在特定浅海、OFDM 或多天线条件下评价。对于近垂直粗糙海面单载波链路，主要到达具有几何结构，散射尾部又可能跨越多个时延抽头；如何在降低估计方差与保留长尾能量之间取得平衡，是本课题的估计问题。

### 2.3 相干水声通信自适应均衡与抗误码传播技术研究现状

Stojanovic 等将载波同步、分数间隔判决反馈均衡和自适应参数更新结合，展示了相干水声接收机的基本结构。[3] NLMS 更新简单，RLS 对相关输入的适应通常更快，但需要更多运算和稳定性控制。Pelekanakis 与 Chitre 比较了稀疏自适应滤波器在水声均衡中的表现；长期浅海试验则说明噪声与环境变化会影响迭代接收机的实际排序。[6,12]

为减轻判决反馈的误差传递，国内研究提出了间接自适应均衡、软迭代信道估计、Kalman 与 Turbo 结合的均衡，以及叠加训练序列的频域迭代接收方法。[17-20] 利用时变相关性的联合信道估计与数据检测也得到研究。[13] 近垂直海试中的软迭代方案还结合了时变横向滤波、相位旋转模型和锁相环。[18] 在运动条件下，共同时间伸缩和非均匀多普勒也需要不同层次的补偿。[14-15] 上述工作为联合接收提供了基础；本课题将聚焦单载波接收机，研究仅依靠导频和接收端可观测量的可靠性控制更新，并量化深衰落、先验失配与受控时变下的误码传播及跟踪代价。

## 第3章 课题研究内容、技术方案与预期创新

### 3.1 总体研究框架与技术路线

研究以单发单收、近垂直、单载波 Gray QPSK 链路为对象。声传播仿真负责提供分级信道数据和物理误差参照，核心研究对象是接收端的导频估计与自适应均衡。工作按“声场到基带转换—导频估计—可靠性控制均衡—独立通信评价”展开：先使各级数据具有统一的相位、时延和功率定义，再在相同数据与导频预算下比较常规方法和拟研究方法。

![课题研究技术路线：分级信道数据、导频估计、自适应均衡与独立评价](../results/visualization/graduate_proposal_figures/research_route.svg)

**图 3-1 课题研究技术路线。** 解析控制和声学样本经跨频相位及时间窗检查后形成接收机输入，导频估计与均衡算法在独立样本上评价；不同数据层级对应不同的物理与统计问题。

### 3.2 研究内容一：宽带声场基带转换与分级数据集构建

声学模型采用 $\exp(-\mathrm i\omega t)$ 时间约定时，设同一源归一化下的原生复响应为 $H_{\mathrm{native}}(f)$。与通信接收机采用的正向时间约定一致的绝对频响以及以直达径时间 $\tau_{\mathrm d}$ 为原点的相对频响定义为

$$
H_{\mathrm{DSP,abs}}(f)=H_{\mathrm{native}}^{*}(f),\qquad
H_{\mathrm{DSP,rel}}(f)=H_{\mathrm{DSP,abs}}(f)
\exp(\mathrm i2\pi f\tau_{\mathrm d}).
$$

在载频 $f_c$ 附近，发射成形滤波器 $G_{\mathrm{tx}}(\nu)$ 和接收匹配滤波器 $G_{\mathrm{rx}}(\nu)$ 与传播频响共同形成等效基带响应

$$
H_{\mathrm{bb}}(\nu)=G_{\mathrm{rx}}(\nu)
H_{\mathrm{DSP,rel}}(f_c+\nu)G_{\mathrm{tx}}(\nu),\qquad |\nu|\le B/2.
$$

对 $N$ 个不重复、等间隔频点 $\nu_k$ 采样，频率间隔为 $\Delta f$。将中心化频率序列按离散傅里叶变换的索引顺序重排后，基带抽头由

$$
h_{\mathrm{bb}}[n]=\frac{1}{N}\sum_{k=0}^{N-1}
H_{\mathrm{bb}}(\nu_k)\exp\!\left(\mathrm i\frac{2\pi kn}{N}\right),
\qquad n=0,\ldots,N-1
$$

给出。相应周期时窗为 $T=1/\Delta f$，离散时间间隔为 $\Delta t=1/(N\Delta f)$。直达与散射路径的有效延迟跨度、成形滤波尾部及保护间隔需落在该时窗内，并以未见频点复响应检查频域重建误差。采用公共 $\tau_{\mathrm d}$，可使不同接收位置与不同声学方法的时延保持可比。

**表 3-1 三级递进信道验证数据集及研究用途**

| 层级 | 信道构成 | 算法研究用途 | 当前基础与后续建设 |
|---|---|---|---|
| 一级：解析基准 | 平面、直达及已知双径的解析复频响 | 验证相位转换、抽头恢复、估计器代数及噪声下渐近行为 | 已有解析频响、时窗和 QPSK 控制；完善统一接收帧 |
| 二级：三维固定粗糙海面 | 经高阶数值参考检查的固定海面离散复响应 | 在真实散射相位、相消频点和数值误差条件下检验估计与均衡机制 | 已有四个独立新增频点和预测频率网格；继续检查带内频点与时窗后形成宽带接收数据 |
| 三级：统计与时变扩展 | 参数化受控粗糙海面实现、连续或明确规定的时变过程 | 检验支撑先验失配、动态跟踪及跨环境泛化 | 按独立海面实现建立开发与测试划分，并记录模型适用域 |

跨层误差传递实验以同一参考频响产生接收波形，保持发射比特、噪声实现及接收机结构不变，只更换估计器或均衡器可见的设计频响，比较由此产生的估计和 BER/PER 差异。第三级中的受控统计和时变数据用于算法鲁棒性研究，其物理适用范围由独立参考与场景定义确定。

### 3.3 研究内容二：几何柔性延迟支撑导频估计

在导频段，令 $\boldsymbol y\in\mathbb C^M$ 为匹配滤波后的接收向量，$\boldsymbol X\in\mathbb C^{M\times L}$ 为由已知导频、成形脉冲及采样时刻构成的卷积矩阵，$\boldsymbol h\in\mathbb C^L$ 为待估计的基带抽头，则

$$
\boldsymbol y=\boldsymbol X\boldsymbol h+\boldsymbol n.
$$

以发射端深度 $z_{\mathrm t}$、接收端深度 $z_{\mathrm r}$、水平偏移 $\rho$ 和声速 $c$ 表示平面控制几何，直达与镜像反射路径长度分别为

$$
L_{\mathrm d}=\sqrt{\rho^2+(z_{\mathrm t}-z_{\mathrm r})^2},\qquad
L_{\mathrm r}=\sqrt{\rho^2+(z_{\mathrm t}+z_{\mathrm r})^2},\qquad
\Delta\tau_0=(L_{\mathrm r}-L_{\mathrm d})/c.
$$

在 $z_{\mathrm t}=1500\,\mathrm m$、$z_{\mathrm r}=5\,\mathrm m$、$\rho=150.5\,\mathrm m$、$c=1500\,\mathrm{m/s}$ 的代表几何下，平面反射相对直达的理论延迟约为 $6.633\,\mathrm{ms}$。对粗糙海面有效反射区 $\Omega_{\mathrm{eff}}$ 中的表面点 $\boldsymbol p$，候选超额延迟写为

$$
\Delta\tau(\boldsymbol p)=
\frac{\|\boldsymbol r_{\mathrm t}-\boldsymbol p\|
      +\|\boldsymbol r_{\mathrm r}-\boldsymbol p\|-L_{\mathrm d}}{c},
\quad \boldsymbol p\in\Omega_{\mathrm{eff}}.
$$

以 $\Omega_{\mathrm{eff}}$ 上 $\Delta\tau(\boldsymbol p)$ 的最小值和最大值形成候选反射时延区间，再与开发集上平均抽头能量的累积分布共同确定核心支撑集 $\mathcal S$；其两侧设置宽度为 $b$ 个抽头的过渡区，以覆盖离栅格时延、滤波旁瓣和弱散射尾部。支撑集和权重在独立测试前固定。定义 $d(\ell,\mathcal S)=\min_{j\in\mathcal S}|\ell-j|$，则可采用平滑递增的权重

$$
w_\ell=w_{\mathrm{in}}+
(w_{\mathrm{out}}-w_{\mathrm{in}})
\left[\min\left(1,\frac{d(\ell,\mathcal S)}{b}\right)\right]^p,
\quad
0<w_{\mathrm{in}}<w_{\mathrm{out}},\ p\ge1,
\qquad
\boldsymbol W_{\mathcal S}=\operatorname{diag}(w_0,\ldots,w_{L-1}).
$$

核心区采用较小惩罚，过渡区逐渐增加，区外保持较大但有限的惩罚。加权估计的目标函数及闭式表达为

$$
\hat{\boldsymbol h}
=\arg\min_{\boldsymbol h}
\bigl\|\boldsymbol y-\boldsymbol X\boldsymbol h\bigr\|_2^2
+\lambda\bigl\|\boldsymbol W_{\mathcal S}\boldsymbol h\bigr\|_2^2,
\qquad
\hat{\boldsymbol h}
=\bigl(\boldsymbol X^{\mathrm H}\boldsymbol X+
\lambda\boldsymbol W_{\mathcal S}^{\mathrm H}\boldsymbol W_{\mathcal S}\bigr)^{-1}
\boldsymbol X^{\mathrm H}\boldsymbol y.
$$

实际数值计算通过线性方程组求解。设 $\boldsymbol A=\boldsymbol X^{\mathrm H}\boldsymbol X+\lambda\boldsymbol W_{\mathcal S}^{\mathrm H}\boldsymbol W_{\mathcal S}$，当噪声零均值、协方差为 $\sigma_n^2\boldsymbol I$ 且抽头模型准确时，其正则化偏差为 $-\lambda\boldsymbol A^{-1}\boldsymbol W_{\mathcal S}^{\mathrm H}\boldsymbol W_{\mathcal S}\boldsymbol h$，估计误差协方差为 $\sigma_n^2\boldsymbol A^{-1}\boldsymbol X^{\mathrm H}\boldsymbol X\boldsymbol A^{-\mathrm H}$。这两项分别刻画支撑惩罚引入的偏差与导频噪声导致的方差。若真实到达落在抽头栅格之间，还应在观测模型中加入离栅格建模残差，单独分析其影响。

实验将使用相同导频预算比较 LS、统一岭正则化、基于开发集协方差的 LMMSE、硬支撑截断和上述柔性支撑方法；分别改变支撑偏移、过渡区宽度、离栅格延迟、散射尾部能量及海面参数。以估计 NMSE、频带复响应误差、支撑外能量、BER/PER 和净速率评价方差下降与遗漏偏差的权衡。

### 3.4 研究内容三：置信度门控自适应均衡

接收端先采用分数间隔前馈滤波器，并以可选 DFE 处理剩余码间干扰。设 $\boldsymbol r_n$ 是每符号多个采样构成的接收向量，$\boldsymbol w_n$ 为前馈系数，$\boldsymbol b_n$ 为反馈系数，先前判决向量为 $\hat{\boldsymbol s}_{n-1}$，则均衡输出为

$$
\tilde s[n]=\boldsymbol w_n^{\mathrm H}\boldsymbol r_n
-\boldsymbol b_n^{\mathrm H}\hat{\boldsymbol s}_{n-1}
=\boldsymbol q_n^{\mathrm H}\boldsymbol u_n,
\quad
\boldsymbol u_n=
\begin{bmatrix}\boldsymbol r_n\\-\hat{\boldsymbol s}_{n-1}\end{bmatrix}.
$$

在导频段，期望符号 $s_{\mathrm{ref}}[n]$ 为已知发射符号；在判决引导阶段，$s_{\mathrm{ref}}[n]=\hat s[n]$ 为最近 QPSK 星座点。更新误差为 $e[n]=s_{\mathrm{ref}}[n]-\tilde s[n]$。固定步长 NLMS 和固定遗忘因子 RLS 构成两条基线。

判决引导阶段的可靠性仅由接收端可观测量计算。第一项是利用导频估计信道和当前暂定符号构造的接收波形先验预测残差 $\rho_e[n]=\|\boldsymbol r_n-\hat{\boldsymbol r}_n\|_2^2/(\hat\sigma_e^2+\epsilon)$；第二项是导频信道估计协方差的归一化迹 $\rho_h=\operatorname{tr}(\boldsymbol C_h)/(\|\hat{\boldsymbol h}\|_2^2+\epsilon)$；第三项是软输出到最近星座点的距离 $\rho_d[n]=|\tilde s[n]-\hat s[n]|^2/(\hat\sigma_s^2+\epsilon)$。其中 $\hat\sigma_e^2$ 和 $\hat\sigma_s^2$ 由导频段噪声与残差估计，$\epsilon>0$ 用于稳定归一化。

令 $a,b,c\ge0$ 为开发集选定的权重，数据符号的综合门控系数定义为

$$
\alpha[n]=\alpha_{\min}+
(1-\alpha_{\min})
\exp\{-a\rho_e[n]-b\rho_h-c\rho_d[n]\},
\qquad 0\le\alpha_{\min}<1.
$$

已知导频取 $\alpha[n]=1$；数据符号的异常残差、较大估计不确定度或较远星座距离均降低其更新权重。取 $\alpha_{\min}=0$ 时，可暂停明显不可靠的判决更新。加权 NLMS 候选为

$$
\boldsymbol q_{n+1}=\boldsymbol q_n+
\frac{\mu\alpha[n]}{\epsilon+\|\boldsymbol u_n\|_2^2}
\boldsymbol u_n e^{*}[n].
$$

为研究跟踪速度与稳定性的折中，还将考察受门控的变遗忘因子 RLS。记 $\bar\rho_e[n]$ 为导频或高置信度符号上的有界残差统计，$\operatorname{clip}(x,l,u)$ 表示截取到 $[l,u]$，则

$$
\lambda_n=\operatorname{clip}
\bigl(\lambda_{\max}-\kappa\alpha[n]\bar\rho_e[n],
\lambda_{\min},\lambda_{\max}\bigr),
\qquad 0<\lambda_{\min}\le\lambda_{\max}\le1,
$$

$$
\boldsymbol k_n=
\frac{\alpha[n]\boldsymbol P_{n-1}\boldsymbol u_n}
{\lambda_n+\alpha[n]\boldsymbol u_n^{\mathrm H}
\boldsymbol P_{n-1}\boldsymbol u_n},
\quad
\boldsymbol q_n=\boldsymbol q_{n-1}+\boldsymbol k_n e^{*}[n],
\quad
\boldsymbol P_n=\lambda_n^{-1}
\bigl(\boldsymbol P_{n-1}-\boldsymbol k_n
\boldsymbol u_n^{\mathrm H}\boldsymbol P_{n-1}\bigr).
$$

当 $\alpha[n]=0$ 时暂停系数和 $\boldsymbol P$ 的更新；对于持续的高置信度失配，通过后续导频重新估计并调整遗忘因子。门控参数、遗忘因子范围和重启条件均由开发集确定，在测试集固定。该机制使低可靠性判决对反馈滤波器的影响受控，同时保留可靠观测下的跟踪能力。

算法评价采用配对设计：同一海面实现、发射比特、导频位置和噪声种子供各方法使用，算法参数只在开发集选择，独立测试按海面实现或连续记录分组。主指标包括 BER、PER、错误突发长度、失配后恢复符号数、估计 NMSE、导频与保护间隔开销及在线运算量。设置“无门控—仅残差门控—残差与估计不确定度—三项联合门控”消融，并在静态、受控相位漂移、路径幅相渐变和突变条件下比较固定 NLMS/RLS 与候选方法。

**表 3-2 单载波接收机拟采用的初始实验参数**

| 参数 | 初始设置 | 研究作用 |
|---|---|
| 调制与载频 | Gray QPSK，6 kHz | 建立相干单载波基准 |
| 符号率与成形 | 2000 symbol/s，根升余弦滚降 0.25 | 名义占用带宽约 2.5 kHz |
| 帧结构 | 128 个导频、1024 个数据符号 | 比较估计质量与导频开销 |
| 接收采样 | 每符号 2 个复采样 | 支持分数间隔均衡 |
| 对照算法 | LS、岭正则化、LMMSE、固定 NLMS/RLS | 确定估计与均衡基线 |

表中参数为初始实验设计；成形滤波跨度、保护间隔及均衡抽头长度将依据频带与信道记忆在开发阶段确定。

### 3.5 拟解决的关键科学问题

（1）**相位一致的跨层转换。** 在有限频率采样下，如何同时保持声学原生复场、公共直达时延原点和成形滤波后的基带响应一致，并识别频带内相消和时窗混叠对接收算法的影响。

（2）**非严格稀疏散射下的估计权衡。** 当近垂直几何给出主要路径位置，而粗糙海面形成离栅格与低能量长尾时，如何选取柔性支撑权重，使有限导频下的方差降低不被支撑遗漏偏差抵消。

（3）**动态条件下的可靠跟踪。** 在深衰落、相位漂移与判决错误相互作用时，如何仅使用接收端可观测量调节均衡更新，同时兼顾突发误码抑制、变化后的恢复速度和在线计算成本。

### 3.6 本课题的预期创新点

（1）建立以通信损失为终端指标的声学复场—基带信道—接收算法跨层误差传递评价框架，使频域重建精度与 BER/PER、净速率的关系能够在同一实验协议下比较。

（2）提出几何引导的柔性延迟支撑导频估计方法，以平滑权重表示核心到达和散射尾部的不确定性，研究导频有限、离栅格到达和海况失配时的偏差—方差规律。

（3）提出融合预测残差、估计不确定度与星座距离的置信度门控自适应均衡方法，在单载波接收机中研究判决反馈错误传播与跟踪速度之间的可调节机制。

## 第4章 前期研究基础与可行性分析

### 4.1 已有软件平台与声学建模基础

课题组已建立 MATLAB 垂直水声信道和 MPSK 通信仿真框架，可分别输出直达、海面反射及总复频响，并保存频率轴与相位参考。现有通信评价环节能够在已知信道条件下完成峰值同步、频域 MMSE 均衡和误码统计，为导频接收机的实现提供了可复用的信号处理基础。前期解析平面和双径控制还完成了原生声学相位与 DSP 相位的转换、离散频率网格与时域往返检查。

二维声场研究采用固定 4 kHz 同源高斯波束，在平面、正弦及三个预先指定的 PM 粗糙海面实现上，对 PE、完整 Kirchhoff 积分与经验边界积分参考进行了配对比较。图 4-1 展示其中一个样本的有限海面嵌入几何及坡度；图 4-2 给出相同样本在 $z=3\,\mathrm m$ 的原 411 点接收线上的原生反射复场幅度和逐点复残差。

![典型 PM 粗糙海面样本的实际嵌入几何与坡度](../results/comparisons/two_dimensional_channel_closeout/20261026_s3a_seed260003_01/figures/F01_PM_geometry_v2.png)

**图 4-1 典型粗糙海面样本的几何形态与坡度。** 上图为原始周期 PM 剖面及平滑嵌入后的有效海面，下图为有效海面坡度；有效海面最大绝对起伏约 $0.431\,\mathrm m$。

![二维 4 kHz 同源高斯波束下的接收线声场与方法复残差](../results/comparisons/two_dimensional_channel_closeout/20261026_s3a_seed260003_01/figures/F03_method_line_v2.png)

**图 4-2 二维粗糙海面接收线的声场及复残差。** 上图为 PM 海面与同源平面条件下的反射复场幅度，下图为 PE、完整 Kirchhoff 积分相对边界积分参考的逐点复残差；接收线包含 411 个原生位置。

**表 4-1 二维固定粗糙海面样本的同源复场误差**

| 固定海面样本 | 接收支持 | 样本点数 | PE 相对 BIE 的复 $L_2$ | FullFK 相对 BIE 的复 $L_2$ |
|---|---|---:|---:|---:|
| 1 | 中心九点 | 9 | 0.109921 | 0.000561 |
| 1 | 接收线 | 411 | 0.127452 | 0.000763 |
| 2 | 中心九点 | 9 | 0.002090 | 0.000267 |
| 2 | 接收线 | 411 | 0.093966 | 0.000400 |
| 3 | 中心九点 | 9 | 0.060506 | 0.000900 |
| 3 | 接收线 | 411 | 0.119206 | 0.000669 |

表中每一行均按该行原生接收支持上的 $\|u_m-u_{\mathrm{BIE}}\|_2/\|u_{\mathrm{BIE}}\|_2$ 计算，九点和接收线使用各自的分母。三个固定样本显示空间支持与海面实现均影响方法差异，为后续信道样本选择和模型误差传递实验提供了定量依据。

在三维代表性几何下，课题组已对一个固定双模粗糙海面开展高阶 Dirichlet 场展开参考与完整 Kirchhoff 积分检查。四个事先确定的新增频点分别为 4125、4750、6562.5 和 7250 Hz。图 4-3 显示反射项在复平面上的位置以及锁定预测与参考的复差。现有 65 点频率序列覆盖 4—8 kHz，间隔为 62.5 Hz，属于重建预测网格；其中的独立高阶数值参考仍以已计算频点为准。带内补充验证和时窗检查是形成宽带通信测试数据的下一步工作。

![三维固定粗糙海面四个离散频点的原生反射复场及预测残差](../results/validation/deepsea3d_65grid_four_20261027/four_point_reflection_complex_v3.png)

**图 4-3 三维固定粗糙海面离散频点的反射复场与预测残差。** 左图以实部和虚部表示四个频点的 Dirichlet 场展开参考与事前锁定预测；右图单独显示预测复差，量级为 $10^{-9}\,\mathrm m^{-1}$。

**表 4-2 三维固定粗糙海面四个独立频点的重构误差**

| 频率/Hz | 频点类型 | 反射复相对误差 $e_R$ | 总场复相对误差 $e_H$ | 相消指标 $C$ |
|---:|---|---:|---:|---:|
| 4125 | 相消压力点 | $3.232\times10^{-5}$ | $2.654\times10^{-4}$ | 17.05 |
| 4750 | 常规点 | $5.932\times10^{-5}$ | $8.566\times10^{-5}$ | 3.00 |
| 6562.5 | 相消压力点 | $1.312\times10^{-5}$ | $7.979\times10^{-5}$ | 12.63 |
| 7250 | 常规点 | $2.753\times10^{-5}$ | $1.329\times10^{-5}$ | 1.00 |

其中 $e_R=|R_{\mathrm{pred}}-R_{\mathrm{ref}}|/|R_{\mathrm{ref}}|$、$e_H=|H_{\mathrm{pred}}-H_{\mathrm{ref}}|/|H_{\mathrm{ref}}|$，$C=(|D|+|R_{\mathrm{ref}}|)/|H_{\mathrm{ref}}|$。相消指标较高时，总场分母变小，因而评估通信算法时还需同时观察直达分量尺度上的误差和接收符号性能。四点结果支持固定海面离散重构的阶段性研究，目标 PM 统计海面的连续宽带参考将按第三级数据建设逐步扩展。

通信基础方面，前期已完成解析平面与双径频响的相位转换、FFT/IFFT 代数和 QPSK 噪声控制。表 4-3 给出 Gray QPSK 在 AWGN 下的独立误码计数，每个信噪比点采用 $5\times10^5$ 个比特，并与解析误码率对照。这些控制确立了噪声功率、符号映射和计数方式，可作为后续导频估计与均衡实验的接收机基线。

**表 4-3 解析通信链路的 QPSK 噪声控制结果**

| $E_b/N_0$/dB | 错误比特数 | 仿真 BER | 解析 BER |
|---:|---:|---:|---:|
| 0 | 39514 | 0.079028 | 0.078650 |
| 4 | 6281 | 0.012562 | 0.012501 |
| 8 | 87 | 0.000174 | 0.000191 |

### 4.2 可行性分析与风险应对方案

**理论模型方面，**近垂直几何给出了直达与主要反射的可解释延迟，导频卷积模型、加权正则化及自适应均衡递推均可在统一的单载波基带框架中表达。解析双径数据可用于检验推导和实现；固定粗糙海面复响应则可用于识别相消与散射尾部对算法的影响。

**平台与数据方面，**已有声场求解、复频响保存、QPSK 噪声控制和已知信道均衡能力，使课题可以先完成接收机基线与算法机制实验，再把新获得的宽带和统计海面数据逐层纳入。二维同源比较提供方法差异的空间参照，三维四频点参考为高阶数值核对提供离散锚点。

**实验风险方面，**对高阶数值参考的计算成本，采用解析基准、固定三维样本和统计扩展的递进安排，并在每层设置频率重建、时间窗与独立样本检查；对延迟先验失配，保留无结构估计及硬支撑估计作为对照；对错误判决和跟踪迟滞，同时评价误码突发、恢复时间与计算成本。这样的实验设计使模型精度、先验条件与算法性能能够分别解释。

## 第5章 论文工作计划与进度安排

论文工作按正式开题后的相对月份组织。各阶段形成可检查的算法、数据和图表，并在进入下一层实验前固定前一层的参数及评价协议。

**表 5-1 论文工作总体日程安排**

| 时间 | 主要研究工作 | 阶段成果 |
|---|---|---|
| 第 1—2 个月 | 梳理水声信道估计与均衡文献，确定单载波链路、导频预算及三级数据协议 | 研究方案、基线定义和实验参数表 |
| 第 3—4 个月 | 完成原生声场到基带抽头的统一转换、解析双径控制及固定三维案例的带内检查 | 相位与时窗一致的信道数据、接收链基准 |
| 第 5—6 个月 | 实现 LS、岭正则化和 LMMSE 对照，研究几何柔性延迟支撑估计 | 估计误差、支撑失配和导频开销分析 |
| 第 7—8 个月 | 建立固定 NLMS/RLS 与分数间隔均衡基线，研究置信度门控更新 | 静态收敛、判决错误传播和受控时变消融结果 |
| 第 9—10 个月 | 开展独立环境样本的配对通信实验和跨层误差传递分析 | BER/PER、净速率、恢复时间与复杂度比较 |
| 第 11—12 个月 | 汇总研究结果，完成论文撰写、复核与答辩准备 | 学位论文及可复核的实验材料 |

## 参考文献

[1] STOJANOVIC M, PREISIG J C. Underwater acoustic communication channels: Propagation models and statistical characterization[J]. IEEE Communications Magazine, 2009, 47(1): 84-89. DOI: 10.1109/MCOM.2009.4752682.

[2] QARABAQI P, STOJANOVIC M. Statistical characterization and computationally efficient modeling of a class of underwater acoustic communication channels[J]. IEEE Journal of Oceanic Engineering, 2013, 38(4): 701-717. DOI: 10.1109/JOE.2013.2278787.

[3] STOJANOVIC M, CATIPOVIC J A, PROAKIS J G. Phase-coherent digital communications for underwater acoustic channels[J]. IEEE Journal of Oceanic Engineering, 1994, 19(1): 100-111. DOI: 10.1109/48.289455.

[4] BERGER C R, ZHOU S, PREISIG J C, et al. Sparse channel estimation for multicarrier underwater acoustic communication: From subspace methods to compressed sensing[J]. IEEE Transactions on Signal Processing, 2010, 58(3): 1708-1721. DOI: 10.1109/TSP.2009.2038424.

[5] STOJANOVIC M. Adaptive channel estimation for underwater acoustic MIMO OFDM systems[C]//2009 IEEE 13th Digital Signal Processing Workshop and 5th IEEE Signal Processing Education Workshop. Marco Island, FL: IEEE, 2009. DOI: 10.1109/DSP.2009.4785909.

[6] PELEKANAKIS K, CHITRE M. Comparison of sparse adaptive filters for underwater acoustic channel equalization/estimation[C]//2010 IEEE International Conference on Communication Systems. Singapore: IEEE, 2010. DOI: 10.1109/ICCS.2010.5686514.

[7] DAVIES J, DINALE J, MITCHELL P, et al. Towards open underwater acoustic channel models & tools for acoustic communication (ACOMMS) waveform-receiver performance assessment[C]//2024 Seventh Underwater Communications and Networking Conference. Sestri Levante: IEEE, 2024. DOI: 10.1109/UComms64662.2024.10847978.

[8] PIERSON W J, MOSKOWITZ L. A proposed spectral form for fully developed wind seas based on the similarity theory of S. A. Kitaigorodskii[J]. Journal of Geophysical Research, 1964, 69(24): 5181-5190. DOI: 10.1029/JZ069i024p05181.

[9] KACHOYAN B J, MACASKILL C. Acoustic scattering from an arbitrarily rough surface[J]. The Journal of the Acoustical Society of America, 1987, 82(5): 1720-1726. DOI: 10.1121/1.395165.

[10] CEBALLOS CARRASCOSA P, STOJANOVIC M. Adaptive channel estimation and data detection for underwater acoustic MIMO-OFDM systems[J]. IEEE Journal of Oceanic Engineering, 2010, 35(3): 635-646. DOI: 10.1109/JOE.2010.2052326.

[11] FENG X, ZHOU M, WANG J, et al. Model-driven deep learning-based estimation for underwater acoustic channels with uncertain sparsity[J]. IEEE Transactions on Wireless Communications, 2024, 23(6): 5710-5725. DOI: 10.1109/TWC.2023.3327995.

[12] OTNES R, EGGEN T. Underwater acoustic communications: Long-term test of turbo equalization in shallow water[J]. IEEE Journal of Oceanic Engineering, 2008, 33(3): 321-334. DOI: 10.1109/JOE.2008.925893.

[13] WU Z, HAN X, WEI L, et al. Joint temporal correlated channel estimation and data detection for underwater acoustic communications[C]//2025 IEEE 15th International Conference on Signal Processing, Communications and Computing. Hong Kong: IEEE, 2025: 1-5. DOI: 10.1109/ICSPCC66825.2025.11194391.

[14] LI B, ZHOU S, STOJANOVIC M, et al. Multicarrier communication over underwater acoustic channels with nonuniform Doppler shifts[J]. IEEE Journal of Oceanic Engineering, 2008, 33(2): 198-209. DOI: 10.1109/JOE.2008.920471.

[15] YANG Y, MA L, LIU S, et al. Underwater acoustic OTFS with nonuniform Doppler shifts: Modeling and off-grid block-sparse channel estimation algorithm[J]. IEEE Journal of Oceanic Engineering, 2026, 51(1): 641-658. DOI: 10.1109/JOE.2025.3599256.

[16] 杨健敏, 王佳惠, 乔钢, 等. 水声通信及网络技术综述[J]. 电子与信息学报, 2024, 46(1): 1-21. DOI: 10.11999/JEIT230424.

[17] 魏昕, 赵力, 邹采荣. 浅海水声通信中的间接自适应均衡算法[J]. 电子与信息学报, 2009, 31(4): 916-919. DOI: 10.3724/SP.J.1146.2007.01979.

[18] 唐怀东, 朱敏, 武岩波. 一种水声通信Turbo均衡中的软迭代信道估计算法[J]. 电子与信息学报, 2013, 35(3): 677-682. DOI: 10.3724/SP.J.1146.2012.00820.

[19] 杨斌斌, 鄢社锋, 章绍晨, 等. 基于Kalman滤波的水声混合双向迭代信道均衡算法[J]. 电子与信息学报, 2022, 44(6): 1879-1886. DOI: 10.11999/JEIT211343.

[20] 杨光, 丁寒雪, 郭庆华, 等. 基于叠加训练序列和低复杂度频域Turbo均衡的时变水声信道估计和均衡[J]. 电子与信息学报, 2021, 43(3): 850-856. DOI: 10.11999/JEIT200315.
