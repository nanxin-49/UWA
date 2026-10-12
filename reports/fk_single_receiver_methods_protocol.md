# 固定单接收点 Kirchhoff：方法与预先规定的试验协议

## 同一离散物理算子

采用现有 `solve_full_kirchhoff_2d` 原生 `exp(-iωt)` 约定：

\[
u_R=-\sum_j q_j\{G_j(n_xu_x+n_zu_z)_j+u_j(n_xG_x+n_zG_z)_j\}.
\]

海面 `(x,eta)`，水侧单位法向 `(-eta',1)/sqrt(1+eta'^2)`；`q`由相邻真实海面点弧长梯形权重构造。完整Hankel Green与源坐标导数保留，不裁剪空间孔径、不截谱、不拟合幅相。

阶段一复用原FastFK的加权SVD：`C(eta,kx)=exp(i(k-kz)eta)`，权重为源谱能量；保留统一rank3、33个Chebyshev高度节点，保留精确carrier `exp(-ik eta)`。原源周期IFFT只用于预计算入射模式及**物理**x/z偏导。海面→标量Rx调用未改动FullFK，完全没有出射低秩或接收线FFT。plan不含当前海面高度求值，可跨realization复用。

阶段二令 `P=G*u`，由乘法法则 `P_x=G*u_x+G_x*u`、`P_z=G*u_z+G_z*u`，积分为 `-sum(q*(nx*P_x+nz*P_z))`。两项同号使该合并严格成立；这里不是把法向或弧长冻结成平面。每个固定x处用精确谱IFFT生成高度节点上的入射场及导数，与精确Hankel Green合成偏导表。表去除高度快相位后存储 `exp(2ikz)*P_x,z`；apply行内重心插值并乘回 `exp(-2ik eta)`。载波移除不引入近似。直接插值物理偏导表，不误将残差表的导数当作物理偏导。无空间插值，无FFT出射卷积，无需预设额外低秩。

2026-10-12契约补强：两类plan绑定实际使用的源字段（kx/kz/coeff/period/n_original/x0/Tx）、频率、声速、网格、高度界和节点数；FastIn还绑定rank，**不绑定Rx**；复合核绑定Rx，**不绑定无关rank**。界内海面/法向和未使用的源metadata变化允许复用。域越界、旧依赖、错误kind和缺失plan必须匹配各自错误标识符。仅支持现有传播谱、周期对齐均匀网格和单个Rx；不是生产API。plan复用正确不代表新Rx误差也通过。

## 配对精度与数值护栏（运行前固定）

- 4/6/8 kHz、c1500、Tx100、Rx(0,3)，既有sigma.3的984点Gaussian谱；海面[-60,60]4097点；界±.5 m、rank3、节点33。它是研究fixture，不替换Tx99安装案例。
- canonical260001与既有Sk生成的11个独立Gaussian Fourier系数海面260002–260012；同12模式PM频带、同C2 plateau42/support50包络，无高度重标定。所有样本对所有方法配对。保存系数而非只保存seed。
- flat0/±.5、A=.02/.45且K=.47正弦控制；每种方法对原始FullFK相对复误差门槛1e-3，另报绝对误差、绝对/相对幅度误差和原生相位差。所有定义有效的误差必须有限。接近零场保留既有1e-14分母保护，同时新增raw相对误差；零声压的raw比例及相位明确undefined，不以0相位宣称准确。
- 精确乘积偏导与原FullFK闭合1e-12；canonical及强正弦33→65节点差1e-5、rank3→5差1e-3、4097→8193积分差1e-3，分别归入插值、压缩、积分离散控制。flat±.5验证节点精确命中。
- SHA256校验既有输入reference；执行前后校验FullFK/BIE与原FastFK基准。两类plan逐项验证height、源coeff/kx/dispersion、频率、声速、源周期/归一化/原点、Tx、网格、高度界、节点、rank、Rx、metadata、surface、kind与missing；允许的复用和重建按1e-12闭合。注入无关异常及无异常，确认拒绝匹配器会失败。

## 相消诊断与有条件冻结

保持原Rx(0,3)主域的51组输入与全部门槛。仅为诊断，在既有强正弦A=.45/K=.47海面，用FullFK扫描Rx x∈[-60,60]、z3并细化|x|≤50的局部低压位置，不改源、不拟合幅相、不纳入加速计时。记录绝对/相对误差、峰值归一化压力、相消比`|sum(I)|/sum(|I|)`；低于峰值1e-3（或绝对1e-14）的near_zero仅是解释标签。补测相同位置rank5、65高度节点及8193积分点，区分压缩、插值及积分敏感性，不替换原rank3/33节点。

相消点的1e-3失败必须单独保存，不用主域PASS覆盖。冻结资格仅指原Rx0研究fixture；不表示任意Rx、相消点或物理模型已通过。CSV/MAT/代码哈希与不变门槛保存于`freeze_manifest.json`；补强前证据保留在本case的`pre_hardening/`。审核与所有未关闭问题见[冻结审核报告](fk_single_receiver_baseline_freeze_report.md)。

## 公平计时与判据

五个配对块，轮换三种方法顺序；每块apply六个不同PM海面并取每海面均摊时间。相同海面生成和法向更新计入重复任务，两种plan均真正复用。保存每块原始值及顺序；重复计时是噪声控制，不计为新独立环境或PM人口可靠性证据。

冷启动：热JIT会话中新建plan加一次apply，FullFK包括重新精确生成当前海面入射场和标量积分。预计算与重复成本分别报告；文件IO、报告、护栏不进入计算计时。交叉点仅由实测中位数估算 `ceil(T_setup/(T_Full-T_apply))`，不是统计置信上界。不清理共享会话、不启动额外MATLAB。

推荐必须同时满足精度及实际收益。比较单次任务和重复realization两种范围；如果当前海面精确入射已由其他流程免费提供，FullFK只需标量Green积分，不能套用包括入射构建的加速倍数。内存报告plan `whos`大小，不声称峰值。没有收益的方案如实保留负面证据，不追加无目的优化。

生产PE、公共接口和原FullFK/BIE/FastFK不变。只运行独立验证工作流，不接入信道/通信链，因此不运行无关生产标量/宽带回归。算法一致性不等于Kirchhoff物理真实性；旧canonical BIE证据不覆盖本轮新随机海面，不能推广到3-D、深海或BER。
