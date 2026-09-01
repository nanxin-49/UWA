# Bellhop 2020 旋转坐标内部粗糙墙技术可行性审查

## 1. 审查结论

**结论：`FEASIBLE_WITH_LIMITS`。**

在以下严格限定下，可以通过较小的、validation-only 的 Bellhop 2020 源码修改完成 proof-of-concept，而不需要重写大部分 Gaussian beam influence，也不需要引入任何 PE/Kirchhoff/SSA 相位屏或附加散射公式：

- 仅使用 2D Bellhop；
- 均匀声速 `c = 1500 m/s`；
- 单条参数化内部墙，只允许一次反射；
- 墙面为 pressure-release/vacuum；
- 不允许海底反射和其他边界反射；
- 反射后使用固定的 orientation-preserving 刚性旋转；
- 场计算先限定到当前验证实际使用的 Cartesian geometric-hat influence；
- influence 只接收坐标变换后的反射分支，不能跨越物理坐标与变换坐标之间的接缝。

在此范围内，物理反射可以直接复用 2D `bellhop.f90` 内部的现有 `Reflect2D`。需要新增的是墙面几何、真实交点事件、一次反射状态和反射后坐标分支管理，而不是新的反射物理。

以下扩展不在本结论覆盖范围内：layered SSP、3D Bellhop、多次墙面反射、动态海面、任意 beam influence 类型、跨接缝的 arrivals/eigenray 搜索。若要求同时支持这些项目，就不能再视为“小型 validation-only 修改”。

## 2. 审查对象和一个容易误改的位置

审查对象是：

`E:/MISC/BELLHOP/AcousticsToolbox_2020/Bellhop`

Bellhop 版本标记为 `2020_11_4`。

源码中存在两个名为 `Reflect2D` 的过程：

1. `bellhop.f90` 约 598–770 行的内部过程，是标准 2D `bellhop.exe` 实际调用的版本；
2. `ReflectMod.f90` 约 1–174 行的模块过程，主要由 Nx2D/3D 路径使用。

本方案必须修改或复用第一个。只修改 `ReflectMod.f90` 不会改变标准 2D Bellhop 的反射行为。

## 3. 当前 2D 调用链

### 3.1 环境与边界读取

- `bellhop.f90:130–135`：读取 `.ati` top boundary 和 `.bty` bottom boundary；
- `ReadEnvironmentBell.f90:154–213`：解释 beam/run type；
- `ReadEnvironmentBell.f90:442–458`：处理 vacuum boundary；
- `bdryMod.f90:28–34`：`BdryPt` 保存边界点、切向、法向、节点切/法向、长度、曲率等状态。

### 3.2 边界切向、法向和曲率

`bdryMod.f90:228–328` 的 `ComputeBdryTangentNormal` 完成现有 top/bottom 几何预处理：

- 约 271 行：由边界图 `z(r)` 计算斜率；
- 274–284 行：归一化切向，并按 top/bottom 选择外法向；
- 289–307 行：构造 curvilinear interpolation 使用的节点切向和法向；
- 309–322 行：由切向角变化计算有符号曲率。

`GetTopSeg`/`GetBotSeg` 位于 `bdryMod.f90:332–376`，通过边界 range 坐标索引分段。它们假设边界是按 range 单调排列的图 `z(r)`，因此不能不加修改地用于旋转后的近竖直墙 `r_w(z)`。

### 3.3 步进、crossing detection 与真实交点

`Step.f90:8–103` 的 `Step2D` 先做 Euler 预测，再用中点状态完成一步。

`Step.f90:107–176` 的 `ReduceStep2D`：

- 约 122 行构造试探终点；
- 124–132 行限制 SSP depth crossing；
- 134–148 行以边界局部直线平面计算 top/bottom 交点并缩短步长；
- 150–166 行限制边界和 SSP 分段 crossing；
- 168 行选择最小允许步长。

`bellhop.f90:512–554` 在完成步进后检查边界 signed distance 是否由正变为非正；curvilinear 模式下在 523–530 和 540–547 行插值局部切向/法向，然后在 532 或 549 行调用 `Reflect2D`。

因此内部墙不能仅在步后检测：必须参加 `ReduceStep2D` 的步长限制，确保状态真正落在墙面交点，而不是在穿墙后的点执行反射。

### 3.4 2D `Reflect2D`

`bellhop.f90:598–770` 的内部 `Reflect2D` 完成：

- 615–616 行：生成相同坐标的 incident/reflected 双节点；
- 618–619 行：计算入射切向和法向分量 `Tg`、`Th`；
- 621–624 行：按

  \[
  \mathbf t_+ = \mathbf t_- - 2(\mathbf t_-\cdot\mathbf n_b)\mathbf n_b
  \]

  计算镜面反射方向；
- 631–637 行：构造入射和反射两侧的 ray-normal 标架；
- 639–651 行：组合边界曲率和声速梯度跳变，形成动态射线修正 `RN`、`RM`；
- 645–648 行：当 `BotTop == 'TOP'` 时翻转部分符号；
- 660–663 行：更新

  \[
  p_+ = p_- + q_- RN,\qquad q_+ = q_-;
  \]

- 671–673 行：vacuum boundary 保持幅度并给 `Phase` 增加 \(\pi\)。

这正是目标方案需要保留的 native Bellhop 反射行为。

### 3.5 receiver influence 与相干场累加

`bellhop.f90:280–306` 在每条 ray trace 后选择 influence。当前验证文件写出的 run type 为 `C *`，第二个字符缺省为 `G`，实际进入 `InfluenceGeoHatCart`。

这里必须区分两个常被统称为“Gaussian”的概念：当前 `.sbp` 实现的是 Gaussian **发射角权重/源指向性**；当前 coherent-field accumulator 使用的是 geometric-hat Cartesian beam influence，并不是 Cerveny Gaussian influence。POC 应保留同一 `.sbp`，同时锁定当前 `G` influence。若改用其他 Bellhop Gaussian/Cerveny beam type，需要另做 range、标架和相位回归，不能由本报告直接视为已通过。

`influence.f90:397–500` 的 `InfluenceGeoHatCart`：

- 407–410 行：建立 `q0` 和初始相位/动态射线状态；
- 415–417 行：从 receiver range 建立初始索引；
- 421–434 行：计算 ray segment 的切向、法向并插值 `q` 和 travel time；
- 436–438 行：根据 `q` 的符号变化加入 caustic phase；
- 440 行：计算 Cartesian beam footprint；
- 451–495 行：在 receiver ranges 上施加 beam influence；
- 625–644 行的 `ApplyContribution`：分派到 arrivals、eigenrays 或 coherent field 累加。

`InfluenceGeoHatCart` 对单段左右传播有一定兼容性，但其索引初始化和分支翻转并不是为“同一数组中途反向”设计的。其他 influence 的正 range 假设更明确，例如：

- `InfluenceCervenyRayCen` 约 95 行跳过 `ir1 >= ir2`；
- `InfluenceCervenyCart` 约 215、225 行包含正向 range 假设；
- SGB influence 约 648–717 行按正向 receiver range 循环。

所以不能把“完整入射分支 + 反射返回分支 + 坐标跳变”直接交给现有 accumulator。

### 3.6 ray 状态与输出

`bellhopMod.f90:32–36` 的 `ray2DPt` 保存：

- `x(2)`、`t(2)`；
- 两组标量动态射线解 `p(2)`、`q(2)`；
- `c`、`Amp`、`Phase`、`tau`；
- top/bottom bounce counts。

标准 2D 状态中没有 `Rfa` 字段，因此本方案不需要处理 `Rfa` 的坐标变换。

`WriteRay.f90:17–47` 会压缩 ray 点并按照物理 top/bottom 保存输出。它不知道内部墙，也不知道两个坐标图之间的接缝。若输出跨接缝的整条轨迹，默认压缩和绘图会产生一条非物理连接线。

## 4. internal rough wall 的几何与符号约定

设墙为按深度参数化的曲线

\[
\mathbf x_w(z)=\begin{bmatrix}r_w(z)\\z\end{bmatrix},
\]

墙左侧 \(r<r_w(z)\) 为水体内部，墙右侧为外部。若墙点按增加的 \(z\) 排列，定义

\[
g=r_w'(z),\qquad a=\sqrt{1+g^2},
\]

则推荐使用

\[
\mathbf t_w=\frac{1}{a}\begin{bmatrix}g\\1\end{bmatrix},\qquad
\mathbf n_w=\frac{1}{a}\begin{bmatrix}1\\-g\end{bmatrix}.
\]

这里 \(\mathbf n_w\) 指向墙外，且在二维旋转矩阵

\[
\mathbf J=\begin{bmatrix}0&-1\\1&0\end{bmatrix}
\]

下满足

\[
\mathbf n_w=-\mathbf J\mathbf t_w.
\]

这与 Bellhop top boundary 的标架方向一致。因此调用现有 `Reflect2D` 时应把内部墙当作 `TOP` 语义，而不是 `BOT`。

沿上述切向方向的有符号曲率为

\[
\kappa_w=-\frac{r_w''(z)}{\left[1+r_w'(z)^2\right]^{3/2}}.
\]

对于

\[
r_w(z)=R_0-\eta(z),
\]

有

\[
\kappa_w=\frac{\eta''(z)}{\left[1+\eta'(z)^2\right]^{3/2}}.
\]

### 4.1 最稳妥的几何复用方式

不建议为 wall 独立重写斜率和曲率公式。推荐将相同的原始粗糙剖面先表示为 Bellhop native top 曲线

\[
\mathbf x_{top}(s)=\begin{bmatrix}s\\\eta(s)\end{bmatrix},
\]

让现有 top-boundary 几何预处理产生 \(\mathbf t_{top}\)、\(\mathbf n_{top}\) 和 \(\kappa_{top}\)，再通过固定的 90° proper rotation

\[
\mathbf Q=\begin{bmatrix}0&-1\\1&0\end{bmatrix}
\]

构造墙：

\[
\begin{aligned}
\mathbf x_w &= \begin{bmatrix}R_0\\0\end{bmatrix}+\mathbf Q\mathbf x_{top},\\
\mathbf t_w &= \mathbf Q\mathbf t_{top},\\
\mathbf n_w &= \mathbf Q\mathbf n_{top},\\
\kappa_w &= \kappa_{top}.
\end{aligned}
\]

因为 \(\det\mathbf Q=+1\)，切向、外法向和有符号曲率的 handedness 全部保留。这是使 native top 与 rotated wall 严格可比的关键。

为此可把 `ComputeBdryTangentNormal` 拆成：

- 保持现有 top/bottom 接口不变的 wrapper；
- 一个显式接收点集、边界类型和插值类型的 geometry core。

validation wall 调用同一个 core，再旋转输出。这样比复制代码更容易确保 Bellhop native 与 rotated-wall 使用完全相同的插值和曲率定义。

### 4.2 是否可以直接调用现有 `Reflect2D`

**可以，但有前置条件。**

- 墙必须使用上述 `TOP` 标架；
- `t`、`n`、`kappa` 必须来自同一个有向参数化；
- boundary halfspace 必须设置成 Bellhop 原生 vacuum；
- 只调用一次 `Reflect2D`；
- 反射完成后才执行坐标变换；
- 不得在坐标变换后再次调用反射函数。

flat/tilted straight wall 的 \(\kappa=0\)，曲率符号风险不存在。curved wall 的主要风险是：调用者若按外法向重新定义了曲率，又把墙标为 `TOP`，就可能与 `Reflect2D:645–648` 的内部符号翻转叠加，造成 double sign error。

POC 应在反射点强制检查：

\[
\|\mathbf t\|\simeq1,\quad
\|\mathbf n\|\simeq1,\quad
\mathbf t\cdot\mathbf n\simeq0,\quad
\mathbf n\cdot(\mathbf J\mathbf t)\simeq-1,
\]

并确认入射侧 `Th` 的符号与 native top 完全一致。若需要将来通用化，最小的安全重构是把 `TOP/BOT` 字符串隐含的 frame sign 提取成显式参数

\[
\sigma=\mathbf n\cdot(\mathbf J\mathbf t)\in\{-1,+1\},
\]

而不是复制一套 wall reflection 公式。

curvilinear 插值时，现有源码对插值后的节点切向/法向没有再次归一化。为了 native/rotated 的严格协变比较，第一阶段应对两者实施完全相同的处理，不应只对 wall 单独归一化；同时应记录 norm defect，作为数值审计指标。

## 5. 墙面交点检测

现有 `GetTopSeg`/`GetBotSeg` 基于 range 索引，不能直接用于 `r_w(z)`。内部墙需要独立的 `GetWallSeg(z)` 或真正的二维线段求交。

对当前 ray 段和第 \(j\) 个墙段，可解

\[
\mathbf x_0+h\mathbf u
=\mathbf w_j+\lambda(\mathbf w_{j+1}-\mathbf w_j),
\]

从所有满足

\[
0<h\le h_{trial},\qquad 0\le\lambda\le1
\]

的候选中取最小正 \(h\)。求得的 \(h\) 必须反馈给 `ReduceStep2D`，使最终 ray point 位于真实交点。

由于 `Step2D` 会调用两次 step reduction，第一次 Euler 预测阶段不能立即把 wall 标记成“已反射”。一次反射状态只能在最终接受的步长和交点确定后锁存。还必须专门处理：

- 顶点命中，避免相邻两段重复触发；
- grazing incidence，设置 `|Th|` 下限；
- 墙端点，避免把有限墙段外推成非物理反射；
- step 恰好从墙面开始时，使用方向性和容差避免零步长再次命中。

## 6. 推荐的反射后固定坐标变换

### 6.1 推荐：绕名义墙中心旋转 π

取固定中心

\[
\mathbf x_0=\begin{bmatrix}R_0\\0\end{bmatrix},
\]

使用

\[
\boxed{
\mathbf x'=\mathbf M\mathbf x+\mathbf b,
\qquad
\mathbf M=-\mathbf I,
\qquad
\mathbf b=2\mathbf x_0
}
\]

即

\[
r'=2R_0-r,\qquad z'=-z.
\]

这是 determinant 为 +1 的 proper rigid rotation，而不是镜像。对于物理反射后满足 \(t_r<0\) 的 ray，

\[
\mathbf t'=\mathbf M\mathbf t=-\mathbf t
\]

使新的 range tangent 为正。POC 应显式拒绝未满足 \(t_r<-\epsilon\) 的高斜率或 grazing case，而不是假设所有墙形都能被该固定变换变成 range-monotone。

物理返回侧 receiver

\[
\mathbf x_R=\begin{bmatrix}R_0-d\\z_R\end{bmatrix}
\]

映射为固定位置

\[
\mathbf x'_R=\begin{bmatrix}R_0+d\\-z_R\end{bmatrix}.
\]

它只取决于物理 receiver 和固定 \(R_0\)，不依赖具体 ray 的墙面交点。因此不会产生 ray-dependent receiver。

对于 flat case，\(R_0=100\,\mathrm m\)、\(d=3\,\mathrm m\)，receiver 的新 range 正好为 103 m。

### 6.2 长度、travel time 和反射系数

因为

\[
\mathbf M^T\mathbf M=\mathbf I,
\]

线段长度保持不变，累计 `tau` 原样保留。该变换在 `Reflect2D` 完成之后执行，因此：

- 不重新计算 reflection coefficient；
- 不再增加 vacuum \(\pi\) 相位；
- `Amp` 和 `Phase` 原样保留；
- 只变换几何位置和方向。

### 6.3 ray tangent、normal、p/q 与 caustic phase

对 proper rotation：

\[
\mathbf x'=\mathbf M\mathbf x+\mathbf b,
\quad
\mathbf t'=\mathbf M\mathbf t,
\quad
\mathbf n_{ray}'=\mathbf M\mathbf n_{ray}.
\]

由于二维 proper rotation 满足

\[
\mathbf J\mathbf M=\mathbf M\mathbf J,
\]

Bellhop 由切向生成的 ray-normal 标架与变换严格协变。

`ray2DPt%p(2)` 和 `%q(2)` 是相对于该随 ray 旋转的横向标架定义的标量动态射线解，不是全局坐标向量。因此 proper rotation 下：

- `p` 不变；
- `q` 不变；
- 不应人为乘以 -1；
- `q` 的过零结构不变；
- caustic/KMAH phase 不变；
- Gaussian beam amplitude 和 phase 不引入额外修正。

标准 2D `ray2DPt` 没有 `Rfa`，所以没有额外的 `Rfa` 变换。若未来使用不同 Bellhop 分支或 3D 状态，则必须重新审查该结论。

### 6.4 为什么不推荐普通镜像

看似自然的变换

\[
\mathbf M=\operatorname{diag}(-1,1),
\qquad
\mathbf b=\begin{bmatrix}2R_0\\0\end{bmatrix}
\]

会把返回 ray 变成正 range，并保持深度不变，但 \(\det\mathbf M=-1\)。此时

\[
\mathbf J\mathbf M=-\mathbf M\mathbf J,
\]

ray-normal 标架发生 orientation reversal。若不处理，`p/q` 的符号语义可能与由 tangent 重建的 normal 不一致；若简单翻转 `q`，又可能在当前 influence 中制造假的 caustic 相位变化。

因此普通镜像不是 POC 的推荐方案。π rotation 避免了这个 handedness 风险。

### 6.5 坐标接缝

π rotation 不会固定粗糙墙上的每一个交点，只固定旋转中心。因此物理入射分支与变换后的反射分支属于两个坐标图，在数组中会出现几何接缝。

该接缝不能被当作真实 ray segment：

- 不得给接缝增加距离或 travel time；
- 不得在接缝上做 beam influence；
- 不得让 ray compression 把它连接成可见物理线段。

推荐在 trace 中保存 `wallBranchStart`，反射完成并旋转之后，只把变换后的反射分支复制/打包给现有 influence。物理入射分支可保留在独立的 validation diagnostics 中。

## 7. receiver influence 能否继续使用

**在限定到 transformed post-wall branch 和当前 `InfluenceGeoHatCart` 时，可以不修改 field accumulator。**

条件是：

1. 传入 influence 的 ray 数组从墙面反射节点开始；
2. 数组中所有后续 segment 的 range 单调增加；
3. mapped receiver ranges 是固定的、升序排列的；
4. `z'=-z` 后若 receiver depth 数组顺序改变，要重排输入并在输出时恢复索引；
5. `tau`、`Amp`、`Phase`、`p/q` 均保留物理反射分支的累计值；
6. `SourceDeclAngle` 仍是原始发射角，`RcvrDeclAngle` 使用变换后的切向角；
7. 不允许 accumulator 跨越坐标接缝。

在均匀声速下，只给 influence 传 post-wall branch 仍能保持完整传播归一化：动态射线状态和 `tau` 已经从源累计到墙后，`q0=c/Dalpha` 仍来自原始 beam fan。`ScalePressure` 使用首 ray point 的声速在此也安全，因为全域 `c` 相同。这个论证不能直接推广到 layered/asymmetric SSP。

`InfluenceGeoHatCart` 会在传入 ray 的首点把内部 caustic phase 计数重新置零。当前 uniform-c 入射段是直线传播，墙前 `q` 不会发生 caustic crossing，所以从 reflected wall node 开始打包是安全的；墙后的 `q` crossing 仍会被正常检测。若未来允许墙前聚焦介质或其他 caustic，必须显式携带墙前 KMAH index，届时“完全不改 accumulator”的结论不再成立。

arrivals/eigenrays 最小 POC 可暂时把一次 internal pressure-release wall bounce 计入现有 `NumTopBnc=1`，因为没有真实 top/bottom 反射。这样不必修改 arrival 文件格式。报告和输出元数据必须明确：这个 top count 在 validation binary 中代表 internal wall。

## 8. 其他隐含的 range-increase 假设

### 8.1 ray tracing 和 termination

- `bellhop.f90:556–568` 使用 `abs(range)>RBOX` 和 `abs(depth)>ZBOX`，本身不要求正 range；
- 约 561 行原有“kill backward ray”逻辑已被注释，因此 tracer 可短暂处理负向 range；
- 但 validation 分支仍应在旋转后检查 range tangent 始终为正，并按 mapped beam box 终止。

### 8.2 边界与 SSP interpolation

- top/bottom 分段查找按 range 单调图工作，不能复用作竖直墙索引；
- SSP 查询按全局 depth 进行，`z'=-z` 后只有均匀或经过严格对称变换的 SSP 才保持物理等价；
- 本审查因此只批准 uniform SSP；
- artificial top/bottom 应放在足够宽的对称、声学匹配区域，避免 π rotation 后触发非物理边界。

### 8.3 receiver indexing 和 beam influence

- 多个 influence 例程假定 receiver ranges 有序且通常等间隔；
- 一些 Cerveny 和 SGB 路径明确跳过反向 range segment；
- 自动 beam 数量和 beam box 的设置会参考最大 receiver range；mapped receivers 必须在读取/初始化阶段就提供；
- 第一阶段只批准 `C *` 当前默认进入的 Cartesian geometric-hat 路径，其他 beam types 需要逐一做旋转协变回归。

### 8.4 arrivals/eigenrays

- arrivals/eigenrays 共用 influence 的 receiver crossing 与 bounce count；
- 只要输入是单调的 post-wall branch，可复用现有 crossing；
- 若要同时输出物理入射分支和 transformed 分支，必须使用独立的 validation diagnostics，不能把接缝交给 eigenray 插值。

### 8.5 ray compression/output

- `WriteRay2D` 只认识 top/bottom boundary；
- 默认 ray output 不能表达坐标图切换；
- 最小实现可只写 transformed post-wall branch，或在 validation module 中另写两个分支；
- 不建议为了 POC 修改正式 `.ray` 格式。

## 9. 最小源码修改路径

| 文件/函数 | 最小修改目的 | 是否改变反射物理 |
|---|---|---|
| 新建 validation-only `InternalWall2DMod.f90` | 读取独立 wall sidecar；保存原始剖面和旋转后几何；墙段索引/求交；一次反射状态；π rotation；post-wall 分支打包和诊断 | 否 |
| `bdryMod.f90` / `ComputeBdryTangentNormal` | 将几何计算拆成可由 top/bottom 和 wall 共用的 core，确保切向、法向、曲率和插值严格一致 | 否 |
| `Step.f90` / `Step2D`, `ReduceStep2D` | 增加可选 wall intersection step limit 和最终 hit 事件；保证 ray 落在真实交点 | 否 |
| `bellhop.f90` / `TraceRay2D` | 初始化 wall；只在未反射时检测；调用现有内部 `Reflect2D(...,'TOP',vacuum,...)`；随后做一次 π rotation；记录 `wallBranchStart`；把 post-wall branch 交给 influence | 复用现有物理 |
| `bellhop.f90` 主循环 | validation 模式下选择 packed post-wall ray；锁定/检查允许的 run/beam type | 否 |
| `WriteRay.f90`，可选 | 仅当必须用正式 R-mode 可视化时识别分支；更推荐 validation module 独立输出，避免改正式格式 | 否 |
| `Makefile` | 加入 wall module 和依赖，生成独立 validation binary | 否 |

不建议修改：

- `influence.f90`：限定为 packed、单调 post-wall branch 时无需修改；
- `bellhopMod.f90`：可用模块/local state 和 `NumTopBnc` 代理避免扩展公开 ray/arrival schema；
- PE marching、Gaussian 源定义、生产通信链：均不应触碰。

wall 输入不应伪装成 `.ati`。建议使用显式 sidecar，例如 `.iw2`，并只由独立命名的 validation binary（例如 `bellhop_iwall_2020.exe`）启用。正式 Bellhop 2020 executable 保持不变。

当前机器审查时未在 PATH 中发现 `gfortran`/`make`。这不改变源码层面的可行性结论，但实际 POC 开工前需要先确认与 Bellhop 2020 兼容的 Fortran 构建链；本次审查没有安装工具链。

## 10. proof-of-concept 验证顺序

### A. flat internal wall

设置

\[
r_w=R_0=100\ \mathrm m,
\]

使用与现有基线完全相同的 `.sbp`、beam fan、beam count、step 和声源设置。只累加 transformed post-wall branch。

必须检查：

- 精确交点 residual；
- 反射前后镜面方向；
- `Phase` 只增加一次 \(\pi\)；
- 反射瞬间 `Amp` 不变；
- \(\kappa=0\) 时 `p/q` 不受曲率 kick；
- 累计 travel time 为 \(103/1500\) s；
- mapped receiver 固定在 103 m；
- 相干场恢复

  \[
  H_{ref}=-P_{BH}(103\ \mathrm m).
  \]

均匀介质、平墙下，总几何路径是 100 m + 3 m。native `Reflect2D` 只加入 pressure-release \(\pi\) 相位，π rotation 不改变 `p/q/tau/Amp/Phase`，所以这个等式应当是强验收门，而不是仅检查“趋势相近”。还要对 step 和 beam count 做收敛扫描。

### B. tilted straight wall

使用斜直线，解析计算：

- ray-wall 交点；
- 单位切向和外法向；
- \(\mathbf t_+=\mathbf t_- -2(\mathbf t_-\cdot\mathbf n)\mathbf n\)；
- \(\kappa=0\)；
- π rotation 后的固定 receiver 位置和正 range tangent。

逐 ray 比较位置、方向、travel time、`Phase` 和 `p/q`，先不依赖相干场恰好掩盖几何错误。

### C. smooth sinusoidal wall

取

\[
r_w(z)=100-\eta(z),
\]

native Bellhop 的 `.ati` top 和 rotated-wall 必须使用同一组原始 \(\eta\) 样点、同一 interpolation type 和同一几何 core。

比较顺序：

1. 对应 ray 的交点；
2. tangent/normal/curvature；
3. reflected tangent；
4. reflection 前后 `p/q`；
5. travel time 和 boundary phase；
6. receiver coherent field。

这个阶段是曲率符号的主验收。应先用低振幅、长波长正弦，再提高曲率并做 profile sampling、step 和 beam-count 收敛。

### D. fixed 1D PM realization

只有 A–C 全部通过后，才加入固定 seed 的一维 PM realization。继续只比较

\[
Bellhop_{native}\longleftrightarrow Bellhop_{rotated-wall},
\]

不运行 PE。必须保存同一原始剖面及其插值设置，并检查：

- 最大 slope/curvature；
- grazing 和墙顶点命中数；
- transformed branch 是否始终 range-monotone；
- profile sampling、ray step、beam fan/count 的联合收敛。

## 11. 主要数学和数值风险

按优先级排序：

1. **曲率符号/标架 double flip。** `Reflect2D` 对 `TOP` 有内部符号处理，wall 几何必须由 native top 标架 proper-rotate 得到。
2. **步后反射而非真实交点。** 若 wall 不参与 `ReduceStep2D`，误差会进入方向、曲率 kick、travel time 和相位。
3. **坐标接缝被误当成传播段。** 会制造巨大假距离、假 travel time 和错误 receiver influence。
4. **improper mirror 反转 ray-normal handedness。** 可能污染 `p/q` 和 caustic phase；采用 π rotation 避免。
5. **固定旋转不能保证任意粗糙墙反射后都为正 range。** 需要 slope/incidence gate，超限 ray 应明确拒绝。
6. **curvilinear 插值后 t/n 非严格单位化。** native 与 rotated 必须同处理，并记录误差。
7. **PM 剖面高频曲率对采样敏感。** 必须在 smooth sinusoid 通过后再进入，并做采样收敛。
8. **现有 influence 类型并非全部支持反向或分支翻转。** POC 仅批准 `InfluenceGeoHatCart` 和 packed post-wall branch。
9. **ray/arrival 输出语义。** 默认文件格式没有 wall bounce 和坐标图切换概念，必须用 validation 元数据明确标注。
10. **负深度和人工边界。** π rotation 产生 `z'=-z`；只在 uniform c 和足够宽、匹配的人工上下边界中安全。

## 12. 停止条件

以下任一情况出现时，应停止扩展，不进入 D，也不修改 PE：

- A 不能在数值收敛后恢复 `-P_BH(103 m)`；
- C 中 native top 与 rotated wall 的 `t/n/kappa` 无法逐点协变；
- 必须修改 `InfluenceGeoHatCart` 的核心 Gaussian beam 公式才能让 flat case 成立；
- 必须对 `p/q` 人为重置、重新拟合或逐 ray 校正相位；
- fixed π rotation 后大量目标 ray 仍非 range-monotone；
- 为表达分支不得不重写 arrivals/eigenray/coherent accumulator 的大部分逻辑。

在这些停止条件之外，当前源码结构支持一个小而隔离的 validation fork，且不要求改变 Bellhop 的核心物理。
