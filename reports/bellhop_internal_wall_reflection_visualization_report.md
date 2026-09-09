# Bellhop tilted / sinusoidal internal-wall 反射可视化

## 结论

本轮没有重新运行 MATLAB、Bellhop 或 PE，也没有修改任何传播、反射或
receiver influence 公式。图件由已经通过的 validation 输出直接后处理得到，目的
是把此前只以表格和诊断量保存的反射过程画清楚。

采用的统一表达分成三个阶段：

```text
Tx -> 真实 wall hit -> 原生 Reflect2D 后的物理反射分支
                         |
                         +-> fixed proper pi rotation
                             -> 独立的正 range 计算分支 -> mapped Rx
```

物理反射图与旋转后计算图刻意放在不同 panel，避免把坐标重映射误画成第二次
物理反射。局部 frame panel 使用相同坐标比例，直接显示 incident/reflected
direction、wall tangent 和 TOP normal。

## 数据来源与代表设置

- tilted：`step=0.05 m`、`5001 beams`，读取
  `results/validation/bellhop_internal_tilted_wall_poc/cases/wall/step_0p05_beams_5001_wall.iwdiag`；墙为
  `r=100+0.005z`。
- sinusoidal：`step=0.05 m`、`5001 beams`、`161 profile samples`，分别读取
  `A=0.25 m` 与 `A=0.50 m` 的既有 `.iwdiag` 和 `.iw3`。墙线直接使用 `.iw3`
  中的原始 sampled profile，没有按公式另行生成。
- 收敛图直接读取既有 `tilted_wall_convergence.csv` 和
  `sinusoidal_wall_convergence.csv`。sinusoidal 图只采用 geometry/frame 与
  native curvature-kick 数据，没有复用此前已确认错误的 receiver-column 场配对。

全局 ray fan 取最接近 `-10/-5/0/5/10 deg` 的已存 ray。紫色 ray 是最接近目标
receiver 的离散 beam，不是重新求解的连续 eigenray；它在 tilted、弱 sinusoidal、
较强 sinusoidal 三个 case 的 receiver-plane depth miss 分别为 `8.43 mm`、
`-6.57 mm`、`8.43 mm`。sinusoidal 局部图另选 `|alpha|<=10 deg` 内
`|kappa|` 最大的已存 ray，以便显示非零曲率对 `p` 的原生 kick。

## 图件

### 一张图内的完整路径叠加

![End-to-end internal-wall overlay](../results/visualization/bellhop_internal_wall_reflection/04_internal_wall_end_to_end_overlay.png)

这是面向直观理解的主图。每个 case 都在同一坐标面板中显示从 Tx 到墙面、原始
反射以及 mapped receiver 的全过程：蓝色实线为入射段，灰色虚线为
`Reflect2D` 输出的原始 backward-range 物理反射分支，绿色实线为 fixed proper
`pi` rotation 后送往 `r'=103 m` receiver 的正-range 分支。紫色点线只连接原始
hit 和其 rotated image，明确表示一次瞬时坐标映射，不代表声线传播、额外路程或
第二次反射。细线 fan 提供周围 ray 的方向背景，加粗线为最接近目标 receiver 的
已存 beam。每个 panel 内的 `wall-to-receiver zoom` 保留同一数据坐标，专门放大
约 `3 m` 的短分支，使 `r=97 m` 灰色原始返回段与 `r'=103 m` 绿色接收段可以直接
辨认，同时仍保留上方约 `100 m` 入射全过程。

### Tilted straight wall

![Tilted internal-wall reflection process](../results/visualization/bellhop_internal_wall_reflection/01_tilted_wall_reflection_process.png)

左图是 Tx 到真实斜墙交点再到 native backward-range 的物理反射；中图是目标
beam 在 `Reflect2D` 处的局部 frame；右图只显示 proper rotation 后隔离出来的
正-range 分支。存储诊断给出 `kappa=0`、phase change 为一次 `pi`、`Amp` 不跳变，
rotation 不改变 `p/q`。

### Smooth sinusoidal walls

![Sinusoidal internal-wall reflection process](../results/visualization/bellhop_internal_wall_reflection/02_sinusoidal_wall_reflection_process.png)

两行分别是弱曲率 `A=0.25 m` 和稍强曲率 `A=0.50 m`。橙红色加粗 ray 是用于
局部 curvature 展示的代表 ray；在相同 hit 区域，`|kappa|` 与原生
`|p_after-p_before|` 随振幅增加，`q`、`Amp` 和一次 pressure-release phase 仍按
已有诊断保持预期行为。右列显示 proper rotation 后各分支均沿正 range 到达
固定 `r'=103 m` receiver plane。

### 已存收敛诊断

![Internal-wall convergence diagnostics](../results/visualization/bellhop_internal_wall_reflection/03_internal_wall_convergence_diagnostics.png)

tilted wall-hit residual 保持在浮点误差量级，mapped branch 的最小 range increment
始终为正；sinusoidal tangent error 随 profile samples 加密下降，而非零 curvature
kick 对同一振幅趋于稳定。该图是已有验证结果的可视化，不是新验收。

## 可复现入口与限制

只读后处理入口为：

```powershell
D:\Anaconda\python.exe scripts/reporting/generate_bellhop_internal_wall_reflection_visuals.py
```

输出目录为 `results/visualization/bellhop_internal_wall_reflection/`，其中
`visualization_manifest.csv` 记录每张图的输入文件、代表 ray 行号及关键数值。

由于环境是 uniform `c=1500 m/s`，ray segment 可由已存 hit 与 incident/reflected
direction 线性重建；本轮没有读取新的 `.ray` 轨迹。全局 panels 为同时容纳约
100 m incident path 与 3 m mapped branch 而使用各自坐标范围，定量判断方向应看
单位比例的局部 frame panel，而不是用屏幕像素比较不同 panels 的倾角。
