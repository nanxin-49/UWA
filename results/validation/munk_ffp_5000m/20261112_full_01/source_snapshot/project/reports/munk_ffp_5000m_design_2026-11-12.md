# 5000 m Munk 独立波动参考：复用审查与事前设计

**当前状态：J01–J07、最终数值/图件审计及交付均已完成，给定模型数值门通过；停止待人工审查。见结果报告与最终checkpoint。**

历史事前/申请状态（保留）：状态：实现与小型资格先导；尚无完整场结果。独立工作树`C:/Users/ASUS/.codex/worktrees/munk_ffp_20261112/Explain`，分支`experiment/munk-ffp-5000m`，固定提交`6c536328b483daa21ced98e01872bb7387ddf03c`。主检出并发任务在审查期间从a37fd75推进至6c53632，未修改/暂存/转移其未提交内容；旧Bellhop实验工作树不修改。

## 已有实现的逐项判断

| 实现 | 代码中实际含义 | 复用方式/限制 |
|---|---|---|
| deepsea3d_point_field | exp(ikR)/(4πR)点源Green，exp(-iωt) | 直接调用作独立均匀基准及已知奇点减除，不重写 |
| deepsea3d_flat | 压力释放平面镜像，含明确路径/相位元数据 | 直接调用作均匀边界基准，不重写 |
| weyl_point_source_reference_vertical | 连续Weyl积分，传播/倏逝换元 | 小型检查复用；现有逐r自适应积分不当作5000m全场批量求解器 |
| deepsea_bie_gauss | Gauss-Legendre节点/权重 | 直接复用，不重写求积节点 |
| solve_full_kirchhoff_2d、solve_deepsea3d_fullfk | 均匀Green核的Kirchhoff边界近似 | 不支持Munk分层垂向Green，不能称FFP |
| deepsea3d_fe_recursion | FE=Field Expansion场展开递归，不是有限元 | 仅复用约定/已有基准，不冒称三维分层有限元 |
| deepsea3d_tfe_modal_plan/solve/operators | 均匀k、横向FFT和垂向Chebyshev的变换场展开 | 未支持5000m Munk、点源jump及流体底；直接替换k²不足以获得资格 |
| AcousticsToolbox_2020/Scooter/scooter.f90 | 真正分层有限元/波数Green求解，含密度/流体半空间 | 保留为候选交叉检查。实源码采用5/12与1/12权重，不能套普通二阶FE色散估计；仍需独立垂向网格收敛/源点插值检查。未启动外部SCOOTER |
| 同工具包TransformG.f90/fieldsco.m | 默认远距Hankel渐近式，fieldsco.m将r=0改1m，H选项个别Hankel项各自奇异 | 不用于零轴重建。使用精确J0逆变换及完整轮廓连接段 |

本机scooter.exe SHA256=`B7AD8C47F65E1EE3B4FAF304BA3C2D0285999CA801CA0EF2663382DC42296483`，未替换版本。主检出及a28b的现行脚本未定位到已保存的Munk垂直时延入口；用户在本轮已明确允许找不到时补写。新积分作为独立验证，不当作波动场结果。

## 数学模型与新实现边界

物理目标是三维点源在轴对称分层介质中的二维切面，不是二维线源。采用exp(-iωt)，g满足`g''+[k(z)^2-kappa^2]g=-delta(z-zs)`（水体密度常数）。精确关系为`G(r,z)=(1/(2π))*integral_0^Inf g(kappa,z)*J0(kappa*r)*kappa dkappa`。J0(0)=1，除真实点源外直接计算零轴。

本轮资格用新验证辅助函数计算分片常速层的精确上下行波，多次反射通过有界反射/透射递推，避免用5000m长域的低阶波相位差代替精确传播。连续Munk通过层厚10/5/2.5m收敛；这是新分层Green实现，复用已有点源/镜像/Weyl/求积，不重复实现这些已完成基准。源点jump、表面Dirichlet和流体底admittance独立验证。未接入src或生产API。

从有限矩形轮廓`0 -> -i*a -> K-i*a -> K`积分，以避开实轴振荡和近极点；所有连接段必须计算，不能只平移横段。底介质被动吸收、出射分支Im(q)>=0，轮廓位于物理奇点下侧。用已有解析点源减除均匀奇异核，积分差值后加回解析Green；在源点用NaN+singularity mask，不人为给有限峰值。a独立变化必须收敛，不作为海水吸收。K增加与面板宽度减半分别测试。

主要保存`G_green`与`pressure_1m=4πG_green`两套复数；TL用后者，与历史1/R点源参考同尺度，不拟合增益、不声称Pa。海底cp1600m/s、rho1.8g/cm³、衰减.8dB/波长是示例假设，海水吸收0；不写为论文原参数。Munk参数仍为c0=1500、axis/scale1300m、epsilon=.00737。

## 资格门与全场交付

- 均匀kernel与解析谱核：复相对误差≤1e-10；均匀空间Green/镜像复L2≤1e-3、TL≤.05dB、相位≤.01rad（远离消零点，另报全支持绝对误差）。
- Munk直达垂向travel integral与分层one-pass累计传播时延：时延误差≤.05/(2π*10000)s；这是直达相位的离散检查，不把含多途的单频场相位误认为唯一时延。
- 非源轴压有限；表面相对边界残差≤1e-10；复场左右对称误差≤1e-12，注明J0偶性是算子性质，不能独自证明物理准确。
- K截断、求积面板宽度、轮廓偏移、Munk层厚必须分别测试。全支持复L2≤1e-3；近消零点不以无限大相对误差掩盖全域绝对误差，记录共同支持与分母。单个门失败不发完整场PASS。
- 全局x=-5:.05:5m、z=0:5:5000m为概览采样；全球5m不用于解析细条纹或错误相位解缠。局部z=4900:.01:5000m，满足lambda/10，源点只排除(0,4950)。保存完整复场和配置。
- 五组PNG/FIG及离线HTML：全场、轴幅相、局部高分辨率、均匀参考、收敛。全场和局部每个接收点都由同一Green积分得到；不镜像填场、不从Bellhop插值。

## 当前执行门

初次MCP附着失败，用户恢复后确认R2025b Update6/PID35640；当前仅1.122GiB系统可用RAM。只做低内存串行小先导，不建parpool，不启第二个MATLAB或SSH。每作业小块处理，剩余RAM动态监测；前期控制/计时预算≤180s，满域生产尚未授权额外预算。

正式多配置全场预计可能超过180s；必须先用先导计时给出可审查预算，再等待明确资源授权。不能将代码可运行、文件存在、轴上非零或图像相似写成完整资格。保留全部失败/未运行门，目标仍为完整5000m波动参考，不缩为均匀基准。

## 原始依据

- 官方SCOOTER描述：https://oalib-acoustics.org/website_resources/AcousticsToolbox/manual/node64.html （本机源码优先）。
- 官方FFP说明：https://oalib-acoustics.org/website_resources/AcousticsToolbox/manual/node63.html 。
- 局部方法审查对照本机Scooter/scooter.f90、TransformG.f90及Matlab/Scooter/fieldsco.m；哈希与引用行保存在本轮输入清单，不修改工具包。
