# M1：100 m 静态信道数值收敛与吸收层资格执行收尾

日期：2026-11-11（Asia/Shanghai）。最终状态：**M1_INCOMPLETE**。

本轮获准执行全部 M1，但现有 MATLAB MCP 首次连接检查返回 `failed to attach to MATLAB session`。未获得实际会话版本、真实计算 PID、日志与长作业状态追踪能力，触发用户指定的执行前置停止门。**5 次无吸收层 PE、10 组 AS 和 3 次吸收层 PE 均未派发；没有获得新的数值合格配置，不建议进入 M2。** 本报告是执行受阻的收尾记录，不是数值通过报告。

结果目录：[20261111_m1_execution_preflight_01](../results/validation/m1_convergence/20261111_m1_execution_preflight_01/)。本轮没有启动/终止 MATLAB 进程，没有使用 SSH、GPU、外部参考求解器或其他物理工况；未合并 main、未推送。

## 1. 基线、封存及检查范围

- 初始分支 `compare`，HEAD=`c90d0ca`，工作区干净；完整 HEAD 见 [checkpoint_preflight.json](../results/validation/m1_convergence/20261111_m1_execution_preflight_01/checkpoint_preflight.json) 中的 `git.head`。
- M0 生产提交 `969127fb2861a5d3b423f63bd7c3c34232284d71`；保留 `M0_SEALED`。已阅读 AGENTS 当前规则、PROJECT_CONTEXT 当前索引、M0 收尾报告、M1 封存预审与执行回执协议，没有读取 cash/。
- 重新计算 M0 清单的 **181/181** 产物 SHA；48 项源码/文档 **42 项一致**，6 项差异恰为封存回执列出的导航漂移：PROJECT_CONTEXT、README、进度 JSON/Markdown、scripts/README、vertical_comm_guide。所有清单内 MATLAB `.m` 文件 SHA 均一致。这是身份复核，没有重新执行 M0 数值验收。
- 详细逐文件结果见 [m0_identity_audit.json](../results/validation/m1_convergence/20261111_m1_execution_preflight_01/m0_identity_audit.json)。本轮读取时的输入 SHA 独立保存在 `inputs_sha256.json`，后续导航更新不回写该历史快照。
- MCP 未连接，**R2025b Update 6 本轮未核实**；历史版本记录不替代当前会话证明。操作系统发现名为 MATLAB 的 PID 7528、28836，但无法确定哪个对应 MCP 计算会话，不将它们当作已确认作业 PID，也未操作这些进程。

## 2. 冻结物理与参数映射

固定 Tx=(0,0,100) m、Rx=(0,0,3) m、c=1500 m/s，Gaussian `exp(-rho^2/(2*sigma^2))`、sigma=.3 m、连续单位峰值，平面 eta=0、R0=-1；直达/反射声程 97/103 m。CPU double、单线程串行、无 bubble/Doppler。此处不是 Tx99 m 的 Bellhop 展示安装，不包含海底边界。

生产入口仍应为 `vertical_channel_model`。真实代码字段为 `xw/yw`、`nx/ny`、`stepz_lamb`、`sponge_ratio`、`alpha_max_np_per_m`。无吸收层 `sponge_ratio=0` 需显式 `validation_allow_extended_sponge_ratio=true` 才在允许范围内。末平面配置为 `save_mode='slice', nout=1`；表面全步诊断需 `surface_wavefield_diagnostics=true`，且 `surface_wavefield_max_z_samples` 至少覆盖实际步数+1，不沿用默认256来抽样代替全步检查。

相位 helper 源码仍明确 `phase_contract_version='2.0.0'`：

\[
P_D=d\,e^{i\omega97/1500},\quad P_R=r\,e^{i\omega103/1500},\quad
H_{DSP}=\overline{P_D+P_R}\,e^{i\omega97/1500}.
\]

即 direct-DSP 分量为 `conj(d)`、`exp(-i*omega*.004)*conj(r)`；已适配 reduced 或 DSP-ready 输入不得重复共轭。执行时须保留 native/physical/DSP 的 direct、reflect、total 与元数据。本轮未产生上述新数组。

吸收层原默认由 `vertical_channel_model_impl.m` 确认为 `.12/.15`，`vertical_wape_propagator.m:1215` 定义：每侧厚度 `Lsponge=.12*W`，起点 `|x|=.5W-Lsponge`，单方向 `alpha_x=.15*((|x|-x_start)/Lsponge)^3`（只在层内、限幅0至.15），y方向相同，最后 `alpha_xy=alpha_x+alpha_y`。因此角点可接近 **.30 Np/m**，.15不是二维总最大值。传播 screen 的实衰减因子为 `exp(-alpha_xy*ds)`，位于现有分裂步乘法中。候选 W=192.1875 m 时每侧厚度23.0625 m、起点73.03125 m。未改吸收系数、比例或有效支持，也未评估其实际影响。

冻结配置和18项矩阵见 `frozen_config.json`、`job_matrix.json`。该 JSON 是待执行需求清单，显式标记为未验证运行配置，不冒充已执行入口或数值结果。

## 3. 逐项执行结果

`未派发`不等于数值失败，也不等于PASS。所有误差、最差位置、全步边缘和耗时均为 NA，原始表见 [metrics.csv](../results/validation/m1_convergence/20261111_m1_execution_preflight_01/metrics.csv)。

| 作业 | f/kHz | W/m | N | dx/m | stepz_lamb | 结果 |
|---|---:|---:|---:|---:|---:|---|
| P4 | 4 | 192.1875 | 984 | .1953125 | .5 | MCP前置门失败，未派发 |
| P6 | 6 | 同上 | 984 | 同上 | .5 | 未派发 |
| P8 | 8 | 同上 | 984 | 同上 | .5 | 未派发 |
| Z4 | 4 | 同上 | 984 | 同上 | .25 | 未派发 |
| Z8 | 8 | 同上 | 984 | 同上 | .25 | 未派发 |
| A-C4 | 4 | 192.1875 | 984 | .1953125 | — | 未派发 |
| A-C6 | 6 | 同上 | 984 | 同上 | — | 未派发 |
| A-C8 | 8 | 同上 | 984 | 同上 | — | 未派发 |
| A-W4 | 4 | 255.859375 | 1310 | .1953125 | — | 未派发，参考边缘资格未知 |
| A-W6 | 6 | 同上 | 1310 | 同上 | — | 同上 |
| A-W8 | 8 | 同上 | 1310 | 同上 | — | 同上 |
| A-G4 | 4 | 192.1875 | 1968 | .09765625 | — | 未派发，参考边缘资格未知 |
| A-G6 | 6 | 同上 | 1968 | 同上 | — | 同上 |
| A-G8 | 8 | 同上 | 1968 | 同上 | — | 同上 |
| A-R4 | 4 | 383.984375 | 1966 | .1953125 | — | 未派发，参考边缘资格未知 |
| SP4 | 4 | 192.1875 | 984 | .1953125 | .5 | 第一阶段未合格，禁止派发 |
| SP6 | 6 | 同上 | 984 | 同上 | .5 | 同上 |
| SP8 | 8 | 同上 | 984 | 同上 | .5 | 同上 |

AS 待执行入口为 `exact_angular_spectrum_one_step_vertical`，必须独立生成源、不调用生产 marching；计划保存97/100/103 m。固定 rho≤.5/1/2 m，共同格点直接取样，不插值、重归一化或逐路径移峰。A-W/A-G/A-R 自身边界门尚未运行，不能因参考互相接近而获资格。

## 4. 四项结论、最差误差及边界门

| 核心问题 | 本轮结论 | 缺失证据 |
|---|---|---|
| 横向窗口192.1875 m是否足够 | 未取得资格 | 三频扩窗对照、A-W/A-R参考自身边缘门 |
| dx=.1953125 m是否足够精细 | 未取得资格 | 三频A-G嵌套网格对照及参考边缘门 |
| .5 lambda纵向步长是否足够 | 未取得资格 | P4/Z4、P8/Z8、PE与同格AS对照 |
| 无吸收边界及原吸收层是否合格 | 均未取得资格 | 每一传播步边缘检查、三频SP对照 |

4/6/8 kHz 每频的 direct、reflect、total 最差复L2、最大幅度差、相位RMS及对应坐标均 **NA/未运行**；全步边缘能量最大值、边界最大幅度及所在步数也均 **NA**。不能将其记为0或“0失败”。待执行门不变：PE对同格AS≤1e-10，接口闭合≤1e-12；中心复L2≤1e-3、幅度≤.05 dB、相位RMS≤.01 rad；全步外围5%能量<1e-4、边界幅度≤-40 dB。相消总场使用直达标度绝对复差，不除近零总场。

6 kHz .25 lambda 未安排直接生产实测。均匀、无吸收、平面常系数下离散传播算子 `U(a)U(b)=U(a+b)` 的解析适用性仍保留，但本轮6 kHz同格AS亦未运行，不能写成实测通过或用解析性质越过全部生产门。原随机面窗口历史PASS未升级为本轮平面三频PASS。

## 5. 资源、日志和预算

用户本轮预算为累计14400 s、单作业3600 s，P4/SP4先导各>900 s停止后续，逐作业至少保留30%启动时可用RAM，磁盘至少2 GiB；大AS需至少1.5 GiB可用RAM。原预审10800 s是旧建议，本轮以14400 s为准；不得通过切换执行方式绕过停止门。

Node 原生OS静态快照：总RAM16,370,106,368 B（约15.25 GiB），可用1,782,767,616 B（约1.66 GiB），E盘可用180,852,674,560 B（约168.43 GiB）。CIM RAM读取被操作系统拒绝，错误与替代测量均已记录，不采信其零值占位。该快照不替代未来逐作业资源门，亦不是峰值实测。

实际 PE/AS 求解耗时 **0 s**；18作业均未启动，无MATLAB作业峰值。MCP连接探测未独立计时，完整附件返回保存在 `mcp_attach_response.json`，不能把它记作MATLAB执行0秒成功；文本预检与包装各有脚本计时。预算未因未执行而记为耗尽，也不能报告完整M1用时。日志是工具错误/壳层回执，**不存在MATLAB diary、MAT原始数据或数值审计**，不存在MATLAB进程退出码；字段明确null/NOT_RUN，不造数据满足清单。

## 6. 推荐配置与下一步

**当前没有本轮取得完整资格的生产推荐配置。** 192.1875 m、984²、dx=.1953125 m、.5 lambda、无吸收层仍只是冻结的待验证候选；原`.12/.15`吸收层没有被判合格或不合格。最终不能标记 `M1_COMPLETE_SPONGE_QUALIFIED` 或 `M1_COMPLETE_NO_SPONGE_ONLY`。

恢复条件是重新连接**现有** MATLAB MCP 会话，核实版本 R2025b Update 6、真实PID、单线程、可靠日志/状态追踪及当前资源；校验本检查点绑定的输入与源码，再从尚未运行的P4开始。已获本任务范围与预算无需逐子步骤重新申请，但本轮工具不可用时不另开进程、不自动远程迁移。按实测P4耗时更新全部剩余预测，不能装入预算则停止。第一阶段全部通过才能派发SP4/SP6/SP8，原吸收层失败时保留合格无吸收层，不调参找PASS。本任务止于M1，M2/CIR、M3、通信、深海及正式S5/S6仍未启动。

可复现的首个连接检查（现有MCP中执行；本轮失败，未实际执行）：

```matlab
fprintf('VERSION=%s\nPWD=%s\nPID=%d\n',version,pwd,feature('getpid'));
[u,s] = memory;
fprintf('AVAILABLE_RAM=%.0f\nTOTAL_RAM=%.0f\nTHREADS=%d\n', ...
    s.PhysicalMemory.Available,s.PhysicalMemory.Total,maxNumCompThreads);
```

通过连接后才允许 `setup_vertical_project` 及符合完整监控合同的M1作业；本轮没有伪造可直接运行的完整M1 runner。独立AS helper 的签名已存在，但其单独调用不能替代生产和全步资格验证。SHA身份审计与文本打包脚本仅用于本次回执，不是数值入口，不覆盖旧run。

## 7. 文件与审计交付

新增本报告、隔离结果目录和18作业未运行指标表；更新 PROJECT_CONTEXT 当前状态、scripts/README 导航和进度JSON，保留旧失败与停止原因。传播、信道、相位合同、公共接口和方法公式均未改变，因此没有标量/宽带通信回归要求所针对的实现改动，也不修改 `vertical_comm_guide.md` 的物理方法。MCP不可用的检查和正式MATLAB看板渲染均记录未执行；按维护文档允许的纯文本JSON/模板包装刷新两个生成视图，不运行任何物理计算。

源/产物SHA清单见结果目录的 `inputs_sha256.json`、`artifacts_sha256.csv`；最终审计与检查点分别为 `delivery_audit.json`、`checkpoint_final.json`。SHA验证只证明身份与交付闭合，不证明数值资格。生成视图引用文件存在性检查单独记录，不升级科学结论。
