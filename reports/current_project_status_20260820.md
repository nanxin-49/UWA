# 项目技术状态报告（2026-08-20）

## 总体判断

项目已经形成一条可运行的公共 PE/WAPE 信道链，以及一条固定环境下的 joint-kstat + cached/adjoint + 解析统计研究链。载波相位修复、精确离散伴随、F=64 二阶统计、U=5/U=8 条件模型和 two-node 通信验证已有正式 PASS 证据。当前主要风险不在这些公式本身，而在公共横向窗口仍沿用兼容默认、随机海面宽带收敛未闭环，以及若干外部物理基准的绝对幅度与孔径问题。

## 已完成

- `direct_dsp` 统一直达与反射的载波参考；标准几何名义相对时延为 4 ms。
- F=65 相位审计、旧 schema 条件模型迁移和公共字段闭合通过。
- uniform CPU-double 单最近网格接收机的 exact discrete adjoint 与投影通过。
- PM 中央裁剪伴随零嵌、dense/FFT `C_H/P_H` 及 F=9/F=64 样本统计通过。
- U=5/U=8 schema 2.x 条件模型、条件库和 two-node full 通信通过。
- Gaussian 直达窗口在 160 m/no-sponge 收敛。
- 固定海面完整反射链在实际 192.1875 m/no-sponge 严格收敛。
- 4 kHz 随机海面集合中，192.1875 m 对 15/15 realization 通过。
- 根目录公共调用名称保留，核心实现已按功能迁入 `src/`。

## 部分完成或仅诊断

- 160.15625 m/no-sponge 的反射接收响应近似收敛，但未通过外围能量严格门限。
- SSA1 只实现第一阶 Dirichlet/微扰极限统计核，不是 NLSSA 或高阶 SSA。
- Li2009 为部分复现，尚未完成横向网格收敛。
- 2k 验证只诊断相位系数近似和最低阶海谱/SSA 趋势，不验证完整 PE 或真实海洋散射。
- 已有传播图册解释了 reduced complex envelope、海面前后场、接收平面、cached/adjoint、C/P、PDP 和 LFM，但它不是新的独立物理基准。

## 未完成与开放问题

- 随机海面 3–5 kHz 候选—256 m 参考的完整成对集合未完成；不能报告 ensemble 最差群时延。
- Bellhop 绝对幅度、Gaussian/点源归一化和无限孔径等价仍开放。
- 公共默认横向窗口/sponge 尚未切换到严格推荐配置。
- adjoint/statistics v1 尚不覆盖 layered、bubble、Doppler、GPU、多接收机或插值接收。
- 未对未验证风速做条件模型插值；U=5/U=8 只是两个已登记节点。

## 默认与推荐配置差异

| 用途 | 当前公共默认 | 当前验证建议 | 说明 |
|---|---|---|---|
| 相位参考 | `direct_dsp` | `direct_dsp` | 已一致 |
| surface model | `kirchhoff_spatial` | 按研究目的显式选择 | 默认未改变 |
| 横向窗口 | 50 m | direct 160 m；完整反射 192.1875 m | 推荐尚未进入默认 |
| sponge | ratio 0.12 / alpha 0.15 | 上述严格窗口使用 off | 小窗口 sponge 会改变接收响应 |
| bubble | off | 仅在专门实验中启用 | 伴随 v1 不覆盖 |

## 接下来建议

1. 在不改公共默认的前提下，完成随机海面宽带成对窗口验证。
2. 独立处理公共默认窗口迁移的计算成本和兼容策略。
3. 关闭 Bellhop 绝对幅度/源归一化问题，再决定外部基准的验收权重。
4. 若继续 SSA/Li2009，先完成横向网格与孔径收敛，不把 2k 诊断升级为物理验证结论。

## 权威证据

- `reports/pe_phase_reference_release_candidate_report.md`
- `reports/adjoint_pe_receiver_projection_feasibility_report.md`
- `reports/pe_reflected_chain_window_validation_report.md`
- `reports/pe_random_surface_window_robustness_report.md`
- `reports/validation_2k_phase_report.md`
- `reports/validation_2k_ocean_spectra_report.md`

目录整理本身不改变任何传播、相位、统计或通信公式，也不重新声明昂贵历史验证的数值有效性。

## 目录迁移后的轻量回归

运行环境：MATLAB R2025b Update 1，CPU double；未重跑昂贵 F=64/512 ensemble、完整 Bellhop 矩阵或 192/256 m 正式窗口验证。

| 检查 | 配置/结果 | 状态 |
|---|---|---|
| setup 与路径遮蔽 | 从根目录和外部目录初始化；6 个公共/核心符号各只有一个定义 | PASS |
| MATLAB 静态检查 | 155 个 `.m` 文件，177 条风格/分析提示，语法错误 0 | PASS |
| scalar direct-only | PE 64²、16 m、20→2 m、6 kHz；反射严格为 0 | PASS |
| scalar direct+reflect | 同网格 Kirchhoff spatial；`H` 分量闭合 `3.47e-18`，参考标量闭合 0 | PASS |
| 公共包装透明性 | `vertical_channel_model` 与 impl 同输入差异 0 | PASS |
| seeded/override/diagnostics | 64² fixed-seed；surface override、端点、接收点和诊断透明性全部通过 | PASS |
| cached/adjoint smoke | PE 32² / PM 64²，4/6/8 kHz；内积 `3.35e-15`，reduced/direct-DSP 投影 `1.04e-15/1.05e-15` | PASS |
| 宽带通信 smoke | F=9、PE 64²、QPSK+AWGN；IFFT `2.96e-16`，分量闭合 `5.55e-17`，无 NaN/Inf | PASS |
| 2k 迁移复现 | 两个脚本从新位置完成；5 个关键 CSV 的 SHA-256 与迁移前逐文件一致 | PASS |
| 根目录清洁 | 仅项目级文档、setup、4 个兼容入口及 Git 配置文件 | PASS |

若干独立 MATLAB 批处理（包括冷启动 setup 和 PE FFT 检查）在全部断言通过并打印上述数值后，于 R2025b 进程退出阶段触发 DDUX/threadpool `std::terminate`；通信与 2k 批处理正常退出。该现象发生在结果完成后的 MATLAB 关闭阶段，未观察到数值或文件损坏，但应作为本机 MATLAB 运行时问题保留记录，不把它伪装成零退出码。
