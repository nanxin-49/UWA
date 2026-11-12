# `src/` 模块索引

### 2026-11-12 活塞P0/P1最终状态

2026-11-12最终封存。活塞P0与4 kHz P1自由场、平面、弱C2正弦均通过，独立MAT/CSV/密度重建审计失败0项，80条FIG曲线核验一致，五组最终PNG已目视检查。弱面20 cm反射复场PE Model-0/BIE相对L2为0.0048932994，FullFK/BIE为2.0947621e-6，均超过五倍经验数值分辨率门。周期截窗初场失败保留，接受结果采用验证专用连续角谱平方根PE与精确平面背景加紧支撑相屏扰动。20 cm可作为该范围主源。后续生产入口与更广工况需要单独检查，本轮止于P1。 [数据与报告](../results/comparisons/piston_2d_single_frequency/20261112_remote_p0p1_01/report.md)。

`baffled_strip_piston_field_vertical`返回同一速度孔径的复压力及解析x/z梯度。验证调度入口为 `run_piston_p0p1_remote_job`，连续PE实现与独立审计位于 `scripts/validation/`，保存数据图件入口为 `render_piston_p0p1_remote`。新增FullFK入射回调保留原Fourier分支并已回归。源定义、偏移载波、相屏表达式、离散矩阵、失败门、资源及局限见最终报告。Gaussian默认和公开WAPE核心沿用原实现。

`src/` 只保存可复用实现，不放演示入口、一次性验证脚本或生成结果。公共调用仍从仓库根目录的兼容包装文件进入。

| 模块 | 责任 | 主要依赖方向 |
|---|---|---|
| `channel/` | 公共信道实现、载波相位参考、CIR 表示转换 | propagation、surface；整理各分量输出 |
| `propagation/` | PE/WAPE marching、Gaussian 初场 | surface、bubble |
| `surface/` | PM 谱、显式 Kirchhoff、joint-kstat、SSA 边界 | 不依赖公共入口 |
| `receiver/` | cached forward、精确伴随、接收投影、PM→PE 统计收缩 | propagation、statistics |
| `statistics/` | joint-kstat 联合模型、条件模型、经验模型和采样 | receiver 输出或已保存统计 |
| `communication/` | taps、MPSK、噪声和 ensemble 评估 | channel/CIR 输出 |
| `bubble/` | 气泡谱、环境和有效介质 | 由 propagation 相位屏调用 |

## 依赖原则

2026-11-12 活塞独立源API：`propagation/baffled_strip_piston_field_vertical(x_m,z_m,source,nquad)`返回峰值Pa及物理x/z梯度Pa/m，严格二维速度孔径/Hankel积分，`source`明确频率、c/rho/v0、总宽与真实源深度。它尚未作为生产默认模式接入；Gaussian生产路径与公共配置未改。当前P0部分通过但PE有限初场表示未过门、远端新进程等待许可，完整资格见[执行记录](../reports/piston_p0p1_remote_execution_2026-11-12.md)。

2026-11-11 M1逐步观察接口：公共`paramsV.validation_step_observer`默认`[]`，可选CPU验证回调`observer(psi_xy,state)`；`state`包含`path`（direct/incident/reflect）、`step`、`n_steps`、`z_m`、`frequency_hz`。回调无返回值，观察每段初场及每一传播步，可记录边缘/接收圆盘与资源并抛异常停止作业；不修改传播算子或相位合同。回调运行时间计入作业预算，保存配置前应清除函数句柄。相关M1入口见`scripts/README.md`。

2026-11-10 相位合同：`channel/apply_pe_channel_phase_reference_vertical` 的原生 reduced 输入经完整包络共轭进入 direct-DSP；`channel/assert_channel_dsp_metadata_vertical` 拦截旧项目合同。条件模型 schema 3、经验模型来源元数据及迁移规则见[指南](../vertical_comm_guide.md)和[M0验收](../reports/m0_phase_contract_closeout_2026-11-10.md)。传播/伴随/海面空间统计算子保持原定义。

```text
root public wrapper
        |
        v
channel -> propagation -> surface
             |
             +---------> bubble

receiver <-> statistics -> communication
```

- `src/` 不依赖 `examples/`、`reports/` 或 `results/`。
- 验证专用 Bellhop、Weyl、Li2009、properness、artifact 和发布元数据辅助函数位于 `scripts/validation/support/`。
- 物理公式与限制以根目录 `vertical_comm_guide.md` 为准；本文件只解释模块边界。
- 添加新模块后同步更新 `setup_vertical_project.m`。不得通过改变当前工作目录解决路径问题。
