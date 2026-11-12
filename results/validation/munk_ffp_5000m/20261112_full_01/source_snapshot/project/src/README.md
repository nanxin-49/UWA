# `src/` 模块索引

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
