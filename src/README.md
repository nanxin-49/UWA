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
