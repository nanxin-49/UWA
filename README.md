# 垂直水声 PE/WAPE 信道与通信项目

本仓库用于研究海底发射机到近海面接收机的垂直水声传播。公共模型以 PE/WAPE 计算直达与海面反射频响，并提供粗糙海面、joint-kstat、接收端精确离散伴随、二阶统计和 MPSK 通信验证工具。

## 快速开始

在 MATLAB 中先初始化一次路径：

```matlab
cd('E:/MISC/CARPE3D_matlab/Explain')
setup_vertical_project;
```

最短调用：

```matlab
paramsV = struct('f0',6000,'enable_wideband',false, ...
    'enable_surface_reflection',true,'show_figures',false);
output = vertical_channel_model(paramsV);
```

也可以运行兼容入口：

```matlab
main_vertical              % 信道演示
explain_main_vertical      % main_vertical 的兼容别名
comm_main_vertical_psk     % MPSK 通信演示
```

这些名称保持不变；实际演示脚本位于 `examples/`，核心实现位于 `src/`。

## 公共入口与目录

| 路径 | 用途 |
|---|---|
| `vertical_channel_model.m` | 稳定公共 API 包装；初始化路径并调用信道实现 |
| `main_vertical.m` / `explain_main_vertical.m` | 信道演示兼容入口 |
| `comm_main_vertical_psk.m` | 通信演示兼容入口 |
| `setup_vertical_project.m` | 幂等加入 `src/`、示例和必要脚本目录，不改变工作目录 |
| `src/` | 可复用生产实现，按 channel/propagation/surface/receiver/statistics/communication/bubble 划分 |
| `examples/` | 演示脚本实际内容 |
| `scripts/` | 验证、比较、实验和报告入口 |
| `reports/` | 一次性验证报告和当前状态报告 |
| `results/` | MAT、CSV、图片、视频及检查点 |
| `docs/history/` | 非权威历史规格与旧说明，仅供追溯 |

模块和依赖方向见 `src/README.md`，脚本命令见 `scripts/README.md`。

## 当前能力与边界

| 项目 | 当前状态 |
|---|---|
| uniform CPU-double PE/WAPE，直达与海面反射 | 可用；公共默认相位参考为 `direct_dsp` |
| Kirchhoff spatial/k-domain 与 joint-kstat | 可用；公共默认 surface model 未改变 |
| cached forward、精确离散伴随、PM 周期 FFT `C/P` | 固定环境、单最近网格点接收机原型已验证 |
| U=5/U=8 条件统计生成器 | 已有 F=64 验证证据；适用于已登记固定条件 |
| bubble、SSA 诊断分支 | 已实现，但适用范围和验证强度不同，见技术指南 |
| layered、Doppler、GPU、多接收机、插值接收的伴随统计链 | 尚未纳入当前原型 |

当前公共窗口参数仍保持兼容默认值，并未自动切换到近期收敛验证的推荐窗口。正式研究运行前请阅读技术指南中的窗口建议和证据边界。

## 文档导航

| 文档 | 主要读者 | 内容 |
|---|---|---|
| `README.md` | 初次进入仓库的人 | 快速入口、目录和导航 |
| `vertical_comm_guide.md` | 研究与开发人员 | 当前物理、接口、结果语义和限制的权威说明 |
| `PROJECT_CONTEXT.md` | 后续模型与维护任务 | 顶部当前索引和完整追加式项目日志 |
| `scripts/README.md` | 运行验证的人 | 脚本入口、环境变量、成本和输出位置 |
| `reports/current_project_status_20260820.md` | 项目审阅者 | 已完成、部分完成、开放问题和推荐配置差异 |
| `reports/` | 验证审阅者 | 单次实验配置、结果和结论 |

历史文件原文保存在 `docs/history/`，不得作为当前 API、默认值或验证状态的唯一依据。

## 结果存放规则

- 验证结果：`results/validation/<case>/`
- 可视化结果：`results/visualization/<case>/`
- 一次性报告：`reports/`
- 大型正式结果不放回根目录，也不由目录整理任务自动删除或重算。

当前技术状态日期为 2026-08-20。目录迁移仅改变文件位置和初始化方式，没有改变 PE marching、相位公式、海面统计公式或公共默认配置。

## PE--Bellhop 展开坐标验证

当前新增的验证入口
`scripts/validation/validate_pe_bellhop_unfolded_flat_gaussian_vertical.m`
采用镜像展开：直达路径对应 Bellhop 97 m，平面压力释放反射路径对应
103 m 的镜像接收端并乘以 `-1`。该入口使用现有生产 Gaussian 源的角谱
指向性生成 `.sbp`，不修改 PE 核心或生产默认窗口；正式输出位于
`results/validation/pe_bellhop_unfolded_flat_gaussian/`。
