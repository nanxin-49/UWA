# Bellhop–BIE 与 PE–Bellhop–BIE 对比图册

本报告导航至独立图册与完整复数数据。实验固定为二维、4 kHz、c=1500 m/s、Tx100/Rx3、Gaussian sigma=0.3 m与匹配Bellhop X线源；三种海面只改变A=0.01/0.05/0.20 m，K=0.10 rad/m，C2平台42 m、支撑50 m。

2026-10-19新增[原始二维声场扩展](pe_bellhop_bie_xz_visualization_report.md)：直接展示强海面u(x,z)的幅值、相位、幅相合成和原始剖面。新扩展保留原生源标度与相位说明；下方原有三案例G比值图册继续用于既定定量对照，原文件未覆盖。

- [中文图册、逐图说明与指标](../results/visualization/pe_bellhop_bie_comparison/20261018_194331_523/comparison_report.md)
- [自包含复数绘图数据 MAT](../results/visualization/pe_bellhop_bie_comparison/20261018_194331_523/plot_data/comparison_data.mat)
- 运行目录：`results/visualization/pe_bellhop_bie_comparison/20261018_194331_523/`。

比较量统一为接收端**反射场的粗糙/平面复数比G**。PE是原始Model-0；BIE使用活动误差分解MAT内的保存派生曲线，没有重新求解或重新进行收敛认证。原R3–R5活动目录缺失，另以`results/archive/`旧工作树保存的原始数值作交叉核对，三个案例PE/BIE与活动副本逐点一致；未访问`cash/`。

Bellhop平面及三个粗糙面各补算5001/10001波束，共8个独立作业。4097点壁面剖面覆盖[-80,80] m，声线步长0.05 m；主图使用10001波束。所有方法共用984点接收网格及Stage-0固定M99的411点入射能量权重。

## 图片目录

|图号|文件|内容|
|---|---|---|
|01|`figures/01_environment`|收发几何、三档海面及比较量|
|BB01|`figures/bellhop_bie/BB01_fields`|Bellhop/BIE的归一幅度与主值相位|
|BB02|`figures/bellhop_bie/BB02_residuals`|Bellhop/BIE幅度比差值dB、Bellhop−BIE相位mrad|
|BB03|`figures/bellhop_bie/BB03_metrics`|固定权重复L2、相位及幅度指标；波束敏感性|
|TR01|`figures/three_way/TR01_fields`|PE/Bellhop/BIE三方幅相曲线|
|TR02|`figures/three_way/TR02_residuals`|三组配对的复差幅度和有符号相位差|
|TR03|`figures/three_way/TR03_metrics`|三组配对指标的数量级与海面高度关系|

每张图均输出2400×1350白底PNG与矢量PDF。主值相位跨越±pi处保存并显示断点；不逐案例调相位，不拟合幅度，不做峰值归一化。相位图使用mrad，1 mrad=10^-3 rad。

## 复现与来源

```matlab
setup_vertical_project;
run_dir = fullfile(project_root, 'results', 'visualization', ...
    'pe_bellhop_bie_comparison', '20261018_194331_523');
generate_pe_bellhop_bie_comparison_visuals(run_dir);
generate_pe_bellhop_bie_comparison_visuals(run_dir, fullfile(run_dir,'rebuilt_figures'));
audit_pe_bellhop_bie_visualization(run_dir, fullfile(run_dir,'rebuilt_figures'));
```

其中`project_root`设为当前仓库绝对路径。纯绘图入口只读取本数据包，不启动声学计算，也不使用MAT内历史工作树路径。CSV保留实部/虚部、固定权重和掩膜、幅相变换、相位断点及配对残差。

专用internal-wall程序从`results/archive/cl-verification-worktree-unmanaged-20260930/validation/bellhop_internal_pm_wall_poc/bin/bellhop_iwall_pm_2020.exe`复用。SHA-256为`81de6a512037d63f06ede0537d5a591fde5e41b4ab9a95ff57f53f66917a5360`，与活动构建清单完全一致，两份overlay源码哈希也一致。用户给出的`E:\MISC\BELLHOP\AcousticsToolbox_2020`是工具箱目录；本轮实际执行的是上述经哈希核对的专用程序，没有使用普通Bellhop代替、安装编译器或重建程序。

`source_data/`保存源MAT、原报告、构建清单、overlay源码、程序、输入/输出和作业检查点。`metrics/`将本轮重算指标、历史报告值、平面控制、波束敏感性及实际进程耗时分开。并行进程耗时不是串行基准，也不据此主张方法加速。

本轮只添加验证与绘图工作流；生产API、冻结求解器和物理公式未改，因此不运行无关的生产信道/通信回归，也不改`vertical_comm_guide.md`。结论只对应三个确定性二维海面，不外推PM、三维或通信性能。流程归因见完整中文图册。

## 实际结果与交付检查（2026-10-18）

|海面|本轮PE–BIE复L2|本轮Bellhop–BIE复L2|5001/10001波束比值变化复L2|
|---|---:|---:|---:|
|弱 A=0.01 m|0.0031575715|8.0919750e-6|3.5113444e-6|
|中 A=0.05 m|0.016089960|3.8895510e-5|4.8939664e-6|
|强 A=0.20 m|0.080908108|9.9760075e-5|1.1121732e-5|

在此固定夹具中，Bellhop与保存BIE的差异均小于原PE相位屏，PE差异随高度增大；不由此推断其他频率、海面或维度的普适优劣。10001波束平面控制与保存Stage-0平面场的M99复差为0。

Bellhop有限角覆盖外有375点平面场为零，粗糙/平面比无法定义；显式以复NaN及`ratio_valid`标记，保留原始SHD，无插值。它们全部位于M99以外，既定411个比较点全部有效，未缩减支持。

- 8/8作业正常退出；壁面交点残差最大约3.54e-12 m、反射相位跳跃误差约7.11e-15 rad，接收坐标误差0；p/q旋转及正向传播诊断通过。
- 7张PNG均2400×1350；7张矢量PDF均单页1152×648 pt。逐张查看PNG及Poppler渲染PDF，标题、图例、单位和中文说明无裁切。
- 仅复制绘图MAT至`rebuild_verification/plot_data/`，使用包内快照脚本重绘，7/7 PNG像素完全一致；指标重算差0。
- 复数CSV往返最大绝对误差4.8937e-15，native表示映射差0，输入快照及程序哈希通过。[审计记录](../results/visualization/pe_bellhop_bie_comparison/20261018_194331_523/metrics/artifact_audit.json)。
- 三个新增MATLAB文件Code Analyzer均无消息；`git diff --check`通过。新增PE求解0次，BIE求解0次，Bellhop作业8个。

完整数据包同时保留`source_data/parameters.json`与`source_data/workflow_code/`。脱离项目路径重绘时，仅需将绘图MAT放在`run_dir/plot_data/comparison_data.mat`，把包内`workflow_code`加入MATLAB路径，再调用绘图入口；本轮独立重绘已采用此方式。
