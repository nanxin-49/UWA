# PE–Bellhop–BIE 原始二维声场可视化

用户于2026-10-19确认真实x–z声场需求并批准最小补算范围。新增工作与原三案例图册分别保存。

运行目录：`results/visualization/pe_bellhop_bie_xz/20261019_115214_348/`。
完整[中文二维场图册](../results/visualization/pe_bellhop_bie_xz/20261019_115214_348/xz_field_report.md)；
[完整复声场MAT](../results/visualization/pe_bellhop_bie_xz/20261019_115214_348/plot_data/xz_data.mat)。

## 固定范围

- 只做强海面A=.20 m、K=.10 rad/m、C2平台42/支撑50 m；4 kHz/c1500、Tx100 m、Gaussian sigma=.3 m、Bellhop X源不变。
- 横向沿用原M99的411个真实接收位置，约[-40.039,40.039] m；深度z=.5:.125:6 m，45层，每种方法18495个场点。
- PE从保存的海面反射初始场传播至新增深度。BIE从保存4480节点边界密度评价Green积分，不重新求解边界方程。
- 仅新增一个10001波束Bellhop网格作业，使用既有哈希核对的internal-wall专用程序、4097点壁面、步长.05 m；不安装工具，不修改冻结求解器或生产代码。

## 原始场与表示边界

绘制u(x,z)的幅值、相位、幅相合成图及1/3/6 m原始幅值剖面。直接使用场数组，不用G或误差指标替代声场。

PE保存约化包络，同时按解析因子exp(+i k0(100+z))恢复载波供波前显示。平面控制恢复后与BIE原生场相对差3.8550e-13，幅值不变。BIE为原生exp(-i omega t)，Bellhop沿用conj(raw SHD)原生表示。

Bellhop与Gaussian场的原始源标度及参考相位不同：旧平面中心幅值比5.7737926、相位差0.77789567 rad。仅记录该诊断，没有将其作为幅相修正系数。PE/BIE幅值共用色标，Bellhop使用显式独立原生源单位；不声称绝对声压Pa、TL或跨方法原生相位误差排序。

所有源数组、边界密度、源文件哈希、原始Bellhop输出、绘图变换保留。原三案例图册及其定量结论不覆盖、不改写。

## 入口

```matlab
setup_vertical_project;
s = prepare_pe_bellhop_bie_xz_data('init');
prepare_pe_bellhop_bie_xz_data('stage',s.run_dir); % 只生成一个作业输入
prepare_pe_bellhop_bie_xz_data('pebie',s.run_dir);
% 独立Bellhop进程完成并保存process_result.json后：
prepare_pe_bellhop_bie_xz_data('finish',s.run_dir);
generate_pe_bellhop_bie_xz_visuals(s.run_dir);
```

`finish`不启动外部进程，避免超时后重复计算。纯绘图入口只依赖新数据包的`plot_data/xz_data.mat`。
本轮是验证/绘图扩展，未修改物理方法或公共API，`vertical_comm_guide.md`不需方法修订；生产信道/通信回归不适用。

## 新增直观图

![原始幅值二维场](../results/visualization/pe_bellhop_bie_xz/20261019_115214_348/figures/XZ01_raw_magnitude.png)

三种方法的幅值主体均为中央较强、向两侧减弱的声束包络。PE/BIE共用幅值范围，Bellhop原生源单位单独给出，不能按跨列颜色直接判断绝对声压比。

![原始相位二维场](../results/visualization/pe_bellhop_bie_xz/20261019_115214_348/figures/XZ02_raw_phase.png)

相位条纹直接展示原始复声场随位置和深度变化的波前结构。颜色在±pi循环，各自源相位约定保留；整体色移不能直接当作粗糙反射模型误差。

![原始幅相合成二维场](../results/visualization/pe_bellhop_bie_xz/20261019_115214_348/figures/XZ03_complex_field.png)

颜色表示相位，明暗表示幅值；中央声束较亮、两侧较暗，波前条纹与声束包络同时可见。圆形图例和幅值上限明确给出映射，没有逐图gamma增强或峰值拟合。

![原始幅值剖面](../results/visualization/pe_bellhop_bie_xz/20261019_115214_348/figures/XZ04_raw_profiles.png)

1/3/6 m真实深度截线显示PE与BIE的幅值主体仍很接近，6 m处局部起伏差别较易辨认。Bellhop原生幅值保留在下排，纵轴从零开始，未除以平面场。

## 完成检查（2026-10-19）

- 每方法411×45=18495个实际场点，全部复数有限；4张2400×1350 PNG及单页1152×648 pt PDF完成。PDF用审核后的PNG整页封装，避免R2025b重复导出时的布局裁切；不标为矢量PDF。
- PE/BIE新3 m截线相对差分别为1.4786e-14、5.2563e-15；Bellhop与原作业逐点一致（相对差0）。接收坐标误差0，壁面交点、压力释放反射相位、p/q旋转及正向传播诊断通过。
- 保存的BIE 12/14 ppw密度在44个新域探针上的场变化1.4232e-8；只作为本轮场评价敏感性，未重做BIE边界求解/全域收敛认证。
- 唯一新增Bellhop作业wall686.2287 s；BIE含截线/探针的场评价约24.9855 s。PE仅从保存起始场评价45个深度，无上行全流程重算。这些耗时不是方法间性能基准。
- 仅复制绘图MAT至独立目录，使用包内快照脚本重绘；4/4 PNG逐像素一致。CSV最大绝对往返误差7.1089e-16，载波和RGB映射公式复核差0。[审计记录](../results/visualization/pe_bellhop_bie_xz/20261019_115214_348/metrics/artifact_audit.json)。
- 新增4个MATLAB文件Code Analyzer无消息；原有生产/冻结求解器未改。README、脚本索引和项目日志已更新；软件流程归因见完整图册。
