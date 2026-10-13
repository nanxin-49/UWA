# 三维 FullFK 单变量振幅适用域扩展

任务日期：2026-10-13。**新增比例0.5、1.5、2.0均获条件性 G1 PASS；1.0复用冻结基准。
未发现可信物理失败或高阶参考失稳。全G1目标域仍OPEN，WP2=false。**

## 固定条件与方法

Tx1500/Rx5、q0、6 kHz、c1500；原dual波矢(.4,0)/(.2,.6) rad/m、phase(.2,.7)，
C2 plateau8/support12 m、原生exp(-i omega t)、单位Green点源完全不变。
唯一物理变量为(.02,.01) m的统一比例alpha，不重新采样相位、海面或归一化。
冻结的FullFK、Dirichlet FE、递推、源和海面helper未修改；比例1.0无新增求解。

按[预设协议](deepsea3d_g1_amplitude_protocol.md)，每个新增点独立完成N12/16、
L24/48、h=.125/.0625、M32/40和真实粗糙面原网格/半网格检查。若16阶不足，
入口可追加冻结允许的20阶；本轮16阶实际均足够，并非默认假定。
仅在0.5、1.5参考护栏全部通过、最大单次计算约30秒后，才安排2.0。

原门槛保持：参考各项数值变化/真实边界残差<1e-4；FK数值变化<1e-3；
可信参考下反射复误差及非零散射增量复误差各<1%。没有幅相拟合或放宽门槛。
原相位差为`arg(P_FK conj(P_ref))`，散射增量为`P_ref-P_flat`；不借平面背景
稀释散射误差。参考不收敛时只能INCONCLUSIVE，不能判FK物理失败。

## 物理比较

|alpha|最大高度m|最大局部坡度|反射复误差|增量复误差|反射相位差rad|增量相位差rad|判定|
|---:|---:|---:|---:|---:|---:|---:|---|
|0.5|.015|.007637|6.189e-5|9.128e-5|-6.148e-5|-3.861e-5|条件PASS|
|1.0|.030|.015273|1.361e-4|1.056e-4|-1.358e-4|-7.102e-5|冻结PASS复用|
|1.5|.045|.022910|2.066e-4|1.175e-4|-2.065e-4|-9.849e-5|条件PASS|
|2.0|.060|.030546|2.565e-4|1.267e-4|-2.557e-4|-1.218e-4|条件PASS|

反射误差和相位偏差在四个点上随振幅增长，反射幅度误差最大仅2.117e-5；
主要差异表现为相位。散射增量误差增长较缓。高度及坡度同时由同一个alpha改变，
不能由此区分哪个是唯一控制量，也不能拟合普遍物理适用公式。
坡度/高度统计使用固定r<=12 m网格；RMS、幅度误差及完整精度在
[summary.csv](../results/validation/deepsea3d_g1_amplitude/summary.csv)。

## 独立参考的实际稳定证据

|新增alpha|粗网格N12→16增量变化|细网格N14→16增量变化|网格变化|域变化|原网格边界残差|半网格残差|
|---:|---:|---:|---:|---:|---:|---:|
|0.5|1.99e-12|4.66e-15|7.60e-7|1.00e-7|1.71e-11|1.94e-7|
|1.5|1.25e-6|2.75e-8|1.33e-6|1.47e-7|2.97e-10|2.16e-7|
|2.0|4.62e-5|1.83e-6|1.51e-6|1.78e-7|3.19e-8|2.44e-7|

细网格的N14→16由同次输出的完整系数作严格截断部分和，不是另一次随机输入。
另外每个点都有单独N12/N16运行及边界残差。2.0的N12原网格残差已升至2.56e-5，
而N16降至3.19e-8；高阶检查确有必要。没有把16阶推广为所有振幅的充分条件。
M32/40边界评价变化记录为0；这是浮点记录，不是严格零余项。
有限Fourier指数评价尾界另存，**不含浮点误差，也不是FE或连续PDE余项界**。

新增点FK自身网格增量变化最大9.32e-9，Gauss/梯形差最大1.73e-7。
最终参考的细阶数/网格/域变化约1e-6量级，比FK约1e-4偏差更小；粗阶数差反映
被舍弃的12阶误差，不能误当最终16阶参考的严格误差界。
当前可信程度与冻结基准一致：有多重数值收敛及真实边界证据的条件性参考，
没有连续算子误差定理或任意C2海面的收敛半径证明。

证据：[guards.csv](../results/validation/deepsea3d_g1_amplitude/guards.csv)、
[全部18个参考运行](../results/validation/deepsea3d_g1_amplitude/reference_runs.csv)、
[逐阶原生复系数](../results/validation/deepsea3d_g1_amplitude/reference_terms.csv)、
[源Taylor及评价尾界](../results/validation/deepsea3d_g1_amplitude/evaluation_checks.csv)。
4项分类回归检查通过，覆盖“参考/数值失败不得判物理FAIL”和“增量不能被背景隐藏”。

## 成本、复现及后续决定

三个新增点保存的FE/FK计算wall分别46.30/49.30/41.75秒，合计**137.35秒**；
CPU合计785.70秒。1.0无新增求解，不重复计算或计入历史137.45秒基准成本。
最大单次参考约29.94秒；最终网格1536²、16个复谱载荷576 MiB，另有临时/forcing
数组，不是进程峰值内存。成本不含IO、几何统计、测试、分析和文档。
[成本记录](../results/validation/deepsea3d_g1_amplitude/cost_summary.json)。

setup后调用`validate_deepsea3d_g1_amplitude('init')`，重复`'next'`至complete=1；
每次最多一个参考解或一个FK检查块，可从checkpoint续跑。audit前运行并保存分类测试：

```matlab
info=setup_vertical_project;
tests=runtests('test_deepsea3d_amplitude_classification');
assert(all([tests.Passed]));
writetable(table(tests),fullfile(info.project_root,'results','validation', ...
    'deepsea3d_g1_amplitude','classification_tests.csv'));
validate_deepsea3d_g1_amplitude('audit');
```

输出`amplitude.mat`仅含数值/字符串struct，并保留冻结基准数值状态；
来源/结果hash和验收JSON在本case目录。新文件仅为编排、分类和报告，无冻结模型变更。
工作流沿用MATLAB技能，归属文献见[冻结基准报告](deepsea3d_g1_finite_report.md)末尾。

**实际观察范围是四个离散振幅点0.5/1/1.5/2，不能认证整个[0.5,2]连续区间。**
未找到本谱/几何下的可信1%失败边界，alpha>2仍未知。若继续振幅方向，可先对2.5
进行低成本N16/20稳定性预检；2.0的粗阶数变化已明显增大，不宜直接大批扫描。
参考不稳则保留INCONCLUSIVE并再评估稳定参考方法。全目标域G1保持OPEN，WP2不自动放行。
