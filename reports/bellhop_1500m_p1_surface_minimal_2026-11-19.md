# 修订P1最小海面反射：数值门通过，新增图件审计失败后停止

2026-11-19。工作树experiment/bellhop-rotated-fig12b。

**历史图例修复已完成；唯一S/4000Hz/C801的物理/数值检查通过，D+S半空间镜像比较也通过。随后新增S图的图例布局审计触发断言，MATLAB外层退出1，统一流程以失败状态停止。未把本轮标记为完整封存通过，也不宣布整个P1通过。**

原自动失败报告另存结果目录automatic_failure_report.md；原MATLAB日志、失败FIG/PNG、checkpoint_final.json和supervisor_receipt.json不修改。本文只整理已保存证据，没有重新执行Bellhop、修改公式或新增审计重跑。

## 1. 先行完成的历史图例修复

原201诊断图及201/401/801对比的三组PNG/FIG已另存新版本，显式绑定图例句柄、关闭AutoUpdate并固定位置，去除data1和裁切。24条曲线从新FIG重开后与保存数据逐项精确一致，4个图例均在画布内，三张PNG已实际查看，6项原图件SHA保持。没有重算201/401/801。

[修复报告](bellhop_D_control_legend_repair_2026-11-19.md)、[修复目录](../results/visualization/bellhop_D_control_legend_repair/20261119_01/)、[目视审计PASS](../results/visualization/bellhop_D_control_legend_repair/20261119_01/visual_audit.json)。修复时本地既有MATLAB曾提示许可证警告，但函数完整返回审计通过；用户随后恢复，本轮S代码静态检查成功。未关闭、清理或重新启动共享会话。

## 2. 冻结配置与唯一新增任务

[执行前冻结计划](bellhop_1500m_p1_surface_minimal_plan_2026-11-19.md)及结果目录frozen_task.json保持。原201束P1失败不重标、不删除；修订候选C束数=801，仅授权S4000/C801，不修改生产默认。

U0：水深1500m，Tx=(0,1490)m，Rx=(x,6)m，x=[0,10,50,100,200]m，是五个独立链路测试位置。c1500m/s，水密度1g/cm³，水吸收0，平面压力释放海面；流体底cp1600m/s、rho1.8g/cm³、损耗0，但选定D/S两分量没有底反射。理想全向3D点源、1m压力参考、源宽N/A。

原S201输入仅将ENV束数行201改为801；SSP、R25D、IWCOV、ROTQ逐字节保持。S展开距离1496m，原对称扇区±8.37618313623413°，步长0.5m，Q间距0.625m，REAL8，正式C ELF SHA `5d87293a27e7fafc3aedefb569f3f35065294f42426503c4881d6eb55e6769e9`，固定`p_abs=-conj(p_native)`及物理时间因子`exp(-iωt)`。D4000/C801只读复用。

## 3. 独立几何与压力释放检查

801束全部恰好一次海面、零底反射，所有17项独立事件几何/系数检查通过。实际第一反射事件与恒声速平面镜像几何逐束比较，未使用示意声线代替数据。

|检查|实测最大误差|门槛|
|---|---:|---:|
|海面物理深度|0m|1e-4m|
|反射点水平坐标|2.81091e-11m|1e-4m|
|反射方向|3.33067e-16rad|1e-6rad|
|展开路径长度|2.01908e-10m|1e-4m|
|映射接收range的路径时延|1.34670e-13s|1e-7s|
|压力释放复系数与−1之差|6.98296e-15|1e-12（双精度日志重建审计）|
|相位跃迁与π之差|7.10543e-15rad|1e-12|
|反射前后幅度差|0|1e-12|

实际求值和缩放坐标REAL8逐位检查、原5点顺序恢复、proper-rotation正反映射通过，S单分量原有19项检查全部通过。完整801条第一事件含入射/反射/旋转方向、p/q、tau、Amp、Phase、线段/参数均保存在all801_first_surface_states.csv。

五个Rx的独立镜面反射点为`x_s=1490x/1496`，`L_S=hypot(x,1496)`。轴上S解析时延997.333333ms，D为989.333333ms，相差8ms。五位置几何CSV是解析参考；本轮没有运行A/R或对五点另作本征声线打靶，不把C事件日志称为ARR或完整RAY。

## 4. S与D+S复声压验收

`p_S,ref=-exp(ikL_S)/L_S`；`p_DS,num=p_D801_saved+p_S,num`；`p_DS,ref=p_D,ref-exp(ikL_S)/L_S`。先复数相干相加，没有相加dB或包络，没有增益/相位拟合。

|场|复数L2 ≤0.001|max复差/E ≤0.002|max有效相位/rad ≤0.01|相位支持|
|---|---:|---:|---:|---|
|S|0.000343625498|0.000406380934|0.000406380935|5/5|
|D+S|0.000392530118|0.000305833057|0.000412921442|3/5：x=50,100,200m|

六项复场门均通过。S包络为|S解析|；D+S独立包络`E=|D解析|+|S解析|`。全部五点保留原复压并参与复数L2与包络归一化复差；不使用近零总场作为逐点相对分母。相位严格按双方模均≥0.01E筛选，x=0和10m不进入D+S相位门，没有删除这些点的复压数据。

|x/m|S复差/E_S|D+S复差/(E_D+E_S)|D+S相位是否评价|
|---:|---:|---:|---|
|0|1.04035e-8|3.65158e-8|否|
|10|0.000344321280|3.42787e-6|否|
|50|0.000406380934|2.35580e-5|是|
|100|0.000393413079|8.77365e-5|是|
|200|0.000390820114|0.000305833057|是|

轴上原始D+S数值约`−2.702581696e-6 + i4.681060091e-6`（1m压力参考单位），独立镜像参考约`−2.702624787e-6 + i4.681083450e-6`，包络约0.001342303645。相消后仍有不同1/L造成的有限余量，不是精确零。

D+S仅标记`U0_D_S_pressure_release_halfspace_control`，不是含海底反射的完整U0波导场、无限阶通信信道或CIR。当前数值证据支持本固定单频最小海面控制，但新增图件审计失败，未完成本轮完整封存；整个P1仍未通过。

## 5. 停止原因与时间/资源

唯一Bellhop PID1831605，真实退出0，native0.824823s，stdout/stderr为空，PRT无FATAL。远端MATLAB R2025b Update1许可检查通过，随后外层真实退出1，原因是execute_bellhop_1500m_s801.m第49行中新S图例位置/AutoUpdate联合断言失败，不是声学求解或上述数值门失败。

失败发生在comparison.mat、所有复压/几何CSV、FIG/PNG已经保存之后，comparison_receipt.json和完整自动数值报告尚未来得及写出，故自动回执把任务记为EXECUTION_STOPPED_NO_QUALIFICATION。原回执保持，不追改为PASS。新S图件保留供人工审查，本轮没有在失败后重跑绘图或算例，也没有改曲率公式。

统一监督器在35.199264s内以失败状态终止，未触发300s总时限；295s动作截止、远端170s、native120s保持。启动可用RAM130826506240bytes，监督采样最低128822185984bytes，高于30%保留39247951872bytes及256MiB门；磁盘和输出门未失败，远端输出3840128bytes。启动RAM≥1.5GiB。没有中断共享会话或其他用户进程。

79项返回SHA通过，375项历史SHA执行前后全部一致。旧201失败、D201/401/801数据、原图与修复图均保留；D未重算。随后仅只读解释已保存CSV并补全本文，不算新声学执行或失败重试。

## 6. 证据与当前边界

[结果目录](../results/validation/bellhop_1500m_applicability/20261119_S4000_C801_01/)、[S及D+S指标](../results/validation/bellhop_1500m_applicability/20261119_S4000_C801_01/returned/S_and_DS_metrics.csv)、[原始五点复压](../results/validation/bellhop_1500m_applicability/20261119_S4000_C801_01/returned/S_and_DS_fields.csv)、[反射几何检查](../results/validation/bellhop_1500m_applicability/20261119_S4000_C801_01/returned/surface_geometry_checks.csv)、[全部801条第一事件](../results/validation/bellhop_1500m_applicability/20261119_S4000_C801_01/returned/all801_first_surface_states.csv)、[原始MAT](../results/validation/bellhop_1500m_applicability/20261119_S4000_C801_01/returned/comparison.mat)、[原退出与监督回执](../results/validation/bellhop_1500m_applicability/20261119_S4000_C801_01/supervisor_receipt.json)。

本轮完成历史图件修复及唯一S数值计算；新增S图件审计失败后停止，不把失败回执抹除。没有B/BS、其他P1、P2–P5、FFP/Munk；没有生产代码/公共接口修改，不需要通信生产回归。没有提交、推送或合并。等待人工审查。
进度视图已通过恢复的本地既有MATLAB MCP刷新（仅文档渲染，无声学计算）。回执保留101项引用检查中的16项历史引用缺失告警（其中1项为旧Gaussian Bellhop checks.csv，其余来自其他工作流），未恢复或搜索归档；本轮S结果与报告引用未新增缺失。最终只读复核：375项历史SHA、79项返回SHA均无变化。详见final_integrity_review.json。
