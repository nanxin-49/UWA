# K1 风速—物理谱截止联合声场敏感性

N7从已审查旧K1任务导入；其余五面在本恢复任务独立完成原PS/BIE/FullFK层级。U6/N12为严格匹配的历史R1锚点。

固定0.10m二维条带、4kHz、seed260001、原C2、phase contract2.0.0、411真实接收点。每个海面使用自己的BIE范数与0.01全线RMS相位阈值。N40不表示完整PM谱已收敛。

| U_PM (m/s) | N | 模型 | 复场L2 (%) | 5倍经验阈值 (%) | 可分辨 |
|---:|---:|---|---:|---:|---|
|6|12|PS|15.923498|0.00018957684|1|
|6|12|FullFK|0.074720825|0.00015644143|1|
|6|7|PS|10.25099|0.00011405036|1|
|6|7|FullFK|0.020832844|0.00010261058|1|
|6|40|PS|26.009111|0.0011602305|1|
|6|40|FullFK|0.18955118|0.00028342228|1|
|5|12|PS|11.383285|0.00012099949|1|
|5|12|FullFK|0.063221441|0.00011371724|1|
|8|12|PS|42.251728|0.00034246637|1|
|8|12|FullFK|0.086026444|0.00028106169|1|
|5|40|PS|19.439822|0.0012602482|1|
|5|40|FullFK|0.18023772|0.00023026395|1|
|8|40|PS|53.845146|0.0008244742|1|
|8|40|FullFK|0.19553116|0.00040518593|1|

可分辨性是经验数值分辨，不是严格真解或工程合格线。未分辨比较如实保留，不自动调网格/门槛。七个条件为同一随机块，411点不当独立重复，不外推总体显著性、宽带/CIR或通信。原始/C2几何与谱差异保存在每案例input/k0_*文件，不把它们当相同带限面。

![PS absolute](figures/K1_PS_absolute.png)

![PS normalized](figures/K1_PS_normalized.png)

![PS amplitude](figures/K1_PS_amplitude.png)

![PS phase](figures/K1_PS_phase.png)

![FullFK absolute](figures/K1_FullFK_absolute.png)

![FullFK normalized](figures/K1_FullFK_normalized.png)

![FullFK amplitude](figures/K1_FullFK_amplitude.png)

![FullFK phase](figures/K1_FullFK_phase.png)

资格由原MODELS/C0–C4封存链、本地独立保存数值审计及新固定画布图件零像素审计组成。旧大画布图件FAIL及原远端退出回执保留，未伪造旧FINAL封存；见local-closeout资格清单。
