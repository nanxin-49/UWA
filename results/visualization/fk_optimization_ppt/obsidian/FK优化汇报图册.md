# FK 优化：汇报可视化

日期：2026-10-08。仅整理已有验证证据，不推进主线研究。

> 使用方式：将本目录整体复制到 Obsidian vault，打开本报告；请保留同级“图片”文件夹。图册按汇报顺序排列，每图下方先给出一句中文概括，再解释图中内容、主要发现及其含义。

面向组会/PPT，用 11 张白底宽屏图讲清“物理图像 → 算法精度与速度 → 物理适用域 → 跨频一致性”。每图独立可用；术语统一为 Incident、Reflected、FullFK、FastFK、BIE。本 Obsidian 副本附带全部 11 张 PNG；同名矢量 PDF 保留在源项目中。

## 数据来源（路径均相对于源项目）

- Prototype（源项目：`reports/spectral_fast_kirchhoff_prototype_report.md`）：`results/validation/spectral_fast_kirchhoff_prototype/` 的 `all_metrics.csv`、`scaling.csv`、`weak_fields.mat`、`strong_height_fields.mat`。
- PM pilot（源项目：`reports/spectral_fast_kirchhoff_pm_pilot_report.md`）：`results/validation/spectral_fast_kirchhoff_pm_pilot/` 的 `field_metrics.csv`、`seed_260001_fields.mat`。
- 物理边界（源项目：`reports/full_kirchhoff_pm_validity_boundary_report.md`）：`results/validation/full_kirchhoff_pm_validity_boundary/` 的 `validity_metrics.csv`、`case_data/cutoff_n12_receiver_fields.csv`。
- 跨频一致性（源项目：`reports/spectral_fast_kirchhoff_multifrequency_consistency_report.md`）：`results/validation/spectral_fast_kirchhoff_multifrequency_consistency/` 的 `dense_metrics.csv`、`bie_anchors.csv`、`phase_group_delay.csv`、主 `reflected_cir.csv`；相位语义以 vertical_comm_guide（源项目：`vertical_comm_guide.md`） 为准。

没有重测速度，没有运行 PE、FK、BIE 或随机 realization。弱海面 Incident 只按原 prototype 的离散角谱公式从保存系数重建；PM Incident 直接读取保存的 exact input。重建剖面保存在源项目可视化输出根目录的两个 `*_incident_source.csv`。`visualization_audit.json` 记录关键数值，`source_manifest.csv` 记录输入 SHA256。

## 图册

阅读时先区分两种比较：**FastFK–FullFK** 检查快速算法是否保留原算法的复数场，属于算法保真；**FullFK–BIE** 检查 FK 物理近似与独立波动方程参考是否一致，属于物理有效性。复数场同时包含幅度和相位，幅度曲线看起来重合并不意味着相位也完全一致。下文的 rank 写作“入射阶数/出射阶数”，表示压缩时保留的基函数数量；热启动速度指准备好计算所需的基与核之后，重复执行场计算的速度。

### 图 1｜项目路线与成果总览

![[图片/ppt_overview_roadmap.png]]

工作链已完成从 FullFK 验证到 FastFK 压缩和固定海面的跨频一致性检查。

上方五个框按研究顺序展开：先确定上行传播到海面的入射场，再以 FullFK 汇总海面各位置的相干反射贡献，随后用 FastFK 压缩这一步计算，最后分别检查 PM 随机海面下的物理适用域和跨频一致性。BIE 是独立的波动方程参考，用来检查 FullFK 的物理近似，不能由 FullFK 与 FastFK 相互重合来替代。

下方三个数字框概括已有成果：确定性 prototype 在 N=4097 时取得约 141 倍热启动加速，PM pilot 的统一 3/8 配置取得约 89–112 倍，固定 canonical 海面的 4–8 kHz 验证采用统一 3/12，最大加权复误差为 2.31×10⁻⁴。这一页说明，当前已得到经过检验的快速压缩实现，同时保留了物理适用条件；可进入小规模 realization 数据阶段是后续允许的范围，不表示本图册已生成这些数据或完成生产集成。

### 图 2｜典型海面与传播几何

![[图片/ppt_geometry_setup.png]]

声波从 Tx 上行至粗糙海面，反射后到达 Rx 平面，箭头仅为传播示意。

左图说明几何关系：深度 z 正向向下，Tx 位于 z=100 m，声波向上传播并照射粗糙海面，反射场再向下到达 z=3 m 的接收平面。蓝线表示海面，紫色虚线表示 Rx 平面。由于 100 m 的深度尺度远大于海面起伏，左图中粗糙度看起来很小，因此右图单独放大了海面和接收平面之间的区域。

右图展示保存的 canonical PM 海面，而非任意画出的波浪；两条竖向虚线标出约 ±39.84 m 的 M99 核心范围，即本轮固定的 99% 入射能量足迹所定义的评价区域。海面在中心区域保留 PM 起伏，在边缘用 C2 包络平滑嵌入有限计算支撑。此图帮助理解“在哪里发生散射、在哪里比较接收场”，但箭头不是计算所得的声线路径，也不能由 Tx 符号推断存在海底反射边界；这组 Tx100/Rx3 研究夹具与 Tx99 的 Bellhop 安装图不同。

### 图 3｜到达海面的入射场

![[图片/ppt_incident_field.png]]

海面处的入射幅度与相位清晰可见，FastFK 压缩输入与精确输入基本重合。

两列分别对应弱确定性海面和 canonical PM 海面，上排显示到达真实海面位置的入射幅度，下排显示展开相位。横轴都是海面横向坐标 x；各列幅度以本案例的精确入射峰值归一化。幅度在中心较强、向两侧减弱，体现照射能量集中于有限区域；PM 案例的相位曲线又叠加了真实海面高度变化引起的结构，不能只根据幅度的平滑外形把入射场视为等相位平面波。

蓝色实线为 FullFK/BIE 使用的精确离散角谱输入，橙色虚线为保存的 FastFK 入射 rank3 近似，二者在幅度和相位上基本重合。全网格入射压缩复误差分别约为 3.67×10⁻⁹ 和 2.88×10⁻⁴，说明这两个案例的低阶入射表示保留了主要场结构。这里检查的是输入压缩，不能把它称为新增的 PE 独立传播验证；共享输入也使后续 FullFK–BIE 比较能够集中检查海面散射算子的差异。

### 图 4｜反射场与接收场对比

![[图片/ppt_reflected_field_comparison.png]]

FullFK 与 FastFK 接收场基本重合，canonical PM 的 BIE 残差用于独立物理检查。

三列依次展示弱确定性海面、强高度海面和 canonical PM 海面的反射接收场。上排为 z=3 m 接收平面上的幅度剖面：蓝色 FullFK 与橙色 FastFK 在三个案例中基本重合，PM 列还加入紫色 BIE。每列均以本列 FullFK 核心峰值为共同幅度基准，因而同一列的算法可以直接比较；不同列分别归一化，不能仅凭峰值都接近 1 推断其绝对声压相同。

下排把不易从幅度中看出的相位差单独展开，单位为毫弧度。橙线表示 FastFK–FullFK 的算法残差，PM 列的紫线表示 FullFK–BIE 的物理残差；各列纵轴尺度不同，读图时应结合刻度。对应 FastFK 加权复误差约为 2.90×10⁻⁹、3.21×10⁻⁵ 和 2.48×10⁻⁴。它们说明 FastFK 在这些案例中保留了接收场的主要复数结构，而 BIE 比较另行支持 canonical PM 在已测核心区域内的物理有效性；这不等于整条接收线或全部 PM 海面都已通过。

### 图 5｜反射前后的场变化

![[图片/ppt_incident_vs_reflected.png]]

同一 PM 海面使平滑入射场在反射后形成更有结构的接收场分布。

上方给出 canonical PM 海面在核心区域的实际起伏，下方左图显示入射到这一真实海面上的幅度，右图显示反射后到达 Rx 平面的幅度。两幅场图共用“精确 Incident 峰值”这一归一化基准，所以可以同时观察反射前后的幅度和空间形态变化，而不只是比较各自归一化后的形状。

入射幅度较为平滑，反射接收场出现多个局部峰谷，体现不同海面位置贡献在接收平面上的相干叠加及空间重分布；右图中 FullFK 和 FastFK 仍基本重合。右侧局部峰值略高于入射归一化峰值并不自动意味着能量放大，因为比较对象位于不同几何位置，且这里没有对整个声场的能流作积分。此图用于建立散射的直观图像，不能解释成二维传播全场图或能量守恒验收。

### 图 6｜Prototype 精度与加速

![[图片/ppt_prototype_accuracy_speed.png]]

确定性案例兼具高精度与加速，N=4097 时热启动速度收益约为 141 倍。

左面板汇总 weak、strong-height、high-K、K070、sigma1 和 Rx2 六个确定性控制，采用保留完整几何的 Prototype-B、统一 3/6 rank 和 padding=2。纵轴是对数尺度的加权复 L2；六个点均低于 10⁻⁴，说明压缩在已测高度、横向粗糙变化、源宽度与接收深度控制中保持精度。中间面板固定入射 rank=4，只改变出射 rank，误差随保留阶数增加而下降，说明出射核的压缩精度需要单独检查，不能仅凭入射场低秩就宣布完整散射算子低秩。

右面板比较已有网格规模下的重复计算耗时：随着表面点数和接收点数 N 增加，FullFK 的计算时间增长更快，FastFK 的 FFT 路径保持较低耗时。N=4097 时热启动约加速 141 倍，计入准备成本后的冷启动约为 14.7 倍，因此重复使用准备结果时收益更明显。这些是当时硬件与会话中的有限规模测量，不能直接外推为渐近复杂度证明或端到端生产系统的加速倍数。

### 图 7｜PM pilot：算法与物理分开看

![[图片/ppt_pm_pilot_algorithm_vs_physics.png]]

统一 3/8 rank 通过算法保真门，但高截止案例仍未通过 FullFK–BIE 物理门。

左面板检查五个固定 PM 输入的算法保真，包括三个系数 realization 和低、高截止控制。冻结的 3/6 配置误差约为 4×10⁻³，超过 10⁻³ 门槛；采用对五个案例共享的 3/8 后，误差降至约 1.9–2.5×10⁻⁴并全部通过。这说明确定性阶段的配置不能直接照搬到 PM 阶段，但通过增加统一的出射 rank，仍能有效压缩同一个 FullFK 算子。

右面板另行以 BIE 检查 FullFK 的物理近似：canonical 的加权复误差约为 8.00×10⁻⁴，低于门槛，高截止案例约为 1.20×10⁻³，高于门槛。两面板的参照对象不同：左图通过回答“算得是否像 FullFK”，右图回答“FullFK 是否足够接近独立波动参考”。因此 pilot 的整体物理门仍未通过，继续提高 FastFK rank 不能消除原 FullFK 模型的物理近似误差。

### 图 8｜PM 条件适用域

![[图片/ppt_pm_validity_boundary.png]]

固定条件下 Kmax=0.706858 通过、0.746128 失败，适用域还取决于幅度与受照结构。

左面板沿同一 master Fourier 系数族逐步增加实际物理截止波数 Kmax，横轴不是数值网格的 Nyquist 上限。曲线显示误差并非处处随截止值单调增加；经过数值护栏检查，Kmax=0.706858 rad/m 的误差约为 8.53×10⁻⁴，而 0.746128 rad/m 约为 1.20×10⁻³，分列 10⁻³ 门槛两侧。端点的小误差棒表示网格、BIE 精化和窗口变化量之和，约 6.8×10⁻⁶，明显小于两端距门槛的差距，支持这两个端点的通过/失败判断；它不是统计置信区间。

右面板把全部 15 个受控输入的物理误差与“受照区域 RMS 曲率/k”并列观察，并突出保持同一 Kmax、只改变新增高波数模态幅度的案例。即使 Kmax 相同，误差也会随幅度和海面结构变化而落在门槛两侧，因此仅用 Kmax 或高度不能给出统一适用规则。受照曲率是有用的候选描述量，但图中没有拟合普适边界，也没有唯一识别误差的因果参数；本结论限于这一个 realization、4 kHz、Tx100/Rx3 和固定核心孔径。

### 图 9｜4–8 kHz 跨频一致性

![[图片/ppt_multifrequency_consistency.png]]

固定 canonical 海面的 4–8 kHz 一致性通过，但 16 ms 周期窗 CIR 仍仅作诊断。

左上面板展示 65 个频点的 FastFK–FullFK 加权复误差，整个 4–8 kHz 频带采用同一 3/12 rank，最大值约为 2.31×10⁻⁴，低于 10⁻³ 门槛。右上面板是 4、6、8 kHz 的独立 FullFK–BIE 锚点，误差约为 8.00×10⁻⁴、5.30×10⁻⁴、4.08×10⁻⁴，说明三个已测锚点的核心场也满足物理门。它们是三处物理参考检查，不应写成对全部 65 点逐一完成 BIE 验证。

左下直接画出 FullFK 与 FastFK 的相对相位差，最大约为 1.74×10⁻⁴ rad；标题同时给出由保存相位数据检查得到的最大群时延差约 70.2 ns，说明压缩对该固定海面的频率相位变化影响很小。右下是在明确 native-to-DSP 符号桥接后得到的 reflected-only CIR，两条幅度曲线基本重合，复 L2 约为 7.17×10⁻⁵。这个结果支持有限频带内两种算法的响应一致性，但 65 点 IFFT 对应 16 ms 周期窗，尚未证明全部长路径尾部无混叠；约 3.69 ms 的采样峰也不等同于约 3.595 ms 的相位群时延。

### 图 10｜统一 rank 的阶段演化

![[图片/ppt_rank_evolution.png]]

统一 rank 从 3/6 增至 3/8 和 3/12，按阶段共享配置而非逐案例或逐频调参。

横轴按三个验证阶段排列，蓝线表示入射 rank 始终取 3，橙线表示出射 rank 从 6 增至 8，再增至 12。入射与出射分别控制不同部分的分离近似，图中强调的是随已测范围扩展而调整共享配置，不是把两个数字相乘或相加后得到某个“完整核 rank”。

下方注释解释调整的证据：3/6 在 PM pilot 不够准确，因此五个 PM 输入共同改用 3/8；扩展至 4–8 kHz 后，3/8 从 6 kHz 起失败，3/10 在 8 kHz 仍略超门槛，有限阶梯中的下一配置 3/12 才通过整个频带。每阶段的失败试验均保留，没有逐案例或逐频选择最优阶数。此图说明配置选择有可追溯的验收依据，但不宣称这些 rank 已达到全局最小，也不保证它们适用于未测海况或频带。

### 图 11｜成果与限制总结

![[图片/ppt_summary_dashboard.png]]

四阶段证据支持条件适用域内的 FastFK 路线，尚不能外推为生产或总体海况结论。

表格把四阶段的目标、代表误差、速度收益和当前判断放在一起。Prototype 回答确定性 FullFK 是否可被准确加速；PM pilot 回答压缩能否延续到五个固定随机海面输入；PM boundary 单独确定 FullFK 的条件物理适用域；multi-frequency 检查同一 canonical 海面跨频后的复数响应与相位一致性。各行误差的参照对象不同，前三个算法配置与物理边界行不能混为一个统一的误差排名。

总体成果是：已测输入上的 FastFK 压缩精度和计算收益成立，且固定海面的跨频一致性通过；仍需保留的限制是，高截止 PM 控制存在物理门失败，适用域属于单 realization 条件结论，CIR 为周期窗诊断。表格下方提出的小规模 realization 数据集是接下来可以单独开展的工作，并未在本图册中执行。汇报结束时可据此明确区分“当前已经建立的证据”与“下一阶段准备做的事”。

推荐完整顺序：**1 → 2 → 3 → 4 → 5 → 6 → 7 → 8 → 9 → 10 → 11**。若限制为十页，可将图 5 作为备选。

若只保留六张：**1、3、4、6、8、9**。依次讲全局、入射、反射、精度/速度、物理边界、宽带一致性；图 8 已包含物理误差，讲述时口头强调它不同于 FastFK–FullFK。

## 阅读边界与归一化

- 本套几何和场图属于二维 Tx100/Rx3 研究夹具；不代表 Tx99 的 Bellhop 安装验证，底部 Tx 符号也不表示加入了海底边界模型。传播箭头只是示意。
- 图 3 显示 exact FullFK/BIE 共享输入与 FastFK rank3 输入；**不是新增的 PE 独立一致性验收**。不把压缩后的入射场冒称 exact 输入。
- 图 3 各案例幅度以自身 exact Incident 峰值归一；图 4 每列共用 FullFK 核心峰值；图 5 前后共用 Incident 峰值；图 9 CIR 共用 FullFK CIR 峰值。没有幅度拟合、相位对齐或峰位移动。
- 场图显示约 ±40 m 的核心孔径，省略的外部点仍在原始文件中。相位差沿用保存 M99 / Full-only −40 dB 掩膜；原相位剖面使用 unwrap，只为展示，不进行定相位。
- 所有 E 为相对加权复 L2；算法误差以 FullFK 为参照，物理误差以 BIE 为参照。rank 图的 outgoing 阶梯固定 incoming=4，不能与最终选定 3/6 混为同一实验。
- 图 8 端点误差棒为三个数值变化量之和，约 6.8e-6，**不是统计置信区间**。其余多数点是 baseline-only；连接线仅指示已采样 continuation 顺序，不给出插值定理、拟合公式或通用阈值。
- 4–8 kHz 一致性仅针对固定 canonical seed260001/C2 海面、65 点和统一 3/12。三个 FullFK–BIE 锚点的物理结论限于冻结 M99 核心。
- CIR 使用已验证 native-to-DSP 全包络符号桥接和反射分量；65 点矩形 IFFT 的 16 ms 周期窗**未认证所有长尾无混叠**。图中约 3.69 ms 的采样峰不等同于约 3.595 ms 的相位群时延。
- 141x 和 89–112x 是已有 kernel warm wall-time 比值，不能等同于端到端生产加速；跨频最低 warm/cold 为 41.94x/9.83x。所有时间依赖当时硬件/会话。

## 复现与检查

`setup_vertical_project; generate_fk_optimization_ppt_visuals();`

入口：scripts/reporting/generate_fk_optimization_ppt_visuals.m（源项目：`scripts/reporting/generate_fk_optimization_ppt_visuals.m`）。自动解析根路径，可从任意 MATLAB cwd 调用；默认重绘本报告输出，不覆盖验证数据。白底、16:9、大字号，PNG 和矢量 PDF 成对导出。Code Analyzer、导出文件尺寸、关键指标复核和人工版面检查见源项目输出 `results/visualization/fk_optimization_ppt/qa_summary.md`。

仅新增报告/可视化入口及导航，生产模型、接口和方法公式均未改变；因此不需要标量信道/宽带通信回归，也不修改 `vertical_comm_guide.md` 的物理方法正文。`PROJECT_CONTEXT.md` 只追加可视化导航与主线暂停说明。

工作流使用 scientific-visualization 与 matlab 技能；归因：Kassis T., Agarwal V., He Y., Patel D., Brueckner A. M. (2026), [Scientific Agent Skills](https://doi.org/10.48550/arXiv.2609.00065)（当前 arXiv 记录已核对）。此引用为工作流来源，不是声学结论依据。
