# M05：PHI_PLUS90 几何图册与声场缺口

预算已追加到 7506.601 s，但当前真实物理 RAM 未过门，用户要求暂时保持资源阻塞。下方 M05a 仍是正式几何图；声学 M05b–e 尚无真实 PHI 复场。预算和 RAM 数值见[执行补记](execution.md)，M05e 将来使用的坐标对称性检查见[输入审计](P1_symmetry_input_audit.md)。

![M05a：BASE_M 过零与 PHI_PLUS90 波峰的高度、坡度、曲率](figures_final/M05a_geometry_only.png)

[可编辑 FIG](figures_final/M05a_geometry_only.fig) · [解析几何 CSV](figures_final/geometry.csv) · [冻结协议和预算门](protocol.md)

## M05a：完整 C2 几何

**看到了什么：** BASE_M 在 x=0 过零，PHI_PLUS90 在 x=0 达 +0.05 m 波峰；中心坡度由 .005 变为 0，曲率由 0 变为 −.0005 m⁻¹。灰线为 `|x|=42/50 m` 的平台/窗口边界；中心接收机投影仅是 x=0，不表示唯一反射点。含窗口全域最大斜率和曲率见[汇总 CSV](figures_final/geometry_summary.csv)，不能用理想平台值替代。

**这些数据支持什么：** 只改变了输入初相位 φ，外侧 `|x|≥50 m` 保持严格平坦，解析几何已核验。这为未来完整 FullFK 解析平面外场补全保留几何前提。

**目前还不能证明什么：** 没有新 PHI 声场，中心 PE–BIE 或 FullFK–BIE 残差的方向、大小、可分辨性均未知。φ 同时改变整段海面，高度、坡度、曲率及传播干涉不能单独归因。

## 已合格 K_HIGH 图像，作为机制背景

下列四图原样复制到[轻量发布包](publication_review_M04/README.md)，只对应 `K_HIGH=(.05,.2,0)` 与 BASE_M，不代表 PHI 场：

![M04b：K_HIGH 的 BIE 参考物理响应](publication_review_M04/M04b_BIE_physical_response.png)

![M04c：K_HIGH 的 PE/FullFK 模型残差](publication_review_M04/M04c_model_error_pair.png)

![M04d：整线、中心与真实九点](publication_review_M04/M04d_line_center_nine.png)

**看到了什么：** K_HIGH 改变 BIE 参考场和模型残差；BASE_M 中心低误差与整线/边缘差异明显。详细固定支持指标见[原图册](../20261023_s2c_b_khigh_02/mechanism_atlas.md)。

**这些数据支持什么：** 已完成的 A、K 受控输入案例说明中心现象与空间支持、幅相贡献相关：保存逐点 CSV 复算的三档 A 情形 PE 整线相位平方占比约 `99.972–99.973%`，K_HIGH 为 `99.961%`。这是复差的代数分解，值得用初相位对照继续检查；这些位置不是独立海况样本。

**目前还不能证明什么：** A、K 和未来 φ 分别改变多项几何描述；现有 K_HIGH 场不能预测 PHI 的模型误差，更不能作局部过零的因果判定。M05b–e 和三参数声场总览因 PHI 参考门失败而未生成。

图形工作流参考：Kassis, T., Agarwal, V., He, Y., Patel, D., & Brueckner, A. M. (2026), *Scientific Agent Skills: A Library of Procedural Knowledge for Research Agents*, [arXiv:2609.00065](https://doi.org/10.48550/arXiv.2609.00065)。此文仅用于程序性绘图/审计方法，不作为声学证据。
