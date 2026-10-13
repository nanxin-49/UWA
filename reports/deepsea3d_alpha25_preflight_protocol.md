# alpha2.5 低成本参考预检协议

仅原dual振幅(.02,.01) m乘2.5，固定Tx1500/Rx5、q0、6 kHz、c1500及全部谱/phase/C2约定。
冻结FE及其他求解器不修改；旧alpha2结果只读取作同配置比较，不重算。
L24 m、h=.125 m（384²点），先N12/N16、M40、真实海面配点及半网格残差。
若12→16增量变化、N16最后两阶变化或任一残差未过原1e-4，则只补N20。
最高20阶；仍不足标记INCONCLUSIVE，不加阶或放宽门。完整G1状态仍INCONCLUSIVE，
即使预检通过也不做本轮的正式网格/域加密或FK比较。

保留原生逐阶复系数、部分和及边界系数范数，检查峰值后衰减、最后两阶及整阶差。
pressure取消条件指数为sum(abs(terms))/abs(sum(terms))，只诊断接收点求和，
不是内部递推或连续PDE条件数。另以低成本诊断回放观测每阶递推贡献范数和/结果范数；
回放必须逐值核对冻结FE的系数范数和原生接收系数，不用它替代物理参考。
记录wall/CPU、谱载荷及MATLAB进程前后内存快照；快照/载荷不冒充进程峰值。

决策仅为PROCEED_FORMAL_G1或INCONCLUSIVE：各收敛/边界量<1e-4且系数未显示
高阶失稳，才建议正式验证。若只是N12不足但N20稳定，不称FE失稳。
若尾部增长、半网格残差或消项使参考不足，停止并评估TFE/其他参考。
全G1目标域OPEN，WP2=false；1%FK物理门未被检验，无FK物理PASS/FAIL结论。

复现：setup后`validate_deepsea3d_alpha25_preflight('init')`，重复`'next'`至complete，
然后`'cancellation'`诊断回放、`'audit'`审核。init不覆盖已有状态，next每次一个参考运行。
结果`results/validation/deepsea3d_alpha25_preflight/`；成本中区分冻结参考和诊断回放。
