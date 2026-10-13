# alpha2.5细网格参考验证预设协议

只沿用冻结alpha2.5 cfg（Tx1500/Rx5、q0、6 kHz、dual/C2），L24，h=.0625。
N16/N20、M40，均检查原网格与半网格真实边界；不得超过20阶，保持全部倏逝模态。
旧h=.125结果只读，不重复求解。不得调用FullFK、扩域、TFE/BIE或WP2。

原参考门保持1e-4：fine N16→20、N20最后两阶、匹配N16/N20的coarse→fine接收
增量变化、最终配点和半网格真实边界残差各自检查。点压力收敛不得代替边界系数
稳定；最后五阶压力幅度/边界系数范数须逐阶检查，若增长须明确INCONCLUSIVE。
相同L下系数范数按网格步长h转换为近似连续L2再比较，避免样本数翻倍被误当增长。
M40有限Fourier指数评价尾项与FE截断误差分开；没有网格以外或PDE严格余项界。

对每个fine N诊断回放内部贡献范数和/结果范数，断言逐阶压力及系数范数与冻结FE一致。
指标只描述消项，不是数学条件数；同时看高阶尾部及真实边界残差，不以单指标放行。
记录参考与观测回放wall/CPU、谱载荷和MATLAB进程前后内存快照，均不冒充峰值内存。

决策只有FINE_REFERENCE_CHECK_PASS或INCONCLUSIVE；alpha2.5物理G1始终INCONCLUSIVE，
全目标域OPEN，WP2=false。通过后可建议下一步域扩展/物理比较，但本轮不执行。
init→两次next→两次observe→audit，结果`results/validation/deepsea3d_alpha25_fine/`。
