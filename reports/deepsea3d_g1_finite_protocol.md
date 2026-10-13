# G1 单案例有限振幅参考协议（2026-10-12）

沿用现有WP0约定和原FullFK；不重跑解析、一阶或旧网格验证。
唯一母案例：Tx(0,0,1500)、Rx(0,0,5)、c1500、6 kHz，原dual同一C2面，
仅把振幅固定为(.02,.01) m；波矢(.4,0)/(.2,.6)、phase(.2,.7)、plateau8/support12。
理论高度界.03 m，k*height界约.754，二阶及以上效应预期不可忽略。
这是一个有限振幅紧支撑二维粗糙面，不是PM海况。

可行性比较：朴素密集3D BIE不启动；half-space/WGF/BEM需要奇异求积、无限边界
处理和新实现。高阶Dirichlet FE不采用Kirchhoff局部边界假设，成本低但可能不稳定。
先运行`validate_deepsea3d_g1_preflight`（L24，h=.25/.125，N6），据自身网格、
真实边界残差和耗时决定是否继续。先前一阶参考不作为有限振幅oracle。

FE采用t*eta同伦，增量u_delta的平面边界系数满足
`b_n=-[t^n]u_flat(t eta)-sum_(j=1)^(n-1)eta^j/j!*D_z^j b_(n-j)`，
`D_z=ifft2(i*kz*fft2)`，kz正实/正虚分支，全部倏逝模态保留，不滤波。
入射/镜像源Taylor系数用解析sqrt、倒数、exp幂级数计算，不用有限差分或FK核。
正/负eta区域都检查实际粗糙面Dirichlet残差：用独立更高M阶传播指数求值，
不只检查自动闭合的N阶多项式。比较M24/32控制边界求值尾项。

预检可行后：N8/10/12/14阶数、h=.125/.0625网格、L24/32/48域扩展，
必要时补充一次更紧配置，保存全部失败，不为通过选择性截断模态或更换海面。
参考数值目标：阶数/网格/域/传播求值变化相对于增量分别<1e-4，真实边界残差<1e-4。
各变化之和只作经验收敛指示，不当严格数学余项界。
若出现不收敛或FE失稳，应保留OPEN，转向TFE/匹配BIE，不启动大批量计算。

待可信参考锁定后，在同一场、频率、海面原生exp(-i omega t)下比较FullFK。
分别报告P_ref及DeltaP=P_ref-P_flat的复误差、幅度误差、原始phase差；无对齐/拟合。
单频条件物理筛查采用1%反射和1%非零增量作为本次新增规则，phase差另列；
它不是已完成band-relative G1，不能自动批准WP2。参考不确定性应明显小于1%。
有限振幅效应比较一阶项与完整高阶和，不重复既有弱极限实验。
同一固定海面附近的一个局部相消频率可作附加敏感性诊断；另做参考自身收敛，
不将单点相消诊断当宽带研究或全球最差点证明。

结果放`results/validation/deepsea3d_g1_finite/`，入口及中文报告在validation/reports。
生产、二维基准、原三维FK/phase/surface helpers保持哈希；不进行被排除研究。

复现（不依赖MATLAB当前目录）：

```matlab
info=setup_vertical_project;
validate_deepsea3d_g1_preflight;
validate_deepsea3d_g1_finite('init');
for job=1:10, validate_deepsea3d_g1_finite('next'); end
validate_deepsea3d_g1_finite('offgrid'); % 独立半网格真实边界残差
validate_deepsea3d_g1_finite('compare');
validate_deepsea3d_g1_finite('cancellation');
tests=runtests('test_deepsea3d_finite_reference');
assert(all([tests.Passed]));
writetable(table(tests),fullfile(info.project_root,'results','validation', ...
    'deepsea3d_g1_finite','tests.csv'));
validate_deepsea3d_g1_finite('audit');
```

`init`重建本case状态；`next`每次只做一个待完成参考，保存checkpoint，可逐次续跑。
cost_summary/provenance记录本次实际保存运行；若重跑，计时不同，不能套用旧hash或成本。
