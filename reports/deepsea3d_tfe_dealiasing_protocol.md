# 三维TFE混叠与展开失效诊断协议

固定冻结PM realization/seed261013、6kHz、q=.10、Hm0=.5父面及C2 support12；
读取prototype.mat，不生成海面，不改旧求解器。新验证入口validate_deepsea3d_tfe_dealiasing(job)。

1. L1含二因子，L2含三因子（eta²H、ex²Wss、eta ex Wxs、eta lapeta Ws）。
   用同一基网格采样的eta/ex/ey/lapeta Fourier插值作配对系数，仅显式去掉两条Nyquist线，
   场保留对称|j|<n/2。两分支使用相同源项、系数和模态空间，区别仅乘积投影及DtN乘积。
2. 3/2填充不能一般消除三因子混叠；2n网格同时计算全部乘积，最后一次投影，
   不在二次子乘积后先投影。制造模态检查FFT归一化、复模态、Nyquist及三因子低频卷积。
   基网格系数只是连续C2几何的数值插值；物理谱和窗口不变，仍须检验系数截断。
3. 在每个N检查点，用独立3n/4n网格重建同一场，重新评价真实连续几何、
   原生真实边界、物理Helmholtz及DtN；保留未投影残差与投影残差，不能用Galerkin正交掩盖高频残差。
   同时比较系数插值的精确高分辨率评价与collocation，隔离混叠及几何系数截断。
4. 先dual_fg与quarter_fg（N16/h.25/M32/L24/b.8），控制通过后原PM配对
   pm_collocation与pm_fg，单次递推到N32并保存N16/24/32，不重复求解三个阶数。
   必要时仅一个half_fg振幅诊断，不认证缩放面。若原目标仍明显失败则不增加N。
5. 禁止support24、FullFK–PM、Monte Carlo和WP2。仅当去混叠有显著收益且尾部递减，
   才考虑一次既有h.125网格的pm_fine_fg；必须通过可用物理内存检查，不靠交换区扩大规模。
   新粗网格阶段总核预算600s，所有新核/独立评价总预算900s；每job保存成本，超预算不再安排新job。
   N≤32，L=24，M32，b=.8；命名工作数组和快照记录，非峰值。
6. 门保持1e-4，ODE残差1e-8；所有阶数/残差/真实边界/透明边界/尾部与控制比较同列。
   阴性结果为PM INCONCLUSIVE、G1 OPEN、WP2=false，不以本轮获取G1为目标。

结果目录results/validation/deepsea3d_tfe_dealiasing/，job拒绝覆盖；
MAT保留N16/24/32完整场，CSV保存压力/残差/尾部，可独立重新评价。

## 对照分支实现纠正

初始pm_collocation体右端未投影掉新生成的Nyquist线，而Galerkin分支已投影；
不能将该初始分支作为同空间配对。原始MAT/CSV原样保留并标记无效配对，
追加pm_collocation_fixed，仍N32且复用16/24/32检查点，不重复Galerkin或扩阶数。
新制造模态断言collocation体右端也满足声明的对称模态空间。
正式收益只由fixed与pm_fg比较，初始错误运行也计入实际成本。

## 停止决策与最终复现

修正配对后原PM残差及高阶接收尾部仍不过门，停止在N32；不运行half_fg或pm_fine_fg。
命名数组上限900MiB，粗网格估计384.75MiB；细网格及4倍评价1539MiB超上限，不分配。
有效复现顺序dual_fg→quarter_fg→pm_fg→pm_collocation_fixed；
probe_deepsea3d_fg_geometry重评保存场，audit_deepsea3d_fg_artifacts核对结果。
初始pm_collocation已退役为只读数据，最终代码不支持重跑旧错误分支。
