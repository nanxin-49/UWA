# 本轮执行记录

现有MATLAB R2025b Update6，cwd始终E:/MISC/CARPE3D_matlab/Explain；未调用任何声学核、FFT源评价、保存密度新点评价或新案例求解。分析运行analysis_IO_seconds=5.0875977，不含全部图形/I-O/审计，不伪报Goal总墙钟。核成本为0，累计1408.660105、余额5791.339895 s。

run_two_dimensional_closeout_s2c_a在第一次1630尝试中发现汇总筛选遗漏reference=BIE条件（同candidate有两种reference），未产出正式指标/图；修正后1642运行成功。失败目录原样收纳execution_attempts/failed_1630，只作为过程记录，不是第二正式run。

图件修订依次纠正sprintf转义、显式对数轴、九点跨案例统一色标、全支持经验RMS参考虚线、x=0辅助线误入legend，以及相位跨界只插NaN分隔符、不构造新相位样本。初始/中间图均保留在execution_attempts，正式仅figures。

初轮审计通过基本CSV/MAT与恒等式；补充基线几何/照射核验时，单行来源CSV被readtable误猜分隔符，改显式逗号。φ的CSV十五位文本舍入导致与JSON精确浮点相等判据过严，按1e-14文本往返容差修正，不修改原声学/几何数组。每个失败类型仅一次相应修正，无反复放宽物理误差门。最终audit/全部通过；没有新物理作业或数值模型调参。

实际查看正式全部11PNG，坐标、单位、wrapped相位断点、近零0、弱中FullFK九点UNRESOLVED、真实9点/411线/局部11点、图例及对数尺度已检查。11FIG重开；保存MAT重绘PNG逐像素最大差0，图册22相对图件链接存在且非空。三个新MATLAB文件Code Analyzer零消息；git diff --check通过，仅有Windows换行提示。

生产求解器/公共API未改，不执行范围外通信回归；无cd、shell MATLAB、cash访问、Bellhop、S2c-B、S3–S6、提交或推送。图册在Codex打开请求已排队；中文Markdown相对PNG/FIG链接可直接预览/打开。