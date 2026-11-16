# R1-BIE runtime.json启动故障修复

仅本地修复/检查；没有SSH连接或C1–C4运行。原失败包、远端session、日志、退出码及历史数值结果均未修改。

根因是AWK的`printf`与输出重定向语法冲突：未括号化的`a>b?a:b`不能可靠表达待打印的最大值。本机Git Bash AWK对原句直接报syntax error；其他实现也可能按重定向解析，不能依赖它产生有效数字。空输出再被原字符串拼接写入JSON，会形成非法配置。

修复为：

- 独立Bash配置helper先验证当前/历史值是非负整数字符串，规范前导零，限制在MATLAB double可精确表示及Bash算术安全的范围；再使用`if ((available>initial))`求最大值。
- `starting_available_bytes`保留本次可用RAM，`initial_available_bytes`保留历史基准，`RAM_baseline_bytes=max(starting,initial)`，不互相替代。
- 系统文件不可读、MemAvailable缺失/重复、单位不符、空值/非数字/超界均报错。首次运行可建立初始文件；恢复时初始文件缺失明确停止。
- 新Python标准库helper通过`json.dump(...,allow_nan=False)`写临时文件，解析/验证字段类型、语义和策略后原子发布，再对最终文件独立回读验证。已有runtime文件拒绝覆盖。
- 任何配置失败以78退出并写本次session的停止原因，发生于MATLAB启动之前；运行期间内存读取失败也明确停止本任务。
- 其他当前启动配置逻辑已检查：旧最大值AWK及runtime JSON拼接均已移除；RSS求和与RAM门的AWK不涉及printf重定向歧义。历史启动器未改写。

生产远端新增明确依赖`python3`（仅标准库）；缺失时停止，不自动安装或换环境。本地Git Bash测试仅在测试进程内把python3映射到现有Anaconda Python，不写入生产包或系统配置。

验证：29项Bash/JSON测试通过，包含大于/等于/小于、0、前导零、空值、非数字、负数、小数、指数、超界、文件缺失、系统读取失败、非法JSON和拒绝覆盖；PowerShell两个脚本解析及5组JSON独立语义检查通过。新包37项SHA、Bash dry-run通过；带7个声学入口哨兵的MATLAB dry-run通过，Code Analyzer 0消息，声学调用/SSH连接0。证据：`results/validation/r1_bie_runtime_fix/`。

新包：`results/comparisons/piston_2d_pm/20261116_r1_bie_manual_runtime_fix_02/`；新远端目录：`/home/misc/cwx/r1_bie_manual_runtime_fix_20261116_02`。不要向旧失败目录覆盖上传。

```powershell
$p = 'E:\MISC\CARPE3D_matlab\Explain\results\comparisons\piston_2d_pm\20261116_r1_bie_manual_runtime_fix_02'
$ctl = Join-Path $p 'code\manage_r1_bie_manual.ps1'
& $ctl -Action UploadAndStart -PackageDir $p
# 后续按需：
& $ctl -Action Status -PackageDir $p
& $ctl -Action Fetch -PackageDir $p
```

除新的runtime输入校验/序列化外，C0导入、C1–C4配方、缓存核心、30%RAM策略、原数值/收敛门、断点恢复及最终独立审计保持不变。远端实际启动链仍待用户手动验证；本轮不宣称C1–C4数值PASS，不提交/推送Git。
