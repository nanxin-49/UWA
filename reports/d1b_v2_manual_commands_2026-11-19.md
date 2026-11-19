# D1b-v2 资源原型手动命令

`02` 已由用户手动运行并通过保存审查。以下命令保留为本次执行记录。 [实测与准入审查](d1b_v2_resource_review_2026-11-19.md)

在配置了 `matlab-ubuntu-new` SSH 别名的 Windows 工作站，以 PowerShell 7 执行。以下为通信统计修复后的 `02` 包。首跑 `01` 的日志与失败保留。`Start` 只派发最多 1,800 s 的原型。 [D1b-v2 方案](d1b_v2_plan_2026-11-19.md) · [首跑与修复](d1b_v2_first_run_repair_2026-11-19.md)

```powershell
$tool = 'E:\MISC\CARPE3D_matlab\Explain\scripts\validation\support\manage_d1b_v2.ps1'
$package = 'E:\MISC\CARPE3D_matlab\Explain\results\validation\d1b_v2\preparation_20261119_02'
$run = 'd1b_v2_20261119_02'

& $tool -Action DryRun -PackageDir $package -TaskName $run
```

看到 `DRY_RUN_PASS` 后上传。

```powershell
& $tool -Action Upload -PackageDir $package -TaskName $run
```

上传和复核完成后，单独执行 Start。

```powershell
& $tool -Action Start -PackageDir $package -TaskName $run
```

Start 返回后可关闭终端。需要时查看状态，出现终态再取回数据。

```powershell
& $tool -Action Status -PackageDir $package -TaskName $run
& $tool -Action Fetch -PackageDir $package -TaskName $run
```

目标目录 `/home/misc/cwx/d1b_v2_20261119_02` 专属。任务目录已存在时 Upload 退出，已有 `output` 时 Start 退出。遇到超时或启动异常，先查看原任务 Status，保留失败目录和日志。终态归档含 `output/exit.json`、`output/checkpoint.json`、每阶段 SHA、`matlab.log`、`supervisor.log`、RSS/压力样本及可用的 `resource_summary.json`。原始 PCT JobStorage 缓存留在对应任务目录。Fetch 验证归档和内部 SHA 后取回到新目录。`PROTOTYPE_COMPLETE_REVIEW_REQUIRED` 表明有原型结果可审阅，完整声场另审。

冻结包 [input_package.tar](../results/validation/d1b_v2/preparation_20261119_02/input_package.tar) 的 SHA-256 为 `e1fd5c04ff5f2a7d8921574aa2a7ad8d961ad0386c34dca853df85193101d9a8`。归档输入清单 SHA-256 为 `c5bcccfa415b0e7494c140635ff348c9b1f5897fe897b3417d4591d59a99b492`。
