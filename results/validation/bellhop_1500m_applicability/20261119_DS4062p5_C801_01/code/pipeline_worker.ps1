param([string]$RunDir,[double]$ActionLimit)
$ErrorActionPreference='Stop'
$clock=[Diagnostics.Stopwatch]::StartNew()
$utf=[Text.UTF8Encoding]::new($false)
$root='C:\Users\ASUS\.codex\worktrees\bellhop_rotated_fig12b\Explain'
$rel='results/validation/bellhop_1500m_applicability/20261119_DS4062p5_C801_01'
$remote='/home/misc/cwx/bellhop_1500m_DS4062p5_C801_20261119_01'
$hostName='matlab-ubuntu-new'
$sshFlags=@('-o','BatchMode=yes','-o','ConnectTimeout=10')
function Stage([string]$name){
    [IO.File]::WriteAllText((Join-Path $RunDir 'pipeline_stage.txt'),$name,$utf)
    Write-Output $name
}
function RunStep([string]$exe,[string[]]$argv,[double]$cap,[string]$tag){
    $cap=[Math]::Min($cap,$ActionLimit-$clock.Elapsed.TotalSeconds-7)
    if($cap -le 0){throw 'WHOLE_PIPELINE_DEADLINE'}
    $info=[Diagnostics.ProcessStartInfo]::new()
    $info.FileName=$exe
    $info.WorkingDirectory=$RunDir
    foreach($arg in $argv){$info.ArgumentList.Add($arg)}
    $info.UseShellExecute=$false
    $info.CreateNoWindow=$true
    $info.RedirectStandardOutput=$true
    $info.RedirectStandardError=$true
    $p=[Diagnostics.Process]::new()
    $p.StartInfo=$info
    $timer=[Diagnostics.Stopwatch]::StartNew()
    $null=$p.Start()
    $stdout=$p.StandardOutput.ReadToEndAsync()
    $stderr=$p.StandardError.ReadToEndAsync()
    $timedout=$false
    while(-not $p.WaitForExit(30)){
        if($timer.Elapsed.TotalSeconds -ge $cap){
            $timedout=$true
            $p.Kill($true)
            $null=$p.WaitForExit(1000)
            break
        }
    }
    foreach($pair in @(@('stdout',$stdout),@('stderr',$stderr))){
        [IO.File]::WriteAllText((Join-Path $RunDir ($tag+'_'+$pair[0]+'.log')),$pair[1].GetAwaiter().GetResult(),$utf)
    }
    [IO.File]::WriteAllText((Join-Path $RunDir ($tag+'_exit.txt')),[string]$p.ExitCode,$utf)
    if($timedout){throw "BOUNDED_STEP_TIMEOUT_$tag"}
    return [int]$p.ExitCode
}
function Verify-History {
    foreach($row in Import-Csv (Join-Path $RunDir 'history_repair/source_sha256.csv')){
        if((Get-FileHash -LiteralPath (Join-Path $root $row.path)).Hash -ne $row.sha256){
            throw "HISTORY_SHA_CHANGED_$($row.path)"
        }
    }
}
function Verify-Manifest([string]$dir,[string]$manifest){
    $count=0
    foreach($line in Get-Content -LiteralPath (Join-Path $dir $manifest)){
        if($line -notmatch '^([a-f0-9]{64})  \./([^:\\]+)$'){throw 'Invalid SHA manifest'}
        $sha=$Matches[1]
        $name=$Matches[2]
        if($name.Split('/') -contains '..'){throw 'Invalid manifest path'}
        if((Get-FileHash -LiteralPath (Join-Path $dir $name)).Hash -ne $sha){throw "SHA mismatch $name"}
        $count++
    }
    if($count -eq 0){throw 'Empty manifest'}
    return $count
}
$problem=''
$inputCount=0
$returnCount=0
$numeric=$null
$metrics=$null
$native=@()
$dest=Join-Path $RunDir 'returned'
try{
    Stage 'FROZEN_INPUT_AND_HISTORY_SHA'
    Verify-History
    $historyAudit=Get-Content (Join-Path $RunDir 'history_repair/visual_review.json') -Raw | ConvertFrom-Json
    if($historyAudit.status -ne 'PASS'){throw 'Historical PNG audit not passed'}
    $inputCount=Verify-Manifest (Join-Path $RunDir 'payload') 'input_SHA256SUMS'
    Stage 'REMOTE_CREATE_AND_TRANSFER'
    if((RunStep 'tar.exe' @('-cf','payload.tar','-C','payload','.') 10 'pack') -ne 0){throw 'Snapshot archive failed'}
    if((RunStep 'ssh.exe' ($sshFlags+@($hostName,"test -d /home/misc/cwx && mkdir -- $remote")) 12 'remote_create') -ne 0){throw 'Unique remote directory failed'}
    if((RunStep 'scp.exe' ($sshFlags+@('payload.tar',"${hostName}:${remote}/payload.tar")) 20 'upload') -ne 0){throw 'Frozen snapshot upload failed'}
    Stage 'REMOTE_START_SOLVE_D_THEN_S'
    $remoteCommand='cd '+$remote+' && tar -xf payload.tar && mv payload.tar ../'+($remote.Split('/')[-1])+'_input.tar && { bash ./execute_budgeted.sh; result=$?; if test ! -s returned.tar && test ! -f remote_abort_reason.txt; then find . -type f ! -name return_SHA256SUMS ! -name returned.tar -print0 | sort -z | xargs -0 sha256sum > return_SHA256SUMS && tar --exclude=./returned.tar -cf returned.tar .; fi; exit "$result"; }'
    $rc=RunStep 'ssh.exe' ($sshFlags+@($hostName,$remoteCommand)) 184 'remote_execute'
    Stage 'RETRIEVE_AND_VERIFY_OUTPUT'
    if((RunStep 'scp.exe' ($sshFlags+@("${hostName}:${remote}/returned.tar",'returned.tar')) 25 'retrieval') -ne 0){throw "Return archive missing, remote exit $rc"}
    New-Item -ItemType Directory -Path $dest | Out-Null
    if((RunStep 'tar.exe' @('-xf','returned.tar','-C',$dest) 10 'extract') -ne 0){throw 'Return extraction failed'}
    $returnCount=Verify-Manifest $dest 'return_SHA256SUMS'
    foreach($c in @('D','S')){
        $file=Join-Path $dest "$c/case_data/process_receipt.json"
        if(Test-Path -LiteralPath $file){$native+=Get-Content -LiteralPath $file -Raw | ConvertFrom-Json}
    }
    if(Test-Path (Join-Path $dest 'DS4062_metrics.csv')){$metrics=Import-Csv (Join-Path $dest 'DS4062_metrics.csv')}
    if($rc -ne 0){throw "Remote executor exit $rc. See returned MATLAB_stderr.log and checkpoints."}
    $numeric=Get-Content (Join-Path $dest 'comparison_receipt.json') -Raw | ConvertFrom-Json
    if($numeric.native_jobs -ne 2 -or @($numeric.native_exit_codes | Where-Object {$_ -ne 0}).Count -or
       $numeric.failed_gates -or $numeric.figure.status -ne 'FIG_AND_PNG_FILE_CHECKS_PASS' -or
       @($native).Count -ne 2){throw 'Numeric or serialized figure audit failed'}
    Stage 'AWAITING_VISUAL_REVIEW'
    $visualFile=Join-Path $RunDir 'visual_review.txt'
    $visualClock=[Diagnostics.Stopwatch]::StartNew()
    while(-not (Test-Path -LiteralPath $visualFile) -and $visualClock.Elapsed.TotalSeconds -lt 55 -and
          $clock.Elapsed.TotalSeconds -lt $ActionLimit-25){Start-Sleep -Milliseconds 150}
    if(-not (Test-Path -LiteralPath $visualFile) -or (Get-Content -LiteralPath $visualFile -Raw).Trim() -ne 'PASS'){throw 'PNG visual audit failed or absent'}
}catch{
    $problem=$_.Exception.Message
    Stage 'AUDIT_FAILED_STOPPED'
}
Stage 'CLOSEOUT'
Verify-History
$status=if($problem){'STOPPED_AT_FAILED_GATE'}else{'U0_D_S_4000_AND_4062P5_CONTROL_PASS_P1_INCOMPLETE'}
$summary=if($problem){"执行停止。$problem"}else{'4000 Hz保存图件补审通过，4062.5 Hz D、S与D+S全部冻结门及图件审计通过。'}
$reportRel='reports/bellhop_1500m_p1_ds4062_2026-11-19.md'
$out=[Collections.Generic.List[string]]::new()
@("# 1500 m P1 修订验证",'',"**$summary**",'','2026-11-19。实验分支 experiment/bellhop-rotated-fig12b。','',
'4000 Hz 图从保存的 S801 comparison.mat 和 CSV 重建。重新打开 FIG 后逐项核对六条曲线、两个图例位置与标签及 AutoUpdate=off。PNG 已目视检查，两处图例均完整。旧失败图、失败回执、原始数值及 SHA 继续保留。4000 Hz 没有新增求解。',
'',
'4062.5 Hz 使用原冻结计划第 2、4 槽，按 D 后 S 串行，均为正式 REAL8 C ELF、801 束。与对应 4000 Hz 验证输入只差频率。U0 等声速、1 m 参考全向三维点源，源宽不适用。水深1500 m，Tx1490 m、Rx6 m，五个独立偏移 0/10/50/100/200 m。步长0.5 m、Q间距0.625 m。保留原扇区及固定相位转换。',
'',
'先做复数相干叠加。参考为独立自由场 D 和负号镜像 S。逐点复差以独立路径包络 E 归一化，相位只在数值和参考的模均达到 0.01E 时评价。冻结门为复L2≤0.001、最大复差/E≤0.002、相位≤0.01 rad。无拟合。','') | ForEach-Object {$out.Add($_)}
if($metrics){
    $out.Add('|场|复L2|最大复差/E|最大有效相位差/rad|有效点|门通过 L2／复差／相位|')
    $out.Add('|---|---:|---:|---:|---:|---|')
    foreach($m in $metrics){$out.Add("|$($m.field)|$($m.complex_L2)|$($m.max_complex_over_E)|$($m.max_phase_rad)|$($m.phase_support)|$($m.L2_pass)／$($m.max_complex_pass)／$($m.phase_pass)|")}
    $out.Add('')
}
if($numeric){
    $out.Add("海面独立事件检查 $($numeric.surface_checks) 项通过。全801束一次海面、零底反射，反射点x最大误差 $($numeric.max_reflection_x_error_m) m，复反射系数最大 |R+1| 为 $($numeric.max_R_plus1)。")
    $out.Add('')
    foreach($a in $numeric.axis){$out.Add("- $($a.frequency_hz) Hz 轴上 |D+S|/E 数值 $($a.axis_numeric_magnitude_over_E)，解析 $($a.axis_analytic_magnitude_over_E)，S相对D相位 $($a.axis_S_minus_D_phase_rad) rad。")}
    $out.Add('')
    $out.Add('轴线两路径时差为8 ms，4062.5 Hz 相比4000 Hz多半个周期，压力释放负号使原相消转为相长。')
    $out.Add('')
}
if($native){$out.Add("原生进程耗时 $($native.elapsed_s -join ' / ') s，退出码 $($native.exit_code -join ' / ')。")}
$out.Add("输入SHA $inputCount 项、返回SHA $returnCount 项。统一 300 s 监督以最终 supervisor_receipt.json 为准。远端原1.5 GiB启动门、256 MiB及30%运行保留、20 GiB磁盘门和2 GiB总输出限制。")
$out.Add('')
$out.Add("证据目录 [本轮结果](../$rel/) 包含修复图、冻结输入、原始输出与日志、审计CSV/MAT、PNG/FIG及SHA清单。")
$out.Add('')
$out.Add('结论限定为 U0 的 D/S 压力释放双频控制。整个P1尚待 B、BS、其余模式及对照任务。等待人工审查。')
[IO.File]::WriteAllLines((Join-Path $root $reportRel),$out,$utf)
$tracker=Join-Path $root 'reports/active_task_tracker.md'
$rows=Get-Content -LiteralPath $tracker
$rows=@($rows | ForEach-Object {
    if($_.StartsWith('| Bellhop 1500m P1 U0 |')){
        "| Bellhop 1500m P1 U0 | 2026-11-19；U0 1500m，4000/4062.5 Hz，801束，5独立位置，3D点源宽N/A | **$status**。$summary [中文报告](bellhop_1500m_p1_ds4062_2026-11-19.md) | 剩余P1任务待审。 |"
    }else{$_}
})
$rows+= "最后动作（2026-11-19）$summary 等待人工审查。"
[IO.File]::WriteAllLines($tracker,$rows,$utf)
$ctx=Join-Path $root 'PROJECT_CONTEXT.md'
$ctxText=[IO.File]::ReadAllText($ctx)
$anchor=$ctxText.IndexOf('### S4000/C801最小海面验证')
if($anchor -lt 0){throw 'Project context heading missing'}
[IO.File]::WriteAllText($ctx,$ctxText.Insert($anchor,"### P1 U0 D/S双频控制（2026-11-19）`n`n$summary [报告]($reportRel)。原4000 Hz图件审计失败记录保留，整个P1待审。`n`n"),$utf)
$progressPath=Join-Path $root 'reports/project_progress.json'
$progress=Get-Content -LiteralPath $progressPath -Raw | ConvertFrom-Json
$b=$progress.workstreams | Where-Object id -eq 'bellhop'
$b.status=if($problem){'stopped'}else{'partial'}
$b.status_label=$status
$b.summary=$summary
$b.meaning='本轮限定1500 m U0、全向3D点源、五个独立位置的D/S及压力释放镜像控制。'
$b.gap='整个P1尚待B/BS、A/R、以及原始REAL8对照等证据。'
$b.next='等待人工审查。'
$b.evidence+=@([pscustomobject]@{label='1500 m P1 D/S双频控制';path=$reportRel;kind='report';revision='';expected_exists=$true})
$progress.events+=@([pscustomobject]@{id='bellhop-u0-ds4062-c801-20261119';date='2026-11-19';title='补审4000 Hz图件与4062.5 Hz修订验证';text=$summary;workstream_id='bellhop'})
$progress.headline="$summary 整个P1尚待其他验证。"
[IO.File]::WriteAllText($progressPath,($progress|ConvertTo-Json -Depth 100),$utf)
if((RunStep 'C:\Program Files\nodejs\node.exe' @((Join-Path $RunDir 'code/pack_progress.cjs'),$root) 8 'progress_pack') -ne 0){throw 'Progress packaging failed'}
$checkpoint=[ordered]@{status=$status;failure=$problem;input_sha_count=$inputCount;returned_sha_count=$returnCount;historical_sha_verified=7;native_receipts_found=$native.Count;native=$native;report=$reportRel;historical_4000_Hz_recomputed=$false;next='manual review only';supervisor='supervisor_receipt.json'}
[IO.File]::WriteAllText((Join-Path $RunDir 'checkpoint_final.json'),($checkpoint|ConvertTo-Json -Depth 15),$utf)
$files=@(Get-ChildItem (Join-Path $RunDir 'history_repair') -File)
if(Test-Path -LiteralPath $dest){$files+=@(Get-ChildItem -LiteralPath $dest -Recurse -File)}
$files+=Get-Item -LiteralPath (Join-Path $RunDir 'checkpoint_final.json'),(Join-Path $root $reportRel)
$files | ForEach-Object {
    [pscustomobject]@{path=$_.FullName.Substring($root.Length+1).Replace('\','/');sha256=(Get-FileHash -LiteralPath $_.FullName).Hash.ToLowerInvariant()}
} | Export-Csv -LiteralPath (Join-Path $RunDir 'delivery_sha256.csv') -NoTypeInformation
Stage $status
if($problem){exit 20}
