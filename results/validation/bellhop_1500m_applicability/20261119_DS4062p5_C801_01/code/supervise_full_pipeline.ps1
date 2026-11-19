param([string]$RunDir,[string]$Worker,[double]$BudgetSeconds=300)
$ErrorActionPreference='Stop'
$clock=[Diagnostics.Stopwatch]::StartNew()
$sealReserve=[Math]::Min([double]5,[double]($BudgetSeconds*.2))
$actionLimit=$BudgetSeconds-$sealReserve
$utf8=[Text.UTF8Encoding]::new($false)
[IO.File]::WriteAllText((Join-Path $RunDir 'supervisor_start_utc.txt'),[DateTimeOffset]::UtcNow.ToString('O'),$utf8)
$info=[Diagnostics.ProcessStartInfo]::new()
$info.FileName=(Get-Process -Id $PID).Path
foreach($arg in @('-NoProfile','-File',$Worker,'-RunDir',$RunDir,'-ActionLimit',"$actionLimit")){$info.ArgumentList.Add($arg)}
$info.UseShellExecute=$false; $info.CreateNoWindow=$true
$info.RedirectStandardOutput=$true; $info.RedirectStandardError=$true
$p=[Diagnostics.Process]::new(); $p.StartInfo=$info
$null=$p.Start(); $stdout=$p.StandardOutput.ReadToEndAsync(); $stderr=$p.StandardError.ReadToEndAsync()
$expired=$false; $lastStage=''
while(-not $p.WaitForExit(30)){
 if($clock.Elapsed.TotalSeconds -ge $actionLimit){$expired=$true; $p.Kill($true); $null=$p.WaitForExit(1000); break}
 $stageFile=Join-Path $RunDir 'pipeline_stage.txt'
 if(Test-Path -LiteralPath $stageFile){$stage=[IO.File]::ReadAllText($stageFile); if($stage -ne $lastStage){Write-Output ("PIPELINE "+$stage); $lastStage=$stage}}
}
if(-not $p.HasExited){throw 'Task worker did not terminate'}
$rc=$p.ExitCode
[IO.File]::WriteAllText((Join-Path $RunDir 'worker_stdout.log'),$stdout.GetAwaiter().GetResult(),$utf8)
[IO.File]::WriteAllText((Join-Path $RunDir 'worker_stderr.log'),$stderr.GetAwaiter().GetResult(),$utf8)
$receipt=[ordered]@{status=$(if($expired){'HARD_DEADLINE_STOP'}elseif($rc -eq 0){'PIPELINE_FINISHED'}else{'PIPELINE_FAILED_STOPPED'});budget_s=$BudgetSeconds;action_cutoff_s=$actionLimit;elapsed_s=$clock.Elapsed.TotalSeconds;worker_exit_code=$rc;expired=$expired;worker_pid=$p.Id;scope='MATLAB startup/native solve/retrieval/local audits/figures/reports/navigation/seal all inside supervised worker';within_budget=($clock.Elapsed.TotalSeconds -lt $BudgetSeconds);new_jobs_after_deadline=0}
[IO.File]::WriteAllText((Join-Path $RunDir 'supervisor_receipt.json'),($receipt|ConvertTo-Json),$utf8)
Write-Output ("SUPERVISOR "+($receipt|ConvertTo-Json -Compress))
if($expired){exit 124}; exit $rc