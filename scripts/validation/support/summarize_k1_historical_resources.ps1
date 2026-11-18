param([Parameter(Mandatory=$true)][string]$OutputDir)
$ErrorActionPreference='Stop'
$repo=[IO.Path]::GetFullPath((Join-Path $PSScriptRoot '../../..')); $out=[IO.Path]::GetFullPath($OutputDir)
if(Test-Path -LiteralPath $out){throw 'New historical evidence directory required'}
New-Item -ItemType Directory -Path $out|Out-Null
$base=Join-Path $repo results/comparisons/piston_2d_pm
$runs=@(@{name='R1';root='r1_bie_manual_c4_policy_v2_20261116_04_retrieved_20261116_173609'},@{name='R2_260002';root='r2_seed260002_retrieved_20261117_01'},@{name='R2_260003';root='r2_seed260003_retrieved_20261117_01'})
$ids=@('bie_coarse','bie_main','bie_self','bie_space','bie_window'); $nodes=@(3200,3840,3840,4480,4480)
$rows=@(); $sources=@(); $sessions=@()
foreach($run in $runs){
 $root=Join-Path $base $run.root
 for($j=0;$j -lt 5;$j++){
  $path=Join-Path $root ('jobs/'+$ids[$j]+'/timings.json')
  if(-not(Test-Path -LiteralPath $path)){$path=Join-Path $root ('origin_evidence/jobs/'+$ids[$j]+'/timings.json')}
  if(-not(Test-Path -LiteralPath $path)){throw ('Missing historical timing '+$path)}
  $t=Get-Content -LiteralPath $path -Raw|ConvertFrom-Json
  $r=Join-Path $root ('control/resources/C'+$j+'.csv')
  if(-not(Test-Path -LiteralPath $r)){$r=Join-Path $root ('origin_evidence/control/resources/C'+$j+'.csv')}
  $rss=$null
  if(Test-Path -LiteralPath $r){$data=@(Import-Csv -LiteralPath $r -Header epoch,available_bytes,RSS_kB); $rss=($data|ForEach-Object {[double]$_.RSS_kB*1024}|Measure-Object -Maximum).Maximum; $sources+=@{path=$r.Substring($repo.Length+1).Replace('\','/');sha256=(Get-FileHash -LiteralPath $r).Hash.ToLowerInvariant()}}
  $call=$t.solver_call_s; $boundary='solver_call_s'
  if($null -eq $call){$call=$t.source_C0_function_s; $boundary='imported_C0_function_s; component timings unavailable'}
  $rows+=[pscustomobject]@{run=$run.name;stage=('C'+$j);nodes=$nodes[$j];geometry_s=$t.geometry_s;cache_prepare_s=$t.cache_prepare_s;base_assembly_s=$t.base_assembly_s;near_replacement_s=$t.near_replacement_s;dense_solve_s=$t.dense_solve_s;receivers_s=$t.receivers_s;boundary_residual_s=$t.boundary_residual_s;solver_call_s=$call;timer_boundary=$boundary;stage_process_group_peak_RSS_bytes=$rss;single_complex_dense_matrix_bytes=16*[double]$nodes[$j]*$nodes[$j];density_reuse=$false}
  $sources+=@{path=$path.Substring($repo.Length+1).Replace('\','/');sha256=(Get-FileHash -LiteralPath $path).Hash.ToLowerInvariant()}
 }
 foreach($d in Get-ChildItem -LiteralPath (Join-Path $root sessions) -Directory){
  $p=Join-Path $d.FullName time.txt; if(-not(Test-Path -LiteralPath $p)){continue}; $text=[IO.File]::ReadAllText($p)
  $max=[regex]::Match($text,'Maximum resident set size \(kbytes\):\s*(\d+)')
  $elapsed=[regex]::Match($text,'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*(\S+)')
  $groupPeak=$null; $resource=Join-Path $d.FullName resources.csv
  if(Test-Path -LiteralPath $resource){$all=@(Import-Csv -LiteralPath $resource -Header epoch,available_bytes,RSS_kB,stage); $groupPeak=($all|ForEach-Object {[double]$_.RSS_kB*1024}|Measure-Object -Maximum).Maximum; $sources+=@{path=$resource.Substring($repo.Length+1).Replace('\','/');sha256=(Get-FileHash -LiteralPath $resource).Hash.ToLowerInvariant()}}
  $sessions+=[pscustomobject]@{run=$run.name;session=$d.Name;time_v_peak_single_process_RSS_bytes=[double]$max.Groups[1].Value*1024;sampled_whole_session_group_peak_RSS_bytes=$groupPeak;elapsed_text=$elapsed.Groups[1].Value;scope='time-v maximum single process; whole session RSS includes model/finalization phases'}
  $sources+=@{path=$p.Substring($repo.Length+1).Replace('\','/');sha256=(Get-FileHash -LiteralPath $p).Hash.ToLowerInvariant()}
 }
}
$rows|Export-Csv -LiteralPath (Join-Path $out BIE_stage_costs.csv) -NoTypeInformation -Encoding utf8
$sessions|Export-Csv -LiteralPath (Join-Path $out sessions.csv) -NoTypeInformation -Encoding utf8
$totals=@(); foreach($run in $runs){$r=@($rows|Where-Object {$_.run -eq $run.name}); $totals+=@{run=$run.name;five_solver_calls_s=($r.solver_call_s|Measure-Object -Sum).Sum;assembly_plus_boundary_s=(($r|ForEach-Object {$_.base_assembly_s+$_.near_replacement_s+$_.boundary_residual_s})|Measure-Object -Sum).Sum;dense_solve_s=($r.dense_solve_s|Measure-Object -Sum).Sum}}
$receipt=[ordered]@{status='HISTORICAL_BIE_RESOURCE_EVIDENCE_EXTRACTED';failed_checks=0;stages=$rows.Count;runs=$totals;largest_nodes=4480;maximum_sampled_BIE_stage_group_RSS_bytes=($rows.stage_process_group_peak_RSS_bytes|Measure-Object -Maximum).Maximum;maximum_whole_session_group_RSS_bytes=($sessions.sampled_whole_session_group_peak_RSS_bytes|Measure-Object -Maximum).Maximum;working_estimate_bytes_per_K1_case=12884901888;K0_RSS_used=$false;memory_scaling='N^2 dense complex matrix plus multiple temporaries; same maximum4480 nodes';compute_scaling='assembly/residual nominal quadratic, dense LU nominal cubic; use measured stage timing rather than a total cubic extrapolation';new_acoustic_calls=0;sources_sha256=$sources}
[IO.File]::WriteAllText((Join-Path $out receipt.json),($receipt|ConvertTo-Json -Depth 15),[Text.UTF8Encoding]::new($false))
[pscustomobject]$receipt|Select-Object status,stages,largest_nodes,maximum_sampled_BIE_stage_group_RSS_bytes,maximum_whole_session_group_RSS_bytes,runs|ConvertTo-Json -Depth 6
