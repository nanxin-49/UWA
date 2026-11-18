param(
 [Parameter(Mandatory=$true)][string]$SourceDir,
 [Parameter(Mandatory=$true)][string]$OutputDir
)
# Read-only retrieval/provenance/CSV review. Does not replace MATLAB MAT/FIG audits.
$ErrorActionPreference='Stop'
$repo=[IO.Path]::GetFullPath((Join-Path $PSScriptRoot '../../..'))
$source=[IO.Path]::GetFullPath($SourceDir)
$out=[IO.Path]::GetFullPath($OutputDir)
$sourceBound=Join-Path $repo 'results/experiments/pm_wind_cutoff'
$outBound=Join-Path $repo 'results/validation/pm_wind_cutoff_retrieval'
if(-not $source.StartsWith($sourceBound+[IO.Path]::DirectorySeparatorChar,[StringComparison]::OrdinalIgnoreCase)){throw 'Unexpected retrieval location'}
if(-not $out.StartsWith($outBound+[IO.Path]::DirectorySeparatorChar,[StringComparison]::OrdinalIgnoreCase)){throw 'Unexpected review location'}
if(Test-Path -LiteralPath $out){throw 'New review directory required'}
New-Item -ItemType Directory -Path $out|Out-Null
$utf8=[Text.UTF8Encoding]::new($false)
function Save-Json($path,$value){[IO.File]::WriteAllText($path,($value|ConvertTo-Json -Depth 50),$utf8)}
function Read-Json($path){Get-Content -LiteralPath $path -Raw|ConvertFrom-Json}
function Sha($path){(Get-FileHash -LiteralPath $path -Algorithm SHA256).Hash.ToLowerInvariant()}
function Check($ok,$message){if(-not $ok){throw $message}}
function Verify-Rows($root,$rows){
 $count=0; $seen=[Collections.Generic.HashSet[string]]::new([StringComparer]::OrdinalIgnoreCase)
 foreach($row in $rows){
  Check ($row.path -notmatch '(^/|:|\\|(^|/)\.\.(/|$)|^cash(/|$))') ('Unsafe path '+$row.path)
  Check ($seen.Add($row.path)) ('Duplicate path '+$row.path)
  $file=Join-Path $root $row.path
  Check (Test-Path -LiteralPath $file -PathType Leaf) ('Missing '+$row.path)
  Check (((Get-Item -LiteralPath $file).Attributes -band [IO.FileAttributes]::ReparsePoint) -eq 0) ('Link '+$row.path)
  Check ((Sha $file) -eq $row.sha256) ('SHA drift '+$row.path)
  $count++
 }
 return $count
}
Start-Transcript -LiteralPath (Join-Path $out review.log)|Out-Null
$timer=[Diagnostics.Stopwatch]::StartNew()
try {
 $transfer=Read-Json (Join-Path $source retrieval/manifest.json)
 $transferCount=Verify-Rows $source $transfer
 $package=Read-Json (Join-Path $source manifest.json)
 $packageCount=Verify-Rows $source $package.files
 $binding=Sha (Join-Path $source manifest.json)
 $prepared=Join-Path $sourceBound preparation_20261117_05
 Check ($binding -eq (Sha (Join-Path $prepared manifest.json))) 'Retrieved package differs from frozen05'
 Check ($binding -eq 'cf48069437260281d34a560ba6c3bcee1047d61acdc18154e8c1d68aeb6fd6c6') 'Unexpected package binding'
 $cfg=Read-Json (Join-Path $source config.json)
 Check ($cfg.kind -eq 'K0_SURFACE_ONLY' -and $cfg.wind_label -eq 'U_PM' -and $cfg.K1_enabled -eq $false) 'Configuration scope'
 Check ($cfg.cases.Count -eq 45 -and $cfg.sample_count -eq 8192 -and $cfg.L_m -eq 160) 'Matrix/grid'
 $complete=Read-Json (Join-Path $source output/checkpoint.json)
 Check ($complete.status -eq 'K0_COMPLETE_SURFACE_ONLY' -and $complete.case_count -eq 45 -and $complete.manifest_sha256 -eq $binding) 'Global checkpoint'
 $globalCount=Verify-Rows (Join-Path $source output) $complete.artifacts_sha256
 $sessions=@(Get-ChildItem -LiteralPath (Join-Path $source sessions) -Directory)
 Check ($sessions.Count -eq 1) 'Review expects exactly one saved session'
 $session=$sessions[0].FullName
 $exit=Read-Json (Join-Path $session exit.json); $runtime=Read-Json (Join-Path $session runtime.json)
 Check ($exit.matlab_exit_code -eq 0 -and $exit.supervisor_exit_code -eq 0 -and $exit.status -eq 'K0_COMPLETE_SURFACE_ONLY' -and -not $exit.stop_reason) 'Real exit receipt'
 Check ($runtime.command[0] -eq '/home/misc/MATLAB/R2025b/bin/matlab' -and $runtime.command[1] -eq '-singleCompThread' -and $runtime.platform -eq 'linux' -and $runtime.concurrency -eq 1) 'Runtime target/serial'
 Check ($exit.elapsed_s -lt $runtime.max_wall_s -and $exit.minimum_available_bytes -ge .3*$runtime.baseline_available_bytes) 'Resource gates'
 $stats=Import-Csv -LiteralPath (Join-Path $source output/summary_attempt_001/statistics.csv)
 Check ($stats.Count -eq 45) 'Summary rows'
 $caseRows=@(); $coeff=@{}; $caseFiles=0; $generatedFiles=0; $summaryMaxDifference=0.0
 $versionLines=@(); $expectedIDs=@(); $gaussianMaxDifference=0.0
 foreach($c in $cfg.cases){
  $id='seed{0}_U{1}_N{2:00}' -f $c.seed,$c.U_PM_mps,$c.N
  $expectedIDs+=$id; $folder=Join-Path $source ('output/'+$id)
  $cp=Read-Json (Join-Path $folder checkpoint.json); $gen=Read-Json (Join-Path $folder generated.json)
  Check ($cp.status -eq 'AUDITED_PASS' -and $cp.manifest_sha256 -eq $binding -and $gen.manifest_sha256 -eq $binding) ('Commit '+$id)
  Check ($cp.case.seed -eq $c.seed -and $cp.case.U_PM_mps -eq $c.U_PM_mps -and $cp.case.N -eq $c.N -and $cp.failed_or_unrun_checks.Count -eq 0) ('Case scope '+$id)
  $caseFiles+=Verify-Rows $folder $cp.artifacts_sha256
  $generatedFiles+=Verify-Rows $folder $gen.files
  $a=Read-Json (Join-Path $folder audit.json); $m=Read-Json (Join-Path $folder metrics.json)
  $config=Read-Json (Join-Path $folder config.json)
  Check ($config.run_config.kind -eq 'K0_SURFACE_ONLY' -and $config.run_config.family -eq 'legacy12_gaussian_extension_v1' -and $config.extension_seed -eq $c.seed+1000000) ('Family '+$id)
  Check ($a.status -eq 'AUDITED_PASS' -and $a.failed_checks -eq 0 -and $a.physical_cutoff_fixed -eq $true) ('Saved audit '+$id)
  $rms=($a.grid_rms_relative|Measure-Object -Maximum).Maximum
  Check ($rms -lt 1e-4 -and $a.grid_max_slope_absolute -lt 1e-4 -and $a.reconstruction_max_abs -lt 1e-11) ('Saved numerical gates '+$id)
  Check ($m.sample_count -eq 8192 -and $m.dx_m -eq .01953125 -and $m.L_m -eq 160) ('Sampling '+$id)
  Check ($cp.matlab_release -eq '2025b' -and $cp.execution_mode -eq 'surface_only_serial') ('Execution '+$id)
  $row=@($stats|Where-Object {$_.seed -eq [string]$c.seed -and $_.U_PM_mps -eq [string]$c.U_PM_mps -and $_.N -eq [string]$c.N})
  Check ($row.Count -eq 1) ('Summary identity '+$id)
  foreach($property in $m.PSObject.Properties){
   $delta=[Math]::Abs([double]$row[0].($property.Name)-[double]$property.Value)
   $summaryMaxDifference=[Math]::Max($summaryMaxDifference,$delta)
   Check ($delta -lt 1e-11) ('Summary/JSON disagreement '+$id+' '+$property.Name)
  }
  $table=@(Import-Csv -LiteralPath (Join-Path $folder coefficients.csv))
  Check ($table.Count -eq $c.N) ('Coefficient count '+$id); $coeff[$id]=$table
  $versionLines+=Get-Content -LiteralPath (Join-Path $folder generation.log) -TotalCount 1
  Check (@(Get-ChildItem -LiteralPath $folder -Filter generation_started.json -File).Count -eq 1) ('Generation record '+$id)
  Check (@(Get-ChildItem -LiteralPath $folder -Filter audit_attempt_*.log -File).Count -eq 1) ('Unexpected audit retries '+$id)
  $caseRows+=[pscustomobject]@{case_id=$id;seed=$c.seed;U_PM_mps=$c.U_PM_mps;N=$c.N;remote_audit_status=$a.status;remote_reconstruction_max_abs=$a.reconstruction_max_abs;remote_grid_RMS_relative_max=$rms;remote_grid_max_slope_absolute=$a.grid_max_slope_absolute;generation_s=$gen.elapsed_s;audit_s=$cp.resource_usage.audit_elapsed_s;local_MAT_reaudit='UNRUN_MCP_UNAVAILABLE'}
 }
 $actualIDs=@(Get-ChildItem -LiteralPath (Join-Path $source output) -Directory -Filter seed*|ForEach-Object {$_.Name})
 Check ($actualIDs.Count -eq 45 -and @(Compare-Object ($actualIDs|Sort-Object) ($expectedIDs|Sort-Object)).Count -eq 0) 'Extra/missing cases'
 $nestedComparisons=0; $windComparisons=0; $anchorComparisons=0
 foreach($seed in @(260001,260002,260003)){
  $anchor=@(Import-Csv -LiteralPath (Join-Path $source "inputs/seed$seed.csv"))
  foreach($u in @(5,6,8)){
   $max=$coeff[('seed{0}_U{1}_N40' -f $seed,$u)]
   foreach($n in @(7,12,19,32,40)){
    $t=$coeff[('seed{0}_U{1}_N{2:00}' -f $seed,$u,$n)]
    for($j=0;$j -lt $n;$j++){
     foreach($p in $t[$j].PSObject.Properties){Check ([double]$p.Value -eq [double]$max[$j].($p.Name)) 'Nested coefficient drift'; $nestedComparisons++}
    }
   }
   $base=$coeff[('seed{0}_U6_N40' -f $seed)]
   for($j=0;$j -lt 40;$j++){
    foreach($field in @('gaussian_cos','gaussian_sin')){
     $delta=[Math]::Abs([double]$max[$j].$field-[double]$base[$j].$field)
     $gaussianMaxDifference=[Math]::Max($gaussianMaxDifference,$delta); Check ($delta -eq 0) 'Across-U Gaussian drift'; $windComparisons++
    }
    Check ([double]$max[$j].Sk_m3 -eq [double]$base[$j].Sk_m3 -or $u -ne 6) 'Sk baseline'
    if($u -ne 6){Check ([double]$max[$j].Sk_m3 -ne [double]$base[$j].Sk_m3) 'Only wind metadata changed'}
   }
  }
  $base=$coeff[('seed{0}_U6_N12' -f $seed)]
  for($j=0;$j -lt 12;$j++){
   Check ([double]$base[$j].cos_depth_coeff_m -eq [double]$anchor[$j].cos_coeff_m) 'Historical cosine drift'
   Check ([double]$base[$j].sin_depth_coeff_m -eq [double]$anchor[$j].sin_coeff_m) 'Historical sine drift'; $anchorComparisons+=2
  }
 }
 $render=Read-Json (Join-Path $source output/summary_attempt_001/render_receipt.json)
 Check ($render.rows -eq 45 -and $render.failed_checks -eq 0 -and $render.reopened_curves -eq 81 -and $render.surface_spectrum_reopened_curves -eq 36) 'Saved FIG data audit'
 $resource=@(Import-Csv -LiteralPath (Join-Path $session resource.tsv) -Delimiter "`t")
 $minimum=($resource|ForEach-Object {[double]$_.MemAvailable_bytes}|Measure-Object -Minimum).Minimum
 $rss=($resource|ForEach-Object {[double]$_.owned_matlab_RSS_bytes}|Measure-Object -Maximum).Maximum
 Check ($minimum -eq $exit.minimum_available_bytes -and $rss -eq $exit.owned_pid_peak_rss_bytes) 'Resource TSV/exit mismatch'
 $caseRows|Export-Csv -LiteralPath (Join-Path $out remote_case_audit_index.csv) -NoTypeInformation -Encoding utf8
 $candidate=@(); foreach($row in $stats){
  if($row.seed -eq '260001' -and (($row.U_PM_mps -eq '6' -and $row.N -in @('7','12','40')) -or ($row.U_PM_mps -in @('5','8') -and $row.N -in @('12','40')))){
   $candidate+=[pscustomobject]@{seed=260001;U_PM_mps=$row.U_PM_mps;N=$row.N;Kmax_rad_per_m=$row.Kmax_rad_per_m;raw_RMS_height_m=$row.raw_rms_height_m;raw_RMS_slope=$row.raw_rms_slope;raw_RMS_second_per_m=$row.raw_rms_second_per_m;C2_max_abs_slope=$row.embedded_max_abs_slope;C2_min_elevation_m=$row.embedded_min_elevation_m;C2_max_elevation_m=$row.embedded_max_elevation_m;sampled_gap_above_auxiliary_zminus2_m=2-[double]$row.embedded_max_elevation_m;sampled_receiver_clearance_m=3+[double]$row.embedded_min_elevation_m;status='K1_PLANNING_ONLY_NO_PROPAGATION_QUALIFICATION'}
  }
 }; Check ($candidate.Count -eq 7) 'K1 candidate count'
 $candidate|Export-Csv -LiteralPath (Join-Path $out k1_candidate_geometry.csv) -NoTypeInformation -Encoding utf8
 # Retain descriptive fixed-U contrasts, without treating 411 points as repeats.
 $contrasts=@(); foreach($seed in @(260001,260002,260003)){foreach($u in @(5,6,8)){
  $a=@($stats|Where-Object {$_.seed -eq [string]$seed -and $_.U_PM_mps -eq [string]$u -and $_.N -eq '12'})[0]
  $b=@($stats|Where-Object {$_.seed -eq [string]$seed -and $_.U_PM_mps -eq [string]$u -and $_.N -eq '40'})[0]
  $contrasts+=[pscustomobject]@{seed=$seed;U_PM_mps=$u;height_ratio_N40_N12=[double]$b.raw_rms_height_m/[double]$a.raw_rms_height_m;slope_ratio_N40_N12=[double]$b.raw_rms_slope/[double]$a.raw_rms_slope;second_derivative_ratio_N40_N12=[double]$b.raw_rms_second_per_m/[double]$a.raw_rms_second_per_m}
 }}
 $contrasts|Export-Csv -LiteralPath (Join-Path $out saved_metric_contrasts.csv) -NoTypeInformation -Encoding utf8
 # Original bytes must still match after all read-only inspections.
 $afterCount=Verify-Rows $source $transfer
 $maxRms=($caseRows.remote_grid_RMS_relative_max|Measure-Object -Maximum).Maximum
 $maxSlope=($caseRows.remote_grid_max_slope_absolute|Measure-Object -Maximum).Maximum
 $maxRecon=($caseRows.remote_reconstruction_max_abs|Measure-Object -Maximum).Maximum
 $receipt=[ordered]@{status='K0_REMOTE_AUDITED_RETRIEVAL_INTEGRITY_PASS_LOCAL_MAT_REAUDIT_UNRUN';source_dir=$source;package_manifest_sha256=$binding;retrieval_manifest_sha256=(Sha (Join-Path $source retrieval/manifest.json));retrieved_tar_sha256=(Sha (Join-Path $source result.tar));retrieved_tar_expected_external_sha256='not independently retained; internal snapshot manifest verified';retrieval_files_before=$transferCount;retrieval_files_after=$afterCount;package_files=$packageCount;global_artifact_references=$globalCount;case_artifact_references=$caseFiles;generated_artifact_references=$generatedFiles;case_count=45;remote_audited_pass=45;remote_reconstruction_max_abs=$maxRecon;remote_grid_RMS_relative_max=$maxRms;remote_grid_max_slope_absolute=$maxSlope;CSV_JSON_max_abs_difference=$summaryMaxDifference;nested_scalar_comparisons=$nestedComparisons;across_U_gaussian_comparisons=$windComparisons;across_U_gaussian_max_difference=$gaussianMaxDifference;original_U6_coefficient_comparisons=$anchorComparisons;matlab_exit_code=$exit.matlab_exit_code;supervisor_exit_code=$exit.supervisor_exit_code;remote_elapsed_s=$exit.elapsed_s;minimum_available_GiB=$minimum/1GB;minimum_available_fraction_of_baseline=$minimum/$runtime.baseline_available_bytes;owned_pid_peak_RSS_GiB=$rss/1GB;resource_samples=$resource.Count;saved_FIG_reopened_curves=117;local_MAT_reaudit='UNRUN_MCP_ATTACH_FAILED';local_FIG_reopen='UNRUN_MCP_ATTACH_FAILED';SSH_connections_this_review=0;new_acoustic_solves=0;failed_checks=0;unrun_checks=@('local MATLAB MAT/CSV/numerical re-audit','local MATLAB FIG reopening');next_allowed_step='K0 saved review only; refine K1 plan, no acoustic execution';elapsed_review_s=$timer.Elapsed.TotalSeconds}
 Save-Json (Join-Path $out receipt.json) $receipt
 Save-Json (Join-Path $out source_file_hashes.json) $transfer
 Write-Output ("RETRIEVAL_PASS files={0} cases=45 nested={1} Gaussian={2} original={3} unchanged_after={4}" -f $transferCount,$nestedComparisons,$windComparisons,$anchorComparisons,$afterCount)
} catch {
 Save-Json (Join-Path $out failure.json) @{status='FAILED_REVIEW';error=$_.ToString();elapsed_s=$timer.Elapsed.TotalSeconds}
 throw
} finally {Stop-Transcript|Out-Null}
