param([Parameter(Mandatory=$true)][string]$SourceDir,[Parameter(Mandatory=$true)][string]$PackageDir,[Parameter(Mandatory=$true)][string]$OutputDir)
$ErrorActionPreference='Stop'
$root=[IO.Path]::GetFullPath($SourceDir); $package=[IO.Path]::GetFullPath($PackageDir)
$dest=[IO.Path]::GetFullPath($OutputDir)
if(Test-Path -LiteralPath $dest){throw 'New D0 review directory required'}
New-Item -ItemType Directory -Path $dest|Out-Null
$utf8=[Text.UTF8Encoding]::new($false)
function Json($path){Get-Content -LiteralPath $path -Raw|ConvertFrom-Json}
function Hash($path){(Get-FileHash -LiteralPath $path -Algorithm SHA256).Hash.ToLowerInvariant()}
function Verify($base,$rows){
 foreach($row in $rows){
  if($row.path -match '(^/|:|\\|(^|/)\.\.(/|$)|^cash(/|$))'){throw 'Unsafe relative path'}
  if((Hash (Join-Path $base $row.path)) -ne $row.sha256){throw ('SHA drift '+$row.path)}
 }
}
$t=[Diagnostics.Stopwatch]::StartNew()
$manifest=Json (Join-Path $root retrieval/manifest.json)
Verify $root $manifest
$input=Json (Join-Path $root manifest.json); Verify $root $input.files
if((Hash (Join-Path $root manifest.json)) -ne (Hash (Join-Path $package manifest.json))){throw 'Frozen package mismatch'}
$receipt=Get-Content -LiteralPath ($root+'_transport/remote_archive_receipt.txt') -Raw
$sha=[regex]::Match($receipt,'ARCHIVE_SHA256=([a-f0-9]{64})').Groups[1].Value
if($sha.Length -ne 64 -or (Hash ($root+'_transport/result.tar')) -ne $sha){throw 'Transport archive SHA mismatch'}
$cp=Json (Join-Path $root output/checkpoint.json)
if($cp.status -ne 'D0_COMPLETE_INCIDENT_ONLY' -or @($cp.failed_or_unrun_checks).Count){throw 'Not a complete D0 checkpoint'}
Verify (Join-Path $root output) $cp.artifacts_sha256
$binding=Hash (Join-Path $root manifest.json)
$cases=@();$all=@()
foreach($depth in @(100,300,600,1000,1490)){
 $folder=Join-Path $root ('output/depth_{0:d4}m' -f $depth)
 $c=Json (Join-Path $folder checkpoint.json)
 if($c.status -ne 'AUDITED_PASS' -or $c.manifest_sha256 -ne $binding){throw 'Case checkpoint mismatch'}
 Verify $folder $c.artifacts_sha256
 $a=Json (Join-Path $folder audit.json)
 if($a.failed_checks -ne 0 -or @($a.checks|Where-Object {-not $_.passed -or $_.value -ge $_.limit}).Count){throw 'Failed numerical gate'}
 foreach($v in $a.checks){$all+=[pscustomobject]@{depth_m=$depth;name=$v.name;value=$v.value;limit=$v.limit;passed=$v.passed}}
 $cases+=[pscustomobject]@{depth_m=$depth;checks=$a.checks.Count;artifacts=$c.artifacts_sha256.Count;elapsed_s=$c.resource_usage.elapsed_s;evanescent_branch_pressure_norm_ratio=$a.branch_pressure_relative}
}
$sessions=@(Get-ChildItem -LiteralPath (Join-Path $root sessions) -Directory)
if($sessions.Count -ne 1){throw 'Expected exactly one D0 session'}
$exit=Json (Join-Path $sessions[0].FullName exit.json)
$runtime=Json (Join-Path $sessions[0].FullName runtime.json)
if($exit.matlab_exit_code -ne 0 -or $exit.supervisor_exit_code -ne 0 -or $exit.status -ne 'D0_COMPLETE_INCIDENT_ONLY'){throw 'Failed real exit status'}
if($runtime.manifest_sha256 -ne $binding -or $runtime.workers -ne 1 -or $runtime.command -notcontains '-singleCompThread'){throw 'Runtime policy mismatch'}
if($exit.minimum_available_bytes -lt .30*$runtime.initial_available_bytes){throw 'RAM gate violated'}
$all|Export-Csv -LiteralPath (Join-Path $dest numerical_gates.csv) -NoTypeInformation -Encoding utf8
$cases|Export-Csv -LiteralPath (Join-Path $dest depth_receipts.csv) -NoTypeInformation -Encoding utf8
$result=[ordered]@{status='RETRIEVAL_SHA_AND_REMOTE_GATES_PASS';manifest_files=$manifest.Count;input_files=$input.files.Count;source_dir=$root;archive_sha256=$sha;manifest_sha256=$binding;checks=$all.Count;failed_checks=0;remote_elapsed_s=$exit.elapsed_s;remote_matlab_exit_code=$exit.matlab_exit_code;remote_supervisor_exit_code=$exit.supervisor_exit_code;minimum_RAM_GiB=($exit.minimum_available_bytes/1GB);minimum_RAM_fraction=($exit.minimum_available_bytes/$runtime.initial_available_bytes);owned_PID_peak_RSS_GiB=($exit.owned_PID_peak_RSS_bytes/1GB);elapsed_s=$t.Elapsed.TotalSeconds;scope='SHA plus saved remote gates only; local MAT/CSV/FIG audit is separate';SSH_calls=0;new_solves=0}
[IO.File]::WriteAllText((Join-Path $dest receipt.json),($result|ConvertTo-Json -Depth 15),$utf8)
$result|ConvertTo-Json -Depth 15
