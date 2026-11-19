param(
 [Parameter(Mandatory=$true)][ValidateSet('DryRun','UploadAndStart','Upload','Start','Status','Fetch')][string]$Action,
 [Parameter(Mandatory=$true)][string]$PackageDir,
 [string]$TaskName,
 [string]$DownloadDir
)
$ErrorActionPreference='Stop'
if($PSVersionTable.PSVersion.Major -lt 7){throw 'Use PowerShell 7 (pwsh)'}
$package=[IO.Path]::GetFullPath($PackageDir)
$cfg=Get-Content -LiteralPath (Join-Path $package config.json) -Raw|ConvertFrom-Json
if(-not $TaskName){if($cfg.recovery){$TaskName=$cfg.recovery.remote_task_name}else{$TaskName='d1_depth_20261118_01'}}
if($TaskName -notmatch '^d1_depth_[A-Za-z0-9_-]+$'){throw 'Invalid D1 task name'}
if($cfg.recovery -and $TaskName -ne $cfg.recovery.remote_task_name){throw 'Recovery requires its isolated task name'}
$alias='matlab-ubuntu-new'; $remote='/home/misc/cwx/'+$TaskName
. (Join-Path $PSScriptRoot invoke_r1_bie_status_stdin.ps1)
function Verify-Package {
 & python -B (Join-Path $package code/d1_package.py) verify $package
 if($LASTEXITCODE -ne 0){throw 'Local frozen input verification failed'}
}
function Invoke-D1Remote([string]$script){
 $r=Invoke-R1StatusStdin -Executable ssh -ArgumentList @('-T','-o','BatchMode=yes',$alias,'bash','-s','--') -InputText ("set -euo pipefail`n"+$script+"`n")
 if($r.exit_code -ne 0){throw ("SSH exit="+$r.exit_code+"`n"+$r.stderr+"`n"+$r.stdout)}
 return $r.stdout
}
function Upload-Package {
 Verify-Package
 $tar=Join-Path $package input_package.tar
 $hash=(Get-Content -LiteralPath (Join-Path $package input_package.tar.sha256) -Raw).Trim().Split(' ')[0]
 if($hash -notmatch '^[a-f0-9]{64}$' -or (Get-FileHash -LiteralPath $tar).Hash.ToLowerInvariant() -ne $hash){throw 'Archive SHA mismatch'}
 Invoke-D1Remote "test ! -e '$remote' && mkdir '$remote'"
 & scp -o BatchMode=yes $tar "${alias}:$remote/input_package.tar"
 if($LASTEXITCODE -ne 0){throw 'Upload failed; partial directory retained; no automatic retry'}
 Invoke-D1Remote @"
cd '$remote'
printf '%s  input_package.tar\n' '$hash' | sha256sum -c -
tar -xf input_package.tar
python3 code/d1_package.py verify '$remote'
"@
}
switch($Action){
 'DryRun'{
  Verify-Package
  $hash=(Get-Content -LiteralPath (Join-Path $package input_package.tar.sha256) -Raw).Trim().Split(' ')[0]
  if((Get-FileHash -LiteralPath (Join-Path $package input_package.tar)).Hash.ToLowerInvariant() -ne $hash){throw 'Archive SHA mismatch'}
  Write-Output "REMOTE_ROOT=$remote"
  Write-Output 'DRY_RUN_PASS: no network, no MATLAB process, no fields; UploadAndStart remains manual.'
 }
 'Upload'{Upload-Package}
 'UploadAndStart'{Upload-Package; Invoke-D1Remote "bash '$remote/code/launch_d1_depth.sh' start '$remote'"}
 'Start'{Verify-Package; Invoke-D1Remote "bash '$remote/code/launch_d1_depth.sh' start '$remote'"}
 'Status'{Invoke-D1Remote "bash '$remote/code/launch_d1_depth.sh' status '$remote'"}
 'Fetch'{
  $text=Invoke-D1Remote "bash '$remote/code/launch_d1_depth.sh' fetch '$remote'"
  $paths=@($text -split "`n"|Where-Object {$_ -like 'ARCHIVE_PATH=*'}); $hashes=@($text -split "`n"|Where-Object {$_ -like 'ARCHIVE_SHA256=*'})
  if($paths.Count -ne 1 -or $hashes.Count -ne 1){throw 'Missing unique archive receipt'}
  $archive=$paths[0].Substring(13).Trim(); $hash=$hashes[0].Substring(15).Trim()
  if($archive -notmatch ('^'+[regex]::Escape($remote)+'/retrieval/result_[0-9_]+\.tar$') -or $hash -notmatch '^[a-f0-9]{64}$'){throw 'Invalid archive identity'}
  if(-not $DownloadDir){$DownloadDir=Join-Path (Split-Path $package -Parent) ($TaskName+'_retrieved_'+(Get-Date -Format yyyyMMdd_HHmmss))}
  $dest=[IO.Path]::GetFullPath($DownloadDir)
  if(Test-Path -LiteralPath $dest){throw 'New retrieval directory required'}
  $transport=$dest+'_transport'
  if(Test-Path -LiteralPath $transport){throw 'Transport directory already exists'}
  New-Item -ItemType Directory -Path $transport|Out-Null
  $local=Join-Path $transport result.tar
  & scp -o BatchMode=yes "${alias}:$archive" $local
  if($LASTEXITCODE -ne 0 -or (Get-FileHash -LiteralPath $local).Hash.ToLowerInvariant() -ne $hash){throw 'Transport SHA failed'}
  [IO.File]::WriteAllText((Join-Path $transport remote_archive_receipt.txt),$text,[Text.UTF8Encoding]::new($false))
  & python -B (Join-Path $package code/d1_package.py) extract $local $dest
  if($LASTEXITCODE -ne 0){throw 'Safe extraction or internal SHA verification failed'}
  Write-Output "FETCH_VERIFIED=$dest"
  Write-Output 'Transport/SHA verified; numerical and visual interpretation awaits Codex saved-data review.'
 }
}
