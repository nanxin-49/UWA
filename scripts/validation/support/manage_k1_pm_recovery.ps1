param(
 [Parameter(Mandatory=$true)][ValidateSet('DryRun','ShowCommands','Upload','Start','Status','Fetch')][string]$Action,
 [Parameter(Mandatory=$true)][string]$PackageDir,
 [string]$TaskName,
 [string]$DownloadDir,
 [string]$SshConfig
)
$ErrorActionPreference='Stop'; $pkg=[IO.Path]::GetFullPath($PackageDir)
$cfg=Get-Content -LiteralPath (Join-Path $pkg config.json) -Raw|ConvertFrom-Json
if(-not $TaskName){$TaskName=$cfg.remote_task_name}
if($TaskName -notmatch '^k1_pm_[A-Za-z0-9_-]+$' -or $TaskName -ne $cfg.remote_task_name){throw 'Use the frozen task name'}
$remote='/home/misc/cwx/'+$TaskName; $alias='matlab-ubuntu-new'
$sshOptions=@('-o','BatchMode=yes','-o','ConnectTimeout=10'); if($SshConfig){$sshOptions=@('-F',$SshConfig)+$sshOptions}
function Preflight {
 & python -B (Join-Path $pkg code/k1_pm_recovery_preflight.py) $pkg
 if($LASTEXITCODE -ne 0){throw 'K1 input/anchor/SHA preflight failed'}
}
function Invoke-K1Ssh([string]$Command){
 & ssh @sshOptions $alias $Command
 if($LASTEXITCODE -ne 0){throw ('SSH failed with exit '+$LASTEXITCODE)}
}
switch($Action){
 'DryRun'{Preflight; Write-Output 'DRY_RUN_ONLY: no network, MATLAB or acoustic solve; live admission not evaluated'}
 'ShowCommands'{Write-Output "REMOTE_ROOT=$remote"; Write-Output "START: bash $remote/code/launch_k1_pm_recovery.sh start $remote"; Write-Output 'One process, one computational thread; only the user calls Start'}
 'Upload'{
  Preflight; $archive=Join-Path $pkg input_package.tar
  $expected=(Get-Content -LiteralPath (Join-Path $pkg input_package.tar.sha256) -Raw).Trim().Split(' ')[0]
  if($expected -notmatch '^[a-f0-9]{64}$' -or (Get-FileHash -LiteralPath $archive).Hash.ToLowerInvariant() -ne $expected){throw 'Input archive SHA mismatch'}
  Invoke-K1Ssh "test ! -e '$remote' && mkdir '$remote'"
  & scp @sshOptions $archive "${alias}:$remote/input_package.tar"
  if($LASTEXITCODE -ne 0){throw 'Upload failed; preserve partial directory'}
  Invoke-K1Ssh "cd '$remote' && printf '%s  input_package.tar\n' '$expected' | sha256sum -c - && tar -xf input_package.tar && python3 -B code/k1_pm_recovery_preflight.py '$remote'"
 }
 'Start'{Preflight; Invoke-K1Ssh "bash '$remote/code/launch_k1_pm_recovery.sh' start '$remote'"}
 'Status'{Invoke-K1Ssh "bash '$remote/code/launch_k1_pm_recovery.sh' status '$remote'"}
 'Fetch'{
  $reply=& ssh @sshOptions $alias "bash '$remote/code/launch_k1_pm_recovery.sh' fetch '$remote'"
  if($LASTEXITCODE -ne 0){throw 'Fetch failed/refused; task may still be active'}
  $paths=@($reply|Where-Object {$_ -like 'ARCHIVE_PATH=*'}); $hashes=@($reply|Where-Object {$_ -like 'ARCHIVE_SHA256=*'})
  if($paths.Count -ne 1 -or $hashes.Count -ne 1){throw 'Missing unique retrieval receipt'}
  $path=$paths[0].Substring(13); $hash=$hashes[0].Substring(15)
  if($path -notmatch ('^'+[regex]::Escape($remote)+'/retrieval/result_[0-9_]+\.tar$') -or $hash -notmatch '^[a-f0-9]{64}$'){throw 'Unsafe remote retrieval identity'}
  if(-not $DownloadDir){$DownloadDir=Join-Path (Split-Path $pkg -Parent) ($TaskName+'_retrieved_'+(Get-Date -Format yyyyMMdd_HHmmss))}
  $dest=[IO.Path]::GetFullPath($DownloadDir); if(Test-Path -LiteralPath $dest){throw 'New retrieval directory required'}
  New-Item -ItemType Directory -Path $dest|Out-Null; $local=Join-Path $dest result.tar
  & scp @sshOptions "${alias}:$path" $local
  if($LASTEXITCODE -ne 0 -or (Get-FileHash -LiteralPath $local).Hash.ToLowerInvariant() -ne $hash){throw 'Retrieval transport SHA mismatch'}
  & python -B (Join-Path $pkg code/k1_pm_recovery_preflight.py) $dest extract $local
  if($LASTEXITCODE -ne 0){throw 'Safe extraction/internal SHA failed'}
  [IO.File]::WriteAllText((Join-Path $dest transport_receipt.json),(@{remote_archive=$path;expected_remote_archive_sha256=$hash;local_archive_sha256=(Get-FileHash -LiteralPath $local).Hash.ToLowerInvariant();transport_pass=$true;numeric_status='not inferred from transport'}|ConvertTo-Json),[Text.UTF8Encoding]::new($false))
  Write-Output "FETCH_SHA_VERIFIED=$dest"
 }
}
