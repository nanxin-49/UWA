param(
 [Parameter(Mandatory=$true)][ValidateSet('ShowCommands','Upload','Start','Status','Fetch')][string]$Action,
 [Parameter(Mandatory=$true)][string]$PackageDir,
 [string]$TaskName='pm_wind_cutoff_k0_20261117_01',
 [string]$DownloadDir
)
$ErrorActionPreference='Stop'
if($TaskName -notmatch '^pm_wind_cutoff_k0_[A-Za-z0-9_-]+$'){throw 'Invalid task name'}
$package=[IO.Path]::GetFullPath($PackageDir); $alias='matlab-ubuntu-new'
$remote='/home/misc/cwx/'+$TaskName
function Preflight {
 & python -B (Join-Path $package code/pm_k0_preflight.py) $package
 if($LASTEXITCODE -ne 0){throw 'Local manifest/config verification failed'}
}
function Invoke-K0Ssh([string]$command){
 & ssh -o BatchMode=yes $alias $command
 if($LASTEXITCODE -ne 0){throw ('SSH exit '+$LASTEXITCODE)}
}
switch($Action){
 'ShowCommands'{
  Write-Output "REMOTE_ROOT=$remote"
  Write-Output "USER_START: bash $remote/code/launch_pm_wind_cutoff_k0.sh start $remote"
  Write-Output 'No network connection made. Start never enables K1.'
 }
 'Upload'{
  Preflight
  $tar=Join-Path $package input_package.tar
  $hash=(Get-Content -LiteralPath (Join-Path $package input_package.tar.sha256) -Raw).Trim().Split(' ')[0]
  if($hash -notmatch '^[0-9a-f]{64}$' -or (Get-FileHash -LiteralPath $tar).Hash.ToLowerInvariant() -ne $hash){throw 'Archive SHA mismatch'}
  Invoke-K0Ssh "test ! -e '$remote' && mkdir '$remote'"
  & scp -o BatchMode=yes $tar "${alias}:$remote/input_package.tar"
  if($LASTEXITCODE -ne 0){throw 'Upload failed; partial task retained'}
  Invoke-K0Ssh "cd '$remote' && printf '%s  input_package.tar\n' '$hash' | sha256sum -c - && tar -xf input_package.tar && python3 code/pm_k0_preflight.py '$remote'"
 }
 'Start'{Preflight; Invoke-K0Ssh "bash '$remote/code/launch_pm_wind_cutoff_k0.sh' start '$remote'"}
 'Status'{Invoke-K0Ssh "bash '$remote/code/launch_pm_wind_cutoff_k0.sh' status '$remote'"}
 'Fetch'{
  $text=& ssh -o BatchMode=yes $alias "bash '$remote/code/launch_pm_wind_cutoff_k0.sh' fetch '$remote'"
  if($LASTEXITCODE -ne 0){throw 'Fetch refused or failed; task may still hold K0 lock'}
  $paths=@($text|Where-Object {$_ -like 'ARCHIVE_PATH=*'}); $hashes=@($text|Where-Object {$_ -like 'ARCHIVE_SHA256=*'})
  if($paths.Count -ne 1 -or $hashes.Count -ne 1){throw 'Missing archive receipt'}
  $archive=$paths[0].Substring(13); $hash=$hashes[0].Substring(15)
  if($archive -notmatch ('^'+[regex]::Escape($remote)+'/retrieval/result_[0-9_]+\.tar$') -or $hash -notmatch '^[a-f0-9]{64}$'){throw 'Invalid archive identity'}
  if(-not $DownloadDir){$DownloadDir=Join-Path (Split-Path $package -Parent) ($TaskName+'_retrieved_'+(Get-Date -Format yyyyMMdd_HHmmss))}
  $dest=[IO.Path]::GetFullPath($DownloadDir)
  if(Test-Path -LiteralPath $dest){throw 'New retrieval directory required'}
  New-Item -ItemType Directory -Path $dest|Out-Null
  $local=Join-Path $dest result.tar
  & scp -o BatchMode=yes "${alias}:$archive" $local
  if($LASTEXITCODE -ne 0 -or (Get-FileHash -LiteralPath $local).Hash.ToLowerInvariant() -ne $hash){throw 'Transport SHA failed'}
  $entries=& tar -tf $local
  if($LASTEXITCODE -ne 0){throw 'Cannot inspect archive'}
  foreach($entry in $entries){if($entry -match '(^/|:|\\|(^|/)\.\.(/|$))'){throw 'Unsafe archive member'}}
  & tar -xf $local -C $dest
  if($LASTEXITCODE -ne 0){throw 'Extraction failed'}
  $manifest=Get-Content -LiteralPath (Join-Path $dest retrieval/manifest.json) -Raw|ConvertFrom-Json
  foreach($row in $manifest){
   if($row.path -match '(^/|:|\\|(^|/)\.\.(/|$))'){throw 'Unsafe manifest path'}
   if((Get-FileHash -LiteralPath (Join-Path $dest $row.path)).Hash.ToLowerInvariant() -ne $row.sha256){throw ('Result SHA mismatch '+$row.path)}
  }
  Write-Output "FETCH_VERIFIED=$dest"
  Write-Output 'Transport verified; inspect exit.json and numerical checkpoints before claiming K0 PASS.'
 }
}
