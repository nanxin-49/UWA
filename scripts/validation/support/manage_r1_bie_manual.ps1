param(
    [ValidateSet('ShowCommands','UploadAndStart','Status','Resume','Fetch')]
    [string]$Action='ShowCommands',
    [Parameter(Mandatory=$true)][string]$PackageDir,
    [string]$DownloadDir
)
$ErrorActionPreference='Stop'
$package=(Resolve-Path -LiteralPath $PackageDir).Path
$manifest=Get-Content -LiteralPath (Join-Path $package 'manifest.json') -Raw | ConvertFrom-Json
$taskName=[string]$manifest.remote_task_name
if($taskName -notmatch '^r1_bie_manual_[A-Za-z0-9_-]+$'){throw 'Invalid task name'}
$remoteRoot='/home/misc/cwx/'+$taskName
$alias='matlab-ubuntu-new'
$launcher=$remoteRoot+'/code/launch_r1_bie_manual.sh'
function Invoke-CheckedSsh([string]$Command){
    & ssh -o BatchMode=yes $alias $Command
    if($LASTEXITCODE -ne 0){throw ('SSH command failed with exit '+$LASTEXITCODE)}
}
function Test-PackageFiles {
    foreach($f in $manifest.files){
        $target=[IO.Path]::GetFullPath((Join-Path $package $f.path))
        if(-not $target.StartsWith($package+[IO.Path]::DirectorySeparatorChar)){throw 'Manifest path escapes package'}
        if((Get-FileHash -LiteralPath $target -Algorithm SHA256).Hash.ToLowerInvariant() -ne $f.sha256){throw ('Package SHA mismatch: '+$f.path)}
    }
}
switch($Action){
    'ShowCommands' {
        Write-Output ('REMOTE_ROOT='+$remoteRoot)
        Write-Output ('START: bash '+$launcher+' start '+$remoteRoot)
        Write-Output ('STATUS: bash '+$launcher+' status '+$remoteRoot)
        Write-Output 'No SSH connection or computation was started.'
    }
    'UploadAndStart' {
        Test-PackageFiles
        $archive=Join-Path $package 'input_package.tar'
        if(-not(Test-Path -LiteralPath $archive -PathType Leaf)){throw 'Missing prepared input_package.tar'}
        $expected=(Get-Content -LiteralPath (Join-Path $package 'input_package.tar.sha256') -Raw).Trim().Split(' ')[0]
        if((Get-FileHash -LiteralPath $archive -Algorithm SHA256).Hash.ToLowerInvariant() -ne $expected){throw 'Input archive SHA mismatch'}
        Invoke-CheckedSsh "test ! -e '$remoteRoot' && mkdir '$remoteRoot'"
        & scp -o BatchMode=yes $archive "${alias}:$remoteRoot/input_package.tar"
        if($LASTEXITCODE -ne 0){throw 'SCP upload failed; remote partial upload retained'}
        Invoke-CheckedSsh "tar -xf '$remoteRoot/input_package.tar' -C '$remoteRoot' && bash '$launcher' start '$remoteRoot'"
    }
    'Status' { Invoke-CheckedSsh "bash '$launcher' status '$remoteRoot'" }
    'Resume' { Invoke-CheckedSsh "bash '$launcher' start '$remoteRoot'" }
    'Fetch' {
        $output=& ssh -o BatchMode=yes $alias "bash '$remoteRoot/code/archive_r1_bie_manual.sh' '$remoteRoot'"
        if($LASTEXITCODE -ne 0){throw 'Remote retrieval snapshot refused or failed'}
        $line=@($output | Where-Object {$_ -like 'ARCHIVE_PATH=*'})
        if($line.Count -ne 1){throw 'Missing archive receipt'}
        $remoteArchive=$line[0].Substring(13)
        if($remoteArchive -notmatch ('^'+[regex]::Escape($remoteRoot)+'/retrieval/result_[A-Za-z0-9_]+\.tar$')){throw 'Unexpected remote archive path'}
        if(-not $DownloadDir){
            $DownloadDir=Join-Path (Split-Path -Parent $package) ($taskName+'_retrieved_'+(Get-Date -Format 'yyyyMMdd_HHmmss'))
        }
        $download=[IO.Path]::GetFullPath($DownloadDir)
        if(Test-Path -LiteralPath $download){throw 'Download directory must not already exist'}
        New-Item -ItemType Directory -Path $download | Out-Null
        $localArchive=Join-Path $download 'result.tar'
        & scp -o BatchMode=yes "${alias}:$remoteArchive" $localArchive
        if($LASTEXITCODE -ne 0){throw 'SCP retrieval failed'}
        $entries=& tar -tf $localArchive
        if($LASTEXITCODE -ne 0){throw 'Cannot inspect archive'}
        foreach($entry in $entries){
            if($entry -match '(^/|:|(^|/)\.\.(/|$)|\\)'){throw 'Unsafe archive path'}
        }
        & tar -xf $localArchive -C $download
        if($LASTEXITCODE -ne 0){throw 'Cannot extract archive'}
        $sha=Get-ChildItem -LiteralPath (Join-Path $download 'retrieval') -Filter 'SHA256SUMS_*.txt' -File
        if(@($sha).Count -ne 1){throw 'Expected one retrieval SHA manifest'}
        $count=0
        foreach($row in Get-Content -LiteralPath $sha.FullName){
            if($row -notmatch '^([a-f0-9]{64})  (.+)$'){throw 'Malformed SHA row'}
            $expected=$Matches[1]; $relative=$Matches[2]
            $target=[IO.Path]::GetFullPath((Join-Path $download $relative))
            if(-not $target.StartsWith($download+[IO.Path]::DirectorySeparatorChar)){throw 'SHA target escapes retrieval directory'}
            if((Get-FileHash -LiteralPath $target -Algorithm SHA256).Hash.ToLowerInvariant() -ne $expected){throw ('Retrieved SHA mismatch '+$relative)}
            $count++
        }
        Write-Output ('RETRIEVED_SHA_PASS files='+$count+' path='+$download)
    }
}
