param(
    [ValidateSet('ShowCommands','UploadAndStart','Status','Continue','Fetch')][string]$Action='ShowCommands',
    [Parameter(Mandatory=$true)][string]$PackageDir,
    [string]$DownloadDir
)
$ErrorActionPreference='Stop'
$package=(Resolve-Path -LiteralPath $PackageDir).Path
$manifest=Get-Content -Raw -LiteralPath (Join-Path $package manifest.json)|ConvertFrom-Json
if($manifest.schema_version -ne 2 -or $manifest.scope -ne 'R1_BIE_C4_CONTINUATION'){throw 'Not a C4 v2 package'}
$taskName=[string]$manifest.remote_task_name
if($taskName -notmatch '^r1_bie_manual_c4_policy_v2_[A-Za-z0-9_-]+$'){throw 'Invalid task name'}
$remoteRoot='/home/misc/cwx/'+$taskName
$alias='matlab-ubuntu-new'
$launcher=$remoteRoot+'/code/launch_r1_bie_c4_continuation.sh'
function Invoke-CheckedSsh([string]$Command){
    & ssh -o BatchMode=yes $alias $Command
    if($LASTEXITCODE -ne 0){throw ('SSH failed with exit '+$LASTEXITCODE)}
}
function Assert-Ready {
    & python (Join-Path $package 'code/r1_bie_c4_preflight.py') ready $package
    if($LASTEXITCODE -ne 0){throw 'Local migration and MATLAB acceptance are not PASS; refusing upload/start'}
}
switch($Action){
    'ShowCommands' {
        Write-Output ('REMOTE_ROOT='+$remoteRoot)
        Write-Output ('START: bash '+$launcher+' start '+$remoteRoot)
        Write-Output ('STATUS: bash '+$launcher+' status '+$remoteRoot)
        Write-Output 'Only C4 may solve. This action makes no network connection.'
    }
    'UploadAndStart' {
        Assert-Ready
        $archive=Join-Path $package input_package.tar
        $expected=(Get-Content -Raw -LiteralPath (Join-Path $package input_package.tar.sha256)).Trim().Split(' ')[0]
        if($expected -notmatch '^[a-f0-9]{64}$' -or (Get-FileHash -LiteralPath $archive -Algorithm SHA256).Hash.ToLowerInvariant() -ne $expected){throw 'Archive SHA mismatch'}
        Invoke-CheckedSsh "test ! -e '$remoteRoot' && mkdir '$remoteRoot'"
        & scp -o BatchMode=yes $archive "${alias}:$remoteRoot/input_package.tar"
        if($LASTEXITCODE -ne 0){throw 'Upload failed; partial remote upload retained'}
        Invoke-CheckedSsh "cd '$remoteRoot' && printf '%s  input_package.tar\n' '$expected' | sha256sum -c - && tar -xf input_package.tar && bash '$launcher' start '$remoteRoot'"
    }
    'Status' {
        # Send a read-only presenter; never patch a sealed/running task package.
        $presenter=Join-Path $PSScriptRoot 'status_r1_bie_c4_session.sh'
        $script=[IO.File]::ReadAllText($presenter).Replace("`r",'')
        $script+="`nr1_c4_status_main " + '$1' + "`n"
        # Keep all quote-rich script bytes off Windows native command lines.
        . (Join-Path $PSScriptRoot 'invoke_r1_bie_status_stdin.ps1')
        $ssh=(Get-Command ssh.exe -CommandType Application).Source
        $result=Invoke-R1StatusStdin -Executable $ssh -ArgumentList @('-o','BatchMode=yes',$alias,('bash -s -- '+$remoteRoot)) -InputText $script
        if($result.stdout){Write-Output $result.stdout.TrimEnd()}
        if($result.stderr){[Console]::Error.Write($result.stderr)}
        if($result.exit_code -ne 0){throw ('Read-only SSH status failed with exit '+$result.exit_code)}
    }
    'Continue' {Assert-Ready; Invoke-CheckedSsh "bash '$launcher' start '$remoteRoot'"}
    'Fetch' {
        $output=& ssh -o BatchMode=yes $alias "bash '$remoteRoot/code/archive_r1_bie_manual.sh' '$remoteRoot'"
        if($LASTEXITCODE -ne 0){throw 'Stopped-task retrieval snapshot failed/refused'}
        $line=@($output|Where-Object {$_ -like 'ARCHIVE_PATH=*'})
        if($line.Count -ne 1){throw 'Missing unique archive receipt'}
        $remoteArchive=$line[0].Substring(13)
        if($remoteArchive -notmatch ('^'+[regex]::Escape($remoteRoot)+'/retrieval/result_[A-Za-z0-9_]+\.tar$')){throw 'Unexpected remote archive path'}
        if(-not $DownloadDir){$DownloadDir=Join-Path (Split-Path -Parent $package) ($taskName+'_retrieved_'+(Get-Date -Format yyyyMMdd_HHmmss))}
        $download=[IO.Path]::GetFullPath($DownloadDir)
        if(Test-Path -LiteralPath $download){throw 'New download directory required'}
        New-Item -ItemType Directory -Path $download|Out-Null
        $localArchive=Join-Path $download result.tar
        & scp -o BatchMode=yes "${alias}:$remoteArchive" $localArchive
        if($LASTEXITCODE -ne 0){throw 'Retrieval SCP failed'}
        $entries=& tar -tf $localArchive
        if($LASTEXITCODE -ne 0){throw 'Cannot inspect retrieved archive'}
        foreach($entry in $entries){if($entry -match '(^/|:|(^|/)\.\.(/|$)|\\)'){throw 'Unsafe archive member'}}
        & tar -xf $localArchive -C $download
        if($LASTEXITCODE -ne 0){throw 'Archive extraction failed'}
        $list=@(Get-ChildItem -LiteralPath (Join-Path $download retrieval) -Filter 'SHA256SUMS_*.txt' -File)
        if($list.Count -ne 1){throw 'Expected one top-level retrieval SHA list'}
        $count=0; $seen=[Collections.Generic.HashSet[string]]::new([StringComparer]::Ordinal)
        foreach($row in Get-Content -LiteralPath $list[0].FullName){
            if($row -notmatch '^([a-f0-9]{64})  (.+)$'){throw 'Malformed retrieval hash'}
            $expected=$Matches[1];$relative=$Matches[2]
            $target=[IO.Path]::GetFullPath((Join-Path $download $relative))
            if(-not $target.StartsWith($download+[IO.Path]::DirectorySeparatorChar) -or -not $seen.Add($relative)){throw 'Unsafe or duplicate retrieval target'}
            if((Get-FileHash -LiteralPath $target -Algorithm SHA256).Hash.ToLowerInvariant() -ne $expected){throw ('Retrieved SHA mismatch '+$relative)}
            $count++
        }
        Write-Output ('RETRIEVED_SHA_PASS files='+$count+' path='+$download)
    }
}
