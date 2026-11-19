param(
    [Parameter(Mandatory)][ValidateSet('DryRun','Upload','Start','Status','Fetch')][string]$Action,
    [Parameter(Mandatory)][string]$PackageDir,
    [ValidatePattern('^d1b_scale_[A-Za-z0-9_-]{1,64}$')][string]$TaskName='d1b_scale_20261119_01',
    [string]$DownloadDir
)
$ErrorActionPreference='Stop'
if($PSVersionTable.PSVersion.Major -lt 7){throw 'Run with PowerShell 7 (pwsh)'}
. (Join-Path $PSScriptRoot 'invoke_r1_bie_status_stdin.ps1')
$package=$ExecutionContext.SessionState.Path.GetUnresolvedProviderPathFromPSPath($PackageDir)
$archive=Join-Path $package 'input_package.tar'
$runner=Join-Path $package 'code/d1b_scale_remote.py'
$sshHost='matlab-ubuntu-new'
$remote='/home/misc/cwx/'+$TaskName

function Assert-File([string]$Path) {
    if(-not (Test-Path -LiteralPath $Path -PathType Leaf)){throw "Required file not found: '$Path'"}
}
function Verify-Package {
    $hashFile=Join-Path $package 'input_package.tar.sha256'
    foreach($file in @($runner,$archive,$hashFile,(Join-Path $package 'manifest.json'),(Join-Path $package 'config.json'))){Assert-File $file}
    & python -B $runner verify $package
    if($LASTEXITCODE -ne 0){throw 'Frozen scale pilot manifest verification failed'}
    $script:inputHash=(Get-Content -LiteralPath $hashFile -Raw).Trim().Split(' ')[0]
    if($inputHash -cnotmatch '^[a-f0-9]{64}$' -or
        (Get-FileHash -Algorithm SHA256 -LiteralPath $archive).Hash.ToLowerInvariant() -cne $inputHash){
        throw 'Frozen scale pilot archive SHA mismatch'
    }
}
function Invoke-Remote([string]$ScriptText) {
    $receipt=Invoke-R1StatusStdin -Executable ssh -ArgumentList @('-T','-o','BatchMode=yes','-o','ConnectTimeout=10',$sshHost,'bash','-s','--') -InputText ("set -euo pipefail`n"+$ScriptText+"`n")
    if($receipt.exit_code -ne 0){throw "SSH exit $($receipt.exit_code)`n$($receipt.stderr)`n$($receipt.stdout)"}
    $receipt.stdout
}

switch($Action){
    'DryRun' {
        Verify-Package
        Write-Output "DRY_RUN_PASS $remote $inputHash"
    }
    'Upload' {
        Verify-Package
        Invoke-Remote "umask 077; test -d /home/misc/cwx; mkdir '$remote'"
        & scp -o BatchMode=yes -o ConnectTimeout=10 $archive "${sshHost}:$remote/input_package.tar"
        if($LASTEXITCODE -ne 0){throw 'Upload failed. Partial directory retained.'}
        Invoke-Remote @"
cd '$remote'
printf '%s  input_package.tar\n' '$inputHash' | sha256sum -c -
tar --no-same-owner --no-same-permissions -xf input_package.tar
printf '%s  input_package.tar\n' '$inputHash' > input_package.tar.sha256
python3 -B code/d1b_scale_remote.py verify '$remote'
"@
    }
    'Start' {
        Verify-Package
        Invoke-Remote @"
cd '$remote'
printf '%s  input_package.tar\n' '$inputHash' | sha256sum -c -
python3 -B code/d1b_scale_remote.py start '$remote'
"@
    }
    'Status' { Invoke-Remote "python3 -B '$remote/code/d1b_scale_remote.py' status '$remote'" }
    'Fetch' {
        Assert-File $runner
        $receipt=Invoke-Remote "python3 -B '$remote/code/d1b_scale_remote.py' fetch '$remote'"
        $paths=@($receipt -split "`n" | Where-Object {$_ -like 'ARCHIVE_PATH=*'})
        $hashes=@($receipt -split "`n" | Where-Object {$_ -like 'ARCHIVE_SHA256=*'})
        if($paths.Count -ne 1 -or $hashes.Count -ne 1){throw 'Expected one fetch receipt'}
        $source=$paths[0].Substring(13).Trim(); $hash=$hashes[0].Substring(15).Trim()
        if($source -cnotmatch ('^'+[regex]::Escape($remote)+'/retrieval/result_[0-9]+\.tar$') -or
            $hash -cnotmatch '^[a-f0-9]{64}$'){throw 'Unexpected fetch archive identity'}
        if(-not $DownloadDir){
            $DownloadDir=Join-Path (Split-Path $package -Parent) ($TaskName+'_retrieved_'+(Get-Date -Format yyyyMMdd_HHmmss))
        }
        $destination=$ExecutionContext.SessionState.Path.GetUnresolvedProviderPathFromPSPath($DownloadDir)
        $transport=$destination+'_transport'
        if((Test-Path -LiteralPath $destination) -or (Test-Path -LiteralPath $transport)){throw 'Use a new retrieval path'}
        New-Item -ItemType Directory -Path $transport | Out-Null
        $file=Join-Path $transport 'result.tar'
        & scp -o BatchMode=yes -o ConnectTimeout=10 "${sshHost}:$source" $file
        if($LASTEXITCODE -ne 0){throw 'Fetch transfer failed'}
        [IO.File]::WriteAllText((Join-Path $transport 'receipt.txt'),$receipt,[Text.UTF8Encoding]::new($false))
        & python -B $runner extract $file $destination $hash
        if($LASTEXITCODE -ne 0){throw 'Fetch SHA or safe extraction failed'}
    }
}
