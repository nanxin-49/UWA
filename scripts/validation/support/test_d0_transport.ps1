param([Parameter(Mandatory=$true)][string]$PackageDir,[Parameter(Mandatory=$true)][string]$OutputDir,[Parameter(Mandatory=$true)][string]$FixtureTar)
$ErrorActionPreference='Stop'
if(Test-Path -LiteralPath $OutputDir){throw 'New mock test directory required'}
New-Item -ItemType Directory -Path $OutputDir|Out-Null
$utf8=[Text.UTF8Encoding]::new($false)
$tool=Join-Path $OutputDir manage_d0_incident_depth.ps1
Copy-Item -LiteralPath (Join-Path $PSScriptRoot manage_d0_incident_depth.ps1) -Destination $tool
# TEST ONLY replacement next to the copied manager. Formal package is untouched.
$mock=@'
function Invoke-R1StatusStdin {
 param($Executable,$ArgumentList,$InputText)
 if($Executable -ne 'ssh' -or $ArgumentList -notcontains 'BatchMode=yes' -or $ArgumentList -notcontains 'matlab-ubuntu-new'){throw 'Unexpected transport'}
 $global:d0Scripts.Add($InputText)
 if($global:d0Fail){return [pscustomobject]@{exit_code=75;stdout='TEST_ONLY failure';stderr='serial lock busy'}}
 if($InputText -match " fetch '"){$stdout="ARCHIVE_PATH=/home/misc/cwx/d0_incident_depth_TEST_ONLY/retrieval/result_20261118_000000.tar`nARCHIVE_SHA256=$global:d0Hash`n"}
 else {$stdout='TEST_ONLY remote mock'}
 return [pscustomobject]@{exit_code=0;stdout=$stdout;stderr=''}
}
'@
[IO.File]::WriteAllText((Join-Path $OutputDir invoke_r1_bie_status_stdin.ps1),$mock,$utf8)
$global:d0Scripts=[Collections.Generic.List[string]]::new(); $global:d0Copies=0; $global:d0Fail=$false
$global:d0Hash=(Get-FileHash -LiteralPath $FixtureTar).Hash.ToLowerInvariant(); $global:d0Fixture=$FixtureTar
function global:scp {
 $Values=$args
 $global:d0Copies++
 if($Values[-2] -like '*:*/retrieval/*'){Copy-Item -LiteralPath $global:d0Fixture -Destination $Values[-1]}
 $global:LASTEXITCODE=0
}
$checks=[Collections.Generic.List[object]]::new()
function Assert-Mock($name,$value){if(-not $value){throw $name};$checks.Add(@{name=$name;passed=$true})}
try {
 & $tool -Action DryRun -PackageDir $PackageDir -TaskName d0_incident_depth_TEST_ONLY
 Assert-Mock 'DryRun zero transport' ($global:d0Scripts.Count -eq 0 -and $global:d0Copies -eq 0)
 & $tool -Action UploadAndStart -PackageDir $PackageDir -TaskName d0_incident_depth_TEST_ONLY
 Assert-Mock 'UploadAndStart exactly three scripts and one copy' ($global:d0Scripts.Count -eq 3 -and $global:d0Copies -eq 1)
 Assert-Mock 'Non-overwrite remote creation' ($global:d0Scripts[0].Contains('test ! -e'))
 Assert-Mock 'Remote archive SHA before extraction' ($global:d0Scripts[1].Contains('sha256sum -c -'))
 Assert-Mock 'One manual start' ($global:d0Scripts[2].Contains(" start '/home/misc/cwx/d0_incident_depth_TEST_ONLY'"))
 & $tool -Action Status -PackageDir $PackageDir -TaskName d0_incident_depth_TEST_ONLY
 Assert-Mock 'Status read-only command' ($global:d0Scripts[3].Contains(" status '"))
 $global:d0Fail=$true; $before=$global:d0Scripts.Count; $caught=$false
 try {& $tool -Action Start -PackageDir $PackageDir -TaskName d0_incident_depth_TEST_ONLY} catch {$caught=$_.Exception.Message.Contains('SSH exit=75')}
 Assert-Mock 'Start failure no retry' ($caught -and $global:d0Scripts.Count -eq $before+1)
 $global:d0Fail=$false
 $dest=Join-Path $OutputDir fetched_partial
 & $tool -Action Fetch -PackageDir $PackageDir -TaskName d0_incident_depth_TEST_ONLY -DownloadDir $dest
 Assert-Mock 'Fetch archive receipt preserved' (Test-Path -LiteralPath ($dest+'_transport/remote_archive_receipt.txt'))
 Assert-Mock 'Fetch verifies and retains failed exit' ((Get-Content -LiteralPath (Join-Path $dest sessions/test/exit.json) -Raw|ConvertFrom-Json).matlab_exit_code -eq 1)
 $before=$global:d0Scripts.Count; $caught=$false
 try {& $tool -Action Status -PackageDir $PackageDir -TaskName 'd0_incident_depth_bad;echo'} catch {$caught=$true}
 Assert-Mock 'Unsafe task rejected before transport' ($caught -and $global:d0Scripts.Count -eq $before)
 [IO.File]::WriteAllText((Join-Path $OutputDir scripts.json),($global:d0Scripts|ConvertTo-Json),$utf8)
 $receipt=@{status='TRANSPORT_MOCK_PASS';checks=$checks;network_calls=0;MATLAB_processes=0;mock_only=$true}
 [IO.File]::WriteAllText((Join-Path $OutputDir receipt.json),($receipt|ConvertTo-Json -Depth 10),$utf8)
 Write-Output ('MOCK_PASS checks='+$checks.Count+' actual_network_calls=0')
} finally {Remove-Item Function:/scp}
