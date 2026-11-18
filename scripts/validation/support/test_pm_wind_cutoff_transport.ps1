param([Parameter(Mandatory=$true)][string]$PackageDir,[Parameter(Mandatory=$true)][string]$OutputDir)
$ErrorActionPreference='Stop'
if(Test-Path -LiteralPath $OutputDir){throw 'New test directory required'}
New-Item -ItemType Directory -Path $OutputDir|Out-Null
Start-Transcript -LiteralPath (Join-Path $OutputDir test.log)|Out-Null
try {
 $manager=Join-Path $PSScriptRoot manage_pm_wind_cutoff_k0.ps1
 $global:pmK0MockCalls=@(); $global:pmK0MockArchive=''; $global:pmK0MockHash=''
 function ssh {
  $command=$args[-1]; $global:pmK0MockCalls+=('ssh '+$command); $global:LASTEXITCODE=0
  if($command -like '* fetch *'){
   Write-Output 'ARCHIVE_PATH=/home/misc/cwx/pm_wind_cutoff_k0_test/retrieval/result_20261117_01.tar'
   Write-Output ('ARCHIVE_SHA256='+$global:pmK0MockHash)
  } else {Write-Output 'MOCK_ONLY_NO_CONNECTION'}
 }
 function scp {
  $global:pmK0MockCalls+=('scp '+($args -join ' ')); $global:LASTEXITCODE=0
  if($args[-2] -like '*:*/retrieval/*.tar'){Copy-Item -LiteralPath $global:pmK0MockArchive -Destination $args[-1]}
 }
 $tests=0
 & $manager -Action ShowCommands -PackageDir $PackageDir -TaskName pm_wind_cutoff_k0_test
 if($global:pmK0MockCalls.Count -ne 0){throw 'ShowCommands attempted network'}; $tests++
 & $manager -Action Upload -PackageDir $PackageDir -TaskName pm_wind_cutoff_k0_test
 if($global:pmK0MockCalls.Count -ne 3){throw 'Upload mock sequence'}; $tests++
 & $manager -Action Start -PackageDir $PackageDir -TaskName pm_wind_cutoff_k0_test
 & $manager -Action Status -PackageDir $PackageDir -TaskName pm_wind_cutoff_k0_test
 if($global:pmK0MockCalls.Count -ne 5){throw 'Start/Status mock sequence'}; $tests+=2
 $source=Join-Path $OutputDir snapshot; New-Item -ItemType Directory -Path $source,(Join-Path $source retrieval)|Out-Null
 [IO.File]::WriteAllText((Join-Path $source payload.txt),'mock result',[Text.UTF8Encoding]::new($false))
 $hash=(Get-FileHash -LiteralPath (Join-Path $source payload.txt)).Hash.ToLowerInvariant()
 $json=ConvertTo-Json -InputObject @(@{path='payload.txt';sha256=$hash})
 [IO.File]::WriteAllText((Join-Path $source retrieval/manifest.json),$json,[Text.UTF8Encoding]::new($false))
 $global:pmK0MockArchive=Join-Path $OutputDir result.tar
 & tar -cf $global:pmK0MockArchive -C $source payload.txt retrieval/manifest.json
 if($LASTEXITCODE -ne 0){throw 'Mock archive build failed'}
 $global:pmK0MockHash=(Get-FileHash -LiteralPath $global:pmK0MockArchive).Hash.ToLowerInvariant()
 & $manager -Action Fetch -PackageDir $PackageDir -TaskName pm_wind_cutoff_k0_test -DownloadDir (Join-Path $OutputDir fetched)
 if(-not (Test-Path -LiteralPath (Join-Path $OutputDir fetched/payload.txt))){throw 'Fetch mock not verified'}; $tests++
 $global:pmK0MockHash='0'*64; $rejected=$false
 try {& $manager -Action Fetch -PackageDir $PackageDir -TaskName pm_wind_cutoff_k0_test -DownloadDir (Join-Path $OutputDir rejected)} catch {$rejected=$true}
 if(-not $rejected){throw 'Bad transport SHA accepted'}; $tests++
 $before=$global:pmK0MockCalls.Count; $rejected=$false
 try {& $manager -Action Start -PackageDir $PackageDir -TaskName '../../elsewhere'} catch {$rejected=$true}
 if(-not $rejected -or $global:pmK0MockCalls.Count -ne $before){throw 'Unsafe task accepted'}; $tests++
 $r=@{status='MOCK_TRANSPORT_PASS';tests=$tests;real_SSH_connections=0;real_MATLAB_launches=0;mock_calls=$global:pmK0MockCalls}
 $r|ConvertTo-Json -Depth 8|Set-Content -LiteralPath (Join-Path $OutputDir receipt.json) -Encoding utf8
 Write-Output "MOCK_TRANSPORT_PASS tests=$tests real_SSH_connections=0"
} finally {Stop-Transcript|Out-Null}
