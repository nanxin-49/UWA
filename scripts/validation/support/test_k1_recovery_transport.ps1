param([Parameter(Mandatory=$true)][string]$PackageDir,[Parameter(Mandatory=$true)][string]$OutputDir)
$ErrorActionPreference='Stop'
if(Test-Path -LiteralPath $OutputDir){throw 'New transport test directory required'}
New-Item -ItemType Directory -Path $OutputDir|Out-Null; Start-Transcript -LiteralPath (Join-Path $OutputDir test.log)|Out-Null
try{
 $testPackage=Join-Path $OutputDir TEST_ONLY_package; Copy-Item -LiteralPath $PackageDir -Destination $testPackage -Recurse
 $archive=Join-Path $testPackage input_package.tar
 if(-not(Test-Path -LiteralPath $archive)){
  & tar -cf $archive -C $testPackage code cases anchor lineage provenance config.json manifest.json
  if($LASTEXITCODE -ne 0){throw 'Fixture archive failed'}
  [IO.File]::WriteAllText((Join-Path $testPackage input_package.tar.sha256),((Get-FileHash -LiteralPath $archive).Hash.ToLowerInvariant()+'  input_package.tar'),[Text.UTF8Encoding]::new($false))
 }
 $manager=Join-Path $testPackage code/manage_k1_pm_recovery.ps1
 $cfg=Get-Content -LiteralPath (Join-Path $testPackage config.json) -Raw|ConvertFrom-Json
 $global:k1MockCalls=@(); $global:k1MockRemote='/home/misc/cwx/'+$cfg.remote_task_name
 $global:k1MockArchive=''; $global:k1MockHash=''; $tests=0
 function ssh {
  $command=$args[-1]; $global:k1MockCalls+=('ssh '+$command); $global:LASTEXITCODE=0
  if($command -like '* fetch *'){
   Write-Output ('ARCHIVE_PATH='+$global:k1MockRemote+'/retrieval/result_20261118_01.tar')
   Write-Output ('ARCHIVE_SHA256='+$global:k1MockHash)
  } else {Write-Output 'MOCK_ONLY_NO_SSH_CONNECTION'}
 }
 function scp {
  $global:k1MockCalls+=('scp '+($args -join ' ')); $global:LASTEXITCODE=0
  if($args[-2] -like '*:*/retrieval/*.tar'){Copy-Item -LiteralPath $global:k1MockArchive -Destination $args[-1]}
 }
 & $manager -Action DryRun -PackageDir $testPackage
 & $manager -Action ShowCommands -PackageDir $testPackage
 if($global:k1MockCalls.Count){throw 'Offline commands attempted network'}; $tests+=2
 & $manager -Action Upload -PackageDir $testPackage -SshConfig C:/TEST_ONLY_CONFIG_NOT_READ
 if($global:k1MockCalls.Count -ne 3){throw 'Upload command sequence'}; $tests++
 & $manager -Action Start -PackageDir $testPackage
 & $manager -Action Status -PackageDir $testPackage
 if($global:k1MockCalls.Count -ne 5){throw 'Start/Status sequence'}; $tests+=2
 $snapshot=Join-Path $OutputDir snapshot; New-Item -ItemType Directory -Path $snapshot,(Join-Path $snapshot retrieval)|Out-Null
 [IO.File]::WriteAllText((Join-Path $snapshot config.json),'TEST_ONLY_TRANSPORT_PAYLOAD',[Text.UTF8Encoding]::new($false))
 $rows=@(@{path='config.json';sha256=(Get-FileHash -LiteralPath (Join-Path $snapshot config.json)).Hash.ToLowerInvariant()})
 [IO.File]::WriteAllText((Join-Path $snapshot retrieval/manifest.json),(ConvertTo-Json -InputObject $rows),[Text.UTF8Encoding]::new($false))
 $global:k1MockArchive=Join-Path $OutputDir result.tar
 & tar -cf $global:k1MockArchive -C $snapshot config.json retrieval/manifest.json
 if($LASTEXITCODE -ne 0){throw 'Mock retrieval tar failed'}
 $global:k1MockHash=(Get-FileHash -LiteralPath $global:k1MockArchive).Hash.ToLowerInvariant()
 & $manager -Action Fetch -PackageDir $testPackage -DownloadDir (Join-Path $OutputDir fetched)
 $receipt=Get-Content -LiteralPath (Join-Path $OutputDir fetched/transport_receipt.json) -Raw|ConvertFrom-Json
 if(-not $receipt.transport_pass -or $receipt.expected_remote_archive_sha256 -ne $global:k1MockHash){throw 'External transport SHA not persisted'}; $tests++
 $global:k1MockHash='0'*64; $bad=$false
 try{& $manager -Action Fetch -PackageDir $testPackage -DownloadDir (Join-Path $OutputDir rejected)}catch{$bad=$true}
 if(-not $bad){throw 'Wrong transport SHA accepted'}; $tests++
 $before=$global:k1MockCalls.Count; $bad=$false
 try{& $manager -Action Start -PackageDir $testPackage -TaskName '../../elsewhere'}catch{$bad=$true}
 if(-not $bad -or $global:k1MockCalls.Count -ne $before){throw 'Invalid task before-network check'}; $tests++
 $result=@{status='K1_RECOVERY_MOCK_TRANSPORT_PASS';tests=$tests;real_SSH_connections=0;MATLAB_launches=0;mock_calls=$global:k1MockCalls}
 $result|ConvertTo-Json -Depth 12|Set-Content -LiteralPath (Join-Path $OutputDir receipt.json) -Encoding utf8
 Write-Output "K1_RECOVERY_MOCK_TRANSPORT_PASS tests=$tests real_SSH_connections=0"
}finally{Stop-Transcript|Out-Null}
