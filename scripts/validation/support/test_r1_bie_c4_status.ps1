param(
    [Parameter(Mandatory=$true)][string]$ProjectRoot,
    [Parameter(Mandatory=$true)][string]$PackageDir,
    [Parameter(Mandatory=$true)][string]$OutputDir,
    [string]$Bash='D:\Git\bin\bash.exe'
)
$ErrorActionPreference='Stop'
if(Test-Path -LiteralPath $OutputDir){throw 'New output directory required'}
New-Item -ItemType Directory -Path $OutputDir|Out-Null
$out=(Resolve-Path -LiteralPath $OutputDir).Path
$helper=Join-Path $ProjectRoot 'scripts/validation/support/status_r1_bie_c4_session.sh'
$manager=Join-Path $ProjectRoot 'scripts/validation/support/manage_r1_bie_c4_continuation.ps1'
$transport=Join-Path $ProjectRoot 'scripts/validation/support/invoke_r1_bie_status_stdin.ps1'
$passed=0
foreach($file in @($manager,$transport,$PSCommandPath)){
    $tokens=$null;$errors=$null
    [Management.Automation.Language.Parser]::ParseFile($file,[ref]$tokens,[ref]$errors)|Out-Null
    if($errors.Count){throw 'PowerShell parsing failed'}
    $passed++
}
. $transport
$sample="apostrophe's `"quotes`" (round) 中文"
[IO.File]::WriteAllText((Join-Path $out sample.json),(@{key=$sample}|ConvertTo-Json),[Text.UTF8Encoding]::new($false))
$jsonPath=(Join-Path $out sample.json).Replace('\','/')
if($jsonPath -match '^([A-Za-z]):/(.*)$'){$jsonPath='/'+$Matches[1].ToLowerInvariant()+'/'+$Matches[2]}
$probe=[IO.File]::ReadAllText($helper).Replace("`r",'')+"`n"+'python3() { python -X utf8 "$@"; }; r1_status_json_field "$1" key'+"`n"
$result=Invoke-R1StatusStdin -Executable $Bash -ArgumentList @('-s','--',$jsonPath) -InputText $probe
[IO.File]::WriteAllText((Join-Path $out native_quote_probe_stdout.log),$result.stdout,[Text.UTF8Encoding]::new($false))
[IO.File]::WriteAllText((Join-Path $out native_quote_probe_stderr.log),$result.stderr,[Text.UTF8Encoding]::new($false))
if($result.exit_code -ne 0 -or $result.stdout.TrimEnd() -cne $sample){throw 'Native stdin quote/unicode round-trip failed'}
$passed++
$large=('# padding to exceed command-line size '+('x'*20000))+"`n"+$probe
$result=Invoke-R1StatusStdin -Executable $Bash -ArgumentList @('-s','--',$jsonPath) -InputText $large
if($result.exit_code -ne 0 -or $result.stdout.TrimEnd() -cne $sample){throw 'Large stdin payload failed'}
$passed++
$result=Invoke-R1StatusStdin -Executable $Bash -ArgumentList @('-s') -InputText "printf 'fixture error\n' >&2; exit 7`n"
if($result.exit_code -ne 7 -or -not $result.stderr.Contains('fixture error')){throw 'Real native exit/stderr lost'}
$passed++
$manifest=Get-Content -Raw -LiteralPath (Join-Path $PackageDir manifest.json)|ConvertFrom-Json
if($manifest.remote_task_name -notmatch '^r1_bie_manual_c4_policy_v2_[A-Za-z0-9_-]+$'){throw 'Unexpected task root'}
$managerSource=[IO.File]::ReadAllText($manager)
if(-not $managerSource.Contains('Invoke-R1StatusStdin -Executable $ssh') -or $managerSource.Contains('$quotedScript') -or $managerSource.Contains('ToBase64String')){throw 'Manager did not select stdin transport'}
$passed++
$receipt=[ordered]@{status='PASS';failed_checks=0;passed_checks=$passed;transport='PowerShell .NET native stdin -> Bash/Python; no script text in command line';apostrophe_double_quote_parentheses_unicode_roundtrip=$true;large_payload_bytes=([Text.Encoding]::UTF8.GetByteCount($large));real_exit_and_stderr_preserved=$true;SSH_connections=0;MATLAB_starts=0}
[IO.File]::WriteAllText((Join-Path $out receipt.json),($receipt|ConvertTo-Json -Depth 5),[Text.UTF8Encoding]::new($false))
Write-Output ('PASS transport_checks='+$passed+' SSH=0 MATLAB=0')