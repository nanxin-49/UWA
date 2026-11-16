param(
    [Parameter(Mandatory=$true)][string]$ProjectRoot,
    [Parameter(Mandatory=$true)][string]$PackageDir,
    [ValidateSet('Prepare','Seal')][string]$Action='Prepare',
    [ValidatePattern('^r1_bie_manual_c4_policy_v2_[A-Za-z0-9_-]+$')]
    [string]$RemoteTaskName='r1_bie_manual_c4_policy_v2_20261116_01'
)
$ErrorActionPreference='Stop'
$project=(Resolve-Path -LiteralPath $ProjectRoot).Path
$package=[IO.Path]::GetFullPath($PackageDir)
if(-not $package.StartsWith($project+[IO.Path]::DirectorySeparatorChar)){throw 'Package must stay in project'}
function Sha([string]$p){(Get-FileHash -LiteralPath $p -Algorithm SHA256).Hash.ToLowerInvariant()}
function WriteJson([string]$p,$v){[IO.File]::WriteAllText($p,($v|ConvertTo-Json -Depth 30),[Text.UTF8Encoding]::new($false))}
if($Action -eq 'Prepare'){
    if(Test-Path -LiteralPath $package){throw 'New package directory required; never replace an existing package'}
    $source=Join-Path $project 'results/comparisons/piston_2d_pm/r1_bie_manual_runtime_fix_20261116_02_retrieved_20261116_144235'
    $pin=Get-Content -Raw -LiteralPath (Join-Path $project 'results/validation/r1_bie_c4_continuation_plan/20261116_static_01/static_checks.json')|ConvertFrom-Json
    $lists=@(Get-ChildItem -LiteralPath (Join-Path $source retrieval) -Filter 'SHA256SUMS_*.txt' -File)
    if($lists.Count -ne 1 -or (Sha $lists[0].FullName) -ne $pin.retrieval_manifest_sha256 -or (Sha (Join-Path $source manifest.json)) -ne $pin.legacy_manifest_sha256){throw 'Origin manifest differs from approved plan evidence'}
    $old=Get-Content -Raw -LiteralPath (Join-Path $source manifest.json)|ConvertFrom-Json
    foreach($f in $old.files){if((Sha (Join-Path $source $f.path)) -ne $f.sha256){throw 'Original package SHA drift'}}
    $records=@(Get-Content -LiteralPath $lists[0].FullName)
    if($records.Count -ne 110){throw 'Expected 110 origin artifacts'}
    $seen=[Collections.Generic.HashSet[string]]::new([StringComparer]::Ordinal)
    foreach($row in $records){
        if($row -notmatch '^([a-f0-9]{64})  (.+)$'){throw 'Malformed origin SHA row'}
        $expected=$Matches[1]; $rel=$Matches[2]
        $path=[IO.Path]::GetFullPath((Join-Path $source $rel))
        if(-not $path.StartsWith($source+[IO.Path]::DirectorySeparatorChar) -or -not $seen.Add($rel)){throw 'Unsafe or duplicate origin path'}
        if((Sha $path) -ne $expected){throw ('Origin drift: '+$rel)}
    }
    $code=@('scripts/validation/run_r1_bie_c4_continuation.m','scripts/validation/validate_r1_bie_c4_preparation.m',
        'scripts/validation/support/migrate_r1_bie_checkpoints.m','scripts/validation/support/r1_bie_continuation_io.m',
        'scripts/validation/support/r1_bie_continuation_gates.m','scripts/validation/support/audit_r1_bie_continuation_final.m',
        'scripts/validation/support/launch_r1_bie_c4_continuation.sh','scripts/validation/support/manage_r1_bie_c4_continuation.ps1',
        'scripts/validation/support/status_r1_bie_c4_session.sh','scripts/validation/support/invoke_r1_bie_status_stdin.ps1',
        'scripts/validation/support/r1_bie_c4_preflight.py','scripts/validation/support/test_r1_bie_c4_package.py','scripts/validation/support/r1_bie_manual_io.m',
        'scripts/validation/support/audit_r1_bie_manual_stage.m','scripts/validation/support/solve_helmholtz_bie_halfplane_vertical.m',
        'scripts/validation/support/solve_helmholtz_bie_halfplane_cached_vertical.m','scripts/validation/support/closeout_s3a_pm_surface.m',
        'src/propagation/baffled_strip_piston_field_vertical.m','scripts/validation/support/piston_independent_angular_spectrum_vertical.m',
        'scripts/comparisons/run_piston_pm_r1_job.m','scripts/validation/audit_piston_pm_r1.m','scripts/reporting/render_piston_pm_r1.m',
        'scripts/validation/support/archive_r1_bie_manual.sh','scripts/validation/support/r1_bie_runtime.sh',
        'scripts/validation/support/r1_bie_runtime_json.py')
    foreach($f in $code){
        $current=Join-Path $project $f; if(-not(Test-Path -LiteralPath $current -PathType Leaf)){throw ('Missing code '+$f)}
        $legacy=Join-Path $source ('code/'+(Split-Path $f -Leaf))
        if(Test-Path -LiteralPath $legacy){if((Sha $current) -ne (Sha $legacy)){throw ('Frozen numerical/runtime dependency changed: '+$f)}}
    }
    New-Item -ItemType Directory -Path $package,(Join-Path $package code),(Join-Path $package origin_evidence),(Join-Path $package control)|Out-Null
    foreach($row in $records){
        $parts=$row -split '  ',2; $rel=$parts[1]; $dst=Join-Path $package ('origin_evidence/'+$rel)
        New-Item -ItemType Directory -Force -Path (Split-Path -Parent $dst)|Out-Null
        Copy-Item -LiteralPath (Join-Path $source $rel) -Destination $dst
        if((Sha $dst) -ne $parts[0]){throw 'Origin copy SHA mismatch'}
    }
    $listRel=[IO.Path]::GetRelativePath($source,$lists[0].FullName).Replace('\','/')
    $listDest=Join-Path $package ('origin_evidence/'+$listRel)
    New-Item -ItemType Directory -Force -Path (Split-Path -Parent $listDest)|Out-Null
    Copy-Item -LiteralPath $lists[0].FullName -Destination $listDest
    foreach($f in $code){Copy-Item -LiteralPath (Join-Path $project $f) -Destination (Join-Path $package code)}
    foreach($f in Get-ChildItem -LiteralPath (Join-Path $package code) -Filter '*.sh' -File){
        # Copied unchanged legacy scripts already use LF; normalize new Bash code only.
        $legacy=Join-Path $source ('code/'+$f.Name)
        if(-not(Test-Path -LiteralPath $legacy)){[IO.File]::WriteAllText($f.FullName,[IO.File]::ReadAllText($f.FullName).Replace("`r`n","`n"),[Text.UTF8Encoding]::new($false))}
    }
    Copy-Item -LiteralPath (Join-Path $source common.mat) -Destination (Join-Path $package common.mat)
    Copy-Item -LiteralPath (Join-Path $source input) -Destination (Join-Path $package input) -Recurse
    Copy-Item -LiteralPath (Join-Path $source 'control/initial_available_bytes.txt') -Destination (Join-Path $package 'control/initial_available_bytes.txt')
    $plan='reports/r1_bie_nonmonotonic_policy_c4_plan_2026-11-16.md'
    Copy-Item -LiteralPath (Join-Path $project $plan) -Destination (Join-Path $package approved_plan.md)
    WriteJson (Join-Path $package migration_sources.json) ([ordered]@{legacy_manifest_sha256=$pin.legacy_manifest_sha256;retrieval_manifest_path=$listRel;retrieval_manifest_sha256=$pin.retrieval_manifest_sha256;source_run=(Split-Path $source -Leaf);original_session_exit_code=1;original_elapsed_s=304;historical_completed_solves=3})
    WriteJson (Join-Path $package approval.json) ([ordered]@{implementation_approved=$true;date='2026-11-16';policy_version='r1-bie-reference-v2';source_plan=$plan;source_plan_sha256=(Sha (Join-Path $project $plan));C4_started_by_agent=$false;user_scope='Approved v2 policy, zero-solve C0-C3 migration, C3 saved-density audit, local validation and package preparation only. No SSH, no C4 start, no Git.'})
    $gates=[ordered]@{sensitivity=1e-4;linear=1e-10;boundary=1e-4;reconstruction=1e-11;geometry=1e-12;source=1e-7;PS=1e-4;PS_step=1e-10;FullFK=1e-5;resolution_factor=5;phase_mask_fraction=.01}
    WriteJson (Join-Path $package policy.json) ([ordered]@{schema_version=2;version='r1-bie-reference-v2';grid_monotonicity_mode='diagnostic';reference_stage=3;pair_names=@('grid_10_12','grid_12_14','self_64_96','window_60_70');gates=$gates;legacy_policy='Strict spatial monotonicity failed; preserve original FAIL and complete diagnostic';approval_path='approval.json'})
    $files=@(Get-ChildItem -LiteralPath (Join-Path $package code),(Join-Path $package origin_evidence),(Join-Path $package input),(Join-Path $package control) -Recurse -File | ForEach-Object {[ordered]@{path=[IO.Path]::GetRelativePath($package,$_.FullName).Replace('\','/');sha256=(Sha $_.FullName)}})
    foreach($name in @('common.mat','policy.json','approval.json','migration_sources.json','approved_plan.md')){$files+=([ordered]@{path=$name;sha256=(Sha (Join-Path $package $name))})}
    WriteJson (Join-Path $package manifest.json) ([ordered]@{schema_version=2;scope='R1_BIE_C4_CONTINUATION';remote_task_name=$RemoteTaskName;allowed_new_solve_stage=4;single_process=$true;single_thread=$true;parpool=$false;RAM_reserve_fraction=.3;automatic_retry=$false;policy_sha256=(Sha (Join-Path $package policy.json));files=$files})
    & python (Join-Path $package 'code/r1_bie_c4_preflight.py') inputs $package
    if($LASTEXITCODE -ne 0){throw 'Input provenance verification failed'}
    Write-Output ('PREPARED_INPUTS_ONLY '+$package+'; migration and MATLAB acceptance still mandatory; no SSH/solve')
} else {
    & python (Join-Path $package 'code/r1_bie_c4_preflight.py') ready $package
    if($LASTEXITCODE -ne 0){throw 'Not ready: migration and local MATLAB acceptance must PASS before sealing/upload/start'}
    foreach($name in @('shell_receipt.json','integrity_receipt.json')){
        $proof=Get-Content -Raw -LiteralPath (Join-Path $package ('preparation_validation/'+$name))|ConvertFrom-Json
        if($proof.status -ne 'PASS' -or $proof.failed_checks -ne 0 -or $proof.new_acoustic_solves -ne 0){throw ('Incomplete local acceptance: '+$name)}
    }
    if(Test-Path -LiteralPath (Join-Path $package input_package.tar)){throw 'Archive already sealed; use a new package for changes'}
    $items=@(Get-ChildItem -LiteralPath $package -File -Recurse | Where-Object {$_.Name -notin @('package.sha256','input_package.tar','input_package.tar.sha256')})
    $shaLines=$items|Sort-Object FullName|ForEach-Object {(Sha $_.FullName)+'  '+[IO.Path]::GetRelativePath($package,$_.FullName).Replace('\','/')}
    [IO.File]::WriteAllText((Join-Path $package package.sha256),([string]::Join("`n",$shaLines)+"`n"),[Text.UTF8Encoding]::new($false))
    $names=@(Get-ChildItem -LiteralPath $package | Select-Object -ExpandProperty Name)
    & tar -cf (Join-Path $package input_package.tar) -C $package @names
    if($LASTEXITCODE -ne 0){throw 'Input archive failed'}
    [IO.File]::WriteAllText((Join-Path $package input_package.tar.sha256),(Sha (Join-Path $package input_package.tar))+'  input_package.tar'+"`n",[Text.UTF8Encoding]::new($false))
    Write-Output ('READY_SEALED '+$package+'; no SSH/solve started')
}
