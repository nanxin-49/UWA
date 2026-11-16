param(
    [Parameter(Mandatory=$true)][string]$ProjectRoot,
    [string]$PackageDir,
    [ValidatePattern('^r1_bie_manual_[A-Za-z0-9_-]+$')]
    [string]$RemoteTaskName='r1_bie_manual_20261116_01'
)
$ErrorActionPreference='Stop'
$project=(Resolve-Path -LiteralPath $ProjectRoot).Path
if(-not $PackageDir){$PackageDir=Join-Path $project 'results/comparisons/piston_2d_pm/20261116_r1_bie_manual_package_01'}
$package=[IO.Path]::GetFullPath($PackageDir)
if(-not $package.StartsWith($project+[IO.Path]::DirectorySeparatorChar)){throw 'Package must stay inside project'}
if(Test-Path -LiteralPath $package){throw 'Package directory must not already exist'}
$parent=Join-Path $project 'results/comparisons/piston_2d_pm/20261116_r1_bie_cache_01'
$parentManifest=Get-Content -LiteralPath (Join-Path $parent 'FINAL_LOCAL_SHA256SUMS.txt')
function Assert-ParentSha([string]$Relative){
    $line=@($parentManifest | Where-Object {$_ -match ('  '+[regex]::Escape($Relative)+'$')})
    if($line.Count -ne 1){throw ('Missing validated parent binding: '+$Relative)}
    $sha=(Get-FileHash -LiteralPath (Join-Path $parent $Relative) -Algorithm SHA256).Hash.ToLowerInvariant()
    if($sha -ne $line[0].Substring(0,64)){throw ('Parent evidence SHA drift: '+$Relative)}
}
foreach($f in @('common.mat','checkpoint.json','jobs/bie_coarse/data.mat','jobs/bie_coarse/points.csv',
    'jobs/bie_coarse/receipt.json','C0_independent_audit/receipt.json','equivalence/C0_full.json','equivalence/metrics.csv')){Assert-ParentSha $f}
foreach($name in @('solve_helmholtz_bie_halfplane_vertical.m','solve_helmholtz_bie_halfplane_cached_vertical.m')){
    Assert-ParentSha ('code/'+$name)
    $current=Join-Path $project ('scripts/validation/support/'+$name)
    if((Get-FileHash -LiteralPath $current -Algorithm SHA256).Hash -ne (Get-FileHash -LiteralPath (Join-Path $parent ('code/'+$name)) -Algorithm SHA256).Hash){throw 'Validated cache code changed'}
}
New-Item -ItemType Directory -Path $package,(Join-Path $package code),(Join-Path $package parent),(Join-Path $package input/gaussian) | Out-Null
$code=@('scripts/validation/run_r1_bie_manual.m','scripts/validation/support/r1_bie_manual_io.m',
    'scripts/validation/support/audit_r1_bie_manual_stage.m','scripts/validation/support/audit_r1_bie_manual_final.m',
    'scripts/validation/support/solve_helmholtz_bie_halfplane_vertical.m',
    'scripts/validation/support/solve_helmholtz_bie_halfplane_cached_vertical.m',
    'scripts/validation/support/closeout_s3a_pm_surface.m','src/propagation/baffled_strip_piston_field_vertical.m',
    'scripts/validation/support/piston_independent_angular_spectrum_vertical.m',
    'scripts/comparisons/run_piston_pm_r1_job.m','scripts/validation/audit_piston_pm_r1.m',
    'scripts/reporting/render_piston_pm_r1.m','scripts/validation/support/launch_r1_bie_manual.sh',
    'scripts/validation/support/archive_r1_bie_manual.sh','scripts/validation/support/manage_r1_bie_manual.ps1')
foreach($f in $code){Copy-Item -LiteralPath (Join-Path $project $f) -Destination (Join-Path $package code)}
Copy-Item -LiteralPath (Join-Path $parent 'common.mat') -Destination (Join-Path $package common.mat)
$map=@{'jobs/bie_coarse/data.mat'='C0_data.mat';'jobs/bie_coarse/points.csv'='C0_points.csv';
    'jobs/bie_coarse/receipt.json'='C0_receipt.json';'C0_independent_audit/receipt.json'='C0_audit_receipt.json';
    'checkpoint.json'='checkpoint.json';'equivalence/C0_full.json'='C0_equivalence.json';'equivalence/metrics.csv'='small_equivalence.csv'}
foreach($f in $map.Keys){Copy-Item -LiteralPath (Join-Path $parent $f) -Destination (Join-Path $package ('parent/'+$map[$f]))}
foreach($name in @('frozen_source.mat','reference_data.mat','seed_260001_coefficients.csv')){
    Assert-ParentSha ('input/'+$name)
    Copy-Item -LiteralPath (Join-Path $parent ('input/'+$name)) -Destination (Join-Path $package input)
}
foreach($f in Get-ChildItem -LiteralPath (Join-Path $parent input/gaussian) -File){
    Assert-ParentSha ('input/gaussian/'+$f.Name)
    Copy-Item -LiteralPath $f.FullName -Destination (Join-Path $package input/gaussian)
}
foreach($f in Get-ChildItem -LiteralPath (Join-Path $package code) -Filter '*.sh' -File){
    [IO.File]::WriteAllText($f.FullName,[IO.File]::ReadAllText($f.FullName).Replace("`r`n","`n"),(New-Object Text.UTF8Encoding($false)))
}
$files=Get-ChildItem -LiteralPath (Join-Path $package code),(Join-Path $package parent),(Join-Path $package input) -File -Recurse | ForEach-Object {
    [pscustomobject]@{path=[IO.Path]::GetRelativePath($package,$_.FullName).Replace('\','/');sha256=(Get-FileHash -LiteralPath $_.FullName -Algorithm SHA256).Hash.ToLowerInvariant()}
}
$files+=([pscustomobject]@{path='common.mat';sha256=(Get-FileHash -LiteralPath (Join-Path $package common.mat) -Algorithm SHA256).Hash.ToLowerInvariant()})
$record=[ordered]@{schema_version=1;scope='R1_BIE_C1_C4_MANUAL';remote_task_name=$RemoteTaskName;
    parent_run_id='20261116_r1_bie_cache_01';parent_historical_process_s=525;
    execution_policy=[ordered]@{single_process=$true;singleCompThread=$true;parpool=$false;RAM_reserve_fraction=.30;
        cumulative_time_cap=$false;cost_prediction_gate=$false;automatic_stage_retry=$false;C0_recompute_allowed=$false;
        stage_sensitivity_gate=1e-4;stage_linear_residual_gate=1e-10;stage_independent_reconstruction_gate=1e-11;
        completion_action='STOP for R1 review; never R2'};
    recipes=@(@(60,12,96),@(60,12,64),@(60,14,96),@(70,12,96));files=@($files)}
[IO.File]::WriteAllText((Join-Path $package manifest.json),($record | ConvertTo-Json -Depth 12),(New-Object Text.UTF8Encoding($false)))
$sha=@($files | ForEach-Object {$_.sha256+'  '+$_.path})
$sha+=((Get-FileHash -LiteralPath (Join-Path $package manifest.json) -Algorithm SHA256).Hash.ToLowerInvariant()+'  manifest.json')
[IO.File]::WriteAllText((Join-Path $package package.sha256),([string]::Join("`n",$sha)+"`n"),(New-Object Text.UTF8Encoding($false)))
& tar -cf (Join-Path $package input_package.tar) -C $package code parent input common.mat manifest.json package.sha256
if($LASTEXITCODE -ne 0){throw 'Cannot create input archive'}
$archiveHash=(Get-FileHash -LiteralPath (Join-Path $package input_package.tar) -Algorithm SHA256).Hash.ToLowerInvariant()
[IO.File]::WriteAllText((Join-Path $package input_package.tar.sha256),$archiveHash+'  input_package.tar'+"`n",(New-Object Text.UTF8Encoding($false)))
Write-Output ('PACKAGE_PREPARED '+$package+'; no SSH or MATLAB started')
