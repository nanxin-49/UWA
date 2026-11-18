param([Parameter(Mandatory=$true)][string]$OutputDir)
$ErrorActionPreference='Stop'
$repo=[IO.Path]::GetFullPath((Join-Path $PSScriptRoot '../../..'))
$dest=[IO.Path]::GetFullPath($OutputDir)
$allowed=Join-Path $repo 'results/experiments/pm_wind_cutoff'
if(-not $dest.StartsWith($allowed+[IO.Path]::DirectorySeparatorChar,[StringComparison]::OrdinalIgnoreCase)){throw 'Use case-specific results/experiments/pm_wind_cutoff directory'}
if(Test-Path -LiteralPath $dest){throw 'New package directory required; never overwrite'}
New-Item -ItemType Directory -Path $dest,(Join-Path $dest 'code'),(Join-Path $dest 'inputs'),(Join-Path $dest 'provenance')|Out-Null
$utf8=[Text.UTF8Encoding]::new($false)
function Write-Json($path,$value){[IO.File]::WriteAllText($path,($value|ConvertTo-Json -Depth 35),$utf8)}
$sourcePaths=@(
 'scripts/experiments/pm_sensitivity_io.m',
 'scripts/experiments/build_pm_wind_cutoff_surface.m',
 'scripts/experiments/run_pm_wind_cutoff_k0.m',
 'scripts/experiments/run_pm_wind_cutoff_k1.m',
 'scripts/validation/audit_pm_wind_cutoff_k0_case.m',
 'scripts/reporting/render_pm_wind_cutoff_k0.m',
 'scripts/validation/support/pm_k0_preflight.py',
 'scripts/validation/support/pm_k0_supervisor.py',
 'scripts/validation/support/launch_pm_wind_cutoff_k0.sh',
 'scripts/validation/support/manage_pm_wind_cutoff_k0.ps1')
$sources=@()
foreach($rel in $sourcePaths){
 $src=Join-Path $repo $rel
 Copy-Item -LiteralPath $src -Destination (Join-Path $dest ('code/'+[IO.Path]::GetFileName($rel)))
 $sources+=@{path=$rel;sha256=(Get-FileHash -LiteralPath $src -Algorithm SHA256).Hash.ToLowerInvariant()}
}
$anchors=@(
 'results/validation/bellhop_internal_pm_fixed_realization/fixed_pm_fourier_coefficients.csv',
 'results/validation/spectral_fast_kirchhoff_pm_pilot/seed_260002_coefficients.csv',
 'results/validation/spectral_fast_kirchhoff_pm_pilot/seed_260003_coefficients.csv')
for($j=0;$j -lt 3;$j++){
 $src=Join-Path $repo $anchors[$j]; $seed=260001+$j
 Copy-Item -LiteralPath $src -Destination (Join-Path $dest "inputs/seed$seed.csv")
 $sources+=@{path=$anchors[$j];sha256=(Get-FileHash -LiteralPath $src -Algorithm SHA256).Hash.ToLowerInvariant()}
 if($j -gt 0){
  $parent=Join-Path $repo "results/validation/pe_bellhop_pm_amplitude_ensemble/seed_${seed}_coefficients.csv"
  if((Get-FileHash -LiteralPath $parent).Hash -ne (Get-FileHash -LiteralPath $src).Hash){throw 'Anchor parent SHA differs'}
  $sources+=@{path=$parent.Substring($repo.Length+1).Replace('\','/');sha256=(Get-FileHash -LiteralPath $parent).Hash.ToLowerInvariant()}
 }
}
$reviewed=@('src/surface/raw_pm_spectrum_grid_vertical.m','src/surface/pm_surface_boundary_model.m',
 'scripts/validation/support/build_paired_pm_gaussian_realization.m',
 'scripts/validation/support/load_fixed_pm_profile_for_pe_bellhop_validation.m',
 'scripts/validation/support/closeout_s3a_pm_surface.m','reports/spectral_fast_kirchhoff_pm_pilot_report.md',
 'reports/r1_bie_v2_final_review_2026-11-16.md','reports/r2_saved_closeout_final_2026-11-17.md')
foreach($rel in $reviewed){$src=Join-Path $repo $rel; $sources+=@{path=$rel;sha256=(Get-FileHash -LiteralPath $src).Hash.ToLowerInvariant()}}
$cases=@(); foreach($seed in 260001..260003){foreach($u in @(5,6,8)){foreach($n in @(7,12,19,32,40)){$cases+=@{seed=$seed;U_PM_mps=$u;N=$n}}}}
$cfg=[ordered]@{schema_version=1;kind='K0_SURFACE_ONLY';family='legacy12_gaussian_extension_v1';wind_label='U_PM';wind_height_m='unknown';wind_averaging_duration_s='unknown';U_PM_mps=@(5,6,8);N=@(7,12,19,32,40);seeds=@(260001,260002,260003);L_m=160;sample_count=8192;coarse_sample_count=4096;C2_plateau_m=42;C2_zero_m=50;normalization='none';cases=$cases;max_wall_s=1800;RAM_reserve_fraction=.30;estimated_incremental_surface_arrays_GiB=.05;minimum_available_RAM_GiB=2;workers=1;threads=1;K1_enabled=$false}
Write-Json (Join-Path $dest config.json) $cfg
$head=& git -C $repo rev-parse HEAD; if($LASTEXITCODE -ne 0){throw 'Git provenance unavailable'}
Write-Json (Join-Path $dest 'provenance/sources.json') @{git_head=$head;sources=$sources;scope='K0 package preparation only; no SSH or acoustic solves';original_PM_files_unchanged=$true}
$files=@(Get-ChildItem -LiteralPath $dest -Recurse -File|Sort-Object FullName|ForEach-Object {@{path=$_.FullName.Substring($dest.Length+1).Replace('\','/');sha256=(Get-FileHash -LiteralPath $_.FullName -Algorithm SHA256).Hash.ToLowerInvariant()}})
Write-Json (Join-Path $dest manifest.json) @{schema_version=1;files=$files}
& python -B (Join-Path $dest code/pm_k0_preflight.py) $dest
if($LASTEXITCODE -ne 0){throw 'Package preflight failed'}
& tar -cf (Join-Path $dest input_package.tar) -C $dest code inputs provenance config.json manifest.json
if($LASTEXITCODE -ne 0){throw 'Archive creation failed'}
$sha=(Get-FileHash -LiteralPath (Join-Path $dest input_package.tar)).Hash.ToLowerInvariant()
[IO.File]::WriteAllText((Join-Path $dest input_package.tar.sha256),"$sha  input_package.tar`n",$utf8)
Write-Output "PACKAGE=$dest"
Write-Output "ARCHIVE_SHA256=$sha"
