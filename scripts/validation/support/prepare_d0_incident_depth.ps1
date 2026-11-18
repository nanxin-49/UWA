param([Parameter(Mandatory=$true)][string]$OutputDir)
$ErrorActionPreference='Stop'
$repo=[IO.Path]::GetFullPath((Join-Path $PSScriptRoot '../../..'))
$dest=[IO.Path]::GetFullPath($OutputDir)
$allowed=Join-Path $repo 'results/experiments/d0_incident_depth'
if(-not $dest.StartsWith($allowed+[IO.Path]::DirectorySeparatorChar,[StringComparison]::OrdinalIgnoreCase)){throw 'Use a new results/experiments/d0_incident_depth directory'}
if(Test-Path -LiteralPath $dest){throw 'Never overwrite a preparation directory'}
$utf8=[Text.UTF8Encoding]::new($false)
function Write-Json($path,$value){[IO.File]::WriteAllText($path,($value|ConvertTo-Json -Depth 50),$utf8)}
function Hash($path){(Get-FileHash -LiteralPath $path -Algorithm SHA256).Hash.ToLowerInvariant()}
New-Item -ItemType Directory -Path $dest,(Join-Path $dest code),(Join-Path $dest inputs/k0),(Join-Path $dest provenance)|Out-Null
$paths=@('src/propagation/baffled_strip_piston_field_vertical.m',
 'src/channel/apply_pe_channel_phase_reference_vertical.m',
 'scripts/validation/support/piston_independent_angular_spectrum_vertical.m',
 'scripts/experiments/pm_sensitivity_io.m',
 'scripts/experiments/d0_incident_diagnostics.m','scripts/experiments/run_d0_incident_depth.m',
 'scripts/validation/support/d0_incident_spectrum_reference.m','scripts/validation/audit_d0_incident_case.m',
 'scripts/reporting/render_d0_incident_depth.m','scripts/validation/support/d0_package.py',
 'scripts/validation/support/d0_supervisor.py','scripts/validation/support/launch_d0_incident_depth.sh',
 'scripts/validation/support/manage_d0_incident_depth.ps1','scripts/validation/support/invoke_r1_bie_status_stdin.ps1')
$sources=@()
foreach($rel in $paths){
 $src=Join-Path $repo $rel
 Copy-Item -LiteralPath $src -Destination (Join-Path $dest ('code/'+[IO.Path]::GetFileName($rel)))
 $sources+=@{path=$rel;sha256=(Hash $src)}
}
$k0='results/experiments/pm_wind_cutoff/pm_wind_cutoff_k0_20261117_05_retrieved_20261117_232303'
$case=$k0+'/output/seed260001_U6_N12'
$cp=Get-Content -LiteralPath (Join-Path $repo ($case+'/checkpoint.json')) -Raw|ConvertFrom-Json
if($cp.status -ne 'AUDITED_PASS' -or $cp.case.seed -ne 260001 -or $cp.case.N -ne 12 -or $cp.case.U_PM_mps -ne 6){throw 'K0 selection mismatch'}
foreach($row in $cp.artifacts_sha256){
 if($row.path -notmatch '^[A-Za-z0-9_.-]+$'){throw 'Unexpected K0 artifact path'}
 $rel=$case+'/'+$row.path; $src=Join-Path $repo $rel
 if((Hash $src) -ne $row.sha256){throw ('K0 SHA drift '+$rel)}
 Copy-Item -LiteralPath $src -Destination (Join-Path $dest ('inputs/k0/'+$row.path))
 $sources+=@{path=$rel;sha256=$row.sha256}
}
Copy-Item -LiteralPath (Join-Path $repo ($case+'/checkpoint.json')) -Destination (Join-Path $dest inputs/k0/checkpoint.json)
$sources+=@{path=$case+'/checkpoint.json';sha256=(Hash (Join-Path $repo ($case+'/checkpoint.json')))}
$reviewed=@('AGENTS.md','reports/pm_wind_cutoff_k0_review_2026-11-18.md',
 'reports/r1_bie_v2_final_review_2026-11-16.md','reports/r2_saved_closeout_final_2026-11-17.md',
 'results/comparisons/piston_2d_single_frequency/20261112_remote_p0p1_01/report.md',
 'scripts/validation/run_piston_p0p1_remote_job.m','scripts/validation/audit_piston_p0p1_remote.m',
 'reports/codex_token_execution_protocol.md')
foreach($rel in $reviewed){$sources+=@{path=$rel;sha256=(Hash (Join-Path $repo $rel))}}
$cfg=[ordered]@{
 schema_version=1;kind='D0_INCIDENT_ONLY';depths_m=@(100,300,600,1000,1490);fine_depths_m=@(100,1490)
 source=@{frequency_hz=4000;c_mps=1500;rho_kgpm3=1025;v0_mps=.001;total_width_m=.1}
 water_depth_context_m=1500;x_source_m=0;time_convention='exp(-i omega t)';phase_contract_version='2.0.0';representation='absolute_physical'
 source_model='2D infinite uniform in-phase strip, ideal rigid rear baffle, upward halfspace only'
 medium_model='uniform lossless water; no SSP, absorption or current'
 x_halfwidth_m=220;dx_m=.05;fine_dx_m=.025;aperture_quadrature_order=64;audit_aperture_order=128
 audit_x_m=@(-220,-200,-150,-100,-50,-20,0,20,50,100,150,200,220)
 PM=@{U_PM_mps=6;N=12;seed=260001;L_m=160;C2_plateau_m=42;C2_zero_m=50;geometry='original saved embedded samples abs(x)<=50; z_actual=-eta; no extension'}
 gates=@{source_relative=1e-7;gradient_relative=1e-6;phase_absolute_Pa=1e-10;csv_relative=1e-12;fft_relative=1e-12;dft_relative=1e-10;common_point_relative=1e-12;sampling_relative=1e-4;geometry_absolute_m=1e-11;max_phase_step_rad=([Math]::PI/4);phase_mask_fraction=1e-4;spectrum_phase_mask_fraction=1e-6}
 workers=1;threads=1;RAM_reserve_fraction=.30;minimum_available_RAM_GiB=4;estimated_working_RAM_GiB=2
 max_wall_s=3600;max_depth_wall_s=900;reflection_enabled=$false;bottom_enabled=$false
 retry_policy='none; every attempted package is one-shot; manual inspection after failure'
 excluded=@('bottom','downward radiation','surface reflection','PS','BIE','FullFK','wideband','CIR','communications','Bellhop','SSP','absorption','currents','shell scattering')
}
Write-Json (Join-Path $dest config.json) $cfg
$head=& git -C $repo rev-parse HEAD
if($LASTEXITCODE -ne 0){throw 'Git provenance failed'}
$dirty=@(& git -C $repo status --short --untracked-files=no)
Write-Json (Join-Path $dest provenance/sources.json) @{git_head=$head;tracked_worktree_status=$dirty;sources=$sources;scope='Preparation only, no SSH and no formal field evaluation';selected_K0_case=$case}
$files=@(Get-ChildItem -LiteralPath $dest -Recurse -File|Sort-Object FullName|ForEach-Object {@{path=$_.FullName.Substring($dest.Length+1).Replace('\','/');sha256=(Hash $_.FullName)}})
Write-Json (Join-Path $dest manifest.json) @{schema_version=1;files=$files}
& python -B (Join-Path $dest code/d0_package.py) verify $dest
if($LASTEXITCODE -ne 0){throw 'D0 preflight failed; retain candidate'}
& tar -cf (Join-Path $dest input_package.tar) -C $dest code inputs provenance config.json manifest.json
if($LASTEXITCODE -ne 0){throw 'D0 archive creation failed'}
$hash=Hash (Join-Path $dest input_package.tar)
[IO.File]::WriteAllText((Join-Path $dest input_package.tar.sha256),"$hash  input_package.tar`n",$utf8)
Write-Output "PACKAGE=$dest"
Write-Output "ARCHIVE_SHA256=$hash"
