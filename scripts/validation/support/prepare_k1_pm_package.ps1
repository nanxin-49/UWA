param([Parameter(Mandatory=$true)][string]$PackageDir,
 [ValidateSet('Prepare','Manifest','Archive')][string]$Action='Prepare')
$ErrorActionPreference='Stop'
$repo=[IO.Path]::GetFullPath((Join-Path $PSScriptRoot '../../..')); $pkg=[IO.Path]::GetFullPath($PackageDir)
$bound=Join-Path $repo results/comparisons/pm_wind_cutoff_k1
if(-not $pkg.StartsWith($bound+[IO.Path]::DirectorySeparatorChar,[StringComparison]::OrdinalIgnoreCase)){throw 'Use isolated K1 comparison directory'}
$utf8=[Text.UTF8Encoding]::new($false)
function Sha($p){(Get-FileHash -LiteralPath $p).Hash.ToLowerInvariant()}
function Json($p,$v){[IO.File]::WriteAllText($p,($v|ConvertTo-Json -Depth 50),$utf8)}
function Rows($root,$paths){@($paths|ForEach-Object {@{path=$_;sha256=(Sha (Join-Path $root $_))}})}
$r1=Join-Path $repo results/comparisons/piston_2d_pm/r1_bie_manual_c4_policy_v2_20261116_04_retrieved_20261116_173609
$r1work=Join-Path $r1 finalization/attempts/attempt_0001/work
$k0=Join-Path $repo results/experiments/pm_wind_cutoff/pm_wind_cutoff_k0_20261117_05_retrieved_20261117_232303
$cases=@(@{case_id='U6_N07';U_PM_mps=6;modes=7},@{case_id='U6_N40';U_PM_mps=6;modes=40},@{case_id='U5_N12';U_PM_mps=5;modes=12},@{case_id='U8_N12';U_PM_mps=8;modes=12},@{case_id='U5_N40';U_PM_mps=5;modes=40},@{case_id='U8_N40';U_PM_mps=8;modes=40})
if($Action -eq 'Prepare'){
 if(Test-Path -LiteralPath $pkg){throw 'New package directory required'}
 foreach($d in @('','code','cases','anchor','provenance','resource_evidence')){New-Item -ItemType Directory -Path (Join-Path $pkg $d)|Out-Null}
 $prior=Join-Path $repo results/comparisons/piston_2d_pm/20261116_r1_v2_review_01/review_receipt.json
 if((Sha $prior) -ne 'eb6d86a8e845ee2785efab4b93470f91d26026ed0d415bdbd808023152b43265'){throw 'R1 formal qualification SHA drift'}
 $r1lists=@(Get-ChildItem -LiteralPath (Join-Path $r1 retrieval) -Filter 'SHA256SUMS_*.txt' -File)
 if($r1lists.Count -ne 1 -or (Sha $r1lists[0].FullName) -ne 'fe3ffb19c6638d8cc6f378d4b5b72ef6c26cf275b95ea103e1d748835a8fd963'){throw 'R1 retrieval provenance drift'}
 $r1count=0; foreach($line in Get-Content -LiteralPath $r1lists[0].FullName){
  if($line -notmatch '^([0-9a-f]{64})  (.+)$'){throw 'R1 SHA row'}; $digest=$Matches[1]; $rel=$Matches[2] -replace '^\./',''
  if($rel -match '(^/|:|(^|/)\.\.(/|$)|^cash/)'){throw 'Unsafe R1 path'}
  if((Sha (Join-Path $r1 $rel)) -ne $digest){throw ('R1 original drift '+$rel)}; $r1count++
 }; if($r1count -ne 287){throw 'R1 retrieval count'}
 $k0list=Get-Content -LiteralPath (Join-Path $k0 retrieval/manifest.json) -Raw|ConvertFrom-Json
 foreach($row in $k0list){if((Sha (Join-Path $k0 $row.path)) -ne $row.sha256){throw ('K0 original drift '+$row.path)}}
 foreach($file in @('common.mat','comparison.mat','model_metrics.csv','receiver_fields_and_errors.csv','independent_final_metrics_audit.json')){
  Copy-Item -LiteralPath (Join-Path $r1work $file) -Destination (Join-Path $pkg ('anchor/'+$file))
 }
 New-Item -ItemType Directory -Path (Join-Path $pkg anchor/audit)|Out-Null
 Copy-Item -LiteralPath (Join-Path $r1work audit/receipt.json) -Destination (Join-Path $pkg anchor/audit/receipt.json)
 Copy-Item -LiteralPath $prior -Destination (Join-Path $pkg anchor/R1_formal_review.json)
 Copy-Item -LiteralPath (Join-Path $r1 finalization/committed.json) -Destination (Join-Path $pkg anchor/R1_final_seal.json)
 Copy-Item -LiteralPath $r1lists[0].FullName -Destination (Join-Path $pkg provenance/R1_retrieval_SHA256SUMS.txt)
 Copy-Item -LiteralPath (Join-Path $k0 output/seed260001_U6_N12/coefficients.csv) -Destination (Join-Path $pkg anchor/K0_U6_N12_coefficients.csv)
 Copy-Item -LiteralPath (Join-Path $k0 output/seed260001_U6_N12/metrics.json) -Destination (Join-Path $pkg anchor/K0_U6_N12_metrics.json)
 Copy-Item -LiteralPath (Join-Path $k0 output/seed260001_U6_N12/spectra.csv) -Destination (Join-Path $pkg anchor/K0_U6_N12_spectra.csv)
 Copy-Item -LiteralPath (Join-Path $r1 input/seed_260001_coefficients.csv) -Destination (Join-Path $pkg anchor/R1_coefficients.csv)
 Copy-Item -LiteralPath (Join-Path $r1 policy.json) -Destination (Join-Path $pkg anchor/policy.json)
 Copy-Item -LiteralPath (Join-Path $r1 jobs/bie_window/points.csv) -Destination (Join-Path $pkg anchor/R1_receiver_origin.csv)
 $points=Get-Content -LiteralPath (Join-Path $r1 jobs/bie_window/points.csv)
 if($points.Count -ne 412){throw 'R1 411 coordinates'}
 $rx=@('x_m,z_m')+@($points[1..411]|ForEach-Object {($_ -split ',')[0..1] -join ','})
 [IO.File]::WriteAllText((Join-Path $pkg anchor/receivers.csv),([string]::Join("`n",$rx)+"`n"),$utf8)
 $newCode=@('scripts/validation/run_k1_pm_case.m','scripts/comparisons/run_k1_piston_pm_job.m',
 'scripts/validation/support/k1_pm_io.m','scripts/validation/prepare_k1_pm_anchor.m',
 'scripts/validation/validate_k1_pm_preparation.m','scripts/reporting/render_k1_pm_stable_saved.m',
 'scripts/reporting/summarize_k1_pm_cases.m','scripts/reporting/render_k1_pm_summary_fields.m','scripts/validation/support/k1_pm_preflight.py',
 'scripts/validation/support/k1_pm_supervisor.py','scripts/validation/support/launch_k1_pm.sh',
 'scripts/validation/support/manage_k1_pm.ps1','scripts/experiments/pm_sensitivity_io.m')
 $frozenCode=@('scripts/validation/support/r1_bie_manual_io.m','scripts/validation/support/r1_bie_continuation_io.m',
 'scripts/validation/support/r1_bie_continuation_gates.m','scripts/validation/support/audit_r1_bie_manual_stage.m',
 'scripts/validation/support/solve_helmholtz_bie_halfplane_vertical.m','scripts/validation/support/solve_helmholtz_bie_halfplane_cached_vertical.m',
 'scripts/validation/support/closeout_s3a_pm_surface.m','scripts/validation/support/closeout_s3a_pm_fullfk.m',
 'scripts/validation/support/fk_single_receiver_incident.m','src/propagation/baffled_strip_piston_field_vertical.m',
 'scripts/validation/support/piston_independent_angular_spectrum_vertical.m','scripts/validation/audit_piston_pm_r1.m',
 'scripts/validation/support/audit_r2_pm_final.m','scripts/validation/support/r2_pm_io.m',
 'scripts/reporting/summarize_r2_pm_fields.m','scripts/validation/audit_r2_final_figures.m')
 $versions=@()
 foreach($rel in ($newCode+$frozenCode)){
  $src=Join-Path $repo $rel; $name=[IO.Path]::GetFileName($rel)
  $baseline=Join-Path $repo ('results/comparisons/piston_2d_pm/r2_preparation_20261117_03/seed260003/code/'+$name)
  if($rel -in $frozenCode -and (Test-Path -LiteralPath $baseline) -and (Sha $src) -ne (Sha $baseline)){throw ('Frozen algorithm/auditor drift '+$rel)}
  Copy-Item -LiteralPath $src -Destination (Join-Path $pkg ('code/'+$name))
  $versions+=@{path=$rel;sha256=(Sha $src);frozen_dependency=($rel -in $frozenCode)}
 }
 foreach($c in $cases){
  $case=Join-Path $pkg ('cases/'+$c.case_id); New-Item -ItemType Directory -Path $case,(Join-Path $case input),(Join-Path $case code)|Out-Null
  Copy-Item -LiteralPath (Join-Path $pkg code) -Destination $case -Recurse -Force
  Copy-Item -LiteralPath (Join-Path $pkg anchor/policy.json) -Destination (Join-Path $case policy.json)
  Copy-Item -LiteralPath (Join-Path $pkg anchor/receivers.csv) -Destination (Join-Path $case input/receivers.csv)
  $original=Join-Path $k0 ('output/seed260001_U'+$c.U_PM_mps+'_N'+('{0:00}' -f $c.modes))
  foreach($pair in @(@('coefficients.csv','k0_coefficients.csv'),@('spectra.csv','k0_spectra.csv'),@('surface.csv','k0_surface.csv'),@('metrics.json','k0_metrics.json'),@('checkpoint.json','k0_checkpoint.json'))){Copy-Item -LiteralPath (Join-Path $original $pair[0]) -Destination (Join-Path $case ('input/'+$pair[1]))}
  # Preserve decimal tokens, relabel physical depth coefficient columns explicitly.
  $lines=Get-Content -LiteralPath (Join-Path $original coefficients.csv); $rows=@('k_rad_per_m,cos_coeff_m,sin_coeff_m,Sk_m3')
  foreach($line in $lines[1..($lines.Count-1)]){$v=$line -split ','; $rows+=(@($v[1],$v[3],$v[4],$v[2]) -join ',')}
  [IO.File]::WriteAllText((Join-Path $case input/coefficients.csv),([string]::Join("`n",$rows)+"`n"),$utf8)
  $frozen=[ordered]@{seed=260001;case_id=$c.case_id;U_PM_mps=$c.U_PM_mps;wind_label='U_PM';modes=$c.modes;family='legacy12_gaussian_extension_v1';
   source=@{frequency_hz=4000;c_mps=1500;rho_kgpm3=1025;v0_mps=.001;total_width_m=.1;z_source_m=100};period_m=160;
   C2_plateau_halfwidth_m=42;C2_support_halfwidth_m=50;receivers=411;receiver_z_m=3;phase_contract_version='2.0.0';
   representation='absolute_physical peak Pa; exp(-i omega t)';phase_mask_fraction_of_BIE_RMS=.01;
   PS_recipes=@(@(128,.05,.1875),@(256,.05,.1875),@(512,.05,.1875),@(1024,.05,.1875),@(1024,.025,.1875),@(1024,.025,.09375));
   FullFK_recipes=@(@(400,50),@(800,50),@(960,60),@(1440,90));BIE_recipes=@(@(60,10,96),@(60,12,96),@(60,12,64),@(60,14,96),@(70,12,96))}
  Json (Join-Path $case input/frozen_config.json) $frozen
 }
 $resources=Join-Path $repo results/validation/k1_preparation
 foreach($r in @('resource_20261118_03/host_stdout.json','resource_20261118_04/matlab_diagnostic.json','resource_20261118_05/host_stdout.json')){
  $safe=($r -replace '/','_'); Copy-Item -LiteralPath (Join-Path $resources $r) -Destination (Join-Path $pkg ('resource_evidence/'+$safe))
 }
 $taskName='k1_pm_'+((Split-Path $pkg -Leaf) -replace '^preparation_','')
 if($taskName -notmatch '^k1_pm_[0-9]{8}_[0-9]{2}$'){throw 'Package leaf must be preparation_YYYYMMDD_NN'}
 Json (Join-Path $pkg config.json) ([ordered]@{schema_version=1;kind='K1_PM_SIX_NEW_SURFACES';remote_task_name=$taskName;cases=$cases;
  execution_mode='A_SERIAL_CASES_SINGLE_COMPUTATIONAL_THREAD';max_processes=1;threads=1;parpool=$false;automatic_retry=$false;
  RAM_reserve_fraction=.3;case_working_bytes=12884901888;task_RSS_cap_bytes=12884901888;minimum_disk_bytes=5368709120;
  minimum_free_inodes=5000;max_case_wall_s=1800;max_total_wall_s=10800;max_summary_wall_s=600;
  admission_memory_PSI_full_avg10_max=5;admission_io_PSI_full_avg10_max=10;runtime_memory_PSI_full_avg10_max=10;
  runtime_io_PSI_full_avg10_max=30;pressure_stop_consecutive_samples=10;load1_max=21;user_manual_start_only=$true})
 Json (Join-Path $pkg provenance/sources.json) @{code=$versions;R1_original_files_checked=$r1count;K0_original_files_checked=$k0list.Count;
  R1_prior_review_sha256=(Sha $prior);R1_retrieval_manifest_sha256=(Sha $r1lists[0].FullName);K0_manifest_sha256=(Sha (Join-Path $k0 manifest.json));
  original_solver_or_numerical_gate_modified=$false;new_acoustic_solves_started=$false;parallel_license_limit='unknown; no additional process seat inferred'}
 Copy-Item -LiteralPath (Join-Path $repo scripts/comparisons/run_r2_piston_pm_job.m) -Destination (Join-Path $pkg provenance/r2_reference_model_job.m)
 Copy-Item -LiteralPath (Join-Path $repo results/visualization/pm_wind_cutoff_review/20261118_01/k1_selected_geometry.csv) -Destination (Join-Path $pkg provenance/k1_selected_geometry.csv)
 Copy-Item -LiteralPath (Join-Path $repo results/validation/k1_preparation/historical_20261118_02/receipt.json) -Destination (Join-Path $pkg resource_evidence/historical_BIE_receipt.json)
 Copy-Item -LiteralPath (Join-Path $repo results/validation/k1_preparation/historical_20261118_02/BIE_stage_costs.csv) -Destination (Join-Path $pkg resource_evidence/historical_BIE_stages.csv)
 Copy-Item -LiteralPath (Join-Path $repo results/validation/pm_wind_cutoff_retrieval/numeric_20261118_01/receipt.json) -Destination (Join-Path $pkg provenance/K0_saved_review_receipt.json)
 Write-Output 'K1_PREPARED_UNSEALED; next local strict anchor verification; no SSH/start'
} elseif($Action -eq 'Manifest'){
 $proof=Join-Path $pkg anchor/anchor_proof.json
 $a=Get-Content -LiteralPath $proof -Raw|ConvertFrom-Json
 if($a.status -ne 'STRICT_R1_K0_ANCHOR_MATCH_PASS' -or $a.failed_checks -ne 0){throw 'Strict anchor proof required'}
 foreach($c in $cases){
  $case=Join-Path $pkg ('cases/'+$c.case_id); Copy-Item -LiteralPath $proof -Destination (Join-Path $case input/r1_anchor_proof.json)
  $files=@(Get-ChildItem -LiteralPath (Join-Path $case code),(Join-Path $case input) -Recurse -File|ForEach-Object {@{path=$_.FullName.Substring($case.Length+1).Replace('\','/');sha256=(Sha $_.FullName)}})
  $files+=@{path='policy.json';sha256=(Sha (Join-Path $case policy.json))}
  Json (Join-Path $case manifest.json) @{schema_version=1;scope='K1_PM_CASE';seed=260001;case_id=$c.case_id;remote_task_name=$c.case_id;
   U_PM_mps=$c.U_PM_mps;N=$c.modes;max_case_processes=1;threads=1;parpool=$false;automatic_retry=$false;RAM_reserve_fraction=.3;
   policy_sha256=(Sha (Join-Path $case policy.json));anchor_proof_sha256=(Sha $proof);files=$files}
 }
 $files=@(Get-ChildItem -LiteralPath $pkg -Recurse -File|Where-Object {$_.Name -notin @('manifest.json','input_package.tar','input_package.tar.sha256')}|ForEach-Object {@{path=$_.FullName.Substring($pkg.Length+1).Replace('\','/');sha256=(Sha $_.FullName)}})
 foreach($c in $cases){$rel='cases/'+$c.case_id+'/manifest.json'; $files+=@{path=$rel;sha256=(Sha (Join-Path $pkg $rel))}}
 Json (Join-Path $pkg manifest.json) @{schema_version=1;scope='K1_PM_MASTER';files=$files;anchor_proof_sha256=(Sha $proof)}
 Write-Output ('K1_MANIFEST_READY files='+$files.Count)
} else {
 & python -B (Join-Path $pkg code/k1_pm_preflight.py) $pkg
 if($LASTEXITCODE -ne 0){throw 'K1 preflight failed'}
 if(Test-Path -LiteralPath (Join-Path $pkg input_package.tar)){throw 'Already archived'}
 & tar -cf (Join-Path $pkg input_package.tar) -C $pkg code cases anchor provenance resource_evidence config.json manifest.json
 if($LASTEXITCODE -ne 0){throw 'Archive failed'}
 $hash=Sha (Join-Path $pkg input_package.tar)
 [IO.File]::WriteAllText((Join-Path $pkg input_package.tar.sha256),"$hash  input_package.tar`n",$utf8)
 Write-Output ('K1_ARCHIVE_SHA256='+$hash)
}
