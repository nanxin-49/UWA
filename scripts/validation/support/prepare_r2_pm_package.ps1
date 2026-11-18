param(
    [Parameter(Mandatory=$true)][string]$ProjectRoot,
    [Parameter(Mandatory=$true)][string]$PackageDir,
    [Parameter(Mandatory=$true)][ValidateSet(260002,260003)][int]$Seed,
    [ValidateSet('Prepare','Seal')][string]$Action='Prepare'
)
$ErrorActionPreference='Stop'
$project=(Resolve-Path -LiteralPath $ProjectRoot).Path
$package=[IO.Path]::GetFullPath($PackageDir)
if(-not $package.StartsWith($project+[IO.Path]::DirectorySeparatorChar)){throw 'Package must stay in project'}
function Sha([string]$p){(Get-FileHash -LiteralPath $p -Algorithm SHA256).Hash.ToLowerInvariant()}
function WriteJson([string]$p,$v){[IO.File]::WriteAllText($p,($v|ConvertTo-Json -Depth 35),[Text.UTF8Encoding]::new($false))}
if($Action -eq 'Prepare'){
    if(Test-Path -LiteralPath $package){throw 'New package directory required'}
    $r1=Join-Path $project 'results/comparisons/piston_2d_pm/r1_bie_manual_c4_policy_v2_20261116_04_retrieved_20261116_173609'
    $r1manifest=Get-Content -Raw -LiteralPath (Join-Path $r1 manifest.json)|ConvertFrom-Json
    # Only exact active inputs/code are read. Never enumerate or restore cash/.
    foreach($f in $r1manifest.files){
        if($f.path -in @('policy.json','input/seed_260001_coefficients.csv') -or $f.path -like 'code/*'){
            if((Sha (Join-Path $r1 $f.path)) -ne $f.sha256){throw ('R1 source drift '+$f.path)}
        }
    }
    New-Item -ItemType Directory -Path $package,(Join-Path $package code),(Join-Path $package input),(Join-Path $package provenance)|Out-Null
    $code=@('scripts/validation/run_r2_pm_seed.m','scripts/comparisons/run_r2_piston_pm_job.m',
        'scripts/validation/validate_r2_pm_preparation.m','scripts/validation/support/r2_pm_io.m',
        'scripts/validation/support/audit_r2_pm_final.m','scripts/reporting/render_r2_pm_seed.m',
        'scripts/reporting/summarize_r2_pm_fields.m','scripts/validation/support/r2_pm_preflight.py',
        'scripts/validation/support/launch_r2_pm_seed.sh','scripts/validation/support/manage_r2_pm_seed.ps1',
        'scripts/validation/support/archive_r2_pm_seed.sh','scripts/validation/support/status_r2_pm_session.sh',
        'scripts/validation/support/r2_pm_parallel.sh',
        'scripts/validation/support/prepare_r2_pm_package.ps1','scripts/validation/support/test_r2_pm_package.py',
        'scripts/validation/support/test_r2_pm_parallel.sh','scripts/validation/support/test_r2_pm_status.sh',
        'scripts/validation/support/invoke_r1_bie_status_stdin.ps1',
        'scripts/validation/support/r1_bie_runtime.sh','scripts/validation/support/r1_bie_runtime_json.py',
        'scripts/validation/support/r1_bie_manual_io.m','scripts/validation/support/r1_bie_continuation_io.m',
        'scripts/validation/support/r1_bie_continuation_gates.m','scripts/validation/support/audit_r1_bie_manual_stage.m',
        'scripts/validation/support/solve_helmholtz_bie_halfplane_vertical.m',
        'scripts/validation/support/solve_helmholtz_bie_halfplane_cached_vertical.m',
        'scripts/validation/support/closeout_s3a_pm_surface.m','scripts/validation/support/closeout_s3a_pm_fullfk.m',
        'scripts/validation/support/fk_single_receiver_incident.m',
        'src/propagation/baffled_strip_piston_field_vertical.m',
        'scripts/validation/support/piston_independent_angular_spectrum_vertical.m','scripts/validation/audit_piston_pm_r1.m')
    $versions=@()
    foreach($rel in $code){
        $current=Join-Path $project $rel
        $old=Join-Path $r1 ('code/'+(Split-Path $rel -Leaf))
        if(-not(Test-Path -LiteralPath $old)){
            $old=Join-Path $project ('results/comparisons/piston_2d_pm/20261115_r1_seed260001_01/code/'+(Split-Path $rel -Leaf))
        }
        if(Test-Path -LiteralPath $old){if((Sha $current) -ne (Sha $old)){throw ('R1 validated dependency differs '+$rel)}}
        Copy-Item -LiteralPath $current -Destination (Join-Path $package code)
        $versions+=([ordered]@{path=$rel;sha256=(Sha $current);R1_byte_identical=(Test-Path -LiteralPath $old)})
    }
    foreach($name in @('run_piston_pm_r1_job.m','render_piston_pm_r1.m','audit_r1_bie_manual_final.m','launch_r1_bie_manual.sh','run_r1_bie_manual.m')){
        $origin=Join-Path $r1 ('code/'+$name)
        if(-not(Test-Path -LiteralPath $origin)){
            $rel='scripts/validation/support/'+$name
            if($name -eq 'run_r1_bie_manual.m'){$rel='scripts/validation/'+$name}
            $origin=Join-Path $project $rel
        }
        Copy-Item -LiteralPath $origin -Destination (Join-Path $package provenance)
    }
    $coeff=Join-Path $project ('results/validation/spectral_fast_kirchhoff_pm_pilot/seed_'+$Seed+'_coefficients.csv')
    $parent=Join-Path $project ('results/validation/pe_bellhop_pm_amplitude_ensemble/seed_'+$Seed+'_coefficients.csv')
    if((Sha $coeff) -ne (Sha $parent)){throw 'Saved seed differs from original PM Gaussian realization'}
    Copy-Item -LiteralPath $coeff -Destination (Join-Path $package input/coefficients.csv)
    Copy-Item -LiteralPath (Join-Path $r1 input/seed_260001_coefficients.csv) -Destination (Join-Path $package input/r1_coefficients.csv)
    # Coordinates are extracted as original decimal tokens, without interpolation.
    $points=Join-Path $r1 'jobs/bie_window/points.csv'
    $lines=Get-Content -LiteralPath $points
    if($lines.Count -ne 412 -or $lines[0] -ne 'x_m,z_m,real_Pa,imag_Pa'){throw 'Unexpected R1 coordinate CSV'}
    $rx=@('x_m,z_m')+@($lines[1..411]|ForEach-Object {($_ -split ',')[0..1] -join ','})
    [IO.File]::WriteAllText((Join-Path $package input/receivers.csv),([string]::Join("`n",$rx)+"`n"),[Text.UTF8Encoding]::new($false))
    $frozen=[ordered]@{seed=$Seed;source=[ordered]@{frequency_hz=4000;c_mps=1500;rho_kgpm3=1025;v0_mps=.001;total_width_m=.1;z_source_m=100};
        period_m=160;modes=12;C2_plateau_halfwidth_m=42;C2_support_halfwidth_m=50;receivers=411;receiver_z_m=3;
        phase_contract_version='2.0.0';representation='absolute_physical peak Pa; exp(-i omega t)';phase_mask_fraction_of_BIE_RMS=.01;
        BIE_recipes=@(@(60,10,96),@(60,12,96),@(60,12,64),@(60,14,96),@(70,12,96));
        PS_recipes=@(@(128,.05,.1875),@(256,.05,.1875),@(512,.05,.1875),@(1024,.05,.1875),@(1024,.025,.1875),@(1024,.025,.09375));
        FullFK_recipes=@(@(400,50),@(800,50),@(960,60),@(1440,90));
        ensemble='existing independent Gaussian cosine/sine coefficients from unchanged PM Sk; no new draw, normalization or phase-only replacement'}
    WriteJson (Join-Path $package input/frozen_config.json) $frozen
    Copy-Item -LiteralPath (Join-Path $r1 policy.json) -Destination (Join-Path $package policy.json)
    $policy=Get-Content -Raw -LiteralPath (Join-Path $package policy.json)|ConvertFrom-Json
    if($policy.version -ne 'r1-bie-reference-v2' -or $policy.grid_monotonicity_mode -ne 'diagnostic'){throw 'Unexpected policy'}
    WriteJson (Join-Path $package provenance/sources.json) ([ordered]@{seed=$Seed;code=$versions;R1_run=(Split-Path $r1 -Leaf);
        R1_manifest_sha256=(Sha (Join-Path $r1 manifest.json));policy_sha256=(Sha (Join-Path $r1 policy.json));
        coefficient_path=[IO.Path]::GetRelativePath($project,$coeff).Replace('\','/');coefficient_sha256=(Sha $coeff);
        coefficient_parent_path=[IO.Path]::GetRelativePath($project,$parent).Replace('\','/');coefficient_parent_sha256=(Sha $parent);
        receiver_origin_path='jobs/bie_window/points.csv';receiver_origin_sha256=(Sha $points);
        task_scope='2026-11-17 user authorizes preparation/static/lightweight/dry-run only; agent must not SSH or start an acoustic solve';
        density_reuse=$false;R1_results_modified=$false;git_commit_push=$false})
    Copy-Item -LiteralPath (Join-Path $project reports/r2_pm_multiseed_plan_2026-11-17.md) -Destination (Join-Path $package plan.md)
    $files=@(Get-ChildItem -LiteralPath (Join-Path $package code),(Join-Path $package input),(Join-Path $package provenance) -Recurse -File|ForEach-Object {
        [ordered]@{path=[IO.Path]::GetRelativePath($package,$_.FullName).Replace('\','/');sha256=(Sha $_.FullName)}})
    foreach($name in @('policy.json','plan.md')){$files+=([ordered]@{path=$name;sha256=(Sha (Join-Path $package $name))})}
    WriteJson (Join-Path $package manifest.json) ([ordered]@{schema_version=1;scope='R2_PM_TWO_SEED_REVIEW';seed=$Seed;
        remote_task_name=('r2_pm_seed'+$Seed+'_20261117_03');single_process=$true;single_thread=$true;parpool=$false;
        RAM_reserve_fraction=.3;automatic_retry=$false;cumulative_time_cap=$false;cost_prediction_gate=$false;
        max_parallel_seed_processes=2;parallel_scope='only seeds 260002 and 260003; one singleCompThread MATLAB process each';
        estimated_working_bytes_per_seed=12884901888;aggregate_admission_bytes=25769803776;
        policy_sha256=(Sha (Join-Path $package policy.json));files=$files})
    Write-Output ('INPUTS_PREPARED '+$package+'; no SSH, MATLAB or solve started')
} else {
    $m=Get-Content -Raw -LiteralPath (Join-Path $package manifest.json)|ConvertFrom-Json
    if($m.seed -ne $Seed){throw 'Wrong seed'}
    & python (Join-Path $package code/r2_pm_preflight.py) $package
    if($LASTEXITCODE -ne 0){throw 'Local package preflight failed'}
    $proof=Get-Content -Raw -LiteralPath (Join-Path $package preparation_validation/static_receipt.json)|ConvertFrom-Json
    if($proof.status -ne 'PASS' -or $proof.failed_checks -ne 0 -or $proof.new_acoustic_solves -ne 0){throw 'Missing local static acceptance'}
    if(Test-Path -LiteralPath (Join-Path $package input_package.tar)){throw 'Already sealed; use new package for changes'}
    $items=@(Get-ChildItem -LiteralPath $package -Recurse -File)
    $rows=@($items|Sort-Object FullName|ForEach-Object {(Sha $_.FullName)+'  '+[IO.Path]::GetRelativePath($package,$_.FullName).Replace('\','/')})
    [IO.File]::WriteAllText((Join-Path $package package.sha256),([string]::Join("`n",$rows)+"`n"),[Text.UTF8Encoding]::new($false))
    $names=@(Get-ChildItem -LiteralPath $package|Select-Object -ExpandProperty Name)
    & tar -cf (Join-Path $package input_package.tar) -C $package @names
    if($LASTEXITCODE -ne 0){throw 'Archive creation failed'}
    [IO.File]::WriteAllText((Join-Path $package input_package.tar.sha256),(Sha (Join-Path $package input_package.tar))+'  input_package.tar'+"`n",[Text.UTF8Encoding]::new($false))
    Write-Output ('SEALED_STATIC_ACCEPTANCE_ONLY '+$package+'; MATLAB dry-run/Analyzer unavailable if recorded; user starts manually')
}
