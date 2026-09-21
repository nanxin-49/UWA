param(
    [string]$OfficialToolboxRoot = '',
    [string]$CompilerRoot = '',
    [string]$OutputDir = ''
)

$ErrorActionPreference = 'Stop'
if ([string]::IsNullOrWhiteSpace($OfficialToolboxRoot)) {
    $OfficialToolboxRoot = $env:BELLHOP_TOOLBOX_ROOT
}
if ([string]::IsNullOrWhiteSpace($OfficialToolboxRoot)) {
    throw 'Set -OfficialToolboxRoot or BELLHOP_TOOLBOX_ROOT to the AcousticsToolbox 2020 source root.'
}
$supportDir = [System.IO.Path]::GetFullPath($PSScriptRoot)
$projectRoot = [System.IO.Path]::GetFullPath((Join-Path $supportDir '..\..\..\..'))
if ([string]::IsNullOrWhiteSpace($OutputDir)) {
    $OutputDir = Join-Path $projectRoot 'results\validation\bellhop_internal_flat_wall_poc'
}
$OutputDir = [System.IO.Path]::GetFullPath($OutputDir)
$expectedRoot = [System.IO.Path]::GetFullPath((Join-Path $projectRoot 'results\validation'))
if (-not $OutputDir.StartsWith($expectedRoot, [System.StringComparison]::OrdinalIgnoreCase)) {
    throw "OutputDir must remain under $expectedRoot"
}

$officialRoot = [System.IO.Path]::GetFullPath($OfficialToolboxRoot)
foreach ($required in @('Bellhop\bellhop.f90', 'Bellhop\Step.f90', 'misc\Makefile')) {
    if (-not (Test-Path -LiteralPath (Join-Path $officialRoot $required) -PathType Leaf)) {
        throw "Missing official Bellhop 2020 source: $required"
    }
}

$buildRoot = Join-Path $OutputDir 'build\AcousticsToolbox_2020'
$binDir = Join-Path $OutputDir 'bin'
if (Test-Path -LiteralPath $buildRoot) {
    $resolvedBuild = [System.IO.Path]::GetFullPath($buildRoot)
    if (-not $resolvedBuild.StartsWith($expectedRoot, [System.StringComparison]::OrdinalIgnoreCase)) {
        throw "Refusing to remove unexpected build path: $resolvedBuild"
    }
    Remove-Item -LiteralPath $resolvedBuild -Recurse -Force
}
New-Item -ItemType Directory -Path $buildRoot, $binDir -Force | Out-Null

Copy-Item -LiteralPath (Join-Path $officialRoot 'misc') -Destination $buildRoot -Recurse
Copy-Item -LiteralPath (Join-Path $officialRoot 'Bellhop') -Destination $buildRoot -Recurse
if (Test-Path -LiteralPath (Join-Path $officialRoot 'tslib')) {
    Copy-Item -LiteralPath (Join-Path $officialRoot 'tslib') -Destination $buildRoot -Recurse
}
Copy-Item -LiteralPath (Join-Path $supportDir 'Step.f90') -Destination (Join-Path $buildRoot 'Bellhop\Step.f90') -Force
Copy-Item -LiteralPath (Join-Path $supportDir 'bellhop.f90') -Destination (Join-Path $buildRoot 'Bellhop\bellhop.f90') -Force

if ([string]::IsNullOrWhiteSpace($CompilerRoot)) {
    $CompilerRoot = Join-Path $OutputDir 'toolchain_tmp\mingw64'
}
$compilerBin = Join-Path ([System.IO.Path]::GetFullPath($CompilerRoot)) 'bin'
$gfortran = Join-Path $compilerBin 'gfortran.exe'
foreach ($tool in @($gfortran, (Join-Path $compilerBin 'ar.exe'))) {
    if (-not (Test-Path -LiteralPath $tool -PathType Leaf)) { throw "Missing compiler tool: $tool" }
}

$oldPath = $env:Path
try {
    $env:Path = "$compilerBin;$oldPath"
    $compileFlags = @('-O2','-std=legacy','-ffast-math','-fno-range-check','-ffree-line-length-none','-static','-I../misc','-I../tslib')
    $miscSources = @('FatalError.f90','MathConstants.f90','pchipMod.f90','AttenMod.f90',
        'PolyMod.f90','RefCoef.f90','beampattern.f90','subtabulate.f90','calculateweights.f90',
        'monotonicMod.f90','SortMod.f90','SourceReceiverPositions.f90','sspMod.f90','RWSHDFile.f90',
        'interpolation.f90','MergeVectorsMod.f90','munk.f90','splinec.f90','norms.f90',
        'cross_products.f90','PekRoot.f90','RootFinderSecantMod.f90','ReadEnvironmentMod.f90')
    $miscObjects = @('FatalError.o','beampattern.o','MathConstants.o','RefCoef.o','SourceReceiverPositions.o',
        'pchipMod.o','AttenMod.o','sspMod.o','RWSHDFile.o','interpolation.o','MergeVectorsMod.o','munk.o',
        'ReadEnvironmentMod.o','SortMod.o','splinec.o','subtabulate.o','calculateweights.o','norms.o',
        'cross_products.o','monotonicMod.o','PolyMod.o','PekRoot.o','RootFinderSecantMod.o')
    Push-Location (Join-Path $buildRoot 'misc')
    try {
        foreach ($source in $miscSources) {
            & $gfortran -c @compileFlags $source
            if ($LASTEXITCODE -ne 0) { throw "Compilation failed: misc/$source" }
        }
        & (Join-Path $compilerBin 'ar.exe') rcs libmisc.a @miscObjects
        if ($LASTEXITCODE -ne 0) { throw 'Archiving libmisc.a failed.' }
    } finally { Pop-Location }

    $bellhopSources = @('bellhopMod.f90','angleMod.f90','ArrMod.f90','bdryMod.f90','sspMod.f90',
        'Cone.f90','ReflectMod.f90','WriteRay.f90','influence.f90','Step.f90','ReadEnvironmentBell.f90',
        'RayNormals.f90','bellhop.f90')
    $bellhopObjects = @('angleMod.o','ArrMod.o','bdryMod.o','bellhopMod.o','sspMod.o','ReflectMod.o',
        'WriteRay.o','influence.o','Step.o','ReadEnvironmentBell.o','bellhop.o')
    Push-Location (Join-Path $buildRoot 'Bellhop')
    try {
        foreach ($source in $bellhopSources) {
            & $gfortran -c @compileFlags $source
            if ($LASTEXITCODE -ne 0) { throw "Compilation failed: Bellhop/$source" }
        }
        & $gfortran -o bellhop.exe @compileFlags @bellhopObjects ..\misc\libmisc.a
        if ($LASTEXITCODE -ne 0) { throw 'Bellhop validation binary link failed.' }
    } finally { Pop-Location }
} finally {
    $env:Path = $oldPath
}

$builtExe = Join-Path $buildRoot 'Bellhop\bellhop.exe'
$validationExe = Join-Path $binDir 'bellhop_iwall_flat_2020.exe'
if (-not (Test-Path -LiteralPath $builtExe -PathType Leaf)) {
    throw "Expected binary was not produced: $builtExe"
}
Copy-Item -LiteralPath $builtExe -Destination $validationExe -Force

$manifest = [ordered]@{
    purpose = 'Bellhop 2020 validation-only flat internal wall POC'
    official_toolbox_root = $officialRoot
    official_bellhop_sha256 = (Get-FileHash -Algorithm SHA256 (Join-Path $officialRoot 'Bellhop\bellhop.f90')).Hash
    official_step_sha256 = (Get-FileHash -Algorithm SHA256 (Join-Path $officialRoot 'Bellhop\Step.f90')).Hash
    overlay_bellhop_sha256 = (Get-FileHash -Algorithm SHA256 (Join-Path $supportDir 'bellhop.f90')).Hash
    overlay_step_sha256 = (Get-FileHash -Algorithm SHA256 (Join-Path $supportDir 'Step.f90')).Hash
    validation_executable = $validationExe
    validation_executable_sha256 = (Get-FileHash -Algorithm SHA256 $validationExe).Hash
    compiler_root = [System.IO.Path]::GetFullPath($CompilerRoot)
    compiler_version = (& $gfortran --version | Select-Object -First 1)
}
$manifest | ConvertTo-Json | Set-Content -LiteralPath (Join-Path $OutputDir 'build_manifest.json') -Encoding utf8
Write-Output $validationExe
