# Shared native stdin transport for read-only status and local regression tests.
# PowerShell 7 / .NET ProcessStartInfo.ArgumentList; no command-text interpolation.
function Invoke-R1StatusStdin {
    param(
        [Parameter(Mandatory=$true)][string]$Executable,
        [Parameter(Mandatory=$true)][string[]]$ArgumentList,
        [Parameter(Mandatory=$true)][string]$InputText
    )
    $start=[Diagnostics.ProcessStartInfo]::new()
    $start.FileName=$Executable
    $start.UseShellExecute=$false
    $start.CreateNoWindow=$true
    $start.RedirectStandardInput=$true
    $start.RedirectStandardOutput=$true
    $start.RedirectStandardError=$true
    $encoding=[Text.UTF8Encoding]::new($false)
    $start.StandardInputEncoding=$encoding
    $start.StandardOutputEncoding=$encoding
    $start.StandardErrorEncoding=$encoding
    foreach($argument in $ArgumentList){$start.ArgumentList.Add($argument)}
    $process=[Diagnostics.Process]::new()
    $process.StartInfo=$start
    try {
        if(-not $process.Start()){throw 'Cannot start read-only status transport'}
        $stdout=$process.StandardOutput.ReadToEndAsync()
        $stderr=$process.StandardError.ReadToEndAsync()
        $process.StandardInput.Write($InputText.Replace("`r",''))
        $process.StandardInput.Close()
        $process.WaitForExit()
        return [pscustomobject]@{exit_code=$process.ExitCode;stdout=$stdout.GetAwaiter().GetResult();stderr=$stderr.GetAwaiter().GetResult()}
    } finally {
        $process.Dispose()
    }
}
