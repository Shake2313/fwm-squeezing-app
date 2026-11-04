param(
    [string]$LeanBin = (Join-Path ([System.IO.Path]::GetTempPath()) 'gabes-lean-tools\lean-4.35.0-windows\bin'),
    [string]$OutputDir = (Join-Path $PSScriptRoot '..\..\..\docs\grand_challenge\lean_pilot_2026_10_28\run')
)
$ErrorActionPreference = 'Stop'
# The host may provide a minimal environment; change only this process.
if ($env:PATHEXT -notmatch '(?i)(^|;)\.EXE(;|$)') { $env:PATHEXT = '.EXE;' + $env:PATHEXT }
[Console]::OutputEncoding = New-Object System.Text.UTF8Encoding($false)
$utf8 = New-Object System.Text.UTF8Encoding($false)
$LeanBin = [System.IO.Path]::GetFullPath($LeanBin)
$OutputDir = [System.IO.Path]::GetFullPath($OutputDir)
New-Item -ItemType Directory -Path $OutputDir -Force | Out-Null
$lake = Join-Path $LeanBin 'lake.exe'
$lean = Join-Path $LeanBin 'lean.exe'
$version = (& $lean --version | Out-String).Trim()
if ($LASTEXITCODE -ne 0 -or $version -notmatch 'version 4\.35\.0,') {
    throw 'Expected the pinned Lean 4.35.0 toolchain.'
}
$checks = New-Object System.Collections.ArrayList
function Invoke-LeanCheck([string]$Name, [string[]]$CommandArgs, [bool]$ExpectFailure) {
    $clock = [Diagnostics.Stopwatch]::StartNew()
    $previousPreference = $ErrorActionPreference
    $ErrorActionPreference = 'Continue'
    $text = (& $lake @CommandArgs 2>&1 | Out-String)
    $code = $LASTEXITCODE
    $ErrorActionPreference = $previousPreference
    $clock.Stop()
    [System.IO.File]::WriteAllText((Join-Path $OutputDir ($Name + '.log')), $text, $utf8)
    $valid = ($code -eq 0)
    if ($ExpectFailure) {
        $valid = ($code -ne 0 -and $text.Contains('error: `grind` failed') -and
                  $text.Contains('[eqc] False propositions') -and
                  $text.Contains('negative/') -and
                  $text -notmatch 'unknown (module|identifier)|failed to import|file not found')
    } else {
        $valid = $valid -and ($text -notmatch 'sorryAx|warning:|error:')
    }
    [void]$checks.Add([ordered]@{
        name = $Name; args = $CommandArgs; exit_code = $code
        seconds = [Math]::Round($clock.Elapsed.TotalSeconds, 4)
        expected_failure = $ExpectFailure; passed = $valid
    })
    if (-not $valid) { throw "Unexpected result in $Name. Inspect its saved log." }
    return $text
}
Push-Location $PSScriptRoot
try {
    $buildOutput = Invoke-LeanCheck 'build' @('build') $false
    $proofOutput = Invoke-LeanCheck 'proof_audit' @('env', 'lean', 'Main.lean') $false
    $source = [System.IO.File]::ReadAllText((Join-Path $PSScriptRoot 'Main.lean'))
    $expected = @([regex]::Matches($source, '(?m)^theorem (\w+)') | ForEach-Object { 'GabesPilot.' + $_.Groups[1].Value })
    $auditMatches = [regex]::Matches($proofOutput, "'(GabesPilot\.\w+)' depends on axioms:\s*\[([^\]]*)\]")
    $audits = @($auditMatches | ForEach-Object {
        $names = @($_.Groups[2].Value.Split(',') | ForEach-Object { $_.Trim() } | Where-Object { $_ })
        $bad = @($names | Where-Object { $_ -notin @('propext','Classical.choice','Quot.sound') })
        if ($bad.Count -gt 0) { throw "Unapproved axioms: $bad" }
        [ordered]@{ theorem = $_.Groups[1].Value; axioms = $names }
    })
    $actual = @($audits | ForEach-Object { $_.theorem })
    if ($expected.Count -ne 13 -or $actual.Count -ne $expected.Count -or
        @(Compare-Object $expected $actual).Count -ne 0) { throw 'Incomplete theorem audit.' }
    $naive = Invoke-LeanCheck 'negative_naive' @('env','lean','negative/NaiveNoise.lean') $true
    $missing = Invoke-LeanCheck 'negative_missing' @('env','lean','negative/MissingNoise.lean') $true
    if (-not $naive.Contains('5 * 8') -or -not $missing.Contains('3 * 8')) {
        throw 'Expected rational counterexamples missing from diagnostics.'
    }
    $hashes = [ordered]@{}
    foreach ($relative in @('Main.lean','negative/NaiveNoise.lean','negative/MissingNoise.lean','lean-toolchain','lakefile.toml','lake-manifest.json','setup.ps1','run_pilot.ps1')) {
        $hashes[$relative] = (Get-FileHash -LiteralPath (Join-Path $PSScriptRoot $relative) -Algorithm SHA256).Hash.ToLowerInvariant()
    }
    $receipt = [ordered]@{
        schema = 'gabes-lean-pilot-v1'
        timestamp_utc = [DateTime]::UtcNow.ToString('o')
        status = 'PASS'; adoption = 'pilot_only_pending_user_decision'
        lean_version = $version; lean_bin = $LeanBin
        checks = @($checks); theorem_count = $expected.Count; axiom_audit = $audits
        source_sha256 = $hashes
        scope = 'Exact rational scalar covariance algebra; no general CP, real-matrix, Python implementation or GC milestone certification.'
    }
    [System.IO.File]::WriteAllText((Join-Path $OutputDir 'receipt.json'), ($receipt | ConvertTo-Json -Depth 10) + [Environment]::NewLine, $utf8)
    Write-Output "PASS: 13 theorem axiom audits; 2 expected mathematical rejections. $OutputDir"
} finally { Pop-Location }
