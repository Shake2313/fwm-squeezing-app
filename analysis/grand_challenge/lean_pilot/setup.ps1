param(
    [string]$InstallRoot = (Join-Path ([System.IO.Path]::GetTempPath()) 'gabes-lean-tools')
)
$ErrorActionPreference = 'Stop'
# The host may provide a minimal environment; change only this process.
if ($env:PATHEXT -notmatch '(?i)(^|;)\.EXE(;|$)') { $env:PATHEXT = '.EXE;' + $env:PATHEXT }
$version = '4.35.0'
$expectedHash = '5fc033e93e1a7bc59c5ebf311103d2c9cd18d4b3eddda4ffa3414e8adb0ef79e'
$url = "https://github.com/leanprover/lean4/releases/download/v$version/lean-$version-windows.zip"
New-Item -ItemType Directory -Path $InstallRoot -Force | Out-Null
$archive = Join-Path $InstallRoot "lean-$version-windows.zip"
$bin = Join-Path $InstallRoot "lean-$version-windows\bin"
if (-not (Test-Path -LiteralPath $archive)) {
    & curl.exe -fL --retry 2 --connect-timeout 20 -o $archive $url
    if ($LASTEXITCODE -ne 0) { throw 'Lean download failed.' }
}
$actualHash = (Get-FileHash -LiteralPath $archive -Algorithm SHA256).Hash.ToLowerInvariant()
if ($actualHash -ne $expectedHash) { throw 'Lean archive SHA-256 mismatch.' }
if (-not (Test-Path -LiteralPath (Join-Path $bin 'lean.exe'))) {
    & tar.exe -xf $archive -C $InstallRoot
    if ($LASTEXITCODE -ne 0) { throw 'Lean extraction failed.' }
}
& (Join-Path $bin 'lean.exe') --version
if ($LASTEXITCODE -ne 0) { throw 'Lean executable failed.' }
Write-Output "Portable Lean ready: $bin"
Write-Output 'No global PATH, CI or production dependencies changed.'
