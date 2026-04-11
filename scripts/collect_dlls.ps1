# collect_dlls.ps1 - Copy required MinGW DLLs to bin directory
# Called from build_windows.bat after successful build

$ErrorActionPreference = "SilentlyContinue"

# Find MinGW bin directory: prefer MINGW_DIR env var, fall back to PATH lookup, then common locations
$mingwBin = $null

if ($env:MINGW_DIR -and (Test-Path "$env:MINGW_DIR\bin")) {
    $mingwBin = "$env:MINGW_DIR\bin"
} else {
    # Try to find g++.exe in PATH (set by setup-mingw action)
    $gpp = Get-Command g++.exe -ErrorAction SilentlyContinue
    if ($gpp) {
        $mingwBin = Split-Path $gpp.Source
    } elseif (Test-Path "C:\msys64\mingw64\bin") {
        $mingwBin = "C:\msys64\mingw64\bin"
    } elseif (Test-Path "C:\mingw64\bin") {
        $mingwBin = "C:\mingw64\bin"
    }
}

if (-not $mingwBin) {
    Write-Host "[WARN] MinGW bin directory not found. Skipping DLL collection."
    exit 0
}

Write-Host "Using MinGW from: $mingwBin"

# DLL patterns required by curcuma
$patterns = @(
    'libgcc_s_seh*',
    'libgfortran*',
    'libgomp*',
    'libquadmath*',
    'libstdc++*',
    'libwinpthread*'
)

# bin directory is at project root / bin (one level up from scripts/)
$binDir = Join-Path $PSScriptRoot "..\bin"
if (-not (Test-Path $binDir)) {
    New-Item -ItemType Directory -Path $binDir -Force | Out-Null
}

$copied = 0
foreach ($pattern in $patterns) {
    $files = Get-ChildItem -Path $mingwBin -Filter $pattern -ErrorAction SilentlyContinue
    foreach ($file in $files) {
        Copy-Item $file.FullName -Destination $binDir -Force
        Write-Host "  Copied: $($file.Name)"
        $copied++
    }
}

if ($copied -eq 0) {
    Write-Host "[WARN] No DLLs found in $mingwBin matching required patterns."
} else {
    Write-Host "Collected $copied DLL(s)."
}