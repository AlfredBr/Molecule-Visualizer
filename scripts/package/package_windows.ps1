# MolVis Windows Release Packaging Script
# Creates a distributable release package

param(
    [string]$Version = "1.0.0",
    [switch]$Build = $false
)

$ErrorActionPreference = "Stop"

Write-Host "=== MolVis Windows Release Packager ===" -ForegroundColor Cyan
Write-Host "Version: $Version" -ForegroundColor Yellow

# Paths
$ScriptDir = $PSScriptRoot
$ProjectRoot = (Resolve-Path "$ScriptDir/../..").Path
$BuildDir = Join-Path $ProjectRoot "build_release"
$OutputDir = Join-Path $ProjectRoot "release"
$ReleaseDir = Join-Path $OutputDir "MolVis-$Version-Windows"

# Clean previous release
if (Test-Path $OutputDir) {
    Write-Host "Cleaning previous release directory..." -ForegroundColor Yellow
    Remove-Item -Recurse -Force $OutputDir
}

# Build if requested
if ($Build) {
    Write-Host "Building MolVis..." -ForegroundColor Yellow
    Push-Location $ProjectRoot
    try {
        & .\build.bat
        if ($LASTEXITCODE -ne 0) {
            throw "Build failed with exit code $LASTEXITCODE"
        }
    }
    finally {
        Pop-Location
    }
}

# Check that molvis.exe exists
$ExePath = Join-Path $ProjectRoot "molvis.exe"
if (-not (Test-Path $ExePath)) {
    Write-Host "Error: molvis.exe not found. Run with -Build flag or build manually first." -ForegroundColor Red
    exit 1
}

# Create release directory structure
Write-Host "Creating release package..." -ForegroundColor Yellow
New-Item -ItemType Directory -Path $ReleaseDir -Force | Out-Null

# Copy main executable
Write-Host "  Copying molvis.exe..." -ForegroundColor Gray
Copy-Item $ExePath $ReleaseDir

# Copy configuration files
$ConfigFiles = @("molvis.ini", "imgui.ini")
foreach ($cfg in $ConfigFiles) {
    $cfgPath = Join-Path $ProjectRoot $cfg
    if (Test-Path $cfgPath) {
        Write-Host "  Copying $cfg..." -ForegroundColor Gray
        Copy-Item $cfgPath $ReleaseDir
    }
}

# Create ZIP archive
Write-Host "Creating ZIP archive..." -ForegroundColor Yellow
$ZipPath = Join-Path $OutputDir "MolVis-$Version-Windows-x64.zip"
Compress-Archive -Path $ReleaseDir -DestinationPath $ZipPath -Force

# Calculate SHA256 checksum
Write-Host "Calculating SHA256 checksum..." -ForegroundColor Yellow
$Hash = (Get-FileHash -Path $ZipPath -Algorithm SHA256).Hash.ToLower()

# Create checksum file
$ChecksumFile = Join-Path $OutputDir "MolVis-$Version-Windows-x64.zip.sha256"
"$Hash  MolVis-$Version-Windows-x64.zip" | Out-File -FilePath $ChecksumFile -Encoding utf8 -NoNewline

# Summary
Write-Host ""
Write-Host "=== Release Package Created ===" -ForegroundColor Green
Write-Host "Directory: $ReleaseDir" -ForegroundColor White
Write-Host "Archive:   $ZipPath" -ForegroundColor White
Write-Host "Checksum:  $ChecksumFile" -ForegroundColor White
Write-Host ""
Write-Host "SHA256: $Hash" -ForegroundColor Cyan
Write-Host ""

# List package contents
Write-Host "Package Contents:" -ForegroundColor Yellow
Get-ChildItem -Path $ReleaseDir -Recurse | ForEach-Object {
    $relativePath = $_.FullName.Substring($ReleaseDir.Length + 1)
    $size = if ($_.PSIsContainer) { "[DIR]" } else { "{0:N0} KB" -f ($_.Length / 1KB) }
    Write-Host "  $relativePath - $size" -ForegroundColor Gray
}

Write-Host ""
Write-Host "To upload to GitHub Release:" -ForegroundColor Yellow
Write-Host "  1. Go to https://github.com/YOUR_USERNAME/MolVis/releases/new" -ForegroundColor Gray
Write-Host "  2. Create a new tag: v$Version" -ForegroundColor Gray
Write-Host "  3. Upload: $ZipPath" -ForegroundColor Gray
Write-Host "  4. Upload: $ChecksumFile" -ForegroundColor Gray
Write-Host ""
