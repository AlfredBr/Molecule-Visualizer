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
$ProjectRoot = $PSScriptRoot
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

# Copy configuration file (default settings)
$IniPath = Join-Path $ProjectRoot "molvis.ini"
if (Test-Path $IniPath) {
    Write-Host "  Copying molvis.ini..." -ForegroundColor Gray
    Copy-Item $IniPath $ReleaseDir
}

# Copy ImGui configuration
$ImGuiIni = Join-Path $ProjectRoot "imgui.ini"
if (Test-Path $ImGuiIni) {
    Write-Host "  Copying imgui.ini..." -ForegroundColor Gray
    Copy-Item $ImGuiIni $ReleaseDir
}

# Copy documentation
Write-Host "  Copying documentation..." -ForegroundColor Gray
Copy-Item (Join-Path $ProjectRoot "README.md") $ReleaseDir
Copy-Item (Join-Path $ProjectRoot "LICENSE") $ReleaseDir
Copy-Item (Join-Path $ProjectRoot "RELEASE_NOTES.md") $ReleaseDir

# Copy docs folder if it exists and has content
$DocsDir = Join-Path $ProjectRoot "docs"
if (Test-Path $DocsDir) {
    $DocsFiles = Get-ChildItem $DocsDir -File
    if ($DocsFiles.Count -gt 0) {
        $ReleaseDocsDir = Join-Path $ReleaseDir "docs"
        New-Item -ItemType Directory -Path $ReleaseDocsDir -Force | Out-Null
        foreach ($file in $DocsFiles) {
            # Skip internal development docs
            if ($file.Name -notmatch "PLAN|TODO|INTERNAL") {
                Copy-Item $file.FullName $ReleaseDocsDir
            }
        }
    }
}

# Find and copy required CUDA runtime DLLs
Write-Host "  Locating CUDA runtime DLLs..." -ForegroundColor Gray
$CudaPath = $env:CUDA_PATH
if ($CudaPath) {
    $CudaDllDir = Join-Path $CudaPath "bin"
    $RequiredDlls = @(
        "cudart64_*.dll"
    )
    
    foreach ($pattern in $RequiredDlls) {
        $dlls = Get-ChildItem -Path $CudaDllDir -Filter $pattern -ErrorAction SilentlyContinue
        foreach ($dll in $dlls) {
            Write-Host "    Copying $($dll.Name)..." -ForegroundColor Gray
            Copy-Item $dll.FullName $ReleaseDir
        }
    }
} else {
    Write-Host "  Warning: CUDA_PATH not set, skipping CUDA DLL copy" -ForegroundColor Yellow
    Write-Host "  Users will need CUDA runtime installed" -ForegroundColor Yellow
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
