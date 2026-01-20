# MolVis Windows Release Script
# Builds, packages, and uploads a release to GitHub
#
# Usage:
#   .\release_windows.ps1 -Version "0.3.1"
#   .\release_windows.ps1 -Version "0.3.1" -Draft        # Create as draft
#   .\release_windows.ps1 -Version "0.3.1" -SkipBuild    # Skip build step
#   .\release_windows.ps1 -Version "0.3.1" -SkipUpload   # Build/package only

param(
    [Parameter(Mandatory=$true)]
    [string]$Version,

    [switch]$Draft = $false,
    [switch]$SkipBuild = $false,
    [switch]$SkipUpload = $false
)

$ErrorActionPreference = "Stop"
$ProjectRoot = $PSScriptRoot

Write-Host ""
Write-Host "╔═══════════════════════════════════════════════════════════╗" -ForegroundColor Cyan
Write-Host "║         MolVis Windows Release Script                     ║" -ForegroundColor Cyan
Write-Host "╚═══════════════════════════════════════════════════════════╝" -ForegroundColor Cyan
Write-Host ""
Write-Host "  Version:     $Version" -ForegroundColor White
Write-Host "  Draft:       $Draft" -ForegroundColor White
Write-Host "  Skip Build:  $SkipBuild" -ForegroundColor White
Write-Host "  Skip Upload: $SkipUpload" -ForegroundColor White
Write-Host ""

# Check for GitHub CLI if we're uploading
if (-not $SkipUpload) {
    $ghPath = Get-Command gh -ErrorAction SilentlyContinue
    if (-not $ghPath) {
        Write-Host "ERROR: GitHub CLI (gh) is not installed." -ForegroundColor Red
        Write-Host ""
        Write-Host "Install it with:" -ForegroundColor Yellow
        Write-Host "  winget install GitHub.cli" -ForegroundColor Gray
        Write-Host "  -or-" -ForegroundColor Gray
        Write-Host "  choco install gh" -ForegroundColor Gray
        Write-Host ""
        Write-Host "Then authenticate with:" -ForegroundColor Yellow
        Write-Host "  gh auth login" -ForegroundColor Gray
        Write-Host ""
        exit 1
    }

    # Check if authenticated
    $authStatus = gh auth status 2>&1
    if ($LASTEXITCODE -ne 0) {
        Write-Host "ERROR: Not authenticated with GitHub CLI." -ForegroundColor Red
        Write-Host "Run: gh auth login" -ForegroundColor Yellow
        exit 1
    }
}

# ============================================================
# Step 1: Build
# ============================================================
if (-not $SkipBuild) {
    Write-Host "═══════════════════════════════════════════════════════════" -ForegroundColor Blue
    Write-Host " Step 1: Building MolVis" -ForegroundColor Blue
    Write-Host "═══════════════════════════════════════════════════════════" -ForegroundColor Blue
    Write-Host ""

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
    Write-Host ""
    Write-Host "Build completed successfully!" -ForegroundColor Green
    Write-Host ""
} else {
    Write-Host "Skipping build step..." -ForegroundColor Yellow
}

# Verify molvis.exe exists
$ExePath = Join-Path $ProjectRoot "molvis.exe"
if (-not (Test-Path $ExePath)) {
    Write-Host "ERROR: molvis.exe not found at $ExePath" -ForegroundColor Red
    Write-Host "Run without -SkipBuild or build manually first." -ForegroundColor Yellow
    exit 1
}

# ============================================================
# Step 2: Package
# ============================================================
Write-Host "═══════════════════════════════════════════════════════════" -ForegroundColor Blue
Write-Host " Step 2: Creating Release Package" -ForegroundColor Blue
Write-Host "═══════════════════════════════════════════════════════════" -ForegroundColor Blue
Write-Host ""

$OutputDir = Join-Path $ProjectRoot "release"
$ReleaseDir = Join-Path $OutputDir "MolVis-$Version-Windows-x64"
$ZipName = "MolVis-$Version-Windows-x64.zip"
$ZipPath = Join-Path $OutputDir $ZipName
$ChecksumName = "$ZipName.sha256"
$ChecksumPath = Join-Path $OutputDir $ChecksumName

# Clean previous release
if (Test-Path $OutputDir) {
    Write-Host "Cleaning previous release directory..." -ForegroundColor Yellow
    Remove-Item -Recurse -Force $OutputDir
}

# Create release directory
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

# Copy documentation
Write-Host "  Copying documentation..." -ForegroundColor Gray
$DocFiles = @("README.md", "LICENSE", "RELEASE_NOTES.md")
foreach ($doc in $DocFiles) {
    $docPath = Join-Path $ProjectRoot $doc
    if (Test-Path $docPath) {
        Copy-Item $docPath $ReleaseDir
    }
}

# Copy docs folder (excluding internal docs)
$DocsDir = Join-Path $ProjectRoot "docs"
if (Test-Path $DocsDir) {
    $DocsFiles = Get-ChildItem $DocsDir -File | Where-Object { $_.Name -notmatch "PLAN|TODO|INTERNAL" }
    if ($DocsFiles.Count -gt 0) {
        $ReleaseDocsDir = Join-Path $ReleaseDir "docs"
        New-Item -ItemType Directory -Path $ReleaseDocsDir -Force | Out-Null
        foreach ($file in $DocsFiles) {
            Write-Host "  Copying docs/$($file.Name)..." -ForegroundColor Gray
            Copy-Item $file.FullName $ReleaseDocsDir
        }
    }
}

# Copy CUDA runtime DLLs
Write-Host "  Locating CUDA runtime DLLs..." -ForegroundColor Gray
$CudaPath = $env:CUDA_PATH
if ($CudaPath) {
    $CudaDllDir = Join-Path $CudaPath "bin"
    $CudaDlls = Get-ChildItem -Path $CudaDllDir -Filter "cudart64_*.dll" -ErrorAction SilentlyContinue
    foreach ($dll in $CudaDlls) {
        Write-Host "    Copying $($dll.Name)..." -ForegroundColor Gray
        Copy-Item $dll.FullName $ReleaseDir
    }
} else {
    Write-Host "  WARNING: CUDA_PATH not set, CUDA DLLs not included" -ForegroundColor Yellow
    Write-Host "  Users will need CUDA runtime installed" -ForegroundColor Yellow
}

# Create ZIP archive
Write-Host "  Creating ZIP archive..." -ForegroundColor Gray
Compress-Archive -Path $ReleaseDir -DestinationPath $ZipPath -Force

# Calculate SHA256 checksum
Write-Host "  Calculating SHA256 checksum..." -ForegroundColor Gray
$Hash = (Get-FileHash -Path $ZipPath -Algorithm SHA256).Hash.ToLower()
"$Hash  $ZipName" | Out-File -FilePath $ChecksumPath -Encoding utf8 -NoNewline

Write-Host ""
Write-Host "Package created successfully!" -ForegroundColor Green
Write-Host "  Archive:  $ZipPath" -ForegroundColor White
Write-Host "  Checksum: $ChecksumPath" -ForegroundColor White
Write-Host "  SHA256:   $Hash" -ForegroundColor Cyan
Write-Host ""

# List package contents
Write-Host "Package Contents:" -ForegroundColor Yellow
Get-ChildItem -Path $ReleaseDir -Recurse | ForEach-Object {
    $relativePath = $_.FullName.Substring($ReleaseDir.Length + 1)
    $size = if ($_.PSIsContainer) { "[DIR]" } else { "{0:N0} KB" -f ($_.Length / 1KB) }
    Write-Host "  $relativePath ($size)" -ForegroundColor Gray
}
Write-Host ""

if ($SkipUpload) {
    Write-Host "Skipping GitHub upload (use without -SkipUpload to publish)." -ForegroundColor Yellow
    Write-Host ""
    Write-Host "Files ready for manual upload:" -ForegroundColor White
    Write-Host "  $ZipPath" -ForegroundColor Gray
    Write-Host "  $ChecksumPath" -ForegroundColor Gray
    exit 0
}

# ============================================================
# Step 3: Create Git Tag
# ============================================================
Write-Host "═══════════════════════════════════════════════════════════" -ForegroundColor Blue
Write-Host " Step 3: Creating Git Tag" -ForegroundColor Blue
Write-Host "═══════════════════════════════════════════════════════════" -ForegroundColor Blue
Write-Host ""

$TagName = "v$Version"

# Check if tag already exists locally
$existingTag = git tag -l $TagName 2>$null
if ($existingTag) {
    Write-Host "Tag $TagName already exists locally." -ForegroundColor Yellow
    $confirm = Read-Host "Delete and recreate? (y/N)"
    if ($confirm -eq 'y' -or $confirm -eq 'Y') {
        git tag -d $TagName 2>$null
        # Also delete remote tag if it exists
        git push origin --delete $TagName 2>$null
    } else {
        Write-Host "Using existing tag $TagName" -ForegroundColor Yellow
    }
}

# Create tag if it doesn't exist
$existingTag = git tag -l $TagName 2>$null
if (-not $existingTag) {
    Write-Host "Creating tag $TagName..." -ForegroundColor Gray
    git tag -a $TagName -m "Release $TagName"
    if ($LASTEXITCODE -ne 0) {
        Write-Host "ERROR: Failed to create tag" -ForegroundColor Red
        exit 1
    }
}

# Push tag to remote
Write-Host "Pushing tag to GitHub..." -ForegroundColor Gray
git push origin $TagName 2>&1
if ($LASTEXITCODE -ne 0) {
    Write-Host "WARNING: Tag may already exist on remote" -ForegroundColor Yellow
}

Write-Host ""

# ============================================================
# Step 4: Create GitHub Release
# ============================================================
Write-Host "═══════════════════════════════════════════════════════════" -ForegroundColor Blue
Write-Host " Step 4: Creating GitHub Release" -ForegroundColor Blue
Write-Host "═══════════════════════════════════════════════════════════" -ForegroundColor Blue
Write-Host ""

# Build release notes body
$ReleaseBody = @"
## MolVis v$Version

GPU-accelerated molecular visualization for Windows (CUDA + DirectX 11).

### Downloads

| Platform | File | Notes |
|----------|------|-------|
| Windows x64 | ``$ZipName`` | Requires NVIDIA GPU with CUDA support |

### Requirements

- Windows 10/11
- NVIDIA GPU with CUDA support (Compute Capability 5.0+)
- CUDA drivers (included in recent NVIDIA Game Ready drivers)

### Installation

1. Download the ZIP file
2. Extract to a folder of your choice
3. Run ``molvis.exe``

### Verify Download

```powershell
(Get-FileHash $ZipName -Algorithm SHA256).Hash.ToLower()
```
Expected: ``$Hash``

---
See [RELEASE_NOTES.md](RELEASE_NOTES.md) for detailed changes.
"@

# Check if release already exists
$existingRelease = gh release view $TagName 2>$null
if ($LASTEXITCODE -eq 0) {
    Write-Host "Release $TagName already exists." -ForegroundColor Yellow
    $confirm = Read-Host "Delete and recreate? (y/N)"
    if ($confirm -eq 'y' -or $confirm -eq 'Y') {
        Write-Host "Deleting existing release..." -ForegroundColor Gray
        gh release delete $TagName --yes 2>$null
    } else {
        Write-Host "Uploading assets to existing release..." -ForegroundColor Yellow
        gh release upload $TagName $ZipPath $ChecksumPath --clobber
        Write-Host ""
        Write-Host "Assets uploaded successfully!" -ForegroundColor Green
        exit 0
    }
}

# Create the release
Write-Host "Creating GitHub release..." -ForegroundColor Gray

$ghArgs = @(
    "release", "create", $TagName,
    $ZipPath,
    $ChecksumPath,
    "--title", "MolVis v$Version",
    "--notes", $ReleaseBody
)

if ($Draft) {
    $ghArgs += "--draft"
}

& gh @ghArgs

if ($LASTEXITCODE -ne 0) {
    Write-Host "ERROR: Failed to create GitHub release" -ForegroundColor Red
    exit 1
}

Write-Host ""
Write-Host "╔═══════════════════════════════════════════════════════════╗" -ForegroundColor Green
Write-Host "║              Release Published Successfully!              ║" -ForegroundColor Green
Write-Host "╚═══════════════════════════════════════════════════════════╝" -ForegroundColor Green
Write-Host ""
Write-Host "  View release at:" -ForegroundColor White
gh release view $TagName --web 2>$null
Write-Host ""
