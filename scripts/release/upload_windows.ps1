# Uploads the Windows release assets to an existing GitHub release using GitHub CLI
# Usage: .\scripts\release\upload_windows.ps1 -Version "0.3.6"

param(
    [Parameter(Mandatory=$true)]
    [string]$Version
)

$ErrorActionPreference = "Stop"

$ScriptDir = $PSScriptRoot
$ProjectRoot = (Resolve-Path "$ScriptDir/../..").Path
$Tag = "v$Version"
$ReleaseDir = Join-Path $ProjectRoot "release"
$Zip = Join-Path $ReleaseDir "MolVis-$Version-Windows-x64.zip"
$Sha = Join-Path $ReleaseDir "MolVis-$Version-Windows-x64.zip.sha256"

if (!(Test-Path $Zip)) {
    Write-Host "ERROR: $Zip not found. Build/package first." -ForegroundColor Red
    exit 1
}
if (!(Test-Path $Sha)) {
    Write-Host "ERROR: $Sha not found. Build/package first." -ForegroundColor Red
    exit 1
}

# Check for GitHub CLI
if (-not (Get-Command gh -ErrorAction SilentlyContinue)) {
    Write-Host "ERROR: GitHub CLI (gh) is not installed." -ForegroundColor Red
    Write-Host "Install with: winget install GitHub.cli" -ForegroundColor Yellow
    exit 1
}

Write-Host "Uploading Windows assets to release $Tag..." -ForegroundColor Cyan

$cmd = @(
    'gh', 'release', 'upload', $Tag,
    $Zip, $Sha,
    '--clobber'
)

$proc = Start-Process -FilePath $cmd[0] -ArgumentList $cmd[1..($cmd.Length-1)] -NoNewWindow -Wait -PassThru
if ($proc.ExitCode -eq 0) {
    Write-Host "Upload complete!" -ForegroundColor Green
} else {
    Write-Host "Upload failed." -ForegroundColor Red
    exit $proc.ExitCode
}
