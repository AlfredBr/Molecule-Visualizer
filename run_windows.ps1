Set-Alias run ./run_windows.ps1

Remove-Item -Force ./molvis.exe -ErrorAction SilentlyContinue

./build.bat

if (Test-Path "./molvis.exe") {
	Write-Host "Launching molvis.exe..."
	Start-Process -FilePath "./molvis.exe"
} else {
	Write-Host "molvis.exe not found. Please build the project first."
}
