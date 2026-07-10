@echo off
REM ============================================================================
REM MolVis Build Script
REM ============================================================================
REM Sets up MSVC environment and runs nmake
REM
REM Usage: build.bat [target]
REM   build.bat          - Build the application
REM   build.bat run      - Build and run
REM   build.bat clean    - Clean build artifacts
REM   build.bat legacy   - Build original CUDA-only version
REM   build.bat help     - Show all targets
REM ============================================================================

REM vcvarsall may place exclamation marks in environment values. Delayed
REM expansion would corrupt those values after returning from the script.
setlocal EnableExtensions DisableDelayedExpansion

echo.
echo  ╔═══════════════════════════════════════╗
echo  ║  MolVis - Molecular Visualization     ║
echo  ║  CUDA + Dear ImGui + DirectX 11       ║
echo  ╚═══════════════════════════════════════╝
echo.

REM Try to find Visual Studio installation
set "VCVARSALL="
set "VCVARS_TOOLSET_ARG="

REM Check Visual Studio 2025 Insiders
if exist "C:\Program Files\Microsoft Visual Studio\18\Insiders\VC\Auxiliary\Build\vcvarsall.bat" (
    set "VCVARSALL=C:\Program Files\Microsoft Visual Studio\18\Insiders\VC\Auxiliary\Build\vcvarsall.bat"
    REM CUDA 13.1 crashes in cudafe++ with the preview 14.51 toolset. Use the
    REM installed stable 14.44 toolset when available.
    if exist "C:\Program Files\Microsoft Visual Studio\18\Insiders\VC\Tools\MSVC\14.44.35207\bin\Hostx64\x64\cl.exe" set "VCVARS_TOOLSET_ARG=-vcvars_ver=14.44"
    echo Found: Visual Studio 2025 Insiders
    goto :found_vs
)

REM Check Visual Studio 2022 Community
if exist "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvarsall.bat" (
    set "VCVARSALL=C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvarsall.bat"
    echo Found: Visual Studio 2022 Community
    goto :found_vs
)

REM Check Visual Studio 2022 Professional
if exist "C:\Program Files\Microsoft Visual Studio\2022\Professional\VC\Auxiliary\Build\vcvarsall.bat" (
    set "VCVARSALL=C:\Program Files\Microsoft Visual Studio\2022\Professional\VC\Auxiliary\Build\vcvarsall.bat"
    echo Found: Visual Studio 2022 Professional
    goto :found_vs
)

REM Check Visual Studio 2022 Enterprise
if exist "C:\Program Files\Microsoft Visual Studio\2022\Enterprise\VC\Auxiliary\Build\vcvarsall.bat" (
    set "VCVARSALL=C:\Program Files\Microsoft Visual Studio\2022\Enterprise\VC\Auxiliary\Build\vcvarsall.bat"
    echo Found: Visual Studio 2022 Enterprise
    goto :found_vs
)

REM Check Visual Studio 2022 Build Tools
if exist "C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvarsall.bat" (
    set "VCVARSALL=C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvarsall.bat"
    echo Found: Visual Studio 2022 Build Tools
    goto :found_vs
)

REM Check Visual Studio 2019
if exist "C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Auxiliary\Build\vcvarsall.bat" (
    set "VCVARSALL=C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Auxiliary\Build\vcvarsall.bat"
    echo Found: Visual Studio 2019 Community
    goto :found_vs
)

echo ERROR: Could not find Visual Studio installation.
echo.
echo CUDA on Windows requires MSVC as the host compiler.
echo Please install Visual Studio 2019, 2022, or later with the
echo "Desktop development with C++" workload.
echo.
echo You can also install just the Build Tools from:
echo https://visualstudio.microsoft.com/downloads/#build-tools-for-visual-studio-2022
exit /b 1

:found_vs
echo.
echo Setting up x64 build environment...
REM A build launched from an IDE may already contain a different VS developer
REM environment. Clear its selection markers so vcvarsall honors our toolset.
set "VSCMD_ARG_TGT_ARCH="
set "VSCMD_VER="
set "VCToolsVersion="
set "VCToolsInstallDir="
call "%VCVARSALL%" x64 %VCVARS_TOOLSET_ARG% >nul 2>&1

if errorlevel 1 (
    echo ERROR: Failed to set up MSVC environment.
    exit /b 1
)

REM IDE terminals can inherit a very long PATH containing several duplicate
REM Visual Studio environments. nvcc prepends more tool directories before it
REM launches cl.exe; on Windows this can overflow the child environment and
REM make an installed compiler appear to be missing. vcvarsall has already set
REM INCLUDE and LIB, so retain only the tools required by this build.
set "PATH=%VCToolsInstallDir%bin\Hostx64\x64;%CUDA_PATH%\bin;%WindowsSdkDir%bin\%WindowsSDKVersion%x64;C:\Windows\System32;C:\Windows"

echo Environment configured successfully.
echo.

where cl.exe >nul 2>&1
if errorlevel 1 (
    echo ERROR: MSVC environment initialization completed, but cl.exe is not in PATH.
    echo Visual Studio installation: %VCVARSALL%
    exit /b 1
)

where nmake.exe >nul 2>&1
if errorlevel 1 (
    echo ERROR: MSVC environment initialization completed, but nmake.exe is not in PATH.
    exit /b 1
)

REM Check for CUDA
where nvcc >nul 2>&1
if errorlevel 1 (
    echo ERROR: nvcc not found in PATH.
    echo Please install CUDA Toolkit and ensure it's in your PATH.
    exit /b 1
)

echo CUDA compiler found:
nvcc --version | findstr /C:"release"
echo.

REM Change to the script's directory
cd /d "%~dp0"

REM Kill any running instance before building
taskkill /F /IM molvis.exe >nul 2>&1

REM Ensure build directory exists
if not exist build mkdir build

REM Run make with the provided target (or 'all' by default)
if "%~1"=="" (
    echo Building MolVis...
    echo.
    nmake /nologo /f Makefile all
) else (
    echo Target: %~1
    echo.
    nmake /nologo /f Makefile %~1
)

if errorlevel 1 (
    echo.
    echo Build failed!
    exit /b 1
)

echo.
echo Build complete!
