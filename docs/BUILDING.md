# Building, Running, Packaging, and Deploying MolVis

This guide covers the complete workflow for building MolVis from source on both macOS and Windows.

## Table of Contents

- [macOS](#macos)
  - [Prerequisites](#macos-prerequisites)
  - [Building](#building-on-macos)
  - [Running](#running-on-macos)
  - [Packaging](#packaging-on-macos)
  - [Code Signing & Notarization](#code-signing--notarization)
  - [Deploying](#deploying-on-macos)
- [Windows](#windows)
  - [Prerequisites](#windows-prerequisites)
  - [Building](#building-on-windows)
  - [Running](#running-on-windows)
  - [Packaging](#packaging-on-windows)
  - [Deploying](#deploying-on-windows)

---

## macOS

### macOS Prerequisites

1. **Xcode Command Line Tools**
   ```bash
   xcode-select --install
   ```

2. **CMake** (3.20 or later)
   ```bash
   brew install cmake
   ```

3. **SDL2**
   ```bash
   brew install sdl2
   ```

4. **Metal Toolchain** (usually included with Xcode)
   ```bash
   # Verify Metal compiler is available
   xcrun -find metal
   ```

### Building on macOS

**Quick build:**
```bash
./build_mac.sh
```

**Manual build:**
```bash
mkdir -p build_mac && cd build_mac
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . --config Release -j$(sysctl -n hw.ncpu)
```

**Clean build:**
```bash
./build_mac.sh clean
```

The built application will be at: `build_mac/bin/MolVis.app`

### Running on macOS

**From build script:**
```bash
./build_mac.sh run
```

**Directly:**
```bash
open build_mac/bin/MolVis.app
# or
./build_mac/bin/MolVis.app/Contents/MacOS/MolVis
```

### Packaging on macOS

Create a distributable DMG file:

```bash
./scripts/package/package_mac.sh
```

The DMG will be created at: `dist/MolVis-<version>-macOS.dmg`

### Code Signing & Notarization

> **Current Status:** The app is built with an ad-hoc signature, which means users must bypass Gatekeeper to run it (right-click → Open, or `xattr -cr MolVis.app`).

#### Why This Matters

macOS Gatekeeper blocks unsigned or improperly signed apps downloaded from the internet. For users to install without workarounds, the app must be:

1. **Signed** with an Apple Developer ID Application certificate
2. **Notarized** by Apple's notary service
3. **Stapled** with the notarization ticket

#### Requirements for Proper Signing

1. **Apple Developer Program membership** ($99/year)
   - Sign up at: https://developer.apple.com/programs/

2. **Developer ID Application certificate**
   - Create at: https://developer.apple.com/account/resources/certificates
   - Select "Developer ID Application" certificate type
   - Download and install in your Keychain

3. **App-specific password** (for notarization)
   - Create at: https://appleid.apple.com/account/manage
   - Under "Security" → "App-Specific Passwords"

4. **Store credentials in Keychain** (recommended):
   ```bash
   xcrun notarytool store-credentials "MolVis-Notarize" \
       --apple-id "your@email.com" \
       --team-id "XXXXXXXXXX" \
       --password "app-specific-password"
   ```

#### Signing and Notarizing

**Option 1: Automated (with credentials configured)**

Set environment variables:
```bash
export APPLE_IDENTITY="Developer ID Application: Your Name (XXXXXXXXXX)"
export APPLE_KEYCHAIN_PROFILE="MolVis-Notarize"
```

Then run the release script:
```bash
./scripts/release/release_mac.sh --build-only
```

**Option 2: Manual signing**

```bash
# Sign the app with hardened runtime
codesign --force --deep --options runtime \
    --sign "Developer ID Application: Your Name (XXXXXXXXXX)" \
    build_mac/bin/MolVis.app

# Create DMG
./scripts/package/package_mac.sh

# Sign the DMG
codesign --force --sign "Developer ID Application: Your Name (XXXXXXXXXX)" \
    dist/MolVis-<version>-macOS.dmg

# Submit for notarization
xcrun notarytool submit dist/MolVis-<version>-macOS.dmg \
    --keychain-profile "MolVis-Notarize" \
    --wait

# Staple the notarization ticket
xcrun stapler staple dist/MolVis-<version>-macOS.dmg
```

#### Verifying Signatures

```bash
# Check code signature
codesign -dv --verbose=4 build_mac/bin/MolVis.app

# Verify signature integrity
codesign --verify --deep --strict build_mac/bin/MolVis.app

# Check Gatekeeper acceptance (should show "accepted")
spctl --assess --type execute --verbose build_mac/bin/MolVis.app

# Check notarization status
xcrun stapler validate dist/MolVis-<version>-macOS.dmg
```

#### Without Developer ID (Current State)

Users downloading the unsigned DMG must bypass Gatekeeper:

**Method 1: Right-click to open**
1. Right-click (or Control-click) on MolVis.app
2. Select "Open" from the context menu
3. Click "Open" in the dialog

**Method 2: Remove quarantine attribute**
```bash
xattr -cr /Applications/MolVis.app
```

**Method 3: System Preferences**
1. Try to open the app (it will be blocked)
2. Go to System Preferences → Security & Privacy → General
3. Click "Open Anyway" next to the MolVis message

### Deploying on macOS

**To GitHub Releases:**
```bash
# Full release (build + package + upload)
./scripts/release/release_mac.sh

# Build and package only
./scripts/release/release_mac.sh --build-only

# Upload existing package
./scripts/release/upload_mac.sh 0.4.0 dist/MolVis-0.4.0-macOS.dmg
```

---

## Windows

### Windows Prerequisites

1. **Visual Studio 2019, 2022, or 2025**
   - Install with "Desktop development with C++" workload
   - Download from: https://visualstudio.microsoft.com/

2. **CUDA Toolkit 12.0+**
   - Download from: https://developer.nvidia.com/cuda-downloads
   - Ensure `CUDA_PATH` environment variable is set

3. **NVIDIA GPU** with Compute Capability 5.0+
   - Required for CUDA acceleration

### Building on Windows

**Quick build (from regular PowerShell or Command Prompt):**
```powershell
.\build.bat
```

When Visual Studio Insiders includes both preview and stable MSVC toolsets, `build.bat` selects the installed 14.44 toolset for CUDA 13.1 compatibility. It also normalizes developer-environment paths, so it is safe to launch from an IDE terminal that already initialized Visual Studio.

**From Developer Command Prompt:**
```powershell
nmake
```

**With environment setup:**
```powershell
# Load Visual Studio environment (dot-source the script)
. .\scripts\setup_env.ps1

# Build
nmake
```

**Clean build:**
```powershell
.\build.bat clean
# or
nmake clean
```

The built executable will be at: `molvis.exe`

### Running on Windows

```powershell
.\molvis.exe
# or
.\run_windows.ps1
```

### Packaging on Windows

Create a distributable ZIP package:

```powershell
.\scripts\package\package_windows.ps1 -Version "0.4.0"

# Build and package in one step
.\scripts\package\package_windows.ps1 -Version "0.4.0" -Build
```

The package will be created at: `release/MolVis-<version>-Windows-x64.zip`

### Deploying on Windows

**To GitHub Releases:**
```powershell
# Full release (build + package + upload)
.\scripts\release\release_windows.ps1 -Version "0.4.0"

# Create as draft
.\scripts\release\release_windows.ps1 -Version "0.4.0" -Draft

# Package only (no upload)
.\scripts\release\release_windows.ps1 -Version "0.4.0" -SkipUpload

# Upload existing package
.\scripts\release\upload_windows.ps1 -Version "0.4.0"
```

---

## CI/CD with GitHub Actions

The repository includes GitHub Actions workflows for automated releases. See `.github/workflows/release.yml`.

Releases are triggered by:
- Pushing a version tag: `git tag v0.4.0 && git push origin v0.4.0`
- Manual workflow dispatch

---

## Troubleshooting

## Database Tests

The molecule registry and external JSON loader have host-only tests that do not require a GPU at runtime. With a configured CMake build, run:

```bash
cmake --build . --target MolVisDatabaseTests
ctest --output-on-failure
```

The tests validate all compiled presets, stable IDs, atom and bond bounds, the semiconductor material models, and representative valid and invalid JSON input.

External molecule authoring is documented in `docs/MOLECULE_FORMAT.md`.

### macOS

**"MolVis.app is damaged and can't be opened"**
- The app was quarantined by macOS. Run: `xattr -cr MolVis.app`

**"Metal compiler not found"**
- Install the Metal toolchain:
  ```bash
  xcodebuild -downloadComponent MetalToolchain
  xcodebuild -runFirstLaunch
  ```

**SDL2 not found**
- Ensure SDL2 is installed: `brew install sdl2`
- Check pkg-config: `pkg-config --libs sdl2`

### Windows

**"CUDA not found" or build errors**
- Verify CUDA is installed: `nvcc --version`
- Check CUDA_PATH: `echo %CUDA_PATH%`

**"nmake not found"**
- Run from Visual Studio Developer Command Prompt, or
- Source the environment: `. .\scripts\setup_env.ps1`

**Missing cudart64_*.dll at runtime**
- Install CUDA Toolkit, or
- Copy the DLL from `%CUDA_PATH%\bin` to the executable directory

---

## Script Reference

| Script | Platform | Description |
|--------|----------|-------------|
| `build_mac.sh` | macOS | Build the application |
| `build.bat` | Windows | Build the application |
| `scripts/setup_env.ps1` | Windows | Set up MSVC environment |
| `scripts/package/package_mac.sh` | macOS | Create DMG package |
| `scripts/package/package_windows.ps1` | Windows | Create ZIP package |
| `scripts/release/release_mac.sh` | macOS | Build, package, and release |
| `scripts/release/release_windows.ps1` | Windows | Build, package, and release |
| `scripts/release/upload_mac.sh` | macOS | Upload to existing release |
| `scripts/release/upload_windows.ps1` | Windows | Upload to existing release |
| `scripts/generate_macos_icon.sh` | macOS | Generate .icns from PNG |
| `scripts/bundle_sdl2.sh` | macOS | Bundle SDL2 into app |
