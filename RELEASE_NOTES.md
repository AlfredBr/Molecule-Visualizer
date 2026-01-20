
# MolVis Release Notes

## v0.3.6 (2026-01-20)

v0.3.6 release.

## v0.3.5 - macOS Distribution Fix (2026-01)

### Bug Fixes
- **Fixed SDL2 dependency** - SDL2 library is now bundled inside the app
- **No Homebrew required** - App now runs on any Mac without needing SDL2 installed
- **Fixed "app is damaged" workaround** - Added clear installation instructions

### Changes
- Added `scripts/bundle_sdl2.sh` to bundle SDL2 into the app bundle
- Updated build process to automatically bundle dependencies
- DMG README now includes Gatekeeper bypass instructions

---

## v0.3.2 - Application Icon Support (2026-01)

### New Features
- **Custom App Icon** - Support for custom application icon from `molvis-icon.png`
- **macOS Icon Generation** - Automatic `.icns` generation from PNG for app bundle and DMG
- **Icon Generation Script** - Added `scripts/generate_macos_icon.sh` for creating macOS icons

### Changes
- Build scripts now auto-generate macOS icon if `molvis-icon.png` exists
- CMake updated to copy icon to app bundle Resources
- Packaging script updated to use the icon for DMG volume

---

## v0.3.1 - Release Automation (2026-01)

### Changes
- **Removed binary from repo** - molvis.exe is now distributed via GitHub Releases only
- **GitHub Actions CI/CD** - Automated Windows builds and releases on tag push
- **Release packaging script** - Added `package_windows.ps1` for local release builds

---

## v0.3.0 - macOS Support (2026-01)

### New Features
- **macOS Port** - Full support for macOS using Metal and SDL2
- **Cross-Platform Build** - CMake-based build system supporting both Windows and macOS
- **Apple Silicon Optimized** - Metal compute shaders optimized for M1/M2/M3/M4 chips

### macOS Implementation
- Metal compute shaders for GPU-accelerated ray tracing (`MoleculeShaders.metal`)
- SDL2 for cross-platform windowing
- Dear ImGui with Metal + SDL2 backends
- Proper macOS app bundle with Info.plist

### Technical Details
- Platform-agnostic molecule types in `molecule_types.h`
- Abstract renderer interface in `renderer_interface.h`
- Metal renderer implementation in `metal_renderer.mm`
- macOS entry point in `main_mac.mm`

### Build Changes
- Added `CMakeLists.txt` for cross-platform CMake builds
- Added `build_mac.sh` for easy macOS building
- Added `platform/macos/Info.plist` for app bundle configuration
- Updated Dear ImGui to latest version with all backends

---

## v0.2.0 - Project Reorganization (2025-01)

### New Features
- **Dear ImGui Integration** - Full GUI with viewport, molecule selector, and controls
- **CUDA + DirectX 11 Rendering** - GPU-accelerated molecular visualization in ImGui viewport
- **Molecule Database** - 10 built-in presets (Water, Methane, Ethanol, Benzene, Caffeine, Aspirin, Dopamine, Glucose, ATP, Random)
- **Interactive Controls** - Mouse drag to rotate, scroll to zoom, auto-rotate toggle

### Project Structure
- Reorganized into proper directory structure:
  - `src/` - Main application source code
  - `src/renderer/` - CUDA rendering engine
  - `src/molecule/` - Molecule data structures and presets
  - `third_party/imgui/` - Dear ImGui library
  - `legacy/` - Original CUDA-only implementation
  - `build/` - Build artifacts
- Updated Makefile with cleaner organization
- Enhanced build.bat with better output and auto-kill of running instances

### Technical Details
- Uses staging buffer approach for CUDA→DX11 texture transfer
- Ray-sphere intersection for atoms, ray-cylinder for bonds
- Phong shading with specular highlights
- CPK color convention for atoms

---

## v0.1.0 - Initial Release

### Features
- CUDA-based molecular rendering
- Ball-and-stick visualization
- CPK coloring convention
- Win32 window management
- Keyboard controls for rotation and molecule selection