
# MolVis Release Notes

## Unreleased

### New Features
- Added 6 new presets to the Sugars & Carbohydrates category: Mannose, Xylose, Trehalose, Cellobiose, Sorbitol, and Glucosamine
- Added a hybrid molecule database: compiled presets remain available while validated JSON files can add or explicitly override records after compilation
- Added reload, user-folder, provenance, and load-error controls to both platform UIs
- Added a Materials & Semiconductors category with finite crystal models for MoS2, WS2, WSe2, silicon, SiGe, 3C-SiC, and wurtzite GaN
- Added support for molybdenum, tungsten, selenium, germanium, and gallium across the database, periodic table, CUDA renderer, and Metal shader

### Developer Improvements
- Added stable molecule IDs, host-only database tests, schema documentation, and an example external molecule
- Saved the database modernization design in `docs/MOLECULE_DATABASE_PLAN.md`
- Fixed Windows builds launched from IDE terminals by selecting the CUDA-compatible MSVC 14.44 toolset and preventing duplicated developer-environment paths from breaking `nvcc`

## v0.4.1 (2026-01-26)

### New Features
- Added **Cyanoacrylate** (Super Glue) molecule to the Plastics category

### Improvements
- **Hide Hydrogen** option now only appears for molecules with more than 2 hydrogen atoms

---

## v0.4.0 (2026-01-25)

### Project Improvements
- **Reorganized scripts directory** - Packaging and release scripts moved to organized subdirectories
  - `scripts/package/` - DMG and ZIP packaging scripts
  - `scripts/release/` - GitHub release and upload scripts
- **Removed unused directories** - Cleaned up empty `platform/windows/` and `include/` directories
- **Added comprehensive build documentation** - New `docs/BUILDING.md` with complete instructions for:
  - Building on macOS and Windows
  - Code signing and notarization for macOS
  - Packaging and deployment workflows
- **Consolidated version numbers** - Aligned version across CMakeLists.txt, Info.plist, and package scripts

---

## v0.3.8 (2026-01-21)

### New Features
- Periodic Table panel in the UI showing a simple grid of supported elements with CPK colors
- Elements present in the currently selected molecule are highlighted for quick reference

### Notes
- The table displays a subset of elements used by MolVis (H, C, N, O, F, Na, Al, Si, P, S, Cl, B, Br, I, Ti, Fe, Cu, Pt) placed at their standard period/group positions

## v0.3.7 (2026-01-21)

### New Features
- **Added Missing Dietary Sugars** - Lactose, Maltose, and Galactose molecule presets
  - **Lactose** (C12H22O11) - Milk sugar disaccharide (glucose + galactose)
  - **Maltose** (C12H22O11) - Malt sugar disaccharide (two glucose units)
  - **Galactose** (C6H12O6) - Milk monosaccharide (glucose epimer)

### Educational Content
- Added detailed descriptions explaining:
  - Lactose's role in milk and lactose intolerance mechanisms
  - Maltose production from grain starch and its use in brewing
  - Galactose's importance for brain development and myelin formation
- All three molecules added to Sugar & Carbohydrate category for easy discovery

### Bug Fixes
- **Fixed Lactose/Galactose Stereochemistry** - Lactose now correctly displays galactose ring with proper C3/C4 stereochemistry, distinguishing it visually from Maltose

---

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
