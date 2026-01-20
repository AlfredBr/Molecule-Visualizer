# MolVis Release Notes

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